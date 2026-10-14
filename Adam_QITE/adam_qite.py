#!/usr/bin/env python
"""An Adam-like step-size controller for variational QITE, built on Route 3's indicators.

Extends `../Attack_Route_03/adaptive_qite.py` (Route 3 of the internal memo
"Adaptive Time-Stepping for Variational QITE", `../Research_Papers/QITE_timestep.pdf`).
Route 3 as implemented there drives delta_tau with a hand-tuned multiplicative ratchet
(GROW=1.25 on a clean accept, SHRINK=0.5 on rejection) -- Assumption 1 of that script
flags both numbers as tuned by eye,
and `../Audits/Audit_Project_Wide/R_Future_Routes.md` (R3) names replacing them as a priority.

The starting observation here: a local truncation error estimate is ALREADY being computed and
thrown away. `energy_check()` forms dE_actual - dE_pred with dE_pred = -(Y.theta_dot)*delta_tau
the exact first-order model of the step, so the difference is O(delta_tau^2) -- an error proxy
at ZERO extra circuits, inside Route 3's "cheap indicators, no extra stages" premise. Today it
only feeds a binary in-band/out-of-band retry test. `../Audits/Audit_Project_Wide/R_Future_Routes.md`
states the project "has no local truncation error estimate at all" and routes R1
(step-doubling, 3 f-evals
per step) through that gap; this gets one for free instead.

Adam is then not decoration. Its ratio m_hat/sqrt(v_hat) is a signal-to-noise ratio on that
estimator, so growth is throttled when the error signal cannot be told apart from its own
fluctuation -- memo Eq. (9)'s question ("when does the estimator stop carrying signal") answered
inside Route 3, with no shot-noise model (R5) and no extra stages (R1/R2).

SCOPE: the Adam moments drive the SCALAR delta_tau only. The update stays
theta <- theta + delta_tau * M^+ Y, so this is still a genuine McLachlan/QITE trajectory and tau
keeps its physical meaning. Per-parameter Adam on the flow vector (a diagonal preconditioner
overwriting M^+) is deliberately NOT implemented -- it would stop being imaginary-time evolution.

Three scenarios, one function each, each reusing the seed of the Route-3 target it is compared
against so theta_0 is identical and the difference is the controller alone:

  E1  target_e1_controller_comparison()  Adam vs Route 3 vs fixed step, noiseless   (SEED_A)
  E2  target_e2_noisy_controller()       Adam vs Route 3 under depolarizing noise   (SEED_D)
  E3  target_e3_stiffness_9param()       Adam vs Route 3 on the rank-6-of-9 ansatz  (SEED_C)

Run with the physics environment:

    ~/micromamba/envs/qml-julia/bin/python adam_qite.py --target all

See `../Notes_on_Project/Adam_Like_Step_Controller_For_Variational_QITE.md` for the prose.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "Attack_Route_03"))

# Imported, not copied: Attack_Route_03 itself was copy-pasted out of the QDP replication and
# has already had to be kept in sync by hand once. adaptive_qite.py has an
# `if __name__ == "__main__"` guard and does only constant construction at module level, so
# importing it is side-effect-free (it also sets matplotlib's Agg backend, which we want).
from adaptive_qite import (  # noqa: E402
    COLORS,
    DE_BAND_ABS,
    DE_BAND_REL,
    DELTA_TAU0,
    DELTA_TAU_MAX,
    DELTA_TAU_MIN,
    GROW,
    MAX_RETRIES,
    MAX_STEPS,
    RCOND_DEFAULT,
    RESIDUAL_REL_TOL,
    SEED_A,
    SEED_C,
    SEED_D,
    SHRINK,
    STIFFNESS_CAP_COEFF,
    TAU_MAX,
    adaptive_qite_run,
    banner,
    energy_and_variance,
    energy_check,
    exact_ground_energy,
    fixed_step_qite_run,
    h2_hamiltonian,
    make_noisy_rho_fn,
    make_rho_fn,
    measure_f,
    route3_indicators,
    trajectory_deviation,
    variational_floor,
)
import matplotlib.pyplot as plt  # noqa: E402  (adaptive_qite already selected the Agg backend)

FIGURES = HERE / "figures"

# -- Adam controller constants ---------------------------------------------------------
# NOT tuned here. BETA1/BETA2/EPS are stock Adam defaults; K_I is Gustafsson's 0.3/(p+1) at
# p=1 (explicit Euler), the value R_Future_Routes.md R3 already proposes. Every other knob
# (GROW, SHRINK, DELTA_TAU*, DE_BAND_*, RESIDUAL_REL_TOL, RCOND_DEFAULT, STIFFNESS_CAP_COEFF,
# MAX_STEPS, MAX_RETRIES, TAU_MAX) is imported from adaptive_qite.py and never redefined, so
# the two controllers differ only in how delta_tau evolves between accepted steps.
BETA1 = 0.9
BETA2 = 0.999
EPS = 1e-8
K_I = 0.15

ASSUMPTIONS: list[str] = [
    "Assumptions 1-8 of ../Attack_Route_03/adaptive_qite.py are inherited unchanged "
    "(finite-difference M and Y at h=1e-6 rather than Hadamard tests; no shot noise, "
    "so 'noisy' means biased-but-"
    "deterministic; the F8/M5 policy of committing a trial step once max_retries or "
    "delta_tau_min is exhausted, with accepted=False recorded; Eq. (7) debiasing applied to M "
    "and Y before the solve with the M1-corrected Y factor q^2+q(1-q)/2). This script changes "
    "the step-size law and nothing else.",

    "Reading |dE_actual - dE_pred| as an O(delta_tau^2) local truncation error proxy is OUR "
    "claim, not the memo's -- the memo never forms an error estimate for Route 3 at all (that "
    "is why R_Future_Routes.md R1 proposes step-doubling to get one). It follows from dE_pred "
    "= -(Y.theta_dot)*delta_tau being the exact first-order model of the Euler step being "
    "taken, so the leading discrepancy is the step's own second-order term. It is a proxy in "
    "the energy, not in theta: it measures the truncation error along the one direction the "
    "controller can actually observe cheaply.",

    "e_hat is normalized by the SAME band energy_check() already uses to reject a step, "
    "band = max(DE_BAND_ABS, DE_BAND_REL*|dE_pred|), so e_hat > 1 is exactly the existing "
    "reject condition. No new tolerance is introduced: the controller reuses the gate's own "
    "tolerance as its setpoint, and the gate itself is left in place unchanged.",

    "The control law delta_tau *= clip(m_e_hat**(-s_e*K_I), SHRINK, GROW) is our synthesis. "
    "The m_e_hat**(-K_I) part is Gustafsson's integral law (R3) with the tolerance folded into "
    "e_hat; the new part is the exponent being scaled by the Adam confidence ratio "
    "s_e = m_e_hat/(sqrt(v_e_hat)+EPS), which lies in (0,1] by Jensen. A consistent error "
    "signal (s_e -> 1) gets the full Gustafsson response; a noise-dominated one (s_e -> 0) "
    "drives the exponent to 0 and freezes delta_tau. GROW and SHRINK survive only as clip "
    "bounds (a trust region), no longer as the step law itself -- which is the specific thing "
    "Assumption 1 of the parent script complains about.",

    "Adam's first-moment EMA already supplies the smoothing a PI controller's I-term provides, "
    "so no separate proportional term k_P (e_{n-1}/e_n)^k_P is included in this version. That "
    "is a YAGNI call, not a claim that it would not help.",

    "The McLachlan residual channel gets its own EMA pair (m_r, v_r) and enters as a growth "
    "CEILING only: growth <- min(growth, (RESIDUAL_REL_TOL/m_r_hat)**(s_r*K_I)) when "
    "m_r_hat exceeds RESIDUAL_REL_TOL, which is <= 1 by construction. This preserves "
    "Assumption 2 of the parent script: r_n does not depend on delta_tau, so it can never "
    "justify rejecting or shrinking the CURRENT step, only withholding growth from the next "
    "one. The same K_I is reused for both channels rather than inventing a second exponent.",

    "The EMAs update once per ACCEPTED step, using the accepted trial's e_hat, and the bias "
    "correction counts accepted steps. A retry cascade therefore does not feed its rejected "
    "trials into the moments -- otherwise a single bad step would be counted several times "
    "and would poison v_e for the following ~1/(1-BETA2) steps.",

    "m_e_hat is floored at EPS before exponentiation, so an exactly-perfect prediction "
    "(dE_actual == dE_pred, which does happen at delta_tau_min) gives maximum growth rather "
    "than a division by zero. The clip to GROW makes the floor's exact value irrelevant.",

    "No constant in this file was adjusted after seeing a result. Where Adam loses to the "
    "Route-3 ratchet that is reported as the finding. The hyperparameter ablation "
    "(BETA1/BETA2 sensitivity, momentum-only vs RMS-only) is NOT in this pass; this fairness "
    "protocol stands in for it, and the ablation remains open work.",
]


# --------------------------------------------------------------------------------------
# The controller
# --------------------------------------------------------------------------------------

class AdamStepController:
    """Adam-style moments on Route 3's two delta_tau-relevant indicators.

    Channel e: the normalized energy local-error proxy, e_hat = |dE_actual - dE_pred|/band.
               Drives growth directly (Gustafsson, exponent gated by confidence).
    Channel r: the relative McLachlan residual ||M theta_dot - Y||/||Y||. delta_tau-
               independent, so it can only withhold growth (parent Assumption 2).
    """

    def __init__(self, beta1=BETA1, beta2=BETA2, eps=EPS, k_i=K_I, grow=GROW, shrink=SHRINK,
                 residual_rel_tol=RESIDUAL_REL_TOL):
        self.beta1, self.beta2, self.eps, self.k_i = beta1, beta2, eps, k_i
        self.grow, self.shrink = grow, shrink
        self.residual_rel_tol = residual_rel_tol
        self.m_e = self.v_e = self.m_r = self.v_r = 0.0
        self.n = 0

    def update(self, e_hat, r_hat):
        """Fold one accepted step's indicators in, return the multiplicative step factor."""
        self.n += 1
        b1, b2 = self.beta1, self.beta2

        self.m_e = b1 * self.m_e + (1.0 - b1) * e_hat
        self.v_e = b2 * self.v_e + (1.0 - b2) * e_hat ** 2
        self.m_r = b1 * self.m_r + (1.0 - b1) * r_hat
        self.v_r = b2 * self.v_r + (1.0 - b2) * r_hat ** 2

        # Adam's bias correction. b1=0 (the degeneracy check) gives c1 = 1, i.e. a no-op,
        # which is right: with no memory there is no initialization bias to correct.
        c1 = 1.0 - b1 ** self.n
        c2 = 1.0 - b2 ** self.n
        m_e_hat = self.m_e / c1 if c1 > 0 else self.m_e
        v_e_hat = self.v_e / c2 if c2 > 0 else self.v_e
        m_r_hat = self.m_r / c1 if c1 > 0 else self.m_r
        v_r_hat = self.v_r / c2 if c2 > 0 else self.v_r

        # Confidence ratios: mean over RMS, in (0, 1] by Jensen. This is the signal-to-noise
        # gate -- it is what makes the controller stop reacting to an error estimate that is
        # mostly its own fluctuation, without ever needing a shot-noise model.
        s_e = m_e_hat / (np.sqrt(v_e_hat) + self.eps)
        s_r = m_r_hat / (np.sqrt(v_r_hat) + self.eps)

        # Gustafsson's integral law with the tolerance already folded into e_hat (setpoint 1),
        # exponent scaled by confidence.
        growth = max(m_e_hat, self.eps) ** (-s_e * self.k_i)

        # Residual channel: ceiling only, never a floor.
        if m_r_hat > self.residual_rel_tol:
            growth = min(growth, (self.residual_rel_tol / m_r_hat) ** (s_r * self.k_i))

        growth = float(np.clip(growth, self.shrink, self.grow))
        return dict(growth=growth, m_e_hat=float(m_e_hat), v_e_hat=float(v_e_hat),
                    s_e=float(s_e), m_r_hat=float(m_r_hat), s_r=float(s_r))


def adam_qite_run(theta0, rho_fn, H, n_params, N, tau_max=TAU_MAX, delta_tau0=DELTA_TAU0,
                  p=0.0, debias=True, rcond=RCOND_DEFAULT,
                  stiffness_cap_coeff=STIFFNESS_CAP_COEFF,
                  delta_tau_min=DELTA_TAU_MIN, delta_tau_max=DELTA_TAU_MAX,
                  max_steps=MAX_STEPS, max_retries=MAX_RETRIES,
                  beta1=BETA1, beta2=BETA2, k_i=K_I):
    """Route 3's measurement, indicators and accept/reject gate, with the GROW/SHRINK ratchet
    replaced by AdamStepController. Returns the same log dict shape as adaptive_qite_run()
    (plus the Adam state), so trajectory_deviation() and the M6 circuit accounting apply
    unchanged.
    """
    theta = np.array(theta0, dtype=float)
    tau = 0.0
    delta_tau = delta_tau0
    E0, _ = energy_and_variance(rho_fn(theta), H)
    ctrl = AdamStepController(beta1=beta1, beta2=beta2, k_i=k_i)

    log = dict(tau=[0.0], E=[E0], delta_tau=[delta_tau0], cond_M=[np.nan],
               sigma_min=[np.nan], cond_spectral=[np.nan],
               residual=[np.nan], dE_pred=[np.nan], dE_exact_pred=[np.nan],
               dE_actual=[np.nan], accepted=[True], retries=[0], reasons=[[]],
               cap_bound=[False], s_eff=[np.nan],
               e_hat=[np.nan], m_e=[np.nan], v_e=[np.nan], s_e=[np.nan],
               m_r=[np.nan], s_r=[np.nan], growth=[np.nan],
               n_evals=0, n_circ=0)

    step = 0
    while tau < tau_max and step < max_steps:
        theta_dot_val, M, Y, rho_n = measure_f(theta, rho_fn, H, N, p, rcond, debias)
        log["n_evals"] += 1
        log["n_circ"] += 2 * N + 1
        E_n, var_n = energy_and_variance(rho_n, H)
        residual, cond_M, sigma_min, cond_spectral = route3_indicators(M, Y, theta_dot_val, rcond)
        r_hat = residual / (np.linalg.norm(Y) + 1e-12)

        # Proactive cond(M) cap, identical to adaptive_qite_run. stiffness_cap_coeff=inf
        # disables it (target E3's third arm) without an inf/inf nan when cond_M is also inf.
        if np.isinf(stiffness_cap_coeff):
            cap = np.inf
        else:
            cap = stiffness_cap_coeff / cond_M if cond_M > 0 else np.inf
        cap_bound = cap < delta_tau and cap < delta_tau_max
        delta_tau = max(delta_tau_min, min(delta_tau, cap, delta_tau_max))

        retries = 0
        reasons_seen = []
        while True:
            theta_trial = theta + delta_tau * theta_dot_val
            E_trial, _ = energy_and_variance(rho_fn(theta_trial), H)
            log["n_circ"] += 1
            dE_pred, dE_exact_pred, dE_actual, reasons = energy_check(
                E_n, var_n, Y, theta_dot_val, delta_tau, E_trial)
            reasons_seen.extend(reasons)
            accepted = not reasons
            if accepted or retries >= max_retries or delta_tau <= delta_tau_min:
                break
            delta_tau = max(delta_tau * SHRINK, delta_tau_min)
            retries += 1

        # The local error proxy, normalized by energy_check's own reject band so that
        # e_hat > 1 is exactly the condition that would have rejected this step.
        band = max(DE_BAND_ABS, DE_BAND_REL * abs(dE_pred))
        e_hat = abs(dE_actual - dE_pred) / band
        state = ctrl.update(e_hat, r_hat)

        theta = theta_trial
        tau += delta_tau
        log["tau"].append(tau)
        log["E"].append(E_trial)
        log["delta_tau"].append(delta_tau)
        log["cond_M"].append(cond_M)
        log["sigma_min"].append(sigma_min)
        log["cond_spectral"].append(cond_spectral)
        log["residual"].append(residual)
        log["dE_pred"].append(dE_pred)
        log["dE_exact_pred"].append(dE_exact_pred)
        log["dE_actual"].append(dE_actual)
        log["accepted"].append(accepted)
        log["retries"].append(retries)
        log["reasons"].append(reasons_seen)
        log["cap_bound"].append(cap_bound)
        log["s_eff"].append(1 + retries)
        log["e_hat"].append(e_hat)
        log["m_e"].append(state["m_e_hat"])
        log["v_e"].append(state["v_e_hat"])
        log["s_e"].append(state["s_e"])
        log["m_r"].append(state["m_r_hat"])
        log["s_r"].append(state["s_r"])
        log["growth"].append(state["growth"])

        delta_tau = float(np.clip(delta_tau * state["growth"], delta_tau_min, delta_tau_max))
        step += 1

    log["exit_reason"] = "tau_max" if tau >= tau_max else "max_steps"
    for key in ("tau", "E", "delta_tau", "cond_M", "sigma_min", "cond_spectral", "residual",
                "dE_pred", "dE_exact_pred", "dE_actual", "retries", "s_eff",
                "e_hat", "m_e", "v_e", "s_e", "m_r", "s_r", "growth"):
        log[key] = np.array(log[key], dtype=float)
    log["accepted"] = np.array(log["accepted"], dtype=bool)
    log["cap_bound"] = np.array(log["cap_bound"], dtype=bool)
    return log


# --------------------------------------------------------------------------------------
# Shared reporting helpers (M6 currency)
# --------------------------------------------------------------------------------------

def circuit_cumsum(run, n_params):
    """Per-step cumulative circuit count in M6's currency: 2N+1 for measure_f's central
    differences, plus one rho_fn call per trial (1 + retries for a gated run, 1 for the
    fixed-step baseline, which has no retry loop).
    """
    n_steps = len(run["tau"]) - 1
    per_eval = 2 * n_params + 1
    if "retries" in run:
        return np.cumsum(per_eval + 1 + run["retries"][1:])
    return np.cumsum(np.full(n_steps, per_eval + 1))


def energy_at_budget(run, n_params, budget):
    """Energy this run had reached once it had spent `budget` circuits."""
    cum = circuit_cumsum(run, n_params)
    i = int(np.searchsorted(cum, budget))
    return float(run["E"][min(i + 1, len(run["E"]) - 1)])


def run_summary(name, run, E_exact, n_params):
    forced = int(np.sum(~run["accepted"][1:])) if "accepted" in run else 0
    retried_ok = (int(np.sum((run["retries"][1:] > 0) & run["accepted"][1:]))
                  if "retries" in run else 0)
    dt = run["delta_tau"][1:] if "delta_tau" in run else None
    dt_str = f"[{dt.min():.4g}, {dt.max():.4g}]" if dt is not None and len(dt) else "fixed"
    print(f"  {name:<22} steps={len(run['tau']) - 1:4d}  circuits={run['n_circ']:6d}  "
          f"final E={run['E'][-1]:.6f}  error={run['E'][-1] - E_exact:+.3e}")
    print(f"  {'':<22} delta_tau range={dt_str}  retried-then-accepted={retried_ok}  "
          f"forced through unaccepted={forced}")


def adam_signal_diagnostic(run, label=""):
    """Why the confidence gate behaves as it does: the distribution of e_hat.

    Adam's second-moment normalization assumes a signal whose fluctuation is roughly
    symmetric about its mean -- true of minibatch gradient noise, which is what it was
    designed for. A truncation-error magnitude is not that: it is non-negative and spiky, so
    sqrt(v_hat) is set by rare large values while m_hat sits near the median. The reported
    mean/median ratio quantifies that skew, and s_e mean shows what it costs the controller.
    """
    e = run["e_hat"][1:]
    g = run["growth"][1:]
    se = run["s_e"][1:]
    finite = e[np.isfinite(e)]
    print(f"  Adam signal diagnostic{label}:")
    print(f"    e_hat   median={np.median(finite):.3e}  mean={finite.mean():.3e}  "
          f"p90={np.percentile(finite, 90):.3e}  p99={np.percentile(finite, 99):.3e}  "
          f"max={finite.max():.3e}")
    print(f"            mean/median={finite.mean() / np.median(finite):.1f} "
          f"(1 = symmetric, >>1 = heavy right tail)   std/mean={finite.std() / finite.mean():.2f}")
    print(f"    s_e     mean={np.nanmean(se):.4f}  ->  effective exponent s_e*K_I = "
          f"{np.nanmean(se) * K_I:.5f}  vs ungated K_I = {K_I}")
    ungated = np.clip(np.maximum(run["m_e"][1:], EPS) ** (-K_I), SHRINK, GROW)
    print(f"    growth  mean={np.nanmean(g):.5f}  frozen (<1.001) on "
          f"{np.mean(g < 1.001):.0%} of steps   (ungated law would average "
          f"{np.nanmean(ungated):.4f})")


# --------------------------------------------------------------------------------------
# Target E1 -- Adam vs the Route-3 ratchet vs a fixed step, noiseless
# --------------------------------------------------------------------------------------

def target_e1_controller_comparison(plots=True):
    banner("E1  Adam-like controller vs Route-3 ratchet vs fixed step (noiseless, 5-param)")
    print("  Same theta_0 (SEED_A), same ansatz, same gate, same cond(M) cap as Target A")
    print("  of ../Attack_Route_03/adaptive_qite.py -- the ONLY difference is how delta_tau")
    print("  evolves between accepted steps. The Route-3 arm below must reproduce that")
    print("  script's Target A exactly (172 steps, final E = -1.122993, 2069 circuits) or")
    print("  the comparison is not measuring what it claims to.\n")

    H = h2_hamiltonian()
    E_exact = exact_ground_energy(H)
    n_params, N = 5, 5
    theta0 = np.random.default_rng(SEED_A).uniform(-0.5, 0.5, size=n_params)
    rho_fn = make_rho_fn(n_params=n_params)

    adam = adam_qite_run(theta0, rho_fn, H, n_params, N)
    route3 = adaptive_qite_run(theta0, rho_fn, H, n_params, N)
    fixed = fixed_step_qite_run(theta0, rho_fn, H, n_params, N, delta_tau=DELTA_TAU0)

    print(f"  Exact H2 ground energy (eigh): {E_exact:.6f}\n")
    run_summary("Adam controller", adam, E_exact, n_params)
    run_summary("Route 3 (ratchet)", route3, E_exact, n_params)
    print(f"  {'fixed step':<22} steps={len(fixed['tau']) - 1:4d}  "
          f"circuits={fixed['n_circ']:6d}  final E={fixed['E'][-1]:.6f}  "
          f"error={fixed['E'][-1] - E_exact:+.3e}")

    # F9/M6: the honest comparison is at matched CIRCUIT count, not matched f-evals.
    budget = min(int(circuit_cumsum(adam, n_params)[-1]),
                 int(circuit_cumsum(route3, n_params)[-1]),
                 int(circuit_cumsum(fixed, n_params)[-1]))
    print(f"\n  [M6] At matched circuit budget ({budget} circuits):")
    print(f"    Adam    E = {energy_at_budget(adam, n_params, budget):.6f}")
    print(f"    Route 3 E = {energy_at_budget(route3, n_params, budget):.6f}")
    print(f"    fixed   E = {energy_at_budget(fixed, n_params, budget):.6f}")

    floor = variational_floor(rho_fn, n_params, H)
    print(f"\n  Variational floor of this ansatz (M6): E_floor = {floor:.6f}")
    print(f"    Adam    final E - floor = {adam['E'][-1] - floor:+.3e}")
    print(f"    Route 3 final E - floor = {route3['E'][-1] - floor:+.3e}")
    print(f"    fixed   final E - floor = {fixed['E'][-1] - floor:+.3e}")

    print(f"\n  Adam state over the run: e_hat median={np.nanmedian(adam['e_hat']):.4g}, "
          f"s_e range=[{np.nanmin(adam['s_e']):.4g}, {np.nanmax(adam['s_e']):.4g}], "
          f"growth range=[{np.nanmin(adam['growth']):.4g}, {np.nanmax(adam['growth']):.4g}]")
    n_res_gated = int(np.sum(adam["m_r"][1:] > RESIDUAL_REL_TOL))
    print(f"    McLachlan channel withheld growth on {n_res_gated} of "
          f"{len(adam['tau']) - 1} steps (m_r_hat > {RESIDUAL_REL_TOL})")

    # Saturation diagnostic. growth = m_e_hat**(-s_e*K_I) is >1 whenever the SMOOTHED
    # normalized error sits below the gate band (m_e_hat < 1), so the law can only shrink on
    # steps the gate would have rejected anyway. If e_hat lives orders of magnitude below its
    # setpoint of 1, the controller is pinned at the GROW clip and delta_tau is really being
    # set by delta_tau_max / the cond(M) cap / the retry loop -- not by Adam. Reported, not
    # tuned away: DE_BAND_ABS/REL are imported from the parent script and fixing them here
    # after seeing this number is exactly the fudge the project's rules forbid.
    n_steps = len(adam["tau"]) - 1
    n_over = int(np.sum(adam["e_hat"][1:] > 1.0))
    n_shrink = int(np.sum(adam["growth"][1:] < 1.0))
    n_pinned = int(np.sum(adam["growth"][1:] >= GROW - 1e-12))
    print(f"    setpoint diagnostic: e_hat > 1 (gate band) on {n_over}/{n_steps} steps; "
          f"growth < 1 on {n_shrink}/{n_steps}; growth pinned at GROW={GROW} on "
          f"{n_pinned}/{n_steps}")
    adam_signal_diagnostic(adam, " (E1, SEED_A)")

    if plots:
        fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2))
        ax = axes[0]
        ax.plot(adam["tau"], adam["E"] - E_exact, "-", color=COLORS[2], lw=1.8, label="Adam")
        ax.plot(route3["tau"], route3["E"] - E_exact, "-", color=COLORS[0], lw=1.6,
                label="Route 3 (ratchet)")
        ax.plot(fixed["tau"], fixed["E"] - E_exact, "--", color=COLORS[1], lw=1.4,
                label=f"fixed ($\\Delta\\tau$={DELTA_TAU0})")
        ax.axhline(floor - E_exact, color="gray", ls=":", lw=1.3, label="variational floor")
        ax.set_yscale("log")
        ax.set_xlabel(r"imaginary time $\tau$")
        ax.set_ylabel(r"$E(\tau) - E_{exact}$")
        ax.set_title("E1: energy convergence")
        ax.legend(frameon=True, fontsize=8)
        ax.grid(alpha=0.3)

        ax = axes[1]
        ax.plot(adam["tau"][1:], adam["delta_tau"][1:], "-", color=COLORS[2], lw=1.5,
                label="Adam")
        ax.plot(route3["tau"][1:], route3["delta_tau"][1:], "-", color=COLORS[0], lw=1.3,
                alpha=0.8, label="Route 3")
        ax.set_yscale("log")
        ax.set_xlabel(r"imaginary time $\tau$")
        ax.set_ylabel(r"$\Delta\tau$")
        ax.set_title("E1: step-size trajectory")
        ax.legend(frameon=True, fontsize=8)
        ax.grid(alpha=0.3)

        ax = axes[2]
        ax.plot(adam["tau"][1:], adam["s_e"][1:], "-", color=COLORS[3], lw=1.5,
                label=r"$s_e$ (energy channel)")
        ax.plot(adam["tau"][1:], adam["s_r"][1:], "-", color=COLORS[4], lw=1.3,
                label=r"$s_r$ (residual channel)")
        ax.set_ylim(0, 1.05)
        ax.set_xlabel(r"imaginary time $\tau$")
        ax.set_ylabel(r"confidence $\hat m/\sqrt{\hat v}$")
        ax.set_title("E1: Adam signal-to-noise ratios")
        ax.legend(frameon=True, fontsize=8)
        ax.grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(FIGURES / "target_e1_controllers.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'target_e1_controllers.png'}")

    return dict(adam=adam, route3=route3, fixed=fixed, E_exact=E_exact, floor=floor)


# --------------------------------------------------------------------------------------
# Target E2 -- the load-bearing test: does EMA smoothing help under noisy dE readings?
# --------------------------------------------------------------------------------------

def target_e2_noisy_controller(plots=True):
    banner("E2  Adam vs Route 3 under global depolarizing noise (5-param, debiased)")
    print("  This is where the EMA either earns its place or does not. The Route-3 ratchet")
    print("  reacts to every individual noisy dE (one bad reading halves delta_tau); the Adam")
    print("  controller sees the same readings through m_e and additionally discounts them by")
    print("  s_e when they stop being self-consistent. Each controller is scored against ITS")
    print("  OWN p=0 run, so this measures noise robustness, not the p=0 gap already in E1.\n")

    H = h2_hamiltonian()
    n_params, N = 5, 5
    theta0 = np.random.default_rng(SEED_D).uniform(-0.5, 0.5, size=n_params)
    rho_fn_clean = make_rho_fn(n_params=n_params)

    base_adam = adam_qite_run(theta0, rho_fn_clean, H, n_params, N)
    base_route3 = adaptive_qite_run(theta0, rho_fn_clean, H, n_params, N)

    # The p=0 references themselves, reported because a stalled baseline changes what the
    # RMS column below even means (trajectory_deviation compares on common support only).
    for label, base in (("Adam", base_adam), ("Route 3", base_route3)):
        print(f"  p=0 reference, {label:<8} {len(base['tau']) - 1:4d} steps, "
              f"final tau = {base['tau'][-1]:.4f}, final E = {base['E'][-1]:.6f}, "
              f"exit = {base['exit_reason']}")
    adam_signal_diagnostic(base_adam, " (E2 p=0 baseline, SEED_D)")
    print()

    p_grid = (1e-4, 1e-3, 1e-2, 5e-2)
    rms_adam, rms_route3, se_mean = [], [], []
    runs_adam = {}
    print(f"  {'p':>8} {'RMS Adam':>13} {'cov%':>6} {'RMS Route3':>13} {'cov%':>6} "
          f"{'mean s_e':>9} {'steps A/R':>11}")
    for p in p_grid:
        rho_fn_noisy = make_noisy_rho_fn(p, n_params=n_params)
        run_a = adam_qite_run(theta0, rho_fn_noisy, H, n_params, N, p=p, debias=True)
        run_r = adaptive_qite_run(theta0, rho_fn_noisy, H, n_params, N, p=p, debias=True)
        dev_a = trajectory_deviation(run_a, base_adam)
        dev_r = trajectory_deviation(run_r, base_route3)
        rms_adam.append(dev_a["rms_common"])
        rms_route3.append(dev_r["rms_common"])
        se_mean.append(float(np.nanmean(run_a["s_e"][1:])))
        runs_adam[p] = run_a
        flags = (" STALLED(A)" if dev_a["stalled"] else "") + (" STALLED(R)" if dev_r["stalled"] else "")
        print(f"  {p:8.4f} {dev_a['rms_common']:13.3e} {dev_a['coverage'] * 100:6.1f} "
              f"{dev_r['rms_common']:13.3e} {dev_r['coverage'] * 100:6.1f} "
              f"{se_mean[-1]:9.4f} "
              f"{len(run_a['tau']) - 1:5d}/{len(run_r['tau']) - 1:<5d}{flags}")

    print("\n  Ratio RMS_Adam / RMS_Route3 (below 1 means Adam tracks its own clean run more")
    print("  closely under the same noise):")
    for p, ra, rr in zip(p_grid, rms_adam, rms_route3):
        print(f"    p={p:<8.4f} ratio = {ra / rr if rr > 0 else np.nan:.4f}")

    if plots:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
        ax = axes[0]
        ax.loglog(p_grid, rms_adam, "-o", color=COLORS[2], ms=4, lw=1.6, label="Adam")
        ax.loglog(p_grid, rms_route3, "-s", color=COLORS[0], ms=4, lw=1.6, label="Route 3")
        ax.set_xlabel("depolarizing probability $p$")
        ax.set_ylabel("RMS energy deviation from own $p=0$ run")
        ax.set_title("E2: noise robustness")
        ax.legend(frameon=True)
        ax.grid(alpha=0.3, which="both")

        ax = axes[1]
        for p, color in zip(p_grid, (COLORS[0], COLORS[1], COLORS[2], COLORS[3])):
            run = runs_adam[p]
            ax.plot(run["tau"][1:], run["s_e"][1:], "-", color=color, lw=1.3,
                    label=f"p={p:g}")
        ax.set_ylim(0, 1.05)
        ax.set_xlabel(r"imaginary time $\tau$")
        ax.set_ylabel(r"$s_e = \hat m_e/\sqrt{\hat v_e}$")
        ax.set_title("E2: error-signal confidence vs. noise")
        ax.legend(frameon=True, fontsize=8)
        ax.grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(FIGURES / "target_e2_noisy.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'target_e2_noisy.png'}")

    return dict(p_grid=p_grid, rms_adam=rms_adam, rms_route3=rms_route3)


# --------------------------------------------------------------------------------------
# Target E3 -- the rank-6-of-9 ansatz: can the Adam law replace the hand-tuned cap?
# --------------------------------------------------------------------------------------

def target_e3_stiffness_9param(plots=True):
    banner("E3  Adam vs Route 3 on the rank-6-of-9 (9-param) ansatz")
    print("  Target C of ../Attack_Route_03/adaptive_qite.py handles this ansatz's rank")
    print("  deficiency with a PROACTIVE cap delta_tau <= STIFFNESS_CAP_COEFF/cond(M),")
    print("  coefficient 2000 picked by hand. The third arm below disables that cap")
    print("  entirely and asks whether the variance-normalized law absorbs the stiffness")
    print("  on its own -- with a Route-3 cap-off control arm, so a difference can be")
    print("  attributed to the law or the cap. Reference (Target C):")
    print("  rcond=1e-8 -> error 8.868e-04 in 108 steps; rcond=1e-2 -> 2.143e-03 in 17.\n")

    H = h2_hamiltonian()
    E_exact = exact_ground_energy(H)
    n_params, N = 9, 9
    theta0 = np.random.default_rng(SEED_C).uniform(-0.5, 0.5, size=n_params)
    rho_fn = make_rho_fn(n_params=n_params)

    results = {}
    print(f"  {'rcond':>8} {'arm':<20} {'final E':>12} {'error':>12} "
          f"{'cond(M) range':>24} {'delta_tau range':>22} {'steps':>6}")
    for rcond in (1e-8, 1e-2):
        # Four arms, not three: "Adam cap-off vs Route-3 cap-on" would change the step law
        # AND the cap at the same time, so a win could not be attributed to either. The
        # Route-3 cap-off arm is the control that separates them.
        arms = {
            "Route 3 (cap on)": adaptive_qite_run(theta0, rho_fn, H, n_params, N, rcond=rcond),
            "Route 3 (cap off)": adaptive_qite_run(theta0, rho_fn, H, n_params, N, rcond=rcond,
                                                   stiffness_cap_coeff=np.inf),
            "Adam (cap on)": adam_qite_run(theta0, rho_fn, H, n_params, N, rcond=rcond),
            "Adam (cap off)": adam_qite_run(theta0, rho_fn, H, n_params, N, rcond=rcond,
                                            stiffness_cap_coeff=np.inf),
        }
        results[rcond] = arms
        for arm_name, run in arms.items():
            cond_vals = run["cond_M"][1:]
            finite = cond_vals[np.isfinite(cond_vals)]
            cond_range = f"[{finite.min():.3g}, {finite.max():.3g}]" if len(finite) else "all inf"
            dt = run["delta_tau"][1:]
            dt_range = f"[{dt.min():.4g}, {dt.max():.4g}]"
            print(f"  {rcond:8.1e} {arm_name:<20} {run['E'][-1]:12.6f} "
                  f"{run['E'][-1] - E_exact:+12.3e} {cond_range:>24} "
                  f"{dt_range:>22} {len(run['tau']) - 1:6d}")

    if plots:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
        styles = {"Route 3 (cap on)": (COLORS[0], "-"), "Route 3 (cap off)": (COLORS[1], "--"),
                  "Adam (cap on)": (COLORS[2], "-"), "Adam (cap off)": (COLORS[3], "--")}
        for arm_name, run in results[1e-8].items():
            color, ls = styles[arm_name]
            axes[0].plot(run["tau"][1:], np.clip(run["cond_M"][1:], 0, 1e8), ls, color=color,
                         lw=1.5, label=arm_name)
            axes[1].plot(run["tau"][1:], run["delta_tau"][1:], ls, color=color, lw=1.5,
                         label=arm_name)
        axes[0].set_yscale("log")
        axes[0].set_xlabel(r"imaginary time $\tau$")
        axes[0].set_ylabel(r"$\mathrm{cond}_F(M)$ (clipped at $10^8$)")
        axes[0].set_title("E3: stiffness, 9-param ansatz (rcond=1e-8)")
        axes[0].legend(frameon=True, fontsize=8)
        axes[0].grid(alpha=0.3)
        axes[1].set_yscale("log")
        axes[1].set_xlabel(r"imaginary time $\tau$")
        axes[1].set_ylabel(r"$\Delta\tau$")
        axes[1].set_title("E3: step size with and without the cap")
        axes[1].legend(frameon=True, fontsize=8)
        axes[1].grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(FIGURES / "target_e3_stiffness.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'target_e3_stiffness.png'}")

    return dict(results=results, E_exact=E_exact)


# --------------------------------------------------------------------------------------

def print_summary():
    banner("Assumptions (ours, on top of ../Attack_Route_03/adaptive_qite.py's)")
    for i, a in enumerate(ASSUMPTIONS, 1):
        print(f"  {i}. {a}\n")


TARGETS = {
    "e1": target_e1_controller_comparison,
    "e2": target_e2_noisy_controller,
    "e3": target_e3_stiffness_9param,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", default="all", choices=[*TARGETS, "all"])
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    FIGURES.mkdir(exist_ok=True)
    plots = not args.no_plots

    chosen = TARGETS.keys() if args.target == "all" else [args.target]
    for name in chosen:
        TARGETS[name](plots=plots)

    if args.target == "all":
        print_summary()


if __name__ == "__main__":
    main()
