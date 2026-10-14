#!/usr/bin/env python
"""P2 of `../Plan_Beyond_Trot_Into_Route_03_and_Adam.md`: give the step controller a real `tol`,
and use tolerance proportionality (the replication's T2' test) as its acceptance test.

Neither shipped controller has a tolerance to sweep. Route 3 ratchets (GROW/SHRINK), and the
Adam controller uses energy_check's reject band as its setpoint, which G4 showed is circular:
it can only ever see in-band errors, so growth < 1 on 0/190 steps. Here the memo's Eq. (17)
PI law is driven by G5's zero-cost estimator e_n = |dE_actual - dE_pred|/2 against an explicit
`tol`, while DE_BAND_* stays exactly what it was: the reject gate.

Two arms, because the plan's pre-registered alpha ~ 1 only follows for one of them:
  eps   error per step       e_n        <= tol   (Delta tau ~ tol^(1/2), global ~ tol^(1/2))
  epus  error per unit step  e_n / dtau <= tol   (Delta tau ~ tol,       global ~ tol)

Nothing in `../../Attack_Route_03/adaptive_qite.py`, `../../Adam_QITE/adam_qite.py` or
`../Proposal_1/p1_singularity.py` is modified; all three are imported. `controlled_qite_run`
is adam_qite_run's loop with the controller made pluggable, and target `ids` proves the copy
is faithful by reproducing adam_qite_run bit for bit through an Adam adapter.

Targets, one function each:

  ids   check_loop_faithful()         the pluggable loop == adam_qite_run, bit for bit (E1)
  ref   p2_reference()                exact-metric DOP853 reference, its self-check, and the
                                      fd-vs-exact ODE gap that floors every error below
  eps   p2_sweep_eps()                PI on error per step, tol = 1e-6 .. 1e-2
  epus  p2_sweep_epus()               PI on error per unit step, tol = 10^-3.5 .. 1
  fit   p2_fit_and_predictions()      alpha per arm, circuit slope, the break test, figures
  ablate p2_first_step_ablation()     POST HOC, not pre-registered: the fixed first step, and
                                      the sweep rerun with it at DELTA_TAU_MIN

Run with the physics environment:

    ~/micromamba/envs/qml-julia/bin/python p2_tol_proportionality.py --target all 2>&1 \
      | tee Output.txt

`--check` runs everything without figures and exits nonzero if any pinned number drifts. No
tolerance flag: floats at PIN_RTOL, integers exact. A drift is a result to record, never a
reason to edit PINNED until it agrees.

See `P2_Tol_Proportionality_Findings.md` for the write-up.
"""

from __future__ import annotations

import argparse
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
QITE = HERE.parent.parent
sys.path.insert(0, str(QITE / "Attack_Route_03"))
sys.path.insert(0, str(QITE / "Adam_QITE"))
sys.path.insert(0, str(HERE.parent / "Proposal_1"))

import adaptive_qite as aq  # noqa: E402  (also selects matplotlib's Agg backend)
import adam_qite as ad  # noqa: E402
import p1_singularity as p1  # noqa: E402  (exact_metric_and_force only)
import matplotlib.pyplot as plt  # noqa: E402

FIGURES = HERE / "figures"
COLORS = aq.COLORS

# -- Frozen before the first controller run --------------------------------------------
N_PARAMS = 5                      # Target A: SEED_A theta_0, full-rank 5-param ansatz
K_I_PI = 0.3 / 2                  # memo Eq. (17), k_I ~ 0.3/(p+1) at p = 1 (explicit Euler)
K_P_PI = 0.4 / 2                  # memo Eq. (17), k_P ~ 0.4/(p+1) at p = 1
E_FLOOR_PI = 1e-16                # e_n floor before the power law; an exact prediction happens
MAX_STEPS_SWEEP = 5000            # ours; the smallest tols need ~1-2k steps to reach tau_max
TAU_REF = 3.6                     # reference horizon; runs end in [3, 3 + DELTA_TAU_MAX)
REF_RTOL, REF_RTOL_CHECK = 1e-11, 1e-13   # exact-metric reference and its self-check
FD_REF_RTOL = 1e-7                # tightest fd-ODE solve that finishes (1e-8 stalls on fd noise)
TOLS_EPS = tuple(10.0 ** (-6 + 0.5 * i) for i in range(9))      # 1e-6 .. 1e-2
TOLS_EPUS = tuple(10.0 ** (-3.5 + 0.5 * i) for i in range(8))   # 10^-3.5 .. 1
FIT_MAX_CAP_FRAC = 0.25           # a run enters the fit only if the cond cap binds on < 25 %
FIT_MAX_BOUND_FRAC = 0.25         # ... Delta tau sits at DELTA_TAU_MIN/MAX on < 25 %
FIT_MIN_FLOOR_RATIO = 10.0        # ... and its error is > 10x the fd-vs-exact reference floor
PIN_RTOL = 1e-9                   # float reproducibility margin on a fixed BLAS, not a fit knob

ASSUMPTIONS: list[str] = [
    "Every physical choice (H2, the 5-param ansatz, SEED_A theta_0, drho at h=1e-6, rcond "
    "1e-6, the cond(M) cap with COEFF 2000, energy_check's gate and DE_BAND_*, GROW, SHRINK, "
    "DELTA_TAU0/MIN/MAX, MAX_RETRIES, TAU_MAX) is imported from adaptive_qite.py unchanged. "
    "Its own assumption list applies.",

    "OURS: two tolerance semantics. The plan's prediction (alpha ~ 1) is pre-registered for "
    "'explicit Euler with local error control'; for error-per-step control (e_n <= tol) with a "
    "p = 1 method, e_n ~ C dtau^2 gives dtau ~ tol^(1/2) and the memo's own E_tot ~ C T dtau^p "
    "gives global ~ tol^(1/2). alpha ~ 1 needs error per unit step (e_n/dtau <= tol). Both are "
    "run; the plan's prediction is scored as written against the per-step arm.",

    "OURS: e_n = |dE_actual - dE_pred| / 2, the estimator G5 validated against step doubling "
    "(ratio 2.060 median over the first 40 E1 steps); the /2 makes it estimate what step "
    "doubling estimates. It is an energy-space proxy, not a theta-space error.",

    "OURS: the PI law is memo Eq. (17), dtau_{n+1} = dtau_n (tol/e_n)^k_I (e_{n-1}/e_n)^k_P, at "
    "k_I = 0.15, k_P = 0.2 (p = 1). The first step is I-only. The factor is clipped to "
    "[SHRINK, GROW] (imported, not retuned) and then to [DELTA_TAU_MIN, DELTA_TAU_MAX] exactly as "
    "adam_qite_run clips. It updates once per accepted step on the surviving trial's e_n, as the "
    "Adam controller does; since tol is now separate from the band, G4's circularity is broken "
    "without feeding it rejected trials (G4 item 2 is not needed and not taken).",

    "OURS: rejection stays energy_check's band, unchanged. There is no reject-on-e_n>tol (the "
    "noiseless form of memo Eq. 14); that alternative is untested here.",

    "OURS: the reference is the EXACT-metric ODE theta_dot = pinv(M_exact, 1e-6) Y_exact "
    "(P1's exact_metric_and_force) solved by DOP853 at rtol 1e-11, atol 1e-14, self-checked "
    "against rtol 1e-13. The fd ODE the controllers integrate cannot be solved that tightly "
    "(DOP853 at rtol 1e-8 stalls on the ~1e-10 fd noise, found while sizing the reference), so "
    "the fd-vs-exact gap is measured by an fd DOP853 solve at rtol 1e-7 and reported as a "
    "floor: errors below 10x it are shown but kept out of the fit.",

    "OURS: the primary error is |E_run(tau_end) - E_ref(tau_end)| at each run's own tau_end; "
    "the secondary is the max over accepted steps of |E_run(tau_n) - E_ref(tau_n)|. The plan's "
    "'final E - E_floor' is printed too but is the wrong quantity for proportionality: at "
    "tau = 3 the flow itself has not converged, so it mixes truncation error with the distance "
    "the exact flow still has to go.",

    "OURS, fixed before the sweep: a run enters the alpha fit iff it exits on tau_max, the cap "
    "binds on < 25 % of its steps, Delta tau sits at DELTA_TAU_MIN/MAX on < 25 %, and its "
    "primary error is > 10x the reference floor. alpha is the least-squares log-log slope over "
    "those runs. The 10x floor criterion was added after sizing the reference and before any "
    "controller run.",

    "OURS: the 'break tol' is the first grid tol (ascending) whose local slope to the previous "
    "point is < alpha_fit/2; the 'cap tol' is the first grid tol with >= 25 % of steps "
    "cap-bound. The plan's break prediction holds iff both exist within one grid point.",

    "No Adam arm is in the sweep: the plan's precondition (the G2 fix, beta1 = beta2) is not met.",
]

# Pinned from the first full run (Output.txt). nan = no fit possible.
PINNED: dict[str, float | int] = {
    'ids/adam_E': -1.122851752157184,
    'ids/adam_steps': 190,
    'ids/adam_circ': 2285,
    'ref/E3': -1.1235788430747347,
    'ref/self_check': 8.251621608224013e-12,
    'ref/floor_all': 3.418269478583369e-08,
    'ref/floor_end': 1.1866818638850418e-10,
    'eps/1e-06/err': 1.1066818965232983e-05,
    'eps/1e-06/steps': 1234,
    'eps/1e-06/circ': 14808,
    'eps/3.16e-06/err': 1.9021556839859954e-05,
    'eps/3.16e-06/steps': 761,
    'eps/3.16e-06/circ': 9132,
    'eps/1e-05/err': 2.524228554090513e-05,
    'eps/1e-05/steps': 497,
    'eps/1e-05/circ': 5964,
    'eps/3.16e-05/err': 2.3316159837749595e-05,
    'eps/3.16e-05/steps': 351,
    'eps/3.16e-05/circ': 4212,
    'eps/0.0001/err': 1.4298887411889893e-05,
    'eps/0.0001/steps': 270,
    'eps/0.0001/circ': 3240,
    'eps/0.000316/err': 5.670167167570739e-06,
    'eps/0.000316/steps': 225,
    'eps/0.000316/circ': 2700,
    'eps/0.001/err': 3.6386614772343506e-05,
    'eps/0.001/steps': 200,
    'eps/0.001/circ': 2400,
    'eps/0.00316/err': 0.00042490165675235225,
    'eps/0.00316/steps': 184,
    'eps/0.00316/circ': 2208,
    'eps/0.01/err': 0.0002543831119601325,
    'eps/0.01/steps': 180,
    'eps/0.01/circ': 2162,
    'epus/0.000316/err': 2.947265042618774e-05,
    'epus/0.000316/steps': 2260,
    'epus/0.000316/circ': 27120,
    'epus/0.001/err': 3.081567639329208e-05,
    'epus/0.001/steps': 822,
    'epus/0.001/circ': 9864,
    'epus/0.00316/err': 2.453248806211228e-05,
    'epus/0.00316/steps': 371,
    'epus/0.00316/circ': 4452,
    'epus/0.01/err': 1.1781226647133991e-05,
    'epus/0.01/steps': 228,
    'epus/0.01/circ': 2736,
    'epus/0.0316/err': 0.00024098331224364777,
    'epus/0.0316/steps': 186,
    'epus/0.0316/circ': 2233,
    'epus/0.1/err': 0.00047658629697866495,
    'epus/0.1/steps': 176,
    'epus/0.1/circ': 2115,
    'epus/0.316/err': 0.0007060704835570863,
    'epus/0.316/steps': 172,
    'epus/0.316/circ': 2068,
    'epus/1/err': 0.0024411607664969193,
    'epus/1/steps': 56,
    'epus/1/circ': 676,
    'fit/eps/alpha': float("nan"),
    'fit/eps/alpha_max': float("nan"),
    'fit/eps/circ': float("nan"),
    'fit/epus/alpha': float("nan"),
    'fit/epus/alpha_max': float("nan"),
    'fit/epus/circ': float("nan"),
    'fit/route3/err': 0.00063180598613366,
    'fit/adam/err': 0.0007277250511912303,
    'pred/1': 0,
    'pred/2': 0,
    'pred/3': 0,
    'pred/4': 0,
    'pred/5': 1,
    'pred/6': 0,
    'fit/E_floor': -1.1246000000000003,
    'abl/eps/alpha_win': 0.4782977470619917,
    'abl/eps/alpha_win_default': 0.26735933546016627,
    'abl/eps/circ': -0.38544745563212884,
    'abl/eps/1e-06/err_win': 0.0024687571204682346,
    'abl/eps/0.01/err_win': 0.08184707469961139,
    'abl/eps/err_first': 0.0005291381899473757,
    'abl/epus/alpha_win': 0.8897683580668071,
    'abl/epus/alpha_win_default': 0.5145184801308876,
    'abl/epus/circ': -0.7653930464344677,
    'abl/epus/0.000316/err_win': 0.002056837691229718,
    'abl/epus/1/err_win': 0.08855690855871301,
    'abl/epus/err_first': 0.0005291381899473757,
}

DRIFTS: list[str] = []
UNPINNED: dict[str, float | int] = {}
N_COMPARED = 0


def compare(label, measured, key, indent=4):
    """Print measured vs pinned and record drift. Never raises."""
    global N_COMPARED
    pad = " " * indent
    if key not in PINNED:
        print(f"{pad}{label:<44s} {measured!r}   (not pinned)")
        UNPINNED[key] = int(measured) if isinstance(measured, (int, np.integer)) else float(measured)
        return
    N_COMPARED += 1
    pinned = PINNED[key]
    if isinstance(pinned, int):
        ok = int(measured) == pinned
    elif np.isnan(pinned):
        ok = bool(np.isnan(measured))    # "no fit possible" is a pinned outcome too
    else:
        ok = abs(measured - pinned) <= PIN_RTOL * abs(pinned)
    print(f"{pad}{label:<44s} {measured!r}   pinned {pinned!r}   {'MATCH' if ok else 'DRIFT'}")
    if not ok:
        DRIFTS.append(f"{key}: measured {measured!r}, pinned {pinned!r}")


def claim(label, ok, indent=4):
    """A yes/no statement that must hold. Recorded as a drift if it fails."""
    global N_COMPARED
    N_COMPARED += 1
    print(f"{' ' * indent}{label:<70s} {'HOLDS' if ok else 'FAILS'}")
    if not ok:
        DRIFTS.append(f"claim failed: {label}")


def prediction(label, ok, key, indent=4):
    """A pre-registered prediction. Its outcome is a RESULT: printed HELD/FAILED and pinned, so
    a failed prediction does not fail --check, but a change in its outcome does."""
    print(f"{' ' * indent}prediction: {label:<58s} {'HELD' if ok else 'FAILED'}")
    compare("  outcome (1 = held)", int(ok), key, indent=indent)


def banner(text):
    aq.banner(text)


def tkey(tol):
    return f"{tol:.3g}"


# --------------------------------------------------------------------------------------
# Controllers and the pluggable loop
# --------------------------------------------------------------------------------------

class PIStepController:
    """Memo Eq. (17) on G5's estimator, against an explicit tol.

    mode 'eps':  e_n = |dE_actual - dE_pred| / 2            (error per step)
    mode 'epus': e_n = |dE_actual - dE_pred| / 2 / dtau     (error per unit step)
    """

    def __init__(self, tol, mode, k_i=K_I_PI, k_p=K_P_PI, grow=aq.GROW, shrink=aq.SHRINK):
        if mode not in ("eps", "epus"):
            raise ValueError(mode)
        self.tol, self.mode, self.k_i, self.k_p = tol, mode, k_i, k_p
        self.grow, self.shrink = grow, shrink
        self.e_prev = None

    def update(self, err_raw, band, r_hat, delta_tau):
        e = err_raw / 2.0
        if self.mode == "epus":
            e = e / delta_tau
        e = max(e, E_FLOOR_PI)
        factor = (self.tol / e) ** self.k_i
        if self.e_prev is not None:
            factor *= (self.e_prev / e) ** self.k_p
        self.e_prev = e
        growth = float(np.clip(factor, self.shrink, self.grow))
        return dict(growth=growth, e=e, raw=factor)


class AdamAdapter:
    """ad.AdamStepController behind the same interface -- used only to prove the loop is
    adam_qite_run's, bit for bit."""

    def __init__(self):
        self.ctrl = ad.AdamStepController()

    def update(self, err_raw, band, r_hat, delta_tau):
        state = self.ctrl.update(err_raw / band, r_hat)
        return dict(growth=state["growth"], e=err_raw / band, raw=state["growth"])


def controlled_qite_run(theta0, rho_fn, H, n_params, N, controller, tau_max=aq.TAU_MAX,
                        delta_tau0=aq.DELTA_TAU0, rcond=aq.RCOND_DEFAULT,
                        stiffness_cap_coeff=aq.STIFFNESS_CAP_COEFF,
                        delta_tau_min=aq.DELTA_TAU_MIN, delta_tau_max=aq.DELTA_TAU_MAX,
                        max_steps=aq.MAX_STEPS, max_retries=aq.MAX_RETRIES):
    """adam_qite_run's loop (noiseless path), line for line, with the step law pluggable.
    Measurement, cond(M) cap, energy_check gate and retry loop are unchanged."""
    theta = np.array(theta0, dtype=float)
    tau = 0.0
    delta_tau = delta_tau0
    E0, _ = aq.energy_and_variance(rho_fn(theta), H)

    log = dict(tau=[0.0], E=[E0], delta_tau=[delta_tau0], cond_M=[np.nan], accepted=[True],
               retries=[0], cap_bound=[False], at_bound=[False], e=[np.nan], raw=[np.nan],
               growth=[np.nan], n_circ=0)

    step = 0
    while tau < tau_max and step < max_steps:
        theta_dot_val, M, Y, rho_n = aq.measure_f(theta, rho_fn, H, N, 0.0, rcond, True)
        log["n_circ"] += 2 * N + 1
        E_n, var_n = aq.energy_and_variance(rho_n, H)
        residual, cond_M, sigma_min, cond_spectral = aq.route3_indicators(M, Y, theta_dot_val,
                                                                          rcond)
        r_hat = residual / (np.linalg.norm(Y) + 1e-12)

        if np.isinf(stiffness_cap_coeff):
            cap = np.inf
        else:
            cap = stiffness_cap_coeff / cond_M if cond_M > 0 else np.inf
        cap_bound = cap < delta_tau and cap < delta_tau_max
        delta_tau = max(delta_tau_min, min(delta_tau, cap, delta_tau_max))
        at_bound = delta_tau >= delta_tau_max or delta_tau <= delta_tau_min

        retries = 0
        while True:
            theta_trial = theta + delta_tau * theta_dot_val
            E_trial, _ = aq.energy_and_variance(rho_fn(theta_trial), H)
            log["n_circ"] += 1
            dE_pred, _, dE_actual, reasons = aq.energy_check(
                E_n, var_n, Y, theta_dot_val, delta_tau, E_trial)
            accepted = not reasons
            if accepted or retries >= max_retries or delta_tau <= delta_tau_min:
                break
            delta_tau = max(delta_tau * aq.SHRINK, delta_tau_min)
            retries += 1

        band = max(aq.DE_BAND_ABS, aq.DE_BAND_REL * abs(dE_pred))
        state = controller.update(abs(dE_actual - dE_pred), band, r_hat, delta_tau)

        theta = theta_trial
        tau += delta_tau
        for key, val in (("tau", tau), ("E", E_trial), ("delta_tau", delta_tau),
                         ("cond_M", cond_M), ("accepted", accepted), ("retries", retries),
                         ("cap_bound", cap_bound), ("at_bound", at_bound), ("e", state["e"]),
                         ("raw", state["raw"]), ("growth", state["growth"])):
            log[key].append(val)

        delta_tau = float(np.clip(delta_tau * state["growth"], delta_tau_min, delta_tau_max))
        step += 1

    log["exit_reason"] = "tau_max" if tau >= tau_max else "max_steps"
    log["theta"] = theta
    for key in ("tau", "E", "delta_tau", "cond_M", "retries", "e", "raw", "growth"):
        log[key] = np.array(log[key], dtype=float)
    for key in ("accepted", "cap_bound", "at_bound"):
        log[key] = np.array(log[key], dtype=bool)
    return log


def _setup():
    H = aq.h2_hamiltonian()
    theta0 = np.random.default_rng(aq.SEED_A).uniform(-0.5, 0.5, size=N_PARAMS)
    return H, theta0, aq.make_rho_fn(n_params=N_PARAMS)


# --------------------------------------------------------------------------------------
# Target ids -- the loop copy is faithful
# --------------------------------------------------------------------------------------

def check_loop_faithful(plots=True):
    banner("ids  the pluggable loop reproduces adam_qite_run bit for bit (E1, SEED_A)")
    H, theta0, rho_fn = _setup()
    ref = ad.adam_qite_run(theta0, rho_fn, H, N_PARAMS, N_PARAMS)
    mine = controlled_qite_run(theta0, rho_fn, H, N_PARAMS, N_PARAMS, AdamAdapter())
    same = (len(ref["tau"]) == len(mine["tau"])
            and np.array_equal(ref["E"], mine["E"])
            and np.array_equal(ref["delta_tau"], mine["delta_tau"])
            and np.array_equal(ref["growth"][1:], mine["growth"][1:])
            and ref["n_circ"] == mine["n_circ"])
    compare("adam_qite_run final E", float(ref["E"][-1]), "ids/adam_E")
    compare("adam_qite_run steps", len(ref["tau"]) - 1, "ids/adam_steps")
    compare("adam_qite_run circuits", int(ref["n_circ"]), "ids/adam_circ")
    claim("E, delta_tau, growth and circuits identical, step by step", same)


# --------------------------------------------------------------------------------------
# Target ref -- the reference solution and its floor
# --------------------------------------------------------------------------------------

def _exact_f(H):
    def f(tau, theta):
        M, Y = p1.exact_metric_and_force(theta, N_PARAMS, H)[:2]
        return np.linalg.pinv(M, rcond=aq.RCOND_DEFAULT) @ Y
    return f


@lru_cache(maxsize=None)
def reference():
    """Exact-metric reference (dense), the self-check difference, and the fd-vs-exact floors."""
    H, theta0, rho_fn = _setup()
    kw = dict(method="DOP853", dense_output=True)
    sol = solve_ivp(_exact_f(H), (0.0, TAU_REF), theta0, rtol=REF_RTOL, atol=REF_RTOL * 1e-3,
                    **kw)
    chk = solve_ivp(_exact_f(H), (0.0, TAU_REF), theta0, rtol=REF_RTOL_CHECK,
                    atol=REF_RTOL_CHECK * 1e-3, **kw)

    def fd_f(tau, theta):
        return aq.measure_f(theta, rho_fn, H, N_PARAMS, 0.0, aq.RCOND_DEFAULT, True)[0]

    fd = solve_ivp(fd_f, (0.0, TAU_REF), theta0, rtol=FD_REF_RTOL, atol=FD_REF_RTOL * 1e-3, **kw)

    grid = np.linspace(0.0, TAU_REF, 3601)

    def energies(s):
        return np.array([aq.energy_and_variance(rho_fn(s.sol(x)), H)[0] for x in grid])

    E_ref, E_chk, E_fd = energies(sol), energies(chk), energies(fd)
    end = grid >= aq.TAU_MAX
    return dict(sol=sol, nfev=(sol.nfev, chk.nfev, fd.nfev),
                self_check=float(np.max(np.abs(E_ref - E_chk))),
                floor_all=float(np.max(np.abs(E_ref - E_fd))),
                floor_end=float(np.max(np.abs(E_ref - E_fd)[end])),
                E_ref_3=float(E_ref[np.searchsorted(grid, aq.TAU_MAX)]))


def E_ref_at(taus):
    H, _, rho_fn = _setup()
    sol = reference()["sol"]
    return np.array([aq.energy_and_variance(rho_fn(sol.sol(t)), H)[0] for t in np.atleast_1d(taus)])


def p2_reference(plots=True):
    banner("ref  exact-metric reference solution, self-check, fd-vs-exact floor")
    r = reference()
    print(f"  DOP853 evaluations: exact rtol {REF_RTOL:g}: {r['nfev'][0]}, exact rtol "
          f"{REF_RTOL_CHECK:g}: {r['nfev'][1]}, fd rtol {FD_REF_RTOL:g}: {r['nfev'][2]}")
    compare("E_ref(tau = 3)", r["E_ref_3"], "ref/E3")
    compare("self-check max |E(1e-11) - E(1e-13)|", r["self_check"], "ref/self_check")
    compare("fd-vs-exact floor, max over tau", r["floor_all"], "ref/floor_all")
    compare("fd-vs-exact floor, tau in [3, 3.6]", r["floor_end"], "ref/floor_end")
    claim("reference self-check < 1e-10 (well below every floor)", r["self_check"] < 1e-10)


# --------------------------------------------------------------------------------------
# Targets eps / epus -- the sweeps
# --------------------------------------------------------------------------------------

def run_metrics(run):
    tau_end = float(run["tau"][-1])
    E_ref_traj = E_ref_at(run["tau"][1:])
    theta_ref = reference()["sol"].sol(tau_end)
    n = len(run["tau"]) - 1
    return dict(
        tau_end=tau_end,
        err=float(abs(run["E"][-1] - E_ref_traj[-1])),
        err_max=float(np.max(np.abs(run["E"][1:] - E_ref_traj))),
        dtheta=float(np.linalg.norm(run["theta"] - theta_ref)),
        steps=n, circ=int(run["n_circ"]), retries=int(np.sum(run["retries"][1:])),
        cap_frac=float(np.mean(run["cap_bound"][1:])),
        bound_frac=float(np.mean(run["at_bound"][1:])),
        grow_clip=int(np.sum(run["raw"][1:] > aq.GROW)),
        shrink_clip=int(np.sum(run["raw"][1:] < aq.SHRINK)),
        growth_lt1=int(np.sum(run["growth"][1:] < 1.0)),
        exit=run["exit_reason"], E_end=float(run["E"][-1]))


@lru_cache(maxsize=None)
def sweep(mode):
    H, theta0, rho_fn = _setup()
    tols = TOLS_EPS if mode == "eps" else TOLS_EPUS
    rows, runs = [], {}
    for tol in tols:
        run = controlled_qite_run(theta0, rho_fn, H, N_PARAMS, N_PARAMS,
                                  PIStepController(tol, mode), max_steps=MAX_STEPS_SWEEP)
        m = run_metrics(run)
        m["tol"] = tol
        m["e_over_tol"] = float(np.median(run["e"][1:]) / tol)
        rows.append(m)
        runs[tol] = run
    return rows, runs


def _in_fit(m, floor):
    return (m["exit"] == "tau_max" and m["cap_frac"] < FIT_MAX_CAP_FRAC
            and m["bound_frac"] < FIT_MAX_BOUND_FRAC and m["err"] > FIT_MIN_FLOOR_RATIO * floor)


def _print_sweep(mode):
    rows, _ = sweep(mode)
    floor = reference()["floor_end"]
    print(f"  {'tol':>9} {'err(tau_end)':>12} {'max err':>10} {'|dtheta|':>9} {'steps':>6} "
          f"{'circ':>6} {'retr':>4} {'cap%':>5} {'bnd%':>5} {'g<1':>4} {'clipG':>5} "
          f"{'med e/tol':>9} {'fit':>3}")
    for m in rows:
        print(f"  {m['tol']:9.3g} {m['err']:12.4e} {m['err_max']:10.3e} {m['dtheta']:9.2e} "
              f"{m['steps']:6d} {m['circ']:6d} {m['retries']:4d} {100 * m['cap_frac']:5.1f} "
              f"{100 * m['bound_frac']:5.1f} {m['growth_lt1']:4d} {m['grow_clip']:5d} "
              f"{m['e_over_tol']:9.3g} {'yes' if _in_fit(m, floor) else 'no':>3}")
    print()
    for m in rows:
        k = f"{mode}/{tkey(m['tol'])}"
        compare(f"tol {tkey(m['tol'])} error at tau_end", m["err"], f"{k}/err")
        compare(f"tol {tkey(m['tol'])} steps", m["steps"], f"{k}/steps")
        compare(f"tol {tkey(m['tol'])} circuits", m["circ"], f"{k}/circ")
    bad = [m["tol"] for m in rows if m["exit"] != "tau_max"]
    if bad:
        print(f"    NOTE runs that hit MAX_STEPS_SWEEP: {bad}")


def p2_sweep_eps(plots=True):
    banner("eps  PI (memo Eq. 17) on error per step, e_n = |dE_actual - dE_pred|/2 <= tol")
    _print_sweep("eps")


def p2_sweep_epus(plots=True):
    banner("epus  PI (memo Eq. 17) on error per unit step, e_n / dtau <= tol")
    _print_sweep("epus")


# --------------------------------------------------------------------------------------
# Target fit -- alpha, predictions, figures
# --------------------------------------------------------------------------------------

def _slope(x, y):
    return float(np.polyfit(np.log10(x), np.log10(y), 1)[0])


def _arm_fit(mode):
    rows, _ = sweep(mode)
    floor = reference()["floor_end"]
    inc = [m for m in rows if _in_fit(m, floor)]
    tol = np.array([m["tol"] for m in rows])
    err = np.array([m["err"] for m in rows])
    local = np.diff(np.log10(err)) / np.diff(np.log10(tol))
    fit = dict(n=len(inc), alpha=np.nan, alpha_max=np.nan, circ_slope=np.nan)
    if len(inc) >= 3:
        fit["alpha"] = _slope([m["tol"] for m in inc], [m["err"] for m in inc])
        fit["alpha_max"] = _slope([m["tol"] for m in inc], [m["err_max"] for m in inc])
        fit["circ_slope"] = _slope([m["tol"] for m in inc], [m["circ"] for m in inc])
    fit["local"] = local
    fit["i_break"] = next((i for i in range(1, len(rows))
                           if np.isfinite(fit["alpha"]) and local[i - 1] < fit["alpha"] / 2), None)
    fit["i_cap"] = next((i for i, m in enumerate(rows) if m["cap_frac"] >= FIT_MAX_CAP_FRAC), None)
    return rows, fit


def _existing_controllers():
    """Route 3 (Target A) and Adam (E1) on the same reference, as work-precision anchors."""
    H, theta0, rho_fn = _setup()
    r3 = aq.adaptive_qite_run(theta0, rho_fn, H, N_PARAMS, N_PARAMS)
    adam = ad.adam_qite_run(theta0, rho_fn, H, N_PARAMS, N_PARAMS)
    out = {}
    for name, run in (("route3", r3), ("adam", adam)):
        tau_end = float(run["tau"][-1])
        out[name] = dict(err=float(abs(run["E"][-1] - E_ref_at(tau_end)[0])),
                         steps=len(run["tau"]) - 1, circ=int(run["n_circ"]))
    return out


def p2_fit_and_predictions(plots=True):
    banner("fit  tolerance proportionality: alpha per arm, circuits, the break, predictions")
    floor = reference()["floor_end"]
    print(f"  fit set: exit on tau_max, cap-bound < {FIT_MAX_CAP_FRAC:.0%}, Delta tau at a bound "
          f"< {FIT_MAX_BOUND_FRAC:.0%}, error > {FIT_MIN_FLOOR_RATIO:g} x floor "
          f"({FIT_MIN_FLOOR_RATIO * floor:.3e})\n")
    fits = {}
    for mode in ("eps", "epus"):
        rows, fit = _arm_fit(mode)
        fits[mode] = (rows, fit)
        print(f"  [{mode}] runs in the fit: {fit['n']} of {len(rows)}"
              + ("" if fit["n"] >= 3 else "  -> NO FIT POSSIBLE (needs >= 3): alpha and circuit "
                 "slope are nan"))
        print(f"    local slopes d log err / d log tol: "
              + " ".join(f"{s:+.2f}" for s in fit["local"]))
        compare("alpha (error at tau_end)", fit["alpha"], f"fit/{mode}/alpha")
        compare("alpha (max error along the run)", fit["alpha_max"], f"fit/{mode}/alpha_max")
        compare("circuit slope d log circ / d log tol", fit["circ_slope"], f"fit/{mode}/circ")
        tols = [m["tol"] for m in rows]
        print(f"    break tol: {tkey(tols[fit['i_break']]) if fit['i_break'] is not None else 'none'}"
              f";  cap tol (>= 25 % cap-bound): "
              f"{tkey(tols[fit['i_cap']]) if fit['i_cap'] is not None else 'none'}\n")

    anchors = _existing_controllers()
    print("  existing controllers on the same reference (work-precision anchors):")
    for name, a in anchors.items():
        print(f"    {name:<7} error at tau_end {a['err']:.4e}   steps {a['steps']}   "
              f"circuits {a['circ']}")
        compare(f"{name} error at tau_end", a["err"], f"fit/{name}/err", indent=6)
    print()

    eps_fit, epus_fit = fits["eps"][1], fits["epus"][1]
    a_eps, a_epus = eps_fit["alpha"], epus_fit["alpha"]
    print("  Predictions (written before the sweep, see the findings note):")
    if not (np.isfinite(a_eps) and np.isfinite(a_epus)):
        print("    NOTE a nan alpha scores FAILED: 1-3, 4 and 6 then fail because no slope could be\n"
              "    fitted under the pre-registered fit criteria, not because one was measured out\n"
              "    of range.")
    prediction("1 (plan, as written) eps alpha in [0.8, 1.2]",
               bool(0.8 <= a_eps <= 1.2), "pred/1")
    prediction("2 eps alpha in [0.4, 0.6]", bool(0.4 <= a_eps <= 0.6), "pred/2")
    prediction("3 epus alpha in [0.8, 1.2]", bool(0.8 <= a_epus <= 1.2), "pred/3")
    ib, ic = eps_fit["i_break"], eps_fit["i_cap"]
    prediction("4 (plan) eps slope breaks where the cap takes over",
               ib is not None and ic is not None and abs(ib - ic) <= 1, "pred/4")
    all_rows = fits["eps"][0] + fits["epus"][0]
    prediction("5 growth < 1 on > 0 steps in every run",
               all(m["growth_lt1"] > 0 for m in all_rows), "pred/5")
    prediction("6 eps circuit slope in [-0.6, -0.4]",
               bool(-0.6 <= eps_fit["circ_slope"] <= -0.4), "pred/6")
    print(f"\n  final E - E_floor (plan's quantity, for the record; E_floor from variational_floor):")
    H, _, rho_fn = _setup()
    E_floor = aq.variational_floor(rho_fn, N_PARAMS, H)
    compare("E_floor", E_floor, "fit/E_floor")
    print(f"    E_ref(3) - E_floor = {reference()['E_ref_3'] - E_floor:.4e}  <- the EXACT flow's "
          f"own distance at tau = 3")
    for mode in ("eps", "epus"):
        rows = fits[mode][0]
        print(f"    {mode:<5} " + "  ".join(f"{m['E_end'] - E_floor:.3e}" for m in rows))

    if plots:
        _figure(fits, anchors, floor)


# --------------------------------------------------------------------------------------
# Target ablate -- POST HOC, not pre-registered: why the error does not track tol
# --------------------------------------------------------------------------------------

ABL_WINDOW = 2.0                  # post hoc: max error over tau <= 2, before the cap binds
ABL_N_ASYM = 3                    # post hoc: slopes over the 3 smallest tols of each arm


@lru_cache(maxsize=None)
def sweep_small_first_step(mode):
    """The same sweep with the first step at DELTA_TAU_MIN instead of the fixed DELTA_TAU0,
    so the controller sets Delta tau from the first step on."""
    H, theta0, rho_fn = _setup()
    tols = TOLS_EPS if mode == "eps" else TOLS_EPUS
    return {tol: controlled_qite_run(theta0, rho_fn, H, N_PARAMS, N_PARAMS,
                                     PIStepController(tol, mode), delta_tau0=aq.DELTA_TAU_MIN,
                                     max_steps=MAX_STEPS_SWEEP)
            for tol in tols}


def _window_metrics(run):
    tau = run["tau"][1:]
    err = np.abs(run["E"][1:] - E_ref_at(tau))
    w = tau <= ABL_WINDOW
    cap_tau = tau[run["cap_bound"][1:]]
    return dict(err_win=float(err[w].max()), tau_peak=float(tau[w][np.argmax(err[w])]),
                err_first=float(err[0]), e_first=float(run["e"][1]),
                first_cap=float(cap_tau.min()) if cap_tau.size else np.inf,
                circ=int(run["n_circ"]), err=float(err[-1]), exit=run["exit_reason"],
                tau=tau, err_traj=err)


def p2_first_step_ablation(plots=True):
    banner("ablate  POST HOC: the fixed first step, and the error before the cap binds")
    print("  Not pre-registered. Added after the fit target: the max error along every run peaks\n"
          "  near tau ~ 1 whatever tol is, and every default run starts with the same fixed\n"
          f"  step DELTA_TAU0 = {aq.DELTA_TAU0:g}. Here: (i) that step's error, (ii) the sweep rerun "
          f"with the\n  first step at DELTA_TAU_MIN = {aq.DELTA_TAU_MIN:g}, scored on the max error over "
          f"tau <= {ABL_WINDOW:g}\n  (the cap never binds there, checked below). Slopes over the "
          f"{ABL_N_ASYM} smallest tols.\n")
    out = {}
    for mode in ("eps", "epus"):
        _, runs = sweep(mode)
        default = {tol: _window_metrics(r) for tol, r in runs.items()}
        small = {tol: _window_metrics(r) for tol, r in sweep_small_first_step(mode).items()}
        out[mode] = (default, small)
        tols = sorted(default)
        firsts = [default[t]["err_first"] for t in tols]
        print(f"  [{mode}] default runs: error after the fixed first step {firsts[0]:.4e}, "
              f"estimator e_1 {default[tols[0]]['e_first']:.4e}")
        claim(f"[{mode}] first-step error identical in every default run (tol-blind)",
              all(f == firsts[0] for f in firsts))
        print(f"    {'tol':>9} | {'default: max err':>16} {'at tau':>6} | {'small first step: max err':>25}"
              f" {'at tau':>6} {'err(tau_end)':>12} {'circ':>6} {'1st cap tau':>11}")
        for t in tols:
            d, s = default[t], small[t]
            print(f"    {t:9.3g} | {d['err_win']:16.4e} {d['tau_peak']:6.2f} | {s['err_win']:25.4e} "
                  f"{s['tau_peak']:6.2f} {s['err']:12.4e} {s['circ']:6d} {s['first_cap']:11.3f}")
        tol_a = np.array(tols)
        for name, res in (("default", default), ("small first step", small)):
            e = np.array([res[t]["err_win"] for t in tols])
            local = np.diff(np.log10(e)) / np.diff(np.log10(tol_a))
            print(f"    {name:<16} local slopes of max err: " + " ".join(f"{x:+.2f}" for x in local))
        e_s = [small[t]["err_win"] for t in tols[:ABL_N_ASYM]]
        c_s = [small[t]["circ"] for t in tols[:ABL_N_ASYM]]
        compare("small first step: slope of max err (asymptotic)",
                _slope(tols[:ABL_N_ASYM], e_s), f"abl/{mode}/alpha_win")
        compare("default: slope of max err (asymptotic)",
                _slope(tols[:ABL_N_ASYM], [default[t]["err_win"] for t in tols[:ABL_N_ASYM]]),
                f"abl/{mode}/alpha_win_default")
        compare("small first step: circuit slope (asymptotic)", _slope(tols[:ABL_N_ASYM], c_s),
                f"abl/{mode}/circ")
        for t in (tols[0], tols[-1]):
            compare(f"small first step tol {tkey(t)} max err", small[t]["err_win"],
                    f"abl/{mode}/{tkey(t)}/err_win")
        compare("default first-step error", firsts[0], f"abl/{mode}/err_first")
        claim(f"[{mode}] the cap never binds inside tau <= {ABL_WINDOW:g} (either sweep)",
              all(r["first_cap"] > ABL_WINDOW for res in (default, small) for r in res.values()))
        claim(f"[{mode}] every small-first-step run exits on tau_max",
              all(r["exit"] == "tau_max" for r in small.values()))
        print()
    if plots:
        _figure_ablation(out)


def _figure_ablation(out):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    arm_color = {"eps": COLORS[0], "epus": COLORS[1]}
    ax = axes[0]
    for mode, (default, small) in out.items():
        tols = np.array(sorted(default))
        for res, mk, mfc, name in ((default, "o", "white", f"{mode}, fixed first step"),
                                   (small, "o", arm_color[mode], f"{mode}, first step at min")):
            ax.plot(tols, [res[t]["err_win"] for t in tols], mk + "-", color=arm_color[mode],
                    mfc=mfc, ms=6, lw=1.1, label=name)
        e0 = small[tols[0]]["err_win"]
        a = 0.5 if mode == "eps" else 1.0
        span = tols[:ABL_N_ASYM + 2]
        ax.plot(span, e0 * (span / tols[0]) ** a, ":" if a == 0.5 else "--", color="gray",
                lw=0.9, label=f"slope {a:g}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("requested tol"); ax.set_ylabel(rf"max $|E - E_{{ref}}|$ over $\tau\leq{ABL_WINDOW:g}$")
    ax.set_title("(a) post hoc: error before the cap, vs tol"); ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = axes[1]
    default, small = out["eps"]
    for tol, ls in ((TOLS_EPS[0], "-"), (TOLS_EPS[4], "--")):
        ax.plot(default[tol]["tau"], default[tol]["err_traj"], ls, color=COLORS[2], lw=1.2,
                label=f"eps tol={tkey(tol)}, fixed first step")
        ax.plot(small[tol]["tau"], small[tol]["err_traj"], ls, color=COLORS[0], lw=1.2,
                label=f"eps tol={tkey(tol)}, first step at min")
    ax.axvline(ABL_WINDOW, color="gray", ls=":", lw=1.0)
    ax.set_yscale("log")
    ax.set_xlabel(r"imaginary time $\tau$"); ax.set_ylabel(r"$|E(\tau) - E_{ref}(\tau)|$")
    ax.set_title("(b) error along the run, eps arm")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout()
    out_path = FIGURES / "p2_first_step_ablation.png"
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    print(f"  wrote {out_path}")


def _figure(fits, anchors, floor):
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.4))
    arm_color = {"eps": COLORS[0], "epus": COLORS[1]}
    arm_label = {"eps": "per step, $e_n\\leq$ tol", "epus": "per unit step, $e_n/\\Delta\\tau\\leq$ tol"}

    ax = axes[0, 0]
    for mode, (rows, fit) in fits.items():
        tol = np.array([m["tol"] for m in rows])
        err = np.array([m["err"] for m in rows])
        inc = np.array([_in_fit(m, floor) for m in rows])
        ax.plot(tol, err, "-", color=arm_color[mode], lw=1.2, alpha=0.6)
        ax.plot(tol[inc], err[inc], "o", color=arm_color[mode], ms=7,
                label=f"{arm_label[mode]}, " + (f"$\\alpha$={fit['alpha']:.2f}"
                                                 if np.isfinite(fit["alpha"]) else "no fit (<3 runs)"))
        ax.plot(tol[~inc], err[~inc], "o", mfc="white", mec=arm_color[mode], ms=7)
        if inc.any():
            t0, e0 = tol[inc][0], err[inc][0]
            for a, ls in ((0.5, ":"), (1.0, "--")):
                span = np.array([tol.min(), tol.max()])
                ax.plot(span, e0 * (span / t0) ** a, ls, color=arm_color[mode], lw=0.9)
    ax.axhline(floor, color="gray", ls="-.", lw=1.0, label="fd-vs-exact floor")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("requested tol"); ax.set_ylabel(r"$|E(\tau_{end}) - E_{ref}(\tau_{end})|$")
    ax.set_title("(a) achieved error vs tol  (dotted $\\alpha$=1/2, dashed $\\alpha$=1)")
    ax.legend(fontsize=7.5, loc="upper left"); ax.grid(alpha=0.3)

    ax = axes[0, 1]
    for mode, (rows, fit) in fits.items():
        ax.plot([m["circ"] for m in rows], [m["err"] for m in rows], "o-",
                color=arm_color[mode], lw=1.2, ms=6, label=arm_label[mode])
    for (name, a), mk, c in zip(anchors.items(), ("s", "D"), (COLORS[2], COLORS[3])):
        ax.plot(a["circ"], a["err"], mk, color=c, ms=9,
                label={"route3": "Route 3 ratchet (Target A)", "adam": "Adam (E1)"}[name])
    ax.axhline(floor, color="gray", ls="-.", lw=1.0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("circuits (M6 currency)"); ax.set_ylabel("error at $\\tau_{end}$")
    ax.set_title("(b) work-precision"); ax.legend(fontsize=7.5); ax.grid(alpha=0.3)

    ax = axes[1, 0]
    for mode, (rows, fit) in fits.items():
        _, runs = sweep(mode)
        tols = [m["tol"] for m in rows]
        for tol, ls in zip((tols[0], tols[len(tols) // 2], tols[-1]), ("-", "--", ":")):
            run = runs[tol]
            ax.plot(run["tau"][1:], run["delta_tau"][1:], ls, color=arm_color[mode], lw=1.3,
                    label=f"{mode} tol={tkey(tol)}")
    ax.set_yscale("log")
    ax.set_xlabel(r"imaginary time $\tau$"); ax.set_ylabel(r"$\Delta\tau$")
    ax.set_title(r"(c) step size along the run"); ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3)

    ax = axes[1, 1]
    for mode, (rows, fit) in fits.items():
        ax.plot([m["tol"] for m in rows], [100 * m["cap_frac"] for m in rows], "o-",
                color=arm_color[mode], lw=1.2, ms=6, label=f"{mode}: cond(M) cap binding")
        ax.plot([m["tol"] for m in rows], [100 * m["bound_frac"] for m in rows], "s--",
                color=arm_color[mode], lw=1.0, ms=5, mfc="white",
                label=f"{mode}: $\\Delta\\tau$ at min/max")
    ax.axhline(100 * FIT_MAX_CAP_FRAC, color="gray", ls=":", lw=1.0)
    ax.set_xscale("log"); ax.set_ylim(-3, 103)
    ax.set_xlabel("requested tol"); ax.set_ylabel("% of steps")
    ax.set_title("(d) what else is setting $\\Delta\\tau$"); ax.legend(fontsize=7.5); ax.grid(alpha=0.3)

    fig.tight_layout()
    out = FIGURES / "p2_tol_proportionality.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"\n  wrote {out}")


# --------------------------------------------------------------------------------------

def print_summary():
    banner("Assumptions (ours, on top of adaptive_qite.py's and adam_qite.py's)")
    for i, a in enumerate(ASSUMPTIONS, 1):
        print(f"  {i}. {a}\n")


TARGETS = {
    "ids": check_loop_faithful,
    "ref": p2_reference,
    "eps": p2_sweep_eps,
    "epus": p2_sweep_epus,
    "fit": p2_fit_and_predictions,
    "ablate": p2_first_step_ablation,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", default="all", choices=[*TARGETS, "all"])
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--check", action="store_true",
                        help="run every target without figures; exit nonzero on any drift")
    args = parser.parse_args()

    FIGURES.mkdir(exist_ok=True)
    plots = not (args.no_plots or args.check)
    chosen = TARGETS.keys() if (args.target == "all" or args.check) else [args.target]
    for name in chosen:
        TARGETS[name](plots=plots)

    if args.target == "all" and not args.check:
        print_summary()
    banner(f"Ledger: {N_COMPARED} pinned values / claims compared, {len(DRIFTS)} drifted or failed")
    for item in DRIFTS:
        print(f"  - {item}")
    if UNPINNED:
        print(f"  {len(UNPINNED)} value(s) not pinned yet:")
        for key, value in UNPINNED.items():
            print(f"    {key!r}: {value!r},")
    return 1 if DRIFTS else 0


if __name__ == "__main__":
    sys.exit(main())
