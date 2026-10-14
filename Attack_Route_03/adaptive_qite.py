#!/usr/bin/env python
"""Implementation of Sections 1-3 (Route 3 only) of the internal memo "Adaptive
Time-Stepping for Variational QITE: An RKF-Type Controller with a Hardware-Noise Floor"
(`../Research_Papers/QITE_timestep.pdf`, 5pp, no arXiv/DOI).

This is an IMPLEMENTATION of a proposed method, not a replication of reported numbers --
the memo contains no numerical claims to score against. Four scenarios, one function each:

  A  target_a_noiseless_convergence()  Stage 1: the ODE + Route 3 controller, p=0
  B  target_b_noisy_debiasing()        Stage 2: global depolarizing noise + Eq.(7) debiasing
  C  target_c_stiffness_9param()       Route 3's cond(M) stiffness proxy on a rank-deficient M
  D  target_d_indicator_noise_floor()  our own extension of Sec.2's noise-floor idea to Route 3

Builds on the H2 (STO-3G, Bravyi-Kitaev) system and the M/Y/cond(M) machinery already
implemented and validated in
`../Replications/practical-limitations-qdp/replicate.py` (Saxena-Shalabi-Kyaw). That
replication found the 5-parameter ansatz's M is ill-conditioned (cond_F ~ 863 at the
paper's own angles) and the 9-parameter ansatz's M is exactly rank-6-of-9 (a 2-qubit pure
state manifold has only 2d-2=6 real dimensions) -- both facts are load-bearing here.

Run with the physics environment:

    ~/micromamba/envs/qml-julia/bin/python adaptive_qite.py --target all

Every threshold/constant the Route-3 controller uses is invented for this script (the memo
gives no numbers for Route 3) and is collected in ASSUMPTIONS, printed at the end of a full
run. See `../Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation.md` for the prose
explanation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pennylane as qml
import scipy.optimize

SEED = 0
# F7/M5: each target gets its own independent rng, so theta_0 no longer depends on which
# other targets ran before it in the same invocation (previously a single module-level `rng`
# was consumed in sequence, so e.g. `--target d` alone drew a different theta_0 than `--target
# d` inside `--target all`, which followed A/B/C's draws first).
SEED_A, SEED_B, SEED_C, SEED_D = SEED, SEED + 1, SEED + 2, SEED + 3

HERE = Path(__file__).resolve().parent
FIGURES = HERE / "figures"

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------

# H2 in the minimal STO-3G basis, Bravyi-Kitaev -- same system and gate assignment as the
# QDP replication, so results are directly comparable.
G = (0.2252, 0.3435, -0.4347, 0.5716, 0.0910, 0.0910)

H_DERIV = 1e-6          # central-difference step for drho, same as the QDP replication
COLORS = ("#0173B2", "#DE8F05", "#029E73", "#CC78BC", "#D55E00")

# Route-1/2 machinery (embedded stages, RKF, step-doubling) is out of scope: Route 3 uses a
# single f-evaluation per step plus two cheap extra observables (E, Var(H)).

# -- Route 3 controller thresholds (all invented here; the memo gives no numbers) ------
RCOND_DEFAULT = 1e-6        # pseudoinverse cutoff for the 5-parameter (full-rank) ansatz
DELTA_TAU0 = 0.05           # initial step
DELTA_TAU_MIN = 1e-4
DELTA_TAU_MAX = 0.5
GROW = 1.25                 # step growth factor on a clean accept
SHRINK = 0.5                # step shrink factor on rejection / high residual
STIFFNESS_CAP_COEFF = 2000.0  # proactive cap: delta_tau <= COEFF / cond(M); chosen so the
                               # cap only bites as cond(M) approaches ~1/RCOND_DEFAULT, where
                               # the pseudoinverse solve itself starts being regularized away
ENERGY_MONOTONE_TOL = 1e-9  # allowed dE > 0 slack (floating point, not a real violation)
DE_BAND_ABS = 1e-3          # |dE_actual - dE_pred| accept band, absolute part
DE_BAND_REL = 0.5           # ... relative part, fraction of |dE_pred|
RESIDUAL_REL_TOL = 0.05     # ||M theta_dot - Y|| / ||Y|| above this freezes growth
MAX_STEPS = 500
MAX_RETRIES = 20

TAU_MAX = 3.0

ASSUMPTIONS: list[str] = [
    "The memo gives no numbers for Route 3's thresholds. RCOND_DEFAULT, DELTA_TAU0/MIN/MAX, "
    "GROW, SHRINK, STIFFNESS_CAP_COEFF, ENERGY_MONOTONE_TOL, DE_BAND_ABS/REL and "
    "RESIDUAL_REL_TOL above are all our choices, tuned only by inspection of target A's "
    "trajectory (never against a paper-reported value -- there isn't one).",
    "Combination rule for the three Route-3 indicators (our synthesis, not in the memo): "
    "cond(M) is applied as a PROACTIVE cap on delta_tau before the trial step (memo: 'cap "
    "delta_tau'); the energy-monotonicity/variance-identity check is a REACTIVE accept/"
    "reject/retry gate, since it is the only one of the three that actually depends on "
    "delta_tau (memo: 'accept the step iff ... reject/shrink on violation'); the McLachlan "
    "residual is a REACTIVE-BUT-NOT-RETRIED caution signal -- since r_n = ||M theta_dot - Y|| "
    "does not depend on delta_tau at all (it is fixed once theta_n is), retrying with a "
    "smaller step cannot change it, so a high relative residual freezes growth for the next "
    "step instead of rejecting the current one (memo: 'a spike signals ... delta_tau should "
    "shrink', read here as 'should not grow').",
    "Eq. (7) debiasing is applied to M and Y individually before the pseudoinverse solve "
    "(M <- M/(1-p)^2N, Y <- Y/(1-p)^N), matching Algorithm 1's literal step order 'measure "
    "M~, Y~; de-bias by (7); solve for k1' and the closed forms M^err=(1-p)^2N M, "
    "Y^err=(1-p)^N Y already validated in the QDP replication (Eqs. B12-B13) -- algebraically "
    "identical to the memo's own f~ = (1-p)^-N f statement when M is full rank.",
    "Eq. (9)'s literal noise floor (delta_tau_floor = (kappa_b eps_res ||theta_dot||/C_tr)^1/p) "
    "is derived for the embedded-pair estimator e_n of Routes 1/2 and is NOT computed here: "
    "Route 3 never forms e_n. Target D builds our own analogous diagnostic instead (energy-"
    "trajectory deviation from the noiseless run, debiased vs biased, vs p) and is explicitly "
    "not a reproduction of Eq. (9).",
    "M, Y are built from finite-difference density-matrix derivatives (drho, h=1e-6) on "
    "default.qubit/default.mixed, exactly the method already used and validated (against the "
    "closed form Eq. B9, to 1e-10) in the QDP replication -- not Hadamard-test circuits. The "
    "'evaluation count' reported for target A therefore counts calls to this shared "
    "measure_f(), not physical circuit shots.",
    "Initial theta_0 is rng.uniform(-0.5, 0.5, size=N_PARAMS) at SEED=0, not the QDP "
    "replication's published Appendix-C angles -- those were chosen there to reproduce a "
    "specific reported cond(M), not as a sensible QITE starting point. rho_0 = |00>, matching "
    "the QDP replication's convention (also unstated by the source papers). Each target uses "
    "its own default_rng(SEED_{A,B,C,D}=0,1,2,3) (F7/M5), so a single target run reproduces "
    "identically whether invoked alone or as part of --target all; targets do NOT share "
    "theta_0 with each other, so any agreement between two targets' numbers is a statement "
    "about the physics, not a shared random draw.",
    "No shot noise is modelled; the global depolarizing channel is applied to the exact "
    "density matrix (deterministic), so 'noisy' here means biased-but-deterministic, not "
    "stochastic. Route 4 (statistical gating) is explicitly out of scope for this task, and "
    "shot noise is the mechanism that would make it necessary.",
    "F8/M5 policy: on exhausting max_retries or reaching delta_tau_min, the trial step is "
    "committed to the trajectory even though the gate still rejects it -- there is nothing "
    "better to do at the minimum step size. log['accepted'] records this (False), so it "
    "stays usable as a quality signal; a run report should read 'forced through unaccepted' "
    "counts, not just 'steps needed a retry'.",
]

# --------------------------------------------------------------------------------------
# Shared layer -- H2 Hamiltonian, ansatz, M/Y/cond(M) (adapted from the QDP replication)
# --------------------------------------------------------------------------------------

I2 = np.eye(2, dtype=complex)
PX = np.array([[0, 1], [1, 0]], dtype=complex)
PY = np.array([[0, -1j], [1j, 0]], dtype=complex)
PZ = np.diag([1, -1]).astype(complex)
PAULIS = (I2, PX, PY, PZ)


def h2_hamiltonian():
    """H = g0 I + g1 Z0 + g2 Z1 + g3 Z0Z1 + g4 Y0Y1 + g5 X0X1."""
    g0, g1, g2, g3, g4, g5 = G
    kron = np.kron
    return (g0 * np.eye(4)
            + g1 * kron(PZ, I2)
            + g2 * kron(I2, PZ)
            + g3 * kron(PZ, PZ)
            + g4 * kron(PY, PY)
            + g5 * kron(PX, PX))


def exact_ground_energy(H):
    return float(np.linalg.eigvalsh(H)[0])


def variational_floor(rho_fn, n_params, H, n_restarts=8, seed=12345):
    """F9/M6: the lowest E(theta) this ansatz can reach, via direct classical minimization --
    the reference line Target A's convergence plot is missing. Without it, a small gap between
    two runs' final energies cannot be told apart from both runs converging onto the same
    ansatz-expressibility wall from slightly different distances.
    """
    def energy(theta):
        rho = rho_fn(theta)
        return float(np.trace(rho @ H).real)

    rng = np.random.default_rng(seed)
    best = np.inf
    for _ in range(n_restarts):
        theta0 = rng.uniform(-np.pi, np.pi, size=n_params)
        res = scipy.optimize.minimize(energy, theta0, method="Nelder-Mead",
                                       options=dict(xatol=1e-10, fatol=1e-12, maxiter=5000))
        best = min(best, res.fun)
    return float(best)


def gate_list(theta, n_params=5):
    """Fig. 4 ansaetze of the source project's QDP paper, same convention as the replication.

    5 params: RX(t1) q0, RX(t2) q1, RY(t3) q0, RY(t4) q1, CRZ(t5) q0->q1
    9 params: RY(t1) q0, RY(t2) q1, RZ(t3) q0, RZ(t4) q1, CRX(t5) q0->q1,
              RY(t6) q0, RY(t7) q1, RZ(t8) q0, RZ(t9) q1
    """
    if n_params == 5:
        return [("RX", theta[0], 0), ("RX", theta[1], 1),
                ("RY", theta[2], 0), ("RY", theta[3], 1),
                ("CRZ", theta[4], [0, 1])]
    if n_params == 9:
        return [("RY", theta[0], 0), ("RY", theta[1], 1),
                ("RZ", theta[2], 0), ("RZ", theta[3], 1),
                ("CRX", theta[4], [0, 1]),
                ("RY", theta[5], 0), ("RY", theta[6], 1),
                ("RZ", theta[7], 0), ("RZ", theta[8], 1)]
    raise ValueError(f"no ansatz with {n_params} parameters")


def prepare(init_state):
    for wire, bit in enumerate(init_state):
        if bit == "1":
            qml.PauliX(wires=wire)


def make_rho_fn(n_params=5, init_state="00", stop=None):
    """rho(theta), noiseless. `stop` truncates the circuit after that many gates (F11/M5:
    restored from the QDP replication so check_eq_b9_self() can run against this script's own
    metric_M/drho, not just the replication's).
    """
    dev = qml.device("default.qubit", wires=2)

    @qml.qnode(dev)
    def circuit(theta):
        prepare(init_state)
        for name, angle, wires in gate_list(theta, n_params)[:stop]:
            getattr(qml, name)(angle, wires=wires)
        return qml.density_matrix(wires=[0, 1])

    return lambda theta: np.asarray(circuit(np.asarray(theta, dtype=float)), dtype=complex)


def global_depolarizing_kraus(p):
    """E(rho) = (1-p) rho + p I/d, Pauli-twirled Kraus operators, d=4."""
    d = 4
    import itertools
    kraus = [np.sqrt(1.0 - p) * np.eye(d, dtype=complex)]
    for pa, pb in itertools.product(PAULIS, repeat=2):
        kraus.append(np.sqrt(p) / d * np.kron(pa, pb))
    return kraus


def make_noisy_rho_fn(p, n_params=5, init_state="00"):
    """rho^err(theta): the ansatz with a global depolarizing channel after every gate."""
    dev = qml.device("default.mixed", wires=2)
    kraus = global_depolarizing_kraus(p)

    @qml.qnode(dev)
    def circuit(theta):
        prepare(init_state)
        for name, angle, wires in gate_list(theta, n_params):
            getattr(qml, name)(angle, wires=wires)
            qml.QubitChannel(kraus, wires=[0, 1])
        return qml.density_matrix(wires=[0, 1])

    return lambda theta: np.asarray(circuit(np.asarray(theta, dtype=float)), dtype=complex)


def drho(rho_fn, theta, h=H_DERIV):
    theta = np.asarray(theta, dtype=float)
    out = []
    for k in range(len(theta)):
        step = np.zeros_like(theta)
        step[k] = h
        out.append((rho_fn(theta + step) - rho_fn(theta - step)) / (2 * h))
    return np.array(out)


def metric_M(d_rho):
    """M_{k,q} = Tr[(d rho/d theta_k)^dagger d rho/d theta_q]."""
    n = len(d_rho)
    M = np.empty((n, n))
    for k in range(n):
        for q in range(n):
            M[k, q] = np.trace(d_rho[k].conj().T @ d_rho[q]).real
    return M


def vector_Y(d_rho, rho, H):
    """Y_k = -Tr[(d rho/d theta_k) {H, rho}].

    F12/M5: the full McLachlan Y_k (Technical_Resume Sec. 2) is
    -Tr[(d_k rho)^dag ({H,rho} - 2<H>rho)], with two simplifications made here:
    (1) the dagger is dropped -- a no-op because d_k rho is Hermitian for this
        (unitary-gate) parametrization;
    (2) the -2<H>rho term is dropped -- exact only because Tr[d_k rho . rho] = 0, which holds
        both for a pure state (Tr[rho^2]=1 is theta-independent) and, separately, under this
        global depolarizing model (Tr[d_k rho] = 0, so the term vanishes even for the mixed
        rho^err). Neither condition is guaranteed for a general noise channel or a
        non-normalized ansatz -- both are on R_Future_Routes' roadmap -- so this term would
        need to be restored there.
    """
    anti = H @ rho + rho @ H
    return np.array([-np.trace(d_rho[k] @ anti).real for k in range(len(d_rho))])


def check_eq_b9_self(verbose=True):
    """F11/M5: Eq. (B9) internal-consistency check, run against THIS script's own
    metric_M/drho/make_rho_fn -- not just the QDP replication's copy, which the note's
    Assumption 5 previously (and only) cited. For a single-qubit rotation gate with generator
    P_k, M_{k,k} = (1 - Tr[P_k rho_before P_k rho_before]) / 2, exact at p=0.
    """
    theta = np.full(5, 0.3)
    rho_fn = make_rho_fn(n_params=5)
    d = drho(rho_fn, theta)
    M = metric_M(d)
    generators = {0: (PX, 0), 1: (PX, 1), 2: (PY, 0), 3: (PY, 1)}
    worst = 0.0
    for k, (pauli, wire) in generators.items():
        rho_before = make_rho_fn(n_params=5, stop=k)(theta)
        P = np.kron(pauli, I2) if wire == 0 else np.kron(I2, pauli)
        predicted = 0.5 * (1.0 - np.trace(P @ rho_before @ P @ rho_before).real)
        worst = max(worst, abs(M[k, k] - predicted))
        if verbose:
            print(f"    M[{k},{k}] = {M[k, k]:.12f}   Eq.(B9) = {predicted:.12f}   "
                  f"|diff| = {abs(M[k, k] - predicted):.2e}")
    assert worst < 1e-8, f"Eq. (B9) violated by {worst:.2e}: the M convention is wrong"
    return worst


def frobenius_cond(M, rcond=None):
    """cond(M) = ||M||_F ||M^-1||_F.

    With rcond given, uses the SAME regularized pseudoinverse the solve actually uses, so this
    is always finite (it reports how ill-conditioned the regularized solve is, which is what
    Route 3's stiffness cap should react to) -- e.g. a truly rank-deficient M with a generous
    rcond gives a moderate, well-behaved cond(M), not infinity. Only the unregularized case
    (rcond=None, an exact inverse) can be genuinely singular, reported as inf.
    """
    norm = np.linalg.norm(M)
    if rcond is None:
        sv = np.linalg.svd(M, compute_uv=False)
        if sv[-1] <= 1e-12 * sv[0]:
            return np.inf
        inv = np.linalg.inv(M)
    else:
        inv = np.linalg.pinv(M, rcond=rcond)
    return norm * np.linalg.norm(inv)


def energy_and_variance(rho, H):
    """E = Tr[rho H], Var_rho(H) = Tr[rho H^2] - E^2 -- the RHS of Eq. (13), dE/dtau=-2Var."""
    E = float(np.trace(rho @ H).real)
    E2 = float(np.trace(rho @ H @ H).real)
    return E, max(E2 - E * E, 0.0)


# --------------------------------------------------------------------------------------
# Stage 1 (Sec. 1) / Stage 2 (Sec. 2) -- the ODE and the debiased noisy measurement
# --------------------------------------------------------------------------------------

def measure_f(theta, rho_fn, H, N, p=0.0, rcond=RCOND_DEFAULT, debias=True):
    """One evaluation of f(theta) = M^+(theta) Y(theta), Eq. (2).

    p, N, debias implement Sec. 2 / Eq. (7). Eq. (7)'s bias factors (1-p)^2N, (1-p)^N assume
    noise on the derivative line only. This script applies the depolarizing channel to the
    whole circuit (make_noisy_rho_fn), so rho^err is mixed and enters both d(rho) and the
    {H, rho} anticommutator in vector_Y. Working that through with q := (1-p)^N:
    M^err = q^2 M (Eq. 7's M factor is exact here too), but
    Y^err = [q^2 + q(1-q)/2] Y_clean, not q Y_clean -- this is the QDP replication's own
    "Theorem 2 is a factor of two pessimistic" result (project-wide audit, F1/M1), applied to
    debias M and Y before the pseudoinverse solve, matching Algorithm 1's step order. p=0
    makes this a no-op, so the same function serves the noiseless case (Stage 1) and the
    noisy case (Stage 2).
    """
    rho = rho_fn(theta)
    d = drho(rho_fn, theta)
    M = metric_M(d)
    Y = vector_Y(d, rho, H)
    if p > 0.0 and debias:
        q = (1.0 - p) ** N
        M = M / q**2
        Y = Y / (q**2 + q * (1.0 - q) / 2.0)
    theta_dot = np.linalg.pinv(M, rcond=rcond) @ Y
    return theta_dot, M, Y, rho


# --------------------------------------------------------------------------------------
# Stage 3 -- Route 3: McLachlan residual, energy monotonicity, cond(M) stiffness proxy
# --------------------------------------------------------------------------------------

def route3_indicators(M, Y, theta_dot_val, rcond):
    """The two static (delta_tau-independent) Route-3 signals, Sec. 3.3 points 1 and 3.

    Also returns sigma_min and the unregularized spectral cond (F10/M4): frobenius_cond is
    computed *through* the pseudoinverse cutoff, so with rcond fixed it is bounded above by
    construction and cannot by itself distinguish M becoming singular from the cutoff simply
    saturating. sigma_min and cond_spectral are not run through rcond at all.
    """
    residual = float(np.linalg.norm(M @ theta_dot_val - Y))
    cond_M = frobenius_cond(M, rcond=rcond)
    sv = np.linalg.svd(M, compute_uv=False)
    sigma_min = float(sv[-1])
    cond_spectral = float(sv[0] / sv[-1]) if sv[-1] > 0 else np.inf
    return residual, cond_M, sigma_min, cond_spectral


def energy_check(E_n, var_n, Y, theta_dot_val, delta_tau, E_trial):
    """Sec. 3.3 point 2, corrected (F5/M4): gate on the *variational* prediction
    dE_var = -Y.theta_dot delta_tau -- the first-order model of the step the ansatz is
    actually taking -- rather than the exact-ITE identity dE_exact = -2 Var_rho(H) delta_tau
    (Eq. 13). The two differ by the ansatz's expressibility deficit, which does not vanish at
    convergence, so gating on Eq. (13) fights the controller's own goal there (see M4 note).
    Eq. (13)'s value is still returned as a diagnostic, not as a gate input.
    """
    dE_var_pred = -float(Y @ theta_dot_val) * delta_tau
    dE_exact_pred = -2.0 * var_n * delta_tau
    dE_actual = E_trial - E_n
    band = max(DE_BAND_ABS, DE_BAND_REL * abs(dE_var_pred))
    reasons = []
    if dE_actual > ENERGY_MONOTONE_TOL:
        reasons.append("energy increased")
    if abs(dE_actual - dE_var_pred) > band:
        reasons.append("dE off linear prediction")
    return dE_var_pred, dE_exact_pred, dE_actual, reasons


def adaptive_qite_run(theta0, rho_fn, H, n_params, N, tau_max=TAU_MAX,
                       delta_tau0=DELTA_TAU0, p=0.0, debias=True, rcond=RCOND_DEFAULT,
                       stiffness_cap_coeff=STIFFNESS_CAP_COEFF,
                       delta_tau_min=DELTA_TAU_MIN, delta_tau_max=DELTA_TAU_MAX,
                       grow=GROW, shrink=SHRINK, max_steps=MAX_STEPS,
                       max_retries=MAX_RETRIES):
    """The concrete Route-3 algorithm: Stage 1 Euler update, gated by Stage 3's three
    indicators, fed by Stage 2's (optionally noisy, debiased) measurement of f.
    """
    theta = np.array(theta0, dtype=float)
    tau = 0.0
    delta_tau = delta_tau0
    E0, _ = energy_and_variance(rho_fn(theta), H)

    # F9/M6: log["n_evals"] only counts measure_f calls, which undercounts the true quantum
    # cost -- each measure_f call is itself 2N+1 circuit (density-matrix) evaluations via
    # central differences, and each retry-loop trial pays one more (rho_fn(theta_trial) for
    # energy_and_variance). log["n_circ"] and log["s_eff"] track the honest count, in the
    # memo's Eq. (10) currency: C(delta_tau) = s_eff * n_circ(N) * n_shots / delta_tau, with
    # n_shots symbolic (no shot-noise model exists yet, see R5).
    log = dict(tau=[0.0], E=[E0], delta_tau=[delta_tau0], cond_M=[np.nan],
               sigma_min=[np.nan], cond_spectral=[np.nan],
               residual=[np.nan], dE_pred=[np.nan], dE_exact_pred=[np.nan],
               dE_actual=[np.nan], accepted=[True], retries=[0], reasons=[[]],
               cap_bound=[False], s_eff=[np.nan], n_evals=0, n_circ=0)

    step = 0
    while tau < tau_max and step < max_steps:
        theta_dot_val, M, Y, rho_n = measure_f(theta, rho_fn, H, N, p, rcond, debias)
        log["n_evals"] += 1
        log["n_circ"] += 2 * N + 1
        E_n, var_n = energy_and_variance(rho_n, H)
        residual, cond_M, sigma_min, cond_spectral = route3_indicators(M, Y, theta_dot_val, rcond)
        rel_residual = residual / (np.linalg.norm(Y) + 1e-12)

        # Proactive cond(M) cap (Sec. 3.3 point 3: "cap delta_tau"). Floored at delta_tau_min
        # so a singular M (cond_M = inf, e.g. the rank-deficient 9-parameter ansatz) throttles
        # down to the smallest allowed step rather than locking delta_tau at exactly zero.
        # cap_bound records whether the cap (rather than delta_tau_max or the incoming
        # delta_tau) was the actual binding constraint this step (F10/M4 attribution).
        cap = stiffness_cap_coeff / cond_M if cond_M > 0 else np.inf
        cap_bound = cap < delta_tau and cap < delta_tau_max
        delta_tau = max(delta_tau_min, min(delta_tau, cap, delta_tau_max))

        retries = 0
        reasons_seen = []
        while True:
            theta_trial = theta + delta_tau * theta_dot_val
            E_trial, _ = energy_and_variance(rho_fn(theta_trial), H)
            log["n_circ"] += 1  # F9/M6: this rho_fn(theta_trial) call is a real circuit cost
            dE_pred, dE_exact_pred, dE_actual, reasons = energy_check(
                E_n, var_n, Y, theta_dot_val, delta_tau, E_trial)
            reasons_seen.extend(reasons)
            accepted = not reasons
            # F8/M5 policy, made explicit: on exhausting max_retries or delta_tau_min, the
            # trial step is taken anyway (there is nothing better to do at the minimum step
            # size) -- but accepted stays False, so log["accepted"] remains a real quality
            # signal and log["reasons"] records every reason seen across all attempts, not
            # just the (possibly still-failing) final one.
            if accepted or retries >= max_retries or delta_tau <= delta_tau_min:
                break
            delta_tau = max(delta_tau * shrink, delta_tau_min)
            retries += 1

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
        log["s_eff"].append(1 + retries)  # Eq. (10)'s effective-stages-per-step factor

        if accepted and rel_residual <= RESIDUAL_REL_TOL:
            delta_tau = min(delta_tau * grow, delta_tau_max)
        elif rel_residual > RESIDUAL_REL_TOL:
            delta_tau = max(delta_tau * shrink, delta_tau_min)
        step += 1

    # F4/M3: record why the loop stopped -- "tau_max" (finished) vs "max_steps" (stalled).
    # A run that exits on max_steps did not fail to track the baseline, it failed to finish;
    # target_b/d must not silently interpolate past a stalled run's last point.
    log["exit_reason"] = "tau_max" if tau >= tau_max else "max_steps"

    for key in ("tau", "E", "delta_tau", "cond_M", "sigma_min", "cond_spectral", "residual",
                "dE_pred", "dE_exact_pred", "dE_actual", "retries", "s_eff"):
        log[key] = np.array(log[key], dtype=float)
    log["accepted"] = np.array(log["accepted"], dtype=bool)
    log["cap_bound"] = np.array(log["cap_bound"], dtype=bool)
    return log


def fixed_step_qite_run(theta0, rho_fn, H, n_params, N, tau_max=TAU_MAX,
                         delta_tau=DELTA_TAU0, p=0.0, debias=True, rcond=RCOND_DEFAULT,
                         max_steps=MAX_STEPS):
    """Plain fixed-Delta-tau Euler baseline: Stage 1 only, no Route 3 control."""
    theta = np.array(theta0, dtype=float)
    tau = 0.0
    E0, _ = energy_and_variance(rho_fn(theta), H)
    tau_hist, E_hist = [0.0], [E0]
    n_evals = 0
    n_circ = 0  # F9/M6: same currency as adaptive_qite_run -- 2N+1 for measure_f, +1 for
    # the rho_fn(theta) call this loop makes every step just to record E_hist
    step = 0
    while tau < tau_max and step < max_steps:
        theta_dot_val, _, _, _ = measure_f(theta, rho_fn, H, N, p, rcond, debias)
        n_evals += 1
        n_circ += 2 * N + 1
        theta = theta + delta_tau * theta_dot_val
        tau += delta_tau
        E, _ = energy_and_variance(rho_fn(theta), H)
        n_circ += 1
        tau_hist.append(tau)
        E_hist.append(E)
        step += 1
    return dict(tau=np.array(tau_hist), E=np.array(E_hist), n_evals=n_evals, n_circ=n_circ)


def trajectory_deviation(run, baseline):
    """F4/M3: RMS energy deviation of `run` from `baseline`, on their COMMON support only --
    never np.interp past a run's own tau range, which silently holds the last E flat and
    scores that extrapolation as if it were tracking data. Also reports coverage and, if the
    run stalled (exit_reason == "max_steps"), flags that explicitly as a distinct outcome
    from "tracked poorly".
    """
    tau_end = min(run["tau"][-1], baseline["tau"][-1])
    mask = baseline["tau"] <= tau_end
    E_interp = np.interp(baseline["tau"][mask], run["tau"], run["E"])
    rms = float(np.sqrt(np.mean((E_interp - baseline["E"][mask]) ** 2)))
    coverage = float(run["tau"][-1] / baseline["tau"][-1])
    stalled = run.get("exit_reason") == "max_steps"
    return dict(rms_common=rms, coverage=coverage, tau_common_end=tau_end, stalled=stalled)


def banner(text):
    print()
    print("=" * 88)
    print(text)
    print("=" * 88)


# --------------------------------------------------------------------------------------
# Target A -- Stage 1: noiseless ODE + Route 3 controller vs. a fixed-step baseline
# --------------------------------------------------------------------------------------

def target_a_noiseless_convergence(plots=True):
    banner("A  Stage 1 (Sec. 1) -- noiseless adaptive Euler QITE, Route 3 controller")

    H = h2_hamiltonian()
    E_exact = exact_ground_energy(H)
    n_params, N = 5, 5
    theta0 = np.random.default_rng(SEED_A).uniform(-0.5, 0.5, size=n_params)
    rho_fn = make_rho_fn(n_params=n_params)

    adaptive = adaptive_qite_run(theta0, rho_fn, H, n_params, N)
    fixed = fixed_step_qite_run(theta0, rho_fn, H, n_params, N, delta_tau=DELTA_TAU0)

    # F8/M5: distinguish "retried and then satisfied the gate" from "forced through anyway"
    # -- adaptive_qite_run always commits theta_trial, even when the retry budget/delta_tau_min
    # is exhausted with accepted=False still set. Only the second is a real cause for concern.
    retried_ok = int(np.sum((adaptive["retries"][1:] > 0) & adaptive["accepted"][1:]))
    forced_through = int(np.sum(~adaptive["accepted"][1:]))
    print(f"  Exact H2 ground energy (eigh):        {E_exact:.6f}")
    print(f"  Adaptive: {len(adaptive['tau']) - 1} steps, {adaptive['n_evals']} f-evals, "
          f"{retried_ok} retried-then-accepted, {forced_through} forced through unaccepted, "
          f"final tau = {adaptive['tau'][-1]:.4f}")
    print(f"    final E = {adaptive['E'][-1]:.6f}  (error {adaptive['E'][-1] - E_exact:+.3e})")
    print(f"    delta_tau range: [{adaptive['delta_tau'].min():.4g}, "
          f"{adaptive['delta_tau'].max():.4g}]")

    # F6/M4: reason histogram -- "energy increased" is a real Eq. (13) monotonicity
    # violation; "dE off linear prediction" is the (now variational) band test. Both were
    # previously folded into a single "zero violations" claim.
    all_reasons = [r for step_reasons in adaptive["reasons"] for r in step_reasons]
    n_increase = all_reasons.count("energy increased")
    n_band = all_reasons.count("dE off linear prediction")
    print(f"    reason histogram (across all retry attempts): "
          f"energy increased={n_increase}, dE off linear prediction={n_band}")

    # F10/M4: attribute the delta_tau collapse to the cond(M) cap vs. the energy band,
    # using the last quarter of the run (where Target A's note reports the collapse).
    tail = max(1, len(adaptive["cap_bound"]) // 4)
    cap_hits_tail = int(np.sum(adaptive["cap_bound"][-tail:]))
    retry_steps_tail = int(np.sum(adaptive["retries"][-tail:] > 0))
    print(f"    last {tail} steps: cond(M) cap binding on {cap_hits_tail}, "
          f"energy-band retry triggered on {retry_steps_tail}")
    print(f"    sigma_min range: [{np.nanmin(adaptive['sigma_min']):.4g}, "
          f"{np.nanmax(adaptive['sigma_min']):.4g}]  "
          f"cond_spectral (unregularized) range: "
          f"[{np.nanmin(adaptive['cond_spectral']):.4g}, "
          f"{np.nanmax(adaptive['cond_spectral']):.4g}]")

    print(f"  Fixed-step (delta_tau={DELTA_TAU0}): {len(fixed['tau']) - 1} steps, "
          f"{fixed['n_evals']} f-evals, final tau = {fixed['tau'][-1]:.4f}")
    print(f"    final E = {fixed['E'][-1]:.6f}  (error {fixed['E'][-1] - E_exact:+.3e})")

    # F9/M6: the honest currency is circuit count, not f-eval count -- adaptive pays
    # 1+retries trial circuits per step on top of each measure_f's 2N+1, which n_evals never
    # counted. Report both the (biased) f-eval comparison, kept for continuity with earlier
    # runs of this note, and the corrected circuit-count comparison.
    n_match_evals = min(adaptive["n_evals"], fixed["n_evals"])
    E_adapt_evals = (adaptive["E"][n_match_evals] if n_match_evals < len(adaptive["E"])
                      else adaptive["E"][-1])
    E_fixed_evals = (fixed["E"][n_match_evals] if n_match_evals < len(fixed["E"])
                      else fixed["E"][-1])
    print(f"\n  [uncorrected] At matched f-eval budget ({n_match_evals}): "
          f"adaptive E = {E_adapt_evals:.6f}, fixed E = {E_fixed_evals:.6f}")
    print(f"  Total circuit count (F9/M6): adaptive={adaptive['n_circ']}, "
          f"fixed={fixed['n_circ']}  (adaptive pays "
          f"{adaptive['n_circ'] - (2 * n_params + 1) * adaptive['n_evals']} extra circuits "
          f"across all retries, never counted by n_evals)")

    # Corrected comparison: walk both step-by-step circuit-count histories and read off E at
    # the smaller run's total circuit count.
    adapt_circ_cum = np.cumsum((2 * n_params + 1) + 1 + adaptive["retries"][1:])
    fixed_circ_cum = np.cumsum(np.full(len(fixed["tau"]) - 1, 2 * n_params + 2))
    n_match_circ = min(int(adapt_circ_cum[-1]), int(fixed_circ_cum[-1]))
    i_adapt = int(np.searchsorted(adapt_circ_cum, n_match_circ))
    i_fixed = int(np.searchsorted(fixed_circ_cum, n_match_circ))
    E_adapt_circ = adaptive["E"][min(i_adapt + 1, len(adaptive["E"]) - 1)]
    E_fixed_circ = fixed["E"][min(i_fixed + 1, len(fixed["E"]) - 1)]
    print(f"  [corrected] At matched circuit budget ({n_match_circ} circuits): "
          f"adaptive E = {E_adapt_circ:.6f}, fixed E = {E_fixed_circ:.6f}")

    floor = variational_floor(rho_fn, n_params, H)
    print(f"\n  Variational floor of this 5-parameter ansatz (F9/M6, scipy.optimize.minimize, "
          f"multi-restart): E_floor = {floor:.6f}  (E_floor - E_exact = {floor - E_exact:+.3e})")
    print(f"    adaptive final E - floor = {adaptive['E'][-1] - floor:+.3e}")
    print(f"    fixed final E - floor    = {fixed['E'][-1] - floor:+.3e}")
    if floor > min(adaptive["E"][-1], fixed["E"][-1]) + 1e-6:
        print("    WARNING: a run landed below the reported floor -- minimizer likely stuck, "
              "add restarts.")

    if plots:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
        ax = axes[0]
        ax.plot(adaptive["tau"], adaptive["E"] - E_exact, "-", color=COLORS[0], lw=1.8,
                label="adaptive (Route 3)")
        ax.plot(fixed["tau"], fixed["E"] - E_exact, "--", color=COLORS[1], lw=1.6,
                label=f"fixed step ($\\Delta\\tau$={DELTA_TAU0})")
        ax.axhline(floor - E_exact, color="gray", ls=":", lw=1.3,
                   label="ansatz variational floor (F9/M6)")
        ax.set_yscale("log")
        ax.set_xlabel(r"imaginary time $\tau$")
        ax.set_ylabel(r"$E(\tau) - E_{exact}$")
        ax.set_title("A: energy convergence")
        ax.legend(frameon=True)
        ax.grid(alpha=0.3)

        ax = axes[1]
        ax.plot(adaptive["tau"][1:], adaptive["delta_tau"][1:], "-o", color=COLORS[0],
                ms=3, lw=1.2)
        ax.set_xlabel(r"imaginary time $\tau$")
        ax.set_ylabel(r"$\Delta\tau$")
        ax.set_title("A: adaptive step-size trajectory")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(FIGURES / "target_a_convergence.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'target_a_convergence.png'}")

    return dict(adaptive=adaptive, fixed=fixed, E_exact=E_exact)


# --------------------------------------------------------------------------------------
# Target B -- Stage 2: global depolarizing noise, with and without Eq. (7) debiasing
# --------------------------------------------------------------------------------------

def target_b_per_step_check(theta0, rho_fn_clean, H, n_params, N, p_values):
    """F2/F3, M2: the per-step check the note has always reported a number for, given a code
    path here. At fixed theta0, for each p: the flow-MAGNITUDE error (not "direction" -- see
    below) biased and debiased against the clean flow, the collinearity residual (predicted
    ~0, since global depolarizing noise only rescales theta_dot, it never rotates it -- F2),
    and the measured rescale c_meas against the closed form c(p) = 1 + (1-q)/(2q), q=(1-p)^N.
    """
    print("  Per-step check at fixed theta_0 (F2/F3, M2) -- 'flow-MAGNITUDE error', not")
    print("  'flow-direction error': global depolarizing noise only rescales theta_dot, it")
    print("  never rotates it (collinearity residual below, predicted ~0).\n")
    theta_dot_clean, _, _, _ = measure_f(theta0, rho_fn_clean, H, N, p=0.0)
    norm_clean = np.linalg.norm(theta_dot_clean)

    print(f"  {'p':>8} {'err debiased':>13} {'err biased':>12} {'collin.resid':>14} "
          f"{'c_meas':>10} {'c(p)':>10}")
    for p in p_values:
        rho_fn_noisy = make_noisy_rho_fn(p, n_params=n_params)
        theta_dot_biased, _, _, _ = measure_f(theta0, rho_fn_noisy, H, N, p=p, debias=False)
        theta_dot_debiased, _, _, _ = measure_f(theta0, rho_fn_noisy, H, N, p=p, debias=True)

        err_debiased = np.linalg.norm(theta_dot_debiased - theta_dot_clean) / norm_clean
        err_biased = np.linalg.norm(theta_dot_biased - theta_dot_clean) / norm_clean
        collin_resid = np.linalg.norm(
            theta_dot_biased / np.linalg.norm(theta_dot_biased) - theta_dot_clean / norm_clean)
        c_meas = float(theta_dot_biased @ theta_dot_clean) / norm_clean**2
        q = (1.0 - p) ** N
        c_pred = 1.0 + (1.0 - q) / (2.0 * q)

        print(f"  {p:8.4f} {err_debiased*100:12.5f}% {err_biased*100:11.5f}% "
              f"{collin_resid:14.3e} {c_meas:10.6f} {c_pred:10.6f}")
    print()


def target_b_noisy_debiasing(plots=True):
    banner("B  Stage 2 (Sec. 2) -- depolarizing noise, Eq. (7) debiasing")

    H = h2_hamiltonian()
    E_exact = exact_ground_energy(H)
    n_params, N = 5, 5
    theta0 = np.random.default_rng(SEED_B).uniform(-0.5, 0.5, size=n_params)
    rho_fn_clean = make_rho_fn(n_params=n_params)

    p_values = (2e-2, 1e-1)

    target_b_per_step_check(theta0, rho_fn_clean, H, n_params, N, p_values)

    baseline = adaptive_qite_run(theta0, rho_fn_clean, H, n_params, N)
    runs = {}
    print(f"  Baseline (p=0) final E = {baseline['E'][-1]:.6f} "
          f"(error {baseline['E'][-1] - E_exact:+.3e})\n")
    print(f"  {'p':>8} {'debias':>8} {'final E':>12} {'error':>12} {'steps':>7} {'evals':>7}")
    for p in p_values:
        rho_fn_noisy = make_noisy_rho_fn(p, n_params=n_params)
        for debias in (True, False):
            run = adaptive_qite_run(theta0, rho_fn_noisy, H, n_params, N, p=p, debias=debias)
            runs[(p, debias)] = run
            print(f"  {p:8.4f} {str(debias):>8} {run['E'][-1]:12.6f} "
                  f"{run['E'][-1] - E_exact:+12.3e} {len(run['tau']) - 1:7d} "
                  f"{run['n_evals']:7d}")

    print("\n  Deviation of each run's trajectory from the noiseless baseline, on their")
    print("  COMMON tau support only (F4/M3 -- no extrapolation past a stalled run):")
    for p in p_values:
        for debias in (True, False):
            run = runs[(p, debias)]
            dev = trajectory_deviation(run, baseline)
            stall_flag = " STALLED (exit on max_steps)" if dev["stalled"] else ""
            print(f"    p={p:<7.4f} debias={str(debias):<5} "
                  f"RMS(common)={dev['rms_common']:.3e}  "
                  f"coverage={dev['coverage']*100:5.1f}%  "
                  f"tau_common_end={dev['tau_common_end']:.3f}{stall_flag}")

    if plots:
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        ax.plot(baseline["tau"], baseline["E"] - E_exact, "-", color="black", lw=2.0,
                label="noiseless")
        styles = {(p_values[0], True): ("-", COLORS[0]), (p_values[0], False): ("--", COLORS[0]),
                  (p_values[1], True): ("-", COLORS[3]), (p_values[1], False): ("--", COLORS[3])}
        for (p, debias), run in runs.items():
            ls, color = styles[(p, debias)]
            label = f"p={p}, {'debiased' if debias else 'biased'}"
            ax.plot(run["tau"], run["E"] - E_exact, ls, color=color, lw=1.5, label=label)
        ax.set_yscale("log")
        ax.set_xlabel(r"imaginary time $\tau$")
        ax.set_ylabel(r"$E(\tau) - E_{exact}$")
        ax.set_title("B: effect of Eq. (7) debiasing under depolarizing noise")
        ax.legend(frameon=True, fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(FIGURES / "target_b_debiasing.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'target_b_debiasing.png'}")

    return dict(baseline=baseline, runs=runs, E_exact=E_exact)


# --------------------------------------------------------------------------------------
# Target C -- Route 3's cond(M) stiffness proxy on the rank-deficient 9-param ansatz
# --------------------------------------------------------------------------------------

def target_c_stiffness_9param(plots=True):
    banner("C  Route 3 point 3 -- cond(M) stiffness proxy, 9-param (rank-6-of-9) ansatz")
    print("  The QDP replication found this ansatz's M is exactly rank 6 of 9 (a 2-qubit")
    print("  pure-state manifold has only 2d-2=6 real dimensions), i.e. cond(M) = infinity")
    print("  by construction unless the pseudoinverse is regularized. This target shows")
    print("  Route 3's stiffness cap responding to that directly.\n")

    H = h2_hamiltonian()
    E_exact = exact_ground_energy(H)
    n_params, N = 9, 9
    theta0 = np.random.default_rng(SEED_C).uniform(-0.5, 0.5, size=n_params)
    rho_fn = make_rho_fn(n_params=n_params)

    print(f"  {'rcond':>10} {'final E':>12} {'error':>12} {'cond(M) range':>28} "
          f"{'delta_tau range':>22} {'steps':>7}")
    results = {}
    for rcond in (1e-8, 1e-2):
        run = adaptive_qite_run(theta0, rho_fn, H, n_params, N, rcond=rcond)
        results[rcond] = run
        cond_vals = run["cond_M"][1:]
        finite = cond_vals[np.isfinite(cond_vals)]
        cond_range = f"[{finite.min():.3g}, {finite.max():.3g}]" if len(finite) else "all inf"
        dt = run["delta_tau"][1:]
        print(f"  {rcond:10.1e} {run['E'][-1]:12.6f} {run['E'][-1] - E_exact:+12.3e} "
              f"{cond_range:>28} [{dt.min():.4g}, {dt.max():.4g}]      {len(run['tau']) - 1:7d}")

    if plots:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
        for rcond, color in zip((1e-8, 1e-2), (COLORS[0], COLORS[3])):
            run = results[rcond]
            cond_vals = np.clip(run["cond_M"][1:], 0, 1e6)
            axes[0].plot(run["tau"][1:], cond_vals, "-", color=color, lw=1.6,
                         label=f"rcond={rcond:.0e}")
            axes[1].plot(run["tau"][1:], run["delta_tau"][1:], "-", color=color, lw=1.6,
                         label=f"rcond={rcond:.0e}")
        axes[0].set_yscale("log")
        axes[0].set_xlabel(r"imaginary time $\tau$")
        axes[0].set_ylabel(r"$\mathrm{cond}_F(M)$ (clipped at $10^6$)")
        axes[0].set_title("C: stiffness proxy, 9-param ansatz")
        axes[0].legend(frameon=True)
        axes[0].grid(alpha=0.3)
        axes[1].set_xlabel(r"imaginary time $\tau$")
        axes[1].set_ylabel(r"$\Delta\tau$")
        axes[1].set_title("C: stiffness cap engaging")
        axes[1].legend(frameon=True)
        axes[1].grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(FIGURES / "target_c_stiffness.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'target_c_stiffness.png'}")

    return dict(results=results, E_exact=E_exact)


# --------------------------------------------------------------------------------------
# Target D -- our extension of Sec. 2's noise-floor idea to Route 3 (NOT Eq. 9)
# --------------------------------------------------------------------------------------

def target_d_indicator_noise_floor(plots=True):
    banner("D  our extension of Sec. 2 to Route 3 -- NOT the memo's Eq. (9)")
    print("  Eq. (9) is defined in terms of the embedded-pair estimator e_n, which Route 3")
    print("  never computes. This target asks the analogous question for Route 3 directly:")
    print("  at what noise level p does even the DEBIASED trajectory stop tracking the")
    print("  noiseless one, and does debiasing actually buy anything at each p?\n")

    H = h2_hamiltonian()
    n_params, N = 5, 5
    theta0 = np.random.default_rng(SEED_D).uniform(-0.5, 0.5, size=n_params)
    rho_fn_clean = make_rho_fn(n_params=n_params)
    baseline = adaptive_qite_run(theta0, rho_fn_clean, H, n_params, N)

    p_grid = (0.0, 1e-4, 1e-3, 1e-2, 5e-2)
    rms_debiased, rms_biased = [], []
    print(f"  {'p':>8} {'RMS debiased':>14} {'cov%':>6} {'RMS biased':>14} {'cov%':>6}")
    for p in p_grid:
        if p == 0.0:
            rms_debiased.append(0.0)
            rms_biased.append(0.0)
            print(f"  {p:8.4f} {0.0:14.3e} {100.0:6.1f} {0.0:14.3e} {100.0:6.1f}")
            continue
        rho_fn_noisy = make_noisy_rho_fn(p, n_params=n_params)
        run_d = adaptive_qite_run(theta0, rho_fn_noisy, H, n_params, N, p=p, debias=True)
        run_b = adaptive_qite_run(theta0, rho_fn_noisy, H, n_params, N, p=p, debias=False)
        dev_d = trajectory_deviation(run_d, baseline)
        dev_b = trajectory_deviation(run_b, baseline)
        rms_debiased.append(dev_d["rms_common"])
        rms_biased.append(dev_b["rms_common"])
        flags = (" STALLED(d)" if dev_d["stalled"] else "") + (" STALLED(b)" if dev_b["stalled"] else "")
        print(f"  {p:8.4f} {dev_d['rms_common']:14.3e} {dev_d['coverage']*100:6.1f} "
              f"{dev_b['rms_common']:14.3e} {dev_b['coverage']*100:6.1f}{flags}")

    if plots:
        fig, ax = plt.subplots(figsize=(6.0, 4.4))
        ax.plot(p_grid, rms_debiased, "-o", color=COLORS[0], lw=1.6, label="debiased")
        ax.plot(p_grid, rms_biased, "-s", color=COLORS[4], lw=1.6, label="biased")
        ax.set_xscale("symlog", linthresh=1e-4)
        ax.set_yscale("log")
        ax.set_xlabel("depolarizing probability $p$")
        ax.set_ylabel("RMS energy deviation from noiseless run")
        ax.set_title("D: our Route-3 noise-floor diagnostic (not Eq. 9)")
        ax.legend(frameon=True)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(FIGURES / "target_d_noise_floor.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'target_d_noise_floor.png'}")

    return dict(p_grid=p_grid, rms_debiased=rms_debiased, rms_biased=rms_biased)


# --------------------------------------------------------------------------------------

def print_summary():
    banner("Assumptions (ours, not the memo's)")
    for i, a in enumerate(ASSUMPTIONS, 1):
        print(f"  {i}. {a}")


TARGETS = {
    "a": target_a_noiseless_convergence,
    "b": target_b_noisy_debiasing,
    "c": target_c_stiffness_9param,
    "d": target_d_indicator_noise_floor,
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
