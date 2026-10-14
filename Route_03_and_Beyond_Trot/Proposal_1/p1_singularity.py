#!/usr/bin/env python
"""P1 of `../Plan_Beyond_Trot_Into_Route_03_and_Adam.md`: classify the singularity the cond(M)
cap responds to in Target C / E3, before deciding anything about the cap.

The question, verbatim from the plan: is the reported cond(M) = 9.66e6 along Target C / E3
(9-parameter, rank-6-of-9 ansatz, rcond = 1e-8) the conditioning INSIDE the physical
6-dimensional range, or finite-difference leakage of the three gauge directions past the
pseudoinverse cutoff? The prediction was written into the plan before this script existed.

Nothing in `../../Attack_Route_03/adaptive_qite.py` or `../../Adam_QITE/adam_qite.py` is
modified. Both are imported (the Adam `sys.path` pattern). Per-step (theta, M, Y) is captured by
temporarily swapping each module's global `measure_f` for a pass-through recorder; the replays
reproducing the shipped logs to the last printed digit is what shows the recorder is inert.

Targets, one function each:

  ids     check_identities()          plan sec. 2 point 2: M = 2 Cov(Gtilde_j, Gtilde_k), and
                                      the Y sign -- exact state-form M, Y vs drho/metric_M/vector_Y
  c       p1_target_c_replay()        Target C, both rconds: full spectrum, what pinv keeps,
                                      rank-6 cond vs the logged cond, parameter vs state speed
  null    p1_null_vectors()           right-singular vectors at the cond peak, by gate
  e3      p1_e3_arms()                the four E3 arms at rcond = 1e-8, same diagnostics
  exact   p1_exact_metric_rerun()     Target C re-run on the exact M, Y (no finite differences)
  a       p1_target_a_control()       Target A (5-param, full rank) as the control
  cap     p1_regime_aware_cap()       the plan's follow-on: cap on cond over the top r = 6 only

Run with the physics environment:

    ~/micromamba/envs/qml-julia/bin/python p1_singularity.py --target all 2>&1 \
      | tee Output.txt

`--check` runs everything without figures and exits nonzero if any pinned number drifts. There
is no tolerance flag: floats are compared at PIN_RTOL, a floating-point reproducibility margin
fixed in this file, and integers (step counts, kept-rank counts) must match exactly. A drift is
a result to record, never a reason to edit PINNED until it agrees.

See `P1_Singularity_Classification_Findings.md` for the write-up.
"""

from __future__ import annotations

import argparse
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
QITE = HERE.parent.parent
sys.path.insert(0, str(QITE / "Attack_Route_03"))
sys.path.insert(0, str(QITE / "Adam_QITE"))

import adaptive_qite as aq  # noqa: E402  (also selects matplotlib's Agg backend)
import adam_qite as ad  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pennylane as qml  # noqa: E402

FIGURES = HERE / "figures"
COLORS = aq.COLORS

# -- Frozen before the first run -------------------------------------------------------
SEED = 0                  # only for the random-theta identity check; trajectories use SEED_A/C
N_IDENTITY_THETAS = 4     # random angles per ansatz in check_identities
R_STRUCT = 6              # structural rank of the 9-param ansatz, 2d-2 with d=4 -- a priori
IDENTITY_TOL = 1e-8       # fd truncation+roundoff is ~1e-10; this is the B9 check's own bound
GAUGE_BOUND = 1e-12       # sigma_7/sigma_1 must sit below this for "no leakage" (cutoff is 1e-8)
PIN_RTOL = 1e-9           # float reproducibility margin on a fixed BLAS, not a fit knob
STIFF_COND = 1e5          # a step "sits in the fold" if its logged cond exceeds this; reporting only

ASSUMPTIONS: list[str] = [
    "Every physical choice (H2 Hamiltonian, ansaetze, drho at h=1e-6, the Route-3 gate, the "
    "cap, the seeds, the Adam law) is imported from adaptive_qite.py / adam_qite.py unchanged. "
    "Their own assumption lists apply.",
    "OURS: the exact metric uses the state form of plan sec. 2, "
    "M_jk = 2(Re<Gt_j Gt_k> - <Gt_j><Gt_k>), Y_k = -i<[Gt_k, H]>, with Gt_k = P_k G_k P_k^dag, "
    "P_k the gates after k, G_k = P/2 for RX/RY/RZ and |1><1| (x) P/2 for CRX/CRZ. Gate "
    "matrices come from qml.matrix on adaptive_qite.gate_list, so wire order and sign "
    "conventions are PennyLane's, not re-derived. Valid for p = 0 (pure states) only; every "
    "run in this script is noiseless.",
    "OURS: the structural rank r = 6 is fixed a priori from 2d-2, never fitted to a spectrum.",
    "OURS: rank_r_cond is reported in the Frobenius form ||M||_F sqrt(sum_{i<=r} sigma_i^-2), "
    "which is what adaptive_qite.frobenius_cond returns when pinv keeps exactly r values, so the "
    "comparison with the logged cond_M is like for like. The spectral form sigma_1/sigma_r is "
    "printed alongside.",
    "OURS: 'state speed' is ||d rho/d tau||_F = sqrt(theta_dot^T M theta_dot) (exact for a pure "
    "state, since M is the Gram matrix of the d_k rho). It is the physical speed of the "
    "trajectory, independent of how the ansatz is parametrized.",
    "OURS: the regime-aware cap (target cap) is implemented by swapping adaptive_qite's global "
    "route3_indicators for one that returns rank_r_cond in the cond_M slot. The cap "
    "coefficient (2000) is untouched, as the plan required.",
]

# Pinned after the first full run (Output.txt). Full precision.
PINNED: dict[str, float | int] = {
    'c/1e-08/error': 0.0008868181396866071,
    'c/1e-08/steps': 108,
    'c/1e-08/cond_max': 9663747.736069854,
    'c/pred_gauge_range': 0,
    'c/0.01/error': 0.0021434047708084325,
    'c/0.01/steps': 17,
    'c/0.01/cond_max': 125.38857061439889,
    'c/0.01/dropped': 11,
    'c/peak_tau': 0.19828322964616446,
    'c/peak_speed_theta': 128.43414782028475,
    'c/peak_speed_state': 0.18622156320885339,
    'null/v6_th8': 0.933728059143997,
    'e3/route3_cap_on/error': 0.0008868181396866071,
    'e3/route3_cap_on/steps': 108,
    'e3/route3_cap_on/fold': 86,
    'e3/route3_cap_off/error': 0.0013979350857165063,
    'e3/route3_cap_off/steps': 22,
    'e3/route3_cap_off/fold': 4,
    'e3/adam_cap_on/error': 0.0008868181396866071,
    'e3/adam_cap_on/steps': 108,
    'e3/adam_cap_on/fold': 86,
    'e3/adam_cap_off/error': 0.000565295558648593,
    'e3/adam_cap_off/steps': 28,
    'e3/adam_cap_off/fold': 4,
    'exact/1e-08/error': 0.000885419437947288,
    'exact/1e-08/steps': 108,
    'exact/1e-08/gap': 0.00011006000211699618,
    'exact/0.01/error': 0.002143404750947875,
    'exact/0.01/steps': 17,
    'exact/0.01/gap': 1.2278327798931343e-09,
    'a/error': 0.022605762509368743,
    'a/steps': 172,
    'a/dropped': 12,
    'a/v5_crz': 0.9991635802633979,
    'a/pop_max': 0.9990163198181179,
    'cap/1e-08/error': 0.000886166932451804,
    'cap/1e-08/steps': 108,
    'cap/0.01/error': 0.003405080644245162,
    'cap/0.01/steps': 105,
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
    else:
        ok = abs(measured - pinned) <= PIN_RTOL * abs(pinned)
    print(f"{pad}{label:<44s} {measured!r}   pinned {pinned!r}   {'MATCH' if ok else 'DRIFT'}")
    if not ok:
        DRIFTS.append(f"{key}: measured {measured!r}, pinned {pinned!r}")


def claim(label, ok, indent=4):
    """A pre-registered yes/no statement. Recorded as a drift if it fails."""
    global N_COMPARED
    N_COMPARED += 1
    print(f"{' ' * indent}{label:<70s} {'HOLDS' if ok else 'FAILS'}")
    if not ok:
        DRIFTS.append(f"claim failed: {label}")


def prediction(label, ok, key, indent=4):
    """A pre-registered prediction from the plan. Its outcome is a RESULT: printed as HELD or
    FAILED and pinned like any other number, so a failed prediction does not fail --check,
    but a change in its outcome does."""
    print(f"{' ' * indent}prediction: {label:<58s} {'HELD' if ok else 'FAILED'}")
    compare("  outcome (1 = held)", int(ok), key, indent=indent)


# --------------------------------------------------------------------------------------
# The four P1 functions
# --------------------------------------------------------------------------------------

_P1 = np.diag([0, 1]).astype(complex)
_HALF = {"RX": aq.PX / 2, "RY": aq.PY / 2, "RZ": aq.PZ / 2}


def _gate_and_generator(name, angle, wires):
    """(U, G) with U = exp(-i angle G), both 4x4 on wire order [0, 1]."""
    U = qml.matrix(getattr(qml, name)(angle, wires=wires), wire_order=[0, 1])
    if name in _HALF:
        g = _HALF[name]
        G = np.kron(g, aq.I2) if wires == 0 else np.kron(aq.I2, g)
    else:  # CRX / CRZ, control 0, target 1
        G = np.kron(_P1, _HALF[name[1:]])
    return np.asarray(U, dtype=complex), G


def exact_metric_and_force(theta, n_params, H, init_state="00"):
    """State-form M and Y from the rotated generators (plan sec. 2 point 2), no differences.

    d_k|psi> = -i Gt_k |psi>, so d_k rho = -i[Gt_k, rho] and, for a pure state,
    M_jk = Tr[d_j rho d_k rho] = 2(Re<Gt_j Gt_k> - <Gt_j><Gt_k>),
    Y_k  = -Tr[d_k rho {H, rho}] = -i<[Gt_k, H]>.
    Returns (M, Y, psi).
    """
    data = [_gate_and_generator(*g) for g in aq.gate_list(theta, n_params)]
    psi = np.zeros(4, dtype=complex)
    psi[int(init_state, 2)] = 1.0
    for U, _ in data:
        psi = U @ psi
    rotated = []
    for k, (_, G) in enumerate(data):
        P = np.eye(4, dtype=complex)
        for U, _ in data[k + 1:]:
            P = U @ P
        rotated.append(P @ G @ P.conj().T)
    V = np.array([G @ psi for G in rotated])          # Gt_k |psi>
    mean = np.array([np.vdot(psi, v).real for v in V])
    M = 2.0 * ((V.conj() @ V.T).real - np.outer(mean, mean))
    Y = np.array([(-1j * (np.vdot(psi, G @ H @ psi) - np.vdot(psi, H @ G @ psi))).real
                  for G in rotated])
    return M, Y, psi


def singular_spectrum_log(M):
    """All singular values, largest first -- not just sigma_min."""
    return np.linalg.svd(M, compute_uv=False)


def null_vector_attribution(M, r):
    """Right-singular vectors r-1 .. n-1 (0-based): the softest physical direction (index r-1)
    and the gauge directions beyond rank r, as rows of loadings on the parameters."""
    _, s, Vt = np.linalg.svd(M)
    return s, Vt[r - 1:]


def rank_r_cond(M, r):
    """cond over the top r singular values: (Frobenius form, spectral form)."""
    s = np.linalg.svd(M, compute_uv=False)
    return (float(np.linalg.norm(M) * np.sqrt(np.sum(1.0 / s[:r] ** 2))),
            float(s[0] / s[r - 1]))


# --------------------------------------------------------------------------------------
# Instrumentation: record (theta, M, Y, rcond) at every measure_f call
# --------------------------------------------------------------------------------------

@contextmanager
def recording(replacement=None):
    """Swap aq.measure_f and ad.measure_f for a recorder (optionally around `replacement`)."""
    original_aq, original_ad = aq.measure_f, ad.measure_f
    inner = replacement or original_aq
    log = []

    def recorder(theta, rho_fn, H, N, p=0.0, rcond=aq.RCOND_DEFAULT, debias=True):
        out = inner(theta, rho_fn, H, N, p, rcond, debias)
        log.append((np.array(theta, dtype=float), out[1].copy(), out[2].copy(), rcond))
        return out

    aq.measure_f = ad.measure_f = recorder
    try:
        yield log
    finally:
        aq.measure_f, ad.measure_f = original_aq, original_ad


def exact_measure_f(theta, rho_fn, H, N, p=0.0, rcond=aq.RCOND_DEFAULT, debias=True):
    """Drop-in measure_f on the exact M, Y. p must be 0 (see ASSUMPTIONS)."""
    assert p == 0.0, "exact_measure_f is pure-state only"
    M, Y, psi = exact_metric_and_force(theta, N, H)
    return np.linalg.pinv(M, rcond=rcond) @ Y, M, Y, np.outer(psi, psi.conj())


def per_step(log, run, r=R_STRUCT, with_exact=True, H=None):
    """Everything P1 measures at each recorded step, as arrays."""
    n = len(log)
    out = dict(tau=run["tau"][:n], cond_logged=run["cond_M"][1:n + 1], kept=np.empty(n, int),
               sv=np.empty((n, len(log[0][1]))), cond_r_F=np.empty(n), cond_r_2=np.empty(n),
               speed_theta=np.empty(n), speed_state=np.empty(n), align_soft=np.empty(n),
               fd_err=np.full(n, np.nan))
    for i, (theta, M, Y, rcond) in enumerate(log):
        s = singular_spectrum_log(M)
        out["sv"][i] = s
        out["kept"][i] = int(np.sum(s > rcond * s[0]))
        out["cond_r_F"][i], out["cond_r_2"][i] = rank_r_cond(M, min(r, len(s)))
        theta_dot = np.linalg.pinv(M, rcond=rcond) @ Y
        out["speed_theta"][i] = np.linalg.norm(theta_dot)
        out["speed_state"][i] = np.sqrt(max(theta_dot @ M @ theta_dot, 0.0))
        soft = np.linalg.svd(M)[2][min(r, len(s)) - 1]
        out["align_soft"][i] = (soft @ theta_dot) ** 2 / (theta_dot @ theta_dot)
        if with_exact:
            M_ex, _, _ = exact_metric_and_force(theta, len(theta), H)
            out["fd_err"][i] = np.abs(M - M_ex).max()
    return out


def banner(text):
    aq.banner(text)


def _gate_names(n_params):
    names = []
    for k, (name, _, wires) in enumerate(aq.gate_list(np.zeros(n_params), n_params)):
        where = "q0q1" if isinstance(wires, list) else f"q{wires}"
        names.append(f"th{k + 1} {name} {where}")
    return names


def _c_setup():
    H = aq.h2_hamiltonian()
    theta0 = np.random.default_rng(aq.SEED_C).uniform(-0.5, 0.5, size=9)
    return H, aq.exact_ground_energy(H), theta0, aq.make_rho_fn(n_params=9)


# --------------------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------------------

def check_identities(plots=True):
    banner("ids  plan sec. 2 point 2 -- exact state-form M, Y vs drho / metric_M / vector_Y")
    H = aq.h2_hamiltonian()
    rng = np.random.default_rng(SEED)
    worst_M = worst_Y = worst_Yflip = 0.0
    for n_params in (5, 9):
        rho_fn = aq.make_rho_fn(n_params=n_params)
        for _ in range(N_IDENTITY_THETAS):
            theta = rng.uniform(-np.pi, np.pi, size=n_params)
            d = aq.drho(rho_fn, theta)
            M_fd, Y_fd = aq.metric_M(d), aq.vector_Y(d, rho_fn(theta), H)
            M_ex, Y_ex, _ = exact_metric_and_force(theta, n_params, H)
            worst_M = max(worst_M, np.abs(M_ex - M_fd).max())
            worst_Y = max(worst_Y, np.abs(Y_ex - Y_fd).max())
            worst_Yflip = max(worst_Yflip, np.abs(Y_ex + Y_fd).max())
    print(f"  {2 * N_IDENTITY_THETAS} random theta (5- and 9-param), default_rng({SEED}):")
    print(f"    max |M_exact - M_fd|  = {worst_M:.3e}")
    print(f"    max |Y_exact - Y_fd|  = {worst_Y:.3e}")
    print(f"    max |Y_exact + Y_fd|  = {worst_Yflip:.3e}   (the opposite sign, for contrast)")
    claim(f"M = 2 Cov(Gt_j, Gt_k) matches metric_M(drho) to {IDENTITY_TOL:g}",
          worst_M < IDENTITY_TOL)
    claim(f"Y = -i<[Gt_k, H]> matches vector_Y to {IDENTITY_TOL:g} (sign as derived)",
          worst_Y < IDENTITY_TOL)

    M0, _, _ = exact_metric_and_force(np.zeros(5), 5, H)
    print(f"\n  5-param ansatz at theta = 0: CRZ row of M = {np.array2string(M0[4], precision=3)}")
    claim("CRZ row is exactly zero at theta = 0 (zero-variance mechanism)",
          np.abs(M0[4]).max() == 0.0)

    theta = np.full(5, 0.3)
    M_ex, _, _ = exact_metric_and_force(theta, 5, H)
    worst_b9 = 0.0
    for k, (pauli, wire) in {0: (aq.PX, 0), 1: (aq.PX, 1), 2: (aq.PY, 0), 3: (aq.PY, 1)}.items():
        rho_before = aq.make_rho_fn(n_params=5, stop=k)(theta)
        P = np.kron(pauli, aq.I2) if wire == 0 else np.kron(aq.I2, pauli)
        worst_b9 = max(worst_b9, abs(M_ex[k, k]
                                     - 0.5 * (1 - np.trace(P @ rho_before @ P @ rho_before).real)))
    print(f"  Eq. (B9) diagonal on the exact M at theta = 0.3: worst |diff| = {worst_b9:.3e}")
    claim("exact M satisfies Eq. (B9) to 1e-12", worst_b9 < 1e-12)
    return dict(worst_M=worst_M, worst_Y=worst_Y)


def p1_target_c_replay(plots=True):
    banner("c  Target C replay, rcond in {1e-8, 1e-2}: what does the logged cond(M) measure?")
    print("  Prediction (plan P1, written before this run): sigma_7..9 from drho sit near")
    print("  1e-10..1e-12, below the 1e-8*sigma_1 cutoff, so rank_r_cond (r = 6) equals the")
    print("  logged cond along E3 and the cap is reading conditioning inside the physical range.\n")
    H, E_exact, theta0, rho_fn = _c_setup()
    results = {}
    for rcond in (1e-8, 1e-2):
        with recording() as log:
            run = aq.adaptive_qite_run(theta0, rho_fn, H, 9, 9, rcond=rcond)
        d = per_step(log, run, H=H)
        results[rcond] = (run, d)
        tag = f"c/{rcond:g}"
        print(f"  rcond = {rcond:g}")
        compare("final error (Output: see note)", run["E"][-1] - E_exact,
                f"{tag}/error")
        compare("steps", len(run["tau"]) - 1, f"{tag}/steps")
        compare("max logged cond_M", float(np.nanmax(d["cond_logged"])), f"{tag}/cond_max")
        gauge = float((d["sv"][:, 6] / d["sv"][:, 0]).max())
        print(f"    kept by pinv, per step:                     {sorted(set(d['kept'].tolist()))}")
        print(f"    max sigma_7/sigma_1 over the run:           {gauge:.3e}   "
              f"(cutoff {rcond:g})")
        print(f"    max |M_fd - M_exact| over the run:          {np.nanmax(d['fd_err']):.3e}")
        ratio = d["cond_r_F"] / d["cond_logged"]
        print(f"    rank-6 cond_F / logged cond_M:              [{ratio.min():.15f}, "
              f"{ratio.max():.15f}]")
        print(f"    max rank-6 spectral cond sigma_1/sigma_6:   {d['cond_r_2'].max():.4e}")
        print(f"    steps with logged cond_M > {STIFF_COND:g}:           "
              f"{int(np.sum(d['cond_logged'] > STIFF_COND))} of {len(d['kept'])}")
        claim(f"sigma_7/sigma_1 < {GAUGE_BOUND:g} on every step (no gauge leakage)",
              gauge < GAUGE_BOUND)
        if rcond == 1e-8:
            claim("pinv keeps exactly 6 singular values on every step",
                  set(d["kept"].tolist()) == {R_STRUCT})
            claim("rank-6 cond_F equals the logged cond_M to 1e-12 on every step",
                  np.abs(ratio - 1).max() < 1e-12)
            prediction("sigma_7/sigma_1 lies in [1e-12, 1e-10]", 1e-12 <= gauge <= 1e-10,
                       "c/pred_gauge_range")
        else:
            print(f"    NOTE rcond=1e-2 also discards physical directions: kept drops to "
                  f"{d['kept'].min()} of 6 on {int(np.sum(d['kept'] < R_STRUCT))} steps")
            compare("steps where pinv keeps < 6", int(np.sum(d["kept"] < R_STRUCT)),
                    f"{tag}/dropped")
        print()

    run, d = results[1e-8]
    i = int(np.nanargmax(d["cond_logged"]))
    print("  rcond = 1e-8 around the cond peak: parameter speed vs state speed")
    print(f"    {'tau':>7} {'cond_M':>10} {'|theta_dot|':>12} {'||drho/dtau||':>14} "
          f"{'frac along v6':>14}")
    for j in sorted({0, i - 30, i - 15, i, i + 15, len(d['tau']) - 1}):
        print(f"    {d['tau'][j]:7.4f} {d['cond_logged'][j]:10.3e} {d['speed_theta'][j]:12.4g} "
              f"{d['speed_state'][j]:14.4g} {d['align_soft'][j]:14.8f}")
    compare("peak tau (start of step)", float(d["tau"][i]), "c/peak_tau")
    compare("|theta_dot| at peak", float(d["speed_theta"][i]), "c/peak_speed_theta")
    compare("||drho/dtau|| at peak", float(d["speed_state"][i]), "c/peak_speed_state")
    s_ratio = d["speed_theta"].max() / d["speed_theta"][0]
    st_ratio = d["speed_state"].max() / d["speed_state"].min()
    print(f"    |theta_dot| varies by x{s_ratio:.1f} over the run; ||drho/dtau|| by "
          f"x{st_ratio:.2f} (max/min)")

    if plots:
        fig, axes = plt.subplots(2, 2, figsize=(11, 7.6), sharex=True)
        tau = d["tau"]
        ax = axes[0, 0]
        rel = d["sv"] / d["sv"][:, :1]
        for k in range(5):
            ax.plot(tau, rel[:, k], "-", color="0.65", lw=1.0,
                    label=r"$\sigma_{1..5}$" if k == 0 else None)
        ax.plot(tau, rel[:, 5], "-", color=COLORS[0], lw=2.0, label=r"$\sigma_6$ (softest physical)")
        for k in range(6, 9):
            ax.plot(tau, rel[:, k], "-", color=COLORS[4], lw=1.0,
                    label=r"$\sigma_{7..9}$ (gauge)" if k == 6 else None)
        ax.axhline(1e-8, color="black", ls="--", lw=1.2, label="pinv cutoff (rcond=1e-8)")
        ax.set_yscale("log")
        ax.set_ylabel(r"$\sigma_i/\sigma_1$ of $M_{fd}$")
        ax.set_title("(a) full spectrum along Target C")
        ax.legend(frameon=True, fontsize=8, loc="center right", bbox_to_anchor=(1.0, 0.47))

        ax = axes[0, 1]
        ax.plot(tau, d["cond_logged"], "-", color=COLORS[0], lw=2.0, label="logged cond$_F(M)$")
        ax.plot(tau, d["cond_r_F"], "--", color=COLORS[1], lw=1.6,
                label="rank-6 cond$_F$ (this script)")
        ax.set_yscale("log")
        ax.set_ylabel("condition number")
        ax.set_title("(b) the cap reads the rank-6 conditioning")
        ax.legend(frameon=True, fontsize=8)

        ax = axes[1, 0]
        ax.plot(tau, d["speed_theta"], "-", color=COLORS[2], lw=2.0)
        ax.set_yscale("log")
        ax.set_ylabel(r"$\|\dot\theta\|$")
        ax.set_title("(c) parameter speed")

        ax = axes[1, 1]
        ax.plot(tau, d["speed_state"], "-", color=COLORS[3], lw=2.0)
        ax.set_yscale("log")
        ax.set_ylabel(r"$\|d\rho/d\tau\|_F=\sqrt{\dot\theta^T M\dot\theta}$")
        ax.set_title("(d) state speed")
        for ax in axes.flat:
            ax.set_xscale("log")
            ax.grid(alpha=0.3)
        for ax in axes[1]:
            ax.set_xlabel(r"imaginary time $\tau$ (start of step)")
        fig.tight_layout()
        fig.savefig(FIGURES / "p1_target_c_spectrum.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'p1_target_c_spectrum.png'}")
    return results


def p1_null_vectors(plots=True):
    banner("null  right-singular vectors at the Target C cond peak (rcond = 1e-8)")
    H, _, theta0, rho_fn = _c_setup()
    with recording() as log:
        run = aq.adaptive_qite_run(theta0, rho_fn, H, 9, 9, rcond=1e-8)
    i = int(np.nanargmax(run["cond_M"][1:]))
    theta, M, Y, _ = log[i]
    s, vecs = null_vector_attribution(M, R_STRUCT)
    names = _gate_names(9)
    print(f"  step {i}, tau = {run['tau'][i]:.4f}, sigma = "
          f"{np.array2string(s, precision=3, max_line_width=120)}\n")
    print(f"  {'':<14}" + "".join(f"{'v' + str(R_STRUCT + j):>9}" for j in range(len(vecs))))
    for k, name in enumerate(names):
        print(f"  {name:<14}" + "".join(f"{v[k]:9.3f}" for v in vecs))
    soft = vecs[0]
    top = np.argsort(-np.abs(soft))[:3]
    print(f"\n  softest physical direction v6: top loadings "
          + ", ".join(f"{names[k]} {soft[k]:+.3f}" for k in top))
    compare("v6 loading on th8 (RZ q0, layer 2)", float(abs(soft[7])), "null/v6_th8")

    _, _, psi = exact_metric_and_force(theta, 9, H)
    z0 = float(np.vdot(psi, np.kron(aq.PZ, aq.I2) @ psi).real)
    print(f"  <Z_0> at the peak = {z0:.4f}; M[th8,th8] = {M[7, 7]:.4e}, "
          f"(1 - <Z_0>^2)/2 = {(1 - z0 ** 2) / 2:.4e}")
    print("  (th8 is RZ on q0 followed only by RZ on q1, so Gt_8 = Z_0/2 and M_88 = 2 Var(Z_0/2).)")
    print(f"  sigma_6 / M[th8,th8] = {s[5] / M[7, 7]:.3e}: the soft direction is a near-linear")
    print("  dependence of th8's tangent on the others, not only a weak row.")
    for j, v in enumerate(vecs[1:], start=R_STRUCT + 1):
        top = np.argsort(-np.abs(v))[:2]
        print(f"  gauge v{j}: largest loadings " +
              ", ".join(f"{names[k]} {v[k]:+.3f}" for k in top))

    if plots:
        fig, ax = plt.subplots(figsize=(8.4, 3.6))
        im = ax.imshow(np.abs(vecs), cmap="Blues", vmin=0, vmax=1, aspect="auto")
        for r_ in range(vecs.shape[0]):
            for c_ in range(vecs.shape[1]):
                val = abs(vecs[r_, c_])
                ax.text(c_, r_, f"{val:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if val > 0.55 else "black")
        ax.set_xticks(range(9), [n.replace(" ", "\n", 1) for n in names], fontsize=7)
        ax.set_yticks(range(vecs.shape[0]),
                      [r"$v_6$ soft" if j == 0 else rf"$v_{{{R_STRUCT + j}}}$ gauge"
                       for j in range(vecs.shape[0])])
        ax.set_title(rf"|loadings| at the cond peak, $\tau$={run['tau'][i]:.3f}")
        fig.colorbar(im, ax=ax, fraction=0.03)
        fig.tight_layout()
        fig.savefig(FIGURES / "p1_null_vectors.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'p1_null_vectors.png'}")
    return dict(s=s, vecs=vecs, z0=z0)


def p1_e3_arms(plots=True):
    banner("e3  the four E3 arms at rcond = 1e-8, same diagnostics")
    H, E_exact, theta0, rho_fn = _c_setup()
    arms = {
        "route3_cap_on": lambda: aq.adaptive_qite_run(theta0, rho_fn, H, 9, 9, rcond=1e-8),
        "route3_cap_off": lambda: aq.adaptive_qite_run(theta0, rho_fn, H, 9, 9, rcond=1e-8,
                                                       stiffness_cap_coeff=np.inf),
        "adam_cap_on": lambda: ad.adam_qite_run(theta0, rho_fn, H, 9, 9, rcond=1e-8),
        "adam_cap_off": lambda: ad.adam_qite_run(theta0, rho_fn, H, 9, 9, rcond=1e-8,
                                                 stiffness_cap_coeff=np.inf),
    }
    print(f"  {'arm':<16} {'error':>12} {'steps':>6} {'kept':>6} {'max s7/s1':>10} "
          f"{'|cond6/logged-1|':>17} {'steps cond>1e5':>15}")
    for name, fn in arms.items():
        with recording() as log:
            run = fn()
        d = per_step(log, run, with_exact=False)
        gauge = float((d["sv"][:, 6] / d["sv"][:, 0]).max())
        dev = float(np.abs(d["cond_r_F"] / d["cond_logged"] - 1).max())
        n_fold = int(np.sum(d["cond_logged"] > STIFF_COND))
        print(f"  {name:<16} {run['E'][-1] - E_exact:12.4e} {len(run['tau']) - 1:6d} "
              f"{str(sorted(set(d['kept'].tolist()))):>6} {gauge:10.2e} {dev:17.2e} "
              f"{n_fold:15d}")
        compare(f"{name} error", run["E"][-1] - E_exact, f"e3/{name}/error", indent=6)
        compare(f"{name} steps", len(run["tau"]) - 1, f"e3/{name}/steps", indent=6)
        compare(f"{name} steps in the fold", n_fold, f"e3/{name}/fold", indent=6)
        claim(f"{name}: no leakage, pinv keeps 6, rank-6 cond = logged",
              gauge < GAUGE_BOUND and set(d["kept"].tolist()) == {R_STRUCT} and dev < 1e-12,
              indent=6)


def p1_exact_metric_rerun(plots=True):
    banner("exact  Target C on the EXACT M, Y (no finite differences)")
    H, E_exact, theta0, rho_fn = _c_setup()
    print(f"  {'rcond':>8} {'arm':<6} {'error':>14} {'steps':>6} {'max cond_M':>12}")
    for rcond in (1e-8, 1e-2):
        fd = aq.adaptive_qite_run(theta0, rho_fn, H, 9, 9, rcond=rcond)
        with recording(replacement=exact_measure_f):
            ex = aq.adaptive_qite_run(theta0, rho_fn, H, 9, 9, rcond=rcond)
        for label, run in (("fd", fd), ("exact", ex)):
            print(f"  {rcond:8g} {label:<6} {run['E'][-1] - E_exact:14.6e} "
                  f"{len(run['tau']) - 1:6d} {np.nanmax(run['cond_M']):12.5e}")
        compare("exact-metric error", ex["E"][-1] - E_exact, f"exact/{rcond:g}/error")
        compare("exact-metric steps", len(ex["tau"]) - 1, f"exact/{rcond:g}/steps")
        # compare on a common tau grid, not by step index (the tau grids differ slightly)
        grid = fd["tau"][fd["tau"] <= min(fd["tau"][-1], ex["tau"][-1])]
        gap = np.abs(np.interp(grid, fd["tau"], fd["E"]) - np.interp(grid, ex["tau"], ex["E"]))
        j = int(np.argmax(gap))
        print(f"    max |E_fd(tau) - E_exact_metric(tau)| = {gap[j]:.3e} at tau = {grid[j]:.4f}; "
              f"at the common end {gap[-1]:.3e}")
        compare("max energy gap fd vs exact metric", float(gap[j]), f"exact/{rcond:g}/gap")
        before = gap[grid < 0.15]
        print(f"    same gap restricted to tau < 0.15 (before the Target C fold): "
              f"{before.max():.3e}")
        print()


def p1_target_a_control(plots=True):
    banner("a  control: Target A (5-param, full rank, rcond = 1e-6)")
    H = aq.h2_hamiltonian()
    E_exact = aq.exact_ground_energy(H)
    theta0 = np.random.default_rng(aq.SEED_A).uniform(-0.5, 0.5, size=5)
    rho_fn = aq.make_rho_fn(n_params=5)
    with recording() as log:
        run = aq.adaptive_qite_run(theta0, rho_fn, H, 5, 5)
    d = per_step(log, run, r=5, H=H)
    compare("final error", run["E"][-1] - E_exact, "a/error")
    compare("steps", len(run["tau"]) - 1, "a/steps")
    dropped = d["kept"] < 5
    print(f"    kept by pinv, per step: {sorted(set(d['kept'].tolist()))}; pinv drops a PHYSICAL "
          f"direction on {int(dropped.sum())} of {len(dropped)} steps"
          + (f", tau in [{d['tau'][dropped].min():.3f}, {d['tau'][dropped].max():.3f}]"
             if dropped.any() else ""))
    compare("steps where pinv keeps < 5", int(dropped.sum()), "a/dropped")
    i = int(np.argmin(d["sv"][:, 4] / d["sv"][:, 0]))
    theta, M, _, _ = log[i]
    s, vecs = null_vector_attribution(M, 5)
    names = _gate_names(5)
    print(f"    softest step tau = {d['tau'][i]:.4f}, sigma_5/sigma_1 = {s[4] / s[0]:.3e} "
          f"(cutoff {aq.RCOND_DEFAULT:g})")
    print("    v5 loadings: " + ", ".join(f"{n} {v:+.3f}" for n, v in zip(names, vecs[0])))
    compare("v5 loading on th5 (CRZ)", float(abs(vecs[0][4])), "a/v5_crz")
    # CRZ is the last gate, so Gt_5 = |1><1|_0 (x) Z_1 / 2 and M_55 = 2 Var of that.
    _, _, psi = exact_metric_and_force(theta, 5, H)
    pops = np.abs(psi) ** 2
    z1 = float(np.vdot(psi, np.kron(aq.I2, aq.PZ) @ psi).real)
    print(f"    basis populations |00>,|01>,|10>,|11> = {np.array2string(pops, precision=5)}")
    print(f"    <Z_1> = {z1:.5f}; M[CRZ,CRZ] = {M[4, 4]:.4e}")
    print("    -> the state is within ~1e-3 of a computational basis state, on which CRZ only")
    print("       adds a phase: its rotated generator has near-zero variance. That is the")
    print("       plan's sec. 2 zero-variance mechanism, developing as the state converges.")
    compare("max basis population at the softest step", float(pops.max()), "a/pop_max")
    print(f"    |theta_dot| max/min x{d['speed_theta'].max() / d['speed_theta'].min():.1f}; "
          f"||drho/dtau|| at the softest step {d['speed_state'][i]:.3e}")


def p1_regime_aware_cap(plots=True):
    banner("cap  the plan's follow-on: cap on cond over the top r = 6 values only")
    print("  Implemented by swapping aq.route3_indicators so its cond_M slot carries the rank-6")
    print("  Frobenius cond. STIFFNESS_CAP_COEFF = 2000 untouched.\n")
    H, E_exact, theta0, rho_fn = _c_setup()
    original = aq.route3_indicators

    def rank6_indicators(M, Y, theta_dot_val, rcond):
        residual, _, sigma_min, cond_spectral = original(M, Y, theta_dot_val, rcond)
        return residual, rank_r_cond(M, R_STRUCT)[0], sigma_min, cond_spectral

    runs = {}
    print(f"  {'rcond':>8} {'cap on':<22} {'error':>12} {'steps':>6} {'cap binding':>12} "
          f"{'min dtau':>10}")
    for rcond in (1e-8, 1e-2):
        runs[(rcond, "frobenius")] = aq.adaptive_qite_run(theta0, rho_fn, H, 9, 9, rcond=rcond)
        aq.route3_indicators = rank6_indicators
        try:
            runs[(rcond, "rank6")] = aq.adaptive_qite_run(theta0, rho_fn, H, 9, 9, rcond=rcond)
        finally:
            aq.route3_indicators = original
        for kind in ("frobenius", "rank6"):
            run = runs[(rcond, kind)]
            label = "cond_F (shipped)" if kind == "frobenius" else "rank-6 cond_F (P1)"
            print(f"  {rcond:8g} {label:<22} {run['E'][-1] - E_exact:12.4e} "
                  f"{len(run['tau']) - 1:6d} {int(np.sum(run['cap_bound'])):12d} "
                  f"{run['delta_tau'][1:].min():10.4g}")
        compare("rank-6 cap error", runs[(rcond, "rank6")]["E"][-1] - E_exact,
                f"cap/{rcond:g}/error", indent=6)
        compare("rank-6 cap steps", len(runs[(rcond, "rank6")]["tau"]) - 1,
                f"cap/{rcond:g}/steps", indent=6)

    if plots:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), sharey=True)
        for ax, rcond in zip(axes, (1e-8, 1e-2)):
            for kind, color, ls in (("frobenius", COLORS[0], "-"), ("rank6", COLORS[1], "--")):
                run = runs[(rcond, kind)]
                ax.plot(run["tau"][1:], run["delta_tau"][1:], ls, color=color, lw=1.6,
                        label="cap on cond$_F$ (shipped)" if kind == "frobenius"
                        else "cap on rank-6 cond$_F$")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(r"imaginary time $\tau$")
            ax.set_title(f"Target C, rcond={rcond:g}")
            ax.grid(alpha=0.3)
        axes[0].set_ylabel(r"$\Delta\tau$")
        axes[1].legend(frameon=True, fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGURES / "p1_regime_aware_cap.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'p1_regime_aware_cap.png'}")


# --------------------------------------------------------------------------------------

def print_summary():
    banner("Assumptions (ours, on top of adaptive_qite.py's and adam_qite.py's)")
    for i, a in enumerate(ASSUMPTIONS, 1):
        print(f"  {i}. {a}\n")


TARGETS = {
    "ids": check_identities,
    "c": p1_target_c_replay,
    "null": p1_null_vectors,
    "e3": p1_e3_arms,
    "exact": p1_exact_metric_rerun,
    "a": p1_target_a_control,
    "cap": p1_regime_aware_cap,
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
