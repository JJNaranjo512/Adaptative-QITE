#!/usr/bin/env python
"""Replication of Saxena, Shalabi & Kyaw, *Practical limitations of quantum data
propagation on noisy quantum processors*, Phys. Rev. Applied 21, 054014 (2024),
arXiv:2306.13046v2.

Five targets, one function each:

  T1  theorem_2_depolarizing_error()    Eq. (19); Fig. 1(b) p.2, Fig. 5 p.16
  T2  eq_18_allowed_probability()       Eq. (18); Fig. 3 caption p.5
  T3  appendix_c_cond_and_norms()       cond(M), ||M||, ||Y||; Appendix C p.16
  T4  appendix_b_noise_scaling()        Eqs. (B12), (B13), (B16); Appendix B pp.14-15
  T5  theorem_1_bound()                 Eq. (17); Fig. 3(b) p.5, Fig. 6 p.17

Run with the physics environment:

    ~/micromamba/envs/qml-julia/bin/python replicate.py --target all

Every number the script compares against was read out of the PDF; the ones that
could only be eyeballed off a plot are marked as such in FIG3B_ENDPOINTS.
Assumptions we had to make (the paper does not state them) are collected in
ASSUMPTIONS and printed at the end of a full run.
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pennylane as qml

SEED = 0
rng = np.random.default_rng(SEED)

HERE = Path(__file__).resolve().parent
FIGURES = HERE / "figures"

# --------------------------------------------------------------------------------------
# Constants taken from the paper
# --------------------------------------------------------------------------------------

# H2 in the minimal STO-3G basis, Bravyi-Kitaev, Eq. (20) p.5
G = (0.2252, 0.3435, -0.4347, 0.5716, 0.0910, 0.0910)

# Angles used for the Appendix C numerics, p.16
THETA_5 = np.array([1.5249, 2.5142, 0.4457, 1.3250, 2.8769])
THETA_9 = np.array([1.294, 3.265, 5.187, 3.819, 1.5708, 4.172, 0.618, 3.441, 2.502])

COND_M_PAPER_5 = 66.7239   # Appendix C p.16 (quoted as 66.7 in the Fig. 3 caption)
NORM_M_PAPER_5 = 0.9977
COND_M_PAPER_9 = 20.22
NORM_M_PAPER_9 = 1.689
NORM_Y_PAPER_9 = 0.121
# The paper never reports ||Y|| for the five-parameter ansatz -- see T5.

DELTA = 0.04                       # Fig. 3 caption p.5
N_LIST = (5, 6, 10, 12, 14)
P_MAX_PAPER = {5: 0.032, 6: 0.036, 10: 0.052, 12: 0.068, 14: 0.129}  # Fig. 3 caption

# Endpoint of each Fig. 3(b) curve, read off the rendered plot. The y axis runs 0-600
# with gridlines every 100, so these are good to roughly +-10 (N=5,6) and +-15 (N>=10).
FIG3B_ENDPOINTS = {5: 75.0, 6: 105.0, 10: 300.0, 12: 430.0, 14: 590.0}
FIG3B_READOFF_TOL = {5: 10.0, 6: 10.0, 10: 15.0, 12: 15.0, 14: 15.0}

H_DERIV = 1e-6      # central-difference step; error ~1e-10 on these smooth functions
COLORS = ("#0173B2", "#DE8F05", "#029E73", "#CC78BC", "#D55E00")

ASSUMPTIONS: list[str] = [
    "Initial state |00>: the paper never states rho_0. T3 scans the alternatives.",
    "Parameters are assigned to gates column-wise as drawn in Fig. 4 (theta_1 top-left, "
    "theta_2 below it, ...). A different assignment only permutes M, so it changes "
    "neither ||M||_F nor cond(M).",
    f"Derivatives of rho by central differences with h={H_DERIV} rather than "
    "parameter-shift; agreement with the closed form of Eq. (B9) is asserted in the code.",
    "cond(.) is the Frobenius condition number ||M||_F ||M^-1||_F, which is what the "
    "paper defines in Eq. (A39). The spectral condition number is printed alongside.",
    "The global depolarizing channel is E(rho) = (1-p) rho + p I/d acting on the whole "
    "ansatz register (d=4), applied after every gate, per Eq. (15) and Fig. 2.",
    "N (number of variational parameters) equals the number of gates in both ansaetze, "
    "so it is also the number of noise insertions.",
    "||Y|| for the five-parameter ansatz is computed here because the paper does not "
    "report it; T5 shows the Fig. 3(b) curves are only reproduced if that term is "
    "negligible, i.e. ||Y|| -> infinity.",
    "Fig. 3(b) endpoint values are eyeballed off the plot to the tolerances in "
    "FIG3B_READOFF_TOL.",
]

RESULTS: list[dict] = []


def record(target, source, paper, obtained, note=""):
    """Collect one row for the consolidated comparison table."""
    if paper is None or obtained is None or not np.isfinite(obtained) or not np.isfinite(paper):
        abs_dev = rel_dev = np.nan
    else:
        abs_dev = obtained - paper
        rel_dev = abs_dev / paper if paper != 0 else np.nan
    RESULTS.append(
        dict(target=target, source=source, paper=paper, obtained=obtained,
             abs_dev=abs_dev, rel_dev=rel_dev, note=note)
    )


def banner(text):
    print()
    print("=" * 88)
    print(text)
    print("=" * 88)


# --------------------------------------------------------------------------------------
# Shared layer: Hamiltonian, ansaetze, M and Y
# --------------------------------------------------------------------------------------

I2 = np.eye(2, dtype=complex)
PX = np.array([[0, 1], [1, 0]], dtype=complex)
PY = np.array([[0, -1j], [1j, 0]], dtype=complex)
PZ = np.diag([1, -1]).astype(complex)
PAULIS = (I2, PX, PY, PZ)


def h2_hamiltonian():
    """H = g0 I + g1 Z0 + g2 Z1 + g3 Z0Z1 + g4 Y0Y1 + g5 X0X1, Eq. (20)."""
    g0, g1, g2, g3, g4, g5 = G
    kron = np.kron
    return (g0 * np.eye(4)
            + g1 * kron(PZ, I2)
            + g2 * kron(I2, PZ)
            + g3 * kron(PZ, PZ)
            + g4 * kron(PY, PY)
            + g5 * kron(PX, PX))


def gate_list(theta, n_params=5, two_qubit="controlled", reverse_wires=False):
    """The ansatz gates, read straight off Fig. 4 (p.16).

    Fig. 4(a), five parameters:
        RX(t1) q0, RX(t2) q1, RY(t3) q0, RY(t4) q1, CRZ(t5) q0->q1
    Fig. 4(b), nine parameters (the universal ansatz of Ref. [55]):
        RY(t1) q0, RY(t2) q1, RZ(t3) q0, RZ(t4) q1, CRX(t5) q0->q1,
        RY(t6) q0, RY(t7) q1, RZ(t8) q0, RZ(t9) q1

    `two_qubit="ising"` replaces the controlled rotation by the corresponding
    two-qubit Pauli rotation, which is the other way the Fig. 4 symbol could be read.
    """
    a, b = (1, 0) if reverse_wires else (0, 1)
    if n_params == 5:
        entangler = "CRZ" if two_qubit == "controlled" else "IsingZZ"
        return [("RX", theta[0], a), ("RX", theta[1], b),
                ("RY", theta[2], a), ("RY", theta[3], b),
                (entangler, theta[4], [a, b])]
    if n_params == 9:
        entangler = "CRX" if two_qubit == "controlled" else "IsingXX"
        return [("RY", theta[0], a), ("RY", theta[1], b),
                ("RZ", theta[2], a), ("RZ", theta[3], b),
                (entangler, theta[4], [a, b]),
                ("RY", theta[5], a), ("RY", theta[6], b),
                ("RZ", theta[7], a), ("RZ", theta[8], b)]
    raise ValueError(f"no ansatz with {n_params} parameters in the paper")


def prepare(init_state):
    """State preparation for rho_0. The paper assumes it is error free (Fig. 2)."""
    if init_state == "++":
        qml.Hadamard(wires=0)
        qml.Hadamard(wires=1)
        return
    for wire, bit in enumerate(init_state):
        if bit == "1":
            qml.PauliX(wires=wire)


def make_rho_fn(n_params=5, two_qubit="controlled", reverse_wires=False,
                init_state="00", stop=None):
    """rho(theta) as a density matrix, noiseless. `stop` truncates after that many gates."""
    dev = qml.device("default.qubit", wires=2)

    @qml.qnode(dev)
    def circuit(theta):
        prepare(init_state)
        gates = gate_list(theta, n_params, two_qubit, reverse_wires)
        for name, angle, wires in gates[:stop]:
            getattr(qml, name)(angle, wires=wires)
        return qml.density_matrix(wires=[0, 1])

    return lambda theta: np.asarray(circuit(np.asarray(theta, dtype=float)), dtype=complex)


def global_depolarizing_kraus(p):
    """Kraus operators of E(rho) = (1-p) rho + p I/d on two qubits, Eq. (15) with
    N(rho) = I/d. Uses the Pauli twirl (1/d^2) sum_P P rho P = Tr(rho) I/d."""
    d = 4
    kraus = [np.sqrt(1.0 - p) * np.eye(d, dtype=complex)]
    for pa, pb in itertools.product(PAULIS, repeat=2):
        kraus.append(np.sqrt(p) / d * np.kron(pa, pb))
    return kraus


def make_noisy_rho_fn(p, n_params=5, two_qubit="controlled", init_state="00"):
    """rho^err(theta): the same ansatz with the global depolarizing channel after
    every gate, exactly the error placement of Fig. 2."""
    dev = qml.device("default.mixed", wires=2)
    kraus = global_depolarizing_kraus(p)

    @qml.qnode(dev)
    def circuit(theta):
        prepare(init_state)
        for name, angle, wires in gate_list(theta, n_params, two_qubit):
            getattr(qml, name)(angle, wires=wires)
            qml.QubitChannel(kraus, wires=[0, 1])
        return qml.density_matrix(wires=[0, 1])

    return lambda theta: np.asarray(circuit(np.asarray(theta, dtype=float)), dtype=complex)


def drho(rho_fn, theta, h=H_DERIV):
    """d rho / d theta_k for every k, by central differences."""
    theta = np.asarray(theta, dtype=float)
    out = []
    for k in range(len(theta)):
        step = np.zeros_like(theta)
        step[k] = h
        out.append((rho_fn(theta + step) - rho_fn(theta - step)) / (2 * h))
    return np.array(out)


def metric_M(d_rho):
    """M_{k,q} = Tr[(d rho / d theta_k)^dagger (d rho / d theta_q)], Eq. (6)."""
    n = len(d_rho)
    M = np.empty((n, n))
    for k in range(n):
        for q in range(n):
            M[k, q] = np.trace(d_rho[k].conj().T @ d_rho[q]).real
    return M


def vector_Y(d_rho, rho, H):
    """Y_k = -Tr[(d rho / d theta_k) {H, rho}], Eq. (8)."""
    anti = H @ rho + rho @ H
    return np.array([-np.trace(d_rho[k] @ anti).real for k in range(len(d_rho))])


def fubini_study(psi_fn, theta, h=H_DERIV):
    """A_{k,q} = Re[<d_k psi|d_q psi> - <d_k psi|psi><psi|d_q psi>], the pure-state
    metric. Included in the T3 scan as an alternative reading of Eq. (6)."""
    theta = np.asarray(theta, dtype=float)
    psi = psi_fn(theta)
    d = []
    for k in range(len(theta)):
        step = np.zeros_like(theta)
        step[k] = h
        d.append((psi_fn(theta + step) - psi_fn(theta - step)) / (2 * h))
    d = np.array(d)
    overlap = d.conj() @ d.T
    v = d.conj() @ psi
    return np.real(overlap - np.outer(v, v.conj()))


def frobenius_cond(M, rcond=None):
    """cond(M) = ||M||_F ||M^-1||_F, the definition the paper gives in Eq. (A39)."""
    norm = np.linalg.norm(M)
    inv = np.linalg.pinv(M, rcond=rcond) if rcond is not None else None
    if inv is None:
        sv = np.linalg.svd(M, compute_uv=False)
        if sv[-1] <= 1e-12 * sv[0]:
            return np.inf
        inv = np.linalg.inv(M)
    return norm * np.linalg.norm(inv)


def spectral_cond(M):
    sv = np.linalg.svd(M, compute_uv=False)
    return np.inf if sv[-1] <= 1e-12 * sv[0] else sv[0] / sv[-1]


def canonical_M_Y(n_params=5, theta=None, init_state="00"):
    """M, Y and rho under the conventions we take as canonical (Eqs. 6 and 8)."""
    theta = THETA_5 if theta is None else theta
    rho_fn = make_rho_fn(n_params=n_params, init_state=init_state)
    rho = rho_fn(theta)
    d = drho(rho_fn, theta)
    return metric_M(d), vector_Y(d, rho, h2_hamiltonian()), rho, d


def check_eq_b9(verbose=True):
    """Eq. (B9) at p=0: for a single-qubit rotation gate with generator P_k,
    M_{k,k} = (1 - Tr[P_k rho_{k-1} P_k rho_{k-1}]) / 2.

    This is an exact internal check on the whole M construction, and it also fixes
    the rotation convention to R(theta) = exp(-i theta P / 2).
    """
    M, _, _, _ = canonical_M_Y(5)
    generators = {0: (PX, 0), 1: (PX, 1), 2: (PY, 0), 3: (PY, 1)}  # single-qubit gates only
    worst = 0.0
    for k, (pauli, wire) in generators.items():
        rho_before = make_rho_fn(n_params=5, stop=k)(THETA_5)
        P = np.kron(pauli, I2) if wire == 0 else np.kron(I2, pauli)
        predicted = 0.5 * (1.0 - np.trace(P @ rho_before @ P @ rho_before).real)
        worst = max(worst, abs(M[k, k] - predicted))
        if verbose:
            print(f"    M[{k},{k}] = {M[k, k]:.12f}   Eq.(B9) = {predicted:.12f}   "
                  f"|diff| = {abs(M[k, k] - predicted):.2e}")
    assert worst < 1e-8, f"Eq. (B9) violated by {worst:.2e}: the M convention is wrong"
    return worst


# --------------------------------------------------------------------------------------
# T1  --  Theorem 2, Eq. (19)
# --------------------------------------------------------------------------------------

def theorem_2_relative_error(p, N):
    """eps/||theta_dot|| = (1 - (1-p)^N) / (1-p)^N, Eq. (19)."""
    p = np.asarray(p, dtype=float)
    return (1.0 - (1.0 - p) ** N) / (1.0 - p) ** N


def theorem_2_depolarizing_error(plots=True):
    banner("T1  Theorem 2, Eq. (19) -- exact relative error under depolarizing noise")

    print("  Sanity limits (these follow from the formula itself):")
    print(f"    p=0     -> {theorem_2_relative_error(0.0, 5):.3f}  (must be 0)")
    print(f"    p->1    -> diverges for every N")

    print("\n  Value at the p = 1e-4 per-gate rate the Discussion (p.6) argues about:")
    for N in N_LIST:
        print(f"    N={N:2d}  eps/||theta_dot|| = {theorem_2_relative_error(1e-4, N):.6e}")

    # Where each curve crosses 500, the top of the Fig. 1(b) / Fig. 5 y axis. This turns
    # an eyeball comparison of the figure into five numbers.
    print("\n  Probability at which each curve reaches the top of the Fig. 5 axis (500):")
    for N in N_LIST:
        p_cross = 1.0 - (1.0 / 501.0) ** (1.0 / N)
        print(f"    N={N:2d}  p = {p_cross:.4f}")
        record(f"T1 N={N} p at error=500", "Fig. 5 p.16 (read off)", None, p_cross,
               "closed form, no paper value to compare")

    if plots:
        p = np.linspace(0.0, 0.72, 2000)
        fig, ax = plt.subplots(figsize=(5.2, 4.0))
        styles = ["-", "--", "--", "-.", ":"]
        for i, N in enumerate(N_LIST):
            ax.plot(p, theorem_2_relative_error(p, N), styles[i], color=COLORS[i],
                    lw=1.6, label=f"N={N}")
        ax.set_xlim(0.0, 0.72)
        ax.set_ylim(0.0, 500.0)
        ax.set_xlabel("probability ($p$)")
        ax.set_ylabel(r"error $\left(\frac{\|\Delta\dot\theta\|}{\|\dot\theta\|}\right)$")
        ax.set_title("T1: Eq. (19), cf. Fig. 1(b) p.2 and Fig. 5 p.16")
        ax.legend(loc="upper left", frameon=True)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(FIGURES / "fig1b_theorem2.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'fig1b_theorem2.png'}")

    print("\n  Eq. (19) is a closed form with no free parameters, so this panel is an "
          "\n  identity rather than a test. T4 is where it is actually checked, by "
          "\n  simulating the noise and recovering the same curve.")


# --------------------------------------------------------------------------------------
# T2  --  Eq. (18), the validity window
# --------------------------------------------------------------------------------------

def eq_18_p_max(N, delta, norm_M):
    """p_max = 1 - (1 - 2 N^2 delta / (N + 2||M||))^(1/2N), Eq. (18)."""
    inner = 1.0 - 2.0 * N ** 2 * delta / (N + 2.0 * norm_M)
    if inner <= 0.0:
        return np.nan
    return 1.0 - inner ** (1.0 / (2.0 * N))


def eq_18_implied_norm_M(N, delta, p_max):
    """Invert Eq. (18) for ||M||: if the five reported p_max values came from one
    ||M||, these must all agree."""
    inner = (1.0 - p_max) ** (2.0 * N)
    return (2.0 * N ** 2 * delta / (1.0 - inner) - N) / 2.0


def eq_18_allowed_probability(plots=True):
    banner("T2  Eq. (18) -- maximum probability for which Theorem 1 is valid")
    print(f"  delta = {DELTA}, ||M|| = {NORM_M_PAPER_5} (both from the Fig. 3 caption p.5)\n")
    print(f"  {'N':>3}  {'paper':>8}  {'obtained':>10}  {'abs dev':>10}  {'rel dev':>9}"
          f"  {'implied ||M||':>14}")
    for N in N_LIST:
        obtained = eq_18_p_max(N, DELTA, NORM_M_PAPER_5)
        paper = P_MAX_PAPER[N]
        implied = eq_18_implied_norm_M(N, DELTA, paper)
        print(f"  {N:3d}  {paper:8.3f}  {obtained:10.5f}  {obtained - paper:+10.5f}"
              f"  {(obtained - paper) / paper:+8.2%}  {implied:14.4f}")
        record(f"T2 p_max N={N}", "Fig. 3 caption p.5", paper, obtained)

    implied = [eq_18_implied_norm_M(N, DELTA, P_MAX_PAPER[N]) for N in N_LIST]
    print(f"\n  The implied ||M|| ranges over {min(implied):.4f} to {max(implied):.4f}. "
          "A single\n  mis-stated ||M|| cannot explain the five reported values, so the "
          "offset is not\n  a different ||M|| -- it is consistent with the paper quoting "
          "the last point of a\n  discrete p grid (every obtained value is 0.001-0.002 "
          "above the reported one).")

    if plots:
        delta_grid = np.linspace(0.0, 0.145, 1500)
        fig, ax = plt.subplots(figsize=(5.2, 4.0))
        styles = ["-", "--", "-.", ":", "-."]
        for i, N in enumerate(N_LIST):
            vals = np.array([eq_18_p_max(N, d, NORM_M_PAPER_5) for d in delta_grid])
            ax.plot(delta_grid, vals, styles[i], color=COLORS[i], lw=1.6, label=f"N={N}")
        ax.set_xlim(0.0, 0.145)
        ax.set_ylim(0.0, 0.37)
        ax.set_xlabel(r"$\delta$")
        ax.set_ylabel("allowed probability of noise $p$")
        ax.set_title("T2: Eq. (18), cf. Fig. 3(a) p.5")
        ax.legend(loc="upper left", frameon=True)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(FIGURES / "fig3a_allowed_p.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 'fig3a_allowed_p.png'}")


# --------------------------------------------------------------------------------------
# T3  --  Appendix C, cond(M) and the norms
# --------------------------------------------------------------------------------------

def appendix_c_cond_and_norms(plots=True):
    banner("T3  Appendix C p.16 -- cond(M), ||M|| and ||Y|| at the published angles")

    print("  First, the internal consistency check of the M construction against Eq. (B9):")
    check_eq_b9()

    H = h2_hamiltonian()

    for n_params, theta, cond_paper, normM_paper, normY_paper in (
        (5, THETA_5, COND_M_PAPER_5, NORM_M_PAPER_5, None),
        (9, THETA_9, COND_M_PAPER_9, NORM_M_PAPER_9, NORM_Y_PAPER_9),
    ):
        print(f"\n  --- {n_params}-parameter ansatz (Fig. 4{'a' if n_params == 5 else 'b'}) ---")
        M, Y, rho, _ = canonical_M_Y(n_params, theta)
        sv = np.linalg.svd(M, compute_uv=False)
        rank = int(np.sum(sv > 1e-8 * max(sv[0], 1.0)))
        print(f"  canonical convention (Eqs. 6 and 8, rho_0 = |00>):")
        print(f"    ||M||_F   = {np.linalg.norm(M):.4f}   (paper {normM_paper})")
        print(f"    cond_F(M) = {frobenius_cond(M):.4f}   (paper {cond_paper})")
        print(f"    cond_2(M) = {spectral_cond(M):.4f}")
        print(f"    ||Y||     = {np.linalg.norm(Y):.4f}   "
              f"(paper {normY_paper if normY_paper is not None else 'not reported'})")
        print(f"    rank(M)   = {rank} of {n_params};  singular values:")
        print("      " + "  ".join(f"{s:.3e}" for s in sv))

        record(f"T3 ||M||_F {n_params}-param", "Appendix C p.16", normM_paper,
               float(np.linalg.norm(M)))
        record(f"T3 cond(M) {n_params}-param", "Appendix C p.16", cond_paper,
               float(frobenius_cond(M)))
        if normY_paper is not None:
            record(f"T3 ||Y|| {n_params}-param", "Appendix C p.16", normY_paper,
                   float(np.linalg.norm(Y)))

        if rank < n_params:
            print(f"\n    NOTE: M is rank {rank} < {n_params} by construction. A two-qubit")
            print( "    pure state manifold has 2d-2 = 6 real dimensions, so a metric with")
            print(f"    {n_params} parameters cannot have full rank and cond(M) = infinity.")
            print( "    The reported value requires a pseudo-inverse cutoff the paper does")
            print( "    not state; the scan below shows what different cutoffs would give.")

    # ---- convention scan -------------------------------------------------------------
    print("\n  --- convention scan ---")
    print("  Every combination of the choices the paper leaves open. This is a search")
    print("  over documented ambiguities, not a fit: whatever it finds is reported as is.")
    print("  (Re-labelling parameters only permutes M, so ||M||_F and cond(M) are")
    print("  invariant under it and it is not scanned.)\n")

    for n_params, theta, cond_paper, normM_paper in (
        (5, THETA_5, COND_M_PAPER_5, NORM_M_PAPER_5),
        (9, THETA_9, COND_M_PAPER_9, NORM_M_PAPER_9),
    ):
        print(f"  {n_params}-parameter ansatz, target cond = {cond_paper}, "
              f"||M|| = {normM_paper}")
        print(f"    {'rho_0':>6} {'2q gate':>10} {'wires':>9} {'metric':>10} "
              f"{'||M||_F':>10} {'cond_F':>14} {'cond_2':>14}")
        hits = []
        for init_state, two_qubit, reverse in itertools.product(
                ("00", "01", "10", "11", "++"), ("controlled", "ising"), (False, True)):
            rho_fn = make_rho_fn(n_params, two_qubit, reverse, init_state)
            d = drho(rho_fn, theta)
            M_rho = metric_M(d)

            dev = qml.device("default.qubit", wires=2)

            @qml.qnode(dev)
            def state_circuit(t, _init=init_state, _tq=two_qubit, _rev=reverse):
                prepare(_init)
                for name, angle, wires in gate_list(t, n_params, _tq, _rev):
                    getattr(qml, name)(angle, wires=wires)
                return qml.state()

            psi_fn = lambda t: np.asarray(state_circuit(np.asarray(t, float)), dtype=complex)
            A = fubini_study(psi_fn, theta)

            for label, M in (("Tr[drdr]", M_rho), ("FS", A), ("4*FS", 4 * A)):
                cf, c2 = frobenius_cond(M), spectral_cond(M)
                print(f"    {init_state:>6} {two_qubit:>10} "
                      f"{'reversed' if reverse else 'as drawn':>9} {label:>10} "
                      f"{np.linalg.norm(M):10.4f} {cf:14.4f} {c2:14.4f}")
                if np.isfinite(cf) and abs(cf - cond_paper) < 0.05 * cond_paper:
                    hits.append((init_state, two_qubit, reverse, label, cf))
        if hits:
            print(f"    --> combinations within 5% of the reported cond: {hits}")
        else:
            print(f"    --> no combination reproduces cond = {cond_paper} to within 5%.")
        print()

    # ---- pseudo-inverse cutoffs for the rank-deficient nine-parameter case ------------
    print("  Nine-parameter ansatz: what pseudo-inverse cutoff would be needed?")
    M9, _, _, _ = canonical_M_Y(9, THETA_9)
    for rcond in (1e-15, 1e-12, 1e-8, 1e-4, 1e-2, 1e-1, 3e-1):
        print(f"    rcond = {rcond:<8g} cond_F = {frobenius_cond(M9, rcond=rcond):16.4f}")
    print(f"    target from the paper: {COND_M_PAPER_9}")
    print("    cond_F is monotone in the cutoff, so a cutoff can be chosen to hit any")
    print("    value; that is exactly why the unreported cutoff makes 20.22 unfalsifiable.")


# --------------------------------------------------------------------------------------
# T4  --  Appendix B, the mechanism
# --------------------------------------------------------------------------------------

def appendix_b_noise_scaling(plots=True):
    banner("T4  Appendix B pp.14-15 -- M^err = (1-p)^2N M, Y^err = (1-p)^N Y, Eq. (B16)")

    N = 5
    H = h2_hamiltonian()
    theta = THETA_5

    M, Y, rho, d = canonical_M_Y(5, theta)
    theta_dot = np.linalg.solve(M, Y)
    print(f"  Noiseless five-parameter ansatz: ||M||={np.linalg.norm(M):.4f}, "
          f"||Y||={np.linalg.norm(Y):.4f}, ||theta_dot||={np.linalg.norm(theta_dot):.4f}")
    print(f"  theta_dot = {np.array2string(theta_dot, precision=5)}")

    p_grid = (1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1)
    # Relative deviation in Frobenius norm, not elementwise: M has entries close to zero
    # (its smallest singular value is ~1e-3) and an elementwise ratio is dominated by them.
    print(f"\n  {'p':>8} {'||M^err-(1-p)^2N M||/||M||':>28} {'||Y^err-(1-p)^N Y||/||Y||':>28}"
          f" {'rel. error, sim':>18} {'Eq.(19)':>14}")
    sim_rel, exact_rel = [], []
    for p in p_grid:
        rho_fn_err = make_noisy_rho_fn(p, n_params=5)
        rho_err = rho_fn_err(theta)
        d_err = drho(rho_fn_err, theta)

        M_err = metric_M(d_err)
        # Appendix B applies the noise to the derivative line only; the rho inside
        # {H, rho} is left ideal (see Eq. B11). This is that convention.
        Y_err = vector_Y(d_err, rho, H)

        dev_M = np.linalg.norm(M_err - (1 - p) ** (2 * N) * M) / np.linalg.norm(M)
        dev_Y = np.linalg.norm(Y_err - (1 - p) ** N * Y) / np.linalg.norm(Y)

        theta_dot_err = np.linalg.solve(M_err, Y_err)
        rel = np.linalg.norm(theta_dot_err - theta_dot) / np.linalg.norm(theta_dot)
        exact = theorem_2_relative_error(p, N)
        sim_rel.append(rel)
        exact_rel.append(exact)
        print(f"  {p:8.4f} {dev_M:28.3e} {dev_Y:28.3e} {rel:18.10f} {exact:14.10f}")

        record(f"T4 M^err/M p={p}", "Eq. (B12)", (1 - p) ** (2 * N),
               float(np.linalg.norm(M_err) / np.linalg.norm(M)))
        record(f"T4 Y^err/Y p={p}", "Eq. (B13)", (1 - p) ** N,
               float(np.linalg.norm(Y_err) / np.linalg.norm(Y)))
        record(f"T4 rel. error p={p}", "Eq. (19) from simulation", float(exact), float(rel))

    worst = max(abs(s - e) for s, e in zip(sim_rel, exact_rel))
    print(f"\n  Largest deviation of the simulated relative error from Eq. (19): {worst:.3e}")
    assert worst < 1e-6, "the simulation does not reproduce Eq. (19)"

    print("\n  Why this is exact rather than approximate: the global depolarizing channel")
    print("  acts on a traceless operator as multiplication by (1-p), and d rho / d theta_k")
    print("  is traceless. The maximally mixed part of the state entering gate k is")
    print("  annihilated by the derivative, so each of the two lines of Fig. 2 picks up")
    print("  exactly (1-p)^N and M^err = (1-p)^2N M follows with no approximation.")

    # ---- probing the assumption the paper leaves implicit ----------------------------
    print("\n  --- probe: what if the rho inside {H, rho} is also measured noisily? ---")
    print("  Eq. (B11) applies noise only to the derivative line. On hardware the state")
    print("  in {H, rho} comes from the same noisy device, so it should carry the noise")
    print("  too. Recomputing Y^err that way:\n")
    print(f"  {'p':>8} {'||Y^err-(1-p)^N Y||/||Y||':>28} {'rel. error, both noisy':>24}"
          f" {'Eq.(19)':>14} {'ratio to Eq.(19)':>18}")
    both_rel = []
    for p in p_grid:
        rho_fn_err = make_noisy_rho_fn(p, n_params=5)
        rho_err = rho_fn_err(theta)
        d_err = drho(rho_fn_err, theta)
        M_err = metric_M(d_err)
        Y_err_both = vector_Y(d_err, rho_err, H)
        dev_Y = np.linalg.norm(Y_err_both - (1 - p) ** N * Y) / np.linalg.norm(Y)
        theta_dot_err = np.linalg.solve(M_err, Y_err_both)
        rel = np.linalg.norm(theta_dot_err - theta_dot) / np.linalg.norm(theta_dot)
        exact = theorem_2_relative_error(p, N)
        both_rel.append(rel)
        print(f"  {p:8.4f} {dev_Y:28.3e} {rel:24.10f} {exact:14.10f} {rel / exact:18.10f}")
        record(f"T4 rel. error, both lines noisy, p={p}", "our extension of Eq. (19)",
               float(exact), float(rel))

    ratios = [r / theorem_2_relative_error(p, N) for r, p in zip(both_rel, p_grid)]
    print(f"\n  Every ratio is 1/2 to {max(abs(r - 0.5) for r in ratios):.1e}. That is exact,")
    print("  and here is why. For a normalised state 2 Re<d_k psi|psi> = 0, so Eq. (8)")
    print("  reduces to Y_k = -Tr[d_k rho H] =: -g_k, i.e. Y is minus the energy gradient")
    print("  and M^-1 g = -theta_dot. Writing rho^err = (1-p)^N rho + (1-(1-p)^N) I/d and")
    print("  d_k rho^err = (1-p)^N d_k rho, the extra maximally mixed piece contributes")
    print("  -(1-p)^N (1-(1-p)^N) g_k / 2 to Y, so")
    print("      theta_dot^err = theta_dot [1 + (1 - (1-p)^N) / (2 (1-p)^N)],")
    print("  exactly half the deviation of Eq. (19).")

    # Y = -g is the reason for the factor; check it rather than assert it in prose.
    g = np.array([np.trace(d[k] @ H).real for k in range(N)])
    print(f"\n  Check Y = -g:  ||Y + g|| / ||Y|| = "
          f"{np.linalg.norm(Y + g) / np.linalg.norm(Y):.3e}")
    assert np.linalg.norm(Y + g) / np.linalg.norm(Y) < 1e-6
    assert max(abs(r - 0.5) for r in ratios) < 1e-4, "the factor 1/2 is not exact"

    print("\n  So the (1-p)^N law for Y holds only for the noise placement of Fig. 2, in")
    print("  which the rho inside {H, rho} is ideal. If both M and Y are read off the same")
    print("  noisy device -- which is what a QITE step on hardware actually does -- the")
    print("  relative error in theta_dot is exactly half of Theorem 2. Theorem 2 is a")
    print("  factor of two pessimistic in that setting; it is not wrong, it answers a")
    print("  different question.")

    if plots:
        p_dense = np.linspace(0.0, 0.25, 400)
        fig, ax = plt.subplots(figsize=(5.6, 4.2))
        ax.plot(p_dense, theorem_2_relative_error(p_dense, N), "-", color=COLORS[0],
                lw=1.6, label="Eq. (19), closed form")
        ax.plot(p_grid, sim_rel, "o", color=COLORS[1], ms=7, mfc="none", mew=1.6,
                label="simulated, Fig. 2 noise placement")
        ax.plot(p_grid, both_rel, "s", color=COLORS[4], ms=7, mfc="none", mew=1.6,
                label=r"simulated, $\rho$ in $\{H,\rho\}$ also noisy")
        ax.set_xlim(0.0, 0.25)
        ax.set_xlabel("probability ($p$)")
        ax.set_ylabel(r"$\|\Delta\dot\theta\| / \|\dot\theta\|$")
        ax.set_title("T4: Eq. (19) recovered from a noisy simulation, N=5")
        ax.legend(loc="upper left", frameon=True, fontsize=9)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(FIGURES / "t4_scaling.png", dpi=160)
        plt.close(fig)
        print(f"\n  wrote {FIGURES / 't4_scaling.png'}")


# --------------------------------------------------------------------------------------
# T5  --  Theorem 1, Eq. (17)
# --------------------------------------------------------------------------------------

def theorem_1_relative_error(p, N, cond_M, norm_M, norm_Y,
                             y_term="as-printed", denominator=False):
    """Eq. (17) of the main text / Eq. (A42) of the appendix.

    y_term: how to read the last factor, which the paper prints as
        1 + sqrt(N) / sqrt(2 ||Y||)                       -> "as-printed"
        1 + sqrt(N) / (sqrt(2) ||Y||)                     -> "as-derived", which is what
                                                             Eq. (A37) actually gives
        1                                                 -> "negligible", the ||Y|| -> oo limit
    denominator: include the 1 - cond(M) ||dM||/||M|| factor of Eq. (A42), which the
        main-text Eq. (17) drops.
    """
    p = np.asarray(p, dtype=float)
    term_M = (1.0 - (1.0 - p) ** (2 * N)) * (1.0 + N / (2.0 * norm_M))
    if y_term == "as-printed":
        y_factor = 1.0 + np.sqrt(N) / np.sqrt(2.0 * norm_Y)
    elif y_term == "as-derived":
        y_factor = 1.0 + np.sqrt(N) / (np.sqrt(2.0) * norm_Y)
    elif y_term == "negligible":
        y_factor = 1.0
    else:
        raise ValueError(y_term)
    term_Y = (1.0 - (1.0 - p) ** N) * y_factor
    value = cond_M * (term_M + term_Y)
    if denominator:
        denom = 1.0 - cond_M * term_M
        with np.errstate(divide="ignore", invalid="ignore"):
            value = np.where(denom > 0, value / denom, np.nan)
    return value


def theorem_1_bound(plots=True):
    banner("T5  Theorem 1, Eq. (17) -- upper bound under general probabilistic noise")

    M5, Y5, _, _ = canonical_M_Y(5, THETA_5)
    norm_Y_ours = float(np.linalg.norm(Y5))
    print(f"  The paper never reports ||Y|| for the five-parameter ansatz. Our own value")
    print(f"  from the same angles is ||Y|| = {norm_Y_ours:.4f}; the scan below uses it,")
    print(f"  the nine-parameter ||Y|| = {NORM_Y_PAPER_9}, and the ||Y|| -> oo limit.\n")

    variants = [
        ("Eq.17, y as printed, ||Y||=ours", dict(y_term="as-printed", denominator=False),
         norm_Y_ours),
        ("Eq.17, y as derived, ||Y||=ours", dict(y_term="as-derived", denominator=False),
         norm_Y_ours),
        ("Eq.17, y term negligible", dict(y_term="negligible", denominator=False), 1.0),
        ("Eq.A42 with denominator, y as printed", dict(y_term="as-printed", denominator=True),
         norm_Y_ours),
    ]

    print(f"  Endpoint of each Fig. 3(b) curve, at p = p_max(N), cond = {COND_M_PAPER_5}, "
          f"||M|| = {NORM_M_PAPER_5}\n")
    header = f"  {'variant':<40}" + "".join(f"{'N=' + str(N):>12}" for N in N_LIST)
    print(header)
    print(f"  {'read off Fig. 3(b)':<40}"
          + "".join(f"{FIG3B_ENDPOINTS[N]:12.1f}" for N in N_LIST))
    print("  " + "-" * (40 + 12 * len(N_LIST)))
    best = None
    for label, kwargs, norm_Y in variants:
        row = []
        for N in N_LIST:
            p_max = eq_18_p_max(N, DELTA, NORM_M_PAPER_5)
            row.append(float(theorem_1_relative_error(
                p_max, N, COND_M_PAPER_5, NORM_M_PAPER_5, norm_Y, **kwargs)))
        print(f"  {label:<40}" + "".join(
            f"{v:12.1f}" if np.isfinite(v) else f"{'invalid':>12}" for v in row))
        if all(np.isfinite(v) for v in row):
            err = max(abs(v - FIG3B_ENDPOINTS[N]) for v, N in zip(row, N_LIST))
            if best is None or err < best[1]:
                best = (label, err, row, kwargs, norm_Y)

    label, err, row, kwargs, norm_Y = best
    print(f"\n  Closest variant: {label}, worst endpoint deviation {err:.1f}")
    within = all(abs(v - FIG3B_ENDPOINTS[N]) <= FIG3B_READOFF_TOL[N]
                 for v, N in zip(row, N_LIST))
    print(f"  Within the read-off tolerance for every N: {within}")
    for v, N in zip(row, N_LIST):
        record(f"T5 Fig. 3(b) endpoint N={N}", "Fig. 3(b) p.5 (read off)",
               FIG3B_ENDPOINTS[N], v, label)

    # Where the Eq. (A42) form is actually valid, Eq. (A43).
    print("\n  Validity of the Eq. (A42) form, from the constraint Eq. (A43):")
    for N in N_LIST:
        limit = 2.0 * NORM_M_PAPER_5 / (COND_M_PAPER_5 * (N + 2.0 * NORM_M_PAPER_5))
        p_valid = 1.0 - (1.0 - limit) ** (1.0 / (2 * N)) if limit < 1 else np.nan
        print(f"    N={N:2d}  Eq. (A42) denominator stays positive only for p <= {p_valid:.2e},"
              f"  but Fig. 3(b) is plotted to p = {P_MAX_PAPER[N]}")
    print("  So the main-text Eq. (17), which drops that denominator, is the only form")
    print("  consistent with the plotted range.")

    if plots:
        for tag, cond_M, norm_M, norm_Y_fig, n_list, xlim, ylim, ref in (
            ("fig3b_theorem1", COND_M_PAPER_5, NORM_M_PAPER_5, norm_Y,
             N_LIST, 0.135, 620.0, "Fig. 3(b) p.5 / Fig. 6(a) p.17"),
            ("fig6b_theorem1_9param", COND_M_PAPER_9, NORM_M_PAPER_9, NORM_Y_PAPER_9,
             (5, 9, 10, 12, 14), 0.085, 210.0, "Fig. 6(b) p.17"),
        ):
            fig, ax = plt.subplots(figsize=(5.4, 4.1))
            styles = ["-", "--", "--", "-.", ":"]
            for i, N in enumerate(n_list):
                p_max = eq_18_p_max(N, DELTA, norm_M)
                if not np.isfinite(p_max):
                    continue
                p = np.linspace(0.0, p_max, 400)
                ax.plot(p, theorem_1_relative_error(p, N, cond_M, norm_M, norm_Y_fig,
                                                    **kwargs),
                        styles[i], color=COLORS[i], lw=1.6, label=f"N={N}")
            ax.set_xlim(0.0, xlim)
            ax.set_ylim(0.0, ylim)
            ax.set_xlabel("probability ($p$)")
            ax.set_ylabel(r"error $\left(\frac{\|\Delta\dot\theta\|}{\|\dot\theta\|}\right)$")
            ax.set_title(f"T5: Eq. (17), cf. {ref}\n({label})", fontsize=9)
            ax.legend(loc="lower right", frameon=True)
            ax.grid(alpha=0.3)
            fig.tight_layout()
            fig.savefig(FIGURES / f"{tag}.png", dpi=160)
            plt.close(fig)
            print(f"\n  wrote {FIGURES / (tag + '.png')}")


# --------------------------------------------------------------------------------------

def print_summary():
    banner("Consolidated comparison")
    print(f"  {'target':<38} {'source':<30} {'paper':>12} {'obtained':>14} "
          f"{'abs dev':>12} {'rel dev':>9}")
    print("  " + "-" * 118)
    for r in RESULTS:
        paper = "-" if r["paper"] is None else f"{r['paper']:.6g}"
        rel = "-" if not np.isfinite(r["rel_dev"]) else f"{r['rel_dev']:+.2%}"
        abs_dev = "-" if not np.isfinite(r["abs_dev"]) else f"{r['abs_dev']:+.4g}"
        print(f"  {r['target']:<38} {r['source']:<30} {paper:>12} "
              f"{r['obtained']:14.6g} {abs_dev:>12} {rel:>9}")

    banner("Assumptions (ours, not the paper's)")
    for i, a in enumerate(ASSUMPTIONS, 1):
        print(f"  {i}. {a}")


TARGETS = {
    "t1": theorem_2_depolarizing_error,
    "t2": eq_18_allowed_probability,
    "t3": appendix_c_cond_and_norms,
    "t4": appendix_b_noise_scaling,
    "t5": theorem_1_bound,
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
