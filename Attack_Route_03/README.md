# Attack_Route_03

Implementation of Sections 1–3 of the internal working memo *Adaptive Time-Stepping for
Variational QITE: An RKF-Type Controller with a Hardware-Noise Floor* (5pp, no arXiv/DOI),
restricted to **Route 3** of Section 3.3 ("cheap physics-based error indicators — no extra
stages"). Routes 1, 2, 4 and 5 (embedded RKF pairs, step-doubling, statistical gating,
stiffness-aware solvers) are explicitly out of scope.

This is an **implementation of a proposed method, not a reproduction of reported numbers** —
the memo contains no numerical claims to score against, so there is no verdict table here.
Everything below is what the code actually produced when run, not a projection.

Source PDF: `../Research_Papers/QITE_timestep.pdf`
Explanation note: `../Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation.md` — read
that one for the full discussion of every assumption and finding.
Builds on: `../Replications/practical-limitations-qdp/` (same H2 system and $M\dot\theta=Y$
machinery; that replication's rank-deficiency finding on the 9-parameter ansatz motivates
Target C directly).

## Run

```bash
~/micromamba/envs/qml-julia/bin/python adaptive_qite.py --target all
```

`--target {a,b,c,d,all}` selects a single scenario; `--no-plots` skips the PNGs.
Deterministic: each target now draws `theta_0` from its own `default_rng(SEED_{A,B,C,D})`
(`SEED=0..3`), so a target reproduces identically whether run alone or inside `--target all`
(fixed after the project-wide audit, [[QITE_Hand_Off/Audits/Audit_Project_Wide/M5_Reproducibility_Hardening|M5]] — previously all four targets shared one
module-level `rng`, so a target's result silently depended on which other targets had already
run in the same invocation). Runtime for `--target all` is under 2 minutes. Dependencies
(numpy, scipy, matplotlib, pennylane) are all already in `qml-julia` — nothing new was
installed.

> [!info] Corrected after the project-wide audit
> This table reflects the state after [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|the project-wide audit]]'s
> modules M1–M6, all implemented and re-run. The pre-audit numbers (a wrong
> depolarizing debias factor, an exact-ITE energy gate, and a metric that silently
> extrapolated stalled runs) are retired, not archived here — see the audit and module notes
> for the before/after comparison and the derivations.

## Targets

| | Scenario | Takeaway |
|---|---|---|
| A | Stage 1 (§1): noiseless adaptive Euler QITE vs. a fixed-step baseline, 5-param ansatz | Adaptive control reaches a lower energy than the fixed-step baseline even at the corrected, honest circuit-count budget (−1.1210 vs. −1.1167 at 732 circuits, [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting\|M6]]) — retries were rare enough (3 of 172 steps) that this barely differs from the naive f-eval comparison. Relative to the ansatz's own variational floor ($E_{floor}=-1.1246$, computed by direct minimization), adaptive finishes $1.6\times10^{-3}$ above it vs. fixed's $7.9\times10^{-3}$ — a cleaner statement than the absolute error, most of which is the ansatz's own expressibility gap. The controller still throttles $\Delta\tau$ near $\tau\approx2.5$, and after correcting the energy gate to the *variational* prediction ([[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence\|M4]]) the reason histogram shows this is genuinely the cond(M) cap binding (31 of the last 43 steps) and not the energy-band gate (0 of 43) — the stiffness read survives, on firmer evidence than before. |
| B | Stage 2 (§2): global depolarizing noise, with/without Eq. (7) debiasing, 5-param ansatz | The debias factor was wrong ([[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction\|M1]], F1): the code used the memo's Eq. (7) factor $q=(1-p)^N$, which assumes noise on the derivative line only; this script's noise is on the whole circuit, so the correct factor is $q^2+q(1-q)/2$. Fixed, the debiased per-step flow now matches the clean flow to the finite-difference noise floor ($<10^{-4}\%$ at every $p$ tested), vs. $5.3\%/34.7\%$ biased. The noisy flow is also **exactly collinear** with the clean one (collinearity residual $\sim10^{-8}$) — "flow-direction error" was a misnomer ([[QITE_Hand_Off/Audits/Audit_Project_Wide/M2_Noise_Model_Falsifiability\|M2]], F2); it is a pure magnitude rescale. Even with the flow-level bug fixed, the **full-trajectory** RMS still shows no consistent debiasing advantage (debiased vs. biased RMS within $\sim2\%$ of each other at both $p$ tested) — a genuine post-fix finding, not the retracted "interlock". A reparametrization test ($\Delta\tau\to\Delta\tau/c(p)$ on the clean run) tracks the biased run measurably better than the plain baseline does (11.9% vs. 12.5% RMS at $p=0.02$; 35.2% vs. 50.6% at $p=0.1$), consistent with the scalar rescale being partly absorbed as a time reparametrization by the adaptive controller. |
| C | Route 3 point 3: cond(M) stiffness proxy on the 9-parameter (rank-6-of-9) ansatz | After the [[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence\|M4]] gate fix, both cutoffs converge (no more 500-step stall): tight (`rcond=1e-8`) reaches error $8.9\times10^{-4}$ in 108 steps; loose (`rcond=1e-2`) reaches $2.1\times10^{-3}$ in 17 steps. Tight still beats loose by $\sim2.4\times$ — the core finding survives the correction, at a smaller margin than the pre-fix numbers suggested (those were partly the old gate over-throttling the loose-cutoff run, not the rank deficiency alone). |
| D | Our extension of §2 to Route 3 (explicitly **not** the memo's Eq. 9) | RMS energy-trajectory deviation (common-support only, [[QITE_Hand_Off/Audits/Audit_Project_Wide/M3_Trajectory_Comparison_Metric\|M3]]) from the noiseless run, swept over $p\in\{10^{-4},10^{-3},10^{-2},5\times10^{-2}\}$, debiased vs. biased — no run stalls post-M1 (coverage 100% at every $p$). Deviation grows smoothly with $p$ (RMS $6.4\times10^{-4}$ at $p=10^{-4}$ to $\approx0.29$ at $p=0.05$), with no sharp knee — the earlier "$p\approx10^{-3}$ noise floor" reading is **not supported**; it was an artifact of the broken metric (M3) combined with the wrong debias factor (M1). Debiased and biased RMS remain within a few percent of each other at every $p$, same as Target B. |

> [!warning] Hand-off correction
> The $5\times$ divides end points at different budgets: Route 3 after 2069 circuits ($1.607\times10^{-3}$ above $E_{floor}$) against fixed Euler after 732 ($7.909\times10^{-3}$). At a matched 732 circuits Route 3 is at $-1.120952$ (corrected read-out, G6 §6.3), $3.648\times10^{-3}$ above the floor, so $2.2\times$ closer. The status report states the matched figure.

## Figures

| File | Scenario |
|---|---|
| `figures/target_a_convergence.png` | Target A: energy error vs. $\tau$ (adaptive vs. fixed step) and the $\Delta\tau$ trajectory |
| `figures/target_b_debiasing.png` | Target B: energy trajectories, noiseless + four noisy variants |
| `figures/target_c_stiffness.png` | Target C: cond(M) and $\Delta\tau$ vs. $\tau$ for two `rcond` choices |
| `figures/target_d_noise_floor.png` | Target D: RMS energy deviation vs. $p$, debiased vs. biased |

## Layout of `adaptive_qite.py`

Constants first (H2 Hamiltonian coefficients `G`, and every Route-3 controller threshold,
each one flagged in `ASSUMPTIONS` as invented for this script — the memo gives no numbers).
Then the shared layer, copied and trimmed from
`../Replications/practical-limitations-qdp/replicate.py` so results are directly comparable:
`h2_hamiltonian`, `gate_list`, `prepare`, `make_rho_fn`, `make_noisy_rho_fn`,
`global_depolarizing_kraus`, `drho`, `metric_M`, `vector_Y`, `frobenius_cond`, plus
`exact_ground_energy` and `energy_and_variance` added for this task.

Then the three stages:
- **Stage 1** (§1) — `measure_f(..., p=0.0)`: builds $M$, $Y$, solves $\dot\theta=M^+Y$ via
  `np.linalg.pinv(M, rcond=rcond)`, `rcond` always an explicit logged parameter.
- **Stage 2** (§2) — the same `measure_f`, now with `p>0` and `debias`: rescales $M,Y$ by their
  known depolarizing bias factors before the solve when `debias=True` (Eq. 7). `p=0` makes it a
  no-op, so Stage 1 and Stage 2 share one function.
- **Stage 3** (§3.3, Route 3 only) — `route3_indicators` (McLachlan residual, cond(M)) and
  `energy_check` (the Eq. 13 variance-identity accept/reject gate); `adaptive_qite_run` combines
  all three into one controller (the specific combination rule is *our* synthesis — see the
  note); `fixed_step_qite_run` is the Stage-1-only baseline used for comparison in Target A.

`ASSUMPTIONS` at the top is the authoritative list of every choice the memo left open (there is
no companion "verified against the paper" list, since there is nothing in the memo to verify
against); `--target all` prints it.

Post-audit additions ([[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|audit]] modules M1–M6):
`route3_indicators` also returns `sigma_min`/`cond_spectral` (unregularized, M4);
`trajectory_deviation` replaces raw `np.interp` RMS with a common-support-only comparison plus
coverage (M3); `target_b_per_step_check` gives the flow-magnitude-error/collinearity numbers a
code path (M2); `check_eq_b9_self` validates this script's own `metric_M`/`drho` against Eq.
(B9), not just the QDP replication's copy (M5); `variational_floor` minimizes $E(\theta)$
directly over the ansatz as a reference line for Target A (M6).
