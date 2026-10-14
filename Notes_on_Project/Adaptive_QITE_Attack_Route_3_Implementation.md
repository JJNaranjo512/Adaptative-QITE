---
title: "Implementation — Adaptive QITE timestep control, Route 3 (Sections 1–3)"
tags:
  - implementation
  - physics
  - quantum-computation
  - adaptative-qite
code: QITE_Hand_Off/Attack_Route_03/adaptive_qite.py
---

# Implementation — Adaptive QITE Timestep Control, Route 3

Code for §1–3 (Route 3 only) of an internal 5-page working memo, *Adaptive Time-Stepping for
Variational QITE: An RKF-Type Controller with a Hardware-Noise Floor* — no arXiv/DOI, source
PDF at `QITE_Hand_Off/Research_Papers/QITE_timestep.pdf`. Per the earlier answer
on this task, there's no formal `Paper_Reading` note for the memo since it has no bibliographic
fields to fill — this note stands in for it and links straight to the PDF.

Code: `QITE_Hand_Off/Attack_Route_03/adaptive_qite.py`, run with

```bash
~/micromamba/envs/qml-julia/bin/python adaptive_qite.py --target all
```

`SEED = 0`, everything below is from one deterministic run. Builds directly on
[[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]] (same H2 system, same $M\dot\theta=Y$ machinery) and
on the finding in [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] that the 9-parameter ansatz's $M$
is exactly rank-6-of-9 — Target C below is that finding's direct consequence for an adaptive
controller.

**This is an implementation of a proposed method, not a reproduction of reported numbers.** The
memo has no numerical claims to score against, so there's no target/verdict table here — just
what the code actually did, plus an honest accounting of two places where the result surprised
me.

> [!warning] Corrected after the project-wide audit — read [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]] first
> The project-wide audit found the debias factor below (Assumption 3) was wrong for this
> script's noise placement, and that the Target B "genuine surprise" I wrote up was that bug,
> not physics — see [[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|M1_Depolarizing_Debias_Correction]]. Fixing it, plus the controller gate
> and metric issues in [[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence|M4_Controller_Near_Convergence]]/[[QITE_Hand_Off/Audits/Audit_Project_Wide/M3_Trajectory_Comparison_Metric|M3_Trajectory_Comparison_Metric]],
> changed the Target A/B/C/D numbers below materially. I'm leaving the original write-up in
> place (it's an honest record of what I actually thought at the time, and *why* I thought it)
> and marking every paragraph that the audit superseded, rather than rewriting history.

## What the memo says (§1–3, notation translated to this vault's conventions)

**§1 — QITE as an ODE.** McLachlan's principle gives $M(\theta)\dot\theta = Y(\theta)$; the
memo frames the imaginary-time flow as the autonomous ODE $\dot\theta = f(\theta) :=
M^+(\theta)Y(\theta)$, with explicit-Euler as the base discretization,
$\theta_{n+1} = \theta_n + \Delta\tau\,\dot\theta_n$. It then sketches embedded RKF pairs
(multiple $f$-evaluations per step, an error estimate from their difference) — that's Routes 1
and 2, out of scope here.

**§2 — the hardware-noise floor.** Global depolarizing noise biases the measured flow by a
known, purely multiplicative factor: $\tilde f = (1-p)^{-N} f$ (Eq. 7), so it "should be removed
by rescaling... before the integrator ever sees the data." What's left after that debiasing is a
smaller, genuinely stochastic residual $\varepsilon_{res}$ (shot noise, imperfect knowledge of
$p$, non-depolarizing components). The memo's own error estimator $e_n$ (Routes 1/2's embedded
pair) then has a noise floor, Eq. (9): below it, $e_n$ is pure noise and the controller
stalls or oscillates.

**§3.3 — Route 3.** Three indicators that need no extra circuit stages, only quantities already
being computed or one cheap extra observable:
1. McLachlan residual $r_n = \|M\dot\theta - Y\|$ — a spike means the linear step model is
   failing, shrink $\Delta\tau$.
2. Energy monotonicity / variance identity, Eq. (13): $dE/d\tau = -2\,\mathrm{Var}_\rho(H) \le
   0$. Predict $\widehat{\Delta E} = -2\,\mathrm{Var}(H)\Delta\tau$, accept the step iff the
   measured $\Delta E$ agrees within a band; reject/shrink on a monotonicity violation.
3. $\mathrm{cond}(M)$ as a cheap stiffness proxy — growing $\mathrm{cond}(M)$ ⇒ cap
   $\Delta\tau$ or switch regularization.

The memo never combines these three into one controller itself — its own worked algorithm (§4)
combines Route 1 + Route 4 instead. Synthesizing the three Route 3 indicators into a single
step-control rule is this task's own construction (see Assumption 2 below).

## What was actually built

`adaptive_qite.py` reuses the H2 (STO-3G, Bravyi-Kitaev) system, both ansätze, and the
$M$/$Y$/$\mathrm{cond}_F(M)$ machinery from `Replications/practical-limitations-qdp/replicate.py`
essentially unchanged, so results here are directly comparable to that replication's numbers.

**Stage 1/2 share one function**, `measure_f(theta, rho_fn, H, N, p=0.0, rcond, debias=True)`:
builds $M$, $Y$ from finite-difference density-matrix derivatives (same method already validated
against the closed form Eq. B9 to 1e-10 in the QDP replication — not Hadamard-test circuits, see
Assumption 5), rescales them by the known depolarizing bias factors when `p>0` and `debias=True`
(applying Eq. 7 to $M$ and $Y$ individually, before the solve — see Assumption 3), then solves
$\dot\theta = M^+Y$ via `np.linalg.pinv(M, rcond=rcond)` with `rcond` always an explicit,
logged argument rather than a silent default. `p=0` makes debiasing a no-op, so this one
function covers both §1 and §2.

**Stage 3** is `route3_indicators` (residual, $\mathrm{cond}_F(M)$) plus `energy_check` (the Eq.
13 predicted-vs-actual gate), combined in `adaptive_qite_run`. The combination rule I landed on,
after the first version locked up on the singular 9-parameter case (see Assumption 2's history
below):

- $\mathrm{cond}(M)$ is a **proactive cap**: before even trying a step,
  $\Delta\tau \leftarrow \min(\Delta\tau,\ \mathrm{COEFF}/\mathrm{cond}(M),\ \Delta\tau_{max})$,
  floored at $\Delta\tau_{min}$ so a literally-infinite cond(M) throttles to the smallest
  allowed step rather than to zero.
- The energy check is the only **delta_tau-dependent** accept/reject/retry gate (it's the only
  one of the three whose signal changes as $\Delta\tau$ shrinks), matching the memo's own
  "accept... reject/shrink" language for point 2.
- The residual is delta_tau-*independent* (it's fixed once $\theta_n$ is, before any trial step
  is even taken), so retrying with a smaller $\Delta\tau$ can't change it. It's implemented as a
  **caution signal**: a high relative residual ($r_n/\|Y\|$ above 5%) freezes growth for the
  *next* step instead of rejecting the current one.

`adaptive_qite_run` ties this into the imaginary-time loop; `fixed_step_qite_run` is a plain
fixed-$\Delta\tau$ Euler baseline (Stage 1 only) used for the Target A comparison.

## Why Route 3 specifically

No extra quantum circuit stages beyond one $f$-evaluation and two cheap classical-postprocessing
observables ($E$, $\mathrm{Var}(H)$) per step — versus Routes 1/2's $O(N^2)$-circuit embedded
stages. It's also the route the memo's own §5 discussion calls "likely the practical entry point
before full embedded pairs are attempted." And its cond(M) indicator ties directly back to this
project's own load-bearing finding (the rank-6-of-9 metric) rather than being an abstract choice.

## Assumptions (all ours, the memo specifies none of this)

1. **Every controller threshold is invented.** `RCOND_DEFAULT`, `DELTA_TAU0/MIN/MAX`, `GROW`,
   `SHRINK`, `STIFFNESS_CAP_COEFF`, `ENERGY_MONOTONE_TOL`, `DE_BAND_ABS/REL`,
   `RESIDUAL_REL_TOL` were tuned only by inspecting Target A's own trajectory — never against a
   paper-reported value, because there isn't one. `STIFFNESS_CAP_COEFF=2000` in particular was
   chosen so the cap only meaningfully bites as $\mathrm{cond}(M)$ approaches
   $1/\texttt{RCOND\_DEFAULT}\,(=10^6)$, i.e. roughly where the pseudoinverse solve itself starts
   being regularized away — a smaller coefficient (I tried 50) made the cap dominate the whole
   5-parameter run and it never reached $\tau_{max}$ within the step budget.
2. **The three-indicator combination rule is my synthesis**, not the memo's (see above): cond(M)
   proactive cap, energy check reactive retry, residual a growth-freeze. The first version I
   wrote treated all three as retry-triggering rejections; on the 9-parameter (cond$=\infty$)
   ansatz this made the proactive cap literally divide by infinity and freeze $\Delta\tau$ at
   exactly zero forever — a real bug, not a finding, fixed by flooring the cap at
   `delta_tau_min`. Worth recording because it's exactly the kind of silent-bad-solve failure
   mode the plan's verification step asked to rule out.
3. ~~**Eq. (7) debiasing is applied to $M,Y$ individually** before the solve
   ($M\leftarrow M/(1-p)^{2N}$, $Y\leftarrow Y/(1-p)^N$), matching Algorithm 1's literal step
   order ("measure $\tilde M,\tilde Y$; de-bias by (7); solve for $k_1$") and the closed forms
   $M^{err}=(1-p)^{2N}M$, $Y^{err}=(1-p)^NY$ already validated to ~1e-10 in the QDP replication —
   algebraically identical to the memo's $\tilde f=(1-p)^{-N}f$ when $M$ is full rank.~~
   **Wrong — corrected after the project-wide audit, [[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|M1_Depolarizing_Debias_Correction]].** Eq. (7)'s
   $Y\leftarrow Y/(1-p)^N$ assumes noise on the derivative line only. This script's noise sits on
   the whole circuit (`make_noisy_rho_fn` applies the channel after every gate), so
   $\rho^{err}$ is mixed and enters both $\partial_k\rho$ and the $\{H,\rho\}$ anticommutator in
   `vector_Y`. Working it through with $q:=(1-p)^N$: $M^{err}=q^2M$ (my $M$ factor was right),
   but $Y^{err}=[q^2+q(1-q)/2]\,Y_{clean}$, not $qY_{clean}$. **This is exactly the "Theorem 2 is
   a factor of two pessimistic" result I already derived in
   [[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]] and never propagated here** — the two weeks it sat
   unused is the actual lesson. Corrected code: `M /= q**2; Y /= q**2 + q*(1-q)/2`.
4. **Eq. (9)'s literal noise floor does not transfer to Route 3** and is not computed here — it's
   defined in terms of the embedded-pair estimator $e_n$, which Route 3 never forms. Target D
   builds an explicitly-labelled substitute instead (below).
5. **$M,Y$ come from finite-difference density-matrix derivatives** (`drho`, $h=10^{-6}$) on
   `default.qubit`/`default.mixed`, the same method already used and validated in the QDP
   replication — not physical Hadamard-test circuits. "Evaluation count" in Target A counts
   calls to `measure_f`, not circuit shots.
6. **$\theta_0 = $ `rng.uniform(-0.5, 0.5, size=N_PARAMS)`** at `SEED=0`, not the QDP
   replication's published Appendix-C angles (those exist to reproduce a specific reported
   cond(M), not as a sensible QITE starting point). $\rho_0=|00\rangle$, matching the QDP
   replication's own unstated convention.
7. **No shot noise.** The depolarizing channel is applied to the exact density matrix
   (deterministic bias, not a stochastic estimator), so "noisy" here means biased-but-exact.
   Route 4 (statistical gating) — explicitly out of scope — is the mechanism that would make
   shot noise necessary to model.

## Results (from the `--target all` run above; see `Attack_Route_03/README.md` for the figure list)

**Target A — noiseless, 5-parameter ansatz.** *Numbers below are post-audit (
[[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence|M4_Controller_Near_Convergence]], [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting|M6_Honest_Cost_Accounting]]); the original "zero
violations of Eq. (13)" framing conflated a real monotonicity check with a mis-specified energy
gate — see M4.* The adaptive run reaches E = −1.122993 (172 steps). At the **corrected, honest
circuit-count budget** (M6 — not raw f-evaluations, which undercounted every retry's extra
circuit) of 732 circuits, adaptive E = −1.120992 vs. the fixed-step baseline's −1.116691 —
adaptive still wins at equal cost, though the margin barely moves from the naive comparison
because only 3 of 172 steps needed a retry in this run. Against the ansatz's own **variational
floor** ($E_{floor}=-1.124600$, direct minimization, M6) — a reference line this note previously
lacked — adaptive finishes $1.6\times10^{-3}$ above it, fixed $7.9\times10^{-3}$ above it: most
of the absolute error against $E_{exact}$ is the ansatz's own expressibility gap, not the
integrator.

$\mathrm{cond}(M)$ still climbs into the $10^5$–$10^6$ range as $\tau\to3$, and $\Delta\tau$
still collapses there — but now on firmer evidence. The energy gate previously used the *exact*
-ITE identity $-2\mathrm{Var}(H)\Delta\tau$ as the accept/reject test; the ansatz's deficit
doesn't vanish at convergence, so that gate fights its own goal near $\tau_{max}$ (see M4). Fixed
to gate on the *variational* prediction $-Y\cdot\dot\theta\,\Delta\tau$ instead (keeping the
exact identity as a diagnostic only), and re-run: the reason histogram over the last quarter of
the run shows the cond(M) cap binding on 31 of 43 steps and the energy band triggering a retry on
0 — the $\Delta\tau$ collapse is genuinely the cond(M) cap, not the gate misfiring. The
unregularized spectral cond (also new, M4) independently reaches $1.4\times10^6$, not bounded by
the `rcond` cutoff the way the regularized $\mathrm{cond}_F$ is — so this reads as a real
stiffness event, on stronger footing than the original write-up had.

**Target B — noisy, 5-parameter ansatz, debiasing.** ~~This is the one genuine surprise.~~ **It
wasn't — it was Assumption 3's wrong bias factor.** Full retraction,
[[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|M1_Depolarizing_Debias_Correction]] and [[QITE_Hand_Off/Audits/Audit_Project_Wide/M2_Noise_Model_Falsifiability|M2_Noise_Model_Falsifiability]]:

The "~5% at $p=0.01$ growing to ~67% at $p=0.2$" per-step number I reported was real (it's what a
wrong-but-internally-consistent bias factor produces), but the code path for it didn't exist
until M2 added `target_b_per_step_check()` — I'd computed it by hand and never committed the
check. Run now: with the *correct* factor $q^2+q(1-q)/2$, the debiased flow matches the clean
flow to the finite-difference noise floor (<0.0001% at every $p$) — Eq. (7)'s corrected form is
exact for this noise model, not just "doing what it says." And the noisy flow direction is
**exactly collinear** with the clean one (collinearity residual $\sim10^{-8}$, $c_{meas}$ matches
the closed form $c(p)=1+(1-q)/(2q)$ to 6 digits) — "flow-direction error" was never the right
name for this; it's a pure magnitude rescale.

I also reused the same theta_0 assumption twice without checking it: Targets B and D draw from
the shared module-level `rng` in sequence, so "reproduced identically across Target B and Target
D" compared two runs with *different* $\theta_0$ — the claim was false as written
([[QITE_Hand_Off/Audits/Audit_Project_Wide/M5_Reproducibility_Hardening|M5_Reproducibility_Hardening]], now each target has its own seed).

What survives, now that the flow-level bug is actually fixed: **the full-trajectory RMS still
shows no consistent debiasing advantage** (debiased vs. biased within ~2% of each other at both
$p=0.02$ and $p=0.1$, common-support metric, [[QITE_Hand_Off/Audits/Audit_Project_Wide/M3_Trajectory_Comparison_Metric|M3_Trajectory_Comparison_Metric]]) — this is a
*real*, post-fix finding, not the residual-freeze interlock I originally invoked (that
explanation is retracted along with the bug it was built to explain). A reparametrization test is
a better account of it: run the noiseless controller with every $\Delta\tau$ divided by $c(p)$,
and it tracks the biased run measurably better than the plain baseline does (RMS 11.9% vs. 12.5%
at $p=0.02$; 35.2% vs. 50.6% at $p=0.1$) — consistent with the scalar flow rescale being partly
absorbed by the adaptive controller as a reparametrization of imaginary time, which is a much
simpler mechanism than an indicator interlock, and a genuine structural limit: nothing homogeneous
in $\dot\theta$ can tell $c\,\dot\theta$ apart from a faster flow.

**Target C — 9-parameter (rank-6-of-9) ansatz, stiffness proxy.** *Numbers below are post-M4* —
the pre-fix numbers (`rcond=1e-8`: 20 steps, error $1.2\times10^{-2}$; `rcond=1e-2`: stalled at
the 500-step cap, error $2.9\times10^{-2}$) turned out to be dominated by the old exact-ITE
energy gate over-throttling both runs, the loose-cutoff one especially. Re-run with the corrected
variational gate: `rcond=1e-8` reaches error $8.9\times10^{-4}$ in 108 steps ($\mathrm{cond}(M)$
ranging 32.9–$9.7\times10^6$); `rcond=1e-2` reaches $2.1\times10^{-3}$ in 17 steps
($\mathrm{cond}(M)$ 19.1–125), no longer stalling at all. **The core finding survives**: tight
still beats loose by $\sim2.4\times$ in final error, so a lower cond(M) reading still doesn't
mean a better-conditioned solve — but the margin is far smaller than the 500-step-stall numbers
suggested, because most of *that* gap was the broken gate, not the rank deficiency alone. Still
exactly the "regularization cutoff has to be promoted to an explicit, tunable part of the method"
point from [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]], just a more honest version of it.

**Target D — our noise-floor extension (not Eq. 9).** *Post-M1/M3.* With the debias bug fixed
and the trajectory metric restricted to common support (no run stalls anymore — coverage 100% at
every $p$), RMS energy deviation from the noiseless run grows smoothly with $p$ (from
$6.4\times10^{-4}$ at $p=10^{-4}$ to $\approx0.29$ at $p=0.05$), with **no sharp knee anywhere in
this range**. ~~If I were picking a Route-3-specific analogue of Eq. (9)'s "floor"... below
about $p=10^{-3}$...~~ **Retracted** — that reading was built on the pre-fix, stall-contaminated
RMS curve (M3) layered on the wrong debias factor (M1); with both corrected there's no floor-like
feature to point at, just smooth growth. Debiasing gives no consistent advantage here either,
same as Target B, and for the same reason (the flow rescale, not a direction error, partly
absorbed by the controller).

## Pointer back to the rank-deficiency finding

[[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]] found the 9-parameter ansatz's $M$ is exactly rank 6
of 9 at the paper's own angles (a 2-qubit pure-state manifold has only $2d-2=6$ real dimensions)
and the 5-parameter ansatz's $M$ is unexpectedly ill-conditioned ($\mathrm{cond}_F\approx863$
against the paper's claimed 66.7). Target C above is that finding turned into a live
demonstration: an adaptive controller running on the rank-deficient ansatz without an explicit,
tuned regularization choice would either lock up (the bug in Assumption 2, before the floor fix)
or silently take a bad step; with cond(M) promoted to an explicit cap and `rcond` an explicit,
logged parameter, it instead visibly throttles down and — with a sensibly chosen cutoff —
recovers to a good energy anyway.

## On how the timestep is been updated

The timestep $\Delta\tau$ is **not** an RKF embedded-error estimate. Route 3 uses a single Euler step plus three physics indicators, combined in `adaptive_qite_run` as a homemade rule (the memo never specifies numbers or this combination).

It starts at $\Delta\tau_0 = 0.05$, and is always clipped to $[10^{-4}, 0.5]$.

Each iteration does three updates, in this order:

**1. Proactive stiffness cap (before the trial)**

After measuring $M$ and solving $M\dot\theta = Y$,
$$\Delta\tau \leftarrow \max\bigl(\Delta\tau_{\min},\; \min(\Delta\tau,\; 2000/\mathrm{cond}(M),\; \Delta\tau_{\max})\bigr)$$
If $M$ is singular, $\mathrm{cond}(M)=\infty$ and this drives $\Delta\tau$ to the floor instead of to zero.

**2. Reactive shrink on energy-gate failure (retry loop)**

The trial is $\theta \leftarrow \theta + \Delta\tau\,\dot\theta$. `energy_check` rejects the step if energy increases (beyond $10^{-9}$) or if $\Delta E$ is outside a band around the **variational** prediction $-Y\cdot\dot\theta\,\Delta\tau$ (absolute $10^{-3}$ or $50\%$ of that prediction). On reject:
$$\Delta\tau \leftarrow \max(\Delta\tau \times 0.5,\; \Delta\tau_{\min})$$
and the trial is redone, up to 20 retries. If retries or $\Delta\tau_{\min}$ are exhausted, the last trial is **committed anyway** (`accepted=False`).

**3. Grow / freeze / shrink for the next step (after commit)**

- If the step was accepted **and** the relative McLachlan residual $\|M\dot\theta - Y\|/\|Y\| \le 0.05$: $\Delta\tau \leftarrow \min(\Delta\tau \times 1.25,\; \Delta\tau_{\max})$.
- If the residual is high: $\Delta\tau \leftarrow \max(\Delta\tau \times 0.5,\; \Delta\tau_{\min})$. Residual does not depend on $\Delta\tau$, so it is **not** retried; it only blocks growth (or shrinks) for the next iteration.
- If accepted but residual is high, the `elif` still shrinks. If rejected but residual is OK, $\Delta\tau$ is left as the (already shrunken) committed value.

The residual is never used to reject the current step; only energy is. Eq. (9)’s noise-floor formula is explicitly not computed here.

### Where does the updating parameters come from?

They are **invented heuristics**, not numbers from the memo. The memo only says things like “cap $\Delta\tau$” and “reject/shrink on violation”; it never gives thresholds. Assumption 1 in `adaptive_qite.py` and the implementation note is explicit: they were **tuned by inspecting Target A’s own trajectory**, never against a reported value.

What each number is, and the (thin) rationale:

| Constant              | Value            | Where it came from                                                                                                                                                                                                                                                                                 |
| --------------------- | ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `DELTA_TAU0`          | 0.05             | Same as the fixed-step baseline; a convenient starting guess                                                                                                                                                                                                                                       |
| `DELTA_TAU_MIN/MAX`   | $10^{-4}$, $0.5$ | Trust region so the controller cannot go to 0 or explode                                                                                                                                                                                                                                           |
| `GROW`                | 1.25             | Classic ODE-controller ratchet (“grow slowly”)                                                                                                                                                                                                                                                     |
| `SHRINK`              | 0.5              | Same family (“halve on reject”)                                                                                                                                                                                                                                                                    |
| `STIFFNESS_CAP_COEFF` | 2000             | Chosen so $\Delta\tau \le 2000/\mathrm{cond}(M)$ only bites when $\mathrm{cond}(M)$ approaches $1/\texttt{RCOND}=10^6$, i.e. where the pseudoinverse is already being regularized. A trial of 50 made the cap dominate the whole 5-param run and it never reached $\tau_{\max}$ in the step budget |
| `RCOND_DEFAULT`       | $10^{-6}$        | Standard `pinv` cutoff, not a paper number                                                                                                                                                                                                                                                         |
| `ENERGY_MONOTONE_TOL` | $10^{-9}$        | FP slack, not a physics tolerance                                                                                                                                                                                                                                                                  |
| `DE_BAND_ABS/REL`     | $10^{-3}$, $0.5$ | Accept band around the variational $\Delta E$ prediction; later notes (Adam QITE / G4) call this band **too loose** as a control setpoint                                                                                                                                                          |
| `RESIDUAL_REL_TOL`    | 0.05             | $5\%$ relative McLachlan residual; freeze-growth threshold, not derived                                                                                                                                                                                                                            |

So $2000/\mathrm{cond}(M)$ is a **dimensional guess** (coefficient / stiffness), and $\times 0.5 / \times 1.25$ is a **hand-tuned multiplicative ratchet**. That is why later work ([[QITE_Hand_Off/Notes_on_Project/Adam_Like_Step_Controller_For_Variational_QITE|Adam_Like_Step_Controller_For_Variational_QITE]]) treats replacing `GROW`/`SHRINK` with a PI/Adam law as the next step — the later controller still **imports** these same thresholds from Route 3 rather than re-deriving them.

The *structure* of the three indicators (cap / energy gate / residual caution) is a synthesis of the memo’s wording. The *numbers* are yours, fitted to one noiseless 5-parameter H2 run.