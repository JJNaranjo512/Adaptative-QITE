---
title: "R — Future routes implementable now"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]]"
---

# R — Future routes implementable now

Companion to [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]]. What can be built next in
`Attack_Route_03/` (or a sibling folder), ordered by value per unit of effort, with what each one
unlocks. Section numbers refer to `Research_Papers/QITE_timestep.pdf`.

> [!note] Precondition
> Modules M1–M6 first. Route 2 will be benchmarked *against* Route 3, so Route 3's numbers have to
> be trustworthy before the comparison means anything. [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting|M6_Honest_Cost_Accounting]] in particular
> supplies the currency all of these are compared in.

## The structural gap these fill

The project currently has **no local truncation error estimate at all.** Route 3 monitors physics
indicators — McLachlan residual, energy monotonicity, $\mathrm{cond}(M)$ — but never estimates the
Euler step's own error $e_n$. Consequently the memo's central argument (Eqs. 8–11: the noise floor
$\Delta\tau_{\text{floor}}$, the optimum $\Delta\tau^\star$, the conclusion that *the controller's
real job on NISQ hardware is to maximize the step*) **cannot be evaluated at all with what is
built**. R1 is the cheapest thing that changes that, which is why it leads.

## R1 — Route 2: step-doubling / Richardson (§3.2, Eq. 12)

**Highest value, lowest effort. Start here.**

Advance once with $\Delta\tau$ and twice with $\Delta\tau/2$ using plain Euler; the difference is the
error estimate,

$$e_n = \frac{\theta^{(\Delta\tau)} - \theta^{(\Delta\tau/2,\times2)}}{2^p - 1}, \qquad p = 1,$$

and the extrapolated combination $2\theta^{(\Delta\tau/2,\times2)} - \theta^{(\Delta\tau)}$ gives a
**second-order update for free**.

Why it is the right first move:

- **Zero change to the measurement pipeline.** It reuses `measure_f` exactly as it stands — the memo
  says so explicitly ("implementable with zero change to the existing QITE measurement pipeline").
- Cost is 3 $f$-evaluations per adapted step, known and constant, so [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting|M6_Honest_Cost_Accounting]]'s
  metric applies directly.
- The memo notes that the two half-steps and the full step share an initial point, so stage-noise
  correlations **partially cancel** in $e_n$, lowering the estimator's noise floor relative to Eq.
  (8). That is a testable claim in its own right.
- It gives the first honest Route-3-vs-embedded-estimator comparison: does the McLachlan residual
  actually track $e_n$? If it does, Route 3's whole premise (cheap indicators substitute for
  expensive stages) is validated; if not, that is the more interesting result.

## R2 — Route 1: embedded pairs (§3.1)

Heun–Euler 2(1) first ($s=2$, $e_n = \tfrac{\Delta\tau}{2}\|k_2-k_1\|$), then Bogacki–Shampine 3(2)
with FSAL (marginal cost 3 evaluations per step). Needs only `measure_f` at a shifted $\theta$ —
mechanically almost the same work as R1, so do them together and compare on Eq. (10).

The memo's own caution is worth keeping in view: hardware noise caps per-step accuracy at
$\varepsilon_{\text{res}}\|\dot\theta\|\Delta\tau$ regardless of integrator order, so **high order
buys nothing below the floor**. Expect R2 to beat R1 modestly at best, and record it honestly if it
does not.

## R3 — PI (Gustafsson) controller, Eq. (17)

$$\Delta\tau_{n+1} = \Delta\tau_n\left(\frac{\text{tol}}{e_n}\right)^{k_I}\left(\frac{e_{n-1}}{e_n}\right)^{k_P},
\qquad k_I \approx \frac{0.3}{p+1},\ k_P \approx \frac{0.4}{p+1}$$

Replaces the invented `GROW = 1.25` / `SHRINK = 0.5` multipliers, which
[[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] Assumption 1 already flags as tuned by eye against
Target A's own trajectory. Implementable the moment R1 or R2 supplies an $e_n$; **needs no shot
noise model**. It also introduces a real `tol` parameter, which Route 3 currently lacks entirely —
and `tol` is what the memo's noise-floor argument constrains.

Do R1 and R3 in the same pass: an error estimator without a principled controller, or a controller
without an estimator, is half of the mechanism.

## R4 — Route 5: regularization coupling (§3.5)

Co-adapt the regularization parameter with the step size. This is the **direct continuation of the
project's own load-bearing finding** — [[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]]'s rank-6-of-9
result, and [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]]'s conclusion that "the cutoff has to be
promoted to an explicit component of the method". Target C demonstrated that picking `rcond` well
matters more than the $\mathrm{cond}(M)$ it produces; R4 is the step from *demonstrating* that to
*handling* it.

Concretely:

1. Swap `np.linalg.pinv(M, rcond=...)` for Tikhonov, $(M + \lambda I)^{-1}Y$, so $\lambda$ is a
   continuous, differentiable knob rather than a discrete SVD cutoff.
2. Minimize the memo's Eq. (8) per-step total over $(\lambda, \Delta\tau)$ jointly, as §3.5 proposes,
   trading the truncation constant $C_{tr}$ against the noise amplification $\mathrm{cond}(M_\lambda)$.

Purely classical, no new quantum cost, and it directly addresses the free parameter the replication
identified as unpinned. §3.5's other two options (Rosenbrock steps from a finite-differenced
Jacobian of recent $\dot\theta$; frozen-$M$ exponential substeps) are also classical and are natural
follow-ons once $\lambda$ is under control.

## R5 — shot noise → Route 4 statistical gating (§3.4, Eqs. 14–17)

Currently blocked by Assumption 7: "No shot noise. The depolarizing channel is applied to the exact
density matrix, so 'noisy' here means biased-but-exact." Route 4's whole premise is that the
accept/reject decision should be a **hypothesis test** against the propagated standard error
$\sigma_e$, not a threshold comparison — which is meaningless without a stochastic estimator.

This is also the only path to evaluating memo Eq. (9) honestly. Target D is explicitly labelled as a
substitute for it, and rightly so, but the substitute is a different quantity: Eq. (9) is a statement
about when the *estimator* stops carrying signal, and it needs both $e_n$ (from R1/R2) and
$\varepsilon_{\text{res}}$ (from a shot model). Larger effort than R1–R4, and it should follow them.

## R6 — scale past 2-qubit H$_2$

Everything so far is one molecule, two qubits, two ansätze. LiH at $N\approx14$ is already the case
[[QITE_Hand_Off/Notes_on_Project/Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors|Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors]] projects failure
for, so it is the natural target. Two things only become visible there:

- Whether the mid-flow stiffness event Target A reports is **generic** or an artifact of a
  5-parameter manifold. As it stands it is $n=1$.
- Where Route 3's $O(1)$ extra-circuit cost actually beats Routes 1/2's $O(N^2)$ per stage. At $N=5$
  the distinction is academic; the memo's argument for Route 3 as "the practical entry point" is an
  argument about large $N$.

## R7 — resolve the Appendix C question with Gaurav

Non-code, cheap, and **open since the replication note was written**.
[[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]] states its own scan space is exhausted (a 60-row sweep over
initial state, gate set, wire order and metric normalization, none within 5 % of the reported
66.7239) and that the remaining unknowns — Appendix C's angles, initial state, gate set, parameter
ordering, and for the 9-parameter case the pseudo-inverse cutoff — have to come from the authors.

It blocks nothing technically, but it is the cheapest open item in the project and it has been
sitting.

## R8 — TITE composition

[[QITE_Hand_Off/Notes_on_Project/Trading_Imaginary_Time_for_Randomness_in_Ground_State_Preparation|Trading_Imaginary_Time_for_Randomness_in_Ground_State_Preparation]] (Arulandu, Martyn & Chuang)
proposes appending randomly-timed real-time evolution after ITE to quadratically suppress trace
distance, allowing $\beta \to \beta/2$. It is implementation-agnostic, so it may compose with the
variational update without touching the $M\dot\theta = Y$ solve.

**Lowest priority.** It is an orthogonal cost lever that does not interact with the step-size
controller at all, and the note currently has no follow-through. Worth keeping on the list precisely
so it is not silently forgotten.

## Suggested order

```
M1 → M4 → M3 → M2 → M5 → M6        (audit modules: make the baseline trustworthy)
      ↓
R1 + R3 together                    (error estimator + principled controller)
      ↓
R2                                  (embedded pairs, compare on Eq. 10)
      ↓
R4                                  (regularization coupling — the project's own open thread)
      ↓
R5, R6                              (shot noise; scale-up)

R7 in parallel, any time — it is an email.
R8 last, or never, depending on where the project goes.
```
