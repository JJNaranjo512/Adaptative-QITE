---
title: "M4 — Controller behaviour near convergence"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]]"
---

# M4 — Controller behaviour near convergence

> [!success] Implemented — vindicated the original Target A claim rather than overturning it, see [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE#6. Implementation status|hub §6]]
> Gate switched to $-Y\cdot\dot\theta\,\Delta\tau$; exact-ITE kept as a diagnostic; reason
> histogram and cap-vs-band attribution added. Re-run: over the last 43 steps of Target A the
> $\mathrm{cond}(M)$ cap binds on 31 and the energy band triggers 0 retries — the $\Delta\tau$
> collapse is genuinely the stiffness cap, not the mis-specified gate. Unregularized
> $\sigma_{min}$/spectral-cond added and independently confirm real ill-conditioning ($1.4\times10^6$).
> Target C's numbers changed the most (see M6/hub §6) — both `rcond` choices now converge.

**Findings:** F5, F6 (high), F10 (medium)  **Effort:** M  **Blocks:** Target A's headline; forces a
re-run of everything.

## What is wrong

### F5 — the accept/reject gate uses the wrong linear model

`energy_check` (`adaptive_qite.py:302-312`) predicts

$$\widehat{\Delta E} = -2\,\mathrm{Var}_\rho(H)\,\Delta\tau$$

and rejects when the measured $\Delta E$ falls outside a band around it. That is memo Eq. (13) — the
identity obeyed by **exact** imaginary-time flow. The *variational* flow on a finite ansatz obeys

$$\frac{dE}{d\tau} = \nabla E \cdot \dot\theta = -Y\cdot\dot\theta = -\|\dot\theta\|_M^2,$$

and the two agree only to the extent that the ansatz can represent the exact flow. Their difference
is precisely the ansatz expressibility deficit — which is exactly what the McLachlan residual is
*already* there to measure, so the gate is double-counting it.

The deficit does not vanish as the run converges. It grows in relative terms. At the 5-parameter
ansatz's converged point, $E - E_{\text{exact}} = 2.59\times10^{-2}$ (with
$E_{\text{exact}} = -1.145599124$, recomputed from `h2_hamiltonian()` during this audit), so
$\mathrm{Var}_\rho(H)$ stays $O(10^{-2})$ while $\Delta E_{\text{actual}} \to 0$ as
$\|Y\|\to0$. Order-of-magnitude at $\Delta\tau = 0.3$ near convergence:
$|\Delta E_{\text{actual}} - \widehat{\Delta E}| \approx 2.5\times10^{-2}$ against a band of
$\max(10^{-3},\,0.5\,|\widehat{\Delta E}|) \approx 1.2\times10^{-2}$ — **permanent rejection.**

So the controller shrinks $\Delta\tau$ near convergence. The memo says the opposite is the point:

> adaptivity is valuable not to make steps small where the flow is fast, but to make them **as
> large as noise permits where the flow is slow** (e.g. near convergence, where $\|Y\|\to0$ and
> $\dot\theta$ flattens). — §3.1

The controller as built is fighting its own design goal in exactly the regime the memo singles out.

**Corollary.** Target A's $\Delta\tau$ collapse at $\tau\approx2.5$ now has two candidate causes —
the $\mathrm{cond}(M)$ cap and this band — and the note attributes it entirely to the first, reading
it as a physical "stiffness event". That attribution is not established.

### F6 — "zero violations of Eq. (13)" is narrower than stated

`energy_check` returns two distinct reason strings: `"energy increased"` (a genuine monotonicity
violation) and `"dE off linear prediction"` (the band test). The note reports "zero violations of
Eq. (13), satisfying the plan's verification criterion" — which is a statement about the first only.
The second is never counted or printed, and by F5 it is likely nonzero and growing along the run.
`log["reasons"]` already stores this per step; nothing reads it back.

### F10 — cond(M) is reported through the regularizer

`frobenius_cond(M, rcond=...)` (`adaptive_qite.py:241-258`) computes
$\|M\|_F\|M^+\|_F$ with the cutoff already applied. That bounds the reported value: with
`RCOND_DEFAULT = 1e-6` and $N=5$, $\|M^+\|_F \le \sqrt{5}/(10^{-6}\sigma_{max})$, so
$\mathrm{cond}_F \lesssim \sqrt{5}\times10^{6}$. "cond(M) climbs past $10^5$ / beyond $10^6$" is
therefore in part a statement about the cutoff saturating, not about $M$ becoming singular.

The three documents also disagree on the number: [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] and
[[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] say "past $10^5$", `Attack_Route_03/README.md` says
"beyond 1e6".

## Why it matters

Target A's claim — adaptive control detects a genuine mid-flow stiffness event and throttles
correctly — is the result that [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §4 has already been
updated to carry, as the answer to its own open question 1. If the $\Delta\tau$ collapse is instead
(or additionally) the energy band misfiring near convergence, that update overstates what was shown.
This is the finding most likely to change a conclusion already written into the system of record.

## Change to make

1. **Use the variational prediction for the gate.** Predict $\widehat{\Delta E} = -Y\cdot\dot\theta\,
   \Delta\tau$ and test the measured $\Delta E$ against *that*. This is a true first-order model of
   the step actually being taken, so its failure means what the memo wants it to mean: the linear
   model of the step is breaking down, shrink $\Delta\tau$.
2. **Keep $-2\mathrm{Var}(H)\Delta\tau$ as a separate diagnostic**, reported but not gating. The
   gap between the two is the ansatz deficit, and tracking it along the trajectory is genuinely
   interesting — it is the same quantity the McLachlan residual probes, from the energy side.
3. **Keep the monotonicity test as its own gate.** $\Delta E > 0$ is a real violation at any
   $\Delta\tau$ and should still reject.
4. **Print the reason histogram** from the existing `log["reasons"]`, split by reason string, and
   report the two counts separately in the note.
5. **Report $\sigma_{min}$ and the unregularized spectral cond alongside the regularized
   $\mathrm{cond}_F$** so a stiffness claim rests on a quantity the cutoff does not bound.
6. **Attribute the $\Delta\tau$ collapse.** With the reason histogram and the cap value logged per
   step, it is directly readable whether the cap or the band drove $\Delta\tau$ down at
   $\tau\approx2.5$. Log which of the two was binding.

## How to verify

- Re-run Target A. Prediction: with the variational gate, $\Delta\tau$ **grows** toward the end of
  the run rather than collapsing, matching the memo's §3.1 expectation. If it still collapses and
  the histogram shows the $\mathrm{cond}(M)$ cap binding, then the stiffness event is real and the
  note's current reading is vindicated — record that outcome either way.
- Check that the reported $\mathrm{cond}_F$ ceiling moves when `RCOND_DEFAULT` changes but
  $\sigma_{min}$ does not. That separates the regularizer's contribution from $M$'s.
- Reconcile the $10^5$ / $10^6$ discrepancy against whichever number the re-run produces.

## Documents to update

- [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] — Target A paragraph, Assumption 1 and 2.
- `Attack_Route_03/README.md` — Target A row, including the "beyond 1e6" figure.
- [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §4, the extension paragraph on cond(M) along a trajectory, and §6.
