---
title: "M6 — Honest cost accounting"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]]"
---

# M6 — Honest cost accounting

> [!success] Implemented and confirmed — headline claim survives, see [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE#6. Implementation status|hub §6]]
> Circuit-count accounting (`n_circ`, `s_eff`) added to both runners; the closed-form check holds
> exactly (12 circuits/step $=2N+2$ for the fixed-step baseline). At the corrected, matched
> circuit budget (732), adaptive still beats fixed: $-1.1210$ vs. $-1.1167$. Variational floor
> added ($E_{floor}=-1.1246$, multi-restart `scipy.optimize.minimize`, plotted as a reference
> line on Target A's figure) — adaptive lands $5\times$ closer to it than fixed
> ($1.6\times10^{-3}$ vs. $7.9\times10^{-3}$).

> [!warning] Hand-off correction
> The $5\times$ divides end points at different budgets: Route 3 after 2069 circuits ($1.607\times10^{-3}$ above $E_{floor}$) against fixed Euler after 732 ($7.909\times10^{-3}$). At a matched 732 circuits Route 3 is at $-1.120952$ (corrected read-out, G6 §6.3), $3.648\times10^{-3}$ above the floor, so $2.2\times$ closer. The status report states the matched figure.

**Findings:** F9 (medium)  **Effort:** M  **Blocks:** any comparison of Route 3 against Routes 1/2
in [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] — this module supplies the currency.

## What is wrong

Target A's headline is that adaptive control beats the fixed-step baseline **at matched cost**:

```python
n_match = min(adaptive["n_evals"], fixed["n_evals"])
E_adapt_at_budget = adaptive["E"][n_match] ...
```

`n_evals` counts calls to `measure_f` only (`adaptive_qite.py:336`). But the adaptive run also pays,
per retry, an extra `rho_fn(theta_trial)` plus an `energy_and_variance` evaluation
(`adaptive_qite.py:349-350`) — and those are quantum measurements on hardware, the very thing the
budget is supposed to be counting. The fixed-step baseline never retries and never pays them. So the
adaptive run's true cost is systematically undercounted, by an amount that grows precisely in the
regime where the note reports the adaptive run doing well (the retry-heavy stretch near
$\tau\approx2.5$).

Separately, the two final energies are $-1.119671$ (adaptive) and $-1.116691$ (fixed), i.e. errors
of $2.59\times10^{-2}$ and $2.89\times10^{-2}$ against $E_{\text{exact}} = -1.145599124$ (recomputed
from `h2_hamiltonian()` during this audit). Both sit close to the **variational floor of the
5-parameter ansatz**, which the script never computes. Without that reference line it is not
possible to say whether the $3\times10^{-3}$ gap between the two runs is the integrator doing better
or two runs converging onto the same wall from slightly different distances.

## Why it matters

"Adaptive wins at equal cost" is Target A's entire claim, and it is the claim that justifies Route 3
as an approach. It is also the number that Routes 1 and 2 will be measured against — and those
routes cost 3 evaluations per step by construction, so a comparison that undercounts Route 3's
retries will be biased against them from the start. Getting the currency right before building R1/R2
is worth more than getting it right afterwards.

## Change to make

1. **Count everything.** Add a counter for every call that would be a circuit on hardware: each
   `measure_f` (which is itself $2N+1$ density-matrix evaluations via central differences), each
   trial-step energy evaluation, each retry. Report the breakdown, not one number.
2. **Adopt the memo's Eq. (10).** The memo already defines the right metric,

   $$C(\Delta\tau) = \frac{s_{\text{eff}}\cdot n_{\text{circ}}(N)\cdot n_{\text{shots}}}{\Delta\tau},$$

   quantum cost per unit imaginary time. Implement it with $n_{\text{shots}}$ symbolic for now
   (there is no shot model yet — see R5 in [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]) and $s_{\text{eff}}$ the measured
   effective stages per accepted step, which for Route 3 is $1 + \text{retries}$. This makes Route 3
   and Routes 1/2 directly comparable on the memo's own terms.
3. **Add the variational floor.** Minimize $E(\theta)$ directly over the 5-parameter ansatz
   (`scipy.optimize.minimize` on the existing `rho_fn`, several restarts) and plot it as a horizontal
   reference on the Target A convergence figure. Then the claim becomes "adaptive reaches the floor
   in fewer evaluations", which is a statement about the integrator, or it becomes "both are at the
   floor and the difference is noise", which is also worth knowing.
4. **Report final $\tau$.** The script prints it but neither the note nor the README records it, so
   it is currently impossible to tell from the documents whether the two runs were compared over the
   same imaginary-time horizon at all.

## How to verify

- The counted cost of the fixed-step run must equal $60 \times (2N+1+1)$ exactly for
  $\Delta\tau=0.05$, $\tau_{max}=3$ — a closed-form check on the counter itself.
- Re-run Target A with the full accounting. If adaptive still wins, the claim stands on a firmer
  basis; if it does not, the note says so. Do not adjust `DELTA_TAU0` for the baseline to recover
  the result.
- The variational floor must be $\le$ both final energies. If either run lands *below* it, the
  minimizer is stuck and needs more restarts.

## Documents to update

- [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] — Target A paragraph, Assumption 5's
  "evaluation count counts calls to `measure_f`" sentence.
- `Attack_Route_03/README.md` — Target A row, the "61 evals" figure.
- [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §6.
