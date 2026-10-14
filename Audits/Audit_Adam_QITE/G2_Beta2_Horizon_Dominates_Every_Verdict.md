---
title: "G2 — the two EMAs use mismatched horizons, and that is what E2 measured"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]"
---

# G2 — Mismatched EMA horizons

**Severity: critical.** Closes G2. Recorded, not fixed.

## The defect

The confidence ratio is the whole idea of the controller (`adam_qite.py:204`):

$$s_e = \frac{\hat m_e}{\sqrt{\hat v_e} + \varepsilon}$$

For this to be the signal-to-noise ratio the design claims, both moments must describe the **same
window** of the signal. The code gives them different ones (`:99-100`):

$$\beta_1 = 0.9 \Rightarrow \text{10-step memory}, \qquad
\beta_2 = 0.999 \Rightarrow \text{1000-step memory}$$

on runs of **172–500 steps**. So $s_e$ divides a ten-step local mean by a global RMS taken over a
window longer than the entire experiment, applied to a signal that falls four orders of magnitude
as the run converges. It is not a signal-to-noise ratio; it is a measure of how far the run has
converged since it started. An early transient spike enters $\hat v_e$ and is never forgotten,
because the run ends before the EMA does.

The same mismatch is what breaks Assumption 4's stated bound. $\hat m/\sqrt{\hat v} \le 1$ is
Jensen's inequality, and Jensen needs **one** weighting. With two, the bound simply does not apply
— which is why the script's own output reports $s_e$ up to 1.015, and why a reversed mismatch
($\beta_1=0.999,\beta_2=0.9$) produces a mean $s_e$ of **258**. One root cause, two symptoms; see
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G7_Latent_Defects|G7_Latent_Defects]] for the second.

These are the correct stock values for what Adam was built for — minibatch SGD over $10^5$–$10^6$
steps, where a 1000-step window is genuinely local. They are out of range by three orders of
magnitude here. The README's fairness protocol ("stock Adam defaults, **not tuned here**") is
right not to fit them to the result; but never *varying* a parameter is not evidence that it does
not matter.

## Evidence

**V9 — is it $\beta_2$'s value, or the mismatch?** `SEED_D`, E2's stalling baseline:

```
  configuration                       steps   tau_end      final E  mean s_e  frozen       exit
  Route 3 (ratchet)                     234    3.0097    -1.123745       n/a     n/a    tau_max
  b1=0.9 b2=0.999                       500    2.8615    -1.123002    0.0893     72%  max_steps
  b1=0.999 b2=0.999  <-- matched        254    3.0269    -1.123757    0.5170      0%    tau_max
  b1=0.99 b2=0.99    <-- matched        259    3.0226    -1.123736    0.4387      0%    tau_max
  b1=0.9 b2=0.9      <-- matched        244    3.0093    -1.123736    0.3938      0%    tau_max
  b1=0.999 b2=0.9                       242    3.0077    -1.123740  258.1220      0%    tau_max
```

**Every matched pair completes. Only the mismatched pair stalls.** The value does not matter —
$\beta=0.9$, $0.99$ and $0.999$ all work, including $\beta_2=0.999$ itself, the very value the
stock configuration uses. Growth freezes on 0 % of steps in every matched run against 72 % in the
shipped one, and $s_e$ recovers from $0.089$ to $0.39$–$0.52$.

This also falsifies the heavy tail as the *sufficient* cause. The signal is equally heavy-tailed in
all six rows — same seed, same trajectory, same $\hat e$ distribution — and with matched horizons it
causes no freezing at all.

**V10 — in M6's matched-circuit currency. Two of the three seeds are shown; see the `SEED_C`
caveat below:**

```
  SEED_A 5-param  (E_floor=-1.124600, matched budget=2069 circuits)
  arm                                           steps   circ     E-floor    E@budget     s_e  froz
  Route 3 (ratchet)                               172   2069  +1.607e-03  +1.607e-03     n/a   n/a
  Adam b1=0.9  b2=0.999  (stock, mismatched)      190   2285  +1.748e-03  +1.860e-03   0.140   49%
  Adam b1=0.9  b2=0.9    (matched)                179   2151  +1.190e-03  +1.800e-03   0.200    0%
  Adam b1=0.99 b2=0.99   (matched)                184   2212  +1.437e-03  +1.841e-03   0.355    0%
  Adam b1=0.999 b2=0.999 (matched)                181   2176  +1.447e-03  +1.803e-03   0.407    0%

  SEED_D 5-param  (E_floor=-1.124600, matched budget=2816 circuits)
  Route 3 (ratchet)                               234   2816  +8.555e-04  +8.555e-04     n/a   n/a
  Adam b1=0.9  b2=0.999  (stock, mismatched)      500   6006  +1.598e-03  +4.165e-03   0.089   72%
  Adam b1=0.9  b2=0.9    (matched)                244   2934  +8.641e-04  +3.495e-03   0.394    0%
  Adam b1=0.99 b2=0.99   (matched)                259   3114  +8.637e-04  +3.942e-03   0.439    0%
  Adam b1=0.999 b2=0.999 (matched)                254   3054  +8.430e-04  +3.893e-03   0.517    0%
```

## What this overturns, and what it does not

**Overturned.** README:66's "Adam's own $p=0$ baseline **stalls**: 500 steps reaching only
$\tau=2.86$, against Route 3's 234 steps reaching $\tau=3.01$" is not a property of the Adam
controller. It is a property of running it with two different EMA horizons. Matched, it finishes in
244 steps at $\tau=3.0093$ within $9\times10^{-6}$ of Route 3. Everything E2 concludes downstream
of that stall inherits the artifact, and [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G3_Noise_Experiment_Cannot_Resolve_Noise|G3_Noise_Experiment_Cannot_Resolve_Noise]] shows the
stall is in fact *what the published 1.024 ratio is measuring*.

**Not overturned — strengthened.** E1's "a wash" verdict survives every configuration tested.
**At matched circuit budget Route 3 wins against all four Adam configurations at both seeds shown
above**, without exception. `SEED_C` was swept too and agreed, but its rows are omitted for the
`rcond` reason below and — the scratchpad being gone — are no longer reproducible; the regenerable
claim is the two-seed one, and that is what `verify_audit.py --target g2` re-derives.
$\beta_1=\beta_2=0.9$'s better *final* energy at `SEED_A`
($+1.190\times10^{-3}$ vs $+1.607\times10^{-3}$) is bought with 4 % more circuits and vanishes
under the matched-budget column — exactly the trap [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting|M6_Honest_Cost_Accounting]] was written to
catch. The README's overall conclusion, that the Adam machinery does not earn its place on these
problems, is therefore **correct**; the reasoning it gives for E2 is not.

**Caveat on the 9-parameter arm.** V10's `SEED_C` rows use `RCOND_DEFAULT = 1e-6`, which is *not*
one of the two `rcond` values E3 actually runs, so they are not comparable to the E3 table and are
omitted here. The `rcond`-resolved version of this sweep is open work.

## Proposed work

1. **Tie the two horizons together.** Either set $\beta_1 = \beta_2$, or derive both from the
   expected step count ($\tau_{\max}/\overline{\Delta\tau}$, or `MAX_STEPS`). Deriving a decay from
   the problem is a principled choice that can be defended in `ASSUMPTIONS`; it is not tuning, and
   it should be made *before* looking at the resulting energies.
2. **Sweep $\beta_1$ and $\beta_2$ jointly across all three seeds** in matched-circuit currency.
   This is the ablation the README already lists as open; G2's finding is that it is not optional,
   because it is what the E2 verdict actually rested on.
3. **Report $\hat m_e$ and $\sqrt{\hat v_e}$ separately** in `adam_signal_diagnostic`, not only
   their ratio. The ratio hides which moment is misbehaving — the reason this needed a sweep to
   find, when the shipped diagnostic prints $s_e$ on every run.
4. **Assert the bound you claim.** If Assumption 4 says $s_e \in (0,1]$, the controller should
   check it and fail loudly, per [[QITE_Hand_Off/Audits/Audit_Project_Wide/M2_Noise_Model_Falsifiability|M2_Noise_Model_Falsifiability]]'s pattern of turning a stated
   assumption into something the script verifies. It would have caught this on the first run.

**Falsifiable prediction, tested.** Before running V9 this note predicted that matching the two
horizons at *any* common value would remove the stall, and that if the stall survived
$\beta_1=\beta_2=0.999$ the heavy tail would be the whole story. It did not survive: matched
$0.999$ completes in 254 steps. The prediction held.
