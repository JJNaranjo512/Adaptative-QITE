---
title: "G3 — E2 has no noise to average and cannot resolve its own signal"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]"
---

# G3 — The noise experiment cannot resolve noise

**Severity: critical.** Closes G3. Recorded, not fixed.

E2 is described in the code as "the load-bearing test: does EMA smoothing help under noisy dE
readings?" (`adam_qite.py:500`) and it returns the suite's strongest verdict — "**The premise
fails.**" Three independent problems mean it cannot support that verdict, or its opposite.

## 3.1 There is no stochastic noise in the experiment

Assumption 1 inherits the parent's Assumption 7 verbatim: "**no shot noise, so 'noisy' means
biased-but-deterministic**". `make_noisy_rho_fn(p)` applies a global depolarizing channel to the
density matrix analytically. For a fixed $\theta$ it returns the same $\rho$ every time.

So the quantity the EMA is supposed to smooth **does not fluctuate**. The premise being tested —
"the ratchet reacts to every individual noisy $\Delta E$; the EMA averages them" (`:506-509`) —
describes a stochastic estimator, and the experiment supplies a deterministic, smooth,
$p$-dependent bias. An exponential moving average cannot remove a bias; averaging a constant
returns the constant.

This is not a subtle point about the noise model — it is stated plainly in the script's own
assumption list, three lines above the target that depends on it being false. Whatever $s_e$ is
measuring in E2, it is not measurement noise; with $p$ fixed and the channel deterministic, the
dispersion in $\hat e$ comes entirely from the trajectory's own curvature as it converges, the
same thing it measures at $p=0$.

**A controller that gates on signal-to-noise cannot be evaluated in a suite with no noise
realizations.** That is R5's gap, and the README is explicit that Adam was meant to answer the
question "with no shot-noise model (R5)". It cannot: the question is only well-posed once R5
exists.

## 3.2 Both controllers have the same noise response

If E2 did measure noise robustness, the two arms would respond differently to $p$. They do not.
From `Output_Pre_Audit.txt`, dividing each RMS by its $p$ (V3):

```
         p   Adam RMS/p     R3 RMS/p
    0.0001       6.5312       6.3761
    0.0010       6.5215       6.3646
    0.0100       6.4059       6.2527
    0.0500       5.8195       5.7688
```

Both deviations are **linear in $p$ with the same coefficient** across three decades — the
signature of a residual first-order bias left over after the Eq. (7) debiasing, propagating
identically through both controllers. Correspondingly the published ratio is constant:

```
  published ratio  : mean=1.0206  spread=0.0159
```

A genuine robustness difference would vary with $p$; a spread of $\pm0.016$ over three decades says
the 2.4 % is a fixed offset, not a response.

## 3.3 The offset is the stalled baseline, and the signal is smaller than it

Each arm is scored against **its own** $p=0$ run (`:517-518`, `:538-539`). Those two references are
not equivalent:

```
  p=0 Adam    : 500 steps, tau_end=2.8615, E=-1.123002, exit=max_steps
  p=0 Route 3 : 234 steps, tau_end=3.0097, E=-1.123745, exit=tau_max
  --> the two baselines differ in tau reach by 0.1482 and in final E by 7.430e-04
```

So the two RMS values are integrals over **different $\tau$ windows against different reference
trajectories**, and their ratio is not a comparison. Worse, at the two smallest $p$ the thing being
measured is smaller than the thing contaminating it: the deviation at $p=10^{-4}$ is
$6.5\times10^{-4}$, while the gap between the two baselines is $7.4\times10^{-4}$. **E2 has no
resolving power at $p \le 10^{-3}$.**

Re-scoring both arms against a single shared baseline on the common window (V3):

```
         p | --- as published (own baseline) ---  | --- shared baseline (Route 3 p=0) ---
           |     RMS Adam      RMS R3     ratio |     RMS Adam      RMS R3     ratio
    0.0001 |    6.531e-04   6.376e-04    1.0243 |    1.010e-03   6.374e-04    1.5850
    0.0010 |    6.521e-03   6.365e-03    1.0246 |    6.627e-03   6.363e-03    1.0414
    0.0100 |    6.406e-02   6.253e-02    1.0245 |    6.265e-02   6.249e-02    1.0026
    0.0500 |    2.910e-01   2.884e-01    1.0088 |    2.805e-01   2.884e-01    0.9727

  published ratio  : mean=1.0206  spread=0.0159
  shared-baseline  : mean=1.1504  spread=0.6123
```

The stable 1.024 dissolves. At the two largest $p$ — the only ones where the signal exceeds the
contamination — **Adam is level or better** (1.003, 0.973). The direction of the published verdict
is not robust to the choice of reference.

And by [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G2_Beta2_Horizon_Dominates_Every_Verdict|G2_Beta2_Horizon_Dominates_Every_Verdict]], the stalled Adam baseline that creates the
offset is itself an artifact of mismatched EMA horizons. With them matched, the Adam baseline
reaches $\tau=3.0093$ and $E=-1.123736$ — a $9\times10^{-6}$ gap to Route 3 instead of
$7.4\times10^{-4}$, i.e. the confound shrinks by two orders of magnitude and E2 would have to be
re-run from scratch.

## Relation to F4 / M3

[[QITE_Hand_Off/Audits/Audit_Project_Wide/M3_Trajectory_Comparison_Metric|M3_Trajectory_Comparison_Metric]] fixed the parent's habit of extrapolating a stalled run past
its own $\tau$ range, and added `exit_reason` and `coverage` precisely so a stall would be visible
rather than silently scored. E2 uses that machinery correctly — `trajectory_deviation` is called as
intended, and the table does print `STALLED(A)`.

The defect is one level up: **M3 made a stall visible within one comparison; E2 then took a ratio
of two comparisons, one of which contains a stall, and reported the ratio.** The flag is printed
two columns away from the number it invalidates, and the "Ratio RMS_Adam / RMS_Route3" block below
it carries no flag at all. The spirit of F4 is reintroduced in a form M3's fix does not catch.

## Proposed work

1. **Do not report E2's verdict** until it is redesigned. The honest current statement is *"the
   premise is untested"*, not *"the premise fails"*.
2. **Score both arms against one reference**, on a window both cover. If a per-controller
   reference is wanted for a different reason, report both columns, as V3 does.
3. **Require resolving power.** A comparison is only meaningful where the measured deviation
   exceeds the baseline-to-baseline gap. At minimum, print that gap next to the RMS column so a
   reader can see when the experiment has nothing to say.
4. **Blocked on R5.** The premise E2 was built for needs shot noise — repeated runs at fixed $p$
   with different shot seeds, and a deviation measured *across realizations*, not against a clean
   trajectory. Until then, an EMA on a deterministic signal has nothing to average, and this is the
   experiment [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] should sequence before the controller is judged.

**Falsifiable prediction.** Under R5 with genuine shot noise at fixed $p$, the ratchet's
$\Delta\tau$ should show variance across realizations that the EMA visibly reduces, and the
published ordering at small $p$ should reverse. If the EMA still does not help once there is real
noise to average, *that* is the result the README currently claims — and it would then be
supportable.
