---
title: "G5 — the free error estimator is real, untested, and blind to stiffness"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]"
---

# G5 — The error estimator is validated, and where it fails

**Severity: medium.** Closes G5. Recorded, not fixed.

This is the one finding in this audit that goes in the project's favour, and it is the one the code
never checked.

## The claim

The module's founding observation (README:14-21, Assumption 2) is that Route 3 already computes a
local truncation error estimate and throws it away: `energy_check` forms
$\Delta E_{\text{actual}} - \Delta E_{\text{pred}}$ with
$\Delta E_{\text{pred}} = -(Y\cdot\dot\theta)\Delta\tau$, so the difference is $O(\Delta\tau^2)$ —
an error proxy at **zero extra circuits**, where [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]'s R1 pays 3 $f$-evals per step
for the same thing.

Assumption 2 is careful to flag this as "OUR claim, not the memo's". It is then used as the
foundation of the entire controller and **never tested against anything**.

## The algebra holds

For this pure-state parametrization, $Y_k = -\partial_k E$. Writing $\rho = |\psi\rangle\langle\psi|$
and $\partial_k\rho = |\partial_k\psi\rangle\langle\psi| + |\psi\rangle\langle\partial_k\psi|$:

$$\mathrm{Tr}[\partial_k\rho\,\{H,\rho\}] = 2\,\mathrm{Re}\,\mathrm{Tr}[\partial_k\rho\,\rho H]
= 2\,\mathrm{Re}\langle\psi|H|\partial_k\psi\rangle + 2\,\mathrm{Re}\big(\langle\partial_k\psi|\psi\rangle E\big)$$

The second term vanishes because $\langle\partial_k\psi|\psi\rangle$ is purely imaginary under
normalization, leaving $\mathrm{Tr}[\partial_k\rho\{H,\rho\}] = \partial_k E$, hence
$Y_k = -\partial_k E$ by `vector_Y`'s sign convention. So

$$\Delta E_{\text{pred}} = -(Y\cdot\dot\theta)\Delta\tau = (\nabla E\cdot\dot\theta)\Delta\tau$$

is the **exact** first-order Taylor term of $E(\theta + \Delta\tau\,\dot\theta)$, and the residual
is second order. The claim is sound. (This also depends on the parent's F12/M5 simplifications to
`vector_Y` — dropping the dagger and the $-2\langle H\rangle\rho$ term — which hold here for the
same pure-state reason and are documented at `adaptive_qite.py:272-287`.)

## It matches a real error estimator

V2 re-runs E1's trajectory and, at each accepted step, additionally performs the R1 step-doubling
that [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] proposes: one full step of $\Delta\tau$ against two half-steps of
$\Delta\tau/2$ with a fresh $f$ at the midpoint, read out in energy.

```
   delta_tau       e_hat   |dEa-dEp|    e_E_true    ratio     e_theta    cond(M)
     0.05000   7.359e-02   4.579e-04   2.575e-04    1.778   7.345e-04        152
     0.12207   1.529e-01   5.144e-03   3.428e-03    1.501   6.173e-03       80.1
     0.19073   1.423e-01   1.399e-02   1.258e-02    1.112   1.835e-02       56.8
     0.29802   2.604e-01   6.409e-02   1.081e-02    5.931   9.809e-02       77.5
     0.23283   9.454e-01   9.595e-02   3.311e-02    2.898   1.498e-01        369
     0.31697   7.544e-01   2.270e-02   1.305e-02    1.740   1.360e+00   3.61e+03
     0.00641   2.769e-01   2.769e-04   1.277e-04    2.168   7.316e-03   2.25e+04
     0.01556   1.772e-02   1.772e-05   1.549e-06   11.439   4.721e-03   2.92e+04
     0.02781   2.693e-02   2.693e-05   1.285e-05    2.096   9.778e-05   5.01e+04
     0.01475   4.644e-03   4.644e-06   2.266e-06    2.050   2.112e-05   1.36e+05
     0.01035   1.919e-03   1.919e-06   9.430e-07    2.035   9.368e-06   1.93e+05
     0.00800   1.008e-03   1.008e-06   4.975e-07    2.027   5.250e-06    2.5e+05

  Pearson  corr(log e_raw, log e_E_true) = 0.9902
  Spearman corr(e_raw, e_E_true)          = 0.9882
  ratio |dEa-dEp| / e_E_true: median=2.060  min=0.350  max=11.439
```

**Probe count, recorded late.** The three aggregate statistics above are computed over the **first 40
accepted steps** of the E1 trajectory, and the twelve rows shown are a subsample of those 40. V2 did
not record the count, which is a defect of the same class as §6.2 in
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility|G6_Attribution_And_Reproducibility]] — the numbers are correct but the note alone does not pin them
down, because the statistics are *monotone in the probe count*: 20 probes gives Pearson 0.9701 and
median 2.157, 60 gives 0.9948 and 2.026, and the full 190-step trajectory gives 0.9983 and 2.005. The
count is now pinned as the default in `verify_audit.py --target g5`, which prints all 40
rows and the sensitivity, so this table is regenerable. The three sensitivity pairs quoted in this
paragraph are themselves in the script's drift ledger — they were printed but
unchecked until an independent pass noticed it, which is the same defect class as a missing row:
a number the ledger does not carry is compared against nothing. The full-trajectory pair needs
`--g5-full-trajectory`, and the script says so where it skips them. The trend is in the estimator's favour: the
more of the converged tail is included, the closer the ratio sits to the predicted 2.

**The ratio settles at 2.03–2.09 over the whole converged tail**, and that is exactly the predicted
value. Step doubling has full-step error $\approx C\Delta\tau^2$ and two-half-step error
$\approx 2C(\Delta\tau/2)^2 = C\Delta\tau^2/2$, so their difference is $C\Delta\tau^2/2$ — half of
what $|\Delta E_{\text{actual}} - \Delta E_{\text{pred}}| \approx C\Delta\tau^2$ measures. A ratio
of 2, recovered to within 3 %.

So the premise is not just defensible, it is **quantitatively confirmed**: the discarded quantity
is the same local truncation error R1 would pay three $f$-evaluations per step to obtain, up to a
factor of 2, at zero extra circuits. This is a genuine result. It is currently asserted in an
`ASSUMPTIONS` string and nowhere demonstrated.

## Where it breaks: stiffness

The ratio is not constant everywhere. It swings from 0.35 to 11.4, and the excursions are not
random — they cluster on the steps with large $\Delta\tau$ and rapidly changing $\mathrm{cond}(M)$
(rows 4–8 above: ratios 5.93, 2.90, 1.74, 2.17, and the 11.44 outlier at
$\mathrm{cond}(M) = 2.9\times10^4$).

The reason is structural. What $\hat e$ measures is
$\tfrac12\Delta\tau^2\,\dot\theta^\top\nabla^2E\,\dot\theta$ — the curvature of $E$ along a
**straight line** in $\theta$. The Euler step's true local error is
$\nabla E\cdot\tfrac12\Delta\tau^2\ddot\theta$ with
$\ddot\theta = (\partial_\theta f)f$ and $f = M^{+}Y$. Expanding $\partial_\theta f$ produces a
$\nabla^2 E$ term — which is what $\hat e$ captures — **plus a term in $\partial_\theta M^{+}$**,
which it does not. In words: $\hat e$ is blind to the metric $M$ changing along the step.

That blind spot is exactly the stiffness the $\mathrm{cond}(M)$ cap exists to handle, which
predicts that the Adam law cannot subsume the cap — and
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility|G6_Attribution_And_Reproducibility]] shows E3's cap-off arms behaving in just that way. The
$\partial_\theta M^{+}$ expansion above is sketched, not derived in full here; it is recorded as
the mechanism the numbers are consistent with, not as a proven identity.

## Proposed work

1. **Promote V2 into the script as a test.** A `check_lte_proxy()` in the style of the parent's
   `check_eq_b9_self()` (`adaptive_qite.py:290`), asserting the ratio-of-2 against step doubling on
   a handful of steps, converts the module's founding assumption into something the code verifies.
   This is the single highest-value change in this audit: it is cheap, it is a real result, and
   right now the strongest claim in the folder is the least supported.
2. **Report the ratio as a diagnostic**, alongside $\mathrm{cond}(M)$. Its departure from 2 is a
   direct measurement of how much of the local error the estimator is missing — a stiffness
   indicator in its own right, and a cheaper one than the cap's.
3. **Write the result up.** "A zero-cost local truncation error estimate that agrees with
   step-doubling to a constant factor" is the kind of finding
   [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] exists to record, and it currently appears nowhere
   outside a docstring. It also partially closes the gap
   [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]] §1 identifies as the project's first structural problem
   ("the project has **no local truncation error estimate at all**"), which would change R1's
   priority in [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]].

**Falsifiable prediction — tested, and it FAILED.** If the blind spot is $\partial_\theta M^{+}$,
then the ratio's departure from 2 should correlate with the *step-to-step change* in
$\mathrm{cond}(M)$ rather than with its level. V2 showed the clustering but never tested that
correlation. `verify_audit.py --target g5` now does, and the prediction does not survive:

```
  Spearman(|log(ratio/2)|, |dlog cond(M)|) = +0.3176   <- change, as predicted
  Spearman(|log(ratio/2)|, cond(M))        = -0.8938   <- level
```

The *level* correlates nearly three times as strongly as the *change*, which is the opposite of
what the mechanism predicts. Before reading that as a refutation, though, the confound has to be
named: $\mathrm{cond}(M)$ rises almost monotonically along the run (Spearman against step index
$= +0.9786$) while the ratio settles onto 2, so step index alone predicts both, and the raw
correlations cannot separate mechanism from convergence. Controlling for it with a rank partial
correlation:

```
  partial(|log(ratio/2)|, |dlog cond(M)| | step) = +0.1523   <- change
  partial(|log(ratio/2)|, cond(M)        | step) = -0.5437   <- level
```

The conclusion holds after the control. **So the prediction is not supported.** What this refutes
is narrow but real: the clustering in the table above is *not* evidence for the
$\partial_\theta M^{+}$ mechanism specifically, because a monotone stiffness level explains the
departures better than the step-to-step metric change does. The algebra in "Where it breaks" is
unaffected — the missing $\partial_\theta M^{+}$ term is still absent from $\hat e$, and that
remains a derivation, not a measurement. What fails is my attempt to *detect* it this way. A probe
that can distinguish the two needs to hold the convergence stage fixed — for instance by comparing
steps at matched $\Delta\tau$ and matched $\mathrm{cond}(M)$ level but different
$|\Delta\log\mathrm{cond}(M)|$, which this trajectory does not sample densely enough to support.
That is open work, and it supersedes item 2 of the proposed work above as the sharper question.

Unlike the three refutations recorded elsewhere in this audit, this prediction was not a
reading-stage guess — it was stated in this note, *after* the evidence was in, and then refuted by
the verification pass. Which is the argument for shipping the script: a note that only asserts a
correlation gets to keep it.
