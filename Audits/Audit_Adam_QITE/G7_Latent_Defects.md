---
title: "G7 — defects that did not fire in this suite but are live in the code"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]"
---

# G7 — Latent defects

**Severity: low.** Closes G7. Recorded, not fixed.

Four defects that are real in the source and that **no experiment in this suite exercises**. They
are collected here rather than spread across G1–G6 so the severity ranking stays honest: none of
them changes a reported number. Two of them would, on a different problem.

## 7.1 The missing `accepted` guard

The parent grows the step only under two conditions (`adaptive_qite.py:493-496`):

```python
if accepted and rel_residual <= RESIDUAL_REL_TOL:
    delta_tau = min(delta_tau * grow, delta_tau_max)
elif rel_residual > RESIDUAL_REL_TOL:
    delta_tau = max(delta_tau * shrink, delta_tau_min)
```

`adam_qite.py:310` applies growth unconditionally:

```python
delta_tau = float(np.clip(delta_tau * state["growth"], delta_tau_min, delta_tau_max))
```

The parent's F8/M5 policy commits a trial step once `max_retries` or `delta_tau_min` is exhausted,
recording `accepted=False` (`:271`). Such a step is the worst step the run produced — and under
`adam_qite.py:310` it is **rewarded with growth**. The same missing guard makes `ctrl.update` fire
at `:283` on a step the gate refused, which contradicts the module's own Assumption 7 verbatim:
"The EMAs update once per **ACCEPTED** step, using the accepted trial's `e_hat`, and the bias
correction counts accepted steps."

**It never fires in this suite.** V4 instruments every Adam run in E2 and E3:

```
  scenario                          steps   forced   retried_ok        exit
  E2 p=0.0001  Adam                   500        0            3   max_steps
  E2 p=0.001   Adam                   500        0            3   max_steps
  E2 p=0.01    Adam                   449        0            3     tau_max
  E2 p=0.05    Adam                   381        0            3     tau_max
  E3 rcond=1e-08 Adam (cap on)        108        0            0     tau_max
  E3 rcond=1e-08 Adam (cap off)        28        0            2     tau_max
  E3 rcond=0.01  Adam (cap on)         14        0            0     tau_max
  E3 rcond=0.01  Adam (cap off)        14        0            0     tau_max

  total forced-through-unaccepted steps across all Adam runs: 0
```

Zero, everywhere, including E1's printed `forced through unaccepted=0`. So this is a correctness
defect with **no effect on any published number**, and my pre-run prediction that it would matter
was wrong — recorded here rather than quietly dropped.

It still matters for two reasons. It makes the code disagree with its own stated assumption, which
is the kind of drift an audit exists to catch before it does fire. And it means README:117's "the
controller is a strict generalization of the ratchet" is **false in general** — the two laws also
diverge on the parent's residual branch, which shrinks on `rel_residual > RESIDUAL_REL_TOL` where
Adam's residual channel can only withhold growth (§7.3). The degeneracy check in
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility|G6_Attribution_And_Reproducibility]] §6.2 confirms the generalization holds *on the trajectories
tested*, which is a weaker statement than the one made.

Related and smaller: the retry loop at `:276` hardcodes the module-level `SHRINK`, where
`adaptive_qite_run` takes `shrink` as a keyword argument. A caller passing `shrink=` to
`adam_qite_run` gets it silently ignored. No caller does.

## 7.2 `s_e` breaks its stated bound, and can overflow

Assumption 4 and README:27 state $s_e = \hat m_e/(\sqrt{\hat v_e}+\varepsilon) \in (0,1]$ "by
Jensen". The saved output reports:

```
  Adam state over run: e_hat median=2.902e-05, s_e range=[5.097e-05, 1.015], growth range=[1, 1.25]
```

$1.015 > 1$. The cause is in [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G2_Beta2_Horizon_Dominates_Every_Verdict|G2_Beta2_Horizon_Dominates_Every_Verdict]] — Jensen requires both
moments under one weighting, and $\beta_1 \ne \beta_2$ gives two. This entry records the
consequences rather than the cause.

The violation is 1.5 % here and harmless. Under a reversed mismatch ($\beta_1=0.999$,
$\beta_2=0.9$) V9 measured a **mean $s_e$ of 258**, and the run emitted:

```
  RuntimeWarning: overflow encountered in scalar power
```

from `adam_qite.py:209`, `growth = max(m_e_hat, self.eps) ** (-s_e * self.k_i)`. With
$\hat m_e \ll 1$ and $s_e$ large, the exponent $-s_e k_I$ is large and positive against a small
base, or large and negative against a base near 1 — either way the power overflows before the clip
at `:215` ever sees it. The clip catches `inf` correctly, so the run completes with the right
answer; but the code is relying on `np.clip(inf, SHRINK, GROW)` rather than on the bound it claims.

**The plots hide it.** `adam_qite.py:485` and `:571` both set `ax.set_ylim(0, 1.05)` on the $s_e$
panel. The one shipped visualization of the quantity is clipped just above the bound the quantity
violates, so a reader of `target_e1_controllers.png` cannot see the excursion that the text output
prints two lines up. This is not deliberate — the limit is the natural choice *if* the bound holds
— but it is why the violation survived into the README.

Fix: assert the bound and fail loudly, per [[QITE_Hand_Off/Audits/Audit_Project_Wide/M2_Noise_Model_Falsifiability|M2_Noise_Model_Falsifiability]]'s pattern; clamp or
log-guard the exponent at `:209`; and set the $s_e$ panel's `ylim` from the data.

## 7.3 The residual channel has never fired

Assumption 6 gives the McLachlan residual its own EMA pair and its own growth ceiling
(`:211-213`). `Output_Pre_Audit.txt`:

```
  McLachlan residual channel withheld growth on 0 of 190 steps (m_r_hat never exceeded 0.05)
```

Zero of 190 in E1, and V11 confirms zero at every band setting down to $10^{-6}$. Half the
controller's state (`m_r`, `v_r`, `s_r`) is maintained on every step, reported in the diagnostic,
and has never influenced a single step of a single experiment. It is a design feature with zero
supporting evidence either way — not shown to help, not shown to be unnecessary.

That is not automatically a defect: `RESIDUAL_REL_TOL = 0.05` is inherited from the parent, where
the same channel *does* fire, and a guard that never trips on well-conditioned problems is doing
its job. But the module's own fairness protocol (Assumption 9) is about not claiming what was not
measured, and the README presents the residual ceiling as part of the design without noting it is
untested. Either exercise it — a stiffer problem, or a run at a tolerance where it engages — or
say plainly that it is inactive in every experiment reported.

## 7.4 Summary of what did not fire

| Defect | Live in source | Fired in E1/E2/E3 | Evidence |
|---|---|---|---|
| growth applied to unaccepted steps | `:310` | no (0 occurrences) | V4 |
| `ctrl.update` on unaccepted steps | `:283` | no (0 occurrences) | V4 |
| `shrink=` kwarg ignored | `:276` | no (no caller passes it) | read |
| $s_e > 1$ | `:204` | **yes**, to 1.015 | saved output |
| exponent overflow | `:209` | no (only under reversed mismatch) | V9 |
| residual channel inactive | `:211-213` | **yes** — 0/190, i.e. it never acts | saved output, V11 |

**Falsifiable prediction.** A problem stiff enough to exhaust `max_retries` — deeper rank
deficiency than E3's rank-6-of-9, or `delta_tau_min` raised — will produce forced-through steps,
and 7.1 will then be visible as growth applied immediately after a rejected step. This audit did
not construct such a case; it is the cheapest way to turn 7.1 from a reading of the source into a
measurement.
