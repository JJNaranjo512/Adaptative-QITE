---
title: "M3 — Trajectory comparison metric"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]]"
---

# M3 — Trajectory comparison metric

> [!success] Implemented and confirmed — see [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE#6. Implementation status|hub §6]]
> `trajectory_deviation()` added (common support, coverage, stall flag via `exit_reason`). The
> synthetic regression test passes: a truncated-but-identical run scores ~0 RMS at partial
> coverage, vs. the old metric's spurious nonzero. No run stalls once [[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|M1]]
> is applied, so Target D's "$p\approx10^{-3}$ floor" reading is retracted outright.

**Findings:** F4 (high)  **Effort:** S  **Blocks:** every quantitative claim in Targets B and D.

## What is wrong

Both Target B (`adaptive_qite.py:509`) and Target D (`adaptive_qite.py:626-627`) compare a noisy
run against the noiseless baseline with

```python
E_interp = np.interp(baseline["tau"], run["tau"], run["E"])
rms = float(np.sqrt(np.mean((E_interp - baseline["E"]) ** 2)))
```

`np.interp` does not extrapolate — outside the sample range it **returns the endpoint value**. So
when a noisy run covers less imaginary time than the baseline, its energy is held flat at the last
value it reached, across the entire uncovered tail, and that flat segment is then scored as if it
were tracking data.

This is not hypothetical for this code. [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] records that
at $p=0.1$ the debiased run "stalls at $\tau\approx1.2$" while the baseline runs to $\tau_{max}=3$.
That means **roughly 60 % of the reported RMS at $p=0.1$ is computed against a constant**, not
against a trajectory.

## Why it matters

The headline of Target B is that the biased run "tracks the noiseless baseline slightly *better* in
RMS energy deviation (0.5749 vs 1.308 at $p=0.1$)". Under this metric, a run that stalls early is
penalised for the flatness of its own extrapolation rather than for its dynamics, and a run that
travels further is rewarded for covering the grid — regardless of how well either tracks. The
comparison is therefore not measuring what it is captioned as measuring.

The same defect propagates into Target D, whose entire output is a sweep of this RMS against $p$,
and hence into the "noise floor around $p\approx10^{-3}$" reading that the note and
[[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] both carry.

Note the interaction with [[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|M1_Depolarizing_Debias_Correction]]: the stall itself is a downstream
consequence of the wrong debias factor changing the residual scale. Fixing F1 may remove the stalls
entirely — but the metric should be fixed regardless, because silently scoring extrapolated data is
a defect even when no run stalls.

## Change to make

Replace the single RMS number with an explicit, three-part comparison:

1. **Common support only.** Restrict to $\tau \in [0, \min(\tau^{\text{run}}_{\text{end}},
   \tau^{\text{base}}_{\text{end}})]$ and compute RMS there. Never interpolate outside a run's own
   range.
2. **Report coverage.** Print $\tau^{\text{run}}_{\text{end}} / \tau^{\text{base}}_{\text{end}}$
   alongside every RMS, so a partial run is visible as partial rather than folded into one scalar.
3. **Treat a stall as its own outcome.** A run that exits on `MAX_STEPS` rather than on
   `tau >= tau_max` did not fail to track — it failed to *finish*, which is a different statement
   about the controller. Add the exit reason to the log and print it. Right now
   `adaptive_qite_run` exits the loop on either condition with no record of which.

A run that stalls should be reported as "stalled at $\tau=1.2$ (40 % coverage), RMS over common
support = X", not as a single RMS that silently mixes the two effects.

## How to verify

- Construct a synthetic case: truncate the baseline itself at $\tau=1.2$ and score it against the
  full baseline. The old metric returns a large nonzero RMS for a trajectory that is *identical*
  wherever it exists; the new one returns ~0 with 40 % coverage. That is the regression test.
- Re-run Targets B and D after [[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|M1_Depolarizing_Debias_Correction]] and check whether any run still
  stalls. If none do, say so — the corrected numbers then stand on their own and this module's
  coverage column is simply 100 % everywhere.

## Documents to update

- [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] — Target B and Target D paragraphs, including the
  0.1225 / 0.1325 / 0.5749 / 1.308 figures and the $p\approx10^{-3}$ floor reading.
- `Attack_Route_03/README.md` — Target B and D rows.
- [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §6, if the noise-floor reading changes.
