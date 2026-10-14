---
title: "Audit — Adaptative_QITE project (project-wide audit)"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
scope: QITE_Hand_Off
---

# Audit — Adaptative_QITE (project-wide)

Full read of every note, script and figure in `QITE_Hand_Off/` (23 files, 2611
lines of prose + code), plus the source memo `Research_Papers/QITE_timestep.pdf` and the two
replication scripts. This note is the hub: **where the project stands**, a **findings register**,
and an index of the **modules** that close each finding. Routes to build next live in
[[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]].

Scope is **physics and code**. Vault-hygiene defects are noted at the bottom but not scheduled.

> [!warning] Headline
> The project already derived the correct hardware depolarizing bias factor — in its own QDP
> replication, recorded as "Theorem 2 is a factor of two pessimistic" — and then did not use it
> in `adaptive_qite.py`. That single defect ([[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|F1]]) is
> sufficient to explain Target B's headline "surprise", and its predicted magnitudes reproduce
> the numbers the Route-3 note reports to two significant figures.

## 1. Stage

| Strand | State |
|---|---|
| Literature | 4 PDFs, 3 `Paper_Reading` notes + 1 implementation note standing in for the un-citable internal memo |
| Replication | **Complete** — [[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP\|Replication_Practical_Limitations_QDP]], 5 targets, honest per-target verdicts |
| Method prototype | **Built and run, not yet trustworthy** — [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation\|Adaptive_QITE_Attack_Route_3_Implementation]], Route 3 of 5, 688 lines, 4 targets, 4 figures |
| Routes 1, 2, 4, 5 | Not started |
| Shot noise | Not modelled (Assumption 7) |
| Systems beyond 2-qubit H$_2$ | None |
| TITE lever | Note written, no follow-through |

The replication is the strongest artifact in the project. It produced a genuinely new analytical
result (the factor of two) and one unresolved *external* question — Appendix C's $\mathrm{cond}(M)$,
where the note itself concludes the scan space is exhausted and the question belongs to Gaurav.
That question was still open when this audit was written.

**Net: end of a first prototype cycle.** Two things characterise the position:

1. The project has **no local truncation error estimate at all**. Route 3 monitors physics
   indicators — residual, energy monotonicity, $\mathrm{cond}(M)$ — but never estimates the Euler
   step's own error. So the memo's central $\Delta\tau^\star$ / noise-floor argument (Eqs. 8–11)
   cannot currently be evaluated, and the controller has no tolerance parameter to speak of.
2. Of the four Route-3 headline results, **A, B and D each rest on at least one defect below.
   C is the one that survives roughly intact** — and C is also the one that ties back to the
   project's own load-bearing rank-deficiency finding.

So the blocking question is not "what to build next" but "are the current numbers reportable".

## 2. Findings register

Line references are to `Attack_Route_03/adaptive_qite.py` unless stated. Throughout, write
$q := (1-p)^N$.

| # | Sev | Finding | Module |
|---|---|---|---|
| F1 | critical | Debias factor contradicts the project's own derived result | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction\|M1_Depolarizing_Debias_Correction]] |
| F2 | critical | "Flow-direction error" is a misnomer — the noisy flow is exactly collinear | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M2_Noise_Model_Falsifiability\|M2_Noise_Model_Falsifiability]] |
| F3 | critical | A reported per-step number has no code path | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M2_Noise_Model_Falsifiability\|M2_Noise_Model_Falsifiability]] |
| F4 | high | RMS trajectory metric silently extrapolates stalled runs | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M3_Trajectory_Comparison_Metric\|M3_Trajectory_Comparison_Metric]] |
| F5 | high | Energy gate uses the exact-ITE model, not the variational one; fights the memo near convergence | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence\|M4_Controller_Near_Convergence]] |
| F6 | high | "Zero violations of Eq. (13)" is narrower than stated | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence\|M4_Controller_Near_Convergence]] |
| F7 | medium | Targets are not independently reproducible (shared RNG) | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M5_Reproducibility_Hardening\|M5_Reproducibility_Hardening]] |
| F8 | medium | Rejected steps are committed to the trajectory anyway | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M5_Reproducibility_Hardening\|M5_Reproducibility_Hardening]] |
| F9 | medium | The matched-budget comparison is not matched | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting\|M6_Honest_Cost_Accounting]] |
| F10 | medium | $\mathrm{cond}(M)$ is reported *through* the regularizer, so `rcond` sets its ceiling | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence\|M4_Controller_Near_Convergence]] |
| F11 | low | Physics layer duplicated with no agreement test | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M5_Reproducibility_Hardening\|M5_Reproducibility_Hardening]] |
| F12 | low | `vector_Y` drops the $-2\langle H\rangle\rho$ term without saying so | [[QITE_Hand_Off/Audits/Audit_Project_Wide/M5_Reproducibility_Hardening\|M5_Reproducibility_Hardening]] |

### Which reported claims are affected

| Claim | Source | Status after audit |
|---|---|---|
| A: "adaptive wins at equal cost, −1.1197 vs −1.1167 at 61 evals" | note + README | **Unsupported as stated** — F9 (retries uncounted), F5/F10 (cause of the $\Delta\tau$ collapse ambiguous) |
| A: "zero violations of Eq. (13)" | note | **Narrower than stated** — F6 |
| B: "debiasing does not help; biased tracks better in RMS" | note + README | **Explanation wrong, and the metric is broken** — F1, F2, F4 |
| B: "debiasing shrinks per-step flow-direction error ~5% → ~67%" | note + README | Number is **right**, wording is **wrong**, and it is **not reproducible from the script** — F1, F2, F3 |
| B: "reproduced identically across Target B and Target D" | note | **False as written** — different $\theta_0$ (F7) |
| C: "tight `rcond` beats loose `rcond` despite a higher reported cond(M)" | note + README | **Survives**, with the caveat that $\mathrm{cond}_F$ is monotone in the cutoff by construction (F10), which the replication note already states |
| D: "noise floor around $p \approx 10^{-3}$" | note + README | **Unsupported** — built on the same broken RMS (F4) and the same wrong debias (F1) |

## 3. Modules

Ordered by the sequence they should be executed in, not by finding number.

| Module | Closes | Effort | Why here in the order |
|---|---|---|---|
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction\|M1_Depolarizing_Debias_Correction]] | F1 | S | One expression; invalidates B and D until re-run |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence\|M4_Controller_Near_Convergence]] | F5, F6, F10 | M | Changes the controller itself, so everything downstream must re-run anyway |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M3_Trajectory_Comparison_Metric\|M3_Trajectory_Comparison_Metric]] | F4 | S | Needed before B/D can be re-read |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M2_Noise_Model_Falsifiability\|M2_Noise_Model_Falsifiability]] | F2, F3 | S | Turns F1/F2 into an assertion the script checks |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M5_Reproducibility_Hardening\|M5_Reproducibility_Hardening]] | F7, F8, F11, F12 | M | Housekeeping; do before adding routes, not after |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting\|M6_Honest_Cost_Accounting]] | F9 | M | Supplies the currency Routes 1/2 will be compared in |

Then build: [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]].

**Sequencing rationale.** Route 2 will be benchmarked *against* Route 3, so Route 3's numbers have
to be trustworthy first. M1 and M4 both force a full re-run of Targets A–D; doing them before M3
and M2 avoids re-running twice.

## 4. Noted, not scheduled

Out of scope by decision, recorded so they are not lost:

- [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §4–§7 contain bare LaTeX outside `$...$` delimiters
  (`\frac`, `\text{cond}(M)`, `\epsilon`, …), so that math does not render in Obsidian. §1–§3 are
  fine; the break starts at the Theorem 1 synthesis.
- All three `Paper_Reading` notes are still `status: reading`, including the QDP paper that has
  been fully replicated and mined for a new result.
- `FileClasses/` defines only `Application` and `Paper_Reading`. The Replication and Implementation
  notes — two note types in active use — have no schema and no template, so `frontmatter_lint.py`
  skips them entirely.
- The technical resume §5 says "H$_2$ ($N=6$)" where every other document uses $N=5$ for the
  simple ansatz.

## 5. What this audit did not check (at time of writing)

- The four PNGs in `Attack_Route_03/figures/` were not regenerated or compared against the numbers
  in the note; the audit read the plotting code, not the images. (Regenerated during
  implementation, §6 below.)
- `replicate.py`'s five targets were not re-run. Its verdicts are taken as given — the audit only
  compared its shared physics functions against Attack_Route_03's copies (F11).
- The memo's Routes 1, 2, 4, 5 were read for [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] but no cost model was validated.

## 6. Implementation status

All six modules were implemented and re-run against `Attack_Route_03/adaptive_qite.py` the same
day the audit was written, in the planned sequence (M1 → M4 → M3 → M2 → M5 → M6). Every
falsifiable prediction in a module note was checked against a real run, not assumed. Detail and
before/after numbers are in each module note; summary:

| Module | Outcome |
|---|---|
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction\|M1_Depolarizing_Debias_Correction]] | **Confirmed exactly.** Corrected code drives the debiased flow error to the finite-difference noise floor (<0.0001%), from the predicted-and-observed 5.3%/34.7% (biased). |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M2_Noise_Model_Falsifiability\|M2_Noise_Model_Falsifiability]] | **Confirmed.** Collinearity residual $\sim10^{-8}$, $c_{meas}$ matches $c(p)$ to 6 digits. Reparametrization test also confirmed the mechanism: the biased run tracks a $\tau/c(p)$-rescaled baseline measurably better than it tracks the plain one. |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M3_Trajectory_Comparison_Metric\|M3_Trajectory_Comparison_Metric]] | **Confirmed** on a synthetic regression (truncated-but-identical run scores ~0, not the old metric's spurious nonzero). No run stalls once M1 is applied — coverage is 100% everywhere post-fix, so the "$p\approx10^{-3}$ noise floor" reading in Target D is retracted, not just caveated. |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence\|M4_Controller_Near_Convergence]] | **Partially confirmed, with a genuine (not predicted) result.** The variational gate does fire far less often than the exact-ITE gate did. But Target A's $\Delta\tau$ collapse turned out to still be the $\mathrm{cond}(M)$ cap (31/43 tail steps), not the old gate misfiring (0/43) — the stiffness reading is *vindicated*, on firmer evidence, not overturned. Target C changed the most: post-fix, both `rcond` choices converge (no more 500-step stall); the core finding survives at a smaller margin ($\sim2.4\times$, not $\sim20$–$100\times$). |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M5_Reproducibility_Hardening\|M5_Reproducibility_Hardening]] | **Confirmed and fixed.** Per-target seeding verified stable across `--target x` vs. `--target all`. `check_eq_b9_self()` passes (worst diff $8.9\times10^{-11}$) against this script's own `metric_M`/`drho`, closing the F11 gap. Forced-accept steps are now distinguished from real retries in the log and the printed summary. |
| [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting\|M6_Honest_Cost_Accounting]] | **Confirmed the headline claim survives, on the corrected currency.** Adaptive still beats fixed at matched circuit cost (not just f-evals): $-1.1210$ vs. $-1.1167$ at 732 circuits. Relative to the ansatz's own variational floor ($E_{floor}=-1.1246$, computed and plotted), adaptive lands $5\times$ closer than fixed ($1.6\times10^{-3}$ vs. $7.9\times10^{-3}$). |

`Attack_Route_03/README.md`, [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] and
[[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §4 and §6 were updated alongside it to carry the
corrected numbers; the pre-audit claims are marked (strikethrough or a dated callout) rather than
silently replaced, per the vault's "read a note before overwriting it" rule. All four figures in
`Attack_Route_03/figures/` were regenerated from the corrected code.
