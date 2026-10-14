---
title: "Report — the Adam-like step controller for variational QITE"
tags:
  - report
  - physics
  - quantum-computation
  - adaptative-qite
code: QITE_Hand_Off/Adam_QITE/adam_qite.py
audit: "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]"
---

# Report — the Adam-like step controller for variational QITE

This report is the entry point to `Adam_QITE/`. It states what the module does, what it measures,
how the timestep is actually updated, and what survives of its original claims after the
Adam audit. It is written for a reader who has not seen the code.

Two things should be said at the top. First, **the module is an original extension, not a
reproduction**: the source memo proposes no controller of this kind, so there is no published
number to score against and no verdict table in the usual replication sense. Second, **the audit
overturned part of the module's own conclusions.** Where `README.md` and this report disagree, this
report is current; the README's verdict cells have been corrected to match, and the reasons are
recorded in the linked module notes.

## 1. What the work is

Variational QITE evolves ansatz parameters along McLachlan's variational principle,
$M\dot\theta = Y$, integrated as explicit Euler:

$$\theta_{n+1} = \theta_n + \Delta\tau\, M^{+}Y$$

The parent implementation, Route 3 (`../Attack_Route_03/adaptive_qite.py`), chooses $\Delta\tau$
with a hand-tuned multiplicative ratchet — $\times 1.25$ on a clean accepted step, $\times 0.5$ on
a rejected one. Assumption 1 of that script flags both constants as tuned by eye, and
[[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R3]] names replacing them as a priority.

This module replaces the ratchet with a control law built on an observation about what Route 3
already computes. Its `energy_check` forms $\Delta E_{\text{actual}} - \Delta E_{\text{pred}}$
where $\Delta E_{\text{pred}} = -(Y\cdot\dot\theta)\Delta\tau$. Because
$Y_k = -\partial_k E$ holds exactly for this pure-state parametrization, $\Delta E_{\text{pred}}$
is the *exact* first-order Taylor term of the step, so the residual is $O(\Delta\tau^2)$ — a local
truncation error proxy obtained at **zero extra circuits**. Route 3 uses it only as a binary
in-band/out-of-band retry test and otherwise discards it. [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]'s R1 pays three
$f$-evaluations per step (step doubling) for the same information.

The error signal is normalized by the band the retry gate already uses, so no new tolerance is
introduced:

$$\hat e_n = \frac{|\Delta E_{\text{actual}} - \Delta E_{\text{pred}}|}
{\max(\texttt{DE\_BAND\_ABS},\ \texttt{DE\_BAND\_REL}\,|\Delta E_{\text{pred}}|)}$$

so $\hat e_n > 1$ *is* the existing reject condition. Adam's two exponential moving averages are
then applied to $\hat e$, and the confidence ratio $s_e = \hat m_e/(\sqrt{\hat v_e}+\varepsilon)$
scales the exponent of a Gustafsson integral step-size law.

**Scope, and why this is still QITE.** The moments drive the **scalar** $\Delta\tau$ only. The
parameter update remains $\theta \leftarrow \theta + \Delta\tau\,M^{+}Y$, so the trajectory is a
genuine McLachlan/QITE trajectory and $\tau$ keeps its physical meaning. Per-parameter Adam on the
flow vector — a diagonal preconditioner overwriting $M^{+}$ — is deliberately **not** implemented,
because it would stop being imaginary-time evolution. This is the single most important design
boundary in the module and it is respected throughout.

## 2. What it measures

Three targets, each reusing the seed of the Route-3 target it is compared against, so $\theta_0$ is
identical across arms. Every shared threshold (`GROW`, `SHRINK`, `DELTA_TAU*`, `DE_BAND_*`,
`RESIDUAL_REL_TOL`, `RCOND_DEFAULT`, `STIFFNESS_CAP_COEFF`, `MAX_STEPS`, `MAX_RETRIES`, `TAU_MAX`)
is imported from the parent and never redefined, so the two controllers differ **only** in how
$\Delta\tau$ evolves between accepted steps.

| Target | Scenario | Seed |
|---|---|---|
| E1 | Adam vs Route-3 ratchet vs fixed step, noiseless, 5-parameter | `SEED_A` |
| E2 | Adam vs Route 3 under global depolarizing noise, $p\in\{10^{-4},10^{-3},10^{-2},5\times10^{-2}\}$, Eq. (7) debiased | `SEED_D` |
| E3 | Rank-6-of-9 (9-parameter) ansatz, `rcond` $\in\{10^{-8},10^{-2}\}$, $\mathrm{cond}(M)$ cap on and off | `SEED_C` |

The system is 2-qubit H$_2$ throughout. All three targets report in the **matched-circuit-budget
currency** of [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting|M6_Honest_Cost_Accounting]]: $2N+1$ circuits per `measure_f` plus one per retry
trial, so a controller that takes cheaper steps is not credited for taking more of them. This
currency is what makes the E1 comparison meaningful, and it is what dissolves one of the module's
apparent wins (§4).

E2 additionally scores each controller **against its own noiseless baseline**, not against a shared
reference, because the two controllers visit different $\tau$ grids. That choice is sound in
principle and is precisely where the experiment breaks in practice (§4).

## 3. How $\Delta\tau$ actually gets updated

This is the part most likely to be misread, so it is given in full. Per accepted step:

1. $\hat e_n$ is formed from the energy-check residual, normalized by the gate band.
2. The two EMAs update, **once per accepted step** (a retry cascade does not feed rejected trials
   into the moments), with bias correction counting accepted steps:
   $m \leftarrow \beta_1 m + (1-\beta_1)\hat e$, $v \leftarrow \beta_2 v + (1-\beta_2)\hat e^2$,
   $\hat m = m/(1-\beta_1^n)$, $\hat v = v/(1-\beta_2^n)$.
3. The confidence ratio $s_e = \hat m_e/(\sqrt{\hat v_e}+\varepsilon)$ is computed.
4. The growth factor is $\texttt{growth} = \max(\hat m_e,\varepsilon)^{-s_e k_I}$ with
   $k_I = 0.15$ — Gustafsson's $0.3/(p+1)$ at $p=1$ for explicit Euler, the value R3 already
   proposes.
5. `growth` is clipped to $[\texttt{SHRINK}, \texttt{GROW}] = [0.5, 1.25]$. `GROW`/`SHRINK` survive
   only as a trust region, no longer as the step law.
6. The McLachlan residual channel has its own EMA pair and enters as a growth **ceiling** only,
   never as a shrink, preserving parent Assumption 2 ($r_n$ is $\Delta\tau$-independent, so it
   cannot justify shrinking the current step).
7. $\Delta\tau \leftarrow \Delta\tau \cdot \texttt{growth}$, then the $\mathrm{cond}(M)$ stiffness
   cap $\texttt{STIFFNESS\_CAP\_COEFF}/\mathrm{cond}(M)$ is applied, and the accept/reject retry
   gate runs.

**The load-bearing point: shrinking never comes from the Adam law.** Because $\hat e > 1$ is the
reject condition, an *accepted* step always has $\hat e \le 1$; the EMA of quantities $\le 1$ is
itself $\le 1$; and $\hat m_e \le 1$ raised to a negative exponent is always $\ge 1$. So the
control law can only ever grow the step. Measured on E1: $\hat e > 1$ on **0 of 190** steps
(max $\hat e = 0.945$), and `growth < 1` on **0 of 190** steps. Every reduction of $\Delta\tau$ in
this module comes from one of two inherited mechanisms — the retry loop's `SHRINK`, or the
$\mathrm{cond}(M)$ cap — neither of which is new work.

This is not a coding error; it follows from reusing a rejection threshold as a control setpoint.
Reusing it was the right call for *not inventing a tolerance*, but a good rejection band and a good
control setpoint are different objects. See [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G4_Controller_Cannot_Shrink|G4]].

## 4. Results as they now stand

The audit ran twelve experiments (V1–V12) and every number below is regenerable — see §8.

| Claim | Original verdict | Status after the audit |
|---|---|---|
| E1: Adam vs the ratchet, noiseless | "A wash." | **Stands, and is better supported than the README argued.** At matched circuit budget Route 3 wins against all four Adam configurations swept, at both seeds shown. See [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G2_Beta2_Horizon_Dominates_Every_Verdict\|G2]]. |
| E2: the EMA should average noisy $\Delta E$ | "The premise fails." | **Withdrawn — the premise is untested, not failed.** The experiment has no stochastic noise for an EMA to average (parent Assumption 7); the RMS ratio is $p$-independent to $\pm0.016$; at $p\le10^{-3}$ the signal is smaller than the gap between the two baselines being differenced. See [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G3_Noise_Experiment_Cannot_Resolve_Noise\|G3]]. |
| E2: "Adam's own $p=0$ baseline stalls" | 500 steps to $\tau=2.86$ vs Route 3's 234 to $\tau=3.01$ | **Withdrawn as an artifact.** The stall is caused by *mismatched EMA horizons*, not by the Adam controller. See below. |
| E3: cap-off win | "The one real win, and it is attributable." | **Re-attributed to the cap, and weakened.** Four of the eight rows are duplicate trajectories, there is no cost-matched column and no matched-mean-step control arm, and $n=1$. See [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility\|G6]]. |
| Assumption 4: $s_e \in (0,1]$ by Jensen | asserted | **False as implemented.** The run itself reports $s_e$ up to $1.015$; a reversed mismatch gives mean $s_e = 258$. Jensen's inequality needs one weighting, and the code uses two. |
| The $O(\Delta\tau^2)$ free-error-estimate premise | asserted, never tested | **Confirmed, quantitatively.** See §5. |

### The EMA horizon mismatch

This is the audit's headline finding. The two decay rates are the stock Adam defaults,
$\beta_1 = 0.9$ and $\beta_2 = 0.999$ — memories of roughly 10 and 1000 steps respectively, on runs
of 172–500 steps. So $s_e$ divides a ten-step local mean by an RMS taken over a window longer than
the entire experiment, applied to a signal that falls four orders of magnitude as the run
converges. It is not a signal-to-noise ratio; it is a measure of how far the run has converged
since it started.

Sweeping the pair on `SEED_D`, E2's stalling baseline:

| configuration | steps | $\tau_{\text{end}}$ | final $E$ | mean $s_e$ | frozen | exit |
|---|---|---|---|---|---|---|
| Route 3 (ratchet) | 234 | 3.0097 | $-1.123745$ | n/a | n/a | `tau_max` |
| $\beta_1{=}0.9$, $\beta_2{=}0.999$ (shipped) | 500 | 2.8615 | $-1.123002$ | 0.0893 | 72 % | `max_steps` |
| $\beta_1{=}\beta_2{=}0.999$ — matched | 254 | 3.0269 | $-1.123757$ | 0.5170 | 0 % | `tau_max` |
| $\beta_1{=}\beta_2{=}0.99$ — matched | 259 | 3.0226 | $-1.123736$ | 0.4387 | 0 % | `tau_max` |
| $\beta_1{=}\beta_2{=}0.9$ — matched | 244 | 3.0093 | $-1.123736$ | 0.3938 | 0 % | `tau_max` |
| $\beta_1{=}0.999$, $\beta_2{=}0.9$ (reversed) | 242 | 3.0077 | $-1.123740$ | 258.1220 | 0 % | `tau_max` |

**Every matched pair completes at $\tau_{\max}$; only the mismatched pair stalls.** The value does
not matter — including $\beta_2 = 0.999$ itself, the very value the shipped configuration uses and
the one originally blamed. Growth freezes on 0 % of steps in every matched run against 72 % in the
shipped one. This also falsifies the heavy tail of $\hat e$ as a *sufficient* cause: the signal is
equally heavy-tailed in all six rows, and with matched horizons it causes no freezing at all.

Note that the 500-step row exits on `MAX_STEPS`, so its $\tau = 2.8615$ is a function of an upstream
constant rather than of convergence.

## 5. The one confirmed positive

The module's founding premise — that `energy_check` discards a usable local truncation error
estimate — is **correct, and is the result worth carrying forward.** It was asserted in an
`ASSUMPTIONS` string and demonstrated nowhere.

Re-running the E1 trajectory and additionally performing R1's step doubling at each accepted step
(one full step of $\Delta\tau$ against two half-steps of $\Delta\tau/2$ with a fresh $f$ at the
midpoint), over the first 40 accepted steps:

- log-domain Pearson correlation between the free estimator and the step-doubling error: **0.9902**
- Spearman correlation: **0.9882**
- ratio $|\Delta E_{\text{actual}} - \Delta E_{\text{pred}}| / e_{\text{true}}$: median **2.060**,
  range $[0.350, 11.439]$

The median of 2 is the *predicted* value, not a fit: step doubling has full-step error
$\approx C\Delta\tau^2$ against two-half-step error $\approx C\Delta\tau^2/2$, so their difference
is $C\Delta\tau^2/2$ — half of what the free residual measures. **A factor of 2, recovered to
within 3 %,** and it tightens as more of the converged tail is included (the full 190-step
trajectory gives Pearson 0.9983 and median 2.005).

So the discarded quantity is the same local truncation error R1 would pay three $f$-evaluations per
step to obtain, up to a known constant factor, at zero extra circuits. This partially closes the gap
[[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]] §1 identifies as the project's first structural problem ("the
project has no local truncation error estimate at all"), and it should change R1's priority in
[[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]].

**Where it is blind, and one refuted attempt to show it.** The estimator captures the curvature of
$E$ along a straight line in $\theta$, but the Euler step's true local error also contains a term in
$\partial_\theta M^{+}$ — the metric changing along the step — which the estimator cannot see. That
is a derivation, and it remains standing. The audit additionally predicted that the ratio's
departure from 2 would track the *step-to-step change* in $\mathrm{cond}(M)$ rather than its level,
and **that prediction failed**: the level correlates far more strongly than the change
($-0.8938$ vs $+0.3176$), and the conclusion survives controlling for the obvious confound that
$\mathrm{cond}(M)$ rises monotonically along the run (rank partial correlations $-0.5437$ vs
$+0.1523$). The blind spot is still argued from the algebra; it has not been *measured*, and the
probe proposed for measuring it needs a different design. See
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot|G5]].

## 6. Flaws

The register from [[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]], one line each. G1 is fixed; the rest are recorded,
not fixed.

| | Severity | Finding |
|---|---|---|
| G1 | — | The sibling import broke when the folder moved; **fixed**, and the fresh run reproduces every prior number, confirming the fix changed nothing but the import. [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G1_Broken_Import_And_Stale_Paths\|→]] |
| G2 | critical | The two EMAs use mismatched horizons (10-step vs 1000-step on 172–500-step runs), which is what E2's stall measured and what breaks the claimed $s_e \in (0,1]$ bound. [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G2_Beta2_Horizon_Dominates_Every_Verdict\|→]] |
| G3 | critical | E2 cannot resolve noise: there is no stochastic channel for an EMA to average, and the reported ratio is measuring the horizon artifact instead. [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G3_Noise_Experiment_Cannot_Resolve_Noise\|→]] |
| G4 | high | The inherited reject band is far too loose to serve as a control setpoint, which is why the law never shrinks. [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G4_Controller_Cannot_Shrink\|→]] |
| G5 | medium | The founding error-estimator premise is *correct* but was never tested by the code; its stiffness blind spot is derived, not measured. [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot\|→]] |
| G6 | high | E3's win is not cleanly attributable (duplicate rows, no cost-matched column, $n=1$); `energy_at_budget` has an off-by-one that credits each arm with an unpaid step; headline numbers were produced in a scratchpad that no longer existed. [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility\|→]] |
| G7 | low | Latent defects: the retry loop hardcodes the module-level `SHRINK` rather than the instance's; the $s_e$ figure axis clips at 1.05 and hides the $1.015$ excursion. [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G7_Latent_Defects\|→]] |

The `energy_at_budget` off-by-one is worth one extra sentence, because it touches a reported result:
corrected, E1's matched-budget gap changes from $+3.5\times10^{-5}$ to $+3.0\times10^{-5}$. The sign
is preserved and the verdict is unchanged — the adaptive-vs-fixed gap is 925 times the correction —
so this is a defect in the instrument, not a retraction of a finding.

## 7. What was not checked

- **Shot noise.** There is none anywhere in the module (R5 in [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]). The E2
  "noise" is a deterministic depolarizing channel, which is why G3 holds.
- **$\beta_1$ alone**, and momentum-only vs RMS-only ablations. The horizon sweep varies the pair;
  it does not isolate the two moments.
- **Systems beyond 2-qubit H$_2$.** Every number in this report is one molecule.
- **Seeds beyond three.** `SEED_A`, `SEED_C`, `SEED_D` only; the V10 matched-budget table covers
  two of them.
- **E3's $n=1$.** The cap-off win is a single configuration at a single seed on a single ansatz.
- **A robust dispersion estimate** (MAD, or a log-domain second moment) in place of
  $\sqrt{\hat v}$ — the change the diagnostic predicts is most likely to matter, since the heavy
  tail of $\hat e$ is real (mean/median 948–2317) even though it is not the cause of the stall.

## 8. Open work, and how to reproduce this

Priority order, from the module notes:

1. **Tie the two EMA horizons together** — set $\beta_1 = \beta_2$, or derive both from the
   expected step count. Deriving a decay from the problem is defensible in `ASSUMPTIONS`; it is not
   tuning, and the choice should be made before looking at the resulting energies.
2. **Promote the step-doubling comparison into a test** in the style of the parent's
   `check_eq_b9_self`, turning the module's founding assumption into something the code verifies.
3. **Assert the bound the code claims.** If Assumption 4 says $s_e \in (0,1]$, the controller should
   check it and fail loudly; that would have caught G2 on the first run.
4. **Separate the setpoint from the reject band** — the band cannot be retuned after the fact
   without becoming a fudge, so this needs a principled setpoint derived independently.
5. **Re-run E3** cost-matched, with a matched-mean-step control arm, at more than one seed.
6. **Design a probe that can actually detect the $\partial_\theta M^{+}$ blind spot**, holding the
   convergence stage fixed — the current trajectory does not sample it densely enough.

To reproduce everything in this report:

```bash
cd QITE_Hand_Off/Adam_QITE
~/micromamba/envs/qml-julia/bin/python adam_qite.py --target all
~/micromamba/envs/qml-julia/bin/python verify_audit.py --target all
```

The first regenerates the shipped suite (~75 s; stdout in `Output_Post_Audit.txt`, with the
pre-audit run preserved in `Output_Pre_Audit.txt` for comparison). The second re-derives the
audit's headline numbers for G2, G5 and G6 from the shipped modules without modifying either
(~28 s; stdout in `Output_Audit_Verify.txt`). It compares **147 recorded values against
the audit notes and reports 0 drift**, prints its own assumptions and provenance, and exits
non-zero if any number moves. Drift is to be recorded as a result, never absorbed by loosening a
tolerance — the script deliberately has no tolerance flag.

That count is a **lower bound on the audit's numeric content, not a full accounting of it**: a
number the script does not carry is compared against nothing, which is a failure mode a drift
ledger cannot see. Two such numbers have already been found — a missing horizon-sweep row, which
crashed, and G5's probe-count sensitivity, which printed silently and was caught only by an
independent check. The script now prints a *What this ledger does NOT cover* section naming the
rest, chief among them G6 §6.3's `rcond` table, which carries the evidence for re-attributing E3
and is not regenerated here.

### The independent check

The verification script checks the audit against the code. It cannot check itself, so a separate
agent was dispatched to do that, with no write access to the vault. Its findings are recorded here
because two of them changed the script; each was re-verified directly before being acted on, since
an agent's report is model output and carries no more authority than any other claim.

| Check | Verdict |
|---|---|
| Does `RECORDED` faithfully transcribe the notes? | **Partially reproduced.** Two recorded numbers were outside the ledger: G5's probe-count sensitivity (three pairs, printed but compared against nothing) and the stock row of G6 §6.2's degeneracy table. Both are now in it — the count rose 140 → 147, and all seven new comparisons match |
| Do the fresh runs still show the recorded symptoms? | **Reproduced**, all five exactly: $s_e$ range $[5.097\times10^{-5},\,1.015]$, `growth < 1` on 0 of 190 steps, $\hat e > 1$ on 0 of 190, the McLachlan residual channel firing on 0 of 190, and E2's ratios 1.0243 / 1.0246 / 1.0245 / 1.0088 |
| Does the fresh suite output match the pre-audit one? | **Reproduced.** No genuine numeric difference between `Output_Pre_Audit.txt` and `Output_Post_Audit.txt` |
| Does anything in `README.md` or the working note contradict the fresh runs? | No remaining contradiction in the corrected `README.md`. The working note still carries the superseded verdicts and now opens with a banner saying so |

Its most useful finding was not a number. It observed that the script's replay-versus-log
assertion is a **plumbing check, not physics validation** — the replay and the shipped run call the
same `measure_f` and `energy_check`, so an error inside either is invisible to it by construction.
That is correct, and the script's assumptions now say it in those words. The check is still
falsifiable (hidden state between steps, or a logged $\Delta\tau$ that is not the committed one,
would both trip it), but the non-circular check in the folder is the step-doubling comparison of
§5, which had a target it could have missed and did not.

The agent was dispatched twice; the first attempt died on a rate limit and returned nothing, so
nothing from it reached these notes.

**Where to go deeper:** [[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]] for the full findings register and the
experiment list; `README.md` for the operator's guide to the code; and
[[QITE_Hand_Off/Notes_on_Project/Adam_Like_Step_Controller_For_Variational_QITE|Adam_Like_Step_Controller_For_Variational_QITE]] for the first-person working note that preceded
the audit — which still carries the E2 and E3 verdicts corrected in §4 above, and has not yet been
updated. [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] does not yet mention this work at all; the §5
result is the part that belongs there.
