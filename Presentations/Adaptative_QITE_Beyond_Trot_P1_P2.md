---
title: "Adaptative QITE — Beyond Trotterization into Route 3: P1 & P2 results"
tags:
  - presentation
  - physics
  - quantum-computation
  - adaptative-qite
---

<!--
Results deck for the two proposals of [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]] that have
been run: P1 (singularity classification) and P2 (tol proportionality). Format: question → prediction → figure → reading →
insight.

The **Slides Extended** plugin is installed in this vault (theme "black", KaTeX math). Open
this note and run "Slides Extended: Start presentation". `---` on its own line = new slide.

Every number here is copied from one of:
  Route_03_and_Beyond_Trot/Proposal_1/P1_Singularity_Classification_Findings.md
  Route_03_and_Beyond_Trot/Proposal_2/P2_Tol_Proportionality_Findings.md
  Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam.md
which in turn trace to Proposal_1/Output.txt and Proposal_2/Output.txt
(51 and 93 pinned values, `--check` exits 0, reruns byte-identical).
Post-hoc results and unrun predictions are labelled as such on the slide.
-->

---
# Beyond Trotterization → Route 3

## P1 & P2: what the cond($M$) cap reads, and whether the controller can hit a `tol`

---

## What this deck covers

Two questions, both opened by the Beyond Trotterization (BT) replication:

1. **P1** — when cond($M$) blows up in Target C / E3 and the cap slams $\Delta\tau$ shut,
   *what* is it actually measuring?
2. **P2** — if the step controller is given an explicit `tol`, does the achieved error
   follow it?

For each: **question** → **prediction (written before the run)** → **figures** →
**conclusions**. Then the common thread and what comes next.

---

## Shared setup

H$_2$, STO-3G, Bravyi–Kitaev → **2 qubits**. McLachlan:

$$M(\theta)\,\dot\theta = Y(\theta), \qquad \dot\theta = M^{+}(\theta)\,Y(\theta)$$

- $M$, $Y$ from finite-difference `drho` at $h=10^{-6}$
- Route 3 safety brake: **cond cap** $\Delta\tau \le 2000/\mathrm{cond}(M)$
- **Target A** — 5-param ansatz, full rank, `rcond` $10^{-6}$
- **Target C / E3** — 9-param ansatz, **rank 6 of 9** by construction ($2d-2$), `rcond`
  $10^{-8}$ or $10^{-2}$

---

## What BT gave us

| BT piece | Used for |
|---|---|
| **Rotated generators** $\tilde G_k=P_kG_kP_k^\dagger$ ⇒ exact $M$, $Y$ | removes fd noise as a suspect (P1) |
| **T6** — boundary-layer vs structural singularity, null-vector attribution | "classify before responding" (P1) |
| **T2′** — achieved error ∝ requested `rtol` | acceptance test for a controller (P2) |

State form of the metric (verified in P1 to $2.594\times10^{-10}$):

$$M_{jk}=2\big(\mathrm{Re}\langle\tilde G_j\tilde G_k\rangle-\langle\tilde G_j\rangle\langle\tilde G_k\rangle\big)$$

---
# Part I
## P1 — what is cond($M$) measuring?

[[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/P1_Singularity_Classification_Findings|P1_Singularity_Classification_Findings]]

---

## P1 — the question

At `rcond` $=10^{-8}$, Target C / E3's cond($M$) peaks at $9.66\times10^{6}$ and the cap pins
$\Delta\tau$. Three suspects:

1. **Gauge leakage (a bug)** — fd error seeping into the 3 structural null directions
   → E3 and Target C contaminated, must be re-run
2. **Real stiffness** — the physics has genuinely fast and slow directions
3. **Chart artefact** — the state moves smoothly, the angles describing it degenerate
   (gimbal lock)

> [!note] Prediction (before the run)
> $\sigma_{7..9}$ sit near $10^{-10}$–$10^{-12}$, below the cutoff → no leakage; rank-6 cond
> matches logged cond; the soft direction is a near-degenerate **gate pair**.

---

## P1 — no leakage

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_target_c_spectrum.png|640]]

- **(a)** All nine $\sigma_i/\sigma_1$. Gauge values (orange) at ~$10^{-16}$–$10^{-18}$,
  **~8 decades under** the cutoff (dashed). Only $\sigma_6$ (blue) dips, at $\tau\approx0.2$
- **(b)** Logged cond vs rank-6 cond: indistinguishable, $|\mathrm{cond}_6/\mathrm{logged}-1|\le1.11\times10^{-15}$

**Insight**: the cap reads conditioning **inside** the 6-dim physical range →
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility|G6]] supported, E3 and Target C **don't need re-running**.

---

## P1 — but it's a chart singularity

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_target_c_spectrum.png|640]]

- **(c)** Parameter speed $\|\dot\theta\|$ spikes to **128.4** at $\tau=0.1983$
- **(d)** State speed $\|d\rho/d\tau\|_F=\sqrt{\dot\theta^\top M\dot\theta}$ is **0.1862** there
  and has no feature at the peak
- Over the run: $\|\dot\theta\|$ varies **×160.8**, $\|d\rho/d\tau\|$ only **×8.83**

**Insight**: stiff in $\theta$, smooth in $\rho$. The *state* isn't doing anything violent —
the *parametrization* is. The capped arms spend **86 of 108** steps crossing this fold.

---

## P1 — which gate is soft?

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_null_vectors.png|720]]

- **$v_6$ (soft)**: essentially **one gate**, th8 = RZ on q0, layer 2 (**0.93**). 99.9998 %
  of $\dot\theta$ points along it at the peak
- **$v_{7..9}$ (gauge)**: the gate pairs live here — th1/th6, th2/th7 — and th4 alone
- Why th8: $\langle Z_0\rangle=0.9944$, so $M_{88}=\tfrac12(1-\langle Z_0\rangle^2)=5.568\times10^{-3}$,
  and its tangent is nearly a combination of the others ($\sigma_6/M_{88}=2.6\times10^{-5}$)

**Insight**: the prediction's "gate pair" was wrong — pairs are the harmless gauge,
already handled by `rcond`. The thing the cap fights is a single RZ near its pole.

---

## P1 — the fold amplifies finite-difference error

Target C rerun on the **exact** $M$, $Y$:

| `rcond`   | arm        | error                   | steps |
| --------- | ---------- | ----------------------- | ----- |
| $10^{-8}$ | fd         | $8.868181\times10^{-4}$ | 108   |
| $10^{-8}$ | exact      | $8.854194\times10^{-4}$ | 108   |
| $10^{-2}$ | fd = exact | $2.143405\times10^{-3}$ | 17    |

fd-vs-exact energy gap at $10^{-8}$: $1.292\times10^{-11}$ before the fold →
$1.101\times10^{-4}$ max after it.

**Insight**: a $10^{-10}$ perturbation of $M$ becomes $10^{-4}$ in the trajectory. Quote
Target C / E3 errors to **2–3 significant figures**. E3's ranking ($5.65$ vs $8.87$ vs
$13.98\times10^{-4}$) is not that close, so it stands.

---

## P1 — two side findings

**`rcond` $=10^{-2}$ hides the singularity.** It keeps only 4–6 values (fewer than 6 on 11 of
17 steps): it throws away *physical* directions too. Logged cond max $125.39$, while the true
spectral $\sigma_1/\sigma_6$ reaches $1.9281\times10^{5}$. "Well-conditioned" because the stiff
direction is projected out.

**Target A silently drops a direction.** The full-rank 5-param ansatz loses CRZ on the
**last 12 of 172** steps ($\tau\in[2.945,\,2.999]$). The state converges to ~$|10\rangle$,
where CRZ only adds a phase → zero-variance row. That's the plan's §2 mechanism appearing
*on its own* as the state converges.

---

## P1 follow-on — a regime-aware cap?

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_regime_aware_cap.png|720]]

Cap reads $2000/\mathrm{cond}_6(M)$ instead of the logged cond. Coefficient unchanged.

| `rcond` | cap | error | steps |
|---|---|---|---|
| $10^{-8}$ | shipped | $8.8682\times10^{-4}$ | 108 |
| $10^{-8}$ | rank-6 | $8.8617\times10^{-4}$ | 108 |
| $10^{-2}$ | shipped | $2.1434\times10^{-3}$ | 17 |
| $10^{-2}$ | rank-6 | $3.4051\times10^{-3}$ | 105 |

**Insight**: identical at $10^{-8}$; at $10^{-2}$ it brakes on stiffness the solve already
discarded (right panel, $\Delta\tau$ dives to $10^{-4}$) and ends **worse**. The cap should read
the conditioning **of the solve actually performed** — which is what it already does.
**Not adopted.**

---

## P1 — scorecard

| Piece of the prediction | Outcome |
|---|---|
| $\sigma_{7..9}$ below the cutoff (no leakage) | ✅ held — max $\sigma_7/\sigma_1=2.021\times10^{-16}$ |
| $\sigma_{7..9}$ at $10^{-10}$–$10^{-12}$ specifically | ❌ failed (favourably) — ~$10^{-16}$: $M$ is a Gram matrix, fd error enters the null space **quadratically** |
| rank-6 cond = logged cond | ✅ held — to $10^{-15}$, all four E3 arms |
| cap reads real stiffness → G6 supported | ✅ held |
| soft direction is a gate **pair** | ❌ failed — one gate, th8 |

---

## P1 — conclusions

> [!success] What P1 decides
> - **No gauge leakage.** G6's re-attribution of E3 to the cap stands
> - **E3 and Target C don't need re-running** on the exact metric
> - Errors at `rcond` $=10^{-8}$: **2–3 significant figures**, no more
> - The stiffness is a **chart singularity** — a fourth kind next to structural /
>   developing / boundary-layer
> - **Don't ship** the regime-aware cap in the plan's form
> - Plan §2's two identities (state form of $M$, sign of $Y$) are **verified numerically**

---
# Part II
## P2 — can the controller hit a requested `tol`?

[[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/P2_Tol_Proportionality_Findings|P2_Tol_Proportionality_Findings]]

---

## P2 — the question

Neither controller had a `tol` to sweep. Route 3 ratchets; Adam uses the reject band as its
setpoint, which [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G4_Controller_Cannot_Shrink|G4]] showed is circular (growth $<1$ on 0 of 190
steps).

**Fix**: memo Eq. (17) PI law at $p=1$, driven by
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot|G5]]'s zero-cost estimator:

$$\Delta\tau_{n+1}=\Delta\tau_n\Big(\frac{\mathrm{tol}}{e_n}\Big)^{0.15}\Big(\frac{e_{n-1}}{e_n}\Big)^{0.2},
\qquad e_n=\tfrac12|\Delta E_{\rm actual}-\Delta E_{\rm pred}|$$

Two tolerance semantics, two theory exponents:

| arm | control | theory $\alpha$ |
|---|---|---|
| **EPS** (per step) | $e_n\le\mathrm{tol}$ | $1/2$ |
| **EPUS** (per unit step) | $e_n/\Delta\tau\le\mathrm{tol}$ | $1$ |

---

## P2 — pre-registered predictions and fit rule

1. (plan, as written) EPS $\alpha\in[0.8,1.2]$ — *expected to fail*: theory says 1/2
2. EPS $\alpha\in[0.4,0.6]$
3. EPUS $\alpha\in[0.8,1.2]$
4. (plan) slope breaks where the cap takes over
5. G4's circularity gone: growth $<1$ in every run
6. EPS circuits $\propto\mathrm{tol}^{-1/2}$

**Metric**: $|E_{\rm run}(\tau_{\rm end})-E_{\rm ref}(\tau_{\rm end})|$ against an exact-metric
DOP853 reference (rtol $10^{-11}$).
**Fit rule**: a run enters only if the cap binds on $<25\%$ of steps (plus 3 other criteria);
need **≥ 3** such runs.

---

## P2 — the headline: no fit is possible

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_tol_proportionality.png|640]]

- **(a)** err($\tau_{\rm end}$) sits near $10^{-5}$ across four decades of `tol`. Local slopes
  wander from $-0.80$ to $+2.62$
- Filled dots = fit-eligible runs: only **2 per arm** (the two smallest tols)

> [!failure] Scored as written
> Predictions 1, 2, 3, 4, 6 **fail** — because **no slope could be fitted**, not because a
> slope landed out of range.

---

## P2 — what is actually setting $\Delta\tau$

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_tol_proportionality.png|640]]

- **(c)** Small tols really do take smaller steps — until $\tau\approx2.1$, where every run
  collapses into the same cap-driven pattern
- **(d)** The cap binds on $\ge25\%$ of steps from the **3rd grid point** in each arm
  (EPS $10^{-5}$, EPUS $3.16\times10^{-3}$) → nothing left to fit

> [!success] Prediction 5 held — G4's circularity is gone
> Growth $<1$ in every run (527 of 1234 steps at EPS $10^{-6}$). Where `tol` binds, median
> $e_n/\mathrm{tol}=0.992$. **The PI law controls $e_n$.** It's the end-point error that
> doesn't follow.

---

## P2 — work-precision

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_tol_proportionality.png|640]]

**(b)**, on the same reference:

| controller | err($\tau_{\rm end}$) | circuits |
|---|---|---|
| Route 3 ratchet (Target A) | $6.3181\times10^{-4}$ | 2069 |
| Adam (E1) | $7.2773\times10^{-4}$ | 2285 |
| PI, EPS $3.16\times10^{-4}$ | $5.67\times10^{-6}$ | 2700 |
| PI, EPUS $10^{-2}$ | $1.18\times10^{-5}$ | 2736 |

> [!warning] Promising, not a claim
> One seed, and err($\tau_{\rm end}$) is **not monotone** in `tol` — part of this may be where
> the contraction happens to leave each run. A second seed first.

---

## P2 — why the end point ignores `tol`

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_first_step_ablation.png|720]]

**(b)** Error along the run, EPS arm:

1. **The first step is tol-blind.** Fixed $\Delta\tau_0=0.05$ injects $5.2914\times10^{-4}$ in
   **every** run ($e_1=229\times$ `tol` at $10^{-6}$) — the jump at $\tau=0.05$
2. **The flow amplifies, then forgets.** Error peaks at $\tau\approx0.94$–$1.15$, then the
   imaginary-time contraction erases 2–3 decades before $\tau=3$

**Insight**: $\tau_{\rm end}$ remembers the late, cap-dominated stretch — not the controller.
(Point 2 is a reading of the figure; no contraction rate was computed.)

---

## P2 — control the first step, and the textbook slopes appear

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_first_step_ablation.png|720]]

> [!warning] Post hoc — added after the fit results, doesn't rescue predictions 1–6
> First step at `DELTA_TAU_MIN` $=10^{-4}$; score = max error over the cap-free $\tau\le2$.

| arm | slope, fixed first step | slope, first step at min | theory |
|---|---|---|---|
| EPS | 0.267 | **0.478** | 1/2 |
| EPUS | 0.515 | **0.890** | 1 |

**Insight (a)**: with the first step controlled, the pre-cap error scales with the right
exponent for each semantics. **G5's estimator does drive error control** — the "negative
result about G5" branch is triggered by the *metric*, not the estimator.

---

## P2 — conclusions

> [!success] What P2 decides
> - **G4 closed on mechanism.** An explicit `tol`, separate from the reject band, gives a
>   controller that shrinks and holds $e_n\approx\mathrm{tol}$. No band retuning
> - **The PI controller works**, with two fixes before it's a real tolerance:
>   1. a **tol-aware first step** (start at `DELTA_TAU_MIN`, or pick $\Delta\tau_0$ from $e_1$)
>   2. use **EPUS** if we want $\alpha\approx1$
> - **T2′ needs a different metric for QITE.** The end point of a contractive gradient flow
>   erases the controller's work; the max error **before the cap binds** shows it
> - Past $\tau\approx2.1$, `tol` barely matters — the cap decides

---

## The common thread

Both experiments end at the same place: **the late regime belongs to the cond cap**, and the
cap is reacting to the **parametrization**.

- P1: the brake fights a chart fold (th8 near its pole) that the state crosses smoothly —
  86 of 108 steps spent there
- P1: Target A's late collapse is the CRZ row fading as the state converges
- P2: past $\tau\approx2.1$ every run, at every `tol`, is cap-limited

**So the lever isn't a better `tol` or a different scalar in the cap — it's controlling the
step in *state* space.**

---

## Next (none of this has been run)

- **P1b — state-space step control.** Detect "$\dot\theta$ along the soft direction while
  $\sqrt{\dot\theta^\top M\dot\theta}$ stays flat" (both already available per step); bound
  $\Delta\tau\,\|d\rho/d\tau\|$ instead, or freeze th8 across the fold.
  *Guess:* Target C drops from 108 steps to ~22–28 without the cap-off error penalty
- **P2b — pre-registered fixed version.** EPUS, tol-aware first step, pre-cap metric, on a
  **second seed** and Target C. *Guess:* slope in $[0.8,1.1]$, and the work-precision lead
  over Route 3 survives
- **P3 — exact testbeds** (BT two-level qubit, imaginary-time qubit): score the memo's
  estimators against the **true** local error, as curves across $\Delta\tau$

---

## Sources

- Plan: [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]]
- P1: [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/P1_Singularity_Classification_Findings|P1_Singularity_Classification_Findings]] — `Proposal_1/p1_singularity.py`,
  `Output.txt`, **51 pinned, 0 drift**
- P2: [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/P2_Tol_Proportionality_Findings|P2_Tol_Proportionality_Findings]] — `Proposal_2/p2_tol_proportionality.py`,
  `Output.txt`, **93 pinned, 0 drift**
- BT replication: [[QITE_Hand_Off/Notes_on_Project/Replication_Beyond_Trotterization|Replication_Beyond_Trotterization]]
- Audits: [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G4_Controller_Cannot_Shrink|G4_Controller_Cannot_Shrink]], [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot|G5_Error_Estimator_Validated_And_Its_Blind_Spot]],
  [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility|G6_Attribution_And_Reproducibility]]

Both scripts import `adaptive_qite.py` / `adam_qite.py` **without modifying them**; reruns are
byte-identical.
