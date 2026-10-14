---
title: "P2 findings — a real tol for the step controller, and tolerance proportionality"
tags:
  - report
  - physics
  - quantum-computation
  - adaptative-qite
sources:
  - "[[QITE_Hand_Off/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]]"
  - "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G4_Controller_Cannot_Shrink|G4_Controller_Cannot_Shrink]]"
  - "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot|G5_Error_Estimator_Validated_And_Its_Blind_Spot]]"
  - "[[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]"
  - "[[QITE_Hand_Off/Adam_QITE/Report_Adam_QITE|Report_Adam_QITE]]"
  - "[[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/P1_Singularity_Classification_Findings|P1_Singularity_Classification_Findings]]"
code: QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/p2_tol_proportionality.py
---

# P2 findings: does the achieved error track a requested `tol`?

This note reports the run of **P2** from [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]]. Neither of
our controllers has a tolerance to sweep. Route 3 ratchets, and Adam uses `energy_check`'s reject
band as its setpoint, which [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G4_Controller_Cannot_Shrink|G4]] showed is circular: growth $<1$ on 0
of 190 steps. So I gave the step law an explicit `tol`. It's the memo's Eq. (17) PI law driven by
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot|G5]]'s zero-cost estimator
$e_n=|\Delta E_{\rm actual}-\Delta E_{\rm pred}|/2$. `DE_BAND_*` stays as the reject gate. Then I
asked the T2′ question: does the achieved error scale like $\mathrm{tol}^\alpha$?

> [!failure] TL;DR
> - **The acceptance test as pre-registered fails: no fit is possible.** The error at $\tau_{\rm end}$
>   against an exact reference doesn't track `tol` in either arm. Local slopes wander from $-0.80$
>   to $+2.62$. Only 2 runs per arm pass the pre-registered fit criteria, which need 3. So 5 of the
>   6 predictions score FAILED. **Read that carefully:** predictions 1–3, 4 and 6 fail because no
>   slope could be fitted at all, not because a slope was measured out of range.
> - **G4's circularity is gone (prediction 5 held).** With `tol` separate from the band, growth
>   $<1$ happens in every run: 527 of 1234 steps at `tol` $=10^{-6}$. Where `tol` sits below the
>   estimator's natural scale, the controller holds median $e_n/\mathrm{tol}$ at 0.99. So the PI
>   law *does* control $e_n$. It's the end-point error that doesn't follow.
> - **Why, post hoc (not pre-registered).** Two things.
>   - The first step is the fixed `DELTA_TAU0` $=0.05$ in every run. It injects a
>     $5.2914\times10^{-4}$ error whatever `tol` is ($e_1=229\times$ `tol` at $10^{-6}$).
>   - The imaginary-time flow first amplifies errors, with a peak at $\tau\approx1$ in every run,
>     then contracts them by 2–3 decades. So $\tau_{\rm end}$ mostly remembers the late,
>     cap-dominated regime, not what the controller did.
> - **With the first step at `DELTA_TAU_MIN`, the pre-cap error tracks `tol` with the textbook
>   exponents.** On max error over $\tau\le2$, the slope over the 3 smallest tols is **0.478**
>   for error per step (theory 1/2) and **0.890** for error per unit step (theory 1). So G5's
>   estimator does drive error control. The plan's "negative result about the G5 premise" is
>   formally triggered by the pre-registered metric, and I don't think it's the right reading.
> - **The plan's $\alpha\approx1$ was wrong for the arm it named.** For error-per-step control of
>   a $p=1$ method, theory says $\alpha=1/2$. The ablation lands on 0.478.
> - **Work-precision at $\tau_{\rm end}$** looks great for PI (e.g. $5.67\times10^{-6}$ at 2700
>   circuits against Route 3's $6.32\times10^{-4}$ at 2069), but the curve isn't monotone in `tol`.
>   I wouldn't claim a 100× win from one seed.

## 1. The predictions, as written before the run

### The plan's, quoted (P2, "Prediction")

> - Global error $\propto\mathrm{tol}^{\alpha}$ with $\alpha\approx1$ (explicit Euler with local
>   error control gives global error roughly proportional to tol) over the decades where neither
>   the variational floor nor the cond cap binds.
> - The `tol` at which the slope breaks marks where the controller stops controlling. I expect the
>   break to coincide with the cap taking over (Target A's last-quarter regime).
> - If there's **no** decade with slope ≈1, then G5's estimator doesn't actually drive error
>   control even though it correlates with the true LTE. That would be a negative result about
>   the G5 premise.

### A flaw I flagged before running

"Local error control" means two different things, and $\alpha\approx1$ only follows for one of them.

- **Error per step (EPS)**, $e_n\le\mathrm{tol}$. With $e_n\approx C\Delta\tau^2$ you get
  $\Delta\tau\propto\mathrm{tol}^{1/2}$. The memo's own $E_{\rm tot}\approx CT\Delta\tau^p$ then
  gives global $\propto\mathrm{tol}^{p/(p+1)}=\mathrm{tol}^{1/2}$.
- **Error per unit step (EPUS)**, $e_n/\Delta\tau\le\mathrm{tol}$. This gives
  $\Delta\tau\propto\mathrm{tol}$ and global $\propto\mathrm{tol}$.

The plan's P2 section describes the per-step form ($e_n$ against `tol`). I didn't edit the plan's
prediction. It's scored as written against the EPS arm, and I ran both arms.

### Mine, written into the plan before the sweep

1. **(plan, as written)** EPS arm: $\alpha\in[0.8,1.2]$. I expected this to **fail**.
2. EPS arm: $\alpha\in[0.4,0.6]$ over the uncapped decades.
3. EPUS arm: $\alpha\in[0.8,1.2]$.
4. **(plan)** The slope breaks where the cap takes over. The break tol is the first grid point with
   local slope $<\alpha_{\rm fit}/2$. The cap tol is the first with $\ge25\%$ of steps cap-bound.
   The prediction holds if the two are within one grid point.
5. G4's circularity is gone: growth $<1$ on $>0$ steps in every run.
6. EPS circuits $\propto\mathrm{tol}^{-1/2}$: slope in $[-0.6,-0.4]$ over the same decades.

### Scorecard

| # | Prediction | Outcome | Why |
|---|---|---|---|
| 1 | (plan) EPS $\alpha\in[0.8,1.2]$ | **failed** | **no fit possible**: 2 of 9 runs pass the fit criteria |
| 2 | EPS $\alpha\in[0.4,0.6]$ | **failed** | **no fit possible** (same 2 of 9) |
| 3 | EPUS $\alpha\in[0.8,1.2]$ | **failed** | **no fit possible**: 2 of 8 |
| 4 | (plan) break where the cap takes over | **failed** | with no $\alpha_{\rm fit}$ there's no break tol. Cap tol: $10^{-5}$ (EPS), $3.16\times10^{-3}$ (EPUS) |
| 5 | growth $<1$ in every run | **held** | min over all 17 runs is 1 step (EPUS `tol` $=1$). 527 steps at EPS $10^{-6}$ |
| 6 | EPS circuit slope $\in[-0.6,-0.4]$ | **failed** | **no fit possible** |

The two points that *do* pass in each arm are the two smallest tols. The EPS local slope between
them is $+0.47$, inside prediction 2's window, and the EPUS slope is $+0.04$. But two points
aren't a fit, and I set the $\ge3$ rule before running. So these count as failures, and the
pinned outcomes are 0.

## 2. Method

Everything is in `p2_tol_proportionality.py`. It **imports** `adaptive_qite.py`, `adam_qite.py`
and P1's `p1_singularity.py`, and changes none of them. So Targets A–D, the Adam 147-check ledger
and P1's 51 pins all stay valid. The plan wanted `PIStepController` inside `adam_qite.py`. Keeping
it in `Proposal_2/` is the same deviation P1 made, for the same reason.

- **The loop.** `controlled_qite_run` is `adam_qite_run`'s noiseless loop, line for line, with
  only the step law pluggable. Measurement, the cond cap, the `energy_check` gate and retries are
  unchanged.
  - Target `ids` proves the copy is faithful: an Adam adapter reproduces `adam_qite_run` on E1 bit
    for bit. That's final $E=-1.122851752157184$, 190 steps, 2285 circuits, with $E$,
    $\Delta\tau$, growth and circuits identical step by step.
- **The controller.** Memo Eq. (17) at $p=1$:
  $\Delta\tau_{n+1}=\Delta\tau_n(\mathrm{tol}/e_n)^{0.15}(e_{n-1}/e_n)^{0.2}$.
  - The first update is I-only.
  - The factor is clipped to `[SHRINK, GROW]`, both imported, not retuned.
  - EPS uses $e_n=|\Delta E_{\rm actual}-\Delta E_{\rm pred}|/2$. EPUS uses $e_n/\Delta\tau$.
- **The reference.** The plan's §5 said to solve the same fd ODE the controllers integrate. That
  turned out to be impossible at the tolerance I needed: DOP853 at rtol $10^{-8}$ stalls on the fd
  noise, and so do LSODA and Radau. So I switched to the **exact-metric** ODE
  $\dot\theta=\mathrm{pinv}(M_{\rm exact},10^{-6})Y_{\rm exact}$, using P1's
  `exact_metric_and_force` with DOP853 at rtol $10^{-11}$.
  - Self-check against rtol $10^{-13}$: $8.25\times10^{-12}$.
  - The fd-vs-exact ODE gap, from an fd solve at rtol $10^{-7}$, is reported as a floor:
    $3.42\times10^{-8}$ max over $\tau$ and $1.19\times10^{-10}$ over $\tau\in[3,3.6]$.
  - Every measured error sits far above both floors, so the reference is not what limits anything
    below.
- **Error metric (pre-registered).** $|E_{\rm run}(\tau_{\rm end})-E_{\rm ref}(\tau_{\rm end})|$
  at each run's own $\tau_{\rm end}$. The max over steps is secondary.
- **Fit criteria (pre-registered).** A run enters the fit if all of these hold:
  - it exits on $\tau_{\max}$;
  - the cap binds on $<25\%$ of its steps;
  - $\Delta\tau$ sits at a bound on $<25\%$;
  - its error is $>10\times$ the floor.
  
  $\alpha$ is then the log-log least-squares slope over those runs, and it needs at least 3.
- **Grids.**
  - EPS: `tol` $=10^{-6}\ldots10^{-2}$ in half decades, 9 runs.
  - EPUS: `tol` $=10^{-3.5}\ldots1$, 8 runs.
  - `MAX_STEPS_SWEEP` $=5000$.
  - Target A: `SEED_A`, 5 params, `rcond` $10^{-6}$, cap on, $\tau_{\max}=3$.
- **Pinning.** Every number here is pinned at full precision in `PINNED`, with `--check` and no
  tolerance knob. "No fit possible" is pinned too, as `nan`.

Targets: `ids`, `ref`, `eps`, `epus`, `fit` (all pre-registered), and `ablate` (**post hoc**, §5).

## 3. Results: the sweeps

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_tol_proportionality.png]]
*Fig. 1.*
- *(a) Error at $\tau_{\rm end}$ vs `tol`. Filled dots are the runs that pass the fit criteria;
  dotted and dashed lines are slope 1/2 and slope 1.*
- *(b) Work-precision against Route 3 (Target A) and Adam (E1).*
- *(c) $\Delta\tau(\tau)$ for three tols per arm.*
- *(d) Fraction of steps where the cond cap binds, and where $\Delta\tau$ sits at a bound.*

### EPS arm, $e_n\le\mathrm{tol}$

| tol | err($\tau_{\rm end}$) | max err | steps | circuits | cap % | growth $<1$ | GROW clips | med $e_n$/tol | in fit |
|---|---|---|---|---|---|---|---|---|---|
| $10^{-6}$ | $1.1067\times10^{-5}$ | $6.977\times10^{-3}$ | 1234 | 14808 | 10.8 | 527 | 131 | 0.992 | yes |
| $3.16\times10^{-6}$ | $1.9022\times10^{-5}$ | $9.066\times10^{-3}$ | 761 | 9132 | 18.4 | 295 | 145 | 0.974 | yes |
| $10^{-5}$ | $2.5242\times10^{-5}$ | $1.291\times10^{-2}$ | 497 | 5964 | 29.0 | 163 | 154 | 0.897 | no |
| $3.16\times10^{-5}$ | $2.3316\times10^{-5}$ | $1.917\times10^{-2}$ | 351 | 4212 | 41.6 | 91 | 164 | 0.775 | no |
| $10^{-4}$ | $1.4299\times10^{-5}$ | $2.918\times10^{-2}$ | 270 | 3240 | 54.4 | 52 | 165 | $3.35\times10^{-3}$ | no |
| $3.16\times10^{-4}$ | $5.6702\times10^{-6}$ | $4.468\times10^{-2}$ | 225 | 2700 | 65.3 | 29 | 167 | $2.61\times10^{-4}$ | no |
| $10^{-3}$ | $3.6387\times10^{-5}$ | $6.479\times10^{-2}$ | 200 | 2400 | 73.5 | 14 | 166 | $4.75\times10^{-5}$ | no |
| $3.16\times10^{-3}$ | $4.2490\times10^{-4}$ | $8.954\times10^{-2}$ | 184 | 2208 | 79.9 | 6 | 163 | $9.5\times10^{-6}$ | no |
| $10^{-2}$ | $2.5438\times10^{-4}$ | $1.130\times10^{-1}$ | 180 | 2162 | 80.0 | 2 | 169 | $2.9\times10^{-6}$ | no |

Local slopes $d\log\mathrm{err}/d\log\mathrm{tol}$: $+0.47,+0.25,-0.07,-0.42,-0.80,+1.61,+2.13,-0.45$.
Retries: 0, except 2 at $10^{-2}$. $\Delta\tau$ never sat at `DELTA_TAU_MIN`/`MAX`.

### EPUS arm, $e_n/\Delta\tau\le\mathrm{tol}$

| tol | err($\tau_{\rm end}$) | max err | steps | circuits | cap % | growth $<1$ | GROW clips | med $e_n$/tol | in fit |
|---|---|---|---|---|---|---|---|---|---|
| $3.16\times10^{-4}$ | $2.9473\times10^{-5}$ | $8.413\times10^{-3}$ | 2260 | 27120 | 6.2 | 1110 | 144 | 1 | yes |
| $10^{-3}$ | $3.0816\times10^{-5}$ | $1.420\times10^{-2}$ | 822 | 9864 | 17.8 | 339 | 162 | 0.985 | yes |
| $3.16\times10^{-3}$ | $2.4532\times10^{-5}$ | $2.751\times10^{-2}$ | 371 | 4452 | 39.6 | 105 | 167 | 0.633 | no |
| $10^{-2}$ | $1.1781\times10^{-5}$ | $5.366\times10^{-2}$ | 228 | 2736 | 64.9 | 30 | 166 | $2.43\times10^{-3}$ | no |
| $3.16\times10^{-2}$ | $2.4098\times10^{-4}$ | $8.966\times10^{-2}$ | 186 | 2233 | 78.0 | 7 | 164 | $4.07\times10^{-4}$ | no |
| $10^{-1}$ | $4.7659\times10^{-4}$ | $1.209\times10^{-1}$ | 176 | 2115 | 81.8 | 2 | 167 | $1.08\times10^{-4}$ | no |
| $0.316$ | $7.0607\times10^{-4}$ | $1.325\times10^{-1}$ | 172 | 2068 | 83.1 | 2 | 169 | $3.23\times10^{-5}$ | no |
| $1$ | $2.4412\times10^{-3}$ | $1.325\times10^{-1}$ | 56 | 676 | 67.9 | 1 | 54 | $1.67\times10^{-4}$ | no |

Local slopes: $+0.04,-0.20,-0.64,+2.62,+0.59,+0.34,+1.08$.

### What the tables say

- **The controller controls $e_n$ where it can.** For the smallest tols, median $e_n$/tol is 0.99
  (EPS $10^{-6}$) and 1 (EPUS $3.16\times10^{-4}$). That's the setpoint behaviour G4 found
  missing in Adam.
- **Above the estimator's natural scale, `tol` doesn't bind at all.** From EPS $10^{-4}$ and EPUS
  $10^{-2}$ upward, median $e_n$/tol drops to $10^{-3}$ and below. $\Delta\tau$ is then set by the
  GROW clip (131–169 hits per run, 54 in the short EPUS `tol` $=1$ run) and by the cap (54–83%
  of steps). Past that point, sweeping `tol` doesn't change what the integrator does much. Steps
  flatten towards 170–270.
- **The cap binds early in `tol`.** It already covers $\ge25\%$ of steps at EPS $10^{-5}$ and
  EPUS $3.16\times10^{-3}$, the third grid point in each arm. That leaves only 2 fit-eligible runs
  per arm. That's the direct reason there's no fit. The cap never binds before $\tau=2$ (claim in
  `ablate`), and in the ablation sweep it first binds at $\tau=2.091$–$2.409$. This is Target A's
  late "developing stiffness" regime. P1 tied its end to the CRZ row fading.
- **The max error along the run does rise with `tol`.** It's monotone in both arms: $6.98\times10^{-3}$
  up to $0.113$ for EPS. But it's shallow (§5), and that's the clue.

### Existing controllers on the same reference

| Controller | err($\tau_{\rm end}$) | steps | circuits |
|---|---|---|---|
| Route 3 ratchet (Target A) | $6.3181\times10^{-4}$ | 172 | 2069 |
| Adam (E1) | $7.2773\times10^{-4}$ | 190 | 2285 |

At comparable cost, several PI runs end much closer to the reference:
- EPS $3.16\times10^{-4}$: $5.67\times10^{-6}$ at 2700 circuits.
- EPUS $10^{-2}$: $1.18\times10^{-5}$ at 2736 circuits.
- EPS $10^{-3}$: $3.64\times10^{-5}$ at 2400 circuits.

But the end-point error isn't monotone in `tol` (see the local slopes), so part of this could be
where the contraction happens to leave each run. I'd want a second seed before calling it a real
improvement.

### The plan's quantity: final $E-E_{\rm floor}$

$E_{\rm floor}=-1.1246000000000003$ (from `variational_floor`). The **exact** flow is still
$E_{\rm ref}(3)-E_{\rm floor}=1.0212\times10^{-3}$ away at $\tau=3$. Every run lands between
$8.3\times10^{-4}$ and $3.45\times10^{-3}$ (full rows in the output). That quantity is dominated
by how far the flow still has to go, and by each run's slightly different $\tau_{\rm end}$, not by
integration error. It's the wrong metric for proportionality, as I suspected in the plan.

## 4. Why the end-point error doesn't track `tol`

Two things I only saw once the sweep was done.

1. **The first step is tol-blind.** `controlled_qite_run` inherits `DELTA_TAU0` $=0.05$ from
   `adam_qite_run`, and the controller only acts *after* step 1. That step's energy error is
   $5.2914\times10^{-4}$ in **every** default run (pinned as a claim). Its estimator reading is
   $e_1=2.2894\times10^{-4}$, which is $229\times$ `tol` at $10^{-6}$. Everything after that is
   controlled, but the injected error is already there.
2. **The flow amplifies, then forgets.** In every run the trajectory error peaks at
   $\tau\approx0.94$–$1.15$ and then decays by 2–3 decades before $\tau=3$ (Fig. 2b). In
   imaginary time the late flow is a contraction towards the ground state of the ansatz. So the
   error at $\tau_{\rm end}$ mostly reflects the last stretch, which is cap-dominated in every run
   past the second grid point. That's why the end-point error sits at $10^{-5}$ across four decades
   of `tol`, and why its sign (and so the local slopes) jumps around.

Point 2 is an interpretation of Fig. 2b, not a separate measurement. I didn't compute a
contraction rate.

## 5. Post hoc: controlling the first step (not pre-registered)

> [!warning] Post hoc
> This target was added **after** the fit target's results, to test the explanation in §4. It
> doesn't rescue predictions 1–6, which stay as scored. The window $\tau\le2$ and "slope over the
> 3 smallest tols" are choices I made after seeing the data. The window is justified by the cap
> never binding there (checked as a claim in both sweeps). In the first-step-at-min sweep the
> first cap-bound step is at $\tau\ge2.091$.

Same sweep, with one change: the first step is `DELTA_TAU_MIN` $=10^{-4}$ instead of 0.05, so the
controller sets $\Delta\tau$ from the start. The score is the max error over $\tau\le2$.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_first_step_ablation.png]]
*Fig. 2.*
- *(a) Max error over $\tau\le2$ vs `tol`. Open markers: fixed first step. Filled: first step at
  minimum. Grey guides show slopes 1/2 and 1.*
- *(b) The error along the run for EPS $10^{-6}$ and $10^{-4}$. The fixed first step shows up as a
  jump at $\tau=0.05$; after the peak at $\tau\approx1$, the flow contracts the error.*

| arm | slope of max err, 3 smallest tols: fixed first step | … first step at min | theory | circuit slope (first step at min) | theory |
|---|---|---|---|---|---|
| EPS | 0.267 | **0.478** | 1/2 | $-0.385$ | $-1/2$ |
| EPUS | 0.515 | **0.890** | 1 | $-0.765$ | $-1$ |

Local slopes of the max error with the first step at min:
- EPS: $+0.48,+0.47,+0.46,+0.43,+0.40,+0.34,+0.26,+0.19$.
- EPUS: $+0.94,+0.84,+0.71,+0.49,+0.24,+0.05,0.00$.

Both roll off gradually above the tols where median $e_n$/tol leaves ~1 in §3. That's where `tol`
stops binding, as it should.

- **So G5's estimator does drive error control.** With the first step controlled, the error in the
  region the cap doesn't touch scales with the textbook exponent for each semantics. The plan's
  negative-result branch ("no decade with slope ≈1 → the G5 premise fails") is triggered on
  paper by the pre-registered metric. But I read it as a failure of the **metric** (the end point
  of a contractive flow) plus a **tol-blind first step**, not of G5.
- **Circuits scale more weakly than theory** ($-0.385$ against $-1/2$). My guess, not tested: the
  last third of every run is cap-limited, so it's a roughly tol-independent step overhead that
  flattens the slope.
- **The end point barely moves.** err($\tau_{\rm end}$) with the first step at min is still
  $2.0\times10^{-5}$–$5.1\times10^{-5}$ across the small tols (table in the output), so point 2 of
  §4 stands independently of point 1.
- **Size of the fixed-step damage.** At EPS $10^{-6}$ the max error drops from $6.977\times10^{-3}$
  to $2.469\times10^{-3}$ just by controlling step 1. At EPUS $3.16\times10^{-4}$ it drops from
  $8.413\times10^{-3}$ to $2.057\times10^{-3}$.

## 6. What this decides

- **G4 (the setpoint question): closed on mechanism.** An explicit `tol`, separate from the reject
  band, gives a controller that shrinks when it should (growth $<1$ in every run) and holds
  $e_n\approx\mathrm{tol}$. No band retuning was needed. G4's item 2 (feeding rejected trials to
  the controller) turned out to be unnecessary.
- **R3 (the PI controller): works as a controller, with two fixes needed before it's a
  tolerance.**
  1. The first step has to be tol-aware: start at `DELTA_TAU_MIN`, or pick $\Delta\tau_0$ from
     $e_1$.
  2. Default to **EPUS** if we want $\alpha\approx1$. EPS gives the $\alpha=1/2$ it should.
- **T2′ as an acceptance test needs a different metric for QITE.** The error at $\tau_{\rm end}$
  of a contractive gradient flow can't show tolerance proportionality, because the flow erases
  it. The max trajectory error before the cap binds does show it. For P3/P5 I'd pre-register
  that metric.
- **The cap still decides the late regime.** Past $\tau\approx2.1$, `tol` barely matters. That's
  P1's chart/zero-variance story again, and P1's "state-space step control" idea is what would
  address it. It's not a `tol` question.
- **G5's blind spot ($\partial_\theta M^+$)** doesn't show up as a failure here. Target A is full
  rank and noiseless, so this run doesn't probe it.

### What I'd run next (not run)

A **P2b** that pre-registers the fixed version: EPUS, first step at `DELTA_TAU_MIN`, max error over
the cap-free window, on a **second seed** plus Target C.
- **Prediction** (a guess): slope in $[0.8,1.1]$ over the 3 smallest tols on both, and the
  work-precision lead over Route 3 survives the second seed at matched circuits.

## 7. Assumptions (ours)

These are the 10 items in `ASSUMPTIONS` in the script, in short:

1. Every physical choice is imported from `adaptive_qite.py` unchanged: H$_2$, the ansatz,
   `SEED_A`, `drho` at $h=10^{-6}$, `rcond` $10^{-6}$, the cap (2000), `energy_check` and
   `DE_BAND_*`, GROW/SHRINK, `DELTA_TAU0/MIN/MAX`, retries, $\tau_{\max}$.
2. Two tolerance semantics (EPS/EPUS). The plan's $\alpha\approx1$ is scored against EPS as
   written.
3. $e_n=|\Delta E_{\rm actual}-\Delta E_{\rm pred}|/2$, the G5 estimator. It's an energy-space
   proxy, not a $\theta$-space error.
4. Eq. (17) at $k_I=0.15$, $k_P=0.2$. The first update is I-only, the factor is clipped to
   `[SHRINK, GROW]`, and it updates once per accepted step, like the Adam controller.
5. Rejection stays `energy_check`'s band. Reject-on-$e_n>$ tol (the noiseless Eq. 14) is untested.
6. The reference is the exact-metric ODE (P1's `exact_metric_and_force`), not the fd ODE, which
   can't be solved tightly enough. The fd-vs-exact gap is reported as a floor.
7. The primary error is at $\tau_{\rm end}$ and the secondary is the max over steps. The plan's
   $E-E_{\rm floor}$ is printed for the record.
8. The fit criteria (§2) were fixed before the sweep. The $10\times$ floor criterion was added after
   sizing the reference and before any controller run.
9. Break tol / cap tol are defined as in prediction 4.
10. No Adam arm: the plan's precondition (the G2 fix, $\beta_1=\beta_2$) isn't met.

Post hoc, for §5 only: the first step at `DELTA_TAU_MIN`, the window $\tau\le2$, and slopes over
the 3 smallest tols.

## 8. Provenance

Every number in this note comes from `Proposal_2/Output.txt`, produced by
`~/micromamba/envs/qml-julia/bin/python p2_tol_proportionality.py --target all`. A second full run
was byte-identical, and `--check` exits 0: 93 pinned values and claims, 0 drifted.

| Numbers | Output section |
|---|---|
| Adam replay $-1.122851752157184$ / 190 / 2285, step-by-step identity | `ids` |
| $E_{\rm ref}(3)$, self-check $8.25\times10^{-12}$, floors $3.42\times10^{-8}$ / $1.19\times10^{-10}$, DOP853 evaluation counts | `ref` |
| EPS table, local slopes | `eps`, `fit` |
| EPUS table, local slopes | `epus`, `fit` |
| fit counts, cap tols, predictions, Route 3 / Adam anchors, $E_{\rm floor}$, $E-E_{\rm floor}$ rows | `fit` |
| first-step error $5.2914\times10^{-4}$, $e_1$, ablation table, slopes 0.478 / 0.890 / $-0.385$ / $-0.765$, first cap $\tau$, peak $\tau$ | `ablate` |
| G4's 0/190, G5's ratio 2.060 | [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G4_Controller_Cannot_Shrink\|G4_Controller_Cannot_Shrink]], [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot\|G5_Error_Estimator_Validated_And_Its_Blind_Spot]] (via the plan) |
