---
title: "Presenter notes — Beyond Trotterization into Route 3: P1 & P2"
tags:
  - presentation
  - physics
  - quantum-computation
  - adaptative-qite
deck: "[[QITE_Hand_Off/Presentations/Adaptative_QITE_Beyond_Trot_P1_P2|Adaptative_QITE_Beyond_Trot_P1_P2]]"
---

# Presenter notes — P1 & P2 deck

Companion to [[QITE_Hand_Off/Presentations/Adaptative_QITE_Beyond_Trot_P1_P2|Adaptative_QITE_Beyond_Trot_P1_P2]]. Written for me coming back to
this cold: nothing is assumed, every term and number on the slides is explained here.

> [!tip] How to use this note
> - §1–§4 are **read-before**: the 60-second version, the backstory, the refresher on the earlier
>   experiments (Targets A–D, E1–E3), the glossary.
> - §5 is **slide by slide**, numbered exactly like the deck (26 slides). Each slide has:
>   **Point** (one line) → **Say** (what to say) → **If you forgot** (background) →
>   **Reading the figure** (figure slides) → **Likely questions** → **Don't overclaim**.
> - §6 is a bank of hard questions. §7 is the numbers I must not misquote.
>
> Sources for everything: [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/P1_Singularity_Classification_Findings|P1_Singularity_Classification_Findings]],
> [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/P2_Tol_Proportionality_Findings|P2_Tol_Proportionality_Findings]], [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]], plus the
> background notes linked where used. No number here is new.

---

## 1. The 60-second version

If I only get one minute, say this:

1. Our Route 3 step controller has a safety brake: when cond($M$) gets big, it shrinks the time step.
2. **P1** asked what that big cond($M$) actually is. Answer: not a bug, not hard physics — it's a
   **coordinate problem**. One angle (th8) stops moving the state, so the solver has to spin it
   very fast, and the brake reacts to that. The state itself moves smoothly the whole time.
3. Consequence: our old E3/Target C results stand, no re-runs needed, quote them to 2–3 sig figs.
4. **P2** gave the controller a real tolerance knob (`tol`) and asked whether the error follows it.
   As pre-registered: **no** — the final error ignores `tol`.
5. But the controller itself works, and the reasons are understood: a fixed first step, and an
   imaginary-time flow that "forgets" errors. Fix the first step and measure before the brake
   kicks in, and the textbook slopes appear (0.478 vs theory 1/2, 0.890 vs theory 1) — post hoc.
6. Common thread: late in every run the brake decides everything. The next lever is to control
   the step in **state space**, not parameter space. Not run yet.

---

## 2. Backstory you need in your head

### Timeline

| Step | What happened |
|---|---|
| 1 | **QDP replication** — reproduced Saxena–Shalabi–Kyaw's noise paper on H$_2$. Found the 9-param ansatz's $M$ is *exactly* rank 6 of 9 ([[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP\|Replication_Practical_Limitations_QDP]]) |
| 2 | **Route 3 controller built** — the adaptive step-size controller ([[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation\|Adaptive_QITE_Attack_Route_3_Implementation]]) |
| 3 | **First audit** (M1–M6) — fixed a wrong debias factor, honest circuit counting, etc. ([[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE\|Audit_Adaptative_QITE]]) |
| 4 | **Adam controller built** — a smarter step law on top of Route 3 (`Adam_QITE/`) |
| 5 | **Adam audit** (G1–G7) — G4, G5, G6 are the ones this deck uses ([[QITE_Hand_Off/Adam_QITE/Report_Adam_QITE\|Report_Adam_QITE]]) |
| 6 | **BT replication** — reproduced the Beyond Trotterization paper ([[QITE_Hand_Off/Notes_on_Project/Replication_Beyond_Trotterization\|Replication_Beyond_Trotterization]]) |
| 7 | **Plan** with 7 proposals (P1–P7) mapping BT ideas onto our controllers; **P1 and P2 run** from it |

### The pieces, in one or two lines each

- **The project.** How far can variational QITE go under hardware noise, and can an adaptive
  step-size controller help? Work with the LG Toronto AI team (Thi Ha Kyaw, Jack Baker, Gaurav
  Saxena).
- **The memo.** An internal 5-page note (`Research_Papers/QITE_timestep.pdf`, no arXiv/DOI) that
  proposes 5 "Routes" to build a step controller. Its core idea: on noisy hardware there's a noise
  floor, so the controller should take steps **as large as noise permits**.
- **Route 3** = the memo's "cheap physics indicators" route. The only route we've actually built.
- **Adam controller** = our own extension: replaces Route 3's crude ×1.25 / ×0.5 step ratchet with
  an Adam-like law. Mostly a wash; its audit found the design problems G2–G6.
- **BT** = *Beyond Trotterization* (Assi, Vogl, Kumari & LeBlanc 2026, PRB 113, 214314,
  arXiv:2511.15124). A **real-time** simulation paper, but its equation of motion
  $g\dot c=-b$ is the same shape as ours, $M\dot\theta=Y$. So its tricks for the matrix transfer.
- **P1–P7** = the 7 proposals in [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]]. Only P1 and P2 have
  been run. P3–P7 are ideas.

---

## 3. Refresher: the earlier experiments (Targets A–D, E1–E3)

P1 and P2 don't start from scratch. They re-open questions left by earlier runs. Here's what
those runs were.

- **Targets A–D** test **Route 3 on its own**. Built at step 2, fixed and re-verified after the
  first audit. Code: `Attack_Route_03/adaptive_qite.py`.
- **E1–E3** compare **Adam vs Route 3**. Built at step 4, audited at step 5. Code:
  `Adam_QITE/adam_qite.py`.
- Each E reuses the seed of a Target, so both start from the **same $\theta_0$**:
  **E1↔A, E2↔D, E3↔C**.
- Every system is H$_2$ on 2 qubits.

> [!summary] Target A: noiseless, 5 parameters, adaptive vs fixed step
> - **What:** Route 3's adaptive $\Delta\tau$ vs a fixed $\Delta\tau=0.05$. Healthy, full-rank
>   5-param ansatz, `rcond` $10^{-6}$, `SEED_A`.
> - **Why:** the basic question: does adapting the step help at all?
> - **Outcome (post-audit): yes.**
>   - At the honest budget of 732 circuits: adaptive $-1.1210$ vs fixed $-1.1167$.
>   - Measured from the ansatz's own best energy ($E_{\rm floor}=-1.1246$), adaptive lands
>     **2.2× closer at a matched 732 circuits**: $3.6\times10^{-3}$ vs $7.9\times10^{-3}$ (the $1.6\times10^{-3}$ / 5× figure is Route 3's end point after 2069 circuits).
>   - Full run: 172 steps, 2069 circuits.
> - **The loose end:** late in the run $\Delta\tau$ collapses.
>   - The cond cap binds on **31 of the last 43** steps.
>   - Spectral cond reaches $1.4\times10^6$.
>   - At the time it was read as "real developing stiffness".
> - **In this deck:**
>   - P1 slide 11 explains that late stiffness: the CRZ gate's row fades as the state converges
>     to ~$|10\rangle$.
>   - P2 runs entirely on Target A's setup.

> [!summary] Target B: depolarizing noise, with and without debiasing
> - **What:** Noise rescales the measured $\dot\theta$ by a known factor. "Debiasing" divides it
>   back out.
> - **Outcome:**
>   - The original "surprise" (debiasing didn't help) was a **bug**: a wrong debias factor.
>     It's retracted.
>   - After the fix, per-step debiasing is exact. The noisy flow is exactly collinear with the
>     clean one.
>   - Over a whole trajectory, debiased and biased still land within ~2 % of each other. The
>     controller absorbs a pure rescale as a time reparametrization.
> - **In this deck:** not used.

> [!summary] Target C: the 9-parameter, rank-6-of-9 ansatz ⭐
> - **What:** Same noiseless H$_2$, but a **9-parameter** ansatz.
>   - A 2-qubit pure state has only 6 real dimensions, so this $M$ is **always** rank 6 of 9.
>     3 directions are dead by construction.
>   - Run at two pseudoinverse cutoffs: **tight** `rcond` $10^{-8}$ and **loose** $10^{-2}$.
>   - Seed `SEED_C`.
> - **Why:** to see how an adaptive controller behaves on a matrix that's singular by
>   construction, and whether the cutoff choice matters.
> - **Outcome (post-audit):**
>
>   | `rcond` | error | steps | cond($M$) range |
>   |---|---|---|---|
>   | $10^{-8}$ (tight) | $8.9\times10^{-4}$ | 108 | 32.9 – $9.7\times10^6$ |
>   | $10^{-2}$ (loose) | $2.1\times10^{-3}$ | 17 | 19.1 – 125 |
>
>   - Tight wins by **~2.4×** in error but takes ~6× more steps.
>   - Tight's cond spikes to ~$10^7$ near $\tau\approx0.2$ and the cap throttles $\Delta\tau$ hard.
>   - Loose never looks stiff.
>   - *Superseded pre-fix numbers, don't quote:* tight error $1.2\times10^{-2}$ in 20 steps;
>     loose stalled at 500 steps with $2.9\times10^{-2}$. Both were artefacts of the old energy
>     gate over-throttling.
> - **The loose end it left: P1's question.** Is that $10^7$ real? Or is it finite-difference noise
>   leaking into the 3 dead directions? And why does loose look so healthy?
> - **In this deck:** P1 slides 6–12 are about this run. Answer: no leakage, it's a chart
>   singularity (th8). And loose looks healthy only because it throws away real directions.

> [!summary] Target D: noise-floor sweep
> - **What:** sweep the depolarizing $p$ and look for a "noise floor" below which shrinking
>   $\Delta\tau$ stops helping.
> - **Outcome:**
>   - The original "$p\approx10^{-3}$ floor" was **retracted**. It sat on a broken comparison
>     metric plus the Target B bug.
>   - After the fix, deviation grows smoothly with $p$, from $6.4\times10^{-4}$ at $p=10^{-4}$ to
>     $\approx0.29$ at $p=0.05$. There's no knee.
> - **In this deck:** not used. The real noise-floor test is planned for P3c.

> [!summary] E1: Adam vs Route 3 vs fixed step, on Target A's setup
> - **What:** the Adam-style step law against Route 3's ×1.25 / ×0.5 ratchet. Noiseless,
>   5 params, `SEED_A`.
> - **Outcome: a wash.**
>   - Adam: $1.748\times10^{-3}$ above the floor, in 190 steps / 2285 circuits.
>   - Route 3: $1.607\times10^{-3}$ above the floor, in 172 / 2069.
>   - Both crush the fixed step.
> - **The audit finding born here, G4:** Adam's law grew the step on every one of its 190 steps,
>   so it never once shrank (growth $<1$ on 0/190). That's circular by construction: it used the
>   reject band as its setpoint, and it only ever saw steps that had already passed the band.
> - **In this deck:**
>   - P2 exists to fix G4 (slide 16).
>   - Route 3's and Adam's E1 points are the green square and pink diamond in P2's
>     work-precision plot (slide 20).

> [!summary] E2: Adam vs Route 3 under noise
> - **Outcome:** the verdict was **withdrawn** at the audit (G3). The depolarizing channel here is
>   deterministic, so there was no random noise for Adam's averaging to average.
> - **In this deck:** not used.

> [!summary] E3: Adam vs Route 3 on Target C's setup ⭐
> - **What:** the 9-param, rank-6-of-9 ansatz (`SEED_C`). There are **4 arms**: each controller
>   with the cond cap **on** and **off**.
>   - The main comparison is at `rcond` $10^{-8}$.
>   - Repeated at $10^{-2}$.
> - **Why:** to see whether a better step law beats Route 3 on the stiff problem, and to separate
>   the law's effect from the cap's (that's what the cap-off arms are for).
> - **Outcome (from the Adam README):**
>
>   | `rcond` | arm | error | steps |
>   |---|---|---|---|
>   | $10^{-8}$ | Route 3, cap on | $8.868\times10^{-4}$ | 108 |
>   | $10^{-8}$ | Route 3, cap off | $1.398\times10^{-3}$ | 22 |
>   | $10^{-8}$ | Adam, cap on | $8.868\times10^{-4}$ | 108 |
>   | $10^{-8}$ | **Adam, cap off** | $\mathbf{5.653\times10^{-4}}$ | **28** |
>   | $10^{-2}$ | Route 3, cap on = off | $2.143\times10^{-3}$ | 17 |
>   | $10^{-2}$ | Adam, cap on = off | $3.314\times10^{-3}$ | 14 |
>
>   - **The cap-on arms are identical.** At cond ~$9.66\times10^6$ the cap forces
>     $\Delta\tau=2000/\mathrm{cond}=2.07\times10^{-4}$. That overrides both step laws the same
>     way. It's a sanity check, not a bug.
>   - **At $10^{-2}$ the cap never binds**, since cond $\le125$. So cap-on = cap-off.
>   - **The headline:** Adam without the cap wins ($5.653\times10^{-4}$ in 28 steps). But removing
>     the cap *alone* makes Route 3 **worse** ($1.398\times10^{-3}$).
> - **The audit's downgrade (G6):** the win is real, but it's **credited to removing the cap, not
>   to Adam**.
>   - 4 of the 8 rows are duplicate trajectories.
>   - There's no circuit-cost column.
>   - There's no matched-mean-step control arm.
>   - It's $n=1$: one seed, one ansatz.
> - **The loose end:** G6's reading only holds if the logged cond($M$) is trustworthy. If fd noise
>   leaked into the dead directions, E3 would be contaminated and need re-running.
> - **In this deck:**
>   - P1 settles it: no leakage, so G6 stands and E3 doesn't need re-running (slides 6, 7, 14).
>   - Quote E3 errors to 2–3 sig figs. The ranking ($5.65$ / $8.87$ / $13.98\times10^{-4}$) still
>     stands (slide 10).

### How the earlier runs feed this deck

| Earlier run | Question it left open | Answered by | Slides |
|---|---|---|---|
| Target C | Is the $10^7$ cond real or fd leakage? Why does loose `rcond` look healthy? | P1 | 6–12 |
| E3 + G6 | Does the cap read real conditioning, so the cap-based attribution holds? | P1 | 6, 7, 10, 14 |
| Target A | What is the late-run stiffness? | P1 (side finding) | 11 |
| E1 + G4 | Can a controller with a real setpoint shrink, and hit a `tol`? | P2 | 16–23 |
| Target A / E1 | Baseline cost and accuracy for work-precision | P2 | 20 |

---

## 4. Glossary (cheat sheet)

### The physics

| Term | Plain words | Why it matters here |
|---|---|---|
| **QITE** | Imaginary-time evolution $e^{-\tau H}$: any state slides down to the ground state as $\tau$ grows. Not unitary, so we fake it with a parametrized circuit | The whole project |
| **McLachlan's principle** | "Pick the $\dot\theta$ that makes the circuit state follow the true imaginary-time flow as closely as possible." Gives a linear system | Produces $M\dot\theta=Y$ |
| $M\dot\theta=Y$ | $M$ = how much each parameter moves the state (a metric). $Y$ = basically minus the energy gradient. $\dot\theta$ = how fast to turn each angle | Every step solves this |
| $\dot\theta$ | Parameter velocity: how fast the circuit angles change per unit $\tau$ | P1's "parameter speed" |
| **explicit Euler** | $\theta_{n+1}=\theta_n+\Delta\tau\,\dot\theta_n$. Order $p=1$: local error per step $\sim\Delta\tau^2$ | The integrator everything uses |
| **H$_2$, 2 qubits** | Hydrogen molecule, STO-3G basis, Bravyi–Kitaev mapping → 2 qubits | The only system so far ($n=1$ problem) |
| **pure 2-qubit state** has $2d-2=6$ real dimensions | ($d=4$.) Any 2-qubit pure state is described by 6 real numbers | Why the 9-param ansatz can only ever use 6 directions |

### The linear algebra

| Term | Plain words | Why it matters here |
|---|---|---|
| **singular values** $\sigma_1\ge\dots\ge\sigma_9$ | How strongly $M$ responds along each of its natural directions. Tiny $\sigma$ = a direction that barely moves the state | P1 Fig. 1(a) plots all 9 |
| **right-singular vectors** $v_i$, **loadings** | The direction in parameter space that goes with $\sigma_i$. Loadings = how much each angle contributes | P1 Fig. 2 heatmap |
| **cond($M$)** | Ratio big/small singular values: how badly a tiny error in $Y$ gets amplified in $\dot\theta$. cond$_F$ = Frobenius version (the logged one), spectral = $\sigma_1/\sigma_{\min}$ | The thing the brake reads |
| **pinv**, **`rcond`** | Pseudoinverse: solves $M\dot\theta=Y$ while *ignoring* directions with $\sigma_i<\texttt{rcond}\times\sigma_1$. `rcond` is that cutoff | $10^{-8}$ = tight (keeps nearly everything), $10^{-2}$ = loose (throws a lot away) |
| **rank 6 of 9** | The 9-param ansatz has 9 angles but the state only has 6 dims, so 3 directions are *always* dead | Structural, not numerical |
| **gauge direction** | One of those 3 dead directions: changing the angles that way doesn't change the state at all (e.g. turn two identical gates opposite ways) | Harmless *if* `pinv` drops them |
| **gauge leakage** | The worry: numerical noise makes a dead direction look slightly alive, so `pinv` keeps it and the solve goes wrong | P1 suspect #1 |
| **chart singularity** | The coordinates (angles) misbehave but the thing they describe is fine. Like gimbal lock, or longitude at the North Pole | P1's main finding |
| **Gram matrix** | A matrix built from inner products of vectors: $M_{jk}=\mathrm{Tr}[\partial_j\rho\,\partial_k\rho]$ | Explains why fd noise enters the null space *squared* |

### How $M$ and $Y$ are computed

| Term | Plain words |
|---|---|
| **finite differences, `drho`** | Approximate $\partial_k\rho$ by nudging angle $k$ by $h=10^{-6}$ each way and subtracting. Cheap, but carries ~$10^{-10}$ error per entry |
| **exact metric** (rotated generators) | BT's trick: the derivative of a gate is "the gate's generator, rotated by the gates before it", $\tilde G_k=P_kG_kP_k^\dagger$. Gives $M$, $Y$ exactly, no $h$. Function `exact_metric_and_force` in P1's script |
| **state form of $M$** | $M_{jk}=2(\mathrm{Re}\langle\tilde G_j\tilde G_k\rangle-\langle\tilde G_j\rangle\langle\tilde G_k\rangle)$ = twice the covariance of the rotated generators. A generator with zero variance in the current state → zero row in $M$ |

### The controllers

| Term | Plain words |
|---|---|
| **Route 3's 3 indicators** | (1) McLachlan residual $\|M\dot\theta-Y\|$ — a caution flag that freezes step growth; (2) energy gate `energy_check` — predict $\Delta E$, accept the step only if the measured $\Delta E$ is within a band (`DE_BAND_*`); (3) **cond cap** |
| **cond cap** ("the brake") | Before a step: $\Delta\tau\le 2000/\mathrm{cond}(M)$. 2000 is `STIFFNESS_CAP_COEFF`, never retuned in P1/P2. At cond $9.66\times10^6$ this forces $\Delta\tau\approx2.07\times10^{-4}$ |
| **GROW / SHRINK** | Route 3's ratchet: ×1.25 after a clean step, ×0.5 after a rejected one. Hand-tuned. In P2 they survive only as clip bounds |
| **Adam controller** | Our step law using Adam-style moving averages of an error estimate. Shares all Route 3 knobs; differs only in how $\Delta\tau$ evolves |
| **PI controller** (memo Eq. 17) | Classic adaptive-step law: grow $\Delta\tau$ when the error is below `tol`, shrink it when above. "I" part pushes to the setpoint, "P" part damps sudden changes. Gains $k_I=0.15$, $k_P=0.2$ = memo's $0.3/(p+1)$, $0.4/(p+1)$ at $p=1$ |
| **`tol`** | The error per step we *ask* for. New in P2 |
| **$e_n$** | The error estimate per step: $\tfrac12|\Delta E_{\rm actual}-\Delta E_{\rm pred}|$. Free (zero extra circuits). Validated by G5 |
| **EPS** vs **EPUS** | Error per step ($e_n\le$ `tol`) vs error per unit step ($e_n/\Delta\tau\le$ `tol`). Theory: global error $\propto$ `tol`$^{1/2}$ for EPS, $\propto$ `tol`$^{1}$ for EPUS |
| **$\alpha$** | The measured exponent: error $\propto$ `tol`$^\alpha$. "Tolerance proportionality" means a clean straight line on log–log |
| **circuits (M6 currency)** | Honest cost: every circuit actually run, including retries. From the first audit ([[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting\|M6_Honest_Cost_Accounting]]) |

### The test scenarios

Short version. The full story of each run, with outcomes, is in **§3**.

| Name | What |
|---|---|
| **Target A** | Noiseless, **5-param** ansatz (full rank), `rcond` $10^{-6}$, `SEED_A`. Route 3: 172 steps, 2069 circuits. Variational floor $E_{\rm floor}=-1.124600$ (best energy this ansatz can reach) |
| **Target C** | Noiseless, **9-param** ansatz (rank 6 of 9), `rcond` $10^{-8}$ or $10^{-2}$. Route 3 at $10^{-8}$: error $8.868\times10^{-4}$ in 108 steps; at $10^{-2}$: $2.143\times10^{-3}$ in 17 |
| **E1** | Adam vs Route 3 on Target A's setup. Adam: 190 steps / 2285 circuits. A wash |
| **E3** | Adam vs Route 3 on Target C's setup, **4 arms** (each controller with cap on / cap off) at `rcond` $10^{-8}$ |

### The audit findings used here (Adam audit)

| | Plain words |
|---|---|
| **[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G4_Controller_Cannot_Shrink\|G4]]** | Adam's step law could **never shrink** the step (0 of 190 steps). Reason: it used the reject band as its target, and it only ever saw steps that had already passed the band. Circular by construction |
| **[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot\|G5]]** | The free error estimate $|\Delta E_{\rm actual}-\Delta E_{\rm pred}|$ is real: it tracks proper step-doubling with log-Pearson 0.9902, ratio 2.060 (predicted 2). Hence the $\tfrac12$. Blind spot: changes in $M^+$ |
| **[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility\|G6]]** | E3's "Adam wins without the cap" result is real but **credited to removing the cap**, not to Adam. And it's one seed ($n=1$) |

### From BT and the method

| | Plain words |
|---|---|
| **T6** | BT replication target: their matrix $g$ is singular at $t=0$ ("boundary layer") vs our structural rank deficiency. Taught us to **classify** a singularity by looking at *which directions* are weak |
| **T2′** | BT replication target: RK45's error tracks its requested `rtol` linearly over nine decades. That's what a working adaptive integrator looks like → P2's test |
| **pre-registered** | Prediction + scoring rule written in the plan *before* running. Scored as written even when it fails |
| **post hoc** | Analysis added *after* seeing results. Allowed, but labelled and can't rescue a failed prediction |
| **pinned / `--check`** | Every reported number is stored at full precision in the script; `--check` re-runs and fails on any change, no tolerance. P1: 51 pinned, P2: 93, both 0 drift |
| **reference solution** | For P2: the "true" trajectory, from the exact-metric ODE solved with a high-accuracy solver (DOP853, rtol $10^{-11}$). Errors are measured against it |

---

## 5. Slide by slide

### Slide 1 — Title

**Point:** two experiments, both born from the BT replication.

**Say:** "This deck covers two experiments. Both take ideas from the Beyond
Trotterization paper and point them at our own step-size controller. P1 asks what our brake is
really reacting to. P2 asks whether the controller can hit an error target we give it."

---

### Slide 2 — What this deck covers

**Point:** set up the two questions and the format.

**Say:** "Same format for each: the question, what I predicted *before* running, the figures, and
what it decides. Predictions are written down first so they can fail — and several did."

> [!info] If you forgot
> "Slams $\Delta\tau$ shut" = the cond cap forces tiny steps. Target C / E3 are the 9-parameter
> runs where cond($M$) spikes to ~$10^7$ (§4 glossary; full story of both in §3).

---

### Slide 3 — Shared setup

**Point:** remind everyone what the controller is solving and what the brake is.

**Say:** "Everything is on H$_2$, two qubits. Each step we solve $M\dot\theta=Y$ with a
pseudoinverse and take an Euler step. Route 3 has a brake: the step can't exceed 2000 over
cond($M$). Two test setups: Target A with 5 parameters, healthy matrix; Target C with 9
parameters, where 3 of the 9 directions are dead by construction."

> [!info] If you forgot
> - **What Targets A and C were, and how they came out:** refresher §3.
> - **Why rank 6 of 9:** a 2-qubit pure state has only $2d-2=6$ real degrees of freedom
>   ($d=4$). Nine angles can't all do something different. Found in the QDP replication.
> - **`rcond` choice:** $10^{-8}$ keeps nearly everything (sees the stiffness), $10^{-2}$ throws a
>   lot away (looks healthy). This matters on slide 11.
> - **fd at $h=10^{-6}$:** see glossary, "finite differences".

> [!question] Likely questions
> - *Why 2000?* It's Route 3's constant (`STIFFNESS_CAP_COEFF`). Neither P1 nor P2 retuned it,
>   on purpose — retuning after seeing results is a fudge.

---

### Slide 4 — What BT gave us

**Point:** BT is real-time and we're imaginary-time, but the matrix is the same kind of object.

**Say:** "BT's equation is the same linear solve as ours, so their matrix tricks carry over.
Three things we took: an exact way to compute $M$ — no finite differences; a way to classify
singularities by looking at which directions are weak; and the tolerance-proportionality test,
which is how you check an adaptive integrator really controls error."

> [!info] If you forgot
> - The formula on the slide = $M$ is twice the covariance of the "rotated generators" in the
>   current state. I derived it by hand in the plan, and P1 checked it numerically against the fd
>   $M$: max difference $2.594\times10^{-10}$ (that's just fd noise). Sign of $Y$ also verified
>   ($3.955\times10^{-10}$; the opposite sign gives 1.589, so it's unambiguous).
> - What does **not** transfer: BT's headline (one parameter set for any initial state) dies in
>   imaginary time — the state-independent force is a trace of a commutator, which is zero. QITE
>   must stay state-based.

---

### Slide 5 — Part I divider

**Say:** "Part one: what is cond($M$) measuring?"

---

### Slide 6 — P1, the question

**Point:** one number, three possible causes, and they need different responses.

**Say:** "The brake reads one number. When it's huge, there are three stories. One: a bug —
finite-difference noise leaking into the dead directions. If so, E3 and Target C are
contaminated and must be re-run. Two: real stiffness in the physics. Three: the physics is fine
but our angles are a bad coordinate system there — like gimbal lock. My prediction was: no
leakage, it's real stiffness inside the 6 live directions, and the weak direction is a pair of
nearly duplicate gates."

> [!info] If you forgot
> - **What Target C and E3 were, and how they came out:** refresher §3. In short: Target C =
>   the 9-param run at tight vs loose `rcond`. E3 = Adam vs Route 3, cap on/off, on that same setup.
> - **Why it matters:** the Adam audit (G6) credited E3's result to the cap. That only holds if
>   the logged cond is trustworthy. And the right *response* differs: smaller steps fix stiffness,
>   but are arguably the wrong medicine for a coordinate artefact.
> - **The prediction's logic:** fd error ~$10^{-10}$ per entry ≪ cutoff $10^{-8}\sigma_1$, so the
>   dead directions should sit around $10^{-10}$–$10^{-12}$ and get dropped.

---

### Slide 7 — P1, no leakage (Fig. 1 a, b)

**Point:** the dead directions are far below the cutoff; the brake reads the real 6-dim matrix.

*Target C / E3 refresher: §3.*

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_target_c_spectrum.png|560]]

**Say:** "Top left: all nine singular values along the run. The three orange ones — the dead
gauge directions — sit around $10^{-16}$ to $10^{-18}$. The dashed line is the cutoff. They're
eight decades below it. So no leakage. Top right: the logged cond versus a cond computed from
only the 6 live directions. They're the same curve, to $10^{-15}$. So the brake was reading real
conditioning. E3 and Target C stand."

> [!example] Reading the figure
> - **(a)** y-axis = $\sigma_i/\sigma_1$, log scale; x-axis = imaginary time $\tau$, log scale.
>   Grey = the 5 strong live directions. **Blue = $\sigma_6$**, the weakest live one: it dips to
>   ~$10^{-7}$ at $\tau\approx0.2$. **Orange = $\sigma_{7..9}$**, the gauge. Black dashed = `pinv`
>   cutoff at `rcond` $10^{-8}$.
> - **(b)** blue solid = logged cond$_F$, orange dashed = rank-6 cond$_F$. Point at the peak
>   (~$10^7$ near $\tau\approx0.2$) and at how the two lines overlap everywhere.
> - The wiggles in orange near $\tau=0.2$ are just eigensolver roundoff — meaningless at $10^{-17}$.

> [!question] Likely questions
> - *Why are they $10^{-16}$, not the $10^{-10}$ you predicted?* $M$ is a Gram matrix. A true
>   dead direction picks up the fd error **squared**: $(10^{-10})^2=10^{-20}$. What's left is
>   the eigensolver's own roundoff, ~$10^{-16}$. The per-entry error really is ~$10^{-10}$
>   ($1.652\times10^{-10}$), it just doesn't land in the null space. Prediction wrong, in the good
>   direction.
> - *Does `pinv` keep exactly 6?* Yes, on every step of Target C and all four E3 arms.

---

### Slide 8 — P1, chart singularity (Fig. 1 c, d)

**Point:** the angles race, the state strolls. It's a coordinate problem.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_target_c_spectrum.png|560]]

**Say:** "Bottom left: how fast the angles move. At $\tau\approx0.2$ it spikes to 128.
Bottom right: how fast the *state* moves. Nothing happens there — 0.186, smooth. Over the whole
run the angle speed varies by a factor 160, the state speed by under 9. So the state crosses this
point gently; it's the parametrization that has to swing one angle very hard. That's a chart
singularity. And the brake spends 86 of 108 steps getting across it."

> [!example] Reading the figure
> - **(c)** green, $\|\dot\theta\|$, log scale. Sharp spike at $\tau=0.198$.
> - **(d)** pink, $\|d\rho/d\tau\|_F=\sqrt{\dot\theta^\top M\dot\theta}$ = state speed. Smooth hump
>   peaking later near $\tau\approx2$, **no feature at 0.2**. Point at the two x-positions side by
>   side.

> [!info] If you forgot
> - State speed formula: because $M$ is the Gram matrix of the $\partial_k\rho$,
>   $\dot\theta^\top M\dot\theta=\|d\rho/d\tau\|^2$ exactly (pure state). Costs nothing extra.
> - "86 of 108" = steps with logged cond $>10^5$ in the capped arms. Cap-off arms: 4.
> - At the peak, 99.9998 % of $\dot\theta$ points along the one soft direction $v_6$.

> [!question] Likely questions
> - *So is the brake wrong?* Not wrong — Euler's error scales with $\|\dot\theta\|$, so cutting
>   $\Delta\tau$ is the correct response for a step in $\theta$. It's just **expensive** for
>   something the state crosses smoothly.
> - *Is this a new kind of singularity?* A fourth kind next to the Technical Resume's structural /
>   developing / boundary-layer — or arguably "developing", correctly located in the chart, not the
>   physics.

---

### Slide 9 — P1, which gate is soft? (Fig. 2)

**Point:** one gate, th8, causes it. The gate pairs are in the harmless directions.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_null_vectors.png|640]]

**Say:** "This is the breakdown of the four weakest directions at the peak, by gate. The soft
live direction, top row, is almost entirely th8 — an RZ on qubit 0 in layer 2, loading 0.93.
Why: qubit 0 is sitting very close to $|0\rangle$, and an RZ on a qubit near $|0\rangle$ only
adds a phase — it barely moves the state. The pairs I predicted are there, but in the bottom three
rows — the dead gauge directions — which `pinv` already handles."

> [!example] Reading the figure
> - Rows = directions ($v_6$ soft, $v_{7,8,9}$ gauge). Columns = the 9 angles, labelled with gate
>   type and qubit. Colour/number = |loading| (0 to 1).
> - $v_6$: th8 = 0.93, smaller bits on th3 (0.23) and th2 (0.17).
> - $v_7$: th1 & th6 (both RY q0, 0.69 each) — two identical gates turning opposite ways.
> - $v_8$: th2 & th7 (both RY q1, 0.64 / 0.61 on the heatmap; $+0.638$ / $-0.609$ in the P1
>   note). $v_9$: th4 alone (RZ q1, 1.00).

> [!info] If you forgot
> - th8's diagonal: $M_{88}=\tfrac12(1-\langle Z_0\rangle^2)$. With $\langle Z_0\rangle=0.9944$ →
>   $5.568\times10^{-3}$. Small but not tiny.
> - The rest of the softness: $\sigma_6/M_{88}=2.6\times10^{-5}$ — th8's effect is *nearly a
>   combination of the other gates'*. Near-$Z$-eigenstate shrinks the row, near-dependence does the
>   rest.

> [!warning] Don't overclaim
> The prediction "the soft direction is a gate *pair*" **failed**. Say so.

---

### Slide 10 — P1, the fold amplifies fd error

**Point:** fd vs exact metric gives the same run, but small errors get magnified → limited sig figs.

**Say:** "I re-ran Target C with the exact $M$ instead of finite differences. Same steps, same
shape, final errors agree to two to three figures. But before the fold, fd and exact differ by
$10^{-11}$; after it, by $10^{-4}$. The fold is an amplifier. So any Target C or E3 error at
`rcond` $10^{-8}$ is good to 2–3 significant figures — and that's enough, because the E3 rankings
are far further apart."

> [!info] If you forgot
> - **E3's four arms and its full table:** refresher §3, E3 card.
> - The gap numbers: $1.292\times10^{-11}$ for $\tau<0.15$, max $1.101\times10^{-4}$ (at
>   $\tau=2.3223$), $6.068\times10^{-6}$ at the end.
> - At `rcond` $10^{-2}$ fd = exact to all printed digits: that solve throws the stiff direction
>   away, so nothing gets amplified (gap max $1.228\times10^{-9}$).
> - E3 errors to compare: $5.65$, $8.87$, $13.98\times10^{-4}$ (Adam cap-off, capped, Route 3
>   cap-off). Differences are tens of percent, not in the 3rd sig fig.

> [!question] Likely questions
> - *So is it a bug?* No — same trajectory shape, step count, and peak cond to 4 figures. It's
>   sensitivity, not leakage.

---

### Slide 11 — P1, two side findings

**Point:** the loose cutoff hides the problem; Target A quietly loses a direction too.

**Say:** "Two things I didn't plan for. First: at `rcond` $10^{-2}$ the run looks healthy — cond
never above 125. But that's because the solve is throwing away real physical directions too; it
keeps as few as 4 of 6. The true conditioning there reaches $1.9\times10^5$. It looks
well-conditioned because it can't see the stiff direction. Second: even the healthy 5-parameter
Target A drops a real direction on its last 12 steps. The state ends almost exactly on
$|10\rangle$, where the CRZ gate only adds a phase, so its row in $M$ goes to zero."

> [!info] If you forgot
> - "Fewer than 6 on 11 of 17 steps" = how often the loose solve discards a live direction.
> - Target A details: the drop happens for $\tau\in[2.945,\,2.999]$; population of $|10\rangle$ is
>   0.99902; CRZ loading on the dropped direction 0.999.
> - This is the "zero-variance mechanism" from plan §2: a gate whose rotated generator has zero
>   variance in the current state contributes a zero row. Before, we'd only seen it at $\theta=0$.

> [!warning] Don't overclaim
> That Target A's old "late stiffness" ($\sigma_{\min}$ down to $3.464\times10^{-7}$) *is* this
> CRZ row fading is **likely, not shown** — I didn't match it step by step.

---

### Slide 12 — P1 follow-on, a regime-aware cap (Fig. 3)

**Point:** tried a "smarter" brake; it's useless at one cutoff and harmful at the other. Rejected.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_regime_aware_cap.png|640]]

**Say:** "The plan said: if no leakage, try a brake that reads the 6-direction cond instead.
Left: at $10^{-8}$ it's identical — of course, the two conds agree to $10^{-15}$. Right: at
$10^{-2}$ it's actively worse. It brakes on a direction the solve already ignored — the step size
dives to $10^{-4}$ — and it ends with a *bigger* error in six times the steps. Lesson: the brake
should read the conditioning of the solve you actually do, which is what it already does."

> [!example] Reading the figure
> - Both panels: $\Delta\tau$ per step vs $\tau$, log–log. Blue = shipped cap, orange dashed =
>   rank-6 cap.
> - **Left:** curves identical, both dip to ~$2\times10^{-4}$ at the fold ($\tau\approx0.2$).
> - **Right:** identical until $\tau\approx0.55$, then orange dives to $10^{-4}$ and crawls back.

> [!info] If you forgot
> - Numbers: $10^{-2}$ shipped $2.1434\times10^{-3}$ in 17 steps; rank-6 $3.4051\times10^{-3}$ in
>   105 (cap binding on 67).
> - At $10^{-8}$ the tiny error difference ($8.8682$ vs $8.8617\times10^{-4}$) is the fold's
>   sensitivity from slide 10, not an improvement.

---

### Slide 13 — P1 scorecard

**Point:** core held, two specifics failed.

**Say:** "Three of five pieces held — the ones that matter. The two that failed: the gauge values
were much smaller than I said, which is the good direction, and the soft direction was one gate,
not a pair."

---

### Slide 14 — P1 conclusions

**Point:** what P1 decides.

**Say:** "So: no leakage, G6's attribution stands, no re-runs. Quote those errors to 2–3 figures.
It's a chart singularity. Don't ship the regime-aware cap. And the two identities from the plan
are now verified, so they can go into the Technical Resume."

> [!info] If you forgot
> The "two identities" = the state form of $M$ and the sign of $Y$ (slide 4). Moving them into
> [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] is still my open decision — not done yet.

---

### Slide 15 — Part II divider

**Say:** "Part two: can the controller hit a tolerance we ask for?"

---

### Slide 16 — P2, the question

**Point:** our controllers had no `tol`; Adam's stand-in setpoint was circular. Give them a real one.

**Say:** "Neither controller had a tolerance knob. Route 3 just ratchets up and down. Adam used the
reject band as its target — and the audit showed that's circular: it only ever sees steps that
already passed, so it can never decide to shrink. So I built a standard PI controller with an
explicit `tol`, fed by the free error estimate G5 validated. Two flavours: error per step, and
error per unit step. Theory says the first should give slope one-half, the second slope one."

> [!info] If you forgot
> - **What E1 was, and where G4 came from:** refresher §3, E1 card.
> - **G4 in one line:** growth $<1$ on 0 of 190 steps — not a tuning issue, an identity.
> - **The formula:** first factor pushes $\Delta\tau$ up if $e_n<$ `tol`, down if above; second
>   factor damps jumps. Exponents 0.15 / 0.2 from the memo's Eq. (17) at $p=1$, not tuned. Factor
>   clipped to [SHRINK, GROW].
> - **Why ½ and 1:** Euler's local error $e\approx C\Delta\tau^2$. EPS ($e\le$ `tol`) →
>   $\Delta\tau\propto$ `tol`$^{1/2}$ → global error $\propto$ `tol`$^{1/2}$. EPUS
>   ($e/\Delta\tau\le$ `tol`) → $\Delta\tau\propto$ `tol` → global $\propto$ `tol`.
> - **The band stays** as the reject gate, unchanged. Setpoint and reject threshold are now
>   separate things — that's the fix for G4 without retuning anything.
> - Ran on Target A's setup (5 params, `SEED_A`, $\tau_{\max}=3$).

---

### Slide 17 — P2, predictions and fit rule

**Point:** what I committed to before running, including a known flaw.

**Say:** "Six predictions. Number one is the plan's, as written — slope ~1 for error-per-step.
I flagged before running that theory says one-half there, but I left it and scored it as written.
The metric: final-energy error versus a high-accuracy reference. And a fit rule: a run only counts
if the brake isn't dominating it, and I need at least three such runs."

> [!info] If you forgot
> - Full fit rule: run exits on $\tau_{\max}$; cap binds on $<25\%$ of steps; $\Delta\tau$ at a
>   bound on $<25\%$; error $>10\times$ the reference floor. Then least-squares slope, ≥3 runs.
> - **Why an exact-metric reference:** the plan wanted to solve the same fd ODE tightly, but every
>   solver (DOP853, LSODA, Radau) stalls on the fd noise. So the reference uses the exact $M$, $Y$.
>   The fd-vs-exact gap is reported as a floor ($3.42\times10^{-8}$ max) — every measured error is
>   far above it.
> - Grids: EPS `tol` $10^{-6}$…$10^{-2}$ in half decades (9 runs); EPUS $10^{-3.5}$…1 (8 runs).

---

### Slide 18 — P2, no fit possible (Fig. 1 a)

**Point:** final error doesn't follow `tol`. Pre-registered test fails.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_tol_proportionality.png|560]]

**Say:** "Top left: achieved final error against requested tolerance. If the controller were
working the way the test assumes, the points would follow the dotted or dashed lines. Instead the
error just sits around $10^{-5}$ across four decades of `tol`. Only the two smallest-`tol` runs in
each arm pass the fit rule — filled dots. Two isn't three, so no fit. Five of six predictions are
scored failed — but because nothing could be fitted, not because a slope came out wrong."

> [!example] Reading the figure
> - **(a)** x = requested `tol`, y = $|E(\tau_{\rm end})-E_{\rm ref}(\tau_{\rm end})|$, log–log.
>   Blue = EPS, orange = EPUS. Filled = fit-eligible, open = not. Dotted = slope ½ guide, dashed =
>   slope 1 guide. Grey dash-dot at the bottom = fd-vs-exact floor (~$10^{-10}$), far below.

> [!warning] Don't overclaim
> - The two eligible EPS points give a local slope of $+0.47$ — inside prediction 2's window. **Do
>   not** present that as a pass. Two points aren't a fit, and the ≥3 rule was set beforehand.
> - Local slopes overall: $-0.80$ to $+2.62$. Jumpy.

---

### Slide 19 — P2, what sets $\Delta\tau$ (Fig. 1 c, d)

**Point:** the brake takes over early; but the PI law *does* control $e_n$ — G4 is gone.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_tol_proportionality.png|560]]

**Say:** "Bottom left: step size along the run. Small tolerances really do take smaller steps —
the controller is listening. But after $\tau$ around 2.1 every run falls into the same
brake-driven pattern. Bottom right: the share of steps where the brake decides. It crosses 25 %
already at the third grid point in each arm — that's why only two runs were eligible. Good news:
prediction 5 held. The controller shrinks when it should in every run — 527 of 1234 steps at the
smallest tol — and where `tol` actually binds, it holds the error estimate right at the target,
median ratio 0.99. The G4 circularity is fixed."

> [!example] Reading the figure
> - **(c)** $\Delta\tau$ vs $\tau$ (linear x). Blue = EPS, orange = EPUS; solid / dashed / dotted =
>   three tols. Point at the solid lines down at $10^{-3}$, the dotted ones up near $10^{-1}$, and
>   all of them converging into the same zig-zag past $\tau\approx2.5$.
> - **(d)** filled circles = % of steps the cond cap binds; open squares = % with $\Delta\tau$ at
>   min/max (always 0). Grey dotted line = the 25 % fit threshold.

> [!info] If you forgot
> - Above `tol` ≈ $10^{-4}$ (EPS) / $10^{-2}$ (EPUS) the estimator is naturally smaller than `tol`,
>   so `tol` doesn't bind at all: steps are set by the GROW clip and the cap. Sweeping `tol` there
>   barely changes anything.
> - The cap never binds before $\tau=2$. The late regime is Target A's old "developing stiffness",
>   which P1 tied to the CRZ row fading.

---

### Slide 20 — P2, work-precision (Fig. 1 b)

**Point:** PI ends much closer to the reference at similar cost — but don't oversell.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_tol_proportionality.png|560]]

**Say:** "Top right: final error against honest circuit cost. Our existing controllers — Route 3
square, Adam diamond — sit at about $6$–$7\times10^{-4}$ for around 2000 circuits. Several PI runs
end around $10^{-5}$ or better for 2400–2700 circuits. That looks like a 100× win, but it's one
seed and the curve isn't monotone in `tol`, so part of it may be luck in where the contraction
leaves each run. I want a second seed before claiming anything."

> [!example] Reading the figure
> - **(b)** x = circuits (M6 currency, includes retries), y = error at $\tau_{\rm end}$. Green square
>   = Route 3 (Target A), pink diamond = Adam (E1). Blue/orange lines = the PI sweeps; the dip near
>   ~2700 circuits is EPS $3.16\times10^{-4}$.

> [!warning] Don't overclaim
> No "100× better" claim. One seed, non-monotone curve.

---

### Slide 21 — P2, why the end point ignores `tol` (Fig. 2 b)

**Point:** two mechanisms: a tol-blind first step and a flow that forgets.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_first_step_ablation.png|640]]

**Say:** "Right panel: the error along the run. Two things. One: the very first step is always
0.05, before the controller ever acts. That injects the same $5.3\times10^{-4}$ error in every
run — 229 times the tolerance at the smallest `tol`. Two: the error grows until $\tau$ around 1,
then the imaginary-time flow squeezes it down by two to three decades. Imaginary time contracts
toward the ground state, so the end point remembers the last stretch — which is brake-dominated —
and forgets what the controller did."

> [!example] Reading the figure
> - **(b)** y = $|E(\tau)-E_{\rm ref}(\tau)|$ vs $\tau$. Green = fixed first step, blue = first step
>   at minimum. Solid = `tol` $10^{-6}$, dashed = $10^{-4}$. Point at the green jump right at
>   $\tau=0.05$, the hump around $\tau\approx1$, then the steady decay. Vertical dotted line at
>   $\tau=2$ = the window used next slide.
> - The sharp dips are where the error crosses zero (sign change) — on a log plot of $|\cdot|$ they
>   look like spikes down.

> [!warning] Don't overclaim
> Point 2 (contraction) is a **reading of the figure**. I didn't compute a contraction rate.

---

### Slide 22 — P2, control the first step (Fig. 2 a)

**Point:** post hoc: with the first step controlled, the textbook slopes appear.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/figures/p2_first_step_ablation.png|640]]

**Say:** "This is post hoc — I added it after seeing the failure, to test my explanation. Same
sweep, one change: the first step starts at the minimum, so the controller is in charge from the
beginning. And I score the max error before $\tau=2$, where the brake never binds. Now the
slopes are 0.478 for error-per-step — theory one-half — and 0.890 for error-per-unit-step —
theory one. So the free estimator does drive error control. The failure was the metric, plus the
first step, not the estimator."

> [!example] Reading the figure
> - **(a)** x = `tol`, y = max error over $\tau\le2$. Open markers = fixed first step, filled =
>   first step at min. Grey dotted / dashed = slope ½ / 1 guides. Point at how the filled blue line
>   lies on the ½ guide at small `tol`, and filled orange on the 1 guide. Both bend flat at large
>   `tol` — that's where `tol` stops binding, as expected.

> [!info] If you forgot
> - Fixed-step damage: at EPS $10^{-6}$, max error drops from $6.977\times10^{-3}$ to
>   $2.469\times10^{-3}$ just by controlling step 1.
> - Circuits scale more weakly than theory ($-0.385$ vs $-\tfrac12$). My guess (untested): the
>   brake-limited last third adds a fixed overhead.
> - Even with step 1 fixed, the **end-point** error still barely moves ($2.0$–$5.1\times10^{-5}$),
>   so the contraction story stands on its own.

> [!warning] Don't overclaim
> Post hoc: window $\tau\le2$ and "3 smallest tols" were chosen after seeing the data. It does
> **not** rescue predictions 1–6. It motivates P2b, which will pre-register this.

---

### Slide 23 — P2 conclusions

**Point:** what P2 decides.

**Say:** "G4 is closed: a separate `tol` gives a controller that shrinks and holds its target,
no band retuning. The PI controller works, but needs two fixes to be a real tolerance: a first
step that respects `tol`, and the per-unit-step form if we want slope 1. And for QITE,
tolerance proportionality has to be measured before the brake binds — the end point of a
contracting flow erases the evidence."

---

### Slide 24 — The common thread

**Point:** both experiments point at the same culprit and the same lever.

**Say:** "Put P1 and P2 together. In both, the late part of the run belongs to the brake, and
the brake is reacting to the parametrization — th8 in Target C, the CRZ row in Target A. So the
next lever isn't a better tolerance or a cleverer number inside the cap. It's controlling the step
in *state* space: how far the state moves, not how far the angles move."

> [!warning] Don't overclaim
> "The lever is state-space control" is my **interpretation** of P1+P2, not something measured.

---

### Slide 25 — Next (not run)

**Point:** three concrete next experiments, each with a guess.

**Say:** "P1b: detect the fold — angle speed along one direction while state speed stays flat;
both are free — and bound how far the state moves per step instead. My guess: Target C goes from
108 steps to about 22–28, without the error penalty that removing the cap costs Route 3. P2b:
rerun P2 the fixed way, pre-registered, on a second seed and Target C. Guess: slope between 0.8
and 1.1. Then P3: tiny exactly-solvable one-qubit testbeds, where we know the true error and can
score every estimator against truth."

> [!info] If you forgot
> - "Cap-off error penalty": in E3, removing the cap made Route 3 worse ($1.398\times10^{-3}$ in
>   22 steps vs $8.868\times10^{-4}$ in 108).
> - P3 closes the memo's Routes 1 and 2 (embedded pairs, step doubling) and the project's "no true
>   local truncation error" gap ([[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]).
> - Any Adam arm in P2–P5 needs the G2 fix first ($\beta_1=\beta_2$,
>   [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G2_Beta2_Horizon_Dominates_Every_Verdict|G2_Beta2_Horizon_Dominates_Every_Verdict]]).

> [!warning] Don't overclaim
> All three predictions on this slide are **guesses**. None of it has been run.

---

### Slide 26 — Sources

**Say:** "Everything is pinned: 51 values for P1, 93 for P2, re-checked with zero drift, and reruns
are byte-identical. Neither script touches the original controller code, so all earlier results
stay valid."

---

## 6. Tough questions bank

> [!question] "Why not just remove the cap?"
> In E3, removing it helped Adam ($5.653\times10^{-4}$ in 28 steps) but hurt Route 3
> ($1.398\times10^{-3}$ in 22). One seed. The cap's response is *correct* for an Euler step in
> $\theta$; it's just costly. P1b is the principled alternative.

> [!question] "Why trust a finite-difference metric at all?"
> P1 checked it against the exact one: entries agree to ~$10^{-10}$, null space is clean to
> $10^{-16}$, and Target C's trajectory is the same up to the fold-amplified $10^{-4}$ gap.

> [!question] "Isn't the post hoc slope cherry-picked?"
> It's labelled post hoc, and it doesn't change the scorecard (still 5/6 failed). The window
> $\tau\le2$ is justified by the cap never binding there — checked in both sweeps. The fix is to
> pre-register it: that's P2b.

> [!question] "You predicted α≈1 and got nothing — isn't P2 just a failure?"
> As pre-registered, yes, and it's reported that way. But it also closed G4 (prediction 5 held),
> and it diagnosed *why* the metric fails, with a clean post hoc check.

> [!question] "Everything is H$_2$ with one seed. Does any of this generalise?"
> Honest answer: unknown. That's R6 / P5 in the plan (bigger systems). P1's mechanism — a gate
> near a $Z$ eigenstate barely moving the state — is generic, but that's an argument, not a result.

> [!question] "What does 'chart singularity' buy us in practice?"
> It says the expensive steps are spent on coordinates, not physics. That suggests either
> state-space step control or locally freezing/reparametrizing th8 across the fold. Both cheap to
> detect with numbers we already compute.

> [!question] "Why didn't you re-run E3 on the exact metric to be safe?"
> I did re-run Target C on it (slide 10): same step count, same peak cond to 4 figures. The
> difference stays below every gap between E3's arms.

> [!question] "Why didn't you put the PI controller in `adam_qite.py` like the plan said?"
> Kept it in `Proposal_2/` and imported everything, so the original modules — and their audited
> results (Targets A–D, Adam's 147 checks, P1's pins) — stay untouched.

---

## 7. Numbers you must not misquote

| Number | Meaning | Source |
|---|---|---|
| $9.66\times10^{6}$ | peak logged cond($M$), Target C / E3, `rcond` $10^{-8}$ | P1 |
| $2.0\times10^{-16}$ | max $\sigma_7/\sigma_1$ (gauge) — no leakage | P1 |
| $10^{-15}$ | agreement of rank-6 cond with logged cond | P1 |
| 128 / 0.186 | $\|\dot\theta\|$ / state speed at the peak, $\tau=0.198$ | P1 |
| ×160.8 / ×8.83 | range of $\|\dot\theta\|$ / state speed over the run | P1 |
| 86 of 108 | capped-arm steps spent in the fold (cond $>10^5$) | P1 |
| th8, 0.934 | the soft gate (RZ q0, layer 2) and its loading | P1 |
| $1.101\times10^{-4}$ | max fd-vs-exact energy gap after the fold | P1 |
| 2–3 sig figs | precision of Target C / E3 errors at $10^{-8}$ | P1 |
| $8.868\times10^{-4}$ / 108 | Target C, `rcond` $10^{-8}$, capped | P1 / Route 3 |
| $2.143\times10^{-3}$ / 17 | Target C, `rcond` $10^{-2}$ | P1 / Route 3 |
| $3.405\times10^{-3}$ / 105 | rank-6 cap at $10^{-2}$ — why it's rejected | P1 |
| 2 of 9, 2 of 8 | fit-eligible runs, EPS / EPUS — why no fit | P2 |
| 527 of 1234 | shrinking steps at EPS $10^{-6}$ — G4 gone | P2 |
| 0.99 | median $e_n$/tol where `tol` binds | P2 |
| $5.29\times10^{-4}$ | error injected by the fixed first step | P2 |
| 0.478 / 0.890 | post hoc slopes, EPS / EPUS (theory ½ / 1) | P2 |
| $5.67\times10^{-6}$ @ 2700 vs $6.32\times10^{-4}$ @ 2069 | PI vs Route 3 work-precision — one seed | P2 |
| 51 / 93 | pinned values in P1 / P2, 0 drift | P1 / P2 |
