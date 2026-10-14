---
title: "P1 findings — classifying the singularity behind Target C / E3's cond(M)"
tags:
  - report
  - physics
  - quantum-computation
  - adaptative-qite
sources:
  - "[[QITE_Hand_Off/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]]"
  - "[[QITE_Hand_Off/Notes_on_Project/Replication_Beyond_Trotterization|Replication_Beyond_Trotterization]]"
  - "[[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]]"
  - "[[QITE_Hand_Off/Adam_QITE/Report_Adam_QITE|Report_Adam_QITE]]"
  - "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility|G6_Attribution_And_Reproducibility]]"
  - "[[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]]"
code: QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/p1_singularity.py
---

# P1 findings: what Target C / E3's cond(M) is actually measuring

This note reports the run of **P1** from [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]]. The question:
at `rcond` $=10^{-8}$, Target C / E3's logged $\mathrm{cond}(M)$ peaks at $9.66\times10^{6}$ and
the cap $2000/\mathrm{cond}(M)$ pins $\Delta\tau$. Is that real stiffness inside the 6-dim
physical range, or finite-difference leakage into the 3 gauge directions of the rank-6-of-9
ansatz?

> [!success] TL;DR
> - **No gauge leakage.** The decisive part of the prediction held. $\sigma_7/\sigma_1$ never
>   goes above $2.0\times10^{-16}$ on any step. `pinv` keeps exactly 6 values on every step of
>   Target C and of all four E3 arms. The rank-6 cond equals the logged cond to $10^{-15}$. So the
>   cap reads conditioning **inside** the physical range, which **supports G6's
>   re-attribution**. E3 and Target C **don't need re-running**.
> - **Two specific parts of the prediction failed.**
>   - The gauge values sit near $10^{-16}$, not $10^{-10}$–$10^{-12}$. $M$ is a Gram matrix, so
>     the fd error enters the null space quadratically. This one fails in the good direction.
>   - The soft direction is essentially **one gate** (th8, RZ on q0 in layer 2, loading 0.934), not
>     a "near-degenerate pair".
> - **New.** The "stiffness" is a **chart singularity**: stiff in $\theta$, smooth in $\rho$.
>   - At the peak ($\tau=0.198$), $\|\dot\theta\|=128$ and 99.9998% of $\dot\theta$ points along
>     that one soft direction.
>   - Meanwhile the state speed $\|d\rho/d\tau\|_F$ is $0.186$ and completely smooth.
>   - The capped arms spend 86 of their 108 steps in this fold ($\mathrm{cond}>10^5$). The
>     cap-off arms spend 4.
> - The fold amplifies the ~$10^{-10}$ fd error. Final errors at `rcond` $=10^{-8}$ are only good
>   to about 2–3 significant figures, which doesn't change any E3 ranking.
> - **The plan's regime-aware cap: don't adopt it.** At $10^{-8}$ it's the same as the shipped
>   cap. At $10^{-2}$ it's harmful ($3.405\times10^{-3}$ in 105 steps against $2.143\times10^{-3}$
>   in 17).
> - **Plan §2 identities.** Both are verified numerically: the state form of $M$ and the sign of
>   $Y$.

## 1. The prediction, as written before the run

Quoted from the plan (P1, "Prediction"):

> The central difference in `drho` has roundoff $\sim\epsilon_{\rm mach}/h\approx10^{-10}$ per
> entry, well below the $10^{-8}\sigma_{\max}$ cutoff. So I expect $\sigma_{7..9}$ from `drho` to
> sit near $10^{-10}$–$10^{-12}$, already discarded by `rcond` $=10^{-8}$. That makes the reported
> $9.66\times10^{6}$ the conditioning **inside** the physical 6-dimensional range, a genuine
> developing stiffness and not gauge leakage. If that holds:
> - `rank_r_cond` matches the logged cond($M$) along E3;
> - the cap in E3 is reading real stiffness, which *supports* G6's re-attribution;
> - the null-vector attribution at the cond peak names the near-degenerate gate pair.

### Scorecard

| Piece of the prediction | Outcome | Evidence |
|---|---|---|
| $\sigma_{7..9}$ below the $10^{-8}\sigma_1$ cutoff (no leakage) | **held** | max $\sigma_7/\sigma_1=2.021\times10^{-16}$ at `rcond` $10^{-8}$; $1.321\times10^{-16}$ at $10^{-2}$ |
| $\sigma_{7..9}$ sit at $10^{-10}$–$10^{-12}$ specifically | **failed** (favourable) | they sit at ~$10^{-16}$–$10^{-19}$, see Fig. 1(a). The pinned prediction outcome is 0 |
| `rank_r_cond` matches logged cond along E3 | **held** | $\lvert\mathrm{cond}_6/\mathrm{logged}-1\rvert\le1.11\times10^{-15}$ on all four arms |
| The cap reads real stiffness, so G6 is supported | **held** | follows from the rows above |
| Attribution names a near-degenerate gate *pair* | **failed** | $v_6$ is dominated by a single gate (th8, 0.934). The gates *do* pair up, but in the gauge vectors, not the soft one |

Why the magnitude was wrong: I reasoned about the error *per entry* of `drho`, ~$10^{-10}$. But
$M_{jk}=\mathrm{Tr}[\partial_j\rho\,\partial_k\rho]$ is a Gram matrix. A true null vector $v$
of the exact $M$ picks up $v^\top M_{\rm fd}v=\|\sum_kv_k\,\delta(\partial_k\rho)\|^2$, which is
**second order** in the fd error, so ~$10^{-20}$ relative to entries of order 1. What survives is
the eigensolver's own roundoff at ~$10^{-16}$. The per-entry error $\max|M_{\rm fd}-M_{\rm exact}|$
really is ~$10^{-10}$ ($1.652\times10^{-10}$ along the run). It just doesn't land in the null
space.

## 2. Method

Everything is in `p1_singularity.py`. It **imports** `adaptive_qite.py` and `adam_qite.py` and
changes neither of them, so Targets A–D and the Adam 147-check ledger stay valid.

- **Recording without changing behaviour.** A context manager wraps the module-global
  `measure_f` in both modules with a pass-through recorder that logs $(\theta, M, Y)$ per step.
  To prove the recorder is inert, every replay reproduces the logged final error and step count
  bit for bit (pinned: $8.868181396866071\times10^{-4}$ / 108, $2.1434047708084325\times10^{-3}$
  / 17, and so on).
- **Exact metric.** `exact_metric_and_force` builds $M$ and $Y$ in the state form of plan §2 from
  rotated generators $\tilde G_k=P_kG_kP_k^\dagger$. Gate matrices come from `qml.matrix` on the
  same `gate_list`, so wire order and sign conventions can't drift.
- **Diagnostics per step.**
  - `singular_spectrum_log`: all 9 singular values.
  - `null_vector_attribution`: right-singular vectors from $v_6$ onwards, with loadings on the
    parameters.
  - `rank_r_cond`: in both the Frobenius form, which is like-for-like with `frobenius_cond`, and
    the spectral form $\sigma_1/\sigma_6$. $r=6$ is fixed from $2d-2$, not fitted.
- **Pinning.** 51 values and claims are pinned at full precision, with `--check` and no tolerance
  knob. Two full runs were byte-identical, and the ledger reports **0 drifted**.

Targets: `ids` (identities), `c` (Target C replay at both rconds), `null` (attribution at the
peak), `e3` (four arms), `exact` (Target C rerun on the exact $M$, $Y$), `a` (5-param control),
`cap` (the regime-aware cap follow-on).

## 3. Plan §2 identities: verified

| Check | Result | Verdict |
|---|---|---|
| $M=2(\mathrm{Re}\langle\tilde G_j\tilde G_k\rangle-\langle\tilde G_j\rangle\langle\tilde G_k\rangle)$ vs `metric_M(drho)`, 8 random $\theta$ (5- and 9-param) | max diff $2.594\times10^{-10}$ | holds |
| $Y_k=-i\langle[\tilde G_k,H]\rangle$ vs `vector_Y`, same $\theta$ | max diff $3.955\times10^{-10}$ (the opposite sign gives $1.589$) | holds, **sign as derived** |
| CRZ row of the 5-param $M$ at $\theta=0$ | `[0. 0. 0. 0. 0.]` exactly | holds (zero-variance mechanism) |
| Eq. (B9) diagonal on the exact $M$ at $\theta=0.3$ | worst diff $7.216\times10^{-16}$ | holds |

This settles the two caveats that plan §2 flagged as "derived here, not checked". Per plan §7
they can now go into [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §2/§4. I haven't done that here.

## 4. Target C and E3: the cap reads the rank-6 conditioning

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_target_c_spectrum.png]]
*Fig. 1. Target C, `rcond` $=10^{-8}$.*
- *(a) All nine $\sigma_i/\sigma_1$ of $M_{\rm fd}$. The gauge values sit ~8 decades under the
  cutoff.*
- *(b) Logged cond vs rank-6 cond. The curves are indistinguishable.*
- *(c) Parameter speed.*
- *(d) State speed.*

| Quantity | `rcond` $=10^{-8}$ | `rcond` $=10^{-2}$ |
|---|---|---|
| final error / steps (replay = logged) | $8.8682\times10^{-4}$ / 108 | $2.1434\times10^{-3}$ / 17 |
| values kept by `pinv` per step | always 6 | 4, 5 or 6 (fewer than 6 on 11 of 17 steps) |
| max $\sigma_7/\sigma_1$ | $2.021\times10^{-16}$ | $1.321\times10^{-16}$ |
| max $\lvert M_{\rm fd}-M_{\rm exact}\rvert$ | $1.652\times10^{-10}$ | $1.715\times10^{-10}$ |
| rank-6 cond$_F$ / logged cond | $[1-10^{-15},\,1+10^{-15}]$ | $[1,\,2593]$ |
| max logged cond | $9.6637\times10^{6}$ | $125.39$ |
| max spectral $\sigma_1/\sigma_6$ | $6.8367\times10^{6}$ | $1.9281\times10^{5}$ |
| steps with logged cond $>10^5$ | 86 of 108 | 0 of 17 |

The $10^{-2}$ column says something the plan didn't anticipate. At that `rcond`, `pinv` throws
away **physical** directions too, not just the gauge: it keeps 4 of 6 at worst. That's why the
logged cond there is small (125 max) while the true rank-6 cond reaches $1.9\times10^5$. The
"well-conditioned" $10^{-2}$ run is well conditioned because it projects the stiff direction out.

### E3, all four arms (`rcond` $=10^{-8}$)

| arm | error | steps | kept | max $\sigma_7/\sigma_1$ | $\lvert\mathrm{cond}_6/\mathrm{logged}-1\rvert$ | steps cond $>10^5$ |
|---|---|---|---|---|---|---|
| Route 3, cap on | $8.8682\times10^{-4}$ | 108 | 6 | $2.02\times10^{-16}$ | $1.11\times10^{-15}$ | 86 |
| Route 3, cap off | $1.3979\times10^{-3}$ | 22 | 6 | $1.42\times10^{-16}$ | $5.55\times10^{-16}$ | 4 |
| Adam, cap on | $8.8682\times10^{-4}$ | 108 | 6 | $2.02\times10^{-16}$ | $1.11\times10^{-15}$ | 86 |
| Adam, cap off | $5.6530\times10^{-4}$ | 28 | 6 | $1.42\times10^{-16}$ | $5.55\times10^{-16}$ | 4 |

All four reproduce the E3 numbers in `Adam_QITE/README.md`. No arm leaks, so the E3 comparison
is between controllers acting on the same, correctly-ranked solve. That's what
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility|G6_Attribution_And_Reproducibility]] needed.

## 5. The finding I didn't predict: it's a chart singularity

At the cond peak (step 43, $\tau=0.1983$) the spectrum is

$\sigma=[0.9991,\ 0.9828,\ 0.1463,\ 0.09246,\ 0.02142,\ \mathbf{1.461\times10^{-7}},\ 1.4\times10^{-16},\ 3.3\times10^{-17},\ 8.8\times10^{-19}]$

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_null_vectors.png]]
*Fig. 2. Loadings of $v_6$ (softest physical direction) and $v_{7..9}$ (gauge) on the nine
parameters at the peak.*

- **The soft direction $v_6$** is basically th8, RZ on q0 in layer 2 (+0.934), with smaller
  pieces on th3 (RZ q0, +0.233) and th2 (RY q1, −0.172).
  - th8 is followed only by an RZ on q1, so $\tilde G_8=Z_0/2$ and $M_{88}=2\,\mathrm{Var}(Z_0/2)=\tfrac12(1-\langle Z_0\rangle^2)$.
  - At the peak $\langle Z_0\rangle=0.9944$, which gives $M_{88}=5.568\times10^{-3}$, matching
    the numerical diagonal to the printed digits.
  - But $\sigma_6/M_{88}=2.6\times10^{-5}$. So the softness isn't just a weak row: th8's tangent
    vector is nearly a linear combination of the others. Q0 sitting close to a $Z$ eigenstate
    shrinks the row, and the near-dependence does the rest.
- **The gauge vectors** are where the gate pairs are:
  - $v_7$ = th1 − th6 (both RY q0, ∓0.693);
  - $v_8$ ≈ th2 − th7 (both RY q1, +0.638 / −0.609);
  - $v_9$ ≈ th4 alone (RZ q1, −0.998).

  That's the T6-style "coincident generators" picture, but it applies to the structural null
  space, which `rcond` already handles fine.

The important part is the speeds (Fig. 1 c/d):

| $\tau$ | logged cond | $\|\dot\theta\|$ | $\|d\rho/d\tau\|_F$ | fraction of $\dot\theta$ along $v_6$ |
|---|---|---|---|---|
| 0.0000 | $1.156\times10^{4}$ | 0.7989 | 0.1428 | 0.351 |
| 0.1880 | $2.492\times10^{6}$ | 35.23 | 0.1826 | 0.99979 |
| 0.1949 | $7.392\times10^{6}$ | 98.8 | 0.1848 | 0.999992 |
| **0.1983** | $\mathbf{9.664\times10^{6}}$ | **128.4** | **0.1862** | **0.9999978** |
| 0.2016 | $7.476\times10^{6}$ | 100.2 | 0.1877 | 0.999989 |
| 2.8223 | $3.285\times10^{1}$ | 0.6401 | 0.257 | 0.276 |

Over the run, $\|\dot\theta\|$ varies by ×160.8 (max/min) and $\|d\rho/d\tau\|_F$ by only
×8.83, with no feature at the peak. The *state* isn't doing anything violent at $\tau\approx0.2$.
The *parametrization* is: to move the state smoothly, the ansatz has to swing th8 fast, because
th8 barely moves the state there.

So this is a coordinate (chart) singularity: stiff in $\theta$, smooth in $\rho$. It's a fourth
kind to put next to the Technical Resume §4's structural / developing / boundary-layer. Or
arguably it's the "developing" kind correctly identified: it develops in the chart, not in the
physics. The cap responds to it by cutting $\Delta\tau$, which is a correct response for the
explicit-Euler step in $\theta$ (Euler error scales with $\|\dot\theta\|$). But it's an
expensive one: the capped arms spend **86 of 108** steps with cond $>10^5$ getting across a fold
that the state itself crosses smoothly.

### The fold amplifies the finite-difference error

Rerunning Target C with the exact $M$, $Y$ in place of `drho`:

| `rcond` | arm | error | steps | max cond |
|---|---|---|---|---|
| $10^{-8}$ | fd | $8.868181\times10^{-4}$ | 108 | $9.66375\times10^{6}$ |
| $10^{-8}$ | exact | $8.854194\times10^{-4}$ | 108 | $9.66367\times10^{6}$ |
| $10^{-2}$ | fd | $2.143405\times10^{-3}$ | 17 | $125.389$ |
| $10^{-2}$ | exact | $2.143405\times10^{-3}$ | 17 | $125.389$ |

On a common $\tau$ grid, the fd-vs-exact energy gap at $10^{-8}$ is:

- $1.292\times10^{-11}$ before the fold ($\tau<0.15$);
- $1.101\times10^{-4}$ at its maximum ($\tau=2.3223$);
- $6.068\times10^{-6}$ at the end.

At $10^{-2}$, where the stiff direction is projected out, the gap stays at
$1.228\times10^{-9}$ max and $1.986\times10^{-11}$ at the end.

So a $10^{-10}$ perturbation of $M$ becomes a $10^{-4}$ trajectory difference *because of the
fold*. It's a sensitivity amplifier, not a leakage bug. Same trajectory shape, same step count,
same peak cond to 4 figures.

Practical consequence: any final error quoted at `rcond` $=10^{-8}$ on Target C is only good to
about 2–3 significant figures. The regime-aware cap below shows the same thing from the other
side: changing the cap input by ~$10^{-15}$ relative moves the error from $8.8682\times10^{-4}$
to $8.8617\times10^{-4}$. None of the E3 rankings are that close ($5.65$ vs $8.87$ vs
$13.98\times10^{-4}$), so they stand.

## 6. Control: Target A (5-param, full rank, `rcond` $=10^{-6}$)

I expected this to be the boring control. It isn't quite.

| Quantity | Value |
|---|---|
| final error / steps (replay = logged) | $2.2606\times10^{-2}$ / 172 |
| `pinv` keeps | 5, except **4 on 12 of 172 steps**, $\tau\in[2.945,\,2.999]$ |
| softest step | $\tau=2.9987$, $\sigma_5/\sigma_1=6.921\times10^{-7}$ (cutoff $10^{-6}$) |
| $v_5$ | th5, CRZ q0q1, loading +0.999 |
| basis populations $\lvert00\rangle,\lvert01\rangle,\lvert10\rangle,\lvert11\rangle$ | $7.53\times10^{-4},\ 1.74\times10^{-7},\ 0.99902,\ 2.30\times10^{-4}$ |
| $\langle Z_1\rangle$; $M_{\rm CRZ,CRZ}$ | $0.99954$; $8.3662\times10^{-4}$ |
| $\|\dot\theta\|$ max/min; $\|d\rho/d\tau\|$ at the softest step | ×116.7; $7.619\times10^{-2}$ |

The full-rank ansatz **silently loses a physical direction** at the end of the run. The state
converges to within ~$10^{-3}$ of $|10\rangle$. On a computational basis state CRZ only adds a
phase, so its rotated generator has near-zero variance, and `pinv` at $10^{-6}$ drops the CRZ
direction on the last 12 steps.

That's plan §2's zero-variance mechanism developing on its own as the state converges, not only
at $\theta=0$. It's also consistent with Target A's last-quarter "developing stiffness" in the plan
($\sigma_{\min}$ down to $3.464\times10^{-7}$) being this CRZ row fading. I haven't matched it
step by step against `Output.txt`, so treat that as likely, not shown.

## 7. Follow-on: the regime-aware cap. Tested, not adopted

The prediction's core held, so per the plan I ran the follow-on: the cap reads $2000/\mathrm{cond}_6(M)$
instead of the logged cond. `STIFFNESS_CAP_COEFF` is untouched. It's implemented by swapping
`aq.route3_indicators` so its `cond_M` slot carries the rank-6 Frobenius cond.

![[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/figures/p1_regime_aware_cap.png]]
*Fig. 3. $\Delta\tau$ per step, shipped cap vs rank-6 cap, at both rconds.*

| `rcond` | cap on | error | steps | cap binding | min $\Delta\tau$ |
|---|---|---|---|---|---|
| $10^{-8}$ | cond$_F$ (shipped) | $8.8682\times10^{-4}$ | 108 | 82 | $2.07\times10^{-4}$ |
| $10^{-8}$ | rank-6 cond$_F$ | $8.8617\times10^{-4}$ | 108 | 82 | $2.07\times10^{-4}$ |
| $10^{-2}$ | cond$_F$ (shipped) | $2.1434\times10^{-3}$ | 17 | 0 | $0.05$ |
| $10^{-2}$ | rank-6 cond$_F$ | $3.4051\times10^{-3}$ | 105 | 67 | $10^{-4}$ |

- **At $10^{-8}$ it's the same cap.** Of course it is: section 4 showed the two conds agree to
  $10^{-15}$. The error difference is the fold's sensitivity from section 5, not a real
  improvement.
- **At $10^{-2}$ it's actively worse.** It caps on stiffness that the regularized solve already
  discarded. The solve never sees the soft direction, but the rank-6 cond still does, so the
  controller crawls through 105 steps (67 binding) and ends with a *larger* error.

Lesson: the cap should read the conditioning **of the solve actually performed**. That's exactly
what `frobenius_cond` over the kept values already does. The "classify before responding" idea
survives, but the response shouldn't be a different scalar in the same cap.

## 8. What this decides

- **G6 is supported.** In E3 the cap is reading real conditioning of the 6-dim physical solve,
  with no gauge contamination. The re-attribution of E3's difference to the cap stands.
- **E3 and Target C don't need re-running on the exact metric.** Their trajectories are the exact
  metric's up to the fold-amplified ~$10^{-4}$ energy gap, which stays below every gap between
  arms.
- **Precision caveat.** Target C / E3 errors at `rcond` $=10^{-8}$ should be quoted to 2–3
  significant figures, not more.
- **Don't ship the regime-aware cap** in the form the plan proposed.
- **Plan §2's two identities are verified.** They're ready for the Technical Resume, pending my
  call on plan §7's open decision.

### What I'd try next (not run)

The real lever is the chart singularity.

- **Detection is cheap.** $\dot\theta$'s fraction along $v_{r}$ goes to ~1 while the state speed
  $\sqrt{\dot\theta^\top M\dot\theta}$ stays flat. Both numbers are already available per step.
- **The response** would be to control the step in *state* space, bounding
  $\Delta\tau\,\|d\rho/d\tau\|$ instead of $\Delta\tau\,\|\dot\theta\|$ via cond. Another option is
  to locally re-parametrize (freeze th8 across the fold).
- **Prediction** (a guess, labelled as one): Target C at `rcond` $=10^{-8}$ drops from 108 steps
  to about the cap-off 22–28 without the cap-off error penalty that Route 3 pays. Worth a slot as
  a P1b before P2, or folded into P5 where the testbed has boundary layers by construction.

## 9. Assumptions (ours)

1. Every physical choice is imported unchanged from `adaptive_qite.py` / `adam_qite.py`: the H$_2$
   Hamiltonian, the ansätze, `drho` at $h=10^{-6}$, the Route-3 gate, the cap, the seeds, and the
   Adam law. Their own assumption lists apply.
2. The exact metric uses plan §2's state form, with $G_k=P/2$ for RX/RY/RZ and
   $|1\rangle\langle1|\otimes P/2$ for CRX/CRZ. Gate matrices come from `qml.matrix`. This is
   valid for $p=0$ only, and every run here is noiseless.
3. $r=6$ is fixed a priori from $2d-2$, never fitted to a spectrum.
4. `rank_r_cond` is reported in the Frobenius form $\|M\|_F\sqrt{\sum_{i\le r}\sigma_i^{-2}}$,
   which is what `frobenius_cond` returns when `pinv` keeps exactly $r$ values. The spectral form
   is printed alongside.
5. "State speed" means $\|d\rho/d\tau\|_F=\sqrt{\dot\theta^\top M\dot\theta}$. This is exact for a
   pure state, since $M$ is the Gram matrix of the $\partial_k\rho$.
6. The regime-aware cap is a swap of the global `route3_indicators`. The coefficient 2000 is
   unchanged, so there was no retuning.

## 10. Provenance

Every number in this note comes from `Proposal_1/Output.txt`, produced by
`~/micromamba/envs/qml-julia/bin/python p1_singularity.py --target all`. Any of them can be
re-checked with `--check`.

| Numbers | Output section |
|---|---|
| identity diffs $2.594\times10^{-10}$, $3.955\times10^{-10}$, $1.589$, $7.216\times10^{-16}$ | `ids` |
| Target C replay, spectrum stats, kept counts, speed table, ×160.8 / ×8.83 | `c` |
| peak spectrum, $v_6$–$v_9$ loadings, $\langle Z_0\rangle$, $M_{88}$, $\sigma_6/M_{88}$ | `null` |
| four-arm E3 table | `e3` |
| exact-metric errors, energy gaps $1.101\times10^{-4}$ / $6.068\times10^{-6}$ / $1.292\times10^{-11}$ | `exact` |
| Target A kept counts, $\tau$ window, populations, $\langle Z_1\rangle$, ×116.7 | `a` |
| regime-aware cap table | `cap` |
| E3 reference numbers ($5.653\times10^{-4}$/28, $1.398\times10^{-3}$/22) | `Adam_QITE/README.md` (E3), reproduced by the `e3` target |
| Target A's $3.464\times10^{-7}$ | `Attack_Route_03/Output.txt`, via the plan's §8 |
