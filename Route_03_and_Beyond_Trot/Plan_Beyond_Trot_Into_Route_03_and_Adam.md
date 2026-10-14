---
title: "Plan — Beyond Trotterization ideas applied to Route 3, Adam_QITE and the memo's attack routes"
tags:
  - plan
  - physics
  - quantum-computation
  - adaptative-qite
sources:
  - "[[QITE_Hand_Off/Notes_on_Project/Beyond_Trotterization_Variational_Product_Formulas_for_Quantum_Simulation|Beyond_Trotterization_Variational_Product_Formulas_for_Quantum_Simulation]]"
  - "[[QITE_Hand_Off/Notes_on_Project/Replication_Beyond_Trotterization|Replication_Beyond_Trotterization]]"
  - "[[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]]"
  - "[[QITE_Hand_Off/Adam_QITE/Report_Adam_QITE|Report_Adam_QITE]]"
  - "[[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]"
---

# Plan — Beyond Trotterization → Route 3, Adam_QITE and the memo routes

How the ideas in [[QITE_Hand_Off/Notes_on_Project/Beyond_Trotterization_Variational_Product_Formulas_for_Quantum_Simulation|Assi, Vogl, Kumari & LeBlanc (2026)]]
and what [[QITE_Hand_Off/Notes_on_Project/Replication_Beyond_Trotterization|our replication]] of it found could change what's
already built in `Attack_Route_03/` ([[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]]) and
`Adam_QITE/` ([[QITE_Hand_Off/Adam_QITE/Report_Adam_QITE|Report_Adam_QITE]]), and where they sit next to the five attack routes
of the memo, `Research_Papers/QITE_timestep.pdf` (section/equation numbers below are the memo's
unless marked "BT").

> [!warning] This is a plan, not a result
> Nothing in this note has been run, **except P1 and P2** (see the callouts under P1
> and P2, [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/P1_Singularity_Classification_Findings|P1_Singularity_Classification_Findings]] and [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/P2_Tol_Proportionality_Findings|P2_Tol_Proportionality_Findings]]). Every proposal carries a **prediction written before the
> experiment**, and a prediction can fail — per the project's audit history that's the point of
> writing it down first. Numbers quoted as *existing* results are cited to their source (table at
> the end); numbers inside a "Prediction" are guesses and are labelled that way. Two identities in
> §2 were derived here by hand; P1 has since **verified both numerically**.

## 0. The one-paragraph version

BT's EOM $\sum_k g_{jk}\dot c_k = -b_j$ is structurally our $M\dot\theta = Y$, so everything BT
learned about *the metric* transfers, and almost nothing about *the target operator* does
(imaginary-time evolution isn't unitary). The transferable pieces are: (i) a clean way to tell a
**boundary-layer** singularity from a **structural** one (T6), which our $\mathrm{cond}(M)$ cap
currently conflates; (ii) **exact** derivatives via rotated generators instead of finite
differences; (iii) **exactly solvable testbeds** where the true local error is known, which is what
the memo's Routes 1, 2 and 4 need to be judged and what the project has never had; (iv) the
**tolerance-proportionality** test (T2′) as a principled replacement for the band-as-setpoint that
broke Adam (G4); and (v) the **order-crossover** lesson (T3″): an error comparison at a single step
size says nothing, you need the whole curve. The highest-value first move is P1 (classify the
singularity in Target C / E3), then P2 (a real `tol`), then P3 (exact testbeds).

## 1. What BT and the replication actually supply

| Piece | Where | What it is | Why it matters here |
|---|---|---|---|
| EOM $g\dot c = -b$ | BT Eq. (3); T2 convention gate | $g_{jk}=\mathrm{Tr}[\tilde O_j\tilde O_k]/D$, $b_j=\mathrm{Tr}[\tilde O_jH]/D$ | Same linear solve as $M\dot\theta=Y$, so same conditioning questions |
| Rotated-generator trick | `replicate.py::rotated_generators`, `g_and_b` | $\partial_jU_a=i\tilde O_jU_a$, $\tilde O_j=P_jO_jP_j^\dagger$ ⇒ exact $g,b$, no finite differences | Our `drho` uses central differences at $h=10^{-6}$; exact derivatives remove one confounder from every $\mathrm{cond}(M)$ reading |
| cond($g$) boundary layer | T6 | cond($g$) $\sim t^{-2}$ at $t\to0$, rank 2 of 3 at $t=0$, null vector $\propto(1,0,-1)$; along the $N=5$ QIM trajectory cond falls from $3.4\times10^{12}$ (BAB) / $2.1\times10^{16}$ (ABAB) at $t_0$ to 1.36 / 4.85 by $t=2$ | A **third kind** of singularity next to QDP's structural rank-6-of-9 and Route 3's developing stiffness — see [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE\|Technical_Resume_Project_Limiting_QITE]] §4, "classify before responding" |
| Insensitivity to the regularization | T6 | $\mathcal E_F(t{=}2)$ moves $\le5.12\times10^{-10}$ relative over $t_0\in[10^{-8},10^{-3}]$, `rcond` $\in[10^{-14},10^{-4}]$ | The test that *proves* a singularity is benign. We have never run it on Route 3 |
| Error ∝ tolerance | T2′ | RK45 floor tracks `rtol` linearly over nine decades | The standard way to check an adaptive integrator actually *controls* error. Our controllers have no `tol` to sweep |
| Order crossover | T3″ | variational $O(t^3)$ vs Ruth $O(t^5)$; crossover at $t=0.376$ ($N=10$); Ruth 5 orders better at $t=10^{-3}$ | A single-step-size benchmark can flip sign. The memo's Route 1 vs Route 2 comparison must be run as curves |
| Ordering parameter Ξ | BT Eq. (65); T5 | sign of the residual difference $R_{AB}-R_{BA}$ picks the better ordering (7/7 once relabelled) | A residual evaluated *before* simulating predicts ansatz quality — the same logic as Route 3's McLachlan residual |
| Exact solutions | BT Eqs. (32)/(33); T2 | two-level $c(t)$ in closed form, satisfies the EOM to $1.8\times10^{-9}$ | A testbed where the true local truncation error of any integrator is known |
| Fast generators | `replicate.py::Generator`, `qim_generators` | `diag`/`prod`/`dense` representations, $N=10$ affordable | Scale-up past 2-qubit H$_2$ (R6) without new infrastructure |
| Cubic angles | BT Eq. (69); T4 | closed-form corrections, error ratio vs TS 0.63 flat up to $\tau=1$ | Only relevant if a real-time block enters (TITE, P7) |

## 2. What does *not* transfer, and why

1. **State independence (BT's headline) is gone in imaginary time, for two separate reasons.**
   A product of unitaries can't approximate the non-unitary $e^{-H\tau}$, so there is no
   operator-level target to vary against. And the natural "infinite-temperature" force vanishes
   identically. For $\partial_k|\psi\rangle = -i\tilde G_k|\psi\rangle$ one gets
   $Y_k = -\partial_kE = -i\langle[\tilde G_k,H]\rangle_\psi$ (sign by `vector_Y`'s convention,
   G5 derivation), and a state-independent version would be $\mathrm{Tr}\big([\tilde G_k,H]\big)/D = 0$.
   Real-time evolution is driven by a *symmetric* correlator (BT's $\mathrm{Tr}[\tilde O_jH]/D$);
   imaginary time is driven by a *commutator*, whose trace is zero. So QITE has to stay
   state-based. *(Derived here, not checked numerically; the trace identity is exact, but the
   sign and factor conventions need confirming against `vector_Y`.)* **Verified (P1):** the
   sign and factor match `vector_Y` to $3.955\times10^{-10}$ over 8 random $\theta$
   ([[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/P1_Singularity_Classification_Findings|P1_Singularity_Classification_Findings]] §3).
2. **Our $M$ is BT's $g$ with the thermal average swapped for a state average.** Same
   derivation: $M_{jk} = \mathrm{Tr}[\partial_j\rho\,\partial_k\rho] =
   2\big(\mathrm{Re}\langle\tilde G_j\tilde G_k\rangle_\psi - \langle\tilde G_j\rangle_\psi\langle\tilde G_k\rangle_\psi\big)$,
   i.e. twice the symmetrized covariance of the rotated generators in $|\psi\rangle$. BT's
   $g_{jk}=\mathrm{Tr}[\tilde O_j\tilde O_k]/D$ is the same object at infinite temperature (and for
   traceless blocks the centering term is zero). This does two things. It explains why the
   rank-deficiency mechanisms line up: repeated generators give identical rows in both. And it
   adds one mechanism BT can't have: **a generator whose rotated form has zero variance in the
   current state gives a zero row.** Example: the 5-parameter H$_2$ ansatz at $\theta=0$ has an
   exactly zero CRZ row, because the control qubit sits in $|0\rangle$. Its diagonal is consistent
   with the Eq. (B9) check already in `check_eq_b9_self`
   ($M_{kk}=\tfrac12(1-\langle P_k\rangle^2)$ for a pure state). *(Derived here; P1 verifies it
   numerically as its first step.)* **Verified (P1):** the state form matches
   `metric_M(drho)` to $2.594\times10^{-10}$, and the CRZ row at $\theta=0$ is exactly zero. P1
   also caught the same mechanism developing *on its own* at the end of Target A
   ([[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/P1_Singularity_Classification_Findings|P1_Singularity_Classification_Findings]] §3, §6).
3. **Energy monotonicity has no real-time analogue.** Route 3's point 2 (Eq. 13) is QITE-only. BT
   has energy *conservation*, which is necessary but not sufficient, so none of BT's diagnostics
   replace our energy gate.
4. **BT is noiseless and its traces are classical.** The memo's central premise, that *stages are
   expensive* ($O(N^2)$ circuits each, §1), only comes back for BT if $g,b$ come from Hadamard
   tests/LCU (BT §II.C, asserted not run) or Krylov partial traces (BT Fig. 14, a partial negative
   result). So BT gives us no evidence about the noise floor (Eqs. 8–9). Anything noise-related
   here has to be injected by us (P3c, P4).

## 3. Proposals

Ordered by value per unit effort. Format for each: **source → lands in → build → prediction →
cost**.

### P1 — Classify the singularity before the $\mathrm{cond}(M)$ cap responds to it

- **Source.** T6's structural-vs-boundary-layer split; BT's null-vector attribution (T6's
  $(1,0,-1)$); the rotated-generator trick; [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §4's three
  regimes.
- **Lands in.** Route 3 Target C; Adam E3 and its re-attribution to the cap
  ([[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility|G6_Attribution_And_Reproducibility]]); memo Route 3 point 3 and Route 5's "regularization
  coupling".
- **The problem.** The cap $2000/\mathrm{cond}(M)$ reads one scalar and can't tell three different
  situations apart:
  - the 9-parameter null space, which is a *gauge* (rank $6 = 2d-2$ by construction, so no step
    size fixes it);
  - Target A's developing stiffness near $\tau\approx2.5$–3 ($\sigma_{\min}$ down to
    $3.464\times10^{-7}$, spectral cond up to $1.445\times10^{6}$, cap binding on 31 of the last 43
    steps);
  - a possible boundary layer.

  In E3 at `rcond` $=10^{-8}$, cond($M$) reaches $9.66\times10^{6}$ and the cap pins
  $\Delta\tau$ for both controllers identically. Removing it helped Adam ($5.653\times10^{-4}$ in
  28 steps) and hurt Route 3 ($1.398\times10^{-3}$ in 22 steps). Nobody has looked at *which
  directions* carry the small singular values.
- **Build** (new functions only; `route3_indicators`'s return signature stays put because
  `adam_qite.py` imports it):
  1. `exact_metric_and_force(theta, n_params, H)`: the state-form rotated-generator $M$ and $Y$ for
     RX/RY/RZ/CRX/CRZ, using $M=2\,\mathrm{Cov}_\psi(\tilde G_j,\tilde G_k)$ from §2. Asserted
     against `metric_M(drho(...))`/`vector_Y` at a few random $\theta$ in the style of
     `check_eq_b9_self`. This is also where §2's identity gets checked.
  2. `singular_spectrum_log(M)`: all singular values per step, not just $\sigma_{\min}$.
  3. `null_vector_attribution(M, r)`: loadings of the right-singular vectors beyond rank $r$ on the
     parameters, T6-style (which gates coincide).
  4. `rank_r_cond(M, r)`: condition number over the top-$r$ singular values, $r=6$ for the
     9-parameter ansatz (known a priori from $2d-2$, **not** fitted).
- **Prediction** (the reasoning is written out so it can be wrong in a specific way). The central
  difference in `drho` has roundoff $\sim\epsilon_{\rm mach}/h\approx10^{-10}$ per entry, well below
  the $10^{-8}\sigma_{\max}$ cutoff. So I expect $\sigma_{7..9}$ from `drho` to sit near
  $10^{-10}$–$10^{-12}$, already discarded by `rcond` $=10^{-8}$. That makes the reported
  $9.66\times10^{6}$ the conditioning **inside** the physical 6-dimensional range, a genuine
  developing stiffness and not gauge leakage. If that holds:
  - `rank_r_cond` matches the logged cond($M$) along E3;
  - the cap in E3 is reading real stiffness, which *supports* G6's re-attribution;
  - the null-vector attribution at the cond peak names the near-degenerate gate pair.

  **If it fails** ($\sigma_7$ above the cutoff anywhere), the E3 trajectory at `rcond` $=10^{-8}$ is
  contaminated by finite-difference leakage into the gauge directions. Then E3 and Target C both
  need re-running on the exact metric before any claim about them stands. Either outcome decides
  something that is currently open.
- **Follow-on, only if the prediction holds.** A **regime-aware cap**: apply
  $\mathrm{COEFF}/\mathrm{cond}_r(M)$ with $r$ the structural rank, and log the regime
  (structural / developing / boundary) per step. That's the concrete version of "classify before
  responding". It doesn't touch `STIFFNESS_CAP_COEFF`, so there's no retuning.
- **Cost.** Small. Needs no new runs beyond replaying Targets A/C with extra logging. Minutes of
  compute.

> [!success] Run: no gauge leakage, and the stiffness is a chart singularity
> Full write-up: [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_1/P1_Singularity_Classification_Findings|P1_Singularity_Classification_Findings]]. Code and output are in `Proposal_1/`:
> 51 pinned values, 0 drift, byte-identical reruns.
> - **The core held.** Max $\sigma_7/\sigma_1=2.0\times10^{-16}$, and `pinv` keeps exactly 6 on
>   every step of Target C and all four E3 arms. The rank-6 cond equals the logged cond to
>   $10^{-15}$. So the cap reads real conditioning inside the physical range: **G6 is supported,
>   and E3/C don't need re-running.**
> - **Two specifics failed.**
>   - The gauge values sit at ~$10^{-16}$, not $10^{-10}$–$10^{-12}$: $M$ is a Gram matrix, so
>     the fd error enters its null space quadratically.
>   - The soft direction is one gate (th8, RZ q0 layer 2, loading 0.934), not a pair. The pairs
>     are in the gauge vectors (th1/th6, th2/th7).
> - **New: it's a chart singularity.** At the peak ($\tau=0.198$), $\|\dot\theta\|=128$ with
>   99.9998% along that soft direction. The state speed $\|d\rho/d\tau\|$ stays smooth: ×8.83
>   over the whole run, against ×160.8 for $\|\dot\theta\|$. The capped arms spend 86 of 108
>   steps crossing it.
> - The fold amplifies the $10^{-10}$ fd error to a $1.1\times10^{-4}$ energy gap against the
>   exact metric. Errors at `rcond` $=10^{-8}$ are good to 2–3 significant figures, and the E3
>   rankings stand.
> - **Target A control.** `pinv` silently drops the physical CRZ direction on the last 12 of 172
>   steps. That's §2's zero-variance mechanism as the state converges to ~$|10\rangle$.
> - **Regime-aware cap: rejected.** At $10^{-8}$ it's identical to the shipped cap. At $10^{-2}$
>   it's worse ($3.405\times10^{-3}$ in 105 steps against $2.143\times10^{-3}$ in 17), because it
>   caps on stiffness the regularized solve already discarded.
> - **Next idea** (not run): detect "θ̇ along the soft direction while the state speed is flat",
>   and step-control in state space instead.

### P2 — A real `tol`, and tolerance proportionality as the acceptance test for the controller

- **Source.** T2′: for a working adaptive integrator the achieved error tracks the requested
  tolerance (RK45: ~linear over nine decades of `rtol`).
- **Lands in.** Adam's finding 2 / [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G4_Controller_Cannot_Shrink|G4_Controller_Cannot_Shrink]] (the reject band is too loose
  to be a setpoint; $\hat e>1$ on 0 of 190 steps, max 0.945); [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] R3 (the PI
  controller, Eq. 17); the memo's "tol ≥ noise floor" rule (§2, Algorithm 1 line 6).
- **Build.**
  - A `PIStepController` next to `AdamStepController` in `adam_qite.py`, implementing Eq. (17)
    with $k_I=0.3/(p+1)$, $k_P=0.4/(p+1)$ at $p=1$.
  - The input is the **unnormalized** $e_n=|\Delta E_{\rm actual}-\Delta E_{\rm pred}|$, the
    estimator G5 validated against step doubling (log-Pearson 0.9902, median ratio 2.060 over the
    first 40 accepted steps). It is divided by 2 per G5's derivation so it estimates the same
    quantity step doubling does.
  - `tol` is a new, explicit argument. `DE_BAND_*` stays as the *reject* gate, unchanged, which
    separates setpoint from band without retuning the band (the fudge G4 rules out).
  - Sweep `tol` over four decades on Target A's $\theta_0$ (`SEED_A`) and record final energy error
    against the variational floor ($E_{\rm floor}=-1.124600$), steps, and M6 circuits.
- **Prediction.**
  - Global error $\propto\mathrm{tol}^{\alpha}$ with $\alpha\approx1$ (explicit Euler with local
    error control gives global error roughly proportional to tol) over the decades where neither
    the variational floor nor the cond cap binds.
  - The `tol` at which the slope breaks marks where the controller stops controlling. I expect the
    break to coincide with the cap taking over (Target A's last-quarter regime).
  - If there's **no** decade with slope ≈1, then G5's estimator doesn't actually drive error
    control even though it correlates with the true LTE. That would be a negative result about
    the G5 premise.
- **Precondition.** None for the PI arm. If an Adam arm is added to the same sweep, fix
  [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G2_Beta2_Horizon_Dominates_Every_Verdict|G2]] first ($\beta_1=\beta_2$), or the sweep measures
  the horizon artifact again.
- **Cost.** Small–medium. About 8 runs × ~20 s each.

> [!failure] Run: the controller controls, but the end-point error doesn't track `tol`
> Full write-up: [[QITE_Hand_Off/Route_03_and_Beyond_Trot/Proposal_2/P2_Tol_Proportionality_Findings|P2_Tol_Proportionality_Findings]]. Code and output are in `Proposal_2/`. Every
> value is pinned, `--check` exits 0, and reruns are byte-identical.
> - **The α above was wrong for the arm it names.** For error-per-step control ($e_n\le$ tol) of a
>   $p=1$ method, global error goes as $\mathrm{tol}^{1/2}$. $\alpha\approx1$ needs error per unit
>   step. I flagged this before running, left the text as written, and ran both arms.
> - **As pre-registered, it fails: no fit is possible.** Against an exact-metric reference,
>   err($\tau_{\rm end}$) sits near $10^{-5}$ across four decades of `tol`, with local slopes from
>   $-0.80$ to $+2.62$. The cap binds on $\ge25\%$ of steps from the third grid point, so only 2
>   runs per arm are fit-eligible. The break prediction fails too, with no $\alpha$ to break from.
> - **G4's circularity is gone.** Growth $<1$ happens in every run (527 of 1234 steps at $10^{-6}$),
>   and median $e_n$/tol is 0.99 where `tol` is below the estimator's scale.
> - **Post hoc, labelled as such:**
>   - The fixed first step ($\Delta\tau_0=0.05$) injects $5.29\times10^{-4}$ regardless of `tol`.
>   - The flow then contracts errors after $\tau\approx1$, so $\tau_{\rm end}$ forgets the
>     controller.
>   - With the first step at `DELTA_TAU_MIN`, max error over the cap-free $\tau\le2$ scales with
>     slope **0.478** (per step) and **0.890** (per unit step). So G5 does drive error control, and
>     the "negative result" branch above is triggered by the metric, not the estimator.
> - **Next** (not run): a P2b that pre-registers the unit-step form, a tol-aware first step, and the
>   pre-cap metric, on a second seed and Target C.

### P3 — Exactly solvable testbeds for the memo's estimators (Routes 1, 2, 4)

- **Source.** BT Eqs. (32)/(33) (exact two-level $c(t)$, confirmed in T2), plus the T3″
  lesson that estimator comparisons have to be curves across step size.
- **Lands in.** [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] R1 (step doubling, Eq. 12), R2 (Heun–Euler, BS3 with FSAL,
  §3.1), R5 (Route 4 gating, Eqs. 14–16), and the memo's noise-floor algebra (Eqs. 8, 9, 11), which
  R_Future_Routes says "cannot be evaluated at all with what is built".
- **(a) Real time, BT two-level model.** $H=h_x\sigma_x+h_z\sigma_z$ at BT's Fig. 2 parameters
  ($h_x=5$, $h_z=2$), ABA ansatz. Integrate $\dot c=-g^{+}b$ with explicit Euler, step doubling,
  Heun–Euler 2(1) and BS3(2)+FSAL at fixed $\Delta t$ on a log grid. The exact $c(t)$ gives the
  **true** local error per step, so each estimator's $e_n$ is scored against truth rather than
  against another estimator. Reuse `two_level_generators`, `g_and_b`, `eqs_32_33_c`. Start at
  $t_0>0$ as the replication does (T6).
- **(b) Imaginary time, our construction on the same qubit.** $e^{-H\tau}|0\rangle$ has a closed
  form, so $E(\tau)$ is exact.
  - A two-parameter RZ·RX ansatz covers the Bloch sphere. The McLachlan flow is then exact ITE
    away from the pole $\theta_{\rm RX}=0$, where the RZ row vanishes: the §2 zero-variance
    mechanism, a toy boundary layer.
  - A three-parameter RX·RZ·RX ansatz is structurally rank ≤ 2 of 3: a toy of the 9-parameter
    gauge problem.
  - This re-tests G5's ratio of 2 against the **true** LTE instead of against step doubling.
  - It is also the first place P1's regime classifier can be checked on a case where the answer
    is known in advance (full-rank / boundary / structural, one each).
  - *(That RZ·RX gives exact ITE is derived here from the manifold dimension, not checked; the
    build asserts it.)*
- **(c) Seeded noise injection.**
  - Perturb $(g,b)$ in (a), or $(M,Y)$ in (b), by $\varepsilon\,\|\cdot\|\,\xi$ with $\xi$ from
    `default_rng(SEED)`. This is a **stand-in** for shot noise, flagged as such in ASSUMPTIONS.
    It is not R5's physical shot model, but it is the first *stochastic* estimator the project
    would have.
  - It makes memo Eq. (8) directly testable: fit $C_{\rm tr}$ from the noiseless run, predict
    $\Delta\tau_{\rm floor}$ from Eq. (9) and $\Delta\tau^\star$ from Eq. (11), then measure both.
  - Route 4's hypothesis-test gate (Eqs. 14–16) can be exercised here because $\sigma_e$ is known
    by construction.
  - It also gives Adam's E2 the stochastic noise its premise always needed
    ([[QITE_Hand_Off/Audits/Audit_Adam_QITE/G3_Noise_Experiment_Cannot_Resolve_Noise|G3]]).
- **Predictions.**
  - (a) In the noiseless case, each estimator's $e_n$ tracks the true LTE with slope $p+1$ on
    log–log axes. Step doubling's ratio to truth is ~1 by construction. Heun beats Euler at equal
    $f$-evaluations only above some $\Delta t$, and there **is** a crossover (T3″ logic).
  - (c) Below $\Delta\tau_{\rm floor}$, $e_n$ flattens from $\propto\Delta\tau^{p+1}$ to
    $\propto\Delta\tau$ (Eq. 8's two terms). The measured floor agrees with Eq. (9) to within a
    factor ~2 if the stages' errors are really independent. It comes out *lower* than Eq. (9) for
    step doubling if the memo's §3.2 claim holds (shared-initial-point correlations cancel), and
    Heun's advantage over Euler disappears below the floor (§3.1: "high order buys nothing below
    the floor").

  Each of these can fail independently. The §3.2 correlation claim is the one I'd bet least on.
- **Cost.** Small. Everything is $D=2$, runs in seconds, and needs no new dependencies.

### P4 — Noise meets the boundary layer: is there an optimal start time $t_0^\star$?

- **Source.** T6 (cond $\sim t^{-2}$, result insensitive to $t_0$ **when noiseless**) combined
  with memo Eq. (6) (relative error in $\dot\theta$ $\lesssim\mathrm{cond}\times$ noise).
- **Lands in.** Route 5 regularization coupling (memo §3.5; R4). It is the noisy version of the
  Technical Resume's prescription "start at $t_0>0$ from the Trotter–Suzuki value rather than
  regularize harder".
- **Build.** On the $N=5$ QIM (T6's setting), inject P3(c)'s seeded noise into $(g,b)$ at
  $\varepsilon\in\{10^{-6},10^{-4},10^{-2}\}$. Sweep $t_0$ over $[10^{-8},10^{-1}]$ and record
  $\mathcal E_F(t{=}2)$. Reuse `qim_generators`, `integrate_c` (swap in a noisy `cdot`) and `E_F`,
  with the BAB and ABAB orderings.
- **Prediction** (schematic, and the exponent is the thing being tested). The noise contribution
  near the start integrates to $\sim\int_{t_0}\varepsilon\,t^{-2}\,dt\sim\varepsilon/t_0$. The
  initial-condition error from the Eq. (10) linear start grows like $t_0^{k}$, with $k$ its order.
  Balancing the two gives $t_0^\star\sim\varepsilon^{1/(k+1)}$, and T6's insensitivity should
  **break** into a U-shaped curve with its minimum moving as $\varepsilon^{1/(k+1)}$.
- **The competing outcome.** The amplified directions near $t=0$ are close to the null vector
  $(1,0,-1)$, which is a near-gauge that barely moves $U_a$. So insensitivity may *survive* the
  noise. That would mean BT's boundary layer is benign under noise too, which is good news for the
  "start at $t_0$" prescription and bad news for the scaling argument above. Either result goes
  into the Technical Resume §4 regime discussion.
- **Cost.** Medium. Each $N=5$ integration takes seconds; the grid is about 3 ε × 15 $t_0$ × 2
  orderings.

### P5 — QITE with a BT-product (Hamiltonian-variational) ansatz on the QIM, $N=4$–10

- **Source.** BT's ansatz Eq. (4), the fast `Generator` machinery, Ξ (Eq. 65), and the
  boundary-layer handling.
- **Lands in.** [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] R6 (scale past H$_2$); Adam E3's $n=1$ problem; G5's open
  "design a probe that can detect the $\partial_\theta M^+$ blind spot"
  ([[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot|G5_Error_Estimator_Validated_And_Its_Blind_Spot]]; [[QITE_Hand_Off/Adam_QITE/Report_Adam_QITE|Report_Adam_QITE]] §8,
  item 6).
- **Build.**
  - Ansatz $|\psi(\theta)\rangle=\prod_j e^{i\theta_jO_j}|\psi_0\rangle$ with $O_j$ alternating QIM
    blocks A/B, in BAB and ABAB orderings, depth 2–4 layers.
  - Exact $M$ and $Y$ in state form (§2), which costs $O(mD)$ vectors instead of BT's $O(D^2)$
    operators. Losing state independence buys back an order of cost.
  - $H$, $E_{\rm exact}$ and the ITE ground truth come straight from `qim_generators`.
  - Run Route 3 (`adaptive_qite_run` generalized to take a `measure_f` callable) and the P2 PI
    controller.
- **Three things only this testbed gives:**
  1. **State-form T6.** Any product $|\psi_0\rangle$ is an eigenstate of one block, so at
     $\theta=0$ that row of $M$ is zero, and repeated generators give coincident rows. The start is
     singular by construction, exactly as in BT. Handle it T6's way: start at small $\tau_0$ and
     **prove** insensitivity by sweeping $\tau_0$ and `rcond`, rather than choosing a regularizer.
  2. **A Ξ analogue for QITE.** Compare Route 3's McLachlan residual $r_0=\|M\dot\theta-Y\|$
     between orderings at $\tau_0$, before evolving. **Prediction:** the lower-$r_0$ ordering
     reaches the lower final energy on at least 6 of 7 pre-registered $(N,h_x,h_z,J)$ points (T5's
     test shape). If not, the residual is a poor ansatz-selection signal, and Route 3's use of it
     as a *caution* signal (not a gate) looks better justified.
  3. **The G5 blind-spot probe.** In a boundary layer cond $\propto\tau^{-2}$ with an
     ordering-dependent prefactor (T6: $3.4\times10^{12}$ vs $2.1\times10^{16}$ in real time; the
     imaginary-time prefactors are unknown). At a matched cond **level**, the two orderings sit at
     different $\tau$, so at matched $\Delta\tau$ their $|\Delta\log\mathrm{cond}|=2\Delta\tau/\tau$
     differs. That is precisely the "matched $\Delta\tau$, matched level, different change"
     sampling G5 said the H$_2$ trajectory couldn't provide. **Prediction:** G5's
     ratio-departure-from-2 correlates with $|\Delta\log\mathrm{cond}|$ at matched level. The H$_2$
     attempt refuted the unmatched version (partial correlations $-0.5437$ level vs $+0.1523$
     change), so this is a second, better-designed shot at the same mechanism, and it may fail
     again.
- **Cost.** Largest item. New `measure_f` path and more targets, but no new dependencies, and $N=10$
  state vectors are cheap.

### P6 — Route 5 hooks (speculative, recorded so they aren't lost)

- **Analytic Jacobians.** Where BT's traces are closed-form (QIM/XXZ, BT App. C), so is
  $\partial_\theta f$. That makes the memo's linearly-implicit Rosenbrock step (§3.5) cheap on
  those models. It only matters once P5 exists.
- **Frozen-layer macro-steps.** BT's stroboscopic construction (solve once for $U_a(\tau)$, then
  repeat it) is the real-time relative of the memo's "freeze $M,Y$ over a macro-step and integrate
  the frozen linear system exactly". The analogy is loose, and its only use is as a guide to which
  macro-step length to test.
- **Problem-shaped Tikhonov.** $(M+\lambda\,g_{\infty})^{-1}Y$ with $g_\infty$ BT's
  infinite-temperature metric in place of $\lambda I$ (§2, point 2). It regularizes along the
  generators' own geometry instead of isotropically. Pure speculation; no prediction yet.

### P7 — TITE composition with a BT propagator (lowest priority, R8)

[[QITE_Hand_Off/Notes_on_Project/Trading_Imaginary_Time_for_Randomness_in_Ground_State_Preparation|TITE]] appends randomly timed
real-time evolution after ITE and needs a propagator $e^{-iHt}$. BT supplies one that is cheaper
at fixed accuracy. **But only above the T3″ crossover**: the four-exponential ansatz is $O(t^3)$,
so for short real-time layers Ruth or plain TS is the better propagator. Any TITE sampling
distribution that puts weight on short $t$ (TITE's Gaussian and Cauchy-like fallbacks do) needs a
per-layer $\tau$ kept above the crossover, with the step count adjusted to match. This doesn't
interact with the step controller at all, which is why it's last.

## 4. Cross-map: memo routes × BT

| Memo route | What BT / the replication supplies | What changes in Route 3 / Adam | Item |
|---|---|---|---|
| **1** Embedded pairs (§3.1) | BT's own integrator *is* a Route 1 method (RK45 = Dormand–Prince 5(4)), and T2′ is its tol-proportionality test; exact two-level solution | Heun–Euler / BS3 scored against true LTE, as curves across $\Delta\tau$ | P3a, P3c |
| **2** Step doubling (§3.2, Eq. 12) | Exact testbeds; the order-crossover warning | G5's ratio-of-2 re-tested against truth, not against another estimator; §3.2's correlation-cancellation claim tested | P3a–c |
| **3** Physics indicators (§3.3) | T6 null-vector attribution; Ξ as a residual-based selector; exact derivatives | Regime-aware cond cap; residual as ansatz-ordering selector; exact $M$/$Y$ replaces `drho` for diagnostics | P1, P5 |
| **4** Statistical gating (§3.4, Eqs. 14–17) | Nothing directly (BT is noiseless) | Seeded-noise stand-in makes $\sigma_e$ known, so Eqs. 14–16 can be exercised; PI law Eq. 17 with a real `tol` | P2, P3c |
| **5** Stiffness-aware (§3.5) | Boundary-layer vs structural classification; insensitivity sweep; analytic traces | $t_0^\star$ under noise; problem-shaped regularization; Rosenbrock on analytic models | P4, P6 |

## 5. What would change, folder by folder

The rules carried over from the replication and the audits:
- **Import, don't copy.** Adam's `sys.path` pattern. `replicate.py` is importable: its `main()`
  is guarded and module level only builds constants.
- **Freeze every constant before the first accuracy run**, and list it in `ASSUMPTIONS`.
- **A `--check` gate with pinned full-precision numbers**, like `T6_PINNED`.
- **Findings pinned in a dated verify script with no tolerance flag**, like
  `verify_audit.py`.

| Folder | Change | Notes |
|---|---|---|
| `Attack_Route_03/adaptive_qite.py` | **Add** `exact_metric_and_force`, `singular_spectrum_log`, `null_vector_attribution`, `rank_r_cond` (P1); make `adaptive_qite_run` accept an optional `measure_f` callable (P5) | Default behaviour byte-identical; Targets A–D re-run to confirm (the M5 reproducibility standard) |
| `Adam_QITE/adam_qite.py` | **Add** `PIStepController` (P2); G2 fix as a separate, prior commit if Adam arms join any new sweep | `AdamStepController` untouched so the audit's 147-comparison ledger still applies |
| `Replications/Beyond_Trot/replicate.py` | **None.** It's a replication and its `--check` gate stays meaningful only if it doesn't move | Imported for `Generator`, `qim_generators`, `two_level_generators`, `g_and_b`, `integrate_c`, `eqs_32_33_c`, `E_F` |
| `Route_03_and_Beyond_Trot/` (this folder) | New `bt_route3.py`, one function per target (`p1_singularity_classification`, `p2_tol_proportionality`, `p3a_two_level_estimators`, `p3b_itebed_qubit`, `p3c_noise_floor`, `p4_t0_star`, `p5_qim_hva_qite`), `--target`, `--check`, `SEED`, `ASSUMPTIONS`; a README; figures | A plan note like this one becomes the entry point, then a report once runs exist |

## 6. Suggested order

```
P1  classify the singularity (Target C / E3 replay, exact metric)   ── DONE: no re-run needed (Proposal_1/)
 │
P2  real tol + PI law, tol-proportionality sweep                    ── DONE: G4 closed; T2' at tau_end fails, pre-cap slopes 0.48/0.89 post hoc (Proposal_2/)
 │
P3a,b  exact testbeds (real-time BT qubit, imaginary-time qubit)    ── R1/R2 scored against truth
 │
P3c + P4  seeded noise: memo Eqs. 8/9/11 measured; t0* under noise  ── first evaluation of the noise floor
 │
P5  BT-product ansatz QITE on the QIM, N = 4..10                    ── R6 scale-up, Ξ analogue, G5 probe, E3 at n > 1
 │
P6  Route 5 hooks, only where P5 shows stiffness worth handling

P7 last, or never — orthogonal to step control.
Precondition for any Adam arm in P2–P5: G2 fix (β1 = β2) first.
```

What each step closes: P1 → the E3 attribution (G6) and Target C's cond reading; P2 → G4, R3; P3 →
R1, R2 and the "no local truncation error estimate" gap; P3c → G3 (E2's missing stochastic noise),
part of R5; P4 → R4; P5 → R6, Adam report §8 item 6, E3's $n=1$.

## 7. Open decisions (for me, before building)

- Is a **seeded Gaussian perturbation of $(M,Y)$** acceptable as the stochastic stand-in for P3c,
  or should R5's proper shot model come first? The stand-in is much cheaper and makes Eq. 9
  testable now, but it can't say anything about shot-noise *structure* (covariances between
  entries).
- **P5's system**: QIM (BT's machinery reused directly) vs LiH (the QDP paper's projected failure
  case). QIM is cheaper and reuses the replication. LiH is closer to the project's stated question.
  Plan assumes QIM first.
- Whether §2's two identities deserve a short standalone note once P1 verifies them. They'd belong
  in [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §2/§4 only **after** verification.

## 8. Provenance of every quoted number

| Number | Meaning | Source |
|---|---|---|
| $3.4\times10^{12}$ / $2.1\times10^{16}$ → 1.36 / 4.85 | cond($g$) at $t_0$ → $t=2$, BAB / ABAB, $N=5$ | [[QITE_Hand_Off/Notes_on_Project/Replication_Beyond_Trotterization\|Replication_Beyond_Trotterization]] T6 |
| $5.12\times10^{-10}$ | max relative shift of $\mathcal E_F(t{=}2)$ over $t_0$, `rcond` sweeps | same, T6 |
| rank 2 of 3, null $(1,0,-1)$, cond $\sim c_1^{-2}$ | two-level metric at $t=0$ | same, T6; BT Eq. (35) |
| floor ∝ `rtol` ($9.9\times10^{-4}$ … $3.3\times10^{-13}$) | T2′ tolerance sweep | same, Discrepancy 5 |
| $t=0.376$; slopes 3.00 / 5.00 | variational vs Ruth crossover, $N=10$ | same, Discrepancy 4 (T3″) |
| 7/7 | Ξ predicts better ordering (relabelled) | same, T5 |
| 0.63 | Eq. (69) / TS error ratio up to $\tau=1$ | same, closing section (T4) |
| $1.8\times10^{-9}$ | Eqs. (32)/(33) satisfy EOM (median) | same, T2 |
| 31 of 43; $\sigma_{\min}\ge3.464\times10^{-7}$; spectral cond $\le1.445\times10^{6}$ | Target A last-quarter cap binding; spectra | `Attack_Route_03/Output.txt` |
| $E_{\rm floor}=-1.124600$ | 5-param variational floor | [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation\|Adaptive_QITE_Attack_Route_3_Implementation]] (M6) |
| $8.868\times10^{-4}$, 108 steps; cond $[32.9,\,9.66\times10^6]$ | Target C / E3, `rcond`$=10^{-8}$, cap on | `Output.txt`; `Adam_QITE/README.md` |
| $2.143\times10^{-3}$, 17 steps; cond $[19.1,125]$ | Target C, `rcond`$=10^{-2}$ | `Output.txt` |
| $5.653\times10^{-4}$ / 28; $1.398\times10^{-3}$ / 22 | E3 cap-off, Adam / Route 3 | `Adam_QITE/README.md` (E3 in full) |
| 0.9902; 2.060; first 40 steps | G5 log-Pearson; median ratio | [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot\|G5_Error_Estimator_Validated_And_Its_Blind_Spot]] |
| $-0.5437$ / $+0.1523$ | G5 rank partial correlations, level / change | same |
| 0 of 190; max 0.945 | $\hat e>1$ count; max $\hat e$ | [[QITE_Hand_Off/Adam_QITE/Report_Adam_QITE\|Report_Adam_QITE]] §3 |
| $k_I\approx0.3/(p+1)$, $k_P\approx0.4/(p+1)$ | PI gains | memo Eq. (17) |
| $h_x=5$, $h_z=2$ | BT two-level parameters | BT Fig. 2; T2 |
