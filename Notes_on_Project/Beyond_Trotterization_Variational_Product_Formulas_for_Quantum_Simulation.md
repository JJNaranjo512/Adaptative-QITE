---
fileClass: Paper_Reading
title: "Beyond Trotterization: Variational Product Formulas for Quantum Simulation"
authors: Ibsal Assi, Michael Vogl, Meenu Kumari, J. P. F. LeBlanc
year: 2026
venue: journal
arxiv_id: 2511.15124
doi: 10.1103/941r-ql3g
status: reading
tags: [paper, physics, quantum-computation, paper-notes]
---

# Beyond Trotterization: Variational Product Formulas for Quantum Simulation

> [!info] Citation
> Ibsal Assi, Michael Vogl, Meenu Kumari, J. P. F. LeBlanc, "*Beyond trotterization: Variational
> product formulas for quantum simulation*," Phys. Rev. B **113**, 214314 (2026) — Editors'
> Suggestion. [arXiv:2511.15124](https://arxiv.org/abs/2511.15124) |
> [DOI: 10.1103/941r-ql3g](https://doi.org/10.1103/941r-ql3g)
>
> Received 2025-12-08, revised 2026-04-20, accepted 2026-05-21, published 2026-06-15.
> Original filename in inbox: `Meenu_PRB_paper.pdf`. Title, authors and DOI read off page 1 of
> the published PDF; arXiv ID (v1 2025-11-19, v2 2026-06-23) confirmed against arxiv.org.
> The filed PDF is the published PRB version, not the preprint.

## TL;DR
> Trotter-Suzuki fixes the exponents of a product formula by algebra; this paper lets a global
> action principle choose them instead, so $U_a(t)=\prod_j e^{ic_{2j}(t)A}e^{ic_{2j+1}(t)B}$ has
> time-dependent variational parameters obeying $\sum_k g_{jk}\dot c_k + iF_j = 0$. Because the
> variation is over the *operator* $U(t)$ rather than a state, one parameter set works for every
> initial state at fixed $H$ — and a 4-exponential variational ansatz beats Ruth's
> 7-exponential 4th-order formula (≈1.5% vs >3% error at $t=1$ on a 10-qubit Ising chain) while
> using about half the gates per layer.

## Why I'm reading this
- Landed in `91_Raw_Files/`; routed here for research relevance only — Juan confirmed that there is no application/supervisor angle, despite Kumari being at IQC
  Waterloo. Read it as an Adaptative_QITE paper.
- The hook is that its equation of motion is the *same object* this project keeps tripping over.
  Their $g_{jk}$ is a quantum geometric tensor and their $\sum_k g_{jk}\dot c_k = -iF_j$ is
  structurally the $M\dot\theta = Y$ solve from [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]]. So
  the conditioning and regularization questions that dominate the adaptive-step programme apply
  here too — and this paper, usefully, does *not* answer them (see §5).

## 1. Context & Motivation
- **Problem addressed:** Trotter-Suzuki decomposition of $U(t)=e^{-iHt}$ needs small time steps
  for accuracy, hence deep circuits, hence runtimes past NISQ coherence times. Higher-order
  Suzuki formulas buy accuracy by adding exponentials per layer, which costs depth back.
- **Why it matters:** Depth is the binding constraint. A product formula that is more accurate
  *per layer* lets you take larger $\tau$ and run fewer layers, cutting total gate count for
  fixed final accuracy.
- **Prior work / what this builds on:** McLachlan-style variational dynamics (VQE adapted to
  dynamics, varQITE) minimize the residual between $\partial_t|\psi(\theta)\rangle$ and
  $-iH|\psi(\theta)\rangle$ — i.e. they vary a *state*. The authors' complaint is that this makes
  accuracy state-dependent, breaks the exact unitary structure needed for error mitigation like
  PEC, and inherits VQA local-minima problems. Ref. [23] is their own prior action principle,
  which this paper builds on. Closely related in machinery to
  [[QITE_Hand_Off/Notes_on_Project/Understanding_Quantum_Imaginary_Time_Evolution_and_its_Variational_Form|Understanding_Quantum_Imaginary_Time_Evolution_and_its_Variational_Form]].

## 2. Core Idea
> Apply the variational principle to the time-evolution *operator* over a whole interval, not to
> a state at an instant. Minimizing $S=\int dt\,\mathrm{Tr}[U^\dagger(i\partial_t U - HU)]$ over
> a product-of-exponentials ansatz turns "what should the Trotter angles be?" into an ODE for
> those angles. The angles come out Hamiltonian-specific rather than universal, which is exactly
> why they beat a fixed algebraic construction of the same depth.

The payoff that a state-based method cannot match: solve the ODE once for a given $H$ and step
$\tau$, and the resulting $U_a(\tau)$ is a fixed circuit that works for *any* initial state. It
stays manifestly unitary, so it is still a legal gate sequence and still compatible with PEC.

## 3. Method / Approach
- **Key assumptions:**
  - $H$ is time-independent, and splits as $H=A+B$ (or $A+B+C$, …) where each block is a sum of
    *mutually commuting* terms — so each $e^{ic A}$ factorizes into parallel one- and two-qubit
    gates.
  - $c_k(0)=0$, which makes $U_a(0)=\mathbb{1}$ and keeps $U_a$ unitary throughout.
  - Traces $g_{jk}$, $F_j$ are computable — the load-bearing assumption, see §5.
- **Formalism / model used:** action principle → Euler-Lagrange EOM for the parameters; ODE
  integrated with RK45. Three test systems: a two-level toy model, the quantum Ising model (QIM)
  with transverse *and* longitudinal field, and the XXZ chain with next-nearest-neighbour
  coupling (nonintegrable).
- **Key equations:**
  - Action and EOM (Eqs. 2, 3):
    $$ S=\int dt\,\mathrm{Tr}\!\left[U^\dagger(i\partial_t U - HU)\right],\qquad
       \sum_k g_{jk}\dot c_k + iF_j = 0 $$
    with the quantum geometric tensor and generalized force
    $$ g_{jk}=\mathrm{Tr}\!\left[\Big(\tfrac{\partial U_a}{\partial c_j}\Big)^{\!\dagger}
       \tfrac{\partial U_a}{\partial c_k}\right],\qquad
       F_j=\mathrm{Tr}\!\left[\Big(\tfrac{\partial U_a}{\partial c_j}\Big)^{\!\dagger} H U_a\right] $$
  - The ansatz (Eq. 4), $M$ controlling depth:
    $$ U_a(t)=\prod_{j=0}^{M} e^{ic_{2j}(t)A}\,e^{ic_{2j+1}(t)B} $$
    Small-$t$ limit gives $c_{2j}=c_{2j+1}=-t/(M+1)+O(t^2)$, recovering Trotter-Suzuki — the
    consistency check that the method contains TS as its linearization (Eq. 10).
  - Workhorse ansätze: 3-exponential $e^{ic_0A}e^{ic_1B}e^{ic_2A}$ (Eq. 17, the analogue of
    symmetric 2nd-order TS) and 4-exponential $e^{ic_0A}e^{ic_1B}e^{ic_2A}e^{ic_3B}$ (Eq. 18).
  - Error metrics (Eqs. 27, 28), normalized so $0\le E_F\le 1$ with $D$ the Hilbert dimension:
    $$ E_F(t)=\tfrac{1}{2\sqrt{D}}\|U(t)-U_a(t)\|_F,\qquad
       E_F^{\mathrm{Strob}}(\tau,n)=\tfrac{1}{2\sqrt{D}}\|U(n\tau)-[U_a(\tau)]^n\|_F $$
    ($2\sqrt{D}$, not $\sqrt{2D}$ — corrected against the printed Eq. 27 p.6.
    $2\sqrt{D}$ is the one that gives $E_F\le1$, since $\|U-U_a\|_F^2\le 4D$ for unitaries.)
  - **Operator ordering matters**, unlike in TS. Residual
    $R=\|i\dot U_a - HU_a\|_F^2$ and ordering parameter $\Xi=\mathrm{sign}(R_{AB}-R_{BA})$
    (Eq. 65): $\Xi=-1$ selects $AB$. For the QIM, $\Xi=\mathrm{sgn}[4N(h_x^2-h_z^2)-J^2(N-1)]$
    (Eq. 68), so in the thermodynamic limit $AB$ is right when $h_x^2<h_z^2+J^2/4$.
  - Closed-form approximate parameters to cubic order — the "ready-to-use" deliverable. For the
    QIM 3-exponential ansatz (Eq. 69):
    $$ c_0=c_2\approx-\frac{t}{2}-\frac{h_x^2}{12}\,\frac{J^2(N-1)+2h_z^2N}{J^2(N-1)+4h_z^2N}t^3,
       \qquad c_1\approx-t+\frac{1}{24}\Big[\Big(1-\frac{1}{N}\Big)\frac{J^2}{2}+h_z^2\Big]t^3 $$
    Both have well-defined $N\to\infty$ limits. These are drop-in replacements for TS rotation
    angles at *zero* extra circuit cost.
- **Notation translation** (their symbol → my convention):

  | Their notation | Meaning | My notation |
  | --- | --- | --- |
  | $g_{jk}$ | Quantum geometric tensor; LHS matrix of the parameter ODE | $M$ (the metric in $M\dot\theta=Y$) |
  | $F_j$ (or $D_j$ in App. A) | Generalized force; RHS driving the parameter flow | $Y$ |
  | $c_j(t)$ | Variational parameters = circuit rotation angles | $\theta$ |
  | $t$ (real time), $\tau$ (stroboscopic step) | $U=e^{-iHt}$ | $\tau$ is imaginary time here — **sign/Wick caution**, their $\tau$ is my $\Delta\tau$, their $t$ has no QITE analogue |
  | $M$ | Number of ansatz layers | — (clashes with my $M$; call it $M_{\text{layers}}$) |
  | $\Xi$ | Ordering parameter, $\mathrm{sign}(R_{AB}-R_{BA})$ | no analogue |
  | $E_F$ | Normalized Frobenius operator error | no analogue (project measures energy error) |

  > [!warning] Notation collision
  > Their $M$ is a layer count; the project's $M$ is the metric matrix. Their metric is $g$.
  > When porting anything from this paper, rename before you regret it.

## 4. Results
All numbers below are read from the published PDF, with the figure/table they came from. Values
quoted off plots are marked as such.

| Target | Value | Source |
|---|---|---|
| Gate count per layer, variational ansatz Eq. (18), $N$ qubits | $R_x=2N$, $R_z=4N-2$, CNOT $=4N-4$ (2 $A$-exp, 2 $B$-exp) | Table I, p.8 |
| Gate count per layer, Ruth Eq. (38) | $R_x=4N$, $R_z=6N-3$, CNOT $=6N-6$ (4 $A$-exp, 3 $B$-exp) | Table I, p.8 |
| Reduction at $N=10$ | $R_x$ −50%, $R_z$ and CNOT ≈ −33% | text, p.8 (consistent with Table I) |
| QIM $E_F$ at $t=1.0$ ($h_x=h_z=J=1$, $N=10$) | variational ≈1.5%, Ruth >3% ("more than twice") | Fig. 4(a) + text, p.8 |
| QIM magnetization rel. error vs symmetric TS, $\tau=0.2,\,0.4$ | variational 4–5× better; gap widens with $\tau$ | Figs. 3(b), 3(d) |
| Two-level model, 3-exponential ansatz | $E_F\sim10^{-7}$ — *exact*, the residual is RK45 integrator noise | Fig. 2(b), $h_x=5$, $h_z=2$ |
| XXZ with NNN, long-time stroboscopic error at $t=100$, $\tau=0.2$ | TS ≈11%, variational ≈5.5% (twofold) | Fig. 11(b), $J_1=2$, $J_2=0.5$, $\Delta_{1,2}=0.2$, $N=5$ |
| Abstract headline | error reduced by "factors of 2 to 5" | abstract |
| Parameter convergence in $N$ | $c_j(\tau)$ converge smoothly to $c_j^\infty$; computed up to **1000 qubits** at $\tau=0.5$ | Figs. 6, 7 |

Why the $N=1000$ result is cheaper than it looks: for the QIM, $g$ and $F$ share a common $2^N$
factor that cancels in Eq. (3), leaving terms growing at most linearly in $N$ — so the *analytic*
route (Appendix C) costs $O(1)$, no matrix exponentials. That is a closed-form special case, not
a general scaling result.

Exact solution for the two-level model (Eqs. 32–33), $\Omega=\sqrt{h_x^2+h_z^2}$:
$$ c_0(t)=c_2(t)=\frac{1}{2h_x}\tan^{-1}\!\Big(\!-\frac{h_x\tan(\Omega t)}{\Omega}\Big),\qquad
   c_1(t)=\frac{1}{h_z}\sin^{-1}\!\Big(\!-\frac{h_z\sin(\Omega t)}{\Omega}\Big) $$
Linearizing gives $c_0=c_2\approx-t/2$, $c_1\approx-t$, the symmetric TS angles.

- **Figures worth remembering:**
  - Fig. 1 — algorithm flowchart: integrate the EOM on a classical machine, hand $c(\tau)$ to the
    circuit. Makes the classical/quantum split explicit.
  - Table I — the only pure-arithmetic claim in the paper; the natural replication target.
  - Fig. 13 — approximate cubic parameters (Eq. 69) vs symmetric TS, QIM $N=5$, $\tau=0.1$:
    the free-lunch claim, same circuit, better angles.
  - Fig. 14 — Krylov partial traces. Honest negative result: 1 basis vector is *worse* than the
    full trace over the whole window; 2 and 5 beat it only over the first few periods.

- **Replication status:** done — see [[QITE_Hand_Off/Notes_on_Project/Replication_Beyond_Trotterization|Replication_Beyond_Trotterization]]
  (code: `Replications/Beyond_Trot/replicate.py`, eight targets T1–T8). Summary: Table I,
  the two-level closed forms Eqs. (31)–(36), Fig. 4(a), Fig. 11(b) and Fig. 13's ratio all
  reproduce. Four things do not, and they are the interesting part:
  1. the p.8 **text** (1.5 % / "over 3 %") disagrees with the p.8 **figure** (0.0212 / 0.0692);
     this replication lands on the figure to ~1 % and misses the prose by 41–129 %;
  2. Eqs. (61)/(65)/(68)/(69) only hold with $A$ and $B$ **interchanged** relative to the
     Eq. (37)/(54) splitting — 0/15 points as printed, 15/15 relabelled;
  3. Fig. 13 plots the **raw** $\lVert\cdot\rVert_F$, not the Eq. (27)-normalized
     $\mathcal{E}_F^{\rm Strob}$ its axis claims (off by $2\sqrt{D}=11.3137$), while Fig. 11(b)
     at the same $D$ does normalize;
  4. the claim on p.8 that the variational ansatz beats Ruth "over the entire time interval
     $[0,2]$" is **false at small $t$** — Ruth is $O(t^5)$ and the four-exponential ansatz is
     $O(t^3)$, so Ruth wins below the crossover at $t = 0.376$ (fitted slopes 5.00 and 3.00).

## 5. Critical Evaluation
- **Strengths:**
  - The state-independence argument is real and is the honest differentiator from McLachlan/
    varQITE methods: one solve, any initial state, exact unitarity retained.
  - Eq. (69) is a genuinely free improvement — same gate count, same circuit, better angles.
    If it holds up, it is the most directly usable output of the paper.
  - They beat a *higher-order* formula (Ruth, 7 exponentials) with 4 exponentials and report
    the ordering dependence honestly instead of quietly picking the best one.
  - The Krylov partial-trace section reports a method that partly *fails*. Good sign.

- **Weaknesses / limitations / assumptions that might not hold:**
  1. **The classical cost is the whole ballgame, and it is never benchmarked.** $g_{jk}$ and
     $F_j$ are traces over a $2^N$-dimensional space. Three escape hatches are offered and none
     is demonstrated end-to-end: (a) analytic closed forms — only derived for the QIM and XXZ
     *with the particular splittings chosen*, Appendix C; (b) MPO contraction at $O(N\chi^3)$ —
     asserted by reference, never run; (c) Hadamard tests + LCU on hardware — asserted, and
     §II.C states outright that hardware implementation "is beyond the scope of the present
     work." So every number in the paper comes from exact classical simulation at $N\le10$
     (mostly $N=5$–$6$), or from the closed forms. The scalability claim rests on the analytic
     special case doing work it cannot do for a general $H$.
  2. **The conditioning of $g$ is never reported — and $g$ is singular at $t=0$.** This is the
     project-relevant defect. Any ansatz that repeats an operator (ABA, ABAB) has
     $T_1(c_1)=\mathrm{Tr}[Ae^{ic_1B}Ae^{-ic_1B}]\to\mathrm{Tr}[A^2]$ as $c_1\to0$, which makes
     rows 0 and 2 of $g$ (Eq. A5) identical at the initial time. Checked numerically on the
     paper's own Eq. (35) at its Fig. 2 parameters ($h_x=5$, $h_z=2$):

     | $c_1$ | $\det g$ | rank | $\mathrm{cond}(g)$ |
     |---|---|---|---|
     | $0$ | $0$ | **2** (of 3) | $\sim10^{17}$ (numerically infinite) |
     | $10^{-6}$ | $4.0\times10^{-8}$ | 3 | $2.5\times10^{11}$ |
     | $10^{-3}$ | $4.0\times10^{-2}$ | 3 | $2.5\times10^{5}$ |
     | $0.2$ | $1.29\times10^{3}$ | 3 | $1.06\times10^{1}$ |

     The null vector at $c_1=0$ is $\propto(1,0,-1)$: the *sum* $c_0+c_2$ is determined at
     $t=0$, the *difference* is not. $\mathrm{cond}(g)$ diverges as $c_1^{-2}$, i.e. as
     $t^{-2}$. The paper's $c_0=c_2$ is therefore a null-space gauge choice, taken silently on
     grounds of symmetry, and RK45 is started from precisely the singular point. It works —
     symmetry picks a consistent branch — but the paper never names the degeneracy, never
     states how the initial solve is regularized, and reports no $\mathrm{cond}(g)$ anywhere.
     For a 2- or 3-parameter ansatz this is benign; for the $2M+1$-parameter ansätze they
     propose for "exact representation", it is the same rank-deficiency wall documented in
     [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §4. *(Verified from the
     printed Eqs. (35) and (A5); reproducible from the paper alone, no simulation needed.)*
     **Now pinned** as target T6 of [[QITE_Hand_Off/Notes_on_Project/Replication_Beyond_Trotterization|Replication_Beyond_Trotterization]] — every value in the
     table above is hard-coded in `Replications/Beyond_Trot/replicate.py`, and
     `replicate.py --check` re-asserts them in about a second and exits nonzero on drift. (The
     gate pins the full-precision computed values, not the 3-s.f. figures printed here; matching
     only the rounded table needs a 2 % band, which is loose enough to hide a real regression.)
     The replication also adds what the table does not show: cond($g$) along the QIM trajectory
     at $N=5$ falls from $3.4\times10^{12}$ (BAB) and $2.1\times10^{16}$ (ABAB) at $t_0$ to 1.36
     and 4.85 by $t=2$, and $\mathcal{E}_F(t{=}2)$ moves by at most $5.12\times10^{-10}$ relative
     across $t_0\in[10^{-8},10^{-3}]$ and `rcond` $\in[10^{-14},10^{-4}]$. So the degeneracy is a boundary layer at
     the initialization point rather than a persistent obstacle — which is the substantive
     difference from the QDP case in [[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]].
  3. **No data or software.** The Data Availability statement says plainly: "There are no
     publicly available research data or software supporting this manuscript." Every figure
     would have to be rebuilt from Eqs. (3), (17), (18) and Appendices A/C to check. For a
     paper whose central claim is a numerical accuracy comparison, that is a real cost.
  4. **Error is measured against a fixed budget, not a fair-cost budget.** Comparisons are
     "per layer at fixed $\tau$" or "at equal exponential count". The argument that larger $\tau$
     compensates for richer layers is made (§II.B) but the paper does not present a clean
     total-gates-to-fixed-$E_F$ curve, which is the number a hardware user actually wants.
  5. **Noiseless throughout.** No decoherence anywhere. The claimed NISQ benefit (fewer gates →
     less noise) is inferred from gate counts, never simulated.
  6. **The Lie-algebraic "exactness" guarantee is not a resource bound.** $M$ is bounded by the
     dimension of the generated Lie algebra, itself bounded by $D^2-1=4^N-1$. Exponential. It
     guarantees a finite $M$ exists, not a useful one.

- **Open questions the paper leaves unanswered:**
  - What is $\mathrm{cond}(g)$ along the trajectory, at ansatz sizes past 4 exponentials?
  - Does the accuracy advantage survive a noise channel, or does the extra classical
    precomputation buy something the hardware immediately destroys?
  - For a general $H$ with no closed-form traces, what is the actual wall-clock crossover where
    the MPO route stops being cheaper than just running more TS layers?

- **What would break this result?** A demonstration that for Hamiltonians without analytic
  traces the ODE solve costs more than the circuit depth it saves — which the paper's own
  restriction to QIM/XXZ quietly leaves open. Also, ill-conditioning of $g$ at larger ansätze
  would make the "add more exponentials" prescription (their stated route to higher accuracy)
  fail exactly where it is needed.

## 6. Connections
- **Relates to:** [[QITE_Hand_Off/Notes_on_Project/Understanding_Quantum_Imaginary_Time_Evolution_and_its_Variational_Form|Understanding_Quantum_Imaginary_Time_Evolution_and_its_Variational_Form]] —
  the closest methodological sibling. Both derive a linear system for parameter velocities from
  a variational principle; that note's varQITE does it for a *state* in imaginary time, this one
  for the *operator* in real time. Reading them together isolates what the state-vs-operator
  choice actually buys, which is the cleanest statement of this paper's contribution.
- **Relates to:** [[QITE_Hand_Off/Notes_on_Project/Trading_Imaginary_Time_for_Randomness_in_Ground_State_Preparation|Trading_Imaginary_Time_for_Randomness_in_Ground_State_Preparation]] — the
  most load-bearing connection. TITE appends *randomly-timed real-time evolution* after ITE, so
  it needs a real-time propagator; this paper makes that propagator cheaper and more accurate at
  fixed depth. The two compose without touching each other's machinery: TITE is agnostic to how
  $e^{-iHt}$ is implemented, and this is a drop-in $e^{-iHt}$.
- **Contradicts / complicates:**
  [[QITE_Hand_Off/Notes_on_Project/Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors|Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors]] — not a
  contradiction so much as a lever on its ceiling. That paper bounds relative error in
  $\dot\theta$ by $\mathrm{cond}(M)$ and the parameter count/depth $N$; fewer gates per layer at
  equal accuracy pushes on the $N$ in those bounds directly. But it cuts both ways: this paper
  is noiseless and reports no $\mathrm{cond}(g)$, so it is silent on exactly the failure mode
  that one documents. See §5 item 2.
- **Follow-up papers to read:** arXiv:2607.29418, a Floquet/time-dependent extension of the same
  framework by overlapping authors (surfaced via the arXiv listing, *not* cited in this PDF —
  TODO: verify it is genuinely the same programme before relying on it).

## 7. My Takeaways
- The transferable idea is not the ansatz, it's the framing: **the regularization/gauge choice in
  a geometric-tensor solve is a free parameter that papers keep making silently.** This one picks
  $c_0=c_2$ by symmetry at a point where $g$ is provably rank-deficient; the project's own audit
  history says the pseudo-inverse cutoff has to be promoted to an explicit, logged component of
  the method. Same lesson, different algorithm.
- Eq. (69) is the piece worth stealing. Closed-form cubic corrections to TS angles, no extra
  gates, well-defined thermodynamic limit. If it reproduces, it is a strictly better default
  than symmetric TS for any QIM real-time sub-block.
- The honest verdict on scalability: this is a small-$N$ / analytically-tractable-Hamiltonian
  result presented as a general method. Worth watching, not yet worth building on.
- Tractable replication targets, in order — see §4 for the target values:
  1. **Table I gate counts.** Pure arithmetic, no simulation, verifies the headline efficiency
     claim and the −50%/−33% figures at $N=10$ in one function.
  2. **Fig. 4(a) at $t=1.0$.** Needs a real $N=10$ QIM propagation plus the RK45 parameter solve,
     but $E_F$ is well-defined and the target (≈1.5% vs >3%) is stated in text, not eyeballed.

## Glossary / New Terms
- **Variational product formula** — a Trotter-like product of exponentials whose exponents are
  time-dependent parameters fixed by a variational principle rather than by algebraic order
  conditions.
- **Quantum geometric tensor (QGT)** — $g_{jk}=\mathrm{Tr}[(\partial_{c_j}U_a)^\dagger
  \partial_{c_k}U_a]$; the metric on the ansatz's parameter manifold. The project's $M$.
- **Ordering parameter $\Xi$** — $\mathrm{sign}(R_{AB}-R_{BA})$; tells you which operator
  ordering to use, computable from static traces before you simulate anything.
- **Ruth's formula** — a 4th-order symmetric TS composition, $U_{\rm Ruth}(t)=U^{(2)}_{\rm TS}(pt)
  U^{(2)}_{\rm TS}(qt)U^{(2)}_{\rm TS}(pt)$ with $p=(2-2^{1/3})^{-1}$, $q=1-2p$. Seven
  exponentials per step; the paper's main comparison baseline.
- **SSS / TSS** — single-step split vs two-step split, the two ways of handling $H=A+B+C$:
  vary all blocks at once (Eq. 19), or nest a second variational decomposition inside
  $e^{ic(B+C)}$ (Eqs. 20–22). SSS is less constrained and wins at long times.
- **Stroboscopic error** — $E_F^{\rm Strob}(\tau,n)$, error of $[U_a(\tau)]^n$ against
  $U(n\tau)$; the quantity that matters for long simulations built from a repeated layer.

---
