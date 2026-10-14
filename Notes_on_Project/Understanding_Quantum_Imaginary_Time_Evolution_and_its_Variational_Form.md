---
fileClass: Paper_Reading
title: "Understanding Quantum Imaginary Time Evolution and its Variational form"
authors: Andreu Anglés-Castillo, Luca Ion, Tanmoy Pandit, Rafael Gomez-Lurbe, Rodrigo Martínez, Miguel Angel Garcia-March
year: 2025
venue: arXiv
arxiv_id: 2510.02015
doi:
status: reading
tags: [paper, physics, quantum-computation, paper-notes]
---

# Understanding Quantum Imaginary Time Evolution and its Variational form

> [!info] Citation
> Anglés-Castillo, A., Ion, L., Pandit, T., Gomez-Lurbe, R., Martínez, R., Garcia-March, M.A.,
> "*Understanding Quantum Imaginary Time Evolution and its Variational form*," arXiv, 2025.
> [arXiv:2510.02015](https://arxiv.org/abs/2510.02015) | DOI: none (preprint only)
>
> Original filename in inbox: `2510.02015v2.pdf`. Metadata (arXiv ID, date, author list)
> extracted from the PDF itself — cross-check against arxiv.org before citing further.

## TL;DR
> A review/tutorial paper that re-derives the original QITE algorithm (Motta et al. 2019) and
> its variational counterpart varQITE (McArdle et al. 2019) step by step, releases a working
> Python implementation of both, and benchmarks them side by side on an 8-qubit transverse-field
> Ising model to compare convergence speed, gate cost, and failure modes.

## Why I'm reading this
- Landed in `91_Raw_Files/` inbox. Routed here because it's a derivation-level companion to
  [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]]: it works out the same McLachlan-derived
  $M\dot\theta = Y$-type linear system the project's noise-sensitivity analysis studies, gives
  explicit pseudocode/circuits for evaluating it, and ships open-source code that could serve as
  a noiseless baseline before layering in Saxena et al.'s noise model.

## 1. Context & Motivation
- **Problem addressed:** Ground-state estimation for a Hamiltonian $H$ is exponentially hard by
  direct diagonalization. QITE (Motta et al.) approximates the non-unitary imaginary-time
  evolution operator with a sequence of unitaries implementable on quantum hardware; varQITE
  (McArdle et al.) instead evolves the parameters of a fixed variational ansatz to track the
  same imaginary-time trajectory. This paper reviews both derivations from scratch, in
  consistent notation, and compares them on one model.
- **Why it matters:** Serves as a citable derivation/tutorial reference and a working code base
  (two Jupyter notebooks, linked in the paper) for exactly the algorithm family the
  `Adaptative_QITE` project analyzes under noise — useful both for pinning down notation and as
  a sanity-check implementation.
- **Prior work / what this builds on:** [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] analyzes
  noise sensitivity of the same varQITE update this paper derives (there via
  [[QITE_Hand_Off/Notes_on_Project/Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors|Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors]], Saxena,
  Shalabi & Kyaw). Directly builds on Motta et al. 2019
  (*Determining eigenstates and thermal states on a quantum computer using quantum imaginary
  time evolution*, Nature Physics) for QITE and McArdle et al. 2019 (*Variational ansatz-based
  quantum simulation of imaginary time evolution*, npj Quantum Information) for varQITE.

## 2. Core Idea
> QITE approximates the (non-unitary, hence non-hardware-implementable) imaginary-time step
> $e^{-\Delta\tau H}$ by finding, at each small time step, the Hermitian generator $A_m$ (written
> as a sum of Pauli strings) whose unitary $e^{-i\Delta\tau A_m}$ best matches the true ITE step
> in norm — this reduces to a linear system for the Pauli coefficients of $A_m$. varQITE instead
> fixes a parameterized ansatz circuit $|\psi(\theta)\rangle$ and uses McLachlan's variational
> principle to derive an ODE for $\theta(\tau)$, evaluated at each step via a different linear
> system solved from quantum-circuit overlap measurements (ancilla-based Hadamard tests for $M$,
> parameter-shift rule for $V$).

## 3. Method / Approach
- **Key assumptions:** Initial state has nonzero overlap with the ground state (Eq. 1). QITE's
  derivation linearizes the unitary difference for small $\Delta\tau \to 0$ (Eq. 8). For QITE
  scaling: the Hamiltonian is decomposed into $T$-local pieces (Trotterization) and each local
  unitary update is further restricted to act on a domain of $D$ qubits around the piece it
  approximates — an explicit approximation, not exact, that trades fidelity for feasibility.
  varQITE's ansatz here is a fixed hardware-efficient circuit: $R_y$-then-$R_z$ rotation layer
  per qubit (the tunable parameters) followed by a CNOT entangling ladder, optionally repeated;
  $2N_{\text{qubits}}$ parameters per layer.
- **Formalism / model used:** McLachlan's variational principle
  $\delta\|(H-E(\theta))|\psi(\theta)\rangle\| = 0$ for varQITE; direct norm-minimization of the
  ITE-step/unitary-step difference for QITE.
- **Key equations:**
  QITE, solved each time step for the Pauli coefficients $a$ of $A_m = \sum_I a_I \sigma_I$:
  $$(S+S^T)a = -b, \qquad S_{IJ} = \langle\phi_{m-1}|\sigma_I\sigma_J|\phi_{m-1}\rangle, \qquad b_I = \frac{i}{c\Delta\tau}\big(\langle\phi_{m-1}|e^{-H\Delta\tau}\sigma_I - \sigma_I e^{-H\Delta\tau}|\phi_{m-1}\rangle\big)$$
  varQITE, solved each time step for $\dot\theta$:
  $$M_{ij}\dot\theta_j = -V_i, \qquad M_{ij} = \mathrm{Re}\!\left[\frac{\partial\langle\psi|}{\partial\theta_i}\frac{\partial|\psi\rangle}{\partial\theta_j} + \frac{\partial|\psi\rangle}{\partial\theta_i}\frac{\partial\langle\psi|}{\partial\theta_j}\bigg|\psi\rangle\langle\psi\bigg|\right],\quad V_i = \mathrm{Re}\!\left[\frac{\partial\langle\psi(\theta)|}{\partial\theta_i}H|\psi(\theta)\rangle\right]$$
  $$\theta_{t+1} = \theta_t - \eta M^{-1}V$$
- **Notation translation** (their symbol → meaning → project's convention, per `03_LG_Project/CLAUDE.md`):
  | Their notation | Meaning | My notation |
  |---|---|---|
  | $\beta$ | imaginary time, $\beta = it$ | $\tau$ ($\tau = it$, same convention) |
  | $M_{ij}\dot\theta_j = -V_i$ | varQITE update linear system | $M\dot\theta = Y$ — **sign flip**: this paper's $V$ is the project's $-Y$; check sign when cross-referencing [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE\|Technical_Resume_Project_Limiting_QITE]] |
  | $N$ (Sec. 2, spin count) | number of qubits/sites in the system | **clashes** with project's $N$ = ansatz parameter count; this paper's parameter count is separately $2N_{\text{qubits}}$ per layer — always disambiguate which $N$ when citing both notes |
  | $D$ | QITE unitary-support domain size (qubits) | new object, no prior analogue in project notes; distinct from $\mathrm{cond}(M)$/$p$ used in the noise-bound paper |

## 4. Results
- **Main result:** On an $N=8$-qubit TFIM ($H = -J\sum Z_iZ_{i+1} + g\sum X_i$, exact $g$ value
  not stated in text — TODO: verify from code/repo), QITE (domain sizes $D=2,4,6$,
  $\Delta\tau=0.25$) and varQITE (2-layer $R_y$-$R_z$+CNOT ansatz, 32 parameters, all initialized
  to $\pi/3$, $\Delta\tau=0.05$) both converge in energy and fidelity to the exact ground state
  (from exact diagonalization) on a comparable time scale — see Fig. 1. No numeric convergence
  values are tabulated in the text; the only quantitative record is the plot itself.
- **Numerical/experimental setup:** TFIM, $N=8$ qubits (up to 20 qubits reported feasible on
  HPC, not shown). QITE: $T=2$-local Hamiltonian pieces, domains $D \in \{2,4,6\}$, initial state
  $|\!\uparrow\rangle^{\otimes N}$. varQITE: EfficientSU2-style ansatz via Qiskit, 2 repetitions
  (32 parameters), all initialized at $\pi/3$. Code released as two Jupyter notebooks at
  `github.com/dark-dryu/understanding-QITE-varQITE` (QITE via SciPy sparse — runs to ~16 qubits
  on a desktop; varQITE via Qiskit's `EfficientSU2`).
- **Figures worth remembering:**
  - Fig. 1 (left) — energy vs. time for exact ITE, QITE ($D=2,4,6$), and varQITE: QITE's
    accuracy improves monotonically with $D$; no numeric values given, read-off only.
  - Fig. 1 (right) — fidelity to exact ground state vs. time for the same runs; similar
    convergence trend to the energy panel.
- **Discussion (Sec. 4.1, qualitative but load-bearing):** QITE's convergence depends directly
  on domain size $D$; measurement cost per time step scales as $\sim 4^D$ Pauli-string
  expectation values, and even $D=2$ needs ≥3 time steps for energy convergence with "very bad"
  fidelity — larger $D$ is quickly infeasible on current hardware (cites Table 1 of Li et al.,
  *Paulihedral*, [10], for concrete gate counts, not reproduced here). varQITE's accuracy instead
  depends on ansatz choice and parameter initialization, not a single tunable knob like $D$; the
  paper explicitly distinguishes convergence-to-local-minimum (insufficient
  expressivity/initialization) from the barren-plateau phenomenon (exponentially vanishing
  gradient *and* variance) — no BP signature was observed in their $N=8$ run, but this does not
  generalize. varQITE only needs the energy expectation value per step, versus QITE's $O(4^D)$
  Pauli-string measurements — a real cost advantage for varQITE at fixed accuracy.

## 5. Critical Evaluation
- **Strengths:** Consistent, worked-through derivation of both QITE and varQITE from the same
  starting point, letting one map notation directly between them; open-source code released
  alongside the paper; useful, non-hand-wavy discussion distinguishing local-minimum convergence
  failure from genuine barren plateaus, a distinction often blurred in the literature.
- **Weaknesses / limitations / assumptions that might not hold:** Entirely noiseless — no
  hardware/gate-error model at all, so it says nothing about the noise-driven condition-number
  blowup [[QITE_Hand_Off/Notes_on_Project/Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors|Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors]]
  documents for the exact same varQITE update. Single system
  size ($N=8$), single random seed for the ansatz's initial parameters (all set to $\pi/3$, no
  averaging over initializations) — no evidence this generalizes to other initializations or
  larger $N$. No numeric table of results, only a plot, so exact reproduction requires reading
  values off Fig. 1 (would need `pdf_text.py --images` or the released code, not the paper text
  alone). Gate-cost comparison for QITE leans on a citation (Paulihedral, [10]) rather than
  restating the numbers.
- **Open questions the paper leaves unanswered:** No systematic ansatz-design recipe for varQITE
  is given — acknowledged directly by the authors. The "natural extension" they propose (use
  QITE with small $D$ to get close to the ground state, then warm-start a variational ansatz
  from that state to avoid barren plateaus) is suggested but not implemented or tested here.
- **What would break this result?** Larger system sizes or deeper/more expressive ansätze, where
  varQITE's local-minimum or barren-plateau failure modes (not observed at $N=8$) could appear;
  for QITE, any domain size $D$ large enough to be practically infeasible on real hardware —
  exactly the noise/gate-count ceiling [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] quantifies.
  Since this paper is noiseless throughout, none of its convergence claims account for the
  hardware error sources that dominate the project's own analysis.

## 6. Connections
- **Relates to:** [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] — this paper derives, in explicit
  step-by-step form, the same $M\dot\theta = Y$-type varQITE update (here $M\dot\theta=-V$,
  note the sign flip above) that the technical resume analyzes for noise sensitivity; a useful
  derivation-level reference for that note's formalism. Also relates to
  [[QITE_Hand_Off/Notes_on_Project/Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors|Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors]]
  (Saxena, Shalabi & Kyaw) — same underlying varQITE update, this paper is the noiseless
  tutorial/derivation companion to that paper's noise-bound analysis. Also relates to
  [[QITE_Hand_Off/Notes_on_Project/Trading_Imaginary_Time_for_Randomness_in_Ground_State_Preparation|Trading_Imaginary_Time_for_Randomness_in_Ground_State_Preparation]] — same ITE-family
  subject area, though that paper targets plain/Trotterized ITE cost scaling in $\beta$ rather
  than the variational $M\dot\theta=Y$ update this one derives.
- **Contradicts / complicates:** [[]] — none identified.
- **Follow-up papers to read:** Motta et al. 2019 (*Determining eigenstates and thermal states on
  a quantum computer using quantum imaginary time evolution*, Nature Physics, ref [8] here — the
  original QITE algorithm) and McArdle et al. 2019 (*Variational ansatz-based quantum simulation
  of imaginary time evolution*, npj Quantum Information, ref [9] — the original varQITE paper),
  both cited as the primary sources this review builds on but not yet filed in the vault.

## 7. My Takeaways
- Two open-source Jupyter notebooks (`github.com/dark-dryu/understanding-QITE-varQITE`) implement
  both QITE (SciPy sparse) and varQITE (Qiskit `EfficientSU2`) — worth checking as a noiseless
  starting point/sanity check before attempting to replicate Saxena et al.'s noise bounds
  numerically.
- Notation clash to remember when cross-referencing project notes: this paper's $N$ is qubit/spin
  count, while the project convention (per `03_LG_Project/CLAUDE.md`) uses $N$ for ansatz
  parameter count. Also watch the sign: $M\dot\theta=-V$ here vs. $M\dot\theta=Y$ elsewhere.
- The local-minimum-vs-barren-plateau distinction (Sec. 4.1) is a clean, citable framing worth
  reusing when writing up the project's own varQITE numerics.

## Glossary / New Terms
- **Domain size $D$ (QITE)** — number of qubits the local approximate unitary $A_m[k]$ is allowed
  to act on, centered on the qubits touched by Hamiltonian piece $h_k$; larger $D$ improves
  fidelity to true ITE at the cost of $O(4^D)$ Pauli-string measurements per step.
- **EfficientSU2** — Qiskit's built-in hardware-efficient parameterized ansatz family, used here
  for the varQITE implementation.
- **Local-minimum convergence vs. barren plateau** — insufficient ansatz expressivity/poor
  initialization can trap optimization at a local minimum without the gradient itself vanishing
  exponentially (which is what defines a true barren plateau); the paper argues these are
  conceptually distinct failure modes often conflated in the literature.

---
