---
fileClass: Paper_Reading
title: "Trading Imaginary Time for Randomness in Ground State Preparation"
authors: Alvan Arulandu, John M. Martyn, Isaac L. Chuang
year: 2026
venue: arXiv
arxiv_id: 2608.00443
doi:
status: reading
tags: [paper, physics, quantum-computation, paper-notes]
---

# Trading Imaginary Time for Randomness in Ground State Preparation

> [!info] Citation
> Alvan Arulandu, John M. Martyn, Isaac L. Chuang, "*Trading Imaginary Time for Randomness in
> Ground State Preparation*," arXiv, 2026. [arXiv:2608.00443](https://arxiv.org/abs/2608.00443) | DOI: TODO: verify
>
> Original filename in inbox: `2608.00443v1.pdf`. Metadata (arXiv ID, date, author list)
> extracted from the PDF itself — cross-check against arxiv.org before citing further.

## TL;DR
> Imaginary-time evolution (ITE) for ground state preparation costs $e^{O(\beta)}$ in
> sample/classical complexity, and its energy accuracy doesn't guarantee accuracy for other
> observables. The authors introduce twirled imaginary-time evolution (TITE): follow ITE with
> real-time evolution for a randomly-drawn duration $t$, which quadratically suppresses the
> trace distance to the ground state ($\epsilon \to O(\epsilon^2)$) at no extra sample/classical
> cost, letting $\beta \to \beta/2$ for the same accuracy.

## Why I'm reading this
- Landed in `91_Raw_Files/` inbox; routed here because it's directly about imaginary-time
  evolution for ground-state preparation and analyzes error/cost scaling in $\beta$ — same
  subject area as [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]]. TODO: confirm with Juan why this
  specific paper was queued up (LG collaboration reading list vs. general interest).

## 1. Context & Motivation
- **Problem addressed:** Standard ITE implementations (Trotterized post-selection, QITE,
  quantum signal processing) all incur $e^{O(\beta)}$ sample complexity and/or classical cost
  because the ITE operator $e^{-\beta H}$ is non-unitary. Additionally, low error in the ground
  state *energy* does not imply low error in other observables — trace distance to the true
  ground state is the quantity that controls worst-case observable error.
- **Why it matters:** The exponential cost in $\beta$ is the practical bottleneck for any ITE
  deployment; a way to reach the same trace-distance accuracy with a smaller $\beta$ would give
  a direct, implementation-agnostic reduction in cost for every ITE variant (Trotterization,
  QITE, QSP-based ITE).
- **Prior work / what this builds on:** Randomized compiling and the Campbell–Hastings mixing
  lemma (quadratic error suppression from averaging over unitaries close to a target). Also
  distinguishes itself from Ref. [54]'s randomized-Trotterization + QITE hybrid (gate-count
  reduction) — TITE instead adds real-time evolution *on top of* any ITE implementation to
  refine the resulting state, so the two approaches are compatible/stackable. Directly relevant
  to the noise/cost-scaling questions raised in [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]],
  though that note's condition-number/McLachlan analysis addresses a different failure mode
  (ill-conditioned $M$ in the variational QITE update) than TITE's (exponential cost in
  Trotterized/post-selected ITE more broadly).

## 2. Core Idea
> Imaginary-time evolution for time $\beta$ leaves the state $O(\epsilon)$-close to the ground
> state, with the leftover error concentrated as a *coherent* superposition over excited states
> carrying relative phases $e^{-i(\lambda_k-\lambda_0)t}$. Applying real-time evolution
> $e^{-iHt}$ for a *randomly drawn* $t$ and averaging over many draws scrambles those relative
> phases toward zero (their expectation is $O(\epsilon)$ by construction of the sampling
> distribution), turning the coherent $O(\epsilon)$ error into an incoherent one that shrinks
> to $O(\epsilon^2)$ — the same quadratic cancellation mechanism as randomized compiling, but
> applied to states via a "mixing lemma for states" rather than to unitaries. Because real-time
> evolution is unitary, this refinement is free in sample complexity: only the additional
> unitary gates cost anything.

## 3. Method / Approach
- **Key assumptions:** Hamiltonian $H$ is gapped and non-degenerate (spectral gap
  $\Delta = \lambda_1-\lambda_0 > 0$, known or lower-bounded); initial state $|\psi\rangle$ has
  non-negligible ground-state overlap $|\langle\psi|\lambda_0\rangle| > c_{\min} > 0$. Result is
  about *algorithmic* error only — doesn't model gate/hardware noise beyond the numerical
  experiments in Sec. IV.
- **Formalism / model used:** Trace distance $d_{tr}(\rho,\sigma) = \frac12\|\rho-\sigma\|_1$ as
  the error metric (bounds worst-case observable error via
  $|\langle\psi|O|\psi\rangle - \langle\lambda_0|O|\lambda_0\rangle| \le 2\|O\|_{op}\, d_{tr}$).
  Mixing lemma for states (their Lemma 2, a state-space analogue of the Campbell–Hastings
  mixing lemma for unitaries) is the core error-suppression tool.
- **Key equations:**
  $$|\psi(\beta)\rangle = \frac{e^{-\beta H}|\psi(0)\rangle}{\|e^{-\beta H}|\psi(0)\rangle\|_2}, \qquad \beta = \frac{\log(1/(c_{\min}\epsilon))}{\Delta} \;\Rightarrow\; d_{tr}(|\psi(\beta)\rangle\langle\psi(\beta)|,\,|\lambda_0\rangle\langle\lambda_0|) \le O(\epsilon)$$
  $$|\psi(\beta,t)\rangle \equiv e^{-iHt}|\psi(\beta)\rangle, \qquad \rho \equiv \mathbb{E}_{t\sim D}\big[|\psi(\beta,t)\rangle\langle\psi(\beta,t)|\big]$$
  $$\text{if } \big|\mathbb{E}_{t\sim D}[e^{-i\omega t}]\big| \le O(\epsilon)\ \ \forall\, \omega \ge \Delta \quad\Rightarrow\quad d_{tr}(\rho,\,|\lambda_0\rangle\langle\lambda_0|) \le O(\epsilon^2)$$
  Optimal twirling distribution (compact support in characteristic-function space, requires
  known $\Delta$):
  $$f_{S_\Delta}(t) \equiv \frac{3\Delta}{8\pi}\,\mathrm{sinc}^4(\Delta t/4)$$
  and two fallback distributions not requiring exact $\Delta$: $C_\beta$ (Cauchy-like, no
  spectral-gap knowledge needed) and $N_{\beta,\Delta}$ (ordinary Gaussian, suboptimal but
  simplest to simulate).
- **Notation translation** (their symbol → my convention):
  | Their notation | Meaning | My notation |
  |---|---|---|
  | $\beta$ | imaginary time | $\tau$ (per [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE\|Technical_Resume_Project_Limiting_QITE]] convention, $\tau = it$) |
  | $\Delta$ | spectral gap $\lambda_1-\lambda_0$ | $\Delta$ (consistent) |
  | $D$, $f_D(t)$ | twirling distribution over real-evolution duration $t$ | — (new object, no prior analogue in project notes) |
  | $N$ | *not used this way here* — paper's $N$ elsewhere is spin-chain size $n$ | avoid confusion with project's $N$ = ansatz parameter count |

## 4. Results
- **Main result (Theorem 3):** TITE achieves $d_{tr}(\rho,|\lambda_0\rangle\langle\lambda_0|) \le O(\epsilon^2)$
  using only the same $\beta = O(\log(1/\epsilon)/\Delta)$ imaginary time that plain ITE needs to
  reach $O(\epsilon)$ — equivalently, for a fixed target trace-distance accuracy, $\beta$ can be
  halved ($\beta \to \beta/2$), giving a **quadratic reduction** in the $e^{O(\beta)}$
  sample/classical cost of any black-box ITE implementation (Trotterization, QITE, QSP-based).
- **Optimality (Appendix A, Lemma 4):** Any randomization scheme using unitaries that stabilize
  the ground state up to phase (which real-time evolution does) cannot suppress trace-distance
  error below $\delta^2$ — so the quadratic suppression TITE achieves is provably optimal for
  this class of method, not just an empirical observation.
- **Numerical/experimental setup:** 1D mixed-field Ising chain, open boundary conditions,
  $H = -J\sum_{i=1}^{n-1} Z_iZ_{i+1} - h_x\sum_j X_j - h_z\sum_j Z_j$, $n=10$ spins,
  $J=1$, $h_x=1.4$, $h_z=0.6$ (non-integrable point). 2nd-order Trotter-Suzuki decomposition,
  $r_{RTE}=100$ / $r_{ITE}=15$ Trotter steps. Noisy runs: Qiskit Aer density-matrix simulator
  with per-gate depolarizing noise $\gamma \in \{0, 10^{-6}, 10^{-5}\}$. Noiseless runs: exact
  sparse-matrix simulation via NumPy, to isolate Trotterization error and choice of $D$ from
  gate noise.
- **Figures worth remembering:**
  - Fig. 5 — noisy circuit simulation ($D=C_\beta$): TITE's trace-distance/energy error still
    quadratically suppressed under gate noise $\gamma$ up to $10^{-5}$, though the advantage
    narrows as $\gamma$ grows (more real-time gates = more noise exposure); spin covariance
    ($C_{56}^{xx}, C_{56}^{zz}$) and magnetization ($m^x, m^z$) converge closer/sooner under
    TITE than ITE, and are comparatively robust to $\gamma$.
  - Fig. 6 — noiseless, isolates Trotterization error: both ITE and TITE plateau beyond a
    critical $\beta^c_{ITE}$ once ITE-Trotterization error dominates; TITE's advantage similarly
    plateaus beyond a critical $\beta^c_{RTE}$ where real-time-evolution Trotterization error
    dominates.
  - Fig. 7 — compares $S_\Delta$, $C_\beta$, $N_{\beta,\Delta}$ (noiseless): all three double the
    trace-distance decay rate ($e^{-\beta\Delta}\to e^{-2\beta\Delta}$); by $\beta=2$, TITE error
    sits ~3 orders of magnitude below ITE, matching what plain ITE would need $\beta\approx4$
    to reach (i.e., quadratically greater sample cost) — $S_\Delta$ best, but margin over
    $C_\beta$/$N_{\beta,\Delta}$ (even with $\Delta$ mis-estimated by $\pm20\%$) is small.
  - Fig. 8 — noiseless, decomposes error into population (diagonal, $\propto e^{-2\beta\Delta}$,
    untouched by real-time evolution) vs. coherent (off-diagonal, $e^{-\beta\Delta}$ for ITE,
    suppressed to population-error level by TITE) — makes the "structure vs. randomness"
    argument (Sec. III F) numerically visible.

## 5. Critical Evaluation
- **Strengths:** Fully implementation-agnostic — proof only assumes the ITE state has been
  prepared to $O(\epsilon)$ trace distance, so it composes with any existing ITE method
  (including the QITE variant analyzed in [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]]) without
  modifying it. Optimality proof (Appendix A) is a genuine strength: rules out chasing
  super-quadratic gains from this class of randomization, so effort is better spent elsewhere
  (e.g. the purification/clock-register idea sketched in Discussion). Numerics test both gate
  noise and Trotterization error as separate, isolated failure modes rather than conflating
  them.
- **Weaknesses / limitations / assumptions that might not hold:** (1) Requires a gapped,
  non-degenerate Hamiltonian with a known or lower-boundable spectral gap $\Delta$ — degenerate
  or near-degenerate spectra (common in frustrated/symmetric systems) need the extension
  mentioned only in passing ("randomization conditioned on gap to target subspace"), not
  worked out. (2) The extra real-time evolution *increases circuit depth*, and Fig. 5 shows the
  quadratic advantage narrows as gate noise increases — the authors state this but don't
  characterize the noise level at which TITE stops being worth it. (3) Real-time evolution cost
  overhead $\mathbb{E}_{t\sim D}[|t|^\alpha]$ scales as $\Delta^{-\alpha}$ for the optimal
  $S_\Delta$ — for small spectral gaps (exactly the regime where ITE is already expensive) the
  real-time evolution itself gets expensive, and the paper doesn't quantify the crossover where
  this eats the quadratic gain. (4) All benchmarks are a single 10-spin non-integrable Ising
  chain — no scaling study across system size, so how the noise-robustness/advantage trade-off
  moves with $n$ is unverified.
- **Open questions the paper leaves unanswered:** How does the TITE-vs-noise crossover depend
  on system size $n$ and spectral gap $\Delta$? What's the concrete recipe for the degenerate
  ground-state-subspace extension? Is the purification-via-clock-register idea (Sec. V)
  actually implementable with near-term hardware, and does it truly beat post-selection cost
  for super-quadratic suppression?
- **What would break this result?** A Hamiltonian with vanishing/unknown spectral gap (breaks
  the twirling-distribution construction and the $\beta$ scaling), an initial state with
  negligible ground-state overlap $c_{\min}$, or hardware where two-qubit gate noise dominates
  enough that the extra real-time-evolution depth erases the trace-distance gain before it
  materializes (per Fig. 5's trend at $\gamma=10^{-5}$).

## 6. Connections
- **Relates to:** [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] — same subject (ITE cost/error
  scaling with $\beta$), complementary failure mode (this paper: exponential-cost ITE
  implementations broadly; resume note: ill-conditioned $M$ in the variational QITE update).
  Also relates to [[QITE_Hand_Off/Notes_on_Project/Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors|Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors]]
  (Saxena, Shalabi & Kyaw, "Practical limitations
  of quantum data propagation on noisy quantum processors") — both analyze error/cost
  limitations of ITE-family algorithms under hardware noise, though from different angles
  (variational-update noise sensitivity there vs. implementation-agnostic cost scaling here).
- **Contradicts / complicates:** [[]] — none identified; TITE is explicitly framed as
  composable with, not competing against, variational/QITE-style approaches.
- **Follow-up papers to read:** [[]] — Ref. [54] (randomized Trotterization + QITE hybrid) and
  Ref. [22] (the QITE algorithm this project's technical resume analyzes) would be the natural
  next reads to check whether TITE's real-time-evolution twirl stacks cleanly on top of the
  McLachlan-variational QITE update characterized there.

## 7. My Takeaways
- The headline number to remember: $\beta \to \beta/2$ for the same trace-distance accuracy,
  i.e. a quadratic reduction in $e^{O(\beta)}$ cost, purely from appending randomly-timed
  real-time evolution — no change to the underlying ITE implementation.
- Worth checking whether this composes with the project's noise-bound analysis: if TITE can be
  layered on top of variational QITE without touching the $M\dot\theta=Y$ update, it's a
  cost-side lever that's orthogonal to the condition-number/noise-floor problem the technical
  resume is about.
- The "structure vs. randomness" framing (Sec. III F, Fig. 8) is a clean way to think about why
  energy converges faster than other observables under plain ITE — coherent vs. incoherent
  error decomposition.

## Glossary / New Terms
- **Twirled imaginary-time evolution (TITE)** — ITE followed by real-time evolution for a
  randomly drawn duration, to quadratically suppress trace distance to the ground state.
- **Temporal twirling** — the randomized real-time evolution step itself, named by analogy to
  Pauli/randomized-compiling "twirls."
- **Mixing lemma for states (Lemma 2)** — state-space generalization of the Campbell–Hastings
  mixing lemma: averaging over states each $O(\epsilon)$-close to a target, with an
  appropriately small first-moment deviation, yields a mixed state $O(\epsilon^2)$-close to the
  target.
- **Characteristic function $\varphi_D(\omega) = \mathbb{E}_{t\sim D}[e^{-i\omega t}]$** — Fourier
  transform of the twirling distribution's density; must decay below $O(\epsilon)$ for all
  $\omega \ge \Delta$ for TITE's error suppression to hold.
- **Structure versus randomness** — additive-combinatorics dichotomy (structured/patterned
  component vs. pseudorandom component) used here to explain why the ITE error term (coherent,
  "structured") is suppressed by randomization while the ground-state population content
  ("structured" itself) is untouched.

---
