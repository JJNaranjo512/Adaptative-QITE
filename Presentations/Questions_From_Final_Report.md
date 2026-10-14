# Questions asked during final report

## 1. Why is that each hardware evaluation of $\dot(\theta)$ costs $O(N^2)$ circuits for N parameters?

Short answer: the metric $M$ is an $N\times N$ matrix, and on hardware each entry needs its own circuit.

Why:
- $M_{jk}=2\big(\mathrm{Re}\langle\tilde G_j\tilde G_k\rangle-\langle\tilde G_j\rangle\langle\tilde G_k\rangle\big)$ (Eq. 6) is a correlation between two generators.
- Hardware can't hand you $\rho$ or $|\psi\rangle$. It only gives expectation values. So each pair $(j,k)$ needs its own measurement, typically a Hadamard test or a pair of parameter-shift circuits.
- $M$ is symmetric, so you need $N(N+1)/2$ entries. That's $O(N^2)$ circuits.
- $Y_k=-\partial_kE$ needs only one gradient per parameter, so it costs $O(N)$.
- $M$ dominates, so one evaluation of $\dot\theta=M^{+}Y$ costs $O(N^2)$.

That's why the step size sets the bill: every step is one full $O(N^2)$ evaluation, and every retry adds more.

Watch out, someone may catch this: Sec. II.E says one evaluation costs $2N+1$ circuits, which is $O(N)$. These two statements don't contradict each other, but they count different things:
- $2N+1$ is the bookkeeping for your simulation. You prepare $\rho(\theta)$ and $\rho(\theta\pm h e_k)$ (that's $2N+1$ states) and then build all of $M$ and $Y$ classically from the full density matrices.
- $O(N^2)$ is the cost on hardware, where you can't read out $\rho$. Full tomography grows exponentially with qubit count.

If asked, say: "The circuit counts in my tables are a simulation currency. They compare controllers fairly with each other, but they underestimate the real hardware cost, which is $O(N^2)$ per evaluation." It would be worth adding that sentence to Sec. II.E.

## 2. Why is that the 9-parameter metric is exactly rank 6 of 9? Why does a two qubit state have $2d - 2 = 6$ real degrees of freedom?

Short answer: $M$ measures how the state moves when you nudge the parameters. A two-qubit pure state can only move in 6 independent directions. However many parameters you add, $M$ can't have rank above 6. With 9 parameters, 3 parameter directions change nothing physical. Those are the gauge modes.

Why 6 degrees of freedom:
- A two-qubit state $|\psi\rangle\in\mathbb C^4$ has $d=4$ complex amplitudes, which is $2d=8$ real numbers.
- Normalization $\langle\psi|\psi\rangle=1$ removes 1.
- Global phase $e^{i\phi}|\psi\rangle$ is physically the same state, which removes 1 more. ($\rho=|\psi\rangle\langle\psi|$ ignores it automatically.)
- That leaves $8-2=6$. In general it's $2d-2$, the real dimension of the space of pure states $\mathbb{CP}^{d-1}$.

Why rank 6:
- $M_{kq}=\mathrm{Tr}[\partial_k\rho,\partial_q\rho]$ is a Gram matrix of the 9 tangent vectors $\partial_k\rho$.
- Those 9 vectors all live in a tangent space of dimension 6. At most 6 of them can be independent, so $\operatorname{rank}M\le 6$.
- "Exactly 6" is what you measured (T3): 6 singular values between 1.27 and 0.043, and 3 at $10^{-17}$–$10^{-18}$, i.e. roundoff. So the ansatz does reach all 6 directions.

The one-liner: "9 knobs, 6 physical directions, so 3 knobs are redundant. $M$ is singular by construction. Any finite $\mathrm{cond}(M)$ for this ansatz, including the paper's 20.22, is really a statement about the pseudo-inverse cutoff."

Likely follow-ups:
- "Is it 6 everywhere?" At most 6 everywhere, and 6 generically. At special points a physical mode can go soft too. At the fold (P1), $\sigma_6$ dips to about $10^{-7}\sigma_1$ but doesn't hit zero. On the 5-parameter ansatz the pseudo-inverse keeps only 4 modes on 12 steps, because of the zero-variance CRZ gate.
- "Noise makes $\rho$ mixed, and mixed states have $d^2-1=15$ dimensions. Doesn't that change the rank?" Not under global depolarizing noise: $\partial_k\rho'=q,\partial_k\rho$, so $M'=q^2M$, which has the same rank.
- "Is 'rank 6' the same thing as 'gauge'?" Yes. The 3 null vectors are parameter combinations that leave $\rho$ unchanged. You identified them in P1: $v_7=$ th1−th6, $v_8\approx$ th2−th7, $v_9\approx$ th4.

Wording to check: Sec. III.C says "a two-qubit pure state has $2d-2=6$ real degrees of freedom" without defining $d$. Add "($d=4$, the Hilbert-space dimension)" so nobody reads $d$ as the number of qubits.

## 3. Whats the context around 5-parameter metric and 9-parameter metric?

Short answer: these are the two circuits from the QDP paper (Saxena–Shalabi–Kyaw, Fig. 4, p. 16) for H₂ on two qubits. The project reuses them so its results are directly comparable to that paper. Each ansatz plays a different role: the 5-parameter one is the easy, under-parametrized test, and the 9-parameter one is the over-parametrized, stiff test.

The system: H₂ in STO-3G, Bravyi–Kitaev mapping, reduced to 2 qubits:
$H=g_0+g_1Z_0+g_2Z_1+g_3Z_0Z_1+g_4Y_0Y_1+g_5X_0X_1$ (QDP Eq. 20).

The circuits (from replicate.py:139, read off QDP Fig. 4):

|                                  | 5-parameter (Fig. 4a)                                                                        | 9-parameter (Fig. 4b)                                                                            |
| -------------------------------- | -------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| Gates                            | RX, RX, RY, RY, then CRZ (q0→q1)                                                             | RY, RY, RZ, RZ, CRX (q0→q1), then RY, RY, RZ, RZ                                                 |
| Params vs. 6 physical directions | 5 < 6, under-parametrized                                                                    | 9 > 6, over-parametrized                                                                         |
| Rank of $M$                      | full (5), generically                                                                        | exactly 6, so 3 gauge modes                                                                      |
| Can it reach the ground state?   | No. Its best energy is $E_{\text{floor}}=-1.124600$, which is $2.1\times10^{-2}$ above exact | Essentially yes. QDP calls it "universal" for two qubits. E3 ends $8.9\times10^{-4}$ above exact |
| Conditioning                     | $\mathrm{cond}_F=863.26$ (paper says 66.72)                                                  | $\mathrm{cond}=\infty$ (paper says 20.22)                                                        |

Where each one is used:
- 5-parameter: Target A (Route 3 is 2.2× closer to the floor than Euler at a matched 732 circuits), Targets B and D (noise), E1 and G2 (Adam vs. ratchet), and P2 (tolerance sweeps). It's the clean benchmark for comparing controllers. That's why everything is measured against $E_{\text{floor}}$ and not the exact energy: no controller can get past what the ansatz can represent.
- 9-parameter: Target C (tight vs. loose rcond), E3, and P1 (fold singularity, $\mathrm{cond}\sim10^7$ at $\tau\approx0.2$). It's the stiffness and degeneracy test.

The one-liner: "The 5-parameter circuit has too few knobs, so it's well defined but can't reach the ground state. The 9-parameter circuit has too many, so it reaches the ground state but its metric is singular by construction. Together they cover both ways an ansatz can be the wrong size."

Likely follow-ups:
- "Why does the 5-parameter one miss the ground state?" It has 5 parameters for a 6-dimensional state space and only one entangler. A 5-dimensional family of states generally doesn't contain the exact ground state.
- "Did you get the circuits exactly right?" The paper leaves the initial state, controlled vs. Ising entangler and wire order unstated. All of these were scanned (the 60-row T3 scan), and none reproduces 66.72.

Gap in the report: the report never lists the gates. Sec. II.D only says "the 5-parameter and 9-parameter ansätze of Ref. [QDP]". A one-line gate list there, or the table above, would save you this question on stage.

## Can you demonstrate me how the term $Y_k$ goes from its trace form to being the variation of energy?

Short answer: for a pure state, $\rho^2=\rho$. That one fact does all the work. It kills the $\langle H\rangle$ term and turns the anticommutator into $\mathrm{Tr}[H,\partial_k\rho]$, which is $\partial_kE$.

Start (Eq. 3):
$$Y_k=-\mathrm{Tr}\left[\partial_k\rho\left({H,\rho}-2\langle H\rangle\rho\right)\right]$$

Step 0: one tool. Differentiate $\rho^2=\rho$:
$$\partial_k\rho,\rho+\rho,\partial_k\rho=\partial_k\rho \qquad (\star)$$

Step 1: the $\langle H\rangle$ term vanishes.
$$\mathrm{Tr}[\partial_k\rho,\rho]=\tfrac12,\partial_k\mathrm{Tr}[\rho^2]=\tfrac12,\partial_k\mathrm{Tr}[\rho]=\tfrac12,\partial_k(1)=0$$
The first equality uses cyclicity: $\partial_k\mathrm{Tr}[\rho^2]=2,\mathrm{Tr}[\partial_k\rho,\rho]$.

Step 2: the anticommutator becomes $\partial_kE$.
$$\mathrm{Tr}[\partial_k\rho,{H,\rho}]=\mathrm{Tr}[\partial_k\rho,H\rho]+\mathrm{Tr}[\partial_k\rho,\rho H]$$
Cycle each trace so that $H$ comes first:
$$=\mathrm{Tr}[H,\rho,\partial_k\rho]+\mathrm{Tr}[H,\partial_k\rho,\rho]=\mathrm{Tr}\big[H(\rho,\partial_k\rho+\partial_k\rho,\rho)\big]\overset{(\star)}{=}\mathrm{Tr}[H,\partial_k\rho]$$
$H$ doesn't depend on $\theta$, so
$$\mathrm{Tr}[H,\partial_k\rho]=\partial_k\mathrm{Tr}[H\rho]=\partial_kE$$

Result:
$$\boxed{Y_k=-\partial_kE}$$
You checked this numerically: $|Y+\nabla E|/|Y|=2.3\times10^{-10}$.

Ket version (this is the route the report's text describes):
With $\partial_k\rho=|\partial_k\psi\rangle\langle\psi|+|\psi\rangle\langle\partial_k\psi|$,
$$\mathrm{Tr}[\partial_k\rho{H,\rho}]=2,\mathrm{Re}\langle\partial_k\psi|H|\psi\rangle+E\big(\langle\psi|\partial_k\psi\rangle+\langle\partial_k\psi|\psi\rangle\big)$$
Normalization makes $\langle\psi|\partial_k\psi\rangle$ purely imaginary, so the bracket is $2,\mathrm{Re}\langle\psi|\partial_k\psi\rangle=0$. The first term is exactly $\partial_kE$. Same answer.

The one-liner: "Purity, $\rho^2=\rho$, kills the mean-energy term and collapses the anticommutator. So for pure states McLachlan's $Y$ is just minus the energy gradient, and VarQITE becomes a natural-gradient descent: $\dot\theta=-M^{+}\nabla E$."

Likely follow-ups:
- "Does it hold under noise?" No. Mixed states have $\rho^2\neq\rho$, so both steps fail. That's why the full trace form is kept for the noise analysis, and why F12 (the dropped term) had to be documented and tested.
- "Then how does Appendix B.1 use $\mathrm{Tr}[\partial_k\rho,\rho]=0$ with noise?" It applies the identity to the ideal pure $\rho$ inside $\rho'=q\rho+(1-q)\mathbb 1/d$, not to $\rho'$ itself. So it's legitimate.
- "Why does this matter for the controller?" Because $Y=-\nabla E$, the gate's prediction $\Delta E_p$ is the exact first-order Taylor term. That makes $|\Delta E_a-\Delta E_p|$ a free second-order error estimate (G5, ratio 2).

## How does the central finite differences serve to compute $d_k \rho$?

Short answer: you nudge one parameter up and down by a tiny $h$, simulate the full density matrix both times, and take the slope:
$$\partial_k\rho(\theta)\approx\frac{\rho(\theta+h,e_k)-\rho(\theta-h,e_k)}{2h},\qquad h=10^{-6}$$
You do that for each of the $N$ parameters, then build $M$ and $Y$ from those matrices by pure linear algebra. The code is drho at Replications/practical-limitations-qdp/replicate.py:221.

Pipeline per evaluation of $\dot\theta$:
1. Simulate $\rho(\theta)$ once, and $\rho(\theta\pm h e_k)$ for each $k$. That's $2N+1$ density matrices, which is where the "$2N+1$ circuits" currency comes from.
2. Get $N$ matrices $\partial_k\rho$ (each $4\times4$) from the formula above.
3. Build $M_{kq}=\mathrm{Tr}[\partial_k\rho,\partial_q\rho]$ and $Y_k=-\mathrm{Tr}[\partial_k\rho({H,\rho}-2\langle H\rangle\rho)]$ classically.
4. Solve $\dot\theta=M^{+}Y$ with the logged rcond.

Why central and not one-sided: Taylor-expand both sides:
$$\rho(\theta\pm h e_k)=\rho\pm h,\partial_k\rho+\tfrac{h^2}{2}\partial_k^2\rho\pm\tfrac{h^3}{6}\partial_k^3\rho+\dots$$
Subtracting cancels the even terms, so the error is $O(h^2)$. A one-sided difference has $O(h)$ error.

Why $h=10^{-6}$: two errors pull in opposite directions:
- truncation error $\sim h^2\approx10^{-12}$, which shrinks with $h$;
- roundoff error $\sim\varepsilon_{\text{mach}}/h\approx10^{-16}/10^{-6}=10^{-10}$, which grows as $h$ shrinks.

At $h=10^{-6}$ roundoff dominates at about $10^{-10}$, which is the size of every validation you report:
- the closed-form diagonal $M_{kk}=\tfrac12(1-\mathrm{Tr}[P_k\rho P_k\rho])$ matches to $8.9\times10^{-11}$;
- the exact state form (P1) matches to $2.594\times10^{-10}$ for $M$ and $3.955\times10^{-10}$ for $Y$;
- the per-entry error on $M$ is $1.652\times10^{-10}$.

Why you used it: it works for any gate, including under noise. You just difference the noisy $\rho$. There are no shift rules to derive and no special cases for CRZ or CRX.

The one-liner: "Wiggle each parameter by $\pm10^{-6}$, take the slope of the density matrix, and build $M$ and $Y$ from the slopes. It's accurate to about $10^{-10}$, I validated it against two closed forms, and it works for any gate and any noise model."

Likely follow-ups:
- "Why not parameter-shift?" Parameter-shift is exact and it's what hardware does. But the controlled rotations (CRZ, CRX) have generator eigenvalues ${0,\pm\tfrac12}$, so they need a four-term rule. Finite differences avoid that and still match to $10^{-10}$. This is the report's App. C answer.
- "Can you do this on hardware?" Not as it stands. Step 1 needs the full $\rho$, which hardware doesn't give you. On hardware you'd measure $M$ and $Y$ entry by entry, at $O(N^2)$ cost, as covered in my first answer.
- "Does the $10^{-10}$ error ever matter?" Only after the fold. There the energy gap between finite-difference and exact runs grows from $1.3\times10^{-11}$ to $1.1\times10^{-4}$, but the final error barely moves ($8.868\times10^{-4}$ vs $8.854\times10^{-4}$).
- "Doesn't finite-difference noise fill the 3 gauge directions?" Barely. $M$ is a Gram matrix, so a true null direction only picks up error at second order ($\sim10^{-20}$). The gauge singular values sit at about $10^{-16}$, which is eigensolver roundoff (P1).