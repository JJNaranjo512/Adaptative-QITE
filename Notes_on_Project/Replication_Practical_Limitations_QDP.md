---
title: "Replication — Practical limitations of quantum data propagation on noisy quantum processors"
tags:
  - replication
  - physics
  - quantum-computation
  - adaptative-qite
paper: "[[QITE_Hand_Off/Notes_on_Project/Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors|Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors]]"
code: QITE_Hand_Off/Replications/practical-limitations-qdp/replicate.py
---

# Replication — Practical Limitations of Quantum Data Propagation

Reproduction of the numerical claims in [[QITE_Hand_Off/Notes_on_Project/Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors|Saxena, Shalabi & Kyaw (2024)]],
Phys. Rev. Applied **21**, 054014, [arXiv:2306.13046v2](https://arxiv.org/abs/2306.13046v2).

Code: `QITE_Hand_Off/Replications/practical-limitations-qdp/replicate.py`, run with

```bash
~/micromamba/envs/qml-julia/bin/python .../replicate.py --target all
```

on the `qml-julia` environment (numpy, scipy, matplotlib, pennylane 0.45.1). Nothing installed.
`SEED = 0`; everything reported below is deterministic. Related: [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]].

## Targets

| Target | Source | Paper value | Obtained | Deviation | Verdict |
|---|---|---|---|---|---|
| T1 Thm 2 law $\frac{1-(1-p)^N}{(1-p)^N}$ | Eq. (19) p.4; Fig. 1(b) p.2; Fig. 5 p.16 | closed form | identical | 0 | reproduced |
| T2 $p_{max}$, $N=5$ | Fig. 3 caption p.5 | 0.032 | 0.03311 | +3.48 % | partially reproduced |
| T2 $p_{max}$, $N=6$ | " | 0.036 | 0.03653 | +1.48 % | partially reproduced |
| T2 $p_{max}$, $N=10$ | " | 0.052 | 0.05349 | +2.86 % | partially reproduced |
| T2 $p_{max}$, $N=12$ | " | 0.068 | 0.06964 | +2.41 % | partially reproduced |
| T2 $p_{max}$, $N=14$ | " | 0.129 | 0.13083 | +1.42 % | partially reproduced |
| T3 $\|M\|_F$, 5-param | Appendix C p.16 | 0.9977 | 0.8793 | −11.9 % | not reproduced |
| T3 cond($M$), 5-param | " | 66.7239 | 863.2577 | +1194 % | **not reproduced** |
| T3 $\|M\|_F$, 9-param | " | 1.689 | 1.6387 | −2.98 % | partially reproduced |
| T3 cond($M$), 9-param | " | 20.22 | $\infty$ (rank 6 of 9) | — | **not reproduced** |
| T3 $\|Y\|$, 9-param | " | 0.121 | 0.1927 | +59.3 % | not reproduced |
| T4 $M^{err} = (1-p)^{2N}M$ | Eq. (B12) p.15 | exact | $\le 4\times10^{-10}$ rel. | ~1e-10 | reproduced |
| T4 $Y^{err} = (1-p)^{N}Y$ | Eq. (B13) p.15 | exact | $\le 1\times10^{-9}$ rel. | ~1e-9 | reproduced |
| T4 $(1-p)^N\dot\theta^{err} = \dot\theta$ | Eq. (B16) p.15 | Eq. (19) | on the curve | max 2.8e-7 | reproduced |
| T5 Fig. 3(b) endpoints | Fig. 3(b) p.5 (read off ±10–15) | 75, 105, 300, 430, 590 | 77.2, 109.7, 295.7, 423.9, 581.7 | ≤ 8.3 abs | partially reproduced |

Figures in `Replications/practical-limitations-qdp/figures/`: `fig1b_theorem2.png`,
`fig3a_allowed_p.png`, `fig3b_theorem1.png`, `fig6b_theorem1_9param.png`, `t4_scaling.png`.

## Assumptions

All of these are **my** choices, not the paper's. The paper does not state any of them.

1. **Initial state $|00\rangle$.** $\rho_0$ is never given. T3 scans $|00\rangle$, $|01\rangle$,
   $|10\rangle$, $|11\rangle$, $|{+}{+}\rangle$; none rescues the reported cond.
2. **Parameter-to-gate assignment column-wise as drawn in Fig. 4** ($\theta_1$ top-left,
   $\theta_2$ below it, …). Another assignment only permutes $M$, so it changes neither
   $\|M\|_F$ nor cond($M$) — this one is harmless.
3. **Central-difference derivatives**, $h = 10^{-6}$, instead of parameter-shift. Accuracy is
   asserted in the code against the closed form of Eq. (B9), agreement ~1e-10.
4. **cond$(\cdot)$ is the Frobenius condition number** $\|M\|_F\|M^{-1}\|_F$, which is what the
   paper defines in Eq. (A39). The spectral cond is printed alongside (624.9 for the 5-param).
5. **The global depolarizing channel is $\mathcal{E}(\rho) = (1-p)\rho + p\,\mathbb{1}/d$ on the
   whole two-qubit register** ($d=4$), applied after every gate, per Eq. (15) and Fig. 2.
   Implemented as an explicit `qml.QubitChannel` with Pauli-twirl Kraus operators.
6. **$N$ = number of variational parameters = number of gates**, so it is also the number of
   noise insertions.
7. **$\|Y\|$ for the five-parameter ansatz is computed by me** (0.3913) because the paper never
   reports it, even though Fig. 3(b) depends on it. See the T5 discrepancy.
8. **Fig. 3(b) endpoints are eyeballed off the plot** to ±10–15 on a 0–600 axis.

## Discrepancies

### T2 — $p_{max}$ is 1.4–3.5 % high (partially reproduced)

Same Eq. (18), same $\delta = 0.04$, same $\|M\| = 0.9977$, and the deviation is one-sided:
every computed value is *above* the reported one, by 0.0005–0.0018 absolute. Inverting Eq. (18)
for the $\|M\|$ each reported value would imply gives 1.1018, 1.0456, 1.0947, 1.0631, 1.0075 —
mutually inconsistent, so a mis-stated $\|M\|$ cannot be the explanation. The pattern (always
low, by less than one grid step) is what you get from quoting the last point of a discrete $p$
grid of spacing ~0.001–0.002 rather than solving Eq. (18) exactly. Not resolvable from the paper.

### T3 — cond($M$) does not reproduce (not reproduced)

This is the substantive failure, and it matters because cond($M$) is the whole mechanism of the
paper's argument.

**Five-parameter ansatz.** Under Eqs. (6)–(8) at the published angles I get cond$_F$ = 863.26
against the reported 66.7239, and $\|M\|_F$ = 0.8793 against 0.9977. cond is scale-free, so no
normalisation convention can move 863 to 66.7. I ran a 60-row scan over the ambiguities the
paper leaves open — $\rho_0 \in \{|00\rangle,|01\rangle,|10\rangle,|11\rangle,|{+}{+}\rangle\}$
× two-qubit gate ∈ {controlled, Ising} × wire order ∈ {as drawn, reversed} × metric ∈
{$\mathrm{Tr}[\partial\rho\partial\rho]$, Fubini–Study, 4·FS}. The cond values cluster at 863.26,
879.59 and 911.08. **No combination lands within 5 % of 66.7239.**

**Nine-parameter ansatz.** $M$ is exactly rank 6, singular values
$[1.267,\ 0.706,\ 0.633,\ 0.383,\ 0.181,\ 0.0432,\ 8.8\!\times\!10^{-17},\ 2.9\!\times\!10^{-17},\ 3.2\!\times\!10^{-18}]$.
This is not numerical noise — a two-qubit *pure* state manifold has $2d-2 = 6$ real dimensions,
so a $9\times 9$ metric built on it **cannot** have full rank. cond($M$) is $\infty$ by
construction, and the reported 20.22 requires a pseudo-inverse cutoff the paper never states.
Sweeping the cutoff: rcond $\in [10^{-15},\,10^{-2}]$ all give 39.41, rcond = 0.1 gives 10.69,
rcond = 0.3 gives 5.66. cond$_F$ is monotone in the cutoff, so *some* cutoff reproduces any
target value — which is exactly why 20.22 as printed is unfalsifiable.

$\|Y\|$ = 0.1927 against 0.121 (+59 %) at the same angles.

The $M$ construction itself is not in doubt: `check_eq_b9()` verifies the diagonal against the
paper's own closed form $M_{kk} = \tfrac12(1-\mathrm{Tr}[P_k\rho_{k-1}P_k\rho_{k-1}])$ (Eq. B9)
to 1e-10 for all four single-qubit gates. So this is a disagreement with Appendix C, not a
convention error on my side — though I cannot rule out that the paper used angles, an ordering
or an initial state different from what Appendix C and Fig. 4 print.

### T5 — Eq. (17) reproduces Fig. 3(b) only if the $\|Y\|$ term is dropped

Four readings evaluated side by side at $p = p_{max}(N)$, against the read-offs
75 / 105 / 300 / 430 / 590:

| variant | N=5 | N=6 | N=10 | N=12 | N=14 |
|---|---|---|---|---|---|
| Eq. 17 as printed, $\|Y\|$ = mine (0.3913) | 103.4 | 146.6 | 396.6 | 575.3 | 824.3 |
| Eq. 17 with $\sqrt{N}/(\sqrt2\|Y\|)$ as *derived* from Eq. (A37) | 119.0 | 168.8 | 457.0 | 665.9 | 969.5 |
| **Eq. 17, $\|Y\|$ term negligible** | **77.2** | **109.7** | **295.7** | **423.9** | **581.7** |
| Eq. A42 with its denominator | invalid | invalid | invalid | invalid | invalid |

Only the third row lands inside the read-off tolerance for every $N$ (worst deviation 8.3).
Since the paper never reports $\|Y\|$ for the five-parameter ansatz, I can't tell whether they
used a large $\|Y\|$ or omitted the term; either way Fig. 3(b) is not reproducible from the
printed inputs alone.

Two paper-internal inconsistencies, evaluated rather than silently resolved:

- Eq. (17) drops the $1 - \mathrm{cond}(M)\|\Delta M\|/\|M\|$ denominator that Eq. (A42)
  carries. Using Eq. (A43), that denominator stays positive only for $p \le 4.28\times10^{-4}$
  ($N=5$) down to $6.68\times10^{-5}$ ($N=14$), while Fig. 3(b) is plotted out to
  $p = 0.032$–$0.129$. So the main-text form is the only one consistent with the plotted range,
  and the appendix bound is invalid over essentially all of it.
- Eq. (A38)'s $\sqrt{N}/\sqrt{2\|Y\|}$ does not follow from Eq. (A37), which gives
  $\sqrt{N}/(\sqrt2\,\|Y\|)$. Both are in the table above; neither matches.

### New result — Theorem 2 is a factor of two pessimistic on hardware (T4 extension)

Appendix B applies noise only to the derivative line and leaves the $\rho$ inside $\{H,\rho\}$
ideal. On hardware both come off the same noisy device. Recomputing with both noisy, the
relative error is **exactly half** of Eq. (19) — every ratio 0.5 to $3.8\times10^{-5}$ across
$p \in [10^{-4}, 0.2]$.

The reason is analytic. For a normalised state $2\,\mathrm{Re}\langle\partial_k\psi|\psi\rangle = 0$,
so Eq. (8) collapses to $Y_k = -\mathrm{Tr}[\partial_k\rho\,H] \equiv -g_k$: $Y$ is minus the
energy gradient, verified numerically as $\|Y+g\|/\|Y\| = 2.3\times10^{-10}$, hence
$M^{-1}g = -\dot\theta$. Writing $\rho^{err} = (1-p)^N\rho + (1-(1-p)^N)\mathbb{1}/d$ with
$\partial_k\rho^{err} = (1-p)^N\partial_k\rho$, the maximally-mixed piece contributes
$-(1-p)^N(1-(1-p)^N)g_k/2$ to $Y$, giving

$$\dot\theta^{err} = \dot\theta\left[1 + \frac{1-(1-p)^N}{2(1-p)^N}\right].$$

Theorem 2 isn't wrong — it answers the question its own Fig. 2 noise placement poses. But the
error a real QITE step incurs under global depolarizing noise is half what Theorem 2 states.

## Verdict

Per target, not global:

- **T1 — reproduced.** Eq. (19) is a closed form with no free parameters; it is an identity, and
  T4 confirms it independently from a noisy circuit simulation.
- **T2 — partially reproduced.** Correct to 1.4–3.5 %, all one-sided; consistent with the paper
  quoting a discrete $p$ grid, but not resolvable from the printed inputs.
- **T3 — not reproduced.** cond($M$) = 863.26 vs 66.7239 for the five-parameter ansatz, with no
  convention in a 60-row scan getting within 5 %; and the nine-parameter $M$ is rank-deficient
  by construction, so its cond is $\infty$ and 20.22 depends on an unstated pseudo-inverse
  cutoff. The $M$ construction is validated against the paper's own Eq. (B9).
- **T4 — reproduced**, to ~1e-10 from a genuine `default.mixed` simulation with a global
  depolarizing channel after every gate, plus the factor-of-½ extension above.
- **T5 — partially reproduced.** Fig. 3(b) is recovered only under the "$\|Y\|$ term negligible"
  reading; the two readings that use a finite $\|Y\|$ overshoot by 35–65 %, and the appendix form
  of the bound is invalid across the plotted range.

### What would settle T3

The angles, the initial state, the gate set and the parameter ordering used for Appendix C, and
— for the nine-parameter case — the pseudo-inverse cutoff. Worth asking Gaurav directly rather
than scanning further; the scan space is exhausted.

### Why this matters for the adaptive-step work

If $M$ is rank-deficient at the ansätze we care about, the QITE step is being taken through a
pseudo-inverse whose cutoff is a free parameter nobody has pinned down. An adaptive
time-stepping controller sitting on top of that is tuning the step size against a $\dot\theta$
whose own regularisation is unspecified — the cutoff has to become an explicit part of the
method, not an implicit `np.linalg.pinv` default. See [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]].
