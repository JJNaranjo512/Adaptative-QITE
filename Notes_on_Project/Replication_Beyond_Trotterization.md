---
title: "Replication — Beyond trotterization: Variational product formulas for quantum simulation"
tags:
  - replication
  - physics
  - quantum-computation
  - adaptative-qite
paper: "[[QITE_Hand_Off/Notes_on_Project/Beyond_Trotterization_Variational_Product_Formulas_for_Quantum_Simulation|Beyond_Trotterization_Variational_Product_Formulas_for_Quantum_Simulation]]"
code: QITE_Hand_Off/Replications/Beyond_Trot/replicate.py
---

# Replication — Beyond Trotterization

Reproduction of the numerical claims in [[QITE_Hand_Off/Notes_on_Project/Beyond_Trotterization_Variational_Product_Formulas_for_Quantum_Simulation|Assi, Vogl, Kumari & LeBlanc (2026)]],
Phys. Rev. B **113**, 214314, [arXiv:2511.15124](https://arxiv.org/abs/2511.15124),
DOI 10.1103/941r-ql3g.

Code: `QITE_Hand_Off/Replications/Beyond_Trot/replicate.py`, run with

```bash
~/micromamba/envs/qml-julia/bin/python .../replicate.py --target all
```

on the `qml-julia` environment (numpy 2.4.3, scipy 1.17.1, matplotlib 3.10.8). Nothing
installed. `SEED = 0`, and in fact no target consumes the RNG — every input is a fixed grid or a
fixed parameter point — two consecutive `--target all` runs of 468 lines differ only in the
three lines that print wall-clock seconds. `--check` re-asserts the T6 pinned numbers in about a second and exits
nonzero on drift. Related: [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]],
[[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]].

The paper's Data Availability statement offers no code or data, so every number below was
rebuilt from the printed equations.

## Targets

| Target | Source | Paper value | Obtained | Deviation | Verdict |
|---|---|---|---|---|---|
| T1 gate counts, $N=2\ldots12$ | Table I p.8 | $R_x=2N$, $R_z=4N-2$, CNOT $=4N-4$; Ruth $4N$, $6N-3$, $6N-6$ | identical, all $N$ | 0 | reproduced |
| T1 $N=10$ reductions | text p.8 | −50 %, −33 %, −33 % | −50 %, −33.33 %, −33.33 % | $5.6\times10^{-17}$ | reproduced |
| T1 circuit vs `expm`, our angles | Fig. 12 p.12 | exact | $2.8\times10^{-16}$ | — | reproduced |
| T1 circuit vs `expm`, Fig. 12 angles | Fig. 12 p.12 | exact | 1.7286 | — | **not reproduced** |
| T2 Eq. (31) vs product ansatz | Eq. (31) p.7 | exact | $1.3\times10^{-15}$ | — | reproduced |
| T2 $g$ vs Eq. (35) | Eq. (35) p.7 | exact | $2.8\times10^{-14}$ | — | reproduced |
| T2 $b$ vs Eq. (36) | Eq. (36) p.7 | exact | $2.8\times10^{-14}$ | — | reproduced |
| T2 Eqs. (32)/(33) satisfy Eq. (34) | Eqs. (32)–(34) p.7 | 0 | $1.8\times10^{-9}$ (median) | — | reproduced |
| T2 $\mathcal{E}_F$ from the closed form | Eqs. (32)/(33) p.7 | 0 | $4.3\times10^{-16}$ | — | reproduced |
| T2′ Fig. 2(b) error floor | Fig. 2(b) p.6, "of order $10^{-7}$" | $10^{-7}$ | $4.0\times10^{-11}$ at our `rtol` | — | **inconclusive** |
| T3 $\mathcal{E}_F(t{=}1)$, Var. ABAB | Fig. 4(a) p.8 (read-off ±0.001) | 0.0212 | 0.021217 | +0.08 % | reproduced |
| T3 $\mathcal{E}_F(t{=}1)$, Var. BABA | " | 0.0212 | 0.021143 | −0.27 % | reproduced |
| T3 $\mathcal{E}_F(t{=}1)$, Ruth ABABABA | " | 0.0692 | 0.068847 | −0.51 % | reproduced |
| T3 $\mathcal{E}_F(t{=}1)$, Ruth BABABAB | " | 0.0514 | 0.052034 | +1.23 % | reproduced |
| T3′ same, against the prose | text p.8 ("1.5 %" / "over 3 %") | 0.015 / 0.030 | 0.0211 / 0.0688 | +41 % / +129 % | **not reproduced** |
| T3 Ruth/var ratio at $t=1$ | text p.8 ("more than twice") | > 2 | 2.461 | — | reproduced |
| T3″ var. beats Ruth on **all** of $[0,2]$ | text p.8 ("over [the] entire time interval") | true | false | — | **not reproduced** |
| T4 Eq. (61) → Eq. (69), $A$ = Ising | Eqs. (61)/(69) p.13 | exact | $3.6\times10^{-16}$ | — | reproduced |
| T4 Eq. (61) → Eq. (69), Eq. (54) labels | " | exact | 2.2334 | — | **not reproduced** |
| T4 Fig. 13 TS at $t=100$, raw $\lVert\cdot\rVert_F$ | Fig. 13 p.14 (read-off ±0.001) | 0.2694 | 0.26951 | +0.04 % | reproduced |
| T4 Fig. 13 Eq. (69) at $t=100$, raw | " | 0.1715 | 0.171821 | +0.19 % | reproduced |
| T4 same two, Eq. (27)-normalized | " | 0.2694 / 0.1715 | 0.023822 / 0.015187 | −91 % | **not reproduced** |
| T5 $\Xi$ two-level | text p.13 | −1 | −1 | 0 | reproduced |
| T5 Eq. (65) = Eq. (68), Eq. (54) labels | Eqs. (65)/(68) p.13 | 15/15 | 0/15 | — | **not reproduced** |
| T5 Eq. (65) = Eq. (68), $A\leftrightarrow B$ | " | 15/15 | 15/15 | 0 | reproduced |
| T5 $\Xi$ predicts the better ordering | Eq. (68) p.13 | 7/7 | 0/7 literal, 7/7 relabelled | — | reproduced **only relabelled** |
| T6 pinned cond($g$) table | vault note §5 item 2 (**ours**) | see below | all reproduce | ≤ 0.3 % | reproduced |
| T7 $\mathcal{E}_F^{\rm Strob}(t{=}100)$ TS | text p.11 (prose, ±0.5 pp) | 0.11 | 0.108968 | −0.94 % | reproduced |
| T7 $\mathcal{E}_F^{\rm Strob}(t{=}100)$ var. | " | 0.055 | 0.055305 | +0.55 % | reproduced |
| T7 TS/var ratio | text p.11 ("twofold") | 2 | 1.9703 | −1.48 % | reproduced |
| T8 median rel-err ratio, $\tau=0.2$ | text p.8 ("factor of 4 or 5") | 4–5 | 4.213 | in band | reproduced |
| T8 median rel-err ratio, $\tau=0.4$ | " | 4–5 | 4.762 | in band | reproduced |
| T8 advantage grows with $\tau$ | text p.8 | true | 1 of 3 statistics | — | **partially reproduced** |

Figures in `Replications/Beyond_Trot/figures/`: `table1_gate_counts.png`, `fig2a_two_level.png`,
`fig2b_two_level.png`, `fig2b_error_floor.png`, `fig4a_qim_error.png`, `fig4b_qim_strob.png`,
`fig3_magnetization.png`, `fig11b_xxz_strob.png`, `fig13_approx_vs_ts.png`, `t4_tau_sweep.png`,
`fig_xi_phase_diagram.png`, `cond_g_trajectory.png`.

## Assumptions

All of these are **my** choices, not the paper's, unless marked otherwise.

1. **Integration starts at $t_0 = 10^{-6}$** with the Eq. (10) linear (Trotter–Suzuki) initial
   condition, because $g$ is exactly rank-deficient at $t = 0$ for any ansatz that repeats an
   operator — proved numerically in T6. The paper is silent on this.
2. **$g\dot c = -b$ is solved with `numpy.linalg.lstsq` at `rcond` $=10^{-10}$** (minimum-norm
   solution). The paper states no regularization at all.
3. **`scipy.solve_ivp(method="RK45")` with `rtol` $=10^{-10}$, `atol` $=10^{-12}$.** The paper
   names RK45 (p.7) but gives no tolerances. This turns out to decide T2′ entirely.
4. **Trace normalization** $g_{jk} = \mathrm{Tr}[\cdot]/D$, $b_j = i\,\mathrm{Tr}[\cdot]/D$, i.e.
   $b = +iF$ in Eq. (3). Not really a free choice — it is *verified* against the paper's own
   printed Eqs. (35)/(36); see the convention gate below.
5. **Branch tracking of Eqs. (32)/(33).** As printed, `arctan`/`arcsin` give a discontinuous
   $c(t)$; continuity requires `atan2` unwrapped.
6. **Gate convention** $R_{x,z}(\theta) = e^{-i\theta\sigma/2}$.
7. **$A\leftrightarrow B$ are interchanged in Eqs. (61), (65), (68), (69)** relative to the
   Eq. (37)/(54) splitting. This is a *determination*, numerically forced in T4/T5, not a
   preference — see Discrepancy 3.
8. **Eq. (37) is an open chain** ($j = 1\ldots N-1$), confirmed on p.7 and by the $(N-1)$ in
   Eqs. (68)/(69).
9. **Read-off bands**: ±0.2 pp on the T3 prose, ±0.5 pp on the T7 prose. The pixel calibrations
   of Figs. 4(a) and 13 are readable to ±0.001, but the band actually applied when scoring them
   is looser — 0.002 for Fig. 4(a), 0.005 for Fig. 13 — to leave room for the calibration itself.
   Every T4 deviation is $\le 3.3\times10^{-4}$, inside the tighter figure either way.
10. **Not mine, but worth recording:** Eq. (27) is printed as
    $\mathcal{E}_F = \lVert U-U_a\rVert_F/(2\sqrt{D})$ — $2\sqrt{D}$, *not* $\sqrt{2D}$. That is
    the form that delivers the $0\le\mathcal{E}_F\le1$ the paper claims, since
    $\lVert U-U_a\rVert_F^2 \le 4D$ for unitaries. The earlier transcription in the paper note
    had this wrong and has been corrected.

### The convention gate

Eq. (3) is complex and varied over real $c_j$, so some Re/Im choice must be made and the paper
never states it. T2 therefore runs *before* any accuracy number is quoted: it builds $g$ and $b$
from the literal trace definitions and asserts them against the printed Eqs. (35)/(36) over a
343-point $(c_0,c_1,c_2)$ grid. Result:

| convention | max $\lvert g - $Eq. (35)$\rvert$ | max $\lvert b - $Eq. (36)$\rvert$ |
|---|---|---|
| $\mathrm{Tr}[\cdot]/D$, $b = +iF$ | $2.8\times10^{-14}$ | $2.1\times10^{-14}$ |
| $\mathrm{Tr}[\cdot]/D$, $b = -iF$ | $2.8\times10^{-14}$ | $5.0\times10^{1}$ |
| $\mathrm{Tr}[\cdot]$, no $1/D$ | $2.5\times10^{1}$ | $2.5\times10^{1}$ |
| $2\,\mathrm{Re}\,\mathrm{Tr}[\cdot]/D$ | $2.5\times10^{1}$ | $2.5\times10^{1}$ |

Only one candidate survives, and the measured imaginary parts are $4\times10^{-16}$, so the
Re/Im question is vacuous for this ansatz class rather than a real choice. Everything downstream
is therefore a disagreement with the paper, not a suspected convention error on my side — the
same role `check_eq_b9()` plays in [[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]].

## Discrepancies

### 1. Fig. 13 is missing its own normalization (T4)

Fig. 13's axis is labelled $\mathcal{E}_F^{\rm Strob}$, but the plotted values match the **raw**
$\lVert U(n\tau) - U_a(\tau)^n\rVert_F$ to four significant figures — 0.26951 against a read-off
of 0.2694, and 0.17182 against 0.1715, ratios of 1.0004 and 1.0019. Applying Eq. (27) as printed
gives 0.023822 and 0.015187, low by exactly $2\sqrt{D} = 11.3137$ at $D = 2^5$.

Fig. 11(b), at the **same** $D = 32$, *does* carry the normalization and reproduces to under 1 %
(T7). So the two figures are inconsistent with each other, and only one of them can be showing
what its axis label says. The physics conclusion of Fig. 13 — Eq. (69) beats symmetric TS at
equal gate count — is a ratio and survives untouched.

### 2. The p.8 text disagrees with the p.8 figure (T3, T3′)

The prose says ≈ 1.5 % for the variational ansatz and "over 3 %" for Ruth at $t = 1$.
Pixel-calibrating the paper's own Fig. 4(a) on the same page gives 0.0212 (both variational
orderings), 0.0692 (Ruth ABABABA) and 0.0514 (Ruth BABABAB). I rendered page 8 at 300 dpi and
confirmed by eye that the variational curves sit at ~0.02 at $t=1$, that the two overlap, and
that Ruth BABABAB lies below Ruth ABABABA — all four qualitative features hold.

This replication lands on the **figure** to within 0.08–1.23 % and therefore misses the
**prose** by +41 % (variational) and +73 % to +129 % (Ruth). The qualitative claim the prose is
making — Ruth's error is "more than twice" the variational one — survives in both readings; the
obtained ratio is 2.461.

### 3. Eqs. (61), (65), (68), (69) use $A$ and $B$ swapped (T4, T5)

Under the Eq. (37)/(54) labels ($A$ = transverse block $h_x\sum S^x$, $B$ = Ising block), the
printed formulas simply do not follow:

| check | Eq. (54) labels | $A$ = Ising block |
|---|---|---|
| Eq. (65) reproduces Eq. (68) | 0/15 points | 15/15 points |
| Eq. (61) reproduces Eq. (69) | off by 2.2334 | $3.6\times10^{-16}$ |
| $\Xi$ predicts the empirically better ordering | 0/7 points | 7/7 points |

One consistent relabelling fixes all three at once, so this is systematic rather than a typo.
It is also the reading consistent with the paper's own choices elsewhere: at $h_x=h_z=J=1$ it
gives Ising-first, i.e. BAB — exactly the ordering Fig. 3 uses and calls best — and for the
two-level model, where the paper names $A = h_x\sigma_x$ explicitly, it gives AB, matching
Fig. 2(a). Since 0/15 and 0/7 are not near-misses, this is a determination rather than a guess.

### 4. The claim that the variational ansatz wins over *all* of $[0,2]$ is false (T3″)

p.8: "The variational ansatz achieved lower error than the Ruth ansatz over the entire time
interval $[0,2]$." It does not, and the reason is structural rather than numerical. Ruth is a
fourth-order formula, so its error is $O(t^5)$; the four-exponential ansatz reduces to
Trotter–Suzuki as $t\to0$ by Eq. (10), so its error is $O(t^3)$. Ruth must therefore win near
the origin. Fitted slopes over the fixed window $t\in[10^{-2},10^{-1}]$ at $N=10$:
$\mathcal{E}_F^{\rm var} \sim t^{2.999}$, $\mathcal{E}_F^{\rm Ruth} \sim t^{4.999}$.

The 60-point log grid brackets the crossover to $0.375 < t < 0.426$, and interpolating the zero
of $\log(\mathcal{E}_F^{\rm var}/\mathcal{E}_F^{\rm Ruth})$ inside that bracket gives
$t = 0.376$ at $N=10$. Below the crossover Ruth is the more accurate of the two,
by up to five orders of magnitude at $t = 10^{-3}$. Fig. 4(a)'s linear axes cannot resolve this
— every curve is visually pinned to zero below $t\approx0.5$ — which is presumably how the claim
survived. The paper's *useful* claim (variational wins in the regime where either error is large
enough to matter, at a third of the gate count) is untouched; the "entire interval" phrasing is
simply too strong.

**Two side-checks behind this section are not pinned.** At $N=5$, where the solve is cheap, I
checked in-session that (i) tightening `rtol` from $10^{-10}$ to $10^{-12}$ leaves both fitted
slopes unchanged to three digits, so they are not integrator noise, and (ii) a 400-point scan of
the bracket reproduces the interpolated crossover (0.381 at $N=5$) to the fourth digit, so the
interpolation is not doing the work. Both were run in a scratch session whose outputs were never
saved, and `replicate.py` runs T3 only at $N=10$ and only at the frozen `rtol` — so these two
statements are testimony, not artifacts, and nothing in the committed script would catch it if
they stopped being true. Folding a cheap `--target t3-n5` variant into the script is the fix.

### 5. Fig. 2(b)'s $10^{-7}$ floor is unfalsifiable as printed (T2′)

The *mechanism* the paper gives is right and reproduces cleanly: the 3-exponential ansatz is
**exact** for the two-level model — substituting the closed forms of Eqs. (32)/(33) gives
$\mathcal{E}_F \le 4.3\times10^{-16}$ — so everything visible in Fig. 2(b) is integrator noise,
just as p.7 says. But the *height* of that floor is set by the RK45 tolerance, which the paper
never states:

| `rtol` / `atol` | $10^{-3}/10^{-6}$ | $10^{-4}/10^{-7}$ | $10^{-5}/10^{-8}$ | $10^{-6}/10^{-9}$ | $10^{-8}/10^{-10}$ | $10^{-10}/10^{-12}$ | $10^{-12}/10^{-14}$ |
|---|---|---|---|---|---|---|---|
| max $\mathcal{E}_F$ | $9.9\times10^{-4}$ | $1.4\times10^{-4}$ | $1.4\times10^{-5}$ | $1.2\times10^{-6}$ | $6.9\times10^{-9}$ | $4.0\times10^{-11}$ | $3.3\times10^{-13}$ |

A floor of $10^{-7}$ brackets to `rtol` between $10^{-8}$ and $10^{-6}$. At the tolerances I
froze before the first accuracy run I get $4.0\times10^{-11}$, four orders below theirs. So the
number $10^{-7}$ measures the authors' integrator settings and nothing about the method — it
cannot be scored either way, which is why T2′ is *inconclusive* rather than *not reproduced*.

### 6. Fig. 12's gate angles need $S = \sigma$, but Eq. (37) says $S = \sigma/2$ (T1)

Building one layer of Eq. (18) at $N=4$ and comparing against the matrix exponential: with
$\alpha=-ch_z$, $\beta=-cJ/2$, $\gamma=-ch_x$ the circuit matches to $2.8\times10^{-16}$; with
Fig. 12's printed $\alpha=-2ch_z$, $\beta=-2cJ$, $\gamma=-2ch_x$ it misses by 1.7286. Gate
*counts*, which is all Table I claims, are unaffected — T1 reproduces exactly for every
$N = 2\ldots12$.

### 7. "The advantage grows with $\tau$" depends on which statistic you pick (T8)

The pointwise ratio of relative errors has no finite mean, because both errors cross zero. Three
statistics:

| statistic | $\tau=0.2$ | $\tau=0.4$ | grows? |
|---|---|---|---|
| median of the ratio | 4.213 | 4.762 | yes |
| ratio of the mean errors | 4.147 | 3.995 | no |
| ratio of the max errors | 3.343 | 3.315 | no |

The "factor of 4 or 5" reproduces on the robust statistic — the median — and it is the only one
in which both values land inside the paper's band. But the accompanying claim that the advantage
becomes more pronounced at larger $\tau$ holds for 1 of 3 statistics, so it is not robust.

### New result — cond($g$) is the paper's unstated regularization (T6)

The paper reports cond($g$) nowhere, yet its whole method inverts $g$. T6 pins what the vault
note §5 previously asserted only in prose. For the two-level Eq. (35) metric, which depends on
$c_1$ alone:

| $c_1$ | $\det g$ | rank | cond($g$) |
|---|---|---|---|
| 0 | 0 | **2 of 3** | $2.2\times10^{17}$ |
| $10^{-6}$ | $4.0\times10^{-8}$ | 3 | $2.50\times10^{11}$ |
| $10^{-3}$ | $4.0\times10^{-2}$ | 3 | $2.50\times10^{5}$ |
| 0.2 | $1.2865\times10^{3}$ | 3 | 10.60 |

The null vector at $c_1=0$ is $\propto(1,0,-1)$ — i.e. at $t=0$ the two $A$-exponentials of an
ABA ansatz are indistinguishable, which is exactly the degeneracy that makes $t_0>0$ mandatory.
cond($g$) $\sim c_1^{-2}$ over $c_1\in[10^{-6},10^{-3}]$, fitted exponent $-2.0000$. Along the
QIM trajectory at $N=5$: cond($g$) $=3.4\times10^{12}$ (BAB) and $2.1\times10^{16}$ (ABAB) at
$t=t_0$, falling to 1.36 and 4.85 by $t=2$.

Encouragingly, the *result* does not depend on the cutoff: varying $t_0$ over
$[10^{-8},10^{-3}]$ and `rcond` over $[10^{-14},10^{-4}]$ moves $\mathcal{E}_F(t{=}2)$ by at most
$5.12\times10^{-10}$ relative. The
ill-conditioning is real but confined to a boundary layer near $t=0$ that the dynamics leaves
behind.

`--check` re-asserts all of this in about a second and exits nonzero on drift. Worth recording
how the gate had to be built, because the obvious version did not work: comparing against the
§5 table's three-significant-figure values forces a ~2 % band, and a 2 % band is wide enough to
hide a real regression — a deliberate 1.55 % perturbation of $\det g$ at $c_1=0.2$ passed
silently. The gate now pins the **full-precision** computed values at a $10^{-9}$ relative
margin (floating-point reproducibility, not a fit knob) and checks agreement with the note's
rounded table separately at $5\times10^{-3}$.

The two directions of that verification rest on different evidence, and the difference is worth
recording. That the **old** gate passed a 1.55 % perturbation silently is an artifact: the run
exists, against the pre-fix script, printing "0 drifted value(s)". That the **new** gate exits 1
on the same perturbation — and on $8\times10^{-6}$, but not on $8\times10^{-12}$ — I verified
interactively, without capturing the output, so it is testimony. The margin it
rests on, `T6_RTOL` $=10^{-9}$, is in the script and inspectable; the demonstration that the gate
fires is not. Anyone relying on it should perturb a pinned value once and watch it fail.

## Verdict

Per target, not global:

- **T1 — reproduced.** Table I is exact for every $N=2\ldots12$, emitted from an actual gate list
  rather than a formula, and the assembled circuit matches `expm` to $2.8\times10^{-16}$. The
  only casualty is Fig. 12's printed angles, which need the opposite spin convention.
- **T2 — reproduced.** Eqs. (31)–(36) all check out to $\sim10^{-14}$, and Eqs. (32)/(33) are an
  exact solution of the equation of motion. This is the paper at its strongest.
- **T2′ — inconclusive.** The mechanism behind the Fig. 2(b) floor is confirmed; its stated value
  is a function of an unreported integrator tolerance and cannot be checked.
- **T3 — reproduced against the figure, not reproduced against the text.** Every curve lands on
  Fig. 4(a) to ~1 %; every curve misses the prose by 41–129 %. The paper is internally
  inconsistent on its own page 8.
- **T3″ — not reproduced.** The "entire interval $[0,2]$" claim fails at small $t$ for an
  order-counting reason: Ruth is $O(t^5)$, the ansatz is $O(t^3)$.
- **T4 — reproduced, but only against the unnormalized norm.** Eq. (69) follows from Eq. (61)
  exactly ($3.6\times10^{-16}$) under the $A\leftrightarrow B$ relabelling, and beats symmetric
  TS at equal gate count. Fig. 13's absolute values reproduce only if Eq. (27) is dropped.
- **T5 — reproduced only under a consistent relabelling.** 0/15 and 0/7 as printed; 15/15 and
  7/7 with $A$ meaning the Ising block throughout.
- **T6 — reproduced.** All pinned values, now backed by a committed script rather than prose.
- **T7 — reproduced.** −0.94 % and +0.55 % against prose values with a ±0.5 pp band, and the
  "twofold improvement" comes out at 1.9703.
- **T8 — partially reproduced.** The "factor of 4 or 5" holds on the median; the claim that it
  grows with $\tau$ is statistic-dependent.

Nothing was tuned. `T0`, `RCOND`, `RTOL` and `ATOL` were frozen before the first accuracy run and
T6 shows the results are insensitive to all of them.

### What would settle the open items

- **Fig. 13 vs Fig. 11(b):** whether Fig. 13 omitted the $1/(2\sqrt{D})$ or Fig. 11(b) applied it
  twice. Only the authors can say; the ratio is unaffected either way.
- **The $A\leftrightarrow B$ swap:** whether Eqs. (61)–(69) were derived with a different
  splitting convention than Eq. (54) declares. The evidence (0/15 → 15/15) is one-sided enough
  that I would treat the relabelled reading as correct and the paper as mislabelled.
- **T2′:** the RK45 `rtol`/`atol` used for Fig. 2(b). One sentence in the paper would fix it.
- **T3 prose vs figure:** which of 1.5 % and 0.0212 the authors intended. Note that 1.5 % is
  roughly what the *stroboscopic* error of Fig. 4(b) reaches early on, so a panel mix-up is a
  plausible explanation — but that is speculation and I have not tested it.

### Why this matters for the adaptive-step work

Two things carry over to the QITE side, both about the metric.

The paper's $g_{jk} = \mathrm{Tr}[(\partial_j U_a)^\dagger \partial_k U_a]$ is the same object as
the project's $M$ in $M\dot\theta = Y$ — a Gram matrix of parameter derivatives — and it is
rank-deficient for the same structural reason: whenever an ansatz repeats a generator, the
derivatives coincide at $t=0$. T6 makes that explicit (rank 2 of 3, null vector $(1,0,-1)$), so
the "cond($M$) wall" of [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §4 is not a QITE artifact; it
is generic to variational product formulas, real and imaginary time alike.

But the encouraging half is the one worth carrying: here the ill-conditioning is a **boundary
layer**, not a persistent obstacle. cond($g$) falls from $10^{12}$–$10^{16}$ at $t_0$ to order 1
by $t=2$, and the final answer is insensitive to `rcond` over ten decades. That is the opposite
of what [[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]] found for the nine-parameter QDP ansatz, where
$M$ is rank-deficient *by construction* at every point of the manifold and the pseudo-inverse
cutoff is a genuine free parameter. The distinction — transient degeneracy at an initialization
point versus a structurally singular metric — is one an adaptive-step controller should be able
to tell apart, and it suggests the right response to a stiff first step is to start at $t_0>0$
from the Trotter–Suzuki value rather than to regularize harder.

Eq. (69) is also still worth stealing on its own terms: closed-form cubic corrections to the TS
angles, zero extra gates, and T4 confirms it beats symmetric TS at every $\tau$ swept up to 1.0
(error ratio 0.63, remarkably flat).

How to act on all this in Route 3, Adam_QITE and the memo's routes (a plan, nothing run yet):
[[QITE_Hand_Off/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]].
