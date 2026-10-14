# Beyond_Trot

Replication of the numerical claims in Assi, Vogl, Kumari & LeBlanc, *Beyond trotterization:
Variational product formulas for quantum simulation*, Phys. Rev. B **113**, 214314 (2026),
[arXiv:2511.15124](https://arxiv.org/abs/2511.15124), DOI 10.1103/941r-ql3g.

Source PDF: `../../Research_Papers/Beyond_Trotterization_Variational_Product_Formulas.pdf`

## Run

```bash
~/micromamba/envs/qml-julia/bin/python replicate.py --target all
```

`--target {t1,…,t8,all}` selects a single target. `--no-plots` skips the PNGs. `--check` runs
only the pinned T6 diagnostic and exits nonzero if any of its hard-coded numbers drift — it
takes about a second and is the cheap regression gate. `--t0` and `--rcond` exist **only** so
that T6 can report sensitivity; they are frozen everywhere else.

Deterministic: `SEED = 0`, and in fact no target consumes the RNG at all — every input is a
fixed grid or a fixed parameter point. Runtime for `--target all` was 39 min on an idle machine
and 72 min on a loaded one — numpy threads out over the available cores, so wall time is not a
stable quantity here. T3 is the single largest target (two RK45 integrations at `N = 10`, 2600
and 2432 right-hand-side evaluations, each RHS being six `1024×1024` conjugations), then T8;
everything else finishes in seconds. Two consecutive `--target all` runs were diffed: all 468
lines are byte-identical except the three that print wall-clock seconds — same RHS-evaluation
counts (2600 and 2432), same digits in every one of the 74 summary rows. Dependencies are numpy 2.4.3, scipy 1.17.1 and matplotlib 3.10.8 —
all already in `qml-julia`, **nothing was installed**.

## Targets

| | Claim | Where the paper reports it | Verdict |
|---|---|---|---|
| T1 | Table I gate counts; $N=10$ reductions of 50 % ($R_x$) and ~33 % ($R_z$, CNOT) | Table I p.8, text p.8 | reproduced (exactly, all $N$) |
| T2 | Two-level closed forms, Eqs. (31)–(36) | Eqs. (31)–(36) p.7 | reproduced (to $3\times10^{-14}$) |
| T2′ | Fig. 2(b) error floor "of order $10^{-7}$", attributed to RK45 | text p.7, Fig. 2(b) p.6 | **mechanism** confirmed, **value** inconclusive (we get $4\times10^{-11}$) |
| T3 | Fig. 4(a) QIM error at $t=1$, $N=10$ | Fig. 4(a) p.8 | reproduced (~1 %) |
| T3′ | same, against the prose "1.5 % … over 3 %" | text p.8 | **not reproduced** (+41 %, +73…129 %) |
| T3″ | "lower error than Ruth over the entire time interval $[0,2]$" | text p.8 | **not reproduced** (fails at small $t$; Ruth is $O(t^5)$, the ansatz $O(t^3)$) |
| T4 | Eq. (69) approximate angles; Fig. 13 at $t=100$ | Eqs. (61)/(69) p.13, Fig. 13 p.14 | reproduced (0.04 %, 0.19 %) — but only against the **unnormalized** norm |
| T5 | Ordering parameter $\Xi$, Eq. (65) → Eq. (68), and $\Xi$ predicts the better ordering | Eqs. (64)/(65)/(68) p.13 | reproduced **only under a consistent A↔B relabelling** (7/7 vs 0/7) |
| T6 | cond($g$) pinned diagnostic | **ours — the paper reports cond($g$) nowhere** | reproduced (all pinned values) |
| T7 | Fig. 11(b) XXZ stroboscopic error: 11 % vs 5.5 % at $t=100$ | text p.11, Fig. 11(b) | reproduced (−0.94 %, +0.55 %) |
| T8 | $\mathcal{M}_z$ relative error "a factor of 4 or 5 times better" | Fig. 3 p.8, text p.8 | partially reproduced (median 4.21 / 4.76; the "grows with $\tau$" half is statistic-dependent) |

## Four internal inconsistencies in the paper

These are the substantive findings, all of them numerically pinned rather than argued:

1. **Fig. 13 is missing the Eq. (27) normalization.** Its axis is labelled
   $\mathcal{E}_F^{\rm Strob}$, but the plotted values match the raw
   $\lVert U-U_a\rVert_F$ to four significant figures (0.26951 vs a read-off of 0.2694;
   0.17182 vs 0.1715), i.e. they are larger than Eq. (27) by exactly $2\sqrt{D}=11.3137$.
   Fig. 11(b) at the **same** $D=32$ does carry the normalization and reproduces to under 1 %,
   so the two figures disagree with each other. The physics conclusion of Fig. 13 is a ratio
   and survives.

2. **The text on p.8 disagrees with Fig. 4(a) on the same page.** The prose says 1.5 % and
   "over 3 %" at $t=1$; pixel-calibrating the paper's own figure gives 0.0212 and 0.0692.
   This replication lands on the figure (0.02114–0.02122, 0.06885, 0.05203) and therefore
   misses the prose. The qualitative claim ("more than twice") holds either way — obtained
   ratio 2.46.

3. **Eqs. (61), (65), (68) and (69) use A and B interchanged relative to the Eq. (37)/(54)
   splitting.** With Eq. (54)'s labels, Eq. (65) reproduces Eq. (68) in **0/15** parameter
   points and the ordering rule predicts the empirically better ordering in **0/7**; with
   "A" meaning the Ising block throughout, both become **15/15** and **7/7**, Eq. (61)
   reproduces Eq. (69) to $3.6\times10^{-16}$, and the rule then agrees with the paper's own
   choice of BAB in Fig. 3 and AB in Fig. 2(a). One relabelling reconciles everything, so it
   is systematic rather than a one-off typo.

4. **"Lower error than Ruth over the entire time interval $[0,2]$" (p.8) is false at small
   $t$.** Ruth is a fourth-order formula, so its error is $O(t^5)$; the four-exponential
   ansatz reduces to Trotter–Suzuki as $t\to0$ by Eq. (10), so its error is $O(t^3)$. Ruth
   must therefore win near the origin, and it does. Fitted slopes over the fixed window
   $t\in[10^{-2},10^{-1}]$: $t^{3.00}$ (variational) and $t^{5.00}$ (Ruth); tightening `rtol`
   by two orders leaves both unchanged to three digits, so neither is integrator noise — though
   that last check was run at $N=5$ in a scratch session and is **not** reproduced by this
   script, which integrates T3 only at the frozen `rtol`. The crossover is bracketed by the
   sampling grid to $0.375<t<0.426$ and interpolates to $t=0.376$; at $t=10^{-3}$ Ruth is five
   orders more accurate. Fig. 4(a)'s linear axes pin every curve to zero below $t\approx0.5$,
   which is presumably how the claim survived. The paper's *useful*
   claim — variational wins where the error is large enough to matter, at a third of the gate
   count — is untouched.

## One underspecification

The Fig. 2(b) floor "of order $10^{-7}$" is not a property of the method. The closed-form
Eqs. (32)/(33) give $E_F\le4.3\times10^{-16}$ — the 3-exponential ansatz is *exact* for the
two-level model — so everything above machine precision in Fig. 2(b) is integrator noise. That
much reproduces, and it is the paper's own explanation. But the *size* of the floor is set by
the RK45 `rtol`/`atol`, which the paper never states: at `rtol=1e-10, atol=1e-12` we get
$3.96\times10^{-11}$, four orders below their $10^{-7}$. Sweeping the tolerance shows the floor
tracks `rtol` essentially linearly:

| `rtol` / `atol` | $10^{-3}/10^{-6}$ | $10^{-4}/10^{-7}$ | $10^{-5}/10^{-8}$ | $10^{-6}/10^{-9}$ | $10^{-8}/10^{-10}$ | $10^{-10}/10^{-12}$ | $10^{-12}/10^{-14}$ |
|---|---|---|---|---|---|---|---|
| max $\mathcal{E}_F$ | $9.9\times10^{-4}$ | $1.4\times10^{-4}$ | $1.4\times10^{-5}$ | $1.2\times10^{-6}$ | $6.9\times10^{-9}$ | $4.0\times10^{-11}$ | $3.3\times10^{-13}$ |

So $10^{-7}$ is a free parameter of their run rather than a result — it brackets to an `rtol`
between $10^{-8}$ and $10^{-6}$, which the paper does not state. T2′ is therefore
*inconclusive*: the claim is unfalsifiable as written, and the only thing Fig. 2(b)'s floor
height measures is the authors' integrator settings.

A fourth inconsistency, smaller: Fig. 12's printed gate angles ($\alpha=-2ch_z$, $\beta=-2cJ$,
$\gamma=-2ch_x$) require $S=\sigma$, while Eq. (37) states $S=\sigma/2$. T1 builds the circuit
both ways and shows ours matches `expm` to $2.8\times10^{-16}$ while Fig. 12's misses by 1.73.
Gate *counts*, which is all Table I claims, are unaffected.

## Figures

| File | Reproduces |
|---|---|
| `figures/table1_gate_counts.png` | Table I p.8 (counts emitted from an actual gate list, not a formula) |
| `figures/fig2a_two_level.png` | Fig. 2(a) p.6 |
| `figures/fig2b_two_level.png` | Fig. 2(b) p.6 |
| `figures/fig2b_error_floor.png` | no paper counterpart — shows the Fig. 2(b) floor is RK45 noise, not ansatz error |
| `figures/fig4a_qim_error.png` | Fig. 4(a) p.8 |
| `figures/fig4b_qim_strob.png` | Fig. 4(b) p.8 |
| `figures/fig3_magnetization.png` | Fig. 3(a)–(d) p.8 |
| `figures/fig11b_xxz_strob.png` | Fig. 11(a) and 11(b) p.11 |
| `figures/fig13_approx_vs_ts.png` | Fig. 13 p.14 — **plotted as the raw Frobenius norm**, which is what the paper's axis actually contains |
| `figures/t4_tau_sweep.png` | no paper counterpart — single-layer error vs $\tau$ for TS, Eq. (69) and full RK45 |
| `figures/fig_xi_phase_diagram.png` | no paper counterpart — the Eq. (68) ordering boundary $h_x^2=h_z^2+J^2/4$ |
| `figures/cond_g_trajectory.png` | no paper counterpart — cond($g$) along the trajectory (T6) |

## Layout of `replicate.py`

Frozen numerical choices first (`T0`, `RCOND`, `RTOL`, `ATOL`), then the constants read out of
the paper, each with its page. Then the spin algebra and a `Generator` class that carries one
block of the splitting in one of three representations — `diag` (the QIM Ising block, where
$e^{icO}$ is a vector), `prod` (the QIM transverse block, a Kronecker power of a $2\times2$),
`dense` (everything else). This is what makes $N=10$ affordable; `_validate_fast_path()`
asserts all of it against `scipy.linalg.expm` at $N=4$ before T3 uses it.

The equations of motion are built twice, deliberately. `g_b_generic()` is the literal
transcription of Eq. (3) — it forms every $\partial U_a/\partial c_j$ densely and takes the
traces, leaving the Re/Im ambiguity visible. `g_and_b()` is the fast route: writing
$\partial_j U_a = i\tilde O_j U_a$ with $\tilde O_j = P_j O_j P_j^\dagger$ collapses Eq. (3)
to $g_{jk}=\mathrm{Tr}[\tilde O_j\tilde O_k]/D$ and $b_j=\mathrm{Tr}[\tilde O_j H]/D$, which
are traces of products of Hermitian matrices and so manifestly real — $U_a$ never enters the
RHS at all. T2 asserts the two against each other and against the paper's printed Eq. (35) and
Eq. (36), which is what settles the convention ($\mathrm{Tr}[\cdot]/D$, $b=+iF$) before any
downstream target runs. The measured imaginary parts are $4\times10^{-16}$, so the Re/Im
question turns out to be vacuous for this ansatz class rather than a real choice.

`ASSUMPTIONS` at the top of the file is the authoritative list of what the paper left open;
`--target all` prints it alongside the consolidated comparison table.

`T6_PINNED` carries two sets of numbers per row and `--check` asserts both. The **regression**
baseline is the full-precision value the script computes, compared at `T6_RTOL = 1e-9` — a
floating-point reproducibility margin, not a fit knob. The vault note's three-significant-figure
table is carried alongside and checked at `T6_NOTE_RTOL = 5e-3`, which is all a rounded prose
table can support. The two are separate on purpose: an earlier version compared only against the
rounded figures, which forced a 2 % band, and a 2 % band is wide enough to swallow a real
regression — a deliberate 1.55 % perturbation of $\det g$ at $c_1=0.2$ passed it silently.
Verified both ways: perturbing a pinned value by 1.55 % or by $8\times10^{-6}$ relative makes
`--check` exit 1; perturbing by $8\times10^{-12}$ does not, which is the intended float margin.
Note the asymmetry in the evidence — the old gate's silent pass was run against the pre-fix
script and its output survives; the new gate's failure was watched interactively and not saved.
If you are about to rely on this gate, perturb one `T6_PINNED` value and confirm it exits 1.

Nothing in this script was tuned to make a plot match. `T0`, `RCOND`, `RTOL` and `ATOL` were
fixed before the first accuracy run; T6 shows that varying `t0` over $[10^{-8},10^{-3}]$ and
`rcond` over $[10^{-14},10^{-4}]$ moves the final $\mathcal{E}_F$ by at most
$5.12\times10^{-10}$ relative, so none of the results below depend on them.
