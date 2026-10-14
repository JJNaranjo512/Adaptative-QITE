# Adam_QITE

An **Adam-like step-size controller** for variational QITE, built on Route 3's indicators.
Extends `../Attack_Route_03/adaptive_qite.py`, which drives $\Delta\tau$ with a hand-tuned multiplicative
ratchet (`GROW = 1.25` on a clean accept, `SHRINK = 0.5` on rejection) — Assumption 1 of that
script flags both numbers as tuned by eye, and [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R3]] names replacing them as a
priority.

This is **our own extension, not part of the memo**. The memo proposes no controller of this
kind; there is nothing to score against, so there is no verdict table. Everything below is what
the code produced, including where it lost.

## The idea

A local truncation error estimate was already being computed and thrown away. `energy_check()`
forms $\Delta E_{\text{actual}} - \Delta E_{\text{pred}}$ with
$\Delta E_{\text{pred}} = -(Y\cdot\dot\theta)\Delta\tau$ the *exact* first-order model of the
Euler step, so the difference is $O(\Delta\tau^2)$ — an error proxy at **zero extra circuits**,
inside Route 3's "cheap indicators, no extra stages" premise. Today it only feeds a binary
in-band/out-of-band retry test. [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] states the project "has no local truncation
error estimate at all" and routes R1 (step-doubling, 3 $f$-evals per step) through that gap.

Adam is then not decoration: $\hat m/\sqrt{\hat v}$ is a signal-to-noise ratio on that estimator,
so growth is throttled when the error signal cannot be told apart from its own fluctuation.

$$\hat e_n = \frac{|\Delta E_{\text{actual}} - \Delta E_{\text{pred}}|}{\max(\texttt{DE\_BAND\_ABS},\ \texttt{DE\_BAND\_REL}\,|\Delta E_{\text{pred}}|)}, \qquad
s_e = \frac{\hat m_e}{\sqrt{\hat v_e}+\varepsilon}, \qquad
\Delta\tau \leftarrow \Delta\tau \cdot \mathrm{clip}\!\left(\hat m_e^{-s_e k_I},\ \texttt{SHRINK},\ \texttt{GROW}\right)$$

$\hat e_n$ is normalized by the band `energy_check()` *already* uses to reject a step, so
$\hat e_n > 1$ is exactly the existing reject condition — no new tolerance is introduced.
`GROW`/`SHRINK` survive only as clip bounds (a trust region), not as the step law. The McLachlan
residual gets its own EMA pair and enters as a growth **ceiling** only, preserving parent
Assumption 2 ($r_n$ is $\Delta\tau$-independent, so it can never justify shrinking the current
step).

**Scope:** the moments drive the **scalar** $\Delta\tau$ only. The update stays
$\theta \leftarrow \theta + \Delta\tau\,M^{+}Y$, so this is still a genuine McLachlan/QITE
trajectory and $\tau$ keeps its physical meaning. Per-parameter Adam on the flow vector (a
diagonal preconditioner overwriting $M^{+}$) is deliberately **not** implemented — it would stop
being imaginary-time evolution.

## Run

```bash
~/micromamba/envs/qml-julia/bin/python adam_qite.py --target all
```

`--target {e1,e2,e3,all}` selects a scenario; `--no-plots` skips the PNGs. Runtime for
`--target all` is about 75 s, of which E2 is roughly 65 s. No new dependencies. Full stdout of
the run these numbers come from is in `Output_Pre_Audit.txt`; `Output_Post_Audit.txt` is a fresh
re-run after the import repair, and reproduces every physics number in it.

```bash
~/micromamba/envs/qml-julia/bin/python verify_audit.py --target all
```

re-derives the audit's headline numbers for G2, G5 and G6 from this module and its parent without
modifying either (~28 s). It has **no tolerance flag** by design: if a number moves, that is a
result to be recorded in `../Audits/Audit_Adam_QITE/`, not a comparison to be loosened.

`adaptive_qite.py` is **imported, not copied** (`sys.path` insert on the sibling
`Attack_Route_03/`) —
`Attack_Route_03` itself was copy-pasted out of the QDP replication and has had to be resynced
by hand once. Every shared threshold (`GROW`, `SHRINK`, `DELTA_TAU*`, `DE_BAND_*`,
`RESIDUAL_REL_TOL`, `RCOND_DEFAULT`, `STIFFNESS_CAP_COEFF`, `MAX_STEPS`, `MAX_RETRIES`,
`TAU_MAX`) is imported and never redefined, so the two controllers differ **only** in how
$\Delta\tau$ evolves between accepted steps. Each target reuses the seed of the Route-3 target it
is compared against (E1↔A, E2↔D, E3↔C), so $\theta_0$ is identical.

## Targets

| | Scenario | Takeaway |
|---|---|---|
| E1 | Adam vs Route-3 ratchet vs fixed step, noiseless, 5-param (`SEED_A`) | **A wash.** Adam finishes $1.748\times10^{-3}$ above the variational floor in 190 steps / 2285 circuits; Route 3 finishes $1.607\times10^{-3}$ above it in 172 steps / 2069 circuits. At the matched 732-circuit budget Adam is marginally ahead ($-1.121027$ vs $-1.120992$), at the end marginally behind. (Those two budget figures are produced by `energy_at_budget`, which [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility\|G6_Attribution_And_Reproducibility]] found credits each arm with one step it has not finished paying for; corrected they are $-1.120983$ vs $-1.120952$. Both arms move the same way, the sign is preserved and **this verdict is unchanged** — the adaptive-vs-fixed gap is 925× the correction.) Both crush the fixed step ($+7.909\times10^{-3}$ above floor). No claim of improvement is supportable here. The Route-3 arm reproduces Target A exactly (172 steps, $E=-1.122993$, 2069 circuits), so the comparison is measuring the step law and nothing else. |
| E2 | Adam vs Route 3 under global depolarizing noise, $p\in\{10^{-4},10^{-3},10^{-2},5\times10^{-2}\}$, debiased (`SEED_D`) | **The premise is untested — verdict withdrawn after the Adam audit.** This was meant to be the load-bearing test (the ratchet reacts to every individual noisy $\Delta E$; the EMA was supposed to average them), and it cannot decide the question: the depolarizing channel here is *deterministic*, so there is no stochastic noise for an EMA to average (parent Assumption 7). The RMS ratio $1.0243/1.0246/1.0245/1.0088$ is $p$-independent to $\pm0.016$, and at $p\le10^{-3}$ the signal is smaller than the gap between the two baselines being differenced — so what the ratio measures is the stall below, not noise rejection. The "Adam's own $p=0$ baseline **stalls**" claim (500 steps to $\tau=2.86$ vs Route 3's 234 to $\tau=3.01$) is **withdrawn as an artifact of mismatched EMA horizons**, not a property of the controller: matched, it finishes in 244 steps at $\tau=3.0093$. See [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G3_Noise_Experiment_Cannot_Resolve_Noise\|G3_Noise_Experiment_Cannot_Resolve_Noise]] and [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G2_Beta2_Horizon_Dominates_Every_Verdict\|G2_Beta2_Horizon_Dominates_Every_Verdict]]. |
| E3 | Rank-6-of-9 (9-param) ansatz, `rcond` $\in\{10^{-8},10^{-2}\}$, four arms (`SEED_C`) | **A win, but not cleanly attributable — verdict downgraded after the Adam audit.** With the $\mathrm{cond}(M)$ cap disabled, Adam reaches error $5.653\times10^{-4}$ in **28 steps**, better than any capped run including Route 3's $8.868\times10^{-4}$ in 108 steps; and removing the cap *alone* makes Route 3 **worse** ($1.398\times10^{-3}$ in 22 steps). Those numbers stand. What does not stand is "attributable": the effect is **re-attributed to the cap**, since four of the eight rows below are duplicate trajectories (the cap binds identically), the "$3.9\times$ fewer steps" framing has **no cost-matched column** and no matched-mean-step control arm, and this is $n=1$ — one seed, one ansatz, one configuration. At the well-conditioned `rcond=1e-2` Adam loses ($3.314\times10^{-3}$ in 14 steps vs $2.143\times10^{-3}$ in 17). See [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility\|G6_Attribution_And_Reproducibility]]. |

### E3 in full

| `rcond` | arm | final $E$ | error | steps |
|---|---|---|---|---|
| $10^{-8}$ | Route 3 (cap on) | $-1.144712$ | $+8.868\times10^{-4}$ | 108 |
| $10^{-8}$ | Route 3 (cap off) | $-1.144201$ | $+1.398\times10^{-3}$ | 22 |
| $10^{-8}$ | Adam (cap on) | $-1.144712$ | $+8.868\times10^{-4}$ | 108 |
| $10^{-8}$ | **Adam (cap off)** | $\mathbf{-1.145034}$ | $\mathbf{+5.653\times10^{-4}}$ | **28** |
| $10^{-2}$ | Route 3 (cap on / off) | $-1.143456$ | $+2.143\times10^{-3}$ | 17 |
| $10^{-2}$ | Adam (cap on / off) | $-1.142285$ | $+3.314\times10^{-3}$ | 14 |

Adam (cap on) is *identical* to Route 3 at `rcond=1e-8` — not a bug: $\mathrm{cond}(M)\sim9.66\times10^6$
makes the cap $2000/\mathrm{cond}(M)=2.07\times10^{-4}$ bind on every step, pinning $\Delta\tau$ and
overriding both step laws the same way. At `rcond=1e-2`, $\mathrm{cond}(M)\le125$ so the cap never
binds and cap-on/cap-off coincide. Both identities are sanity checks that passed.

## Why Adam under-performs: the diagnostic

`adam_signal_diagnostic()` prints this, at both seeds:

| | E1 (`SEED_A`) | E2 baseline (`SEED_D`) |
|---|---|---|
| $\hat e$ median | $2.902\times10^{-5}$ | $1.261\times10^{-5}$ |
| $\hat e$ mean | $2.751\times10^{-2}$ | $2.920\times10^{-2}$ |
| $\hat e$ mean/median | **947.9** | **2316.6** |
| mean $s_e$ | 0.1401 | 0.0893 |
| effective exponent $s_e k_I$ | 0.02102 | 0.01339 |
| growth frozen ($<1.001$) | 49 % of steps | 72 % of steps |
| ungated law would average | 1.2467 | 1.2439 |

Two independent findings, both about the design rather than the tuning:

1. **Adam's variance normalization misreads a heavy-tailed signal.** $\hat e$ is non-negative and
   spiky — a mean/median ratio of 948–2317 — so $\sqrt{\hat v}$ is set by rare large values while
   $\hat m$ sits near the median, and $s_e$ collapses to $\approx0.09$. That throttles the
   Gustafsson exponent $11\times$ below $k_I=0.15$ and freezes growth on up to 72 % of steps,
   which is exactly why E2's baseline needs 500 steps. Adam's second moment was designed for
   minibatch gradient noise, which is roughly symmetric about its mean; a truncation-error
   *magnitude* is not, and the ratio does not mean the same thing on it.
2. **The inherited band is far too loose to be a controller setpoint.** $\hat e > 1$ on
   **0 of 190** steps (max $\hat e = 0.945$), so $\hat m_e < 1$ always, so
   $\hat m_e^{-s_e k_I} > 1$ always: the law **never shrinks** (growth $< 1$ on 0/190 steps) and
   all shrinking is delegated to the retry gate. Reusing `energy_check`'s reject band as the
   setpoint was the right call for *not inventing a tolerance*, but a good rejection threshold
   and a good control setpoint turn out to be different objects. Not fixed here — `DE_BAND_*` are
   imported from the parent script, and retuning them after seeing this is exactly the fudge the
   project's rules forbid.

The controller is a strict generalization of the ratchet, which the degeneracy check confirms:
at $\beta_1=\beta_2=0,\ k_I=10$ it collapses to **exactly** Route 3 (172 steps,
$E=-1.122993$, growth pinned at `GROW` on 100 % of steps). Memoryless with the stock $k_I$ gives
172 steps and $E=-1.122928$ (95 % pinned). So the moments are genuinely doing the work in the
default configuration — they are just, on this problem, doing it in the wrong direction.

These four numbers are no longer hand-run: `verify_audit.py --target g6-2` regenerates
them and checks the collapse **element-wise** against the Route-3 trajectory (identity holds to
$0$ in both $E$ and $\Delta\tau$). The $k_I=10$ arm overflows `m_e_hat ** (-s_e * k_i)` to `inf`
on 10 steps, which `np.clip` pulls back to `GROW` — harmless, and precisely why the collapse is
exact rather than approximate.

Two caveats on the diagnostic table above, both from the audit. The `mean s_e` and
`frozen`-percentage columns are **artifacts of the mismatched EMA horizons**, not properties of the
controller: with $\beta_1=\beta_2$ at any common value, $s_e$ recovers to $0.39$–$0.52$ and growth
freezes on **0 %** of steps. And finding 1's heavy-tail explanation, while the tail itself is real
(mean/median 948–2317), is **not the cause** of the freezing — see
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G2_Beta2_Horizon_Dominates_Every_Verdict|G2_Beta2_Horizon_Dominates_Every_Verdict]]. Finding 2 stands unchanged and is the more important
of the two.

## Figures

| File | Scenario |
|---|---|
| `figures/target_e1_controllers.png` | E1: energy convergence, $\Delta\tau$ trajectories, and the $s_e$/$s_r$ confidence ratios |
| `figures/target_e2_noisy.png` | E2: RMS deviation vs. $p$ for both controllers, and $s_e$ vs. $\tau$ per noise level |
| `figures/target_e3_stiffness.png` | E3: $\mathrm{cond}(M)$ and $\Delta\tau$ vs. $\tau$ for all four arms at `rcond=1e-8` |

## Layout of `adam_qite.py`

Adam constants (`BETA1`, `BETA2`, `EPS`, `K_I`) and `ASSUMPTIONS` first, then
`AdamStepController` (the two EMA channels and the control law), then `adam_qite_run` — Route 3's
measurement, indicators and accept/reject gate verbatim, with only the step law swapped, returning
the same log dict shape as `adaptive_qite_run` so `trajectory_deviation` and the M6 circuit
accounting apply unchanged. Then `circuit_cumsum`/`energy_at_budget`/`run_summary` (M6 currency),
`adam_signal_diagnostic`, and one function per target.

## Fairness protocol

Since the hyperparameter ablation is not in this pass, this stands in for it and is in
`ASSUMPTIONS` verbatim:

- $\beta_1=0.9$, $\beta_2=0.999$, $\varepsilon=10^{-8}$ — **stock Adam defaults, not tuned here**.
- $k_I=0.15$ — Gustafsson's $0.3/(p+1)$ at $p=1$ (explicit Euler), the value [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R3]]
  already proposes. **Not tuned here.**
- Every shared knob imported from `adaptive_qite.py`, never redefined.
- EMAs update once per *accepted* step, so a retry cascade does not feed rejected trials into the
  moments; bias correction counts accepted steps.
- **No constant was adjusted after seeing a result.** Where Adam loses, that is the finding.

## Open

> [!important] Read the report first
> [[QITE_Hand_Off/Adam_QITE/Report_Adam_QITE|Report_Adam_QITE]] is this folder's entry point: what the module does, how
> $\Delta\tau$ is actually updated, and which of the claims above survived the Adam audit.
> The full findings register is [[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]] (G1–G7). Most of its headline
> numbers are regenerated by `verify_audit.py --target all` — 147 comparisons, 0 drift,
> stdout in `Output_Audit_Verify.txt`. That script also prints what it does *not* cover;
> read that section before treating the count as a clean bill of health.
>
> Three things this README got wrong and now says correctly: E2's verdict, E3's attribution, and
> the diagnostic's causal story. One thing it understated: the $O(\Delta\tau^2)$ free-error-estimate
> premise in "The idea" above is **not just defensible but quantitatively confirmed** — it tracks
> R1's step-doubling estimator with log-Pearson $0.9902$ and the predicted ratio of 2 to within 3 %,
> at zero extra circuits. That is the result worth carrying forward, and it is the one the code
> never checked. See [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot|G5_Error_Estimator_Validated_And_Its_Blind_Spot]].
>
> Also note that Assumption 4's claimed bound $s_e \in (0,1]$ is **false as implemented**: the run
> reports $s_e$ up to $1.015$. Jensen's inequality needs one weighting; the code uses two.

The $\beta_1/\beta_2$ sensitivity sweep **has since been run** (audit V9/V10, two seeds, four
configurations, in matched-circuit currency) and it did not go the way this section predicted: the
fix that matters is **tying the two horizons together** ($\beta_1 = \beta_2$ at any common value, or
both derived from the expected step count), not replacing $\sqrt{\hat v}$ with a robust dispersion
estimate. The heavy tail is real but it is not the whole problem, so MAD-in-place-of-RMS drops down
the list. Still genuinely open: the momentum-only vs. RMS-only ablation, the sweep resolved against
both `rcond` values, and a setpoint separated from the reject band — the band cannot be retuned
after the fact without becoming exactly the fudge the project's rules forbid.

Separately, E3's cap-off win deserves checking on the 5-parameter ansatz and at other seeds before
it is treated as general — it is currently $n=1$, and it now also needs a cost-matched column and a
matched-mean-step control arm before "attributable" can be claimed again.
