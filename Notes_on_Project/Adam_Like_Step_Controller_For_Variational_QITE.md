---
title: "Implementation — an Adam-like step-size controller for variational QITE"
tags:
  - implementation
  - physics
  - quantum-computation
  - adaptative-qite
code: QITE_Hand_Off/Adam_QITE/adam_qite.py
---

# Implementation — an Adam-like Step-Size Controller for Variational QITE

> [!warning] Partly superseded by the Adam audit
> Five claims below did not survive [[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]. **The prose underneath is left
> exactly as written** — it is the record of what I thought before the Adam audit, and that record is worth
> more intact than patched. Read it with these corrections in hand;
> [[QITE_Hand_Off/Adam_QITE/Report_Adam_QITE|Report_Adam_QITE]] is the current statement.
>
> - **`:86` — $s_e \in (0,1]$.** False as implemented. The bound needs *one* EMA weighting for both
>   moments; $\beta_1 = 0.9$ against $\beta_2 = 0.999$ breaks the Jensen step, and the run reports
>   $s_e = 1.015$.
> - **`:143` — "the premise fails" (E2).** Withdrawn. E2 carries no stochastic noise for an EMA to
>   average, so it cannot resolve the effect it was built to test. The premise is *untested*.
> - **`:147` — the $p=0$ baseline "stalls" at 500 steps.** An artifact of the same horizon mismatch,
>   not a property of the controller: every matched pair finishes, including
>   $\beta_1=\beta_2=0.999$ at 254 steps and $\beta_1=\beta_2=0.9$ at 244.
> - **`:150` — "the one real win" (E3).** Downgraded. Four of the eight rows are duplicate
>   trajectories, there is no cost-matched column and no matched-mean-step control arm, and $n=1$.
> - **`:209` — "Finding 2 also partly explains E3's win".** The narrower claim made here is the
>   closer one, but the win attributes to the inherited $\mathrm{cond}(M)$ cap, not to Adam.
>
> What *strengthened* instead: the founding premise that `energy_check` already computes a free
> $O(\Delta\tau^2)$ local error estimate. It now has a measurement — it tracks R1 step-doubling to
> the predicted factor of 2, at zero extra circuits.

My own extension of [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Route 3]], not anything the
memo proposes. Nothing to score against, so no verdict table — just what the code did, including
the two places it lost.

Code: `QITE_Hand_Off/Adam_QITE/adam_qite.py`, run with

```bash
~/micromamba/envs/qml-julia/bin/python adam_qite.py --target all
```

Every number below is from the single deterministic run saved in that folder's
`Output_Pre_Audit.txt`. ~3 min, no new dependencies.

## Where this came from

Route 3 as I built it drives $\Delta\tau$ with a ratchet: multiply by `GROW = 1.25` on a clean
accept, by `SHRINK = 0.5` on a rejection. Assumption 1 of that script admits both numbers were
tuned by eye against Target A's own trajectory, and [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R3]] flags replacing them
with a Gustafsson PI controller as a priority. [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] also says the sharper thing:
the project *"has no local truncation error estimate at all"*, and routes R1 (step-doubling, 3
$f$-evals per step) through that gap.

Then I noticed I was already computing one and throwing it away. `energy_check()` forms
$\Delta E_{\text{actual}} - \Delta E_{\text{pred}}$, and after the
[[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence|M4]] fix $\Delta E_{\text{pred}} = -(Y\cdot\dot\theta)\Delta\tau$
is the *exact* first-order model of the Euler step I'm taking. So the leading discrepancy is the
step's own second-order term:

$$\left|\Delta E_{\text{actual}} - \Delta E_{\text{pred}}\right| = O(\Delta\tau^2)$$

A local truncation error proxy for **zero extra circuits** — no embedded stages, no
step-doubling, entirely inside Route 3's "cheap indicators, no extra stages" premise. It's a
proxy *in the energy*, not in $\theta$: it sees the truncation error only along the one direction
I can watch cheaply. But R1 needs 3 $f$-evals per step to get its $e_n$ and this needs none, so
it seemed worth pushing on. Today all it does is feed a binary in-band/out-of-band retry test.

That's what makes the Adam framing more than decoration. Adam's $\hat m/\sqrt{\hat v}$ is a
signal-to-noise ratio, so pointing it at $\hat e$ should throttle growth exactly when the error
estimate stops being distinguishable from its own fluctuation — which is memo Eq. (9)'s question
("when does the estimator stop carrying signal"), answerable inside Route 3 without R5's
shot-noise model.

## The controller

$$\hat e_n = \frac{|\Delta E_{\text{actual}} - \Delta E_{\text{pred}}|}{\max(\texttt{DE\_BAND\_ABS},\ \texttt{DE\_BAND\_REL}|\Delta E_{\text{pred}}|)},
\qquad \hat r_n = \frac{\|M\dot\theta - Y\|}{\|Y\|}$$

Two EMA channels with Adam's bias correction, $\beta_1 = 0.9$, $\beta_2 = 0.999$, then

$$s_e = \frac{\hat m_e}{\sqrt{\hat v_e}+\varepsilon} \in (0,1], \qquad
\Delta\tau \leftarrow \Delta\tau \cdot \mathrm{clip}\left(\hat m_e^{\,-s_e k_I},\ \texttt{SHRINK},\ \texttt{GROW}\right)$$

$\hat m_e^{-k_I}$ on its own is just Gustafsson's integral law with the tolerance folded into
$\hat e$ (setpoint 1). The new part is the exponent being **scaled by the confidence** $s_e$:
consistent error signal $\Rightarrow$ full response, noise-dominated $\Rightarrow$ exponent
$\to 0$ and $\Delta\tau$ freezes.

Three deliberate choices:

- **$\hat e$ is normalized by the band `energy_check` already uses to reject a step.** So
  $\hat e > 1$ *is* the existing reject condition, and I introduce no new tolerance — the
  controller reuses the gate's own threshold as its setpoint. (This turns out to be the bad
  choice. See below.)
- **`GROW`/`SHRINK` get demoted from step law to clip bounds** — a trust region. That was the
  point: they're the specific numbers Assumption 1 complains about.
- **The residual channel can only withhold growth**, never shrink, preserving Assumption 2 of the
  parent script. $\hat r_n$ doesn't depend on $\Delta\tau$, so it can't justify rejecting the
  current step. Same $k_I$ reused rather than inventing a second exponent.

**Scope:** the moments drive the *scalar* $\Delta\tau$ only. The update stays
$\theta \leftarrow \theta + \Delta\tau M^{+}Y$, so this is still a real McLachlan/QITE trajectory
and $\tau$ still means imaginary time. I deliberately did **not** do per-parameter Adam on the
flow vector — a diagonal preconditioner would overwrite the $M^{+}$ natural-gradient geometry and
it would stop being imaginary-time evolution. That's a different (possibly interesting) object,
not this one.

## Keeping the comparison honest

`adaptive_qite.py` is **imported, not copied** — `Attack_Route_03` was itself copy-pasted out of
[[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]] and I've already had to resync it by hand once. Every
shared threshold is imported and never redefined, so the two controllers differ *only* in how
$\Delta\tau$ evolves between accepted steps. Each target reuses the seed of the Route-3 target it
is compared against (E1↔A, E2↔D, E3↔C), so $\theta_0$ is identical.

Since I skipped the hyperparameter ablation this pass, a fairness protocol stands in: $\beta_1$,
$\beta_2$, $\varepsilon$ are stock Adam defaults; $k_I = 0.15$ is Gustafsson's $0.3/(p+1)$ at
$p=1$ straight out of [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R3]]; nothing was tuned here, and **no constant was
adjusted after seeing a result**.

Two checks I'd have caught nothing without:

- The Route-3 arm inside E1 reproduces Target A **exactly** — 172 steps, $E = -1.122993$, 2069
  circuits. If it hadn't, the seed/state sharing would be wrong and nothing downstream would mean
  anything.
- Degeneracy: at $\beta_1 = \beta_2 = 0,\ k_I = 10$ the controller collapses to *exactly* the
  ratchet (172 steps, $E = -1.122993$, growth pinned at `GROW` on 100% of steps). So it's a
  strict generalization and the moments really are what's doing the work in the default config.

## Results

**E1, noiseless, 5-param — a wash.** Adam lands $1.748\times10^{-3}$ above the ansatz's
variational floor in 190 steps / 2285 circuits; Route 3 lands $1.607\times10^{-3}$ above it in
172 steps / 2069. At the matched 732-circuit budget ([[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting|M6]] currency)
Adam is marginally ahead, $-1.121027$ vs $-1.120992$; at the end marginally behind. Both are far
better than fixed-step ($+7.909\times10^{-3}$). I can't claim an improvement from this.

**E2, under depolarizing noise — the premise fails.** This was the load-bearing test. The whole
argument for an EMA was that the ratchet overreacts to each individual noisy $\Delta E$ while the
moving average shouldn't. Adam is **worse at every $p$**: RMS ratios $1.0243$, $1.0246$,
$1.0245$, $1.0088$ at $p = 10^{-4},10^{-3},10^{-2},5\times10^{-2}$. And its own $p=0$ baseline
*stalls* — 500 steps reaching only $\tau = 2.86$, versus Route 3's 234 steps reaching
$\tau = 3.01$.

**E3, the rank-6-of-9 ansatz — the one real win.** With the $\mathrm{cond}(M)$ cap disabled, Adam
reaches error $5.653\times10^{-4}$ in **28 steps** — better than *any* capped run, including
Route 3's $8.868\times10^{-4}$ in 108, at $3.9\times$ fewer steps.

| `rcond` | arm | final $E$ | error | steps |
|---|---|---|---|---|
| $10^{-8}$ | Route 3 (cap on) | $-1.144712$ | $+8.868\times10^{-4}$ | 108 |
| $10^{-8}$ | Route 3 (cap off) | $-1.144201$ | $+1.398\times10^{-3}$ | 22 |
| $10^{-8}$ | Adam (cap on) | $-1.144712$ | $+8.868\times10^{-4}$ | 108 |
| $10^{-8}$ | **Adam (cap off)** | $\mathbf{-1.145034}$ | $\mathbf{+5.653\times10^{-4}}$ | **28** |
| $10^{-2}$ | Route 3 (cap on / off) | $-1.143456$ | $+2.143\times10^{-3}$ | 17 |
| $10^{-2}$ | Adam (cap on / off) | $-1.142285$ | $+3.314\times10^{-3}$ | 14 |

I nearly shipped this with three arms and no "Route 3 (cap off)". That would have been wrong:
"Adam cap-off beats Route-3 cap-on" changes the step law *and* the cap at once, and I couldn't
have attributed the win to either. With the control arm in, it's clean — removing the cap
**alone** makes Route 3 *worse* ($1.398\times10^{-3}$), so the improvement really is the step
law. Note also that at `rcond=1e-2` Adam loses ($3.314\times10^{-3}$ in 14 steps vs
$2.143\times10^{-3}$ in 17), so this is specifically about the ill-conditioned regime.

Two identities in that table are sanity checks that passed, not bugs. Adam (cap on) $\equiv$
Route 3 at `rcond=1e-8` because $\mathrm{cond}(M)\sim 9.66\times10^{6}$ makes the cap
$2000/\mathrm{cond}(M) = 2.07\times10^{-4}$ bind on *every* step, pinning $\Delta\tau$ and
overriding both step laws identically. At `rcond=1e-2`, $\mathrm{cond}(M) \le 125$, the cap never
binds, and cap-on/cap-off coincide.

## Why it under-performs — the part worth keeping

`adam_signal_diagnostic()` prints this at both seeds:

| | E1 (`SEED_A`) | E2 baseline (`SEED_D`) |
|---|---|---|
| $\hat e$ median | $2.902\times10^{-5}$ | $1.261\times10^{-5}$ |
| $\hat e$ mean | $2.751\times10^{-2}$ | $2.920\times10^{-2}$ |
| $\hat e$ mean/median | **947.9** | **2316.6** |
| mean $s_e$ | 0.1401 | 0.0893 |
| effective exponent $s_e k_I$ | 0.02102 | 0.01339 |
| growth frozen ($<1.001$) | 49% of steps | 72% of steps |
| ungated law would average | 1.2467 | 1.2439 |

**Finding 1 — Adam's variance normalization misreads a heavy-tailed signal.** $\hat e$ is
non-negative and spiky: mean/median of 948 to 2317. So $\sqrt{\hat v}$ is set by rare large
values while $\hat m$ sits near the median, and $s_e$ collapses to $\approx 0.09$. That throttles
the Gustafsson exponent **11× below** $k_I = 0.15$ and freezes growth on up to 72% of steps —
which is exactly why E2's baseline needs 500 steps and stalls at $\tau = 2.86$. Adam's second
moment was designed for minibatch gradient noise, which is roughly symmetric about its mean. A
truncation-error *magnitude* is not, and $\hat m/\sqrt{\hat v}$ does not mean the same thing on
it. This is the generalizable result: it's an argument about the shape of the signal, not about
this Hamiltonian or this ansatz.

**Finding 2 — a good rejection threshold is not a good control setpoint.** $\hat e > 1$ on
**0 of 190** steps (max $\hat e = 0.945$). So $\hat m_e < 1$ always, so $\hat m_e^{-s_ek_I} > 1$
always: the law **never shrinks** — growth $< 1$ on 0/190 steps — and all shrinking is delegated
to the retry gate. Normalizing by `energy_check`'s band was the right instinct for *not inventing
a tolerance*, and I'd make the same call again, but the band is ~5 orders of magnitude looser
than where $\hat e$ actually lives, so as a setpoint it never binds. I did **not** retune
`DE_BAND_*` after finding this — they're imported from the parent script, and adjusting them now
to make my own method look better is precisely the fudge the project rules forbid.

Finding 2 also partly explains E3's win, which is worth being honest about: on the rank-deficient
ansatz "never shrink, freeze when unsure" happens to be a *better* policy than the ratchet's
aggressive halving, because the stiffness is structural (rank 6 of 9) rather than transient — the
step doesn't need to shrink, it needs to stop growing. That's a narrower claim than "Adam helps".

## What I'd do next

The diagnostic points at one change: replace $\sqrt{\hat v}$ with a **robust dispersion estimate**
— median absolute deviation, or a second moment taken in the log domain — since the heavy tail is
the entire problem. That keeps the confidence-gating idea, which I still think is right, and drops
the assumption about signal shape that isn't.

Second, E3's cap-off win is $n=1$: one ansatz, one seed. It needs the 5-parameter ansatz and other
seeds before it's more than suggestive.

Third, the $\beta_1/\beta_2$ sensitivity and momentum-only-vs-RMS-only ablation I skipped this
pass. The fairness protocol above stands in, but it's not the same thing.

Nothing here changes [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] yet — the local-error-proxy
observation is the piece most likely to earn a place there (it's a cheaper route to what
[[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R1]] wants), but I'd rather confirm it against a real step-doubling $e_n$ first.
That comparison — does $|\Delta E_{\text{actual}} - \Delta E_{\text{pred}}|$ actually track R1's
$e_n$? — is now the obvious next experiment, and it's the same question [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]] poses
about whether the McLachlan residual tracks $e_n$.
