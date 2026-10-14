---
title: "G4 — the control law is structurally incapable of shrinking"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]"
---

# G4 — The controller cannot shrink

**Severity: high.** Closes G4. Recorded, not fixed.

## The observation the README makes

`Output_Pre_Audit.txt`, E1:

```
    setpoint diagnostic: e_hat > 1 (gate band) on 0/190 steps; growth < 1 on 0/190;
                         growth pinned at GROW=1.25 on 10/190
```

README:108-115 reports this and attributes it to the inherited band being "far too loose to be a
controller setpoint", concluding that "a good rejection threshold and a good control setpoint turn
out to be different objects". The observation is right and the conclusion is right as far as it
goes. The attribution is wrong, and the difference matters because it changes what the fix is.

## The actual mechanism

It is a circularity, not a magnitude. Three facts compose:

1. $\hat e$ is normalized by `energy_check`'s own reject band (`:281-282`), so **$\hat e > 1$ is
   by definition the reject condition** (Assumption 3 says so explicitly, as a virtue).
2. The retry loop (`:266-277`) shrinks $\Delta\tau$ by `SHRINK` and re-tries **until the step is
   in band**, and only then exits.
3. `ctrl.update` is called *after* the retry loop (`:283`), on the surviving trial.

So the controller only ever observes $\hat e \le 1$. Hence $\hat m_e < 1$ always, hence
$\hat m_e^{-s_e k_I} > 1$ always, hence `growth` $\ge 1$ always. The law is not merely
mis-scaled — **it is one-sided by construction**, and the `growth < 1 on 0/190` line is not a
symptom of tuning but an identity.

**Tightening `DE_BAND_*` does not fix it.** A tighter band tightens the gate in lockstep, so the
controller again sees only in-band errors. Any tolerance used simultaneously as the rejection
threshold and the control setpoint produces this. Assumption 3 presents the reuse as introducing
"no new tolerance", which is true and is good discipline — but it is precisely what makes the
controller blind on one side.

V11 sweeps the band over three orders of magnitude (E1, `SEED_A`, both the gate in `energy_check`
and the normalization in `adam_qite.py:281` moved together):

```
   DE_BAND_ABS  DE_BAND_REL  steps   growth<1   e_hat>1  retries  median e_hat      final E
       1.0e-03         0.50    190      0/190     0/190        5     2.902e-05    -1.122852
       5.0e-04         0.50    191      0/191     0/191        5     5.822e-05    -1.122854
       1.0e-04         0.50    193      0/193     0/193        5     2.945e-04    -1.122857
       1.0e-05         0.50    192      0/192     0/192        5     2.996e-03    -1.122859
       1.0e-06         0.05    269      0/269     0/269        5     4.707e-02    -1.123546
```

`growth < 1` is **0/N in every configuration**, and so is $\hat e > 1$. A 1000$\times$ tighter
absolute band and a 10$\times$ tighter relative band do not produce a single shrinking step. The
README's "band too loose" reading is refuted directly.

One refinement to the prediction this note made before running it: the shrinking does not migrate
into the gate either — the retry count is pinned at 5 throughout. What actually happens is that
$\hat m_e$ rises toward 1, so `growth` $= \hat m_e^{-s_e k_I}$ is squeezed toward 1 **from above**,
and the run simply takes more steps (190 → 269). The law's reachable range is $[1, \texttt{GROW}]$
and tightening the band only compresses it, never reflects it.

## What actually sets $\Delta\tau$, then

| Direction | Mechanism | Origin |
|---|---|---|
| down | retry loop, $\times\,\texttt{SHRINK} = 0.5$ per rejection (`:276`) | Route 3's, verbatim |
| down | proactive cap $\texttt{COEFF}/\mathrm{cond}(M)$ (`:255-262`) | Route 3's, verbatim |
| up | $\hat m_e^{-s_e k_I}$, clipped at `GROW` (`:209-215`) | the new part |

So the Adam machinery is a **growth modulator bolted onto an unchanged Route-3 shrink mechanism**.
Assumption 4 claims "`GROW` and `SHRINK` survive only as clip bounds (a trust region), no longer as
the step law itself — which is the specific thing Assumption 1 of the parent script complains
about". That is half true at best:

- `SHRINK` is not a clip bound. It is the live shrink law, running in the retry loop at `:276`,
  and it is the **only** thing that shrinks the step.
- The `SHRINK` end of the clip at `:215` is unreachable, since `growth` $\ge 1$ always. Dead code
  in every run.
- `GROW` still binds: 10/190 E1 steps, and 100 % of steps in the degeneracy configuration.

The stated motivation for the whole module — retiring the two hand-tuned ratchet constants flagged
by the parent's Assumption 1 — is therefore not achieved. Both constants are still load-bearing.

## A secondary consequence: the setpoint is unreachable near convergence

`band = max(DE_BAND_ABS, DE_BAND_REL * |dE_pred|)` with `DE_BAND_ABS = 1e-3`. As the run converges,
$\Delta\tau$ and $|\Delta E_{\text{pred}}|$ both shrink, so the relative term falls away and the
band **floors at a constant $10^{-3}$** while the true local error keeps falling as
$O(\Delta\tau^2)$. Hence the E1 median $\hat e = 2.9\times10^{-5}$: the normalized error is not
near its setpoint of 1, it is five orders of magnitude below it and falling. An integral controller
whose error signal is pinned far below setpoint saturates at maximum growth, which is what the
`GROW`-pinning count records.

This dovetails with [[QITE_Hand_Off/Audits/Audit_Project_Wide/M4_Controller_Near_Convergence|M4_Controller_Near_Convergence]], which established for the parent that near
convergence the $\mathrm{cond}(M)$ cap binds on 31 of 43 steps while the energy band triggers zero
retries — the $\Delta\tau$ collapse there is stiffness, not the gate. The same regime is what
neutralizes the Adam law: the gate has nothing to say, the error signal sits far below setpoint,
so the law asks for maximum growth on every step and the cap overrules it. That is the mechanism
behind E3's cap-on rows being identical to Route 3's, in
[[QITE_Hand_Off/Audits/Audit_Adam_QITE/G6_Attribution_And_Reproducibility|G6_Attribution_And_Reproducibility]] §6.1.

## Proposed work

1. **Separate the gate tolerance from the control setpoint.** They are different objects with
   different jobs — the gate answers "is this step safe", the setpoint answers "is this step the
   right size". Introducing a control tolerance `TOL_CTRL` is not the fudge Assumption 3 avoids,
   provided it is declared and not fitted to a result.
2. **Feed the controller the rejected trials.** The information the law is missing lives in
   exactly the steps the retry loop discards. Updating the EMAs with the *first* trial's $\hat e$
   rather than the surviving one would let the law see $\hat e > 1$ and shrink on its own. Note
   this trades against Assumption 7's reasoning (a retry cascade poisoning $\hat v_e$), so it needs
   the horizon question in [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G2_Beta2_Horizon_Dominates_Every_Verdict|G2_Beta2_Horizon_Dominates_Every_Verdict]] settled first.
3. **Make `DE_BAND_ABS` scale-aware**, or drive the controller on the un-normalized
   $|\Delta E_{\text{actual}} - \Delta E_{\text{pred}}|$ against a relative target, so the setpoint
   remains meaningful as the run converges.
4. **Assert the identity.** `growth < 1` is currently impossible; the script should say so, or
   check it. As it stands the diagnostic prints `growth < 1 on 0/190` as though it were a
   measurement.

**Falsifiable prediction, tested.** Before running V11 this note predicted that tightening
`DE_BAND_ABS` would leave `growth < 1` at 0/N while pushing the shrinking further into the gate.
The first half held exactly (0/N at every band, down to $10^{-6}$); the second half was wrong — the
retry count never moved, and the growth range compressed toward 1 instead. The correction is in
§"The actual mechanism" above. The circularity claim survives; the mechanism for what the band
*does* control was mis-stated and has been fixed.
