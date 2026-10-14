---
title: "M5 — Reproducibility hardening"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]]"
---

# M5 — Reproducibility hardening

> [!success] Implemented and confirmed — see [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE#6. Implementation status|hub §6]]
> Each target now seeds its own `default_rng(SEED_{A,B,C,D})`; verified `--target b` alone
> matches `--target b` inside `--target all`. `stop` restored to `make_rho_fn`; new
> `check_eq_b9_self()` passes (worst diff $8.9\times10^{-11}$) against this script's own
> `metric_M`/`drho`. `log["accepted"]`/log["reasons"]` now distinguish forced-through steps
> from real retries; `vector_Y`'s two simplifications documented in its docstring.

**Findings:** F7, F8 (medium), F11, F12 (low)  **Effort:** M  **Do before** adding Routes 1/2/4/5,
not after — every one of them will inherit this layer.

## F7 — targets are not independently reproducible

`rng` is module-level (`adaptive_qite.py:44`) and consumed in sequence by the four targets. Target A
draws 5 values, B draws 5, C draws 9, D draws 5. So `--target d` on its own draws a *different*
$\theta_0$ than `--target d` inside `--target all`.

Two consequences:

- `Attack_Route_03/README.md` says "Deterministic: `SEED = 0`, the only randomness is the initial
  $\theta_0$ draw." That is true only of the exact `--target all` invocation. The per-target flag
  the README documents in the line above silently changes the physics.
- [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] defends the Target B result as "Not a bug:
  reproduced identically across Target B and Target D below." Targets B and D use **different
  initial angles** (D's draw follows C's nine). The claim as written is false; whatever agreement
  was observed was between two different trajectories, which is a weaker — though still meaningful —
  form of corroboration.

**Fix:** give each target its own `np.random.default_rng(SEED)`, or better, a per-target seed
derived from a named constant (`SEED_A`, `SEED_B`, …) so the draws are stable, independent, and
documented. Then re-check whether B and D actually do agree; if they now share $\theta_0$, the
"reproduced identically" claim becomes true and can stay.

## F8 — rejected steps are committed to the trajectory

`adaptive_qite_run` (`adaptive_qite.py:348-358`):

```python
    if accepted or retries >= max_retries or delta_tau <= delta_tau_min:
        break
    ...
theta = theta_trial
```

The retry loop can exit with `accepted = False` — on exhausting `MAX_RETRIES = 20`, or on hitting
`DELTA_TAU_MIN` — and the trial step is then taken regardless. A step logged `accepted=False` still
enters the trajectory.

That may well be the right policy (there is nothing better to do at the minimum step size), but it
is undocumented, and it makes `log["accepted"]` unusable as a quality signal: `n_rejected` in Target
A counts steps that were taken anyway. Decide the policy explicitly, state it in `ASSUMPTIONS`, and
distinguish "rejected then retried successfully" from "rejected and taken anyway" in the log. The
second is the one that should alarm a reader.

## F11 — the physics layer is duplicated with no agreement test

`adaptive_qite.py` re-implements the shared layer from
`Replications/practical-limitations-qdp/replicate.py`. AST body comparison, run during this audit:

| Function | Status |
|---|---|
| `drho`, `metric_M`, `vector_Y`, `h2_hamiltonian` | byte-identical |
| `gate_list`, `make_rho_fn`, `make_noisy_rho_fn`, `prepare`, `frobenius_cond`, `global_depolarizing_kraus` | trimmed copies — the scan knobs `two_qubit`, `reverse_wires`, `stop`, `'++'` were dropped |

The trimming is benign for the physics. But `make_rho_fn` lost the `stop` argument, and `stop` is
what `check_eq_b9()` uses to build the partial-circuit states for the Eq. (B9) diagonal check. So
**the validation the note cites — "the same method already validated against the closed form Eq. B9
to 1e-10" — cannot be run against Attack_Route_03's own copy.** It was validated against the other
copy, which happens to be identical today and has nothing enforcing that it stays so.

**Fix:** extract the shared layer into one module both scripts import (e.g.
`Replications/practical-limitations-qdp/qite_core.py`, or a small `common/` alongside), or, if
keeping two copies is preferred for the "self-contained replication" property, add a test to
`adaptive_qite.py` that runs `check_eq_b9()` against its *own* `metric_M`/`drho`. Either closes the
gap; the second is cheaper and preserves the current folder structure.

## F12 — `vector_Y` drops a term without saying so

[[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §2 defines

$$Y_k := -\mathrm{Tr}\big[(\partial_k\rho)^\dagger(\{H,\rho\} - 2\langle H\rangle\rho)\big],$$

calling the $-2\langle H\rangle\rho$ piece "crucial, to maintain normalization and faithful
trajectory tracking". `vector_Y` (`adaptive_qite.py:235-238`, and identically in `replicate.py:242`)
implements only the $\{H,\rho\}$ part.

It is numerically harmless **here**: $\mathrm{Tr}[\partial_k\rho\cdot\rho] = \tfrac12\partial_k
\mathrm{Tr}[\rho^2] = 0$ for a pure state, and it also vanishes under this depolarizing model since
$\mathrm{Tr}[\partial_k\rho] = 0$. But it is undocumented, and it will *not* vanish for a general
channel, a non-purity-preserving noise model, or a non-normalized parametrization — all of which are
on the roadmap in [[QITE_Hand_Off/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]. The dagger in the same definition is likewise a no-op only
because $\partial_k\rho$ is Hermitian.

**Fix:** docstring only. State both simplifications and the conditions under which they hold. No
code change.

Also in passing: `n_params` is an accepted-but-unused parameter of both `adaptive_qite_run` and
`fixed_step_qite_run`.

## How to verify

- `--target b` alone and `--target all` produce the same Target B numbers (F7).
- The Eq. (B9) check passes when run from inside `adaptive_qite.py` (F11).
- A deliberate forced rejection at `DELTA_TAU_MIN` shows up in the log as distinct from a
  successful retry (F8).

## Documents to update

- `Attack_Route_03/README.md` — the determinism claim in the Run section.
- [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] — the "reproduced identically across Target B and
  Target D" sentence, and Assumption 5's citation of the Eq. (B9) validation.
