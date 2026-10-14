"""Regenerate the four load-bearing numerical findings of the Adam audit on `Adam_QITE/`.

Why this file exists
--------------------
The audit in `../Audits/Audit_Adam_QITE/` was written against twelve scratchpad experiments (V1-V12)
run on the `qml-julia` environment. Those scratchpads no longer exist. Its strongest results --
G2's EMA-horizon mismatch and G5's validation of the free error estimator -- therefore had no
reproducible source in the repository, which is *exactly* the defect the audit itself raises as
G6 section 6.2 ("a correct claim with no code path", a recurrence of the project-wide audit's F3). This
script closes that hole for the four findings that rest on numbers rather than on reading.

It imports `adam_qite` and `adaptive_qite` and modifies neither. Every physical choice --
Hamiltonian, ansatz, finite-difference M and Y, Eq. (7) debiasing, the pseudoinverse solve, the
accept/reject gate, the cond(M) cap, the seeds -- is imported, never redefined. What this script
adds is measurement: it re-runs the shipped code and prints each audit number next to the value
measured today, with the note and section it came from.

What it regenerates
-------------------
  g2     Matched vs mismatched EMA horizons (G2, V9 and V10). The stock (beta1=0.9,
         beta2=0.999) pair stalls at MAX_STEPS; every matched pair finishes at tau_max,
         including beta1=beta2=0.999, the value the first draft of G2 wrongly blamed.
  g5     The O(delta_tau^2) "free error estimate" against R1 step-doubling (G5, V2). The one
         positive finding in the audit, and the one that appears nowhere outside a docstring.
  g6-2   The ratchet degeneracy: beta1=beta2=0, k_I=10 collapses to Route 3 exactly (G6, V5).
  g6-3   `energy_at_budget`'s one-step over-credit, synthetic and on the real E1 arms (G6, V12).

How to run it
-------------
    cd QITE_Hand_Off/Adam_QITE
    ~/micromamba/envs/qml-julia/bin/python verify_audit.py --target all 2>&1 \
      | tee Output_Audit_Verify.txt

`2>&1` is deliberate: the k_I=10 arm overflows a scalar power (harmlessly -- see g6-2) and that
belongs in the artifact rather than swallowed. About 25 s for `--target all`. This script writes
no files of its own and draws no figures.

A drift is a result, not a failure
----------------------------------
`compare()` never raises. It prints MATCH or DRIFT and records the latter in a ledger; the script
exits non-zero if anything drifted. The correct response to a drift is to record the new value and
the environment difference in `../Audits/Audit_Adam_QITE/` and update the note -- never to adjust this
script until it agrees. That rule is the vault's, and it is why there is no tolerance flag here:
tolerances are derived from the recorded numbers' own printed precision and cannot be loosened
from the command line.

The three `check_*` functions are different in kind. They assert hard, because a failure there
means the script is no longer measuring what it claims to measure.

Sources
-------
  [[Audit_Adam_QITE]]                        the hub, findings register in section 4
  [[G2_Beta2_Horizon_Dominates_Every_Verdict]]          target g2
  [[G5_Error_Estimator_Validated_And_Its_Blind_Spot]]   target g5
  [[G6_Attribution_And_Reproducibility]]                targets g6-2 and g6-3
  ../Audits/Audit_Project_Wide/R_Future_Routes.md               R1, the step-doubling route g5 tests
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy.stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "Attack_Route_03"))  # same bootstrap as adam_qite.py:51

from adaptive_qite import (  # noqa: E402
    DE_BAND_ABS, DE_BAND_REL, DELTA_TAU0, DELTA_TAU_MAX, DELTA_TAU_MIN, GROW, H_DERIV,
    MAX_RETRIES, MAX_STEPS, RCOND_DEFAULT, RESIDUAL_REL_TOL, SEED_A, SEED_C, SEED_D, SHRINK,
    STIFFNESS_CAP_COEFF, TAU_MAX, adaptive_qite_run, banner, energy_and_variance, energy_check,
    fixed_step_qite_run, h2_hamiltonian, make_rho_fn, measure_f, variational_floor)

sys.path.insert(0, str(HERE))  # adam_qite lives here; import it by name, do not exec it
from adam_qite import (  # noqa: E402
    BETA1, BETA2, EPS, K_I, adam_qite_run, circuit_cumsum, energy_at_budget)

N_QUBITS = 2  # N in Eq. (7)'s bias factors; the 2-qubit H2 ansatz, as every target uses

ASSUMPTIONS = [
    "Inherits Assumptions 1-8 of adaptive_qite.py and 1-9 of adam_qite.py by construction. "
    "Nothing physical is redefined here; this script only re-runs the shipped code and measures.",

    "G5 needs theta_n along the E1 trajectory, which adam_qite_run does not log. Rather than "
    "re-implement its step loop -- which would verify the audit against a *copy* of the code, "
    "the very defect the audit is about -- theta_n is replayed from the log's own delta_tau plus "
    "measure_f. This is valid only because (i) the committed step uses the final post-cap, "
    "post-retry delta_tau, which is exactly what log['delta_tau'][n] records, and (ii) nothing "
    "but theta and delta_tau carries between steps. check_replay_matches_shipped_log() asserts "
    "both, so a future edit breaking either trips an assertion instead of quietly moving a number.",

    "WHAT THE REPLAY CHECK IS AND IS NOT. It is falsifiable -- a post-hoc correction to theta, a "
    "logged delta_tau that is the pre-cap rather than the committed value, or any hidden state "
    "carried between steps would all make the worst deviation exceed the tolerance. But it is a "
    "plumbing check, not physics validation: both the replay and the shipped run call the same "
    "measure_f and energy_check, so an error inside either is invisible to it by construction. It "
    "licenses target g5's replay premise and nothing more. The non-circular check in this script "
    "is the step-doubling comparison in target g5, which puts the module's error estimate against "
    "an independently computed one; that one has a target it could have missed, and does not.",

    "THE G5 PROBE COUNT OF 40 IS OURS, RECOVERED BY INFERENCE. The audit records twelve printed "
    "rows and five aggregate statistics but never says how many steps V2 probed, and the "
    "scratchpad is gone. 40 accepted steps is the prefix that reproduces all five statistics to "
    "their recorded digits. The statistics are monotone in the count (see the sensitivity line "
    "printed by target g5), so this parameter is load-bearing: it is pinned as the default, not "
    "chosen per run. G5 has been amended to record it.",

    "The rule by which the audit selected twelve of those 40 rows to print cannot be recovered. "
    "This script prints all 40 instead of inventing a selection rule; the audit's twelve are a "
    "subset and can be checked against these by eye.",

    "'Converged tail' in G5's '2.03-2.09' claim is our operational definition: probed steps 25 "
    "through n_probes. The interval is prose-rounded in the note, so the measured [min, max] is "
    "printed rather than a pass/fail on the prose.",

    "e_E_true is |E_full - E_two_half| and e_theta is ||theta_full - theta_two_half||_2. The "
    "half-step midpoint gets a FRESH measure_f -- that is the circuit cost R1 pays and the "
    "reason the estimator under test is free by comparison. The step-doubling arm is a probe "
    "only: it never feeds back into the replayed trajectory.",

    "The corrected energy_at_budget cannot be imported (it is the defect under test) and must "
    "not be patched into adam_qite.py. energy_at_budget_right() is defined locally and imports "
    "circuit_cumsum, so only the single indexing decision is duplicated and the cost model "
    "cannot drift from the shipped one. It is a demonstration, not a fix.",

    "G2's V10 table covers SEED_A and SEED_D. Its prose said 'all three seeds'; SEED_C was swept "
    "but its rows are omitted for the rcond reason G2 gives, and are no longer reproducible. "
    "This script regenerates the two seeds the note actually shows and says so, rather than "
    "fabricating a third row. G2 has been amended to match.",

    "variational_floor carries its own seed=12345 and does not depend on theta0, so the audit's "
    "two identical E_floor = -1.124600 rows are ONE measurement, not two. Computed once here and "
    "printed for both seeds with that stated.",

    "Runs are memoized on (kind, seed, n_params, kwargs), which is safe because each run is a "
    "deterministic function of its inputs. check_cache_determinism() re-runs one configuration "
    "and asserts identity rather than taking that on trust.",

    "adam_signal_diagnostic is deliberately NOT reused for the sweeps: it closes over the "
    "module-level K_I, so its 'ungated k_I' line would mislabel every arm with a swept k_i. "
    "arm_stats() computes only the columns the audit reports. This non-reusability is G2's "
    "proposed-work item 3.",
]

DRIFTS: list[tuple[str, str, str, str]] = []
N_COMPARED = 0

# --------------------------------------------------------------------------------------
# Recorded numbers. Every audit value lives here once, as a string at the audit's own
# printed precision, with the note and section it came from. No finding function below
# hardcodes a target.
# --------------------------------------------------------------------------------------

G2_V9 = "G2 section 'Evidence', table V9"
G2_V10 = "G2 section 'Evidence', table V10"
G5_V2 = "G5 section 'It matches a real error estimator', V2"
G5_SENS = "G5 section 'It matches a real error estimator', paragraph 'Probe count, recorded late'"
G6_V5 = "G6 section 6.2, V5"
G6_V12 = "G6 section 6.3, V12"

RECORDED: dict[str, dict] = {
    # --- G2 / V9: SEED_D horizon sweep ------------------------------------------------
    "g2.v9.route3": dict(steps=234, tau="3.0097", E="-1.123745", exit="tau_max", source=G2_V9),
    "g2.v9.0.9_0.999": dict(steps=500, tau="2.8615", E="-1.123002", s_e="0.0893",
                            frozen=72, exit="max_steps", source=G2_V9),
    "g2.v9.0.999_0.999": dict(steps=254, tau="3.0269", E="-1.123757", s_e="0.5170",
                              frozen=0, exit="tau_max", source=G2_V9),
    "g2.v9.0.99_0.99": dict(steps=259, tau="3.0226", E="-1.123736", s_e="0.4387",
                            frozen=0, exit="tau_max", source=G2_V9),
    "g2.v9.0.9_0.9": dict(steps=244, tau="3.0093", E="-1.123736", s_e="0.3938",
                          frozen=0, exit="tau_max", source=G2_V9),
    "g2.v9.0.999_0.9": dict(steps=242, tau="3.0077", E="-1.123740", s_e="258.1220",
                            frozen=0, exit="tau_max", source=G2_V9),

    # --- G2 / V10: matched-circuit-budget currency, SEED_A then SEED_D ----------------
    "g2.v10.A.route3": dict(steps=172, circ=2069, dE="+1.607e-03", at_budget="+1.607e-03",
                            source=G2_V10),
    "g2.v10.A.0.9_0.999": dict(steps=190, circ=2285, dE="+1.748e-03", at_budget="+1.860e-03",
                               s_e="0.140", frozen=49, source=G2_V10),
    "g2.v10.A.0.9_0.9": dict(steps=179, circ=2151, dE="+1.190e-03", at_budget="+1.800e-03",
                             s_e="0.200", frozen=0, source=G2_V10),
    "g2.v10.A.0.99_0.99": dict(steps=184, circ=2212, dE="+1.437e-03", at_budget="+1.841e-03",
                               s_e="0.355", frozen=0, source=G2_V10),
    "g2.v10.A.0.999_0.999": dict(steps=181, circ=2176, dE="+1.447e-03", at_budget="+1.803e-03",
                                 s_e="0.407", frozen=0, source=G2_V10),
    "g2.v10.D.route3": dict(steps=234, circ=2816, dE="+8.555e-04", at_budget="+8.555e-04",
                            source=G2_V10),
    "g2.v10.D.0.9_0.999": dict(steps=500, circ=6006, dE="+1.598e-03", at_budget="+4.165e-03",
                               s_e="0.089", frozen=72, source=G2_V10),
    "g2.v10.D.0.9_0.9": dict(steps=244, circ=2934, dE="+8.641e-04", at_budget="+3.495e-03",
                             s_e="0.394", frozen=0, source=G2_V10),
    "g2.v10.D.0.99_0.99": dict(steps=259, circ=3114, dE="+8.637e-04", at_budget="+3.942e-03",
                               s_e="0.439", frozen=0, source=G2_V10),
    "g2.v10.D.0.999_0.999": dict(steps=254, circ=3054, dE="+8.430e-04", at_budget="+3.893e-03",
                                 s_e="0.517", frozen=0, source=G2_V10),
    "g2.v10.budget.A": dict(value=2069, source=G2_V10),
    "g2.v10.budget.D": dict(value=2816, source=G2_V10),
    "g2.v10.floor": dict(value="-1.124600", source=G2_V10),

    # --- G5 / V2: the error estimator against R1 step-doubling ------------------------
    "g5.pearson_log": dict(value="0.9902", source=G5_V2),
    "g5.spearman": dict(value="0.9882", source=G5_V2),
    "g5.ratio_median": dict(value="2.060", source=G5_V2),
    "g5.ratio_min": dict(value="0.350", source=G5_V2),
    "g5.ratio_max": dict(value="11.439", source=G5_V2),

    # The probe-count sensitivity G5 quotes to argue that 40 is load-bearing. These were
    # printed but unchecked until a later pass; an independent check found them outside the
    # ledger, which is the same defect class as the missing V9 row above -- a number in the
    # note that RECORDED does not carry is compared against nothing.
    "g5.sens.20.pearson_log": dict(value="0.9701", source=G5_SENS),
    "g5.sens.20.ratio_median": dict(value="2.157", source=G5_SENS),
    "g5.sens.60.pearson_log": dict(value="0.9948", source=G5_SENS),
    "g5.sens.60.ratio_median": dict(value="2.026", source=G5_SENS),
    # Full trajectory: only reachable under --g5-full-trajectory, and skipped loudly otherwise.
    "g5.sens.full.pearson_log": dict(value="0.9983", source=G5_SENS),
    "g5.sens.full.ratio_median": dict(value="2.005", source=G5_SENS),

    # --- G6 section 6.2: the ratchet degeneracy ---------------------------------------
    "g6_2.route3": dict(steps=172, E="-1.122993", source=G6_V5),
    "g6_2.k10": dict(steps=172, E="-1.122993", pinned=100, source=G6_V5),
    "g6_2.k015": dict(steps=172, E="-1.122928", pinned=95, source=G6_V5),
    # The stock row of the same table. Printed but uncompared until a later pass.
    "g6_2.stock": dict(steps=190, E="-1.122852", pinned=5, source=G6_V5),

    # --- G6 section 6.3: energy_at_budget over-credits by one step --------------------
    "g6_3.budget": dict(value=732, source=G6_V12),
    "g6_3.adam": dict(circ=2285, pct=32, afford=60, correct="-1.120983", reported="-1.121027",
                      source=G6_V12),
    "g6_3.route3": dict(circ=2069, pct=35, afford=60, correct="-1.120952",
                        reported="-1.120992", source=G6_V12),
    "g6_3.fixed": dict(circ=732, pct=100, afford=61, correct="-1.116691",
                       reported="-1.116691", source=G6_V12),
}

# The synthetic table of G6 section 6.3: (budget, steps affordable, correct E, reported E).
# n_params = 5 gives per_eval + 1 = 12 circuits per step, so cum = [12, 24, 36, 48, 60].
RECORDED_G6_3_SYNTHETIC = [
    (12, 1, -1.0, -1.0),
    (20, 1, -1.0, -2.0),
    (24, 2, -2.0, -2.0),
    (30, 2, -2.0, -3.0),
    (47, 3, -3.0, -4.0),
    (48, 4, -4.0, -4.0),
]

# --------------------------------------------------------------------------------------
# Infrastructure
# --------------------------------------------------------------------------------------


def print_provenance():
    """What was actually imported, and every upstream constant these numbers depend on.

    G1's failure mode was a stale path silently resolving to a stale __pycache__, so the
    resolved __file__ and mtime of both modules are printed rather than assumed. An edit to an
    upstream constant then shows up as a one-line diff in the artifact instead of as drift
    nobody can explain.
    """
    banner("Provenance -- what this run measured against")
    import adam_qite
    import adaptive_qite
    print("  modules imported (resolved path, mtime):")
    for mod in (adaptive_qite, adam_qite):
        path = Path(mod.__file__).resolve()
        mtime = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
        print(f"    {mod.__name__:16s} {path}")
        print(f"    {'':16s} mtime {mtime:%Y-%m-%d %H:%M:%S} UTC")
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE,
                             capture_output=True, text=True, timeout=10).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", str(HERE)], cwd=HERE,
                               capture_output=True, text=True, timeout=10).stdout.strip()
        print(f"\n  git HEAD {sha or '(unknown)'}"
              f"{'  (Adam_QITE/ has uncommitted changes)' if dirty else ''}")
    except Exception as exc:                                     # noqa: BLE001
        print(f"\n  git HEAD unavailable: {exc!r}")

    import matplotlib
    import pennylane
    import scipy
    print(f"\n  python {sys.version.split()[0]}   numpy {np.__version__}   "
          f"scipy {scipy.__version__}   pennylane {pennylane.__version__}   "
          f"matplotlib {matplotlib.__version__}")

    print("\n  upstream constants these numbers depend on:")
    print(f"    seeds        SEED_A={SEED_A} SEED_C={SEED_C} SEED_D={SEED_D}")
    print(f"    integration  TAU_MAX={TAU_MAX} MAX_STEPS={MAX_STEPS} MAX_RETRIES={MAX_RETRIES}")
    print(f"    step bounds  DELTA_TAU0={DELTA_TAU0} DELTA_TAU_MIN={DELTA_TAU_MIN} "
          f"DELTA_TAU_MAX={DELTA_TAU_MAX} GROW={GROW} SHRINK={SHRINK}")
    print(f"    gate         DE_BAND_ABS={DE_BAND_ABS} DE_BAND_REL={DE_BAND_REL} "
          f"RESIDUAL_REL_TOL={RESIDUAL_REL_TOL}")
    print(f"    linalg       RCOND_DEFAULT={RCOND_DEFAULT} "
          f"STIFFNESS_CAP_COEFF={STIFFNESS_CAP_COEFF} H_DERIV={H_DERIV}")
    print(f"    adam         BETA1={BETA1} BETA2={BETA2} EPS={EPS} K_I={K_I}")


def _tolerance_from(recorded: str) -> float:
    """Half a unit in the recorded value's last printed decimal place.

    '-1.123002' -> 5e-7, '3.0097' -> 5e-5, '2.060' -> 5e-4, '+1.607e-03' -> 5e-7. The audit's
    own precision sets the bar; there is deliberately no way to widen it from the CLI.
    """
    s = recorded.strip().lstrip("+-")
    if "e" in s.lower():
        mantissa, _, exponent = s.lower().partition("e")
        decimals = len(mantissa.partition(".")[2])
        return 0.5 * 10.0 ** (int(exponent) - decimals)
    return 0.5 * 10.0 ** (-len(s.partition(".")[2]))


def compare(label, measured, key, field="value", indent=4):
    """Print measured against recorded and record a drift. Never raises.

    Ints (step counts, percentages, circuit counts) must match exactly -- they are integers
    produced by threshold crossings, so a single flipped accept/reject decision is a real
    finding and should not be absorbed by a tolerance. Floats get _tolerance_from().
    """
    global N_COMPARED
    entry = RECORDED[key]
    recorded = entry[field]
    source = entry.get("source", "")
    N_COMPARED += 1
    pad = " " * indent

    if isinstance(recorded, int):
        ok = int(measured) == recorded
        line = (f"{pad}{label:<26s} recorded={recorded:>10d}  measured={int(measured):>10d}  "
                f"delta={int(measured) - recorded:+d}")
    elif isinstance(recorded, str) and not recorded.lstrip("+-")[:1].isdigit():
        ok = str(measured) == recorded              # exit_reason and other string fields
        line = f"{pad}{label:<26s} recorded={recorded:>10s}  measured={str(measured):>10s}"
    else:
        tol = _tolerance_from(recorded)
        delta = float(measured) - float(recorded)
        ok = abs(delta) <= tol
        line = (f"{pad}{label:<26s} recorded={recorded:>10s}  measured={float(measured):>10.4e}  "
                f"delta={delta:+.1e}  (tol {tol:.1e})")

    print(f"{line}  {'MATCH' if ok else 'DRIFT'}")
    if not ok:
        DRIFTS.append((f"{key}/{field}", label, str(recorded), str(measured)))
        print(f"{pad}  ^ recorded in: {source}")
    return ok


_CACHE: dict[tuple, dict] = {}


def cached_run(kind, seed, n_params=5, **kwargs):
    """Run and memoize. `kind` selects which shipped runner to call; nothing is reimplemented."""
    key = (kind, seed, n_params, tuple(sorted(kwargs.items())))
    if key in _CACHE:
        return _CACHE[key]

    H = h2_hamiltonian()
    rho_fn = make_rho_fn(n_params=n_params)
    theta0 = np.random.default_rng(seed).uniform(-0.5, 0.5, size=n_params)
    runner = {"adam": adam_qite_run, "route3": adaptive_qite_run,
              "fixed": fixed_step_qite_run}[kind]

    t0 = time.perf_counter()
    # k_I = 10 overflows m_e_hat ** (-s_e * k_i) to inf, which np.clip then pulls back to GROW.
    # Numerically harmless -- and it is precisely why the degeneracy in g6-2 is exact. It is
    # COUNTED AND REPORTED rather than silenced: an overflow is a real numerical event, and an
    # artifact that hides it would misrepresent what the shipped controller does at that k_i.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with np.errstate(over="warn"):
            run = runner(theta0, rho_fn, H, n_params, N_QUBITS, **kwargs)
    dt = time.perf_counter() - t0
    overflows = [w for w in caught if "overflow" in str(w.message).lower()]

    shown = " ".join(f"{k}={v}" for k, v in sorted(kwargs.items())) or "(defaults)"
    steps = len(run["tau"]) - 1
    print(f"    [ran] {kind:<7s} seed={seed} n={n_params} {shown:<44s} "
          f"-> {steps:>3d} steps, {dt:4.1f} s")
    if overflows:
        print(f"           ^ {len(overflows)} overflow RuntimeWarning(s) in the growth power law"
              f", clipped to GROW: {str(overflows[0].message).strip()}")
    for w in (w for w in caught if w not in overflows):
        print(f"           ^ {w.category.__name__}: {str(w.message).strip()}")
    _CACHE[key] = run
    return run


def arm_stats(run, n_params=5):
    """Only the columns the audit reports. See ASSUMPTIONS on why adam_signal_diagnostic is
    not reused here."""
    steps = len(run["tau"]) - 1
    out = dict(steps=steps, tau_end=float(run["tau"][-1]), E_final=float(run["E"][-1]),
               exit=run.get("exit_reason", "n/a"),
               circ_total=int(circuit_cumsum(run, n_params)[-1]))
    if "s_e" in run:
        growth = run["growth"][1:]
        out["mean_s_e"] = float(np.mean(run["s_e"][1:]))
        out["frozen_pct"] = int(round(100.0 * np.mean(growth < 1.001)))
        out["pinned_pct"] = int(round(100.0 * np.mean(growth >= GROW - 1e-12)))
    else:
        out["mean_s_e"] = out["frozen_pct"] = out["pinned_pct"] = None
    return out


def energy_at_budget_right(run, n_params, budget):
    """`energy_at_budget` with the indexing corrected. A demonstration, not a patch.

    The shipped version takes searchsorted's default side='left' and then indexes E[i + 1],
    which credits the run with a step it had not finished paying for. cum[j] is the total spent
    *after* step j+1, so the number of affordable steps is searchsorted(cum, budget, 'right')
    and the energy reached is E[that]. Only this one line differs; circuit_cumsum is imported
    so the cost model cannot drift from the shipped one.
    """
    cum = circuit_cumsum(run, n_params)
    affordable = int(np.searchsorted(cum, budget, side="right"))
    return float(run["E"][min(affordable, len(run["E"]) - 1)]), affordable


# --------------------------------------------------------------------------------------
# Structural self-checks. These assert: a failure means the script has stopped measuring
# what it claims to. Same role as the parent's check_eq_b9_self().
# --------------------------------------------------------------------------------------


def check_replay_matches_shipped_log(run, seed, n_params=5, tol=1e-12):
    """Assert theta_n can be replayed from log['delta_tau'] alone, so g5 need not copy the
    step loop. Returns the worst deviation over E, dE_pred and dE_actual."""
    H = h2_hamiltonian()
    rho_fn = make_rho_fn(n_params=n_params)
    theta = np.random.default_rng(seed).uniform(-0.5, 0.5, size=n_params)
    n_steps = len(run["tau"]) - 1
    worst = 0.0
    for k in range(1, n_steps + 1):
        f, _M, Y, rho_n = measure_f(theta, rho_fn, H, N_QUBITS)
        E_n, var_n = energy_and_variance(rho_n, H)
        h = float(run["delta_tau"][k])
        theta_trial = theta + h * f
        E_trial, _ = energy_and_variance(rho_fn(theta_trial), H)
        dE_pred, _dE_exact, dE_actual, _reasons = energy_check(
            E_n, var_n, Y, f, h, E_trial)
        worst = max(worst,
                    abs(E_trial - float(run["E"][k])),
                    abs(dE_pred - float(run["dE_pred"][k])),
                    abs(dE_actual - float(run["dE_actual"][k])))
        theta = theta_trial
    print(f"    check_replay_matches_shipped_log: worst |replay - log| over {n_steps} steps "
          f"= {worst:.2e}  (tol {tol:.0e})")
    assert worst < tol, (
        f"replaying theta_n from log['delta_tau'] no longer reproduces the shipped log "
        f"(worst deviation {worst:.3e}). One of the two invariants in ASSUMPTION 2 has been "
        f"broken by an edit to adam_qite.py; target g5 cannot be trusted until it is fixed.")
    return worst


def check_cache_determinism():
    """Assert a re-run of the same configuration is bit-identical, since cached_run memoizes."""
    H = h2_hamiltonian()
    rho_fn = make_rho_fn(n_params=5)
    theta0 = np.random.default_rng(SEED_A).uniform(-0.5, 0.5, size=5)
    a = adaptive_qite_run(theta0, rho_fn, H, 5, N_QUBITS)
    b = adaptive_qite_run(theta0, rho_fn, H, 5, N_QUBITS)
    worst = float(np.max(np.abs(a["E"] - b["E"])))
    print(f"    check_cache_determinism: worst |E_run1 - E_run2| = {worst:.2e}  (tol 0)")
    assert worst == 0.0, (
        f"two identical runs disagree by {worst:.3e}; memoizing runs is unsafe and every "
        f"cached number below is suspect")
    return worst


def check_synthetic_cost_model(synth, n_params=5):
    """Assert the fabricated log in g6-3 really does cost 12 circuits a step under the shipped
    cost model, so the synthetic demonstration cannot drift if that model changes."""
    cum = circuit_cumsum(synth, n_params)
    expected = np.array([12.0, 24.0, 36.0, 48.0, 60.0])
    print(f"    check_synthetic_cost_model: circuit_cumsum = {list(cum.astype(int))}, "
          f"expected {list(expected.astype(int))}")
    assert np.array_equal(cum, expected), (
        f"the synthetic log no longer costs 12 circuits per step under circuit_cumsum "
        f"(got {list(cum)}); the six recorded rows below are not comparable")
    return 0.0


# --------------------------------------------------------------------------------------
# Finding g2 -- the EMA horizon mismatch
# --------------------------------------------------------------------------------------

HORIZONS = [(0.9, 0.999), (0.999, 0.999), (0.99, 0.99), (0.9, 0.9)]


def finding_g2_matched_ema_horizons():
    banner("g2 -- matched vs mismatched EMA horizons  [[G2_Beta2_Horizon_Dominates_Every_Verdict]]")
    print("""
  beta1 = 0.9 averages the local-error proxy over ~10 steps; beta2 = 0.999 averages its square
  over ~1000, on runs of 172-500 steps. The confidence ratio s_e = m_e_hat / sqrt(v_e_hat)
  divides a ten-step local mean by an RMS that has barely moved off its zero initialization, so
  Jensen's (0, 1] bound -- which needs ONE weighting, not two -- does not apply. This is what
  E2's "the premise fails" verdict was actually measuring.

  V9 sweeps the pair at SEED_D. The test is structural: does the stall track the *mismatch*
  rather than beta2's value?
""".rstrip())

    print("\n  V9 -- SEED_D, 5-param, all else at defaults:\n")
    route3_d = cached_run("route3", SEED_D)
    arms_d = {hz: cached_run("adam", SEED_D, beta1=hz[0], beta2=hz[1]) for hz in HORIZONS}
    reversed_d = cached_run("adam", SEED_D, beta1=0.999, beta2=0.9)

    header = (f"\n  {'configuration':<36s}{'steps':>7s}{'tau_end':>10s}{'final E':>13s}"
              f"{'mean s_e':>11s}{'frozen':>8s}{'exit':>11s}")
    print(header)
    rows = [("Route 3 (ratchet)", route3_d, "g2.v9.route3")]
    for (b1, b2) in HORIZONS:
        tag = "  <-- matched" if b1 == b2 else ""
        rows.append((f"b1={b1} b2={b2}{tag}", arms_d[(b1, b2)], f"g2.v9.{b1}_{b2}"))
    rows.append(("b1=0.999 b2=0.9  (reversed)", reversed_d, "g2.v9.0.999_0.9"))

    stats_d = {}
    for label, run, key in rows:
        st = arm_stats(run)
        stats_d[key] = st
        s_e = "n/a" if st["mean_s_e"] is None else f"{st['mean_s_e']:.4f}"
        froz = "n/a" if st["frozen_pct"] is None else f"{st['frozen_pct']}%"
        print(f"  {label:<36s}{st['steps']:>7d}{st['tau_end']:>10.4f}{st['E_final']:>13.6f}"
              f"{s_e:>11s}{froz:>8s}{st['exit']:>11s}")

    print("\n  against the recorded V9 table:")
    for _label, _run, key in rows:
        st = stats_d[key]
        print(f"    {key}:")
        compare("steps", st["steps"], key, "steps", indent=6)
        compare("tau_end", st["tau_end"], key, "tau", indent=6)
        compare("final E", st["E_final"], key, "E", indent=6)
        compare("exit_reason", st["exit"], key, "exit", indent=6)
        if "s_e" in RECORDED[key]:
            compare("mean s_e", st["mean_s_e"], key, "s_e", indent=6)
            compare("frozen %", st["frozen_pct"], key, "frozen", indent=6)

    print("\n  NOTE: the b1=0.9 b2=0.999 row exits on max_steps, so its tau_end is set by "
          "MAX_STEPS\n        (currently %d), not by convergence. It moves by design if that "
          "constant is edited." % MAX_STEPS)

    # --- V10: the same question in M6's matched-circuit currency ----------------------
    print("\n\n  V10 -- the same arms in M6's matched-circuit currency, two seeds.")
    print("  E_floor is one measurement, not two: variational_floor carries its own seed=12345")
    print("  and does not depend on theta0 (see ASSUMPTIONS).\n")

    H = h2_hamiltonian()
    floor = variational_floor(make_rho_fn(n_params=5), 5, H)
    print(f"    E_floor = {floor:.6f}")
    compare("E_floor", floor, "g2.v10.floor", indent=4)

    claim_budget = {}
    for seed, tag in ((SEED_A, "A"), (SEED_D, "D")):
        print(f"\n  SEED_{tag} 5-param:\n")
        route3 = cached_run("route3", seed)
        arms = {hz: cached_run("adam", seed, beta1=hz[0], beta2=hz[1]) for hz in HORIZONS}

        all_arms = [("Route 3 (ratchet)", route3, f"g2.v10.{tag}.route3")]
        for (b1, b2) in HORIZONS:
            note = "(stock, mismatched)" if (b1, b2) == (0.9, 0.999) else "(matched)"
            all_arms.append((f"Adam b1={b1} b2={b2} {note}", arms[(b1, b2)],
                             f"g2.v10.{tag}.{b1}_{b2}"))

        totals = {lbl: int(circuit_cumsum(r, 5)[-1]) for lbl, r, _k in all_arms}
        budget = min(totals.values())
        setter = min(totals, key=lambda k: totals[k])
        print(f"    matched budget = {budget} circuits, set by: {setter}")
        compare("matched budget", budget, f"g2.v10.budget.{tag}", indent=4)

        print(f"\n  {'arm':<40s}{'steps':>7s}{'circ':>7s}{'E-floor':>12s}"
              f"{'@budget(ship)':>15s}{'@budget(corr)':>15s}{'s_e':>8s}{'froz':>7s}")
        measured = {}
        for label, run, key in all_arms:
            st = arm_stats(run)
            e_ship = energy_at_budget(run, 5, budget)
            e_corr, _afford = energy_at_budget_right(run, 5, budget)
            s_e = "n/a" if st["mean_s_e"] is None else f"{st['mean_s_e']:.3f}"
            froz = "n/a" if st["frozen_pct"] is None else f"{st['frozen_pct']}%"
            print(f"  {label:<40s}{st['steps']:>7d}{st['circ_total']:>7d}"
                  f"{st['E_final'] - floor:>+12.3e}{e_ship - floor:>+15.3e}"
                  f"{e_corr - floor:>+15.3e}{s_e:>8s}{froz:>7s}")
            measured[key] = dict(st=st, e_ship=e_ship, e_corr=e_corr)

        print("\n  against the recorded V10 table:")
        for _label, _run, key in all_arms:
            m = measured[key]
            print(f"    {key}:")
            compare("steps", m["st"]["steps"], key, "steps", indent=6)
            compare("circuits", m["st"]["circ_total"], key, "circ", indent=6)
            compare("E - floor", m["st"]["E_final"] - floor, key, "dE", indent=6)
            # V10's E@budget column was produced by a scratchpad whose reduction convention is
            # unknown, so both are offered and the matching one is named.
            rec = RECORDED[key]["at_budget"]
            tol = _tolerance_from(rec)
            hit = [n for n, v in (("shipped", m["e_ship"]), ("corrected", m["e_corr"]))
                   if abs((v - floor) - float(rec)) <= tol]
            compare("E@budget (shipped)", m["e_ship"] - floor, key, "at_budget", indent=6)
            print(f"      ^ conventions matching the recorded value: "
                  f"{', '.join(hit) if hit else 'NEITHER'}"
                  f"{'  (they agree at this budget -- luck, not design)' if len(hit) == 2 else ''}")
            if "s_e" in rec if isinstance(rec, dict) else "s_e" in RECORDED[key]:
                compare("mean s_e", m["st"]["mean_s_e"], key, "s_e", indent=6)
                compare("frozen %", m["st"]["frozen_pct"], key, "frozen", indent=6)

        r3 = measured[f"g2.v10.{tag}.route3"]["e_corr"]
        beaten = sum(1 for (b1, b2) in HORIZONS
                     if measured[f"g2.v10.{tag}.{b1}_{b2}"]["e_corr"] > r3)
        claim_budget[tag] = (beaten, len(HORIZONS))

    # --- the two structural claims G2 actually rests on -------------------------------
    print("\n\n  The claims, derived rather than quoted:\n")
    matched_ok = sum(1 for (b1, b2) in HORIZONS
                     if b1 == b2 and arms_d[(b1, b2)]["exit_reason"] == "tau_max")
    n_matched = sum(1 for (b1, b2) in HORIZONS if b1 == b2)
    mism = [arms_d[(0.9, 0.999)], reversed_d]
    mism_stalled = sum(1 for r in mism if r["exit_reason"] == "max_steps")
    print("  CLAIM (G2): the stall tracks horizon MISMATCH, not beta2's value.")
    print(f"    matched pairs reaching tau_max:      {matched_ok}/{n_matched}")
    print(f"    b1=0.999 b2=0.999 (the blamed value) reaches tau_max: "
          f"{arms_d[(0.999, 0.999)]['exit_reason'] == 'tau_max'}")
    print(f"    stock mismatched pair stalls on max_steps: "
          f"{arms_d[(0.9, 0.999)]['exit_reason'] == 'max_steps'}")
    holds = matched_ok == n_matched and arms_d[(0.9, 0.999)]["exit_reason"] == "max_steps"
    print(f"    -> {'HOLDS' if holds else 'DOES NOT HOLD'}")
    print(f"    (the reversed pair b1=0.999 b2=0.9 finishes too, but with mean s_e "
          f"{arm_stats(reversed_d)['mean_s_e']:.1f} -- the same mismatch, symptom inverted, "
          f"and\n     a flat refutation of Assumption 4's claimed bound s_e in (0, 1]; "
          f"see G7 section 7.2.)")

    print("\n  CLAIM (G2): at matched circuit budget Route 3 beats every Adam configuration.")
    for tag, (beaten, total) in claim_budget.items():
        print(f"    SEED_{tag}: {beaten}/{total} Adam arms worse at the matched budget")
    all_beaten = all(b == t for b, t in claim_budget.values())
    print(f"    -> {'HOLDS' if all_beaten else 'DOES NOT HOLD'} on the two seeds V10 shows.")
    print("    G2's prose originally said 'all three seeds'; SEED_C was swept but its rows are")
    print("    omitted for the rcond reason G2 gives and are no longer reproducible. G2 has")
    print("    been amended. The rcond-resolved SEED_C sweep remains open work.")
    return dict(matched_complete=holds, route3_wins=all_beaten)


# --------------------------------------------------------------------------------------
# Finding g5 -- the free error estimator against R1 step-doubling
# --------------------------------------------------------------------------------------


def _probe_step_doubling(run, seed, n_probes, n_params=5):
    """Replay the trajectory; at each step also take two half-steps and compare."""
    H = h2_hamiltonian()
    rho_fn = make_rho_fn(n_params=n_params)
    theta = np.random.default_rng(seed).uniform(-0.5, 0.5, size=n_params)
    n_steps = len(run["tau"]) - 1
    n_probes = min(n_probes, n_steps)
    rows = []
    for k in range(1, n_probes + 1):
        f, _M, Y, rho_n = measure_f(theta, rho_fn, H, N_QUBITS)
        E_n, var_n = energy_and_variance(rho_n, H)
        h = float(run["delta_tau"][k])

        theta_full = theta + h * f
        E_full, _ = energy_and_variance(rho_fn(theta_full), H)
        dE_pred, _dE_exact, dE_actual, _reasons = energy_check(E_n, var_n, Y, f, h, E_full)
        e_raw = abs(dE_actual - dE_pred)

        # R1: two half-steps, with a FRESH f at the midpoint. This is the cost R1 pays.
        theta_half = theta + 0.5 * h * f
        f2, _M2, _Y2, _rho2 = measure_f(theta_half, rho_fn, H, N_QUBITS)
        theta_two = theta_half + 0.5 * h * f2
        E_two, _ = energy_and_variance(rho_fn(theta_two), H)

        e_E_true = abs(E_full - E_two)
        rows.append(dict(
            step=k, delta_tau=h, e_hat=float(run["e_hat"][k]), e_raw=e_raw,
            e_E_true=e_E_true, ratio=e_raw / e_E_true if e_E_true > 0 else np.nan,
            e_theta=float(np.linalg.norm(theta_full - theta_two)),
            cond_M=float(run["cond_M"][k]), accepted=bool(run["accepted"][k])))
        theta = theta_full                                  # never the two-half state
    return rows


def _aggregate(rows):
    e_raw = np.array([r["e_raw"] for r in rows])
    e_true = np.array([r["e_E_true"] for r in rows])
    ratio = np.array([r["ratio"] for r in rows])
    ok = (e_raw > 0) & (e_true > 0)
    return dict(
        pearson_log=float(scipy.stats.pearsonr(np.log(e_raw[ok]), np.log(e_true[ok]))[0]),
        spearman=float(scipy.stats.spearmanr(e_raw, e_true)[0]),
        median=float(np.median(ratio)), min=float(np.min(ratio)), max=float(np.max(ratio)))


def finding_g5_error_estimator_vs_step_doubling(n_probes=40, full_trajectory=False):
    banner("g5 -- the free error estimator vs R1 step-doubling  "
           "[[G5_Error_Estimator_Validated_And_Its_Blind_Spot]]")
    print("""
  Y_k = -dE/dtheta_k exactly, so dE_pred = -Y.theta_dot delta_tau is the exact first-order
  Taylor term and |dE_actual - dE_pred| is a genuine O(delta_tau^2) local-error proxy -- already
  computed by energy_check on every step, then discarded. The README asserts this premise and
  never tests it. V2 tests it against the step-doubling estimator R1 proposes, which costs a
  whole extra measure_f per step.

  Step-doubling algebra predicts the ratio exactly: the full step errs by ~C dtau^2 and two half
  steps by ~2C(dtau/2)^2 = C dtau^2 / 2, so their difference is C dtau^2 / 2 -- half of what
  |dE_actual - dE_pred| ~ C dtau^2 measures. The prediction is a constant ratio of 2.
""".rstrip())

    run = cached_run("adam", SEED_A)              # the E1 trajectory
    print()
    check_replay_matches_shipped_log(run, SEED_A)

    if n_probes != 40:
        print(f"\n  WARNING: n_probes={n_probes} is not the audit-matching 40. The recorded G5")
        print("  statistics do not apply at this count and the comparisons below are")
        print("  informational only (they are monotone in the count -- see the sensitivity line).")

    rows = _probe_step_doubling(run, SEED_A, n_probes)
    print(f"\n  All {len(rows)} probed steps (the audit printed twelve of these; its selection")
    print("  rule is unrecoverable, so all are shown -- see ASSUMPTIONS):\n")
    print(f"  {'step':>5s}{'delta_tau':>11s}{'e_hat':>12s}{'|dEa-dEp|':>12s}{'e_E_true':>12s}"
          f"{'ratio':>9s}{'e_theta':>12s}{'cond(M)':>11s}{'acc':>5s}")
    for r in rows:
        print(f"  {r['step']:>5d}{r['delta_tau']:>11.5f}{r['e_hat']:>12.3e}{r['e_raw']:>12.3e}"
              f"{r['e_E_true']:>12.3e}{r['ratio']:>9.3f}{r['e_theta']:>12.3e}"
              f"{r['cond_M']:>11.3g}{'y' if r['accepted'] else 'n':>5s}")

    agg = _aggregate(rows)
    print(f"\n  Pearson  corr(log e_raw, log e_E_true) = {agg['pearson_log']:.4f}")
    print(f"  Spearman corr(e_raw, e_E_true)         = {agg['spearman']:.4f}")
    print(f"  ratio |dEa-dEp| / e_E_true: median={agg['median']:.3f} "
          f"min={agg['min']:.3f} max={agg['max']:.3f}")

    print("\n  against the recorded V2 aggregates:")
    compare("Pearson (log)", agg["pearson_log"], "g5.pearson_log")
    compare("Spearman", agg["spearman"], "g5.spearman")
    compare("ratio median", agg["median"], "g5.ratio_median")
    compare("ratio min", agg["min"], "g5.ratio_min")
    compare("ratio max", agg["max"], "g5.ratio_max")

    tail = [r["ratio"] for r in rows if r["step"] >= 25]
    if tail:
        print(f"\n  Converged tail (probed steps 25-{rows[-1]['step']}): "
              f"ratio in [{min(tail):.3f}, {max(tail):.3f}], median {np.median(tail):.3f}")
        print("    G5 records 'settles at 2.03-2.09'; that interval is prose-rounded, and the")
        print("    tail minimum in its own printed table is 2.027. The measured window is shown")
        print("    rather than scored against the prose.")
        print(f"    Predicted by step-doubling algebra: exactly 2. Tail median is "
              f"{100 * abs(np.median(tail) / 2 - 1):.1f}% off.")

    # G5's falsifiable prediction, which V2 stated but never tested.
    lr = np.array([abs(np.log(r["ratio"] / 2.0)) for r in rows])
    cond = np.array([r["cond_M"] for r in rows])
    dlogcond = np.abs(np.diff(np.log(cond), prepend=np.log(cond[0])))
    idx = np.array([r["step"] for r in rows], dtype=float)
    rho_change = float(scipy.stats.spearmanr(lr, dlogcond)[0])
    rho_level = float(scipy.stats.spearmanr(lr, cond)[0])
    rho_index = float(scipy.stats.spearmanr(lr, idx)[0])
    rho_cond_index = float(scipy.stats.spearmanr(cond, idx)[0])
    print("\n  G5's own falsifiable prediction, now tested (no recorded target -- V2 stated the")
    print("  clustering but never computed the correlation). If the blind spot is d_theta M+,")
    print("  the ratio's departure from 2 should track the step-to-step CHANGE in cond(M)")
    print("  rather than its level:")
    print(f"    Spearman(|log(ratio/2)|, |dlog cond(M)|) = {rho_change:+.4f}   <- change")
    print(f"    Spearman(|log(ratio/2)|, cond(M))        = {rho_level:+.4f}   <- level")
    print("\n  Both are confounded, and the confound has to be named before either is read as a")
    print("  mechanism: cond(M) rises almost monotonically along the run while the ratio")
    print("  converges to 2, so step index alone predicts both.")
    print(f"    Spearman(|log(ratio/2)|, step index)     = {rho_index:+.4f}   <- the confound")
    print(f"    Spearman(cond(M), step index)            = {rho_cond_index:+.4f}")
    # Partial correlations, controlling for step index by ranking out the confound.
    def _partial(a, b, c):
        ra, rb, rc = (scipy.stats.rankdata(v) for v in (a, b, c))
        rac = scipy.stats.pearsonr(ra, rc)[0]
        rbc = scipy.stats.pearsonr(rb, rc)[0]
        rab = scipy.stats.pearsonr(ra, rb)[0]
        denom = np.sqrt(max((1 - rac**2) * (1 - rbc**2), 1e-30))
        return float((rab - rac * rbc) / denom)
    p_change = _partial(lr, dlogcond, idx)
    p_level = _partial(lr, cond, idx)
    print("\n  Controlling for step index (rank partial correlation):")
    print(f"    partial(|log(ratio/2)|, |dlog cond(M)| | step) = {p_change:+.4f}   <- change")
    print(f"    partial(|log(ratio/2)|, cond(M)        | step) = {p_level:+.4f}   <- level")
    if abs(p_change) > abs(p_level):
        verdict = ("SUPPORTED once the step-index confound is removed: the CHANGE in cond(M) "
                   "carries\n       the residual signal, which is what the d_theta M+ blind spot "
                   "predicts.")
    else:
        verdict = ("NOT SUPPORTED: with step index controlled, the LEVEL of cond(M) still "
                   "explains at\n       least as much as its change, so the clustering G5 "
                   "observed is not evidence for the\n       d_theta M+ mechanism specifically. "
                   "G5's proposed probe needs a different design.")
    print(f"    -> {verdict}")
    print("    Either way this is a NEW result: the raw correlations alone cannot decide it, and")
    print("    reporting them without the confound would have been the wrong answer.")

    print(f"\n  Provenance: n_probes={n_probes} was recovered by inference, not recorded by the")
    print("  audit (see ASSUMPTIONS). The statistics are monotone in the count, so this matters.")
    print("  G5 quotes the sensitivity to make that argument, so it is checked like any other")
    print("  recorded number -- it was printed but unchecked until a later pass:\n")
    for n in (20, 60):
        if n <= len(run["tau"]) - 1:
            a = _aggregate(rows[:n] if n <= len(rows)
                           else _probe_step_doubling(run, SEED_A, n))
            compare(f"prefix {n} Pearson (log)", a["pearson_log"], f"g5.sens.{n}.pearson_log")
            compare(f"prefix {n} ratio median", a["median"], f"g5.sens.{n}.ratio_median")
    print("\n  The trend favours the estimator: the more of the converged tail is included, the")
    print("  closer the ratio sits to the predicted 2.")

    if full_trajectory:
        n_all = len(run["tau"]) - 1
        print(f"\n  --g5-full-trajectory: all {n_all} accepted steps. This is the strongest form")
        print("  of the physics claim, and the end point of the sensitivity trend G5 quotes.\n")
        rows_all = _probe_step_doubling(run, SEED_A, n_all)
        a = _aggregate(rows_all)
        compare("full Pearson (log)", a["pearson_log"], "g5.sens.full.pearson_log")
        compare("full ratio median", a["median"], "g5.sens.full.ratio_median")
        tail_all = [r["ratio"] for r in rows_all if r["step"] >= 90]
        print(f"    Pearson {a['pearson_log']:.4f}  Spearman {a['spearman']:.4f}  "
              f"median {a['median']:.4f}  range [{a['min']:.3f}, {a['max']:.3f}]")
        if tail_all:
            print(f"    tail (steps 90+): median {np.median(tail_all):.4f}, "
                  f"range [{min(tail_all):.4f}, {max(tail_all):.4f}] "
                  f"-- {100 * abs(np.median(tail_all) / 2 - 1):.2f}% off the predicted 2")
    else:
        print("\n  NOT CHECKED in this run: g5.sens.full.pearson_log (0.9983) and")
        print("  g5.sens.full.ratio_median (2.005), the third point of the sensitivity trend.")
        print("  They need all 190 steps probed; pass --g5-full-trajectory to include them.")
        print("  The ledger total below is 2 lower than it would otherwise be.")
    return agg


# --------------------------------------------------------------------------------------
# Finding g6-2 -- the ratchet degeneracy
# --------------------------------------------------------------------------------------


def finding_g6_2_ratchet_degeneracy():
    banner("g6-2 -- the ratchet degeneracy  [[G6_Attribution_And_Reproducibility]] section 6.2")
    print("""
  G6 section 6.2 records that beta1 = beta2 = 0 with a large k_I collapses the Adam controller
  onto Route 3's plain ratchet exactly, and notes this claim had no code path -- the same defect
  as the project-wide audit's F3. It does now.

  Why it is exact, in one line: with beta1 = 0 the EMA keeps only the current sample, and the
  bias correction c1 = 1 - 0**n = 1 is a no-op, so m_e_hat = e_hat and s_e = 1. growth then
  reduces to e_hat ** -k_I, which for k_I large and e_hat < 1 (always -- see G4) overflows past
  GROW and is clipped to it on every step. Clipped-to-GROW on every accepted step IS the ratchet.
""".rstrip())

    print()
    check_cache_determinism()
    print()
    route3 = cached_run("route3", SEED_A)
    k10 = cached_run("adam", SEED_A, beta1=0.0, beta2=0.0, k_i=10.0)
    k015 = cached_run("adam", SEED_A, beta1=0.0, beta2=0.0, k_i=0.15)
    stock = cached_run("adam", SEED_A)

    print(f"\n  {'configuration':<40s}{'steps':>7s}{'final E':>13s}{'pinned@GROW':>14s}")
    rows = [("Route 3 (ratchet)", route3, "g6_2.route3"),
            ("Adam b1=b2=0, k_I=10  (degenerate)", k10, "g6_2.k10"),
            ("Adam b1=b2=0, k_I=0.15", k015, "g6_2.k015"),
            ("Adam stock b1=0.9 b2=0.999", stock, "g6_2.stock")]
    for label, run, _key in rows:
        st = arm_stats(run)
        pinned = "n/a" if st["pinned_pct"] is None else f"{st['pinned_pct']}%"
        print(f"  {label:<40s}{st['steps']:>7d}{st['E_final']:>13.6f}{pinned:>14s}")

    print("\n  against the recorded section 6.2 numbers:")
    for _label, run, key in rows:
        st = arm_stats(run)
        print(f"    {key}:")
        compare("steps", st["steps"], key, "steps", indent=6)
        compare("final E", st["E_final"], key, "E", indent=6)
        if "pinned" in RECORDED[key]:
            compare("pinned at GROW %", st["pinned_pct"], key, "pinned", indent=6)

    # This is the one physics claim in the audit that is an identity, so it is asserted.
    print("\n  IDENTITY: b1=b2=0, k_I=10 must reproduce Route 3 step for step, not merely match")
    print("  its two reported scalars. Comparing the whole trajectory:\n")
    n = min(len(route3["E"]), len(k10["E"]))
    worst_E = float(np.max(np.abs(route3["E"][:n] - k10["E"][:n])))
    worst_dt = float(np.max(np.abs(route3["delta_tau"][:n] - k10["delta_tau"][:n])))
    print(f"    steps             Route 3 = {len(route3['tau']) - 1}   "
          f"degenerate = {len(k10['tau']) - 1}")
    print(f"    max |E_n difference|        over {n} logged points = {worst_E:.2e}")
    print(f"    max |delta_tau_n difference| over {n} logged points = {worst_dt:.2e}")
    assert len(route3["tau"]) == len(k10["tau"]), (
        "the degenerate configuration no longer takes the same number of steps as Route 3; "
        "G6 section 6.2's collapse claim has been broken by an edit")
    assert worst_E < 1e-12 and worst_dt < 1e-12, (
        f"the degenerate configuration no longer reproduces Route 3 exactly "
        f"(worst |dE| = {worst_E:.3e}, worst |d delta_tau| = {worst_dt:.3e})")
    print("    -> IDENTITY HOLDS, element-wise. This is stronger than the two scalars the note")
    print("       records, and it cost nothing since both logs were already in hand.")

    print("\n  NOTE: at k_I = 10, m_e_hat ** (-s_e * k_i) overflows to inf and np.clip pulls it")
    print("  back to GROW = %.2f (adam_qite.py:212 and :218). Numerically harmless, and it is"
          % GROW)
    print("  precisely why the collapse is exact. numpy's overflow warnings are NOT suppressed:")
    print("  cached_run captures and counts them, so the '^ N overflow RuntimeWarning(s)' line")
    print("  above the k_I=10 arm is the measured count, not a silenced event.")
    return dict(worst_E=worst_E, worst_delta_tau=worst_dt)


# --------------------------------------------------------------------------------------
# Finding g6-3 -- energy_at_budget's one-step over-credit
# --------------------------------------------------------------------------------------


def finding_g6_3_energy_at_budget_off_by_one():
    banner("g6-3 -- energy_at_budget over-credits by one step  "
           "[[G6_Attribution_And_Reproducibility]] section 6.3")
    print("""
  energy_at_budget (adam_qite.py:342) takes searchsorted's default side='left' and then indexes
  E[i + 1]. cum[j] is the total already spent after step j+1, so the correct count of affordable
  steps is searchsorted(cum, budget, side='right') and the energy reached is E[that]. The
  shipped version credits each arm with one step it had not finished paying for.

  This matters because energy_at_budget is the instrument that decides E1's verdict.
""".rstrip())

    # --- Part A: synthetic, where the right answer is checkable by hand ---------------
    print("\n  Part A -- synthetic run, 5 steps, no retries. n_params=5 gives 2N+1+1 = 12")
    print("  circuits a step, so cum = [12, 24, 36, 48, 60] and every answer is checkable.\n")
    synth = dict(tau=np.arange(6.0), E=-np.arange(6.0), retries=np.zeros(6))
    check_synthetic_cost_model(synth)

    print(f"\n  {'budget':>8s}{'affordable':>12s}{'correct E':>12s}{'reported E':>12s}{'':>4s}")
    for budget, rec_afford, rec_correct, rec_reported in RECORDED_G6_3_SYNTHETIC:
        e_corr, afford = energy_at_budget_right(synth, 5, budget)
        e_ship = energy_at_budget(synth, 5, budget)
        flag = "ok" if abs(e_ship - e_corr) < 1e-12 else "OFF BY ONE STEP"
        print(f"  {budget:>8d}{afford:>12d}{e_corr:>12.1f}{e_ship:>12.1f}   {flag}")
        for label, meas, rec in (("affordable", afford, rec_afford),
                                 ("correct E", e_corr, rec_correct),
                                 ("reported E", e_ship, rec_reported)):
            ok = (int(meas) == rec) if isinstance(rec, int) else abs(meas - rec) < 1e-12
            global N_COMPARED
            N_COMPARED += 1
            if not ok:
                DRIFTS.append((f"g6_3.synthetic.budget{budget}", label, str(rec), str(meas)))
                print(f"      DRIFT {label}: recorded={rec} measured={meas}   "
                      f"(recorded in {RECORDED['g6_3.budget']['source']})")
    print("    all six rows compared against the recorded synthetic table: 6 rows x 3 fields")
    print("    = 18 comparisons, counted in the drift ledger but printed only on drift, so a")
    print("    silent table above is the pass and the ledger total stays reconcilable.")

    # --- Part B: the real E1 arms ----------------------------------------------------
    print("\n\n  Part B -- the real E1 arms at SEED_A, in M6's currency.\n")
    adam = cached_run("adam", SEED_A)
    route3 = cached_run("route3", SEED_A)
    fixed = cached_run("fixed", SEED_A)

    arms = [("Adam", adam, "g6_3.adam"), ("Route 3", route3, "g6_3.route3"),
            ("fixed", fixed, "g6_3.fixed")]
    totals = {lbl: int(circuit_cumsum(r, 5)[-1]) for lbl, r, _k in arms}
    budget = min(totals.values())
    setter = min(totals, key=lambda k: totals[k])
    print(f"    matched budget = {budget} circuits, set by the {setter} arm "
          f"(this is E1's reduction, verbatim)")
    compare("matched budget", budget, "g6_3.budget", indent=4)

    print(f"\n  {'arm':<10s}{'total':>8s}{'% of run':>10s}{'afford':>8s}"
          f"{'correct E':>13s}{'reported E':>13s}{'':>4s}")
    got = {}
    for label, run, key in arms:
        total = totals[label]
        e_corr, afford = energy_at_budget_right(run, 5, budget)
        e_ship = energy_at_budget(run, 5, budget)
        pct = int(round(100.0 * budget / total))
        flag = "ok" if abs(e_ship - e_corr) < 1e-12 else "OVER-CREDITED"
        print(f"  {label:<10s}{total:>8d}{pct:>9d}%{afford:>8d}{e_corr:>13.6f}"
              f"{e_ship:>13.6f}   {flag}")
        got[key] = dict(total=total, pct=pct, afford=afford, corr=e_corr, ship=e_ship)

    print("\n  against the recorded section 6.3 table:")
    for _label, _run, key in arms:
        g = got[key]
        print(f"    {key}:")
        compare("total circuits", g["total"], key, "circ", indent=6)
        compare("% of run", g["pct"], key, "pct", indent=6)
        compare("steps affordable", g["afford"], key, "afford", indent=6)
        compare("correct E", g["corr"], key, "correct", indent=6)
        compare("reported E", g["ship"], key, "reported", indent=6)

    # --- the two consequences, computed rather than quoted ----------------------------
    a, r = got["g6_3.adam"], got["g6_3.route3"]
    print("\n\n  Consequence 1 -- does the off-by-one change E1's verdict?\n")
    m_ship = r["ship"] - a["ship"]
    m_corr = r["corr"] - a["corr"]
    print(f"    Route 3 minus Adam at budget, as reported = {m_ship:+.1e}")
    print(f"    Route 3 minus Adam at budget, corrected   = {m_corr:+.1e}")
    same = np.sign(m_ship) == np.sign(m_corr)
    print(f"    sign preserved: {same} -> E1's verdict is {'unchanged' if same else 'CHANGED'}")
    gap = abs(got["g6_3.fixed"]["corr"] - r["corr"])
    print(f"    adaptive-vs-fixed gap = {gap:.1e}, {gap / max(abs(m_corr - m_ship), 1e-30):.0f}x")
    print("    the size of the correction -- the adaptive-beats-fixed conclusion is not")
    print("    threatened by this defect. The defect is real; its blast radius is small.")

    print("\n  Consequence 2 -- `min` is the wrong reduction, which matters more.\n")
    print(f"    {budget} circuits is {a['pct']}% of Adam's run and {r['pct']}% of Route 3's.")
    print("    The budget is set by the handicapped fixed-step arm, so the two adaptive arms")
    print("    are compared one third of the way through both -- not at convergence.")
    adaptive_budget = min(a["total"], r["total"])
    ea, _ = energy_at_budget_right(adam, 5, adaptive_budget)
    er, _ = energy_at_budget_right(route3, 5, adaptive_budget)
    print(f"\n    Budget over the ADAPTIVE arms only = {adaptive_budget} circuits:")
    print(f"      Adam    E = {ea:.6f}")
    print(f"      Route 3 E = {er:.6f}   <- {'Route 3' if er < ea else 'Adam'} ahead")
    print("\n    Both fixes are demonstrated here and NEITHER is applied to adam_qite.py. This")
    print("    script's job is to make the defect reproducible against the shipped code.")
    return dict(margin_reported=m_ship, margin_corrected=m_corr, verdict_unchanged=bool(same))


# --------------------------------------------------------------------------------------
# Tail
# --------------------------------------------------------------------------------------

TARGETS = {
    "g2": finding_g2_matched_ema_horizons,
    "g5": finding_g5_error_estimator_vs_step_doubling,
    "g6-2": finding_g6_2_ratchet_degeneracy,
    "g6-3": finding_g6_3_energy_at_budget_off_by_one,
}


def print_assumptions():
    banner("Assumptions of THIS script -- choices that are ours, not the audit's")
    for i, line in enumerate(ASSUMPTIONS, 1):
        print(f"  {i:>2d}. {line}\n")


NOT_COVERED = [
    "G6 section 6.3's rcond table (four arms x two rcond values, E to six digits). The cap-off "
    "arms at rcond=1e-2 are the evidence for re-attributing E3's win to the cond(M) cap, which "
    "is the audit's third refutation -- and it is the largest block of audit numbers this "
    "script does not regenerate. Reproducing it needs a second sweep dimension; not written.",

    "The s_e maximum of 1.015, which falsifies the claimed s_e in (0,1] bound. It is measured, "
    "but by the shipped suite: it lives in Output_Post_Audit.txt (adam_qite.py --target all), "
    "not in RECORDED. So the audit's first refutation is reproducible from the repo but NOT "
    "from this script's ledger, and a drift in it would not appear in the count below.",

    "E1's and E2's headline numbers generally (190 steps / 2285 circuits / E = -1.122852, the "
    "E2 depolarizing ratios 1.0243-1.0088). Same reason: the shipped suite prints them and this "
    "script only re-checks the four findings whose evidence came from a scratchpad. Anything "
    "adam_qite.py already regenerates on its own was left to it deliberately.",

    "G5's '2.03-2.09 tail' and G6 section 6.1's prose intervals. These are prose-rounded, so a "
    "pass/fail against them would be a test of the rounding. The measured interval is printed "
    "instead; the known near-miss (tail min ~2.027, just outside the prose) is reported, not fixed.",

    "G7 section 7.2's figure clipping (set_ylim(0, 1.05) hides the s_e excursion). It is a "
    "plotting defect with no number to drift.",

    "Everything outside the audit's numeric claims: shot noise, beta1 sweeps, systems beyond "
    "2-qubit H2, seeds beyond three, E3's n=1. Listed in the report as what was not checked.",
]


def print_coverage():
    banner("What this ledger does NOT cover")
    print("\n  The drift count below is a lower bound on the audit's numeric content, not a full")
    print("  accounting of it. A number absent from RECORDED is compared against nothing, and")
    print("  two such numbers have already been found this way -- the V9 0.99/0.99 row (which")
    print("  crashed) and G5's probe-count sensitivity (which printed silently). So the gaps are")
    print("  named here rather than left to be inferred from the key list:\n")
    for i, line in enumerate(NOT_COVERED, 1):
        print(f"  {i:>2d}. {line}\n")


def print_drift_ledger():
    banner("Drift ledger")
    if not DRIFTS:
        print(f"\n  0 of {N_COMPARED} recorded numbers drifted. Adam audit reproduces.\n")
        return 0
    print(f"\n  {len(DRIFTS)} of {N_COMPARED} recorded numbers drifted:\n")
    for key, label, rec, meas in DRIFTS:
        print(f"    {key:<34s} {label:<22s} recorded={rec:<14s} measured={meas}")
    print("""
  THESE ARE RESULTS, NOT FAILURES. Record the new values and the environment difference in
  ../Audits/Audit_Adam_QITE/, update the note, and do NOT adjust this script until it agrees --
  that is the vault's rule and the reason this script has no tolerance flag.

  Reading the ledger: a drifted step count with intact energies is an environment note, not a
  physics result. Step counts are integers set by threshold crossings, so a library change in
  the 15th digit can flip one accept/reject decision and cascade. Energies sit on a converged
  plateau and are far more robust. Check the provenance header above before concluding anything.
""".rstrip())
    return 1


def print_recorded_ledger():
    banner("Recorded numbers this script checks, and where each comes from")
    for key in sorted(RECORDED):
        entry = dict(RECORDED[key])
        source = entry.pop("source", "")
        fields = "  ".join(f"{k}={v}" for k, v in entry.items())
        print(f"  {key:<28s} {fields}")
        print(f"  {'':<28s} <- {source}")
    print(f"\n  plus {len(RECORDED_G6_3_SYNTHETIC)} rows of G6 section 6.3's synthetic table.")


def main():
    ap = argparse.ArgumentParser(
        description="Regenerate the numerical findings of the Adam audit on Adam_QITE/.")
    ap.add_argument("--target", default="all", choices=[*TARGETS, "all"])
    ap.add_argument("--g5-probes", type=int, default=40,
                    help="accepted steps probed by target g5 (default 40, the audit-matching "
                         "count recovered by inference; see ASSUMPTIONS)")
    ap.add_argument("--g5-full-trajectory", action="store_true",
                    help="also probe every accepted step, to expose the probe-count dependence")
    ap.add_argument("--list-targets", action="store_true",
                    help="print the recorded-number ledger and exit, running no physics")
    args = ap.parse_args()

    if args.list_targets:
        print_recorded_ledger()
        return 0

    t0 = time.perf_counter()
    print_provenance()

    names = list(TARGETS) if args.target == "all" else [args.target]
    for name in names:
        if name == "g5":
            TARGETS[name](n_probes=args.g5_probes,
                          full_trajectory=args.g5_full_trajectory)
        else:
            TARGETS[name]()

    print_assumptions()
    if args.target == "all":
        print_coverage()
    code = print_drift_ledger()
    print(f"  total wall time {time.perf_counter() - t0:.1f} s")
    return code


if __name__ == "__main__":
    sys.exit(main())
