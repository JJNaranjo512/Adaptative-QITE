---
title: "G1 — The script does not run from its own directory"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]"
---

# G1 — Broken import and stale paths

**Severity: critical.** Closes G1. **Fixed in this pass.**

## The defect

`adam_qite.py:50-57`, as written:

```python
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
...
from adaptive_qite import (COLORS, DE_BAND_ABS, ...)
```

Commit `3aff9cc move: Adam QITE up one level` moved the folder from
`Adaptative_QITE/Attack_Route_03/Adam_QITE/` to `Adaptative_QITE/Adam_QITE/`. `HERE.parent` is now
`Adaptative_QITE/`, but `adaptive_qite.py` lives in `Adaptative_QITE/Attack_Route_03/`. Observed:

```
$ ~/micromamba/envs/qml-julia/bin/python adam_qite.py --target all
ModuleNotFoundError: No module named 'adaptive_qite'
```

The `__pycache__/adam_qite.cpython-311.pyc` left in the folder is from before the move, which is
why nothing flagged it at the time.

Every `../`-relative reference moved with it and none were updated. The module docstring, the
`ASSUMPTIONS` strings, several `print` calls inside the targets, `Adam_QITE/README.md` and the
`code:` frontmatter field of [[QITE_Hand_Off/Notes_on_Project/Adam_Like_Step_Controller_For_Variational_QITE|Adam_Like_Step_Controller_For_Variational_QITE]] all still describe
the pre-move layout. `R_Future_Routes.md` additionally moved into `Audits/Audit_Project_Wide/` when the
project-wide audit was filed, so `../../Audits/R_Future_Routes.md` was wrong on two counts.

## Why it matters more than a broken path usually would

The vault's rules require that a reported number be reproducible from the source. As found, **no
number in the README or in [[QITE_Hand_Off/Notes_on_Project/Adam_Like_Step_Controller_For_Variational_QITE|Adam_Like_Step_Controller_For_Variational_QITE]] could be regenerated
at all** — the script could not start. A reader checking Juan's work would have hit an import
error before reaching any physics. This is the same class as the project-wide audit's F3 ("a reported per-step
number has no code path"), arrived at by a different route, and it is the reason G1 is sequenced
first: nothing else in this audit could be verified until it was worked around.

## The fix applied

`adam_qite.py:51` — point the insert at the sibling directory that actually holds the parent:

```python
sys.path.insert(0, str(HERE.parent / "Attack_Route_03"))
```

and the prose paths, in `adam_qite.py` (docstring, `ASSUMPTIONS`, target `print`s),
`Adam_QITE/README.md`, and the note's `code:` field:

| Was | Now |
|---|---|
| `../adaptive_qite.py` | `../Attack_Route_03/adaptive_qite.py` |
| `../../Research_Papers/QITE_timestep.pdf` | `../Research_Papers/QITE_timestep.pdf` |
| `../../Audits/R_Future_Routes.md` | `../Audits/Audit_Project_Wide/R_Future_Routes.md` |
| `../../Notes_on_Project/Adam_Like_...md` | `../Notes_on_Project/Adam_Like_...md` |
| `03_LG_Project/.../Attack_Route_03/Adam_QITE/adam_qite.py` | `03_LG_Project/.../Adam_QITE/adam_qite.py` |

`Output_Pre_Audit.txt` keeps its old paths deliberately — it is the record of a run that really
did happen in the old location.

## Verification

Run from the script's own directory, no `PYTHONPATH`:

```bash
cd QITE_Hand_Off/Adam_QITE
~/micromamba/envs/qml-julia/bin/python adam_qite.py --target all --no-plots
```

Exit code 0. Diffed against `Output_Pre_Audit.txt` with the path strings normalized, **every
numeric line is identical** — E1 190 steps / 2285 circuits / $E=-1.122852$, the Route-3 arm
172 / 2069 / $-1.122993$, the E1 diagnostics, E2's four RMS rows, E3's eight arms, the assumption
block. So the move broke the import and **nothing else**; the saved output remains a valid record
of the physics.

Three cosmetic differences remain in the diff and are expected:

- the `wrote …/figures/*.png` lines, absent under `--no-plots`;
- the corrected paths themselves, which appear in the E1 and E3 banners;
- the E1 and E3 banner prose, **rewrapped**. The longer path pushed three source lines past the
  file's 100-column convention, so the surrounding `print` calls were reflowed. The words are
  unchanged and the line breaks moved. Verbatim prose agreement with `Output_Pre_Audit.txt` was
  already impossible once the paths changed, so nothing was lost by reflowing.

Runtime for `--target all` is **73 s** measured (E2 dominates at ~65 s; E1 8 s, E3 6 s run
separately), not the ~3 minutes the README claimed. Minor, but it is a reported number and it was
wrong; the README now says ~75 s.
