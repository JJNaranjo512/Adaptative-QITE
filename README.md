---
title: "Start here — Adaptative_QITE hand-off"
tags:
  - handoff
  - physics
  - quantum-computation
  - adaptative-qite
---

# Start here — Adaptative_QITE hand-off

This folder is a self-contained snapshot of the `Adaptative_QITE` project: everything needed to
understand what was done, check that it is right, and continue it. The
order of the work is given by **phases** (below), and every file is reached either by an Obsidian
link or by its path relative to `Adaptative-QITE/`.

> [!note] Link convention in this folder
> Links are written with the path from the folder name, e.g.
> `[[Adaptative-QITE/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical resume]]`.
> Figure embeds use the same form.

## The project in one paragraph

Variational QITE integrates the linear system $M(\theta)\,\dot\theta = Y(\theta)$ from McLachlan's
principle. Each evaluation of $\dot\theta$ costs $O(N^2)$ circuits for $N$ parameters, so the step
size $\Delta\tau$ sets the circuit bill directly, and the quality of the solve depends on
$\mathrm{cond}(M)$ and on the pseudo-inverse cutoff `rcond`. The project starts from an internal
memo ([[Adaptative-QITE/Research_Papers/QITE_timestep.pdf|QITE_timestep.pdf]], no arXiv/DOI) that treats this as an ODE for an
RKF-type controller with a hardware-noise floor, and asks how far variational QITE can go under
noise and whether adaptive time-stepping helps. All numerical work is on H$_2$ in STO-3G on two
qubits (5- and 9-parameter ansätze). Work done with the LG Electronics Toronto AI team (Thi Ha
Kyaw, Jack Baker, Gaurav Saxena).

## Read in this order

1. **The report** — `Presentations/Adaptative_QITE_Status_Report.pdf` (source:
   `Presentations/Adaptative_QITE_Status_Report.tex`). Five phases in paper format, a claims
   ledger (Table XVII), derivations (Appendix B) and reviewer Q&A (Appendix C). Anticipated questions
   with longer answers: [[Adaptative-QITE/Presentations/Questions_From_Final_Report|Questions_From_Final_Report]].
2. **The system of record** — [[Adaptative-QITE/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]]: framework, error bounds,
   benchmark, roadmap, open research questions.
3. **The phase notes**, in the order of the table below.
4. **The two audits** — [[Adaptative-QITE/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]] (project-wide, findings F1–F12, fix modules M1–M6, routes
   in [[Adaptative-QITE/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]) and [[Adaptative-QITE/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]] (findings G1–G7).
5. **The plan and the P1/P2 results** — [[Adaptative-QITE/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]], then P1 and P2 findings, then the
   deck [[Adaptative-QITE/Presentations/Adaptative_QITE_Beyond_Trot_P1_P2|Adaptative_QITE_Beyond_Trot_P1_P2]] and its
   [[Adaptative-QITE/Presentations/Presentator_Notes_Beyond_Trot_P1_P2|Presentator_Notes_Beyond_Trot_P1_P2]] (glossary, slide-by-slide script, tough-questions bank).

## The five phases

Phase numbers are the report's. Phase II (Adam) was written on top of Phase III (Route 3): `adam_qite.py` imports
`adaptive_qite.py`, so read III before II if you want to follow the code.

| Phase | Question | Folder | Entry note | Script | Outcome (report Table I) |
|---|---|---|---|---|---|
| I | Does the Saxena–Shalabi–Kyaw (QDP) noise analysis reproduce? | `Replications/practical-limitations-qdp/` | [[Adaptative-QITE/Notes_on_Project/Replication_Practical_Limitations_QDP\|Replication_Practical_Limitations_QDP]] | `replicate.py` | T1, T4 reproduced; T2, T5 partially; T3 not. Factor-of-two correction to Theorem 2 |
| II | Does Adam-style memory improve the step law? | `Adam_QITE/`, `Audits/Audit_Adam_QITE/` | [[Adaptative-QITE/Adam_QITE/Report_Adam_QITE\|Report_Adam_QITE]] | `adam_qite.py`, `verify_audit.py` | Ratchet wins at matched budget; E2 withdrawn; E3 re-attributed to the cap; zero-cost error estimator validated; 147 values pinned, 0 drift |
| III | Route 3: physics-indicator controller (Targets A–D) | `Attack_Route_03/`, `Audits/Audit_Project_Wide/` | [[Adaptative-QITE/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation\|Adaptive_QITE_Attack_Route_3_Implementation]] | `adaptive_qite.py` | 2.2× closer to the floor than fixed-step Euler at a matched 732 circuits; debias factor corrected; noise-floor claim retracted |
| IV | Do the figures of Assi–Vogl–Kumari–LeBlanc (BT) reproduce? | `Replications/Beyond_Trot/` | [[Adaptative-QITE/Notes_on_Project/Replication_Beyond_Trotterization\|Replication_Beyond_Trotterization]] | `replicate.py` (`--check`) | Figures reproduce; four textual claims refuted; boundary-layer target T6 pinned |
| V | P1: what is the Route-3 stiffness? P2: does `tol` control the error? | `Route_03_and_Beyond_Trot/` | [[Adaptative-QITE/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam\|Plan_Beyond_Trot_Into_Route_03_and_Adam]] | `Proposal_1/p1_singularity.py`, `Proposal_2/p2_tol_proportionality.py` | Chart singularity, no gauge leakage (51 pinned); no tolerance fit possible (93 pinned) |

Per-folder READMEs (results and run instructions): `Replications/practical-limitations-qdp/README.md`,
`Replications/Beyond_Trot/README.md`, `Attack_Route_03/README.md`, `Adam_QITE/README.md`.
Paper notes: [[Adaptative-QITE/Notes_on_Project/Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors|Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors]],
[[Adaptative-QITE/Notes_on_Project/Beyond_Trotterization_Variational_Product_Formulas_for_Quantum_Simulation|Beyond_Trotterization_Variational_Product_Formulas_for_Quantum_Simulation]],
[[Adaptative-QITE/Notes_on_Project/Understanding_Quantum_Imaginary_Time_Evolution_and_its_Variational_Form|Understanding_Quantum_Imaginary_Time_Evolution_and_its_Variational_Form]],
[[Adaptative-QITE/Notes_on_Project/Trading_Imaginary_Time_for_Randomness_in_Ground_State_Preparation|Trading_Imaginary_Time_for_Randomness_in_Ground_State_Preparation]]. PDFs are in `Research_Papers/`.

## Folder map

```
Adaptative-QITE/
  Start_Here_Adaptative-QITE.md        this note
  Presentations/                     report (.tex/.pdf), P1/P2 deck, presenter notes, Q&A bank
  Notes_on_Project/                  technical resume, paper notes, replication and implementation notes
  Research_Papers/                   source PDFs (incl. the internal memo)
  Replications/                      practical-limitations-qdp/, Beyond_Trot/   (replicate.py, README, figures/)
  Attack_Route_03/                   Route 3 controller           (adaptive_qite.py, README, Output.txt, figures/)
  Adam_QITE/                         Adam-like controller         (adam_qite.py, verify_audit.py, README, Report_Adam_QITE.md, Output_*.txt, figures/)
  Audits/                            Audit_Project_Wide/ (hub, M1–M6, R_Future_Routes), Audit_Adam_QITE/ (hub, G1–G7)
  Route_03_and_Beyond_Trot/          Plan, Proposal_1/ (P1), Proposal_2/ (P2)  (script, Output.txt, findings, figures/)
  Environment/quantum-environment.yml   conda spec of the `qml-julia` environment
```

**Keep the relative layout.** Scripts resolve paths from `__file__` and import each other across
sibling folders (`adam_qite.py` imports `Attack_Route_03/adaptive_qite.py`; P1 and P2 import both),
and the report's `\graphicspath` reads figures from the `figures/` folders in place.

## How to run

Python 3.11 with PennyLane 0.45.1, NumPy and SciPy, in the `qml-julia` environment
(`Environment/quantum-environment.yml`). Call the interpreter by absolute path; shell activation
does not work non-interactively:

```bash
~/micromamba/envs/qml-julia/bin/python Adaptative-QITE/Adam_QITE/adam_qite.py --target all
~/micromamba/envs/qml-julia/bin/python Adaptative-QITE/Adam_QITE/verify_audit.py --target all
~/micromamba/envs/qml-julia/bin/python Adaptative-QITE/Route_03_and_Beyond_Trot/Proposal_2/p2_tol_proportionality.py --check
```

Every script takes `--target` and `--no-plots`; the replication, P1 and P2 scripts also take
`--check`, which recomputes every pinned value and fails on any change (there is no tolerance
argument). The `Output*.txt` files next to each script are the verbatim stdout of its shipped run;
they are the provenance for the numbers quoted in the notes. Rebuild the report with
`latexmk -pdf Adaptative_QITE_Status_Report.tex` from `Presentations/`.

## Rules the project follows

- **One function per target.** Every figure, table entry or quoted number has its own function, and
  every quoted value is pinned in a verification script.
- **Verdicts use three words:** *reproduced*, *partially reproduced*, *not reproduced*. A documented
  mismatch is a result; no parameter is tuned until a plot matches.
- **Predictions are written before the run** (from Phase V on) and scored as held or failed.
- Controllers are compared in a common currency of **circuits**: one evaluation of $f$ with $N$
  parameters costs $2N+1$ circuits and each retried trial costs one more.
- Constants of the baseline step law (Route 3) were chosen in this work, not taken from the memo.

## Where the project stands

The report's claims ledger (Table XVII) is the authoritative list of what is established,
retracted, withdrawn or unsupported. In short:

- **Established:** Theorem 1 and noisy-flow scaling of the QDP paper; Theorem 2 overstates the
  hardware velocity error by exactly 2×; the 9-parameter metric is structurally rank 6 of 9; Route 3
  beats fixed-step Euler at 732 circuits; the Route-3 stiffness is a chart singularity, not leakage
  or noise; the PI law with explicit `tol` does control the local-error estimate.
- **Not reproduced / unresolved on the paper's side:** the QDP condition numbers (T3, not reproduced) and the
  Fig. 3(b) curves (T5, partially reproduced: one reading only). No error on the replication side is found by the 60-row convention scan.
- **Not supported / withdrawn / retracted:** Adam-style memory as an improvement (not supported),
  Adam's noise robustness E2 (withdrawn), a noise floor near $p\approx10^{-3}$ and the per-step
  debiasing "surprise" (retracted), end-point error scaling as $\mathrm{tol}^\alpha$ (not fittable).

## Open items — where to continue

From the plan ([[Adaptative-QITE/Route_03_and_Beyond_Trot/Plan_Beyond_Trot_Into_Route_03_and_Adam|Plan_Beyond_Trot_Into_Route_03_and_Adam]], §6–§7) and the deck's closing slide:

- **P1b** — state-space step control: detect $\dot\theta$ growing along the soft direction while
  $\sqrt{\dot\theta^\top M\dot\theta}$ stays flat, and bound $\Delta\tau\,\|d\rho/d\tau\|$ instead.
- **P2b** — the pre-registered fixed version of P2 (EPUS, tolerance-aware first step, pre-cap
  metric) on a **second seed** and Target C.
- **P3** — exact testbeds (BT two-level qubit, imaginary-time qubit) to score the memo's estimators
  against the true local error; then **P3c + P4** (seeded noise, optimal start time), **P5** (BT-product
  ansatz QITE on the QIM, $N=4$–$10$), **P6**, **P7**.
- **Precondition for any Adam arm in P2–P5:** fix G2 first (one EMA weighting, $\beta_1=\beta_2$).
- **Open decisions** (plan §7): whether a seeded Gaussian perturbation of $(M,Y)$ is an acceptable
  stand-in for stochastic noise in P3c; whether P5 starts on QIM or LiH (the plan assumes QIM first);
  whether the two derived identities of plan §2 get a standalone note.
- **Recorded, not all fixed:** findings G1–G7; see the remediation status in
  [[Adaptative-QITE/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]. Also R7 in [[Adaptative-QITE/Audits/Audit_Project_Wide/R_Future_Routes|R_Future_Routes]]: ask Gaurav Saxena (QDP paper
  author, on the team) about the Appendix C condition numbers.
- **Known limits** (report, limitations list): every threshold of the baseline is a chosen
  constant; G2 uses two seeds and E3 is a single run; a few checks (Target B reparametrization
  test, some Phase IV side checks) were run once and are not pinned.
