# practical-limitations-qdp

Replication of the numerical claims in Saxena, Shalabi & Kyaw, *Practical limitations of
quantum data propagation on noisy quantum processors*, Phys. Rev. Applied **21**, 054014
(2024), [arXiv:2306.13046v2](https://arxiv.org/abs/2306.13046v2).

Source PDF: `../../Research_Papers/1_Practical_limitations_quantum_data_propagation_noisy_quantum_processors.pdf`
Paper note: `Practical_Limitations_of_Quantum_Data_Propagation_on_Noisy_Quantum_Processors.md`
Results note: `Replication_Practical_Limitations_QDP.md` — read that one for the verdicts.

## Run

```bash
~/micromamba/envs/qml-julia/bin/python replicate.py --target all
```

`--target {t1,t2,t3,t4,t5,all}` selects a single target; `--no-plots` skips the PNGs.
Deterministic: `SEED = 0`, and only the optional random-angle sanity sweeps use the RNG at all.
Runtime is under a minute. Dependencies are numpy, scipy, matplotlib and pennylane 0.45.1
(`default.qubit` and `default.mixed`) — all already in `qml-julia`, nothing was installed.

## Targets

| | Claim | Where the paper reports it | Verdict |
|---|---|---|---|
| T1 | Thm 2, $\epsilon/\|\dot\theta\| = \frac{1-(1-p)^N}{(1-p)^N}$ | Eq. (19) p.4, Fig. 1(b) p.2, Fig. 5 p.16 | reproduced |
| T2 | max allowed $p$ at $\delta=0.04$ | Fig. 3 caption p.5, Eq. (18) | partially reproduced (+1.4 % to +3.5 %) |
| T3 | cond(M), $\|M\|$, $\|Y\|$ for both ansätze | Appendix C p.16 | **not reproduced** |
| T4 | $M^{err}=(1-p)^{2N}M$, $Y^{err}=(1-p)^{N}Y$, $(1-p)^N\dot\theta^{err}=\dot\theta$ | Eqs. (B12), (B13), (B16) pp.14–15 | reproduced (~1e-10) |
| T5 | Thm 1 bound, Eq. (17) | Fig. 3(b) p.5, Fig. 6 p.17 | partially reproduced (one reading only) |

## Figures

| File | Reproduces |
|---|---|
| `figures/fig1b_theorem2.png` | Fig. 1(b) p.2 / Fig. 5 p.16 |
| `figures/fig3a_allowed_p.png` | Fig. 3(a) p.5 |
| `figures/fig3b_theorem1.png` | Fig. 3(b) p.5 |
| `figures/fig6b_theorem1_9param.png` | Fig. 6 p.17 |
| `figures/t4_scaling.png` | no paper counterpart — Eq. (19) recovered from the noisy simulation |

## Layout of `replicate.py`

Constants first, every one traceable to a page of the paper (`G` from Eq. 20, `THETA_5` and
`THETA_9` from Appendix C, the Appendix C targets, the Fig. 3(b) read-offs and their
tolerances). Then a shared layer — `h2_hamiltonian`, `gate_list`, `make_rho_fn`,
`make_noisy_rho_fn`, `global_depolarizing_kraus`, `drho`, `metric_M`, `vector_Y`,
`frobenius_cond` — then one function per target, then `print_summary`.

`check_eq_b9()` asserts the construction of $M$ against the paper's own closed form
$M_{kk} = \tfrac12(1-\mathrm{Tr}[P_k\rho_{k-1}P_k\rho_{k-1}])$ (Eq. B9) to 1e-10 for all four
single-qubit gates. That check is what makes the T3 disagreement a finding rather than a
suspected convention error on our side.

`ASSUMPTIONS` at the top of the file is the authoritative list of choices the paper left open;
`--target all` prints it. Nothing in this script was tuned to make a plot match.
