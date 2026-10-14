---
title: "M1 — Depolarizing debias correction"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]]"
---

# M1 — Depolarizing debias correction

> [!success] Implemented and confirmed exactly — see [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE#6. Implementation status|hub §6]]
> Corrected `measure_f` drives the debiased flow error to the finite-difference noise floor
> (<0.0001% at every $p$ tested) from the predicted 5.3%/34.7% biased — matching this note's
> table to the digit. Targets B and D re-run; no stalls post-fix.

**Findings:** F1 (critical)  **Effort:** S (one expression + a re-run)  **Blocks:** Targets B and
D are not reportable until this lands; [[QITE_Hand_Off/Audits/Audit_Project_Wide/M2_Noise_Model_Falsifiability|M2_Noise_Model_Falsifiability]] and
[[QITE_Hand_Off/Audits/Audit_Project_Wide/M3_Trajectory_Comparison_Metric|M3_Trajectory_Comparison_Metric]] depend on it.

## What is wrong

`measure_f` (`Attack_Route_03/adaptive_qite.py:284-286`) debiases with

```python
if p > 0.0 and debias:
    M = M / (1.0 - p) ** (2 * N)
    Y = Y / (1.0 - p) ** N
```

following Eq. (7) of the memo, $\tilde f = (1-p)^{-N} f$. That is the right factor *for the noise
placement the QDP paper's Appendix B assumes* — noise on the derivative line only, with the $\rho$
inside $\{H,\rho\}$ left ideal.

It is not the noise placement this script uses. `make_noisy_rho_fn` puts a global depolarizing
channel after every gate and returns $\rho^{err}$, and that same $\rho^{err}$ is then fed to *both*
`drho` and `vector_Y`. Write $q := (1-p)^N$ and
$\rho^{err} = q\,\rho + (1-q)\,\mathbb{1}/d$, so $\partial_k\rho^{err} = q\,\partial_k\rho$ and
$\mathrm{Tr}[\partial_k\rho] = 0$. Then

$$M^{err} = q^2 M \qquad\text{(the } M \text{ rescale in the code is exact)}$$

$$Y^{err}_k = -\mathrm{Tr}\big[\partial_k\rho^{err}\{H,\rho^{err}\}\big]
 = -\Big[q^2 + \tfrac{q(1-q)}{2}\Big] g_k
 = \Big[q^2 + \tfrac{q(1-q)}{2}\Big] Y_k$$

using $Y_k = -\mathrm{Tr}[\partial_k\rho\,H] \equiv -g_k$ for a normalised pure state — the identity
[[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]] already verified numerically to $2.3\times10^{-10}$.

So $Y^{err} \neq q\,Y$, and dividing by $q$ leaves

$$\dot\theta_{\text{debiased}} = \Big[q + \tfrac{1-q}{2}\Big]\,\dot\theta_{\text{clean}} \neq \dot\theta_{\text{clean}}.$$

**This is the project's own result, not a new one.** [[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]]
derived exactly this under "New result — Theorem 2 is a factor of two pessimistic on hardware", and
[[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §4 records it. It was simply never propagated into
`Attack_Route_03`.

## Why it matters

The mis-scaling is what Target B measures and reports as a physical surprise. Predicted relative
magnitude error of the flow, from the algebra above ($N=5$):

| $p$ | debiased $\lvert c-1\rvert$ | biased $\lvert c-1\rvert$ | reduction from debiasing |
|---|---|---|---|
| 0.01 | 2.450 % | 2.577 % | **4.901 %** |
| 0.02 | 4.804 % | 5.315 % | 9.608 % |
| 0.10 | 20.475 % | 34.675 % | 40.951 % |
| 0.20 | 33.616 % | 102.588 % | **67.232 %** |

The note reports "~5% at $p=0.01$ growing to ~67% at $p=0.2$". Both endpoints are reproduced by a
model in which *nothing physical is happening* — the flow is simply being scaled by the wrong
constant. In fact the reduction has a closed form: the debiased and biased magnitude errors are
$(1-q)/2$ and $(1-q)/(2q)$, so

$$\text{reduction} \;=\; 1 - \frac{(1-q)/2}{(1-q)/(2q)} \;=\; 1 - q \;=\; 1-(1-p)^N$$

exactly, independent of the Hamiltonian, the ansatz and $\theta$. The "effect" Target B reports as
a property of debiasing is a two-line identity in $p$ and $N$.

## Change to make

In `measure_f`, use the bias factors that match this script's own noise model:

```python
if p > 0.0 and debias:
    q = (1.0 - p) ** N
    M = M / q**2
    Y = Y / (q**2 + q * (1.0 - q) / 2.0)
```

With this, $\dot\theta$ is **exact** for any $p$ in this noise model — not approximately corrected.
Keep the old expression available behind a flag (`debias="eq7"` vs `debias="hardware"`) so the two
can be compared directly; the difference between them *is* the factor-of-two finding, and showing
it is worth a figure.

Update `ASSUMPTIONS[2]` (`adaptive_qite.py:97-101`), which currently asserts the Eq. (7) form is
"algebraically identical to the memo's $\tilde f = (1-p)^{-N} f$ when $M$ is full rank". That
equivalence is true, but both sides are the wrong factor for this noise placement, and the
assumption should say so.

## How to verify

1. **Exactness.** At a fixed $\theta$ and $p \in \{10^{-3}, 10^{-2}, 0.1, 0.2\}$, assert
   $\|\dot\theta_{\text{debiased}} - \dot\theta_{\text{clean}}\| / \|\dot\theta_{\text{clean}}\| < 10^{-8}$.
   The old expression fails this at $p=0.2$ by ~34 %.
2. **Re-run Targets B and D.** Prediction: the corrected debiased run collapses onto the noiseless
   baseline up to the noisy *energy readout* alone — because debiasing fixes the ODE's right-hand
   side but the energy is still measured on $\rho^{err}$, which no rescaling touches. That residual
   gap is the honest content of Target B, and it should now be the *only* thing left in it.
3. If the corrected debiased run still tracks worse than the biased one, F1 is not the explanation
   and the module says so in the note — do not tune to make it come out the expected way.

## Documents to update

- [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] — the Target B paragraph, Assumption 3.
- `Attack_Route_03/README.md` — the Target B row.
- [[QITE_Hand_Off/Notes_on_Project/Technical_Resume_Project_Limiting_QITE|Technical_Resume_Project_Limiting_QITE]] §6, the "Adaptive step-size control (implemented)"
  bullet, if the re-run changes the picture.
