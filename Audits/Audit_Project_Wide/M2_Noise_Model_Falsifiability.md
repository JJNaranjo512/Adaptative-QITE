---
title: "M2 — Noise model falsifiability"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]]"
---

# M2 — Noise model falsifiability

> [!success] Implemented and confirmed — see [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE#6. Implementation status|hub §6]]
> `target_b_per_step_check()` added; collinearity residual $\sim10^{-8}$, $c_{meas}$ matches
> $c(p)$ to 6 digits. The reparametrization test (verification step 3) also confirmed the
> mechanism: the biased run tracks a $\tau/c(p)$-rescaled clean baseline better than the plain
> one. "Flow-direction error" retired in favor of "flow-magnitude error" throughout.

**Findings:** F2, F3 (both critical)  **Effort:** S  **Depends on:**
[[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|M1_Depolarizing_Debias_Correction]]

## What is wrong

**F3 — a reported number has no code path.** [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]]
(Target B) and `Attack_Route_03/README.md` both report that, checked at a fixed $\theta$,
$\|\dot\theta_{\text{debiased}} - \dot\theta_{\text{clean}}\|$ is smaller than the biased
equivalent, "by ~5% at $p=0.01$ growing to ~67% at $p=0.2$".
`target_b_noisy_debiasing()` (`adaptive_qite.py:480-536`) computes no such quantity. It builds a
noiseless baseline, runs four full trajectories, and prints final energies and interpolated RMS
deviations — nothing else. Running the committed script does not produce those numbers.

The numbers are *correct* ([[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|M1_Depolarizing_Debias_Correction]] reproduces both endpoints
analytically), which makes this a reproducibility defect rather than a factual one. But per the
project's own rules — every reported number traceable to code that produces it — it should not
stand.

**F2 — "flow-direction error" is a misnomer.** Since $M^{err} = q^2 M$ and $Y^{err} \propto Y$,

$$\dot\theta^{err} = c(p)\,\dot\theta_{\text{clean}}, \qquad c(p) = 1 + \frac{1-q}{2q},$$

**exactly collinear** with the clean flow. Global depolarizing noise, in this model, does not rotate
the flow at all — it only rescales it. Both the note and the README describe the quantity as a
direction error.

## Why it matters

The distinction is not cosmetic. A pure scalar rescale of $\dot\theta$ is *partly absorbed* by an
adaptive step-size controller, because $\dot\theta \to c\,\dot\theta$ with $\Delta\tau \to
\Delta\tau/c$ is the same trajectory in $\theta$ — a reparametrization of imaginary time. That is a
far simpler and more likely explanation for "debiasing doesn't visibly help the trajectory" than
the residual-freeze interlock the note currently invokes, and it is testable.

It also bounds what Route 3 can be expected to do about depolarizing noise: nothing the controller
does can distinguish $c\,\dot\theta$ from a genuinely faster flow, so this noise channel is
invisible to *every* step-size indicator that is homogeneous in $\dot\theta$. That is a real
structural statement about the method, and it is currently obscured by the wording.

## Change to make

Add a `target_b_per_step_check()` (or a clearly-named section inside Target B) that, at fixed
$\theta_0$ and over the $p$ grid, computes and prints:

- $\|\dot\theta^{err} - \dot\theta_{\text{clean}}\| / \|\dot\theta_{\text{clean}}\|$, biased and debiased —
  the number the note already reports;
- the **collinearity residual**
  $\left\| \frac{\dot\theta^{err}}{\|\dot\theta^{err}\|} - \frac{\dot\theta_{\text{clean}}}{\|\dot\theta_{\text{clean}}\|} \right\|$,
  which this analysis predicts is $\approx 0$ to machine precision;
- the measured scale $c_{\text{meas}} = \dot\theta^{err}\!\cdot\!\dot\theta_{\text{clean}} / \|\dot\theta_{\text{clean}}\|^2$
  against the predicted $c(p) = 1 + (1-q)/(2q)$.

Then replace "flow-direction error" with "flow-magnitude error" wherever it appears, and state the
collinearity explicitly as a property of this noise channel.

## How to verify

1. Collinearity residual $< 10^{-10}$ across the whole $p$ grid. If it is not, the noise channel is
   doing something the algebra in [[QITE_Hand_Off/Audits/Audit_Project_Wide/M1_Depolarizing_Debias_Correction|M1_Depolarizing_Debias_Correction]] does not capture, and F1's
   fix is incomplete — that would be a genuine finding, worth chasing.
2. $|c_{\text{meas}} - c(p)| / c(p) < 10^{-8}$.
3. **The reparametrization test**, which is the point of this module: run the biased controller at
   $p$ and the noiseless controller with every $\Delta\tau$ divided by $c(p)$. If the two
   trajectories agree in $\theta$, the "biased tracks better" effect is a time-reparametrization
   artifact and Target B's conclusion has to be rewritten around that, not around the residual
   freeze.

## Documents to update

- [[QITE_Hand_Off/Notes_on_Project/Adaptive_QITE_Attack_Route_3_Implementation|Adaptive_QITE_Attack_Route_3_Implementation]] — Target B, and the "genuine surprise" framing.
- `Attack_Route_03/README.md` — the Target B row.
