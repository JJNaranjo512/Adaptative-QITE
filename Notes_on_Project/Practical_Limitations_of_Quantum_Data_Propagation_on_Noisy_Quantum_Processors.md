---
fileClass: Paper_Reading
title: Practical limitations of quantum data propagation on noisy quantum processors
authors: Gaurav Saxena, Ahmed Shalabi, Thi Ha Kyaw
year: "2024"
venue: journal
arxiv_id: https://arxiv.org/abs/2306.13046v2
doi: https://doi.org/10.1103/PhysRevApplied.21.054014
status: reading
tags:
  - paper
  - physics
  - quantum-computation
  - paper-notes
---

# Practical limitations of quantum data propagation on noisy quantum processors

> [!info] Citation
> Saxena, G., Shalabi, A. & Kyaw, T., "*Practical limitations of quantum data propagation on noisy quantum processors*", Phys. Rev. Appl., 2024. [arXiv:ID]([https://arxiv.org/abs/](https://arxiv.org/abs/2306.13046v2)) | [DOI]([https://doi.org/](https://doi.org/10.1103/PhysRevApplied.21.054014))

## TL;DR
> They take QITE (hybrid algorithm), and evaluate what's the effect of global Markovian error sources in NISQ computers when computing $M$ and $Y$. Turns out, the larger the size of $\theta$, $N$, produces errors to quickly grow or diverge. 
> To implement QITE, we would need quantum hardware with really low error rates.

## Why I'm reading this
- To come up to speed on the new project with LG Electronics Toronto AI Team: *Thi Ha Kyaw, Jack Baker, and Gaurav Saxena*. For the looks of it, we are going to try to include an adaptative-step algorithm.

## 1. Context & Motivation
- **Problem addressed:** how does the error (Markovian) grow when implementing QITE. 
- **Why it matters:** This type of error, especially depolarizing channels, are present in NISQ computers, such as superconducting ones. 

## 2. Core Idea
> Variational quantum algorithms uses ansatz to find ground states of Hamiltonians. 

- 

## 3. Method / Approach
- **Key assumptions:**
- **Formalism / model used:**
- **Key equations:**
  $$ $$
- **Notation translation** (their symbol → my convention):
  | Their notation | Meaning | My notation |
  |---|---|---|

## 4. Results
- **Main result:**
- **Numerical/experimental setup** (if applicable): 
- **Figures worth remembering:**
  - Fig. X — 
- **Replication:** [[QITE_Hand_Off/Notes_on_Project/Replication_Practical_Limitations_QDP|Replication_Practical_Limitations_QDP]] — Theorem 2 (Eq. 19) and the
  Appendix B scalings reproduce exactly; Eq. (18) is off by 1.4–3.5 %; **Appendix C's
  cond(M) = 66.7239 and 20.22 do not reproduce** (the 9-parameter $M$ is rank 6 of 9 by
  construction, so its cond is infinite); Fig. 3(b) only reproduces if the $\|Y\|$ term of
  Eq. (17) is dropped.

## 5. Critical Evaluation
- **Strengths:**
- **Weaknesses / limitations / assumptions that might not hold:**
- **Open questions the paper leaves unanswered:**
- **What would break this result?**

## 6. Connections
- **Relates to:** [[]]
- **Contradicts / complicates:** [[]]
- **Follow-up papers to read:** [[]]

## 7. My Takeaways
- 

## Glossary / New Terms
- **Term** — definition

---