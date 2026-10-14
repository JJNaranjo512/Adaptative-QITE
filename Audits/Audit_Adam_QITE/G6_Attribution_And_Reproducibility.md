---
title: "G6 — E3's win is not attributable, and two reported numbers have no code path"
tags:
  - audit
  - physics
  - quantum-computation
  - adaptative-qite
audit: "[[QITE_Hand_Off/Audits/Audit_Adam_QITE/Audit_Adam_QITE|Audit_Adam_QITE]]"
---

# G6 — Attribution and reproducibility

**Severity: medium.** Closes G6. Recorded, not fixed.

Three separable problems, all about whether a reported number means what its caption says.

## 6.1 E3's "one real win" is not attributable

README:67 calls E3 "**The one real win, and it is attributable**", and the reasoning given is
sound in form: the Route-3 cap-off control arm exists precisely so a difference can be assigned to
the step law rather than to the cap. That is the right instinct and it is better than E1 or E2
manage. It still does not carry the claim, for four reasons.

**The four arms are two arms.** From `Output_Pre_Audit.txt`:

| `rcond` | arm | final $E$ | steps |
|---|---|---|---|
| $10^{-8}$ | Route 3 cap on | $-1.144712$ | 108 |
| $10^{-8}$ | **Adam cap on** | $-1.144712$ | **108** |
| $10^{-2}$ | Route 3 cap on / off | $-1.143456$ | 17 |
| $10^{-2}$ | Adam cap on / off | $-1.142285$ | 14 |

At $\mathrm{rcond}=10^{-8}$ the cap binds on every step, so the two cap-on arms are the same
trajectory to six digits and the Adam law contributes nothing. At $\mathrm{rcond}=10^{-2}$,
$\mathrm{cond}(M) \le 125$ gives a cap of $2000/125 = 16 \gg \Delta\tau_{\max}$, so it never binds
and cap-on equals cap-off. README:80-83 identifies both degeneracies and correctly calls them
sanity checks that passed. But it then reports a table of eight rows as though it were eight
measurements. **Four of the eight rows are duplicates**, and the entire claim rests on the single
cap-off row at $\mathrm{rcond}=10^{-8}$.

**The comparison is not cost-matched.** The win is $5.653\times10^{-4}$ in 28 steps against Route 3
cap-off's $1.398\times10^{-3}$ in 22 steps. Adam took **27 % more steps**, and E3 has no circuit
column and no `energy_at_budget` call — the M6 currency E1 uses is simply absent from this target.
A step-law change that produces smaller average steps will buy accuracy with circuits, and this
table cannot tell that apart from the law being better. [[QITE_Hand_Off/Audits/Audit_Project_Wide/M6_Honest_Cost_Accounting|M6_Honest_Cost_Accounting]] exists to
stop exactly this, and E1 applies it; E3 does not.

**There is no matched-mean-step control arm.** The clean control for "is it the law or is it the
step size" is Route 3 cap-off re-run at a $\Delta\tau_{\max}$ that reproduces Adam's mean step. It
is a one-line arm and it is not there.

**$n=1$.** One seed (`SEED_C`), one ansatz, one rank deficiency. README:158-160 already flags this
("E3's cap-off win deserves checking on the 5-parameter ansatz and at other seeds before it is
treated as general"), which is the correct caveat — it is just not compatible with calling the
result attributable in the summary table two sections earlier.

V4 adds one more line on the cap-off arms: they run 28 and 22 steps with **zero** forced-through
steps and 2 retries. Nothing pathological is happening; they are simply short runs on a stiff
problem, which is where [[QITE_Hand_Off/Audits/Audit_Adam_QITE/G5_Error_Estimator_Validated_And_Its_Blind_Spot|G5_Error_Estimator_Validated_And_Its_Blind_Spot]]'s
$\partial_\theta M^{+}$ blind spot predicts $\hat e$ is least trustworthy. That the one favourable
result sits exactly in the regime where the error estimator is least reliable is a reason to want
the seed sweep before the result travels.

## 6.2 The degeneracy check has no code path

README:117-121 reports that the controller is a strict generalization of the ratchet, confirmed by
two runs:

> $\beta_1=\beta_2=0$, $k_I=10$ collapses to exactly Route 3 (172 steps, $E=-1.122993$);
> memoryless with stock $k_I$: 172 steps, $E=-1.122928$, 95 % pinned.

Neither number appears in `adam_qite.py` and neither appears in `Output_Pre_Audit.txt`. There is
no target, no flag and no test that produces them. A reader cannot regenerate them from the
repository.

V5 runs the check externally and **the numbers are exactly right**:

```
  configuration                steps      final E   pinned@GROW
  Route 3 (ratchet)              172    -1.122993           n/a
  b1=b2=0, k_I=10                172    -1.122993          100%
  b1=b2=0, k_I=0.15              172    -1.122928           95%
  b1=0.9 b2=0.999 (stock)        190    -1.122852            5%
```

So this is not a wrong claim — it is a **correct claim with no code path**, which is the same
defect [[QITE_Hand_Off/Audits/Audit_Project_Wide/Audit_Adaptative_QITE|Audit_Adaptative_QITE]] recorded as F3. The fix is to promote V5 into a
`check_ratchet_degeneracy()` alongside the parent's `check_eq_b9_self()`, so the generalization
claim is asserted by the script rather than by a README paragraph. It is worth doing precisely
because the claim is true: it is the cleanest structural statement the module makes, and right now
it is the least defensible one on paper.

## 6.3 `energy_at_budget` over-credits by one step

`adam_qite.py:339-343`:

```python
def energy_at_budget(run, n_params, budget):
    cum = circuit_cumsum(run, n_params)
    i = int(np.searchsorted(cum, budget))
    return float(run["E"][min(i + 1, len(run["E"]) - 1)])
```

`cum[k]` is the cost after $k+1$ steps. With the default `side="left"`, a budget that lands
*strictly inside* step $k+1$ returns `E[i+1]` — the energy the run reached only *after* spending
more than the budget. On an exact boundary it is correct. V8, on a synthetic run of 5 steps at 12
circuits each with $E=-\text{step}$:

```
 budget  steps affordable  correct E  reported E
     12                 1       -1.0        -1.0   ok
     20                 1       -1.0        -2.0   off by one step
     24                 2       -2.0        -2.0   ok
     30                 2       -2.0        -3.0   off by one step
     47                 3       -3.0        -4.0   off by one step
     48                 4       -4.0        -4.0   ok
```

**It fires in E1.** The matched budget is `min` over the three arms, so it is set by whichever arm
is cheapest — the fixed-step one — and lands on a boundary for that arm only (V12):

```
  matched budget = 732 circuits  (min over the three arms)
  Adam     total= 2285 ( 32% of run)  affordable= 60 steps  correct E=-1.120983  reported E=-1.121027  OVER-CREDITED by 1 step
  Route 3  total= 2069 ( 35% of run)  affordable= 60 steps  correct E=-1.120952  reported E=-1.120992  OVER-CREDITED by 1 step
  fixed    total=  732 (100% of run)  affordable= 61 steps  correct E=-1.116691  reported E=-1.116691  ok (lands on a boundary)
```

Two consequences, in increasing order of importance:

1. **The off-by-one does not change E1's verdict.** Both adaptive arms gain one free step, so the
   Adam-minus-Route-3 margin moves from $3.5\times10^{-5}$ to $3.1\times10^{-5}$ and keeps its
   sign. The adaptive-vs-fixed margin *is* inflated, since only the fixed arm is scored honestly —
   but not nearly enough to threaten a gap of $4\times10^{-3}$. Report it as a defect, not as a
   retraction.
2. **The budget is set by the weakest arm.** 732 circuits is 32 % of Adam's run and 35 % of Route
   3's. The headline "at the matched circuit budget Adam is marginally ahead" therefore compares
   the two adaptive controllers **one third of the way through both runs**, where neither has
   converged and the ordering has no reason to be stable — and indeed README:65 notes it reverses
   by the end. This is the larger problem: `min` is the wrong reduction when one arm is a
   deliberately handicapped baseline. Matching the two *adaptive* arms to each other (2069
   circuits, as V7/V10 do) is the comparison the verdict is actually about.

The fix is one line:

```python
i = int(np.searchsorted(cum, budget, side="right"))
return float(run["E"][min(i, len(run["E"]) - 1)])
```

and separately, an E1 budget taken over the arms being compared rather than over all three.

## Proposed work

1. **Add circuit accounting to E3.** Same `circuit_cumsum`/`energy_at_budget` treatment E1 gets;
   without it the cap-off row cannot be read.
2. **Add the matched-mean-step Route-3 cap-off arm**, and drop the four duplicate rows or mark
   them as identities in the table itself rather than in prose below it.
3. **Sweep E3 over seeds and over the 5-parameter ansatz** before the win is cited anywhere.
   README:158-160 already lists this; G6's contribution is that it blocks the word "attributable".
4. **Promote the degeneracy check into the script** as `check_ratchet_degeneracy()`.
5. **Fix `energy_at_budget`** as above, and re-derive E1's budget from the arms under comparison.
6. **Record the correction, not a retraction.** E1's ordering survives; say so when the fix lands,
   so the changed digits are not mistaken for a changed result.

**Falsifiable prediction.** If E3's cap-off win is the step law, it survives cost-matching and
appears at other seeds. If it is the smaller mean step, a Route-3 cap-off arm with
$\Delta\tau_{\max}$ tuned to Adam's mean step will match it within noise. G6 asserts only that the
current table cannot distinguish these; it does not predict which way it goes.
