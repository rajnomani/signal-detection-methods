# Signal Detection Methods: PRR, ROR, IC and EBGM on FAERS

A comparison of four disproportionality methods for pharmacovigilance signal detection, run on 14 years of FDA Adverse Event Reporting System (FAERS) data. The worked example is gabapentin and pregabalin.

## Summary

- **Data:** 56 FAERS quarterly ASCII extracts, 2012 Q4 to 2026 Q3 (20.96 million reports).
- **After cleaning:** 17.99 million unique cases (latest version of each case, FDA-deleted cases removed).
- **Analysis set:** 17,730,954 cases (primary suspect drugs only, cases with 30 or fewer events), 3,391 drugs, 1,191,181 drug-event pairs.
- **Methods:** PRR, ROR, IC (BCPNN-style), EBGM.
- **Main finding:** the methods agree on strong signals and differ in how many weak ones they flag. EBGM is the most conservative: it flags 24.4% of pairs, against 37.8% for IC, 42.8% for PRR and 46.1% for ROR.

## Results

Pairs flagged as a signal by each method:

| Method | Signal rule | Pairs flagged | Share of pairs |
|---|---|---|---|
| PRR | PRR ≥ 2, n ≥ 3, chi-square ≥ 4 | 509,367 | 42.8% |
| ROR | lower 95% CI > 1, n ≥ 3 | 549,272 | 46.1% |
| IC | IC025 > 0 | 450,527 | 37.8% |
| EBGM | EB05 ≥ 2 | 290,933 | 24.4% |
| All four | all rules met | 290,845 | 24.4% |

Almost every EBGM signal is also flagged by the other three methods.

![Signals by method](outputs/figures/signals_by_method.png)

Selected terms for gabapentin and pregabalin:

| Drug | Event | Cases | PRR | ROR | IC025 | EBGM | EB05 |
|---|---|---|---|---|---|---|---|
| Gabapentin | Respiratory