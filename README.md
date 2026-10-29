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
| Gabapentin | Respiratory depression | 243 | 9.10 | 9.14 | 2.92 | 8.77 | 7.88 |
| Pregabalin | Respiratory depression | 256 | 6.69 | 6.71 | 2.48 | 6.46 | 5.82 |
| Gabapentin | Sedation | 478 | 8.00 | 8.07 | 2.81 | 7.79 | 7.22 |
| Pregabalin | Sedation | 337 | 3.90 | 3.92 | 1.76 | 3.83 | 3.50 |
| Gabapentin | Somnolence | 2,565 | 5.50 | 5.72 | 2.37 | 5.42 | 5.24 |
| Pregabalin | Somnolence | 3,724 | 5.60 | 5.84 | 2.40 | 5.49 | 5.34 |

For high-count pairs the four methods give nearly the same score. They diverge for low-count pairs, where EBGM shrinks the score toward 1.

![PRR vs EBGM](outputs/figures/prr_vs_ebgm.png)
![Respiratory depression and sedation](outputs/figures/resp_sedation.png)

Full tables: `outputs/gabapentin_pregabalin_all_methods.csv`.

## Methods

**Data cleaning**
- Kept the latest version of each case (by case ID, then case version and FDA date).
- Removed cases listed in the FDA's quarterly deleted-case files (51,351 cases).
- Standardised drug names (active ingredient where given, otherwise the reported drug name, upper-cased) and reaction terms (MedDRA preferred terms, upper-cased).

**Filters for the main analysis**
- Primary suspect (PS) drugs only.
- Cases with more than 30 reported events excluded.
- Drugs with fewer than 100 cases and pairs with fewer than 3 cases excluded.

**Why these filters.** A first run counted all suspect drugs and had no limit on events per case. It produced implausible top signals (for example appendicolith and C-reactive protein for gabapentin and pregabalin). Cases behind those pairs had on average 54 to 75 drugs and 38 to 78 events, against 1.5 drugs and 2.9 events for the average case. After filtering, the top signals were clinically plausible (for example myoclonus and asterixis for gabapentin, withdrawal and substance use disorder terms for pregabalin). The unfiltered results are kept in `outputs/signals_unfiltered_all_suspects.csv` for comparison. Respiratory depression scores fell after filtering (gabapentin PRR from 11.27 to 9.10) but stayed clear signals.

**Disproportionality methods** (for each drug-event pair, from the 2x2 table of a, b, c, d)
- **PRR** with chi-square.
- **ROR** with a 95% confidence interval.
- **IC** as log2((a + 0.5) / (expected + 0.5)), with an approximate lower bound (IC025).
- **EBGM**: a two-component gamma mixture prior fitted by maximum likelihood on a random sample of 300,000 pairs (conditioned on at least 3 cases per pair), then posterior EBGM and EB05 for every pair.

## Limitations

- FAERS is a spontaneous reporting system. Scores show reporting disproportion, not incidence or risk, and do not show that one drug is more dangerous than another.
- Reporting is affected by media attention, litigation and report type. Duplicate reports remain possible despite case versioning.
- Confounding by indication is present (for example "drug ineffective for unapproved indication" and "diabetic neuropathy" appear as signals).
- Some high-ranking terms rest on small counts (about 10 to 30 cases) and likely reflect small clusters of reports. EBGM shrinks these only moderately here.
- Drug names are not mapped to a standard vocabulary, so spelling variants of the same drug are counted separately.
- The EBGM code is my own implementation, not validated against the FDA's reference software. Treat its numbers as a comparison tool.
- The 30-events-per-case cutoff is a judgement call. Changing it changes the results.

## Reproducing

The data is not in this repository. Download the ASCII files for each quarter from the FDA's FAERS Quarterly Data Files page and unzip them into `data/raw/`. Copy each quarter's `DELETE*.txt` file (from its `Deleted` folder) into `data/deleted_lists/`, then run:

```
pip install -r requirements.txt
python src/build_dataset.py   # slim per-quarter parquet files
python src/combine.py         # one file per table
python src/clean.py           # deduplicate, remove deleted cases
python src/counts.py          # 2x2 counts per drug-event pair
python src/signals.py         # PRR, ROR, IC
python src/ebgm.py            # EBGM and four-method comparison
python src/figures.py         # figures
```

`src/load_data.py` is an early exploration script and `src/diagnose.py` is the check that found the multi-drug report problem.

## Data source

FDA Adverse Event Reporting System (FAERS) Quarterly Data Extract Files, public data from the U.S. Food and Drug Administration. This project is for learning and portfolio purposes and is not medical advice.

## License

MIT