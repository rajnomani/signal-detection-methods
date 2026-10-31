"""Analysis 2: fairer comparison of PRR, ROR, IC and EBGM.

Part A (matched thresholds): the four methods normally use rules of different
strictness (PRR >= 2 is a point estimate; ROR, IC and EBGM use lower limits).
Here each method is judged by the SAME kind of rule, in two versions:
    point estimate >= 2        PRR, ROR, EBGM        (IC is on a log scale, left out)
    lower limit above "no association"
                               PRR (prr_lo > 1), ROR (ror_lo > 1),
                               IC (ic025 > 0), EBGM (eb05 > 1)
Part B (top-k overlap): take the k highest-scoring pairs for each score
(k = 100, 1,000, 10,000) and count how many pairs two scores share.
Ties are broken the same way every time: more cases first, then drug name, then event name.
Pairs with a missing score are ranked last.

Reads data/signals/ (read only). Writes 4 small CSV files to outputs/.
Run from the project folder:   python src/matched.py
"""
import operator
import sys
from itertools import combinations
from pathlib import Path

import pandas as pd

from agreement import kappa, load, make_flags, readme_check

OUT = Path("outputs")
KS = [100, 1000, 10000]
SCORES = ["prr", "ror", "ic", "ebgm", "eb05"]

OPS = {">=": operator.ge, ">": operator.gt}
FAMILIES = {
    "point estimate >= 2": {
        "PRR": ("prr", ">=", 2), "ROR": ("ror", ">=", 2), "EBGM": ("ebgm", ">=", 2)},
    "lower limit above no association": {
        "PRR": ("prr_lo", ">", 1), "ROR": ("ror_lo", ">", 1),
        "IC": ("ic025", ">", 0), "EBGM": ("eb05", ">", 1)},
}


def family_flags(df, rules):
    return pd.DataFrame({m: OPS[op](df[col], t) for m, (col, op, t) in rules.items()})


def part_a(df, out_dir):
    n = len(df)
    counts, pairs = [], []
    for fam, rules in FAMILIES.items():
        fl = family_flags(df, rules)
        for m in fl.columns:
            c = int(fl[m].sum())
            counts.append({"rule": fam, "method": m, "pairs_flagged": c,
                           "share_percent": round(100 * c / n, 2)})
        for a_, b_ in combinations(fl.columns, 2):
            x, y = fl[a_].to_numpy(), fl[b_].to_numpy()
            both = int((x & y).sum())
            only_a, only_b = int((x & ~y).sum()), int((~x & y).sum())
            pairs.append({"rule": fam, "method_a": a_, "method_b": b_, "both": both,
                          "only_a": only_a, "only_b": only_b,
                          "neither": n - both - only_a - only_b,
                          "jaccard": round(both / (both + only_a + only_b), 4),
                          "kappa": round(kappa(x, y), 4)})
    pd.DataFrame(counts).to_csv(out_dir / "matched_flag_counts.csv", index=False)
    pd.DataFrame(pairs).to_csv(out_dir / "matched_pairwise.csv", index=False)


def part_b(df, out_dir):
    tops, cases = {}, []
    for s in SCORES:
        order = df.sort_values([s, "a", "drug", "pt"], ascending=[False, False, True, True],
                               na_position="last").index
        for k in KS:
            idx = order[:k]
            tops[(s, k)] = set(idx)
            cases.append({"k": k, "score": s,
                          "median_cases": round(float(df.loc[idx, "a"].median()), 1)})
    rows = []
    for k in KS:
        for a_, b_ in combinations(SCORES, 2):
            shared = len(tops[(a_, k)] & tops[(b_, k)])
            rows.append({"k": k, "score_a": a_, "score_b": b_, "shared": shared,
                         "share_percent": round(100 * shared / k, 1)})
    pd.DataFrame(rows).to_csv(out_dir / "topk_overlap.csv", index=False)
    pd.DataFrame(cases).to_csv(out_dir / "topk_cases.csv", index=False)


def analyse(df, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    part_a(df, out_dir)
    part_b(df, out_dir)


def main():
    skip = "--skip-readme-check" in sys.argv   # for testing only
    df = load()
    print(f"{len(df):,} drug-event pairs loaded")
    if not readme_check(df, make_flags(df)) and not skip:
        print("\nSTOP: I cannot reproduce the README counts, so nothing was written.")
        sys.exit(1)
    analyse(df, OUT)
    print("\nWrote 4 files to outputs/: matched_flag_counts.csv, matched_pairwise.csv,")
    print("topk_overlap.csv, topk_cases.csv")
    for f in ["matched_flag_counts", "matched_pairwise", "topk_overlap", "topk_cases"]:
        print("\n" + f)
        print(pd.read_csv(OUT / f"{f}.csv").to_string(index=False))


if __name__ == "__main__":
    main()