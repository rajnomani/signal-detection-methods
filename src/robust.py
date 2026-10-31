"""Analysis 3: two checks on the top-k lists from Analysis 2.

Part A (ties): are there missing, infinite or tied scores near the cut-off of the top lists?
    If many pairs share exactly the same score at rank k, the tie-break rule (more cases first,
    then drug name, then event name) decides who gets in.
Part B (at least 10 cases): repeat the top-k lists using only pairs with 10 or more cases,
    so that pairs with very few cases cannot fill the top of the lists.

Reads data/signals/ (read only). Writes 3 small CSV files to outputs/.
Run from the project folder:   python src/robust.py
"""
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from agreement import load, make_flags, readme_check

OUT = Path("outputs")
KS = [100, 1000, 10000]
SCORES = ["prr", "ror", "ic", "ebgm", "eb05"]
MIN_CASES = 10


def part_a(df, out_dir):
    rows = []
    for s in SCORES:
        col = df[s]
        missing = int(col.isna().sum())
        infinite = int(np.isinf(col).sum())
        valid = col.dropna().sort_values(ascending=False).to_numpy()
        for k in KS:
            v = valid[k - 1]
            above = int((col > v).sum())
            equal = int((col == v).sum())
            rows.append({"score": s, "missing": missing, "infinite": infinite, "k": k,
                         "value_at_k": round(float(v), 4), "pairs_above": above,
                         "pairs_equal": equal,
                         "tie_straddles_cutoff": "yes" if above + equal > k else "no"})
    pd.DataFrame(rows).to_csv(out_dir / "ties_check.csv", index=False)


def part_b(df, out_dir):
    sub = df[df["a"] >= MIN_CASES]
    tops, cases = {}, []
    for s in SCORES:
        order = sub.sort_values([s, "a", "drug", "pt"], ascending=[False, False, True, True],
                                na_position="last").index
        for k in KS:
            idx = order[:k]
            tops[(s, k)] = set(idx)
            cases.append({"k": k, "score": s,
                          "median_cases": round(float(sub.loc[idx, "a"].median()), 1)})
    rows = []
    for k in KS:
        for a_, b_ in combinations(SCORES, 2):
            shared = len(tops[(a_, k)] & tops[(b_, k)])
            rows.append({"k": k, "score_a": a_, "score_b": b_, "shared": shared,
                         "share_percent": round(100 * shared / k, 1)})
    pd.DataFrame(rows).to_csv(out_dir / "topk_min10_overlap.csv", index=False)
    pd.DataFrame(cases).to_csv(out_dir / "topk_min10_cases.csv", index=False)
    return len(sub)


def analyse(df, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    part_a(df, out_dir)
    return part_b(df, out_dir)


def main():
    skip = "--skip-readme-check" in sys.argv   # for testing only
    df = load()
    print(f"{len(df):,} drug-event pairs loaded")
    if not readme_check(df, make_flags(df)) and not skip:
        print("\nSTOP: I cannot reproduce the README counts, so nothing was written.")
        sys.exit(1)
    n = analyse(df, OUT)
    print(f"\n{n:,} pairs have at least {MIN_CASES} cases (used for Part B)")
    print("Wrote 3 files to outputs/: ties_check.csv, topk_min10_overlap.csv, topk_min10_cases.csv")
    for f in ["ties_check", "topk_min10_overlap", "topk_min10_cases"]:
        print("\n" + f)
        print(pd.read_csv(OUT / f"{f}.csv").to_string(index=False))


if __name__ == "__main__":
    main()