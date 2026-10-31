"""How much do PRR, ROR, IC and EBGM agree on the same drug-event pairs?

Reads data/signals/signals.parquet and data/signals/ebgm.parquet (read only).
Writes five small CSV files into outputs/. Run it from the project folder:

    python src/agreement.py

Signal rules are the same as in src/ebgm.py:
    PRR   prr >= 2 and a >= 3 and chi2 >= 4
    ROR   ror_lo > 1 and a >= 3
    IC    ic025 > 0
    EBGM  eb05 >= 2
"""
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

SIG = Path("data/signals")
OUT = Path("outputs")

# Counts printed in the repository README. The script stops if it cannot reproduce them.
README = {"PRR": 509367, "ROR": 549272, "IC": 450527, "EBGM": 290933,
          "all four": 290845, "pairs": 1191181}

BANDS = ["3-4", "5-9", "10-29", "30-99", "100+"]


def load():
    import duckdb
    con = duckdb.connect()
    return con.execute(f"""
        SELECT sg.*, bb.ebgm, bb.eb05
        FROM read_parquet('{SIG}/signals.parquet') sg
        JOIN read_parquet('{SIG}/ebgm.parquet') bb USING (drug, pt)
    """).df()


def make_flags(df):
    return pd.DataFrame({
        "PRR": (df["prr"] >= 2) & (df["a"] >= 3) & (df["chi2"] >= 4),
        "ROR": (df["ror_lo"] > 1) & (df["a"] >= 3),
        "IC": df["ic025"] > 0,
        "EBGM": df["eb05"] >= 2,
    })


def readme_check(df, flags):
    got = {k: int(flags[k].sum()) for k in flags}
    got["all four"] = int(flags.all(axis=1).sum())
    got["pairs"] = len(df)
    ok = True
    for k, want in README.items():
        good = got[k] == want
        ok = ok and good
        print(("PASS" if good else "FAIL"), f"{k}: I get {got[k]:,}, README says {want:,}")
    return ok


def kappa(x, y):
    n = len(x)
    po = float((x == y).mean())
    pa, pb = float(x.mean()), float(y.mean())
    pe = pa * pb + (1 - pa) * (1 - pb)
    return (po - pe) / (1 - pe)


def analyse(df, flags, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    n = len(df)
    names = list(flags.columns)

    # 1. how many pairs each method flags
    rows = [{"method": m, "pairs_flagged": int(flags[m].sum()),
             "share_percent": round(100 * flags[m].sum() / n, 2)} for m in names]
    allfour = int(flags.all(axis=1).sum())
    rows.append({"method": "all four", "pairs_flagged": allfour,
                 "share_percent": round(100 * allfour / n, 2)})
    pd.DataFrame(rows).to_csv(out_dir / "agreement_flag_counts.csv", index=False)

    # 2. how many methods flag each pair
    k = flags.sum(axis=1)
    nm = k.value_counts().reindex(range(5), fill_value=0).sort_index()
    pd.DataFrame({"methods_flagging": nm.index, "pairs": nm.values,
                  "share_percent": (100 * nm.values / n).round(2)}
                 ).to_csv(out_dir / "agreement_n_methods.csv", index=False)

    # 3. each pair of methods: overlap, Jaccard, Cohen's kappa
    rows = []
    for a_, b_ in combinations(names, 2):
        x, y = flags[a_].to_numpy(), flags[b_].to_numpy()
        both = int((x & y).sum())
        only_a = int((x & ~y).sum())
        only_b = int((~x & y).sum())
        neither = n - both - only_a - only_b
        rows.append({"method_a": a_, "method_b": b_, "both": both, "only_a": only_a,
                     "only_b": only_b, "neither": neither,
                     "jaccard": round(both / (both + only_a + only_b), 4),
                     "kappa": round(kappa(x, y), 4)})
    pd.DataFrame(rows).to_csv(out_dir / "agreement_pairwise.csv", index=False)

    # 4. by number of cases in the pair (where the methods should differ most)
    band = pd.cut(df["a"], bins=[0, 4, 9, 29, 99, np.inf], labels=BANDS)
    ratio = (df["ebgm"] / df["prr"]).replace([np.inf, -np.inf], np.nan)
    rows = []
    for b in BANDS:
        m = (band == b).to_numpy()
        row = {"cases_in_pair": b, "pairs": int(m.sum())}
        for name in names:
            c = int(flags.loc[m, name].sum())
            row[f"{name}_flagged"] = c
            row[f"{name}_percent"] = round(100 * c / m.sum(), 2) if m.sum() else 0
        row["median_ebgm_over_prr"] = round(float(ratio[m].median()), 3) if m.sum() else ""
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_dir / "agreement_by_count_band.csv", index=False)

    # 5. do the scores put pairs in the same order? (Spearman rank correlation)
    cols = ["prr", "ror", "ic", "ebgm", "eb05"]
    sp = df[cols].corr(method="spearman").round(4)
    sp.insert(0, "score", sp.index)
    sp.to_csv(out_dir / "agreement_spearman.csv", index=False)


def main():
    skip = "--skip-readme-check" in sys.argv   # for testing only
    df = load()
    flags = make_flags(df)
    print(f"{len(df):,} drug-event pairs loaded")
    if not readme_check(df, flags) and not skip:
        print("\nSTOP: I cannot reproduce the README counts, so nothing was written.")
        sys.exit(1)
    analyse(df, flags, OUT)
    print("\nWrote 5 files to outputs/: agreement_flag_counts.csv, agreement_n_methods.csv,")
    print("agreement_pairwise.csv, agreement_by_count_band.csv, agreement_spearman.csv")
    for f in ["flag_counts", "n_methods", "pairwise", "by_count_band", "spearman"]:
        print("\n" + f)
        print(pd.read_csv(OUT / f"agreement_{f}.csv").to_string(index=False))


if __name__ == "__main__":
    main()