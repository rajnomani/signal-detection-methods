"""Independent recount of the agreement results.

This does NOT use src/agreement.py or pandas. It recounts everything with plain SQL
(DuckDB) and compares with the CSV files, which it reads with Python's csv module.
Run from the project folder:   python src/verify_agreement.py
"""
import csv
import sys
from pathlib import Path

SIG = Path("data/signals")
OUT = Path("outputs")
NAMES = ["PRR", "ROR", "IC", "EBGM"]
README = {"PRR": 509367, "ROR": 549272, "IC": 450527, "EBGM": 290933,
          "all four": 290845, "pairs": 1191181}

ALLM = f"""
    SELECT sg.*, bb.ebgm, bb.eb05
    FROM read_parquet('{SIG}/signals.parquet') sg
    JOIN read_parquet('{SIG}/ebgm.parquet') bb USING (drug, pt)
"""
FLAGS = """
    WITH f AS (
      SELECT a,
        CASE WHEN prr >= 2 AND a >= 3 AND chi2 >= 4 THEN 1 ELSE 0 END AS f_prr,
        CASE WHEN ror_lo > 1 AND a >= 3 THEN 1 ELSE 0 END AS f_ror,
        CASE WHEN ic025 > 0 THEN 1 ELSE 0 END AS f_ic,
        CASE WHEN eb05 >= 2 THEN 1 ELSE 0 END AS f_ebgm
      FROM allm)
"""
BAND = """CASE WHEN a <= 4 THEN '3-4' WHEN a <= 9 THEN '5-9' WHEN a <= 29 THEN '10-29'
          WHEN a <= 99 THEN '30-99' ELSE '100+' END"""

results = []


def check(ok, text):
    results.append(ok)
    print(("PASS" if ok else "FAIL"), text)


def read(name):
    with open(OUT / f"agreement_{name}.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def run_checks(q):
    """q(sql) -> list of rows. The SQL needs a table called allm."""
    n, *sums, allfour = q(FLAGS + "SELECT COUNT(*), SUM(f_prr), SUM(f_ror), SUM(f_ic), SUM(f_ebgm), "
                                  "SUM(f_prr*f_ror*f_ic*f_ebgm) FROM f")[0]
    my = dict(zip(NAMES, sums))
    my["all four"], my["pairs"] = allfour, n

    # 1. my own counts against the README
    for k, want in README.items():
        check(int(my[k]) == want, f"README count for {k}: I recount {int(my[k]):,}, README {want:,}")

    # 2. flag counts file
    got = {r["method"]: r for r in read("flag_counts")}
    for k in NAMES + ["all four"]:
        check(int(got[k]["pairs_flagged"]) == int(my[k]), f"flag_counts file, {k}")
        pct = round(100 * int(my[k]) / int(n), 2)
        check(abs(float(got[k]["share_percent"]) - pct) < 0.006, f"flag_counts percentage, {k}")

    # 3. number of methods flagging each pair
    rows = dict(q(FLAGS + "SELECT f_prr+f_ror+f_ic+f_ebgm AS k, COUNT(*) FROM f GROUP BY k"))
    file = {int(r["methods_flagging"]): int(r["pairs"]) for r in read("n_methods")}
    check(all(int(rows.get(k, 0)) == file[k] for k in range(5)), "n_methods file (0 to 4 methods)")
    check(sum(file.values()) == int(n), "n_methods adds up to the number of pairs")

    # 4. pairwise overlaps
    pw = {(r["method_a"], r["method_b"]): r for r in read("pairwise")}
    ok = True
    for (a_, b_), r in pw.items():
        both = q(FLAGS + f"SELECT SUM(f_{a_.lower()}*f_{b_.lower()}) FROM f")[0][0]
        only_a = int(my[a_]) - int(both)
        only_b = int(my[b_]) - int(both)
        neither = int(n) - int(both) - only_a - only_b
        ok = ok and (int(r["both"]), int(r["only_a"]), int(r["only_b"]), int(r["neither"])) == \
            (int(both), only_a, only_b, neither)
    check(ok and len(pw) == 6, "pairwise file: both, only A, only B, neither for all 6 pairs")

    # 5. by number of cases in the pair
    bands = {r[0]: r for r in q(FLAGS + f"SELECT {BAND} AS band, COUNT(*), SUM(f_prr), SUM(f_ror), "
                                         f"SUM(f_ic), SUM(f_ebgm) FROM f GROUP BY band")}
    ok = True
    for r in read("by_count_band"):
        b = bands[r["cases_in_pair"]]
        ok = ok and int(r["pairs"]) == int(b[1]) and all(
            int(r[f"{m}_flagged"]) == int(b[2 + i]) for i, m in enumerate(NAMES))
    check(ok and len(bands) == 5, "by_count_band file: pairs and flagged counts in all 5 bands")
    check(sum(int(b[1]) for b in bands.values()) == int(n), "bands add up to the number of pairs")


def main():
    import duckdb
    con = duckdb.connect()
    con.execute(f"CREATE TEMP TABLE allm AS {ALLM}")
    run_checks(lambda sql: con.execute(sql).fetchall())
    bad = results.count(False)
    print("\n" + ("ALL %d CHECKS PASS" % len(results) if not bad else f"{bad} of {len(results)} CHECKS FAILED"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()