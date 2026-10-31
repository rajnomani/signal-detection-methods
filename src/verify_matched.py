"""Independent recount of Analysis 2 (matched thresholds and top-k overlap).

This does NOT use src/matched.py, src/agreement.py or pandas. It recounts with plain SQL
(DuckDB) and compares with the CSV files, read with Python's csv module.
Run from the project folder:   python src/verify_matched.py
"""
import csv
import statistics
import sys
from itertools import combinations
from pathlib import Path

from verify_agreement import ALLM, README

OUT = Path("outputs")
KS = [100, 1000, 10000]
SCORES = ["prr", "ror", "ic", "ebgm", "eb05"]
# (method, column, operator, limit) written out again, separately from matched.py
FAMILIES = {
    "point estimate >= 2": [("PRR", "prr", ">=", 2), ("ROR", "ror", ">=", 2),
                            ("EBGM", "ebgm", ">=", 2)],
    "lower limit above no association": [("PRR", "prr_lo", ">", 1), ("ROR", "ror_lo", ">", 1),
                                         ("IC", "ic025", ">", 0), ("EBGM", "eb05", ">", 1)],
}
results = []


def check(ok, text):
    results.append(ok)
    print(("PASS" if ok else "FAIL"), text)


def read(name):
    with open(OUT / f"{name}.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def flag_sql(col, op, lim):
    # a missing or not-a-number value never counts as a signal
    return (f"CASE WHEN {col} IS NOT NULL AND NOT isnan({col}) AND {col} {op} {lim} "
            f"THEN 1 ELSE 0 END")


def run_checks(q):
    """q(sql) -> list of rows. The SQL needs a table called allm."""
    n = q("SELECT COUNT(*) FROM allm")[0][0]
    check(int(n) == README["pairs"], f"pair count: I recount {int(n):,}, README {README['pairs']:,}")

    # Part A, counts
    fc = {(r["rule"], r["method"]): r for r in read("matched_flag_counts")}
    check(len(fc) == 7, "matched_flag_counts has 3 + 4 rows")
    mine = {}
    for fam, rules in FAMILIES.items():
        sel = ", ".join(f"SUM({flag_sql(c, o, l)})" for _, c, o, l in rules)
        sums = q(f"SELECT {sel} FROM allm")[0]
        for (m, *_), s in zip(rules, sums):
            mine[(fam, m)] = int(s)
            r = fc[(fam, m)]
            check(int(r["pairs_flagged"]) == int(s), f"matched count, {fam}, {m}")
            check(abs(float(r["share_percent"]) - round(100 * int(s) / int(n), 2)) < 0.006,
                  f"matched percentage, {fam}, {m}")

    # Part A, pairwise overlaps
    pw = {(r["rule"], r["method_a"], r["method_b"]): r for r in read("matched_pairwise")}
    ok = True
    for fam, rules in FAMILIES.items():
        for (ma, ca, oa, la), (mb, cb, ob, lb) in combinations(rules, 2):
            both = int(q(f"SELECT SUM(({flag_sql(ca, oa, la)}) * ({flag_sql(cb, ob, lb)})) "
                         f"FROM allm")[0][0])
            only_a, only_b = mine[(fam, ma)] - both, mine[(fam, mb)] - both
            neither = int(n) - both - only_a - only_b
            r = pw[(fam, ma, mb)]
            ok = ok and (int(r["both"]), int(r["only_a"]), int(r["only_b"]), int(r["neither"])) \
                == (both, only_a, only_b, neither)
    check(ok and len(pw) == 3 + 6, "matched_pairwise: both, only A, only B, neither for all 9 pairs")

    # Part B, rank every pair by every score (own tie-break rule: more cases, drug, event)
    ranks = ", ".join(
        f"ROW_NUMBER() OVER (ORDER BY ({s} IS NULL OR isnan({s})), {s} DESC, a DESC, drug, pt) "
        f"AS rk_{s}" for s in SCORES)
    q(f"CREATE TEMP TABLE ranked AS SELECT drug, pt, a, {ranks} FROM allm")
    ov = {(int(r["k"]), r["score_a"], r["score_b"]): r for r in read("topk_overlap")}
    ok = True
    for k in KS:
        for sa, sb in combinations(SCORES, 2):
            shared = int(q(f"SELECT COUNT(*) FROM ranked WHERE rk_{sa} <= {k} AND rk_{sb} <= {k}")[0][0])
            r = ov[(k, sa, sb)]
            ok = ok and int(r["shared"]) == shared and \
                abs(float(r["share_percent"]) - round(100 * shared / k, 1)) < 0.06
    check(ok and len(ov) == 30, "topk_overlap: shared pairs for all 30 combinations")
    cs = {(int(r["k"]), r["score"]): float(r["median_cases"]) for r in read("topk_cases")}
    ok = True
    for k in KS:
        for s in SCORES:
            vals = [row[0] for row in q(f"SELECT a FROM ranked WHERE rk_{s} <= {k}")]
            ok = ok and abs(statistics.median(vals) - cs[(k, s)]) < 0.06
    check(ok and len(cs) == 15, "topk_cases: median number of cases for all 15 combinations")


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