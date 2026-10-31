"""Independent recount of Analysis 3 (ties check and top-k lists for pairs with 10+ cases).

This does NOT use src/robust.py, src/agreement.py or pandas. It recounts with plain SQL
(DuckDB) and compares with the CSV files, read with Python's csv module.
Run from the project folder:   python src/verify_robust.py
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
results = []


def check(ok, text):
    results.append(ok)
    print(("PASS" if ok else "FAIL"), text)


def read(name):
    with open(OUT / f"{name}.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def run_checks(q):
    """q(sql) -> list of rows. The SQL needs a table called allm."""
    n = q("SELECT COUNT(*) FROM allm")[0][0]
    check(int(n) == README["pairs"], f"pair count: I recount {int(n):,}, README {README['pairs']:,}")

    # Part A: missing, infinite, value at rank k, pairs above, pairs equal
    rows = {(r["score"], int(r["k"])): r for r in read("ties_check")}
    ok = True
    for s in SCORES:
        good = f"{s} IS NOT NULL AND NOT isnan({s})"
        missing = int(q(f"SELECT COUNT(*) FROM allm WHERE NOT ({good})")[0][0])
        infinite = int(q(f"SELECT COUNT(*) FROM allm WHERE ({good}) AND isinf({s})")[0][0])
        for k in KS:
            r = rows[(s, k)]
            v = q(f"SELECT {s} FROM allm WHERE {good} ORDER BY {s} DESC LIMIT 1 OFFSET {k - 1}")[0][0]
            if v == float("inf"):
                above = 0
                equal = int(q(f"SELECT COUNT(*) FROM allm WHERE ({good}) AND isinf({s}) AND {s} > 0")[0][0])
            else:
                above = int(q(f"SELECT COUNT(*) FROM allm WHERE ({good}) AND {s} > {v!r}")[0][0])
                equal = int(q(f"SELECT COUNT(*) FROM allm WHERE ({good}) AND {s} = {v!r}")[0][0])
            strad = "yes" if above + equal > k else "no"
            ok = ok and int(r["missing"]) == missing and int(r["infinite"]) == infinite \
                and (float(r["value_at_k"]) == v or abs(float(r["value_at_k"]) - v) < 6e-5) \
                and int(r["pairs_above"]) == above and int(r["pairs_equal"]) == equal \
                and r["tie_straddles_cutoff"] == strad
    check(ok and len(rows) == 15, "ties_check: all 15 rows (missing, infinite, value at k, above, equal)")

    # Part B: only pairs with 10 or more cases
    m = int(q("SELECT COUNT(*) FROM allm WHERE a >= 10")[0][0])
    bands = {r["cases_in_pair"]: int(r["pairs"]) for r in read("agreement_by_count_band")}
    check(m == bands["10-29"] + bands["30-99"] + bands["100+"],
          f"pairs with 10+ cases: I recount {m:,}, agreement file bands add to "
          f"{bands['10-29'] + bands['30-99'] + bands['100+']:,}")
    ranks = ", ".join(
        f"ROW_NUMBER() OVER (ORDER BY ({s} IS NULL OR isnan({s})), {s} DESC, a DESC, drug, pt) "
        f"AS rk_{s}" for s in SCORES)
    q(f"CREATE TEMP TABLE ranked10 AS SELECT drug, pt, a, {ranks} FROM allm WHERE a >= 10")
    ov = {(int(r["k"]), r["score_a"], r["score_b"]): r for r in read("topk_min10_overlap")}
    ok = True
    for k in KS:
        for sa, sb in combinations(SCORES, 2):
            shared = int(q(f"SELECT COUNT(*) FROM ranked10 WHERE rk_{sa} <= {k} AND rk_{sb} <= {k}")[0][0])
            r = ov[(k, sa, sb)]
            ok = ok and int(r["shared"]) == shared and \
                abs(float(r["share_percent"]) - round(100 * shared / k, 1)) < 0.06
    check(ok and len(ov) == 30, "topk_min10_overlap: shared pairs for all 30 combinations")
    cs = {(int(r["k"]), r["score"]): float(r["median_cases"]) for r in read("topk_min10_cases")}
    ok = True
    for k in KS:
        for s in SCORES:
            vals = [row[0] for row in q(f"SELECT a FROM ranked10 WHERE rk_{s} <= {k}")]
            ok = ok and abs(statistics.median(vals) - cs[(k, s)]) < 0.06
    check(ok and len(cs) == 15, "topk_min10_cases: median number of cases for all 15 combinations")


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