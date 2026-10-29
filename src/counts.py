import duckdb
from pathlib import Path

CLEAN = Path("data/clean")
OUT = Path("data/signals")
TMP = Path("data/tmp")
OUT.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)

MIN_DRUG_CASES = 100   # drugs with fewer cases are dropped
MIN_PAIR_CASES = 3     # drug-event pairs with fewer cases are dropped

con = duckdb.connect()
con.execute("PRAGMA memory_limit='4GB'")
con.execute(f"PRAGMA temp_directory='{TMP}'")

drug_src = f"read_parquet('{CLEAN}/drug.parquet')"
reac_src = f"read_parquet('{CLEAN}/reac.parquet')"

# 1. One row per case-drug, only for cases that also have a reaction
print("Step 1: case-drug table ...")
con.execute(f"""
    COPY (
        SELECT DISTINCT primaryid, drug FROM {drug_src}
        WHERE primaryid IN (SELECT primaryid FROM {reac_src})
    ) TO '{OUT}/case_drug.parquet' (FORMAT PARQUET)
""")

# 2. One row per case-event, only for cases that also have a suspect drug
print("Step 2: case-event table ...")
con.execute(f"""
    COPY (
        SELECT primaryid, pt FROM {reac_src}
        WHERE primaryid IN (SELECT primaryid FROM {drug_src})
    ) TO '{OUT}/case_event.parquet' (FORMAT PARQUET)
""")

cd = f"read_parquet('{OUT}/case_drug.parquet')"
ce = f"read_parquet('{OUT}/case_event.parquet')"

# 3. Total number of cases (N)
N = con.execute(f"SELECT COUNT(DISTINCT primaryid) FROM {cd}").fetchone()[0]
print(f"Total cases N = {N:,}")

# 4. Cases per drug and per event
print("Step 3: drug and event totals ...")
con.execute(f"""
    COPY (
        SELECT drug, COUNT(*) AS n_drug FROM {cd}
        GROUP BY drug HAVING COUNT(*) >= {MIN_DRUG_CASES}
    ) TO '{OUT}/drug_counts.parquet' (FORMAT PARQUET)
""")
con.execute(f"""
    COPY (
        SELECT pt, COUNT(*) AS n_event FROM {ce} GROUP BY pt
    ) TO '{OUT}/event_counts.parquet' (FORMAT PARQUET)
""")

dc = f"read_parquet('{OUT}/drug_counts.parquet')"
ec = f"read_parquet('{OUT}/event_counts.parquet')"
print("Drugs kept:", con.execute(f"SELECT COUNT(*) FROM {dc}").fetchone()[0])

# 5. Cases per drug-event pair (this is the slow step)
print("Step 4: drug-event pairs (this can take several minutes) ...")
con.execute(f"""
    COPY (
        SELECT d.drug, e.pt, COUNT(*) AS a
        FROM {cd} d
        JOIN {dc} k USING (drug)
        JOIN {ce} e USING (primaryid)
        GROUP BY d.drug, e.pt
        HAVING COUNT(*) >= {MIN_PAIR_CASES}
    ) TO '{OUT}/pair_counts.parquet' (FORMAT PARQUET)
""")

# 6. Final a, b, c, d table
print("Step 5: building a, b, c, d table ...")
con.execute(f"""
    COPY (
        SELECT p.drug, p.pt,
               p.a,
               k.n_drug - p.a AS b,
               e.n_event - p.a AS c,
               {N} - k.n_drug - e.n_event + p.a AS d,
               {N} AS n_total
        FROM read_parquet('{OUT}/pair_counts.parquet') p
        JOIN {dc} k USING (drug)
        JOIN {ec} e USING (pt)
    ) TO '{OUT}/counts.parquet' (FORMAT PARQUET)
""")

counts = f"read_parquet('{OUT}/counts.parquet')"
print("Drug-event pairs:", con.execute(f"SELECT COUNT(*) FROM {counts}").fetchone()[0])
print("Rows with a negative cell (should be 0):",
      con.execute(f"SELECT COUNT(*) FROM {counts} WHERE b < 0 OR c < 0 OR d < 0").fetchone()[0])

print("\nCheck: top events for GABAPENTIN and PREGABALIN")
print(con.execute(f"""
    SELECT drug, pt, a, b, c, d FROM {counts}
    WHERE drug IN ('GABAPENTIN', 'PREGABALIN')
    ORDER BY a DESC LIMIT 6
""").df().to_string(index=False))
print("Done.")