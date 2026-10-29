import duckdb
from pathlib import Path

COMB = Path("data/combined")
OUT = Path("data/clean")
OUT.mkdir(parents=True, exist_ok=True)

con = duckdb.connect()
con.execute("PRAGMA memory_limit='4GB'")

demo = f"read_parquet('{COMB}/all_demo.parquet')"
drug = f"read_parquet('{COMB}/all_drug.parquet')"
reac = f"read_parquet('{COMB}/all_reac.parquet')"

# 1. Keep only the latest version of each case
con.execute(f"""
    CREATE TEMP TABLE keep AS
    SELECT primaryid FROM (
        SELECT primaryid,
               ROW_NUMBER() OVER (
                   PARTITION BY caseid
                   ORDER BY TRY_CAST(caseversion AS INTEGER) DESC NULLS LAST,
                            fda_dt DESC NULLS LAST,
                            primaryid DESC
               ) AS rn
        FROM {demo}
    ) WHERE rn = 1
""")

n_before = con.execute(f"SELECT COUNT(*) FROM {demo}").fetchone()[0]
n_after = con.execute("SELECT COUNT(*) FROM keep").fetchone()[0]
print(f"DEMO: {n_before:,} reports -> {n_after:,} unique cases")

# 2. Clean demographics
con.execute(f"""
    COPY (
        SELECT d.* FROM {demo} d JOIN keep USING (primaryid)
    ) TO '{OUT}/demo.parquet' (FORMAT PARQUET)
""")

# 3. Clean drugs: suspect drugs only, standardized names
con.execute(f"""
    COPY (
        SELECT DISTINCT primaryid, drug, role_cod FROM (
            SELECT g.primaryid, g.role_cod,
                   UPPER(TRIM(COALESCE(NULLIF(TRIM(g.prod_ai), ''), g.drugname))) AS drug
            FROM {drug} g JOIN keep USING (primaryid)
            WHERE g.role_cod IN ('PS', 'SS')
        ) WHERE drug IS NOT NULL AND drug <> ''
    ) TO '{OUT}/drug.parquet' (FORMAT PARQUET)
""")

# 4. Clean reactions
con.execute(f"""
    COPY (
        SELECT DISTINCT primaryid, UPPER(TRIM(pt)) AS pt
        FROM {reac} r JOIN keep USING (primaryid)
        WHERE pt IS NOT NULL AND TRIM(pt) <> ''
    ) TO '{OUT}/reac.parquet' (FORMAT PARQUET)
""")

for t in ["demo", "drug", "reac"]:
    n = con.execute(f"SELECT COUNT(*) FROM read_parquet('{OUT}/{t}.parquet')").fetchone()[0]
    print(f"Wrote {OUT}/{t}.parquet  ({n:,} rows)")

print("Done.")