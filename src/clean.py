import duckdb
from pathlib import Path

COMB = Path("data/combined")
RAW = Path("data/deleted_lists")
OUT = Path("data/clean")
TMP = Path("data/tmp")
OUT.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)

con = duckdb.connect()
con.execute("PRAGMA memory_limit='4GB'")
con.execute(f"PRAGMA temp_directory='{TMP}'")

demo = f"read_parquet('{COMB}/all_demo.parquet')"
drug = f"read_parquet('{COMB}/all_drug.parquet')"
reac = f"read_parquet('{COMB}/all_reac.parquet')"

# 1. Load the FDA lists of deleted cases
files = [str(f) for f in RAW.rglob("*.txt") if f.name.upper().startswith("DELETE")]
print(f"Found {len(files)} deleted-case lists")
if files:
    con.execute(f"""
        CREATE TEMP TABLE deleted AS
        SELECT DISTINCT TRIM(caseid) AS caseid
        FROM read_csv({files!r}, header=false, all_varchar=true, delim='$',
                      columns={{'caseid': 'VARCHAR'}}, ignore_errors=true)
        WHERE regexp_matches(TRIM(caseid), '^[0-9]+$')
    """)
else:
    con.execute("CREATE TEMP TABLE deleted (caseid VARCHAR)")

n_del = con.execute(f"""
    SELECT COUNT(DISTINCT caseid) FROM {demo}
    WHERE caseid IN (SELECT caseid FROM deleted)
""").fetchone()[0]
print(f"Deleted cases found in the data: {n_del:,}")

# 2. Keep only the latest version of each non-deleted case
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
        WHERE caseid NOT IN (SELECT caseid FROM deleted)
    ) WHERE rn = 1
""")

n_before = con.execute(f"SELECT COUNT(*) FROM {demo}").fetchone()[0]
n_after = con.execute("SELECT COUNT(*) FROM keep").fetchone()[0]
print(f"DEMO: {n_before:,} reports -> {n_after:,} unique cases")

# 3. Clean demographics (one row per primaryid)
con.execute(f"""
    COPY (
        SELECT * EXCLUDE (rn) FROM (
            SELECT d.*, ROW_NUMBER() OVER (
                PARTITION BY d.primaryid ORDER BY d.quarter DESC
            ) AS rn
            FROM {demo} d JOIN keep USING (primaryid)
        ) WHERE rn = 1
    ) TO '{OUT}/demo.parquet' (FORMAT PARQUET)
""")

# 4. Clean drugs: suspect drugs only, standardized names
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

# 5. Clean reactions
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