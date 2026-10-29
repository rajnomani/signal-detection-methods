import duckdb
from pathlib import Path

PROC = Path("data/processed")
OUT = Path("data/combined")
OUT.mkdir(parents=True, exist_ok=True)

con = duckdb.connect()
con.execute("PRAGMA memory_limit='4GB'")

for table in ["demo", "drug", "reac"]:
    out = OUT / f"all_{table}.parquet"
    con.execute(f"""
        COPY (
            SELECT * FROM read_parquet('{PROC}/*_{table}.parquet', union_by_name=true)
        ) TO '{out}' (FORMAT PARQUET)
    """)
    n = con.execute(f"SELECT COUNT(*) FROM read_parquet('{out}')").fetchone()[0]
    print(f"Wrote {out}  ({n:,} rows)")

print("Done.")