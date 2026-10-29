import pandas as pd
from pathlib import Path

RAW_DIR = Path("data/raw")
OUT_DIR = Path("data/processed")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Only the columns we need (lowercase)
KEEP = {
    "DEMO": {"primaryid", "caseid", "caseversion", "age", "age_cod",
             "sex", "event_dt", "fda_dt", "rept_cod", "occr_country"},
    "DRUG": {"primaryid", "caseid", "drug_seq", "role_cod",
             "drugname", "prod_ai"},
    "REAC": {"primaryid", "caseid", "pt"},
}


def read_table(path, table):
    df = pd.read_csv(
        path,
        sep="$",
        encoding="latin-1",
        dtype=str,
        usecols=lambda c: c.strip().lower() in KEEP[table],
        low_memory=False,
    )
    df.columns = [c.strip().lower() for c in df.columns]
    return df


def find_files(table):
    """Return {quarter: path} for every file of this table type."""
    files = {}
    for f in RAW_DIR.rglob("*.txt"):
        name = f.stem.upper()
        if name.startswith(table):
            files[name[len(table):].replace("_NEW", "")] = f   # e.g. '23Q1'
    return files


demo_files = find_files("DEMO")
drug_files = find_files("DRUG")
reac_files = find_files("REAC")

for quarter in sorted(demo_files):
    if quarter not in drug_files or quarter not in reac_files:
        print(f"Skipping {quarter}: missing a DRUG or REAC file")
        continue
    print(f"Processing {quarter} ...")
    for table, files in [("DEMO", demo_files), ("DRUG", drug_files), ("REAC", reac_files)]:
        df = read_table(files[quarter], table)
        df["quarter"] = quarter
        df.to_parquet(OUT_DIR / f"{quarter}_{table.lower()}.parquet", index=False)
        print(f"   {table}: {len(df):,} rows")

print("Done.")