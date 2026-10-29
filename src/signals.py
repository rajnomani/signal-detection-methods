import duckdb
from pathlib import Path

SIG = Path("data/signals")
OUT = Path("outputs")
OUT.mkdir(exist_ok=True)

con = duckdb.connect()
con.execute("PRAGMA memory_limit='4GB'")
con.execute("PRAGMA temp_directory='data/tmp'")

con.execute(f"""
    CREATE TEMP TABLE base AS
    SELECT drug, pt, a, b, c, d,
           CAST(a AS DOUBLE) AS fa, CAST(b AS DOUBLE) AS fb,
           CAST(c AS DOUBLE) AS fc, CAST(d AS DOUBLE) AS fd,
           CAST(n_total AS DOUBLE) AS fn
    FROM read_parquet('{SIG}/counts.parquet')
""")

con.execute("""
    CREATE TEMP TABLE m AS
    SELECT drug, pt, a, b, c, d,
           (fa / (fa + fb)) / NULLIF(fc / (fc + fd), 0)           AS prr,
           1/fa - 1/(fa + fb) + 1/NULLIF(fc, 0) - 1/(fc + fd)     AS prr_var,
           fn * power(fa * fd - fb * fc, 2)
              / ((fa + fb) * (fc + fd) * (fa + fc) * (fb + fd))   AS chi2,
           (fa * fd) / NULLIF(fb * fc, 0)                         AS ror,
           1/fa + 1/NULLIF(fb, 0) + 1/NULLIF(fc, 0) + 1/fd        AS ror_var,
           log2((fa + 0.5) / (((fa + fb) * (fa + fc) / fn) + 0.5)) AS ic,
           fa
    FROM base
""")

con.execute(f"""
    COPY (
        SELECT drug, pt, a, b, c, d,
               prr,
               exp(ln(prr) - 1.96 * sqrt(prr_var)) AS prr_lo,
               chi2,
               ror,
               exp(ln(ror) - 1.96 * sqrt(ror_var)) AS ror_lo,
               ic,
               ic - 3.3 * power(fa + 0.5, -0.5) - 2 * power(fa + 0.5, -1.5) AS ic025
        FROM m
    ) TO '{SIG}/signals.parquet' (FORMAT PARQUET)
""")

s = f"read_parquet('{SIG}/signals.parquet')"

print("Total drug-event pairs:", con.execute(f"SELECT COUNT(*) FROM {s}").fetchone()[0])
print("\nPairs flagged as a signal by each method:")
print(con.execute(f"""
    SELECT
      SUM(CASE WHEN prr >= 2 AND a >= 3 AND chi2 >= 4 THEN 1 ELSE 0 END) AS prr_signals,
      SUM(CASE WHEN ror_lo > 1 AND a >= 3 THEN 1 ELSE 0 END)             AS ror_signals,
      SUM(CASE WHEN ic025 > 0 THEN 1 ELSE 0 END)                         AS ic_signals,
      SUM(CASE WHEN prr >= 2 AND a >= 3 AND chi2 >= 4
                AND ror_lo > 1 AND ic025 > 0 THEN 1 ELSE 0 END)          AS all_three
    FROM {s}
""").df().to_string(index=False))

print("\nTop 10 gabapentin events by IC025 (a >= 10):")
print(con.execute(f"""
    SELECT pt, a, round(prr, 2) AS prr, round(ror, 2) AS ror, round(ic025, 2) AS ic025
    FROM {s} WHERE drug = 'GABAPENTIN' AND a >= 10
    ORDER BY ic025 DESC LIMIT 10
""").df().to_string(index=False))

print("\nTop 10 pregabalin events by IC025 (a >= 10):")
print(con.execute(f"""
    SELECT pt, a, round(prr, 2) AS prr, round(ror, 2) AS ror, round(ic025, 2) AS ic025
    FROM {s} WHERE drug = 'PREGABALIN' AND a >= 10
    ORDER BY ic025 DESC LIMIT 10
""").df().to_string(index=False))

print("\nRespiratory depression and sedation, both drugs:")
print(con.execute(f"""
    SELECT drug, pt, a, round(prr, 2) AS prr, round(prr_lo, 2) AS prr_lo,
           round(ror, 2) AS ror, round(ror_lo, 2) AS ror_lo, round(ic025, 2) AS ic025
    FROM {s}
    WHERE drug IN ('GABAPENTIN', 'PREGABALIN')
      AND pt IN ('RESPIRATORY DEPRESSION', 'SEDATION', 'SOMNOLENCE')
    ORDER BY pt, drug
""").df().to_string(index=False))

con.execute(f"""
    COPY (
        SELECT * FROM {s} WHERE drug IN ('GABAPENTIN', 'PREGABALIN') ORDER BY drug, ic025 DESC
    ) TO '{OUT}/gabapentin_pregabalin_signals.csv' (HEADER, DELIMITER ',')
""")
print("\nSaved outputs/gabapentin_pregabalin_signals.csv")
print("Done.")