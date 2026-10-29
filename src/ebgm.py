import duckdb
import numpy as np
from pathlib import Path
from scipy.optimize import minimize
from scipy.special import gammaln, gammainc, digamma

SIG = Path("data/signals")
OUT = Path("outputs")

MIN_A = 3              # counts.py only kept pairs with at least this many cases
FIT_SAMPLE = 300_000   # pairs used to fit the prior (None = use all pairs)

con = duckdb.connect()
df = con.execute(f"""
    SELECT drug, pt, a,
           CAST(a + b AS DOUBLE) * (a + c) / n_total AS expected
    FROM read_parquet('{SIG}/counts.parquet')
""").df()

n = df["a"].to_numpy(dtype=float)
E = df["expected"].to_numpy(dtype=float)
print(f"Pairs: {len(n):,}")


def log_nb(n, e, alpha, beta):
    """log P(n) when n ~ Poisson(lambda * e) and lambda ~ Gamma(alpha, rate=beta)."""
    return (gammaln(alpha + n) - gammaln(alpha) - gammaln(n + 1)
            + alpha * np.log(beta / (beta + e)) + n * np.log(e / (beta + e)))


def log_mix(n, e, p):
    a1, b1, a2, b2, w = p
    return np.logaddexp(np.log(w) + log_nb(n, e, a1, b1),
                        np.log1p(-w) + log_nb(n, e, a2, b2))


def unpack(x):
    return np.exp(x[0]), np.exp(x[1]), np.exp(x[2]), np.exp(x[3]), 1 / (1 + np.exp(-x[4]))


def nll(x, n, e):
    """Negative log-likelihood, conditioned on n >= MIN_A (we only kept those pairs)."""
    p = unpack(x)
    ll = log_mix(n, e, p)
    p_low = sum(np.exp(log_mix(np.full_like(n, k), e, p)) for k in range(MIN_A))
    p_low = np.clip(p_low, 0, 1 - 1e-12)
    return -(ll - np.log1p(-p_low)).sum()


# 1. Fit the prior
rng = np.random.default_rng(42)
if FIT_SAMPLE is None or FIT_SAMPLE >= len(n):
    idx = np.arange(len(n))
else:
    idx = rng.choice(len(n), size=FIT_SAMPLE, replace=False)

print("Fitting the prior (this can take a few minutes) ...")
x0 = np.array([np.log(0.2), np.log(0.1), np.log(2.0), np.log(4.0), np.log(0.5)])
bounds = [(-8, 8)] * 4 + [(-6, 6)]
res = minimize(nll, x0, args=(n[idx], E[idx]), method="L-BFGS-B",
               bounds=bounds, options={"maxiter": 300})
a1, b1, a2, b2, w = unpack(res.x)
print(f"Fit converged: {res.success}")
print(f"alpha1={a1:.3f} beta1={b1:.3f} alpha2={a2:.3f} beta2={b2:.3f} w={w:.3f}")

# 2. Posterior for every pair
l1 = np.log(w) + log_nb(n, E, a1, b1)
l2 = np.log1p(-w) + log_nb(n, E, a2, b2)
Q = np.exp(l1 - np.logaddexp(l1, l2))

ebgm = np.exp(Q * (digamma(a1 + n) - np.log(b1 + E))
              + (1 - Q) * (digamma(a2 + n) - np.log(b2 + E)))

print("Calculating EB05 for all pairs ...")
lo = np.full(len(n), 1e-6)
hi = np.full(len(n), 1e4)
for _ in range(45):
    mid = np.sqrt(lo * hi)
    cdf = Q * gammainc(a1 + n, (b1 + E) * mid) + (1 - Q) * gammainc(a2 + n, (b2 + E) * mid)
    too_high = cdf > 0.05
    hi = np.where(too_high, mid, hi)
    lo = np.where(too_high, lo, mid)
eb05 = np.sqrt(lo * hi)

df["ebgm"] = ebgm
df["eb05"] = eb05
df.to_parquet(SIG / "ebgm.parquet", index=False)

# 3. Compare all four methods
con.execute(f"""
    CREATE TEMP TABLE allm AS
    SELECT sg.*, bb.ebgm, bb.eb05
    FROM read_parquet('{SIG}/signals.parquet') sg
    JOIN read_parquet('{SIG}/ebgm.parquet') bb USING (drug, pt)
""")

print("\nPairs flagged as a signal by each method:")
print(con.execute("""
    SELECT
      SUM(CASE WHEN prr >= 2 AND a >= 3 AND chi2 >= 4 THEN 1 ELSE 0 END) AS prr_signals,
      SUM(CASE WHEN ror_lo > 1 AND a >= 3 THEN 1 ELSE 0 END)             AS ror_signals,
      SUM(CASE WHEN ic025 > 0 THEN 1 ELSE 0 END)                         AS ic_signals,
      SUM(CASE WHEN eb05 >= 2 THEN 1 ELSE 0 END)                         AS ebgm_signals,
      SUM(CASE WHEN prr >= 2 AND a >= 3 AND chi2 >= 4 AND ror_lo > 1
                AND ic025 > 0 AND eb05 >= 2 THEN 1 ELSE 0 END)           AS all_four
    FROM allm
""").df().to_string(index=False))

for drug in ["GABAPENTIN", "PREGABALIN"]:
    print(f"\nTop 10 {drug.lower()} events by EB05 (a >= 10):")
    print(con.execute(f"""
        SELECT pt, a, round(ebgm, 2) AS ebgm, round(eb05, 2) AS eb05,
               round(prr, 2) AS prr, round(ic025, 2) AS ic025
        FROM allm WHERE drug = '{drug}' AND a >= 10
        ORDER BY eb05 DESC LIMIT 10
    """).df().to_string(index=False))

print("\nRespiratory depression and sedation, all four methods:")
print(con.execute("""
    SELECT drug, pt, a, round(prr, 2) AS prr, round(ror, 2) AS ror,
           round(ic025, 2) AS ic025, round(ebgm, 2) AS ebgm, round(eb05, 2) AS eb05
    FROM allm
    WHERE drug IN ('GABAPENTIN', 'PREGABALIN')
      AND pt IN ('RESPIRATORY DEPRESSION', 'SEDATION', 'SOMNOLENCE')
    ORDER BY pt, drug
""").df().to_string(index=False))

con.execute(f"""
    COPY (
        SELECT * FROM allm WHERE drug IN ('GABAPENTIN', 'PREGABALIN')
        ORDER BY drug, eb05 DESC
    ) TO '{OUT}/gabapentin_pregabalin_all_methods.csv' (HEADER, DELIMITER ',')
""")
print("\nSaved outputs/gabapentin_pregabalin_all_methods.csv")
print("Done.")