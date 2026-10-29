import duckdb
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from pathlib import Path

SIG = Path("data/signals")
OUT = Path("outputs")
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# ---- Figure 1: how many pairs each method flags ----
con = duckdb.connect()
con.execute(f"""
    CREATE TEMP TABLE allm AS
    SELECT sg.*, bb.ebgm, bb.eb05
    FROM read_parquet('{SIG}/signals.parquet') sg
    JOIN read_parquet('{SIG}/ebgm.parquet') bb USING (drug, pt)
""")
r = con.execute("""
    SELECT COUNT(*) AS total,
      SUM(CASE WHEN prr >= 2 AND a >= 3 AND chi2 >= 4 THEN 1 ELSE 0 END) AS prr_n,
      SUM(CASE WHEN ror_lo > 1 AND a >= 3 THEN 1 ELSE 0 END)             AS ror_n,
      SUM(CASE WHEN ic025 > 0 THEN 1 ELSE 0 END)                         AS ic_n,
      SUM(CASE WHEN eb05 >= 2 THEN 1 ELSE 0 END)                         AS ebgm_n,
      SUM(CASE WHEN prr >= 2 AND a >= 3 AND chi2 >= 4 AND ror_lo > 1
                AND ic025 > 0 AND eb05 >= 2 THEN 1 ELSE 0 END)           AS all_n
    FROM allm
""").df().iloc[0]

total = float(r["total"])
labels = ["PRR", "ROR", "IC", "EBGM", "All four"]
vals = [float(r["prr_n"]), float(r["ror_n"]), float(r["ic_n"]),
        float(r["ebgm_n"]), float(r["all_n"])]

fig, ax = plt.subplots(figsize=(7.5, 4.8))
bars = ax.bar(labels, vals)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,.0f}\n({v / total:.1%})",
            ha="center", va="bottom", fontsize=9)
ax.set_ylim(0, max(vals) * 1.2)
ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
ax.set_ylabel("Drug-event pairs flagged")
ax.set_title(f"Signals flagged by each method ({total:,.0f} pairs, FAERS 2012Q4-2026Q3)")
fig.savefig(FIG / "signals_by_method.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# ---- Figure 2: PRR vs EBGM (gabapentin and pregabalin pairs) ----
d = pd.read_csv(OUT / "gabapentin_pregabalin_all_methods.csv").dropna(subset=["prr", "ebgm"])
fig, ax = plt.subplots(figsize=(6.8, 6))
sc = ax.scatter(d["prr"], d["ebgm"], c=np.log10(d["a"]), s=8, alpha=0.6, cmap="viridis")
lo = min(d["prr"].min(), d["ebgm"].min()) * 0.8
hi = max(d["prr"].max(), d["ebgm"].max()) * 1.2
ax.plot([lo, hi], [lo, hi], "k--", lw=1)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel("PRR")
ax.set_ylabel("EBGM")
ax.set_title("EBGM shrinks scores for low-count pairs\n(gabapentin and pregabalin)")
cb = fig.colorbar(sc, ax=ax)
cb.set_label("log10(number of cases)")
fig.savefig(FIG / "prr_vs_ebgm.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# ---- Figure 3: respiratory depression and sedation ----
terms = ["RESPIRATORY DEPRESSION", "SEDATION", "SOMNOLENCE"]
t = d[d["pt"].isin(terms)]
x = np.arange(len(terms))
w = 0.35
fig, ax = plt.subplots(figsize=(7.5, 4.8))
for i, drug in enumerate(["GABAPENTIN", "PREGABALIN"]):
    s = t[t["drug"] == drug].set_index("pt").reindex(terms)
    lower = (s["ebgm"] - s["eb05"]).values
    ax.bar(x + (i - 0.5) * w, s["ebgm"].values, w,
           yerr=[lower, np.zeros(len(terms))], capsize=4, label=drug.title())
ax.axhline(2, ls="--", lw=1, color="gray")
ax.text(len(terms) - 0.5, 2.1, "EB05 = 2 signal threshold", ha="right", fontsize=8, color="gray")
ax.set_xticks(x)
ax.set_xticklabels([s_.title() for s_ in terms])
ax.set_ylabel("EBGM")
ax.set_title("EBGM by drug (whisker = lower bound, EB05)")
ax.legend()
fig.savefig(FIG / "resp_sedation.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print("Saved 3 figures in outputs/figures/")