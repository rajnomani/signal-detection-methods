import duckdb

con = duckdb.connect()
con.execute("PRAGMA memory_limit='4GB'")
cd = "read_parquet('data/signals/case_drug.parquet')"
ce = "read_parquet('data/signals/case_event.parquet')"

con.execute(f"CREATE TEMP TABLE nd AS SELECT primaryid, COUNT(*) AS n_drugs FROM {cd} GROUP BY primaryid")
con.execute(f"CREATE TEMP TABLE ne AS SELECT primaryid, COUNT(*) AS n_events FROM {ce} GROUP BY primaryid")

print("All cases, average per case:")
print(con.execute("""
    SELECT round(AVG(n_drugs), 2) AS avg_drugs, round(AVG(n_events), 2) AS avg_events
    FROM nd JOIN ne USING (primaryid)
""").df().to_string(index=False))

for drug, pt in [("GABAPENTIN", "APPENDICOLITH"),
                 ("PREGABALIN", "C-REACTIVE PROTEIN"),
                 ("GABAPENTIN", "SEDATION")]:
    print(f"\n{drug} + {pt}:")
    print(con.execute(f"""
        SELECT COUNT(*) AS cases,
               round(AVG(nd.n_drugs), 1) AS avg_drugs,
               round(AVG(ne.n_events), 1) AS avg_events,
               MAX(ne.n_events) AS max_events
        FROM {cd} d
        JOIN {ce} e USING (primaryid)
        JOIN nd USING (primaryid)
        JOIN ne USING (primaryid)
        WHERE d.drug = '{drug}' AND e.pt = '{pt}'
    """).df().to_string(index=False))