"""Run Analysis 2 twice, check the two runs are identical, then verify it.
Run from the project folder:   python src/run_matched.py
"""
import hashlib
import subprocess
import sys
from pathlib import Path

FILES = [f"outputs/{n}.csv" for n in
         ["matched_flag_counts", "matched_pairwise", "topk_overlap", "topk_cases"]]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def run(script):
    r = subprocess.run([sys.executable, script], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        sys.exit(1)
    return r.stdout


print("Run 1 of 2 ...")
run("src/matched.py")
first = {f: sha(f) for f in FILES}
print("Run 2 of 2 ...")
run("src/matched.py")
second = {f: sha(f) for f in FILES}
same = True
print()
for f in FILES:
    ok = first[f] == second[f]
    same = same and ok
    print("identical" if ok else "DIFFERENT", first[f][:12], second[f][:12], f)
print("\nAll files identical across both runs." if same else "\nFILES DIFFER BETWEEN RUNS.")
print("\n--- independent recount ---")
r = subprocess.run([sys.executable, "src/verify_matched.py"])
sys.exit(0 if same and r.returncode == 0 else 1)