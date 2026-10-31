"""Run the agreement analysis twice, check the two runs are identical, then verify it.
Run from the project folder:   python src/run_agreement.py
"""
import hashlib
import subprocess
import sys
from pathlib import Path

FILES = [f"outputs/agreement_{n}.csv" for n in
         ["flag_counts", "n_methods", "pairwise", "by_count_band", "spearman"]]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def run(script):
    r = subprocess.run([sys.executable, script], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        sys.exit(1)
    return r.stdout


print("Run 1 of 2 ...")
run("src/agreement.py")
first = {f: sha(f) for f in FILES}
print("Run 2 of 2 ...")
run("src/agreement.py")
second = {f: sha(f) for f in FILES}
same = True
print()
for f in FILES:
    ok = first[f] == second[f]
    same = same and ok
    print("identical" if ok else "DIFFERENT", first[f][:12], second[f][:12], f)
print("\nAll files identical across both runs." if same else "\nFILES DIFFER BETWEEN RUNS.")
print("\n--- independent recount ---")
r = subprocess.run([sys.executable, "src/verify_agreement.py"])
sys.exit(0 if same and r.returncode == 0 else 1)