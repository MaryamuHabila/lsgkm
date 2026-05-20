#!/usr/bin/env python3

import subprocess, pandas as pd, os

FIMO = "/home/s2451842/miniconda3/envs/meme/bin/fimo"
DB   = "/home/s2451842/miniconda3/envs/meme/.../JASPAR2018_CORE_non-redundant.meme"

REF = "snp_ref_clean.fa"
ALT = "snp_alt_clean.fa"

def run_fimo(fa):
    cmd = [FIMO,"--thresh","1e-4","--text",DB,fa]
    res = subprocess.run(cmd, capture_output=True, text=True)

    rows=[]
    for l in res.stdout.splitlines():
        if l.startswith("#"): continue
        p=l.split("\t")
        if len(p)<8: continue
        rows.append({
            "tf":p[1],
            "seq":p[2],
            "start":int(p[3]),
            "end":int(p[4]),
            "score":float(p[6])
        })
    return pd.DataFrame(rows)

ref = run_fimo(REF)
alt = run_fimo(ALT)

ref.to_csv("fimo_ref.tsv", sep="\t")
alt.to_csv("fimo_alt.tsv", sep="\t")

print("FIMO done")

