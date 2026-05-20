#!/usr/bin/env python3

import pandas as pd
import subprocess, os

WORKDIR = "/home/s2451842/chirp_chip_gkmsvmR/lsgkm_vsmc_nfcore/variants"
MODEL   = "/home/s2451842/chirp_chip_gkmsvmR/lsgkm_vsmc_nfcore/models/H3K27me3_VSMC.model.txt"
GKMPREDICT = "/home/s2451842/chirp_chip_gkmsvmR/lsgkm/src/gkmpredict"

REF = f"{WORKDIR}/snp_ref_windows.fa"
ALT = f"{WORKDIR}/snp_alt_windows.fa"

ref_sc = f"{WORKDIR}/ref_scores.txt"
alt_sc = f"{WORKDIR}/alt_scores.txt"

subprocess.run([GKMPREDICT, REF, MODEL, ref_sc, "-T", "4"], check=True)
subprocess.run([GKMPREDICT, ALT, MODEL, alt_sc, "-T", "4"], check=True)

def read_scores(f):
    d = {}
    for l in open(f):
        p = l.strip().split("\t")
        if len(p)>=2:
            d[p[0]] = float(p[1])
    return d

r = read_scores(ref_sc)
a = read_scores(alt_sc)

rows=[]
for k in r:
    if k+"_alt" in a:
        delta = a[k+"_alt"] - r[k]
        rows.append((k, r[k], a[k+"_alt"], delta))

df = pd.DataFrame(rows, columns=["id","ref","alt","delta"])
df["abs_delta"] = df["delta"].abs()
df.sort_values("abs_delta", ascending=False).to_csv(
    f"{WORKDIR}/variant_prioritisation_results.tsv", sep="\t", index=False
)

print(df.head(10))
``
