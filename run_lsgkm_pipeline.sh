#!/bin/bash
set -euo pipefail

echo "======================================================"
echo " LS-GKM H3K27me3 VSMC FINAL PIPELINE"
echo " $(date)"
echo "======================================================"

# ---------------- CONFIG ----------------
HERE="$(pwd)"

BOOLEAN_TXT="${HERE}/H3K27me3.consensus_peaks.boolean.txt"
CONSENSUS_BED="${HERE}/H3K27me3.consensus_peaks.bed"

REF_FA="/home/s2451842/GTF/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_FAI="${REF_FA}.fai"
BLACKLIST="/home/s2451842/GTF/hg38.blacklist.bed"

LSGKM="/home/s2451842/chirp_chip_gkmsvmR/lsgkm/src"
WORKDIR="/home/s2451842/chirp_chip_gkmsvmR/lsgkm_vsmc_nfcore"

mkdir -p ${WORKDIR}/{tmp,peaks,fasta/split,models,scores}

echo "Working directory: ${WORKDIR}"

# ---------------- STEP 0 ----------------
echo "STEP 0: Normalize blacklist"
awk '{gsub(/^chr/,"",$1); print $0}' $BLACKLIST | sort -k1,1 -k2,2n > ${WORKDIR}/tmp/blacklist.bed

# ---------------- STEP 1 ----------------
echo "STEP 1: Select peaks ≥2 replicates"
python3 <<EOF
import csv

kept=0
with open("${BOOLEAN_TXT}") as f, open("${WORKDIR}/tmp/rep2plus.bed","w") as o:
    r=csv.DictReader(f, delimiter="\t")
    for row in r:
        if int(row["num_samples"]) >= 2:
            o.write(f"{row['chr']}\t{row['start']}\t{row['end']}\n")
            kept+=1
print("Peaks:", kept)
EOF

# ---------------- STEP 2 ----------------
echo "STEP 2: Center to 200bp"
python3 <<EOF
half=100
with open("${WORKDIR}/tmp/rep2plus.bed") as f, open("${WORKDIR}/peaks/pos_raw.bed","w") as o:
    for i,line in enumerate(f,1):
        c,s,e=line.split()
        s,e=int(s),int(e)
        mid=(s+e)//2
        start=max(0,mid-half)
        end=start+200
        o.write(f"{c}\t{start}\t{end}\tpeak_{i}\n")
EOF

# ---------------- STEP 3 ----------------
echo "STEP 3: Blacklist filter"
bedtools intersect -v -a ${WORKDIR}/peaks/pos_raw.bed \
    -b ${WORKDIR}/tmp/blacklist.bed > ${WORKDIR}/peaks/pos.bed

# ---------------- STEP 4 ----------------
echo "STEP 4: Positive FASTA"
bedtools getfasta -fi ${REF_FA} -bed ${WORKDIR}/peaks/pos.bed \
    -fo ${WORKDIR}/fasta/pos.fa -name

# ---------------- STEP 5 ----------------
echo "STEP 5: Generate negatives"

awk '{if($1!="MT" && $1!="Y" && $1 !~ /_/) print $1"\t0\t"$2}' ${REF_FAI} \
    | sort -k1,1 -k2,2n > ${WORKDIR}/tmp/genome.bed

bedtools subtract -a ${WORKDIR}/tmp/genome.bed \
    -b ${CONSENSUS_BED} > ${WORKDIR}/tmp/bg.bed

bedtools subtract -a ${WORKDIR}/tmp/bg.bed \
    -b ${WORKDIR}/tmp/blacklist.bed > ${WORKDIR}/tmp/bg_clean.bed

bedtools makewindows -b ${WORKDIR}/tmp/bg_clean.bed -w 200 \
    > ${WORKDIR}/tmp/neg.bed

shuf -n 4000 ${WORKDIR}/tmp/neg.bed > ${WORKDIR}/tmp/neg_sample.bed

bedtools getfasta -fi ${REF_FA} -bed ${WORKDIR}/tmp/neg_sample.bed \
    -fo ${WORKDIR}/fasta/neg.fa

# ---------------- STEP 6 ----------------
echo "STEP 6: Clean FASTA"

python3 <<EOF
from Bio import SeqIO
import random

pos=[r for r in SeqIO.parse("${WORKDIR}/fasta/pos.fa","fasta") if "N" not in str(r.seq)]
neg=[r for r in SeqIO.parse("${WORKDIR}/fasta/neg.fa","fasta") if "N" not in str(r.seq)]

neg=neg[:len(pos)]

SeqIO.write(pos,"${WORKDIR}/fasta/pos_clean.fa","fasta")
SeqIO.write(neg,"${WORKDIR}/fasta/neg_clean.fa","fasta")
EOF

# ---------------- STEP 7 ----------------
echo "STEP 7: Split"
cd ${WORKDIR}/fasta
split -l 337 pos_clean.fa

# ---------------- STEP 8 ----------------
echo "STEP 8: Train model"
${LSGKM}/gkmtrain -l 10 -k 6 -d 3 \
    pos_clean.fa neg_clean.fa \
    ${WORKDIR}/models/H3K27me3_VSMC

# ---------------- STEP 9 ----------------
echo "STEP 9: Predict"
${LSGKM}/gkmpredict pos_clean.fa \
    ${WORKDIR}/models/H3K27me3_VSMC.model.txt \
    ${WORKDIR}/scores/pos.txt

echo "DONE"
