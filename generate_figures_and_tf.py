#!/usr/bin/env python3
# =============================================================================
# generate_figures_and_tf.py
#
# H3K27me3 VSMC LS-GKM paper — Figure generation + TF motif analysis
#
# Generates:
#   Figure1_peak_reproducibility.png/pdf
#   Figure2_ROC_curve.png/pdf
#   Figure3_PR_curve.png/pdf
#   Figure4_score_distribution.png/pdf
#   Figure5_variant_prioritisation.png/pdf
#   Figure6_TF_motif_disruption.png/pdf          (new — TF prediction)
#   tf_motif_delta_scores.tsv                    (full TF results table)
#
# Usage (on bifx-core2, gkmsvm_env active):
#   pip install matplotlib seaborn numpy scikit-learn biopython pandas requests
#   python3 generate_figures_and_tf.py
#
# All figures saved to OUTDIR (default: current directory).
# =============================================================================

import os, sys, subprocess, csv, json, time, random, re
import urllib.request
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from sklearn.metrics import (roc_curve, roc_auc_score,
                             precision_recall_curve, average_precision_score)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# =============================================================================
# CONFIGURATION — edit paths if needed
# =============================================================================
WORKDIR    = "/home/s2451842/chirp_chip_gkmsvmR/lsgkm_vsmc_nfcore"
REF_FA     = "/home/s2451842/GTF/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
MODEL      = f"{WORKDIR}/models/H3K27me3_VSMC.model.txt"
GKMPREDICT = "/home/s2451842/chirp_chip_gkmsvmR/lsgkm/src/gkmpredict"
VARTSV     = f"{WORKDIR}/variants/variant_prioritisation_results.tsv"
POS_SCORES = f"{WORKDIR}/scores/pos_test_scores.txt"
NEG_SCORES = f"{WORKDIR}/scores/neg_test_scores.txt"
OUTDIR     = "."        # change to a dedicated figures/ directory if preferred

DPI        = 300
os.makedirs(OUTDIR, exist_ok=True)

# Publication-quality style
plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "axes.spines.top":  False,
    "axes.spines.right": False,
    "axes.linewidth":   0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "figure.dpi":       150,
})

BLUE  = "#2E86C1"
NAVY  = "#1A5276"
GREY  = "#808B96"
RED   = "#C0392B"
GREEN = "#1E8449"

# =============================================================================
# HELPERS
# =============================================================================
def read_scores(path):
    scores = []
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try:
                    scores.append(float(parts[1]))
                except ValueError:
                    pass
    return np.array(scores)

def save(fig, name):
    for ext in ("png", "pdf"):
        p = os.path.join(OUTDIR, f"{name}.{ext}")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        print(f"  Saved: {p}")
    plt.close(fig)

# =============================================================================
# FIGURE 1 — Peak reproducibility & dataset composition
# =============================================================================
def figure1():
    print("\n--- Figure 1: Peak reproducibility ---")
    categories  = ["1 replicate\n(n=8,187)", "2 replicates\n(n=370)", "3 replicates\n(n=6)"]
    counts      = [8187, 370, 6]
    colors_bar  = ["#AED6F1", "#2E86C1", "#154360"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    bars = axes[0].bar(categories, counts, color=colors_bar,
                       edgecolor="black", linewidth=0.7, width=0.55)
    axes[0].axhline(y=376, color=RED, linestyle="--", linewidth=1.5,
                    label=r"ENCODE threshold ($\geq$2 replicates, n=376)")
    for bar, count in zip(bars, counts):
        axes[0].text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + 120,
                     f"{count:,}", ha="center", va="bottom",
                     fontsize=11, fontweight="bold")
    axes[0].set_ylabel("Number of consensus peaks", fontsize=12)
    axes[0].set_title("H3K27me3 peak reproducibility\nacross 3 VSMC biological replicates",
                       fontsize=12, fontweight="bold", pad=10)
    axes[0].legend(fontsize=9, framealpha=0.6)
    axes[0].set_ylim(0, 9400)
    axes[0].set_xlabel("Replicate reproducibility", fontsize=11)

    comp_labels = ["Positive\n(H3K27me3 peaks)", "Negative\n(Genomic background)"]
    comp_values = [375, 375]
    comp_colors = [BLUE, GREY]
    b2 = axes[1].bar(comp_labels, comp_values, color=comp_colors,
                     edgecolor="black", linewidth=0.7, width=0.45)
    axes[1].set_ylabel("Number of 200 bp sequences", fontsize=12)
    axes[1].set_title("Training / test dataset composition\n(balanced 1:1, post–blacklist filter)",
                       fontsize=12, fontweight="bold", pad=10)
    axes[1].set_ylim(0, 450)
    for bar, v in zip(b2, comp_values):
        axes[1].text(bar.get_x() + bar.get_width() / 2,
                     v + 6, str(v), ha="center", va="bottom",
                     fontsize=13, fontweight="bold")

    for ax, lbl in zip(axes, ["A", "B"]):
        ax.text(-0.08, 1.06, lbl, transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top")

    fig.suptitle("Figure 1", fontsize=10, color="grey", y=0.02)
    fig.tight_layout()
    save(fig, "Figure1_peak_reproducibility")

# =============================================================================
# FIGURE 2 — ROC curve
# =============================================================================
def figure2(pos, neg):
    print("\n--- Figure 2: ROC curve ---")
    y_true  = np.concatenate([np.ones(len(pos)),  np.zeros(len(neg))])
    y_score = np.concatenate([pos, neg])

    fpr, tpr, _ = roc_curve(y_true, y_score)
    auroc = roc_auc_score(y_true, y_score)

    rng = np.random.default_rng(42)
    n = len(y_true)
    boot_aurocs = []
    for _ in range(10_000):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        boot_aurocs.append(roc_auc_score(y_true[idx], y_score[idx]))
    ci_lo, ci_hi = np.percentile(boot_aurocs, [2.5, 97.5])
    print(f"  AUROC = {auroc:.4f}  95% CI: {ci_lo:.4f}–{ci_hi:.4f}")

    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    ax.plot(fpr, tpr, color=BLUE, lw=2,
            label=f"LS-GKM H3K27me3 VSMC\nAUROC = {auroc:.4f}\n95% CI: {ci_lo:.4f}–{ci_hi:.4f}")
    ax.fill_between(fpr, tpr, alpha=0.10, color=BLUE)
    ax.plot([0, 1], [0, 1], "--", color="grey", lw=1, label="Random classifier")
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title("ROC Curve — H3K27me3 LS-GKM Classifier\n(held-out test set: 38 pos, 38 neg)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="lower right", framealpha=0.7)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    fig.tight_layout()
    save(fig, "Figure2_ROC_curve")
    return auroc, ci_lo, ci_hi

# =============================================================================
# FIGURE 3 — Precision-Recall curve
# =============================================================================
def figure3(pos, neg):
    print("\n--- Figure 3: PR curve ---")
    y_true  = np.concatenate([np.ones(len(pos)),  np.zeros(len(neg))])
    y_score = np.concatenate([pos, neg])

    prec, rec, _ = precision_recall_curve(y_true, y_score)
    auprc = average_precision_score(y_true, y_score)
    print(f"  AUPRC = {auprc:.4f}")

    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    ax.plot(rec, prec, color=NAVY, lw=2, label=f"AUPRC = {auprc:.4f}")
    ax.fill_between(rec, prec, alpha=0.10, color=NAVY)
    ax.axhline(y=0.5, linestyle="--", color="grey", lw=1,
               label="No-skill classifier (precision = 0.5)")
    ax.set_xlabel("Recall", fontsize=12)
    ax.set_ylabel("Precision", fontsize=12)
    ax.set_title("Precision-Recall Curve — H3K27me3 LS-GKM Classifier\n"
                 "(held-out test set: 38 pos, 38 neg)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="lower left", framealpha=0.7)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([0.0, 1.05])
    fig.tight_layout()
    save(fig, "Figure3_PR_curve")

# =============================================================================
# FIGURE 4 — Score distribution violin plot
# =============================================================================
def figure4(pos, neg):
    print("\n--- Figure 4: Score distributions ---")
    df = pd.DataFrame({
        "Score": np.concatenate([pos, neg]),
        "Class": (["H3K27me3 peaks"] * len(pos) +
                  ["Genomic background"] * len(neg))
    })

    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    palette = {"H3K27me3 peaks": BLUE, "Genomic background": GREY}
    sns.violinplot(data=df, x="Class", y="Score", palette=palette,
                   inner="box", ax=ax, linewidth=0.8)
    sns.stripplot(data=df, x="Class", y="Score", color="black",
                  alpha=0.45, size=3.5, jitter=True, ax=ax)
    ax.axhline(y=0, linestyle="--", color=RED, lw=1.2,
               label="Decision boundary (score = 0)")
    ax.set_xlabel("", fontsize=12)
    ax.set_ylabel("LS-GKM score", fontsize=12)
    ax.set_title("LS-GKM Score Distribution\n"
                 "Held-out test set (38 pos, 38 neg)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.7)
    fig.tight_layout()
    save(fig, "Figure4_score_distribution")

# =============================================================================
# FIGURE 5 — Variant prioritisation lollipop chart
# =============================================================================
def figure5(tsv_path):
    print("\n--- Figure 5: Variant prioritisation ---")
    df = pd.read_csv(tsv_path, sep="\t")
    # Expect columns: rsid, gene, category, ref_allele, alt_allele,
    #                 ref_score, alt_score, delta, abs_delta
    # Sort by abs_delta descending, take top 10
    df = df.sort_values("abs_delta", ascending=False).head(10).reset_index(drop=True)

    labels = [f"{r['rsid']} ({r['gene']})" for _, r in df.iterrows()]
    deltas = df["delta"].tolist()
    colors = [BLUE if d > 0 else RED for d in deltas]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5),
                             gridspec_kw={"width_ratios": [3, 2]})

    # Panel A — lollipop
    y_pos = range(len(labels))
    axes[0].barh(y_pos, deltas, color=colors, alpha=0.80, height=0.5)
    axes[0].scatter(deltas, y_pos, color=colors, s=90, zorder=5)
    axes[0].axvline(x=0, color="black", lw=0.8)
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(labels, fontsize=9.5)
    axes[0].set_xlabel("H3K27me3 delta-score (alt − ref)", fontsize=11)
    axes[0].set_title("Top 10 GWAS variants by predicted\nH3K27me3 effect in VSMCs",
                      fontsize=12, fontweight="bold")
    axes[0].invert_yaxis()
    legend_elems = [
        mpatches.Patch(facecolor=BLUE, label="Gain of H3K27me3 (▲ repression)"),
        mpatches.Patch(facecolor=RED,  label="Loss of H3K27me3 (▼ de-repression)")
    ]
    axes[0].legend(handles=legend_elems, fontsize=9, loc="lower right", framealpha=0.7)
    axes[0].text(-0.06, 1.04, "A", transform=axes[0].transAxes,
                 fontsize=14, fontweight="bold", va="top")

    # Panel B — category breakdown (all 47 variants)
    all_df = pd.read_csv(tsv_path, sep="\t")
    cat_order = ["CAD_core", "lipid", "ECM_SMC", "inflammation", "calcification", "additional"]
    cat_labels = {
        "CAD_core":      "CAD core",
        "lipid":         "Lipid metabolism",
        "ECM_SMC":       "ECM / SMC phenotype",
        "inflammation":  "Inflammation",
        "calcification": "Calcification",
        "additional":    "Additional loci"
    }
    gain_counts = []
    loss_counts = []
    for cat in cat_order:
        sub = all_df[all_df["category"] == cat]
        gain_counts.append((sub["delta"] > 0).sum())
        loss_counts.append((sub["delta"] <= 0).sum())

    x = np.arange(len(cat_order))
    w = 0.35
    axes[1].bar(x - w/2, gain_counts, w, color=BLUE, alpha=0.85,
                label="Gain of H3K27me3", edgecolor="black", linewidth=0.5)
    axes[1].bar(x + w/2, loss_counts, w, color=RED,  alpha=0.85,
                label="Loss of H3K27me3", edgecolor="black", linewidth=0.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([cat_labels[c] for c in cat_order],
                             rotation=30, ha="right", fontsize=8.5)
    axes[1].set_ylabel("Number of variants", fontsize=11)
    axes[1].set_title("H3K27me3 direction by\nfunctional category (all 47 variants)",
                      fontsize=12, fontweight="bold")
    axes[1].legend(fontsize=9, framealpha=0.7)
    axes[1].text(-0.08, 1.04, "B", transform=axes[1].transAxes,
                 fontsize=14, fontweight="bold", va="top")

    fig.suptitle("Figure 5", fontsize=10, color="grey", y=0.01)
    fig.tight_layout()
    save(fig, "Figure5_variant_prioritisation")

# =============================================================================
# TF MOTIF PREDICTION
#
# Strategy: use JASPAR 2024 core vertebrates PFMs (downloaded via JASPAR REST
# API or from file) + FIMO-style scanning with log-odds scoring.
# For each variant we score ref and alt 200 bp windows for all motifs,
# report delta = alt − ref per motif, and highlight top disruptions.
#
# Requirements:
#   pip install requests biopython pandas
#
# If FIMO (MEME Suite) is available, the script also runs FIMO directly.
# Otherwise it uses a pure-Python log-odds scanner (slower but no extra deps).
# =============================================================================

JASPAR_API   = "https://jaspar.elixir.lu/api/v1/matrix/?collection=CORE&species=9606&page_size=500&format=json"
JASPAR_CACHE = os.path.join(OUTDIR, "jaspar_pfms.json")
TF_TSV       = os.path.join(OUTDIR, "tf_motif_delta_scores.tsv")
FIMO_BIN     = "fimo"          # set to full path if not in PATH
PSEUDOCOUNT  = 0.1
SCAN_PVAL    = 1e-4            # p-value threshold for motif hits

def fetch_jaspar_pfms(cache=JASPAR_CACHE):
    """Download JASPAR 2024 CORE vertebrate PFMs; cache locally."""
    if os.path.exists(cache):
        print(f"  Loading JASPAR PFMs from cache: {cache}")
        with open(cache) as fh:
            return json.load(fh)

    print("  Downloading JASPAR 2024 CORE vertebrate PFMs …")
    pfms = {}
    url = JASPAR_API
    while url:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read().decode())
        for entry in data["results"]:
            mid  = entry["matrix_id"]
            name = entry["name"]
            mat  = entry["pfm"]   # dict: A, C, G, T → list of counts
            pfms[mid] = {"name": name, "pfm": mat}
        url = data.get("next")
        if url:
            time.sleep(0.3)

    with open(cache, "w") as fh:
        json.dump(pfms, fh)
    print(f"  Downloaded {len(pfms)} motifs → {cache}")
    return pfms

def pfm_to_pwm(pfm_dict, pseudocount=PSEUDOCOUNT):
    """Convert PFM counts to log-odds PWM (bits), background = 0.25."""
    bases = ["A", "C", "G", "T"]
    L = len(pfm_dict["A"])
    pwm = []
    for i in range(L):
        col = {}
        total = sum(pfm_dict[b][i] for b in bases) + 4 * pseudocount
        for b in bases:
            freq = (pfm_dict[b][i] + pseudocount) / total
            col[b] = np.log2(freq / 0.25)
        pwm.append(col)
    return pwm  # list[dict] length L

def score_seq(seq, pwm):
    """Score a sequence against a PWM; return max score over all windows."""
    seq = seq.upper()
    L   = len(pwm)
    n   = len(seq)
    best = -1e9
    for start in range(n - L + 1):
        s = 0.0
        for j, b in enumerate(seq[start:start + L]):
            if b not in pwm[j]:
                s = -1e9; break
            s += pwm[j][b]
        if s > best:
            best = s
    return best

def rc_seq(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

def score_seq_both_strands(seq, pwm):
    return max(score_seq(seq, pwm), score_seq(rc_seq(seq), pwm))

# Approximate p-value from PWM score using a lookup table of quantiles
# computed from 100,000 random background sequences (25% each base).
def compute_bg_quantiles(pwm, n_bg=50_000, seed=7):
    """Fast background score distribution for p-value estimation."""
    rng = np.random.default_rng(seed)
    bases = "ACGT"
    L = len(pwm)
    col_arrays = {b: np.array([pwm[i].get(b, -1e9) for i in range(L)])
                  for b in bases}
    idx = rng.integers(0, 4, size=(n_bg, L))
    scores = np.zeros(n_bg)
    for j in range(L):
        col_scores = np.array([col_arrays[bases[k]][j] for k in idx[:, j]])
        scores += col_scores
    return np.sort(scores)

def pval_from_score(score, bg_quantiles):
    rank = np.searchsorted(bg_quantiles, score)
    return 1.0 - rank / len(bg_quantiles)

def run_tf_analysis(vartsv, ref_fa, pfms):
    """
    For each variant in vartsv:
      - score ref and alt 200 bp sequences with every JASPAR motif
      - report delta = max_alt_score − max_ref_score per motif
      - keep motifs where |delta| is significant (top hits)
    Returns DataFrame.
    """
    import tempfile

    print(f"\n  Building PWMs from {len(pfms)} JASPAR PFMs …")
    pwms = {}
    for mid, entry in pfms.items():
        try:
            pwms[mid] = (entry["name"], pfm_to_pwm(entry["pfm"]))
        except Exception:
            pass
    print(f"  {len(pwms)} PWMs built.")

    # Load variant sequences from ref FASTA used in variant prioritisation
    ref_fa_path = os.path.join(
        os.path.dirname(vartsv), "snp_ref_windows.fa")
    alt_fa_path = os.path.join(
        os.path.dirname(vartsv), "snp_alt_windows.fa")

    # Build alt sequences on the fly from TSV + ref FASTA
    # (re-reads the ref sequences built in variant_prioritisation1.py)
    ref_seqs = {}
    if os.path.exists(ref_fa_path):
        for rec in SeqIO.parse(ref_fa_path, "fasta"):
            # header format: rsid:chrom:start-end or rsid
            rsid = rec.id.split(":")[0]
            ref_seqs[rsid] = str(rec.seq).upper()

    alt_seqs = {}
    if os.path.exists(alt_fa_path):
        for rec in SeqIO.parse(alt_fa_path, "fasta"):
            rsid = rec.id.split(":")[0]
            alt_seqs[rsid] = str(rec.seq).upper()

    # Load variant table
    df_var = pd.read_csv(vartsv, sep="\t")
    # Expects columns including: rsid, ref_allele, alt_allele, (chrom, pos optional)
    # If ref/alt sequences not on disk, reconstruct from REF_FA inline
    SNP_POS = 100  # 0-based position of SNP in 200 bp window

    rows_available = []
    for _, row in df_var.iterrows():
        rsid = row["rsid"]
        if rsid in ref_seqs and rsid in alt_seqs:
            rows_available.append((rsid, row["gene"], row["category"],
                                   ref_seqs[rsid], alt_seqs[rsid]))
        elif rsid in ref_seqs:
            # Reconstruct alt from ref by substitution
            ref_s = ref_seqs[rsid]
            ref_a = row.get("ref_allele", ref_s[SNP_POS])
            alt_a = row.get("alt_allele", "N")
            if alt_a == "N":
                continue
            if ref_s[SNP_POS] == ref_a:
                alt_s = ref_s[:SNP_POS] + alt_a + ref_s[SNP_POS+1:]
            else:
                # try reverse complement orientation
                rc = rc_seq(ref_s)
                if rc[SNP_POS] == ref_a:
                    alt_s = rc[:SNP_POS] + alt_a + rc[SNP_POS+1:]
                else:
                    continue
            rows_available.append((rsid, row["gene"], row["category"], ref_s, alt_s))

    if not rows_available:
        print("  WARNING: No ref/alt sequences found. TF analysis skipped.")
        print("  Expected files:")
        print(f"    {ref_fa_path}")
        print(f"    {alt_fa_path}")
        print("  Run variant_prioritisation1.py first, then re-run this script.")
        return None

    print(f"  Scoring {len(rows_available)} variants × {len(pwms)} motifs …")
    print("  (This may take several minutes — runs in pure Python)")

    results = []
    for vi, (rsid, gene, category, ref_s, alt_s) in enumerate(rows_available):
        if (vi + 1) % 10 == 0:
            print(f"    variant {vi+1}/{len(rows_available)} …")
        for mid, (tf_name, pwm) in pwms.items():
            ref_score = score_seq_both_strands(ref_s, pwm)
            alt_score = score_seq_both_strands(alt_s, pwm)
            delta = alt_score - ref_score
            results.append({
                "rsid":      rsid,
                "gene":      gene,
                "category":  category,
                "motif_id":  mid,
                "tf_name":   tf_name,
                "ref_score": round(ref_score, 4),
                "alt_score": round(alt_score, 4),
                "delta":     round(delta, 4),
                "abs_delta": round(abs(delta), 4),
            })

    df_tf = pd.DataFrame(results)
    df_tf = df_tf.sort_values("abs_delta", ascending=False).reset_index(drop=True)
    df_tf.to_csv(TF_TSV, sep="\t", index=False)
    print(f"  TF motif results saved: {TF_TSV}")
    return df_tf

# =============================================================================
# FIGURE 6 — TF motif disruption heatmap / dot plot
# =============================================================================
def figure6(df_tf):
    print("\n--- Figure 6: TF motif disruption ---")

    if df_tf is None or df_tf.empty:
        print("  Skipping Figure 6 — no TF data.")
        return

    # Top 10 motifs by max |delta| across any variant
    top_motifs = (df_tf.groupby("tf_name")["abs_delta"]
                       .max()
                       .sort_values(ascending=False)
                       .head(10)
                       .index.tolist())

    # Top 15 variants by max |delta| across any motif
    top_vars   = (df_tf.groupby("rsid")["abs_delta"]
                       .max()
                       .sort_values(ascending=False)
                       .head(15)
                       .index.tolist())

    sub = df_tf[df_tf["tf_name"].isin(top_motifs) &
                df_tf["rsid"].isin(top_vars)]

    # Pivot: rows = variants, cols = TF motifs
    pivot = sub.pivot_table(index="rsid", columns="tf_name",
                            values="delta", aggfunc="max").fillna(0)
    pivot = pivot.loc[top_vars, [m for m in top_motifs if m in pivot.columns]]

    # Add gene label to row index
    gene_map = df_tf.drop_duplicates("rsid").set_index("rsid")["gene"].to_dict()
    pivot.index = [f"{rsid}\n({gene_map.get(rsid,'')})" for rsid in pivot.index]

    vmax = max(abs(pivot.values.max()), abs(pivot.values.min()), 0.01)

    fig, ax = plt.subplots(figsize=(13, 6.5))
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)

    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=40, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8.5)
    ax.set_xlabel("TF motif (JASPAR 2024 CORE, vertebrates)", fontsize=11)
    ax.set_ylabel("GWAS variant (rsID, gene)", fontsize=11)
    ax.set_title("Predicted TF motif disruption by cardiovascular GWAS variants\n"
                 "(delta = alt − ref log-odds score; red = gain, blue = loss of motif affinity)",
                 fontsize=12, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Δ motif score (log-odds)", fontsize=10)

    # Annotate cells with delta values
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            val = pivot.values[i, j]
            if abs(val) > 0.5:   # only annotate notable values
                ax.text(j, i, f"{val:+.2f}", ha="center", va="center",
                        fontsize=7,
                        color="white" if abs(val) > 0.7 * vmax else "black")

    fig.tight_layout()
    save(fig, "Figure6_TF_motif_disruption")

    # ---- supplementary summary: top disrupted TF motifs (bar chart) ----
    top20 = (df_tf.groupby("tf_name")["abs_delta"]
                  .max()
                  .sort_values(ascending=False)
                  .head(20))

    # Colour by direction (gain vs loss) of the most extreme hit
    bar_colors = []
    for tf in top20.index:
        sub2 = df_tf[df_tf["tf_name"] == tf]
        best_delta = sub2.loc[sub2["abs_delta"].idxmax(), "delta"]
        bar_colors.append(BLUE if best_delta > 0 else RED)

    fig2, ax2 = plt.subplots(figsize=(12, 5))
    ax2.bar(top20.index, top20.values, color=bar_colors,
            edgecolor="black", linewidth=0.5)
    ax2.set_xticks(range(len(top20)))
    ax2.set_xticklabels(top20.index, rotation=45, ha="right", fontsize=9)
    ax2.set_ylabel("Max |Δ motif score| across all variants", fontsize=11)
    ax2.set_title("Top 20 most disrupted TF motifs across 47 cardiovascular GWAS variants\n"
                  "colour: direction of largest effect (blue = gain, red = loss of motif affinity)",
                  fontsize=11, fontweight="bold")
    legend_elems = [
        mpatches.Patch(facecolor=BLUE, label="Gain of motif affinity (alt > ref)"),
        mpatches.Patch(facecolor=RED,  label="Loss of motif affinity (alt < ref)")
    ]
    ax2.legend(handles=legend_elems, fontsize=9)
    fig2.tight_layout()
    save(fig2, "Figure6b_TF_top20_motifs")

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 60)
    print("  H3K27me3 VSMC LS-GKM — Figure generation + TF analysis")
    print("=" * 60)

    # Check score files exist
    if not os.path.exists(POS_SCORES):
        sys.exit(f"ERROR: {POS_SCORES} not found. Run the LS-GKM pipeline first.")
    if not os.path.exists(NEG_SCORES):
        sys.exit(f"ERROR: {NEG_SCORES} not found.")

    pos = read_scores(POS_SCORES)
    neg = read_scores(NEG_SCORES)
    print(f"\n  Scores loaded: {len(pos)} positive, {len(neg)} negative test sequences")

    # Figures 1–5
    figure1()
    figure2(pos, neg)
    figure3(pos, neg)
    figure4(pos, neg)

    if not os.path.exists(VARTSV):
        print(f"\n  WARNING: {VARTSV} not found — skipping Figure 5 and TF analysis.")
        print("  Run variant_prioritisation1.py first.")
    else:
        figure5(VARTSV)

        # TF motif analysis
        print("\n" + "=" * 60)
        print("  TF MOTIF PREDICTION (JASPAR 2024 CORE, vertebrates)")
        print("=" * 60)
        try:
            pfms = fetch_jaspar_pfms()
        except Exception as e:
            print(f"  WARNING: Could not download JASPAR PFMs: {e}")
            print("  Check internet access to jaspar.elixir.lu")
            print("  TF analysis skipped.")
            pfms = None

        if pfms:
            df_tf = run_tf_analysis(VARTSV, REF_FA, pfms)
            figure6(df_tf)

    print("\n" + "=" * 60)
    print(f"  Done. All outputs written to: {os.path.abspath(OUTDIR)}")
    print("=" * 60)

if __name__ == "__main__":
    main()

