#!/usr/bin/env bash
###############################################################################
#  chr15_haplotype_sharing.sh
#
#  Minimal pipeline behind the chromosome-15 haplotype-sharing (introgression)
#  heatmap: which wild haplotype is the closest SV-sharing partner of each
#  cultivated haplotype in each 1-Mb window along chr15.
#
#  Steps
#    1  per-haplotype SV calls -> length filter -> merged presence/absence matrix
#    2  non-overlapping genomic windows from the genome index
#    3  pairwise Jaccard similarity (cultivated haplotype x wild haplotype)
#    4  best-matching donor per cultivated haplotype per window
#    5  row ordering (subgroup blocks, cosine distance + Ward linkage) and heatmap
#
#  Inputs (see README)
#    data/sv_matrix.tsv   CHROM POS ID SVTYPE SVLEN <hap1> <hap2> ...   (0/1)
#    data/metadata.tsv    Sample_ID  Group(domestic|wild)  [Subgroup]
#    data/genome.fai      samtools faidx index of the reference
#
#  Usage
#    bash chr15_haplotype_sharing.sh              # run missing steps
#    CHR=15 WIN=1000000 MIN_SV=5 bash chr15_haplotype_sharing.sh
#
#  Only existing outputs are re-used: delete results/ to start over.
###############################################################################
set -euo pipefail

## ------------------------------------------------------------------ settings
DATA=${DATA:-data}
OUT=${OUT:-results}
CHR=${CHR:-15}                 # chromosome to display (matches the matrix CHROM field)
WIN=${WIN:-1000000}            # window size (bp)
MIN_SV=${MIN_SV:-5}            # windows with fewer SVs in one haplotype are skipped
MIN_SVLEN=${MIN_SVLEN:-50}     # minimum SV length (bp)
MERGE_DIST=${MERGE_DIST:-1000} # SURVIVOR merging distance (bp)
THREADS=${THREADS:-8}
PYTHON=${PYTHON:-python3}

MATRIX=$DATA/sv_matrix.tsv
META=$DATA/metadata.tsv
FAI=$DATA/genome.fai
WINDOWS=$OUT/windows.bed
PW=$OUT/pairwise_jaccard.tsv
BEST=$OUT/best_donor.tsv
PDF=$OUT/chr15_haplotype_sharing.pdf

mkdir -p "$OUT"
log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH" >&2; exit 1; }; }

## ------------------------------------------------- 1. SV presence/absence matrix
if [ ! -s "$MATRIX" ]; then
  log "step 1: SV calls -> length filter -> merged matrix"
  need bcftools; need SURVIVOR
  ls "$DATA"/vcf/*.vcf > "$OUT/vcf_list.txt"
  SURVIVOR merge "$OUT/vcf_list.txt" "$MERGE_DIST" 1 1 1 0 "$MIN_SVLEN" "$OUT/merged.vcf" >/dev/null
  # GT -> binary presence (0/0, ./. and . are absence; anything else is presence)
  bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN[\t%GT]\n' "$OUT/merged.vcf" |
    awk -F'\t' 'BEGIN{OFS="\t"}
      {for(i=6;i<=NF;i++) $i=($i~/^(0\/0|\.\/\.|\.|\|)$/)?0:1; print}' > "$MATRIX"
  log "     matrix: $(($(wc -l < "$MATRIX")-1)) SV x $(head -1 "$MATRIX" | awk -F'\t' '{print NF-5}') haplotypes"
fi

## ------------------------------------------------------------- 2. window index
if [ ! -s "$WINDOWS" ]; then
  log "step 2: ${WIN}-bp windows from $FAI"
  awk -v w="$WIN" 'BEGIN{OFS="\t"}
    {for(s=0;s<$2;s+=w){e=s+w; if(e>$2)e=$2; print $1,s,e}}' "$FAI" > "$WINDOWS"
fi

## --------------------------------------------- 3. pairwise Jaccard per window
if [ ! -s "$PW" ]; then
  log "step 3: pairwise Jaccard (min $MIN_SV SVs per haplotype per window)"
  "$PYTHON" - "$MATRIX" "$META" "$WINDOWS" "$PW" "$CHR" "$MIN_SV" <<'PY'
import sys
import re
import numpy as np
import pandas as pd

mat, meta_f, win_f, out, chrom, min_sv = sys.argv[1:7]
min_sv = int(min_sv)

meta = pd.read_csv(meta_f, sep="\t")
meta["Sample_ID"] = meta["Sample_ID"].astype(str).str.strip()
meta["Group"] = meta["Group"].astype(str).str.strip().str.lower()
grp = dict(zip(meta["Sample_ID"], meta["Group"]))

M = pd.read_csv(mat, sep="\t", dtype=str)
haps = [c for c in M.columns[5:] if c not in ("", "nan")]
X = (M[haps].fillna("0").values != "0").astype(np.int8)      # SV present / absent
CHROM, POS = M["CHROM"].values, M["POS"].astype(np.int64).values
col = {h: i for i, h in enumerate(haps)}

def kind(h):
    return grp.get(h.rsplit(".hap", 1)[0], None)

dom = [h for h in haps if kind(h) == "domestic"]
wild = [h for h in haps if kind(h) == "wild"]

W = pd.read_csv(win_f, sep="\t", header=None, names=["chrom", "start", "end"])
norm = lambda s: (re.search(r"(\d+)$", str(s)) or [None, str(s)])[1]   # "chr15"/"GDT2T_hap1_chr15" -> "15"
W = W[W.chrom.map(norm) == norm(chrom)].reset_index(drop=True)
if W.empty:
    sys.exit(f"no windows for chromosome {chrom} in {win_f}")
CHROM_N = pd.Series(CHROM).map(norm).values

rows = []
for _, w in W.iterrows():
    sel = (CHROM_N == norm(w.chrom)) & (POS >= w.start) & (POS < w.end)
    if not sel.any():
        continue
    sub = X[sel]
    a_all = sub[:, [col[h] for h in dom]] == 1
    b_all = sub[:, [col[h] for h in wild]] == 1
    for i, d in enumerate(dom):
        a = a_all[:, i]
        if a.sum() < min_sv:
            continue
        for j, wd in enumerate(wild):
            b = b_all[:, j]
            if b.sum() < min_sv:
                continue
            union = np.sum(a | b)
            if union:
                rows.append((w.chrom, int(w.start), d, wd,
                             round(float(np.sum(a & b)) / union, 6),
                             int(a.sum()), int(b.sum())))

cols = ["Chrom", "Window_Start", "Domestic_Haplotype", "Wild_Donor",
        "Jaccard_Similarity", "n_SV_domestic", "n_SV_wild"]
pd.DataFrame(rows, columns=cols).to_csv(out, sep="\t", index=False)
print(f"    {len(rows)} cultivated x wild pairs written to {out}")
PY
fi

## ------------------------------------------- 4. best-matching donor per window
if [ ! -s "$BEST" ]; then
  log "step 4: best-matching wild donor per cultivated haplotype and window"
  "$PYTHON" - "$PW" "$BEST" <<'PY'
import sys
import pandas as pd
pw = pd.read_csv(sys.argv[1], sep="\t")
i = pw.groupby(["Domestic_Haplotype", "Window_Start"])["Jaccard_Similarity"].idxmax()
pw.loc[i].sort_values(["Window_Start", "Domestic_Haplotype"]).to_csv(
    sys.argv[2], sep="\t", index=False)
PY
fi

## ---------------------------------------------------- 5. clustering + heatmap
if [ ! -s "$PDF" ]; then
  log "step 5: clustering and heatmap"
  "$PYTHON" - "$BEST" "$META" "$PDF" "$CHR" <<'PY'
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

best_f, meta_f, pdf, chrom = sys.argv[1:5]
best = pd.read_csv(best_f, sep="\t")
meta = pd.read_csv(meta_f, sep="\t")
meta["Sample_ID"] = meta["Sample_ID"].astype(str).str.strip()
sub = meta["Subgroup"] if "Subgroup" in meta else meta["Group"]
sub = sub.astype(str).str.strip()
sg = dict(zip(meta["Sample_ID"], sub))

def sample(h):
    return h.rsplit(".hap", 1)[0]

doms = sorted(best.Domestic_Haplotype.unique())
wins = sorted(best.Window_Start.unique())
donors = sorted(best.Wild_Donor.unique())

J = best.pivot_table(index="Domestic_Haplotype", columns="Window_Start",
                     values="Jaccard_Similarity", aggfunc="max").reindex(index=doms, columns=wins)
D = best.pivot_table(index="Domestic_Haplotype", columns="Window_Start",
                     values="Wild_Donor", aggfunc="first").reindex(index=doms, columns=wins)

# feature vector: donor identity weighted by Jaccard, one block per window
di = {d: k for k, d in enumerate(donors)}
feat = np.zeros((len(doms), len(wins) * len(donors)))
for i in range(len(doms)):
    for j in range(len(wins)):
        d = D.iat[i, j]
        if isinstance(d, str):
            feat[i, j * len(donors) + di[d]] = J.iat[i, j]

# rows: cultivated subgroups in metadata order, clustered within subgroup
order, seen = [], []
for h in doms:
    g = sg.get(sample(h), "NA")
    if g not in seen:
        seen.append(g)
for g in seen:
    idx = [i for i, h in enumerate(doms) if sg.get(sample(h), "NA") == g]
    if len(idx) > 1:
        d = np.nan_to_num(pdist(feat[idx], metric="cosine"), nan=0.0, posinf=0.0)
        idx = [idx[k] for k in leaves_list(linkage(d, method="ward"))]
    order += idx
rows = [doms[i] for i in order]

# donor subgroup colours (expandable palette, grey = no match)
pal = plt.get_cmap("tab20").colors
groups = list(dict.fromkeys([sg.get(sample(d), "NA") for d in donors]))
cmap = {g: pal[k % 20] for k, g in enumerate(groups)}
NOFILL = (0.75, 0.75, 0.75)

n_r, n_c = len(rows), len(wins)
fig, ax = plt.subplots(figsize=(11, max(4.0, n_r * 0.20)))
for p, h in enumerate(rows):
    y = n_r - 1 - p
    for j, w in enumerate(wins):
        d = D.loc[h, w]
        if not isinstance(d, str):
            ax.add_patch(Rectangle((j, y), 1, 1, facecolor=NOFILL,
                                   edgecolor="white", linewidth=0.2))
            continue
        v = J.loc[h, w]
        ax.add_patch(Rectangle((j, y), 1, 1, facecolor=cmap[sg.get(sample(d), "NA")],
                               edgecolor="white", linewidth=0.2, alpha=max(v, 0.10)))
ax.set_xlim(0, n_c)
ax.set_ylim(0, n_r)
step = max(1, n_c // 10)
ax.set_xticks(np.arange(0, n_c, step) + 0.5)
ax.set_xticklabels([f"{int(wins[t]) / 1e6:g}" for t in range(0, n_c, step)], fontsize=7)
ax.set_yticks(np.arange(n_r) + 0.5)
ax.set_yticklabels([rows[n_r - 1 - k] for k in range(n_r)], fontsize=5)
ax.set_xlabel(f"Chromosome {chrom} position (Mb)", fontsize=8, labelpad=6)
ax.set_ylabel("Cultivated haplotype", fontsize=8)
ax.tick_params(length=2, width=0.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

handles = [Patch(facecolor=cmap[g], edgecolor="none", label=str(g)) for g in groups]
handles.append(Patch(facecolor=NOFILL, edgecolor="grey", linewidth=0.3, label="No match"))
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01),
           ncol=min(len(handles), 8), frameon=False, fontsize=6, handletextpad=0.3)
sm = ScalarMappable(norm=Normalize(0, 1), cmap=plt.cm.Greys)
fig.colorbar(sm, ax=ax, fraction=0.015, pad=0.02).set_label("Jaccard similarity", fontsize=7)
fig.tight_layout()
fig.savefig(pdf, bbox_inches="tight")
fig.savefig(pdf.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
print(f"    figure: {pdf}")
PY
fi

log "done"
