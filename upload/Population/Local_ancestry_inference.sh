import pandas as pd
import numpy as np

# =========================
# 1. 读取VCF
# =========================

cols = [
"CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
"genome.hap1","genome.hap2"]

df = pd.read_csv("merged.vcf", sep="\t", comment="#", names=cols)

# =========================
# 2. GT → allele
# =========================

def gt_to_allele(gt, ref, alt):
    if gt in ["0|0","0/0"]:
        return ref
    elif gt in ["1|1","1/1"]:
        return alt
    elif gt in ["1|0","0|1"]:
        return "HET"
    else:
        return None

def match(a, b):
    if a is None or b is None:
        return 0
    return 1 if a == b  else 0

# =========================
# 3. 50kb window
# =========================

df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
df = df.dropna(subset=["POS"])
df["POS"] = df["POS"].astype(int)

df["window"] = (df["POS"] // 50000) * 50000

# =========================
# 4. 样本（每个hap单独算）
# =========================

samples = ["genome"]
haps = ["hap1","hap2"]

# =========================
# 5. similarity计算
# =========================

results = []

for sample in samples:
    for hap in haps:

        sample_name = f"{sample}.{hap}"

        for (chrom, w), g in df.groupby(["CHROM","window"]):

            ancestry1_score = 0
            ancestry2_score = 0
            n = 0

            for _, row in g.iterrows():

                s = gt_to_allele(row[sample_name], row["REF"], row["ALT"])

                if s is None:
                    continue

                # ancestry1
                ancestry1_1 = gt_to_allele(row["ancestry1.hap1"], row["REF"], row["ALT"])
                ancestry1_1 = gt_to_allele(row["ancestry1.hap2"], row["REF"], row["ALT"])

                # ancestry2
                ancestry2_1 = gt_to_allele(row["ancestry2.hap1"], row["REF"], row["ALT"])
                ancestry2_2 = gt_to_allele(row["ancestry2.hap2"], row["REF"], row["ALT"])

                # similarity（取最相似hap）
                ancestry1_score += max(match(s, ancestry1_1), match(s, ancestry1_2))
                ancestry2_score += max(match(s, ancestry2_1), match(s, ancestry2_2))

                n += 1

            if n == 0:
                continue

            msr_sim = ancestry1_score / n
            ra16_sim = ancestry2_score / n

            results.append([
                sample_name, chrom, w,
                ancestry1_sim, ancestry2_sim
            ])

# =========================
# 6. 转DataFrame
# =========================

res = pd.DataFrame(results,
    columns=["sample","chrom","window","ancestry1_sim","ancestry2_sim"])

# =========================
# 7. ancestry判定（关键）
# =========================

def assign(row, threshold=0):

    diff = abs(row.ancestry2_sim - row.ancestry1_sim)

    if diff < threshold:
        return "unclear"
    elif row.ancestry2_sim > row.ancestry1_sim:
        return "ancestry2"
    else:
        return "ancestry1"

res["ancestry"] = res.apply(assign, axis=1)

# =========================
# 8. confidence（可选）
# =========================

res["confidence"] = abs(res["ancestry1_sim"] - res["ancestry2_sim"])

# =========================
# 9. genome比例
# =========================

summary = (
    res.groupby("sample")["ancestry"]
    .value_counts(normalize=True)
    .unstack()
    .fillna(0)
)

# =========================
# 10. 输出
# =========================

res.to_csv("window_ancestry.tsv", sep="\t", index=False)
summary.to_csv("ancestry_summary.tsv", sep="\t")

print("\n=== Genome-wide ancestry proportion ===")
print(summary)
