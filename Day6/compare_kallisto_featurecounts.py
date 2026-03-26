#!/usr/bin/env python3
import argparse
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def strip_version(identifier: str) -> str:
    return re.sub(r"\.\d+$", "", str(identifier))


def parse_gtf_transcript_to_gene(gtf_path: str) -> pd.DataFrame:
    mapping = []
    with open(gtf_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attributes = parts[8]
            tx_match = re.search(r'transcript_id "([^"]+)"', attributes)
            gene_match = re.search(r'gene_id "([^"]+)"', attributes)
            if tx_match and gene_match:
                transcript_id = strip_version(tx_match.group(1))
                gene_id = strip_version(gene_match.group(1))
                mapping.append((transcript_id, gene_id))

    df_map = pd.DataFrame(mapping, columns=["transcript_id", "gene_id"]).drop_duplicates()
    if df_map.empty:
        raise ValueError("Nie znaleziono mapowania transcript_id -> gene_id w pliku GTF.")
    return df_map


def load_kallisto_gene_counts(kallisto_tsv: str, df_map: pd.DataFrame) -> pd.DataFrame:
    df_k = pd.read_csv(kallisto_tsv, sep="\t")
    if "target_id" not in df_k.columns or "est_counts" not in df_k.columns:
        raise ValueError("abundance.tsv musi zawierać kolumny target_id oraz est_counts.")

    df_k["transcript_id"] = df_k["target_id"].map(strip_version)
    df_k = df_k.merge(df_map, on="transcript_id", how="left")
    df_k = df_k.dropna(subset=["gene_id"])

    df_gene = (
        df_k.groupby("gene_id", as_index=False)["est_counts"]
        .sum()
        .rename(columns={"est_counts": "kallisto_count"})
    )
    return df_gene


def load_featurecounts_gene_counts(featurecounts_txt: str) -> pd.DataFrame:
    df_fc = pd.read_csv(featurecounts_txt, sep="\t", comment="#")
    if df_fc.shape[1] < 7:
        raise ValueError("Plik featureCounts ma nieoczekiwany format.")

    sample_col = df_fc.columns[-1]
    df_gene = df_fc[["Geneid", sample_col]].copy()
    df_gene["Geneid"] = df_gene["Geneid"].map(strip_version)
    df_gene = df_gene.rename(columns={"Geneid": "gene_id", sample_col: "featurecounts_count"})
    return df_gene


def main():
    parser = argparse.ArgumentParser(
        description="Porownanie Kallisto (agregacja do genu) vs featureCounts."
    )
    parser.add_argument("--gtf", required=True, help="Sciezka do pliku GTF")
    parser.add_argument("--kallisto", required=True, help="Sciezka do kallisto abundance.tsv")
    parser.add_argument("--featurecounts", required=True, help="Sciezka do pliku featureCounts")
    parser.add_argument(
        "--out-prefix",
        default="Day6/data/reads/kallisto_vs_featurecounts",
        help="Prefiks plikow wyjsciowych",
    )
    args = parser.parse_args()

    df_map = parse_gtf_transcript_to_gene(args.gtf)
    df_k_gene = load_kallisto_gene_counts(args.kallisto, df_map)
    df_fc_gene = load_featurecounts_gene_counts(args.featurecounts)

    df_merge = df_fc_gene.merge(df_k_gene, on="gene_id", how="inner")
    if df_merge.empty:
        raise ValueError("Brak wspolnych gene_id po merge Kallisto i featureCounts.")

    rho, pval = spearmanr(df_merge["featurecounts_count"], df_merge["kallisto_count"])

    df_merge["log2_featurecounts"] = np.log2(df_merge["featurecounts_count"] + 1)
    df_merge["log2_kallisto"] = np.log2(df_merge["kallisto_count"] + 1)

    merged_path = f"{args.out_prefix}.merged.tsv"
    stats_path = f"{args.out_prefix}.stats.txt"
    plot_path = f"{args.out_prefix}.scatter.png"

    df_merge.to_csv(merged_path, sep="\t", index=False)

    with open(stats_path, "w", encoding="utf-8") as handle:
        handle.write(f"n_common_genes\t{len(df_merge)}\n")
        handle.write(f"spearman_rho\t{rho}\n")
        handle.write(f"spearman_pvalue\t{pval}\n")

    plt.figure(figsize=(7, 7))
    plt.scatter(df_merge["log2_featurecounts"], df_merge["log2_kallisto"], s=8, alpha=0.6)
    plt.xlabel("log2(featureCounts + 1)")
    plt.ylabel("log2(Kallisto est_counts + 1)")
    plt.title(f"Kallisto vs featureCounts (Spearman rho={rho:.3f})")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150)

    print(f"Merged table: {merged_path}")
    print(f"Stats: {stats_path}")
    print(f"Plot: {plot_path}")


if __name__ == "__main__":
    main()
