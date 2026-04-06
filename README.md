# Genome-Wide ATAC-seq Accessibility Analysis of Ohnolog Pairs in *Cyprinus carpio*

## Overview

This repository contains the code and supplementary data for a genome-wide quantitative analysis of chromatin accessibility asymmetry between ohnolog pairs in the allotetraploid common carp (*Cyprinus carpio*), extending prior work from an MSc project examining regulatory divergence during late zygotic development.

The original MSc analysis ([Gibson, 2024](https://github.com/rossjamesgibson21-bit/MSc-Project-Multi-Omic-Analysis-Allotetraploid-Cyprinid)) characterised promoter-proximal epigenetic asymmetry across a curated set of 75 ohnolog pairs using semi-quantitative visual scoring of ATAC-seq and histone modification tracks. That analysis found non-significant ATAC accessibility differences between subgenomes across all expression bias categories — a result that was statistically appropriate given the analytical design, which required visible ATAC signal at both subgenome copies for locus inclusion.

This project scales the ATAC analysis to the full complement of 7,497 expressed ohnolog pairs (baseMean ≥ 150) using a fully quantitative, automated pipeline. The genome-wide analysis reveals that the MSc non-significant result and the genome-wide significant result are complementary rather than contradictory — they reflect two distinct components of accessibility asymmetry that operate at different scales and through different mechanisms.

---

## Biological Question

> To what extent does promoter-proximal chromatin accessibility asymmetry between ohnolog copies correlate with expression bias across the full 7,497-pair ohnolog complement at Late Somite stage?

The curated 75-pair analysis found that ATAC accessibility differences between subgenomes were non-significant, while histone modification asymmetry (H3K27ac, H3K4me3) tracked expression bias robustly in SubgB-dominant pairs. Scaling the ATAC analysis to the full ohnolog complement allows us to test whether:

1. ATAC asymmetry is genuinely absent as a correlate of expression bias, or whether additional signal emerges at genome-wide scale
2. The direction of accessibility asymmetry is consistent with the subgenome dominance patterns identified in the histone mark analysis
3. Accessibility asymmetry varies across expression bias magnitude categories (Moderate, Strong, Extreme)

---

## Key Results

### Genome-Wide Accessibility Asymmetry

ATAC accessibility asymmetry significantly tracks expression bias direction across all 7,497 ohnolog pairs at Late Somite stage (Spearman ρ = 0.115, p = 1.35×10⁻²³) and Pre-Hatch stage (ρ = 0.105, p = 5.60×10⁻²⁰). The direction of asymmetry is consistent with subgenome dominance in both categories:

| Bias Category | n | Mean ATAC Asymmetry | Direction |
|---------------|---|---------------------|-----------|
| SubgA_dominant | 1,851 | +0.305 | SubgA more accessible |
| Balanced | 3,677 | −0.015 | Near-symmetric |
| SubgB_dominant | 1,969 | −0.444 | SubgB more accessible |

All pairwise Mann-Whitney comparisons survive BH-FDR correction (q < 0.05), with rank-biserial effect sizes of r = −0.164 (SubgA vs SubgB dominant), r = 0.104 (SubgB vs Balanced), and r = −0.068 (SubgA vs Balanced).

### Two Components of Accessibility Asymmetry

A restricted analysis separating pairs by ATAC detection status reveals two distinct components:

**1. Binary accessibility divergence** (one copy has no detectable peak): 1,299 pairs carry a disproportionate share of the expression bias correlation. These represent complete promoter accessibility divergence between subgenomes — the strongest form of regulatory asymmetry detectable from peak data alone.

**2. Graded accessibility divergence** (both copies have detectable peaks at different levels): 5,907 pairs show a weaker but still significant correlation with expression bias (ρ = 0.080, p = 6.28×10⁻¹⁰). Mean asymmetry values compress substantially in this subset (SubgB_dominant mean = −0.119 vs −0.444 in full dataset).

### Relationship to MSc Findings

The MSc curated analysis was statistically correct for what it measured. Locus selection criteria required visible ATAC signal at both subgenome copies, which effectively sampled from the graded divergence subset where the signal is genuinely weak. The non-significant MSc result reflected both the smaller sample size and the genuine biological distinction between binary and graded accessibility divergence — not simply insufficient power.

Together the two analyses characterise a layered regulatory hierarchy in *C. carpio* ohnologs:

- **Histone modification asymmetry** (H3K27ac, H3K4me3): strongest correlate of active expression bias, detectable at high effect size (r = 0.65–0.71) in curated pairs with visible signal at both copies
- **Binary accessibility divergence**: strong signal genome-wide, represents complete promoter accessibility loss in one subgenome copy
- **Graded accessibility divergence**: weaker signal in pairs where both copies are accessible, consistent with the MSc non-significant ATAC result

---

## Data

### Input Files

| File | Description | Source |
|------|-------------|--------|
| `DegListWithExNum_NEW.tsv` | DESeq2 expression data for 9,581 ohnolog pairs across Late Somite and Pre-Hatch stages, including Ensembl gene IDs, scaffold coordinates, exon counts, L2FC and padj values | Provided by collaborating wet lab (A. Jimenez-Gonzalez) |
| `LateSomATAC.narrowPeak` | ATAC-seq narrow peak calls for Late Somite stage, aligned to cypCar4 (WAG4.0) assembly | Provided by collaborating wet lab |

> **Note**: Raw data files are not included in this repository as they derive from unpublished work by the collaborating researcher. The pipeline is fully reproducible given these input files.

### Genome Assembly

All coordinates use the *Cyprinus carpio carpio* WAG4.0 assembly (UCSC: cypCar4; NCBI: GCA_905221575.1), with scaffold identifiers in CAJNDQ format. Gene coordinates are retrieved programmatically from the Ensembl REST API using ENSCCRG Ensembl gene IDs.

---

## Pipeline

### Dependencies

```
python >= 3.9
pandas >= 2.1.4
pybedtools >= 0.10.0
scipy
numpy
requests
bedtools >= 2.31.1
openpyxl
statsmodels
```

Install Python dependencies:

```bash
pip install pandas pybedtools scipy numpy requests openpyxl statsmodels
```

Bedtools must be installed separately and available on your PATH:

```bash
bedtools --version
```

### Running the Pipeline

1. Clone this repository
2. Place input data files in the expected locations (or update paths in the notebook header)
3. Open `CarpOhnolog_ATAC_GenomeWide.ipynb` in Jupyter and run all cells

On first run, the notebook queries the Ensembl REST API for TSS coordinates for all ENSCCRG IDs. Results are cached locally in `atac_analysis/ensembl_coords_cache.json` — subsequent runs use the cache and are substantially faster.

### Pipeline Steps

```
1.  Load DEG list → filter to baseMean ≥ 150 (7,497 pairs)
2.  Query Ensembl REST API → TSS coordinates for all ENSCCRG IDs (batched, cached)
3.  Build promoter window BED files → ±2kb around TSS for each subgenome copy
4.  Intersect with ATAC narrowPeak → pybedtools left-outer join
5.  Quantify accessibility → max fold enrichment per promoter window
6.  Calculate ATAC asymmetry → log2(SubgA / SubgB fold enrichment + pseudocount)
7.  Correlate asymmetry with L2FC → Spearman correlation by stage
8.  Mann-Whitney U tests by bias category with BH-FDR correction and rank-biserial r
9.  Restricted analysis → binary vs graded accessibility divergence components
10. Sensitivity analysis → window size stability (±1kb, ±2kb, ±5kb)
```

---

## Outputs

All outputs are written to `~/Desktop/atac_analysis/` by default (configurable in notebook header).

| File | Description |
|------|-------------|
| `ohnolog_atac_asymmetry_results.csv` | Per-pair results: ATAC fold enrichment per subgenome, asymmetry score, L2FC, bias category |
| `atac_asymmetry_summary_by_bias.csv` | Mean and median ATAC asymmetry by expression bias category |
| `promoters_SubgA.bed` | Promoter window BED file for SubgA copies |
| `promoters_SubgB.bed` | Promoter window BED file for SubgB copies |
| `ensembl_coords_cache.json` | Cached Ensembl REST API coordinate results |
| `outputs/*.png` | Figures: scatter, violin, magnitude boxplot, mean enrichment bar, regulatory comparison |

---

## Key Analytical Choices

**Promoter window**: ±2kb around TSS. TSS-centred windows provide a consistent and reproducible promoter definition across all 7,497 pairs. Sensitivity analyses at ±1kb and ±5kb are included in Section 10.3 of the notebook.

**Quantification**: Maximum fold enrichment within the promoter window where multiple peaks overlap. Captures the strongest accessibility signal at each locus.

**Asymmetry metric**: log2(SubgA fold enrichment + ε) / (SubgB fold enrichment + ε), where ε = 0.1 is a pseudocount. Positive values indicate greater SubgA accessibility; negative values indicate greater SubgB accessibility. Directly comparable to the L2FC expression bias metric.

**Statistical framework**: Spearman correlation for continuous association with L2FC; Mann-Whitney U with BH-FDR correction and rank-biserial r for bias category comparisons.

---

## Relationship to Prior Work

This analysis directly extends:

> Gibson, R.J. (2024). *Multi-omic Analysis of an Allotetraploid Cyprinid Reveals Regulatory Divergence of Ohnologs During Late Zygotic Development*. MSc dissertation, University of Birmingham.
> Repository: [MSc-Project-Multi-Omic-Analysis-Allotetraploid-Cyprinid](https://github.com/rossjamesgibson21-bit/MSc-Project-Multi-Omic-Analysis-Allotetraploid-Cyprinid)

The findings from this extended analysis are intended to complement ongoing wet lab work by A. Jimenez-Gonzalez (University of Birmingham) and may contribute to a collaborative publication.

---

## Author

Ross Gibson  
MSc Bioinformatics, University of Birmingham  
GTCS Registered Biology Teacher  
Edinburgh, Scotland

---

## License

Code is released under the MIT License. Data files are not included; please contact the author regarding data access.
