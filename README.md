# Genome-Wide ATAC-seq Accessibility Analysis of Ohnolog Pairs in *Cyprinus carpio*

## Overview

This repository contains the code and supplementary data for a genome-wide quantitative analysis of chromatin accessibility asymmetry between ohnolog pairs in the allotetraploid common carp (*Cyprinus carpio*), extending prior work from an MSc project examining regulatory divergence during late zygotic development.

The original MSc analysis ([Gibson, 2024](https://github.com/rossjamesgibson21-bit/MSc-Project-Multi-Omic-Analysis-Allotetraploid-Cyprinid)) characterised promoter-proximal epigenetic asymmetry across a curated set of 75 ohnolog pairs using semi-quantitative visual scoring of ATAC-seq and histone modification tracks. That analysis identified significant activation-mark asymmetry in SubgB-dominant pairs but non-significant ATAC accessibility differences across all bias categories — a finding that may reflect the limited power of the curated subset and the semi-quantitative methodology rather than a genuine absence of accessibility divergence.

This project addresses that limitation directly by applying a fully quantitative, automated pipeline to the complete set of 7,497 expressed ohnolog pairs (baseMean ≥ 150), using ATAC-seq narrow peak data and programmatic bedtools intersection to characterise accessibility asymmetry at scale.

---

## Biological Question

> To what extent does promoter-proximal chromatin accessibility asymmetry between ohnolog copies correlate with expression bias at the genome-wide scale in *C. carpio* Late Somite stage embryos?

The curated 75-pair analysis suggested that ATAC accessibility differences between subgenomes were non-significant, while histone modification asymmetry (H3K27ac, H3K4me3) tracked expression bias robustly in SubgB-dominant pairs. Scaling the ATAC analysis to the full ohnolog complement allows us to test whether:

1. ATAC asymmetry is genuinely absent as a correlate of expression bias, or whether it emerges at genome-wide scale with sufficient statistical power
2. The direction of accessibility asymmetry is consistent with the subgenome dominance patterns identified in the histone mark analysis
3. Accessibility asymmetry varies across expression bias magnitude categories (Moderate, Strong, Extreme)

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
```

Install Python dependencies:

```bash
pip install pandas pybedtools scipy numpy requests openpyxl
```

Bedtools must be installed separately and available on your PATH:

```bash
bedtools --version
```

### Running the Pipeline

1. Clone this repository
2. Place input data files in the expected locations (or update paths in the script header)
3. Run the pipeline:

```bash
python3 atac_pipeline.py
```

On first run, the script queries the Ensembl REST API for TSS coordinates for all ENSCCRG IDs. Results are cached locally in `atac_analysis/ensembl_coords_cache.json` — subsequent runs use the cache and are substantially faster.

### Pipeline Steps

```
1. Load DEG list → filter to baseMean ≥ 150 (7,497 pairs)
2. Query Ensembl REST API → TSS coordinates for all ENSCCRG IDs (batched, cached)
3. Build promoter window BED files → ±2kb around TSS for each subgenome copy
4. Intersect with ATAC narrowPeak → pybedtools left-outer join
5. Quantify accessibility → max fold enrichment per promoter window
6. Calculate ATAC asymmetry → log2(SubgA / SubgB fold enrichment + pseudocount)
7. Correlate asymmetry with L2FC → Spearman correlation, Mann-Whitney U by bias category
```

---

## Outputs

All outputs are written to `~/Desktop/atac_analysis/` by default (configurable in script header).

| File | Description |
|------|-------------|
| `ohnolog_atac_asymmetry_results.csv` | Per-pair results: ATAC fold enrichment per subgenome, asymmetry score, L2FC, bias category |
| `atac_asymmetry_summary_by_bias.csv` | Mean and median ATAC asymmetry by expression bias category |
| `promoters_SubgA.bed` | Promoter window BED file for SubgA copies |
| `promoters_SubgB.bed` | Promoter window BED file for SubgB copies |
| `ensembl_coords_cache.json` | Cached Ensembl REST API coordinate results |

---

## Key Analytical Choices

**Promoter window**: ±2kb around TSS. This is wider than the ±5kb PADRE-anchored window used in the original MSc analysis, but without PADRE coordinates available for all 7,497 pairs, TSS-centred windows provide a consistent and reproducible promoter definition. Sensitivity analyses with ±1kb and ±5kb windows are straightforward to implement by modifying `PROMOTER_WINDOW` in the script.

**Quantification**: Maximum fold enrichment within the promoter window, where multiple peaks overlap. This captures the strongest accessibility signal rather than cumulative peak density, which is more interpretable in a promoter context.

**Asymmetry metric**: log2(SubgA fold enrichment + ε) / (SubgB fold enrichment + ε), where ε = 0.1 is a pseudocount to handle zero values. Positive values indicate greater SubgA accessibility; negative values indicate greater SubgB accessibility. This is directly comparable to the L2FC metric (log2(SubgA / SubgB expression)).

**Limitations**: TSS coordinates are derived from Ensembl annotations which may not perfectly match the cypCar4 assembly used for ATAC alignment. Coordinate discrepancies between assemblies could introduce noise in the intersection step. Additionally, without the original PADRE BED file, promoter-proximal regulatory elements cannot be anchored to ATAC summits as in the original analysis.

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
