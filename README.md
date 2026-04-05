# Multi-Omic Analysis of Ohnolog Regulatory Divergence in Allotetraploid *Cyprinus carpio*

## Overview

This repository contains the analysis code, processed data, and figures from an MSc Bioinformatics thesis project (University of Birmingham, 2024–2025) investigating how duplicated gene pairs (ohnologs) diverge in their regulatory control following whole-genome duplication (WGD) in common carp (*Cyprinus carpio*).

Common carp is an allotetraploid species whose genome arose from the merger of two ancestral diploid genomes approximately 12.4 million years ago. This created thousands of ohnolog pairs — duplicate genes derived from the same ancestral gene, one on each subgenome (A and B). A central question in polyploid genome evolution is whether these duplicates maintain equivalent expression and regulation, or whether one copy becomes dominant while the other degrades or specialises. This project addresses that question by integrating three independent data types across two developmental stages.

## Biological Question

After whole-genome duplication, do ohnolog pairs show symmetric or asymmetric regulatory divergence? Specifically:

- Do ohnolog copies on subgenomes A and B show balanced or biased expression during embryonic development?
- Is expression bias correlated with differences in chromatin accessibility (ATAC-seq) and histone modification (ChIP-seq) at promoter regions?
- Are functionally coherent gene categories (e.g. ribosome biogenesis, signalling, transcriptional regulation) enriched among ohnologs with particular divergence trajectories?

## Data Types

The analysis integrates three omic layers from curated ohnolog pairs across **Late Somite** and **Pre-Hatch** developmental stages:

| Data Type | What It Measures | Application in This Project |
|-----------|-----------------|---------------------------|
| **RNA-seq** | mRNA expression levels | Quantifying expression bias between ohnolog A and B copies |
| **ATAC-seq** | Chromatin accessibility | Comparing promoter-level openness between ohnolog pairs |
| **Histone ChIP-seq** | Histone modification marks | Assessing epigenetic state at ohnolog regulatory regions |

## Analytical Approach

1. **Expression bias classification**: Ohnolog pairs were classified into divergence categories based on the magnitude and direction of expression bias between subgenome A and B copies, using DESeq2-derived log2 fold-changes and Wilcoxon signed-rank tests. Categories range from *Balanced* (no significant bias) through *Moderate* to *Extreme* subgenome dominance.

2. **Chromatin accessibility comparison**: ATAC-seq signal at promoter regions of ohnolog pairs was quantified and compared between A and B copies to assess whether expression bias is reflected in differential chromatin openness.

3. **Histone mark integration**: ChIP-seq data for activating and repressive histone marks were integrated to determine whether divergent ohnologs show corresponding epigenetic asymmetry.

4. **Functional enrichment and network analysis**: Gene Ontology enrichment and protein–protein interaction (PPI) network analysis (via STRING/zebrafish ortholog mapping) were performed for each divergence class to identify functional themes among biased ohnologs.

## Repository Structure

```
├── README.md
├── notebooks/
│   ├── 01_DE_Analysis.ipynb              # Differential expression and bias classification
│   ├── 02_Ohnolog_Regulatory_Analysis.ipynb  # ATAC-seq and ChIP-seq integration
│   └── 03_Enrichment_and_Networks.ipynb  # GO enrichment and PPI network analysis
├── data/
│   ├── ohnolog_pairs/
│   │   ├── filtered_ohnologs_basemean_150.csv
│   │   └── filtered_ohnologs_Lsom_L2FC_2_to_7.csv
│   ├── divergence_classes/
│   │   ├── Balanced_LSom.csv
│   │   ├── Balanced_PHatch.csv
│   │   ├── Moderate_SubgA_Dom_LSom.csv
│   │   ├── Moderate_SubgA_Dom_PHatch.csv
│   │   ├── Moderate_SubgB_Dom_LSom.csv
│   │   ├── Moderate_SubgB_Dom_PHatch.csv
│   │   ├── Strong_SubgA_Dom_LSom.csv
│   │   ├── Strong_SubgA_Dom_PHatch.csv
│   │   ├── Strong_SubgB_Dom_LSom.csv
│   │   ├── Strong_SubgB_Dom_PHatch.csv
│   │   ├── Extreme_SubgA_Dom_LSom.csv
│   │   ├── Extreme_SubgA_Dom_PHatch.csv
│   │   ├── Extreme_SubgB_Dom_LSom.csv
│   │   ├── Extreme_SubgB_Dom_PHatch.csv
│   │   ├── No_Significant_Bias_LSom.csv
│   │   └── No_Significant_Bias_PHatch.csv
│   ├── zebrafish_mappings/
│   │   └── Zebrafish_mapping.tsv
│   ├── intersection_results/
│   │   ├── Extreme_SubgA_Dom_LSom_intersections.csv
│   │   ├── Extreme_SubgA_Dom_PHatch_intersections.csv
│   │   ├── Extreme_SubgB_Dom_LSom_intersections.csv
│   │   └── Extreme_SubgB_Dom_PHatch_intersections.csv
│   └── summary/
│       ├── Combined_wilcoxon_summary.csv
│       ├── DEGList_EXNum_Compile_PEAK.xlsx
│       ├── DEGList_EXNum_MEAN_PHatch.tsv
│       └── DegListWithExNum.xlsx
├── figures/
│   ├── igv_tracks/
│   │   ├── dnmt1_chr2_20kb.png
│   │   ├── dnmt1_chr4_yscale.png
│   │   ├── mib2_chr18_5kb.png
│   │   ├── mib2_chr24_20kb.png
│   │   ├── tgfb2_chr3_20kb.png
│   │   └── tgfb2_chr48_20kb.png
│   ├── ppi_networks/
│   │   ├── Extreme-SubgA_Dom_LSom_PPI_Main_Network.png
│   │   ├── Extreme-SubgA_Dom_PHatch_Main_PPI_Network.png
│   │   ├── Extreme-SubgB_Dom_LSom_Main_PPI_Network.png
│   │   └── Extreme-SubgB_Dom_PreHatch_Main_PPI_Network.png
│   ├── enrichment/
│   │   └── Enriched_Pathways.png
│   ├── expression/
│   │   ├── Correlation_Matrix_Task2_RE_and_Ohnolog_Analysis.png
│   │   └── Comparison_of_TF_Families.png
│   └── ribosome_biogenesis/
│       ├── Balanced_LSom_Ribosome_Biogenesis_Network.png
│       ├── SubgA_Dom_LSom_Ribosome_Biogenesis_Network.png
│       └── SubgB_Dom_LSom_Ribosome_Biogenesis_Network.png
└── results/
    └── zebrafish_id_lists/
        ├── Balanced_LSom_zebrafish_IDs.txt
        ├── Balanced_PHatch_zebrafish_IDs.txt
        ├── Extreme_SubgA_Dom_LSom_zebrafish_IDs.txt
        ├── Extreme_SubgA_Dom_PHatch_zebrafish_IDs.txt
        ├── Extreme_SubgB_Dom_LSom_zebrafish_IDs.txt
        ├── Extreme_SubgB_Dom_PHatch_zebrafish_IDs.txt
        ├── Moderate_SubgA_Dom_LSom_zebrafish_IDs.txt
        ├── Moderate_SubgA_Dom_PHatch_zebrafish_IDs.txt
        ├── Moderate_SubgB_Dom_LSom_zebrafish_IDs.txt
        ├── Moderate_SubgB_Dom_PHatch_zebrafish_IDs.txt
        ├── No_Significant_Bias_LSom_zebrafish_IDs.txt
        ├── No_Significant_Bias_PHatch_zebrafish_IDs.txt
        ├── Strong_SubgA_Dom_LSom_zebrafish_IDs.txt
        ├── Strong_SubgA_Dom_PHatch_zebrafish_IDs.txt
        ├── Strong_SubgB_Dom_LSom_zebrafish_IDs.txt
        └── Strong_SubgB_Dom_PHatch_zebrafish_IDs.txt
```

## Key Findings

- Ohnolog pairs show a spectrum of expression divergence, from fully balanced to extreme subgenome dominance, with the distribution shifting between Late Somite and Pre-Hatch stages.
- Extreme subgenome-B-dominant ohnologs at Late Somite stage are enriched for ribosome biogenesis and translational machinery, suggesting subgenome-specific control of core biosynthetic processes during early development.
- PPI network analysis of divergence classes reveals functionally coherent modules, with distinct biological processes enriched in A-dominant versus B-dominant ohnolog sets.
- IGV visualisation of candidate ohnologs (*dnmt1*, *mib2*, *tgfb2*) across homeologous chromosomes shows differential ATAC-seq peak profiles consistent with the expression bias observed in the RNA-seq data.

## Tools and Dependencies

**Languages**: Python 3, R

**Key Python libraries**: pandas, NumPy, SciPy, matplotlib, seaborn, Jupyter

**Key R packages**: DESeq2, GenomicRanges, clusterProfiler, ggplot2

**Bioinformatics tools**: STAR (alignment), MACS2 (ATAC-seq peak calling), deepTools (signal visualisation), samtools, bedtools, featureCounts

**External resources**: STRING database (PPI networks via zebrafish orthologs), Gene Ontology, IGV (genome browser)

## Reference Genome

All analyses use the *Cyprinus carpio* chromosome-level genome assembly with subgenome A and B assignments as described in Xu et al. (2019).

## Zebrafish Ortholog Mapping

Because functional annotation and PPI databases are more complete for zebrafish (*Danio rerio*) than for carp, ohnolog gene IDs were mapped to zebrafish orthologs for GO enrichment and STRING network analysis. The mapping file (`data/zebrafish_mappings/Zebrafish_mapping.tsv`) contains these correspondences.

## Status

This project formed the basis of an MSc thesis completed in October 2025. Manuscript preparation is in progress, with a planned extension to genome-wide ATAC-seq analysis beyond the initial curated ohnolog set.

## Author

**Ross Gibson** — MSc Bioinformatics, University of Birmingham
- Email: rossjamesgibson21@gmail.com
- Location: Edinburgh, UK

## Licence

This repository is shared for academic and portfolio purposes. Please contact the author before reusing data or code in publications.
