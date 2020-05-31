### Disclaimer

This set of scripts was developed and updated over a period of time. It is shared for reference only and by no means intended to be used as is. The input files are the result of the previous steps in the pipeline for processing NGS data. The pipeline developed by Fengyuan Hu, Sri Deevi and Olga Shamardina.

The scripts and the input can certainly be improved, optimised and structured better.

### Citation

The CNV-calling method, plot description and definitions can be found in the paper:

* Simeoni I, Shamardina O, Deevi SV, et al. GRID – Genomics of Rare Immune Disorders: a highly sensitive and specific diagnostic gene panel for patients with primary immunodeficiencies. 2019, bioRxiv. [doi:10.1101/431544](https://doi.org/10.1101/431544)

Please cite this paper if you adopt ideas from this repository.


### Make CNV plots

First, we calculate the average coverage for each sample:

```sh
Rscript mean_coverage.R /path/to/conf.R
```

Also we precalculate percentiles for each gene:

```sh
for gene in gene1 gene2 gene3; do Rscript gene_stats.R /path/to/conf.R $gene; done
```

The script to make the plots have three arguments: configuration file, sample ID and gene name:

```sh
Rscript plotCNV.R /path/to/conf.R S1 G1
```

We will normally call it in a loop or in `parallel` for all identified CNVs in each sample.


The script `plotCNV-dense-pdf.R` is a modified version for long genes. It's not used routinely.
