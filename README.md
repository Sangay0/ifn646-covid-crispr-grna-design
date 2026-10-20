# IFN646: COVID-19 RNA-seq, CRISPRi guide design and variant analysis

Group project (IFN646) by Amy O'Toole, Sangay Norbu and Abhishek Chacko.

The project looks for genes disrupted by SARS-CoV-2 infection, designs CRISPRi guide RNAs for one of them (FEV), and checks the guides against human genetic variants.

## What each task does

| Task | What it does | Where |
|---|---|---|
| 1 | Differential expression of SARS-CoV-2 vs mock samples with DESeq2 and limma-voom; consensus genes | `task1/` (R scripts), `results/` (outputs) |
| 2 | Designs CRISPRi guides in exon 1 of FEV with CRISPOR and filters them | `Task2_gRNA_design.ipynb`, `task2/` |
| 3 | Screens the guides and their off-target sites against sampled 1000 Genomes variants (chr 1 and 2) | `task3/` |

## Folder layout

```
README.md
gse159717_rnaseq_deseq_5dpi_counts_raw.tsv   RNA-seq raw counts (Task 1 input)
meta.csv                                     sample table (Task 1 input)
task1/        R scripts for Task 1
results/      Task 1 outputs (gene tables, plots)
Task2_gRNA_design.ipynb                      Task 2 notebook
task2/        Task 2 outputs: guides, off-targets, CRISPOR tables, figures, requirements.txt
task3/        Task 3 notebook and outputs
```

## Data not included

Too large for GitHub, so download them and put them in the places shown:

- GRCh38 genome `GCF_000001405.26_GRCh38_genomic.fna.gz` and annotation `GCF_000001405.26_GRCh38_genomic.gff.gz` from NCBI: repository root (needed for Task 2).
- `sampled_chr1.vcf.zip` and `sampled_chr2.vcf.zip` (1000 Genomes subset, 7 individuals): `task3/input/` (needed for Task 3).

## How to run

Run everything from the repository root.

1. **Task 1:** R 4.5.2 with DESeq2, limma, edgeR, apeglm, ashr, pheatmap and ggplot2. Run a script from `task1/` (for example `source("task1/task1_stepwise.R")`). Results are written to `results/`.
2. **Task 2:** Python 3.10. Run `pip install -r task2/requirements.txt`, then run all cells of `Task2_gRNA_design.ipynb`. The CRISPOR results are saved in `task2/`, so reruns give the same numbers without internet. To remake the report figures, run `python task2/make_report_figures.py`.
3. **Task 3:** See `task3/README.md`.

## Main results

- **Task 2:** FEV was chosen (up-regulated, DESeq2 log2 fold change 3.8). Two recommended guides (285forw, 165forw) and one backup (259forw).
- **Task 3:** No variant overlaps the three target sites. One candidate off-target site of 259forw overlaps a C→G variant. Only 4 of 22 candidate off-target sites could be checked (VCFs for chr 1 and 2 only).

