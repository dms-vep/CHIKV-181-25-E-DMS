# CHIKV intrahost variant analysis

A `Snakemake` pipeline that fetches CHIKV sequencing runs from the Sequence Read Archive (SRA) and identifies viral variation in each sample.

Analysis by Will Hannon

## Organization

The contents of `intrahost/` are organized as follows:

```bash
.
├── configuration # <------ Pipeline configuration
├── data # <--------------- Accessions and metadata 
├── environment.yaml
├── README.md
├── results # <------------ Sequences, alignments, variants, etc...
├── run_analysis.bash
└── workflow # <----------- Code to run the analysis
```

The `Snakemake` pipeline (and associated code) that creates the `Nextstrain` tree is in the [`workflow/`](workflow/) directory:

```bash
workflow/
├── notebooks # <---- Jupyter notebooks
├── profiles # <----- Cluster resource configuration
├── scripts # <------ Python scripts
└── Snakefile # <---- Snakemake workflow
```
