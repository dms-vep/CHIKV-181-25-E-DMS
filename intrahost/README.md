# CHIKV intrahost variant analysis

A `Snakemake` pipeline that fetches CHIKV sequencing runs from the Sequence Read Archive (SRA) and identifies viral variation in each sample relative to the DMS library strain [MW473668](https://www.ncbi.nlm.nih.gov/nuccore/MW473668).

Analysis by Will Hannon

## Organization

The contents of `intrahost/` are organized as follows:

```bash
.
├── configuration # <------ Pipeline configuration
├── data # <--------------- Accessions and metadata 
├── environment.yaml # <------- Conda environment
├── README.md
├── results # <------------ Sequences, alignments, variants, etc...
├── run_analysis.bash
└── workflow # <----------- Code to run the analysis
```

The `Snakemake` pipeline (and associated code) that identifies viral variants is located in the [`workflow/`](workflow/) directory:

```bash
workflow/
├── notebooks # <---- Jupyter notebooks
├── profiles # <----- Cluster resource configuration
├── rules # <----- Snakemake rule files
└── Snakefile # <---- Snakemake workflow
```

## Identifying CHIKV Sequencing Runs

The NCBI tabulates the taxonomic breakdown of each sequencing run in the SRA with the [SRA Taxonomy Tool](https://www.ncbi.nlm.nih.gov/sra/docs/sra-taxonomy-analysis-tool/). This information is [stored on the cloud](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-taxonomy-analysis-table/) and can be queried with Google's BigQuery UI. We used this information to compile a list of SRA accessions for all sequencing runs that contain a minimum number of reads from Chikungunya virus. As of 06/05/25, this consisted of 9606 runs with some number of reads containing CHIKV kmers.

We used Google's BigQuery UI by logging into the [Google Cloud Console](https://console.cloud.google.com/welcome) and selecting `Run a query in BigQuery`. This should take you to something called 'BigQuery Studio' where you'll be prompted to log in and create a new project. Once you've created a project, select `SQL Query` under 'Create new', which will open an empty text file in the browser window. Enter the following SQL command per [these instructions from the NCBI](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-taxonomy-analysis-table/) to get the SRA accessions for all runs containing Chikungunya virus reads:

```sql
SELECT * FROM `nih-sra-datastore.sra_tax_analysis_tool.tax_analysis` WHERE name = 'Chikungunya virus' and total_count > 1
```

The [resulting table](data/SRAQueryResults.csv) is in `data/` directory. We fetched metadata associated with each SRA accession using [NCBI Entrez](https://biopython.org/docs/1.76/api/Bio.Entrez.html) implemented in BioPython. The resulting metadata for each run can be found [here](data/SRA_Runs.csv).

## Variant Detection Pipeline

We used a `Snakemake` pipeline to download publicly deposited sequencing runs from the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) and detect viral variants. Of the 9606 sequencing runs with detectable CHIKV sequences (as of 06/05/25), we downloaded runs with at least 10,000 CHIKV kmer hits——roughly 100X coverage over the ~11Kb CHIKV genome——resulting in 1976 total sequencing runs.We used `parallel-fastq-dump` to download and split read pairs into separate files. Runs where treated differently depending on whether they were from paired-end or single-end layouts. We used `fastp` to trim adaptors, low quality bases, and repetitive stretches and filtered the CHIKV reads from each run using `BBduk` with a kmer length of 31 and a hamming distance of 1. Filtered and trimmed reads were aligned to the sequence of the DMS library strain ([MW473668](https://www.ncbi.nlm.nih.gov/nuccore/MW473668)) using `minimap2` with the preset chosen according to the sequencing technology (Illumina, Oxford Nanopore, or PacBio). We identified viral variants relative to the library sequence using `iVar` with a minimum frequency of 0.05 and a minimum base quality score of 20. The viral variants were further filtered to sites covered by greater than 100 reads and `iVar` p-values less than 0.05. 

## Results

We split the final set of filtered and annotated variants into three main files (all in [`results/summary`](results/summary)):
- [all_variants.csv](results/summary/all_variants.csv): All of the coding variants in the CHIKV E protein. 
- [fixed_variants.csv](results/summary/fixed_variants.csv): Only variants that are 'fixed' in that sample (frequency > 1 - minimum frequency).
- [minor_variants.csv](results/summary/minor_variants.csv): Only minor variants (the minor codon at a site whose frequency is less than 0.5).