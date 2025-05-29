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

## Identifying CHIKV Sequencing Runs

The NCBI tabulates the taxonomic breakdown of each sequencing run in the SRA with the [SRA Taxonomy Tool](https://www.ncbi.nlm.nih.gov/sra/docs/sra-taxonomy-analysis-tool/). This information is [stored on the cloud](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-taxonomy-analysis-table/) and can be queried with Google's BigQuery UI. We used this information to compile a list of SRA accessions for all sequencing runs that contain a minimum number of reads from Chikungunya virus. 

To do this, we used Google's BigQuery UI by logging into the [Google Cloud Console](https://console.cloud.google.com/welcome) and selecting `Run a query in BigQuery`. This should take you to something called 'BigQuery Studio' where you'll be prompted to log in and create a new project. Once you've created a project, select `SQL Query` under 'Create new', which will open an empty text file in the browser window. Enter the following SQL command per [these instructions from the NCBI](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-taxonomy-analysis-table/) to get the SRA accessions for all runs containing Chikungunya virus reads:

```sql
SELECT * FROM `nih-sra-datastore.sra_tax_analysis_tool.tax_analysis` WHERE name = 'Chikungunya virus' and total_count > 1
```

The [resulting table](data/SRAQueryResults.csv) is in `data/` directory.