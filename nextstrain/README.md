# CHIKV E Nextstrain Build

A `Snakemake` pipeline that builds an interactive [Nextstrain](https://nextstrain.org/) tree of CHIKV sequences (just the E protein). Sequences are annotated by collection date, geographic location, and host. The interactive tree can be visualized using [Auspice](https://docs.nextstrain.org/projects/auspice/en/stable/).

Analysis by Will Hannon, adapted from Caleb Carr's [RABV-G pipeline](https://github.com/dms-vep/RABV_Pasteur_G_DMS/tree/main/non-pipeline_analyses/RABV_nextstrain).

## Accessions

To download the accessions, go to [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and click *Search by virus*. In the *Search by virus name or taxonomy* box, enter **Chikungunya virus, taxid:37124** and hit enter. Then click the  *Download* option, select *Accession List* and *Nucleotide* options and hit *Next*. On the next page, select *Download All Records* and hit *Next*. On the next page, select *Accession with version* and click *Download*. Sequences are downloaded from the list of accessions because more information is extracted from the genbank file during the download process. The [current accession list](configuration/sequences.acc) was downloaded on **April 4th, 2025**.


## Rooting

The tree is rooted with [O'nyong'nyong virus (ONNV)](https://www.ncbi.nlm.nih.gov/nuccore/NC_075006.1).

## Notes

There's a reasonably up-to-date tree of the entire CHIKV genome here: https://github.com/ViennaRNA/CHIKV. I can probably grab the accessions from this tree and re-build it with the E protein only. This could be a fast alternative to get up and running before building a more comprehensive tree.