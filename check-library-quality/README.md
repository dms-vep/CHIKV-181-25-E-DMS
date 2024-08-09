# Check Library Quality

*This is a miniature version of `dms-vep-pipeline-3` modified to deal with PacBio samples lacking barcodes. Essentially, this pipeline includes every rule in [build_variants.smk](../dms-vep-pipeline-3/build_variants.smk).*

## Overview

The goal of this analysis is to assess the quality of the CHIKV E protein libraries from GenScript. There are two separate libraries covering different segments of the E protein:

- **E3E2**: 2:488
- **6KE1**: 489:988

For each library, we have PacBio sequencing from three stages of the experiment: 

1. `Genscript_Product`: The raw pool of variants from GenScript
2. `Barcode_PCR`: Variants after adding barcodes with PCR
3. `Plasmid_Pool`: Barcoded variants digested from the plasmid pool after liquid culture

There are no barcodes in the products from GenScript. Also, although there are barcodes in the Library Pool, the plasmids were digested too close to the barcode and the barcodes are soft-clipped during alignment.

## Analysis

To run the pipeline, navigate to this directory:

```
cd check-library-quality
```

Activate the `dms-vep-pipeline-3` conda environment:

```
conda activate dms-vep-pipeline-3
```

And run the snakemake pipeline:

```
snakemake --cores 6 -j 6
```

