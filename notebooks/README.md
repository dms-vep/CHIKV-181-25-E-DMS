# Analysis Notebooks

This directory contains custom analysis notebooks that will may be added to the pipeline at a later date.

## Notebooks

- [`compare-cell-entry`](/notebooks/compare-cell-entry.ipynb): Explores the differences between functional selections (cell entry) in three different cell lines; one expressing the known receptor **Mxra8**, one expressing the entry factor **TIM-1**, and **C636** mosquito cells expression an unknown receptor.

- [`plot-data-on-structures`](/notebooks/plot-data-on-structures.ipynb): Display the functional selections and difference between cell types on the CHIKV E + Mxra8 crystal structure with `dms-viz`.

## `dms-viz` visualizations

Stop codons and InDels are filtered out. Mutations with `times_seen < 2` are filtered out as well.

- [Functional scores on  CHIKV E 6JO8 monomer](./dms-viz/output/CHIKV_6JO8_functional_scores.json)
- [Functional scores on CHIKV E 6JO8 trimer](./dms-viz/output/CHIKV_6JO8_trimer_functional_scores.json)
- [Functional scores on CHIKV E 6NK6 trimer](./dms-viz/output/CHIKV_6NK6_functional_scores.json)
- [Functional difference on CHIKV E 6JO8 monomer](./dms-viz/output/CHIKV_6JO8_functional_diff.json)
- [Functional difference on CHIKV E 6JO8 trimer](./dms-viz/output/CHIKV_6JO8_trimer_functional_diff.json)
- [Functional difference on CHIKV E 6NK6 trimer](./dms-viz/output/CHIKV_6NK6_functional_diff.json)

