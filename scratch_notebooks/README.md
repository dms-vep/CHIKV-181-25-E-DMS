# Scratch Notebooks

This directory contains custom analysis notebooks that will be added to the pipeline at a later date.

## Notebooks

- [average-functional-effects](average-functional-effects.ipynb): Combines functional effects of mutations in 181/25 E in three different cell lines; 293T cells expressing the known receptor **Mxra8**, 293T cells expressing the entry factor **TIM1**, and **C636** mosquito cells. 

- [`compare-cell-entry`](compare-cell-entry.ipynb): Explores the differences between functional selections (cell entry) in three different cell lines; one expressing the known receptor **Mxra8**, one expressing the entry factor **TIM-1**, and **C636** mosquito cells expression an unknown receptor.

- [`plot-data-on-structures`](plot-data-on-structures.ipynb): Display the functional selections and difference between cell types on the CHIKV E + Mxra8 crystal structure with `dms-viz`.

## `dms-viz` visualizations

Stop codons and InDels are filtered out. Mutations with `times_seen <= 2` are filtered out as well. We might want to consider different filters.

There are two papers that elucidate the structure of CHIKV E in complex with Mxra8——[this paper](https://www.sciencedirect.com/science/article/pii/S0092867419303940?via%3Dihub) and [this paper](https://www.cell.com/cell/pdf/S0092-8674(19)30392-7.pdf).

For several reasons, I think the structure from [Micheal Diamond and Daved Fremont's paper](https://www.cell.com/cell/pdf/S0092-8674(19)30392-7.pdf) is better suited for viewing our data. They used a combination of X-ray crystallography, Cryo-EM of viral particles, and computational reconstruction to get structures of 'mature' CHIKV E  in complex with mouse Mxra8 on [VLPs](https://www.rcsb.org/structure/6NK6) and [infectious particles](https://www.rcsb.org/structure/6NK7). It’s easy to orient yourself visually with these structures because they contain the transmembrane domain and capsid; The models reflect the T=4 symmetry of CHIKV E, making it easy to show all 4 binding ‘sites’ with Mxra8; There are structures with (from infectious particles) and without (from the VLP) E3 retention.

I plotted the functional scores and the difference in functional scores between cell types on the [VLP](https://www.rcsb.org/structure/6NK6) and [infectious particle](https://www.rcsb.org/structure/6NK7) structures.

To illustrate the three types of contacts that Mxra8 makes with CHIKV E, I made `dms-viz` visualizations showing only the wrapped, interspike, and intraspike interactions respectively:

- Wrapped [Functional Scores](dms-viz/output/CHIKV_VLP_wrapped_monomer_functional_scores.json) and [Differences](dms-viz/output/CHIKV_VLP_wrapped_monomer_functional_differences.json)
- Intraspike [Functional Scores](dms-viz/output/CHIKV_VLP_intraspike_monomer_functional_scores.json) and [Differences](dms-viz/output/CHIKV_VLP_intraspike_monomer_functional_differences.json)
- Interspike [Functional Scores](dms-viz/output/CHIKV_VLP_interspike_monomer_functional_scores.json) and [Differences](dms-viz/output/CHIKV_VLP_interspike_monomer_functional_differences.json)

Additionally, since E3 retention changes the number of sites that can bind to Mxra8 with high affinity, I plotted the data on a structure without E3 (the VLP) and with E3 (the infectious particle).

- Without E3 [Functional Scores](dms-viz/output/CHIKV_VLP_full_functional_scores.json) and [Differences](dms-viz/output/CHIKV_VLP_full_functional_differences.json)
- With E3 [Functional Scores](dms-viz/output/CHIKV_infectious_full_functional_scores.json) and [Differences](     dms-viz/output/CHIKV_infectious_full_functional_differences.json)
