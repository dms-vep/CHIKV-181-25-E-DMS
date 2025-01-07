# Analysis Notebooks

This directory contains custom analysis notebooks that will may be added to the pipeline at a later date.

## Notebooks

- [`compare-cell-entry`](/notebooks/compare-cell-entry.ipynb): Explores the differences between functional selections (cell entry) in three different cell lines; one expressing the known receptor **Mxra8**, one expressing the entry factor **TIM-1**, and **C636** mosquito cells expression an unknown receptor.

- [`plot-data-on-structures`](/notebooks/plot-data-on-structures.ipynb): Display the functional selections and difference between cell types on the CHIKV E + Mxra8 crystal structure with `dms-viz`.

## `dms-viz` visualizations

Stop codons and InDels are filtered out. Mutations with `times_seen < 2` are filtered out as well.

- [Functional selections on CHIKV E monomer](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV_181-25_E_DMS%2Frefs%2Fheads%2Fplot-individual-mutation-effects%2Fnotebooks%2Fdms-viz%2Foutput%2FCHIKV_6JO8_functional_scores.json%3Ftoken%3DGHSAT0AAAAAACQIDWSPVA72B5HCCJ63EBPSZ34ZEMQ&pe=MXRA8&fi=%257B%2522n_selections%2522%253A1%252C%2522times_seen%2522%253A2%257D&ce=%255B%2522MXRA8%2522%255D&sa=true)
- [Functional selections on CHIKV E timer](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV_181-25_E_DMS%2Frefs%2Fheads%2Fplot-individual-mutation-effects%2Fnotebooks%2Fdms-viz%2Foutput%2FCHIKV_6JO8_trimer_functional_scores.json%3Ftoken%3DGHSAT0AAAAAACQIDWSOR6JP6XVEFNM23WQCZ34ZE4Q&pe=MXRA8&fi=%257B%2522n_selections%2522%253A1%252C%2522times_seen%2522%253A2%257D&sa=true&ce=%255B%2522MXRA8%2522%255D)
- [Functional differences on CHIKV E monomer](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV_181-25_E_DMS%2Frefs%2Fheads%2Fplot-individual-mutation-effects%2Fnotebooks%2Fdms-viz%2Foutput%2FCHIKV_6JO8_functional_diff.json%3Ftoken%3DGHSAT0AAAAAACQIDWSOF3MJLHYDJA3HSNWQZ34ZDUA&pe=C636_v_MXRA8&ce=%255B%2522C636_v_MXRA8%2522%255D&sa=true)
- [Functional differences on CHIKV E trimer](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV_181-25_E_DMS%2Frefs%2Fheads%2Fplot-individual-mutation-effects%2Fnotebooks%2Fdms-viz%2Foutput%2FCHIKV_6JO8_trimer_functional_diff.json%3Ftoken%3DGHSAT0AAAAAACQIDWSO7E2XDLUOICIZD466Z34ZEXA&pe=C636_v_MXRA8&ce=%255B%2522C636_v_MXRA8%2522%255D&sa=true)