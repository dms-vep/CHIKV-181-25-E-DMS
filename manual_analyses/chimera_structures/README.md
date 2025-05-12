# Generate file for ChimeraX visualization of functional effects on structure

Analysis by Xiaohui Ju, helped by Brendan Larsen. This notebook is to generate `.defattr` file that can be read by ChimeraX for visualization of functional effects on structure.

## Organization

The contents of `chimera_structures/` are organized as follows:

```bash
├── prep_293T_Mxra8_entry_data_for_chimera # <---- Jupyter notebook
├── README.md
├── results # <------------ input file, output file for ChimeraX
```

### Data

The notebook reads file ([`entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding_annotated_site_means.csv`](../../results/annotated_summary_csvs/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding_annotated_site_means.csv)) and extract desired columns to generate ([`293T-Mxra8_entry_func_effects_mean.csv`](../../manual_analyses/chimera_structures/results/293T-Mxra8_entry_func_effects_mean.csv)), which is the input file to create ([`293T-Mxra8_entry_func_effects.defattr`](../../manual_analyses/chimera_structures/results/293T-Mxra8_entry_func_effects.defattr)) for ChimeraX visualization of functional effects in 293T-Mxra8 cells on structure (PDB: `6nk7`). Similarly, `.defattr` file for visualization of functional effects in 293T-TIM1 fand C636 cells on structure can be generated easily if needed.
