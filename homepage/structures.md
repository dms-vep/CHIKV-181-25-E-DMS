---
aside: false
---

# Visualizations of deep mutational scanning on envelope protein structures
This provides access to the effects of mutations mapped onto three-dimensional structures of the envelope proteins, as rendered by [dms-viz](https://dms-viz.github.io/dms-viz-docs/).

It uses two structures, both from [this paper](https://pubmed.ncbi.nlm.nih.gov/31080061/), that are cryo-EM structures of the CHIKV envelope proteins bound to mouse Mxra8.
In both structures, there are four envelope heteromers (E2 / E1 or E3 / E2 / E1) in the asymmetric unit even though envelop actually forms a trimer. Note that [as described in the paper](https://pubmed.ncbi.nlm.nih.gov/31080061/), Mxra8 contacts E2 / E1 in three different modes: wrapped, intraspike, and interspike. The structures are:
  - [PDB ID 6nk7](https://www.rcsb.org/structure/6NK7): E2 / E1 on CHIKV virus-like-particles bound to four copies of Mxra8
  - [PDB ID 6nk6](https://www.rcsb.org/structure/6NK6): E3 / E2 / E1 on CHIKV virions bound to one copy of Mxra8.

Note that all of the [dms-viz](https://dms-viz.github.io/dms-viz-docs/) visualizations are interactive, and allow you to click on sites in the line plot to show them on the envelope proteins structure, as well as various other interactive options.
In all of the structures below, the envelope proteins are initialy gray while Mxra8 is shown in yellow; clicking on a site in the line plot will color it on the envelope protein structure based on its measured effect:
  
[[toc]]

## Effects of mutations on cell entry
Below is a visualization of the effects of mutations on cell entry on the trimer in `6nk7`.
Click on the *Sidebar* and *Chart Options* to select which cell line to show entry for; click on the line plot to highlight sites on the structure.
To view it in a standalone page (which may be easier to interact with), click [here](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_on_6nk7_trimer%2Fcell_entry_on_6nk7_trimer.json&ce=%255B%2522entry%2520in%2520293T_Mxra8%2520cells%2522%255D&bc=%23d2cd32).

<iframe src="https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_on_6nk7_trimer%2Fcell_entry_on_6nk7_trimer.json&ce=%255B%2522entry%2520in%2520293T_Mxra8%2520cells%2522%255D&bc=%23d2cd32" width="100%" height="700px"></iframe>

For additional visualizations, see:
  - [Cell entry on the `6nk7` asymmetric unit](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_on_6nk7_asymmetric_unit%2Fcell_entry_on_6nk7_asymmetric_unit.json&ce=%255B%2522entry%2520in%2520293T_Mxra8%2520cells%2522%255D&bc=%23c5bf2b)
  - [Cell entry on the `6nk6` asymmetric unit](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_on_6nk6_asymmetric_unit%2Fcell_entry_on_6nk6_asymmetric_unit.json&ce=%255B%2522entry%2520in%2520293T_Mxra8%2520cells%2522%255D&bc=%23c5bf2b)
  - [Cell entry on the `6nk6` structure showing wrapped contact mode](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_on_6nk6_wrapped%2Fcell_entry_on_6nk6_wrapped.json&bc=%23d1c32e&ce=%255B%2522entry%2520in%2520293T_Mxra8%2520cells%2522%255D)
  - [Cell entry on the `6nk6` structure showing the intraspike contact mode](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_on_6nk6_intraspike%2Fcell_entry_on_6nk6_intraspike.json&bc=%23c8d331&ce=%255B%2522entry%2520in%2520293T_Mxra8%2520cells%2522%255D)
  - [Cell entry on the `6nk6` structure showing the interspike contact mode](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_on_6nk6_interspike%2Fcell_entry_on_6nk6_interspike.json&ce=%255B%2522entry%2520in%2520293T_Mxra8%2520cells%2522%255D&bc=%23cdc82d)

## Differences in effects of mutations on entry in different cells
Below is a visualization of the **differences** in effects of mutations on entry in two different cell lines on the trimer in `6nk7`.
Click on the *Sidebar* and *Chart Options* to select which cell pair to show the differences for; click on the line plot to highlight sites on the structure.
To view it in a standalone page (which may be easier to interact with), click [here](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_diffs_on_6nk7_trimer%2Fcell_entry_diffs_on_6nk7_trimer.json&ce=%255B%2522293T_Mxra8%2520minus%2520293T_TIM1%2522%255D&bc=%23cbcd32).

<iframe src="https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_diffs_on_6nk7_trimer%2Fcell_entry_diffs_on_6nk7_trimer.json&ce=%255B%2522293T_Mxra8%2520minus%2520293T_TIM1%2522%255D&bc=%23cbcd32" width="100%" height="700px"></iframe>

For additional visualizations, see:
  - [Cell entry differences on `6nk7` asymmetric unit](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_diffs_on_6nk7_asymmetric_unit%2Fcell_entry_diffs_on_6nk7_asymmetric_unit.json&ce=%255B%2522293T_Mxra8%2520minus%2520293T_TIM1%2522%255D&bc=%23c3c534)
  - [Cell entry differences on `6nk6` asymmetric unit](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fcell_entry_diffs_on_6nk6_asymmetric_unit%2Fcell_entry_diffs_on_6nk6_asymmetric_unit.json&ce=%255B%2522293T_Mxra8%2520minus%2520293T_TIM1%2522%255D&bc=%23c3c534)

## Effects of mutations on Mxra8 binding
Below is a visualization of the effects of mutations on Mxra8 binding on the trimer in `6nk7`.
Click on the *Sidebar* and *Chart Options* tho select whether to show binding to mouse Mxra8 or human Mxra8.
To view it on a standalone page (which may be easier to interact with), click [here](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fmxra8_binding_on_6nk7_trimer%2Fmxra8_binding_on_6nk7_trimer.json&pe=binding+to+mouse+Mxra8&fi=%257B%2522entry_in_293T_Mxra8_cells%2522%253A-4%257D&ce=%255B%2522binding%2520to%2520mouse%2520Mxra8%2522%255D&bc=%23c1c33c).

<iframe src="https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fmxra8_binding_on_6nk7_trimer%2Fmxra8_binding_on_6nk7_trimer.json&pe=binding+to+mouse+Mxra8&fi=%257B%2522entry_in_293T_Mxra8_cells%2522%253A-4%257D&ce=%255B%2522binding%2520to%2520mouse%2520Mxra8%2522%255D&bc=%23c1c33c" width="100%" height="700px"></iframe>

For additional visualizations, see:

  - [Mxra8 binding on `6nk7` asymmetric unit](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fmxra8_binding_on_6nk7_asymmetric_unit%2Fmxra8_binding_on_6nk7_asymmetric_unit.json&pe=binding+to+mouse+Mxra8&fi=%257B%2522entry_in_293T_Mxra8_cells%2522%253A-4%257D&ce=%255B%2522binding%2520to%2520mouse%2520Mxra8%2522%255D&bc=%23c1c33c)
  - [Mxra8 binding on `6nk6` asymmetric unit](https://dms-viz.github.io/v0/?data=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-vep%2FCHIKV-181-25-E-DMS%2Frefs%2Fheads%2Fmain%2Fresults%2Fdms-viz%2Fmxra8_binding_on_6nk6_asymmetric_unit%2Fmxra8_binding_on_6nk6_asymmetric_unit.json&pe=binding+to+mouse+Mxra8&fi=%257B%2522entry_in_293T_Mxra8_cells%2522%253A-4%257D&ce=%255B%2522binding%2520to%2520mouse%2520Mxra8%2522%255D&bc=%23c1c33c)
