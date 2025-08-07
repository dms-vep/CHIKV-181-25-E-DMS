---
aside: false
---

# Mutation effects on cell entry
This page shows via interactive plots the effects of mutations to the CHIKV envelope proteins on entry into 293T-Mxra8, 293T-TIM1, and C6/36 cells.
Specifically:

[[toc]]

## Heatmap of effects on entry in 293T-Mxra8 cells
In the interactive heatmap, more negative values (red) indicate the mutation impairs entry, values of zero (white) indicate a mutation has no effect on entry, and positive values (blue) indicates a mutation improves entry.
The `x` in each column indicates the amino-acid identity in the parental 181/25 CHIKV strain.
Light gray squares indicate the effect of mutation was not reliably measured in our experiments.
You can use the zoom bar at top to zoom into specific regions, and you can mouseover points for details about specific mutations.
The lineplot at top shows the mean effect of mutations at each site.
Below the plot are various optons to filter data or otherwise adjust the plot.

Click [here](htmls/293T-Mxra8_entry_func_effects.html){target="_self"} for a standalone version of this plot

<Figure caption="Effects of mutations on entry in 293T-Mxra8 cells">
    <Altair :showShadow="true" :spec-url="'htmls/293T-Mxra8_entry_func_effects.html'"></Altair>
</Figure>

## Heatmap of effects on entry in 293T-TIM1 cells
This interactive heatmap is like the [one explained above](/cell_entry.html#heatmap-of-effects-on-entry-in-293t-mxra8-cells) except it shows effects of mutations in 293T-TIM1 cells.

Click [here](htmls/293T-TIM1_entry_func_effects.html){target="_self"} for a standalone version of this plot

<Figure caption="Effects of mutations on entry in 293T-TIM1 cells">
    <Altair :showShadow="true" :spec-url="'htmls/293T-TIM1_entry_func_effects.html'"></Altair>
</Figure>

## Heatmap of effects on entry in C6/36 cells
This interactive heatmap is like the [one explained above](/cell_entry.html#heatmap-of-effects-on-entry-in-293t-mxra8-cells) except it shows effects of mutations in 293T-TIM1 cells.

Click [here](htmls/C636_entry_func_effects.html){target="_self"} for a standalone version of this plot

<Figure caption="Effects of mutations on entry in C6/36 cells">
    <Altair :showShadow="true" :spec-url="'htmls/C636_entry_func_effects.html'"></Altair>
</Figure>

## Single large heatmap of effects on entry in all three cells
The interactive heatmap below shows the measured effects of mutations on entry in all three cells in a single plot, alongside the effects of mutations on [binding to mouse Mxra8](/mxra8_binding).
Note that to compare effects across cells we recommend instead using the interactive plot that shows the [differences among cells](/cell_entry_diffs).

Click [here](htmls/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding_overlaid.html){target="_self"} for a standalone version of the below plot.

<Figure caption="Effects of mutations on entry in all three cells">
    <Altair :showShadow="true" :spec-url="'htmls/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding_overlaid.html'"></Altair>
</Figure>

## Effects of mutations projected on protein structure
**Still need to add this section**

## Numerical data in CSV format
For the numerical values of the mutation effects on cell entry plotted above, see the following files:

  - [Merged CSV of effects of mutations on entry in all three cells](https://github.com/dms-vep/CHIKV-181-25-E-DMS/blob/main/results/summaries/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding.csv): the measurements in this CSV are filtered to only high-confidence measurements. Unless you understand the the QC filtering in detail, we recommend you use this CSV.
  - Effects of mutations on entry in each cell **without** the filtering for high-confidence measurements applied:
    - [non-filtered effects for 293T-Mxra8 cells](https://github.com/dms-vep/CHIKV-181-25-E-DMS/blob/main/results/func_effects/averages/293T-Mxra8_entry_func_effects.csv)
    - [non-filtered effects for 293T-TIM1 cells](https://github.com/dms-vep/CHIKV-181-25-E-DMS/blob/main/results/func_effects/averages/293T-TIM1_entry_func_effects.csv)
    - [non-filtered effects for C6/36 cells](https://github.com/dms-vep/CHIKV-181-25-E-DMS/blob/main/results/func_effects/averages/C636_entry_func_effects.csv)
