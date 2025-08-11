---
aside: false
---

# Compare how mutations affect entry in different cells
This page has interactive plots that make it easy to compare how mutations to the CHIKV envelope proteins affect entry in different cells.
The cells are 293T-Mxra8 cells (human cells expressing the Mxra8 receptor), 293T-TIM1 cells (human cells expressing the TIM1 attachment factor), and C6/36 cells (mosquito cells).
Specifically:

[[toc]]

## Differences in effects across cells
The interactive plot below provides the best way to compare mutation effects on entry across the cells.
The line plot shows the difference in effects of mutations at each site between a *comparator* and *reference* cell line (use the dropdowns below the plot to choose which cell line is the *comparator* and *reference*).
If you click on sites on that line plot, then the heatmaps below the line plot will zoom to show that cite in the center, and you can mouseover the heatmaps to see the effect of the mutation in each cell.
The scatter plot to the right of the heatmaps shows the effect of each mutation in each of the two cells being compared.
In the heatmap, gray indicates a mutation that was not measured and `x` indicates the amino acid in the parental 181/25 strain.
You can zoom in using the zoom bar below the line plot.
The options below the chart also allow you to show difference metrics for how the difference between mutation effects is calculated for the line plot, as well as the floor assigned to the mutation effects in the plots.

(Click [here](htmls/compare_cell_entry_site_zoom.html){target="_self"} for a standalone version of the plot.)

<Figure caption="Zoomable plot of differences in mutation effects on entry across cells">
    <Altair :showShadow="true" :spec-url="'htmls/compare_cell_entry_site_zoom.html'"></Altair>
</Figure>

## Scatter plot of mutation effects in different cells
The scatter plots below show the effects of each mutation in each pair of cells.
You can mouseover the points for details about each mutation, and click on the colors in the key to the right of the plot to only show mutations in different proteins (E3, E2, 6K, E1).

(Click [here](htmls/compare_cell_entry_scatter.html){target="_self"} for a standalone version of the plot.)

<Figure caption="Scatter plot of mutation effects on entry in different cells">
    <Altair :showShadow="true" :spec-url="'htmls/compare_cell_entry_scatter.html'"></Altair>
</Figure>

## Numerical data in CSV format
For the numerical values of the differences in mutation effects on entry across cells, see the following files:
  - [Difference in effects of mutations at each site between each pair of cells](https://github.com/dms-vep/CHIKV-181-25-E-DMS/blob/main/results/compare_cell_entry/site_diffs.csv): CSV reporting site-level differences using three different metrics to quantify the differences.
  - [Differences in effects of each mutation between each pair of cells](https://github.com/dms-vep/CHIKV-181-25-E-DMS/blob/main/results/compare_cell_entry/site_diffs.csv): CSV with effect of each mutation on entry in each cell as well as the differences in these effects across cells. Note that the differences are computed after flooring the per-cell effects, since differences in effects between very negative (highly deleterious) mutations may not be meaningful.
