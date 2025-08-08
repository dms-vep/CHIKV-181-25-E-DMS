---
aside: false
---

# Mutation effects on Mxra8 binding
This page has interactive plots showing the effects of mutations to the CHIKV envelope proteins on Mxra8 binding, as assessed by neutralization of pseudovirus infection of 293T-Mxra8 cells by soluble mouse or human Mxra8.
Note that these measurements may be a bit more noisy than the cell entry measurements, and that the measurements for human Mxra8 and especially noisy and may be less reliable than those for mouse Mxra8 (see the analysis notebooks for [mouse](notebooks/avg_escape_receptor_affinity_mouse_Mxra8.html){target="_self"} and [human](notebooks/avg_escape_receptor_affinity_human_Mxra8.html){target="_self"} Mxra8 binding, respectively, to see the greater noise for the human Mxra8 binding).
Due to this noise, it can be helpful to assess individual measurements by mousing over them in the interactive heatmaps for [mouse Mxra8](mxra8_binding.html#heatmap-of-binding-to-mouse-mxra8){target="_self"} or [human Mxra8](mxra8_binding.html#heatmap-of-binding-to-human-mxra8){target="_self"} and looking at the per-library measurements to see if they are similar.

You can access the data as follows:
[[toc]]

## Heatmap showing binding to mouse and human Mxra8 alongside entry in 293T-Mxra8 cells
The heatmaps below show entry in 293T-Mxra8 cells (293T cells expressing human Mxra8) alongside the effects of mutations on binding to mouse and human Mxra8.
You can use the zoom bar at top to zoom to specific site ranges, and mouseover the heatmap for details on specific mutations.
In the heatmap, negative values indicate worse entry or binding, and the `x` indicates the parental amino-acid in the CHIKV 181/25 envelope proteins.
Light gray squares indicate mutations that were not reliably measured in the deep mutational scanning.
Dark gray squares for the binding plots indicate mutations that are too deleterious for cell entry to reliably measure their effect on Mxra8 binding.
The slider at the bottom allows you to set this threshold for what is considered "too deleterious for cell entry" to report a binding effect.

(Click [here](htmls/binding_mouse_vs_human_Mxra8_overlaid.html){target="_self"} for a standalone version of the plot.)

<Figure caption="Interactive heatmap showing entry in 293T cells expressing human Mxra8 alongside binding to mouse and human Mxra8">
    <Altair :showShadow="true" :spec-url="'htmls/binding_mouse_vs_human_Mxra8_overlaid.html'"></Altair>
</Figure>

## Heatmap of binding to mouse Mxra8
Effects of mutations on binding to mouse Mxra8.
The line plot shows the total effect mutations at a site have on binding (use options below heatmap to select site summary statistic), and you can use the scale bar to zoom on specific sites.
Light gray squares indicate mutations that were not reliably measured in the deep mutational scanning.
Dark gray squares for the binding plots indicate mutations that are too deleterious for cell entry to reliably measure their effect on Mxra8 binding.

When you mouseover points you get the measurement in each of the two library replicates, which can be helpful for assessing confidence in measurements.

(Click [here](htmls/mouse_Mxra8_mut_effect.html){target="_self"} for a standalone version of the plot.)

<Figure caption="Interactive heatmap of effects of mutations on binding to mouse Mxra8">
    <Altair :showShadow="true" :spec-url="'htmls/mouse_Mxra8_mut_effect.html'"></Altair>
</Figure>

## Heatmap of binding to human Mxra8
This plot is like the one shown immediately above but for the effects of mutations on binding to human rather than mouse Mxra8.

(Click [here](htmls/human_Mxra8_mut_effect.html){target="_self"} for a standalone version of the plot.)

<Figure caption="Interactive heatmap of effects of mutations on binding to human Mxra8">
    <Altair :showShadow="true" :spec-url="'htmls/human_Mxra8_mut_effect.html'"></Altair>
</Figure>

## Effects of mutations at each site on binding
This plot shows the total effects of mutations at each site on binding to human and mouse Mxra8.
Mouse over lines for details.

(Click [here](htmls/mxra8_site_chart.html){target="_self"} for a standalone version of the plot.)

<Figure caption="Total effects of mutations at each site on Mxra8 binding">
    <Altair :showShadow="true" :spec-url="'htmls/mxra8_site_chart.html'"></Altair>
</Figure>

## Scatter plots of mutation effects on binding to mouse versus human Mxra8

This plot shows the measured effect of each mutation on binding to mouse versus human Mxra8 (recall there is some noise in these data):

<Figure caption="Effects of mutations on binding to human versus mouse Mxra8">
    <Altair :showShadow="true" :spec-url="'htmls/mxra8_mut_binding_corr.html'"></Altair>
</Figure>

This plot shows the total effects of mutation at each site on binding to mouse versus human Mxra8:

<Figure caption="Effects of all mutations at each site on binding to human versus mouse Mxra8">
    <Altair :showShadow="true" :spec-url="'htmls/mxra8_site_binding_corr.html'"></Altair>
</Figure>

## Numerical data in CSV format
The numerical data are available as follows:
  - [Effects of mutations on binding to human and mouse Mxra8](https://github.com/dms-vep/CHIKV-181-25-E-DMS/blob/main/results/summaries/binding_mouse_vs_human_Mxra8.csv): this CSV is filtered for just high-quality measurements, but note that you may also want to filter for mutations with reasonable cell entry as assessed by the *entry in 293T_Mxra8 cells* (the heatmaps above apply a filter of -4).
  - [Total effects of mutations at each site on binding to mouse or human Mxra8](https://github.com/dms-vep/CHIKV-181-25-E-DMS/blob/main/results/compare_human_mouse_mxra8/site_binding.csv)
