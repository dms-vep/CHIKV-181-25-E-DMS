# Visualize CHIKV Functional Effects with `dms-viz`

Here's some quick code to visualize the CHIKV functional effects in various conditions (TIM-1, Mxra8, and C636 cells).

## Configure `dms-viz`

The python package `configure-dms-viz` provides a command that takes in your data and outputs a visualization file that you can upload to [`dms-viz`](https://dms-viz.github.io/v0/).

To use `configure-dms-viz`, make a conda environment:

```bash
conda env create -f data/dms-viz/environment.yml
conda activate dms-viz
```

### All Functional Data

This command makes a visualization with the average functional effects of mutations in Mxra8, TIM-1, and C636 cells. I combined all three datasets in this notebook: [data/dms-viz/combine-datasets.ipynb](data/dms-viz/combine-datasets.ipynb). I'm coloring only chains 'A' and 'B' and excluding the other chains in the model except for one Mxra8 chain ('O').

After activating the conda environment, copy and paste this to the command line:

```bash
configure-dms-viz format \
  --input data/dms-viz/input/all_func_effects.csv \
  --sitemap data/dms-viz/sitemap/CHIKV_6JO8_sitemap.csv \
  --output data/dms-viz/output/CHIKV_6JO8_functional_scores.json \
  --name "CHIKV Func." \
  --metric "effect" \
  --metric-name "Functional Effect" \
  --exclude-amino-acids "*, -" \
  --included-chains "A B" \
  --excluded-chains "C E D F M N" \
  --condition "condition" \
  --condition-name "Cell Line" \
  --tooltip-cols "{'times_seen': '# Obsv'}" \
  --filter-cols "{'times_seen': 'Times Seen'}" \
  --filter-limits "{'times_seen': [0, 2, 10]}" \
  --structure "6JO8"
```