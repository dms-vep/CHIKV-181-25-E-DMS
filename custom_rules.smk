"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""


rule compare_cell_entry:
    """Compare entry across cell lines."""
    input:
        nb="notebooks/compare_cell_entry.ipynb",
        mut_effects_csv="results/summaries/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding.csv",
    output:
        nb="results/notebooks/compare_cell_entry.ipynb",
        site_diffs_csv="results/compare_cell_entry/site_diffs.csv",
        mut_scatter_chart="results/compare_cell_entry/compare_cell_entry_scatter.html",
        site_zoom_chart="results/compare_cell_entry/compare_cell_entry_site_zoom.html",
    params:
        params_yaml=lambda wc, input: yaml_str(
            {
                # cells and their names in input CSV file
                "cells": {
                    "293T-Mxra8": "293T_Mxra8", "C6/36": "C636", "293T-TIM1": "293T_TIM1",
                },
                # for calculating differences and display, floor mutation effects at this
                "floor_mut_effects": -5,
            }
        ),
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/compare_cell_entry.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p mut_effects_csv {input.mut_effects_csv} \
            -p site_diffs_csv {output.site_diffs_csv} \
            -p mut_scatter_chart {output.mut_scatter_chart} \
            -p site_zoom_chart {output.site_zoom_chart} \
            -y "{params.params_yaml}" \
            &> {log}
        """

docs["Compare entry among cells"] = {
    "Final plots": {
        "Scatter plots of mutation effects in different cells": rules.compare_cell_entry.output.mut_scatter_chart,
        "Zoomable site chart of differences": rules.compare_cell_entry.output.site_zoom_chart,
    },
    "Analysis notebooks": {
        "Notebook making plots comparing entry": rules.compare_cell_entry.output.nb,
    },
    "Data files": {
        "site-differences in entry effects": rules.compare_cell_entry.output.site_diffs_csv,
    }
}

# Files (Jupyter notebooks, HTML plots, or CSVs) that you want included in
# the HTML docs should be added to the nested dict `docs` or `other_target_files`
