"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""

# Compare human and mouse Mxra8 binding -------------------------------------------------

rule compare_human_mouse_mxra8_binding:
    """Compare human versus mouse Mxra8 binding."""
    input:
        nb="notebooks/compare_human_mouse_mxra8_binding.ipynb",
        entry_293T_human_Mxra8="results/func_effects/averages/293T-Mxra8_entry_func_effects.csv",
        binding_human_Mxra8="results/receptor_affinity/averages/human_Mxra8_mut_effect.csv",
        binding_mouse_Mxra8="results/receptor_affinity/averages/mouse_Mxra8_mut_effect.csv",
        addtl_annotations="data/addtl_site_annotations.csv",
    output:
        nb="results/notebooks/compare_human_mouse_mxra8_binding.ipynb",
        corr_chart_html="results/compare_human_mouse_mxra8_binding.html",
    params:
        params_yaml=lambda wc: yaml_str(
            {
                "min_entry_293T_human_Mxra8": -3,
                "min_entry_293T_human_Mxra8_std": 2.25,
                "min_mouse_Mxra8_binding_std": 2.25,
                "min_human_Mxra8_binding_std": 2.5,
                "min_times_seen": 2,
            }
        ),
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    log:
        "results/logs/compare_human_mouse_mxra8_binding.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p entry_293T_human_Mxra8 {input.entry_293T_human_Mxra8} \
            -p binding_human_Mxra8 {input.binding_human_Mxra8} \
            -p binding_mouse_Mxra8 {input.binding_mouse_Mxra8} \
            -p addtl_annotations {input.addtl_annotations} \
            -p corr_chart_html {output.corr_chart_html} \
            -y "{params.params_yaml}" \
            &> {log}
        """

docs["Compare binding to human vs mouse Mxra8"] = {
    "scatter plot": rules.compare_human_mouse_mxra8_binding.output.corr_chart_html,
}


# Compare Mxra8 binding to entry --------------------------------------------------------

rule compare_mxra8_binding_to_entry:
    """Compare Mxra8 binding to cell entry."""
    input:
        data_csv="results/summaries/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding.csv",
        nb="notebooks/compare_mxra8_binding_to_entry.ipynb",
    output:
        corr_chart_html="results/compare_mxra8_binding_to_entry/compare_mxra8_binding_to_entry.html",
        nb="results/notebooks/compare_mxra8_binding_to_entry.ipynb",
    params:
        params_yaml=lambda wc: yaml_str(
            {
                "mut_effects_floor": -5,  # floor mut effects on entry at this value
                "min_293T_Mxra8_entry": -4,  # only consider binding for mutations w effects on 293T-Mxra8 entry >= this
                "cells": ["293T_Mxra8", "C636", "293T_TIM1"],
            }
        ),
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/compare_mxra8_binding_to_entry.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p data_csv {input.data_csv} \
            -p corr_chart_html {output.corr_chart_html} \
            -y "{params.params_yaml}" \
            &> {log}
        """

docs["Compare Mxra8 binding to cell entry"] = {
    "scatter plot": rules.compare_mxra8_binding_to_entry.output.corr_chart_html,
}

# Compare entry across cells ------------------------------------------------------------

rule cell_entry_mut_diffs:
    """Get CSV with differences in entry for mutations between cells after flooring."""
    input:
        mut_effects_csv="results/summaries/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding.csv",
        nb="notebooks/cell_entry_mut_diffs.ipynb",
    output:
        mut_diffs_csv="results/compare_cell_entry/mut_diffs.csv",
        nb="results/notebooks/cell_entry_mut_diffs.ipynb",
    params:
        params_yaml=lambda wc: yaml_str(
            {
                # cells names in input CSV file
                "cells": ["293T_Mxra8", "C636", "293T_TIM1"],
                # for calculating differences and display, floor mutation effects at this
                "floor_mut_effects": -5,
            }
        ),
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/cell_entry_mut_diffs.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p mut_effects_csv {input.mut_effects_csv} \
            -p mut_diffs_csv {output.mut_diffs_csv} \
            -y "{params.params_yaml}" \
            &> {log}
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
        params_yaml=lambda wc: yaml_str(
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
        "mutation-differences in entry effects (after flooring negative values)":
            rules.cell_entry_mut_diffs.output.mut_diffs_csv,
    }
}


# Configure dms-viz JSONs ---------------------------------------------------------------

# read configuration for `configure_dms_viz`
with open("data/dms_viz_config.yaml") as f:
    dms_viz_config = yaml.YAML(typ="safe", pure=True).load(f)

rule configure_dms_viz:
    """Configure a JSON for `dms-viz`."""
    input:
        data_csv=lambda wc: dms_viz_config[wc.viz_name]["data_csv"],
        sitemap_csv=lambda wc: dms_viz_config[wc.viz_name]["sitemap_csv"],
        nb="notebooks/configure_dms_viz.ipynb",
    output:
        dms_viz_json="results/dms-viz/{viz_name}/{viz_name}.json",
        pdb_file="results/dms-viz/{viz_name}/{viz_name}.pdb",
        input_data_csv="results/dms-viz/{viz_name}/{viz_name}_data.csv",
        input_sitemap_csv="results/dms-viz/{viz_name}/{viz_name}_sitemap.csv",
        nb="results/notebooks/configure_dms_viz_{viz_name}.ipynb",
    params:
        params_yaml=lambda wc: yaml_str(
            {
                key: dms_viz_config[wc.viz_name][key]
                for key in [
                    "pdb_id",
                    "pdb_type",
                    "name",
                    "melt_condition_metric_cols",
                    "metric",
                    "opt_params",
                ]
            }
        ),
    conda:
        "envs/dms-viz.yml"
    log:
        "results/logs/configure_dms_viz_{viz_name}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p data_csv {input.data_csv} \
            -p sitemap_csv {input.sitemap_csv} \
            -p dms_viz_json {output.dms_viz_json} \
            -p pdb_file {output.pdb_file} \
            -p input_data_csv {output.input_data_csv} \
            -p input_sitemap_csv {output.input_sitemap_csv} \
            -y "{params.params_yaml}" \
            &> {log}
        """

docs["dms-viz visualizations"] = {
    "dms-viz JSON files": {
        viz_name: rules.configure_dms_viz.output.dms_viz_json.format(viz_name=viz_name)
        for viz_name in dms_viz_config
    },
    "Notebooks prepping dms-viz JSONs": {
        viz_name: rules.configure_dms_viz.output.nb.format(viz_name=viz_name)
        for viz_name in dms_viz_config
    },
}

# Make row-wrapped heatmaps -------------------------------------------------------------

# read configuration for wrapped heatmaps
with open("data/wrapped_heatmap_config.yaml") as f:
    wrapped_heatmap_config = yaml.YAML(typ="safe", pure=True).load(f)


rule wrapped_heatmap:
    """Make row-wrapped heatmaps."""
    input:
        data_csv=lambda wc: wrapped_heatmap_config[wc.wrapped_hm]["data_csv"],
    output:
        chart_html="results/wrapped_heatmaps/{wrapped_hm}_wrapped_heatmap.html",
    params:
        params_dict=lambda wc: wrapped_heatmap_config[wc.wrapped_hm]
    log:
        notebook="results/notebooks/wrapped_heatmap_{wrapped_hm}.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    notebook:
        "notebooks/wrapped_heatmap.py.ipynb"

docs["Row-wrapped heatmaps"] = {
    "Heatmap HTMLs" : {
        wrapped_hm: rules.wrapped_heatmap.output.chart_html.format(wrapped_hm=wrapped_hm)
        for wrapped_hm in wrapped_heatmap_config
    }
}


# Make annotated CSVs ------------------------------------------------------------------

rule annotated_summary_csvs:
    """Make annotated summary CSVs with addtl annotations."""
    input:
        data="results/summaries/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding.csv",
        addtl_annotations="data/addtl_site_annotations.csv",
    output:
        mut="results/annotated_summary_csvs/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding_annotated.csv",
        site_mean="results/annotated_summary_csvs/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding_annotated_site_means.csv",
    log:
        "results/logs/annotated_summary_csvs.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    script:
        "scripts/annotated_summary_csvs.py"

docs["additional data CSVs"] = {
    "CSVs": {
        "annotated entry/binding mutation CSV": rules.annotated_summary_csvs.output.mut,
        "annotated entry/binding site mean CSV": rules.annotated_summary_csvs.output.site_mean,
    }
}


# Make some paper figures ---------------------------------------------------------------

rule paper_figures:
    """Make some paper figures."""
    input:
        **{
            f"func_scores_{s}": rules.func_scores.output.func_scores.format(selection=s)
            for s in func_scores
        },
        **{
            f"func_effects_{c}_{s}": f"results/func_effects/by_selection/{s}_func_effects.csv"
            for c in func_effects_config["avg_func_effects"]
            for s in func_effects_config["avg_func_effects"][c]["selections"]
        },
        codon_variants="results/variants/codon_variants.csv",
        annotated_mut_summary=rules.annotated_summary_csvs.output.mut,
        annotated_site_summary=rules.annotated_summary_csvs.output.site_mean,
        nb="notebooks/paper_figures.ipynb",
    output:
        nb="results/notebooks/paper_figures.ipynb",
    log:
        "results/logs/paper_figures.txt",
    params:
        params_yaml=lambda _, input: yaml_str(
            {
                "params": dict(input.items()),
                "min_times_seen": 2,
                "cell_entry_clip_lower": -6,
            }
        ),
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    shell:
        'papermill {input.nb} {output.nb} -y "{params.params_yaml}" &> {log}'

docs["Custom paper figures"] = {"Notebook w figures" : rules.paper_figures.output.nb}
