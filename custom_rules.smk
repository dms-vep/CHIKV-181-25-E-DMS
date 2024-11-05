"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""

rule format_func_effects_dms_viz:
    """Format the data for functional effects for input into dms-viz"""
    input:
        site_map="data/site_numbering_map.csv",
        functional_data="results/func_effect_shifts/averages/{library}_comparison_shifts.csv",
        input_pdb_file="data/PDBs/6jo8.pdb",
    output:
        output_json_file_name = os.path.join("results/dms-viz/CHIKV181E_{library}_functional_effect_shifts.json"),
    params:
        env_chains = 'A C E B D F',
        name="effect_shift",
    log:
        os.path.join("results/logs/", "dms-viz_file_{library}_functional_effect_shifts.txt"),
    conda:
        "dms-viz.yml"
    shell:
        """
        configure-dms-viz format \
            --name {params.name} \
            --input {input.functional_data} \
            --metric  "shift" \
            --metric-name "Functional Effect Shift" \
            --exclude-amino-acids "*, -" \
            --structure {input.input_pdb_file} \
            --sitemap {input.site_map} \
            --output {output.output_json_file_name} \
            --included-chains "{params.env_chains}" \
            --tooltip-cols "{{'times_seen': '# Obsv'}}" \
            --filter-cols "{{'times_seen': 'Times Seen'}}" \
            --filter-limits "{{'times_seen': [0, 2, 10]}}" \
            &> {log}
        """

# Files (Jupyter notebooks, HTML plots, or CSVs) that you want included in
# the HTML docs should be added to the nested dict `docs`:
other_target_files.append(os.path.join("results/dms-viz/CHIKV181E_E3E2_functional_effect_shifts.json"))
other_target_files.append(os.path.join("results/dms-viz/CHIKV181E_6KE1_functional_effect_shifts.json"))