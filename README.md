# Pseudovirus deep mutational scanning of how mutations to the CHIKV envelope proteins affect entry in various cells
Study led by Xiaohui Ju in the [Bloom lab](https://jbloomlab.org/).

See [https://dms-vep.org/CHIKV-181-25-E-DMS/](https://dms-vep.org/CHIKV-181-25-E-DMS/) for the interactive HTML results of this pipeline.

## Organization of this repo

### `dms-vep-pipeline-3` submodule

Most of the analysis is done by the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), which was added as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to this pipeline via:

    git submodule add https://github.com/dms-vep/dms-vep-pipeline-3

This added the file [.gitmodules](.gitmodules) and the submodule [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), which was then committed to the repo.
Note that if you want a specific commit or tag of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) or to update to a new commit, follow the [steps here](https://stackoverflow.com/a/10916398), basically:

    cd dms-vep-pipeline-3
    git checkout <commit>

and then `cd ../` back to the top-level directory, and add and commit the updated `dms-vep-pipeline-3` submodule.
You can also make changes to the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) that you commit back to that repo.

### Code and configuration
The [snakemake](https://snakemake.readthedocs.io/) pipeline itself is run by `dms-vep-pipeline-3/Snakefile` which reads its configuration from [config.yaml](config.yaml).
The [conda](https://docs.conda.io/) environment used by the pipeline is that specified in the `environment.yml` file in [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3).

### Custom analyses run as part of pipeline
In addition to the core functionality of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), some additional custom analyses are run by the pipeline.
The notebooks for these analyses are in [./notebooks/](notebooks), and the code running them is in [custom_rules.smk](custom_rules.smk).

#### `dms-viz` visualization JSONs
These custom analyses include mappings for [dms-viz](https://dms-viz.github.io/) visualizations.
The configuration for those visualizations is specified in [./data/dms_viz_config.yaml](data/dms_viz_config.yaml), and uses the site maps that map protein sites to the PDB sites as defined in [data/pdb_sitemaps/](data/pdb_sitemaps/).

### [Nextstrain](https://nextstrain.org/) phylogeny
A `Snakemake` pipeline in [`./nextstrain`](./nextstrain/) builds an interactive phylogeny of the CHIKV structural polyprotein (Capisd, E3, E2, 6K, E1) using publicly available sequences from NCBI virus.
The final Auspice JSON produced by that pipeline is copied to to [auspice/CHIKV-181-25-E-DMS.json](auspice/CHIKV-181-25-E-DMS.json) where it can be viewed via  [Nextstrain community build](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html) at [nexstrain.org/community/dms-vep/CHIKV-181-25-E-DMS](nexstrain.org/community/dms-vep/CHIKV-181-25-E-DMS).

Note that this pipeline in [`./nextstrain`](./nextstrain/) is not run by the top-level `snakemake` file but must be run separately as described within that subdirectory.

#### Row-wrapped heatmaps
These custom analyses include making row-wrapped heatmaps that are more sized for paper figures.
The configuration for those visualizations is specified in [data/wrapped_heatmap_config.yaml](data/wrapped_heatmap_config.yaml), and HTMLs of the heatmaps are shown in the auto-rendered documentation.

#### manual_analyses
These custom analyses include making files for ChimeraX visualization of functional effects, Mxra8 binding on structure, analyzing cell entry, Mxra8 binding and validations.

#### Chimera_structures
These custom analyses include making files for ChimeraX visualization of functional effects on structure.
The file generated that can be read by ChimeraX is [`293T-Mxra8_entry_func_effects.defattr`](manual_analyses/chimera_structures/results/293T-Mxra8_entry_func_effects.defattr).

### Data
Input data utilized by the pipeline are located in [./data/](data). 

### Results and documentation
The results of running the pipeline are placed in [./results/](results).
Due to space, only some results are tracked. For those that are not, see the [.gitignore](.gitignore) document.

The pipeline builds HTML documentation for the pipeline in [./docs/](docs). These docs are rendered for viewing at [https://dms-vep.org/CHIKV-181-25-E-DMS/](https://dms-vep.org/CHIKV-181-25-E-DMS/).

## Running the pipeline
To run the pipeline, build the conda environment `dms-vep-pipeline-3` in the `environment.yml` file of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), activate it, and run [snakemake](https://snakemake.readthedocs.io/), such as:

    conda activate dms-vep-pipeline-3
    snakemake -j 16 --software-deployment-methd conda -s dms-vep-pipeline-3/Snakefile

To run on the Hutch cluster via [slurm](https://slurm.schedmd.com/), you can run the file [run_Hutch_cluster.bash](run_Hutch_cluster.bash):

    sbatch run_Hutch_cluster.bash
