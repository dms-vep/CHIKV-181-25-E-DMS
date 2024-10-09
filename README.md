# CHIKV 181/25 Envelope DMS

Experiments by Xiaohui Ju with analysis help by Will Hannon


Deep mutational scanning of the CHIKV envelope proteins from 181/25 strarin using a barcoded pseudotyped lentiviral platform
Study by Xiaohui Ju and Jesse Bloom.

This repo contains data and analyses from deep mutational scanning experiments on the CHIKV E protein. All experiments were performed on the 181/25 strain.


Organization of this repo
dms-vep-pipeline-3 submodule
Most of the analysis is done by the dms-vep-pipeline-3, which was added as a git submodule to this pipeline via:

git submodule add https://github.com/dms-vep/dms-vep-pipeline-3
This added the file .gitmodules and the submodule dms-vep-pipeline-3, which was then committed to the repo. Note that if you want a specific commit or tag of dms-vep-pipeline-3 or to update to a new commit, follow the steps here, basically:

cd dms-vep-pipeline-3
git checkout <commit>
and then cd ../ back to the top-level directory, and add and commit the updated dms-vep-pipeline-3 submodule. You can also make changes to the dms-vep-pipeline-3 that you commit back to that repo.

Code and configuration
The snakemake pipeline itself is run by dms-vep-pipeline-3/Snakefile which reads its configuration from config.yaml. The conda environment used by the pipeline is that specified in the environment.yml file in dms-vep-pipeline-3.

Data
Input data utilized by the pipeline are located in ./data/.

Results and documentation
The results of running the pipeline are placed in ./results/. Due to space, only some results are tracked. For those that are not, see the .gitignore document.

The pipeline builds HTML documentation for the pipeline in ./docs/. These docs are rendered for viewing at https://dms-vep.org/RABV_Pasteur_G_DMS/ as stated above.

Non-pipeline analyses
All other non-pipeline analyses are contained in ./scratch_notebook/. The notebooks in this directory are not part of the main pipeline but have been used to generate files used as input for the pipeline.

Running the pipeline
To run the pipeline, build the conda environment dms-vep-pipeline-3 in the environment.yml file of dms-vep-pipeline-3, activate it, and run snakemake, such as:

conda activate dms-vep-pipeline-3
snakemake -j 32 --use-conda -s dms-vep-pipeline-3/Snakefile
To run on the Hutch cluster via slurm, you can run the file run_Hutch_cluster.bash:

sbatch -c 32 run_Hutch_cluster.bash

