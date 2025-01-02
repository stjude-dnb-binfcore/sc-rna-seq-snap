# How to run the container for the Single cell RNA Seq Snap workflow (ScRNASeqSnap)

We have generated Dockerfile and Definition file that contain all tools, packages, and dependencies necessary to run the code and analyses modules in the `sc-rna-seq-snap` repository. These are customized for `Rstudio/R v4.4.0` and `Seurat v4.4.0`.


## To use the container in an R interactive session on HPC:

From an interactive node on HPC, the user can open their session. Please modify memory and resources as needed for the analysis module to run.
```
bsub -P hpcf_interactive -J hpcf_interactive -n 2 -q standard -R "rusage[mem=16G]" -Is "bash"
```

## Load specific version of Singularity

Please note that a version of Singularity is installed by default on all the cluster nodes at St Jude HPC. 
Otherwise the user needs to ensure and load Singularity module by running the following on HPC:
```
module load singularity/4.1.1
```


1. Pull the singularity container from the `sc-rna-seq-snap` root_dir
```
singularity pull docker://achronistjude/snap_rstudio_4.4.0_seurat_4.4.0:latest
```

2. Start the singularity container

a. To run from the terminal
```
bash run-terminal.sh
```

Then you may navigate to your module of interest, `./sc-rna-seq-snap/analyses/<module_of_interest>`. For example:
```
cd ./sc-rna-seq-snap/analyses/upstream-analysis
bash run-upstream-analysis.sh
```


b. To run from Rstudio
```
bash run-rstudio.sh
```

The `run-rstudio.sh` is running at `IP_ADDR:PORT`. When RStudio launches, please click "Session" -> "Restart R" (at the RStudio web session). For St Jude users, we advice to disconnect from CloudFlare WARP as this might lead to unstable behavior while on VPN.

Again, the user can navigate to their module of interest and explore/run their analyses.


If you encounter issues during this step related to RStudio Server and specifically to an invalid secure cookie error. This might be an issue with how the secure cookie is being handled during an HTTP request. In this case, please check if the following directories have been generated and if so, remove them:
```
rm -r .cache/
```
```
rm -r .config/
```
```
rm -r .local/
```
```
rm -r rstudio-container-tmp/
```

These folders cache history and user info. Then, kill the interactive session, start a new one, and hopefully, it works! 🎉


3. Build container (if needed)

If the user does not have access to the `snap_rstudio_4.4.0_seurat_4.4.0_latest.sif`, they can build their own. 
User can rename the `.sif` file, if they want to (not needed). Run the following from the `./run-container` dir:
```
singularity build snap_rstudio_4.4.0_seurat_4.4.0_latest.sif rstudio_r_4.4.0_seurat_4.4.0.def
```

Then, the user can start the container as explained in the step (3).


## To use the container outside of HPC and singularity:

1. Pull the docker container from the `sc-rna-seq-snap` root_dir
```
docker pull docker://achroni/snap_rstudio_4.4.0_seurat_4.4.0:latest
```

2. Start the docker container

To run from the terminal
```
docker run --platform linux/amd64 --name review -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/sc-rna-seq-snap docker://achroni/snap_rstudio_4.4.0_seurat_4.4.0:latest
```

```
docker container start review
```

```
docker exec -ti review bash
```

Navigate to your module of interest:
```
cd ./sc-rna-seq-snap/analyses/upstream-analysis
bash run-upstream-analysis.sh
```


## Authors

Antonia Chroni, PhD ([@AntoniaChroni](https://github.com/AntoniaChroni)), 
Walid Abu Al-Afia ([@walidabualafia](https://github.com/walidabualafia)), and
John Grime ([@JohnGrime](https://github.com/JohnGrime))




## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
