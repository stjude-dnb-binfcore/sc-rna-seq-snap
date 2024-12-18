<p align="center">
  <img src="figures/img/SCRNA_Logo_Primary.png" alt="ScRNASeqSnap repository logo" width="560px" />
</p>
<p align="center">
  <a href="https://www.repostatus.org/#active"><img src="https://www.repostatus.org/badges/latest/active.svg?style=for-the-badge" alt="The project has reached a stable, usable state and is being actively developed." /></a>

</p>

# Single cell RNA Seq Snap workflow (ScRNASeqSnap)

This repository contains tools and workflows for analyzing single cell and single nuclei RNA (sc/snRNA) data from 10X sequencing technology. 

There are a Singularity and Docker image that contain all tools, packages, and dependencies necessary to run the code and analyses modules in the `sc-rna-seq-snap` repository. For more information on how to run the container and package versioning, please see [./run-container/README.md](https://github.com/AntoniaChroni/sc-rna-seq-snap/blob/main/run-container/README.md).

The `sc-rna-seq-snap` repository is an initiative of the [Bioinformatics Core at the Department of Developmental Neurobiology at the St. Jude Children's Research Hospital](https://www.stjude.org/research/departments/developmental-neurobiology/shared-resources/bioinformatic-core.html).


⚠️ 🚧 🚧 ⚠️
The repo is currently under development and code review process. Analysis modules in the main branch are ready to be used.


## To access the code in this repository:

### Clone option
1. Clone the repository
```
git clone https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap.git
```

### Fork option
1. Fork the repository on your own account from the main page of the `stjude-dnb-binfcore/sc-rna-seq-snap` by clicking the “Fork” button

<img width="650" alt="how-to-fork-repo-1" src="https://github.com/user-attachments/assets/1fc0a459-2c8c-4d2e-ab6b-6abaafae963e">


2. Change the name if you like, but probably not; click “Create fork”

<img width="650" alt="how-to-fork-repo-2" src="https://github.com/user-attachments/assets/914a3db5-6e87-41fb-baf2-a50ffdb2a7c0">


3. Enjoy your new project repo!

<img width="650" alt="how-to-fork-repo-3" src="https://github.com/user-attachments/assets/073abb78-3993-4527-a574-859fd3046d39">


## To run the code in this repository:

1. Replace the `project_parameters.Config.yaml` with your file paths and parameters.

2. Navigate to an analysis module and run the shell script of interest:
```
cd ./sc-rna-seq-snap/analyses/<module_of_interest>
```

**If you have forked the repo, you will need to do the following steps before running the script of interest. If you have cloned the repo, you can skip this.**

You need to do `sync fork` of your project repo at GitHub before running a module, if your branch is behind the main branch of the `stjude-dnb-binfcore/sc-rna-seq-snap:main`. This will update the main branch of your project repo with the new code and modules (if any). This will add code and not break any analyses already run in your project repo. 

Then navigate to your `./sc-rna-seq-snap` project repo and ensure you are at the `main` branch (in case you have featured branches). If not, you will need to `git checkout` to the main branch.
```
git branch
git checkout main
```

Finally, `git pull` to get the most updated changes and code in your project repo. Please be mindful of any changes in files in your project repo that you have done, e.g., `project_parameters.Config.yaml`. You will need to commit or stash (or restore) the changes to the yaml before completing the pull.
```
git pull
```


### Below is the main directory structure listing the analyses and data files used in this repository

```
├── analyses
|  ├── cell-contamination-removal-analysis
|  ├── cell-types-annotation
|  ├── cellranger-analysis
|  ├── cluster-cell-calling
|  ├── fastqc-analysis
|  ├── integrative-analysis
|  ├── README.md
|  └── upstream-analysis
├── figures
├── LICENSE
├── project_parameters.Config.yaml
├── README.md
├── run-container
├── run-rstudio.sh
├── run-terminal.sh
└── SECURITY.md
```

## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
