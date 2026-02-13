<p align="center";">
  <img src="figures/img/SCRNA_Logo_Primary.png" alt="ScRNASeqSnap repository logo" width="560px" />
</p>
<p align="center";">
  <a href="https://www.repostatus.org/#active">
    <img src="https://www.repostatus.org/badges/latest/active.svg?style=for-the-badge" alt="The project has reached a stable, usable state and is being actively developed." />
  </a>
  <a href="https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap">
    <img src="https://img.shields.io/badge/version-1.0.0.beta.7-brightgreen" alt="Version" />
  </a>
</p>

# Single cell RNA Seq Snap workflow (ScRNASeqSnap)

Snap is a comprehensive suite of tools and workflows for analyzing single-cell and single-nucleus RNA (sc/snRNA) data from 10X Genomics sequencing technology supporting human, mouse, and dual genome cohorts. Snap is an initiative of the [Bioinformatics Core](https://www.stjude.org/research/departments/developmental-neurobiology/shared-resources/bioinformatic-core.html) at the Department of Developmental Neurobiology at the St. Jude Children's Research Hospital.


## Table of Contents

1. [Getting Started](#getting-started)
2. [Installation](#installation)
3. [Tutorial and Documentation](#tutorial-and-documentation)
4. [Preparing project metadata](#preparing-project-metadata)
5. [How to Use the Repository](#how-to-use-the-repository)
   - [Accessing the Code](#accessing-the-code)
   - [Running the Code](#running-the-code)
6. [Requesting Resources from the HPCF Cluster](#requesting-resources-from-the-hpcf-cluster)



## Getting Started

### Installation

To begin using the Snap pipeline, follow the instructions below to set up the environment and run the code. A pre-built [Docker image](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/blob/main/run-container/README.md) is available for easy setup, containing all the necessary tools, packages, and dependencies to seamlessly run the code and analysis modules. 

### Tutorial and Documentation

For a step-by-step guide on how to access the code, run the analysis, and request memory from the HPCF cluster, refer to the current README file or the [Snap wiki page](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/wiki). Training sessions can also be provided upon request for St Jude users.

> [Tutorial and documentation for the snap pipeline](https://github.com/stjude-dnb-binfcore/trainings/blob/main/courses/sc-rna-seq-snap-repo/tutorial/snap-tutorial-docs). (legacy): This material reflects a previous version of the Snap pipeline and is no longer actively maintained or supported. It is retained for reference purposes only.



### Preparing project metadata

The pipeline requires a TSV file containing essential metadata for cohort analysis. The file must be named `project_metadata.tsv`. It can include one or more samples, as long as it contains at least the following columns in this exact order: `ID`, `SAMPLE`, and `FASTQ`. The `ID` column must contain unique values. The `SAMPLE` column must contain the `seq_submission_code` along with the ID, e.g., `seq_submission_code1_sample1` or the corresponding library name. The `FASTQ` column must contain the file path to the fastq files. For samples with top-ups or multiple technical replicates, list all associated library names and FASTQ file paths in the same row, using commas to separate each path. Additional metadata columns can be added and arranged as needed by the user (though not required).

The pipeline can be run on both single-sample and multi-sample cohorts.

The file can be stored anywhere, but its filepath must be specified in the `project_parameters.Config.yaml` file.

For user convenience, an example [project_metadata.tsv](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/blob/main/data/project_metadata) file is provided.


### How to Use the Repository

#### Accessing the Code

We recommend that users fork the `sc-rna-seq-snap` repository and then clone their forked repository to their local machine. Team members should use the [stjude-dnb-binfcore](https://github.com/stjude-dnb-binfcore) account, while others can use their preferred GitHub account. We welcome collaborations, so please feel free to reach out if you're interested in being added to the `stjude-dnb-binfcore` account.

1. Fork the repository

Navigate to the main page of the stjude-dnb-binfcore/sc-rna-seq-snap repository and click the "Fork" button.

<img width="650" alt="how-to-fork-repo-1" src="https://github.com/user-attachments/assets/1fc0a459-2c8c-4d2e-ab6b-6abaafae963e">


2. Create Your Fork

You can change the name of the forked repository (optional - unless you will use it for multiple projects). Click "Create fork" to proceed.


<img width="650" alt="how-to-fork-repo-2" src="https://github.com/user-attachments/assets/914a3db5-6e87-41fb-baf2-a50ffdb2a7c0">


3. Enjoy your new project repo!

<img width="650" alt="how-to-fork-repo-3" src="https://github.com/user-attachments/assets/073abb78-3993-4527-a574-859fd3046d39">


4. Clone Your Fork

Once you have created the fork, clone it to your local machine:

```
git clone https://github.com/<FORK_NAME>.git
```

#### Running the Code

1. Configure Your Parameters

Replace the project_parameters.Config.yaml file with your own file paths and parameters.

2. Navigate to an Analysis Module

Change to the relevant directory and run the desired shell script:

```
cd ./sc-rna-seq-snap/analyses/<module_of_interest>
```

3. Sync Your Fork

User needs to ensure that the main branch of the forked repository is always up to date with `stjude-dnb-binfcore/sc-rna-seq-snap:main`. 

If your fork is behind the main repository (`stjude-dnb-binfcore/sc-rna-seq-snap:main`), sync it to ensure you have the latest updates. This will update the main branch of your project repo with the new code and modules (if any). This will add code and not break any analyses already run in your project repo. 

When syncing your forked repository with the main repository, please be cautious of any changes made to the following files, as they are typically modified and specified for project data analysis:

   - `project_parameters.Config.yaml`

Before pulling the latest changes, stash any modifications you have made to these files. This ensures that you won't accidentally overwrite your changes when syncing with the main repository. 

Some useful git commands:

```
git branch
git checkout main
git config pull.rebase false

git status
git add project_parameters.Config.yaml
git commit -m "Update yaml"
```

Finally, `git pull` to get the most updated changes and code in your project repo. Please be mindful of any local changes in files in your project repo that you have done, e.g., `project_parameters.Config.yaml`. You will need to commit or stash (or restore) the changes to the yaml before completing the pull.

```
git pull
```

### Requesting CPU and Memory Resources

While we provide estimates for the computational resources required (based on 8 samples with approximately 50,000 cells), users may need to adjust memory settings based on cohort size and analysis requirements.

Important Considerations:

  - Adjust memory requests according to the size of your cohort and specific analysis needs.
  - For St. Jude users:
    - Refer to the [Introduction to the HPCF cluster](https://wiki.stjude.org/display/HPCF/Introduction+to+the+HPCF+cluster#IntroductiontotheHPCFcluster-queuesQueues:) for detailed guidance.
    - If you require more than 1 TB of memory, use the `large_mem` queue to ensure proper resource allocation.
  

### Below is the main directory structure listing the analyses and data files used in this repository

```
├── analyses
|  ├── cell-contamination-removal-analysis
|  ├── cell-types-annotation
|  ├── cellranger-analysis
|  ├── clone-phylogeny-analysis
|  ├── cluster-cell-calling
|  ├── de-go-analysis
|  ├── fastqc-analysis
|  ├── integrative-analysis
|  ├── project-updates
|  ├── README.md
|  ├── rshiny-app
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

*These tools and pipelines have been developed by the Bioinformatics core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
