FROM rocker/rstudio:4.4.0

#
# Ensure virtual env on PATH so we're using the right venv binaries
#
ENV LC_ALL=C
ENV VIRTUAL_ENV=/opt/container-venv
ENV PATH=/opt/container-venv/bin:$PATH

#
# Ensure core system software is up to date and ready. Need certs for Cloudflare!
#

RUN apt-get -y update && apt-get -y install ca-certificates gdebi-core python3-pip python3-venv

#
# Install Cloudflare certificates so Docker build isn't interrupted by
# certificate errors. Combination of certificates should be valid until 2029.
#

COPY "./cloudflare_certs/original-certificate.pem" /usr/local/share/ca-certificates/2024-certificate.crt
COPY "./cloudflare_certs/SJ-2025-2029-CF-WARP-certificate.pem" /usr/local/share/ca-certificates/2025-2029-certificate.crt
RUN update-ca-certificates



#
# Set up virtual env; we need this to keep the existing numpy/scipy/pandas
# versions cohesive; otherwise, we get mistmatches with e.g. multiqc in the
# analysis pipelines. After activating the virtual env, we should have:
#
# (container-venv) Apptainer> pip3 list
# Package         Version
# --------------- -----------
# numpy           2.2.1
# pandas          2.2.3
# scipy           1.14.1
#
# Note: need to "source" the virtual env activate script, not just run it.
#

RUN python3 -m venv /opt/container-venv
RUN . /opt/container-venv/bin/activate

# Avoids BiocManager errors in Rhtslib; "fatal error: lzma.h: No such file or directory" etc
RUN apt-get -y install liblzma-dev libbz2-dev

# Avoids BiocManager errors in bluster; "libglpk.so.40: cannot open shared object file: No such file or directory"
RUN apt-get -y install libglpk-dev

# Avoids BiocManager errors in infercnv: "libjags.so.4: cannot open shared object file: No such file or directory"
RUN apt-get -y install r-cran-rjags

# Avoids errors in "Rscript -e 'devtools::install_github("welch-lab/RcppPlanc")'"
RUN apt-get -y install cmake libhdf5-dev

# Needed for pipeline stages
RUN apt-get install -y openjdk-8-jre-headless fastqc multiqc curl

# Needed for the cluster-cell-calling stage etc. Pandas as per Antonia's request.
RUN pip3 install pandas leidenalg


#
#
# Install texlive full version
#
ENV TEXLIVE_INSTALL_NO_CONTEXT_CACHE=1 \
    NOPERLDOC=1

# Install TexLive prerequisites
RUN apt-get update && \
    apt-get install -y wget unzip tar \
    make fontconfig perl openjdk-8-jre libgetopt-long-descriptive-perl \
    libdigest-perl-md5-perl libncurses5 \
    python3-pygments && \
    rm -rf /var/lib/apt/lists/*

# Install TexLive to the container
RUN wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz && \
    tar xzf install-tl-unx.tar.gz && rm install-tl-unx.tar.gz && \
    cd install-tl* && \
    echo "selected_scheme scheme-full" > install.profile && \
    echo "tlpdbopt_install_docfiles 0" >> install.profile && \
    echo "tlpdbopt_install_srcfiles 0" >> install.profile && \
    echo "tlpdbopt_autobackup 0" >> install.profile && \
    echo "tlpdbopt_sys_bin /usr/bin" >> install.profile && \
    ./install-tl -profile install.profile && cd .. && rm -rf install-tl*

# Symlink TexLive binaries to PATH
RUN /usr/local/texlive/2024/bin/x86_64-linux/tlmgr path add
#


#
# Use latest R sources to avoid e.g. "matrixStats" version too low for "MatrixGenerics"
#

RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest'))" >> /usr/local/lib/R/etc/Rprofile.site

#
# Install required R package managers
#


RUN Rscript -e 'install.packages("usethis")'
RUN Rscript -e 'install.packages(c("BiocManager", "devtools", "remotes"))'

# Debugging step: Verify installations
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) stop("BiocManager not installed!")'
RUN Rscript -e 'if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools not installed!")'
RUN Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) stop("remotes not installed!")'


#
# This is a fairly time consuming section; lots built from source.
#

RUN Rscript -e 'BiocManager::install(c("celldex"))'
RUN Rscript -e 'BiocManager::install(c("miQC"))'
RUN Rscript -e 'BiocManager::install(c("infercnv"))'
RUN Rscript -e 'BiocManager::install(c("scater"))'
RUN Rscript -e 'BiocManager::install(c("scDblFinder"))'
RUN Rscript -e 'BiocManager::install(c("SingleCellExperiment"))'
RUN Rscript -e 'BiocManager::install(c("SingleR"))'

#
# Some github-based R modules; RcppPlanc needs previously installed cmake
# and libhdf5-dev
#

RUN Rscript -e 'devtools::install_github("igordot/scooter")'
RUN Rscript -e 'devtools::install_github("welch-lab/RcppPlanc")'
RUN Rscript -e 'devtools::install_github("kharchenkolab/numbat")'
RUN Rscript -e 'devtools::install_github("welch-lab/liger")'

#
# Removed "grid" from installs as triggered the warning message:
# "package 'grid' is a base package, and should not be updated"
# "grid" was then subsequently ignored by the installer anyway.
#

RUN Rscript -e 'install.packages(c("clustree", "cowplot", "data.table"))'
RUN Rscript -e 'install.packages(c("flexmix", "flextable", "forcats"))'
RUN Rscript -e 'install.packages(c("fs", "future", "GGally"))'
RUN Rscript -e 'install.packages(c("ggh4x", "ggplot2", "ggpmisc"))'
RUN Rscript -e 'install.packages(c("ggrepel", "ggthemes"))'
RUN Rscript -e 'install.packages(c("harmony", "igraph", "irlba"))'
RUN Rscript -e 'install.packages(c("knitr", "optparse"))'
RUN Rscript -e 'install.packages(c("patchwork", "purrr", "RColorBrewer"))'
RUN Rscript -e 'install.packages(c("remotes", "reshape2", "rliger"))'
RUN Rscript -e 'install.packages(c("rlist", "R.utils", "SeuratObject"))'
RUN Rscript -e 'install.packages(c("shiny", "SoupX", "stringr"))'
RUN Rscript -e 'install.packages(c("tidytext", "tidyverse", "tinytex", "yaml"))'

# Prevent errors in dependencies for the following
RUN apt-get -y update

#
# Was Rscript -e 'remotes::install_version("Seurat", "4.4.0", repos = c("https://packagemanager.posit.co/cran/__linux__/jammy/latest"))'
# but we specified the latest CRAN repo earlier so repo URL is redundant.
# Default version of Seurat is higher than 4.4.0 and is incompatible with
# input object format so install specific version. We also need to fix the
# SeuratWrappers, as the existing version is incompatible with 4.4.0:
#
# https://github.com/satijalab/seurat-wrappers/issues/169
#
# We need to install the Seurat wrappers first, as otherwise they will
# override the install of Seurat 4.4.0 with a later version of Seurat.
#
RUN Rscript -e 'remotes::install_github("satijalab/seurat-wrappers@community-vignette")'
RUN Rscript -e 'remotes::install_version("Seurat", version = "4.4.0")'

#
# As per Antonia's request via Teams, 12/16/2024 09:46.
#

#RUN Rscript -e 'install.packages("Signac")'
RUN Rscript -e 'devtools::install_github("stuart-lab/signac", ref = "develop")'


#
# Install CellRanger v8.0.1
#
# You will need to download the desired version from here: https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions
# Save on the same folder where Dockerfile is located
# You should also replace the link with the one from the most recent download as the url expires after 1-2 days
# Set environment variables
ENV PATH=/opt/cellranger-8.0.1/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/cellranger-8.0.1/lib

RUN curl -o cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1735637913&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=g0mi3a~7RGhqh51x1p6qL0rJ7gvx3~NVWlpW5pg4E7D0xWTISNpa1kmb6AS~qxAiP2mLwDBqjfH9NLmzelBmTx02xRZ3YfsKikcPPRn9KyB8D5Ki2lrIUZkChkfuLdDkoxAfmBOBv~kLvsuefV~2k8-B1WRNCPcoxryhdbbK3-b~FkZgnP~Fe0h5rI~5757jV2r1I6yq3xTs0ayw0WDyL49apdI0bIi4E7-lCPohbnGop0Z~G1Y~QMTR3TN98e8mEKWfmZjboSgG7triA8jDY3g1hbSI5pQh8yzHfww3fTcFN8P1ksjX0NH~~UWNeepSJwPA11MfoQ3cUq3SLvazuA__"
RUN tar -xvzf cellranger-8.0.1.tar.gz -C /opt
#


# Run virtual end on container startup- adjust as needed.
ENTRYPOINT . /opt/container-venv/bin/activate

