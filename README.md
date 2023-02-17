# cellXplore

`cellXplore` is a user-friendly interactive visualisation tool that allows researchers with no prior coding knowledge to explore cellular interactions inferred from single cell data. The basis of `cellXplore` is built within the `cellxgeneVIP` framework (https://github.com/interactivereport/cellxgene_VIP), an open source project that provides click and point functionality to interrogate single cell datasets. Although `cellxgeneVIP` is a powerful resource it doesn't address exploring cell-cell interactions (CCI), an increasingly popular downstream analysis step in scRNA-seq. 

`cellXplore` provides a shared web-based platform bringing together widely-used existing cell-cell interaction packages, allowing users to develop customisable analysis pipelines and interpret results with interactive data visualisations. `cellXplore` requires ••a fully pre-processed object containing ligand-receptor information•• from their single-cell data that can be generated from supported packages such as `CellPhoneDB`, `CellChat` and `NicheNetR`. Users then upload their pre-processed dataset, visualise and filter through results to find interesting CCI’s, streamlining the analysis with no complex bioinformatics necessary. 

# How to install cellXplore

Note: These installation instructions have been taken from `cellxgeneVIP` and contain altered `.yml` files containing the relevant additional packages required to run `cellXplore`.

## 1. Install `Anaconda` if it is not available on any server you may be using (https://docs.anaconda.com/anaconda/install/linux/)
``` bash
bash ~/Downloads/Anaconda3-2020.02-Linux-x86_64.sh
```

## 2. Create and enable the `cellXplore` conda environment from the command line by running the following lines of code
``` bash
git clone https://github.com/iii-cell-atlas/cellxgene_VIP.git
cd cellxgene_VIP

source <path to Anaconda3>/etc/profile.d/conda.sh (Default: /opt/anaconda3/etc/profile.d/conda.sh)
conda config --set channel_priority flexible
conda env create -n <env name, such as: VIP> -f VIP.yml (system-wide R) or VIP_conda_R.yml (local R under conda, no root privilege needed)

For Mac Users, conda env create -n <env name, such as: VIP> -f VIP.macOS.yml

conda activate <env name, such as: VIP>
or
source activate <env name>
```
## 3. Install the main `cellxgene` framework by running config.sh in "cellxgene_VIP" directory
```bash
./config.sh
For Mac User, ./config.macOS.sh
```
## 4. Install the required `R` packages
```bash
export LIBARROW_MINIMAL=false
#  Ensure that the right instance of `R` is used. e.g. system-wide: /bin/R or /usr/bin/R ; local R under conda: ~/.conda/envs/VIP_conda_R/bin/R
#  You can check this by running the following line of code in the terminal.
which R
#  The version of a `Bioconductor` package is controlled by `BiocManager`, whose version is provided

R -q -e 'if(!require(devtools)) install.packages("devtools",repos = "http://cran.us.r-project.org")'
R -q -e 'if(!require(Cairo)) devtools::install_version("Cairo",version="1.5-12",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(foreign)) devtools::install_version("foreign",version="0.8-76",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(ggpubr)) devtools::install_version("ggpubr",version="0.3.0",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(ggrastr)) devtools::install_version("ggrastr",version="0.1.9",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(arrow)) devtools::install_version("arrow",version="2.0.0",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(Seurat)) devtools::install_version("Seurat",version="3.2.3",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(rmarkdown)) devtools::install_version("rmarkdown",version="2.5",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(tidyverse)) devtools::install_version("tidyverse",version="1.3.0",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(viridis)) devtools::install_version("viridis",version="0.5.1",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(hexbin)) devtools::install_version("hexbin",version="1.28.2",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(ggforce)) devtools::install_version("ggforce",version="0.3.3",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(RcppRoll)) devtools::install_version("RcppRoll",version="0.3.0",repos = "https://cran.r-project.org")'
R -q -e 'if(!require(fastmatch)) devtools::install_version("fastmatch",version="1.1-3",repos = "https://cran.r-project.org")'
R -q -e 'if(!require(BiocManager)) devtools::install_version("BiocManager",version="1.30.10",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(fgsea)) BiocManager::install("fgsea")'
R -q -e 'if(!require(rtracklayer)) BiocManager::install("rtracklayer")'
R -q -e 'if(!require(rjson)) devtools::install_version("rjson",version="0.2.20",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(ComplexHeatmap)) BiocManager::install("ComplexHeatmap")'

# These should be already installed as dependencies of above packages
R -q -e 'if(!require(dbplyr)) devtools::install_version("dbplyr",version="1.0.2",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(RColorBrewer)) devtools::install_version("RColorBrewer",version="1.1-2",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(glue)) devtools::install_version("glue",version="1.4.2",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(gridExtra)) devtools::install_version("gridExtra",version="2.3",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(ggrepel)) devtools::install_version("ggrepel",version="0.8.2",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(MASS)) devtools::install_version("MASS",version="7.3-51.6",repos = "https://cran.us.r-project.org")'
R -q -e 'if(!require(data.table)) devtools::install_version("data.table",version="1.13.0",repos = "https://cran.us.r-project.org")'
```
## 5. Run `cellXplore` by specifiying a `.h5ad` file storing scRNA-seq data along with a host and a port, use "ps" to find used ports to spare, see https://chanzuckerberg.github.io/cellxgene/posts/launch for details.
```bash
ps -ef | grep cellxgene
Rscript -e 'reticulate::py_config()'
# Run the following command if the output of the above command doesn't point to the `Python` in your env.
export RETICULATE_PYTHON=`which python`
cellxgene launch --host <xxx> --port <xxx> --disable-annotations --verbose <h5ad file>
```
## 6. From a web browser (••Chrome is preferred••, Version 87.0.4280.88 or 87.0.4280.141 is used), access http(s)://host:port

You should be able to see this in Console of Chrome Developer Tools (Mac: Option+⌘+J, Windows/Linux: Shift+Ctrl+J) if everything is right.
![VIP_ready](https://user-images.githubusercontent.com/29576524/92059839-46482d00-ed60-11ea-8890-8e1b513a1656.png)

*Note: While spinning up the `cellXplore` from a HPC, do **NOT** use qlogin. **ssh directly to the server**.*



# Future Development

Currently we are working towards developing `cellXplore` to handle integrated scRNA-seq data with spatial data to build a comprehensive view of the interactome in its native context. 
