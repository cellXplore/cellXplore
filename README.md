# cellXplore

`cellXplore` is a user-friendly interactive visualisation tool that allows researchers with no prior coding knowledge to explore cellular interactions inferred from single cell data. The basis of `cellXplore` is built within the `cellxgeneVIP` framework (https://github.com/interactivereport/cellxgene_VIP), an open source project that provides click and point functionality to interrogate single cell datasets. Although `cellxgeneVIP` is a powerful resource it doesn't address exploring cell-cell interactions (CCI), an increasingly popular downstream analysis step in scRNA-seq. 

`cellXplore` provides a shared web-based platform bringing together widely-used existing cell-cell interaction packages, allowing users to develop customisable analysis pipelines and interpret results with interactive data visualisations. `cellXplore` requires ••a fully pre-processed object containing ligand-receptor information•• from their single-cell data that can be generated from supported packages such as `CellPhoneDB`, `CellChat` and `NicheNetR`. Users then upload their pre-processed dataset, visualise and filter through results to find interesting CCI’s, streamlining the analysis with no complex bioinformatics necessary. 

# How to install cellXplore

Note: These installation instructions have been taken from `cellxgeneVIP` and contain altered `.yml` files containing the relevant additional packages required to run `cellXplore`.

## 1. Install `Anaconda` if it is not available on any server you may be using (https://github.com/conda-forge/miniforge#mambaforge)
``` bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

## 2. Create and enable the `cellXplore` conda environment from the command line by running the following lines of code
``` bash
git clone https://github.com/cellXplore/cellXplore.git
cd cellXplore

source <path to Anaconda3>/etc/profile.d/conda.sh (Default: /opt/anaconda3/etc/profile.d/conda.sh)
conda config --set channel_priority flexible
conda env create -n cellXplore -f cellXplore.yml (local R under conda, no root privilege needed)
mamba env update -f cellXplore_r_dependencies.yml --name cellXplore

conda activate cellXplore
or
source activate cellXplore
```
## 3. Install the main `cellxgene` framework by running config.sh in "cellxgene_VIP" directory
```bash
./config.sh
For Mac User, ./config.macOS.sh
```
## 4. Run `cellXplore` by specifiying a `.h5ad` file storing scRNA-seq data along with a host and a port, use "ps" to find used ports to spare, see https://chanzuckerberg.github.io/cellxgene/posts/launch for details.
```bash
ps -ef | grep cellxgene
Rscript -e 'reticulate::py_config()'
# Run the following command if the output of the above command doesn't point to the `Python` in your env.
export RETICULATE_PYTHON=`which python`
cellxgene launch --host <xxx> --port <xxx> --disable-annotations --verbose <h5ad file>
```
## 5. From a web browser (••Chrome is preferred••, Version 87.0.4280.88 or 87.0.4280.141 is used), access http(s)://host:port

You should be able to see this in Console of Chrome Developer Tools (Mac: Option+⌘+J, Windows/Linux: Shift+Ctrl+J) if everything is right.
![VIP_ready](https://user-images.githubusercontent.com/29576524/92059839-46482d00-ed60-11ea-8890-8e1b513a1656.png)

*Note: While spinning up the `cellXplore` from a HPC, do **NOT** use qlogin. **ssh directly to the server**.*



# Future Development

Currently we are working towards developing `cellXplore` to handle integrated scRNA-seq data with spatial data to build a comprehensive view of the interactome in its native context. 
