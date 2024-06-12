## conda (recommended)
```
conda install -c conda-forge -c bioconda syngap
```
## manually
```
cd ~/code  # or any directory of your choice
git clone git://github.com/yanyew/SynGAP.git
cd ~/code/SynGAP
conda env create -f SynGAP.environment.yaml -c conda-forge -c bioconda
export PATH=~/code/SynGAP:$PATH
```
## Docker image
```
docker pull yanyew/syngap:1.1.0
docker run -it yanyew/syngap:1.1.0
conda activate syngap # activate the conda environment for SynGAP
```
