# SelectK  

Snakemake workflow was implemented to perform the subsampling and clustering steps

For now, three steps are performed.

* Selection of optimal **resolution**
* Selection of optimal number of **PCs** (principle components) 
* **algorithm** for modularity optimization

on server:

```bash
git clone https://github.com/DayuHu2021/SelectK

conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake

source activate snakemake

# R4.0.3, make sure you load R after source activate conda environment
module load R


R
>install.package("Seurat")



```





```bash
# dry run
snakemake -np 

# real run
snakemake -s snakemakefile(here is SelectK.py)


```
