# SACSflow

The goal of SelectK is to assess the stability of single-cell clusters by subsampling and clustering the cells, providing visualization method for comparing clusters and selecting the optimum suggested parameter for clustering.

For now, three steps are performed.

* Selection of optimal **resolution**
* Selection of optimal number of **PCs** (principle components) 
* **algorithm** for modularity optimization



**Upstream analysis** relies on Snakemake, Snakemake workflow was implemented to perform the subsampling and clustering steps

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



**Downstream analysis** relies on script skVisualization

the whole process  is as follows.

Step1,performing subsampling 100 times, and each time 80% of the original dataset is randomly selected as a subset. 

Step2,performing coclustering and constructing the matrix of copolymerization frequency, copolymerization frequency reflects the similarity between cells. The similarity distance is defined as 1 minus the value of copolymerization frequency.

Step3, calculating silhouette score using the similarity distance.

Step4, assessing the stability of single-cell clusters. Cluster whose silhouette score is higher than 0.7 is regarded as stable cluster. 

Step5,selecting the best resolution. The clustering results with stable cluster proportion less than 50% are regarded as invalid. In the results meeting the above constraints, the parameter with the highest resolution is selected.

The workflow is:



<img src="/Users/kabuda/Library/Application Support/typora-user-images/image-20210528101040926.png" alt="image-20210528101040926"  />

 **Datasets**

1. The dataset of 3k Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics can be found at https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

2. The smart-seq PBMC dataset can be found at https://singlecell.broadinstitute.org/single_cell/study/SCP424/

3. A mixture control dataset can be found at https://github.com/LuyiTian/sc_mixology
