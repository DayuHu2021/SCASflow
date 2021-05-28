#########
##2021-3-11
###Dayu,Hu
##SelectK with snakemake

pcs = {"10","20","30","40","50"}
resolutions = {"0.1","0.5","1","1.5","2","4","8","16"}

rule all:
    input:
        expand("results/ClusterSilhouette_resolution_{res}_PC_{pc}.rds",res=resolutions, pc = pcs)

rule SubSample:
    input:
        "rawdata/sample.rds"
    output:
        "results/SubSample_resolution_{res}_PC_{pc}.rds"
    script: "scripts/SubSample.R"

rule CoClustering:
    input:
        "results/SubSample_resolution_{res}_PC_{pc}.rds"
    output:
        "results/CoClustering_resolution_{res}_PC_{pc}.rds"
    script: "scripts/CoClustering.R"

rule CellSilhouette:
    input:
        "results/CoClustering_resolution_{res}_PC_{pc}.rds"
    output:
        "results/CellSilhouette_resolution_{res}_PC_{pc}.rds"
    script: "scripts/CellSilhouette.R"

rule ClusterSilhouette:
    input:
        "results/CellSilhouette_resolution_{res}_PC_{pc}.rds"
    output:
        "results/ClusterSilhouette_resolution_{res}_PC_{pc}.rds"
    script: "scripts/ClusterSilhouette.R"




