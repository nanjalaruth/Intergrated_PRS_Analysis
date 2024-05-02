# Intergrated_PRS_Analysis (Nextflow pipeline)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

The pipeline downloads PGS Score files from the PGS Catalogue based on PGS ids for different traits.
It then computes PGS scores using PLINK an combines the scores using ElasticNet. 
It then calculates a prediction score and performs association analysis on lipid traits.

## Installation 

1. Nextflow
```
wget -qO- https://get.nextflow.io | bash
```

## Running the pipeline
The pipeline does not require installation as `NextFlow` will automatically fetch it from `GitHub`.

### Own data
Start running your own analysis either by using flags as shown below:

 Run your own analysis by modifying the conf/test.config file to suit the path to your data location and then run the command as below:
 
 ```
 nextflow run nanjalaruth/hla_typing_using_HLA-LA -profile slurm -c <path to your edited config file> -resume
 ```
    
## To run the updated version of this pipeline, run:

 ```
 nextflow pull nanjalaruth/hla_typing_using_HLA-LA
 ```

## Support
I track open tasks using github's [issues](https://github.com/nanjalaruth/Intergrated_PRS_Analysis/issues)
