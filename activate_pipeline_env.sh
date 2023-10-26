#!/bin/bash 

module load bcftools/1.12 # this includes python version python3/3.9.2
source /g/data/pq84/python-3.9.2/bin/activate
module load java/jdk-8.40 
module load R/4.1.0 
export R_LIBS_USER="/g/data/pq84/R"
export PATH=$PATH:/g/data/pq84/bin/plink2/