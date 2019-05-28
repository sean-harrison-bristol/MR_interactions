#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=168:00:00

cd $PBS_O_WORKDIR
module add languages/R-3.4.1-ATLAS
Rscript meta_sensitivity_v8.r 267639401



