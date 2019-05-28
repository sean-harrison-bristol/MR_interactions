#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=168:00:00

cd $PBS_O_WORKDIR
module add languages/R-3.4.1-ATLAS
Rscript factorial_mr_v4_h_vi.r 8376154



