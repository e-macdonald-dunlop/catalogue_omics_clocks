#!/bin/bash


# Grid Engine options (lines prefixed with #$)
#$ -N dexa
#$ -cwd                  
#$ -l h_rt=08:00:00 
#$ -l h_vmem=16G
#$ -V
##$ -t 1-200
##$ -tc 20
#$ -hold_jid chunk
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem
# Initialise the environment modules
. /etc/profile.d/modules.sh
 
# Load Python
#module load R
#module load igmm/apps/plink/1.07  
module load igmm/apps/R/3.6.0
 #which plink
# Run the program


mkdir -p p04_clock_500_iterations

cd p04_clock_500_iterations

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a04_clock_500_iterations.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a04_clock_500_iterations_log.log

cd ../
