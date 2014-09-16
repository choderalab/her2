#!/bin/tcsh
#  Batch script for submitting serial job to MSKCC Torque/Moab cluster.
#
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=36:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
##PBS -k oe
#
# specify queue
#PBS -q gpu
#
# nodes: number of 32-hyperthread nodes
#   ppn: how many cores per node to use (1 through 32)
#       (you are always charged for the entire node)
#PBS -l nodes=1:ppn=1:gpus=1:exclusive
#
# export all my environment variables to the job
##PBS -V
#
# job name (default = name of script file)
#PBS -N setup-her2-siegetank
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
##PBS -o /cbio/jclab/home/chodera/ocores/output

cd "$PBS_O_WORKDIR"

#setenv LD_LIBRARY_PATH /usr/local/cuda-6.0/lib64/

date
hostname

python setup_mutations.py

date

