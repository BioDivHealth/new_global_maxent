#!/bin/bash -l

#$ -S /bin/bash

# wallclock time (format hours:minutes:seconds).
#$ -l h_rt=47:10:0

# RAM.
#$ -l mem=4G

# TMPDIR space (default is 10 GB)
#$ -l tmpfs=10G

# name of the job.
#$ -N R_SDMS_1

# Set the working directory to somewhere in your scratch space.
# Replace "<your_UCL_id>" with your UCL user ID
#$ -wd /home/ucbtdw0/Scratch/sdms/legion_output

# Run the application.
cd $TMPDIR

module unload compilers
module unload mpi
module load r/new


Rscript $HOME/Scratch/sdms/code3.R

