#!/bin/bash
#PBS -V
#PBS -q cpu
#PBS -j n
#PBS -l select=1:ncpus=16:mpiprocs=16:mem=32gb
#PBS -N chembl_pre
#PBS -m a

module load mpi/mpich-x86_64
cd $PBS_O_WORKDIR
hostname

last=$(date "+%s")

source ~/.bashrc
conda activate qvina

python run_qvina-w.py -c ../input/config.yaml > log.txt

current=$(date "+%s")
s=$((current - last))
echo "Elapsed Time: ${s} sec."