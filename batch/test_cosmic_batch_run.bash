#!/bin/bash

# job name
#SBATCH -J NNBAR_simulation
#SBATCH -N 2
#SBATCH --tasks-per-node=10
 
# job time, change for what your job requires
#SBATCH -t 1:00:00

# filenames stdout and stderr - customise, include %j
#SBATCH -o ./log/process_%j.out
#SBATCH -e ./log/process_%j.err

# write this script to stdout-file - useful for scripting errors
cat $0

#SBATCH --mail-user=sze-chun.yiu@fysik.su.se
#SBATCH --mail-type=ALL

module load anaconda3/4.4.0
module load CMake

module load GCCcore/11.2.0  
module load Qt5/5.15.2

module load GCC/11.2.0 intel-compilers/2021.4.0 
module load Boost/1.77.0

module load OpenMPI/4.1.1
module load Arrow/6.0.0

source activate nnbar_env
source /home/scyiu/nnbar/geant4-MT/install/bin/geant4.sh

# Copying the executable onto the local disks of the nodes
srun -n $SLURM_NNODES -N $SLURM_NNODES cp -r * $SNIC_TMP

# Copy the input file onto the headnode - if your MPI program
# reads from all tasks, you need to do the above srun construct
# again
#cp -p input.dat $SNIC_TMP

# change to local disk and start the mpi executable
cd $SNIC_TMP
./nnbar-calo-sim -m ./macro/cosmic_macro/test_macro/test_muon.mac -t 20

# Copy result files back - example assumes only task zero writes
# if in your application result files are written on all nodes
# you need to initiate a copy on each node via srun
ls
#cp -rT * $SLURM_SUBMIT_DIR/output/
rsync -av --delete ./output/* $SLURM_SUBMIT_DIR/output/
