#!/bin/bash

# job name
#SBATCH -J NNBAR_simulation_valgrind_check
#SBATCH -N 1
#SBATCH --tasks-per-node=20
#SBATCH --mail-user=sze-chun.yiu@fysik.su.se
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=4800

# job time, change for what your job requires
#SBATCH -t 12:00:00

# filenames stdout and stderr - customise, include %j
#SBATCH -o ./log/process_%j.out
#SBATCH -e ./log/process_%j.err

# write this script to stdout-file - useful for scripting errors
cat $0

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

export PATH="/home/scyiu/valgrind-3.17.0/build/bin:$PATH"
valgrind --tool=memcheck --leak-check=yes ./nnbar-calo-sim -m ./macro/cosmic_macro/cosmic_run_main.mac -t 2
