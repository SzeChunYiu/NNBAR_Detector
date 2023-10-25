#!/bin/bash

# job name
#SBATCH -p hep
#SBATCH -A hep2023-1-3
#SBATCH -J HIBEAM_neutron5
#SBATCH -N 1
#SBATCH --tasks-per-node=20

# job time, change for what your job requires
#SBATCH -t 120:00:00

# filenames stdout and stderr - customise, include %j
#SBATCH -o ./log/cosmic_neutron_5_%j.out
#SBATCH -e ./log/cosmic_neutron_5_%j.err

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
source /projects/hep/fs12/nnbar/software/geant4-MT/install/bin/geant4.sh

./nnbar-calo-sim -m ./macro/cosmic_macro/cosmic_neutron/run_5.mac -t 20
