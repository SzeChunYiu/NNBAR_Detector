# Building simulation from source
cmake -DGeant4_DIR=/home/scyiu/nnbar/geant4-MT/install/lib/Geant4-10.7.2 -DCMAKE_INSTALL_PREFIX=<path-to-build> -DCRY_PATH=/home/scyiu/nnbar/cry_v1.7/ -DCMAKE_MODULE_PATH=/sw/easybuild/software/Arrow/6.0.0-foss-2021b/lib/cmake/arrow/ -DVERSION_MCPL=1 -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON <path-to-source>

Examples:

cmake -DGeant4_DIR=/home/scyiu/nnbar/geant4-MT/install/lib/Geant4-10.7.2 -DMCPL_BUILD=1 -DTARGET_BUILD=0 -DCMAKE_INSTALL_PREFIX=. -DCMAKE_MODULE_PATH=/sw/easybuild/software/Arrow/6.0.0-foss-2021b/lib/cmake/arrow/ -DVERSION_MCPL=1 -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON ../nnbar-full-det-sim-source

# Profiling the GEANT4 simulation
valgrind --tool=callgrind ./nnbar-calo-sim -m ./macro/test.mac 
./gprof2dot.py -f callgrind callgrind.out.x | dot -Tsvg -o output.svg
