# Calorimeter simulation for NNbar/HIBEAM experiment


<h3>How to Build</h3>

<h5>Within the directory that contains nnbar-calo-sim/, run </h5>
  
$ mkdir nnbar-calo-sim-build <br>
$ cd nnbar-calo-sim-build  <br>
$ cmake -DGeant4_DIR=path_to_Geant4_installation/lib/Geant4-10.0.0/ ../nnbar-calo-sim <br>
$ make <br>              

<h5>If you have shared library problems, you made need to run </h5>

$ source path_to_Geant4_installation/bin/geant4.sh 

<h3>How to Run</h3>

<h5>Interactive Mode</h5>

$ ./nnbar-calo-sim

<h5>Macro Mode</h5>

$ ./nnbar-calo-sim -m nnbar-calo-sim.in
