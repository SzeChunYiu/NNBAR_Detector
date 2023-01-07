# Beampipe + Detector simulation for NNbar/HIBEAM experiment

<h3>Software Dependencies</h3>
- GEANT4.10.07.p02 or older
- ROOT (Version TBA)
- boost library

<h3>How to Build</h3>

<h5>Within the directory that contains nnbar-calo-sim/, run </h5>
  
$ mkdir nnbar-full-det <br>
$ cd nnbar-full-det  <br>
$ sudo cmake -DGeant4_DIR=**\<path to GEANT4 install folder\>**/lib/Geant4-10.7.2 -DVERSION_SHIELD=0 -DCRY_BUILD=0 -DVERSION_MCPL=0 **\<path to source\>** <br>
$ sudo make install -j4 <br>              

<h3> Settings of the simulation </h3>
1. You can set the particle source by -DCRY_BUILD and -DVERSION_MCPL <br>
If -DCRY_BUILD=1 (while we must set -DVERSION_MCPL=0),particles will be generated from CRY library.<br>
If -DVERSION_MCPL=1 (while we must set -DCRY_BUILD=0),particles will be generated from a specified MCPL file. <br>

2. You can switch the cosmic ray shielding geometry on and off by -DVERSION_SHIELD <br>
If -DVERSION_SHIELD=0 (by default), the shield will be turned off <br>
If -DVERSION_SHIELD=1, the shield will be constructed in the simulation <br>

