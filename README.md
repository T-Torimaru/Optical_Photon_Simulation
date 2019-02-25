version : 10.00

@ kekcc

## Procedure from cloning to creating a root file
1. git clone

``` 
$ git clone https://github.com/T-Torimaru/Optical_Photon_Simulation.git
```
2. execute cmake and make

```
$ pwd
   ~/yourworkingdirectory
$ ls
   GEANT_Optical_Photon
$ mkdir build
$ cd build
$ cmake -DGeant4_DIR=/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-11/geant4/10.01/lib64/Geant4-10.1.0 ~/Optical_Photon_Simulation
$ make -j4
```

3. execute the file with macro

```
$ ./OpNovice -m ../Optical_Photon_Simulation/optPhoton.mac
```

## Contents

* optPhoton.mac in Optical_Photon_Simulation directory, not in build directory

  * You can choose energy spectrum of Sr90 or a fixed energy beam and define an initial state of events

* OpNoviceDetectorConstruction

  * geometry
  * LY : Scintillation Efficiency
    * Higher LY values make calculation time be longer because each of scintillation light tracks is calculated.
  * PDE of MPPC
  * surface state
  * visualization attributes (color...)
  * have a sensitivity to MPPC

* OpNoviceEventAction

  * Energy deposits at the scintillator and the trigger counter are filled into a root file

* OpNovicePhysicsList

  * Define physics processes

* OpNovicePrimaryGeneratorAction

  * Define initial states of Events, for instance, beam energy, the number of events, incident position, and etc.
  * optPhoton.mac is substitute for it.

* OpNoviceRunAction

  * Create root file and branches
  * Start and stop a run

* OpNoviceStackingAction

  * count excited scintillation lights

* OpNoviceSteppingAction

  * You can get information of each of steps
    * momentum, position, energy, step length, and etc...