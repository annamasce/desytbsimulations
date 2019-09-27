# DesyTBSimulation

Geant4 simulation of the SciFi telescope employed in the DESY Testbeam in October 2019

$cd ~/elena_tests/examples/
$source setup.sh
$cd ~/elena_tests/examples/advanced/
$git clone https://gitlab.cern.ch/mstramag/desytbsimulation.git newname
$cd ~/elena_tests/examples/advanced/newname/
$mkdir builds
$cd builds/
$cmake ../
$make -j8
$./amsEcal mips.mac

$mv  Run_5000_withSatHit.root  Run_5XXX_withSatHit.root
$mv  mips.root  mips_5XXX.root

These two files are:
Run_5000_withSatHit.root = the data in our data format to be analysed by our analysis
mips.root = a summary of the simulation output

$if you want to see what is happening
$./amsEcal 
and in the command line at the bottom of the window /control/execute mips.mac
(the default is 50000 events, set in ~/elena_tests/examples/advanced/newname/mips.mac, if you want you can change it when you run the online version to have few events, you can change it there, but after you have to cmake ../ and make again)

The configuration of the beam is in :~/elena_tests/examples/advanced/newname/src/PrimaryGeneratorAction.cc lines 102-108.

The shifts of the planes are configurable in : ~/elena_tests/examples/advanced/newname/src/DetectorConstruction.cc lines 333-339
here you can as well change the different materials (Pb or Vacuum if you want the absorber or not) and dimensions (for W or Pb)


