//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "DetectorConstruction.hh"
#include "globals.hh"

#include <vector>
#include <map>

class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:  
    EventAction(DetectorConstruction*,HistoManager*, PrimaryGeneratorAction*);
   ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void SumDeStep(G4int, G4int,G4int, G4int, G4int, G4double,  G4bool, G4String, G4double, G4double);
	
	void WriteFibers(const G4Event*);
			         	    

  private:  
    DetectorConstruction*   detector;
    PrimaryGeneratorAction* primary;
	
	G4int nbOfplanes, nbOfModules, nbOfLayers,nbOffibres, kLayerMax;     
    std::vector<G4double>   EtotLayer;
    std::vector<G4double>   EvisLayer;
    unsigned long long TriggerIDE;
    unsigned long long TriggerCounterE;

    std::vector<uint32_t> PlaneCodeE;
    std::vector<uint32_t> PlaneNumberE;
    std::vector<uint32_t> Board_IPE;
    std::vector<uint32_t> Board_IDE;
    std::vector<uint32_t> STiC_IDE;
    std::vector<uint32_t> Ch_IDE;
    std::vector<uint32_t> Ch_PositionE; // Position on a layer

    std::vector<unsigned long long> AmpE;
    std::vector<unsigned long long> Hit_TimeE;
    std::vector<unsigned long long> Fine_TimeE;
    std::vector<unsigned long long> Trig_TimeE;

    std::vector<double> trigtimes;
    std::vector<double> steptimes;
    std::vector<double> steptracks;
    std::vector<double> Trig_RealTimeE;
    std::vector<double> Hit_RealTimeE; 
	std::vector<double> Track_IDE;
    G4double EtotCalor;
	G4double EvisCalor;
    G4double HvisSum;
    HistoManager*           fHistoManager; 	
	std::map<G4int, G4double> EvisFiber;
    std::map<G4int, G4double> HvisFiber;
    std::map<G4int, G4double> EvisLayerm;
    std::map<G4int, G4double> HvisLayerm;
    std::map<G4int, G4double> Hittimes;
    std::map<G4int, G4double> Trkids;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
