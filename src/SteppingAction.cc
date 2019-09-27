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

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
:G4UserSteppingAction(),detector(det),eventAct(evt)
{
  first = true;
  lvol_world = lvol_module = lvol_plane = lvol_layer = lvol_trigger = lvol_fiber = 0;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step )
{ 
 //some initialisation
 // 
 if (first) {
   // G4cout<<"stepping first"<<G4endl;
   lvol_world  = detector->GetLvolWorld();
   lvol_module = detector->GetLvolModule();
   lvol_plane  = detector->GetLvolPlane();
   lvol_layer  = detector->GetLvolLayer();
   lvol_fiber  = detector->GetLvolFiber();     
   lvol_trigger= detector->GetLvolTrigger();
   first = false;   
 }
 
 //if no edep, return
 //
 G4double edep = step->GetTotalEnergyDeposit();
 ///if (edep == 0.) return;
 
 //locate point in geometry
 //
 G4int iModule = 0;
 G4int iPlane  = 0;
 G4int iLayer  = 0;
 G4int iFiber  = 0;
 G4int iModule2 = 0;
 G4int iPlane2  = 0;
 G4int iLayer2  = 0;
 G4int iFiber2  = 0;
 G4int iTrigger  = 0;
 G4bool samefibre =0;
   
 G4TouchableHandle touch1 = step->GetPreStepPoint()->GetTouchableHandle(); 
 G4LogicalVolume* lvol = touch1->GetVolume()->GetLogicalVolume();
 G4Track* track = step->GetTrack();
 if (!track)return;
 if (track== NULL)return;
 G4String particle = track->GetDefinition()->GetParticleName();
//G4cout<<"particle: "<<track->GetDefinition()->GetParticleName()<<G4endl;
 if (!lvol)return;
 else if(lvol== NULL)return;
 else if (lvol == lvol_world) return;
 else if (lvol == lvol_module) { iModule = touch1->GetCopyNumber(0);}
 else if (lvol == lvol_plane)  { iPlane  = touch1->GetCopyNumber(0);
                                 iModule = touch1->GetCopyNumber(1);}
 else if (lvol == lvol_layer)  { iLayer  = touch1->GetCopyNumber(0);
                                 iPlane  = touch1->GetCopyNumber(1);
                                 iModule = touch1->GetCopyNumber(2);}
 else if (lvol == lvol_fiber)  { iFiber  = touch1->GetCopyNumber(0);
                                 iLayer  = touch1->GetCopyNumber(1);
                                 iPlane  = touch1->GetCopyNumber(2);
                                 iModule = touch1->GetCopyNumber(3);}
 else if (lvol == lvol_trigger){ iTrigger  =1;}
 if(!step->GetPostStepPoint())return;
 else if(step->GetPostStepPoint()==NULL)return;
 
 G4TouchableHandle touch2 = step->GetPostStepPoint()->GetTouchableHandle(); 
 if(!touch2)return;
 if(!touch2->GetVolume())return;
 if(touch2->GetVolume()==NULL)return;
 G4LogicalVolume* lvol2 = touch2->GetVolume()->GetLogicalVolume();
 //G4cout<<"particle: "<<track->GetDefinition()->GetParticleName()<<G4endl;
 if (!lvol2)return;
 else if(lvol2== NULL)return;
 if (lvol2 == lvol_module) { iModule2 = touch2->GetCopyNumber(0);}
 else if (lvol2 == lvol_plane)  { iPlane2  = touch2->GetCopyNumber(0);
                                 iModule2 = touch2->GetCopyNumber(1);}
 else if (lvol2 == lvol_layer)  { iLayer2  = touch2->GetCopyNumber(0);
                                 iPlane2  = touch2->GetCopyNumber(1);
                                 iModule2 = touch2->GetCopyNumber(2);}
 else if (lvol2 == lvol_fiber)  { iFiber2  = touch2->GetCopyNumber(0);
                                 iLayer2  = touch2->GetCopyNumber(1);
                                 iPlane2  = touch2->GetCopyNumber(2);
                                 iModule2 = touch2->GetCopyNumber(3);}
 G4double time=0;
 if(edep>0){time= step->GetPreStepPoint()->GetGlobalTime();}
 // sum edep
 //
 //G4cout<<"Module, Plane, Layer, Fibre "<<iModule<<" "<<iPlane<<" "<<iLayer<<" "<<iFiber<<G4endl;
 //G4cout<<"Module2, Plane2, Layer2, Fibre2 "<<iModule2<<" "<<iPlane2<<" "<<iLayer2<<" "<<iFiber2<<G4endl;
 if(iModule==iModule2 && iPlane==iPlane2 && iLayer ==iLayer2 && iFiber==iFiber2)samefibre=1;
 eventAct->SumDeStep(iModule, iPlane, iLayer, iFiber, iTrigger, edep, samefibre, particle,time);         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SteppingAction::BirksAttenuation(const G4Step* aStep)
{
 //Example of Birk attenuation law in organic scintillators.
 //adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
 //
 G4Material* material = aStep->GetTrack()->GetMaterial();
 G4double birk1       = material->GetIonisation()->GetBirksConstant();
 G4double destep      = aStep->GetTotalEnergyDeposit();
 G4double stepl       = aStep->GetStepLength();  
 G4double charge      = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
 //
 G4double response = destep;
 if (birk1*destep*stepl*charge != 0.)
   {
     response = destep/(1. + birk1*destep/stepl);
   }
 return response;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

