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
/// \file electromagnetic/TestEm1/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include <TTree.h>
#include <TFile.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("amsEcal")
{
fRootfile = 0;
fNtuple1 = 0;
TriggerID = 0;
TriggerCounter = 0;
PlaneCode.clear();
PlaneNumber.clear();
Board_IP.clear(); 
Board_ID.clear(); 
STiC_ID.clear();  
Ch_ID.clear();    
Ch_Position.clear();
Amp.clear(); 
Hit_Time.clear(); 
Fine_Time.clear();
Trig_Time.clear();
Trig_RealTime.clear(); 
Hit_RealTime.clear(); 
Track_ID.clear();

Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
if(fRootfile) delete fRootfile;
delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  //fRootFile = new TFile("SPILL_DATA","RECREATE");
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);     // start histogram numbering from 1
  analysisManager->SetActivation(true);    // enable inactivation of histograms
TriggerID = 0;
TriggerCounter = 0;
PlaneCode.clear();
PlaneNumber.clear();
Board_IP.clear(); 
Board_ID.clear(); 
STiC_ID.clear();  
Ch_ID.clear();    
Ch_Position.clear();
Amp.clear(); 
Hit_Time.clear(); 
Fine_Time.clear();
Trig_Time.clear();
Trig_RealTime.clear(); 
Hit_RealTime.clear(); 
Track_ID.clear();
 // Define histograms start values
  const G4int kMaxHisto = 34;
  const G4String id[] = {"1", "2", "3" , "4", "5", "6", "7"," 8", "9",  "10", "11" , "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23" , "24", "25", "26", "27", "28", "29", "1", "2", "3", "4", "5"};
  const G4String title[] = 
				{ "total Etot in Ecal",       //1
                  "total Evis in Ecal",       //2
                  "Etot profile",             //3
                  "Evis profile",             //4
                  "Evis per fiber",           //5
                  "Energy per channel in layer 1", //6
                  "Energy per channel in layer 2", //7
                  "Energy per channel in layer 3", //8
                  "Energy per channel in layer 4", //9
                  "Energy per channel in layer 5", //10
                  "Energy per channel in layer 6", //11
                  "Energy per channel in layer 7", //12
                  "Energy per channel in layer 8", //13
                  "Energy per channel in layer 9", //14
                  "Energy per channel in layer 10", //15
                  "Energy per channel in layer 11", //16
                  "Energy per channel in layer 12",  //17
                  "Hit per channel in layer 1", //18
                  "Hit per channel in layer 2", //19
                  "Hit per channel in layer 3", //20
                  "Hit per channel in layer 4", //21
                  "Hit per channel in layer 5", //22
                  "Hit per channel in layer 6", //23
                  "Hit per channel in layer 7", //24
                  "Hit per channel in layer 8", //25
                  "Hit per channel in layer 9", //26
                  "Hit per channel in layer 10", //27
                  "Hit per channel in layer 11", //28
                  "Hit per channel in layer 12", //29
                  "Hit distribution fibre",         //30
                  "Hit distribution layer",   //31
                  "Energy distribution in fibres",      //32
                  "Energy distribution in layers",      //33
                  "Hits distribution in layers"       //34
                 };
				 
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto-5; k++) {
  G4int ih = analysisManager->CreateH1(title[k], title[k], nbins, vmin, vmax); analysisManager->SetH1Activation(ih, false);
}
for (G4int k=0; k<5; k++){
     //G4cout<<"histo num: "<<k<<" is: "<<title[kMaxHisto-4+k]<<G4endl;
     G4int ih = analysisManager->CreateH2(title[kMaxHisto-5+k], title[kMaxHisto-5+k], nbins, vmin, vmax,nbins, vmin, vmax);
     analysisManager->SetH2Activation(ih, false);
     //G4cout<<"histo num: "<<ih<<" is: "<<title[kMaxHisto-4+k]<<G4endl;
     }
fNtuple1 = new TTree("data", "data tree");
    fNtuple1->Branch("Eventnum" ,&TriggerCounter);
    fNtuple1->Branch("TriggerID"      , &TriggerID); 
    fNtuple1->Branch("PlaneCode"  ,&PlaneCode);
    fNtuple1->Branch("PlaneNumber",&PlaneNumber);
    fNtuple1->Branch("BoardIP"    ,&Board_IP);
    fNtuple1->Branch("BoardID"    ,&Board_ID); 
    fNtuple1->Branch("STiCID"     ,&STiC_ID);
    fNtuple1->Branch("ChID"       ,&Ch_ID);
    fNtuple1->Branch("ChPosition" ,&Ch_Position);  

    fNtuple1->Branch("TOT",&Amp);
    fNtuple1->Branch("HitTime",&Hit_Time);
    fNtuple1->Branch("FineTime",&Fine_Time);
    fNtuple1->Branch("TriggerTime",&Trig_Time);

    fNtuple1->Branch("HitRealTime",&Hit_RealTime);
    fNtuple1->Branch("TriggerRealTime",&Trig_RealTime);
    fNtuple1->Branch("TrackID",&Track_ID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(unsigned long long TriggerIDE, unsigned long long TriggerCounterE, std::vector<uint32_t> PlaneCodeE,      std::vector<uint32_t> PlaneNumberE, std::vector<uint32_t> Board_IPE, std::vector<uint32_t> Board_IDE, std::vector<uint32_t> STiC_IDE, std::      vector<uint32_t> Ch_IDE, std::vector<uint32_t> Ch_PositionE, std::vector<unsigned long long> AmpE, std::vector<unsigned long long>               Hit_TimeE, std::vector<unsigned long long> Fine_TimeE, std::vector<unsigned long long> Trig_TimeE, std::vector<double> Trig_RealTimeE, std::vector<double> Hit_RealTimeE, std::vector<double> Track_IDE)
{
  TriggerID     = TriggerIDE;
  TriggerCounter= TriggerCounterE;
  PlaneCode     = PlaneCodeE;
  PlaneNumber   = PlaneNumberE;
  Board_IP      = Board_IPE;
  Board_ID      = Board_IDE;
  STiC_ID       = STiC_IDE;
  Ch_ID         = Ch_IDE;
  Ch_Position   = Ch_PositionE;
  Amp           = AmpE;
  Hit_Time      = Hit_TimeE;
  Fine_Time     = Fine_TimeE;
  Trig_Time     = Trig_TimeE;
  Trig_RealTime = Trig_RealTimeE;
  Hit_RealTime  = Hit_RealTimeE;
  Track_ID = Track_IDE;
  if (fNtuple1) fNtuple1->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::WriteTuple()
{
  fRootfile = new TFile("Run_5000_withSatHit.root","RECREATE");
  fNtuple1->Write();       // Writing the tuple to the file
  fRootfile->Write();
  fRootfile->Close();
  return; 
}


