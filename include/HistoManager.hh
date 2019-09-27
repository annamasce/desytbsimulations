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
/// \file electromagnetic/TestEm1/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include <vector>
#include "g4root.hh"
////#include "g4xml.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class TTree;   
class TFile;
class HistoManager
{
  public:
   HistoManager();
  ~HistoManager();
   void Book();
   void WriteTuple();
   void FillNtuple(unsigned long long TriggerIDE, unsigned long long TriggerCounterE,  std::vector<uint32_t> PlaneCodeE, std::vector<uint32_t> PlaneNumberE, std::vector<uint32_t> Board_IPE, std::vector<uint32_t> Board_IDE, std::vector<uint32_t> STiC_IDE, std::vector<uint32_t> Ch_IDE, std::vector<uint32_t> Ch_PositionE, std::vector<unsigned long long> AmpE, std::vector<unsigned long long> Hit_TimeE, std::vector<unsigned long long> Fine_TimeE, std::vector<unsigned long long> Trig_TimeE, std::vector<double> Trig_RealTimeE, std::vector<double> Hit_RealTimeE); 

  private:
    TFile*   fRootfile;
    TTree*   fNtuple1;
    G4String fFileName;
    unsigned long long TriggerID;
    unsigned long long TriggerCounter;
    std::vector<uint32_t> PlaneCode;
    std::vector<uint32_t> PlaneNumber;
    std::vector<uint32_t> Board_IP; 
    std::vector<uint32_t> Board_ID; 
    std::vector<uint32_t> STiC_ID;  
    std::vector<uint32_t> Ch_ID;    
    std::vector<uint32_t> Ch_Position;
    std::vector<unsigned long long> Amp; 
    std::vector<unsigned long long> Hit_Time; 
    std::vector<unsigned long long> Fine_Time;
    std::vector<unsigned long long> Trig_Time;
    std::vector<double> Trig_RealTime; 
    std::vector<double> Hit_RealTime; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

