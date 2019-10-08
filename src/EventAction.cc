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

#include "EventAction.hh"

#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det, HistoManager* histo, PrimaryGeneratorAction* prim)
:detector(det), fHistoManager(histo), primary(prim)
{
  nbOfModules = detector->GetNbModules();	 	
  nbOfLayers  = detector->GetNbLayers();
  nbOfplanes  = detector->GetNbPlanes();
  nbOffibres  = detector->GetNbFibers();
  kLayerMax = nbOfModules*nbOfLayers*nbOfplanes + 1;
  //kLayerMax = nbOfModules*nbOfLayers+1;
  EtotCalor = EvisCalor = 0.;
  //TriggerIDE=numev;
  //TriggerCounterE=numev;
     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  EtotLayer.resize(kLayerMax);
  EvisLayer.resize(kLayerMax);			
  for (G4int k=0; k<kLayerMax; k++) {
    EtotLayer[k] = EvisLayer[k] = 0.0;
  }
  EtotCalor = EvisCalor = HvisSum = 0.;
  EvisFiber.clear();
  HvisFiber.clear();
  EvisLayerm.clear();
  HvisLayerm.clear();
  Hittimes.clear();
//layers.clear();
TriggerIDE = 0;
TriggerCounterE = 0;
PlaneCodeE.clear();
PlaneNumberE.clear();
Board_IPE.clear(); 
Board_IDE.clear(); 
STiC_IDE.clear();  
Ch_IDE.clear();    
Ch_PositionE.clear();
AmpE.clear(); 
Hit_TimeE.clear(); 
Fine_TimeE.clear();
Trig_TimeE.clear();
Trig_RealTimeE.clear(); 
Hit_RealTimeE.clear(); 
steptimes.clear();
trigtimes.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumDeStep(G4int iModule, G4int iPlane, G4int iLayer, G4int iFiber, G4int iTrigger, G4double deStep, G4bool samefibre , G4String particle, G4double time)
{
  if(iTrigger>0){trigtimes.push_back(time);}
  if (iModule > 0) EtotCalor += deStep;
  if(iModule==4 || iModule==5)iModule=iModule-1;
  if(iModule==7 || iModule==8)iModule=iModule-2;
  G4int kLayer = 0; G4int kFiber = 0;
  if (iLayer > 0) {
    //kLayer = (iModule-1)*nbOfLayers + iLayer; 
    //G4cout<<"Module, Plane, Layer "<<iModule<<" "<<iPlane<<" "<<iLayer<<G4endl;
    kLayer = (iModule-1)*nbOfLayers*nbOfplanes+ (iPlane-1)*nbOfLayers + iLayer;
	//G4cout<<"kLayer "<<kLayer<<G4endl;
    EtotLayer[kLayer] += deStep;
  }
  
  if (iFiber > 0) {
    //G4cout<<"Fibre >0"<<G4endl;
	EvisLayer[kLayer] += deStep;
    //G4cout<<"stepping"<<G4endl;
	EvisCalor += deStep;
    //kFiber = 1000*kLayer + iFiber;     
    //G4cout<<"Fibre "<<iFiber<<G4endl;
    kFiber = (iModule-1)*nbOfLayers*nbOfplanes*nbOffibres+ (iPlane-1)*nbOfLayers*nbOffibres+(iLayer-1)*nbOffibres+iFiber;
	//G4cout<<"kFibre "<<kFiber<<" "<<deStep<<G4endl;
    //EvisFiber[kFiber] += deStep;
    HvisSum+=deStep;
    steptimes.push_back(time);
    //HvisFiber[kFiber] +=1;
    //EvisFiber[kFiber] += deStep;
    if(!samefibre){
        if((particle!= "gamma" &&  HvisSum>0.1) || (particle== "gamma" && HvisSum>0.1)){
        HvisFiber[kFiber] +=1.; 
        EvisFiber[kFiber] += HvisSum;
        EvisLayerm[kLayer] += HvisSum;
        HvisLayerm[kLayer] += 1.;
        if(steptimes.size()>0)Hittimes[kFiber]=steptimes.at(0);
        }
        HvisSum = 0.;
        steptimes.clear();
        }
  }
  else{HvisSum = 0.;}
//}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  //TString planame[12]={"X1","Y1","X2","Y2","X3","Y3","X4","Y4","X5","Y5","X6","Y6"};
  //pass informations to Run
  //
  //G4cout<<"EndOfEventAction start"<<G4endl;
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  	
  for (G4int k=0; k<kLayerMax; k++) {
     run->SumEvents_1(k,EtotLayer[k],EvisLayer[k]);   
  }
  TriggerIDE=run->GetNev();
  TriggerCounterE=run->GetNev();
    
  for(int l=0;l<kLayerMax-1;l++){
  if(trigtimes.size()>0){
  Trig_RealTimeE.push_back(trigtimes.at(0));
  Trig_TimeE.push_back((unsigned long long)(trigtimes.at(0)*0.17/4));}
  else{
  Trig_RealTimeE.push_back(0);
  Trig_TimeE.push_back(0);
  }
  }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
  analysisManager->FillH1(1,EtotCalor);
  analysisManager->FillH1(2,EvisCalor);
    
  G4double Ebeam = primary->GetParticleGun()->GetParticleEnergy();
  G4double Eleak = Ebeam - EtotCalor;
  run->SumEvents_2(EtotCalor,EvisCalor,Eleak);
  

  std::map<G4int,G4double>::iterator lit;
  for (lit = EvisLayerm.begin(); lit != EvisLayerm.end(); lit++) {
  G4int kLayer = lit->first;
  G4double Evis = lit->second;
  analysisManager->FillH2(4,kLayer+0.5,Evis);
  }
  std::map<G4int,G4double>::iterator hlit;
  for (hlit = HvisLayerm.begin(); hlit != HvisLayerm.end(); hlit++) {
  G4int kLayer = hlit->first;
  G4double Hvis = hlit->second;
  analysisManager->FillH2(2,kLayer+0.5,Hvis);
  }



  std::map<G4int,G4double>::iterator it;         
  for (it = EvisFiber.begin(); it != EvisFiber.end(); it++) {
     G4int kFiber = it->first;
	 G4int iFiber = kFiber;
     G4double Evis = it->second;
	 analysisManager->FillH1(5,iFiber+0.5,Evis);
     G4int h=(int)((kFiber-1)/512)+6;
     G4int f=(kFiber-1)%512;
     analysisManager->FillH1(h,f+0.5,Evis);
     analysisManager->FillH2(3,iFiber+0.5,Evis);
     unsigned long long tot=(unsigned long long)(Evis*255/20.);
     AmpE.push_back(tot);
}
  std::map<G4int,G4double>::iterator hit;
  for (hit = HvisFiber.begin(); hit != HvisFiber.end(); hit++) {
     G4int kFiber = hit->first;
	 G4int iFiber = kFiber;
     G4double Hvis = hit->second;
     G4int h=(int)((kFiber-1)/512)+18;
     G4int f=(kFiber-1)%512;
     analysisManager->FillH1(h,f+0.5,Hvis);
     analysisManager->FillH2(5,f+0.5,h-17,Hvis);
     analysisManager->FillH2(1,iFiber+0.5,Hvis);
//     layers.push_back(h-17);
     PlaneCodeE.push_back((h-17)%2);
     //PlaneNumberE.push_back(h-17);
     PlaneNumberE.push_back((int)((h-18)/2)+1);
     Board_IPE.push_back(h-17);
     Board_IDE.push_back(h-17);
     STiC_IDE.push_back((int)((f)/128)*2+((f)%2));
     Ch_IDE.push_back(((f-128*((int)(f/128)))-((f-128*((int)(f/128)))%2))/2);
     Ch_PositionE.push_back(f);
     //Hit_TimeE.push_back(0.);
     //Fine_TimeE.push_back(0.);
     //Hit_RealTimeE.push_back();
     //analysisManager->FillNtupleDColumn(1, 0, h-17);
}
std::map<G4int,G4double>::iterator thit;
for (thit = Hittimes.begin(); thit != Hittimes.end(); thit++) {
  G4int kFiber = thit->first;
  G4int iFiber = kFiber;
  G4double tim = thit->second;
  G4int f=(kFiber-1)%512;
        Hit_RealTimeE.push_back(tim);
        double timetmp=tim*0.68-(int)(tim*0.68);
        unsigned long long hittmp;
        if(timetmp>0) hittmp= (unsigned long long)(tim*0.68);
        else hittmp= (unsigned long long)(tim*0.68)-1;
        Hit_TimeE.push_back(hittmp);
        Fine_TimeE.push_back((unsigned long long)((tim*0.68-hittmp)*32));
}
     //analysisManager->FillNtupleDColumn(1, 0, layers);
     //analysisManager->AddNtupleRow();
     fHistoManager->FillNtuple(TriggerIDE, TriggerCounterE, PlaneCodeE, PlaneNumberE, Board_IPE, Board_IDE, STiC_IDE, Ch_IDE, Ch_PositionE, AmpE, Hit_TimeE, Fine_TimeE, Trig_TimeE, Trig_RealTimeE, Hit_RealTimeE);
    
  //write fired fibers on a file
  //
  //// WriteFibers(evt); 
  //G4cout<<"EndOfEventAction ends"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include <fstream>

void EventAction::WriteFibers(const G4Event* evt)
{
  // event is appended on a file
  //
  //G4cout<<"WriteFibers starts"<<G4endl;

  G4String name = G4AnalysisManager::Instance()->GetFileName();
  G4String fileName = name + ".fibers.ascii";
  
  std::ofstream File(fileName, std::ios::app);
  std::ios::fmtflags mode = File.flags();  
  File.setf( std::ios::scientific, std::ios::floatfield );
  G4int prec = File.precision(3);
    
  //write event number  
  //
  File << evt->GetEventID() << G4endl;
  
  //gun particle informations
  //
  G4ParticleGun* gun = primary->GetParticleGun();
  G4double ekin = gun->GetParticleEnergy();
  G4ThreeVector direction = gun->GetParticleMomentumDirection();
  G4ThreeVector position  = gun->GetParticlePosition();
  File << ekin << " " << direction << " " << position << G4endl;  
  
  //write fibers
  //
  File << EvisFiber.size() << G4endl;
  //
  std::map<G4int,G4double>::iterator it;         
  for (it = EvisFiber.begin(); it != EvisFiber.end(); it++) {
     G4int kFiber = it->first;
     G4double Evis = it->second;
     File << " " << std::setw(7) << kFiber << " "<< std::setw(10) << Evis
            << G4endl;
  }
           
  File << G4endl;
    
  // restaure default formats
  File.setf(mode,std::ios::floatfield);
  File.precision(prec);         
  //G4cout<<"WriteFibers ends"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

