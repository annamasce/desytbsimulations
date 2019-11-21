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
  Trkids.clear();
  Particleids.clear();
  Masses.clear();
  Pxs.clear();
  Pys.clear();
  Pzs.clear();
  Pts.clear();
  Thetas.clear();
  Phis.clear();
  Energies.clear();
  Momenta.clear();
//layers.clear();
TriggerIDE = 0;
TriggerCounterE = 0;
PlaneCodeE.clear();
PlaneNumberE.clear();
Board_IPE.clear(); 
Track_IDE.clear();
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
ParticleIDE.clear();
MassE.clear();
PxE.clear();
PyE.clear();
PzE.clear();
PtE.clear();
ThetaE.clear();
PhiE.clear();
EnergyE.clear();
MomentumE.clear();
steptimes.clear();
steptracks.clear();
stepparticleids.clear();
stepmasses.clear();
steppxs.clear();
steppys.clear();
steppzs.clear();
steppts.clear();
stepthetas.clear();
stepphis.clear();
stepenergies.clear();
stepmomenta.clear();
trigtimes.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumDeStep(G4int iModule, G4int iPlane, G4int iLayer, G4int iFiber, G4int iTrigger, G4double deStep, G4bool samefibre , G4String particle, G4double time, G4double trackid, G4double mass, G4double px, G4double py, G4double pz, G4double pt, G4double theta, G4double phi, G4double energy, G4double momentum)
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

  //KEY_HERE!
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
    steptracks.push_back(trackid);
    stepparticleids.push_back(particle);
    stepmasses.push_back(mass);
    steppxs.push_back(px);
    steppys.push_back(py);
    steppzs.push_back(pz);
    steppts.push_back(pt);
    stepthetas.push_back(theta);
    stepphis.push_back(phi);
    stepenergies.push_back(energy);
    stepmomenta.push_back(momentum);
    //HvisFiber[kFiber] +=1;
    //EvisFiber[kFiber] += deStep;
    if(!samefibre){
        if((particle!= "gamma" &&  HvisSum>0.1) || (particle== "gamma" && HvisSum>0.1)){
        HvisFiber[kFiber] +=1.; 
        EvisFiber[kFiber] += HvisSum;
        EvisLayerm[kLayer] += HvisSum;
        HvisLayerm[kLayer] += 1.;
        //Trkids[kFiber]=trackid;
        if(steptimes.size()>0) Hittimes[kFiber]=steptimes.at(0);
        if(steptracks.size()>0) Trkids[kFiber]=steptracks.at(0);
        if(stepparticleids.size()>0){
            Particleids[kFiber]=stepparticleids.at(0);
            G4cout << "stepparticleids.at(0) ----> " << stepparticleids.at(0) << G4endl;
        }
        if(stepmasses.size()>0) Masses[kFiber]=stepmasses.at(0);
        if(steppxs.size()>0) Pxs[kFiber]=steppxs.at(0);
        if(steppys.size()>0) Pys[kFiber]=steppys.at(0);
        if(steppzs.size()>0) Pzs[kFiber]=steppzs.at(0);
        if(steppts.size()>0) Pts[kFiber]=steppts.at(0);
        if(stepthetas.size()>0) Thetas[kFiber]=stepthetas.at(0);
        if(stepphis.size()>0) Phis[kFiber]=stepphis.at(0);
        if(stepenergies.size()>0) Energies[kFiber]=stepenergies.at(0);
        if(stepmomenta.size()>0) Momenta[kFiber]=stepmomenta.at(0);
        }
        HvisSum = 0.;
        steptimes.clear();
        steptracks.clear();
        stepparticleids.clear();
        stepmasses.clear();
        steppxs.clear();
        steppys.clear();
        steppzs.clear();
        steppts.clear();
        stepthetas.clear();
        stepphis.clear();
        stepenergies.clear();
        stepmomenta.clear();
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
  //G4int kFiber = thit->first;
  //G4int iFiber = kFiber;
  G4double tim = thit->second;
  //G4int f=(kFiber-1)%512;
        Hit_RealTimeE.push_back(tim);
        double timetmp=tim*0.68-(int)(tim*0.68);
        unsigned long long hittmp;
        if(timetmp>0) hittmp= (unsigned long long)(tim*0.68);
        else hittmp= (unsigned long long)(tim*0.68)-1;
        Hit_TimeE.push_back(hittmp);
        Fine_TimeE.push_back((unsigned long long)((tim*0.68-hittmp)*32));
}
std::map<G4int,G4double>::iterator trid;
for (trid = Trkids.begin(); trid != Trkids.end(); trid++) {
  //G4int kFiber = trid->first;
  G4double tid = trid->second;
  Track_IDE.push_back(tid);
}

std::map<G4int,G4double>::iterator mhit;
for (mhit = Masses.begin(); mhit != Masses.end(); mhit++) {
  //G4int kFiber = mhit->first;
  G4double hitmass = mhit->second;
  MassE.push_back(hitmass);
}

std::map<G4int,G4String>::iterator parthit;
for (parthit = Particleids.begin(); parthit != Particleids.end(); parthit++) {
  //G4int kFiber = parthit->first;
  G4String hitparticle = parthit->second;
  G4cout << "hitparticle ----> " << hitparticle << G4endl;
  ParticleIDE.push_back(hitparticle);
}

std::map<G4int,G4double>::iterator pxhit;
for (pxhit = Pxs.begin(); pxhit != Pxs.end(); pxhit++) {
  //G4int kFiber = pxhit->first;
  G4double hitpx = pxhit->second;
  PxE.push_back(hitpx);
}

std::map<G4int,G4double>::iterator pyhit;
for (pyhit = Pys.begin(); pyhit != Pys.end(); pyhit++) {
  //G4int kFiber = pyhit->first;
  G4double hitpy = pyhit->second;
  PyE.push_back(hitpy);
}

std::map<G4int,G4double>::iterator pzhit;
for (pzhit = Pzs.begin(); pzhit != Pzs.end(); pzhit++) {
  //G4int kFiber = pzhit->first;
  G4double hitpz = pzhit->second;
  PzE.push_back(hitpz);
}

std::map<G4int,G4double>::iterator pthit;
for (pthit = Pts.begin(); pthit != Pts.end(); pthit++) {
  //G4int kFiber = pthit->first;
  G4double hitpt = pthit->second;
  PtE.push_back(hitpt);
}

std::map<G4int,G4double>::iterator thetahit;
for (thetahit = Thetas.begin(); thetahit != Thetas.end(); thetahit++) {
  //G4int kFiber = thetahit->first;
  G4double hittheta = thetahit->second;
  ThetaE.push_back(hittheta);
}

std::map<G4int,G4double>::iterator phihit;
for (phihit = Phis.begin(); phihit != Phis.end(); phihit++) {
  //G4int kFiber = phihit->first;
  G4double hitphi = phihit->second;
  PhiE.push_back(hitphi);
}

std::map<G4int,G4double>::iterator enhit;
for (enhit = Energies.begin(); enhit != Energies.end(); enhit++) {
  //G4int kFiber = enhit->first;
  G4double hitenergy = enhit->second;
  EnergyE.push_back(hitenergy);
}

std::map<G4int,G4double>::iterator momhit;
for (momhit = Momenta.begin(); momhit != Momenta.end(); momhit++) {
  //G4int kFiber = momhit->first;
  G4double hitmomentum = momhit->second;
  MomentumE.push_back(hitmomentum);
}

//G4cout<<Hittimes.size()<<" "<<Trkids.size()<<G4endl;
     //analysisManager->FillNtupleDColumn(1, 0, layers);
     //analysisManager->AddNtupleRow();

     //KEY_HERE!
     fHistoManager->FillNtuple(TriggerIDE, TriggerCounterE, PlaneCodeE, PlaneNumberE, Board_IPE, Board_IDE, STiC_IDE, Ch_IDE, Ch_PositionE, AmpE, Hit_TimeE, Fine_TimeE, Trig_TimeE, Trig_RealTimeE, Hit_RealTimeE, Track_IDE, ParticleIDE, MassE, PxE, PyE, PzE, PtE, ThetaE, PhiE, EnergyE, MomentumE);
    
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

