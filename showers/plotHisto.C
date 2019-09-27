{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("run2.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");

  TH1D* hist[29];
  for(int i=1;i<30;i++){
  TString num="";
  num.Form("%d",i);
  hist[i] = (TH1D*)f.Get(num);
  hist[i]->Draw("HIST");
  }
  TH2D* hist2[5];
  for(int i=1;i<6;i++){
  TString num="";
  num.Form("%d",i);
  hist2[i] = (TH2D*)f.Get(i);
  hist2[i]->Draw("colz");
  }
}  
