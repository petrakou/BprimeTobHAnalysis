#include "../../BprimeTobHAnalysis/interface/histCummulative.h"

int main (int argc, char *argv[]) {

  TH1D* hgaus = new TH1D ("hgaus", "gaus", 1000, -50, 50) ; 
  hgaus->Sumw2() ; 
  hgaus->FillRandom("gaus", 10000) ; 
  char optUp[10] = "Up" ;
  TGraphAsymmErrors* gerfUp = getCummulative(hgaus, optUp) ; 
  char optDown[10] = "Down" ;
  TGraphAsymmErrors* gerfDown = getCummulative(hgaus, optDown) ; 

  TFile * fout = new TFile("file.root","RECREATE") ; 
  fout->cd() ; 
  hgaus->Write();
  gerfUp->Write();
  gerfDown->Write();
  fout->Close() ;

  return 0 ; 

}
