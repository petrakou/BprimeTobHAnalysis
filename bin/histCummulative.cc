#include "BprimeTobHAnalysisv1/BprimeTobHAnalysis/interface/histCummulative.h"

int main (int argc, char *argv[]) {

  TH1D* hgaus = new TH1D ("hgaus", "gaus", 1000, -50, 50) ; 
  hgaus->Sumw2() ; 
  hgaus->FillRandom("gaus", 10000) ; 
  TGraphAsymmErrors* gerf = getCummulative(hgaus) ; 

  TFile * fout = new TFile("file.root","RECREATE") ; 
  fout->cd() ; 
  hgaus->Write();
  gerf->Write();
  fout->Close() ;

  return 0 ; 

}
