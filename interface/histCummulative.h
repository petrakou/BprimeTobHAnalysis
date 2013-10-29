#include <TH1F.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <cmath> 

using namespace std ; 

template <class Type> 
TGraphAsymmErrors* getCummulative (Type histp) {

  Type hcummulp = (Type)histp->Clone("hcummul") ; 
  hcummulp->Sumw2() ; 
  hcummulp->Reset(0) ; 

  double val(0) ; 
  double err(0) ; 
  for (int ibin = 0; ibin <= histp->GetNbinsX()+1; ++ibin) {
    val += histp->GetBinContent(ibin) ; 
    err = (err*err) + (histp->GetBinError(ibin)*histp->GetBinError(ibin)) ; 
    err = sqrt(err) ; 
    hcummulp->SetBinContent(ibin, val) ; 
    hcummulp->SetBinError(ibin, err) ; 
  }

  Type htot = (Type)histp->Clone("htot") ; 
  htot->Sumw2() ; 
  for (int ibin = 0; ibin <= histp->GetNbinsX()+1; ++ibin) {
    htot->SetBinContent(ibin, hcummulp->GetBinContent(hcummulp->GetNbinsX()+1)) ; 
    htot->SetBinError(ibin, hcummulp->GetBinError(hcummulp->GetNbinsX()+1)) ; 
  }

  TGraphAsymmErrors* geff = new TGraphAsymmErrors(hcummulp, htot, "cl=0.683 b(1,1) mode") ; 
  geff->SetName("geff") ; 

  return geff ; 
}


