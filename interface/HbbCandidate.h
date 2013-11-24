#ifndef BPRIMETOBHANALYSIS_INTERFACE_HBBCANDIDATE_H
#define BPRIMETOBHANALYSIS_INTERFACE_HBBCANDIDATE_H

#include "../../BprimeTobHAnalysis/interface/Particle.h" 

class HbbCandidate {

  public:

  HbbCandidate (const TLorentzVector& p4Higgs, const TLorentzVector& p4b0, const TLorentzVector& p4b1) ;
  ~HbbCandidate () {}; 

  TLorentzVector p4Higgs() { return Higgs.p4() ; } ; 
  TLorentzVector p4Dau0() { return Dau0.p4() ; } ; 
  TLorentzVector p4Dau1() { return Dau1.p4() ; } ; 

  private:

  Particle Higgs ; 
  Particle Dau0 ; 
  Particle Dau1 ; 

} ; 

#endif 

#ifdef BPRIMETOBHANALYSIS_INTERFACE_HBBCANDIDATE_H
HbbCandidate::HbbCandidate (const TLorentzVector& p4Higgs, const TLorentzVector& p4b0, const TLorentzVector& p4b1) : 
  Higgs(p4Higgs),
  Dau0(p4b0),
  Dau1(p4b1) 
{
}
#endif 
