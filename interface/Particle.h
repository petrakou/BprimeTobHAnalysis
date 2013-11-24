#ifndef BPRIMETOBHANALYSIS_INTERFACE_H
#define BPRIMETOBHANALYSIS_INTERFACE_H

//#include "DataFormats/Candidate/interface/Candidate.h"
#include <TLorentzVector.h>

class Particle { 

  public : 
    Particle (const TLorentzVector& p4) ; 
    ~Particle () {} ;

  TLorentzVector p4() { return m_p4 ; } ; 

  private : 
    TLorentzVector m_p4 ; 

} ;

#endif 

#ifdef BPRIMETOBHANALYSIS_INTERFACE_H
Particle::Particle (const TLorentzVector& p4) :
  m_p4(p4)  
{
}
#endif 
