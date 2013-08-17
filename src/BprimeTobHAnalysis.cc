// -*- C++ -*-
//
// Package:    BprimeTobHAnalysis
// Class:      BprimeTobHAnalysis
// 
/**\class BprimeTobHAnalysis BprimeTobHAnalysis.cc Bprime_kit/BprimeTobHAnalysis/src/BprimeTobHAnalysis.cc

 Description: 
 Analyzer class for Bprime -> b Higgs studies 
 - National Taiwan University - 

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eleni Petrakou,27 2-020,+41227674870,
//         Created:  Tue Jul 16 19:48:47 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TH1.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Bprime_kit/BprimeTobH/interface/format.h"
#include "Bprime_kit/BprimeTobH/interface/TriggerBooking.h"
#include "Bprime_kit/BprimeTobH/interface/Njettiness.hh"
#include "Bprime_kit/BprimeTobH/interface/Nsubjettiness.hh"

//
// class declaration
//

class BprimeTobHAnalysis : public edm::EDAnalyzer {
   public:
      explicit BprimeTobHAnalysis(const edm::ParameterSet&);
      ~BprimeTobHAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      int                maxEvents;
      std::string	 inFile;
      //const std::vector<std::string>  inFile;
      //bool               debug;

      TChain*            chain;
      std::string        outFile;
      TFile              *newfile;
      //TTree		 *newtree;

      EvtInfoBranches    EvtInfo;
      VertexInfoBranches VtxInfo;
      LepInfoBranches    LepInfo;
      JetInfoBranches    FatJetInfo;
      JetInfoBranches    SubJetInfo;
      JetInfoBranches    JetInfo;
      JetInfoBranches    GenJetInfo;
      GenInfoBranches    GenInfo;

      TH1I*              h_FatJets_Index;
      TH1I*              h_FatJets_NTracks;
      TH1D*              h_FatJets_Et;
      TH1D*              h_FatJets_Pt;
      TH1D*              h_FatJets_Eta;
      TH1D*              h_FatJets_Phi;
      TH1D*              h_FatJets_Energy;
      TH1D*              h_FatJets_Px;
      TH1D*              h_FatJets_Py;
      TH1D*              h_FatJets_Pz;
      TH1D*              h_FatJets_Mass;
      TH1D*              h_FatJets_Area;
      TH1I*              h_FatJets_JetCharge;
      TH1I*              h_FatJets_NConstituents;
      TH1D*              h_FatJets_NCH;
      TH1D*              h_FatJets_CEF;
      TH1D*              h_FatJets_NHF;
      TH1D*              h_FatJets_NEF;
      TH1D*              h_FatJets_CHF;
      TH1D*              h_FatJets_QGTagsMLP;
      TH1D*              h_FatJets_QGTagsLikelihood;
      TH1I*              h_FatJets_JetIDLOOSE;
      TH1I*              h_FatJets_JetIDTIGHT;
      TH1D*              h_FatJets_PtCorrRaw;
      TH1D*              h_FatJets_PtCorrL2;
      TH1D*              h_FatJets_PtCorrL3;
      TH1D*              h_FatJets_PtCorrL7g;
      TH1D*              h_FatJets_PtCorrL7uds;
      TH1D*              h_FatJets_PtCorrL7c;
      TH1D*              h_FatJets_PtCorrL7b;
      TH1D*              h_FatJets_JetBProbBJetTags;
      TH1D*              h_FatJets_JetProbBJetTags; 
      TH1D*              h_FatJets_TrackCountHiPurBJetTags; 
      TH1D*              h_FatJets_CombinedSVBJetTags; 
      TH1D*              h_FatJets_CombinedSVMVABJetTags; 
      TH1D*              h_FatJets_SoftElecByIP3dBJetTags; 
      TH1D*              h_FatJets_SoftElecByPtBJetTags;
      TH1D*              h_FatJets_SoftMuonBJetTags;
      TH1D*              h_FatJets_SoftMuonByIP3dBJetTags; 
      TH1D*              h_FatJets_SoftMuonByPtBJetTags;
      TH1D*              h_FatJets_DoubleSVHighEffBJetTags; 
      TH1D*              h_FatJets_GenJetPt; 
      TH1D*              h_FatJets_GenJetEta;
      TH1D*              h_FatJets_GenJetPhi;
      TH1I*              h_FatJets_GenMCTag;
      TH1D*              h_FatJets_GenPt; // for partons 
      TH1D*              h_FatJets_GenEta;
      TH1D*              h_FatJets_GenPhi;
      TH1D*              h_FatJets_GenPdgID;
      TH1D*              h_FatJets_GenFlavor;
      // Specific to Fat Jets: 
      TH1D*              h_FatJets_EtPruned;
      TH1D*              h_FatJets_PtPruned;
      TH1D*              h_FatJets_EtaPruned;
      TH1D*              h_FatJets_PhiPruned;
      TH1D*              h_FatJets_EnergyPruned;
      TH1D*              h_FatJets_PxPruned; 
      TH1D*              h_FatJets_PyPruned; 
      TH1D*              h_FatJets_PzPruned; 
      TH1D*              h_FatJets_MassPruned; 
      TH1D*              h_FatJets_AreaPruned; 
      TH1I*              h_FatJets_Jet_SubJet1Idx;
      TH1I*              h_FatJets_Jet_SubJet2Idx;
      TH1D*              h_FatJets_tau1;
      TH1D*              h_FatJets_tau2;
      TH1D*              h_FatJets_tau3;


      TH1D*              h_SubJets_Pt;
// -- to be filled......... 
      // Specific to Subjets:
      TH1I*              h_SubJets_Jet_FatJetIdx;


};

//
// constructors and destructor
//
BprimeTobHAnalysis::BprimeTobHAnalysis(const edm::ParameterSet& iConfig)
{
  maxEvents = iConfig.getParameter<int>("MaxEvents");
  inFile = iConfig.getParameter<std::string>("InputFiles");
  outFile = iConfig.getUntrackedParameter<std::string>("OutputFile");
}


BprimeTobHAnalysis::~BprimeTobHAnalysis()
{
  delete chain;
  newfile->Close();
  delete newfile;
}


// ------------ method called once each job just before starting event loop  ------------
void BprimeTobHAnalysis::beginJob()
{

  chain = new TChain("ntuple/tree");
  chain->Add(inFile.c_str());

  EvtInfo.Register(chain);
  VtxInfo.Register(chain);
  LepInfo.Register(chain);
  GenInfo.Register(chain);
  FatJetInfo.Register(chain,"FatJetInfo");
  SubJetInfo.Register(chain,"SubJetInfo");
  JetInfo.Register(chain,"JetInfo");
  GenJetInfo.Register(chain,"GenJetInfo");

  if(maxEvents<0 || maxEvents>chain->GetEntries())
     maxEvents = chain->GetEntries();

   // example from ExcitedQuark 
   //eSelector = new eventSelector(SelectionParameters,EvtInfo,LepInfo,JetInfo,VtxInfo);
   //cutLevels = eSelector->getCutLevels();

   newfile = new TFile(outFile.c_str(),"recreate");
   newfile->cd();

   h_FatJets_Pt = new TH1D("h_FatJets_Pt","FatJets pT",50,0.,2000.);
// -- to be filled.........
   h_SubJets_Pt = new TH1D("h_SubJets_Pt","SubJets pT",50,0.,2000.);

   h_FatJets_Pt->Sumw2();
   h_SubJets_Pt->Sumw2();

}


// ------------ method called for each event  ------------
void
BprimeTobHAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   if(chain == 0) return;

   cout << endl << "Starting analysis loop\n";
   for(int entry=0; entry<maxEvents; entry++) {
	if((entry%100)==0)
		printf("\nEvent %i of %i.\n", entry, maxEvents);
     	chain->GetEntry(entry);

cout << endl << "Fat size " << FatJetInfo.Size; // Test 
	for (int ii=0; ii < FatJetInfo.Size; ++ii) {
		h_FatJets_Pt->Fill(FatJetInfo.Pt[ii]);
// -- and so on......... 
	}

cout << endl << "Sub size " << SubJetInfo.Size; // Test
        for (int ii=0; ii < SubJetInfo.Size; ++ii) {
                h_SubJets_Pt->Fill(SubJetInfo.Pt[ii]);
        }

        for (int ii=0; ii < GenInfo.Size; ++ii) {
                if (GenInfo.Status[ii] == 3 && (TMath::Abs(GenInfo.PdgID[ii])<=6 || GenInfo.PdgID[ii]==22) ) {
                      //h_ptJets->Fill(GenInfo.Pt[ii]);
                }
        }

   } // entry loop 

}


// ------------ method called once each job just after ending the event loop  ------------
void 
BprimeTobHAnalysis::endJob() 
{
    newfile->mkdir("BprimeAnalysis");
    newfile->cd("BprimeAnalysis");

    h_FatJets_Pt->Write();
    h_SubJets_Pt->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BprimeTobHAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BprimeTobHAnalysis);
