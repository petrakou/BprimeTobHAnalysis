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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "BprimebHAnalysis/BprimeTobH/interface/format.h"
#include "BprimebHAnalysis/BprimeTobH/interface/TriggerBooking.h"
#include "BprimebHAnalysis/BprimeTobH/interface/Njettiness.hh"
#include "BprimebHAnalysis/BprimeTobH/interface/Nsubjettiness.hh"

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

		//// Configurables 
		int              maxEvents_;
		std::string	 inFile_;
		//const std::vector<std::string>  inFile_;
		//bool               debug;

		TChain*            chain_;
		//std::string        outFile_;
		//TTree		 *newtree;
		
                double fatJetPtMin_ ; 
		double fatJetPtMax_ ; 
		double fatJetAbsEtaMax_ ;
		double fatJetMassMin_ ;
		double fatJetMassMax_ ; 
		double fatJetPrunedMassMin_ ;
		double fatJetPrunedMassMax_ ; 
		double fatJetTau2ByTau1Min_ ; 
		double subjet1CSVDiscMin_ ; 
		double subjet1CSVDiscMax_ ; 
		double subjet2CSVDiscMin_ ; 
		double subjet2CSVDiscMax_ ; 

		EvtInfoBranches    EvtInfo;
		VertexInfoBranches VtxInfo;
		GenInfoBranches    GenInfo;
		JetInfoBranches    GenJetInfo;
		JetInfoBranches    JetInfo;
		JetInfoBranches    FatJetInfo;
		JetInfoBranches    SubJetInfo;
		LepInfoBranches    LepInfo;

		edm::Service<TFileService> fs; 

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
		TH1D*              h_FatJets_tau2ByTau1; 
		TH1D*              h_FatJets_tau3ByTau2; 
		TH1D*              h_FatJets_tau3ByTau1; 


		TH1D*              h_SubJets_Pt;
		// Specific to Subjets:
		TH1I*              h_SubJets_Jet_FatJetIdx;

		TH1D* h_SubJet1_Pt ; 
		TH1D* h_SubJet1_Eta ;
		TH1D* h_SubJet1_Mass ; 
                TH1D* h_SubJet1_CombinedSVBJetTags ; 

		TH1D* h_SubJet2_Pt ; 
		TH1D* h_SubJet2_Eta ;
		TH1D* h_SubJet2_Mass ; 
                TH1D* h_SubJet2_CombinedSVBJetTags ; 

		TH1D* h_HiggsJet_Pt ; 
		TH1D* h_HiggsJet_Eta ;
		TH1D* h_HiggsJet_Mass ; 

};

//
// constructors and destructor
//
BprimeTobHAnalysis::BprimeTobHAnalysis(const edm::ParameterSet& iConfig) : 
	maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
	inFile_(iConfig.getParameter<std::string>("InputFiles")),
	//outFile_(iConfig.getUntrackedParameter<std::string>("OutputFile")), 
	fatJetPtMin_(iConfig.getParameter<double>("FatJetPtMin")),
	fatJetPtMax_(iConfig.getParameter<double>("FatJetPtMax")),
	fatJetAbsEtaMax_(iConfig.getParameter<double>("FatJetAbsEtaMax")),
	fatJetMassMin_(iConfig.getParameter<double>("FatJetMassMin")),
	fatJetMassMax_(iConfig.getParameter<double>("FatJetMassMax")), 
	fatJetPrunedMassMin_(iConfig.getParameter<double>("FatJetPrunedMassMin")),
	fatJetPrunedMassMax_(iConfig.getParameter<double>("FatJetPrunedMassMax")),
	fatJetTau2ByTau1Min_(iConfig.getParameter<double>("FatJetTau2ByTau1Min")),
	subjet1CSVDiscMin_(iConfig.getParameter<double>("Subjet1CSVDiscMin")),
	subjet1CSVDiscMax_(iConfig.getParameter<double>("Subjet1CSVDiscMax")),
	subjet2CSVDiscMin_(iConfig.getParameter<double>("Subjet2CSVDiscMin")),
	subjet2CSVDiscMax_(iConfig.getParameter<double>("Subjet2CSVDiscMax")) 
{ 
}


BprimeTobHAnalysis::~BprimeTobHAnalysis() { 
	delete chain_;
	//newfile_->Close();
	//delete newfile_;
}


// ------------ method called once each job just before starting event loop  ------------
void BprimeTobHAnalysis::beginJob() { 

	chain_ = new TChain("ntuple/tree");
	chain_->Add(inFile_.c_str());

	EvtInfo.Register(chain_);
	VtxInfo.Register(chain_);
	GenInfo.Register(chain_);
	GenJetInfo.Register(chain_,"GenJetInfo");
	JetInfo.Register(chain_,"JetInfo");
	FatJetInfo.Register(chain_,"FatJetInfo");
	SubJetInfo.Register(chain_,"SubJetInfo");
	JetInfo.Register(chain_,"JetInfo");
	LepInfo.Register(chain_);

	if(maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

	// example from ExcitedQuark 
	//eSelector = new eventSelector(SelectionParameters,EvtInfo,LepInfo,JetInfo,VtxInfo);
	//cutLevels = eSelector->getCutLevels();

	h_FatJets_Pt                 = fs->make<TH1D>("h_FatJets_Pt"                 ,"Fat jet p_{T} [GeV]"       ,100 ,0.  ,2000.);
	h_FatJets_Eta                = fs->make<TH1D>("h_FatJets_Eta"                ,"Fat jet #eta"              ,50  ,-4. ,4.   );
	h_FatJets_Mass               = fs->make<TH1D>("h_FatJets_Mass"               ,"Fat jet mass [GeV]"        ,100 ,0.  ,2000.); 
	h_FatJets_MassPruned         = fs->make<TH1D>("h_FatJets_MassPruned"         ,"Fat jet pruned mass [GeV]" ,100 ,0.  ,2000.); 
	h_FatJets_tau2ByTau1         = fs->make<TH1D>("h_FatJets_tau2ByTau1"         ,"Fat jet #tau2/#tau1"       ,20  ,0.  ,1.   ); 
	h_FatJets_tau3ByTau2         = fs->make<TH1D>("h_FatJets_tau3ByTau2"         ,"Fat jet #tau2/#tau1"       ,20  ,0.  ,1.   ); 
	h_FatJets_tau3ByTau1         = fs->make<TH1D>("h_FatJets_tau3ByTau1"         ,"Fat jet #tau2/#tau1"       ,20  ,0.  ,1.   ); 
	h_FatJets_CombinedSVBJetTags = fs->make<TH1D>("h_FatJets_CombinedSVBJetTags" ,"Fat jet CSV discriminator" ,20  ,0.  ,1.   ); 

	h_SubJet1_Pt                 = fs->make<TH1D>("h_SubJet1_Pt"                 ,"SubJet1 p_{T} [GeV]"       ,100 ,0.  ,2000.);
	h_SubJet1_Eta                = fs->make<TH1D>("h_SubJet1_Eta"                ,"SubJet1 #eta"              ,50  ,-4. ,4.   );
	h_SubJet1_Mass               = fs->make<TH1D>("h_SubJet1_Mass"               ,"SubJet1 mass [GeV]"        ,100 ,0.  ,2000.); 
	h_SubJet1_CombinedSVBJetTags = fs->make<TH1D>("h_SubJet1_CombinedSVBJetTags" ,"SubJet1 CSV discriminator" ,20  ,0.  ,1.   ); 

	h_SubJet2_Pt                 = fs->make<TH1D>("h_SubJet2_Pt"                 ,"SubJet2 p_{T} [GeV]"       ,100 ,0.  ,2000.);
	h_SubJet2_Eta                = fs->make<TH1D>("h_SubJet2_Eta"                ,"SubJet2 #eta"              ,50  ,-4. ,4.   );
	h_SubJet2_Mass               = fs->make<TH1D>("h_SubJet2_Mass"               ,"SubJet2 mass [GeV]"        ,100 ,0.  ,2000.); 
	h_SubJet2_CombinedSVBJetTags = fs->make<TH1D>("h_SubJet2_CombinedSVBJetTags" ,"SubJet2 CSV discriminator" ,20  ,0.  ,1.   ); 

	h_HiggsJet_Pt                = fs->make<TH1D>("h_HiggsJet_Pt"                ,"Higgs jet p_{T} [GeV]"     ,100 ,0.  ,2000.);
	h_HiggsJet_Eta               = fs->make<TH1D>("h_HiggsJet_Eta"               ,"Higgs jet #eta"            ,50  ,-4. ,4.   );
	h_HiggsJet_Mass              = fs->make<TH1D>("h_HiggsJet_Mass"              ,"Higgs jet mass [GeV]"      ,100 ,0.  ,200. ); 

	h_FatJets_Pt                 -> Sumw2() ; 
	h_FatJets_Eta                -> Sumw2() ; 
	h_FatJets_Mass               -> Sumw2() ; 
	h_FatJets_MassPruned         -> Sumw2() ; 
	h_FatJets_tau2ByTau1         -> Sumw2() ; 
	h_FatJets_tau3ByTau2         -> Sumw2() ; 
	h_FatJets_tau3ByTau1         -> Sumw2() ; 
	h_FatJets_CombinedSVBJetTags -> Sumw2() ; 

	h_SubJet1_Pt                 -> Sumw2() ; 
	h_SubJet1_Eta                -> Sumw2() ; 
	h_SubJet1_Mass               -> Sumw2() ; 
	h_SubJet1_CombinedSVBJetTags -> Sumw2() ; 

	h_SubJet2_Pt                 -> Sumw2() ; 
	h_SubJet2_Eta                -> Sumw2() ; 
	h_SubJet2_Mass               -> Sumw2() ; 
	h_SubJet2_CombinedSVBJetTags -> Sumw2() ; 

	h_HiggsJet_Pt                -> Sumw2() ;
	h_HiggsJet_Eta               -> Sumw2() ;
	h_HiggsJet_Mass              -> Sumw2() ;

	return ;  

}


// ------------ method called for each event  ------------
void BprimeTobHAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
	using namespace edm;
	using namespace std;

	if(chain_ == 0) return;

	edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

	for(int entry=0; entry<maxEvents_; entry++) {

		if((entry%100)==0) printf("\nEvent %i of %i.\n", entry, maxEvents_);

		chain_->GetEntry(entry);

		//DM edm::LogInfo("FatJetSize") << "Fat size " << FatJetInfo.Size; // Test 
		for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

			h_FatJets_Pt->Fill(FatJetInfo.Pt[ifatjet]);

			//// Fat jet selection
			if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) 
				continue; //// apply jet pT cut
			if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
				continue; //// apply jet eta cut
			if ( FatJetInfo.JetIDLOOSE[ifatjet]==0 ) 
				continue; //// apply loose jet ID
			if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
					|| FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) continue; //// apply pruned jet mass cut

				TLorentzVector fatjet_p4;
				fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
						FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

				h_FatJets_Pt->Fill(fatjet_p4.Pt()) ;  
				h_FatJets_Eta->Fill(fatjet_p4.Eta()) ; 
				h_FatJets_Mass->Fill(fatjet_p4.Mag()) ; 
				h_FatJets_MassPruned->Fill(FatJetInfo.MassPruned[ifatjet]) ; 
				h_FatJets_CombinedSVBJetTags->Fill(FatJetInfo.CombinedSVBJetTags[ifatjet]) ; 
				h_FatJets_tau2ByTau1->Fill(FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet]) ; 
				h_FatJets_tau3ByTau2->Fill(FatJetInfo.tau3[ifatjet]/FatJetInfo.tau2[ifatjet]) ; 
				h_FatJets_tau3ByTau1->Fill(FatJetInfo.tau3[ifatjet]/FatJetInfo.tau1[ifatjet]) ; 

				//// Get subjets of fat jets 

				int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
				int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

				edm::LogInfo("SubJetIndices") << " Subjet1 : " << iSubJet1 << " Subjet2 : " << iSubJet2 ; 

				if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) 
					continue; //// skip fat jets for which one of the subjets has pT=0

				TLorentzVector subjet1_p4, subjet2_p4;
				subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], 
						SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
				subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], 
						SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

				double subjet_dR = subjet1_p4.DeltaR(subjet2_p4);
				double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
				double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
				double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

				if( subjet_dyphi < (FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]) ) 
					continue; //// skip infrared unsafe configurations

				h_SubJet1_Pt->Fill(subjet1_p4.Pt()) ;  
				h_SubJet1_Eta->Fill(subjet1_p4.Eta()) ; 
				h_SubJet1_Mass->Fill(subjet1_p4.Mag()) ; 
				h_SubJet1_CombinedSVBJetTags->Fill(SubJetInfo.CombinedSVBJetTags[iSubJet1]) ; 

				h_SubJet2_Pt->Fill(subjet2_p4.Pt()) ;  
				h_SubJet2_Eta->Fill(subjet2_p4.Eta()) ; 
				h_SubJet2_Mass->Fill(subjet2_p4.Mag()) ; 
				h_SubJet2_CombinedSVBJetTags->Fill(SubJetInfo.CombinedSVBJetTags[iSubJet2]) ; 

				//// Higgs tagging
				if (fatjet_p4.Mag() > fatJetMassMin_ 
						&& fatjet_p4.Mag() < fatJetMassMax_ 
						&& FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Min_ 
						&& SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
						&& SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
						&& SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
						&& SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 
					h_HiggsJet_Pt   -> Fill(fatjet_p4.Pt()) ; 
					h_HiggsJet_Eta  -> Fill(fatjet_p4.Eta()) ;
					h_HiggsJet_Mass -> Fill(fatjet_p4.Mag()) ;
				}

		} //// Loop over fat jets 

		//DM edm::LogInfo("SubJetSize") << "Sub size " << SubJetInfo.Size; // Test
		//DM for (int isubjet=0; isubjet < SubJetInfo.Size; ++isubjet) {
		//DM 	h_SubJets_Pt->Fill(SubJetInfo.Pt[isubjet]);
		//DM } //// Loop over subjets  

		//DM for (int ii=0; ii < GenInfo.Size; ++ii) {
		//DM 	if (GenInfo.Status[ii] == 3 && (TMath::Abs(GenInfo.PdgID[ii])<=6 || GenInfo.PdgID[ii]==22) ) {
		//DM 		//h_ptJets->Fill(GenInfo.Pt[ii]);
		//DM 	}
		//DM }

	} //// entry loop 

}


// ------------ method called once each job just after ending the event loop  ------------
void BprimeTobHAnalysis::endJob() { 
	//newfile_->mkdir("BprimeAnalysis");
	//newfile_->cd("BprimeAnalysis");

	//h_FatJets_Pt->Write();
	//h_SubJets_Pt->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BprimeTobHAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BprimeTobHAnalysis);
