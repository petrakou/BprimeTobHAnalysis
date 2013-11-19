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
// Second Author:    Devdatta Majumder 
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
#include <vector>
#include <map>

// Root headers 
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TEfficiency.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "../../BprimeTobH/interface/format.h"
#include "../../BprimeTobH/interface/TriggerBooking.h"
#include "../../BprimeTobH/interface/Njettiness.hh"
#include "../../BprimeTobH/interface/Nsubjettiness.hh"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "../../BprimeTobHAnalysis/interface/JetID.h"

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

    void CreateHistos(const TString&) ; 
    void AddHisto(const TString&, const TString&, const TString&, const int&, const double&, const double&) ; 
    template <class Type>
      void FillHisto(const TString& name, const Type value, const double weight);

    edm::LumiReWeighting LumiWeights_; 

    // ----------member data ---------------------------

    //// Configurables 

    int                             maxEvents_; 
    const int                       reportEvery_; 
    const std::string               inputTTree_;
    const std::vector<std::string>  inputFiles_;
    const std::vector<int>          hltPaths_; 
    const int                       doPUReweighting_ ;
    const std::string               file_PUDistMC_ ;
    const std::string               file_PUDistData_ ;
    const std::string               hist_PUDistMC_ ;
    const std::string               hist_PUDistData_ ;

    const double jetPtMin_ ; 
    const double jetPtMax_ ; 
    const double jetAbsEtaMax_ ;
    const double bjetPtMin_ ; 
    const double fatJetPtMin_ ; 
    const double fatJetPtMax_ ; 
    const double fatJetAbsEtaMax_ ;
    const double fatJetMassMin_ ;
    const double fatJetMassMax_ ; 
    const double fatJetPrunedMassMin_ ;
    const double fatJetPrunedMassMax_ ; 
    const double fatJetTau2ByTau1Max_ ; 
    const double subjet1CSVDiscMin_ ; 
    const double subjet1CSVDiscMax_ ; 
    const double subjet2CSVDiscMin_ ; 
    const double subjet2CSVDiscMax_ ; 
    const double HTMin_ ; 
    const double HTMax_ ; 

    TChain*            chain_;

    EvtInfoBranches    EvtInfo;
    VertexInfoBranches VtxInfo;
    GenInfoBranches    GenInfo;
    JetInfoBranches    GenJetInfo;
    JetInfoBranches    JetInfo;
    JetInfoBranches    FatJetInfo;
    JetInfoBranches    SubJetInfo;
    LepInfoBranches    LepInfo;

    edm::Service<TFileService> fs; 

    bool isData_ ; 
    double evtwt_ ; 
    double puweight_ ; 

    TH1D* h_cutflow ; 

    std::map<TString, TH1D*> hmap_1d ;  

};

//
// constructors and destructor
//
BprimeTobHAnalysis::BprimeTobHAnalysis(const edm::ParameterSet& iConfig) : 
  maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
  reportEvery_(iConfig.getParameter<int>("ReportEvery")),
  inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
  inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")),
  hltPaths_(iConfig.getParameter<std::vector<int> >("HLTPaths")),
  doPUReweighting_(iConfig.getParameter<bool>("DoPUReweighting")), 
  file_PUDistMC_(iConfig.getParameter<std::string>("File_PUDistMC")),
  file_PUDistData_(iConfig.getParameter<std::string>("File_PUDistData")),
  hist_PUDistMC_(iConfig.getParameter<std::string>("Hist_PUDistMC")),
  hist_PUDistData_(iConfig.getParameter<std::string>("Hist_PUDistData")),
  jetPtMin_(iConfig.getParameter<double>("JetPtMin")),
  jetPtMax_(iConfig.getParameter<double>("JetPtMax")),
  jetAbsEtaMax_(iConfig.getParameter<double>("JetAbsEtaMax")),
  bjetPtMin_(iConfig.getParameter<double>("BJetPtMin")),
  fatJetPtMin_(iConfig.getParameter<double>("FatJetPtMin")),
  fatJetPtMax_(iConfig.getParameter<double>("FatJetPtMax")),
  fatJetAbsEtaMax_(iConfig.getParameter<double>("FatJetAbsEtaMax")),
  fatJetMassMin_(iConfig.getParameter<double>("FatJetMassMin")),
  fatJetMassMax_(iConfig.getParameter<double>("FatJetMassMax")), 
  fatJetPrunedMassMin_(iConfig.getParameter<double>("FatJetPrunedMassMin")),
  fatJetPrunedMassMax_(iConfig.getParameter<double>("FatJetPrunedMassMax")),
  fatJetTau2ByTau1Max_(iConfig.getParameter<double>("FatJetTau2ByTau1Max")),
  subjet1CSVDiscMin_(iConfig.getParameter<double>("Subjet1CSVDiscMin")),
  subjet1CSVDiscMax_(iConfig.getParameter<double>("Subjet1CSVDiscMax")),
  subjet2CSVDiscMin_(iConfig.getParameter<double>("Subjet2CSVDiscMin")),
  subjet2CSVDiscMax_(iConfig.getParameter<double>("Subjet2CSVDiscMax")),
  HTMin_(iConfig.getParameter<double>("HTMin")), 
  HTMax_(iConfig.getParameter<double>("HTMax")),
  isData_(0),
  evtwt_(1), 
  puweight_(1)  
{ 

  if (doPUReweighting_) LumiWeights_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_) ;

}


BprimeTobHAnalysis::~BprimeTobHAnalysis() { 
  delete chain_;
}

// ------------ method called once each job just before starting event loop  ------------
void BprimeTobHAnalysis::beginJob() { 

  chain_ = new TChain(inputTTree_.c_str());

  for(unsigned i=0; i<inputFiles_.size(); ++i) {
    chain_->Add(inputFiles_.at(i).c_str());

    TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
    f->Close();
  }

  EvtInfo.Register(chain_);
  VtxInfo.Register(chain_);
  GenInfo.Register(chain_);
  GenJetInfo.Register(chain_,"GenJetInfo");
  JetInfo.Register(chain_,"JetInfo");
  FatJetInfo.Register(chain_,"FatJetInfo");
  SubJetInfo.Register(chain_,"SubJetInfo");
  LepInfo.Register(chain_);

  if(maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

  // example from ExcitedQuark 
  //eSelector = new eventSelector(SelectionParameters,EvtInfo,LepInfo,JetInfo,VtxInfo);
  //cutLevels = eSelector->getCutLevels();

  h_cutflow                    = fs->make<TH1D>("h_cutflow"                   ,"Cut flow"                   ,20  ,0.  ,20.  ); 
  h_cutflow -> Sumw2() ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(1,"AllEvents") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(2,"TriggerSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(3,"FatJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(4,"HiggsJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(5,"JetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(6,"BJetsSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(7,"HTSel") ; 

  for (int ii = 1; ii <= 7; ++ii) 
    CreateHistos(h_cutflow->GetXaxis()->GetBinLabel(ii)) ; 

  return ;  

}

void BprimeTobHAnalysis::CreateHistos(const TString& cutname) {

  AddHisto(cutname ,"_nPVtx_NoPUWt"               ,"N(PV), No PU weight"       ,50     ,-0.5     ,49.5    ) ; 
  AddHisto(cutname ,"_nPVtx_PUWt"                 ,"N(PV)"                     ,50     ,-0.5     ,49.5    ) ; 
  AddHisto(cutname ,"_nJets"                      ,"N(AK5 jets)"               ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_nBJets"                     ,"N(b jets)"                 ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_nFatJets"                   ,"N(fat jets)"               ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_nHJets"                     ,"N(Higgs jets)"             ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_FatJets_Pt"                 ,"p_{T}(fat jets)"           ,1000   ,0.       ,1000.   ) ; 
  AddHisto(cutname ,"_FatJets_Eta"                ,"#eta(fat jets)"            ,50     ,-4.      ,4.      ) ; 
  AddHisto(cutname ,"_FatJets_Mass"               ,"Fat jet mass [GeV]"        ,100    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_FatJets_MassPruned"         ,"Fat jet pruned mass [GeV]" ,100    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_FatJets_tau2ByTau1"         ,"Fat jet #tau2/#tau1"       ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_tau3ByTau2"         ,"Fat jet #tau2/#tau1"       ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_tau3ByTau1"         ,"Fat jet #tau2/#tau1"       ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_CombinedSVBJetTags" ,"Fat jet CSV discriminator" ,20     ,0.       ,1.      ) ; 

  AddHisto(cutname ,"_SubJet1_Pt"                 ,"SubJet1 p_{T} [GeV]"       ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_SubJet1_Eta"                ,"SubJet1 #eta"              ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_SubJet1_Mass"               ,"SubJet1 mass [GeV]"        ,100    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_SubJet1_CombinedSVBJetTags" ,"SubJet1 CSV discriminator" ,20     ,0.       ,1.      ) ; 

  AddHisto(cutname ,"_SubJet2_Pt"                 ,"SubJet2 p_{T} [GeV]"       ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_SubJet2_Eta"                ,"SubJet2 #eta"              ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_SubJet2_Mass"               ,"SubJet2 mass [GeV]"        ,100    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_SubJet2_CombinedSVBJetTags" ,"SubJet2 CSV discriminator" ,20     ,0.       ,1.      ) ; 

  AddHisto(cutname ,"_HiggsJet_Pt"                ,"p_{T} (Higgs jet)[GeV]"    ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_HiggsJet_Eta"               ,"#eta (Higgs jet)"          ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_HiggsJet_Mass"              ,"Mass (Higgs jet) [GeV]"    ,100    ,0.       ,200.    ) ; 

  AddHisto(cutname ,"_BJet_Pt"                    ,"p_{T} (b jet)[GeV]"        ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_BJet_Eta"                   ,"#eta (b jet)"              ,50     ,-4.      ,4.      ) ;

  AddHisto(cutname ,"_HT"                         ,"H_{T}[GeV]"                 ,200   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_bprimePt"                   ,"b' p_{T} [GeV]"             ,100   ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_bprimeMass"                 ,"b' mass [GeV]"              ,40    ,0.       ,2000.   ) ;

  return ; 

}

void BprimeTobHAnalysis::AddHisto(const TString& cutname, const TString& histname, const TString& histtitle, const int& nbins, const double& min, const double& max) { 

  TH1D* h1d ; 
  h1d = fs->make<TH1D>(cutname+histname, cutname+histtitle, nbins, min, max);  
  h1d -> Sumw2() ; 
  hmap_1d[cutname+histname] = h1d ; 

  //TH1D* h1d_mc ; 
  //h1d_mc = fs->make<TH1D>(cutname+histname+"_mc", cutname+histtitle, nbins, min, max);  
  //h1d_mc -> Sumw2() ; 
  //hmap_1d[cutname+histname+"_mc"] = h1d_mc ; 

  //TH1D* h1d_data ; 
  //h1d_data = fs->make<TH1D>(cutname+histname+"_data", cutname+histtitle, nbins, min, max);  
  //h1d_data -> Sumw2() ; 
  //hmap_1d[cutname+histname+"_data"] = h1d_data ; 

  return ; 

}

template <class Type>
void BprimeTobHAnalysis::FillHisto(const TString& name, const Type value, const double weight){
  hmap_1d[name]->Fill(double(value),weight);
  //if (!isData_) hmap_1d[name+"_mc"]->Fill(double(value),weight);
  //else hmap_1d[name+"_data"]->Fill(double(value),weight); 

  return ; 

}

// ------------ method called for each event  ------------
void BprimeTobHAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
  using namespace edm;
  using namespace std;

  if(chain_ == 0) return;

  JetID jetIDTight(JetID::FIRSTDATA,JetID::TIGHT, JetInfo) ; 
  JetID fatjetIDLoose(JetID::FIRSTDATA,JetID::LOOSE, FatJetInfo) ; 
  pat::strbitset retak5 = jetIDTight.getBitTemplate() ;
  pat::strbitset retca8 = fatjetIDLoose.getBitTemplate() ;

  ofstream fout("Evt_NoJets.txt") ; 
  if ( isData_ ) {
    fout << "EvtInfo.RunNo " << " EvtInfo.LumiNo " << " EvtInfo.EvtNo " << std::endl ;
  }

  edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

  for(int entry=0; entry<maxEvents_; entry++) {

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_ ; 

    //// Event variables 
    std::vector<TLorentzVector>fatJets ; 
    std::vector<TLorentzVector>higgsJets ; 
    std::vector<TLorentzVector>jets ; 
    std::vector<TLorentzVector>bJets ; 
    //std::vector<std::pair<int,TLorentzVector> > higgsJets ; 
    //std::vector<std::pair<int,TLorentzVector> > jets ; 
    //std::vector<std::pair<int,TLorentzVector> > bJets ; 
    std::vector<TLorentzVector>bprimes ; 

    bool passHLT(false) ; 

    int nGoodVtxs(0) ;
    int njets(0) ; 
    double HT(0) ; 

    chain_->GetEntry(entry);

    isData_   = EvtInfo.McFlag ? 0 : 1; 
    if ( !isData_ ) evtwt_    = EvtInfo.Weight ; 
    if ( doPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]) ; 

    nGoodVtxs = 0 ;
    /**\ Select good vertices */
    for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) {
      if (   VtxInfo.Type[iVtx]==1
          && VtxInfo.isFake[iVtx]==false
          && VtxInfo.Ndof[iVtx]>4
          && VtxInfo.Rho[iVtx]<2.
          && VtxInfo.z[iVtx]<24.) { ++nGoodVtxs ; }
    }
    if (nGoodVtxs < 1)  { edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex " ; continue ; }

    FillHisto(TString("AllEvents")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt_) ; 
    FillHisto(TString("AllEvents")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt_*puweight_) ; 
    FillHisto(TString("AllEvents")+TString("_nJets"), JetInfo.Size, evtwt_*puweight_) ; 
    FillHisto(TString("AllEvents")+TString("_nFatJets"), FatJetInfo.Size, evtwt_*puweight_) ; 
    h_cutflow -> Fill("AllEvents", 1) ; 

    for ( std::vector<int>::const_iterator ihlt = hltPaths_.begin();
        ihlt != hltPaths_.end(); ++ihlt ) { 
      if (EvtInfo.TrgBook[*ihlt] == 1) { 
        passHLT = true ; 
        break ; 
      }
      else passHLT = false ; 
    }

    if ( !passHLT ) continue ; 

    if ( isData_ ) {
      if ( JetInfo.Size == 0 ) fout << EvtInfo.RunNo << " " << EvtInfo.LumiNo << " " << EvtInfo.EvtNo << std::endl ; 
    }

    FillHisto(TString("TriggerSel")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt_) ; 
    FillHisto(TString("TriggerSel")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt_*puweight_) ; 
    FillHisto(TString("TriggerSel")+TString("_nJets"), JetInfo.Size, evtwt_*puweight_) ; 
    FillHisto(TString("TriggerSel")+TString("_nFatJets"), FatJetInfo.Size, evtwt_*puweight_) ; 
    h_cutflow -> Fill("TriggerSel", 1) ; 

    evtwt_ *= puweight_ ; 

    TLorentzVector higgs_p4 ; 
    for (int igen=0; igen < GenInfo.Size; ++igen) {

      if ( GenInfo.Status[igen] == 3 
          && TMath::Abs(GenInfo.PdgID[igen])==25 
          //&& GenInfo.nDa[igen] == 2 
          //&& TMath::Abs(GenInfo.PdgID[GenInfo.Da1[igen]])==5 
          //&& TMath::Abs(GenInfo.PdgID[GenInfo.Da2[igen]])==5 
         ) { 

        higgs_p4.SetPtEtaPhiM(GenInfo.Pt[igen], GenInfo.Eta[igen], GenInfo.Phi[igen], GenInfo.Mass[igen]) ; 
        TLorentzVector fatjet_p4;
        bool matched ; 

        //h_HiggsPt ->Fill(higgs_p4.Pt()) ; 

        for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) { 

          if ( fabs(FatJetInfo.Eta[ifatjet]) > 1.5 ) continue ; 

          fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
              FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

          if (higgs_p4.DeltaR(fatjet_p4) < 0.5) {
            matched = true ; 
            //h_HiggsPtMatchedJet ->Fill(higgs_p4.Pt()) ; 
            break ; 
          }
          else 
            matched = false ; 

        } //// Loop over fat jets 

        //teff_HiggsJetMatch -> Fill(matched, GenInfo.Pt[igen]) ; 

      } //// Get Higgs boson 

    } //// Loop over all gen particles 

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //Fix h_FatJets_Pt->Fill(FatJetInfo.Pt[ifatjet]);

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ 
          || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) continue; //// apply jet pT cut
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) continue; //// apply jet eta cut
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) continue; //// apply pruned jet mass cut 
      retca8.set(false);
      if ( fatjetIDLoose(ifatjet,retca8) == 0 ) continue; //// apply loose jet ID

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

      FillHisto(TString("TriggerSel")+TString("_FatJets_Pt")                 ,fatjet_p4.Pt() ,evtwt_)  ;  
      FillHisto(TString("TriggerSel")+TString("_FatJets_Eta")                ,fatjet_p4.Eta() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_Mass")               ,fatjet_p4.Mag() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_MassPruned")         ,FatJetInfo.MassPruned[ifatjet] ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_CombinedSVBJetTags") ,FatJetInfo.CombinedSVBJetTags[ifatjet] ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_tau2ByTau1")         ,FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_tau3ByTau2")         ,FatJetInfo.tau3[ifatjet]/FatJetInfo.tau2[ifatjet] ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_tau3ByTau1")         ,FatJetInfo.tau3[ifatjet]/FatJetInfo.tau1[ifatjet] ,evtwt_)  ; 

      //// Get subjets of fat jets 

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

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

      FillHisto(TString("TriggerSel")+TString("_SubJet1_Pt") ,subjet1_p4.Pt() ,evtwt_)  ;  
      FillHisto(TString("TriggerSel")+TString("_SubJet1_Eta") ,subjet1_p4.Eta() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet1_Mass") ,subjet1_p4.Mag() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet1_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet1] ,evtwt_)  ; 

      FillHisto(TString("TriggerSel")+TString("_SubJet2_Pt") ,subjet2_p4.Pt() ,evtwt_)  ;  
      FillHisto(TString("TriggerSel")+TString("_SubJet2_Eta") ,subjet2_p4.Eta() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet2_Mass") ,subjet2_p4.Mag() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet2_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet2] ,evtwt_)  ; 

      //// Selecting fat jets 
      if (fatjet_p4.Mag() > fatJetMassMin_ 
          && fatjet_p4.Mag() < fatJetMassMax_ 
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_) {

        fatJets.push_back(fatjet_p4) ; 

        //// Higgs tagging
        if ( SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
            && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
            && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
            && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 
          FillHisto(TString("TriggerSel")+TString("_HiggsJet_Pt")   ,fatjet_p4.Pt() ,evtwt_)  ; 
          FillHisto(TString("TriggerSel")+TString("_HiggsJet_Eta")  ,fatjet_p4.Eta() ,evtwt_)  ;
          FillHisto(TString("TriggerSel")+TString("_HiggsJet_Mass") ,fatjet_p4.Mag() ,evtwt_)  ;

          higgsJets.push_back(fatjet_p4) ; 

        } //// Higgs tagging 

      } //// Selecting fat jets 

    } //// Loop over fat jets 

    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      retak5.set(false);
      if ( jetIDTight(ijet,retak5) == 0 ) continue; 

      ++njets ; 

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

      bool isJetNotHiggs(false) ; 
      for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) {
          isJetNotHiggs = false ; 
          break ; 
        }
        else {
          isJetNotHiggs = true ; 
        } 
      } 
      if (!isJetNotHiggs) continue ; //// Higgs-b jet disambiguation  

      jets.push_back(jet_p4) ; 

      if (JetInfo.Pt[ijet] > bjetPtMin_ && JetInfo.CombinedSVBJetTags[ijet] > 0.679) {

        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
            JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

        bJets.push_back(bjet_p4) ; 

      } //// Select b-tagged AK5 jets 

    } //// Loop over AK5 jets 

    if (fatJets.size() >= 1) {
      FillHisto(TString("FatJetSel")+TString("_nJets"), njets, evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_nBJets"), bJets.size(), evtwt_) ; 
      for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) { 
        FillHisto(TString("FatJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
        FillHisto(TString("FatJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
      }
      h_cutflow -> Fill("FatJetSel", 1) ; 
    }

    if (higgsJets.size() >= 1) {

      FillHisto(TString("HiggsJetSel")+TString("_nJets"), njets, evtwt_) ; 
      FillHisto(TString("HiggsJetSel")+TString("_nBJets"), bJets.size(), evtwt_) ; 
      for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) { 
        FillHisto(TString("HiggsJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
      }
      h_cutflow -> Fill("HiggsJetSel", 1) ; 

      if (bJets.size() >= 2 ) { 

        FillHisto(TString("BJetsSel")+TString("_nJets"), JetInfo.Size, evtwt_) ; 
        h_cutflow -> Fill("BJetsSel", 1) ;

        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) { 
          HT += ihig->Pt() ; 
        }
        for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) { 
          HT += ib->Pt() ; 
        }

        if (HT < HTMin_ || HT > HTMax_) continue ; 

        FillHisto(TString("HTSel")+TString("_nJets"), njets, evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_nBJets"), bJets.size(), evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_nHJets"), higgsJets.size(), evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_HT") ,HT ,evtwt_)  ; 
        if ( !isData_ ) h_cutflow -> Fill("HTSel", 1) ; 

        //// Reconstruct b' candidates
        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) { 
          const TLorentzVector* closestB_p4 ;
          double deltaR(TMath::Pi()) ; 
          for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) { 
            if ( ihig->DeltaR(*ib) < deltaR) {
              deltaR = ihig->DeltaR(*ib) ; 
              closestB_p4 = &(*ib) ; 
            }
          }
          if (deltaR < TMath::Pi()) {
            bprimes.push_back(*ihig + *closestB_p4) ; 
            FillHisto(TString("HTSel")+TString("_bprimePt") ,(*ihig + *closestB_p4).Pt() ,evtwt_)  ; 
            FillHisto(TString("HTSel")+TString("_bprimeMass") ,(*ihig + *closestB_p4).Mag() ,evtwt_)  ;
          }
        } //// Reconstruct b' candidates

      } //// If at least two b-jets 
    } //// If at least one Higgs jet 

  } //// entry loop 

  fout.close() ; 

}

// ------------ method called once each job just after ending the event loop  ------------
void BprimeTobHAnalysis::endJob() { 
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
