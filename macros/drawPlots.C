#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <THStack.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <TPaveText.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include "CMSstyle.C"
#include "help.C"

using namespace std;

//TString filename = "Final_histograms_BprimebH_PURewt_CorrJetId.root" ; 
TString filename = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_11_BpbH/src/BprimeTobHAnalysisv1/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_01_PURewt_CorrJetId/Final_histograms_BprimebH.root" ; 

Double_t Lint = 19700.0 ; 
TString title1 = "CMS Preliminary, 19.7/fb at #sqrt{s} = 8 TeV";
TString datacaption = "Data"; 

TString dir4plots ="BprimeTobH_01_PURewt_CorrJetId" ;

TString formata = ".pdf";
TString formatb = ".png";
TString formatc = ".C";

bool web = 0;
bool setSampleName = 1;

void DrawAll () ; 
void DrawStacked(TString name, TString histotitle, bool log, bool doData, bool fExtNorm=false, int nRebin=1, bool setXRange=false, double rangeXLow=0., double rangeXHigh=0.);

void drawPlots () {

  gROOT->SetBatch(kTRUE);
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  CMSstyle() ; 

  TString action = "mkdir -p " + dir4plots;
  system(action);

  DrawAll () ; 

  return ; 

}

void DrawAll () {

   DrawStacked("h_cutflow" ,"" ,1 ,1 ,0 ,1 ,1 ,0 ,7); 

   DrawStacked("TriggerSel_nPVtx_NoPUWt" ,"N(PV), No PU weight" ,0 ,1 ,0 ,1 ,1 ,0 ,50); 
   DrawStacked("TriggerSel_nPVtx_PUWt" ,"N(PV)" ,0 ,1 ,0 ,1 ,1 ,0 ,50); 
   DrawStacked("TriggerSel_nFatJets" ,"N(CA8 jets)" ,1 ,1 ,0 ,1 ,1 ,0 ,5); 
   DrawStacked("TriggerSel_nJets" ,"N(AK5 jets)" ,1 ,1 ,0 ,1 ,1 ,0 ,5); 
   DrawStacked("TriggerSel_FatJets_Pt" ,"p_{T}(CA8 jets) (GeV)" ,1 ,1 ,0 ,4 ,1 ,0 ,1000); 
   DrawStacked("TriggerSel_SubJet1_Pt" ,"p_{T}(leading subjet) (GeV)" ,1 ,1 ,0 ,1 ,1 ,0 ,1000); 
   DrawStacked("TriggerSel_SubJet2_Pt" ,"p_{T}(subleading subjet) (GeV)" ,1 ,1 ,0 ,1 ,1 ,0 ,1000); 
   DrawStacked("TriggerSel_FatJets_Mass" ,"M(CA8 jets) (GeV)" ,1 ,1 ,0 ,1 ,1 ,0 ,1000); 
   DrawStacked("TriggerSel_FatJets_MassPruned" ,"pruned mass(CA8 jets) (GeV)" ,1 ,1 ,0 ,1 ,1 ,0 ,1000); 
   DrawStacked("TriggerSel_FatJets_tau2ByTau1" ,"#tau_{2}/#tau_{1}(Pruned CA8 jets)" ,1 ,1 ,0 ,1 ,1 ,0 ,1); 
   DrawStacked("TriggerSel_FatJets_tau3ByTau1" ,"#tau_{3}/#tau_{1}(Pruned CA8 jets)" ,1 ,1 ,0 ,1 ,1 ,0 ,1); 
   DrawStacked("TriggerSel_FatJets_tau3ByTau2" ,"#tau_{3}/#tau_{2}(Pruned CA8 jets)" ,1 ,1 ,0 ,1 ,1 ,0 ,1); 
   DrawStacked("TriggerSel_SubJet1_CombinedSVBJetTags" ,"CSV discriminator(leading subjet)" ,1 ,1 ,0 ,1 ,1 ,0 ,1); 
   DrawStacked("TriggerSel_SubJet2_CombinedSVBJetTags" ,"CSV discriminator(subleading subjet)" ,1 ,1 ,0 ,1 ,1 ,0 ,1); 
   DrawStacked("TriggerSel_SubJet1_Mass" ,"M(leading subjet) (GeV)" ,1 ,1 ,0 ,1 ,1 ,0 ,600); 
   DrawStacked("TriggerSel_SubJet2_Mass" ,"M(subleading subjet) (GeV)" ,1 ,1 ,0 ,1 ,1 ,0 ,300); 

   DrawStacked("FatJetSel_nJets" ,"N(AK5 jets)" ,1 ,1 ,0 ,1 ,1 ,0 ,15); 
   DrawStacked("FatJetSel_nBJets" ,"N(b-tagged AK5 jets)" ,1 ,1 ,0 ,1 ,1 ,0 ,5); 
   DrawStacked("FatJetSel_BJet_Pt" ,"p_{T}(b-tagged AK5 jets) (GeV)" ,1 ,1 ,0 ,1 ,1 ,0 ,1000); 
   DrawStacked("FatJetSel_BJet_Eta" ,"#eta(b-tagged AK5 jets) (GeV)" ,1 ,1 ,0 ,1 ,1 , -3, 3); 

   DrawStacked("HiggsJetSel_nJets" ,"N(AK5 jets)" ,0 ,1 ,0 ,1 ,1 ,0 ,15); 
   DrawStacked("HiggsJetSel_nBJets" ,"N(b-tagged AK5 jets)" ,0 ,1 ,0 ,1 ,1 ,0 ,5); 
   DrawStacked("HiggsJetSel_BJet_Pt" ,"p_{T}(b-tagged AK5 jets) (GeV)" ,1 ,1 ,0 ,1 ,1 ,0 ,1000); 
   DrawStacked("HiggsJetSel_BJet_Eta" ,"#eta(b-tagged AK5 jets) (GeV)" ,1 ,1 ,0 ,1 ,1 , -3, 3); 

   DrawStacked("BJetsSel_nJets" ,"N(AK5 jets)" ,0 ,1 ,0 ,1 ,1 ,0 ,15); 

   DrawStacked("HTSel_nJets" ,"N(AK5 jets)" ,0 ,0 ,0 ,1 ,1 ,0 ,15); 
   DrawStacked("HTSel_nBJets" ,"N(b-tagged AK5 jets)" ,0 ,0 ,0 ,1 ,1 ,0 ,5); 
   DrawStacked("HTSel_nHJets" ,"N(Higgs-tagged CA8 jets)" ,0 ,0 ,0 ,1 ,1 ,0 ,5); 
   DrawStacked("HTSel_HT" ,"HT (GeV)" ,0 ,0 ,0 ,4 ,1 ,750 ,2050); 
   DrawStacked("HTSel_HT" ,"HT (GeV)" ,1 ,0 ,0 ,4 ,1 ,750 ,2050); 

   return ; 

}

void DrawStacked(TString name,
    TString histotitle,
    bool log,
    bool doData,
    bool fExtNorm,
    int nRebin,
    bool setXRange,
    double rangeXLow,
    double rangeXHigh) {

  TH1D* hist_ttjets ;
  TH1D* hist_qcd    ;
  TH1D* hist_sig0   ;
  TH1D* hist_sig1   ;
  TH1D* hist_sig2   ;
  TH1D* hist_data   ;

  TFile *myFile  = TFile::Open(filename,"READ") ;
  myFile->cd();

  hist_qcd              = (TH1D*)myFile->Get("QCD__"+name);
  hist_ttjets           = (TH1D*)myFile->Get("TTJets__"+name);
  hist_sig0             = (TH1D*)myFile->Get("BprimeBprimeToBHBHinc_M-500__"+name);
  hist_sig1             = (TH1D*)myFile->Get("BprimeBprimeToBHBHinc_M-800__"+name);
  hist_sig2             = (TH1D*)myFile->Get("BprimeBprimeToBHBHinc_M-1000__"+name);
  if (doData) hist_data = (TH1D*)myFile->Get("DATA__"+name);

  fix(hist_qcd   )           ; 
  fix(hist_ttjets)           ; 
  fix(hist_sig0  )           ; 
  fix(hist_sig1  )           ; 
  fix(hist_sig2  )           ; 
  if (doData) fix(hist_data) ; 

  if (nRebin > 1) {
    hist_ttjets -> Rebin(nRebin) ;
    hist_qcd    -> Rebin(nRebin) ;
    hist_sig0   -> Rebin(nRebin) ;
    hist_sig1   -> Rebin(nRebin) ;
    hist_sig2   -> Rebin(nRebin) ;
    if (doData) hist_data   -> Rebin(nRebin) ;
  }

  TH1D* hist_bkg = (TH1D*) hist_ttjets->Clone();
  hist_bkg->Add(hist_qcd) ; 

  if ( name.Contains("nPVtx") ) { 
    double scalef = hist_data->Integral()/hist_bkg->Integral() ; 
    hist_bkg->Scale(scalef) ; 
  }

  beautify(hist_qcd   ,42 ,1001 ,1)        ; 
  beautify(hist_ttjets,38 ,1001 ,1)        ; 
  beautify(hist_bkg   ,0  ,0    ,0)        ; 
  beautify(hist_sig0  ,43 ,1001 ,1)        ; 
  beautify(hist_sig1  ,46 ,1001 ,1)        ; 
  beautify(hist_sig2  ,49 ,1001 ,1)        ; 
  if ( name.Contains("nPVtx") || name.Contains("h_cutflow") ) {
    hist_bkg->SetFillStyle(1001);
    hist_bkg->SetFillColor(46);
  }
  else {
    hist_bkg->SetFillStyle(3254);
    hist_bkg->SetFillColor(12);
  }
  if (doData) {
    beautify(hist_data ,1 ,0 ,1) ; 
    hist_data->SetMarkerStyle(20);
    hist_data->SetMarkerSize(0.75);
    hist_data->SetLineWidth(2);
  }

  THStack *stack = new THStack("stack","");
  stack->Add(hist_ttjets) ;
  stack->Add(hist_qcd) ; 

  TH1D *hist_mcUnc, *hist_ratio ; 
  if (doData) {
    hist_mcUnc = (TH1D*)hist_bkg->Clone("hist_mcUnc") ; 
    hist_mcUnc->Sumw2() ; 
    hist_mcUnc->SetTitle("Stat. uncertainty on MC yield");
    hist_mcUnc->Divide(hist_bkg);

    hist_ratio = (TH1D*) hist_data->Clone("hist_ratio");
    hist_ratio->Sumw2() ; 
    hist_ratio->SetTitle("Data/MC");
    hist_ratio->Divide(hist_bkg);
  }

  TCanvas* c1 = new TCanvas();
  c1->cd();
  c1->SetBorderMode(0);
  c1->SetFrameBorderMode(0);

  TPad* pad0 = new TPad("pad0", "",0,0.30,1,.92) ; 
  pad0->Draw();
  pad0->cd();
  pad0->SetBorderMode(0);
  pad0->SetFrameBorderMode(0);
  beautifyTopPad(pad0) ; 
  pad0->SetLogy(log) ; 

  if (!log) {
    hist_bkg->SetMaximum( doData ? hist_data->GetMaximum()*2.1 : hist_bkg->GetMaximum()*2.1) ;
    hist_bkg->SetMinimum(0.) ; 
  }
  else {
    if (name.Contains("tau") || name.Contains("nFatJets") 
        || name.Contains("nJets") || name.Contains("TriggerSel_FatJets_Pt")
        || name.Contains("TriggerSel_SubJet1_Pt")) 
      hist_bkg->SetMaximum( doData ? hist_data->GetMaximum()*5000000 : hist_bkg->GetMaximum()*5000000) ;
    else 
      hist_bkg->SetMaximum( doData ? hist_data->GetMaximum()*2000 : hist_bkg->GetMaximum()*2000) ;
    hist_bkg->SetMinimum(0.1) ; 
  }

  TAxis* ax = hist_bkg->GetXaxis() ; 
  TAxis* ay = hist_bkg->GetYaxis() ; 
  beautifyAxis(ax) ; 
  beautifyAxis(ay) ; 
  ax->SetTitle(histotitle);
  if (doData) {
    ax->SetLabelSize( 0.0 );
    ax->SetTitleSize( 0.0 );
  }
  ay->SetTitle("Entries");
  hist_bkg->SetTitleOffset(0.83,"Y");

  if (setXRange) {
    if (rangeXLow == rangeXHigh) std::cout << "Error: X-axis low and high ranges have same value\n" ;
    else {
      hist_bkg->GetXaxis()->SetRangeUser(rangeXLow, rangeXHigh) ;
    }
  }

  if ( name.Contains("nPVtx") || name.Contains("h_cutflow") ) {
    hist_bkg->Draw("HIST");
    if (doData) hist_data->Draw("SAMEE");
  }
  else {
    hist_bkg->Draw("hist");
    stack->Draw("histSAME");
    if (doData) hist_data->Draw("SAMEE1");
    hist_bkg->Draw("samee2");
  }

  pad0->RedrawAxis();

  TPad *overlay ;
  if (name.Contains("nPVtx")) {
    overlay = new TPad("overlay","",0,0,1,1);
    overlay->SetFillStyle(4000);
    overlay->SetFillColor(0);
    overlay->SetFrameFillStyle(4000);
    overlay->Draw();
    overlay->cd();
    overlay->SetLogy();
    overlay->RedrawAxis() ;
  }

  int move_legend=0;
  TLegend *leg ;
  if (move_legend==1) {
    leg =  new TLegend(0.1,0.57,0.40,.92,NULL,"brNDC");
  }
  else {
    leg = new TLegend(0.57,0.62,0.895,0.93,NULL,"brNDC");
  }
  leg->SetBorderSize(1);
  leg->SetTextFont(132);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);

  if ( name.Contains("nPVtx") || name.Contains("h_cutflow") ) {
    leg->AddEntry(hist_bkg             , "All backgrounds", "f") ; 
  }
  else {
    std::cout << " name contains " << name << std::endl ; 
    leg->AddEntry(hist_ttjets          , "t#bar{t}+jets"          ,"f");
    leg->AddEntry(hist_qcd             , "Non-t#bar{t} multijets" ,"f");
    leg->AddEntry(hist_bkg             , "Bkg. error (stat.)"     ,"f");
  }
  if (doData) leg->AddEntry(hist_data, datacaption              ,"pl");

  leg->Draw();

  pad0->Modified();

  c1->cd();

  if (doData) {
    TPad* pad1 = new TPad("pad1", "",0,0,1,.25) ; 
    pad1->Draw() ;  
    pad1->cd() ;  
    pad1->SetGridx() ; 
    pad1->SetGridy() ; 

    hist_ratio->SetMarkerStyle(20);
    hist_ratio->SetMarkerSize(0.75);
    hist_ratio->SetLineWidth(2);

    hist_mcUnc->SetMarkerStyle(0);
    hist_mcUnc->SetMarkerSize(0);
    hist_mcUnc->SetLineWidth(0);
    hist_mcUnc->SetFillStyle(1001);
    hist_mcUnc->SetFillColor(kYellow);

    hist_mcUnc->GetYaxis()->SetTitle("Data/MC");
    hist_mcUnc->SetTitleOffset(0.9,"X");
    hist_mcUnc->SetTitleOffset(0.31,"Y");
    hist_mcUnc->GetXaxis()->SetTitle(histotitle);
    hist_mcUnc->GetYaxis()->SetNdivisions( 505 );

    TAxis* ax1 = hist_mcUnc->GetXaxis();
    TAxis* ay1 = hist_mcUnc->GetYaxis();

    beautifyBottomPad(pad1,ax1,ay1) ; 

    if (setXRange) {
      if (rangeXLow == rangeXHigh) std::cout << "Error: X-axis low and high ranges have same value\n" ;
      else {
        hist_mcUnc->GetXaxis()->SetRangeUser(rangeXLow, rangeXHigh) ;
      }
    }

    hist_mcUnc->SetMinimum(0.0);
    hist_mcUnc->SetMaximum(2.6);
    hist_mcUnc->Draw("E2");
    hist_ratio->Draw("SAMEE1");

    pad1->Modified();
  }

  c1->cd();

  char temp[100];
  TPaveText *plotlabel0 = new TPaveText(0.17,0.925,0.37,.95,"NDC");
  plotlabel0->SetTextColor(kBlack);
  plotlabel0->SetFillColor(kWhite);
  plotlabel0->SetBorderSize(0);
  plotlabel0->SetTextAlign(12);
  plotlabel0->SetTextSize(0.045);
  sprintf(temp, "%.1f", Lint/1000);
  plotlabel0->AddText( (string("CMS Preliminary 2012, ")  
        + string("L = ") 
        + temp 
        + string("/fb")).c_str()); 
  TPaveText *plotlabel1 = new TPaveText(0.77,0.925,0.90,.95,"NDC");
  plotlabel1->SetTextColor(kBlack);
  plotlabel1->SetFillColor(kWhite);
  plotlabel1->SetBorderSize(0);
  plotlabel1->SetTextAlign(12);
  plotlabel1->SetTextSize(0.045);
  plotlabel1->AddText("#sqrt{s} = 8 TeV") ; 
  TPaveText *plotlabel2 = new TPaveText(0.17,0.925,0.37,0.95,"NDC");
  plotlabel2->SetTextColor(kBlack);
  plotlabel2->SetFillColor(kWhite);
  plotlabel2->SetBorderSize(0);
  plotlabel2->SetTextAlign(12);
  plotlabel2->SetTextSize(0.045);
  sprintf(temp, "%.1f", Lint);
  plotlabel2->AddText((string("#int#font[12]{L}dt = ") + temp + string(" fb^{ -1}")).c_str()); 
  plotlabel0->Draw() ; 
  plotlabel1->Draw() ; 

  c1->Modified();
  c1->cd();
  c1->SetSelected(c1) ;

  TString name_plot=name+"_Linear"+formata;
  if(log) name_plot=name+"_Log"+formata;
  c1->SaveAs(dir4plots+"/"+name_plot);
  name_plot=name+"_Linear"+formatb;
  if(log) name_plot=name+"_Log"+formatb;
  c1->SaveAs(dir4plots+"/"+name_plot);
  name_plot=name+"_Linear"+formatc;
  if(log) name_plot=name+"_Log"+formatc;

  if (log && web) {  
    pad0 ->cd();
    pad0->SetLogy(false);
    c1->cd();
    c1->SaveAs(dir4plots+"/"+name+"_Linear"+formata);
  }

  return ; 

}

