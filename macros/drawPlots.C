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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include "CMSstyle.C"
#include "help.C"

using namespace std;

TString filename = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_11_BpbH/src/BprimeTobHAnalysisv1/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_01_PURewt_CorrJetId/Final_histograms_BprimebH.root" ; 


TString title1 = "CMS Preliminary, 19.8 fb^{-1} at #sqrt{s} = 8 TeV";
TString datacaption = "Data"; 

TString dir4plots ="BprimeTobH_01_PURewt_CorrJetId" ;

TString formata = ".pdf";
TString formatb = ".png";
TString formatc = ".C";

bool web = 0;
bool setSampleName = 1;

void Draw(TString name, TString histotitle, bool log);
void DrawAll () ; 
void DrawStacked(TString name, TString histotitle, bool log, bool doData, bool fExtNorm=false, int nRebin=1, bool setXRange=false, double rangeXLow=0., double rangeXHigh=0.);

void drawPlots () {

  gROOT->SetBatch(kTRUE);
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TString action = "mkdir -p " + dir4plots;
  system(action);

  DrawAll () ; 

  return ; 

}

void DrawAll () {

   DrawStacked("TriggerSel_nPVtx_NoPUWt" ,"N(PV), No PU weight" ,1 ,1 ,0 ,1 ,1 ,0 ,50); 
   DrawStacked("TriggerSel_nPVtx_PUWt" ,"N(PV)" ,1 ,1 ,0 ,1 ,1 ,0 ,50); 
   //DrawStacked("HiggsJet_Pt" ,"p_{T}(Higgs jets)" ,1 ,1 ,0 ,1 ,1 ,0 ,1000); 
   //DrawStacked("nBJets" ,"N(b jets)" ,1 ,1 ,0 ,1 ,1 ,0 ,5); 

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

  //hist_qcd              = (TH1D*)myFile->Get("QCD__"+name+"_mc");
  //hist_ttjets           = (TH1D*)myFile->Get("TTJets__"+name+"_mc");
  //hist_sig0             = (TH1D*)myFile->Get("BprimeBprimeToBHBHinc_M-500__"+name+"_mc");
  //hist_sig1             = (TH1D*)myFile->Get("BprimeBprimeToBHBHinc_M-800__"+name+"_mc");
  //hist_sig2             = (TH1D*)myFile->Get("BprimeBprimeToBHBHinc_M-1000__"+name+"_mc");
  //if (doData) hist_data = (TH1D*)myFile->Get("DATA__"+name+"_data");

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
  hist_bkg->SetFillStyle(3004);
  hist_bkg->SetFillColor(12);
  if (doData) {
    beautify(hist_data ,1 ,0 ,1) ; 
    hist_data->SetMarkerStyle(20);
    hist_data->SetMarkerSize(0.75);
    hist_data->SetLineWidth(2);
  }

  THStack *stack = new THStack("stack","");
  stack->Add(hist_ttjets) ;
  stack->Add(hist_qcd) ; 

  TH1D *hist_ratio, *hist_pull ; 
  if (doData) {
    hist_ratio = (TH1D*) hist_data->Clone();
    hist_ratio->SetName("hist_ratio");
    hist_ratio->SetTitle("");
    hist_ratio->Divide(hist_bkg);
  }

  TCanvas *c1 = new TCanvas("c1", "c1",441,159,782,552);
  c1->Range(0,0,1,1);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);

  TPad* pad0 = new TPad("pad0", "pad0",0,0.25,1.0,1.00);
  pad0 ->Draw();
  pad0 ->cd();
  pad0->SetFillColor(0);
  pad0->SetBorderMode(0);
  pad0->SetBorderSize(2);
  pad0->SetFrameBorderMode(0);
  pad0->SetTopMargin(0.065);

  pad0->SetLogy(log);

  if (!log) {
    hist_bkg->SetMaximum( doData ? hist_data->GetMaximum()*2.1 : hist_bkg->GetMaximum()*2.1) ;
     hist_bkg->SetMinimum(0.) ; 
  }
  else {
    hist_bkg->SetMaximum( doData ? hist_data->GetMaximum()*100 : hist_bkg->GetMaximum()*100) ;
     hist_bkg->SetMinimum(0.1) ; 
  }

  hist_bkg->GetXaxis()->SetTitle(name);
  hist_bkg->GetYaxis()->SetTitle("Entries");
  hist_bkg->SetTitleOffset(0.81,"Y");
  hist_bkg->GetYaxis()->SetLabelSize( 0.05 );
  hist_bkg->GetYaxis()->SetTitleSize( 0.06 );

  if (setXRange) {
    if (rangeXLow == rangeXHigh) std::cout << "Error: X-axis low and high ranges have same value\n" ;
    else {
      hist_bkg->GetXaxis()->SetRangeUser(rangeXLow, rangeXHigh) ;
    }
  }

  if ( name.Contains("hist_bkg") ) {
    hist_bkg->Draw("HIST");
    if (doData) hist_data->Draw("SAMEE");
  }
  else {
    hist_bkg->Draw("hist");
    stack->Draw("histSAME");
    if (doData) hist_data->Draw("SAMEE");
    hist_bkg->Draw("samee2");
  }

  pad0->RedrawAxis();

  int move_legend=0;
  TLegend *leg ;
  if (move_legend==1) {
    leg =  new TLegend(0.1,0.55,0.40,.90,NULL,"brNDC");
  }
  else {
    leg = new TLegend(0.555,0.55,0.855,0.90,NULL,"brNDC");
  }
  leg->SetBorderSize(1);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  if ( !name.Contains("nPVtx") ) {
    leg->AddEntry(hist_ttjets          , "t#bar{t}+jets"          ,"f");
    leg->AddEntry(hist_qcd             , "Non-t#bar{t} multijets" ,"f");
  }
  leg->AddEntry(hist_bkg             , "Bkg. error (stat.)"     ,"f");
  if (doData) leg->AddEntry(hist_data, datacaption              ,"pl");

  leg->Draw();

  TLatex *   tex0 = new TLatex(0.42,1.00,"CMS Preliminary, 19.8 fb^{-1} at #sqrt{s} = 8 TeV");
  tex0->SetNDC();
  tex0->SetTextAlign(13);
  tex0->SetTextFont(42);
  tex0->SetTextFont(62);
  tex0->SetTextSize(0.055);
  tex0->SetLineWidth(2);
  tex0->Draw();

  if (setSampleName) {
    TString sample = "";
    TLatex *tex1 = new TLatex(0.14,0.88,sample);
    tex1->SetNDC();
    tex1->SetTextAlign(13);
    tex1->SetTextFont(42);
    tex1->SetTextFont(62);
    tex1->SetTextSize(0.055);
    tex1->SetLineWidth(2);
    tex1->Draw();
  }

  pad0->Modified();

  c1->cd();

  if (doData) {
    TPad* pad1 = new TPad("pad1", "pad1",0,0.,1.0,0.32);
    pad1->Draw();
    pad1->cd();
    pad1->SetFillColor(0);
    pad1->SetBorderMode(0);
    pad1->SetBorderSize(2);
    pad1->SetGridy();
    pad1->SetBottomMargin(0.31);
    pad1->SetFrameBorderMode(0);

    hist_ratio->SetMarkerStyle(20);
    hist_ratio->SetMarkerSize(0.75);
    hist_ratio->SetLineWidth(2);

    hist_ratio->GetYaxis()->SetTitle("Data/MC");
    hist_ratio->SetTitleOffset(0.9,"X");
    hist_ratio->SetTitleOffset(0.31,"Y");
    hist_ratio->GetXaxis()->SetTitle(histotitle);
    hist_ratio->GetYaxis()->SetNdivisions( 505 );

    double labelsizex=0.12;
    double labelsizey=0.12;
    double titlesizex=0.15;
    double titlesizey=0.14;

    hist_ratio->GetXaxis()->SetLabelSize( labelsizex );
    hist_ratio->GetXaxis()->SetTitleSize( titlesizex );
    hist_ratio->GetYaxis()->SetLabelSize( labelsizey );
    hist_ratio->GetYaxis()->SetTitleSize( titlesizey );

    if (setXRange) {
      if (rangeXLow == rangeXHigh) std::cout << "Error: X-axis low and high ranges have same value\n" ;
      else {
        hist_ratio->GetXaxis()->SetRangeUser(rangeXLow, rangeXHigh) ;
      }
    }

    hist_ratio->SetMinimum(0.0);
    hist_ratio->SetMaximum(2.6);
    hist_ratio->Draw("E1X0");

    pad1->Modified();
  }

  c1->cd();
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
