import FWCore.ParameterSet.Config as cms

process = cms.Process("BprimebH")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) ) # Leave it this way. 

process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("file.root")
)

from Bprime_kit.BprimeTobHAnalysis.bprimetobhanalysis_cfi import *

process.BprimebH = cms.EDAnalyzer('BprimeTobHAnalysis',
    InputFiles        = cms.string("/tmp/petrakou/NtuplesBprimeTobH_BprimeBprimeToBHTWinc_M-1500_TuneZ2star_8TeV-madgraph_BprimeTobH_3_1_xWa.root"),
    MaxEvents         = cms.int32(100),
    OutputFile       = cms.untracked.string("file.root")
)

process.p = cms.Path(process.BprimebH)

