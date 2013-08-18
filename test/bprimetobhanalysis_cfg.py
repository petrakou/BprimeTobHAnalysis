import FWCore.ParameterSet.Config as cms

process = cms.Process("BprimebH")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cout = cms.untracked.PSet(
#		threshold = cms.untracked.string('INFO'), 
#		) 
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) ) # Leave it this way. 

process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("file.root")
)

from BprimebHAnalysis.BprimeTobHAnalysis.bprimetobhanalysis_cfi import *

process.BprimebH = cms.EDAnalyzer('BprimeTobHAnalysis',
		InputFiles = cms.string("dcache:/pnfs/cms/WAX/11/store/user/devdatta/NtuplesBprimeTobH/BprimeBprimeToBHTWinc_M-1500_TuneZ2star_8TeV-madgraph/BprimeTobH_18_1_ghb.root"),
		MaxEvents  = cms.int32(-1),
		FatJetPtMin = cms.double(150.),
		FatJetPtMax = cms.double(1.E6),
		FatJetAbsEtaMax = cms.double(2.5),
                FatJetMassMin = cms.double(100.),
		FatJetMassMax = cms.double(150),
		FatJetPrunedMassMin = cms.double(75.),
		FatJetPrunedMassMax = cms.double(1.E6),
		FatJetTau2ByTau1Min = cms.double(0.5),
		Subjet1CSVDiscMin = cms.double(0.679),
		Subjet1CSVDiscMax = cms.double(1.),
		Subjet2CSVDiscMin = cms.double(0.679),
		Subjet2CSVDiscMax = cms.double(1.),
		) 

process.p = cms.Path(process.BprimebH)

