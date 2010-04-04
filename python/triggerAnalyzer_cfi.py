import FWCore.ParameterSet.Config as cms


triggerAnalyzer = cms.untracked.PSet(
    triggerResults = cms.untracked.InputTag("triggerResults"),
    pathNames      = cms.untracked.vstring('HLT_Jet180','HLT_DiJetAve130','HLT_MET50','HLT_Mu9'),
    debugTrig      = cms.untracked.int32(0)
    )
