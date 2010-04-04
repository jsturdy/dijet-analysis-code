import FWCore.ParameterSet.Config as cms


metAnalyzer = cms.untracked.PSet(
    doMCMET   = cms.untracked.bool(True),
    doCaloMET = cms.untracked.bool(True),
    doPfMET   = cms.untracked.bool(True),
    doTcMET   = cms.untracked.bool(True),
    
    genMETTag = cms.untracked.InputTag("genMetCalo"),
    metTag    = cms.untracked.InputTag("layer1METsAK5"),
    pfmetTag  = cms.untracked.InputTag("layer1METsPF"),
    tcmetTag  = cms.untracked.InputTag("layer1METsTC"),

    debugMET = cms.untracked.int32(0)
    )
