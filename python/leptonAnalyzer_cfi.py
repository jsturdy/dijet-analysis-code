import FWCore.ParameterSet.Config as cms


leptonAnalyzer = cms.untracked.PSet(
    elecMaxEta = cms.untracked.double(2.4),
    elecMaxEt  = cms.untracked.double(15.),
    elecMinEt  = cms.untracked.double(5.),
    elecRelIso = cms.untracked.double(0.5),

    muonMaxEta = cms.untracked.double(2.),
    muonMaxEt  = cms.untracked.double(10.),
    muonMinEt  = cms.untracked.double(5.),
    muonRelIso = cms.untracked.double(0.1),
    
    #tauMaxEta = cms.untracked.double(2.),
    #tauMaxEt  = cms.untracked.double(10.),
    #tauMinEt  = cms.untracked.double(5.),
    #tauRelIso = cms.untracked.double(0.1),
    
    doMCLeps  = cms.untracked.bool(True),
    debugLeps = cms.untracked.int32(0),
    genLepTag = cms.untracked.InputTag("genParticles"),
    
    elecTag   = cms.untracked.InputTag("cleanLayer1Electrons"),
    pfelecTag = cms.untracked.InputTag("pfLayer1Electrons"),
    muonTag   = cms.untracked.InputTag("cleanLayer1Muons"),
    pfmuonTag = cms.untracked.InputTag("pfLayer1Muons")
    #tauTag   = cms.untracked.InputTag("cleanLayer1Taus"),
    #pftauTag = cms.untracked.InputTag("pfLayer1Taus")
)
                 
