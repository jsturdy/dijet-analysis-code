import FWCore.ParameterSet.Config as cms


photonAnalyzer = cms.untracked.PSet(
    photMaxEta = cms.untracked.double(12.4),
    photMaxEt  = cms.untracked.double(10000.),
    photMinEt  = cms.untracked.double(5.),
    photRelIso = cms.untracked.double(1.5),
    
    doMCPhots  = cms.untracked.bool(True),
    debugPhots = cms.untracked.int32(0),
    genPhotTag = cms.untracked.InputTag("genParticles"),
    
    photTag   = cms.untracked.InputTag("cleanLayer1Photons"),
    pfphotTag = cms.untracked.InputTag("pfLayer1Photons")
)
                 
