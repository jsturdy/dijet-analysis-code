import FWCore.ParameterSet.Config as cms


jetAnalyzer = cms.untracked.PSet(

    jetMaxEta = cms.untracked.double(3.0),
    jetMinPt  = cms.untracked.double(30.),
    jetMaxEMF = cms.untracked.double(0.95),
    jetMinEMF = cms.untracked.double(0.),

    selJetMaxEta = cms.untracked.vdouble(2.5,3.),
    selJetMinPt  = cms.untracked.vdouble(50.,50.),
    selJetMaxEMF = cms.untracked.vdouble(0.95,0.95),
    selJetMinEMF = cms.untracked.vdouble(0.05,0.05),
    doMCJets     = cms.untracked.bool(True),
    genJetTag    = cms.untracked.InputTag("ak5GenJets"),

    usePFJets    = cms.untracked.bool(True),
    useJPTJets   = cms.untracked.bool(True),
    useCaloJets  = cms.untracked.bool(True),
    useTrackJets = cms.untracked.bool(True),
    pfJetTag     = cms.untracked.InputTag("cleanLayer1JetsAK5PF"),
    jptJetTag    = cms.untracked.InputTag("cleanLayer1JetsAK5JPT"),
    caloJetTag   = cms.untracked.InputTag("cleanLayer1JetsAK5"),
    trackJetTag  = cms.untracked.InputTag("cleanLayer1JetsAK5Track"),

    #htTag        = cms.InputTag("htTag"),
    #mhtTag       = cms.InputTag("mhtTag"),

    debugJets = cms.untracked.int32(0)

    )
