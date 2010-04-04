import FWCore.ParameterSet.Config as cms


trackAnalyzer = cms.untracked.PSet(
    doMCTracks = cms.untracked.bool(True),
    trackTag   = cms.untracked.InputTag("generalTracks"),
    debugTrack = cms.untracked.int32(0)
    )
