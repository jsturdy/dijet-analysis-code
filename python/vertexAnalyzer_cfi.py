import FWCore.ParameterSet.Config as cms


vertexAnalyzer = cms.untracked.PSet(
    minNVtx    = cms.untracked.int32(1),
    minVtxTrks = cms.untracked.int32(3),
    minVtxNdof = cms.untracked.int32(4),
    maxVtxChi2 = cms.untracked.double(999),
    maxVtxZ    = cms.untracked.double(15.),

    vtxTag = cms.untracked.InputTag("offlinePrimaryVertices"),
    debugVtx = cms.untracked.int32(0)
    )
