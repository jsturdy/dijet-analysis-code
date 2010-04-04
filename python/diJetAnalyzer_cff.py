import FWCore.ParameterSet.Config as cms


from JSturdy.DiJetAnalysis.diJetAnalyzer_cfi import *
doDiJetAnalysis = cms.Sequence(
    diJetAnalyzer
    )
