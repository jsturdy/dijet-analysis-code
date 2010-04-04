import FWCore.ParameterSet.Config as cms


from JSturdy.DiJetAnalysis.jetAnalyzer_cfi import jetAnalyzer
from JSturdy.DiJetAnalysis.metAnalyzer_cfi import metAnalyzer
from JSturdy.DiJetAnalysis.leptonAnalyzer_cfi import leptonAnalyzer
from JSturdy.DiJetAnalysis.trackAnalyzer_cfi import trackAnalyzer
from JSturdy.DiJetAnalysis.triggerAnalyzer_cfi import triggerAnalyzer
from JSturdy.DiJetAnalysis.vertexAnalyzer_cfi import vertexAnalyzer
from JSturdy.DiJetAnalysis.photonAnalyzer_cfi import photonAnalyzer
#from JSturdy.DiJetAnalysis.Analyzer_cfi import 

diJetAnalyzer = cms.EDAnalyzer("DiJetAnalysis",
    jetParameters        = cms.untracked.PSet(jetAnalyzer.clone()),
   #hemisphereParameters = cms.untracked.PSet(hemisphereAnalyzer.clone()),
    leptonParameters     = cms.untracked.PSet(leptonAnalyzer.clone()),
    photonParameters     = cms.untracked.PSet(photonAnalyzer.clone()),
    metParameters        = cms.untracked.PSet(metAnalyzer.clone()),
    vertexParameters     = cms.untracked.PSet(vertexAnalyzer.clone()),
    trackParameters      = cms.untracked.PSet(trackAnalyzer.clone()),
    triggerParameters    = cms.untracked.PSet(triggerAnalyzer.clone()),
    debugDiJets = cms.untracked.int32(0)
)
