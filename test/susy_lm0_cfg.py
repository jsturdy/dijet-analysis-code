import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

# Include PAT Layer 0 & 1 if not running on pattified data
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V11::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# CaloTowerConstituentsMap needed for Electron/Photon-Jet cleaning
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

# Cross-cleaner setup
process.load("SusyAnalysis.PatCrossCleaner.patCrossCleaner_cfi")

process.patcrosscleaner.patMuons      = cms.InputTag("allLayer1Muons")
process.patcrosscleaner.patElectrons  = cms.InputTag("allLayer1Electrons")
process.patcrosscleaner.patPhotons    = cms.InputTag("allLayer1Photons")
#process.patcrosscleaner.patJets       = cms.InputTag("allLayer1JetsPF")
process.patcrosscleaner.patJets       = cms.InputTag("allLayer1JetsSC5")
#process.patcrosscleaner.patMets       = cms.InputTag("allLayer1METstcMET")
#process.patcrosscleaner.patMets       = cms.InputTag("allLayer1METsPF")
process.patcrosscleaner.patMets       = cms.InputTag("allLayer1METsSC5")
# Switch on/off some components
process.patcrosscleaner.doMuonJetCC        = True
process.patcrosscleaner.doElectronJetCC    = True
process.patcrosscleaner.doPhotonJetCC      = True
process.patcrosscleaner.doElectronPhotonCC = True
# Change the jet energy corrections
process.patcrosscleaner.L1JetCorrector      = 'none'
process.patcrosscleaner.L2JetCorrector      = 'L2RelativeJetCorrectorIC5Calo'
process.patcrosscleaner.L3JetCorrector      = 'L3AbsoluteJetCorrectorIC5Calo'
process.patcrosscleaner.L4JetCorrector      = 'none'
process.patcrosscleaner.L5udsJetCorrector   = 'none'
process.patcrosscleaner.L5gluonJetCorrector = 'none'
process.patcrosscleaner.L5cJetCorrector     = 'none'
process.patcrosscleaner.L5bJetCorrector     = 'none'
#process.patcrosscleaner.L6udsJetCorrector   = 'none'
#process.patcrosscleaner.L6gluonJetCorrector = 'none'
#process.patcrosscleaner.L6cJetCorrector     = 'none'
#process.patcrosscleaner.L6bJetCorrector     = 'none'
process.patcrosscleaner.L6JetCorrector      = 'none'
process.patcrosscleaner.L7udsJetCorrector   = 'L7PartonJetCorrectorIC5qJet'
process.patcrosscleaner.L7gluonJetCorrector = 'L7PartonJetCorrectorIC5gJet'
process.patcrosscleaner.L7cJetCorrector     = 'L7PartonJetCorrectorIC5cJet'
process.patcrosscleaner.L7bJetCorrector     = 'L7PartonJetCorrectorIC5bJet'

#Parameters for electron-jet cross-cleaning
#process.patcrosscleaner.ElectronJetCrossCleaning.SusyAnalyzerCleaning = True
process.patcrosscleaner.ElectronJetCrossCleaning.deltaR_min        = 0.5
process.patcrosscleaner.ElectronJetCrossCleaning.SharedEtoJetE     = 0.7
#process.patcrosscleaner.ElectronJetCrossCleaning.SharedEForNIsoEle = -1.
#process.patcrosscleaner.ElectronJetCrossCleaning.IsolationKey  = 'TrackerIso'
process.patcrosscleaner.ElectronJetCrossCleaning.IsolationKey = cms.string('CombRelIso')
#process.patcrosscleaner.ElectronJetCrossCleaning.IsoValueCut   = 1.
process.patcrosscleaner.ElectronJetCrossCleaning.IsoValueCut = 0.2
process.patcrosscleaner.ElectronJetCrossCleaning.ecalIsoWeight = 1.
process.patcrosscleaner.ElectronJetCrossCleaning.hcalIsoWeight = 1.
process.patcrosscleaner.ElectronJetCrossCleaning.trkIsoWeight = 1.

process.patcrosscleaner.ElectronJetCrossCleaning.ElectronID   = 'eidRobustLoose'
# Parameters for photon-jet cross-cleaning
process.patcrosscleaner.PhotonJetCrossCleaning.deltaR_min   = 0.5
process.patcrosscleaner.PhotonJetCrossCleaning.IsoValueCut  = 0.3
process.patcrosscleaner.PhotonJetCrossCleaning.IsolationKey = 'CaloIso'
process.patcrosscleaner.PhotonJetCrossCleaning.PhotonID = 'TightPhoton'
# Parameters for muon-jet cross-cleaning
process.patcrosscleaner.MuonJetCrossCleaning.deltaR_min   = 0.2
#process.patcrosscleaner.MuonJetCrossCleaning.caloIso_max  = 10.0
#process.patcrosscleaner.MuonJetCrossCleaning.trackIso_max = 10.0
process.patcrosscleaner.MuonJetCrossCleaning.useCombRelIso = True # ntuple 3 + // V+jets recommodation
process.patcrosscleaner.MuonJetCrossCleaning.IsoValueCut = 0.2
process.patcrosscleaner.MuonJetCrossCleaning.ecalIsoWeight = 1.
process.patcrosscleaner.MuonJetCrossCleaning.hcalIsoWeight = 1.
process.patcrosscleaner.MuonJetCrossCleaning.trkIsoWeight = 1.

#process.patcrosscleaner.MuonJetCrossCleaning.MuonID = 'TM2DCompatibilityTight'
process.patcrosscleaner.MuonJetCrossCleaning.MuonID = 'AllGlobalMuons'

process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
      "SUSY_LM0_SC5MET_PATtified.root"
      )
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      "file:/uscms/home/sturdy07/working/LM0_patLayer1.root"
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.jaredsusy = cms.EDAnalyzer("JaredSusyAnalyzer",
#process.jaredsusy = cms.EDFilter("JaredSusyAnalyzer",
     UsePfjet = cms.bool(True),
     pfjetTag = cms.InputTag("allLayer1JetsIC5PF"),
     jptTag   = cms.InputTag("allLayer1JetsSC5"),
     jetTag   = cms.InputTag("allLayer1JetsSC5") ,

     genTag     = cms.InputTag("genParticles"),
     jetMaxEta  = cms.double(3.0),
     jetMinPt   = cms.double(30.0),

     vtxTag  = cms.InputTag("offlinePrimaryVertices"),
     tauTag  = cms.InputTag("allLayer1Taus"),
     elecTag = cms.InputTag("allLayer1Electrons"),
     photTag = cms.InputTag("allLayer1Photons"),
     muonTag = cms.InputTag("allLayer1Muons"),

     tcmetTag = cms.InputTag("allLayer1METstcMET"),
     #pfmetTag = cms.InputTag("allLayer1METsPF"),
     metTag   = cms.InputTag("allLayer1METsSC5"),
     
     ccelecTag   = cms.InputTag("patcrosscleaner:ccElectrons"),
     ccJptTag    = cms.InputTag("patcrosscleaner:ccJets"),
     ccpfjetTag  = cms.InputTag("patcrosscleaner:ccJets"),
     ccjetTag    = cms.InputTag("patcrosscleaner:ccJets"),
     ccmuonTag   = cms.InputTag("patcrosscleaner:ccMuons"),
     #ccpfmetTag  = cms.InputTag("patcrosscleaner:ccMETsPF"),
     #cctcmetTag  = cms.InputTag("patcrosscleaner:ccMETsTC"),
     ccmetTag    = cms.InputTag("patcrosscleaner:ccMETs"),
     ccphotonTag = cms.InputTag("patcrosscleaner:ccPhotons"),

    pathNames = cms.vstring('HLT_Jet180','HLT_DiJetAve130','HLT_MET50','HLT_Mu9') ,    
                   
    eventWeight = cms.double(19.67249),
  
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
)

#process.genParticles.abortOnUnknownPDGCode = False

process.selectedLayer2Hemispheres = cms.EDProducer("PATHemisphereProducer",
    #patpfJets    = cms.InputTag("allLayer1JetsIC5PF"),
    patJets      = cms.InputTag("allLayer1JetsSC5"),
    patMuons     = cms.InputTag("allLayer1Muons"),
    patElectrons = cms.InputTag("allLayer1Electrons"),
    patPhotons   = cms.InputTag("allLayer1Photons"),
    patTaus      = cms.InputTag("allLayer1Taus") ,
                                                   
    #pattcMets = cms.InputTag("allLayer1METstcMET"),
    #patpfMets = cms.InputTag("allLayer1METsPF"),
    patMets   = cms.InputTag("allLayer1METsSC5"),
                                                   
    maxElectronEta = cms.double(5.0),
    maxTauEta = cms.double(-1.0),
    maxPhotonEta = cms.double(5.0), 
    maxMuonEta = cms.double(5.0),
    maxJetEta = cms.double(5.0),
                                                   
    minMuonEt = cms.double(1000000.0),#comment out muons                                        
    minTauEt = cms.double(1000000.0),#comment out taus
    minPhotonEt = cms.double(100000.0),#comment out photons
    minElectronEt = cms.double(1000000.0),#comment out electrons
    minJetEt = cms.double(30.0),
                                                   
    combinationMethod = cms.int32(3),
    seedMethod = cms.int32(3)
  
)


process.MessageLogger.categories.extend(['JaredSusyEvent'])
process.MessageLogger.cerr.threshold = 'WARN'
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.p = cms.Path(process.patcrosscleaner*process.selectedLayer2Hemispheres*process.jaredsusy)
#process.p = cms.Path(process.patLayer0*process.patLayer1*process.patcrosscleaner*process.selectedLayer2Hemispheres*process.jaredsusy)
#process.p = cms.Path(process.patLayer0*process.patLayer1*process.patcrosscleaner*process.selectedLayer2Hemispheres*process.jaredsusy)


## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)
