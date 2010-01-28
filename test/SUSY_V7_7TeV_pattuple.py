#
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 33X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV7
#

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#-- Meta data to be logged in DBS ---------------------------------------------
process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.23 $'),
            name = cms.untracked.string('$Source: /local/projects/CMSSW/rep/CMSSW/PhysicsTools/Configuration/test/SUSY_pattuple_cfg.py,v $'),
            annotation = cms.untracked.string('SUSY pattuple definition')
        )

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
        limit = cms.untracked.int32(-1),
            reportEvery = cms.untracked.int32(1)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
process.GlobalTag.globaltag = 'MC_3XY_V15::All' #Data: GR09_P_V8_34X , MC: MC_3XY_V15
#process.GlobalTag.globaltag = 'STARTUP31X_V2::All' #Data: GR09_P_V8_34X , MC: MC_3XY_V15

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands
#Apply SUSYPAT: Parameters are: mcInfo, HLT menu, Jet energy corrections, MC version ('31x' or '31xReReco332')
addDefaultSUSYPAT(process,True,'HLT','Summer09_7TeV','31x')
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
############################## END SUSYPAT specifics ####################################


#-- Output module configuration -----------------------------------------------
process.out.fileName = 'SUSYPAT_7TeV.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )

# Analyzer

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
        "7TeV_PATtified.root"
      )
)
# source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/mc/Summer09/LM0/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0005/2055D8E0-83BD-DE11-8EF7-001E0B1CA8A0.root'
    )
                            )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.jaredsusy = cms.EDAnalyzer("JaredSusyAnalyzer",
#process.jaredsusy = cms.EDFilter("JaredSusyAnalyzer",
     UsePfjet = cms.bool(True),
     pfjetTag = cms.InputTag("cleanLayer1JetsAK5PF"),
     jptTag   = cms.InputTag("cleanLayer1JetsAK5JPT"),
     jetTag   = cms.InputTag("cleanLayer1JetsAK5") ,

     genTag     = cms.InputTag("genParticles"),
     doMCData   = cms.bool(True),
     jetMaxEta  = cms.double(5.0),
     jetMinPt   = cms.double(30.0),

     vtxTag  = cms.InputTag("offlinePrimaryVertices"),
     elecTag = cms.InputTag("cleanLayer1Electrons"),
     muonTag = cms.InputTag("cleanLayer1Muons"),

     pfelecTag = cms.InputTag("pfLayer1Electrons"),
     pfmuonTag = cms.InputTag("pfLayer1Muons"),

     tcmetTag = cms.InputTag("layer1METsTC"),
     pfmetTag = cms.InputTag("pfLayer1METs"),
     metTag   = cms.InputTag("layer1METsAK5"),
     #mhtTag   = cms.InputTag("layer1MHTsAK5"),
     

     pathNames = cms.vstring('HLT_Jet180','HLT_DiJetAve130','HLT_MET50','HLT_Mu9') ,    
                   
     eventWeight = cms.double(19.67249),
  
     triggerResults = cms.InputTag("TriggerResults","","HLT"),
)

process.selectedLayer2Hemispheres = cms.EDProducer("PATHemisphereProducer",
    patpfJets    = cms.InputTag("cleanLayer1JetsAK5PF"),
    patJets      = cms.InputTag("cleanLayer1JetsAK5"),
    patMuons     = cms.InputTag("cleanLayer1Muons"),
    patElectrons = cms.InputTag("cleanLayer1Electrons"),
    patPhotons   = cms.InputTag("cleanLayer1Photons"),
    patTaus      = cms.InputTag("cleanLayer1Taus") ,
                                                   
    patMets   = cms.InputTag("allLayer1METsAK5"),
                                                   
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
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#-- Execution path ------------------------------------------------------------


# Full path
process.p = cms.Path( process.seqSUSYDefaultSequence
                      * process.selectedLayer2Hemispheres*process.jaredsusy
                      )

