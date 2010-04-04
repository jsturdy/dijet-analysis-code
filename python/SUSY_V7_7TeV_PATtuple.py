#
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 35X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV8
#

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#-- Meta data to be logged in DBS ---------------------------------------------
process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.26 $'),
            name = cms.untracked.string('$Source: /local/projects/CMSSW/rep/CMSSW/PhysicsTools/Configuration/test/SUSY_pattuple_cfg.py,v $'),
            annotation = cms.untracked.string('SUSY pattuple definition')
        )

#-- Message Logger ------------------------------------------------------------
#process.MessageLogger.categories.append('PATSummaryTables')
#process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1),
#            reportEvery = cms.untracked.int32(1000)
#            )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#-- Input Source --------------------------------------------------------------
process.source.fileNames = [
    #'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/n/nmohr/QCDDiJet_Pt380to470_MC_31X_V9_ReReco332.root'
    #'rfio://?svcClass=cmscafuser&path=/castor/cern.ch/cms/store/caf/user/fronga/V6production/PYTHIA6_SUSY_LM0_sftsht_10TeV_cff_py_RAW2DIGI_RECO_1.root'
    #'/store/mc/Summer09/QCD_Pt30/ALCARECO/MC_31X_V3_7TeV_StreamHcalCalDijets-v1/0000/FA10893E-66AF-DE11-B45C-001A4BA939F2.root'
    #'/store/mc/Summer09/Zprime_M500GeV_W5GeV-madgraph/GEN-SIM-RECO/MC_31X_V3_7TeV-v3/0000/F8CB67AA-DFEB-DE11-ACEB-001A644EAF10.root'
    '/store/relval/CMSSW_3_4_1/RelValTTbar/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/CE62D4D8-85ED-DE11-8BD2-000423D9853C.root'
    #'/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/SD_InterestingEvents-Dec19thSkim_341_v1/0006/E8C2D2A4-B9ED-DE11-A4E2-0026189438A5.root'     
    ]
process.maxEvents.input = 100
# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

##run on 31X mc samples
#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
#run33xOn31xMC(process)

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
process.GlobalTag.globaltag = 'MC_3XY_V18::All'

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands
#Apply SUSYPAT, parameters are: mcInfo, HLT menu, Jet energy corrections, JetCollections
addDefaultSUSYPAT(process,True,'HLT','Summer09_7TeV_ReReco332','31xReReco332',['AK5PF','AK5JPT','AK5Track'])
#SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
############################### END SUSYPAT specifics ####################################
#
#
##-- Output module configuration -----------------------------------------------
#process.out.fileName = 'SUSYPAT_7TeV.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS
#
## Custom settings
#process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
#process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
#process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
#process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )
del process.outpath
# Analyzer section
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
        #"QCD_7TeV_PATtified.root"
        #"ZPrime_7TeV_PATtified.root"
        "RelValTTbar_7TeV_PATtified.root"
    )
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#fixed met sequence
process.load("JetMETCorrections.Type1MET.MuonMETValueMapProducer_cff")
process.load("JetMETCorrections.Type1MET.MuonTCMETValueMapProducer_cff")
process.load("JetMETCorrections.Type1MET.MetMuonCorrections_cff")
process.load("RecoMET.METProducers.TCMET_cfi")
process.fixedMETSequence = cms.Sequence(
    process.muonMETValueMapProducer *
        process.muonTCMETValueMapProducer *
        process.corMetGlobalMuons *
        process.tcMet
)


#from JSturdy.DiJetAnalysis.diJetAnalyzer_cff import doDiJetAnalysis
process.load("JSturdy.DiJetAnalysis.diJetAnalyzer_cff")


#-- Execution path ------------------------------------------------------------
# Full path
process.p = cms.Path(
    process.fixedMETSequence *
    process.seqSUSYDefaultSequence *
    process.doDiJetAnalysis)
#-- Dump config ------------------------------------------------------------
file = open('SusyPAT_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
