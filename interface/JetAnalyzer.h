#ifndef JETANALYZER
#define JETANALYZER

// System include files
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

// ROOT includes
#include <TNtuple.h>

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// SUSY include files
//#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
//#include "SusyAnalysis/EventSelector/interface/SelectorSequence.h"

//#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

//#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

//#include "UserCode/AnalysisTools/test/ALPGENParticleId.cc"

//#include "PhysicsTools/Utilities/interface/deltaPhi.h"
//#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Hemisphere.h"


//
// Class declaration
//

class JetAnalyzer {
 public:
  JetAnalyzer(const edm::ParameterSet&);
  ~JetAnalyzer();
  
 private:
  //*** CMSSW interface
  /// Called once per job, at start
  void beginJob(const edm::EventSetup&) ;
  /// Called for each event
  bool filter(const edm::Event&, const edm::EventSetup&);
  /// Called once per job, at end
  void endJob();
  
  /// Print a summary of counts for all selectors
  //  virtual void printSummary(void);
  // Print an HLT trigger report
  void printHLTreport(void); // georgia
  
  //*** Plotting
  /// Define all plots
  void initTuple();
  void bookHealthPlots();
  /// Fill all plots for an event
  void fillTuple( const edm::Event& ) {mJetData->Fill(); }
  //void fillHealthPlots( const edm::Event& ) {mJetData->Fill(); }

  //  virtual bool filter(const edm::Event& evt,const edm::EventSetup& iSetup );


 private:

  bool matchJetsByCaloTowers( const pat::Jet&, const pat::Jet& );

  // Configuration parameters
  double jetMaxEta_, jetMinPt_;  /// for preselection cuts on jets and to calculate HT and MHT
  bool doMCData_;
  bool usePfJets_;
  bool useJPTJets_;
  bool useCaloJets_;

  // Data tags
  edm::InputTag pfjetTag_;
  edm::InputTag jptTag_;
  edm::InputTag jetTag_;
  edm::InputTag genTag_;
  edm::InputTag mhtTag_;
  edm::InputTag htTag_;
    
  // Plots
  TNtuple* ntuple_;      /// Will contain all the selector information we want to keep
  TTree * mAllData;      /// Will contain the additional di-jet specific data
  TTree * mJetData;

  //PF Jets
  int    m_NPFJets;
  double m_PFHt;
  double m_PFMHx;
  double m_PFMHy;
  double m_PFMHt;
  double m_PFJetEta[50];
  double m_PFJetPhi[50];
  double m_PFJetE[50];
  double m_PFJetEt[50];
  double m_PFJetPx[50];
  double m_PFJetPy[50];
  double m_PFJetPz[50];
  double m_PFJetPt[50];
  double m_PFJetCharge[50];
  double m_PFJetFem[50];

  //Calo Jets no JPT corrections
  int    m_NCaloJets;
  double m_CaloHt;
  double m_CaloMHx;
  double m_CaloMHy;
  double m_CaloMHt;
  double m_CaloJetEt[50];
  double m_CaloJetPt[50];
  double m_CaloJetPx[50];
  double m_CaloJetPy[50];
  double m_CaloJetPz[50];
  double m_CaloJetE[50];
  double m_CaloJetEta[50];
  double m_CaloJetPhi[50];
  double m_CaloJetFem[50];
  int    m_CaloJetHemi[50];

  double m_CaloJetMCCorrFactor[50];
  double m_CaloJetJPTCorrFactor[50];

  // track info:
  int    m_JetTrackNo[50];
  double m_JetTrackPhi[50];
  double m_JetTrackPhiWeighted[50];
  double m_JetTrackPt[50];

  double m_ccJetMCCorrFactor[50];
  double m_ccJetJPTCorrFactor[50];

  bool    m_ccJetAssoc[50];
  double  m_ccJetAssoc_E[50];
  double  m_ccJetAssoc_px[50];
  double  m_ccJetAssoc_py[50];
  double  m_ccJetAssoc_pz[50];

  int    m_JetPartonId[50];
  int    m_JetPartonMother[50];
  int    m_JetPartonFlavour[50];
  double m_JetPartonPx[50];
  double m_JetPartonPy[50];
  double m_JetPartonPz[50];
  double m_JetPartonEt[50];
  double m_JetPartonEnergy[50];
  double m_JetPartonPhi[50];
  double m_JetPartonEta[50];

  float m_JetsBTag_TkCountHighEff[50];
  float m_JetsBTag_SimpleSecVtx[50];
  float m_JetsBTag_CombSecVtx[50];

  //Generator level information
  double m_GenHt;
  double m_GenMHx;
  double m_GenMHy;
  double m_GenMHt;
  double m_GenJetEt[50];
  double m_GenJetPt[50];
  double m_GenJetE[50];
  double m_GenJetPx[50];
  double m_GenJetPy[50];
  double m_GenJetPz[50];
  double m_GenJetEta[50];
  double m_GenJetPhi[50];


  //JPT corrected Calo Jets
  int    m_NJPTJets;
  double m_JPTHt;
  double m_JPTMHx;
  double m_JPTMHy;
  double m_JPTMHt;
  double m_JPTJetEt[50];
  double m_JPTJetPt[50];
  double m_JPTJetPx[50];
  double m_JPTJetPy[50];
  double m_JPTJetPz[50];
  double m_JPTJetE[50];
  double m_JPTJetEta[50];
  double m_JPTJetPhi[50];
  double m_JPTJetFem[50];
  int    m_JPTJetHemi[50];
  int    m_JPTJetPartonFlavour[50];

  double m_JPTJetMCCorrFactor[50];
  double m_JPTJetJPTCorrFactor[50];


  std::string outputFileName_;

  double localPi;
  //unsigned int *mSelectorResults;

};
#endif
