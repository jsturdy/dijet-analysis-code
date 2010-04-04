#ifndef PHOTONANALYZER
#define PHOTONANALYZER

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

#include "DataFormats/PatCandidates/interface/Photon.h"

//#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"

//#include "UserCode/AnalysisTools/test/ALPGENParticleId.cc"


//
// Class declaration
//
class PhotonAnalyzer {
 public:
  PhotonAnalyzer(const edm::ParameterSet&, TTree*);
  ~PhotonAnalyzer();
  
  //*** CMSSW interface
  /// Called once per job, at start
  void beginJob(const edm::EventSetup&) ;
  /// Called for each event
  //void analyze(const edm::Event&, const edm::EventSetup&);
  bool filter(const edm::Event& evt,const edm::EventSetup& iSetup );
  /// Called once per job, at end
  void endJob();
  
  
  //*** Plotting
  /// Define all plots
  void initTuple();
  //void bookHealthPlots();
  /// Fill all plots for an event
  //void fillPlots( const edm::Event&, const SelectorDecisions& ) {mPhotonData->Fill()};
  //void fillTuple( const edm::Event& ) {mPhotonData->Fill(); }
  //void fillHealthPlots( const edm::Event& ) {mPhotonData->Fill(); }
  
  

private:
  
  //configuration parameters
  edm::ParameterSet photonParams;

  edm::InputTag photTag_;
  edm::InputTag ccphotTag_;
  edm::InputTag pfphotTag_;

  edm::InputTag genTag_;


  double photMaxEta_, photMaxEt_, photMinEt_, photRelIso_;  /// for prelection cuts on photons
  bool   doMCData_;                 /// switch to turn off generator level information
  int    debug_;
    
  // Plots
  TTree * mPhotonData;      /// Will contain the additional di-jet specific data

  // Variables
  int    m_Nphot;
  double m_PhotE[50];
  double m_PhotEt[50];
  double m_PhotPt[50];
  double m_PhotPx[50];
  double m_PhotPy[50];
  double m_PhotPz[50];
  double m_PhotEta[50];
  double m_PhotPhi[50];

  double m_PhotTrkIso[50];
  double m_PhotECalIso[50];
  double m_PhotHCalIso[50];
  double m_PhotAllIso[50];

  //bool m_ccPhotAssoc[50];
  bool m_PhotLooseEM[50];
  bool m_PhotLoosePhoton[50];
  bool m_PhotTightPhoton[50];

  double m_GenPhotPdgId[50];
  double m_GenPhotMother[50];
  double m_GenPhotPx[50];
  double m_GenPhotPy[50];
  double m_GenPhotPz[50];
  double m_GenPhotPt[50];
  double m_GenPhotEt[50];
  double m_GenPhotE[50];
  //

  int genPhotLength;
  int genPhotIds[100];
  int genPhotRefs[100];
  int genPhotStatus[100];
  float genPhotE[100];
  float genPhotPx[100];
  float genPhotPy[100];
  float genPhotPz[100];


  bool init_;                          // vectors initialised or not


  std::string outputFileName_;

  //input from .cfg

  double localPi;

};

#endif
