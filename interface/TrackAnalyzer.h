#ifndef TRACKANALYZER
#define TRACKANALYZER

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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//
// Class declaration
//
class TrackAnalyzer {
 public:
  TrackAnalyzer(const edm::ParameterSet&, TTree*);
  ~TrackAnalyzer();
  
  //*** CMSSW interface
  /// Called once per job, at start
  void beginJob(const edm::EventSetup&);
  /// Called for each event
  //virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool filter(const edm::Event& evt,const edm::EventSetup& iSetup );
  /// Called once per job, at end
  void endJob();

  //*** Plotting
  /// Define all plots
  void initTuple();
  //void bookHealthPlots();
  /// Fill all plots for an event
  //void fillTuple( ) {mTrackData->Fill(); }
  //void fillTuple( ) {mAllData->Fill(); }
  //void fillHealthPlots( ) {mTrackData->Fill(); }



 private:

  // Data tags
  edm::ParameterSet trackParams;
  edm::InputTag  trackTag_;
  bool doMCData_;
  int  debug_;
  //TTree * mAllData;        /// Will contain all the data
  TTree * mTrackData;      /// Will contain the additional track parameters

  bool track_result;

  // track info:
  int    m_JetTrackNo[50];
  double m_JetTrackPhi[50];
  double m_JetTrackPhiWeighted[50];
  double m_JetTrackPt[50];

  double m_MPTPhi;
  double m_MPTPx;
  double m_MPTPy;
  double m_MPTPz;

  double localPi;
};

#endif
