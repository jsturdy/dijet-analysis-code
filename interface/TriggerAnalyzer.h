#ifndef TRIGGERANALYZER
#define TRIGGERANALYZER

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

//#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Framework/interface/TriggerNames.h"

//
// Class declaration
//
class TriggerAnalyzer {
 public:
  TriggerAnalyzer(const edm::ParameterSet&);
  ~TriggerAnalyzer();
  
 private:
  //*** CMSSW interface
  /// Called once per job, at start
  void beginJob(const edm::EventSetup&) ;
  /// Called for each event
  //virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool filter(const edm::Event& evt,const edm::EventSetup& iSetup );

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
  void fillTuple( const edm::Event& ) {mTriggerData->Fill(); }
  //void fillHealthPlots( const edm::Event& ) {mTriggerData->Fill(); }



private:

  bool doMCData_;                 /// switch to turn off generator level information

  // Plots
  TNtuple* ntuple_;      /// Will contain all the selector information we want to keep
  TTree * mAllData;      /// Will contain the additional di-jet specific data
  TTree * mTriggerData;      /// Will contain the additional di-jet specific data

  int m_nHLT;
  int m_HLTArray[200];

  bool m_HLT1JET;
  bool m_HLT2JET;
  bool m_HLT1MET1HT;
  bool m_HLT1Muon;

  bool trigger_result;

  // Data tags
  edm::InputTag triggerResults_; 
  std::vector<std::string> pathNames_;

  edm::TriggerNames triggerNames_;     // TriggerNames class

  unsigned int  nEvents_;              // number of events processed

  unsigned int  nWasRun_;              // # where at least one HLT was run
  unsigned int  nAccept_;              // # of accepted events
  unsigned int  nErrors_;              // # where at least one HLT had error
  std::vector<unsigned int> hlWasRun_; // # where HLT[i] was run
  std::vector<unsigned int> hlAccept_; // # of events accepted by HLT[i]
  std::vector<unsigned int> hlErrors_; // # of events with error in HLT[i]
  bool init_;                          // vectors initialised or not

  double localPi;
  unsigned int *mSelectorResults;

};

#endif
