#ifndef DIJETANALYSIS
#define DIJETANALYSIS

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
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "JSturdy/DiJetAnalysis/interface/JetAnalyzer.h"
#include "JSturdy/DiJetAnalysis/interface/METAnalyzer.h"
//#include "JSturdy/DiJetAnalysis/interface/HemisphereAnalyzer.h"
#include "JSturdy/DiJetAnalysis/interface/LeptonAnalyzer.h"
#include "JSturdy/DiJetAnalysis/interface/VertexAnalyzer.h"
#include "JSturdy/DiJetAnalysis/interface/PhotonAnalyzer.h"
#include "JSturdy/DiJetAnalysis/interface/TrackAnalyzer.h"
#include "JSturdy/DiJetAnalysis/interface/TriggerAnalyzer.h"


//
// Class declaration
//
class DiJetAnalysis : public edm::EDAnalyzer {
public:
  explicit DiJetAnalysis(const edm::ParameterSet&);
  ~DiJetAnalysis();
  
private:
  //*** CMSSW interface
  /// Called once per job, at start
  virtual void beginJob(const edm::EventSetup&) ;
  /// Called for each event
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  /// Called once per job, at end
  virtual void endJob();

  /// Print a summary of counts for all selectors
  virtual void printSummary(void);
 
  //*** Plotting
  //virtual void fillPlots( const edm::Event&, const SelectorDecisions& );
  void initPlots();

  //  virtual bool filter(const edm::Event& evt,const edm::EventSetup& iSetup );


private:

  // Event information
  int run_, event_;

  // Counters
  unsigned int nrEventTotalRaw_;          ///< Raw number of events (+1 at each event)
  double nrEventTotalWeighted_;           ///< Weighted #(events)
  std::vector<float> nrEventSelected_;    ///< Selected #(events) for each module
  std::vector<float> nrEventAllButOne_;   ///< All-but-one selected #(events) for each module
  std::vector<float> nrEventCumulative_;  ///< Cumulative selected #(events) for each module
    
  // Plots
  TNtuple* ntuple_;      /// Will contain all the selector information we want to keep
  TTree * mAllData;      /// Will contain the additional di-jet specific data

  float* variables_;     ///< Container for the tree variables (from selectors)
  bool*  decisions_;     ///< Container for all selector decisions
  bool   globalDec_;     ///< Global decision for event

  int m_Run;
  int m_Event;

  int    mGlobalDecision;

  unsigned int  nEvents_;              // number of events processed

  bool init_;                          // vectors initialised or not

  int debug_;
  
  std::string outputFileName_;

  double localPi;
  unsigned int *mSelectorResults;

  JetAnalyzer        * jetinfo;
  METAnalyzer        * metinfo;
  //HemisphereAnalyzer * heminfo;
  PhotonAnalyzer     * photons;
  LeptonAnalyzer     * leptons;
  VertexAnalyzer     * vertex;
  TrackAnalyzer      * tracks;
  TriggerAnalyzer    * triggers;

};
#endif
