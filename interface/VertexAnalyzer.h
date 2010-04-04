#ifndef VERTEXANALYZER
#define VERTEXANALYZER

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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


//
// Class declaration
//
class VertexAnalyzer {
 public:
  VertexAnalyzer(const edm::ParameterSet&, TTree*);
  ~VertexAnalyzer();
  
  //*** CMSSW interface
  /// Called once per job, at start
  void beginJob(const edm::EventSetup&) ;
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
  //void fillTuple(){mVertexData->Fill(); }
  //void fillTuple( ) {mAllData->Fill(); }
  //void fillHealthPlots(){mVertexData->Fill(); }



 private:
  
  //configuration parameters
  edm::ParameterSet vertexParams;
  edm::InputTag vtxTag_;
  double _minNVtx, _minVtxTrks, _minVtxNdof, _maxVtxChi2, _maxVtxZ;   /// for primary vertex selection
  int    debug_;
  //
  
  // Counters
  unsigned int nrEventTotalRaw_;          ///< Raw number of events (+1 at each event)
    
  // Plots
  TTree * mVertexData;    //Will contain the data passing the vertex selection
  //TTree * mAllData;       //Will contain the preselection data

  bool   vertexDecision;
  
  int    m_nVtx;
  int    m_VtxNTrks[10];
  double m_VtxChi2[5];
  double m_VtxNdof[5];
  double m_VtxIsValid[5];
  double m_VtxSumTrkPt[10];
  double m_VtxNormalizedChi2[5];
  double m_VtxX[5];
  double m_VtxY[5];
  double m_VtxZ[5];
  double m_VtxdX[5];
  double m_VtxdY[5];
  double m_VtxdZ[5];

  bool init_;                          // vectors initialised or not

};

#endif















