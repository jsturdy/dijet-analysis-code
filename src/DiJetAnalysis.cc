
// -*- C++ -*-
//
// Package:    DiJetAnalysis
// Class:      DiJetAnalysis
// 
/**\class DiJetAnalysis DiJetAnalysis.cc JSturdy/DiJetAnalysis/src/DiJetAnalysis.cc

Description: Variable collector/ntupler for SUSY search with Jets + MET

Implementation:Uses the EventSelector interface for event selection and TFileService for plotting.

*/
//
// Original Author:  Markus Stoye, (modified by Jared Sturdy from SusyDiJetAnalysis)
//         Created:  Mon Feb 18 15:40:44 CET 2008
// $Id: DiJetAnalysis.cpp,v 1.3 2010/01/28 23:21:08 sturdy Exp $
//
//
#include "JSturdy/DiJetAnalysis/interface/DiJetAnalysis.h"

#include <TMath.h>
#include <sstream>

using namespace std;
//using namespace reco;
using namespace edm;

//________________________________________________________________________________________
DiJetAnalysis::DiJetAnalysis(const edm::ParameterSet& iConfig)
{ 

  //default parameters
 
  //  edm::LogInfo("SusyDiJetAnalysis") << "Global event weight set to " << eventWeight_;
    
  jetinfo  = new JetAnalyzer(iConfig.getParameter<ParameterSet>("jetAnalyzer"));
  //heminfo  = new HemisphereAnalyzer(iConfig.getParameter<ParameterSet>("hemisphereAnalyzer"));
  leptons  = new LeptonAnalyzer(iConfig.getParameter<ParameterSet>("leptonAnalyzer"));
  metinfo  = new METAnalyzer(iConfig.getParameter<ParameterSet>("metAnalyzer"));
  vertex   = new VertexAnalyzer(iConfig.getParameter<ParameterSet>("vertexAnalyzer"));
  //tracks   = new TrackAnalyzer(iConfig.getParameter<ParameterSet>("trackAnalyzer"));
  //triggers = new TriggerAnalyzer(iConfig.getParameter<ParameterSet>("triggerAnalyzer"));


  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
  initPlots();
}


//________________________________________________________________________________________
DiJetAnalysis::~DiJetAnalysis() {}


//________________________________________________________________________________________
// Method called to for each event
void
DiJetAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  //bool preselection = false;
  edm::LogVerbatim("DiJetAnalysis") << " Start  " << std::endl;

  std::ostringstream dbg;

  mAllData->Fill();
}

//________________________________________________________________________________________
void 
DiJetAnalysis::beginJob(const edm::EventSetup&) {
  //setup job before events
}

//________________________________________________________________________________________
void 
DiJetAnalysis::endJob() {
  //cleanup after all events
}

void
DiJetAnalysis::printSummary( void ) {

  // prints an HLT report -- associates trigger bits with trigger names (prints #events fired the trigger etc)
//  const unsigned int n(pathNames_.size());
//  std::cout << "\n";
//  std::cout << "HLT-Report " << "---------- Event  Summary ------------\n";
//  std::cout << "HLT-Report"
//	    << " Events total = " << nEvents_
//	    << " wasrun = " << nWasRun_
//	    << " passed = " << nAccept_
//	    << " errors = " << nErrors_
//	    << "\n";
//
//  std::cout << endl;
//  std::cout << "HLT-Report " << "---------- HLTrig Summary ------------\n";
//  std::cout << "HLT-Report "
//	    << right << setw(10) << "HLT  Bit#" << " "
//	    << right << setw(10) << "WasRun" << " "
//	    << right << setw(10) << "Passed" << " "
//	    << right << setw(10) << "Errors" << " "
//	    << "Name" << "\n";
//
//  if (init_) {
//    for (unsigned int i=0; i!=n; ++i) {
//      std::cout << "HLT-Report "
//		<< right << setw(10) << i << " "
//		<< right << setw(10) << hlWasRun_[i] << " "
//		<< right << setw(10) << hlAccept_[i] << " "
//		<< right << setw(10) << hlErrors_[i] << " "
//		<< pathNames_[i] << "\n";
//    }
//  } else {
//    std::cout << "HLT-Report - No HL TriggerResults found!" << endl;
//  }
//  
//  std::cout << endl;
//  std::cout << "HLT-Report end!" << endl;
//  std::cout << endl;

}


//________________________________________________________________________________________
void
DiJetAnalysis::initPlots() {

  std::ostringstream variables; // Container for all variables
  
  // 1. Event variables
  variables << "weight:process";
  
  // Register this ntuple
  edm::Service<TFileService> fs;

  // Now we add some additional ones for the dijet analysis
  //mPreselection = fs->make<TTree>( "Preselection", "data after preselection" );
  mAllData = fs->make<TTree>( "AllData", "data after cuts" );
  mAllData->SetAutoSave(10);

    
  edm::LogInfo("DiJetAnalysis") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________

// Define this as a plug-in

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DiJetAnalysis);
