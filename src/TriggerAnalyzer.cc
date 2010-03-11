
// -*- C++ -*-
//
// Package:    SusyDiJetAnalysis
// Class:      TriggerAnalyzer
// 
/**\class TriggerAnalyzer TriggerAnalyzer.cc JSturdy/DiJetAnalysis/src/TriggerAnalyzer.cc

Description: Variable collector/ntupler for SUSY search with Jets + MET

Implementation:Uses the EventSelector interface for event selection and TFileService for plotting.

*/
//
// Original Author:  Markus Stoye, (modified by Jared Sturdy from SusyDiJetAnalysis)
//         Created:  Mon Feb 18 15:40:44 CET 2008
// $Id: TriggerAnalyzer.cpp,v 1.3 2010/01/28 23:21:08 sturdy Exp $
//
//
//#include "SusyAnalysis/EventSelector/interface/BJetEventSelector.h"
#include "JSturdy/DiJetAnalysis/interface/TriggerAnalyzer.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#define _USE_MATH_DEFINES
//#include <math.h>
#include <TMath.h>
#include <sstream>
using namespace std;
using namespace reco;
using namespace edm;

//________________________________________________________________________________________
TriggerAnalyzer::TriggerAnalyzer(const edm::ParameterSet& iConfig)
{ 

  if (iConfig.exists("doMCData"))  doMCData_  = iConfig.getParameter<bool>("doMCData");
 
  // trigger stuff
  triggerResults_ = iConfig.getParameter<edm::InputTag>("triggerResults");
  // trigger path names
  pathNames_ = iConfig.getParameter< std::vector<std::string> >("pathNames");

  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
  initTuple();
}


//________________________________________________________________________________________
TriggerAnalyzer::~TriggerAnalyzer() {}


//________________________________________________________________________________________
// Method called to for each event
bool
TriggerAnalyzer::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  //bool preselection = false;
  edm::LogVerbatim("TriggerEvent") << " Start  " << std::endl;

  std::ostringstream dbg;

  trigger_result = false;

  //Trigger results
  edm::LogVerbatim("TriggerEvent") << " Trigger decision  " << std::endl;

  //get the trigger decision
  m_HLT1JET=false;
  m_HLT2JET=false;
  m_HLT1MET1HT=false;

  // Get the trigger results and check validity
  edm::Handle<edm::TriggerResults> hltHandle;
  iEvent.getByLabel(triggerResults_, hltHandle);
  if ( !hltHandle.isValid() ) {
    edm::LogWarning("HLTEventSelector") << "No trigger results for InputTag " << triggerResults_;
    return false;
  }

  //  std::cout << " get results " << std::endl;

  // Get results
  edm::TriggerNames trgNames;
  trgNames.init(*hltHandle);
  unsigned int trgSize = trgNames.size();
  
  // Example for OR of all specified triggers

  edm::LogWarning("HLTEventSelector") << " triggers " << trgNames.size() << std::endl;

  // GEORGIA
  if (!hltHandle.isValid()) {
    // triggerExists = false;
    std::cout << "HLTriggerResult Not Valid!" << endl;
  }
  else {  
    if (hltHandle->wasrun()) nWasRun_++;
    const bool accept(hltHandle->accept());
    LogDebug("") << "HL TriggerResults decision: " << accept;
    if (accept) ++nAccept_;
    if (hltHandle->error() ) nErrors_++;
  }
  if (!init_) {
    init_=true;
    triggerNames_.init(*hltHandle);
    pathNames_=triggerNames_.triggerNames();
    const unsigned int n(pathNames_.size());
    hlWasRun_.resize(n);
    hlAccept_.resize(n);
    hlErrors_.resize(n);
    for (unsigned int i=0; i!=n; ++i) {
      hlWasRun_[i]=0;
      hlAccept_[i]=0;
      hlErrors_[i]=0;
    }}

  // decision for each HL algorithm
  const unsigned int n(pathNames_.size());
  for (unsigned int i=0; i!=n; ++i) {
    if (hltHandle->wasrun(i)) hlWasRun_[i]++;
    if (hltHandle->accept(i)) hlAccept_[i]++;
    if (hltHandle->error(i) ) hlErrors_[i]++;
  }
  
  m_nHLT=static_cast<int>(n);
  for(unsigned int i=0; i!=n; ++i) {
    m_HLTArray[i]=hltHandle->accept(i);
  }

  //looping over list of trig path names
  for ( std::vector<std::string>::const_iterator i=pathNames_.begin();
	i!=pathNames_.end(); ++i ) {
    // Get index
 
    unsigned int index = trgNames.triggerIndex(*i);
    if ( index==trgSize ) {
      edm::LogWarning("HLTEventSelector") << "Unknown trigger name " << *i;
      continue;
    }
    if ( hltHandle->accept(index) ) {
      LogDebug("HLTEventSelector") << "Event selected by " << *i ;
      std::string trigName = *i;
      if (trigName == "HLT_Jet180")      m_HLT1JET=true;
      if (trigName == "HLT_DiJetAve130") m_HLT2JET=true;
      if (trigName == "HLT_MET50")       m_HLT1MET1HT=true;
      if (trigName == "HLT_Mu9")         m_HLT1Muon=true; 
      
    } 
  }

  // Fill the tree only if all preselection conditions are met
  //if(preselection) //mPreselection->Fill();
  return trigger_result;
  //mTriggerData->Fill();
}


//________________________________________________________________________________________
void 
TriggerAnalyzer::beginJob(const edm::EventSetup&) {}

//________________________________________________________________________________________
void 
TriggerAnalyzer::endJob() {

}

void
TriggerAnalyzer::printHLTreport( void ) {

  // prints an HLT report -- associates trigger bits with trigger names (prints #events fired the trigger etc)
  const unsigned int n(pathNames_.size());
  std::cout << "\n";
  std::cout << "HLT-Report " << "---------- Event  Summary ------------\n";
  std::cout << "HLT-Report"
	    << " Events total = " << nEvents_
	    << " wasrun = " << nWasRun_
	    << " passed = " << nAccept_
	    << " errors = " << nErrors_
	    << "\n";

  std::cout << endl;
  std::cout << "HLT-Report " << "---------- HLTrig Summary ------------\n";
  std::cout << "HLT-Report "
	    << right << setw(10) << "HLT  Bit#" << " "
	    << right << setw(10) << "WasRun" << " "
	    << right << setw(10) << "Passed" << " "
	    << right << setw(10) << "Errors" << " "
	    << "Name" << "\n";

  if (init_) {
    for (unsigned int i=0; i!=n; ++i) {
      std::cout << "HLT-Report "
		<< right << setw(10) << i << " "
		<< right << setw(10) << hlWasRun_[i] << " "
		<< right << setw(10) << hlAccept_[i] << " "
		<< right << setw(10) << hlErrors_[i] << " "
		<< pathNames_[i] << "\n";
    }
  } else {
    std::cout << "HLT-Report - No HL TriggerResults found!" << endl;
  }
  
  std::cout << endl;
  std::cout << "HLT-Report end!" << endl;
  std::cout << endl;

}
// end GEORGIA


//________________________________________________________________________________________
void
TriggerAnalyzer::initTuple() {

  std::ostringstream variables; // Container for all variables
  
  // 1. Event variables
  variables << "weight:process";
  
  // Register this ntuple
  edm::Service<TFileService> fs;

  // Now we add some additional ones for the dijet analysis
  //mPreselection = fs->make<TTree>( "Preselection", "data after preselection" );
  mTriggerData = fs->make<TTree>( "TriggerData", "data after cuts" );
  mTriggerData->SetAutoSave(10);

  // Add the branches
  //mTriggerData->Branch("run",   &Run,   "run/int");
  //mTriggerData->Branch("event", &Event, "event/int");

  mTriggerData->Branch("nHLT",    &m_nHLT,     "nHLT/I");
  mTriggerData->Branch("HLTArray", m_HLTArray, "HLTArray[nHLT]/I");

  //Trigger information
  mTriggerData->Branch("HLT1JET",    &m_HLT1JET,    "HLT1JET/bool");
  mTriggerData->Branch("HLT2JET",    &m_HLT2JET,    "HLT2JET/bool");
  mTriggerData->Branch("HLT1MET1HT", &m_HLT1MET1HT, "HLT1MET1HT/bool");
  mTriggerData->Branch("HLT1MUON",   &m_HLT1Muon,   "HLT1MUON/bool");

  edm::LogInfo("TriggerEvent") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(TriggerAnalyzer);
