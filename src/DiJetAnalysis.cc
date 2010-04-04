
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
// $Id: DiJetAnalysis.cc,v 1.2 2010/03/29 11:19:36 sturdy Exp $
//
//
#include "JSturdy/DiJetAnalysis/interface/DiJetAnalysis.h"

#include <TMath.h>
#include <sstream>

//________________________________________________________________________________________
DiJetAnalysis::DiJetAnalysis(const edm::ParameterSet& iConfig)
{ 

  //default parameters
  if (iConfig.exists("debugDiJets"))     debug_   = iConfig.getUntrackedParameter<int>("debugDiJets");
 
  //  edm::LogInfo("SusyDiJetAnalysis") << "Global event weight set to " << eventWeight_;
  edm::Service<TFileService> fs;

  mAllData = fs->make<TTree>( "AllData", "data after preselection" );
  mAllData->SetAutoSave(10);
    
  jetinfo  = new JetAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("jetParameters"), mAllData);
  metinfo  = new METAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("metParameters"), mAllData);
  photons  = new PhotonAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("photonParameters"), mAllData);
  leptons  = new LeptonAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("leptonParameters"), mAllData);
  vertex   = new VertexAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("vertexParameters"), mAllData);
  tracks   = new TrackAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("trackParameters"), mAllData);
  triggers = new TriggerAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("triggerParameters"), mAllData);
  //heminfo  = new HemisphereAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("hemisphereParameters"), mAllData));


  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
  //initPlots();
}


//________________________________________________________________________________________
DiJetAnalysis::~DiJetAnalysis() {}


//________________________________________________________________________________________
// Method called to for each event
void
DiJetAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;
  using namespace edm;

  bool preselection = false;
  edm::LogVerbatim("DiJetAnalysis") << " Start  " << std::endl;

  std::ostringstream dbg;

  bool myjetresult     = jetinfo->filter(iEvent, iSetup);
  bool mymetresult     = metinfo->filter(iEvent, iSetup);
  bool myphotonresult  = photons->filter(iEvent, iSetup);
  bool myleptonresult  = leptons->filter(iEvent, iSetup);
  bool myvertexresult  = vertex->filter(iEvent, iSetup);
  bool mytrackresult   = tracks->filter(iEvent, iSetup);
  bool mytriggerresult = triggers->filter(iEvent, iSetup);
  //bool myhemresult     = heminfo->filter(iEvent, iSetup);

  preselection = myjetresult;// && myleptonresult;// && mymetresult && myvertexresult && mytriggerresult;
  
  //if (preselection) {
  mAllData->Fill();
  //}
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

}


//________________________________________________________________________________________
void
DiJetAnalysis::initPlots() {

  std::ostringstream variables; // Container for all variables
  
  // 1. Event variables
  variables << "weight:process";
  
  // Register this ntuple

    
  edm::LogInfo("DiJetAnalysis") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________

// Define this as a plug-in

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DiJetAnalysis);
