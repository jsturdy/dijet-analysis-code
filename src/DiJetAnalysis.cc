
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
// $Id: DiJetAnalysis.cc,v 1.3 2010/04/04 00:04:01 sturdy Exp $
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
 

  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
  initPlots();
    
  jetinfo  = new JetAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("jetParameters"), mAllData);
  metinfo  = new METAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("metParameters"), mAllData);
  photons  = new PhotonAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("photonParameters"), mAllData);
  leptons  = new LeptonAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("leptonParameters"), mAllData);
  vertex   = new VertexAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("vertexParameters"), mAllData);
  tracks   = new TrackAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("trackParameters"), mAllData);
  triggers = new TriggerAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("triggerParameters"), mAllData);
  //heminfo  = new HemisphereAnalyzer(iConfig.getUntrackedParameter<edm::ParameterSet>("hemisphereParameters"), mAllData));

  //Setup counters for filters
  passJets[0]        = 0;
  passMET[0]         = 0;
  passLeptons[0]     = 0;
  passPhotons[0]     = 0;
  passVertex[0]      = 0;
  passTracks[0]      = 0;
  passTriggers[0]    = 0;
  //passHemispheres[0] = 0;
  passJets[1]        = 0;
  passMET[1]         = 0;
  passLeptons[1]     = 0;
  passPhotons[1]     = 0;
  passVertex[1]      = 0;
  passTracks[1]      = 0;
  passTriggers[1]    = 0;
  //passHemispheres[1] = 0;

  nEvents_ = 0;

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

  m_Run   = iEvent.id().run();
  m_Event = iEvent.id().event();

  //Run filters

  bool myjetresult     = jetinfo->filter(iEvent, iSetup);
  bool mymetresult     = metinfo->filter(iEvent, iSetup);
  bool myphotonresult  = photons->filter(iEvent, iSetup);
  bool myleptonresult  = leptons->filter(iEvent, iSetup);
  bool myvertexresult  = vertex->filter(iEvent, iSetup);
  bool mytrackresult   = tracks->filter(iEvent, iSetup);
  bool mytriggerresult = triggers->filter(iEvent, iSetup);
  //bool myhemresult     = heminfo->filter(iEvent, iSetup);

  if (myjetresult)     ++passJets[0];
  if (mymetresult)     ++passMET[0];
  if (myleptonresult)  ++passLeptons[0];
  if (myphotonresult)  ++passPhotons[0];
  if (myvertexresult)  ++passVertex[0];
  if (mytrackresult)   ++passTracks[0];
  if (mytriggerresult) ++passTriggers[0];
  //if (myhemresult)   ++passHemispheres[0];

  if ( mymetresult && myleptonresult && myphotonresult
    && myvertexresult && mytrackresult && mytriggerresult ) ++passJets[1];
  if ( myjetresult && myleptonresult && myphotonresult
    && myvertexresult && mytrackresult && mytriggerresult ) ++passMET[1];
  if ( myjetresult && mymetresult && myphotonresult
    && myvertexresult && mytrackresult && mytriggerresult ) ++passLeptons[1];
  if ( myjetresult && mymetresult && myleptonresult
    && myvertexresult && mytrackresult && mytriggerresult ) ++passPhotons[1];
  if ( myjetresult && mymetresult && myleptonresult && myphotonresult
     && mytrackresult && mytriggerresult )                  ++passVertex[1];
  if ( myjetresult && mymetresult && myleptonresult && myphotonresult
    && myvertexresult && mytriggerresult )                  ++passTracks[1];
  if ( myjetresult && mymetresult && myleptonresult && myphotonresult
    && myvertexresult && mytrackresult )                    ++passTriggers[1];
  //if ( myjetresult && mymetresult && myleptonresult && myphotonresult)\
  //&&(myvertexresult && mytrackresult && mytriggerresult ) ++passHemispheres[1];
  
  ++nEvents_;

  preselection = myjetresult;// && myleptonresult;// && mymetresult && myvertexresult && mytriggerresult;
  
  //if (preselection) {
  mAllData->Fill();
  //}
}

//________________________________________________________________________________________
void 
DiJetAnalysis::beginJob() {
  //setup job before events
}

//________________________________________________________________________________________
void 
DiJetAnalysis::endJob() {
  //cleanup after all events
  printSummary();
}

void
DiJetAnalysis::printSummary( void ) {
  printf("============================Summary of DiJetAnalyzer============================\n");
  printf("= Total number of events filtered:     %2d                                     =\n",nEvents_);
  printf("= Events passing the jet filter:       %2d                                     =\n",passJets[0]);
  printf("= Events passing the MET filter:       %2d                                     =\n",passMET[0]);
  printf("= Events passing the lepton veto:      %2d                                     =\n",passLeptons[0]);
  printf("= Events passing the photon veto:      %2d                                     =\n",passPhotons[0]);
  printf("= Events passing the vertex filter:    %2d                                     =\n",passVertex[0]);
  printf("= Events passing the track filter:     %2d                                     =\n",passTracks[0]);
  printf("= Events passing the trigger filter:   %2d                                     =\n",passTriggers[0]);
  //printf("= Events passing the hemisphere filter:%2d                                     =\n",passHemispheres[0]);
  printf("=                                                                              =\n");
  printf("=                                                                              =\n");
  printf("=                                                                              =\n");
  printf("=                                                                              =\n");
  printf("=                                                                              =\n");
  printf("= Events passing all but the jet filter:       %2d                             =\n",passJets[1]);
  printf("= Events passing all but the MET filter:       %2d                             =\n",passMET[1]);
  printf("= Events passing all but the Lepton veto:      %2d                             =\n",passLeptons[1]);
  printf("= Events passing all but the photon veto:      %2d                             =\n",passPhotons[1]);
  printf("= Events passing all but the vertex filter:    %2d                             =\n",passVertex[1]);
  printf("= Events passing all but the track filter:     %2d                             =\n",passTracks[1]);
  printf("= Events passing all but the trigger filter:   %2d                             =\n",passTriggers[1]);
  //printf("= Events passing all but the hemisphere filter:%2d                             =\n",passHemispheres[1]);
  printf("================================================================================\n");
}


//________________________________________________________________________________________
void
DiJetAnalysis::initPlots() {

  std::ostringstream variables; // Container for all variables
  
  // 1. Event variables
  variables << "weight:process";
  
  // Register this ntuple
  edm::Service<TFileService> fs;

  mAllData = fs->make<TTree>( "AllData", "data after preselection" );
  mAllData->SetAutoSave(10);

  mAllData->Branch("Run",    &m_Run,     "Run/int");
  mAllData->Branch("Event",  &m_Event,   "Event/int");


    
  edm::LogInfo("DiJetAnalysis") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________

// Define this as a plug-in

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DiJetAnalysis);
