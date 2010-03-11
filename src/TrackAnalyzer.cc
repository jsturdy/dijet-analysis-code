
// -*- C++ -*-
//
// Package:    DiJetAnalysis
// Class:      TrackAnalyzer
// 
/**\class TrackAnalyzer TrackAnalyzer.cc JSturdy/DiJetAnalysis/src/TrackAnalyzer.cc

Description: Variable collector/ntupler for SUSY search with Jets + MET

Implementation:Uses the EventSelector interface for event selection and TFileService for plotting.

*/
//
// Original Author:  Markus Stoye, (modified by Jared Sturdy from SusyDiJetAnalysis)
//         Created:  Mon Feb 18 15:40:44 CET 2008
// $Id: TrackAnalyzer.cpp,v 1.3 2010/01/28 23:21:08 sturdy Exp $
//
//
#include "JSturdy/DiJetAnalysis/interface/TrackAnalyzer.h"
#include <TMath.h>
using namespace std;
using namespace reco;
using namespace edm;


//________________________________________________________________________________________
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig)
{ 

  doMCData_ = true;
  // Say something about event weights
  if (iConfig.exists("doMCData"))
    doMCData_  = iConfig.getParameter<bool>("doMCData");
  if (doMCData_)
    if (iConfig.exists("genTag"))
      genTag_  = iConfig.getParameter<edm::InputTag>("genTag");
 
  //edm::LogInfo("TrackTest") << "Global event weight set to " << eventWeight_;
    
  //jptTag_    = iConfig.getParameter<edm::InputTag>("jptTag");

  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
  initTuple();
}


//________________________________________________________________________________________
TrackAnalyzer::~TrackAnalyzer() {}


//________________________________________________________________________________________
// Method called to for each event
bool
TrackAnalyzer::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  using namespace edm;

  track_result = false;
  //bool preselection = false;
  edm::LogVerbatim("TrackEvent") << " Start  " << std::endl;

  std::ostringstream dbg;

  
  //get tracks
  
  //    reco::TrackCollection myTracks;
  // iEvent.getByLabel("generalTracks",myTracks);
  edm::Handle<View<reco::Track> >  myTracks;
  iEvent.getByLabel("generalTracks",myTracks);
  double ptMax_ = 500;
  math::XYZTLorentzVector totalP3;
  for(View<reco::Track>::const_iterator elem = myTracks->begin(); 
      elem != myTracks->end(); ++elem) {
    
    if (!(elem->quality(reco::TrackBase::highPurity))) continue;
    
    double elemPt = elem->pt();
    
    if ( elemPt > ptMax_) continue;
    
    math::XYZTLorentzVector p3(elem->px(),elem->py(),elem->pz(),elem->p());
    totalP3 -= p3;
    
  }
  
  m_MPTPhi= totalP3.phi();
  m_MPTPx = totalP3.px();
  m_MPTPy = totalP3.py();
  m_MPTPz = totalP3.pz();
   

  // Fill the tree only if all preselection conditions are met
  //if(preselection) //mPreselection->Fill();
  return track_result;
  //mTrackData->Fill();
}

//________________________________________________________________________________________
void 
TrackAnalyzer::beginJob(const edm::EventSetup&) {}

//________________________________________________________________________________________
void 
TrackAnalyzer::endJob() {

}


//________________________________________________________________________________________
void
TrackAnalyzer::initTuple() {

  std::ostringstream variables; // Container for all variables
  
  // 1. Event variables
  variables << "weight:process";
  
  // Register this ntuple
  edm::Service<TFileService> fs;

  // Now we add some additional ones for the dijet analysis
  mTrackData = fs->make<TTree>( "TrackData", "data after cuts" );
  mTrackData->SetAutoSave(10);

  //jet correction factors
  //mTrackData->Branch("Jet_MCcorrFactor",m_JetMCCorrFactor,"m_JetMCCorrFactor[Njets]/double");
  //mTrackData->Branch("Jet_JPTcorrFactor",m_JetJPTCorrFactor,"m_JetJPTCorrFactor[Njets]/double");
  
  mTrackData->Branch("JetTrackPt",m_JetTrackPt,"JetTrackPt[Njets]/double"); 
  mTrackData->Branch("JetTrackPhi",m_JetTrackPhi,"JetTrackPhi[Njets]/double"); 
  mTrackData->Branch("JetTrackPhiWeighted",m_JetTrackPhiWeighted,"JetTrackPhiWeighted[Njets]/double"); 
  mTrackData->Branch("JetTrackNo",m_JetTrackNo,"JetTrackNo[Njets]/int");

  // MPT Markus 
  mTrackData->Branch("MPTPhi" ,& m_MPTPhi ,"MPTPhi/double");
  mTrackData->Branch("MPTPx" ,& m_MPTPx ,"MPTPx/double");
  mTrackData->Branch("MPTPy" ,& m_MPTPy ,"MPTPy/double");
  mTrackData->Branch("MPTPz" ,& m_MPTPz ,"MPTPz/double");
    
  edm::LogInfo("TrackEvent") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(TrackAnalyzer);
