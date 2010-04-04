
// -*- C++ -*-
//
// Package:    DiJetAnalysis
// Class:      TrackAnalyzer
// 
/**\class TrackAnalyzer TrackAnalyzer.cc JSturdy/DiJetAnalysis/src/TrackAnalyzer.cc

Description: Collects variables related to tracks


*/
//
// Original Author:  Jared Sturdy (from SusyDiJetAnalysis)
//         Created:  Fri Jan 29 16:10:31 PDT 2010
// $Id: TrackAnalyzer.cc,v 1.2 2010/03/29 11:18:22 sturdy Exp $
//
//
#include "JSturdy/DiJetAnalysis/interface/TrackAnalyzer.h"
#include <TMath.h>

//________________________________________________________________________________________
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& pset, TTree* tmpAllData)
{ 
  mTrackData = tmpAllData;
  trackParams = pset;
  doMCData_ = true;
  debug_    =0;

  if (trackParams.exists("debugTracks"))     debug_   = trackParams.getUntrackedParameter<int>("debugTracks");
  if (trackParams.exists("doMCTracks"))
    doMCData_  = trackParams.getUntrackedParameter<bool>("doMCTracks");
  if (doMCData_)
    if (trackParams.exists("trackTag"))
      trackTag_  = trackParams.getUntrackedParameter<edm::InputTag>("trackTag");
 
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
  using namespace reco;
  using namespace edm;

  track_result = true;
  edm::LogVerbatim("TrackEvent") << " Start  " << std::endl;

  std::ostringstream dbg;

  
  //get tracks
  
  edm::Handle<View<reco::Track> >  myTracks;
  iEvent.getByLabel(trackTag_,myTracks);
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
   

  return track_result;
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
  

  mTrackData->Branch("MPTPhi", &m_MPTPhi, "MPTPhi/double");
  mTrackData->Branch("MPTPx",  &m_MPTPx,  "MPTPx/double");
  mTrackData->Branch("MPTPy",  &m_MPTPy,  "MPTPy/double");
  mTrackData->Branch("MPTPz",  &m_MPTPz,  "MPTPz/double");

  edm::LogInfo("TrackEvent") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(TrackAnalyzer);
