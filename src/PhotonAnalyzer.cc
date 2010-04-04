
// -*- C++ -*-
//
// Package:    DiJetAnalysis
// Class:      PhotonAnalyzer
// 
/**\class PhotonAnalyzer PhotonAnalyzer.cc JSturdy/DiJetAnalysis/src/PhotonAnalyzer.cc

Description: Variable collector/ntupler for SUSY search with Jets + MET

*/
//
// Original Author:  Jared Sturdy
//         Created:  Fri Jan 29 16:10:31 PDT 2010
// $Id: PhotonAnalyzer.cc,v 1.1 2010/04/04 00:04:01 sturdy Exp $
//
//


//#include "SusyAnalysis/EventSelector/interface/BJetEventSelector.h"
#include "JSturdy/DiJetAnalysis/interface/PhotonAnalyzer.h"
//#define _USE_MATH_DEFINES
//#include <math.h>
#include <TMath.h>
#include <sstream>

//________________________________________________________________________________________
PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& pset, TTree* tmpAllData)
{ 

  mPhotonData = tmpAllData;
  photonParams = pset;

  //defaults
  photMaxEta_ = 2.5;
  photMaxEt_  = 50.;
  photMinEt_  = 5.;
  photRelIso_ = 0.5;

  doMCData_ = true;
  debug_ = 0;
  // Read in parameters from the config file
  if (photonParams.exists("photMaxEta")) photMaxEta_ = photonParams.getUntrackedParameter<double>("photMaxEta");
  if (photonParams.exists("photMaxEt"))  photMaxEt_  = photonParams.getUntrackedParameter<double>("photMaxEt");
  if (photonParams.exists("photMinEt"))  photMinEt_  = photonParams.getUntrackedParameter<double>("photMinEt");
  if (photonParams.exists("photRelIso")) photRelIso_ = photonParams.getUntrackedParameter<double>("photRelIso");

  if (photonParams.exists("doMCPhots"))   doMCData_   = photonParams.getUntrackedParameter<bool>("doMCPhots");
  if (doMCData_) 
    if (photonParams.exists("genPhotTag")) genTag_  = photonParams.getUntrackedParameter<edm::InputTag>("genPhotTag");
  if (photonParams.exists("debugPhots"))   debug_   = photonParams.getUntrackedParameter<int>("debugPhots");
 
  // get the data tags
  photTag_   = photonParams.getUntrackedParameter<edm::InputTag>("photTag");
  pfphotTag_ = photonParams.getUntrackedParameter<edm::InputTag>("pfphotTag");

  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
  initTuple();
}


//________________________________________________________________________________________
PhotonAnalyzer::~PhotonAnalyzer() {}


//________________________________________________________________________________________
// Method called to for each event
//void
//PhotonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
bool
PhotonAnalyzer::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  using namespace reco;
  using namespace edm;

  //bool preselection = false;
  bool photon_result = true;
  edm::LogVerbatim("PhotonEvent") << " Start  " << std::endl;

  std::ostringstream dbg;
  
  // GEN INFO do only if running on MC data
  if(doMCData_) {
    
    Handle<reco::GenParticleCollection>  genParticles;
    iEvent.getByLabel(genTag_, genParticles);   
    
    int pcount=0;
    for( size_t i = 0; i < genParticles->size(); ++ i ) {
      const reco::Candidate& pCand = (*genParticles)[ i ];
      
      int st = pCand.status();  
      
      if (st==3) {
  	int status = 3;
      } else { // store photons of status 1 
  	if ( (abs(pCand.pdgId()) == 22) ) {
  	  
  	  genPhotIds[pcount]    = pCand.pdgId();
  	  genPhotStatus[pcount] = pCand.status();
  	  genPhotE[pcount]      = pCand.energy();
  	  genPhotPx[pcount]     = pCand.px();
  	  genPhotPy[pcount]     = pCand.py();
  	  genPhotPz[pcount]     = pCand.pz();
  	  
  	  if (pCand.numberOfMothers() > 0 ) { 
  	    const reco::Candidate * mom = pCand.mother();
  	    while (mom->pdgId() == pCand.pdgId()) { mom = mom->mother(); }
  	    
  	    for( size_t j = 0; j < i; ++ j ) {
  	      const Candidate * ref = &((*genParticles)[j]);
  	      if (ref == mom) { genPhotRefs[pcount] = ref->pdgId(); }
  	      //if (ref == mom) { genPhotRefs[pcount] = j; }
  	    }  
  	  } else { genPhotRefs[pcount]=-999;}
  	  pcount++;
  	}
      }
    }
    genPhotLength = pcount;
  }
  

  /*
   *Get the information on all the photons
   *
   */

  // get the photons
  edm::Handle< std::vector<pat::Photon> > photHandle;
  iEvent.getByLabel(photTag_, photHandle);
  if ( !photHandle.isValid() ) {
    edm::LogWarning("PhotonEvent") << "No Photon results for InputTag " << photTag_;
    return false;
  }

  
  edm::LogVerbatim("PhotonEvent") << " start reading in photons " << std::endl;
  // Add the photons
  m_Nphot = photHandle->size();
  int ph = 0;
  if ( m_Nphot > 50 ) {
    m_Nphot = 50;
    if (debug_) printf("Photon/%d      E     Et    Pt    Px    Py    Pz    Eta    Phi\n",m_Nphot);
  }
  for (int i=0;i<m_Nphot;i++) {
    if ( ((*photHandle)[i].pt() > photMinEt_) && !((*photHandle)[i].eta() > photMaxEta_) ) {
      edm::LogVerbatim("PhotonEvent") << " looping over good photons " << std::endl;      
      m_PhotE[i]   = (*photHandle)[i].energy();
      m_PhotEt[i]  = (*photHandle)[i].et();
      m_PhotPt[i]  = (*photHandle)[i].pt();
      m_PhotPx[i]  = (*photHandle)[i].momentum().X();
      m_PhotPy[i]  = (*photHandle)[i].momentum().Y();
      m_PhotPz[i]  = (*photHandle)[i].momentum().Z();
      m_PhotEta[i] = (*photHandle)[i].eta();
      m_PhotPhi[i] = (*photHandle)[i].phi();
      if (debug_) printf("%6d   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f\n", \
	     i,m_PhotE[i],m_PhotEt[i],m_PhotPt[i],m_PhotPx[i],m_PhotPy[i],m_PhotPz[i],m_PhotEta[i],m_PhotPhi[i]);

      m_PhotTrkIso[i]  = (*photHandle)[i].trackIso();
      m_PhotECalIso[i] = (*photHandle)[i].ecalIso();
      m_PhotHCalIso[i] = (*photHandle)[i].hcalIso();
      m_PhotAllIso[i]  = (*photHandle)[i].caloIso();

      //m_PhotLooseEM[i] = (*photHandle)[i].photonID("PhotonCutBasedIDLooseEM");
      m_PhotLoosePhoton[i] = (*photHandle)[i].photonID("PhotonCutBasedIDLoose");
      m_PhotTightPhoton[i] = (*photHandle)[i].photonID("PhotonCutBasedIDTight");
      
      // GenPhoton info
      if (doMCData_) {
      	//reco::Particle* part = const_cast<reco::Particle*>( (*photHandle)[i].genPhoton() );
      	const reco::Candidate* candPhot = (*photHandle)[i].genPhoton();
      	if (debug_) printf("      GenPhoton      E     Et    Pt    Px    Py    Pz    PdgId    Mother\n");
      	if ( candPhot ) {
      	  m_GenPhotPdgId[i] = candPhot->pdgId();
      	  m_GenPhotPx[i]    = candPhot->px();
      	  m_GenPhotPy[i]    = candPhot->py();
      	  m_GenPhotPz[i]    = candPhot->pz();
      	  m_GenPhotPt[i]    = candPhot->pt();
      	  m_GenPhotEt[i]    = candPhot->et();
      	  m_GenPhotE[i]     = candPhot->energy();
      	  const reco::Candidate* photMother = candPhot->mother();
      	  if ( photMother ) {
      	    while (photMother->pdgId() == candPhot->pdgId()) photMother = photMother->mother();
      	    if ( photMother ) {
      	      m_GenPhotMother[i] = photMother->pdgId();
      	      //if ( cand->mother()->pdgId() ==  cand->pdgId()) 
      	      //  {
      	      //	m_GenPhotMother[i] = cand->mother()->mother()->pdgId();
      	      //  }
      	    }
      	  }
      	}
      	else {
      	  m_GenPhotPdgId[i]  = -999;
      	  m_GenPhotMother[i] = -999;
      	  m_GenPhotPx[i]     = -999;
      	  m_GenPhotPy[i]     = -999;
      	  m_GenPhotPz[i]     = -999;
      	  m_GenPhotPt[i]     = -999;
      	  m_GenPhotEt[i]     = -999;
      	  m_GenPhotE[i]      = -999;
      	}
      	if (debug_) printf("      %6d   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f\n", \
      	       i,m_GenPhotE[i],m_GenPhotEt[i],m_GenPhotPt[i],m_GenPhotPx[i],m_GenPhotPy[i],m_GenPhotPz[i],m_GenPhotPdgId[i],m_GenPhotMother[i]);
      }
      ++ph;
    }
  } // loop over pat::Photons
  m_Nphot = ph;
  
  
  // Fill the tree only if all preselection conditions are met
  //if(preselection) //mPreselection->Fill();
  //mPhotonData->Fill();
  return photon_result;
}

//________________________________________________________________________________________
void 
PhotonAnalyzer::beginJob(const edm::EventSetup&) {}

//________________________________________________________________________________________
void 
PhotonAnalyzer::endJob() {

}

//________________________________________________________________________________________
void
PhotonAnalyzer::initTuple() {

  std::ostringstream variables; // Container for all variables
  
  // 1. Event variables
  variables << "weight:process";
  
  // Add the branches

  //add photons
  mPhotonData->Branch("PhotN",   &m_Nphot,   "PhotN/int");  
  mPhotonData->Branch("PhotE",    m_PhotE,   "PhotE[PhotN]/double");
  mPhotonData->Branch("PhotEt",   m_PhotEt,  "PhotEt[PhotN]/double");
  mPhotonData->Branch("PhotPt",   m_PhotPt,  "PhotPt[PhotN]/double");
  mPhotonData->Branch("PhotPx",   m_PhotPx,  "PhotPx[PhotN]/double");
  mPhotonData->Branch("PhotPy",   m_PhotPy,  "PhotPy[PhotN]/double");
  mPhotonData->Branch("PhotPz",   m_PhotPz,  "PhotPz[PhotN]/double");
  mPhotonData->Branch("PhotEta",  m_PhotEta, "PhotEta[PhotN]/double");
  mPhotonData->Branch("PhotPhi",  m_PhotPhi, "PhotPhi[PhotN]/double");

  mPhotonData->Branch("PhotTrkIso",  m_PhotTrkIso,  "PhotTrkIso[PhotN]/double");
  mPhotonData->Branch("PhotECalIso", m_PhotECalIso, "PhotECalIso[PhotN]/double");
  mPhotonData->Branch("PhotHCalIso", m_PhotHCalIso, "PhotHCalIso[PhotN]/double");
  mPhotonData->Branch("PhotAllIso",  m_PhotAllIso,  "PhotAllIso[PhotN]/double");

  //mPhotonData->Branch("Phot_isccPhotAssoc", m_ccPhotAssoc,     "ccPhotAssoc[PhotN]/bool");
  //mPhotonData->Branch("PhotLooseEM",        m_PhotLooseEM,     "PhotLooseEM[PhotN]/bool");
  mPhotonData->Branch("PhotLoosePhoton",    m_PhotLoosePhoton, "PhotLoosePhoton[PhotN]/bool");
  mPhotonData->Branch("PhotTightPhoton",    m_PhotTightPhoton, "PhotTightPhoton[PhotN]/bool");

  if (doMCData_) {
    //from reco::candidate
    mPhotonData->Branch("PhotGenPdgId",  m_GenPhotPdgId,  "PhotGenPdgId[PhotN]/double");
    mPhotonData->Branch("PhotGenMother", m_GenPhotMother, "PhotGenMother[PhotN]/double");
    mPhotonData->Branch("PhotGenPx",     m_GenPhotPx,     "PhotGenPx[PhotN]/double");
    mPhotonData->Branch("PhotGenPy",     m_GenPhotPy,     "PhotGenPy[PhotN]/double");
    mPhotonData->Branch("PhotGenPz",     m_GenPhotPz,     "PhotGenPz[PhotN]/double");
    mPhotonData->Branch("PhotGenPt",     m_GenPhotPt,     "PhotGenPt[PhotN]/double");
    mPhotonData->Branch("PhotGenEt",     m_GenPhotEt,     "PhotGenEt[PhotN]/double");
    mPhotonData->Branch("PhotGenE",      m_GenPhotE,      "PhotGenE[PhotN]/double");
    //from genParticles
    mPhotonData->Branch("genPhotN",     &genPhotLength, "genPhotN/int");
    mPhotonData->Branch("genPhotId",     genPhotIds,    "genPhotIds[genPhotN]/int");
    mPhotonData->Branch("genPhotMother", genPhotRefs,   "genPhotRefs[genPhotN]/int");
    mPhotonData->Branch("genPhotStatus", genPhotStatus, "genPhotStatus[genPhotN]/int");
    mPhotonData->Branch("genPhotE",      genPhotE,      "genPhotE[genPhotN]/float");
    mPhotonData->Branch("genPhotPx",     genPhotPx,     "genPhotPx[genPhotN]/float");
    mPhotonData->Branch("genPhotPy",     genPhotPy,     "genPhotPy[genPhotN]/float");
    mPhotonData->Branch("genPhotPz",     genPhotPz,     "genPhotPz[genPhotN]/float");

  }

  edm::LogInfo("PhotonEvent") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/PluginManager/interface/ModuleDef.h"
//#include "FWCore/Framework/interface/ModuleFactory.h"
//
//DEFINE_EDM_PLUGIN(PhotonAnalyzer);
