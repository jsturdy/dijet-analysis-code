
// -*- C++ -*-
//
// Package:    DiJetAnalysis
// Class:      LeptonAnalyzer
// 
/**\class LeptonAnalyzer LeptonAnalyzer.cc JSturdy/DiJetAnalysis/src/LeptonAnalyzer.cc

Description: Variable collector/ntupler for SUSY search with Jets + MET

*/
//
// Original Author:  Jared Sturdy
//         Created:  Fri Jan 29 16:10:31 PDT 2010
// $Id: LeptonAnalyzer.cc,v 1.2 2010/03/29 11:19:36 sturdy Exp $
//
//


//#include "SusyAnalysis/EventSelector/interface/BJetEventSelector.h"
#include "JSturdy/DiJetAnalysis/interface/LeptonAnalyzer.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#define _USE_MATH_DEFINES
//#include <math.h>
#include <TMath.h>
#include <sstream>

//________________________________________________________________________________________
LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& pset, TTree* tmpAllData)
{ 

  mLeptonData = tmpAllData;
  leptonParams = pset;

  //defaults
  elecMaxEta_ = 3.0;
  elecMaxEt_  = 15.;
  elecMinEt_  = 5.;
  elecRelIso_ = 0.5;

  muonMaxEta_ = 3.0;
  muonMaxEt_  = 25.;
  muonMinEt_  = 5.;
  muonRelIso_ = 0.1;

  //tauMaxEta_ = 12.;
  //tauMaxEt_  = 50.;
  //tauMinEt_  = 5.;
  //tauRelIso_ = 0.5;

  doMCData_ = true;
  debug_    = 0;

  // Read in parameters from the config file
  if (leptonParams.exists("elecMaxEta")) elecMaxEta_ = leptonParams.getUntrackedParameter<double>("elecMaxEta");
  if (leptonParams.exists("elecMaxEt"))  elecMaxEt_  = leptonParams.getUntrackedParameter<double>("elecMaxEt");
  if (leptonParams.exists("elecMinEt"))  elecMinEt_  = leptonParams.getUntrackedParameter<double>("elecMinEt");
  if (leptonParams.exists("elecRelIso")) elecRelIso_ = leptonParams.getUntrackedParameter<double>("elecRelIso");

  if (leptonParams.exists("muonMaxEta")) muonMaxEta_ = leptonParams.getUntrackedParameter<double>("muonMaxEta");
  if (leptonParams.exists("muonMaxEt"))  muonMaxEt_  = leptonParams.getUntrackedParameter<double>("muonMaxEt");
  if (leptonParams.exists("muonMinEt"))  muonMinEt_  = leptonParams.getUntrackedParameter<double>("muonMinEt");
  if (leptonParams.exists("muonRelIso")) muonRelIso_ = leptonParams.getUntrackedParameter<double>("muonRelIso");

  //if (leptonParams.exists("tauMaxEta")) tauMaxEta_ = leptonParams.getUntrackedParameter<double>("tauMaxEta");
  //if (leptonParams.exists("tauMaxEt"))  tauMaxEt_  = leptonParams.getUntrackedParameter<double>("tauMaxEt");
  //if (leptonParams.exists("tauMinEt"))  tauMinEt_  = leptonParams.getUntrackedParameter<double>("tauMinEt");
  //if (leptonParams.exists("tauRelIso")) tauRelIso_ = leptonParams.getUntrackedParameter<double>("tauRelIso");

  if (leptonParams.exists("doMCLeps"))   doMCData_   = leptonParams.getUntrackedParameter<bool>("doMCLeps");
  if (doMCData_) 
    if (leptonParams.exists("genLepTag"))    
      genTag_   = leptonParams.getUntrackedParameter<edm::InputTag>("genLepTag");
  if (leptonParams.exists("debugLeps"))   debug_   = leptonParams.getUntrackedParameter<int>("debugLeps");
 
  // get the data tags
  elecTag_   = leptonParams.getUntrackedParameter<edm::InputTag>("elecTag");
  pfelecTag_ = leptonParams.getUntrackedParameter<edm::InputTag>("pfelecTag");

  muonTag_   = leptonParams.getUntrackedParameter<edm::InputTag>("muonTag");
  pfmuonTag_ = leptonParams.getUntrackedParameter<edm::InputTag>("pfmuonTag");

  //tauTag_   = leptonParams.getUntrackedParameter<edm::InputTag>("tauTag");
  //pftauTag_ = leptonParams.getUntrackedParameter<edm::InputTag>("pftauTag");

  localPi = acos(-1.0);

  // Initialise ntuple branches
  initTuple();
}


//________________________________________________________________________________________
LeptonAnalyzer::~LeptonAnalyzer() {}


//________________________________________________________________________________________
// Method called to for each event
//void
//LeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
bool
LeptonAnalyzer::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;
  using namespace edm;

  m_elecVeto      = false;
  m_muonVeto      = false;
  bool lepton_result = true;

  edm::LogVerbatim("LeptonEvent") << " Start  " << std::endl;

  std::ostringstream dbg;

  // GEN INFO do only if running on MC data
  if(doMCData_) {
    //get pthat of process
    m_Pthat = -999.;
    
    Handle<double> genEventScale;
    iEvent.getByLabel( "genEventScale", genEventScale );
    if ( genEventScale.isValid() ) m_Pthat = *genEventScale;
    
    Handle<reco::GenParticleCollection>  genParticles;
    iEvent.getByLabel(genTag_, genParticles);   
    
    int count=0; int lcount=0; ; int pcount=0; // int tcount=0;
    printf("Status   genpart/%d   PdgId   genE   genPx   genPy   genPz   genMother\n",genParticles->size());
    for( size_t i = 0; i < genParticles->size(); ++ i ) {
      const reco::Candidate& pCand = (*genParticles)[ i ];

      int st = pCand.status();  

      //get status 3 particles
      if (st==3) {
	genIds[count]    = pCand.pdgId();
	genStatus[count] = pCand.status();
	genE[count]      = pCand.energy();
	genPx[count]     = pCand.px();
	genPy[count]     = pCand.py();
	genPz[count]     = pCand.pz();
      
	if (pCand.numberOfMothers() > 0 ) { 
	  const reco::Candidate * mom = pCand.mother();
	  while (mom->pdgId() == pCand.pdgId()) {mom = mom->mother(); }
	  
	  for( size_t j = 0; j < i; ++ j ) {
	    const Candidate * ref = &((*genParticles)[j]);
	    if (ref == mom) { genRefs[count] = ref->pdgId(); } //return mother's pdgId
	    //if (ref == mom) { genRefs[count] = j; } //point to particle that is reference
	  }  
	} else { genRefs[count]=-999;}
	if (debug_>1)  printf("%2d       %4d       %4d       %4.2f     %4.2f    %4.2f    %4.2f     %4d\n", \
	       genStatus[count],count,genIds[count],genE[count],genPx[count],genPy[count],genPz[count],genRefs[count]);
	++count;
      }
      else { // store also electrons or muons of status 1 
	//if ( (abs(pCand.pdgId()) == 11) || (abs(pCand.pdgId()) == 13) ) {
	if ( (abs(pCand.pdgId()) == 11) || (abs(pCand.pdgId()) == 13) || (abs(pCand.pdgId()) == 22) ) {
	  
	  genLepIds[lcount]    = pCand.pdgId();
	  genLepStatus[lcount] = pCand.status();
	  genLepE[lcount]      = pCand.energy();
	  genLepPx[lcount]     = pCand.px();
	  genLepPy[lcount]     = pCand.py();
	  genLepPz[lcount]     = pCand.pz();
	  
	  if (pCand.numberOfMothers() > 0 ) { 
	    const reco::Candidate * mom = pCand.mother();
	    while (mom->pdgId() == pCand.pdgId()) { mom = mom->mother(); }
	    
	    for( size_t j = 0; j < i; ++ j ) {
	      const reco::Candidate * ref = &((*genParticles)[j]);
	      if (ref == mom) { genLepRefs[lcount] = ref->pdgId(); }
	      //if (ref == mom) { genLepRefs[lcount] = j; }
	    }  
	  } else { genLepRefs[lcount]=-999;}
	  if (debug_>1)printf("%2d         %4d         %4d         %4.2f       %4.2f    %4.2f    %4.2f     %4d\n", \
		 genLepStatus[lcount],lcount,genLepIds[lcount],genLepE[lcount],genLepPx[lcount],genLepPy[lcount],genLepPz[lcount],genLepRefs[lcount]);
	  ++lcount;
	}
	//if ( (abs(pCand.pdgId()) == 22) ) {
	//  //  
	//  genPhotIds[pcount]    = pCand.pdgId();
	//  genPhotStatus[pcount] = pCand.status();
	//  genPhotE[pcount]      = pCand.energy();
	//  genPhotPx[pcount]     = pCand.px();
	//  genPhotPy[pcount]     = pCand.py();
	//  genPhotPz[pcount]     = pCand.pz();
	//  
	//  if (pCand.numberOfMothers() > 0 ) { 
	//    const reco::Candidate * mom = pCand.mother();
	//    while (mom->pdgId() == pCand.pdgId()) { mom = mom->mother(); }
	//    
	//    for( size_t j = 0; j < i; ++ j ) {
	//      const reco::Candidate * ref = &((*genParticles)[j]);
	//      if (ref == mom) { genPhotRefs[pcount] = ref->pdgId(); }
	//      //if (ref == mom) { genPhotRefs[pcount] = j; }
	//    }  
	//  } else { genPhotRefs[pcount]=-999;}
	//  if (debug_>1)printf("%2d         %4d         %4d         %4.2f       %4.2f    %4.2f    %4.2f     %4d\n", \
	//	 genPhotStatus[pcount],pcount,genPhotIds[pcount],genPhotE[pcount],genPhotPx[pcount],genPhotPy[pcount],genPhotPz[pcount],genPhotRefs[pcount]);
	//  ++pcount;
	//  /*
	//  //  genPhotIds[pcount]    = pCand.pdgId();
	//  //  genPhotStatus[pcount] = pCand.status();
	//  //  genPhotE[pcount]      = pCand.energy();
	//  //  genPhotPx[pcount]     = pCand.px();
	//  //  genPhotPy[pcount]     = pCand.py();
	//  //  genPhotPz[pcount]     = pCand.pz();
	//  //  
	//  //  if (pCand.numberOfMothers() > 0 ) { 
	//  //    const reco::Candidate * mom = pCand.mother();
	//  //    while (mom->pdgId() == pCand.pdgId()) { mom = mom->mother(); }
	//  //    
	//  //    for( size_t j = 0; j < i; ++ j ) {
	//  //      const Candidate * ref = &((*genParticles)[j]);
	//  //      if (ref == mom) { genPhotRefs[pcount] = ref->pdgId(); }
	//  //      //if (ref == mom) { genPhotRefs[pcount] = j; }
	//  //    }  
	//  //  } else { genPhotRefs[pcount]=-999;}
	//  //  pcount++;
	//  //  printf("%2d         2.2f       %2.2f       %2.2f       %2.2f    %2.2f    %2.2f     %2.2f", \
	//  //	 pcount,genIds[pcount],genStatus[pcount],genE[pcount],genPx[pcount],genPy[pcount],genPz[pcount],genMother[pcount]);
	//  //}
	//  */
	//}
      }
    }
    length = count;
    genLepLength = lcount;
    //genPhotLength = pcount;
  }
  
  /*
   *Get the information on all the electrons
   *
   */

  // get the electrons
  edm::Handle< std::vector<pat::Electron> > elecHandle;
  iEvent.getByLabel(elecTag_, elecHandle);
  if ( !elecHandle.isValid() ) {
    edm::LogWarning("LeptonEvent") << "No Electron results for InputTag " << elecTag_;
    return false;
  }

  
  edm::LogVerbatim("LeptonEvent") << " start reading in electrons " << std::endl;
  // Add the electrons
  m_Nelec = elecHandle->size();
  if ( m_Nelec > 50 ) m_Nelec = 50;
  int el = 0;
  if (debug_) printf("Elec/%d      E     Et    Pt    Px    Py    Pz    Eta    Phi\n",m_Nelec);
  for (int i=0;i<m_Nelec;i++){
    if ( ((*elecHandle)[i].pt()>elecMinEt_) && !((*elecHandle)[i].eta()>elecMaxEta_) ) {
      edm::LogVerbatim("LeptonEvent") << " looping over good electrons " << std::endl;
      m_ElecE[i]      = (*elecHandle)[i].energy();
      m_ElecEt[i]     = (*elecHandle)[i].et();
      m_ElecPt[i]     = (*elecHandle)[i].pt();
      m_ElecPx[i]     = (*elecHandle)[i].momentum().X();
      m_ElecPy[i]     = (*elecHandle)[i].momentum().Y();
      m_ElecPz[i]     = (*elecHandle)[i].momentum().Z();
      m_ElecEta[i]    = (*elecHandle)[i].eta();
      m_ElecPhi[i]    = (*elecHandle)[i].phi();
      m_ElecCharge[i] = (*elecHandle)[i].charge();
      if (debug_) printf("%6d   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f\n", \
	     i,m_ElecE[i],m_ElecEt[i],m_ElecPt[i],m_ElecPx[i],m_ElecPy[i],m_ElecPz[i],m_ElecEta[i],m_ElecPhi[i]);
      
      m_ElecTrkIso[i]  = (*elecHandle)[i].trackIso();
      m_ElecECalIso[i] = (*elecHandle)[i].ecalIso();
      m_ElecHCalIso[i] = (*elecHandle)[i].hcalIso() ;
      m_ElecAllIso[i]  = (*elecHandle)[i].caloIso() ;
      
      //m_ElecECalIsoDeposit[i]  = (*elecHandle)[i].ecalIsoDeposit()->candEnergy() ;
      //m_ElecHCalIsoDeposit[i]  = (*elecHandle)[i].hcalIsoDeposit()->candEnergy() ;
      
      m_ElecIdLoose[i]    = (*elecHandle)[i].electronID("eidLoose");
      m_ElecIdTight[i]    = (*elecHandle)[i].electronID("eidTight");
      m_ElecIdRobLoose[i] = (*elecHandle)[i].electronID("eidRobustLoose");
      m_ElecIdRobTight[i] = (*elecHandle)[i].electronID("eidRobustTight"); 
      
      m_ElecCaloEnergy[i] = (*elecHandle)[i].caloEnergy();
      m_ElecHOverE[i]     = (*elecHandle)[i].hadronicOverEm();
      m_ElecVx[i]         = (*elecHandle)[i].vx();
      m_ElecVy[i]         = (*elecHandle)[i].vy();
      m_ElecVz[i]         = (*elecHandle)[i].vz();
      
      edm::LogVerbatim("LeptonEvent") << " before asking for gsfTrack " << std::endl;
      m_ElecD0[i]               = (*elecHandle)[i].gsfTrack()->d0();
      m_ElecDz[i]               = (*elecHandle)[i].gsfTrack()->dz();
      m_ElecChargeMode[i]       = (*elecHandle)[i].gsfTrack()->chargeMode();	
      m_ElecPtTrkMode[i]        = (*elecHandle)[i].gsfTrack()->ptMode();
      m_ElecQOverPErrTrkMode[i] = (*elecHandle)[i].gsfTrack()->qoverpModeError();
      m_ElecCharge[i]           = (*elecHandle)[i].gsfTrack()->charge();
      m_ElecPtTrk[i]            = (*elecHandle)[i].gsfTrack()->pt();
      m_ElecQOverPErrTrk[i]     = (*elecHandle)[i].gsfTrack()->qoverpError();
      m_ElecNormChi2[i]         = (*elecHandle)[i].gsfTrack()->normalizedChi2();
      m_ElecLostHits[i]         = (*elecHandle)[i].gsfTrack()->lost();
      m_ElecValidHits[i]        = (*elecHandle)[i].gsfTrack()->found();
      
      edm::LogVerbatim("LeptonEvent") << " before asking for trackMomentumAtVtx " << std::endl;
    
      m_ElecEtaTrk[i] = (*elecHandle)[i].trackMomentumAtVtx().Eta();
      m_ElecPhiTrk[i] = (*elecHandle)[i].trackMomentumAtVtx().Phi();
      
      // Added protection statement, against missing SuperCluster collection in 2_1_X PatLayer1 samples
      try { 
	m_ElecWidthClusterEta[i] = (*elecHandle)[i].superCluster()->etaWidth();
	m_ElecWidthClusterPhi[i] = (*elecHandle)[i].superCluster()->phiWidth();
      } catch ( const cms::Exception& e ) {
	m_ElecWidthClusterEta[i]=-999.;
	m_ElecWidthClusterPhi[i]=-999.;
	std::stringstream ss;
	ss << " cms::Exception caught!"
	   << " Invalid edm::Ref<reco::SuperCluster> returned from pat::Electron!" 
	   << std::endl 
	   << " Setting ClusterEta and ClusterPhi to -999.!" 
	   << std::endl 
	   << " Output from cms::Exception::what():"
	   << std::endl 
	   << e.what();
	edm::LogWarning("LeptonEvent") << ss.str();
      }
      
      m_ElecPinTrk[i] = sqrt((*elecHandle)[i].trackMomentumAtVtx().Mag2());
      m_ElecPoutTrk[i] = sqrt((*elecHandle)[i].trackMomentumOut().Mag2());
      
      //get associated gen particle information
      if (doMCData_) {
      	const reco::Candidate* candElec = (*elecHandle)[i].genLepton();
      	if ( candElec ) {
      	  m_GenElecPdgId[i] = candElec->pdgId();
      	  m_GenElecPx[i]    = candElec->px();
      	  m_GenElecPy[i]    = candElec->py();
      	  m_GenElecPz[i]    = candElec->pz();
      	  m_GenElecPt[i]    = candElec->pt();
      	  m_GenElecEt[i]    = candElec->et();
      	  m_GenElecE[i]     = candElec->energy();
        
      	  const reco::Candidate* elecMother = candElec->mother();
      	  if( elecMother ) {
      	    while (elecMother->pdgId() == candElec->pdgId()) elecMother = elecMother->mother();
      	    if ( elecMother ) {
      	      m_GenElecMother[i] = (*elecHandle)[i].genLepton()->mother()->pdgId();
      	      //if ( (*elecHandle)[i].genLepton()->mother()->pdgId() ==  (*elecHandle)[i].genLepton()->pdgId()) 
      	      //  {
      	      //	m_GenElecMother[i] = (*elecHandle)[i].genLepton()->mother()->mother()->pdgId();
      	      //  }
      	    }
      	  }
      	}
      	else {
      	  m_GenElecPdgId[i]  = -999.;
      	  m_GenElecMother[i] = -999.;
      	  m_GenElecPx[i]     = -999.;
      	  m_GenElecPy[i]     = -999.;
      	  m_GenElecPz[i]     = -999.;
      	  m_GenElecPt[i]     = -999.;
      	  m_GenElecEt[i]     = -999.;
      	  m_GenElecE[i]      = -999.;
      	}
      }
      double elecIsoReq = (m_ElecTrkIso[i]+m_ElecECalIso[i]+m_ElecHCalIso[i])/m_ElecPt[i];
      if ( elecIsoReq  > elecRelIso_) m_elecVeto = m_elecVeto || true;
      if ( m_ElecPt[i] > elecMaxEt_ ) m_elecVeto = m_elecVeto || true;
      ++el;
    }
  }//end loop over Electrons
  m_Nelec = el;

  /*
   * get the muons
   *
   */

  edm::Handle< std::vector<pat::Muon> > muonHandle;
  iEvent.getByLabel(muonTag_, muonHandle);
  if ( !muonHandle.isValid() ) {
    edm::LogWarning("LeptonEvent") << "No Muon results for InputTag " << muonTag_;
    return false;
  }
  
  edm::LogVerbatim("LeptonEvent") << " start reading in muons " << std::endl;


  // Add the muons
  m_Nmuon= muonHandle->size();
  if ( m_Nmuon > 50 ) m_Nmuon = 50;
  int mu = 0;
  if (debug_) printf("Muon/%d      E     Et    Pt    Px    Py    Pz    Eta    Phi\n",m_Nmuon);
  for (int i=0;i<m_Nmuon;i++){
    if ( ((*muonHandle)[i].pt()>muonMinEt_) && !((*muonHandle)[i].eta()>muonMaxEta_) ) {
      edm::LogVerbatim("LeptonEvent") << " looping over good muons " << std::endl;      
      m_MuonPt[i]  = (*muonHandle)[i].pt();
      m_MuonE[i]   = (*muonHandle)[i].energy();
      m_MuonEt[i]  = (*muonHandle)[i].et();
      m_MuonPx[i]  = (*muonHandle)[i].momentum().X();
      m_MuonPy[i]  = (*muonHandle)[i].momentum().Y();
      m_MuonPz[i]  = (*muonHandle)[i].momentum().Z();
      m_MuonEta[i] = (*muonHandle)[i].eta();
      m_MuonPhi[i] = (*muonHandle)[i].phi();
      if (debug_) printf("%6d   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f   %2.2f\n", \
	     i,m_MuonE[i],m_MuonEt[i],m_MuonPt[i],m_MuonPx[i],m_MuonPy[i],m_MuonPz[i],m_MuonEta[i],m_MuonPhi[i]);

      //Muon isolation variables
      m_MuonTrkIso[i]   = (*muonHandle)[i].trackIso();
      m_MuonCharge[i]   = (*muonHandle)[i].charge();
      m_MuonECalIso[i]  = (*muonHandle)[i].ecalIso();
      m_MuonHCalIso[i]  = (*muonHandle)[i].hcalIso() ;
      m_MuonAllIso[i]   = (*muonHandle)[i].caloIso() ;

      //Muon classification variables
      m_MuonIsGlobal[i]               = muon::isGoodMuon((*muonHandle)[i], muon::AllGlobalMuons);
      m_MuonIsStandAlone[i]           = muon::isGoodMuon((*muonHandle)[i], muon::AllStandAloneMuons);
      m_MuonIsTracker[i]              = muon::isGoodMuon((*muonHandle)[i], muon::AllTrackerMuons);
      m_MuonIsGlobalTight[i]          = muon::isGoodMuon((*muonHandle)[i], muon::GlobalMuonPromptTight);
      m_MuonIsTMLastStationLoose[i]   = muon::isGoodMuon((*muonHandle)[i], muon::TMLastStationLoose);
      m_MuonTMLastStationTight[i]     = muon::isGoodMuon((*muonHandle)[i], muon::TMLastStationTight);
      m_MuonTM2DCompatibilityLoose[i] = muon::isGoodMuon((*muonHandle)[i], muon::TM2DCompatibilityLoose);
      m_MuonTM2DCompatibilityTight[i] = muon::isGoodMuon((*muonHandle)[i], muon::TM2DCompatibilityTight);
    
      m_MuonECalIsoDeposit[i]  = (*muonHandle)[i].ecalIsoDeposit()->candEnergy() ;
      m_MuonHCalIsoDeposit[i]  = (*muonHandle)[i].hcalIsoDeposit()->candEnergy() ;

      //Muon Vertex information
      // Vertex info is stored only for GlobalMuons (combined muons)
      if((*muonHandle)[i].isGlobalMuon() && (*muonHandle)[i].combinedMuon().isNonnull()){ 

	m_MuonCombChi2[i] = (*muonHandle)[i].combinedMuon()->chi2();
	m_MuonCombNdof[i] = (*muonHandle)[i].combinedMuon()->ndof();

	m_MuonCombVx[i] = (*muonHandle)[i].combinedMuon()->vx();
	m_MuonCombVy[i] = (*muonHandle)[i].combinedMuon()->vy();
	m_MuonCombVz[i] = (*muonHandle)[i].combinedMuon()->vz();
	m_MuonCombD0[i] = (*muonHandle)[i].combinedMuon()->d0();
	m_MuonCombDz[i] = (*muonHandle)[i].combinedMuon()->dz();

      } else {
	m_MuonCombVx[i] = 999.;
	m_MuonCombVy[i] = 999.;
	m_MuonCombVz[i] = 999.;
	m_MuonCombD0[i] = 999.;
	m_MuonCombDz[i] = 999.;
      }

      //Standalone muon information
      if((*muonHandle)[i].isStandAloneMuon() && (*muonHandle)[i].standAloneMuon().isNonnull()){
	m_MuonStandValidHits[i]   = (*muonHandle)[i].standAloneMuon()->found();
	m_MuonStandLostHits[i]    = (*muonHandle)[i].standAloneMuon()->lost();
	m_MuonStandPt[i]          = (*muonHandle)[i].standAloneMuon()->pt();
	m_MuonStandPz[i]          = (*muonHandle)[i].standAloneMuon()->pz();
	m_MuonStandP[i]           = (*muonHandle)[i].standAloneMuon()->p();
	m_MuonStandEta[i]         = (*muonHandle)[i].standAloneMuon()->eta();
	m_MuonStandPhi[i]         = (*muonHandle)[i].standAloneMuon()->phi();
	m_MuonStandChi[i]         = (*muonHandle)[i].standAloneMuon()->chi2();
	m_MuonStandCharge[i]      = (*muonHandle)[i].standAloneMuon()->charge();
	m_MuonStandQOverPError[i] = (*muonHandle)[i].standAloneMuon()->qoverpError();
      } 
      else{
	m_MuonStandValidHits[i]   = 999.;
	m_MuonStandLostHits[i]    = 999.;
	m_MuonStandPt[i]          = 999.;
	m_MuonStandPz[i]          = 999.;
	m_MuonStandP[i]           = 999.;
	m_MuonStandEta[i]         = 999.;
	m_MuonStandPhi[i]         = 999.;
	m_MuonStandChi[i]         = 999.;
	m_MuonStandCharge[i]      = 999.;
	m_MuonStandQOverPError[i] = 999.;
      }

      //Muon tracking information
      if((*muonHandle)[i].isTrackerMuon() && (*muonHandle)[i].track().isNonnull()){
	m_MuonTrkChiNorm[i]     = (*muonHandle)[i].track()->normalizedChi2();
	m_MuonTrkValidHits[i]   = (*muonHandle)[i].track()->found();
	m_MuonTrkLostHits[i]    = (*muonHandle)[i].track()->lost();
	m_MuonTrkD0[i]          = (*muonHandle)[i].track()->d0();
	m_MuonTrkPt[i]          = (*muonHandle)[i].track()->pt();
	m_MuonTrkPz[i]          = (*muonHandle)[i].track()->pz();
	m_MuonTrkP[i]           = (*muonHandle)[i].track()->p();
	m_MuonTrkEta[i]         = (*muonHandle)[i].track()->eta();
	m_MuonTrkPhi[i]         = (*muonHandle)[i].track()->phi();
	m_MuonTrkChi[i]         = (*muonHandle)[i].track()->chi2();
	m_MuonTrkCharge[i]      = (*muonHandle)[i].track()->charge();
	m_MuonTrkQOverPError[i] = (*muonHandle)[i].track()->qoverpError();
	//  m_MuonTrkOuterZ[i]=(*muonHandle)[i].track()->outerZ();
	//  m_MuonTrkOuterR[i]=(*muonHandle)[i].track()->outerRadius();

      }
      else{
	m_MuonTrkChiNorm[i]    = 999.;
	m_MuonTrkValidHits[i]  = 999.;
	m_MuonTrkLostHits[i]   = 999.;
	m_MuonTrkPt[i]         = 999.;
	m_MuonTrkPz[i]         = 999.;
	m_MuonTrkP[i]          = 999.;
	m_MuonTrkEta[i]        = 999.;
	m_MuonTrkPhi[i]        = 999.;
	m_MuonTrkChi[i]        = 999.;
	m_MuonTrkCharge[i]     = 999.;
	m_MuonTrkQOverPError[i]= 999.;
	m_MuonTrkOuterZ[i]     = 999.;
	m_MuonTrkOuterR[i]     = 999.;
      }

      //Muon gen particle association variables
      if (doMCData_) {
      	const reco::Candidate* candMuon = (*muonHandle)[i].genLepton();
      	if ( candMuon ) {
      	  m_GenMuonPdgId[i] = candMuon->pdgId();
      	  m_GenMuonPx[i]    = candMuon->px();
      	  m_GenMuonPy[i]    = candMuon->py();
      	  m_GenMuonPz[i]    = candMuon->pz();
      	  m_GenMuonPt[i]    = candMuon->pt();
      	  m_GenMuonEt[i]    = candMuon->et();
      	  m_GenMuonE[i]     = candMuon->energy();
      	  
      	  const reco::Candidate* muonMother = candMuon->mother();
      	  if( muonMother ) {
      	    while (muonMother->pdgId() == candMuon->pdgId()) muonMother = muonMother->mother();
      	    if ( muonMother ) {
      	      m_GenMuonMother[i] = (*muonHandle)[i].genLepton()->mother()->pdgId();
      	      //if ( (*muonHandle)[i].genLepton()->mother()->pdgId() ==  (*muonHandle)[i].genLepton()->pdgId()) 
      	      //  {
      	      //	m_GenMuonMother[i] = (*muonHandle)[i].genLepton()->mother()->mother()->pdgId();
      	      //  }
      	    }
      	  }
      	}
      	
      	else{
      	  m_GenMuonPdgId[i]  = -999.;
      	  m_GenMuonMother[i] = -999.;
      	  m_GenMuonPx[i]     = -999.;
      	  m_GenMuonPy[i]     = -999.;
      	  m_GenMuonPz[i]     = -999.;
      	  m_GenMuonPt[i]     = -999.;
      	  m_GenMuonEt[i]     = -999.;
      	  m_GenMuonE[i]      = -999.;
      	}
      }
      double muonIsoReq = (m_MuonTrkIso[i]+m_MuonECalIso[i]+m_MuonHCalIso[i])/m_MuonPt[i];
      if ( muonIsoReq  > muonRelIso_) m_muonVeto = m_muonVeto || true;
      if ( m_MuonPt[i] > muonMaxEt_)  m_muonVeto = m_muonVeto || true;
      ++mu;
    }
  }// end loop over muons
  m_Nmuon = mu;
  // return true when none of the events have leptons above threshold
  lepton_result = !(m_elecVeto || m_muonVeto);
  //mLeptonData->Fill();
  return lepton_result;
}

//________________________________________________________________________________________
void 
LeptonAnalyzer::beginJob(const edm::EventSetup&) {}

//________________________________________________________________________________________
void 
LeptonAnalyzer::endJob() {

}

//________________________________________________________________________________________
void
LeptonAnalyzer::initTuple() {

  std::ostringstream variables; // Container for all variables
  
  // 1. Event variables
  variables << "weight:process";


  // Add the branches

  //add electrons
  mLeptonData->Branch("ElecVeto", &m_elecVeto, "ElecVeto/bool");
  //General electron information
  mLeptonData->Branch("ElecN",     &m_Nelec,      "ElecN/int");  
  mLeptonData->Branch("ElecE",      m_ElecE,      "ElecE[ElecN]/double");
  mLeptonData->Branch("ElecEt",     m_ElecEt,     "ElecEt[ElecN]/double");
  mLeptonData->Branch("ElecPt",     m_ElecPt,     "ElecPt[ElecN]/double");
  mLeptonData->Branch("ElecPx",     m_ElecPx,     "ElecPx[ElecN]/double");
  mLeptonData->Branch("ElecPy",     m_ElecPy,     "ElecPy[ElecN]/double");
  mLeptonData->Branch("ElecPz",     m_ElecPz,     "ElecPz[ElecN]/double");
  mLeptonData->Branch("ElecEta",    m_ElecEta,    "ElecEta[ElecN]/double");
  mLeptonData->Branch("ElecPhi",    m_ElecPhi,    "ElecPhi[ElecN]/double");
  mLeptonData->Branch("ElecCharge", m_ElecCharge, "ElecCharge[ElecN]/double");
  mLeptonData->Branch("ElecHOverE", m_ElecHOverE, "ElecHOverE[ElecN]/double");

  //Isolation and tracking variables
  mLeptonData->Branch("ElecTrkIso",     m_ElecTrkIso,   "ElecTrkIso[ElecN]/double");
  mLeptonData->Branch("ElecECalIso",    m_ElecECalIso,  "ElecECalIso[ElecN]/double");
  mLeptonData->Branch("ElecHCalIso",    m_ElecHCalIso,  "ElecHCalIso[ElecN]/double");
  mLeptonData->Branch("ElecAllIso",     m_ElecAllIso,   "ElecAllIso[ElecN]/double");
  mLeptonData->Branch("ElecTrkChiNorm", m_ElecNormChi2, "ElecTrkChiNorm[ElecN]/double");
  //mLeptonData->Branch("NIsoelec",      &m_NIsoelec,     "NIsoelec/int");  

  mLeptonData->Branch("ElecECalIsoDeposit", m_ElecECalIsoDeposit, "ElecECalIsoDeposit[ElecN]/double");
  mLeptonData->Branch("ElecHCalIsoDeposit", m_ElecHCalIsoDeposit, "ElecHCalIsoDeposit[ElecN]/double");

  //Electron identification values
  mLeptonData->Branch("ElecIdLoose",    m_ElecIdLoose,    "ElecIdLoose [ElecN]/double");
  mLeptonData->Branch("ElecIdTight",    m_ElecIdTight,    "ElecIdTight [ElecN]/double");
  mLeptonData->Branch("ElecIdRobLoose", m_ElecIdRobLoose, "ElecIdRobLoose [ElecN]/double");
  mLeptonData->Branch("ElecIdRobTight", m_ElecIdRobTight, "ElecIdRobTight [ElecN]/double");
  mLeptonData->Branch("ElecChargeMode", m_ElecChargeMode, "ElecChargeMode [ElecN]/double");
  mLeptonData->Branch("ElecPtMode",     m_ElecPtTrkMode,  "ElecPtMode [ElecN]/double");


  //Electron vertex information
  mLeptonData->Branch("ElecVx",     m_ElecVx,     "ElecVx[ElecN]/double");
  mLeptonData->Branch("ElecVy",     m_ElecVy,     "ElecVy[ElecN]/double");
  mLeptonData->Branch("ElecVz",     m_ElecVz,     "ElecVz[ElecN]/double");
  mLeptonData->Branch("ElecD0",     m_ElecD0,     "ElecD0[ElecN]/double");
  mLeptonData->Branch("ElecDz",     m_ElecDz,     "ElecDz[ElecN]/double");
  mLeptonData->Branch("ElecPtTrk",  m_ElecPtTrk,  "ElecPtTrk[ElecN]/double");

  //Additonal electron detector information
  //Electron tracking information
  mLeptonData->Branch("ElecQOverPErrTrkMode", m_ElecQOverPErrTrkMode, "ElecQOverPErrTrkMode [ElecN]/double");
  mLeptonData->Branch("ElecCaloEnergy",       m_ElecCaloEnergy,       "ElecCaloEnergy[ElecN]/double");
  mLeptonData->Branch("ElecQOverPErrTrk",     m_ElecQOverPErrTrk,     "ElecQOverPErrTrk[ElecN]/double");
  mLeptonData->Branch("ElecPinTrk",           m_ElecPinTrk,           "ElecPinTrk[ElecN]/double");
  mLeptonData->Branch("ElecPoutTrk",          m_ElecPoutTrk,          "ElecPoutTrk[ElecN]/double"); 
  mLeptonData->Branch("ElecLostHits",         m_ElecLostHits,         "ElecLostHits[ElecN]/double"); 
  mLeptonData->Branch("ElecValidHits",        m_ElecValidHits,        "ElecValidHits[ElecN]/double"); 
  mLeptonData->Branch("ElecNCluster",         m_ElecNCluster,         "ElecNCluster[ElecN]/double"); 
  mLeptonData->Branch("ElecEtaTrk",           m_ElecEtaTrk,           "ElecEtaTrk[ElecN]/double"); 
  mLeptonData->Branch("ElecPhiTrk",           m_ElecPhiTrk,           "ElecPhiTrk[ElecN]/double"); 
  mLeptonData->Branch("ElecWidthClusterEta",  m_ElecWidthClusterEta,  "ElecWidthClusterEta[ElecN]/double"); 
  mLeptonData->Branch("ElecWidthClusterPhi",  m_ElecWidthClusterPhi,  "ElecWidthClusterPhi[ElecN]/double"); 

  if (doMCData_) {
    //Generator level information stored in the electron object
    mLeptonData->Branch("ElecGenPdgId",  m_GenElecPdgId,  "ElecGenPdgId[ElecN]/double");
    mLeptonData->Branch("ElecGenMother", m_GenElecMother, "ElecGenMother[ElecN]/double");
    mLeptonData->Branch("ElecGenPx",     m_GenElecPx,     "ElecGenPx[ElecN]/double");
    mLeptonData->Branch("ElecGenPy",     m_GenElecPy,     "ElecGenPy[ElecN]/double");
    mLeptonData->Branch("ElecGenPz",     m_GenElecPz,     "ElecGenPz[ElecN]/double");
    mLeptonData->Branch("ElecGenPt",     m_GenElecPt,     "ElecGenPt[ElecN]/double");
    mLeptonData->Branch("ElecGenEt",     m_GenElecEt,     "ElecGenEt[ElecN]/double");
    mLeptonData->Branch("ElecGenE",      m_GenElecE,      "ElecGenE[ElecN]/double");
  }

  //add muons
  mLeptonData->Branch("MuonVeto", &m_muonVeto, "MuonVeto/bool");
  //General kinematic variables related to muons
  mLeptonData->Branch("MuonN",         &m_Nmuon,          "MuonN/int");  
  mLeptonData->Branch("MuonE",          m_MuonE,          "MuonE[MuonN]/double");
  mLeptonData->Branch("MuonEt",         m_MuonEt,         "MuonEt[MuonN]/double");
  mLeptonData->Branch("MuonPt",         m_MuonPt,         "MuonPt[MuonN]/double");
  mLeptonData->Branch("MuonPx",         m_MuonPx,         "MuonPx[MuonN]/double");
  mLeptonData->Branch("MuonPy",         m_MuonPy,         "MuonPy[MuonN]/double");
  mLeptonData->Branch("MuonPz",         m_MuonPz,         "MuonPz[MuonN]/double");
  mLeptonData->Branch("MuonEta",        m_MuonEta,        "MuonEta[MuonN]/double");
  mLeptonData->Branch("MuonPhi",        m_MuonPhi,        "MuonPhi[MuonN]/double");
  mLeptonData->Branch("MuonCharge",     m_MuonCharge,     "MuonCharge[MuonN]/double");

  //Muon isolation variables
  //mLeptonData->Branch("NIsomuon",      &m_NIsomuon,       "NIsomuon/int");  
  mLeptonData->Branch("MuonTrkIso",     m_MuonTrkIso,     "MuonTrkIso[MuonN]/double");
  mLeptonData->Branch("MuonECalIso",    m_MuonECalIso,    "MuonECalIso[MuonN]/double");
  mLeptonData->Branch("MuonHCalIso",    m_MuonHCalIso,    "MuonHCalIso[MuonN]/double");
  mLeptonData->Branch("MuonAllIso",     m_MuonAllIso,     "MuonAllIso[MuonN]/double");
  mLeptonData->Branch("MuonTrkChiNorm", m_MuonTrkChiNorm, "MuonTrkChiNorm[MuonN]/double");

  mLeptonData->Branch("MuonECalIsoDeposit", m_MuonECalIsoDeposit, "MuonECalIsoDeposit[MuonN]/double");
  mLeptonData->Branch("MuonHCalIsoDeposit", m_MuonHCalIsoDeposit, "MuonHCalIsoDeposit[MuonN]/double");

  //Muon calorimeter type
  mLeptonData->Branch("MuonIsGlobal",                 m_MuonIsGlobal,               "MuonIsGlobal[MuonN]/bool");
  mLeptonData->Branch("MuonIsStandAlone",             m_MuonIsStandAlone,           "MuonIsStandAlone[MuonN]/bool");
  mLeptonData->Branch("MuonIsGlobalTight",            m_MuonIsGlobalTight,          "MuonIsGlobalTight[MuonN]/bool");
  mLeptonData->Branch("MuonIsTMLastStationLoose",     m_MuonIsTMLastStationLoose,   "MuonIsTMLastStationLoose[MuonN]/bool");
  mLeptonData->Branch("MuonIsTracker",                m_MuonIsTracker,              "MuonIsTracker[MuonN]/bool");
  mLeptonData->Branch("MuonIsTMLastStationTight",     m_MuonTMLastStationTight,     "MuonIsTMLastStationTight[MuonN]/bool");
  mLeptonData->Branch("MuonIsTM2DCompatibilityLoose", m_MuonTM2DCompatibilityLoose, "MuonIsTM2DCompatibilityLoose[MuonN]/bool");
  mLeptonData->Branch("MuonIsTM2DCompatibilityTight", m_MuonTM2DCompatibilityTight, "MuonIsTM2DCompatibilityTight[MuonN]/bool");

  //  mLeptonData->Branch("MuonId", m_MuonId, "MuonId[MuonN]/double");
  mLeptonData->Branch("MuonCombChi2", m_MuonCombChi2, "MuonCombChi2[MuonN]/double");
  mLeptonData->Branch("MuonCombNdof", m_MuonCombNdof, "MuonCombNdof[MuonN]/double");
  mLeptonData->Branch("MuonCombVx",   m_MuonCombVx,   "MuonCombVx[MuonN]/double");
  mLeptonData->Branch("MuonCombVy",   m_MuonCombVy,   "MuonCombVy[MuonN]/double");
  mLeptonData->Branch("MuonCombVz",   m_MuonCombVz,   "MuonCombVz[MuonN]/double");
  mLeptonData->Branch("MuonCombD0",   m_MuonCombD0,   "MuonCombD0[MuonN]/double");
  mLeptonData->Branch("MuonCombDz",   m_MuonCombDz,   "MuonCombDz[MuonN]/double");

  //Muon tracking information
  mLeptonData->Branch("MuonStandValidHits",   m_MuonStandValidHits,   "MuonStandValidHits[MuonN]/double");
  mLeptonData->Branch("MuonStandLostHits",    m_MuonStandLostHits,    "MuonStandLostHits[MuonN]/double");
  mLeptonData->Branch("MuonStandPt",          m_MuonStandPt,          "MuonStandPt[MuonN]/double");
  mLeptonData->Branch("MuonStandPz",          m_MuonStandPz,          "MuonStandPz[MuonN]/double");
  mLeptonData->Branch("MuonStandP",           m_MuonStandP,           "MuonStandP[MuonN]/double");
  mLeptonData->Branch("MuonStandEta",         m_MuonStandEta,         "MuonStandEta[MuonN]/double");
  mLeptonData->Branch("MuonStandPhi",         m_MuonStandPhi,         "MuonStandPhi[MuonN]/double");
  mLeptonData->Branch("MuonStandCharge",      m_MuonStandCharge,      "MuonStandCharge[MuonN]/double");
  mLeptonData->Branch("MuonStandChi",         m_MuonStandChi,         "MuonStandChi[MuonN]/double");
  mLeptonData->Branch("MuonStandQOverPError", m_MuonStandQOverPError, "MuonStandQOverPError[MuonN]/double");

  mLeptonData->Branch("MuonTrkValidHits",   m_MuonTrkValidHits,   "MuonTrkValidHits[MuonN]/double");
  mLeptonData->Branch("MuonTrkLostHits",    m_MuonTrkLostHits,    "MuonTrkLostHits[MuonN]/double");
  mLeptonData->Branch("MuonTrkD0",          m_MuonTrkD0,          "MuonTrkD0[MuonN]/double");
  mLeptonData->Branch("MuonTrkPt",          m_MuonTrkPt,          "MuonTrkPt[MuonN]/double");
  mLeptonData->Branch("MuonTrkPz",          m_MuonTrkPz,          "MuonTrkPz[MuonN]/double");
  mLeptonData->Branch("MuonTrkP",           m_MuonTrkP,           "MuonTrkP[MuonN]/double");
  mLeptonData->Branch("MuonTrkEta",         m_MuonTrkEta,         "MuonTrkEta[MuonN]/double");
  mLeptonData->Branch("MuonTrkPhi",         m_MuonTrkPhi,         "MuonTrkPhi[MuonN]/double");
  mLeptonData->Branch("MuonTrkCharge",      m_MuonTrkCharge,      "MuonTrkCharge[MuonN]/double");
  mLeptonData->Branch("MuonTrkChi",         m_MuonTrkChi,         "MuonTrkChi[MuonN]/double");
  mLeptonData->Branch("MuonTrkQOverPError", m_MuonTrkQOverPError, "MuonTrkQOverPError[MuonN]/double"); 
  mLeptonData->Branch("MuonTrkOuterZ",      m_MuonTrkOuterZ,      "MuonOuterZ[MuonN]/double");
  mLeptonData->Branch("MuonTrkOuterR",      m_MuonTrkOuterR,      "MuonOuterR[MuonN]/double");

  //Generator level muon information
  if (doMCData_) {
    mLeptonData->Branch("MuonGenPdgId",  m_GenMuonPdgId,  "MuonGenPdgId[MuonN]/double");
    mLeptonData->Branch("MuonGenMother", m_GenMuonMother, "MuonGenMother[MuonN]/double");
    mLeptonData->Branch("MuonGenPx",     m_GenMuonPx,     "MuonGenPx[MuonN]/double");
    mLeptonData->Branch("MuonGenPy",     m_GenMuonPy,     "MuonGenPy[MuonN]/double");
    mLeptonData->Branch("MuonGenPz",     m_GenMuonPz,     "MuonGenPz[MuonN]/double");
    mLeptonData->Branch("MuonGenPt",     m_GenMuonPt,     "MuonGenPt[MuonN]/double");
    mLeptonData->Branch("MuonGenEt",     m_GenMuonEt,     "MuonGenEt[MuonN]/double");
    mLeptonData->Branch("MuonGenE",      m_GenMuonE,      "MuonGenE[MuonN]/double");
    
    //generator leptons (electrons and muons)
    mLeptonData->Branch("genN",     &length,    "genN/int");
    mLeptonData->Branch("genid",     genIds,    "genIds[genN]/int");
    mLeptonData->Branch("genMother", genRefs,   "genRefs[genN]/int"); 
    mLeptonData->Branch("genE",      genE,      "genE[genN]/float");
    mLeptonData->Branch("genPx",     genPx,     "genPx[genN]/float");
    mLeptonData->Branch("genPy",     genPy,     "genPy[genN]/float");
    mLeptonData->Branch("genPz",     genPz,     "genPz[genN]/float");
    
    //generator leptons status (electrons and muons)
    mLeptonData->Branch("genLepN",     &genLepLength, "genLepN/int");
    mLeptonData->Branch("genLepId",     genLepIds,    "genLepIds[genLepN]/int");
    mLeptonData->Branch("genLepMother", genLepRefs,   "genLepRefs[genLepN]/int");
    mLeptonData->Branch("genLepStatus", genLepStatus, "genLepStatus[genLepN]/int");
    mLeptonData->Branch("genLepE",      genLepE,      "genLepE[genLepN]/float");
    mLeptonData->Branch("genLepPx",     genLepPx,     "genLepPx[genLepN]/float");
    mLeptonData->Branch("genLepPy",     genLepPy,     "genLepPy[genLepN]/float");
    mLeptonData->Branch("genLepPz",     genLepPz,     "genLepPz[genLepN]/float");

    ////generator photon status 1
    //mLeptonData->Branch("genPhotN",     &genPhotLength, "genPhotN/int");
    //mLeptonData->Branch("genPhotId",     genPhotIds,    "genPhotIds[genPhotN]/int");
    //mLeptonData->Branch("genPhotMother", genPhotRefs,   "genPhotRefs[genPhotN]/int");
    //mLeptonData->Branch("genPhotStatus", genPhotStatus, "genPhotStatus[genPhotN]/int");
    //mLeptonData->Branch("genPhotE",      genPhotE,      "genPhotE[genPhotN]/float");
    //mLeptonData->Branch("genPhotPx",     genPhotPx,     "genPhotPx[genPhotN]/float");
    //mLeptonData->Branch("genPhotPy",     genPhotPy,     "genPhotPy[genPhotN]/float");
    //mLeptonData->Branch("genPhotPz",     genPhotPz,     "genPhotPz[genPhotN]/float");

    //mLeptonData->Branch("AlpPtScale", &m_AlpPtScale,"AlpPtScale/double");
    //mLeptonData->Branch("AlpIdTest", &m_AlpIdTest ,"AlpIdTest/int");
    mLeptonData->Branch("pthat", &m_Pthat, "pthat/double");
  }    

  edm::LogInfo("LeptonEvent") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/PluginManager/interface/ModuleDef.h"
//#include "FWCore/Framework/interface/ModuleFactory.h"
//
//DEFINE_EDM_PLUGIN(LeptonAnalyzer);
