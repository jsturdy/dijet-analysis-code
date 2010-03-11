
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
// $Id: VertexAnalyzer.cpp,v 1.1 2010/01/29 16:10:31 sturdy Exp $
//
//


//#include "SusyAnalysis/EventSelector/interface/BJetEventSelector.h"
#include "JSturdy/DiJetAnalysis/interface/LeptonAnalyzer.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#define _USE_MATH_DEFINES
//#include <math.h>
#include <TMath.h>
#include <sstream>
using namespace std;
using namespace reco;
using namespace edm;

//________________________________________________________________________________________
LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& pset)
{ 

  leptonParams = pset;

  //defaults
  elecMaxEta_ = 12.;
  elecMinEt_  = 0.;
  muonMaxEta_ = 12.;
  muonMinEt_  = 0.;
  doMCData_ = true;

  // Read in parameters from the config file
  if (leptonParams.exists("elecMaxEta")) elecMaxEta_ = leptonParams.getParameter<double>("elecMaxEta");
  if (leptonParams.exists("elecMinEt"))  elecMinEt_  = leptonParams.getParameter<double>("elecMinEt");

  if (leptonParams.exists("muonMaxEta")) muonMaxEta_ = leptonParams.getParameter<double>("muonMaxEta");
  if (leptonParams.exists("muonMinEt"))  muonMinEt_  = leptonParams.getParameter<double>("muonMinEt");

  //if (leptonParams.exists("tauMaxEta")) tauMaxEta_ = leptonParams.getParameter<double>("tauMaxEta");
  //if (leptonParams.exists("tauMinEt"))  tauMinEt_  = leptonParams.getParameter<double>("tauMinEt");

  if (leptonParams.exists("doMCData"))   doMCData_   = leptonParams.getParameter<bool>("doMCData");
  if (doMCData_) if (leptonParams.exists("genTag"))    genTag_    = leptonParams.getParameter<edm::InputTag>("genTag");
 
  // get the data tags
  elecTag_   = leptonParams.getParameter<edm::InputTag>("elecTag");
  pfelecTag_ = leptonParams.getParameter<edm::InputTag>("pfelecTag");

  muonTag_   = leptonParams.getParameter<edm::InputTag>("muonTag");
  pfmuonTag_ = leptonParams.getParameter<edm::InputTag>("pfmuonTag");

  //tauTag_   = leptonParams.getParameter<edm::InputTag>("tauTag");
  //pftauTag_ = leptonParams.getParameter<edm::InputTag>("pftauTag");

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
  using namespace edm;

  //bool preselection = false;
  bool lepton_result = true;
  edm::LogVerbatim("LeptonEvent") << " Start  " << std::endl;

  std::ostringstream dbg;

  //get pthat of process
  m_Pthat = -999.;
  
  Handle<double> genEventScale;
  iEvent.getByLabel( "genEventScale", genEventScale );
  if ( genEventScale.isValid() ) m_Pthat = *genEventScale;
  
  // GEN INFO do only if running on MC data
  if(doMCData_) {
    //ALPGENParticleId myALPGENParticleId;
    //
    //mTempAlpIdTest = myALPGENParticleId.AplGenParID(iEvent,genTag_);
    //mTempAlpPtScale = myALPGENParticleId.getPt();

    Handle<reco::GenParticleCollection>  genParticles;
    iEvent.getByLabel(genTag_, genParticles);   

    int count=0; int lcount=0;// int tcount=0;
    //length=genParticles->size();
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
	} else { genRefs[count]=-1;}

	count++;
      } else { // store also electrons or muons of status 1 
	if ( (abs(pCand.pdgId()) == 11) || (abs(pCand.pdgId()) == 13)  ) {
	  //      if ( (abs(pCand.pdgId()) == 11) || (abs(pCand.pdgId()) == 13) || pCand.pdgId() == 22 ) {

	  genLepIds[lcount]    = pCand.pdgId();
	  genLepStatus[lcount] = pCand.status();
	  genLepE[lcount]      = pCand.energy();
	  genLepPx[lcount]     = pCand.px();
	  genLepPy[lcount]     = pCand.py();
	  genLepPz[lcount]     = pCand.pz();
	
	  if (pCand.numberOfMothers() > 0 ) { 
	    const reco::Candidate * mom = pCand.mother();
	    while (mom->pdgId() == pCand.pdgId()) { mom = mom->mother(); }
	    //if (mom->pdgId() == pCand.pdgId()) { mom = mom->mother(); }

	    for( size_t j = 0; j < i; ++ j ) {
	      const Candidate * ref = &((*genParticles)[j]);
	      if (ref == mom) { genLepRefs[lcount] = ref->pdgId(); }
	      //if (ref == mom) { genLepRefs[lcount] = j; }
	    }  
	  } else { genLepRefs[lcount]=-1;}
	  lcount++;
	}
      }
    }
    length = count;
    genLepLength = lcount;
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

  
  edm::LogVerbatim("LeptonEvent") << " start reading in electrons " << endl;
  // Add the electrons
  m_Nelec= elecHandle->size();
  if ( m_Nelec > 50 ) m_Nelec = 50;
  for (int i=0;i<m_Nelec;i++){
    edm::LogVerbatim("LeptonEvent") << " looping over electrons " << endl;
    m_ElecE[i]      = (*elecHandle)[i].energy();
    m_ElecEt[i]     = (*elecHandle)[i].et();
    m_ElecPt[i]     = (*elecHandle)[i].pt();
    m_ElecPx[i]     = (*elecHandle)[i].momentum().X();
    m_ElecPy[i]     = (*elecHandle)[i].momentum().Y();
    m_ElecPz[i]     = (*elecHandle)[i].momentum().Z();
    m_ElecEta[i]    = (*elecHandle)[i].eta();
    m_ElecPhi[i]    = (*elecHandle)[i].phi();
    m_ElecCharge[i] = (*elecHandle)[i].charge();

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

    edm::LogVerbatim("LeptonEvent") << " before asking for gsfTrack " << endl;
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
    
    edm::LogVerbatim("LeptonEvent") << " before asking for trackMomentumAtVtx " << endl;
    
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

    if (&(*(*elecHandle)[i].genLepton())!=0){
      m_GenElecPdgId[i] = (*elecHandle)[i].genLepton()->pdgId();
      m_GenElecPx[i]    = (*elecHandle)[i].genLepton()->px();
      m_GenElecPy[i]    = (*elecHandle)[i].genLepton()->py();
      m_GenElecPz[i]    = (*elecHandle)[i].genLepton()->pz();
      m_GenElecPt[i]    = (*elecHandle)[i].genLepton()->pt();
      m_GenElecEt[i]    = (*elecHandle)[i].genLepton()->et();
      m_GenElecE[i]     = (*elecHandle)[i].genLepton()->energy();

      if(&(*(*elecHandle)[i].genLepton()->mother())!=0){
	m_GenElecMother[i] = (*elecHandle)[i].genLepton()->mother()->pdgId();
	if ( (*elecHandle)[i].genLepton()->mother()->pdgId() ==  (*elecHandle)[i].genLepton()->pdgId()) 
	  {
	    m_GenElecMother[i] = (*elecHandle)[i].genLepton()->mother()->mother()->pdgId();
	  }
      }
    }
    else {
      m_GenElecPdgId[i]  = 999.;
      m_GenElecMother[i] = 999.;
      m_GenElecPx[i]     = 999.;
      m_GenElecPy[i]     = 999.;
      m_GenElecPz[i]     = 999.;
      m_GenElecPt[i]     = 999.;
      m_GenElecEt[i]     = 999.;
      m_GenElecE[i]      = 999.;
    }

  }//end loop over Electrons


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
  
  edm::LogVerbatim("LeptonEvent") << " start reading in muons " << endl;


  // Add the muons
  m_Nmuon= muonHandle->size();
  if ( m_Nmuon > 50 ) m_Nmuon = 50;
  for (int i=0;i<m_Nmuon;i++){
   
    m_MuonPt[i]  = (*muonHandle)[i].pt();
    m_MuonE[i]   = (*muonHandle)[i].energy();
    m_MuonEt[i]  = (*muonHandle)[i].et();
    m_MuonPx[i]  = (*muonHandle)[i].momentum().X();
    m_MuonPy[i]  = (*muonHandle)[i].momentum().Y();
    m_MuonPz[i]  = (*muonHandle)[i].momentum().Z();
    m_MuonEta[i] = (*muonHandle)[i].eta();
    m_MuonPhi[i] = (*muonHandle)[i].phi();
    
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
    if (&(*(*muonHandle)[i].genLepton())!=0){
      m_GenMuonPdgId[i] = (*muonHandle)[i].genLepton()->pdgId();
      m_GenMuonPx[i]    = (*muonHandle)[i].genLepton()->px();
      m_GenMuonPy[i]    = (*muonHandle)[i].genLepton()->py();
      m_GenMuonPz[i]    = (*muonHandle)[i].genLepton()->pz();
      m_GenMuonPt[i]    = (*muonHandle)[i].genLepton()->pt();
      m_GenMuonEt[i]    = (*muonHandle)[i].genLepton()->et();
      m_GenMuonE[i]     = (*muonHandle)[i].genLepton()->energy();

      /*
       *figure out how best to make sure we get to the mother
       */
//      reco::Candidate *muMother = (*muonHandle)[i].genLepton()->mother();
//      if (muMother!=0) {
//	m_GenMuonMother[i] = muMother->pdgId();
//	if ( muMother->pdgId() ==  (*muonHandle)[i].genLepton()->pdgId()) 
//	  {
//	    m_GenMuonMother[i] = muMother->mother()->pdgId();
//	  }
//	//while ( muMother->pdgId() ==  (*muonHandle)[i].genLepton()->pdgId()) {muMother = muMother->mother(); }
//      }
      if (&(*(*muonHandle)[i].genLepton()->mother())!=0) {
	m_GenMuonMother[i]=(*muonHandle)[i].genLepton()->mother()->pdgId();
	if ( (*muonHandle)[i].genLepton()->mother()->pdgId() ==  (*muonHandle)[i].genLepton()->pdgId()) 
	  {
	    m_GenMuonMother[i]= (*muonHandle)[i].genLepton()->mother()->mother()->pdgId();
	  }
      }	else {
	m_GenMuonMother[i] = 999.;
      }
    }
    else{
      m_GenMuonPdgId[i]  = 999.;
      m_GenMuonMother[i] = 999.;
      m_GenMuonPx[i]     = 999.;
      m_GenMuonPy[i]     = 999.;
      m_GenMuonPz[i]     = 999.;
      m_GenMuonPt[i]     = 999.;
      m_GenMuonEt[i]     = 999.;
      m_GenMuonE[i]      = 999.;
    }

  }


  // Fill the tree only if all preselection conditions are met
  //if(preselection) //mPreselection->Fill();
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
  
  // Register this ntuple
  edm::Service<TFileService> fs;

  // Now we add some additional ones for the dijet analysis
  //mPreselection = fs->make<TTree>( "Preselection", "data after preselection" );
  mLeptonData = fs->make<TTree>( "LeptonData", "data after cuts" );
  mLeptonData->SetAutoSave(10);

  // Add the branches

  //add electrons
  //General electron information
  mLeptonData->Branch("Nelec",     &m_Nelec,      "Nelec/int");  
  mLeptonData->Branch("ElecE",      m_ElecE,      "ElecE[Nelec]/double");
  mLeptonData->Branch("ElecEt",     m_ElecEt,     "ElecEt[Nelec]/double");
  mLeptonData->Branch("ElecPt",     m_ElecPt,     "ElecPt[Nelec]/double");
  mLeptonData->Branch("ElecPx",     m_ElecPx,     "ElecPx[Nelec]/double");
  mLeptonData->Branch("ElecPy",     m_ElecPy,     "ElecPy[Nelec]/double");
  mLeptonData->Branch("ElecPz",     m_ElecPz,     "ElecPz[Nelec]/double");
  mLeptonData->Branch("ElecEta",    m_ElecEta,    "ElecEta[Nelec]/double");
  mLeptonData->Branch("ElecPhi",    m_ElecPhi,    "ElecPhi[Nelec]/double");
  mLeptonData->Branch("ElecCharge", m_ElecCharge, "ElecCharge[Nelec]/double");
  mLeptonData->Branch("ElecHOverE", m_ElecHOverE, "ElecHOverE[Nelec]/double");

  //Isolation and tracking variables
  mLeptonData->Branch("ElecTrkIso",     m_ElecTrkIso,   "ElecTrkIso[Nelec]/double");
  mLeptonData->Branch("ElecECalIso",    m_ElecECalIso,  "ElecECalIso[Nelec]/double");
  mLeptonData->Branch("ElecHCalIso",    m_ElecHCalIso,  "ElecHCalIso[Nelec]/double");
  mLeptonData->Branch("ElecAllIso",     m_ElecAllIso,   "ElecAllIso[Nelec]/double");
  mLeptonData->Branch("ElecTrkChiNorm", m_ElecNormChi2, "ElecTrkChiNorm[Nelec]/double");
  //mLeptonData->Branch("NIsoelec",      &m_NIsoelec,     "NIsoelec/int");  

  mLeptonData->Branch("ElecECalIsoDeposit", m_ElecECalIsoDeposit, "ElecECalIsoDeposit[Nelec]/double");
  mLeptonData->Branch("ElecHCalIsoDeposit", m_ElecHCalIsoDeposit, "ElecHCalIsoDeposit[Nelec]/double");

  //Electron identification values
  mLeptonData->Branch("ElecIdLoose",    m_ElecIdLoose,    "ElecIdLoose [Nelec]/double");
  mLeptonData->Branch("ElecIdTight",    m_ElecIdTight,    "ElecIdTight [Nelec]/double");
  mLeptonData->Branch("ElecIdRobLoose", m_ElecIdRobLoose, "ElecIdRobLoose [Nelec]/double");
  mLeptonData->Branch("ElecIdRobTight", m_ElecIdRobTight, "ElecIdRobTight [Nelec]/double");
  mLeptonData->Branch("ElecChargeMode", m_ElecChargeMode, "ElecChargeMode [Nelec]/double");
  mLeptonData->Branch("ElecPtMode",     m_ElecPtTrkMode,  "ElecPtMode [Nelec]/double");


  //Electron vertex information
  mLeptonData->Branch("ElecVx",     m_ElecVx,     "ElecVx[Nelec]/double");
  mLeptonData->Branch("ElecVy",     m_ElecVy,     "ElecVy[Nelec]/double");
  mLeptonData->Branch("ElecVz",     m_ElecVz,     "ElecVz[Nelec]/double");
  mLeptonData->Branch("ElecD0",     m_ElecD0,     "ElecD0[Nelec]/double");
  mLeptonData->Branch("ElecDz",     m_ElecDz,     "ElecDz[Nelec]/double");
  mLeptonData->Branch("ElecPtTrk",  m_ElecPtTrk,  "ElecPtTrk[Nelec]/double");

  //Additonal electron detector information
  //Electron tracking information
  mLeptonData->Branch("ElecQOverPErrTrkMode", m_ElecQOverPErrTrkMode, "ElecQOverPErrTrkMode [Nelec]/double");
  mLeptonData->Branch("ElecCaloEnergy",       m_ElecCaloEnergy,       "ElecCaloEnergy[Nelec]/double");
  mLeptonData->Branch("ElecQOverPErrTrk",     m_ElecQOverPErrTrk,     "ElecQOverPErrTrk[Nelec]/double");
  mLeptonData->Branch("ElecPinTrk",           m_ElecPinTrk,           "ElecPinTrk[Nelec]/double");
  mLeptonData->Branch("ElecPoutTrk",          m_ElecPoutTrk,          "ElecPoutTrk[Nelec]/double"); 
  mLeptonData->Branch("ElecLostHits",         m_ElecLostHits,         "ElecLostHits[Nelec]/double"); 
  mLeptonData->Branch("ElecValidHits",        m_ElecValidHits,        "ElecValidHits[Nelec]/double"); 
  mLeptonData->Branch("ElecNCluster",         m_ElecNCluster,         "ElecNCluster[Nelec]/double"); 
  mLeptonData->Branch("ElecEtaTrk",           m_ElecEtaTrk,           "ElecEtaTrk[Nelec]/double"); 
  mLeptonData->Branch("ElecPhiTrk",           m_ElecPhiTrk,           "ElecPhiTrk[Nelec]/double"); 
  mLeptonData->Branch("ElecWidthClusterEta",  m_ElecWidthClusterEta,  "ElecWidthClusterEta[Nelec]/double"); 
  mLeptonData->Branch("ElecWidthClusterPhi",  m_ElecWidthClusterPhi,  "ElecWidthClusterPhi[Nelec]/double"); 

  //Generator level information stored in the electron object
  mLeptonData->Branch("ElecGenPdgId",  m_GenElecPdgId,  "ElecGenPdgId[Nelec]/double");
  mLeptonData->Branch("ElecGenMother", m_GenElecMother, "ElecGenMother[Nelec]/double");
  mLeptonData->Branch("ElecGenPx",     m_GenElecPx,     "ElecGenPx[Nelec]/double");
  mLeptonData->Branch("ElecGenPy",     m_GenElecPy,     "ElecGenPy[Nelec]/double");
  mLeptonData->Branch("ElecGenPz",     m_GenElecPz,     "ElecGenPz[Nelec]/double");

  //add muons
  //General kinematic variables related to muons
  mLeptonData->Branch("Nmuon",         &m_Nmuon,          "Nmuon/int");  
  mLeptonData->Branch("MuonE",          m_MuonE,          "MuonE[Nmuon]/double");
  mLeptonData->Branch("MuonEt",         m_MuonEt,         "MuonEt[Nmuon]/double");
  mLeptonData->Branch("MuonPt",         m_MuonPt,         "MuonPt[Nmuon]/double");
  mLeptonData->Branch("MuonPx",         m_MuonPx,         "MuonPx[Nmuon]/double");
  mLeptonData->Branch("MuonPy",         m_MuonPy,         "MuonPy[Nmuon]/double");
  mLeptonData->Branch("MuonPz",         m_MuonPz,         "MuonPz[Nmuon]/double");
  mLeptonData->Branch("MuonEta",        m_MuonEta,        "MuonEta[Nmuon]/double");
  mLeptonData->Branch("MuonPhi",        m_MuonPhi,        "MuonPhi[Nmuon]/double");
  mLeptonData->Branch("MuonCharge",     m_MuonCharge,     "MuonCharge[Nmuon]/double");

  //Muon isolation variables
  //mLeptonData->Branch("NIsomuon",      &m_NIsomuon,       "NIsomuon/int");  
  mLeptonData->Branch("MuonTrkIso",     m_MuonTrkIso,     "MuonTrkIso[Nmuon]/double");
  mLeptonData->Branch("MuonECalIso",    m_MuonECalIso,    "MuonECalIso[Nmuon]/double");
  mLeptonData->Branch("MuonHCalIso",    m_MuonHCalIso,    "MuonHCalIso[Nmuon]/double");
  mLeptonData->Branch("MuonAllIso",     m_MuonAllIso,     "MuonAllIso[Nmuon]/double");
  mLeptonData->Branch("MuonTrkChiNorm", m_MuonTrkChiNorm, "MuonTrkChiNorm[Nmuon]/double");

  mLeptonData->Branch("MuonECalIsoDeposit", m_MuonECalIsoDeposit,"MuonECalIsoDeposit[Nmuon]/double");
  mLeptonData->Branch("MuonHCalIsoDeposit", m_MuonHCalIsoDeposit,"MuonHCalIsoDeposit[Nmuon]/double");

  //Muon calorimeter type
  mLeptonData->Branch("MuonIsGlobal",                m_MuonIsGlobal,              "m_MuonIsGlobal[Nmuon]/bool");
  mLeptonData->Branch("MuonIsStandAlone",            m_MuonIsStandAlone,          "m_MuonIsStandAlone[Nmuon]/bool");
  mLeptonData->Branch("MuonIsGlobalTight",           m_MuonIsGlobalTight,         "m_MuonIsGlobalTight[Nmuon]/bool");
  mLeptonData->Branch("MuonIsTMLastStationLoose",    m_MuonIsTMLastStationLoose,  "m_MuonIsTMLastStationLoose[Nmuon]/bool");
  mLeptonData->Branch("MuonIsTracker",               m_MuonIsTracker,             "m_MuonIsTracker[Nmuon]/bool");
  mLeptonData->Branch("MuonIsTMLastStationTight",    m_MuonTMLastStationTight,    "MuonIsTMLastStationTight[Nmuon]/bool");
  mLeptonData->Branch("MuonIsTM2DCompatibilityLoose",m_MuonTM2DCompatibilityLoose,"MuonIsTM2DCompatibilityLoose[Nmuon]/bool");
  mLeptonData->Branch("MuonIsTM2DCompatibilityTight",m_MuonTM2DCompatibilityTight,"MuonIsTM2DCompatibilityTight[Nmuon]/bool");

  //  mLeptonData->Branch("MuonId",m_MuonId,"m_MuonId[Nmuon]/double");
  mLeptonData->Branch("MuonCombChi2",m_MuonCombChi2,"m_MuonCombChi2[Nmuon]/double");
  mLeptonData->Branch("MuonCombNdof",m_MuonCombNdof,"m_MuonCombNdof[Nmuon]/double");
  mLeptonData->Branch("MuonCombVx",  m_MuonCombVx,  "m_MuonCombVx[Nmuon]/double");
  mLeptonData->Branch("MuonCombVy",  m_MuonCombVy,  "m_MuonCombVy[Nmuon]/double");
  mLeptonData->Branch("MuonCombVz",  m_MuonCombVz,  "m_MuonCombVz[Nmuon]/double");
  mLeptonData->Branch("MuonCombD0",  m_MuonCombD0,  "m_MuonCombD0[Nmuon]/double");
  mLeptonData->Branch("MuonCombDz",  m_MuonCombDz,  "m_MuonCombDz[Nmuon]/double");

  //Muon tracking information
  mLeptonData->Branch("MuonStandValidHits",  m_MuonStandValidHits,  "m_MuonStandValidHits[Nmuon]/double");
  mLeptonData->Branch("MuonStandLostHits",   m_MuonStandLostHits,   "m_MuonStandLostHits[Nmuon]/double");
  mLeptonData->Branch("MuonStandPt",         m_MuonStandPt,         "m_MuonStandPt[Nmuon]/double");
  mLeptonData->Branch("MuonStandPz",         m_MuonStandPz,         "m_MuonStandPz[Nmuon]/double");
  mLeptonData->Branch("MuonStandP",          m_MuonStandP,          "m_MuonStandP[Nmuon]/double");
  mLeptonData->Branch("MuonStandEta",        m_MuonStandEta,        "m_MuonStandEta[Nmuon]/double");
  mLeptonData->Branch("MuonStandPhi",        m_MuonStandPhi,        "m_MuonStandPhi[Nmuon]/double");
  mLeptonData->Branch("MuonStandCharge",     m_MuonStandCharge,     "m_MuonStandCharge[Nmuon]/double");
  mLeptonData->Branch("MuonStandChi",        m_MuonStandChi,        "m_MuonStandChi[Nmuon]/double");
  mLeptonData->Branch("MuonStandQOverPError",m_MuonStandQOverPError,"m_MuonStandQOverPError[Nmuon]/double");

  mLeptonData->Branch("MuonTrkValidHits",   m_MuonTrkValidHits,   "m_MuonTrkValidHits[Nmuon]/double");
  mLeptonData->Branch("MuonTrkLostHits",    m_MuonTrkLostHits,    "m_MuonTrkLostHits[Nmuon]/double");
  mLeptonData->Branch("MuonTrkD0",          m_MuonTrkD0,          "m_MuonTrkD0[Nmuon]/double");
  mLeptonData->Branch("MuonTrkPt",          m_MuonTrkPt,          "m_MuonTrkPt[Nmuon]/double");
  mLeptonData->Branch("MuonTrkPz",          m_MuonTrkPz,          "m_MuonTrkPz[Nmuon]/double");
  mLeptonData->Branch("MuonTrkP",           m_MuonTrkP,           "m_MuonTrkP[Nmuon]/double");
  mLeptonData->Branch("MuonTrkEta",         m_MuonTrkEta,         "m_MuonTrkEta[Nmuon]/double");
  mLeptonData->Branch("MuonTrkPhi",         m_MuonTrkPhi,         "m_MuonTrkPhi[Nmuon]/double");
  mLeptonData->Branch("MuonTrkCharge",      m_MuonTrkCharge,      "m_MuonTrkCharge[Nmuon]/double");
  mLeptonData->Branch("MuonTrkChi",         m_MuonTrkChi,         "m_MuonTrkChi[Nmuon]/double");
  mLeptonData->Branch("MuonTrkQOverPError", m_MuonTrkQOverPError, "m_MuonTrkQOverPError[Nmuon]/double"); 
  mLeptonData->Branch("MuonTrkOuterZ",      m_MuonTrkOuterZ,      "m_MuonOuterZ[Nmuon]/double");
  mLeptonData->Branch("MuonTrkOuterR",      m_MuonTrkOuterR,      "m_MuonOuterR[Nmuon]/double");

  //Generator level muon information
  mLeptonData->Branch("MuonGenPdgId",  m_GenMuonPdgId,  "MuonGenPdgId[Nmuon]/double");
  mLeptonData->Branch("MuonGenMother", m_GenMuonMother, "MuonGenMother[Nmuon]/double");
  mLeptonData->Branch("MuonGenPx",     m_GenMuonPx,     "MuonGenPx[Nmuon]/double");
  mLeptonData->Branch("MuonGenPy",     m_GenMuonPy,     "MuonGenPy[Nmuon]/double");
  mLeptonData->Branch("MuonGenPz",     m_GenMuonPz,     "MuonGenPz[Nmuon]/double");

  //generator leptons (electrons and muons)
  mLeptonData->Branch("genN",     &length,    "genN/int");
  mLeptonData->Branch("genid",     genIds,    "genIds[genN]/int");
  mLeptonData->Branch("genMother", genRefs,   "genRefs[genN]/int");
  mLeptonData->Branch("genE",      genE,      "genE[genN]/float");
  mLeptonData->Branch("genPx",     genPx,     "genPx[genN]/float");
  mLeptonData->Branch("genPy",     genPy,     "genPy[genN]/float");
  mLeptonData->Branch("genPz",     genPz,     "genPz[genN]/float");
  mLeptonData->Branch("genStatus", genStatus, "genStatus[genN]/int");

  //generator leptons status (electrons and muons)
  mLeptonData->Branch("genLepN",     &genLepLength, "genLepN/int");
  mLeptonData->Branch("genLepId",     genLepIds,    "genLepIds[genLepN]/int");
  mLeptonData->Branch("genLepMother", genLepRefs,   "genLepRefs[genLepN]/int");
  mLeptonData->Branch("genLepE",      genLepE,      "genLepE[genLepN]/float");
  mLeptonData->Branch("genLepPx",     genLepPx,     "genLepPx[genLepN]/float");
  mLeptonData->Branch("genLepPy",     genLepPy,     "genLepPy[genLepN]/float");
  mLeptonData->Branch("genLepPz",     genLepPz,     "genLepPz[genLepN]/float");
  mLeptonData->Branch("genLepStatus", genLepStatus, "genLepStatus[genLepN]/int");

  //mLeptonData->Branch("AlpPtScale", &mTempAlpPtScale,"mTempAlpPtScale/double");
  //mLeptonData->Branch("AlpIdTest", &mTempAlpIdTest ,"AlpIdTest/int");
  mLeptonData->Branch("pthat",&m_Pthat,"pthat/double");
    
  edm::LogInfo("LeptonEvent") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/PluginManager/interface/ModuleDef.h"
//#include "FWCore/Framework/interface/ModuleFactory.h"
//
//DEFINE_EDM_PLUGIN(LeptonAnalyzer);
