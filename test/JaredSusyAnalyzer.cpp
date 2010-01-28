
// -*- C++ -*-
//
// Package:    JaredSusyAnalyzer
// Class:      JaredSusyAnalyzer
// 
/**\class JaredSusyAnalyzer JaredSusyAnalyzer.cc JSturdy/AnalysisTools/plugins/JaredSusyAnalyzer.cc

Description: Variable collector/ntupler for SUSY search with Jets + MET

Implementation:Uses the EventSelector interface for event selection and TFileService for plotting.

*/
//
// Original Author:  Markus Stoye, (modified by Jared Sturdy from SusyDiJetAnalysis)
//         Created:  Mon Feb 18 15:40:44 CET 2008
// $Id: JaredSusyAnalyzer.cpp,v 1.39 2009/05/31 12:46:17 pioppi Exp $
//
//
//#include "SusyAnalysis/EventSelector/interface/BJetEventSelector.h"
#include "UserCode/AnalysisTools/test/JaredSusyAnalyzer.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#define _USE_MATH_DEFINES
//#include <math.h>
#include <TMath.h>
#include <sstream>
using namespace std;
using namespace reco;
using namespace edm;

//________________________________________________________________________________________
JaredSusyAnalyzer::JaredSusyAnalyzer(const edm::ParameterSet& iConfig):
  eventWeight_( iConfig.getParameter<double>("eventWeight") ),
  nrEventTotalRaw_(0), nrEventTotalWeighted_(0.0),
  pathNames_(0), nEvents_(0), nWasRun_(0), nAccept_(0), nErrors_(0), hlWasRun_(0), hlAccept_(0), hlErrors_(0), init_(false) //georgia
{ 

  //defaults
  jetMaxEta_ = 3.0;
  jetMinPt_ = 30;
  doMCData_ = true;
  // Say something about event weights
  if (iConfig.exists("jetMaxEta")) jetMaxEta_ = iConfig.getParameter<double>("jetMaxEta");
  if (iConfig.exists("jetMinPt"))  jetMinPt_  = iConfig.getParameter<double>("jetMinPt");
  if (iConfig.exists("doMCData"))  doMCData_  = iConfig.getParameter<bool>("doMCData");
  if (doMCData_) if (iConfig.exists("genTag"))    genTag_    = iConfig.getParameter<edm::InputTag>("genTag");
 
  edm::LogInfo("JaredSusyTest") << "Global event weight set to " << eventWeight_;
    
  // get the data tags
  elecTag_   = iConfig.getParameter<edm::InputTag>("elecTag");
  muonTag_   = iConfig.getParameter<edm::InputTag>("muonTag");
  pfelecTag_ = iConfig.getParameter<edm::InputTag>("pfelecTag");
  pfmuonTag_ = iConfig.getParameter<edm::InputTag>("pfmuonTag");
  vtxTag_    = iConfig.getParameter<edm::InputTag>("vtxTag"); 

  //MET
  tcmetTag_ = iConfig.getParameter<edm::InputTag>("tcmetTag");
  pfmetTag_ = iConfig.getParameter<edm::InputTag>("pfmetTag");
  metTag_   = iConfig.getParameter<edm::InputTag>("metTag");
  //mhtTag_   = iConfig.getParameter<edm::InputTag>("mhtTag");

  //Jets
  usePfjets_ = iConfig.getParameter<bool>("UsePfjet");
  pfjetTag_  = iConfig.getParameter<edm::InputTag>("pfjetTag");
  jetTag_    = iConfig.getParameter<edm::InputTag>("jetTag");
  jptTag_    = iConfig.getParameter<edm::InputTag>("jptTag");

  // trigger stuff
  triggerResults_ = iConfig.getParameter<edm::InputTag>("triggerResults");
  // trigger path names
  pathNames_ = iConfig.getParameter< std::vector<std::string> >("pathNames");

  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
  initPlots();
}


//________________________________________________________________________________________
JaredSusyAnalyzer::~JaredSusyAnalyzer() {}


//________________________________________________________________________________________
// Method called to for each event
void
JaredSusyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  //bool preselection = false;
  edm::LogVerbatim("JaredSusyEvent") << " Start  " << std::endl;

  std::ostringstream dbg;

  run_   = iEvent.id().run();
  event_ = iEvent.id().event();

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
      const reco::Candidate& p = (*genParticles)[ i ];

      int st = p.status();  

      if (st==3) {
	ids[count] = p.pdgId(); genStatus[count]=p.status();
	genE[count]=p.energy(); genPx[count]=p.px(); genPy[count]=p.py(); genPz[count]=p.pz();
      
	if (p.numberOfMothers() > 0 ) { 
	  const reco::Candidate * mom = p.mother();

	  for( size_t j = 0; j < i; ++ j ) {
	    const Candidate * ref = &((*genParticles)[j]);
	    if (ref == mom) { refs[count]= j; }
	  }  
	} else { refs[count]=-1;}

	count++;
      } else { // store also electrons or muons or photons of status 1 
	if ( (abs(p.pdgId()) == 11) || (abs(p.pdgId()) == 13)  ) {
	  //      if ( (abs(p.pdgId()) == 11) || (abs(p.pdgId()) == 13) || p.pdgId() == 22 ) {

	  genLepIds[lcount] = p.pdgId(); genLepStatus[lcount]=p.status();
	  genLepE[lcount]=p.energy(); genLepPx[lcount]=p.px(); genLepPy[lcount]=p.py(); genLepPz[lcount]=p.pz();
	
	  if (p.numberOfMothers() > 0 ) { 
	    const reco::Candidate * mom = p.mother();
	    if (mom->pdgId() == p.pdgId()) { mom = mom->mother(); }

	    for( size_t j = 0; j < i; ++ j ) {
	      const Candidate * ref = &((*genParticles)[j]);
	      if (ref == mom) { genLepRefs[lcount]=ref->pdgId(); }
	    }  
	  } else { genLepRefs[lcount]=-1;}
	  lcount++;
	}
      }
    }
    length=count; genLepLength=lcount;
  }  

  //Trigger results
  edm::LogVerbatim("JaredSusyEvent") << " Trigger decision  " << std::endl;

  //get the trigger decision
  mTempTreeHLT1JET=false;
  mTempTreeHLT2JET=false;
  mTempTreeHLT1MET1HT=false;

  // Get the trigger results and check validity
  edm::Handle<edm::TriggerResults> hltHandle;
  iEvent.getByLabel(triggerResults_, hltHandle);
  if ( !hltHandle.isValid() ) {
    edm::LogWarning("HLTEventSelector") << "No trigger results for InputTag " << triggerResults_;
    return;
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
  
  nHLT=static_cast<int>(n);
  for(unsigned int i=0; i!=n; ++i) {
    HLTArray[i]=hltHandle->accept(i);
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
      if (trigName == "HLT_Jet180")      mTempTreeHLT1JET=true;
      if (trigName == "HLT_DiJetAve130") mTempTreeHLT2JET=true;
      if (trigName == "HLT_MET50")       mTempTreeHLT1MET1HT=true;
      if (trigName == "HLT_Mu9")         mTempTreeHLT1Muon=true; 
      
    } 
  }

 
  mTempTreeEventWeight =eventWeight_;

  mTempTreeRun = run_;
  mTempTreeEvent = event_;

  // get the Vertex collection

  LogDebug("JaredSusyEvent") << "Vertex results for InputTag" << vtxTag_;
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vtxTag_, vertices);
  if ( !vertices.isValid() ) {
    LogDebug("JaredSusyEvent") << "No Vertex results for InputTag" << vtxTag_;
    return;
  } 

  mTempTreenVtx = (*vertices).size();

  for (int i=0; i< (int) (*vertices).size(); i++){  
    //    int indPrim=0;
    const reco::Vertex* pVertex = &(*vertices)[i];
    mTempTreeVtxNormalizedChi2[i] = pVertex->normalizedChi2();
    mTempTreeVtxChi2[i] = pVertex->chi2();
    mTempTreeVtxNdof[i] = pVertex->ndof();
    mTempTreeVtxX[i]    = pVertex->x();
    mTempTreeVtxY[i]    = pVertex->y();
    mTempTreeVtxZ[i]    = pVertex->z();
    mTempTreeVtxdX[i]   = pVertex->xError();
    mTempTreeVtxdY[i]   = pVertex->yError();
    mTempTreeVtxdZ[i]   = pVertex->zError();
  } 

  // get the electrons
  edm::Handle< std::vector<pat::Electron> > elecHandle;
  iEvent.getByLabel(elecTag_, elecHandle);
  if ( !elecHandle.isValid() ) {
    edm::LogWarning("JaredSusyEvent") << "No Electron results for InputTag " << elecTag_;
    return;
  }

  
  edm::LogVerbatim("JaredSusyEvent") << " start reading in electrons " << endl;
  // Add the electrons
  mTempTreeNelec= elecHandle->size();
  if ( mTempTreeNelec > 50 ) mTempTreeNelec = 50;
  for (int i=0;i<mTempTreeNelec;i++){
    mTempTreeElecE[i]      = (*elecHandle)[i].energy();
    mTempTreeElecEt[i]     = (*elecHandle)[i].et();
    mTempTreeElecPt[i]     = (*elecHandle)[i].pt();
    mTempTreeElecPx[i]     = (*elecHandle)[i].momentum().X();
    mTempTreeElecPy[i]     = (*elecHandle)[i].momentum().Y();
    mTempTreeElecPz[i]     = (*elecHandle)[i].momentum().Z();
    mTempTreeElecEta[i]    = (*elecHandle)[i].eta();
    mTempTreeElecPhi[i]    = (*elecHandle)[i].phi();
    mTempTreeElecTrkIso[i] = (*elecHandle)[i].trackIso();

    mTempTreeElecECalIso[i] = (*elecHandle)[i].ecalIso();
    mTempTreeElecHCalIso[i] = (*elecHandle)[i].hcalIso() ;
    mTempTreeElecAllIso[i]  = (*elecHandle)[i].caloIso() ;
    mTempTreeElecCharge[i]  = (*elecHandle)[i].charge();

    //mTempTreeElecECalIsoDeposit[i]  = (*elecHandle)[i].ecalIsoDeposit()->candEnergy() ;
    //mTempTreeElecHCalIsoDeposit[i]  = (*elecHandle)[i].hcalIsoDeposit()->candEnergy() ;

    mTempTreeElecIdLoose[i]    = (*elecHandle)[i].electronID("eidLoose");
    edm::LogVerbatim("JaredSusyEvent") << " looping over electrons " << endl;
    mTempTreeElecIdTight[i]    = (*elecHandle)[i].electronID("eidTight");
    mTempTreeElecIdRobLoose[i] = (*elecHandle)[i].electronID("eidRobustLoose");
    mTempTreeElecIdRobTight[i] = (*elecHandle)[i].electronID("eidRobustTight"); 

    mTempTreeElecCaloEnergy[i] = (*elecHandle)[i].caloEnergy();
    mTempTreeElecHOverE[i]     = (*elecHandle)[i].hadronicOverEm();
    mTempTreeElecVx[i]         = (*elecHandle)[i].vx();
    mTempTreeElecVy[i]         = (*elecHandle)[i].vy();
    mTempTreeElecVz[i]         = (*elecHandle)[i].vz();

    edm::LogVerbatim("JaredSusyEvent") << " before asking for gsfTrack " << endl;
    mTempTreeElecD0[i]               = (*elecHandle)[i].gsfTrack()->d0();
    mTempTreeElecDz[i]               = (*elecHandle)[i].gsfTrack()->dz();
    mTempTreeElecChargeMode[i]       = (*elecHandle)[i].gsfTrack()->chargeMode();	
    mTempTreeElecPtTrkMode[i]        = (*elecHandle)[i].gsfTrack()->ptMode();
    mTempTreeElecQOverPErrTrkMode[i] = (*elecHandle)[i].gsfTrack()->qoverpModeError();
    mTempTreeElecCharge[i]           = (*elecHandle)[i].gsfTrack()->charge();
    mTempTreeElecPtTrk[i]            = (*elecHandle)[i].gsfTrack()->pt();
    mTempTreeElecQOverPErrTrk[i]     = (*elecHandle)[i].gsfTrack()->qoverpError();
    mTempTreeElecNormChi2[i]         = (*elecHandle)[i].gsfTrack()->normalizedChi2();
    mTempTreeElecLostHits[i]         = (*elecHandle)[i].gsfTrack()->lost();
    mTempTreeElecValidHits[i]        = (*elecHandle)[i].gsfTrack()->found();
    
    edm::LogVerbatim("JaredSusyEvent") << " before asking for trackMomentumAtVtx " << endl;
    
    //try {
    //  mTempTreeElecNCluster[i] = (*elecHandle)[i].numberOfClusters();
    //} catch ( const cms::Exception& e ) {
    //  mTempTreeElecNCluster[i] = -999;
    //  std::stringstream ss;
    //  ss << " cms::Exception caught!"
    //	 << " Invalid edm::Ref<reco::SuperCluster> returned from pat::Electron!" 
    //	 << std::endl 
    //	 << " Setting numberOfClusters to -999!"
    //	 << std::endl 
    //	 << " Output from cms::Exception::what():"
    //	 << std::endl 
    //	 << e.what();
    //  edm::LogWarning("JaredSusyEvent") << ss.str();
    //}
    mTempTreeElecEtaTrk[i] = (*elecHandle)[i].trackMomentumAtVtx().Eta();
    mTempTreeElecPhiTrk[i] = (*elecHandle)[i].trackMomentumAtVtx().Phi();

    // Added protection statement, against missing SuperCluster collection in 2_1_X PatLayer1 samples
    try { 
      mTempTreeElecWidthClusterEta[i] = (*elecHandle)[i].superCluster()->etaWidth();
      mTempTreeElecWidthClusterPhi[i] = (*elecHandle)[i].superCluster()->phiWidth();
    } catch ( const cms::Exception& e ) {
      mTempTreeElecWidthClusterEta[i]=-999.;
      mTempTreeElecWidthClusterPhi[i]=-999.;
      std::stringstream ss;
      ss << " cms::Exception caught!"
	 << " Invalid edm::Ref<reco::SuperCluster> returned from pat::Electron!" 
	 << std::endl 
	 << " Setting ClusterEta and ClusterPhi to -999.!" 
	 << std::endl 
	 << " Output from cms::Exception::what():"
	 << std::endl 
	 << e.what();
      edm::LogWarning("JaredSusyEvent") << ss.str();
    }
    
    mTempTreeElecPinTrk[i] = sqrt((*elecHandle)[i].trackMomentumAtVtx().Mag2());
    mTempTreeElecPoutTrk[i] = sqrt((*elecHandle)[i].trackMomentumOut().Mag2());

    if (&(*(*elecHandle)[i].genLepton())!=0){
      mTempTreeGenElecPdgId[i] = (*elecHandle)[i].genLepton()->pdgId();
      mTempTreeGenElecPx[i]    = (*elecHandle)[i].genLepton()->px();
      mTempTreeGenElecPy[i]    = (*elecHandle)[i].genLepton()->py();
      mTempTreeGenElecPz[i]    = (*elecHandle)[i].genLepton()->pz();
      if(&(*(*elecHandle)[i].genLepton()->mother())!=0){
	mTempTreeGenElecMother[i] = (*elecHandle)[i].genLepton()->mother()->pdgId();
	if ( (*elecHandle)[i].genLepton()->mother()->pdgId() ==  (*elecHandle)[i].genLepton()->pdgId()) 
	  {
	    mTempTreeGenElecMother[i] = (*elecHandle)[i].genLepton()->mother()->mother()->pdgId();
	  }
      }
    }
    else {
      mTempTreeGenElecPdgId[i]  = 999.;
      mTempTreeGenElecPx[i]     = 999.;
      mTempTreeGenElecPy[i]     = 999.;
      mTempTreeGenElecPz[i]     = 999.;
      mTempTreeGenElecMother[i] = 999.;
    }

  }//end loop over Electrons

  // get the muons
  edm::Handle< std::vector<pat::Muon> > muonHandle;
  iEvent.getByLabel(muonTag_, muonHandle);
  if ( !muonHandle.isValid() ) {
    edm::LogWarning("JaredSusyEvent") << "No Muon results for InputTag " << muonTag_;
    return;
  }
  
  edm::LogVerbatim("JaredSusyEvent") << " start reading in muons " << endl;


  // Add the muons
  mTempTreeNmuon= muonHandle->size();
  if ( mTempTreeNmuon > 50 ) mTempTreeNmuon = 50;
  for (int i=0;i<mTempTreeNmuon;i++){
   
    mTempTreeMuonPt[i]                     = (*muonHandle)[i].pt();
    mTempTreeMuonE[i]                      = (*muonHandle)[i].energy();
    mTempTreeMuonEt[i]                     = (*muonHandle)[i].et();
    mTempTreeMuonPx[i]                     = (*muonHandle)[i].momentum().X();
    mTempTreeMuonPy[i]                     = (*muonHandle)[i].momentum().Y();
    mTempTreeMuonPz[i]                     = (*muonHandle)[i].momentum().Z();
    mTempTreeMuonEta[i]                    = (*muonHandle)[i].eta();
    mTempTreeMuonPhi[i]                    = (*muonHandle)[i].phi();
    mTempTreeMuonTrkIso[i]                 = (*muonHandle)[i].trackIso();
    mTempTreeMuonCharge[i]                 = (*muonHandle)[i].charge();
    mTempTreeMuonECalIso[i]                = (*muonHandle)[i].ecalIso();
    mTempTreeMuonHCalIso[i]                = (*muonHandle)[i].hcalIso() ;
    mTempTreeMuonAllIso[i]                 = (*muonHandle)[i].caloIso() ;

    //_muon_isGood_All = muon::isGoodMuon(*the_muon, muon::All);
    //_muon_isGood_TrackerMuonArbitrated = muon::isGoodMuon(*the_muon, muon::TrackerMuonArbitrated);
    //_muon_isGood_AllArbitrated = muon::isGoodMuon(*the_muon, muon::AllArbitrated);
    //_muon_isGood_TMOneStationLoose = muon::isGoodMuon(*the_muon, muon::TMOneStationLoose);
    //_muon_isGood_TMOneStationTight = muon::isGoodMuon(*the_muon, muon::TMOneStationTight);
    //_muon_isGood_TMLastStationOptimizedLowPtLoose = muon::isGoodMuon(*the_muon, muon::TMLastStationOptimizedLowPtLoose);
    //_muon_isGood_TMLastStationOptimizedLowPtTight = muon::isGoodMuon(*the_muon, muon::TMLastStationOptimizedLowPtTight);
    
    //mTempTreeMuonIsGlobal[i]               = (*muonHandle)[i].isGlobalMuon();
    mTempTreeMuonIsGlobal[i] = muon::isGoodMuon((*muonHandle)[i], muon::AllGlobalMuons);
    //mTempTreeMuonIsStandAlone[i]           = (*muonHandle)[i].isStandAloneMuon();
    mTempTreeMuonIsStandAlone[i] = muon::isGoodMuon((*muonHandle)[i], muon::AllStandAloneMuons);
    //mTempTreeMuonIsTracker[i]              = (*muonHandle)[i].isTrackerMuon();
    mTempTreeMuonIsTracker[i] = muon::isGoodMuon((*muonHandle)[i], muon::AllTrackerMuons);
    //mTempTreeMuonIsGlobalTight[i]          = (*muonHandle)[i].isGood(pat::Muon::SelectionType(6));
    mTempTreeMuonIsGlobalTight[i] = muon::isGoodMuon((*muonHandle)[i], muon::GlobalMuonPromptTight);
    //mTempTreeMuonIsTMLastStationLoose[i]   = (*muonHandle)[i].isGood(pat::Muon::SelectionType(7));
    mTempTreeMuonIsTMLastStationLoose[i] = muon::isGoodMuon((*muonHandle)[i], muon::TMLastStationLoose);
    //mTempTreeMuonTMLastStationTight[i]     = (*muonHandle)[i].isGood(pat::Muon::SelectionType(8));
    mTempTreeMuonTMLastStationTight[i] = muon::isGoodMuon((*muonHandle)[i], muon::TMLastStationTight);
    //mTempTreeMuonTM2DCompatibilityLoose[i] = (*muonHandle)[i].isGood(pat::Muon::SelectionType(9));
    mTempTreeMuonTM2DCompatibilityLoose[i] = muon::isGoodMuon((*muonHandle)[i], muon::TM2DCompatibilityLoose);
    //mTempTreeMuonTM2DCompatibilityTight[i] = (*muonHandle)[i].isGood(pat::Muon::SelectionType(10));
    mTempTreeMuonTM2DCompatibilityTight[i] = muon::isGoodMuon((*muonHandle)[i], muon::TM2DCompatibilityTight);
    
    mTempTreeMuonECalIsoDeposit[i]  = (*muonHandle)[i].ecalIsoDeposit()->candEnergy() ;
    mTempTreeMuonHCalIsoDeposit[i]  = (*muonHandle)[i].hcalIsoDeposit()->candEnergy() ;

    // Vertex info is stored only for GlobalMuons (combined muons)
    if((*muonHandle)[i].isGlobalMuon() && (*muonHandle)[i].combinedMuon().isNonnull()){ 

      mTempTreeMuonCombChi2[i] = (*muonHandle)[i].combinedMuon()->chi2();
      mTempTreeMuonCombNdof[i] = (*muonHandle)[i].combinedMuon()->ndof();

      mTempTreeMuonCombVx[i] = (*muonHandle)[i].combinedMuon()->vx();
      mTempTreeMuonCombVy[i] = (*muonHandle)[i].combinedMuon()->vy();
      mTempTreeMuonCombVz[i] = (*muonHandle)[i].combinedMuon()->vz();
      mTempTreeMuonCombD0[i] = (*muonHandle)[i].combinedMuon()->d0();
      mTempTreeMuonCombDz[i] = (*muonHandle)[i].combinedMuon()->dz();

    } else {
      mTempTreeMuonCombVx[i] = 999.;
      mTempTreeMuonCombVy[i] = 999.;
      mTempTreeMuonCombVz[i] = 999.;
      mTempTreeMuonCombD0[i] = 999.;
      mTempTreeMuonCombDz[i] = 999.;
    }

    if((*muonHandle)[i].isStandAloneMuon() && (*muonHandle)[i].standAloneMuon().isNonnull()){
      mTempTreeMuonStandValidHits[i]   = (*muonHandle)[i].standAloneMuon()->found();
      mTempTreeMuonStandLostHits[i]    = (*muonHandle)[i].standAloneMuon()->lost();
      mTempTreeMuonStandPt[i]          = (*muonHandle)[i].standAloneMuon()->pt();
      mTempTreeMuonStandPz[i]          = (*muonHandle)[i].standAloneMuon()->pz();
      mTempTreeMuonStandP[i]           = (*muonHandle)[i].standAloneMuon()->p();
      mTempTreeMuonStandEta[i]         = (*muonHandle)[i].standAloneMuon()->eta();
      mTempTreeMuonStandPhi[i]         = (*muonHandle)[i].standAloneMuon()->phi();
      mTempTreeMuonStandChi[i]         = (*muonHandle)[i].standAloneMuon()->chi2();
      mTempTreeMuonStandCharge[i]      = (*muonHandle)[i].standAloneMuon()->charge();
      mTempTreeMuonStandQOverPError[i] = (*muonHandle)[i].standAloneMuon()->qoverpError();
    } 
    else{
      mTempTreeMuonStandValidHits[i]   = 999.;
      mTempTreeMuonStandLostHits[i]    = 999.;
      mTempTreeMuonStandPt[i]          = 999.;
      mTempTreeMuonStandPz[i]          = 999.;
      mTempTreeMuonStandP[i]           = 999.;
      mTempTreeMuonStandEta[i]         = 999.;
      mTempTreeMuonStandPhi[i]         = 999.;
      mTempTreeMuonStandChi[i]         = 999.;
      mTempTreeMuonStandCharge[i]      = 999.;
      mTempTreeMuonStandQOverPError[i] = 999.;
    }

    if((*muonHandle)[i].isTrackerMuon() && (*muonHandle)[i].track().isNonnull()){
      mTempTreeMuonTrkChiNorm[i]     = (*muonHandle)[i].track()->normalizedChi2();
      mTempTreeMuonTrkValidHits[i]   = (*muonHandle)[i].track()->found();
      mTempTreeMuonTrkLostHits[i]    = (*muonHandle)[i].track()->lost();
      mTempTreeMuonTrkD0[i]          = (*muonHandle)[i].track()->d0();
      mTempTreeMuonTrkPt[i]          = (*muonHandle)[i].track()->pt();
      mTempTreeMuonTrkPz[i]          = (*muonHandle)[i].track()->pz();
      mTempTreeMuonTrkP[i]           = (*muonHandle)[i].track()->p();
      mTempTreeMuonTrkEta[i]         = (*muonHandle)[i].track()->eta();
      mTempTreeMuonTrkPhi[i]         = (*muonHandle)[i].track()->phi();
      mTempTreeMuonTrkChi[i]         = (*muonHandle)[i].track()->chi2();
      mTempTreeMuonTrkCharge[i]      = (*muonHandle)[i].track()->charge();
      mTempTreeMuonTrkQOverPError[i] = (*muonHandle)[i].track()->qoverpError();
      //  mTempTreeMuonTrkOuterZ[i]=(*muonHandle)[i].track()->outerZ();
      //  mTempTreeMuonTrkOuterR[i]=(*muonHandle)[i].track()->outerRadius();

    }
    else{
      mTempTreeMuonTrkChiNorm[i]    = 999.;
      mTempTreeMuonTrkValidHits[i]  = 999.;
      mTempTreeMuonTrkLostHits[i]   = 999.;
      mTempTreeMuonTrkPt[i]         = 999.;
      mTempTreeMuonTrkPz[i]         = 999.;
      mTempTreeMuonTrkP[i]          = 999.;
      mTempTreeMuonTrkEta[i]        = 999.;
      mTempTreeMuonTrkPhi[i]        = 999.;
      mTempTreeMuonTrkChi[i]        = 999.;
      mTempTreeMuonTrkCharge[i]     = 999.;
      mTempTreeMuonTrkQOverPError[i]= 999.;
      mTempTreeMuonTrkOuterZ[i]     = 999.;
      mTempTreeMuonTrkOuterR[i]     = 999.;
    }

    if (&(*(*muonHandle)[i].genLepton())!=0){
      mTempTreeGenMuonPdgId[i] = (*muonHandle)[i].genLepton()->pdgId();
      mTempTreeGenMuonPx[i]    = (*muonHandle)[i].genLepton()->px();
      mTempTreeGenMuonPy[i]    = (*muonHandle)[i].genLepton()->py();
      mTempTreeGenMuonPz[i]    = (*muonHandle)[i].genLepton()->pz();
      if (&(*(*muonHandle)[i].genLepton()->mother())!=0) {
	mTempTreeGenMuonMother[i] = (*muonHandle)[i].genLepton()->mother()->pdgId();
	if ( (*muonHandle)[i].genLepton()->mother()->pdgId() ==  (*muonHandle)[i].genLepton()->pdgId()) 
	  {
	    mTempTreeGenMuonMother[i] = (*muonHandle)[i].genLepton()->mother()->mother()->pdgId();
	  }
      } else {
	mTempTreeGenMuonMother[i] = 999.;
      }
    }
    else{
      mTempTreeGenMuonPdgId[i]  = 999.;
      mTempTreeGenMuonMother[i] = 999.;
      mTempTreeGenMuonPx[i]     = 999.;
      mTempTreeGenMuonPy[i]     = 999.;
      mTempTreeGenMuonPz[i]     = 999.;
    }

  }

  //get pthat of process
  mTempTreePthat = -999.;
  
  Handle<double> genEventScale;
  iEvent.getByLabel( "genEventScale", genEventScale );
  if ( genEventScale.isValid() ) mTempTreePthat = *genEventScale;
  
  //get tracks
  
  //    reco::TrackCollection myTracks;
  // iEvent.getByLabel("generalTracks",myTracks);
  edm::Handle<View<reco::Track> >  myTracks;
  iEvent.getByLabel("generalTracks",myTracks);
  double ptMax_ = 500;
  math::XYZTLorentzVector totalP3;
  for(View<reco::Track>::const_iterator elem = myTracks->begin(); elem
	!= myTracks->end(); ++elem) {
    
    if (!(elem->quality(reco::TrackBase::highPurity))) continue;
    
    double elemPt = elem->pt();
    
    if ( elemPt > ptMax_) continue;
    
    math::XYZTLorentzVector p3(elem->px(),elem->py(),elem->pz(),elem->p());
    totalP3 -= p3;
    
  }
  
  mTempTreeMPTPhi= totalP3.phi();
  mTempTreeMPTPx = totalP3.px();
  mTempTreeMPTPy = totalP3.py();
  mTempTreeMPTPz = totalP3.pz();
   
  //benedetta : PFjets
  edm::Handle< std::vector<pat::Jet> > pfjetHandle;

  double pfsumpx = 0;
  double pfsumpy = 0;
  double pfsumpt = 0;

  if (usePfjets_) {
    iEvent.getByLabel(pfjetTag_, pfjetHandle);
    if(!pfjetHandle.isValid()){
      edm::LogWarning("JaredSusyEvent") << "No PFjet results for InputTag "<< pfjetTag_;
      return;
    }
    
    mTempTreeNPFjet = pfjetHandle->size();
    if( mTempTreeNPFjet > 50) mTempTreeNPFjet=50;
    for(int pf=0; pf< mTempTreeNPFjet; pf++){
      if ((*pfjetHandle)[pf].pt() > 10.) {
	mTempTreePFjetEta[pf]    = (*pfjetHandle)[pf].eta();
	mTempTreePFjetPhi[pf]    = (*pfjetHandle)[pf].phi();
	mTempTreePFjetE[pf]      = (*pfjetHandle)[pf].energy();
	mTempTreePFjetPx[pf]     = (*pfjetHandle)[pf].px();
	mTempTreePFjetPy[pf]     = (*pfjetHandle)[pf].py();
	mTempTreePFjetPz[pf]     = (*pfjetHandle)[pf].pz();
	mTempTreePFjetPt[pf]     = (*pfjetHandle)[pf].pt();
	mTempTreePFjetCharge[pf] = (*pfjetHandle)[pf].charge();
	
	if ((*pfjetHandle)[pf].isCaloJet())
	  mTempTreePFjetFem[pf] = (*pfjetHandle)[pf].emEnergyFraction();
	if ((*pfjetHandle)[pf].isPFJet())
	  mTempTreePFjetFem[pf] = (*pfjetHandle)[pf].neutralEmEnergyFraction()+
	    (*pfjetHandle)[pf].chargedEmEnergyFraction();
	
	
	if ((*pfjetHandle)[pf].pt()>jetMinPt_) {
	  if (fabs((*pfjetHandle)[pf].eta())<jetMaxEta_) pfsumpt+=(*pfjetHandle)[pf].pt();
	  pfsumpx += (*pfjetHandle)[pf].px();
	  pfsumpy += (*pfjetHandle)[pf].py();}
      }
    }
  }
  else mTempTreeNPFjet = 0;

  // get the calo jets
  edm::Handle< std::vector<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetTag_, jetHandle);
  if ( !jetHandle.isValid() ) {
    edm::LogWarning("JaredSusyEvent") << "No Jet results for InputTag " << jetTag_;
    return;
  }

  // get the JPT-corrected pat::Jets
  edm::Handle< std::vector<pat::Jet> > jptHandle;
  iEvent.getByLabel(jptTag_, jptHandle);
  if ( !jptHandle.isValid() ) {
    edm::LogWarning("JaredSusyEvent") << "No JetCorrFactor results for InputTag " << jptTag_;
    std::cout << "No JetCorrFactor results for InputTag " << jptTag_ << std::endl;
    return;
  }

  //get number of jets
  mTempTreeNjets = jetHandle->size();

  // Add the jets
  int i=0;
  double jetsumpx = 0;
  double jetsumpy = 0;
  double jetsumpt = 0;
  double gensumpx = 0;
  double gensumpy = 0;
  double gensumpt = 0;

  if ( mTempTreeNjets >50 ) mTempTreeNjets = 50;
  for (int k=0;k<mTempTreeNjets;k++){
    const pat::Jet& uncorrJet = ((*jetHandle)[k].isCaloJet())? (*jetHandle)[k].correctedJet("RAW"): (*jetHandle)[k];
    
    if ( (*jetHandle)[k].pt() > 10. ) {

      jetsumpt += (*jetHandle)[k].pt();
      jetsumpx += (*jetHandle)[k].momentum().X();
      jetsumpy += (*jetHandle)[k].momentum().Y();
      
      if((*jetHandle)[k].genJet()!= 0) {
	gensumpt += (*jetHandle)[k].genJet()->pt();
	gensumpx += (*jetHandle)[k].genJet()->momentum().X();
	gensumpy += (*jetHandle)[k].genJet()->momentum().Y();}

      if ((*jetHandle)[k].pt()>jetMinPt_) {
	if (fabs((*jetHandle)[k].eta())<jetMaxEta_) {

	  for ( uint16_t n = 0; n < ( jptHandle->size() > 50 ? 50 : jptHandle->size() ); n++ ) {
	    if ( matchJetsByCaloTowers( (*jptHandle)[n], (*jetHandle)[k] ) ) {
	      pat::Jet jpt( (*jptHandle)[n] ); // no corrections by default
	      mTempTreeJetJPTCorrFactor[i] = (jpt.isCaloJet()) ? ( jpt.energy() / uncorrJet.energy() ) : -1 ;
	    }
	  }
	
	  const reco::TrackRefVector & mrTracksInJet = (*jetHandle)[k].associatedTracks();
	
	  mTempTreeJetTrackPt[k]          = 0;
	  mTempTreeJetTrackPhi[k]         = 0;
	  mTempTreeJetTrackPhiWeighted[k] = 0;
	  mTempTreeJetTrackNo[k]          = 0;
	
	  float JetPhi = (*jetHandle)[k].phi();
	
	  for (reco::TrackRefVector::iterator aIter = mrTracksInJet.begin();aIter!= mrTracksInJet.end();aIter++)
	    {
	      mTempTreeJetTrackPt[k] += (*aIter)->pt();
	      float myPhi = (*aIter)->phi();
	      if( JetPhi > 2. ) {
		if(myPhi<0) myPhi = myPhi + 2*TMath::Pi();
	      }
	      if( JetPhi < -2. ) {
		if(myPhi>0) myPhi = myPhi - 2*TMath::Pi();
	      }
	      mTempTreeJetTrackPhiWeighted[k] += (*aIter)->pt()*myPhi;
	      mTempTreeJetTrackPhi[k]         += myPhi;
	      mTempTreeJetTrackNo[k]++;
	    
	    }

	  mTempTreeJetTrackPhiWeighted[k] = mTempTreeJetTrackPhiWeighted[k]/ mTempTreeJetTrackPt[k];
	  mTempTreeJetTrackPhi[k]         = mTempTreeJetTrackPhi[k]/float(mTempTreeJetTrackNo[k]);
 
	  mTempTreeJetsPt[i]  = (*jetHandle)[k].pt();
	  mTempTreeJetsE[i]   = (*jetHandle)[k].energy();
	  mTempTreeJetsEt[i]  = (*jetHandle)[k].et();
	  mTempTreeJetsPx[i]  = (*jetHandle)[k].momentum().X();
	  mTempTreeJetsPy[i]  = (*jetHandle)[k].momentum().Y();
	  mTempTreeJetsPz[i]  = (*jetHandle)[k].momentum().Z();
	  mTempTreeJetsEta[i] = (*jetHandle)[k].eta();
	  mTempTreeJetsPhi[i] = (*jetHandle)[k].phi();
	
	  if ((*jetHandle)[k].isCaloJet())
	    mTempTreeJetsFem[i] = (*jetHandle)[k].emEnergyFraction();
	  if ((*jetHandle)[k].isPFJet())
	    mTempTreeJetsFem[i] = (*jetHandle)[k].neutralEmEnergyFraction()+
	      (*jetHandle)[k].chargedEmEnergyFraction();
      
	  mTempTreeJetsBTag_TkCountHighEff[i] = (*jetHandle)[k].bDiscriminator("trackCountingHighEffBJetTags");
	  mTempTreeJetsBTag_SimpleSecVtx[i]   = (*jetHandle)[k].bDiscriminator("simpleSecondaryVertexBJetTags");
	  mTempTreeJetsBTag_CombSecVtx[i]     = (*jetHandle)[k].bDiscriminator("combinedSecondaryVertexBJetTags");
	  mTempTreeJetPartonFlavour[i]        = (*jetHandle)[k].partonFlavour();
	
	  if((*jetHandle)[k].genJet()!= 0) {
	    mTempTreeGenJetsPt[i]  = (*jetHandle)[k].genJet()->pt();
	    mTempTreeGenJetsE[i]   = (*jetHandle)[k].genJet()->energy();
	    mTempTreeGenJetsEt[i]  = (*jetHandle)[k].genJet()->et();
	    mTempTreeGenJetsPx[i]  = (*jetHandle)[k].genJet()->momentum().X();
	    mTempTreeGenJetsPy[i]  = (*jetHandle)[k].genJet()->momentum().Y();
	    mTempTreeGenJetsPz[i]  = (*jetHandle)[k].genJet()->momentum().z();
	    mTempTreeGenJetsEta[i] = (*jetHandle)[k].genJet()->eta();
	    mTempTreeGenJetsPhi[i] = (*jetHandle)[k].genJet()->phi();
	  }
	  else {
	    mTempTreeGenJetsPt[i]  = -999;
	    mTempTreeGenJetsE[i]   = -999;
	    mTempTreeGenJetsEt[i]  = -999;
	    mTempTreeGenJetsPx[i]  = -999;
	    mTempTreeGenJetsPy[i]  = -999;
	    mTempTreeGenJetsPz[i]  = -999;
	    mTempTreeGenJetsEta[i] = -999;
	    mTempTreeGenJetsPhi[i] = -999;
	  }
	  
	  if((*jetHandle)[k].genParton() != 0){
	    mTempTreeJetPartonId[i]     = (*jetHandle)[k].genParton()->pdgId();
	    mTempTreeJetPartonPx[i]     = (*jetHandle)[k].genParton()->px();
	    mTempTreeJetPartonPy[i]     = (*jetHandle)[k].genParton()->py();
	    mTempTreeJetPartonPz[i]     = (*jetHandle)[k].genParton()->pz();
	    mTempTreeJetPartonEt[i]     = (*jetHandle)[k].genParton()->et();
	    mTempTreeJetPartonPhi[i]    = (*jetHandle)[k].genParton()->phi();
	    mTempTreeJetPartonEta[i]    = (*jetHandle)[k].genParton()->eta();
	    mTempTreeJetPartonEnergy[i] = (*jetHandle)[k].genParton()->energy();
	    mTempTreeJetPartonMother[i] = (*jetHandle)[k].genParton()->mother()->pdgId();
	  }
	  else{
	    mTempTreeJetPartonId[i]     = -999;
	    mTempTreeJetPartonPx[i]     = -999;
	    mTempTreeJetPartonPy[i]     = -999;
	    mTempTreeJetPartonPz[i]     = -999;
	    mTempTreeJetPartonEt[i]     = -999;
	    mTempTreeJetPartonPhi[i]    = -999;
	    mTempTreeJetPartonEta[i]    = -999;
	    mTempTreeJetPartonEnergy[i] = -999;
	    mTempTreeJetPartonMother[i] = -999;
	  }
	
	  // Add the JPT corrs
	  int mTempTreeNjptjets = jptHandle->size();
	  if ( mTempTreeNjptjets > 50 ) mTempTreeNjptjets = 50;
	  for ( int m = 0; m < mTempTreeNjptjets; m++ ) {
	    if( (*jptHandle)[m].originalObjectRef() == (*jetHandle)[k].originalObjectRef() ) {
	      pat::Jet jet = ((*jptHandle)[m].isCaloJet()) ? (*jptHandle)[m].correctedJet("RAW") : (*jptHandle)[m];
	    }
	  }
	  i++;
	}
      }
    }
  }
  
  mTempTreeNjets  = i;
  mTempTreeHt     = jetsumpt;
  mTempTreeMHt    = sqrt(jetsumpx*jetsumpx+jetsumpy*jetsumpy);
  mTempTreePFHt   = pfsumpt;
  mTempTreePFMHt  = sqrt(pfsumpx*pfsumpx+pfsumpy*pfsumpy);
  mTempTreeGenHt  = gensumpt;
  mTempTreeGenMHt = sqrt(gensumpx*gensumpx+gensumpy*gensumpy);

  // Get the hemispheres
  Handle< edm::View<pat::Hemisphere> > hemisphereHandle;
  iEvent.getByLabel("selectedLayer2Hemispheres", hemisphereHandle);
  if ( !hemisphereHandle.isValid() ) {
    edm::LogWarning("JaredSusyEvent") << "No Hemisphere results for InputTag ";
    return;
  }
  const edm::View<pat::Hemisphere>& hemispheres = (*hemisphereHandle); // For simplicity...
  
  mTempTreeNhemispheres = 2;
  for (unsigned int i=0;i <  hemispheres.size() ;i++){
    mTempTreeHemispheresPt[i]  = hemispheres[i].pt();
    mTempTreeHemispheresE[i]   = hemispheres[i].energy();
    mTempTreeHemispheresEt[i]  = hemispheres[i].et();
    mTempTreeHemispheresPx[i]  = hemispheres[i].momentum().X();
    mTempTreeHemispheresPy[i]  = hemispheres[i].momentum().Y();
    mTempTreeHemispheresPz[i]  = hemispheres[i].momentum().Z();
    mTempTreeHemispheresEta[i] = hemispheres[i].eta();
    mTempTreeHemispheresPhi[i] = hemispheres[i].phi();

    for(unsigned int dau = 0; dau < hemispheres[i].numberOfDaughters();dau++){
      for (int k=0;k<mTempTreeNjets;k++){
	mTempTreeJetsHemi[k] = -1;
	if(  hemispheres[i].daughter(dau)->phi() >= mTempTreeJetsPhi[k]-0.0001 &&  hemispheres[i].daughter(dau)->phi() <= mTempTreeJetsPhi[k]+0.0001 )  mTempTreeJetsHemi[k] = i;
      }
    }
  }
  mTempTreeHemispheresdPhi = fabs(reco::deltaPhi(hemispheres[0].phi(),hemispheres[1].phi()));

  //
  // get the PF MET result
  //
  edm::Handle< std::vector<pat::MET> > pfmetHandle;
  iEvent.getByLabel(pfmetTag_, pfmetHandle);
  if ( !pfmetHandle.isValid() ) {
    edm::LogWarning("METEventSelector") << "No Met results for InputTag " << pfmetTag_;
    return;
  }
  //
  // sanity check on collection
  //
  if ( pfmetHandle->size()!=1 ) {
    edm::LogWarning("METEventSelector") << "MET collection size is "
					<< pfmetHandle->size() << " instead of 1";
    return;
  }
 
  if(pfmetHandle->front().genMET()!=NULL) {
    const reco::GenMET* myGenMet = pfmetHandle->front().genMET();
    mTempTreepfMET_Gen[0] = myGenMet->px();
    mTempTreepfMET_Gen[1] = myGenMet->py();
    mTempTreepfMET_Gen[2] = myGenMet->pz();
  }
  else{
    mTempTreepfMET_Gen[0] = -99999999;
    mTempTreepfMET_Gen[1] = -99999999;
    mTempTreepfMET_Gen[2] = -99999999;
  }

  // Do the MET save for full corr no cc pfMET
  mTempTreepfMET_Fullcorr_nocc[0]           = pfmetHandle->front().momentum().X();
  mTempTreepfMET_Fullcorr_nocc[1]           = pfmetHandle->front().momentum().Y();
  mTempTreepfMET_Fullcorr_nocc[2]           = pfmetHandle->front().momentum().z();
  mTempTreepfMETphi_Fullcorr_nocc           = pfmetHandle->front().phi();
  mTempTreepfMET_Fullcorr_nocc_significance = pfmetHandle->front().mEtSig();

  // Do the MET save for no corr no cc pfMET
  mTempTreepfMET_Nocorr_nocc[0] = pfmetHandle->front().corEx();//uncorr to bare bones
  mTempTreepfMET_Nocorr_nocc[1] = pfmetHandle->front().corEy(pat::MET::UncorrectionType(0));//uncorr to bare bones
  mTempTreepfMETphi_Nocorr_nocc = pfmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(0));//uncorr to bare bones

  // Do the MET save for muon corr no cc pfMET
  mTempTreepfMET_Muoncorr_nocc[0] = pfmetHandle->front().corEx(pat::MET::UncorrectionType(1));//uncorr for JEC
  mTempTreepfMET_Muoncorr_nocc[1] = pfmetHandle->front().corEy(pat::MET::UncorrectionType(1));//uncorr for JEC 
  mTempTreepfMETphi_Muoncorr_nocc = pfmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(1));//uncorr for JEC

  // Do the MET save for JEC corr no cc pfMET
  mTempTreepfMET_JECcorr_nocc[0] = pfmetHandle->front().corEx(pat::MET::UncorrectionType(2));//uncorr for muons
  mTempTreepfMET_JECcorr_nocc[1] = pfmetHandle->front().corEy(pat::MET::UncorrectionType(2));//uncorr for muons
  mTempTreepfMETphi_JECcorr_nocc = pfmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(2));//uncorr for muons


  //
  // get the tcMET result
  //
  edm::Handle< std::vector<pat::MET> > tcmetHandle;
  iEvent.getByLabel(tcmetTag_, tcmetHandle);
  if ( !tcmetHandle.isValid() ) {
    edm::LogWarning("METEventSelector") << "No Met results for InputTag " << tcmetTag_;
    return;
  }
  //
  // sanity check on collection
  //
  if ( tcmetHandle->size()!=1 ) {
    edm::LogWarning("METEventSelector") << "tcMET collection size is "
					<< tcmetHandle->size() << " instead of 1";
    return;
  }
 
  if(tcmetHandle->front().genMET()!=NULL) {
    const reco::GenMET* myGenMet = tcmetHandle->front().genMET();
    mTempTreetcMET_Gen[0] = myGenMet->px();
    mTempTreetcMET_Gen[1] = myGenMet->py();
    mTempTreetcMET_Gen[2] = myGenMet->pz();
  }
  else{
    mTempTreetcMET_Gen[0] = -99999999;
    mTempTreetcMET_Gen[1] = -99999999;
    mTempTreetcMET_Gen[2] = -99999999;
  }

  // Do the MET save for full corr no cc tcMET
  mTempTreetcMET_Fullcorr_nocc[0]           = tcmetHandle->front().momentum().X();
  mTempTreetcMET_Fullcorr_nocc[1]           = tcmetHandle->front().momentum().Y();
  mTempTreetcMET_Fullcorr_nocc[2]           = tcmetHandle->front().momentum().z();
  mTempTreetcMETphi_Fullcorr_nocc           = tcmetHandle->front().phi();
  mTempTreetcMET_Fullcorr_nocc_significance = tcmetHandle->front().mEtSig();

  // Do the MET save for no corr no cc tcMET
  mTempTreetcMET_Nocorr_nocc[0] = tcmetHandle->front().corEx();//uncorr to bare bones
  mTempTreetcMET_Nocorr_nocc[1] = tcmetHandle->front().corEy(pat::MET::UncorrectionType(0));//uncorr to bare bones
  mTempTreetcMETphi_Nocorr_nocc = tcmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(0));//uncorr to bare bones

  // Do the MET save for muon corr no cc tcMET
  mTempTreetcMET_Muoncorr_nocc[0] = tcmetHandle->front().corEx(pat::MET::UncorrectionType(1));//uncorr for JEC
  mTempTreetcMET_Muoncorr_nocc[1] = tcmetHandle->front().corEy(pat::MET::UncorrectionType(1));//uncorr for JEC 
  mTempTreetcMETphi_Muoncorr_nocc = tcmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(1));//uncorr for JEC

  // Do the MET save for JEC corr no cc tcMET
  mTempTreetcMET_JECcorr_nocc[0] = tcmetHandle->front().corEx(pat::MET::UncorrectionType(2));//uncorr for muons
  mTempTreetcMET_JECcorr_nocc[1] = tcmetHandle->front().corEy(pat::MET::UncorrectionType(2));//uncorr for muons
  mTempTreetcMETphi_JECcorr_nocc = tcmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(2));//uncorr for muons


  //
  // get the calo (AK5) MET result
  //
  edm::Handle< std::vector<pat::MET> > metHandle;
  iEvent.getByLabel(metTag_, metHandle);
  if ( !metHandle.isValid() ) {
    edm::LogWarning("METEventSelector") << "No Met results for InputTag " << metTag_;
    return;
  }
  //
  // sanity check on collection
  //
  if ( metHandle->size()!=1 ) {
    edm::LogWarning("METEventSelector") << "MET collection size is "
					<< metHandle->size() << " instead of 1";
    return;
  }
 
  if(metHandle->front().genMET()!=NULL) {
    const reco::GenMET* myGenMet = metHandle->front().genMET();
    mTempTreeMET_Gen[0] = myGenMet->px();
    mTempTreeMET_Gen[1] = myGenMet->py();
    mTempTreeMET_Gen[2] = myGenMet->pz();
  }
  else{
    mTempTreeMET_Gen[0] = -99999999;
    mTempTreeMET_Gen[1] = -99999999;
    mTempTreeMET_Gen[2] = -99999999;
  }

  // Do the MET save for full corr no cc MET
  mTempTreeMET_Fullcorr_nocc[0]           = metHandle->front().momentum().X();
  mTempTreeMET_Fullcorr_nocc[1]           = metHandle->front().momentum().Y();
  mTempTreeMET_Fullcorr_nocc[2]           = metHandle->front().momentum().z();
  mTempTreeMETphi_Fullcorr_nocc           = metHandle->front().phi();
  mTempTreeMET_Fullcorr_nocc_significance = metHandle->front().mEtSig();

  // Do the MET save for no corr no cc MET
  mTempTreeMET_Nocorr_nocc[0] = metHandle->front().corEx();//uncorr to bare bones
  mTempTreeMET_Nocorr_nocc[1] = metHandle->front().corEy(pat::MET::UncorrectionType(0));//uncorr to bare bones
  mTempTreeMETphi_Nocorr_nocc = metHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(0));//uncorr to bare bones

  // Do the MET save for muon corr no cc MET
  mTempTreeMET_Muoncorr_nocc[0] = metHandle->front().corEx(pat::MET::UncorrectionType(1));//uncorr for JEC
  mTempTreeMET_Muoncorr_nocc[1] = metHandle->front().corEy(pat::MET::UncorrectionType(1));//uncorr for JEC 
  mTempTreeMETphi_Muoncorr_nocc = metHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(1));//uncorr for JEC

  // Do the MET save for JEC corr no cc MET
  mTempTreeMET_JECcorr_nocc[0] = metHandle->front().corEx(pat::MET::UncorrectionType(2));//uncorr for muons
  mTempTreeMET_JECcorr_nocc[1] = metHandle->front().corEy(pat::MET::UncorrectionType(2));//uncorr for muons
  mTempTreeMETphi_JECcorr_nocc = metHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(2));//uncorr for muons

  //
  // sanity check on collection
  //
  nUncorrMET = 2;
  nFullMET   = 3;

  // Fill the tree only if all preselection conditions are met
  //if(preselection) //mPreselection->Fill();
  mAllData->Fill();
}

//________________________________________________________________________________________
bool 
JaredSusyAnalyzer::matchJetsByCaloTowers( const pat::Jet& jet1,
					  const pat::Jet& jet2 ) {
  
  std::vector< edm::Ptr<CaloTower> > towers1 = jet1.getCaloConstituents();
  std::vector< edm::Ptr<CaloTower> > towers2 = jet2.getCaloConstituents();
  
  if ( towers1.empty() || 
       towers2.empty() || 
       towers1.size() != towers2.size() ) { return false; }
  
  std::vector< edm::Ptr<CaloTower> >::const_iterator ii = towers1.begin();
  std::vector< edm::Ptr<CaloTower> >::const_iterator jj = towers1.end();
  for ( ; ii != jj; ++ii ) {
    std::vector< edm::Ptr<CaloTower> >::const_iterator iii = towers2.begin();
    std::vector< edm::Ptr<CaloTower> >::const_iterator jjj = towers2.end();
    for ( ; iii != jjj; ++iii ) { if ( *iii == *ii ) { break; } }
    if ( iii == towers2.end() ) { return false; }
  }

  return true;

}

//________________________________________________________________________________________
void 
JaredSusyAnalyzer::beginJob(const edm::EventSetup&) {}

//________________________________________________________________________________________
void 
JaredSusyAnalyzer::endJob() {

}

void
JaredSusyAnalyzer::printHLTreport( void ) {

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
JaredSusyAnalyzer::initPlots() {

  std::ostringstream variables; // Container for all variables
  
  // 1. Event variables
  variables << "weight:process";
  
  // Register this ntuple
  edm::Service<TFileService> fs;

  // Now we add some additional ones for the dijet analysis
  //mPreselection = fs->make<TTree>( "Preselection", "data after preselection" );
  mAllData = fs->make<TTree>( "AllData", "data after cuts" );
  mAllData->SetAutoSave(10);

  // Add the branches
  mAllData->Branch("run",&mTempTreeRun,"run/int");
  mAllData->Branch("event",&mTempTreeEvent,"event/int");

  mAllData->Branch("nHLT",&nHLT,"nHLT/I");
  mAllData->Branch("HLTArray",HLTArray,"HLTArray[nHLT]/I");


  //Trigger information
  mAllData->Branch("HLT1JET",&mTempTreeHLT1JET,"HLT1JET/bool");
  mAllData->Branch("HLT2JET",&mTempTreeHLT2JET,"HLT2JET/bool");
  mAllData->Branch("HLT1MET1HT",&mTempTreeHLT1MET1HT,"HLT1MET1HT/bool");
  mAllData->Branch("HLT1MUON",&mTempTreeHLT1Muon,"HLT1MUON/bool");

  //general MET information
  mAllData->Branch("nFullMET",  &nFullMET,  "nFullMET/int");
  mAllData->Branch("nUncorrMET",&nUncorrMET,"nUncorrMET/int");

  //pfMET information
  mAllData->Branch("pfMET_fullcorr_nocc",   mTempTreepfMET_Fullcorr_nocc,    "mTempTreepfMET_Fullcorr_nocc[nFullMET]/double");
  mAllData->Branch("pfMETphi_fullcorr_nocc",&mTempTreepfMETphi_Fullcorr_nocc,"mTempTreepfMETphi_Fullcorr_nocc/double");

  mAllData->Branch("pfMET_nocorr_nocc",  mTempTreepfMET_Nocorr_nocc,  "mTempTreepfMET_Nocorr_nocc[nUncorrMET]/double");
  mAllData->Branch("pfMET_muoncorr_nocc",mTempTreepfMET_Muoncorr_nocc,"mTempTreepfMET_Muoncorr_nocc[nUncorrMET]/double");
  mAllData->Branch("pfMET_jeccorr_nocc", mTempTreepfMET_JECcorr_nocc, "mTempTreepfMET_JECcorr_nocc[nUncorrMET]/double");

  mAllData->Branch("pfMETphi_fullcorr_nocc",&mTempTreepfMETphi_Fullcorr_nocc,"mTempTreepfMETphi_Fullcorr_nocc/double");
  mAllData->Branch("pfMETphi_muoncorr_nocc",&mTempTreepfMETphi_Muoncorr_nocc,"mTempTreepfMETphi_Muoncorr_nocc/double");
  mAllData->Branch("pfMETphi_jeccorr_nocc", &mTempTreepfMETphi_JECcorr_nocc, "mTempTreepfMETphi_JECcorr_nocc/double");
  mAllData->Branch("pfMET_Fullcorr_nocc_significance",&mTempTreepfMET_Fullcorr_nocc_significance,"mTempTreepfMET_Fullcorr_nocc_significance/double");
  mAllData->Branch("GenpfMET",&mTempTreepfMET_Gen,"mTempTreepfMET_Gen[3]/double",6400);

  //tcMET information
  mAllData->Branch("tcMET_fullcorr_nocc",   mTempTreetcMET_Fullcorr_nocc,    "mTempTreetcMET_Fullcorr_nocc[nFullMET]/double");
  mAllData->Branch("tcMETphi_fullcorr_nocc",&mTempTreetcMETphi_Fullcorr_nocc,"mTempTreetcMETphi_Fullcorr_nocc/double");

  mAllData->Branch("tcMET_nocorr_nocc",  mTempTreetcMET_Nocorr_nocc,  "mTempTreetcMET_Nocorr_nocc[nUncorrMET]/double");
  mAllData->Branch("tcMET_muoncorr_nocc",mTempTreetcMET_Muoncorr_nocc,"mTempTreetcMET_Muoncorr_nocc[nUncorrMET]/double");
  mAllData->Branch("tcMET_jeccorr_nocc", mTempTreetcMET_JECcorr_nocc, "mTempTreetcMET_JECcorr_nocc[nUncorrMET]/double");

  mAllData->Branch("tcMETphi_fullcorr_nocc",&mTempTreetcMETphi_Fullcorr_nocc,"mTempTreetcMETphi_Fullcorr_nocc/double");
  mAllData->Branch("tcMETphi_muoncorr_nocc",&mTempTreetcMETphi_Muoncorr_nocc,"mTempTreetcMETphi_Muoncorr_nocc/double");
  mAllData->Branch("tcMETphi_jeccorr_nocc", &mTempTreetcMETphi_JECcorr_nocc, "mTempTreetcMETphi_JECcorr_nocc/double");
  mAllData->Branch("tcMET_Fullcorr_nocc_significance",&mTempTreetcMET_Fullcorr_nocc_significance,"mTempTreetcMET_Fullcorr_nocc_significance/double");
  mAllData->Branch("GentcMET",&mTempTreetcMET_Gen,"mTempTreetcMET_Gen[3]/double",6400);

  //calo MET information
  mAllData->Branch("MET_fullcorr_nocc",   mTempTreeMET_Fullcorr_nocc,    "mTempTreeMET_Fullcorr_nocc[nFullMET]/double");
  mAllData->Branch("METphi_fullcorr_nocc",&mTempTreeMETphi_Fullcorr_nocc,"mTempTreeMETphi_Fullcorr_nocc/double");

  mAllData->Branch("MET_nocorr_nocc",  mTempTreeMET_Nocorr_nocc,  "mTempTreeMET_Nocorr_nocc[nUncorrMET]/double");
  mAllData->Branch("MET_muoncorr_nocc",mTempTreeMET_Muoncorr_nocc,"mTempTreeMET_Muoncorr_nocc[nUncorrMET]/double");
  mAllData->Branch("MET_jeccorr_nocc", mTempTreeMET_JECcorr_nocc, "mTempTreeMET_JECcorr_nocc[nUncorrMET]/double");

  mAllData->Branch("METphi_fullcorr_nocc",&mTempTreeMETphi_Fullcorr_nocc,"mTempTreeMETphi_Fullcorr_nocc/double");
  mAllData->Branch("METphi_muoncorr_nocc",&mTempTreeMETphi_Muoncorr_nocc,"mTempTreeMETphi_Muoncorr_nocc/double");
  mAllData->Branch("METphi_jeccorr_nocc", &mTempTreeMETphi_JECcorr_nocc, "mTempTreeMETphi_JECcorr_nocc/double");
  mAllData->Branch("MET_Fullcorr_nocc_significance",&mTempTreeMET_Fullcorr_nocc_significance,"mTempTreeMET_Fullcorr_nocc_significance/double");
  mAllData->Branch("GenMET",&mTempTreeMET_Gen,"mTempTreeMET_Gen[3]/double",6400);

  //cc MET information
  mAllData->Branch("MET_fullcorr_cc",     mTempTreeMET_Fullcorr_cc,      "mTempTreeMET_Fullcorr_cc[nFullMET]/double");
  mAllData->Branch("METphi_fullcorr_cc",  &mTempTreeMETphi_Fullcorr_cc,  "mTempTreeMETphi_Fullcorr_cc/double");

  //other information
  mAllData->Branch("evtWeight",&mTempTreeEventWeight,"evtWeight/double");
  mAllData->Branch("procID",&mTempTreeProcID,"procID/int");
  mAllData->Branch("pthat",&mTempTreePthat,"pthat/double");

  mAllData->Branch("nVtx",      &mTempTreenVtx,   "nVtx/int");
  mAllData->Branch("VertexChi2",mTempTreeVtxChi2,"VertexChi2[nVtx]/double");
  mAllData->Branch("VertexNdof",mTempTreeVtxNdof,"VertexNdof[nVtx]/double");
  mAllData->Branch("VertexNormalizedChi2",mTempTreeVtxNormalizedChi2,"VertexNormalizedChi2[nVtx]/double");
  mAllData->Branch("VertexNormalizedChi2",mTempTreeVtxNormalizedChi2,"VertexNormalizedChi2[nVtx]/double");

  mAllData->Branch("VertexX", mTempTreeVtxX, "VertexX[nVtx]/double");
  mAllData->Branch("VertexY", mTempTreeVtxY, "VertexY[nVtx]/double");
  mAllData->Branch("VertexZ", mTempTreeVtxZ, "VertexZ[nVtx]/double");
  mAllData->Branch("VertexdX",mTempTreeVtxdX,"VertexdX[nVtx]/double");
  mAllData->Branch("VertexdY",mTempTreeVtxdY,"VertexdY[nVtx]/double");
  mAllData->Branch("VertexdZ",mTempTreeVtxdZ,"VertexdZ[nVtx]/double");

  //add hemispheres
  mAllData->Branch("Nhemispheres" ,&mTempTreeNhemispheres ,"Nhemispheres/int");  
  mAllData->Branch("HemisphereE" ,mTempTreeHemispheresE ,"HemisphereE[Nhemispheres]/double");
  mAllData->Branch("HemisphereEt",mTempTreeHemispheresEt,"HemisphereEt[Nhemispheres]/double");
  mAllData->Branch("Hemispherept",mTempTreeHemispheresPt,"Hemispherept[Nhemispheres]/double");
  mAllData->Branch("Hemispherepx",mTempTreeHemispheresPx,"Hemispherepx[Nhemispheres]/double");
  mAllData->Branch("Hemispherepy",mTempTreeHemispheresPy,"Hemispherepy[Nhemispheres]/double");
  mAllData->Branch("Hemispherepz",mTempTreeHemispheresPz,"Hemispherepz[Nhemispheres]/double");
  mAllData->Branch("Hemisphereeta",mTempTreeHemispheresEta,"Hemisphereeta[Nhemispheres]/double");
  mAllData->Branch("Hemispherephi",mTempTreeHemispheresPhi,"Hemispherephi[Nhemispheres]/double");
  mAllData->Branch("Hemispheredphi",&mTempTreeHemispheresdPhi,"Hemispheredphi/double");

  //add uncorrected non-cc jets
  mAllData->Branch("Njets",&mTempTreeNjets ,"Njets/int");  
  mAllData->Branch("Ht"   ,&mTempTreeHt,"Ht/double");
  mAllData->Branch("MHt"  ,&mTempTreeMHt,"MHt/double");
  mAllData->Branch("JetE" ,mTempTreeJetsE ,"JetE[Njets]/double");
  mAllData->Branch("JetEt",mTempTreeJetsEt,"JetEt[Njets]/double");
  mAllData->Branch("Jetpt",mTempTreeJetsPt,"Jetpt[Njets]/double");
  mAllData->Branch("Jetpx",mTempTreeJetsPx,"Jetpx[Njets]/double");
  mAllData->Branch("Jetpy",mTempTreeJetsPy,"Jetpy[Njets]/double");
  mAllData->Branch("Jetpz",mTempTreeJetsPz,"Jetpz[Njets]/double");
  mAllData->Branch("Jeteta",mTempTreeJetsEta,"Jeteta[Njets]/double");
  mAllData->Branch("Jetphi",mTempTreeJetsPhi,"Jetphi[Njets]/double");
  mAllData->Branch("JetFem",mTempTreeJetsFem,"JetFem[Njets]/double");
  mAllData->Branch("JetHemi",mTempTreeJetsHemi,"JetHemi[Njets]/int");

  //jet correction factors
  mAllData->Branch("Jet_MCcorrFactor",mTempTreeJetMCCorrFactor,"mTempTreeJetMCCorrFactor[Njets]/double");
  mAllData->Branch("Jet_JPTcorrFactor",mTempTreeJetJPTCorrFactor,"mTempTreeJetJPTCorrFactor[Njets]/double");
  
  //b-tagging information
  mAllData->Branch("JetBTag_TkCountHighEff",mTempTreeJetsBTag_TkCountHighEff,"JetBTag_TkCountHighEff[Njets]/float");
  mAllData->Branch("JetBTag_SimpleSecVtx",  mTempTreeJetsBTag_SimpleSecVtx,  "JetBTag_SimpleSecVtx[Njets]/float");
  mAllData->Branch("JetBTag_CombSecVtx",    mTempTreeJetsBTag_CombSecVtx,    "JetBTag_CombSecVtx[Njets]/float");

  //information about associated gen jets
  mAllData->Branch("GenHt"   ,&mTempTreeGenHt ,"GenHt/double");
  mAllData->Branch("GenMHt"  ,&mTempTreeGenMHt ,"GenMHt/double");
  mAllData->Branch("GenJetE" ,mTempTreeGenJetsE ,"GenJetE[Njets]/double");
  mAllData->Branch("GenJetEt",mTempTreeGenJetsEt,"GenJetEt[Njets]/double");
  mAllData->Branch("GenJetpt",mTempTreeGenJetsPt,"GenJetpt[Njets]/double");
  mAllData->Branch("GenJetpx",mTempTreeGenJetsPx,"GenJetpx[Njets]/double");
  mAllData->Branch("GenJetpy",mTempTreeGenJetsPy,"GenJetpy[Njets]/double");
  mAllData->Branch("GenJetpz",mTempTreeGenJetsPz,"GenJetpz[Njets]/double");
  mAllData->Branch("GenJeteta",mTempTreeGenJetsEta,"GenJeteta[Njets]/double");
  mAllData->Branch("GenJetphi",mTempTreeGenJetsPhi,"GenJetphi[Njets]/double");
  
  //information about associated partons
  mAllData->Branch("JetPartonId" ,mTempTreeJetPartonId ,"JetPartonId[Njets]/int"); 
  mAllData->Branch("JetPartonMother" ,mTempTreeJetPartonMother ,"JetPartonMother[Njets]/int"); 
  mAllData->Branch("JetPartonPx", mTempTreeJetPartonPx ,"JetPartonPx[Njets]/double"); 
  mAllData->Branch("JetPartonPy" ,mTempTreeJetPartonPy ,"JetPartonPy[Njets]/double"); 
  mAllData->Branch("JetPartonPz" ,mTempTreeJetPartonPz ,"JetPartonPz[Njets]/double"); 
  mAllData->Branch("JetPartonEt" ,mTempTreeJetPartonEt ,"JetPartonEt[Njets]/double"); 
  mAllData->Branch("JetPartonE" ,mTempTreeJetPartonEnergy ,"JetPartonE[Njets]/double"); 
  mAllData->Branch("JetPartonPhi" ,mTempTreeJetPartonPhi ,"JetPartonPhi[Njets]/double"); 
  mAllData->Branch("JetPartonEta" ,mTempTreeJetPartonEta ,"JetPartonEta[Njets]/double"); 
  mAllData->Branch("JetPartonFlavour",mTempTreeJetPartonFlavour,"JetPartonFlavour[Njets]/int");
  mAllData->Branch("JetTrackPt",mTempTreeJetTrackPt,"JetTrackPt[Njets]/double"); 
  mAllData->Branch("JetTrackPhi",mTempTreeJetTrackPhi,"JetTrackPhi[Njets]/double"); 
  mAllData->Branch("JetTrackPhiWeighted",mTempTreeJetTrackPhiWeighted,"JetTrackPhiWeighted[Njets]/double"); 
  mAllData->Branch("JetTrackNo",mTempTreeJetTrackNo,"JetTrackNo[Njets]/int");

  //add electrons
  mAllData->Branch("Nelec" , &mTempTreeNelec, "Nelec/int");  
  mAllData->Branch("ElecE" , mTempTreeElecE , "ElecE[Nelec]/double");
  mAllData->Branch("ElecEt", mTempTreeElecEt, "ElecEt[Nelec]/double");
  mAllData->Branch("Elecpt", mTempTreeElecPt, "Elecpt[Nelec]/double");
  mAllData->Branch("Elecpx", mTempTreeElecPx, "Elecpx[Nelec]/double");
  mAllData->Branch("Elecpy", mTempTreeElecPy, "Elecpy[Nelec]/double");
  mAllData->Branch("Elecpz", mTempTreeElecPz, "Elecpz[Nelec]/double");
  mAllData->Branch("Eleceta",mTempTreeElecEta,"Eleceta[Nelec]/double");
  mAllData->Branch("Elecphi",mTempTreeElecPhi,"Elecphi[Nelec]/double");

  mAllData->Branch("ElecCharge",    mTempTreeElecCharge,   "ElecCharge[Nelec]/double");
  mAllData->Branch("ElecTrkIso",    mTempTreeElecTrkIso,   "ElecTrkIso[Nelec]/double");
  mAllData->Branch("ElecECalIso",   mTempTreeElecECalIso,  "ElecECalIso[Nelec]/double");
  mAllData->Branch("ElecHCalIso",   mTempTreeElecHCalIso,  "ElecHCalIso[Nelec]/double");
  mAllData->Branch("ElecAllIso",    mTempTreeElecAllIso,   "ElecAllIso[Nelec]/double");
  mAllData->Branch("ElecTrkChiNorm",mTempTreeElecNormChi2 ,"ElecTrkChiNorm[Nelec]/double");

  mAllData->Branch("ElecECalIsoDeposit", mTempTreeElecECalIsoDeposit,"ElecECalIsoDeposit[Nelec]/double");
  mAllData->Branch("ElecHCalIsoDeposit", mTempTreeElecHCalIsoDeposit,"ElecHCalIsoDeposit[Nelec]/double");

  mAllData->Branch("ElecIdLoose",   mTempTreeElecIdLoose,   "ElecIdLoose [Nelec]/double");
  mAllData->Branch("ElecIdTight",   mTempTreeElecIdTight,   "ElecIdTight [Nelec]/double");
  mAllData->Branch("ElecIdRobLoose",mTempTreeElecIdRobLoose,"ElecIdRobLoose [Nelec]/double");
  mAllData->Branch("ElecIdRobTight",mTempTreeElecIdRobTight,"ElecIdRobTight [Nelec]/double");
  mAllData->Branch("ElecChargeMode",mTempTreeElecChargeMode,"ElecChargeMode [Nelec]/double");
  mAllData->Branch("ElecPtMode",    mTempTreeElecPtTrkMode, "ElecPtMode [Nelec]/double");

  mAllData->Branch("ElecQOverPErrTrkMode",mTempTreeElecQOverPErrTrkMode,"ElecQOverPErrTrkMode [Nelec]/double");
  mAllData->Branch("ElecCaloEnergy",      mTempTreeElecCaloEnergy,      "ElecCaloEnergy[Nelec]/double");

  mAllData->Branch("ElecHOverE",mTempTreeElecHOverE,"ElecHOverE[Nelec]/double");
  mAllData->Branch("ElecVx",    mTempTreeElecVx,    "ElecVx[Nelec]/double");
  mAllData->Branch("ElecVy",    mTempTreeElecVy,    "ElecVy[Nelec]/double");
  mAllData->Branch("ElecVz",    mTempTreeElecVz,    "ElecVz[Nelec]/double");
  mAllData->Branch("ElecD0",    mTempTreeElecD0,    "ElecD0[Nelec]/double");
  mAllData->Branch("ElecDz",    mTempTreeElecDz,    "ElecDz[Nelec]/double");
  mAllData->Branch("ElecPtTrk", mTempTreeElecPtTrk, "ElecPtTrk[Nelec]/double");

  mAllData->Branch("ElecQOverPErrTrk", mTempTreeElecQOverPErrTrk,"ElecQOverPErrTrk[Nelec]/double");
  mAllData->Branch("ElecPinTrk",       mTempTreeElecPinTrk,      "ElecPinTrk[Nelec]/double");
  mAllData->Branch("ElecPoutTrk",      mTempTreeElecPoutTrk,     "ElecPoutTrk[Nelec]/double"); 
  mAllData->Branch("ElecLostHits",     mTempTreeElecLostHits,    "ElecLostHits[Nelec]/double"); 
  mAllData->Branch("ElecValidHits",    mTempTreeElecValidHits,   "ElecValidHits[Nelec]/double"); 
  mAllData->Branch("ElecNCluster",     mTempTreeElecNCluster,    "ElecNCluster[Nelec]/double"); 
  mAllData->Branch("ElecEtaTrk",       mTempTreeElecEtaTrk,      "ElecEtaTrk[Nelec]/double"); 
  mAllData->Branch("ElecPhiTrk",       mTempTreeElecPhiTrk,      "ElecPhiTrk[Nelec]/double"); 

  mAllData->Branch("ElecWidthClusterEta",mTempTreeElecWidthClusterEta,"ElecWidthClusterEta[Nelec]/double"); 
  mAllData->Branch("ElecWidthClusterPhi",mTempTreeElecWidthClusterPhi,"ElecWidthClusterPhi[Nelec]/double"); 

  mAllData->Branch("ElecGenPdgId", mTempTreeGenElecPdgId, "ElecGenPdgId[Nelec]/double");
  mAllData->Branch("ElecGenMother",mTempTreeGenElecMother,"ElecGenMother[Nelec]/double");
  mAllData->Branch("ElecGenPx",    mTempTreeGenElecPx,    "ElecGenPx[Nelec]/double");
  mAllData->Branch("ElecGenPy",    mTempTreeGenElecPy,    "ElecGenPy[Nelec]/double");
  mAllData->Branch("ElecGenPz",    mTempTreeGenElecPz,    "ElecGenPz[Nelec]/double");

  //add muons
  mAllData->Branch("Nmuon",         &mTempTreeNmuon,        "Nmuon/int");  
  mAllData->Branch("MuonE",         mTempTreeMuonE,         "MuonE[Nmuon]/double");
  mAllData->Branch("MuonEt",        mTempTreeMuonEt,        "MuonEt[Nmuon]/double");
  mAllData->Branch("Muonpt",        mTempTreeMuonPt,        "Muonpt[Nmuon]/double");
  mAllData->Branch("Muonpx",        mTempTreeMuonPx,        "Muonpx[Nmuon]/double");
  mAllData->Branch("Muonpy",        mTempTreeMuonPy,        "Muonpy[Nmuon]/double");
  mAllData->Branch("Muonpz",        mTempTreeMuonPz,        "Muonpz[Nmuon]/double");
  mAllData->Branch("Muoneta",       mTempTreeMuonEta,       "Muoneta[Nmuon]/double");
  mAllData->Branch("Muonphi",       mTempTreeMuonPhi,       "Muonphi[Nmuon]/double");
  mAllData->Branch("MuonCharge",    mTempTreeMuonCharge,    "MuonCharge[Nmuon]/double");
  mAllData->Branch("MuonTrkIso",    mTempTreeMuonTrkIso,    "MuonTrkIso[Nmuon]/double");
  mAllData->Branch("MuonECalIso",   mTempTreeMuonECalIso,   "MuonECalIso[Nmuon]/double");
  mAllData->Branch("MuonHCalIso",   mTempTreeMuonHCalIso,   "MuonHCalIso[Nmuon]/double");
  mAllData->Branch("MuonAllIso",    mTempTreeMuonAllIso,    "MuonAllIso[Nmuon]/double");
  mAllData->Branch("MuonTrkChiNorm",mTempTreeMuonTrkChiNorm,"MuonTrkChiNorm[Nmuon]/double");

  mAllData->Branch("MuonECalIsoDeposit", mTempTreeMuonECalIsoDeposit,"MuonECalIsoDeposit[Nmuon]/double");
  mAllData->Branch("MuonHCalIsoDeposit", mTempTreeMuonHCalIsoDeposit,"MuonHCalIsoDeposit[Nmuon]/double");

  mAllData->Branch("MuonIsGlobal",                mTempTreeMuonIsGlobal,              "mTempTreeMuonIsGlobal[Nmuon]/bool");
  mAllData->Branch("MuonIsStandAlone",            mTempTreeMuonIsStandAlone,          "mTempTreeMuonIsStandAlone[Nmuon]/bool");
  mAllData->Branch("MuonIsGlobalTight",           mTempTreeMuonIsGlobalTight,         "mTempTreeMuonIsGlobalTight[Nmuon]/bool");
  mAllData->Branch("MuonIsTMLastStationLoose",    mTempTreeMuonIsTMLastStationLoose,  "mTempTreeMuonIsTMLastStationLoose[Nmuon]/bool");
  mAllData->Branch("MuonIsTracker",               mTempTreeMuonIsTracker,             "mTempTreeMuonIsTracker[Nmuon]/bool");
  mAllData->Branch("MuonIsTMLastStationTight",    mTempTreeMuonTMLastStationTight,    "MuonIsTMLastStationTight[Nmuon]/bool");
  mAllData->Branch("MuonIsTM2DCompatibilityLoose",mTempTreeMuonTM2DCompatibilityLoose,"MuonIsTM2DCompatibilityLoose[Nmuon]/bool");
  mAllData->Branch("MuonIsTM2DCompatibilityTight",mTempTreeMuonTM2DCompatibilityTight,"MuonIsTM2DCompatibilityTight[Nmuon]/bool");

  //  mAllData->Branch("MuonId",mTempTreeMuonId,"mTempTreeMuonId[Nmuon]/double");
  mAllData->Branch("MuonCombChi2",mTempTreeMuonCombChi2,"mTempTreeMuonCombChi2[Nmuon]/double");
  mAllData->Branch("MuonCombNdof",mTempTreeMuonCombNdof,"mTempTreeMuonCombNdof[Nmuon]/double");
  mAllData->Branch("MuonCombVx",  mTempTreeMuonCombVx,  "mTempTreeMuonCombVx[Nmuon]/double");
  mAllData->Branch("MuonCombVy",  mTempTreeMuonCombVy,  "mTempTreeMuonCombVy[Nmuon]/double");
  mAllData->Branch("MuonCombVz",  mTempTreeMuonCombVz,  "mTempTreeMuonCombVz[Nmuon]/double");
  mAllData->Branch("MuonCombD0",  mTempTreeMuonCombD0,  "mTempTreeMuonCombD0[Nmuon]/double");
  mAllData->Branch("MuonCombDz",  mTempTreeMuonCombDz,  "mTempTreeMuonCombDz[Nmuon]/double");

  mAllData->Branch("MuonStandValidHits",  mTempTreeMuonStandValidHits,  "mTempTreeMuonStandValidHits[Nmuon]/double");
  mAllData->Branch("MuonStandLostHits",   mTempTreeMuonStandLostHits,   "mTempTreeMuonStandLostHits[Nmuon]/double");
  mAllData->Branch("MuonStandPt",         mTempTreeMuonStandPt,         "mTempTreeMuonStandPt[Nmuon]/double");
  mAllData->Branch("MuonStandPz",         mTempTreeMuonStandPz,         "mTempTreeMuonStandPz[Nmuon]/double");
  mAllData->Branch("MuonStandP",          mTempTreeMuonStandP,          "mTempTreeMuonStandP[Nmuon]/double");
  mAllData->Branch("MuonStandEta",        mTempTreeMuonStandEta,        "mTempTreeMuonStandEta[Nmuon]/double");
  mAllData->Branch("MuonStandPhi",        mTempTreeMuonStandPhi,        "mTempTreeMuonStandPhi[Nmuon]/double");
  mAllData->Branch("MuonStandCharge",     mTempTreeMuonStandCharge,     "mTempTreeMuonStandCharge[Nmuon]/double");
  mAllData->Branch("MuonStandChi",        mTempTreeMuonStandChi,        "mTempTreeMuonStandChi[Nmuon]/double");
  mAllData->Branch("MuonStandQOverPError",mTempTreeMuonStandQOverPError,"mTempTreeMuonStandQOverPError[Nmuon]/double");

  mAllData->Branch("MuonTrkValidHits",  mTempTreeMuonTrkValidHits,  "mTempTreeMuonTrkValidHits[Nmuon]/double");
  mAllData->Branch("MuonTrkLostHits",   mTempTreeMuonTrkLostHits,   "mTempTreeMuonTrkLostHits[Nmuon]/double");
  mAllData->Branch("MuonTrkD0",         mTempTreeMuonTrkD0,         "mTempTreeMuonTrkD0[Nmuon]/double");
  mAllData->Branch("MuonTrkPt",         mTempTreeMuonTrkPt,         "mTempTreeMuonTrkPt[Nmuon]/double");
  mAllData->Branch("MuonTrkPz",         mTempTreeMuonTrkPz,         "mTempTreeMuonTrkPz[Nmuon]/double");
  mAllData->Branch("MuonTrkP",          mTempTreeMuonTrkP,          "mTempTreeMuonTrkP[Nmuon]/double");
  mAllData->Branch("MuonTrkEta",        mTempTreeMuonTrkEta,        "mTempTreeMuonTrkEta[Nmuon]/double");
  mAllData->Branch("MuonTrkPhi",        mTempTreeMuonTrkPhi,        "mTempTreeMuonTrkPhi[Nmuon]/double");
  mAllData->Branch("MuonTrkCharge",     mTempTreeMuonTrkCharge,     "mTempTreeMuonTrkCharge[Nmuon]/double");
  mAllData->Branch("MuonTrkChi",        mTempTreeMuonTrkChi,        "mTempTreeMuonTrkChi[Nmuon]/double");
  mAllData->Branch("MuonTrkQOverPError",mTempTreeMuonTrkQOverPError,"mTempTreeMuonTrkQOverPError[Nmuon]/double"); 
  mAllData->Branch("MuonTrkOuterZ",     mTempTreeMuonTrkOuterZ,     "mTempTreeMuonOuterZ[Nmuon]/double");
  mAllData->Branch("MuonTrkOuterR",     mTempTreeMuonTrkOuterR,     "mTempTreeMuonOuterR[Nmuon]/double");

  mAllData->Branch("MuonGenMother",mTempTreeGenMuonMother,"MuonGenMother[Nmuon]/double");
  mAllData->Branch("MuonGenPx",    mTempTreeGenMuonPx,    "MuonGenPx[Nmuon]/double");
  mAllData->Branch("MuonGenPy",    mTempTreeGenMuonPy,    "MuonGenPy[Nmuon]/double");
  mAllData->Branch("MuonGenPz",    mTempTreeGenMuonPz,    "MuonGenPz[Nmuon]/double");

  // PFjets
  mAllData->Branch("NPFjet" ,&mTempTreeNPFjet,"NPFjet/int");
  mAllData->Branch("PFHt"   ,&mTempTreePFHt,"PFHt/double");
  mAllData->Branch("PFMHt"  ,&mTempTreePFMHt,"PFMHt/double");
  mAllData->Branch("PFjetEta",&mTempTreePFjetEta,"PFjetEta[NPFjet]/double");
  mAllData->Branch("PFjetPhi",&mTempTreePFjetPhi,"PFjetPhi[NPFjet]/double");
  mAllData->Branch("PFjetE",&mTempTreePFjetE,"PFjetE[NPFjet]/double");
  mAllData->Branch("PFjetPx",&mTempTreePFjetPx,"PFjetPx[NPFjet]/double");
  mAllData->Branch("PFjetPy",&mTempTreePFjetPy,"PFjetPy[NPFjet]/double");
  mAllData->Branch("PFjetPz",&mTempTreePFjetPz,"PFjetPz[NPFjet]/double");
  mAllData->Branch("PFjetPt",&mTempTreePFjetPt,"PFjetPt[NPFjet]/double");
  mAllData->Branch("PFjetCharge",&mTempTreePFjetCharge,"PFjetCharge[NPFjet]/double");
  mAllData->Branch("PFjetFem",&mTempTreePFjetFem,"PFjetFem[NPFjet]/double");
  
  mAllData->Branch("genN",&length,"genN/int");
  mAllData->Branch("genid",ids,"ids[genN]/int");
  mAllData->Branch("genMother",refs,"refs[genN]/int");
  mAllData->Branch("genE",genE,"genE[genN]/float");
  mAllData->Branch("genPx",genPx,"genPx[genN]/float");
  mAllData->Branch("genPy",genPy,"genPy[genN]/float");
  mAllData->Branch("genPz",genPz,"genPz[genN]/float");
  mAllData->Branch("genStatus",genStatus,"genStatus[genN]/int");

  mAllData->Branch("genLepN",&genLepLength,"genLepN/int");
  mAllData->Branch("genLepId",genLepIds,"genLepIds[genLepN]/int");
  mAllData->Branch("genLepMother",genLepRefs,"genLepRefs[genLepN]/int");
  mAllData->Branch("genLepE",genLepE,"genLepE[genLepN]/float");
  mAllData->Branch("genLepPx",genLepPx,"genLepPx[genLepN]/float");
  mAllData->Branch("genLepPy",genLepPy,"genLepPy[genLepN]/float");
  mAllData->Branch("genLepPz",genLepPz,"genLepPz[genLepN]/float");
  mAllData->Branch("genLepStatus",genLepStatus,"genLepStatus[genLepN]/int");

  //mAllData->Branch("AlpPtScale" ,&mTempAlpPtScale,"mTempAlpPtScale/double");
  //mAllData->Branch("AlpIdTest" ,&mTempAlpIdTest ,"AlpIdTest/int");
  
  // MPT Markus 
  mAllData->Branch("MPTPhi" ,& mTempTreeMPTPhi ,"MPTPhi/double");
  mAllData->Branch("MPTPx" ,& mTempTreeMPTPx ,"MPTPx/double");
  mAllData->Branch("MPTPy" ,& mTempTreeMPTPy ,"MPTPy/double");
  mAllData->Branch("MPTPz" ,& mTempTreeMPTPz ,"MPTPz/double");
    
  edm::LogInfo("JaredSusyEvent") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JaredSusyAnalyzer);
