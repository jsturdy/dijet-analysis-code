
// -*- C++ -*-
//
// Package:    DiJetAnalysis
// Class:      JetAnalyzer
// 
/**\class JetAnalyzer JetAnalyzer.cc JSturdy/DiJetAnalysis/src/JetAnalyzer.cc

Description: Collects variables related to jets, performs dijet preselection
             Energy of jets = (50,50,30...), |eta|<2.5, 0.05<emfrac<0.95
             If successful, it stores the variables and returns the value of the check

*/
//
// Original Author:  Jared Sturdy
//         Created:  Fri Jan 29 16:10:31 PDT 2010
// $Id: JetAnalyzer.cpp,v 1.1 2010/01/29 16:10:31 sturdy Exp $
//
//

#include "JSturdy/DiJetAnalysis/interface/JetAnalyzer.h"

#include <TMath.h>
#include <sstream>
using namespace std;
using namespace reco;
using namespace edm;

//________________________________________________________________________________________
JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig)
{ 

  //defaults
  usePfJets_   = true;
  useJPTJets_  = true;
  useCaloJets_ = true;
  doMCData_    = true;

  jetMaxEta_ = 3.0;
  jetMinPt_ = 30;
  doMCData_ = true;
  // Preselection parameters
  if (iConfig.exists("jetMaxEta")) jetMaxEta_ = iConfig.getParameter<double>("jetMaxEta");
  if (iConfig.exists("jetMinPt"))  jetMinPt_  = iConfig.getParameter<double>("jetMinPt");
  if (iConfig.exists("doMCData"))  doMCData_  = iConfig.getParameter<bool>("doMCData");
  if (doMCData_) 
    if (iConfig.exists("genTag"))    genTag_    = iConfig.getParameter<edm::InputTag>("genTag");
 

  // get the data tags
  usePfJets_   = iConfig.getParameter<bool>("UsePfjet");
  useJPTJets_  = iConfig.getParameter<bool>("UseJPTjet");
  useCaloJets_ = iConfig.getParameter<bool>("UseCalojet");
  if (usePfJets_)   pfjetTag_  = iConfig.getParameter<edm::InputTag>("pfjetTag");
  if (useJPTJets_)  jptTag_    = iConfig.getParameter<edm::InputTag>("jptTag");
  if (useCaloJets_) jetTag_    = iConfig.getParameter<edm::InputTag>("jetTag");
  //If possible to retrieve these values rather than calculate them myself
  htTag_     = iConfig.getParameter<edm::InputTag>("htTag");
  mhtTag_    = iConfig.getParameter<edm::InputTag>("mhtTag");


  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
  initTuple();
  //bookHealthPlots();
}


//________________________________________________________________________________________
JetAnalyzer::~JetAnalyzer() {}


//________________________________________________________________________________________
// Method called to for each event
bool
JetAnalyzer::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  bool jetPreselection = false;
  bool jet_result = true;
  edm::LogVerbatim("DiJetEvent::JetAnalyzer") << " Start  " << std::endl;

  /*****************************
   * PFJets
   * Store the information on the particle flow jets
   * 
   *****************************/
  if (usePfJets_) {
    edm::Handle< std::vector<pat::Jet> > pfjetHandle;
    iEvent.getByLabel(pfjetTag_, pfjetHandle);
    if(!pfjetHandle.isValid()){
      edm::LogWarning("DiJetEvent::JetAnalyzer") << "No PFJet results for InputTag "<< pfjetTag_;
      return false;
    }

    int i = 0;
    double pfsumpx = 0;
    double pfsumpy = 0;
    double pfsumpt = 0;
        
    m_NPFJets = pfjetHandle->size();
    if( m_NPFJets > 50) m_NPFJets = 50;
    for(int pf=0; pf < m_NPFJets; pf++){
      if ((*pfjetHandle)[pf].pt() > jetMinPt_ ) {
	if (fabs((*pfjetHandle)[pf].eta()) < jetMaxEta_) {
	  m_PFJetEta[pf]    = (*pfjetHandle)[pf].eta();
	  m_PFJetPhi[pf]    = (*pfjetHandle)[pf].phi();
	  m_PFJetE[pf]      = (*pfjetHandle)[pf].energy();
	  m_PFJetEt[pf]     = (*pfjetHandle)[pf].et();
	  m_PFJetPx[pf]     = (*pfjetHandle)[pf].px();
	  m_PFJetPy[pf]     = (*pfjetHandle)[pf].py();
	  m_PFJetPz[pf]     = (*pfjetHandle)[pf].pz();
	  m_PFJetPt[pf]     = (*pfjetHandle)[pf].pt();
	  m_PFJetCharge[pf] = (*pfjetHandle)[pf].charge();
	  
	  if ((*pfjetHandle)[pf].isCaloJet())
	    m_PFJetFem[pf] = (*pfjetHandle)[pf].emEnergyFraction();
	  if ((*pfjetHandle)[pf].isPFJet())
	    m_PFJetFem[pf] = (*pfjetHandle)[pf].neutralEmEnergyFraction()+
	      (*pfjetHandle)[pf].chargedEmEnergyFraction();
	  
	  //values that contribute to HT and MHT
	  pfsumpt += (*pfjetHandle)[pf].pt();
	  pfsumpx += (*pfjetHandle)[pf].px();
	  pfsumpy += (*pfjetHandle)[pf].py();
	  //Increment for jets passing preselection
	  i++;
	}
      }
    }
    m_PFHt   = pfsumpt;
    m_PFMHx  = -pfsumpx;
    m_PFMHy  = -pfsumpy;
    m_PFMHt  = -sqrt(pfsumpx*pfsumpx+pfsumpy*pfsumpy);

    m_NPFJets = i;
  }
  else m_NPFJets = 0;

  /*****************************
   * Calo Jets
   * Store the jet information based on 
   * AK5 (AK7,SC5,SC7,IC5,SC7) calo jets
   * 
   *****************************/
  if (useCaloJets_) {
    edm::Handle< std::vector<pat::Jet> > jetHandle;
    iEvent.getByLabel(jetTag_, jetHandle);
    if ( !jetHandle.isValid() ) {
      edm::LogWarning("DiJetEvent::JetAnalyzer") << "No Jet results for InputTag " << jetTag_;
      return false;
    }

    // get the JPT-corrected pat::Jets
    edm::Handle< std::vector<pat::Jet> > jptHandle;
    iEvent.getByLabel(jptTag_, jptHandle);
    if ( !jptHandle.isValid() ) {
      edm::LogWarning("DiJetEvent::JetAnalyzer") << "No JetCorrFactor results for InputTag " << jptTag_;
      std::cout << "No JetCorrFactor results for InputTag " << jptTag_ << std::endl;
      return false;
    }

    //get number of jets
    m_NCaloJets = jetHandle->size();

    // Add the jets
    int i = 0;
    double jetsumpx = 0;
    double jetsumpy = 0;
    double jetsumpt = 0;

    double jptsumpx = 0;
    double jptsumpy = 0;
    double jptsumpt = 0;

    double gensumpx = 0;
    double gensumpy = 0;
    double gensumpt = 0;

    if ( m_NCaloJets >50 ) m_NCaloJets = 50;
    for (int k=0;k<m_NCaloJets;k++){
      const pat::Jet& uncorrJet = ((*jetHandle)[k].isCaloJet())? (*jetHandle)[k].correctedJet("RAW"): (*jetHandle)[k];
    
      if ((*jetHandle)[k].pt() > jetMinPt_) {
	if (fabs((*jetHandle)[k].eta()) < jetMaxEta_) {
	
	  jetsumpt += (*jetHandle)[k].pt();
	  jetsumpx += (*jetHandle)[k].momentum().X();
	  jetsumpy += (*jetHandle)[k].momentum().Y();
	
	  if((*jetHandle)[k].genJet()!= 0) {
	    gensumpt += (*jetHandle)[k].genJet()->pt();
	    gensumpx += (*jetHandle)[k].genJet()->momentum().X();
	    gensumpy += (*jetHandle)[k].genJet()->momentum().Y();}
	
	
	  for ( uint16_t n = 0; n < ( jptHandle->size() > 50 ? 50 : jptHandle->size() ); n++ ) {
	    if ( matchJetsByCaloTowers( (*jptHandle)[n], (*jetHandle)[k] ) ) {
	      pat::Jet jpt( (*jptHandle)[n] ); // no corrections by default
	      m_CaloJetJPTCorrFactor[i] = (jpt.isCaloJet()) ? ( jpt.energy() / uncorrJet.energy() ) : -1 ;
	    }
	  }
	
	  const reco::TrackRefVector & mrTracksInJet = (*jetHandle)[k].associatedTracks();
	
	  m_JetTrackPt[k]          = 0;
	  m_JetTrackPhi[k]         = 0;
	  m_JetTrackPhiWeighted[k] = 0;
	  m_JetTrackNo[k]          = 0;
	
	  float JetPhi = (*jetHandle)[k].phi();
	
	  for (reco::TrackRefVector::iterator aIter = mrTracksInJet.begin();aIter!= mrTracksInJet.end();aIter++)
	    {
	      m_JetTrackPt[k] += (*aIter)->pt();
	      float myPhi = (*aIter)->phi();
	      if( JetPhi > 2. ) {
		if(myPhi<0) myPhi = myPhi + 2*TMath::Pi();
	      }
	      if( JetPhi < -2. ) {
		if(myPhi>0) myPhi = myPhi - 2*TMath::Pi();
	      }
	      m_JetTrackPhiWeighted[k] += (*aIter)->pt()*myPhi;
	      m_JetTrackPhi[k]         += myPhi;
	      m_JetTrackNo[k]++;
	    
	    }
	
	  m_JetTrackPhiWeighted[k] = m_JetTrackPhiWeighted[k]/ m_JetTrackPt[k];
	  m_JetTrackPhi[k]         = m_JetTrackPhi[k]/float(m_JetTrackNo[k]);
	
	  m_CaloJetPt[i]  = (*jetHandle)[k].pt();
	  m_CaloJetE[i]   = (*jetHandle)[k].energy();
	  m_CaloJetEt[i]  = (*jetHandle)[k].et();
	  m_CaloJetPx[i]  = (*jetHandle)[k].momentum().X();
	  m_CaloJetPy[i]  = (*jetHandle)[k].momentum().Y();
	  m_CaloJetPz[i]  = (*jetHandle)[k].momentum().Z();
	  m_CaloJetEta[i] = (*jetHandle)[k].eta();
	  m_CaloJetPhi[i] = (*jetHandle)[k].phi();
	
	  if ((*jetHandle)[k].isCaloJet())
	    m_CaloJetFem[i] = (*jetHandle)[k].emEnergyFraction();
	  if ((*jetHandle)[k].isPFJet())
	    m_CaloJetFem[i] = (*jetHandle)[k].neutralEmEnergyFraction()+
	      (*jetHandle)[k].chargedEmEnergyFraction();
	
	  //m_JetBTag_TkCountHighEff[i] = (*jetHandle)[k].bDiscriminator("trackCountingHighEffBJetTags");
	  //m_JetBTag_SimpleSecVtx[i]   = (*jetHandle)[k].bDiscriminator("simpleSecondaryVertexBJetTags");
	  //m_JetBTag_CombSecVtx[i]     = (*jetHandle)[k].bDiscriminator("combinedSecondaryVertexBJetTags");
	  m_JetPartonFlavour[i]   = (*jetHandle)[k].partonFlavour();
	
	  if((*jetHandle)[k].genJet()!= 0) {
	    m_GenJetPt[i]  = (*jetHandle)[k].genJet()->pt();
	    m_GenJetE[i]   = (*jetHandle)[k].genJet()->energy();
	    m_GenJetEt[i]  = (*jetHandle)[k].genJet()->et();
	    m_GenJetPx[i]  = (*jetHandle)[k].genJet()->momentum().X();
	    m_GenJetPy[i]  = (*jetHandle)[k].genJet()->momentum().Y();
	    m_GenJetPz[i]  = (*jetHandle)[k].genJet()->momentum().z();
	    m_GenJetEta[i] = (*jetHandle)[k].genJet()->eta();
	    m_GenJetPhi[i] = (*jetHandle)[k].genJet()->phi();
	  }
	  else {
	    m_GenJetPt[i]  = -999;
	    m_GenJetE[i]   = -999;
	    m_GenJetEt[i]  = -999;
	    m_GenJetPx[i]  = -999;
	    m_GenJetPy[i]  = -999;
	    m_GenJetPz[i]  = -999;
	    m_GenJetEta[i] = -999;
	    m_GenJetPhi[i] = -999;
	  }
	  
	  if((*jetHandle)[k].genParton() != 0){
	    m_JetPartonId[i]     = (*jetHandle)[k].genParton()->pdgId();
	    m_JetPartonPx[i]     = (*jetHandle)[k].genParton()->px();
	    m_JetPartonPy[i]     = (*jetHandle)[k].genParton()->py();
	    m_JetPartonPz[i]     = (*jetHandle)[k].genParton()->pz();
	    m_JetPartonEt[i]     = (*jetHandle)[k].genParton()->et();
	    m_JetPartonPhi[i]    = (*jetHandle)[k].genParton()->phi();
	    m_JetPartonEta[i]    = (*jetHandle)[k].genParton()->eta();
	    m_JetPartonEnergy[i] = (*jetHandle)[k].genParton()->energy();
	    m_JetPartonMother[i] = (*jetHandle)[k].genParton()->mother()->pdgId();
	  }
	  else{
	    m_JetPartonId[i]     = -999;
	    m_JetPartonPx[i]     = -999;
	    m_JetPartonPy[i]     = -999;
	    m_JetPartonPz[i]     = -999;
	    m_JetPartonEt[i]     = -999;
	    m_JetPartonPhi[i]    = -999;
	    m_JetPartonEta[i]    = -999;
	    m_JetPartonEnergy[i] = -999;
	    m_JetPartonMother[i] = -999;
	  }
	
	  // Add the JPT corrs
	  int m_NJPTJets = jptHandle->size();
	  if ( m_NJPTJets > 50 ) m_NJPTJets = 50;
	  for ( int m = 0; m < m_NJPTJets; m++ ) {
	    if( (*jptHandle)[m].originalObjectRef() == (*jetHandle)[k].originalObjectRef() ) {
	      pat::Jet jet = ((*jptHandle)[m].isCaloJet()) ? (*jptHandle)[m].correctedJet("RAW") : (*jptHandle)[m];
	    }
	  }
	  i++;
	}
      }
    }
  
    m_NCaloJets  = i;
    m_CaloHt     = jetsumpt;
    m_CaloMHx    = -jetsumpx;
    m_CaloMHy    = -jetsumpy;
    m_CaloMHt    = -sqrt(jetsumpx*jetsumpx+jetsumpy*jetsumpy);
    
    m_JPTHt   = jptsumpt;
    m_JPTMHx  = -jptsumpx;
    m_JPTMHy  = -jptsumpy;
    m_JPTMHt  = -sqrt(jptsumpx*jptsumpx+jptsumpy*jptsumpy);
    
    m_GenHt  = gensumpt;
    m_GenMHx = -gensumpx;
    m_GenMHy = -gensumpy;
    m_GenMHt = -sqrt(gensumpx*gensumpx+gensumpy*gensumpy);
  }

  return jet_result;
  //mJetData->Fill();
}

//________________________________________________________________________________________
bool 
JetAnalyzer::matchJetsByCaloTowers( const pat::Jet& jet1,
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
JetAnalyzer::beginJob(const edm::EventSetup&) {}

//________________________________________________________________________________________
void 
JetAnalyzer::endJob() {

}

//________________________________________________________________________________________
void
JetAnalyzer::initTuple() {

  std::ostringstream variables; // Container for all variables
  
  // Register this ntuple
  edm::Service<TFileService> fs;

  mJetData = fs->make<TTree>( "JetData", "data after cuts" );
  mJetData->SetAutoSave(10);

  //add uncorrected non-cc particle flow jets
  mJetData->Branch("NPFJets",   &m_NPFJets,   "NPFJets/int");  
  mJetData->Branch("PFHt",      &m_PFHt,      "PFHt/double");
  mJetData->Branch("PFMHx",     &m_PFMHx,     "PFMHx/double");
  mJetData->Branch("PFMHy",     &m_PFMHy,     "PFMHy/double");
  mJetData->Branch("PFMHt",     &m_PFMHt,     "PFMHt/double");

  mJetData->Branch("PFJetE",       m_PFJetE,       "PFJetE[NPFJets]/double");
  mJetData->Branch("PFJetEt",      m_PFJetEt,      "PFJetEt[NPFJets]/double");
  mJetData->Branch("PFJetpt",      m_PFJetPt,      "PFJetpt[NPFJets]/double");
  mJetData->Branch("PFJetpx",      m_PFJetPx,      "PFJetpx[NPFJets]/double");
  mJetData->Branch("PFJetpy",      m_PFJetPy,      "PFJetpy[NPFJets]/double");
  mJetData->Branch("PFJetpz",      m_PFJetPz,      "PFJetpz[NPFJets]/double");
  mJetData->Branch("PFJeteta",     m_PFJetEta,     "PFJeteta[NPFJets]/double");
  mJetData->Branch("PFJetphi",     m_PFJetPhi,     "PFJetphi[NPFJets]/double");
  mJetData->Branch("PFJetFem",     m_PFJetFem,     "PFJetFem[NPFJets]/double");
  mJetData->Branch("PFJetCharge",  m_PFJetCharge,  "PFJetCharge[NPFJets]/double");
  //mJetData->Branch("PFJetHemi",    m_PFJetHemi,    "PFJetHemi[NPFJets]/int");


  //add uncorrected non-cc jets
  mJetData->Branch("NCaloJets",   &m_NCaloJets,   "NCaloJets/int");  
  mJetData->Branch("CaloHt",      &m_CaloHt,      "CaloHt/double");
  mJetData->Branch("CaloMHx",     &m_CaloMHx,     "CaloMHx/double");
  mJetData->Branch("CaloMHy",     &m_CaloMHy,     "CaloMHy/double");
  mJetData->Branch("CaloMHt",     &m_CaloMHt,     "CaloMHt/double");

  mJetData->Branch("CaloJetE",    m_CaloJetE,    "CaloJetE[NCaloJets]/double");
  mJetData->Branch("CaloJetEt",   m_CaloJetEt,   "CaloJetEt[NCaloJets]/double");
  mJetData->Branch("CaloJetpt",   m_CaloJetPt,   "CaloJetpt[NCaloJets]/double");
  mJetData->Branch("CaloJetpx",   m_CaloJetPx,   "CaloJetpx[NCaloJets]/double");
  mJetData->Branch("CaloJetpy",   m_CaloJetPy,   "CaloJetpy[NCaloJets]/double");
  mJetData->Branch("CaloJetpz",   m_CaloJetPz,   "CaloJetpz[NCaloJets]/double");
  mJetData->Branch("CaloJeteta",  m_CaloJetEta,  "CaloJeteta[NCaloJets]/double");
  mJetData->Branch("CaloJetphi",  m_CaloJetPhi,  "CaloJetphi[NCaloJets]/double");
  mJetData->Branch("CaloJetFem",  m_CaloJetFem,  "CaloJetFem[NCaloJets]/double");
  //mJetData->Branch("CaloJetHemi", m_CaloJetHemi, "CaloJetHemi[NCaloJets]/int");

  //jet correction factors
  //mJetData->Branch("CaloJetMCcorrFactor",  m_CaloJetMCCorrFactor,  "m_CaloJetMCCorrFactor[NCaloJETs]/double");
  mJetData->Branch("CaloJetJPTcorrFactor", m_CaloJetJPTCorrFactor, "m_CaloJetJPTCorrFactor[NCaloJETs]/double");
  //b-tagging information
  //mJetData->Branch("JetBTag_TkCountHighEff", m_JetBTag_TkCountHighEff, "JetBTag_TkCountHighEff[Njets]/float");
  //mJetData->Branch("JetBTag_SimpleSecVtx",   m_JetBTag_SimpleSecVtx,   "JetBTag_SimpleSecVtx[Njets]/float");
  //mJetData->Branch("JetBTag_CombSecVtx",     m_JetBTag_CombSecVtx,     "JetBTag_CombSecVtx[Njets]/float");

  //information about associated gen jets
  mJetData->Branch("GenHt",     &m_GenHt,     "GenHt/double");
  mJetData->Branch("GenMHt",    &m_GenMHt,    "GenMHt/double");
  mJetData->Branch("GenJetE" ,  m_GenJetE,   "GenJetE[NCaloJets]/double");
  mJetData->Branch("GenJetEt",  m_GenJetEt,  "GenJetEt[NCaloJets]/double");
  mJetData->Branch("GenJetpt",  m_GenJetPt,  "GenJetpt[NCaloJets]/double");
  mJetData->Branch("GenJetpx",  m_GenJetPx,  "GenJetpx[NCaloJets]/double");
  mJetData->Branch("GenJetpy",  m_GenJetPy,  "GenJetpy[NCaloJets]/double");
  mJetData->Branch("GenJetpz",  m_GenJetPz,  "GenJetpz[NCaloJets]/double");
  mJetData->Branch("GenJeteta", m_GenJetEta, "GenJeteta[NCaloJets]/double");
  mJetData->Branch("GenJetphi", m_GenJetPhi, "GenJetphi[NCaloJets]/double");
  
  //information about associated partons
  mJetData->Branch("JetPartonId",         m_JetPartonId,         "JetPartonId[NCaloJets]/int"); 
  mJetData->Branch("JetPartonMother",     m_JetPartonMother,     "JetPartonMother[NCaloJets]/int"); 
  mJetData->Branch("JetPartonPx",         m_JetPartonPx,         "JetPartonPx[NCaloJets]/double"); 
  mJetData->Branch("JetPartonPy",         m_JetPartonPy,         "JetPartonPy[NCaloJets]/double"); 
  mJetData->Branch("JetPartonPz",         m_JetPartonPz,         "JetPartonPz[NCaloJets]/double"); 
  mJetData->Branch("JetPartonEt",         m_JetPartonEt,         "JetPartonEt[NCaloJets]/double"); 
  mJetData->Branch("JetPartonE" ,         m_JetPartonEnergy,     "JetPartonE[NCaloJets]/double"); 
  mJetData->Branch("JetPartonPhi",        m_JetPartonPhi,        "JetPartonPhi[NCaloJets]/double"); 
  mJetData->Branch("JetPartonEta",        m_JetPartonEta,        "JetPartonEta[NCaloJets]/double"); 
  mJetData->Branch("JetPartonFlavour",    m_JetPartonFlavour,    "JetPartonFlavour[NCaloJets]/int");

  //information about associated tracks
  mJetData->Branch("JetTrackPt",          m_JetTrackPt,          "JetTrackPt[NCaloJets]/double"); 
  mJetData->Branch("JetTrackPhi",         m_JetTrackPhi,         "JetTrackPhi[NCaloJets]/double"); 
  mJetData->Branch("JetTrackPhiWeighted", m_JetTrackPhiWeighted, "JetTrackPhiWeighted[NCaloJets]/double"); 
  mJetData->Branch("JetTrackNo",          m_JetTrackNo,          "JetTrackNo[NCaloJets]/int");


  ////add JPT corrected calo jets
  //mJetData->Branch("NJPTjets",   &m_NJPTjets,   "NJPTjets/int");  
  //mJetData->Branch("JPTHt",      &m_JPTHt,      "JPTHt/double");
  //mJetData->Branch("JPTMHx",     &m_JPTMHx,     "JPTMHx/double");
  //mJetData->Branch("JPTMHy",     &m_JPTMHy,     "JPTMHy/double");
  //mJetData->Branch("JPTMHt",     &m_JPTMHt,     "JPTMHt/double");
  //
  //mJetData->Branch("JPTJetE",    m_JPTJetE,    "JPTJetE[NJPTJets]/double");
  //mJetData->Branch("JPTJetEt",   m_JPTJetEt,   "JPTJetEt[NJPTJets]/double");
  //mJetData->Branch("JPTJetpt",   m_JPTJetPt,   "JPTJetpt[NJPTJets]/double");
  //mJetData->Branch("JPTJetpx",   m_JPTJetPx,   "JPTJetpx[NJPTJets]/double");
  //mJetData->Branch("JPTJetpy",   m_JPTJetPy,   "JPTJetpy[NJPTJets]/double");
  //mJetData->Branch("JPTJetpz",   m_JPTJetPz,   "JPTJetpz[NJPTJets]/double");
  //mJetData->Branch("JPTJeteta",  m_JPTJetEta,  "JPTJeteta[NJPTJets]/double");
  //mJetData->Branch("JPTJetphi",  m_JPTJetPhi,  "JPTJetphi[NJPTJets]/double");
  //mJetData->Branch("JPTJetFem",  m_JPTJetFem,  "JPTJetFem[NJPTJets]/double");
  //mJetData->Branch("JPTJetHemi", m_JPTJetHemi, "JPTJetHemi[NJPTJets]/int");

  //jet correction factors
  //mJetData->Branch("JPTJetMCcorrFactor",  m_JPTJetMCCorrFactor,  "m_JPTJetMCCorrFactor[NJPTJets]/double");
  //mJetData->Branch("JPTJetJPTcorrFactor", m_JPTJetJPTCorrFactor, "m_JPTJetJPTCorrFactor[NJPTJets]/double");


  
  edm::LogInfo("DiJetEvent::JetAnalyzer") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_EDM_PLUGIN(JetAnalyzer);
