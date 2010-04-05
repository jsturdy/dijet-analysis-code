
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
// $Id: JetAnalyzer.cc,v 1.3 2010/04/04 00:04:01 sturdy Exp $
//
//

#include "JSturdy/DiJetAnalysis/interface/JetAnalyzer.h"

#include <TMath.h>
#include <sstream>

//________________________________________________________________________________________
JetAnalyzer::JetAnalyzer(const edm::ParameterSet& pset, TTree* tmpAllData)
{ 
  mJetData = tmpAllData;
  jetParams = pset;
  //defaults
  usePFJets_    = true;
  useJPTJets_   = true;
  useCaloJets_  = true;
  useTrackJets_ = true;
  doMCData_     = false;
  debug_        = 0;
  // Preselection parameters
  //Basic jet information, minimum number, eta and pt requirement for all jets
  minNJets_  = 1;
  jetMaxEta_ = 5.0;
  jetMinPt_  = 30.;
  jetMaxEMF_ = 0.99;
  jetMinEMF_ = 0.0;
  if (jetParams.exists("debugJets"))     debug_   = jetParams.getUntrackedParameter<int>("debugJets");
  //if (jetParams.exists("minNJets"))  minNJets_  = jetParams.getUntrackedParameter<int>("minNJets");
  if (jetParams.exists("jetMaxEta")) jetMaxEta_ = jetParams.getUntrackedParameter<double >("jetMaxEta");
  if (jetParams.exists("jetMinPt"))  jetMinPt_  = jetParams.getUntrackedParameter<double >("jetMinPt");
  if (jetParams.exists("jetMaxEMF")) jetMaxEMF_ = jetParams.getUntrackedParameter<double >("jetMaxEMF");
  if (jetParams.exists("jetMinEMF")) jetMinEMF_ = jetParams.getUntrackedParameter<double >("jetMinEMF");
  //for (int nj = 0; nj < selJetMaxEta_.size(); ++nj) {
  //  printf("jet %2d, eta max %2.2f, min pt %2.2f, max emf %2.2f, min emf %2.2f\n",nj,selJetMaxEta_.at(nj), selJetMinPt_.at(nj),selJetMaxEMF_.at(nj), selJetMinEMF_.at(nj));
  //}

  //Individual jet requirements
  selJetMaxEta_ = jetParams.getUntrackedParameter<std::vector<double > >("selJetMaxEta");
  selJetMinPt_  = jetParams.getUntrackedParameter<std::vector<double > >("selJetMinPt");
  selJetMaxEMF_ = jetParams.getUntrackedParameter<std::vector<double > >("selJetMaxEMF");
  selJetMinEMF_ = jetParams.getUntrackedParameter<std::vector<double > >("selJetMinEMF");
  if (debug_) {
    std::cout<<"size of dijet vector "<<selJetMaxEta_.size()<<std::endl;
    for (int nj = 0; nj < selJetMaxEta_.size(); ++nj) {
      printf("jet %2d, eta max %2.2f, min pt %2.2f, max emf %2.2f, min emf %2.2f\n",nj,selJetMaxEta_.at(nj), selJetMinPt_.at(nj),selJetMaxEMF_.at(nj), selJetMinEMF_.at(nj));
    }
  }
  
  doMCData_ = true;
  if (jetParams.exists("doMCJets"))     doMCData_     = jetParams.getUntrackedParameter<bool>("doMCJets");
  if (doMCData_) 
    if (jetParams.exists("genJetTag"))    genJetTag_    = jetParams.getUntrackedParameter<edm::InputTag>("genJetTag");
 
  // get the data tags
  usePFJets_    = jetParams.getUntrackedParameter<bool>("usePFJets");
  useJPTJets_   = jetParams.getUntrackedParameter<bool>("useJPTJets");
  useCaloJets_  = jetParams.getUntrackedParameter<bool>("useCaloJets");
  useTrackJets_ = jetParams.getUntrackedParameter<bool>("useTrackJets");
  if (usePFJets_)    pfJetTag_    = jetParams.getUntrackedParameter<edm::InputTag>("pfJetTag");
  if (useJPTJets_)   jptJetTag_   = jetParams.getUntrackedParameter<edm::InputTag>("jptJetTag");
  if (useCaloJets_)  caloJetTag_  = jetParams.getUntrackedParameter<edm::InputTag>("caloJetTag");
  if (useTrackJets_) trackJetTag_ = jetParams.getUntrackedParameter<edm::InputTag>("trackJetTag");
  //If possible to retrieve these values rather than calculate them myself
  //htTag_     = jetParams.getUntrackedParameter<edm::InputTag>("htTag");
  mhtTag_    = jetParams.getUntrackedParameter<edm::InputTag>("mhtTag");


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
  using namespace reco;
  using namespace edm;

  m_pfJetPreselection    = false;
  m_jptJetPreselection   = false;
  m_caloJetPreselection  = false;
  m_trackJetPreselection = false;
  bool jet_result = true;
  edm::LogVerbatim("DiJetEvent::JetAnalyzer") << " Start  " << std::endl;

  /*****************************
   * PFJets
   * Store the information on the particle flow jets
   * 
   *****************************/
  if (usePFJets_) {
    edm::Handle< std::vector<pat::Jet> > pfjetHandle;
    iEvent.getByLabel(pfJetTag_, pfjetHandle);
    if(!pfjetHandle.isValid()){
      edm::LogWarning("DiJetEvent::JetAnalyzer") << "No PFJet results for InputTag "<< pfJetTag_;
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

  //determine preselection requirement based on pf jets
  if (m_NPFJets>selJetMaxEta_.size()) {
    m_pfJetPreselection = true;
    for (unsigned int pfj = 0; pfj < selJetMaxEta_.size(); ++pfj) {
      if (m_PFJetEta[pfj] > selJetMaxEta_.at(pfj)) m_pfJetPreselection = false;
      if (m_PFJetPt[pfj]  < selJetMinPt_.at(pfj))  m_pfJetPreselection = false;
      if (m_PFJetFem[pfj] > selJetMaxEMF_.at(pfj)) m_pfJetPreselection = false;
      if (m_PFJetFem[pfj] < selJetMinEMF_.at(pfj)) m_pfJetPreselection = false;
    }
  }

  /*****************************
   * Track Jets
   * Store the jet information based on 
   * AK5 (AK7,SC5,SC7,IC5,SC7) track jets
   * 
   *****************************/
  if (useTrackJets_) {
    edm::Handle< std::vector<pat::Jet> > trackjetHandle;
    iEvent.getByLabel(trackJetTag_, trackjetHandle);
    if(!trackjetHandle.isValid()){
      edm::LogWarning("DiJetEvent::JetAnalyzer") << "No TrackJet results for InputTag "<< trackJetTag_;
      return false;
    }

    int i = 0;
    double tracksumpx = 0;
    double tracksumpy = 0;
    double tracksumpt = 0;
        
    m_NTrackJets = trackjetHandle->size();
    if( m_NTrackJets > 50) m_NTrackJets = 50;
    for(int track=0; track < m_NTrackJets; track++){
      if ((*trackjetHandle)[track].pt() > jetMinPt_ ) {
	if (fabs((*trackjetHandle)[track].eta()) < jetMaxEta_) {
	  m_TrackJetEta[track]    = (*trackjetHandle)[track].eta();
	  m_TrackJetPhi[track]    = (*trackjetHandle)[track].phi();
	  m_TrackJetE[track]      = (*trackjetHandle)[track].energy();
	  m_TrackJetEt[track]     = (*trackjetHandle)[track].et();
	  m_TrackJetPx[track]     = (*trackjetHandle)[track].px();
	  m_TrackJetPy[track]     = (*trackjetHandle)[track].py();
	  m_TrackJetPz[track]     = (*trackjetHandle)[track].pz();
	  m_TrackJetPt[track]     = (*trackjetHandle)[track].pt();
	  m_TrackJetCharge[track] = (*trackjetHandle)[track].charge();
	  
	  if ((*trackjetHandle)[track].isCaloJet())
	    m_TrackJetFem[track] = (*trackjetHandle)[track].emEnergyFraction();
	  if ((*trackjetHandle)[track].isPFJet())
	    m_TrackJetFem[track] = (*trackjetHandle)[track].neutralEmEnergyFraction()+
	      (*trackjetHandle)[track].chargedEmEnergyFraction();
	  
	  //values that contribute to HT and MHT
	  tracksumpt += (*trackjetHandle)[track].pt();
	  tracksumpx += (*trackjetHandle)[track].px();
	  tracksumpy += (*trackjetHandle)[track].py();
	  //Increment for jets passing preselection
	  i++;
	}
      }
    }
    m_TrackHt   = tracksumpt;
    m_TrackMHx  = -tracksumpx;
    m_TrackMHy  = -tracksumpy;
    m_TrackMHt  = -sqrt(tracksumpx*tracksumpx+tracksumpy*tracksumpy);

    m_NTrackJets = i;
  }
  else m_NTrackJets = 0;
  
  //determine preselection requirement based on pf jets
  if (m_NTrackJets>selJetMaxEta_.size()) {
    m_trackJetPreselection = true;
    for (unsigned int trackj = 0; trackj < selJetMaxEta_.size(); ++trackj) {
      if (m_TrackJetEta[trackj] > selJetMaxEta_.at(trackj)) m_trackJetPreselection = false;
      if (m_TrackJetPt[trackj]  < selJetMinPt_.at(trackj))  m_trackJetPreselection = false;
      if (m_TrackJetFem[trackj] > selJetMaxEMF_.at(trackj)) m_trackJetPreselection = false;
      if (m_TrackJetFem[trackj] < selJetMinEMF_.at(trackj)) m_trackJetPreselection = false;
    }
  }


  /*****************************
   * Calo Jets
   * Store the jet information based on 
   * AK5 (AK7,SC5,SC7,IC5,SC7) calo jets
   * 
   *****************************/
  if (useCaloJets_) {
    edm::Handle< std::vector<pat::Jet> > jetHandle;
    iEvent.getByLabel(caloJetTag_, jetHandle);
    if ( !jetHandle.isValid() ) {
      edm::LogWarning("DiJetEvent::JetAnalyzer") << "No Jet results for InputTag " << caloJetTag_;
      return false;
    }

    // get the JPT-corrected pat::Jets
    edm::Handle< std::vector<pat::Jet> > jptHandle;
    iEvent.getByLabel(jptJetTag_, jptHandle);
    if ( !jptHandle.isValid() ) {
      edm::LogWarning("DiJetEvent::JetAnalyzer") << "No JetCorrFactor results for InputTag " << jptJetTag_;
      std::cout << "No JetCorrFactor results for InputTag " << jptJetTag_ << std::endl;
      return false;
    }

    //get number of jets
    m_NCaloJets = jetHandle->size();

    // Add the jets
    int i = 0;
    double jetsumpx = 0;
    double jetsumpy = 0;
    double jetsumpt = 0;

    double gensumpx = 0;
    double gensumpy = 0;
    double gensumpt = 0;

    if ( m_NCaloJets >50 ) m_NCaloJets = 50;
    for (int k=0;k<m_NCaloJets;k++){
      const pat::Jet& uncorrJet = ((*jetHandle)[k].isCaloJet())? (*jetHandle)[k].correctedJet("RAW"): (*jetHandle)[k];
    
      if ((*jetHandle)[k].pt() > jetMinPt_) {
	if (fabs((*jetHandle)[k].eta()) < jetMaxEta_) {
	  //if ((*jetHandle)[k].emEnergyFraction() > jetMinEMF_) {
	  //if ((*jetHandle)[k].jetID().fHPD < jetMaxfHPD_) {
	  //if ((*jetHandle)[k].jetID().fRBX < jetMaxfRBX_) {
	  //if ((*jetHandle)[k].jetID().n90Hits > jetMinN90_) {
	
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
	  
	  m_CaloJetTrackPt[k]          = 0;
	  m_CaloJetTrackPhi[k]         = 0;
	  m_CaloJetTrackPhiWeighted[k] = 0;
	  m_CaloJetTrackNo[k]          = 0;
	
	  float JetPhi = (*jetHandle)[k].phi();
	
	  for (reco::TrackRefVector::iterator aIter = mrTracksInJet.begin();aIter!= mrTracksInJet.end();aIter++)
	    {
	      m_CaloJetTrackPt[k] += (*aIter)->pt();
	      float myPhi = (*aIter)->phi();
	      if( JetPhi > 2. ) {
		if(myPhi<0) myPhi = myPhi + 2*TMath::Pi();
	      }
	      if( JetPhi < -2. ) {
		if(myPhi>0) myPhi = myPhi - 2*TMath::Pi();
	      }
	      m_CaloJetTrackPhiWeighted[k] += (*aIter)->pt()*myPhi;
	      m_CaloJetTrackPhi[k]         += myPhi;
	      m_CaloJetTrackNo[k]++;
	    
	    }
	
	  m_CaloJetTrackPhiWeighted[k] = m_CaloJetTrackPhiWeighted[k]/ m_CaloJetTrackPt[k];
	  m_CaloJetTrackPhi[k]         = m_CaloJetTrackPhi[k]/float(m_CaloJetTrackNo[k]);
	
	  m_CaloJetPt[i]   = (*jetHandle)[k].pt();
	  m_CaloJetE[i]    = (*jetHandle)[k].energy();
	  m_CaloJetEt[i]   = (*jetHandle)[k].et();
	  m_CaloJetPx[i]   = (*jetHandle)[k].momentum().X();
	  m_CaloJetPy[i]   = (*jetHandle)[k].momentum().Y();
	  m_CaloJetPz[i]   = (*jetHandle)[k].momentum().Z();
	  m_CaloJetEta[i]  = (*jetHandle)[k].eta();
	  m_CaloJetPhi[i]  = (*jetHandle)[k].phi();
	  m_CaloJetN90[i]  = (*jetHandle)[k].jetID().n90Hits;
	  m_CaloJetfHPD[i] = (*jetHandle)[k].jetID().fHPD;
	  m_CaloJetfRBX[i] = (*jetHandle)[k].jetID().fRBX;
	
	  if ((*jetHandle)[k].isCaloJet())
	    m_CaloJetFem[i] = (*jetHandle)[k].emEnergyFraction();
	  if ((*jetHandle)[k].isPFJet())
	    m_CaloJetFem[i] = (*jetHandle)[k].neutralEmEnergyFraction()+
	      (*jetHandle)[k].chargedEmEnergyFraction();
	
	  //m_CaloJetBTag_TkCountHighEff[i] = (*jetHandle)[k].bDiscriminator("trackCountingHighEffBJetTags");
	  //m_CaloJetBTag_SimpleSecVtx[i]   = (*jetHandle)[k].bDiscriminator("simpleSecondaryVertexBJetTags");
	  //m_CaloJetBTag_CombSecVtx[i]     = (*jetHandle)[k].bDiscriminator("combinedSecondaryVertexBJetTags");
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
	  i++;
	}
      }
    }
  
    m_NCaloJets  = i;
    m_CaloHt     = jetsumpt;
    m_CaloMHx    = -jetsumpx;
    m_CaloMHy    = -jetsumpy;
    m_CaloMHt    = -sqrt(jetsumpx*jetsumpx+jetsumpy*jetsumpy);
    
    m_GenHt  = gensumpt;
    m_GenMHx = -gensumpx;
    m_GenMHy = -gensumpy;
    m_GenMHt = -sqrt(gensumpx*gensumpx+gensumpy*gensumpy);
  }

  //determine preselection requirement based on calo jets
  if (m_NCaloJets>selJetMaxEta_.size()) {
    m_caloJetPreselection = true;
    for (unsigned int caloj = 0; caloj < selJetMaxEta_.size(); ++caloj) {
      if (m_CaloJetEta[caloj] > selJetMaxEta_.at(caloj)) m_caloJetPreselection = false;
      if (m_CaloJetPt[caloj]  < selJetMinPt_.at(caloj))  m_caloJetPreselection = false;
      if (m_CaloJetFem[caloj] > selJetMaxEMF_.at(caloj)) m_caloJetPreselection = false;
      if (m_CaloJetFem[caloj] < selJetMinEMF_.at(caloj)) m_caloJetPreselection = false;
    }
  }


  /*****************************
   * JPT Jets
   * Store the jet information based on 
   * AK5 (AK7,SC5,SC7,IC5,SC7) JPT corrected calo jets
   * 
   *****************************/
  if (useJPTJets_) {
    edm::Handle< std::vector<pat::Jet> > jptHandle;
    iEvent.getByLabel(jptJetTag_, jptHandle);
    if ( !jptHandle.isValid() ) {
      edm::LogWarning("DiJetEvent::JetAnalyzer") << "No Jet results for InputTag " << caloJetTag_;
      return false;
    }
    
    //get number of jets
    m_NJPTJets = jptHandle->size();
    
    // Add the jets
    int i = 0;
    double jptsumpx = 0;
    double jptsumpy = 0;
    double jptsumpt = 0;
    
    if ( m_NJPTJets >50 ) m_NJPTJets = 50;
    for (int k=0;k<m_NJPTJets;k++){
      
      if ((*jptHandle)[k].pt() > jetMinPt_) {
	if (fabs((*jptHandle)[k].eta()) < jetMaxEta_) {
	  //if ((*jptHandle)[k].emEnergyFraction() > jetMinEMF_) {
	  //if ((*jptHandle)[k].jetID().fHPD < jetMaxfHPD_) {
	  //if ((*jptHandle)[k].jetID().fRBX < jetMaxfRBX_) {
	  //if ((*jptHandle)[k].jetID().n90Hits > jetMinN90_) {
	  
	  jptsumpt += (*jptHandle)[k].pt();
	  jptsumpx += (*jptHandle)[k].momentum().X();
	  jptsumpy += (*jptHandle)[k].momentum().Y();
	  
	  const reco::TrackRefVector & mrTracksInJet = (*jptHandle)[k].associatedTracks();
	  
	  m_JPTJetTrackPt[k]          = 0;
	  m_JPTJetTrackPhi[k]         = 0;
	  m_JPTJetTrackPhiWeighted[k] = 0;
	  m_JPTJetTrackNo[k]          = 0;
	
	  float JetPhi = (*jptHandle)[k].phi();
	
	  for (reco::TrackRefVector::iterator aIter = mrTracksInJet.begin();aIter!= mrTracksInJet.end();aIter++)
	    {
	      m_JPTJetTrackPt[k] += (*aIter)->pt();
	      float myPhi = (*aIter)->phi();
	      if( JetPhi > 2. ) {
		if(myPhi<0) myPhi = myPhi + 2*TMath::Pi();
	      }
	      if( JetPhi < -2. ) {
		if(myPhi>0) myPhi = myPhi - 2*TMath::Pi();
	      }
	      m_JPTJetTrackPhiWeighted[k] += (*aIter)->pt()*myPhi;
	      m_JPTJetTrackPhi[k]         += myPhi;
	      m_JPTJetTrackNo[k]++;
	    
	    }
	
	  m_JPTJetTrackPhiWeighted[k] = m_JPTJetTrackPhiWeighted[k]/ m_JPTJetTrackPt[k];
	  m_JPTJetTrackPhi[k]         = m_JPTJetTrackPhi[k]/float(m_JPTJetTrackNo[k]);
	  
	  
	  m_JPTJetPt[i]   = (*jptHandle)[k].pt();
	  m_JPTJetE[i]    = (*jptHandle)[k].energy();
	  m_JPTJetEt[i]   = (*jptHandle)[k].et();
	  m_JPTJetPx[i]   = (*jptHandle)[k].momentum().X();
	  m_JPTJetPy[i]   = (*jptHandle)[k].momentum().Y();
	  m_JPTJetPz[i]   = (*jptHandle)[k].momentum().Z();
	  m_JPTJetEta[i]  = (*jptHandle)[k].eta();
	  m_JPTJetPhi[i]  = (*jptHandle)[k].phi();
	  m_JPTJetN90[i]  = (*jptHandle)[k].jetID().n90Hits;
	  m_JPTJetfHPD[i] = (*jptHandle)[k].jetID().fHPD;
	  m_JPTJetfRBX[i] = (*jptHandle)[k].jetID().fRBX;
	  
	  if ((*jptHandle)[k].isCaloJet())
	    m_JPTJetFem[i] = (*jptHandle)[k].emEnergyFraction();
	  if ((*jptHandle)[k].isPFJet())
	    m_JPTJetFem[i] = (*jptHandle)[k].neutralEmEnergyFraction()+
	      (*jptHandle)[k].chargedEmEnergyFraction();
	  
	  //m_JPTJetBTag_TkCountHighEff[i] = (*jptHandle)[k].bDiscriminator("trackCountingHighEffBJetTags");
	  //m_JPTJetBTag_SimpleSecVtx[i]   = (*jptHandle)[k].bDiscriminator("simpleSecondaryVertexBJetTags");
	  //m_JPTJetBTag_CombSecVtx[i]     = (*jptHandle)[k].bDiscriminator("combinedSecondaryVertexBJetTags");
	  i++;
	}
      }
    }
    
    m_NJPTJets  = i;
    
    m_JPTHt   = jptsumpt;
    m_JPTMHx  = -jptsumpx;
    m_JPTMHy  = -jptsumpy;
    m_JPTMHt  = -sqrt(jptsumpx*jptsumpx+jptsumpy*jptsumpy);
  }
  
  //determine preselection requirement based on jpt calo jets
  if (m_NJPTJets>selJetMaxEta_.size()) {
    m_jptJetPreselection = true;
    for (unsigned int jptj = 0; jptj < selJetMaxEta_.size(); ++jptj) {
      if (m_JPTJetEta[jptj] > selJetMaxEta_.at(jptj)) m_jptJetPreselection = false;
      if (m_JPTJetPt[jptj]  < selJetMinPt_.at(jptj))  m_jptJetPreselection = false;
      if (m_JPTJetFem[jptj] > selJetMaxEMF_.at(jptj)) m_jptJetPreselection = false;
      if (m_JPTJetFem[jptj] < selJetMinEMF_.at(jptj)) m_jptJetPreselection = false;
    }
  }


  jet_result = m_pfJetPreselection || m_trackJetPreselection || m_caloJetPreselection || m_jptJetPreselection;
  //mJetData->Fill();
  return jet_result;
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
JetAnalyzer::beginJob() {}

//________________________________________________________________________________________
void 
JetAnalyzer::endJob() {

}

//________________________________________________________________________________________
void
JetAnalyzer::initTuple() {

  std::ostringstream variables; // Container for all variables
  
  /*
  // Register this ntuple
  edm::Service<TFileService> fs;

  mJetData = fs->make<TTree>( "JetData", "data after cuts" );
  mJetData->SetAutoSave(10);
  */

  if (usePFJets_) {
    //add uncorrected non-cc particle flow jets
    mJetData->Branch("PFNJets",   &m_NPFJets,   "PFNJets/int");  
    mJetData->Branch("PFHt",      &m_PFHt,      "PFHt/double");
    mJetData->Branch("PFMHx",     &m_PFMHx,     "PFMHx/double");
    mJetData->Branch("PFMHy",     &m_PFMHy,     "PFMHy/double");
    mJetData->Branch("PFMHt",     &m_PFMHt,     "PFMHt/double");
    
    mJetData->Branch("PFJetE",       m_PFJetE,       "PFJetE[PFNJets]/double");
    mJetData->Branch("PFJetEt",      m_PFJetEt,      "PFJetEt[PFNJets]/double");
    mJetData->Branch("PFJetpt",      m_PFJetPt,      "PFJetpt[PFNJets]/double");
    mJetData->Branch("PFJetpx",      m_PFJetPx,      "PFJetpx[PFNJets]/double");
    mJetData->Branch("PFJetpy",      m_PFJetPy,      "PFJetpy[PFNJets]/double");
    mJetData->Branch("PFJetpz",      m_PFJetPz,      "PFJetpz[PFNJets]/double");
    mJetData->Branch("PFJeteta",     m_PFJetEta,     "PFJeteta[PFNJets]/double");
    mJetData->Branch("PFJetphi",     m_PFJetPhi,     "PFJetphi[PFNJets]/double");
    mJetData->Branch("PFJetFem",     m_PFJetFem,     "PFJetFem[PFNJets]/double");
    mJetData->Branch("PFJetCharge",  m_PFJetCharge,  "PFJetCharge[PFNJets]/double");
    //mJetData->Branch("PFJetHemi",    m_PFJetHemi,    "PFJetHemi[PFNJets]/int");
    mJetData->Branch("PFJetPreselection", &m_pfJetPreselection, "PFJetPreselection/bool");
  }

  if (useTrackJets_) {
    //add uncorrected non-cc track jets
    mJetData->Branch("TrackNJets",   &m_NTrackJets,   "TrackNJets/int");  
    mJetData->Branch("TrackHt",      &m_TrackHt,      "TrackHt/double");
    mJetData->Branch("TrackMHx",     &m_TrackMHx,     "TrackMHx/double");
    mJetData->Branch("TrackMHy",     &m_TrackMHy,     "TrackMHy/double");
    mJetData->Branch("TrackMHt",     &m_TrackMHt,     "TrackMHt/double");
    
    mJetData->Branch("TrackJetE",       m_TrackJetE,       "TrackJetE[TrackNJets]/double");
    mJetData->Branch("TrackJetEt",      m_TrackJetEt,      "TrackJetEt[TrackNJets]/double");
    mJetData->Branch("TrackJetpt",      m_TrackJetPt,      "TrackJetpt[TrackNJets]/double");
    mJetData->Branch("TrackJetpx",      m_TrackJetPx,      "TrackJetpx[TrackNJets]/double");
    mJetData->Branch("TrackJetpy",      m_TrackJetPy,      "TrackJetpy[TrackNJets]/double");
    mJetData->Branch("TrackJetpz",      m_TrackJetPz,      "TrackJetpz[TrackNJets]/double");
    mJetData->Branch("TrackJeteta",     m_TrackJetEta,     "TrackJeteta[TrackNJets]/double");
    mJetData->Branch("TrackJetphi",     m_TrackJetPhi,     "TrackJetphi[TrackNJets]/double");
    mJetData->Branch("TrackJetFem",     m_TrackJetFem,     "TrackJetFem[TrackNJets]/double");
    mJetData->Branch("TrackJetCharge",  m_TrackJetCharge,  "TrackJetCharge[TrackNJets]/double");
    //mJetData->Branch("TrackJetHemi",    m_TrackJetHemi,    "TrackJetHemi[TrackNJets]/int");
    mJetData->Branch("TrackJetPreselection", &m_trackJetPreselection, "TrackJetPreselection/bool");
  }

  if (useCaloJets_) {
    //add uncorrected non-cc jets
    mJetData->Branch("CaloNJets",   &m_NCaloJets,   "CaloNJets/int");  
    mJetData->Branch("CaloHt",      &m_CaloHt,      "CaloHt/double");
    mJetData->Branch("CaloMHx",     &m_CaloMHx,     "CaloMHx/double");
    mJetData->Branch("CaloMHy",     &m_CaloMHy,     "CaloMHy/double");
    mJetData->Branch("CaloMHt",     &m_CaloMHt,     "CaloMHt/double");
    
    mJetData->Branch("CaloJetE",    m_CaloJetE,    "CaloJetE[CaloNJets]/double");
    mJetData->Branch("CaloJetEt",   m_CaloJetEt,   "CaloJetEt[CaloNJets]/double");
    mJetData->Branch("CaloJetpt",   m_CaloJetPt,   "CaloJetpt[CaloNJets]/double");
    mJetData->Branch("CaloJetpx",   m_CaloJetPx,   "CaloJetpx[CaloNJets]/double");
    mJetData->Branch("CaloJetpy",   m_CaloJetPy,   "CaloJetpy[CaloNJets]/double");
    mJetData->Branch("CaloJetpz",   m_CaloJetPz,   "CaloJetpz[CaloNJets]/double");
    mJetData->Branch("CaloJeteta",  m_CaloJetEta,  "CaloJeteta[CaloNJets]/double");
    mJetData->Branch("CaloJetphi",  m_CaloJetPhi,  "CaloJetphi[CaloNJets]/double");
    mJetData->Branch("CaloJetFem",  m_CaloJetFem,  "CaloJetFem[CaloNJets]/double");
    mJetData->Branch("CaloJetfHPD", m_CaloJetfHPD, "CaloJetfHPD[CaloNJets]/double");
    mJetData->Branch("CaloJetfRBX", m_CaloJetfRBX, "CaloJetfRBX[CaloNJets]/double");
    mJetData->Branch("CaloJetn90",  m_CaloJetN90,  "CaloJetn90[CaloNJets]/double");
    //mJetData->Branch("CaloJetHemi", m_CaloJetHemi, "CaloJetHemi[CaloNJets]/int");
    mJetData->Branch("CaloJetPreselection", &m_caloJetPreselection, "CaloJetPreselection/bool");

    //jet correction factors
    //mJetData->Branch("CaloJetMCcorrFactor",  m_CaloJetMCCorrFactor,  "m_CaloJetMCCorrFactor[JPTNJets]/double");
    mJetData->Branch("CaloJetJPTcorrFactor", m_CaloJetJPTCorrFactor, "m_CaloJetJPTCorrFactor[CaloNJets]/double");
    //b-tagging information
    //mJetData->Branch("CaloJetBTag_TkCountHighEff", m_CaloJetBTag_TkCountHighEff, "CaloJetBTag_TkCountHighEff[JPTNJetsJ/float");
    //mJetData->Branch("CaloJetBTag_SimpleSecVtx",   m_CaloJetBTag_SimpleSecVtx,   "CaloJetBTag_SimpleSecVtx[JPTNJets]/float");
    //mJetData->Branch("CaloJetBTag_CombSecVtx",     m_CaloJetBTag_CombSecVtx,     "CaloJetBTag_CombSecVtx[JPTNJets]/float");

    //information about associated gen jets
    mJetData->Branch("GenHt",     &m_GenHt,     "GenHt/double");
    mJetData->Branch("GenMHt",    &m_GenMHt,    "GenMHt/double");
    mJetData->Branch("GenJetE" ,  m_GenJetE,   "GenJetE[CaloNJets]/double");
    mJetData->Branch("GenJetEt",  m_GenJetEt,  "GenJetEt[CaloNJets]/double");
    mJetData->Branch("GenJetpt",  m_GenJetPt,  "GenJetpt[CaloNJets]/double");
    mJetData->Branch("GenJetpx",  m_GenJetPx,  "GenJetpx[CaloNJets]/double");
    mJetData->Branch("GenJetpy",  m_GenJetPy,  "GenJetpy[CaloNJets]/double");
    mJetData->Branch("GenJetpz",  m_GenJetPz,  "GenJetpz[CaloNJets]/double");
    mJetData->Branch("GenJeteta", m_GenJetEta, "GenJeteta[CaloNJets]/double");
    mJetData->Branch("GenJetphi", m_GenJetPhi, "GenJetphi[CaloNJets]/double");
    
    //information about associated partons
    mJetData->Branch("JetPartonId",         m_JetPartonId,         "JetPartonId[CaloNJets]/int"); 
    mJetData->Branch("JetPartonMother",     m_JetPartonMother,     "JetPartonMother[CaloNJets]/int"); 
    mJetData->Branch("JetPartonPx",         m_JetPartonPx,         "JetPartonPx[CaloNJets]/double"); 
    mJetData->Branch("JetPartonPy",         m_JetPartonPy,         "JetPartonPy[CaloNJets]/double"); 
    mJetData->Branch("JetPartonPz",         m_JetPartonPz,         "JetPartonPz[CaloNJets]/double"); 
    mJetData->Branch("JetPartonEt",         m_JetPartonEt,         "JetPartonEt[CaloNJets]/double"); 
    mJetData->Branch("JetPartonE" ,         m_JetPartonEnergy,     "JetPartonE[CaloNJets]/double"); 
    mJetData->Branch("JetPartonPhi",        m_JetPartonPhi,        "JetPartonPhi[CaloNJets]/double"); 
    mJetData->Branch("JetPartonEta",        m_JetPartonEta,        "JetPartonEta[CaloNJets]/double"); 
    mJetData->Branch("JetPartonFlavour",    m_JetPartonFlavour,    "JetPartonFlavour[CaloNJets]/int");
    
    //information about associated tracks
    mJetData->Branch("CaloJetTrackPt",          m_CaloJetTrackPt,          "CaloJetTrackPt[CaloNJets]/double"); 
    mJetData->Branch("CaloJetTrackPhi",         m_CaloJetTrackPhi,         "CaloJetTrackPhi[CaloNJets]/double"); 
    mJetData->Branch("CaloJetTrackPhiWeighted", m_CaloJetTrackPhiWeighted, "CaloJetTrackPhiWeighted[CaloNJets]/double"); 
    mJetData->Branch("CaloJetTrackNo",          m_CaloJetTrackNo,          "CaloJetTrackNo[CaloNJets]/int");
  }

  if (useJPTJets_) {
    //add jpt corrected non-cc jets
    mJetData->Branch("JPTNJets",   &m_NJPTJets,   "JPTNJets/int");  
    mJetData->Branch("JPTHt",      &m_JPTHt,      "JPTHt/double");
    mJetData->Branch("JPTMHx",     &m_JPTMHx,     "JPTMHx/double");
    mJetData->Branch("JPTMHy",     &m_JPTMHy,     "JPTMHy/double");
    mJetData->Branch("JPTMHt",     &m_JPTMHt,     "JPTMHt/double");
    
    mJetData->Branch("JPTJetE",    m_JPTJetE,    "JPTJetE[JPTNJets]/double");
    mJetData->Branch("JPTJetEt",   m_JPTJetEt,   "JPTJetEt[JPTNJets]/double");
    mJetData->Branch("JPTJetpt",   m_JPTJetPt,   "JPTJetpt[JPTNJets]/double");
    mJetData->Branch("JPTJetpx",   m_JPTJetPx,   "JPTJetpx[JPTNJets]/double");
    mJetData->Branch("JPTJetpy",   m_JPTJetPy,   "JPTJetpy[JPTNJets]/double");
    mJetData->Branch("JPTJetpz",   m_JPTJetPz,   "JPTJetpz[JPTNJets]/double");
    mJetData->Branch("JPTJeteta",  m_JPTJetEta,  "JPTJeteta[JPTNJets]/double");
    mJetData->Branch("JPTJetphi",  m_JPTJetPhi,  "JPTJetphi[JPTNJets]/double");
    mJetData->Branch("JPTJetFem",  m_JPTJetFem,  "JPTJetFem[JPTNJets]/double");
    mJetData->Branch("JPTJetfHPD", m_JPTJetfHPD, "JPTJetfHPD[JPTNJets]/double");
    mJetData->Branch("JPTJetfRBX", m_JPTJetfRBX, "JPTJetfRBX[JPTNJets]/double");
    mJetData->Branch("JPTJetn90",  m_JPTJetN90,  "JPTJetn90[JPTNJets]/double");
    //mJetData->Branch("JPTJetHemi", m_JPTJetHemi, "JPTJetHemi[JPTNJets]/int");
    mJetData->Branch("JPTJetPreselection", &m_jptJetPreselection, "JPTJetPreselection/bool");
    
    //information about associated tracks
    mJetData->Branch("JPTJetTrackPt",          m_JPTJetTrackPt,          "JPTJetTrackPt[JPTNJets]/double"); 
    mJetData->Branch("JPTJetTrackPhi",         m_JPTJetTrackPhi,         "JPTJetTrackPhi[JPTNJets]/double"); 
    mJetData->Branch("JPTJetTrackPhiWeighted", m_JPTJetTrackPhiWeighted, "JPTJetTrackPhiWeighted[JPTNJets]/double"); 
    mJetData->Branch("JPTJetTrackNo",          m_JPTJetTrackNo,          "JPTJetTrackNo[JPTNJets]/int");
  }
  
  edm::LogInfo("DiJetEvent::JetAnalyzer") << "Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_EDM_PLUGIN(JetAnalyzer);
