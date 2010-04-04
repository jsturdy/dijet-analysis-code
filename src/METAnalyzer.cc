
// -*- C++ -*-
//
// Package:    DiJetAnalysis
// Class:      METAnalyzer
// 
/**\class METAnalyzer METAnalyzer.cc JSturdy/DiJetAnalysis/src/METAnalyzer.cc


*/
//
// Original Author:  Jared Sturdy
//         Created:  Tue Feb 2 12:11:44 PDT 2010
// $Id: METAnalyzer.cc,v 1.2 2010/03/29 11:19:36 sturdy Exp $
//
//
#include "JSturdy/DiJetAnalysis/interface/METAnalyzer.h"

#include <TMath.h>

//________________________________________________________________________________________
METAnalyzer::METAnalyzer(const edm::ParameterSet& pset, TTree* tmpAllData)
{ 
  metParams = pset;
  mMETData = tmpAllData;

  _doMCData  = true;
  _doCaloMET = true;
  _doPfMET   = true;
  _doTcMET   = true;
  debug_     = 0;

  if (metParams.exists("debugMET"))     debug_   = metParams.getUntrackedParameter<int>("debugMET");

  // get the data tags
  if (metParams.exists("doMCMET"))   _doMCData  = metParams.getUntrackedParameter<bool>("doMCMET");
  if (metParams.exists("doCaloMET")) _doCaloMET = metParams.getUntrackedParameter<bool>("doCaloMET");
  if (metParams.exists("doPfMET"))   _doPfMET   = metParams.getUntrackedParameter<bool>("doPfMET");
  if (metParams.exists("doTcMET"))   _doTcMET   = metParams.getUntrackedParameter<bool>("doTcMET");

  if (_doMCData)  genTag_   = metParams.getUntrackedParameter<edm::InputTag>("genMETTag");
  if (_doCaloMET) metTag_   = metParams.getUntrackedParameter<edm::InputTag>("metTag");
  if (_doPfMET)   pfmetTag_ = metParams.getUntrackedParameter<edm::InputTag>("pfmetTag");
  if (_doTcMET)   tcmetTag_ = metParams.getUntrackedParameter<edm::InputTag>("tcmetTag");
  //mhtTag_   = metParams.getUntrackedParameter<edm::InputTag>("mhtTag");

  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
  initTuple();
}


//________________________________________________________________________________________
METAnalyzer::~METAnalyzer() {}


//________________________________________________________________________________________
// Method called to for each event
bool
METAnalyzer::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;
  using namespace edm;

  met_result = true;
  edm::LogVerbatim("DiJetEvent::METAnalyzer") << " Start  " << std::endl;

  //**********************************
  // get the PF MET result
  //**********************************
  if (_doPfMET) {
    edm::Handle< std::vector<pat::MET> > pfmetHandle;
    iEvent.getByLabel(pfmetTag_, pfmetHandle);
    if ( !pfmetHandle.isValid() ) {
      edm::LogWarning("METEventSelector") << "No Met results for InputTag " << pfmetTag_;
      return false;
    }
    
    // sanity check on collection
    if ( pfmetHandle->size()!=1 ) {
      edm::LogWarning("METEventSelector") << "MET collection size is "
					  << pfmetHandle->size() << " instead of 1";
      return false;
    }
    
    if(pfmetHandle->front().genMET()!=NULL) {
      const reco::GenMET* myGenMet = pfmetHandle->front().genMET();
      m_pfMET_Gen[0] = myGenMet->px();
      m_pfMET_Gen[1] = myGenMet->py();
      m_pfMET_Gen[2] = myGenMet->pz();
    }
    else{
      m_pfMET_Gen[0] = -99999999;
      m_pfMET_Gen[1] = -99999999;
      m_pfMET_Gen[2] = -99999999;
    }
    
    // Do the MET save for full corr no cc pfMET
    m_pfMET_Fullcorr_nocc[0]           = pfmetHandle->front().momentum().X();
    m_pfMET_Fullcorr_nocc[1]           = pfmetHandle->front().momentum().Y();
    m_pfMET_Fullcorr_nocc[2]           = pfmetHandle->front().momentum().z();
    m_pfMETphi_Fullcorr_nocc           = pfmetHandle->front().phi();
    m_pfMETsumEt_Fullcorr_nocc         = pfmetHandle->front().sumEt();
    m_pfMETsignificance_Fullcorr_nocc  = pfmetHandle->front().mEtSig();
    
    // Do the MET save for no corr no cc pfMET
    m_pfMET_Nocorr_nocc[0]   = pfmetHandle->front().corEx();                                       //uncorr to bare bones
    m_pfMET_Nocorr_nocc[1]   = pfmetHandle->front().corEy(pat::MET::UncorrectionType(0));          //uncorr to bare bones
    m_pfMETpt_Nocorr_nocc    = pfmetHandle->front().uncorrectedPt(pat::MET::UncorrectionType(0));  //uncorr to bare bones
    m_pfMETphi_Nocorr_nocc   = pfmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(0)); //uncorr to bare bones
    m_pfMETsumEt_Nocorr_nocc = pfmetHandle->front().corSumEt(pat::MET::UncorrectionType(0));       //uncorr to bare bones
    
    // Do the MET save for muon corr no cc pfMET
    // i.e., remove the JEC corrections
    m_pfMET_Muoncorr_nocc[0]   = pfmetHandle->front().corEx(pat::MET::UncorrectionType(1));          //uncorr for JEC
    m_pfMET_Muoncorr_nocc[1]   = pfmetHandle->front().corEy(pat::MET::UncorrectionType(1));          //uncorr for JEC 
    m_pfMETpt_Muoncorr_nocc    = pfmetHandle->front().uncorrectedPt(pat::MET::UncorrectionType(1));  //uncorr for JEC
    m_pfMETphi_Muoncorr_nocc   = pfmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(1)); //uncorr for JEC
    m_pfMETsumEt_Muoncorr_nocc = pfmetHandle->front().corSumEt(pat::MET::UncorrectionType(1));       //uncorr for JEC
    
    // Do the MET save for JEC corr no cc pfMET
    // i.e., remove the muon corrections
    m_pfMET_JECcorr_nocc[0]   = pfmetHandle->front().corEx(pat::MET::UncorrectionType(2));          //uncorr for muons
    m_pfMET_JECcorr_nocc[1]   = pfmetHandle->front().corEy(pat::MET::UncorrectionType(2));          //uncorr for muons
    m_pfMETpt_JECcorr_nocc    = pfmetHandle->front().uncorrectedPt(pat::MET::UncorrectionType(2));  //uncorr for muons
    m_pfMETphi_JECcorr_nocc   = pfmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(2)); //uncorr for muons
    m_pfMETsumEt_JECcorr_nocc = pfmetHandle->front().corSumEt(pat::MET::UncorrectionType(2));       //uncorr for muons
  }

  //**********************************
  // get the tcMET result
  //**********************************
  if (_doTcMET) {
    edm::Handle< std::vector<pat::MET> > tcmetHandle;
    iEvent.getByLabel(tcmetTag_, tcmetHandle);
    if ( !tcmetHandle.isValid() ) {
      edm::LogWarning("METEventSelector") << "No Met results for InputTag " << tcmetTag_;
      return false;
    }
    
    // sanity check on collection
    if ( tcmetHandle->size()!=1 ) {
      edm::LogWarning("METEventSelector") << "tcMET collection size is "
					  << tcmetHandle->size() << " instead of 1";
      return false;
    }
    
    if(tcmetHandle->front().genMET()!=NULL) {
      const reco::GenMET* myGenMet = tcmetHandle->front().genMET();
      m_tcMET_Gen[0] = myGenMet->px();
      m_tcMET_Gen[1] = myGenMet->py();
      m_tcMET_Gen[2] = myGenMet->pz();
    }
    else{
      m_tcMET_Gen[0] = -99999999;
      m_tcMET_Gen[1] = -99999999;
      m_tcMET_Gen[2] = -99999999;
    }
    
    // Do the MET save for full corr no cc tcMET
    m_tcMET_Fullcorr_nocc[0]           = tcmetHandle->front().momentum().X();
    m_tcMET_Fullcorr_nocc[1]           = tcmetHandle->front().momentum().Y();
    m_tcMET_Fullcorr_nocc[2]           = tcmetHandle->front().momentum().z();
    m_tcMETphi_Fullcorr_nocc           = tcmetHandle->front().phi();
    m_tcMETsumEt_Fullcorr_nocc         = tcmetHandle->front().sumEt();
    m_tcMETsignificance_Fullcorr_nocc  = tcmetHandle->front().mEtSig();
    
    // Do the MET save for no corr no cc tcMET
    m_tcMET_Nocorr_nocc[0]   = tcmetHandle->front().corEx();                                       //uncorr to bare bones
    m_tcMET_Nocorr_nocc[1]   = tcmetHandle->front().corEy(pat::MET::UncorrectionType(0));          //uncorr to bare bones
    m_tcMETpt_Nocorr_nocc    = tcmetHandle->front().uncorrectedPt(pat::MET::UncorrectionType(0));  //uncorr to bare bones
    m_tcMETphi_Nocorr_nocc   = tcmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(0)); //uncorr to bare bones
    m_tcMETsumEt_Nocorr_nocc = tcmetHandle->front().corSumEt(pat::MET::UncorrectionType(0));       //uncorr to bare bones
    
    // Do the MET save for muon corr no cc tcMET
    // i.e., remove the JEC corrections
    m_tcMET_Muoncorr_nocc[0]   = tcmetHandle->front().corEx(pat::MET::UncorrectionType(1));          //uncorr for JEC
    m_tcMET_Muoncorr_nocc[1]   = tcmetHandle->front().corEy(pat::MET::UncorrectionType(1));          //uncorr for JEC 
    m_tcMETpt_Muoncorr_nocc    = tcmetHandle->front().uncorrectedPt(pat::MET::UncorrectionType(1));  //uncorr for JEC
    m_tcMETphi_Muoncorr_nocc   = tcmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(1)); //uncorr for JEC
    m_tcMETsumEt_Muoncorr_nocc = tcmetHandle->front().corSumEt(pat::MET::UncorrectionType(1));       //uncorr for JEC
    
    // Do the MET save for JEC corr no cc tcMET
    // i.e., remove the muon corrections
    m_tcMET_JECcorr_nocc[0]   = tcmetHandle->front().corEx(pat::MET::UncorrectionType(2));          //uncorr for muons
    m_tcMET_JECcorr_nocc[1]   = tcmetHandle->front().corEy(pat::MET::UncorrectionType(2));          //uncorr for muons
    m_tcMETpt_JECcorr_nocc    = tcmetHandle->front().uncorrectedPt(pat::MET::UncorrectionType(2));  //uncorr for muons
    m_tcMETphi_JECcorr_nocc   = tcmetHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(2)); //uncorr for muons
    m_tcMETsumEt_JECcorr_nocc = tcmetHandle->front().corSumEt(pat::MET::UncorrectionType(2));       //uncorr for muons
  }
  
  //**********************************
  // get the calo (AK5) MET result
  //**********************************
  if (_doCaloMET) {
    edm::Handle< std::vector<pat::MET> > metHandle;
    iEvent.getByLabel(metTag_, metHandle);
    if ( !metHandle.isValid() ) {
      edm::LogWarning("METEventSelector") << "No Met results for InputTag " << metTag_;
      return false;
    }

    // sanity check on collection
    if ( metHandle->size()!=1 ) {
      edm::LogWarning("METEventSelector") << "MET collection size is "
					  << metHandle->size() << " instead of 1";
      return false;
    }
 
    if(metHandle->front().genMET()!=NULL) {
      const reco::GenMET* myGenMet = metHandle->front().genMET();
      m_MET_Gen[0] = myGenMet->px();
      m_MET_Gen[1] = myGenMet->py();
      m_MET_Gen[2] = myGenMet->pz();
    }
    else{
      m_MET_Gen[0] = -99999999;
      m_MET_Gen[1] = -99999999;
      m_MET_Gen[2] = -99999999;
    }

    // Do the MET save for full corr no cc MET
    m_MET_Fullcorr_nocc[0]           = metHandle->front().momentum().X();
    m_MET_Fullcorr_nocc[1]           = metHandle->front().momentum().Y();
    m_MET_Fullcorr_nocc[2]           = metHandle->front().momentum().z();
    m_METphi_Fullcorr_nocc           = metHandle->front().phi();
    m_METsumEt_Fullcorr_nocc         = metHandle->front().sumEt();
    m_METsignificance_Fullcorr_nocc  = metHandle->front().mEtSig();

    // Do the MET save for no corr no cc MET
    m_MET_Nocorr_nocc[0]   = metHandle->front().corEx();                                       //uncorr to bare bones
    m_MET_Nocorr_nocc[1]   = metHandle->front().corEy(pat::MET::UncorrectionType(0));          //uncorr to bare bones
    m_METpt_Nocorr_nocc    = metHandle->front().uncorrectedPt(pat::MET::UncorrectionType(0));  //uncorr to bare bones
    m_METphi_Nocorr_nocc   = metHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(0)); //uncorr to bare bones
    m_METsumEt_Nocorr_nocc = metHandle->front().corSumEt(pat::MET::UncorrectionType(0));       //uncorr to bare bones

    // Do the MET save for muon corr no cc MET
    // i.e., remove JEC corrections
    m_MET_Muoncorr_nocc[0]   = metHandle->front().corEx(pat::MET::UncorrectionType(1));          //uncorr for JEC
    m_MET_Muoncorr_nocc[1]   = metHandle->front().corEy(pat::MET::UncorrectionType(1));          //uncorr for JEC 
    m_METpt_Muoncorr_nocc    = metHandle->front().uncorrectedPt(pat::MET::UncorrectionType(1));  //uncorr for JEC
    m_METphi_Muoncorr_nocc   = metHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(1)); //uncorr for JEC
    m_METsumEt_Muoncorr_nocc = metHandle->front().corSumEt(pat::MET::UncorrectionType(1));       //uncorr for JEC

    // Do the MET save for JEC corr no cc MET
    // i.e., remove muon corrections
    m_MET_JECcorr_nocc[0]   = metHandle->front().corEx(pat::MET::UncorrectionType(2));          //uncorr for muons
    m_MET_JECcorr_nocc[1]   = metHandle->front().corEy(pat::MET::UncorrectionType(2));          //uncorr for muons
    m_METpt_JECcorr_nocc    = metHandle->front().uncorrectedPt(pat::MET::UncorrectionType(2));  //uncorr for muons
    m_METphi_JECcorr_nocc   = metHandle->front().uncorrectedPhi(pat::MET::UncorrectionType(2)); //uncorr for muons
    m_METsumEt_JECcorr_nocc = metHandle->front().corSumEt(pat::MET::UncorrectionType(2));       //uncorr for muons
  }

  //
  // sanity check on collection
  //
  nUncorrMET = 2;
  nFullMET   = 3;

  // Fill the tree only if all preselection conditions are met
  //if(preselection) //mPreselection->Fill();
  //mMETData->Fill();
  return met_result;
}

//________________________________________________________________________________________

//________________________________________________________________________________________
void 
METAnalyzer::beginJob(const edm::EventSetup&) {}

//________________________________________________________________________________________
void 
METAnalyzer::endJob() {

}


//________________________________________________________________________________________
void
METAnalyzer::initTuple() {

  std::ostringstream variables; // Container for all variables
  

  // Add the branches
  //general MET information
  mMETData->Branch("nFullMET",   &nFullMET,   "nFullMET/int");
  mMETData->Branch("nUncorrMET", &nUncorrMET, "nUncorrMET/int");

  //pfMET information
  if (_doPfMET) {
    mMETData->Branch("pfMET_fullcorr_nocc",              m_pfMET_Fullcorr_nocc,             "m_pfMET_Fullcorr_nocc[nFullMET]/double");
    mMETData->Branch("pfMETphi_fullcorr_nocc",          &m_pfMETphi_Fullcorr_nocc,          "m_pfMETphi_Fullcorr_nocc/double");
    mMETData->Branch("pfMETsumEt_fullcorr_nocc",        &m_pfMETsumEt_Fullcorr_nocc,        "m_pfMETsumEt_Fullcorr_nocc/double");
    mMETData->Branch("pfMETsignificance_Fullcorr_nocc", &m_pfMETsignificance_Fullcorr_nocc, "m_pfMETsignificance_Fullcorr_nocc/double");
    
    mMETData->Branch("pfMET_nocorr_nocc",       m_pfMET_Nocorr_nocc,      "m_pfMET_Nocorr_nocc[nUncorrMET]/double");
    mMETData->Branch("pfMETphi_nocorr_nocc",   &m_pfMETphi_Nocorr_nocc,   "m_pfMETphi_Nocorr_nocc/double");
    mMETData->Branch("pfMETpt_nocorr_nocc",    &m_pfMETpt_Nocorr_nocc,    "m_pfMETpt_Nocorr_nocc/double");
    mMETData->Branch("pfMETsumEt_nocorr_nocc", &m_pfMETsumEt_Nocorr_nocc, "m_pfMETsumEt_Nocorr_nocc/double");
    
    
    mMETData->Branch("pfMET_muoncorr_nocc",       m_pfMET_Muoncorr_nocc,      "m_pfMET_Muoncorr_nocc[nUncorrMET]/double");
    mMETData->Branch("pfMETphi_muoncorr_nocc",   &m_pfMETphi_Muoncorr_nocc,   "m_pfMETphi_Muoncorr_nocc/double");
    mMETData->Branch("pfMETpt_muoncorr_nocc",    &m_pfMETpt_Muoncorr_nocc,    "m_pfMETpt_Muoncorr_nocc/double");
    mMETData->Branch("pfMETsumEt_muoncorr_nocc", &m_pfMETsumEt_Muoncorr_nocc, "m_pfMETsumEt_Muoncorr_nocc/double");
    
    mMETData->Branch("pfMET_jeccorr_nocc",        m_pfMET_JECcorr_nocc,       "m_pfMET_JECcorr_nocc[nUncorrMET]/double");
    mMETData->Branch("pfMETphi_jeccorr_nocc",    &m_pfMETphi_JECcorr_nocc,    "m_pfMETphi_JECcorr_nocc/double");
    mMETData->Branch("pfMETpt_jeccorr_nocc",     &m_pfMETpt_JECcorr_nocc,     "m_pfMETpt_JECcorr_nocc/double");
    mMETData->Branch("pfMETsumEt_jeccorr_nocc",  &m_pfMETsumEt_JECcorr_nocc,  "m_pfMETsumEt_JECcorr_nocc/double");
    
    mMETData->Branch("GenpfMET", &m_pfMET_Gen, "m_pfMET_Gen[3]/double",6400);
  }

  //tcMET information
  if (_doTcMET) {
    mMETData->Branch("tcMET_fullcorr_nocc",              m_tcMET_Fullcorr_nocc,             "m_tcMET_Fullcorr_nocc[nFullMET]/double");
    mMETData->Branch("tcMETphi_fullcorr_nocc",          &m_tcMETphi_Fullcorr_nocc,          "m_tcMETphi_Fullcorr_nocc/double");
    mMETData->Branch("tcMETsumEt_fullcorr_nocc",        &m_tcMETsumEt_Fullcorr_nocc,        "m_tcMETsumEt_Fullcorr_nocc/double");
    mMETData->Branch("tcMETsignificance_Fullcorr_nocc", &m_tcMETsignificance_Fullcorr_nocc, "m_tcMETsignificance_Fullcorr_nocc/double");
    
    mMETData->Branch("tcMET_nocorr_nocc",       m_tcMET_Nocorr_nocc,      "m_tcMET_Nocorr_nocc[nUncorrMET]/double");
    mMETData->Branch("tcMETpt_nocorr_nocc",    &m_tcMETpt_Nocorr_nocc,    "m_tcMETpt_Nocorr_nocc/double");
    mMETData->Branch("tcMETphi_nocorr_nocc",   &m_tcMETphi_Nocorr_nocc,   "m_tcMETphi_Nocorr_nocc/double");
    mMETData->Branch("tcMETsumEt_nocorr_nocc", &m_tcMETsumEt_Nocorr_nocc, "m_tcMETsumEt_Nocorr_nocc/double");
    
    mMETData->Branch("tcMET_muoncorr_nocc",       m_tcMET_Muoncorr_nocc,      "m_tcMET_Muoncorr_nocc[nUncorrMET]/double");
    mMETData->Branch("tcMETpt_muoncorr_nocc",    &m_tcMETpt_Muoncorr_nocc,    "m_tcMETpt_Muoncorr_nocc/double");
    mMETData->Branch("tcMETphi_muoncorr_nocc",   &m_tcMETphi_Muoncorr_nocc,   "m_tcMETphi_Muoncorr_nocc/double");
    mMETData->Branch("tcMETsumEt_muoncorr_nocc", &m_tcMETsumEt_Muoncorr_nocc, "m_tcMETsumEt_Muoncorr_nocc/double");
    
    mMETData->Branch("tcMET_jeccorr_nocc",        m_tcMET_JECcorr_nocc,       "m_tcMET_JECcorr_nocc[nUncorrMET]/double");
    mMETData->Branch("tcMETpt_jeccorr_nocc",     &m_tcMETpt_JECcorr_nocc,     "m_tcMETpt_JECcorr_nocc/double");
    mMETData->Branch("tcMETphi_jeccorr_nocc",    &m_tcMETphi_JECcorr_nocc,    "m_tcMETphi_JECcorr_nocc/double");
    mMETData->Branch("tcMETsumEt_jeccorr_nocc",  &m_tcMETsumEt_JECcorr_nocc,  "m_tcMETsumEt_JECcorr_nocc/double");
    
    mMETData->Branch("GentcMET", &m_tcMET_Gen, "m_tcMET_Gen[3]/double",6400);
  }

  //calo MET information
  if (_doCaloMET) {
    mMETData->Branch("MET_fullcorr_nocc",              m_MET_Fullcorr_nocc,             "m_MET_Fullcorr_nocc[nFullMET]/double");
    mMETData->Branch("METphi_fullcorr_nocc",          &m_METphi_Fullcorr_nocc,          "m_METphi_Fullcorr_nocc/double");
    mMETData->Branch("METsumEt_fullcorr_nocc",        &m_METsumEt_Fullcorr_nocc,        "m_METsumEt_Fullcorr_nocc/double");
    mMETData->Branch("METsignificance_Fullcorr_nocc", &m_METsignificance_Fullcorr_nocc, "m_METsignificance_Fullcorr_nocc/double");
    
    mMETData->Branch("MET_nocorr_nocc",       m_MET_Nocorr_nocc,      "m_MET_Nocorr_nocc[nUncorrMET]/double");
    mMETData->Branch("METpt_nocorr_nocc",    &m_METpt_Nocorr_nocc,    "m_METpt_Nocorr_nocc/double");
    mMETData->Branch("METphi_nocorr_nocc",   &m_METphi_Nocorr_nocc,   "m_METphi_Nocorr_nocc/double");
    mMETData->Branch("METsumEt_nocorr_nocc", &m_METsumEt_Nocorr_nocc, "m_METsumEt_Nocorr_nocc/double");
    
    mMETData->Branch("MET_muoncorr_nocc",       m_MET_Muoncorr_nocc,      "m_MET_Muoncorr_nocc[nUncorrMET]/double");
    mMETData->Branch("METpt_muoncorr_nocc",    &m_METpt_Muoncorr_nocc,    "m_METpt_Muoncorr_nocc/double");
    mMETData->Branch("METphi_muoncorr_nocc",   &m_METphi_Muoncorr_nocc,   "m_METphi_Muoncorr_nocc/double");
    mMETData->Branch("METsumEt_muoncorr_nocc", &m_METsumEt_Muoncorr_nocc, "m_METsumEt_Muoncorr_nocc/double");
    
    mMETData->Branch("MET_jeccorr_nocc",       m_MET_JECcorr_nocc,       "m_MET_JECcorr_nocc[nUncorrMET]/double");
    mMETData->Branch("METpt_jeccorr_nocc",    &m_METpt_JECcorr_nocc,     "m_METpt_JECcorr_nocc/double");
    mMETData->Branch("METphi_jeccorr_nocc",   &m_METphi_JECcorr_nocc,    "m_METphi_JECcorr_nocc/double");
    mMETData->Branch("METsumEt_jeccorr_nocc",  &m_METsumEt_JECcorr_nocc, "m_METsumEt_JECcorr_nocc/double");
    
    mMETData->Branch("GenMET", &m_MET_Gen, "m_MET_Gen[3]/double",6400);
  }

  //cc MET information
//  mMETData->Branch("MET_fullcorr_cc",              m_MET_Fullcorr_cc,             "m_MET_Fullcorr_cc[nFullMET]/double");
//  mMETData->Branch("METphi_fullcorr_cc",          &m_METphi_Fullcorr_cc,          "m_METphi_Fullcorr_cc/double");
//  mMETData->Branch("METsignificance_fullcorr_cc", &m_METsignificance_Fullcorr_cc, "m_METsignificance_Fullcorr_cc/double");
//  
//  mMETData->Branch("MPTPhi", &m_MPTPhi, "MPTPhi/double");
//  mMETData->Branch("MPTPx",  &m_MPTPx,  "MPTPx/double");
//  mMETData->Branch("MPTPy",  &m_MPTPy,  "MPTPy/double");
//  mMETData->Branch("MPTPz",  &m_MPTPz,  "MPTPz/double");
  
  edm::LogInfo("DiJetEvent::METAnalyzer") << "MET Ntuple variables " << variables.str();
  
}


//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(METAnalyzer);
