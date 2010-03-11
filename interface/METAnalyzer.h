#ifndef METANALYZER
#define METANALYZER

// System include files
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

// ROOT includes
#include <TNtuple.h>

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// SUSY include files
//#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
//#include "SusyAnalysis/EventSelector/interface/SelectorSequence.h"

#include "DataFormats/PatCandidates/interface/MET.h"


//
// Class declaration
//
class METAnalyzer {
 public:
  METAnalyzer(const edm::ParameterSet&);
  ~METAnalyzer();
  
 private:
  //*** CMSSW interface
  /// Called once per job, at start
  void beginJob(const edm::EventSetup&) ;
  /// Called for each event
  //virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool filter(const edm::Event& evt,const edm::EventSetup& iSetup );
  /// Called once per job, at end
  void endJob();

  //*** Plotting
  /// Define all plots
  void initTuple();
  //void bookHealthPlots();
  /// Fill all plots for an event
  void fillTuple( const edm::Event& ) {mMETData->Fill();}
  //void fillHealthPlots( const edm::Event& ) {mMETData->Fill();}
  
 private:
  
  // configuration parameters
  //met tags
  edm::InputTag tcmetTag_;
  edm::InputTag pfmetTag_;
  edm::InputTag metTag_;
  edm::InputTag genTag_;
  //edm::InputTag mhtTag_;
  bool _doMCData;
  bool _doCaloMET;
  bool _doPfMET;
  bool _doTcMET;

  bool met_result;       /// result of the met cut

  TTree * mMETData;      /// Will contain the additional di-jet specific data

  // Generated MET
  double m_pfMET_Gen[3];
  double m_tcMET_Gen[3];
  double m_MET_Gen[3];

  int nFullMET;
 
  // Do the MET save for non cc pfMET
  double m_pfMET_Fullcorr_nocc[3];
  double m_pfMETphi_Fullcorr_nocc;
  double m_pfMETsumEt_Fullcorr_nocc;
  double m_pfMETsignificance_Fullcorr_nocc;

  double m_pfMET_Nocorr_nocc[2];
  double m_pfMETpt_Nocorr_nocc;
  double m_pfMETphi_Nocorr_nocc;
  double m_pfMETsumEt_Nocorr_nocc;
  double m_pfMETsignificance_Nocorr_nocc;

  double m_pfMET_Muoncorr_nocc[2];
  double m_pfMETpt_Muoncorr_nocc;
  double m_pfMETphi_Muoncorr_nocc;
  double m_pfMETsumEt_Muoncorr_nocc;
  double m_pfMETsignificance_Muoncorr_nocc;

  double m_pfMET_JECcorr_nocc[2];
  double m_pfMETpt_JECcorr_nocc;
  double m_pfMETphi_JECcorr_nocc;
  double m_pfMETsumEt_JECcorr_nocc;
  double m_pfMETsignificance_JECcorr_nocc;

  // Do the MET save for non cc tcMET
  double m_tcMET_Fullcorr_nocc[3];
  double m_tcMETphi_Fullcorr_nocc;
  double m_tcMETsumEt_Fullcorr_nocc;
  double m_tcMETsignificance_Fullcorr_nocc;

  double m_tcMET_Nocorr_nocc[2];
  double m_tcMETpt_Nocorr_nocc;
  double m_tcMETphi_Nocorr_nocc;
  double m_tcMETsumEt_Nocorr_nocc;
  double m_tcMETsignificance_Nocorr_nocc;

  double m_tcMET_Muoncorr_nocc[2];
  double m_tcMETpt_Muoncorr_nocc;
  double m_tcMETphi_Muoncorr_nocc;
  double m_tcMETsumEt_Muoncorr_nocc;
  double m_tcMETsignificance_Muoncorr_nocc;

  double m_tcMET_JECcorr_nocc[2];
  double m_tcMETpt_JECcorr_nocc;
  double m_tcMETphi_JECcorr_nocc;
  double m_tcMETsumEt_JECcorr_nocc;
  double m_tcMETsignificance_JECcorr_nocc;

  // Do the MET save for non cc calo MET
  double m_MET_Fullcorr_nocc[3];
  double m_METphi_Fullcorr_nocc;
  double m_METsumEt_Fullcorr_nocc;
  double m_METsignificance_Fullcorr_nocc;

  double m_MET_Nocorr_nocc[2];
  double m_METpt_Nocorr_nocc;
  double m_METphi_Nocorr_nocc;
  double m_METsumEt_Nocorr_nocc;
  double m_METsignificance_Nocorr_nocc;

  double m_MET_Muoncorr_nocc[2];
  double m_METpt_Muoncorr_nocc;
  double m_METphi_Muoncorr_nocc;
  double m_METsumEt_Muoncorr_nocc;
  double m_METsignificance_Muoncorr_nocc;

  double m_MET_JECcorr_nocc[2];
  double m_METpt_JECcorr_nocc;
  double m_METphi_JECcorr_nocc;
  double m_METsumEt_JECcorr_nocc;
  double m_METsignificance_JECcorr_nocc;

  int nUncorrMET;

  //// Do the MET save for cc MET
  //double m_MET_Fullcorr_cc[3];
  //double m_METphi_Fullcorr_cc;
  //double m_METsumEt_Fullcorr_cc;
  //double m_METsignificance_Fullcorr_cc;
  //
  //double m_MET_Nocorr_cc[2];
  //double m_METpt_Nocorr_cc;
  //double m_METphi_Nocorr_cc;
  //double m_METsumEt_Nocorr_cc;
  //double m_METsignificance_Nocorr_cc;
  //
  //double m_MET_Muoncorr_cc[2];
  //double m_METpt_Muoncorr_cc;
  //double m_METphi_Muoncorr_cc;
  //double m_METsumEt_Muoncorr_cc;
  //double m_METsignificance_Muoncorr_cc;
  //
  //double m_MET_JECcorr_cc[2];
  //double m_METpt_JECcorr_cc;
  //double m_METphi_JECcorr_cc;
  //double m_METsumEt_JECcorr_cc;
  //double m_METsignificance_JECcorr_cc;

  double localPi;
};

#endif
