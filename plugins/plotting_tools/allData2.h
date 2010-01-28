//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun  8 08:26:41 2009 by ROOT version 5.18/00a
// from TTree allData/data after cuts
// found on file: lm1_NT3.root
//////////////////////////////////////////////////////////

#ifndef allData2_h
#define allData2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include "TLorentzVector.h"

class allData2 {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  UInt_t          run;
  UInt_t          event;

  Int_t           nHLT;
  Int_t           HLTArray[159];
  UChar_t         HLT1JET;
  UChar_t         HLT2JET;
  UChar_t         HLT1MET1HT;
  UChar_t         HLT1MUON;


  UInt_t          genN;
  UInt_t          genid[50];//1000
  UInt_t          genMother[50];//1000
  Float_t         genE[50];//1000
  Float_t         genPx[50];//1000
  Float_t         genPy[50];//1000
  Float_t         genPz[50];//1000
  UInt_t          genStatus[50];//1000

  UInt_t          genLepN;
  UInt_t          genLepId[50];//100
  UInt_t          genLepMother[50];//100
  Float_t         genLepE[50];//100
  Float_t         genLepPx[50];//100
  Float_t         genLepPy[50];//100
  Float_t         genLepPz[50];//100
  UInt_t          genLepStatus[50];//100


  UInt_t          nVtx;
  Double_t        VtxChi2[5];
  Double_t        VtxNdof[5];
  Double_t        VtxNormalizedChi2[5];
  Double_t        VtxX[5];
  Double_t        VtxY[5];
  Double_t        VtxZ[5];
  Double_t        VtxdX[5]; 
  Double_t        VtxdY[5]; 
  Double_t        VtxdZ[5];

  UInt_t          nFullMET;
  Double_t        MET_fullcorr_nocc[3];
  Double_t        MET_fullcorr_cc[3];
  Double_t        METphi_fullcorr_nocc;
  UInt_t          nUncorrMET;
  Double_t        MET_nocorr_nocc[2];
  //Double_t        METphi_nocorr_nocc;
  Double_t        MET_jeccorr_nocc[2];
  Double_t        METphi_jeccorr_nocc;
  Double_t        MET_muoncorr_nocc[2];
  Double_t        METphi_muoncorr_nocc;
  Double_t        MET_Fullcorr_nocc_significance;
  Double_t        GenMET[3];

  Double_t        evtWeight;
  UInt_t          procID;
  Double_t        pthat;


  UInt_t          Nhemispheres;
  Double_t        HemisphereE[2];
  Double_t        HemisphereEt[2];
  Double_t        Hemispherept[2];
  Double_t        Hemispherepx[2];
  Double_t        Hemispherepy[2];
  Double_t        Hemispherepz[2];
  Double_t        Hemisphereeta[2];
  Double_t        Hemispherephi[2];
  Double_t        Hemispheredphi[2];

  UInt_t          Njets;
  Double_t        Ht;
  Double_t        MHt;
  Double_t        JetE[50];
  Double_t        JetEt[50];
  Double_t        Jetpt[50];
  Double_t        Jetpx[50];
  Double_t        Jetpy[50];
  Double_t        Jetpz[50];
  Double_t        Jeteta[50];
  Double_t        Jetphi[50];
  Double_t        JetFem[50];
  UInt_t          JetHemi[50];
  Double_t        Jet_MCcorrFactor[50];
  Double_t        Jet_JPTcorrFactor[50];
  Float_t         JetsBTag_TkCountHighEff[50];
  Float_t         JetsBTag_SimpleSecVtx[50];
  Float_t         JetsBTag_CombSecVtx[50];
  UChar_t         Jet_isccJetAssoc[50];
  Double_t        Jet_ccJetE[50];
  Double_t        Jet_ccJetpx[50];
  Double_t        Jet_ccJetpy[50];
  Double_t        Jet_ccJetpz[50];
  UInt_t          NPFjet;
  Double_t        PFHt;
  Double_t        PFMHt;
  Double_t        PFjetEta[50];
  Double_t        PFjetPhi[50];
  Double_t        PFjetE[50];
  Double_t        PFjetPx[50];
  Double_t        PFjetPy[50];
  Double_t        PFjetPz[50];
  Double_t        PFjetPt[50];
  Double_t        PFjetCharge[50];


  Double_t        GenHt;
  Double_t        GenMHt;
  Double_t        GenJetE[50];
  Double_t        GenJetEt[50];
  Double_t        GenJetpt[50];
  Double_t        GenJetpx[50];
  Double_t        GenJetpy[50];
  Double_t        GenJetpz[50];
  Double_t        GenJeteta[50];
  Double_t        GenJetphi[50];

  UInt_t          JetPartonId[50];
  UInt_t          JetPartonMother[50];
  Double_t        JetPartonPx[50];
  Double_t        JetPartonPy[50];
  Double_t        JetPartonPz[50];
  Double_t        JetPartonEt[50];
  Double_t        JetPartonE[50];
  Double_t        JetPartonPhi[50];
  Double_t        JetPartonEta[50];
  UInt_t          JetPartonFlavour[50];
  Double_t        JetTrackPt[50];
  Double_t        JetTrackPhi[50];
  Double_t        JetTrackPhiWeighted[50];
  UInt_t          JetTrackNo[50];

  Double_t        GenJetsEt[50];
  Double_t        GenJetsPt[50];
  Double_t        GenJetsE[50];
  Double_t        GenJetsPx[50];
  Double_t        GenJetsPy[50];
  Double_t        GenJetsPz[50];
  Double_t        GenJetsEta[50];
  Double_t        GenJetsPhi[50];
	          
  UInt_t          Nphot;
  Double_t        PhotEt[50];
  Double_t        PhotPt[50];
  Double_t        PhotPx[50];
  Double_t        PhotPy[50];
  Double_t        PhotPz[50];
  Double_t        PhotE[50];
  Double_t        PhotEta[50];
  Double_t        PhotPhi[50];
  Double_t        PhotTrkIso[50];
  Double_t        PhotECalIso[50];
  Double_t        PhotHCalIso[50];
  Double_t        PhotAllIso[50];

  Bool_t          Phot_isccPhotAssoc[50];
  Bool_t          PhotLooseEM[50];
  Bool_t          PhotLoosePhoton[50];
  Bool_t          PhotTightPhoton[50];

  Double_t        GenPhotPdgId[50];
  Double_t        GenPhotMother[50];
  Double_t        GenPhotPx[50];
  Double_t        GenPhotPy[50];
  Double_t        GenPhotPz[50];
	          
  UInt_t          Nelec;
  Double_t        ElecEt[50];
  Double_t        ElecPt[50];
  Double_t        ElecPx[50];
  Double_t        ElecPy[50];
  Double_t        ElecPz[50];
  Double_t        ElecE[50];
  Double_t        ElecEta[50];
  Double_t        ElecPhi[50];
  Double_t        ElecTrkIso[50];
  Double_t        ElecECalIso[50];
  Double_t        ElecHCalIso[50];
  Double_t        ElecAllIso[50];
  Double_t        ElecTrkChiNorm[50];
  Double_t        ElecCharge[50];
  //MICHELE       
  Double_t        ElecIdLoose[50];
  Double_t        ElecIdTight[50];
  Double_t        ElecIdRobLoose[50];
  Double_t        ElecIdRobTight[50];
  Double_t        ElecChargeMode[50];
  Double_t        ElecPtTrkMode[50];
  Double_t        ElecQOverPErrTrkMode[50];
  Double_t        GenElecPdgId[50];
  Double_t        GenElecMother[50];
  Double_t        GenElecPx[50];
  Double_t        GenElecPy[50];
  Double_t        GenElecPz[50];
  Double_t        ElecCaloEnergy[50];
  Double_t        ElecHOverE[50];
  Double_t        ElecVx[50];
  Double_t        ElecVy[50];
  Double_t        ElecVz[50];
  Double_t        ElecD0[50];
  Double_t        ElecDz[50];
  Double_t        ElecPtTrk[50];
  Double_t        ElecQOverPErrTrk[50];
  Double_t        ElecLostHits[50];
  Double_t        ElecValidHits[50];
  Double_t        ElecNCluster[50];
  Double_t        ElecEtaTrk[50];
  Double_t        ElecPhiTrk[50];
  Double_t        ElecWidthClusterEta[50];
  Double_t        ElecWidthClusterPhi[50];
  Double_t        ElecPinTrk[50];
  Double_t        ElecPoutTrk[50];
  Double_t        ElecNormChi2[50];
  Bool_t          Elec_isccElecAssoc[50];

  Double_t        ElecECalIsoDeposit[50];
  Double_t        ElecHCalIsoDeposit[50];

  UInt_t          Nmuon;
  Double_t        MuonEt[50];
  Double_t        MuonPt[50];
  Double_t        MuonPx[50];
  Double_t        MuonPy[50];
  Double_t        MuonPz[50];
  Double_t        MuonE[50];
  Double_t        MuonEta[50];
  Double_t        MuonPhi[50];
  Double_t        MuonTrkIso[50];
  Double_t        MuonECalIso[50];
  Double_t        MuonHCalIso[50];
  Double_t        MuonAllIso[50];
  Double_t        MuonTrkChiNorm[50];
  Double_t        MuonCharge[50];
  Bool_t          MuonIsGlobal[50];
  Bool_t          MuonIsStandAlone[50];
  Bool_t          MuonIsTracker[50]; 
  Bool_t          MuonIsGlobalTight[50];
  Bool_t          MuonIsTMLastStationLoose[50];
  Bool_t          MuonTMLastStationTight[50];
  Bool_t          MuonTM2DCompatibilityLoose[50];
  Bool_t          MuonTM2DCompatibilityTight[50];
  Bool_t          Muon_isccMuonAssoc[50];

  Double_t        MuonECalIsoDeposit[50];
  Double_t        MuonHCalIsoDeposit[50];
  	          
  Double_t        MuonCombChi2[50];
  Double_t        MuonCombNdof[50];
  Double_t        MuonTrkD0[50];
  	          
  //MICHELE       
  Double_t        MuonId[50];
  Double_t        MuonCombVx[50];
  Double_t        MuonCombVy[50];
  Double_t        MuonCombVz[50];
  Double_t        MuonCombD0[50];
  Double_t        MuonCombDz[50];
	          
  Double_t        MuonStandValidHits[50];
  Double_t        MuonStandLostHits[50];
  Double_t        MuonStandPt[50];
  Double_t        MuonStandPz[50];
  Double_t        MuonStandP[50];
  Double_t        MuonStandEta[50];
  Double_t        MuonStandPhi[50];
  Double_t        MuonStandChi[50];
  Double_t        MuonStandCharge[50];
  Double_t        MuonStandQOverPError[50];
	          
  Double_t        MuonTrkValidHits[50];
  Double_t        MuonTrkLostHits[50];
  Double_t        MuonTrkPt[50];
  Double_t        MuonTrkPz[50];
  Double_t        MuonTrkP[50];
  Double_t        MuonTrkEta[50];
  Double_t        MuonTrkPhi[50];
  Double_t        MuonTrkChi[50];
  Double_t        MuonTrkCharge[50];
  Double_t        MuonTrkQOverPError[50];
  Double_t        MuonTrkOuterZ[50];
  Double_t        MuonTrkOuterR[50];
	          
  Double_t        GenMuonPdgId[50];
  Double_t        GenMuonMother[50];
  Double_t        GenMuonPx[50];
  Double_t        GenMuonPy[50];
  Double_t        GenMuonPz[50];

  Double_t        MuonPairMass;
  UInt_t          MuonPairIndex[2];

  Double_t        AlpPtScale;
  UInt_t          AlpIdTest;

  Double_t        MPTPhi;
  Double_t        MPTPx;
  Double_t        MPTPy;
  Double_t        MPTPz;

  // List of branches
  TBranch        *b_run;
  TBranch        *b_event;

  TBranch        *b_nHLT;
  TBranch        *b_HLTArray;
  TBranch        *b_HLT1JET;
  TBranch        *b_HLT2JET;
  TBranch        *b_HLT1MET1HT;
  TBranch        *b_HLT1MUON;

  TBranch        *b_genN;
  TBranch        *b_genid;
  TBranch        *b_genMother;
  TBranch        *b_genE;
  TBranch        *b_genPx;
  TBranch        *b_genPy;
  TBranch        *b_genPz;
  TBranch        *b_genStatus;
  
  TBranch        *b_genLepN;
  TBranch        *b_genLepId;
  TBranch        *b_genLepMother;
  TBranch        *b_genLepE;
  TBranch        *b_genLepPx;
  TBranch        *b_genLepPy;
  TBranch        *b_genLepPz;
  TBranch        *b_genLepStatus;

  TBranch        *b_nFullMET;
  TBranch        *b_MET_fullcorr_nocc;
  TBranch        *b_METphi_fullcorr_nocc;
  TBranch        *b_MET_fullcorr_cc;
  TBranch        *b_nUncorrMET;
  TBranch        *b_MET_nocorr_nocc;
  //TBranch        *b_METphi_nocorr_nocc;
  TBranch        *b_MET_jeccorr_nocc;
  TBranch        *b_METphi_jeccorr_nocc;
  TBranch        *b_MET_muoncorr_nocc;
  TBranch        *b_METphi_muoncorr_nocc;
  TBranch        *b_MET_Fullcorr_nocc_significance;
  TBranch        *b_MET_Gen;

  TBranch        *b_evtWeight;
  TBranch        *b_procID;
  TBranch        *b_pthat;

  TBranch        *b_nVtx;
  TBranch        *b_VertexChi2;
  TBranch        *b_VertexNdof;
  TBranch        *b_VertexNormalizedChi2;
  TBranch        *b_VertexX;
  TBranch        *b_VertexY;
  TBranch        *b_VertexZ;
  TBranch        *b_VertexdX;
  TBranch        *b_VertexdY;
  TBranch        *b_VertexdZ;

  TBranch        *b_Nhemispheres;
  TBranch        *b_HemisphereE;
  TBranch        *b_HemisphereEt;
  TBranch        *b_Hemispherept;
  TBranch        *b_Hemispherepx;
  TBranch        *b_Hemispherepy;
  TBranch        *b_Hemispherepz;
  TBranch        *b_Hemisphereeta;
  TBranch        *b_Hemispherephi;
  TBranch        *b_Hemispheredphi;

  TBranch        *b_Njets;
  TBranch        *b_Ht;
  TBranch        *b_MHt;
  TBranch        *b_JetE;
  TBranch        *b_JetEt;
  TBranch        *b_Jetpt;
  TBranch        *b_Jetpx;
  TBranch        *b_Jetpy;
  TBranch        *b_Jetpz;
  TBranch        *b_Jeteta;
  TBranch        *b_Jetphi;
  TBranch        *b_JetFem;
  TBranch        *b_JetHemi;
  TBranch        *b_Jet_MCcorrFactor;
  TBranch        *b_Jet_JPTcorrFactor;
  TBranch        *b_JetBTag_TkCountHighEff;
  TBranch        *b_JetBTag_SimpleSecVtx;
  TBranch        *b_JetBTag_CombSecVtx;
  TBranch        *b_Jet_isccJetAssoc;
  TBranch        *b_Jet_ccJetE;
  TBranch        *b_Jet_ccJetpx;
  TBranch        *b_Jet_ccJetpy;
  TBranch        *b_Jet_ccJetpz;

  TBranch        *b_GenHt;
  TBranch        *b_GenMHt;
  TBranch        *b_GenJetE;
  TBranch        *b_GenJetEt;
  TBranch        *b_GenJetpt;
  TBranch        *b_GenJetpx;
  TBranch        *b_GenJetpy;
  TBranch        *b_GenJetpz;
  TBranch        *b_GenJeteta;
  TBranch        *b_GenJetphi;

  TBranch        *b_NPFjet;
  TBranch        *b_PFHt;
  TBranch        *b_PFMHt;
  TBranch        *b_PFjetEta;
  TBranch        *b_PFjetPhi;
  TBranch        *b_PFjetE;
  TBranch        *b_PFjetPx;
  TBranch        *b_PFjetPy;
  TBranch        *b_PFjetPz;
  TBranch        *b_PFjetPt;
  TBranch        *b_PFjetCharge;

  TBranch        *b_JetPartonId;
  TBranch        *b_JetPartonMother;
  TBranch        *b_JetPartonPx;
  TBranch        *b_JetPartonPy;
  TBranch        *b_JetPartonPz;
  TBranch        *b_JetPartonEt;
  TBranch        *b_JetPartonE;
  TBranch        *b_JetPartonPhi;
  TBranch        *b_JetPartonEta;
  TBranch        *b_JetPartonFlavour;
  TBranch        *b_JetTrackPt;
  TBranch        *b_JetTrackPhi;
  TBranch        *b_JetTrackPhiWeighted;
  TBranch        *b_JetTrackNo;

  //add photons
  TBranch        *b_Nphot;
  TBranch        *b_PhotE;
  TBranch        *b_PhotEt;
  TBranch        *b_Photpt;
  TBranch        *b_Photpx;
  TBranch        *b_Photpy;
  TBranch        *b_Photpz;
  TBranch        *b_Photeta;
  TBranch        *b_Photphi;
  TBranch        *b_PhotTrkIso;
  TBranch        *b_PhotECalIso;
  TBranch        *b_PhotHCalIso;
  TBranch        *b_PhotAllIso;
  TBranch        *b_Phot_isccPhotAssoc;
  TBranch        *b_PhotLooseEM;
  TBranch        *b_PhotLoosePhoton;
  TBranch        *b_PhotTightPhoton;
  TBranch        *b_PhotGenPdgId;
  TBranch        *b_PhotGenMother;
  TBranch        *b_PhotGenPx;
  TBranch        *b_PhotGenPy;
  TBranch        *b_PhotGenPz;
 
  //add electrons
  TBranch        *b_Nelec;
  TBranch        *b_ElecE;
  TBranch        *b_ElecEt;
  TBranch        *b_Elecpt;
  TBranch        *b_Elecpx;
  TBranch        *b_Elecpy;
  TBranch        *b_Elecpz;
  TBranch        *b_Eleceta;
  TBranch        *b_Elecphi;
  TBranch        *b_ElecCharge;
  TBranch        *b_ElecTrkIso;
  TBranch        *b_ElecECalIso;
  TBranch        *b_ElecHCalIso;
  TBranch        *b_ElecAllIso;
  TBranch        *b_ElecTrkChiNorm;
  TBranch        *b_ElecECalIsoDeposit;
  TBranch        *b_ElecHCalIsoDeposit;
  TBranch        *b_ElecIdLoose;
  TBranch        *b_ElecIdTight;
  TBranch        *b_ElecIdRobLoose;
  TBranch        *b_ElecIdRobTight;
  TBranch        *b_ElecChargeMode;
  TBranch        *b_ElecPtMode;
  TBranch        *b_ElecQOverPErrTrkMode;
  TBranch        *b_ElecCaloEnergy;
  TBranch        *b_ElecHOverE;
  TBranch        *b_ElecVx;
  TBranch        *b_ElecVy;
  TBranch        *b_ElecVz;
  TBranch        *b_ElecD0;
  TBranch        *b_ElecDz;
  TBranch        *b_ElecPtTrk;
  TBranch        *b_ElecQOverPErrTrk;
  TBranch        *b_ElecPinTrk;
  TBranch        *b_ElecPoutTrk;
  TBranch        *b_ElecLostHits;
  TBranch        *b_ElecValidHits;
  TBranch        *b_ElecNCluster;
  TBranch        *b_ElecEtaTrk;
  TBranch        *b_ElecPhiTrk;
  TBranch        *b_ElecWidthClusterEta;
  TBranch        *b_ElecWidthClusterPhi;
  TBranch        *b_ElecGenPdgId;
  TBranch        *b_ElecGenMother;
  TBranch        *b_ElecGenPx;
  TBranch        *b_ElecGenPy;
  TBranch        *b_ElecGenPz;
  TBranch        *b_Elec_isccElecAssoc;

  //add muons
  TBranch        *b_Nmuon;
  TBranch        *b_MuonE;
  TBranch        *b_MuonEt;
  TBranch        *b_Muonpt;
  TBranch        *b_Muonpx;
  TBranch        *b_Muonpy;
  TBranch        *b_Muonpz;
  TBranch        *b_Muoneta;
  TBranch        *b_Muonphi;
  TBranch        *b_MuonCharge;
  TBranch        *b_MuonTrkIso;
  TBranch        *b_MuonECalIso;
  TBranch        *b_MuonHCalIso;
  TBranch        *b_MuonAllIso;
  TBranch        *b_MuonTrkChiNorm;
  TBranch        *b_MuonECalIsoDeposit;
  TBranch        *b_MuonHCalIsoDeposit;
  TBranch        *b_MuonIsGlobal;
  TBranch        *b_MuonIsStandAlone;
  TBranch        *b_MuonIsGlobalTight;
  TBranch        *b_MuonIsTMLastStationLoose;
  TBranch        *b_MuonIsTracker;
  TBranch        *b_MuonIsTMLastStationTight;
  TBranch        *b_MuonIsTM2DCompatibilityLoose;
  TBranch        *b_MuonIsTM2DCompatibilityTight;
  //  TBranch        *b_MuonId;
  TBranch        *b_MuonCombChi2;
  TBranch        *b_MuonCombNdof;
  TBranch        *b_MuonCombVx;
  TBranch        *b_MuonCombVy;
  TBranch        *b_MuonCombVz;
  TBranch        *b_MuonCombD0;
  TBranch        *b_MuonCombDz;
  TBranch        *b_MuonStandValidHits;
  TBranch        *b_MuonStandLostHits;
  TBranch        *b_MuonStandPt;
  TBranch        *b_MuonStandPz;
  TBranch        *b_MuonStandP;
  TBranch        *b_MuonStandEta;
  TBranch        *b_MuonStandPhi;
  TBranch        *b_MuonStandCharge;
  TBranch        *b_MuonStandChi;
  TBranch        *b_MuonStandQOverPError;
  TBranch        *b_MuonTrkValidHits;
  TBranch        *b_MuonTrkLostHits;
  TBranch        *b_MuonTrkD0;
  TBranch        *b_MuonTrkPt;
  TBranch        *b_MuonTrkPz;
  TBranch        *b_MuonTrkP;
  TBranch        *b_MuonTrkEta;
  TBranch        *b_MuonTrkPhi;
  TBranch        *b_MuonTrkCharge;
  TBranch        *b_MuonTrkChi;
  TBranch        *b_MuonTrkQOverPError;
  TBranch        *b_MuonTrkOuterZ;
  TBranch        *b_MuonTrkOuterR;
  TBranch        *b_MuonGenMother;
  TBranch        *b_MuonGenPx;
  TBranch        *b_MuonGenPy;
  TBranch        *b_MuonGenPz;
  TBranch        *b_Muon_isccMuonAssoc;

  TBranch        *b_AlpPtScale;
  TBranch        *b_AlpIdTest;

  TBranch        *b_MPTPhi;
  TBranch        *b_MPTPx;
  TBranch        *b_MPTPy;
  TBranch        *b_MPTPz;

  allData2(TTree *tree=0, std::string outfilename="outfile.root", double=1, double=1, double=1);
  virtual ~allData2();

  std::pair<double, double> relpt_mindr(const std::vector<TLorentzVector>& vJetCollection, const TLorentzVector& vLepton);

  bool isPrompt(long pdgid);

  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  Double_t luminosity_, cross_section_, efficiency_;
  std::string outfilename_;
  std::string infilename_;

  Double_t cut_jet1et, cut_jet2et, cut_jet1eta, cut_jet2eta, cut_jet1phi, cut_jet2phi;
  Double_t cut_jet1emfrac[2], cut_jet2emfrac[2], cut_jet12dphi, cut_jet1metdphi;
  Double_t cut_jet2metdphi, cut_jetemfrac[2], cut_ht, cut_mht, cut_met, cut_meff;
  Double_t cut_alljetet, cut_alljeteta, cut_minjetet, cut_njet;
  Double_t cut_elecet, cut_muonet;
  //Double_t cut_nelec, cut_nmuon;
};

#endif

#ifdef allData2_cxx
allData2::allData2(TTree *tree, std::string outfilename, double luminosity, double cross_section, double efficiency)
{
  outfilename_   = outfilename;
  luminosity_    = luminosity;
  cross_section_ = cross_section;
  efficiency_    = efficiency;

  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SUSY_LM0_PATtified.root");
    if (!f) {
      f = new TFile("SUSY_LM0_PATtified.root");
      if (f) std::cout << "file "<<"SUSY_LM0_PATtified.root"<<" opened\n";
    }
    tree = (TTree*)gDirectory->Get("dijet/allData;");
    if (tree) std::cout << "tree opened\n";

  }
  Init(tree);
}

allData2::~allData2()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t allData2::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t allData2::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void allData2::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (!tree) {std::cout << "tree pointer null\n"; return; }
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("run",   &run,   &b_run);
  fChain->SetBranchAddress("event", &event, &b_event);

  fChain->SetBranchAddress("nHLT",       &nHLT,       &b_nHLT);
  fChain->SetBranchAddress("HLTArray",   HLTArray,    &b_HLTArray);
  fChain->SetBranchAddress("HLT1JET",    &HLT1JET,    &b_HLT1JET);
  fChain->SetBranchAddress("HLT2JET",    &HLT2JET,    &b_HLT2JET);
  fChain->SetBranchAddress("HLT1MET1HT", &HLT1MET1HT, &b_HLT1MET1HT);
  fChain->SetBranchAddress("HLT1MUON",   &HLT1MUON,   &b_HLT1MUON);

  fChain->SetBranchAddress("genN",      &genN,     &b_genN);
  fChain->SetBranchAddress("genid",     genid,     &b_genid);
  fChain->SetBranchAddress("genMother", genMother, &b_genMother);
  fChain->SetBranchAddress("genE",      genE,      &b_genE);
  fChain->SetBranchAddress("genPx",     genPx,     &b_genPx);
  fChain->SetBranchAddress("genPy",     genPy,     &b_genPy);
  fChain->SetBranchAddress("genPz",     genPz,     &b_genPz);
  fChain->SetBranchAddress("genStatus", genStatus, &b_genStatus);

  fChain->SetBranchAddress("genLepN",     &genLepN,    &b_genLepN);
  fChain->SetBranchAddress("genLepId",    genLepId,    &b_genLepId);
  fChain->SetBranchAddress("genLepMother",genLepMother,&b_genLepMother);
  fChain->SetBranchAddress("genLepE",     genLepE,     &b_genLepE);
  fChain->SetBranchAddress("genLepPx",    genLepPx,    &b_genLepPx);
  fChain->SetBranchAddress("genLepPy",    genLepPy,    &b_genLepPy);
  fChain->SetBranchAddress("genLepPz",    genLepPz,    &b_genLepPz);
  fChain->SetBranchAddress("genLepStatus",genLepStatus,&b_genLepStatus);
  
  //vertex information
  fChain->SetBranchAddress("nVtx",                &nVtx,             &b_nVtx);
  fChain->SetBranchAddress("VertexChi2",          VtxChi2,           &b_VertexChi2);
  fChain->SetBranchAddress("VertexNdof",          VtxNdof,           &b_VertexNdof);
  fChain->SetBranchAddress("VertexNormalizedChi2",VtxNormalizedChi2, &b_VertexNormalizedChi2);
  fChain->SetBranchAddress("VertexX",             VtxX,              &b_VertexX);
  fChain->SetBranchAddress("VertexY",             VtxY,              &b_VertexY);
  fChain->SetBranchAddress("VertexZ",             VtxZ,              &b_VertexZ);
  fChain->SetBranchAddress("VertexdX",            VtxdX,             &b_VertexdX);
  fChain->SetBranchAddress("VertexdY",            VtxdY,             &b_VertexdY);
  fChain->SetBranchAddress("VertexdZ",            VtxdZ,             &b_VertexdZ);

  fChain->SetBranchAddress("nFullMET",             &nFullMET,             &b_nFullMET);
  fChain->SetBranchAddress("MET_fullcorr_nocc",    MET_fullcorr_nocc,     &b_MET_fullcorr_nocc);
  fChain->SetBranchAddress("METphi_fullcorr_nocc", &METphi_fullcorr_nocc, &b_METphi_fullcorr_nocc);
  fChain->SetBranchAddress("MET_fullcorr_cc",      MET_fullcorr_cc,       &b_MET_fullcorr_cc);
  fChain->SetBranchAddress("nUncorrMET",           &nUncorrMET,           &b_nUncorrMET);
  fChain->SetBranchAddress("MET_nocorr_nocc",      MET_nocorr_nocc,       &b_MET_nocorr_nocc);
  //fChain->SetBranchAddress("METphi_nocorr_nocc", &METphi_nocorr_nocc, &b_METphi_nocorr_nocc);
  fChain->SetBranchAddress("MET_muoncorr_nocc",    MET_muoncorr_nocc,     &b_MET_muoncorr_nocc);
  fChain->SetBranchAddress("METphi_muoncorr_nocc", &METphi_muoncorr_nocc, &b_METphi_muoncorr_nocc);
  fChain->SetBranchAddress("MET_jeccorr_nocc",     MET_jeccorr_nocc,      &b_MET_jeccorr_nocc);
  fChain->SetBranchAddress("METphi_jeccorr_nocc",  &METphi_jeccorr_nocc,  &b_METphi_jeccorr_nocc);
  fChain->SetBranchAddress("GenMET",               GenMET,                &b_MET_Gen);
  fChain->SetBranchAddress("MET_Fullcorr_nocc_significance", &MET_Fullcorr_nocc_significance, &b_MET_Fullcorr_nocc_significance);

  fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
  fChain->SetBranchAddress("procID",    &procID,    &b_procID);
  fChain->SetBranchAddress("pthat",     &pthat,     &b_pthat);

  fChain->SetBranchAddress("Nhemispheres",  &Nhemispheres, &b_Nhemispheres);
  fChain->SetBranchAddress("HemisphereE",   HemisphereE,   &b_HemisphereE);
  fChain->SetBranchAddress("HemisphereEt",  HemisphereEt,  &b_HemisphereEt);
  fChain->SetBranchAddress("Hemispherept",  Hemispherept,  &b_Hemispherept);
  fChain->SetBranchAddress("Hemispherepx",  Hemispherepx,  &b_Hemispherepx);
  fChain->SetBranchAddress("Hemispherepy",  Hemispherepy,  &b_Hemispherepy);
  fChain->SetBranchAddress("Hemispherepz",  Hemispherepz,  &b_Hemispherepz);
  fChain->SetBranchAddress("Hemisphereeta", Hemisphereeta, &b_Hemisphereeta);
  fChain->SetBranchAddress("Hemispherephi", Hemispherephi, &b_Hemispherephi);
  fChain->SetBranchAddress("Hemispheredphi",Hemispheredphi,&b_Hemispheredphi);

  fChain->SetBranchAddress("Njets",   &Njets,  &b_Njets);
  fChain->SetBranchAddress("Ht",      &Ht,     &b_Ht);
  fChain->SetBranchAddress("MHt",     &MHt,    &b_MHt);
  fChain->SetBranchAddress("JetE",    JetE,    &b_JetE);
  fChain->SetBranchAddress("JetEt",   JetEt,   &b_JetEt);
  fChain->SetBranchAddress("Jetpt",   Jetpt,   &b_Jetpt);
  fChain->SetBranchAddress("Jetpx",   Jetpx,   &b_Jetpx);
  fChain->SetBranchAddress("Jetpy",   Jetpy,   &b_Jetpy);
  fChain->SetBranchAddress("Jetpz",   Jetpz,   &b_Jetpz);
  fChain->SetBranchAddress("Jeteta",  Jeteta,  &b_Jeteta);
  fChain->SetBranchAddress("Jetphi",  Jetphi,  &b_Jetphi);
  fChain->SetBranchAddress("JetFem",  JetFem,  &b_JetFem);
  fChain->SetBranchAddress("JetHemi", JetHemi, &b_JetHemi);
  fChain->SetBranchAddress("Jet_MCcorrFactor",  Jet_MCcorrFactor,  &b_Jet_MCcorrFactor);
  fChain->SetBranchAddress("Jet_JPTcorrFactor", Jet_JPTcorrFactor, &b_Jet_JPTcorrFactor);
  fChain->SetBranchAddress("JetBTag_TkCountHighEff",JetsBTag_TkCountHighEff,&b_JetBTag_TkCountHighEff);
  fChain->SetBranchAddress("JetBTag_SimpleSecVtx",  JetsBTag_SimpleSecVtx,  &b_JetBTag_SimpleSecVtx);
  fChain->SetBranchAddress("JetBTag_CombSecVtx",    JetsBTag_CombSecVtx,    &b_JetBTag_CombSecVtx);
  fChain->SetBranchAddress("Jet_isccJetAssoc",  Jet_isccJetAssoc,  &b_Jet_isccJetAssoc);
  fChain->SetBranchAddress("Jet_ccJetE",  Jet_ccJetE,  &b_Jet_ccJetE);
  fChain->SetBranchAddress("Jet_ccJetpx", Jet_ccJetpx, &b_Jet_ccJetpx);
  fChain->SetBranchAddress("Jet_ccJetpy", Jet_ccJetpy, &b_Jet_ccJetpy);
  fChain->SetBranchAddress("Jet_ccJetpz", Jet_ccJetpz, &b_Jet_ccJetpz);

  fChain->SetBranchAddress("GenHt",    &GenHt,     &b_GenHt);
  fChain->SetBranchAddress("GenHt",    &GenMHt,    &b_GenMHt);
  fChain->SetBranchAddress("GenJetE",  GenJetE,    &b_GenJetE);
  fChain->SetBranchAddress("GenJetEt", GenJetEt,   &b_GenJetEt);
  fChain->SetBranchAddress("GenJetpt", GenJetpt,   &b_GenJetpt);
  fChain->SetBranchAddress("GenJetpx", GenJetpx,   &b_GenJetpx);
  fChain->SetBranchAddress("GenJetpy", GenJetpy,   &b_GenJetpy);
  fChain->SetBranchAddress("GenJetpz", GenJetpz,   &b_GenJetpz);
  fChain->SetBranchAddress("GenJeteta", GenJeteta, &b_GenJeteta);
  fChain->SetBranchAddress("GenJetphi", GenJetphi, &b_GenJetphi);

  fChain->SetBranchAddress("JetPartonId",         JetPartonId,         &b_JetPartonId);
  fChain->SetBranchAddress("JetPartonMother",     JetPartonMother,     &b_JetPartonMother);
  fChain->SetBranchAddress("JetPartonPx",         JetPartonPx,         &b_JetPartonPx);
  fChain->SetBranchAddress("JetPartonPy",         JetPartonPy,         &b_JetPartonPy);
  fChain->SetBranchAddress("JetPartonPz",         JetPartonPz,         &b_JetPartonPz);
  fChain->SetBranchAddress("JetPartonEt",         JetPartonEt,         &b_JetPartonEt);
  fChain->SetBranchAddress("JetPartonE",          JetPartonE,          &b_JetPartonE);
  fChain->SetBranchAddress("JetPartonPhi",        JetPartonPhi,        &b_JetPartonPhi);
  fChain->SetBranchAddress("JetPartonEta",        JetPartonEta,        &b_JetPartonEta);
  fChain->SetBranchAddress("JetPartonFlavour",    JetPartonFlavour,    &b_JetPartonFlavour);
  fChain->SetBranchAddress("JetTrackPt",          JetTrackPt,          &b_JetTrackPt);
  fChain->SetBranchAddress("JetTrackPhi",         JetTrackPhi,         &b_JetTrackPhi);
  fChain->SetBranchAddress("JetTrackPhiWeighted", JetTrackPhiWeighted, &b_JetTrackPhiWeighted);
  fChain->SetBranchAddress("JetTrackNo",          JetTrackNo,          &b_JetTrackNo);

  fChain->SetBranchAddress("Nphot",      &Nphot,     &b_Nphot);
  fChain->SetBranchAddress("PhotE",      PhotE,      &b_PhotE);
  fChain->SetBranchAddress("PhotEt",     PhotEt,     &b_PhotEt);
  fChain->SetBranchAddress("Photpt",     PhotPt,     &b_Photpt);
  fChain->SetBranchAddress("Photpx",     PhotPx,     &b_Photpx);
  fChain->SetBranchAddress("Photpy",     PhotPy,     &b_Photpy);
  fChain->SetBranchAddress("Photpz",     PhotPz,     &b_Photpz);
  fChain->SetBranchAddress("Photeta",    PhotEta,    &b_Photeta);
  fChain->SetBranchAddress("Photphi",    PhotPhi,    &b_Photphi);
  fChain->SetBranchAddress("PhotTrkIso", PhotTrkIso, &b_PhotTrkIso);
  fChain->SetBranchAddress("PhotECalIso",PhotECalIso,&b_PhotECalIso);
  fChain->SetBranchAddress("PhotHCalIso",PhotHCalIso,&b_PhotHCalIso);
  fChain->SetBranchAddress("PhotAllIso", PhotAllIso, &b_PhotAllIso);

  fChain->SetBranchAddress("Phot_isccPhotAssoc",Phot_isccPhotAssoc,    &b_Phot_isccPhotAssoc);
  fChain->SetBranchAddress("PhotLooseEM",       PhotLooseEM,           &b_PhotLooseEM);
  fChain->SetBranchAddress("PhotLoosePhoton",   PhotLoosePhoton,       &b_PhotLoosePhoton);
  fChain->SetBranchAddress("PhotTightPhoton",   PhotTightPhoton,       &b_PhotTightPhoton);
  fChain->SetBranchAddress("PhotGenPdgId",      GenPhotPdgId,          &b_PhotGenPdgId);
  fChain->SetBranchAddress("PhotGenMother",     GenPhotMother,         &b_PhotGenMother);
  fChain->SetBranchAddress("PhotGenPx",         GenPhotPx,             &b_PhotGenPx);
  fChain->SetBranchAddress("PhotGenPy",         GenPhotPy,             &b_PhotGenPy);
  fChain->SetBranchAddress("PhotGenPz",         GenPhotPz,             &b_PhotGenPz);
 
  //add electrons
  fChain->SetBranchAddress("Nelec",  &Nelec, &b_Nelec);
  fChain->SetBranchAddress("ElecE",  ElecE,  &b_ElecE);
  fChain->SetBranchAddress("ElecEt", ElecEt, &b_ElecEt);
  fChain->SetBranchAddress("Elecpt", ElecPt, &b_Elecpt);
  fChain->SetBranchAddress("Elecpx", ElecPx, &b_Elecpx);
  fChain->SetBranchAddress("Elecpy", ElecPy, &b_Elecpy);
  fChain->SetBranchAddress("Elecpz", ElecPz, &b_Elecpz);
  fChain->SetBranchAddress("Eleceta",ElecEta,&b_Eleceta);
  fChain->SetBranchAddress("Elecphi",ElecPhi,&b_Elecphi);

  fChain->SetBranchAddress("ElecCharge",    ElecCharge,   &b_ElecCharge);
  fChain->SetBranchAddress("ElecTrkIso",    ElecTrkIso,   &b_ElecTrkIso);
  fChain->SetBranchAddress("ElecECalIso",   ElecECalIso,  &b_ElecECalIso);
  fChain->SetBranchAddress("ElecHCalIso",   ElecHCalIso,  &b_ElecHCalIso);
  fChain->SetBranchAddress("ElecAllIso",    ElecAllIso,   &b_ElecAllIso);
  fChain->SetBranchAddress("ElecTrkChiNorm",ElecNormChi2 ,&b_ElecTrkChiNorm);

  fChain->SetBranchAddress("ElecECalIsoDeposit", ElecECalIsoDeposit,&b_ElecECalIsoDeposit);
  fChain->SetBranchAddress("ElecHCalIsoDeposit", ElecHCalIsoDeposit,&b_ElecHCalIsoDeposit);

  //MICHELE
  fChain->SetBranchAddress("ElecIdLoose",   ElecIdLoose,   &b_ElecIdLoose );
  fChain->SetBranchAddress("ElecIdTight",   ElecIdTight,   &b_ElecIdTight );
  fChain->SetBranchAddress("ElecIdRobLoose",ElecIdRobLoose,&b_ElecIdRobLoose );
  fChain->SetBranchAddress("ElecIdRobTight",ElecIdRobTight,&b_ElecIdRobTight );
  fChain->SetBranchAddress("ElecChargeMode",ElecChargeMode,&b_ElecChargeMode );
  fChain->SetBranchAddress("ElecPtMode",    ElecPtTrkMode, &b_ElecPtMode );

  fChain->SetBranchAddress("ElecQOverPErrTrkMode",ElecQOverPErrTrkMode,&b_ElecQOverPErrTrkMode );
  fChain->SetBranchAddress("ElecCaloEnergy",      ElecCaloEnergy,      &b_ElecCaloEnergy);

  fChain->SetBranchAddress("ElecHOverE",ElecHOverE,&b_ElecHOverE);
  fChain->SetBranchAddress("ElecVx",    ElecVx,    &b_ElecVx);
  fChain->SetBranchAddress("ElecVy",    ElecVy,    &b_ElecVy);
  fChain->SetBranchAddress("ElecVz",    ElecVz,    &b_ElecVz);
  fChain->SetBranchAddress("ElecD0",    ElecD0,    &b_ElecD0);
  fChain->SetBranchAddress("ElecDz",    ElecDz,    &b_ElecDz);
  fChain->SetBranchAddress("ElecPtTrk", ElecPtTrk, &b_ElecPtTrk);

  fChain->SetBranchAddress("ElecQOverPErrTrk", ElecQOverPErrTrk,&b_ElecQOverPErrTrk);
  fChain->SetBranchAddress("ElecPinTrk",       ElecPinTrk,      &b_ElecPinTrk);
  fChain->SetBranchAddress("ElecPoutTrk",      ElecPoutTrk,     &b_ElecPoutTrk); 
  fChain->SetBranchAddress("ElecLostHits",     ElecLostHits,    &b_ElecLostHits); 
  fChain->SetBranchAddress("ElecValidHits",    ElecValidHits,   &b_ElecValidHits); 
  fChain->SetBranchAddress("ElecNCluster",     ElecNCluster,    &b_ElecNCluster); 
  fChain->SetBranchAddress("ElecEtaTrk",       ElecEtaTrk,      &b_ElecEtaTrk); 
  fChain->SetBranchAddress("ElecPhiTrk",       ElecPhiTrk,      &b_ElecPhiTrk); 

  fChain->SetBranchAddress("ElecWidthClusterEta",ElecWidthClusterEta,&b_ElecWidthClusterEta); 
  fChain->SetBranchAddress("ElecWidthClusterPhi",ElecWidthClusterPhi,&b_ElecWidthClusterPhi); 

  fChain->SetBranchAddress("ElecGenPdgId", GenElecPdgId, &b_ElecGenPdgId);
  fChain->SetBranchAddress("ElecGenMother",GenElecMother,&b_ElecGenMother);
  fChain->SetBranchAddress("ElecGenPx",    GenElecPx,    &b_ElecGenPx);
  fChain->SetBranchAddress("ElecGenPy",    GenElecPy,    &b_ElecGenPy);
  fChain->SetBranchAddress("ElecGenPz",    GenElecPz,    &b_ElecGenPz);
  //PIOPPI
  fChain->SetBranchAddress("Elec_isccElecAssoc",Elec_isccElecAssoc,&b_Elec_isccElecAssoc);

  //add muons
  fChain->SetBranchAddress("Nmuon",         &Nmuon,        &b_Nmuon);  
  fChain->SetBranchAddress("MuonE",         MuonE,         &b_MuonE);
  fChain->SetBranchAddress("MuonEt",        MuonEt,        &b_MuonEt);
  fChain->SetBranchAddress("Muonpt",        MuonPt,        &b_Muonpt);
  fChain->SetBranchAddress("Muonpx",        MuonPx,        &b_Muonpx);
  fChain->SetBranchAddress("Muonpy",        MuonPy,        &b_Muonpy);
  fChain->SetBranchAddress("Muonpz",        MuonPz,        &b_Muonpz);
  fChain->SetBranchAddress("Muoneta",       MuonEta,       &b_Muoneta);
  fChain->SetBranchAddress("Muonphi",       MuonPhi,       &b_Muonphi);
  fChain->SetBranchAddress("MuonCharge",    MuonCharge,    &b_MuonCharge);
  fChain->SetBranchAddress("MuonTrkIso",    MuonTrkIso,    &b_MuonTrkIso);
  fChain->SetBranchAddress("MuonECalIso",   MuonECalIso,   &b_MuonECalIso);
  fChain->SetBranchAddress("MuonHCalIso",   MuonHCalIso,   &b_MuonHCalIso);
  fChain->SetBranchAddress("MuonAllIso",    MuonAllIso,    &b_MuonAllIso);
  fChain->SetBranchAddress("MuonTrkChiNorm",MuonTrkChiNorm,&b_MuonTrkChiNorm);

  fChain->SetBranchAddress("MuonECalIsoDeposit", MuonECalIsoDeposit,&b_MuonECalIsoDeposit);
  fChain->SetBranchAddress("MuonHCalIsoDeposit", MuonHCalIsoDeposit,&b_MuonHCalIsoDeposit);

  fChain->SetBranchAddress("MuonIsGlobal",                MuonIsGlobal,              &b_MuonIsGlobal);
  fChain->SetBranchAddress("MuonIsStandAlone",            MuonIsStandAlone,          &b_MuonIsStandAlone);
  fChain->SetBranchAddress("MuonIsGlobalTight",           MuonIsGlobalTight,         &b_MuonIsGlobalTight);
  fChain->SetBranchAddress("MuonIsTMLastStationLoose",    MuonIsTMLastStationLoose,  &b_MuonIsTMLastStationLoose);
  fChain->SetBranchAddress("MuonIsTracker",               MuonIsTracker,             &b_MuonIsTracker);
  fChain->SetBranchAddress("MuonIsTMLastStationTight",    MuonTMLastStationTight,    &b_MuonIsTMLastStationTight);
  fChain->SetBranchAddress("MuonIsTM2DCompatibilityLoose",MuonTM2DCompatibilityLoose,&b_MuonIsTM2DCompatibilityLoose);
  fChain->SetBranchAddress("MuonIsTM2DCompatibilityTight",MuonTM2DCompatibilityTight,&b_MuonIsTM2DCompatibilityTight);

  //MICHELE
  //  fChain->SetBranchAddress("MuonId",MuonId,&b_MuonId);
  fChain->SetBranchAddress("MuonCombChi2",MuonCombChi2,&b_MuonCombChi2);
  fChain->SetBranchAddress("MuonCombNdof",MuonCombNdof,&b_MuonCombNdof);
  fChain->SetBranchAddress("MuonCombVx",  MuonCombVx,  &b_MuonCombVx);
  fChain->SetBranchAddress("MuonCombVy",  MuonCombVy,  &b_MuonCombVy);
  fChain->SetBranchAddress("MuonCombVz",  MuonCombVz,  &b_MuonCombVz);
  fChain->SetBranchAddress("MuonCombD0",  MuonCombD0,  &b_MuonCombD0);
  fChain->SetBranchAddress("MuonCombDz",  MuonCombDz,  &b_MuonCombDz);

  fChain->SetBranchAddress("MuonStandValidHits",  MuonStandValidHits,  &b_MuonStandValidHits);
  fChain->SetBranchAddress("MuonStandLostHits",   MuonStandLostHits,   &b_MuonStandLostHits);
  fChain->SetBranchAddress("MuonStandPt",         MuonStandPt,         &b_MuonStandPt);
  fChain->SetBranchAddress("MuonStandPz",         MuonStandPz,         &b_MuonStandPz);
  fChain->SetBranchAddress("MuonStandP",          MuonStandP,          &b_MuonStandP);
  fChain->SetBranchAddress("MuonStandEta",        MuonStandEta,        &b_MuonStandEta);
  fChain->SetBranchAddress("MuonStandPhi",        MuonStandPhi,        &b_MuonStandPhi);
  fChain->SetBranchAddress("MuonStandCharge",     MuonStandCharge,     &b_MuonStandCharge);
  fChain->SetBranchAddress("MuonStandChi",        MuonStandChi,        &b_MuonStandChi);
  fChain->SetBranchAddress("MuonStandQOverPError",MuonStandQOverPError,&b_MuonStandQOverPError);

  fChain->SetBranchAddress("MuonTrkValidHits",  MuonTrkValidHits,  &b_MuonTrkValidHits);
  fChain->SetBranchAddress("MuonTrkLostHits",   MuonTrkLostHits,   &b_MuonTrkLostHits);
  fChain->SetBranchAddress("MuonTrkD0",         MuonTrkD0,         &b_MuonTrkD0);
  fChain->SetBranchAddress("MuonTrkPt",         MuonTrkPt,         &b_MuonTrkPt);
  fChain->SetBranchAddress("MuonTrkPz",         MuonTrkPz,         &b_MuonTrkPz);
  fChain->SetBranchAddress("MuonTrkP",          MuonTrkP,          &b_MuonTrkP);
  fChain->SetBranchAddress("MuonTrkEta",        MuonTrkEta,        &b_MuonTrkEta);
  fChain->SetBranchAddress("MuonTrkPhi",        MuonTrkPhi,        &b_MuonTrkPhi);
  fChain->SetBranchAddress("MuonTrkCharge",     MuonTrkCharge,     &b_MuonTrkCharge);
  fChain->SetBranchAddress("MuonTrkChi",        MuonTrkChi,        &b_MuonTrkChi);
  fChain->SetBranchAddress("MuonTrkQOverPError",MuonTrkQOverPError,&b_MuonTrkQOverPError); 
  fChain->SetBranchAddress("MuonTrkOuterZ",     MuonTrkOuterZ,     &b_MuonTrkOuterZ);
  fChain->SetBranchAddress("MuonTrkOuterR",     MuonTrkOuterR,     &b_MuonTrkOuterR);

  fChain->SetBranchAddress("MuonGenMother",GenMuonMother,&b_MuonGenMother);
  fChain->SetBranchAddress("MuonGenPx",    GenMuonPx,    &b_MuonGenPx);
  fChain->SetBranchAddress("MuonGenPy",    GenMuonPy,    &b_MuonGenPy);
  fChain->SetBranchAddress("MuonGenPz",    GenMuonPz,    &b_MuonGenPz);

  //PIOPPI
  fChain->SetBranchAddress("Muon_isccMuonAssoc",Muon_isccMuonAssoc,&b_Muon_isccMuonAssoc);

  //benedetta : PFjets
  fChain->SetBranchAddress("NPFjet",     &NPFjet,     &b_NPFjet);
  fChain->SetBranchAddress("PFHt",       &PFHt,       &b_PFHt);
  fChain->SetBranchAddress("PFMHt",      &PFMHt,      &b_PFMHt);
  fChain->SetBranchAddress("PFjetEta",   &PFjetEta,   &b_PFjetEta);
  fChain->SetBranchAddress("PFjetPhi",   &PFjetPhi,   &b_PFjetPhi);
  fChain->SetBranchAddress("PFjetE",     &PFjetE,     &b_PFjetE);
  fChain->SetBranchAddress("PFjetPx",    &PFjetPx,    &b_PFjetPx);
  fChain->SetBranchAddress("PFjetPy",    &PFjetPy,    &b_PFjetPy);
  fChain->SetBranchAddress("PFjetPz",    &PFjetPz,    &b_PFjetPz);
  fChain->SetBranchAddress("PFjetPt",    &PFjetPt,    &b_PFjetPt);
  fChain->SetBranchAddress("PFjetCharge",&PFjetCharge,&b_PFjetCharge);

  fChain->SetBranchAddress("AlpPtScale", &AlpPtScale, &b_AlpPtScale);
  fChain->SetBranchAddress("AlpIdTest",  &AlpIdTest,  &b_AlpIdTest);

  fChain->SetBranchAddress("MPTPhi", &MPTPhi, &b_MPTPhi);
  fChain->SetBranchAddress("MPTPx",  &MPTPx,  &b_MPTPx);
  fChain->SetBranchAddress("MPTPy",  &MPTPy,  &b_MPTPy);
  fChain->SetBranchAddress("MPTPz",  &MPTPz,  &b_MPTPz);
  Notify();
}

Bool_t allData2::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void allData2::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

std::pair<double, double> allData2::relpt_mindr(const std::vector<TLorentzVector>& vJetCollection,
                                               const TLorentzVector& vLepton) {
  double relpt = 0.;
  double mindr = 99.;

  for (unsigned int ijet = 0; ijet < vJetCollection.size(); ++ijet) {
    if ( vJetCollection[ijet].DeltaR(vLepton) < mindr ) {
      if (vLepton.P() > 0.01 && vJetCollection[ijet].P() > 0.01) 
        relpt = ( (vLepton - vJetCollection[ijet]*vJetCollection[ijet].Dot(vLepton)) * (1./(vLepton.P()*vJetCollection[ijet].P())) ).P();
      mindr = vJetCollection[ijet].DeltaR(vLepton);
    }
  } // jets
 
  std::pair<double,double> thepair(relpt, mindr);
  return thepair;
}

bool allData2::isPrompt(long pdgid) {
  pdgid++;
  //  return indexToPromptBool[ abs(pdgToIndex[pdgid]) ];
  return true;  
}


#endif // #ifdef allData2_cxx

