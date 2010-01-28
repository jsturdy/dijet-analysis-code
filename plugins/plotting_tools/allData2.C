#define allData2_cxx
#include "allData2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <math.h>
#include <algorithm>

#define NUMHISTOS 23

void allData2::Loop() {
  if (fChain == 0) return;

  TFile *file = new TFile(outfilename_.c_str(),"RECREATE");
  file->cd();

  Long64_t nentries = fChain->GetEntries();
  std::cout << "GetEntries() = " << fChain->GetEntries() << std::endl;

  fChain->SetBranchStatus("*",0);  // disable all branches

  //vertex information
  fChain->SetBranchStatus("nVtx",1);
  fChain->SetBranchStatus("VertexChi2",1);
  fChain->SetBranchStatus("VertexNdof",1);
  fChain->SetBranchStatus("VertexNormalizedChi2",1);
  fChain->SetBranchStatus("VertexX",1);
  fChain->SetBranchStatus("VertexY",1);
  fChain->SetBranchStatus("VertexZ",1);
  fChain->SetBranchStatus("VertexdX",1);
  fChain->SetBranchStatus("VertexdY",1);
  fChain->SetBranchStatus("VertexdZ",1);

  //jet information
  fChain->SetBranchStatus("Njets",1);
  fChain->SetBranchStatus("Jetpx",1);
  fChain->SetBranchStatus("Jetpy",1);
  fChain->SetBranchStatus("Jetpz",1);
  fChain->SetBranchStatus("JetE",1);
  fChain->SetBranchStatus("JetEt",1);
  fChain->SetBranchStatus("Jetpt",1);
  fChain->SetBranchStatus("Jetphi",1);
  fChain->SetBranchStatus("Jeteta",1);
  fChain->SetBranchStatus("JetFem",1);
  fChain->SetBranchStatus("Jet_isccJetAssoc",1);
  fChain->SetBranchStatus("Jet_ccJetE",1);
  fChain->SetBranchStatus("Jet_ccJetpx",1);
  fChain->SetBranchStatus("Jet_ccJetpy",1);
  fChain->SetBranchStatus("Jet_ccJetpz",1);
  fChain->SetBranchStatus("Jet_MCcorrFactor",1);
  fChain->SetBranchStatus("Jet_JPTcorrFactor",1);
  fChain->SetBranchStatus("JetBTag_TkCountHighEff",1);
  fChain->SetBranchStatus("JetBTag_SimpleSecVtx",1);
  fChain->SetBranchStatus("JetBTag_CombSecVtx",1);
  fChain->SetBranchStatus("Ht",1);
  fChain->SetBranchStatus("MHt",1);

  fChain->SetBranchStatus("nFullMET",1);
  fChain->SetBranchStatus("nUncorrMET",1);
  fChain->SetBranchStatus("MET_fullcorr_cc",1);
  fChain->SetBranchStatus("MET_fullcorr_nocc",1);
  fChain->SetBranchStatus("MET_jeccorr_nocc",1);
  fChain->SetBranchStatus("MET_muoncorr_nocc",1);
  fChain->SetBranchStatus("MET_nocorr_nocc",1);
  fChain->SetBranchStatus("METphi_fullcorr_nocc",1);
  fChain->SetBranchStatus("METphi_jeccorr_nocc",1);
  fChain->SetBranchStatus("METphi_muoncorr_nocc",1);
  fChain->SetBranchStatus("MPTPhi",1);
  fChain->SetBranchStatus("MPTPx",1);
  fChain->SetBranchStatus("MPTPy",1);
  fChain->SetBranchStatus("MPTPz",1);

  //pf jet information
  fChain->SetBranchStatus("NPFjet",1);
  fChain->SetBranchStatus("PFHt",1);
  fChain->SetBranchStatus("PFMHt",1);
  fChain->SetBranchStatus("PFjetEta",1);
  fChain->SetBranchStatus("PFjetPhi",1);
  fChain->SetBranchStatus("PFjetE",1);
  fChain->SetBranchStatus("PFjetPx",1);
  fChain->SetBranchStatus("PFjetPy",1);
  fChain->SetBranchStatus("PFjetPz",1);
  fChain->SetBranchStatus("PFjetPt",1);
  fChain->SetBranchStatus("PFjetCharge",1);

  //electron information
  fChain->SetBranchStatus("Nelec",1);
  fChain->SetBranchStatus("ElecE",1);
  fChain->SetBranchStatus("ElecEt",1);
  fChain->SetBranchStatus("Elecpt",1);
  fChain->SetBranchStatus("Elecpx",1);
  fChain->SetBranchStatus("Elecpy",1);
  fChain->SetBranchStatus("Elecpz",1);
  fChain->SetBranchStatus("Eleceta",1);
  fChain->SetBranchStatus("Elecphi",1);
  fChain->SetBranchStatus("ElecCharge",1);

  //muon information
  fChain->SetBranchStatus("Nmuon",1);
  fChain->SetBranchStatus("MuonE",1);
  fChain->SetBranchStatus("MuonEt",1);
  fChain->SetBranchStatus("Muonpt",1);
  fChain->SetBranchStatus("Muonpx",1);
  fChain->SetBranchStatus("Muonpy",1);
  fChain->SetBranchStatus("Muonpz",1);
  fChain->SetBranchStatus("Muoneta",1);
  fChain->SetBranchStatus("Muonphi",1);
  fChain->SetBranchStatus("MuonCharge",1);

  fChain->SetBranchStatus("AlpPtScale",1);
  fChain->SetBranchStatus("AlpIdTest",1);


  double localpi  = acos(-1);
  TH1F *h_selections[2];
  TH1F *h_Nelec[4], *h_Nmuon[4];

  TH1F *h_Njets[4][2], *h_jet1eta[4], *h_jet2eta[4];
  TH1F *h_elecEta[4], *h_muonEta[4];

  TH1F *h_jet1emfrac[4], *h_jet2emfrac[4], *h_jetFem[4];
  TH1F *h_jet12dphi[4], *h_jet1metdphi[4], *h_jet2metdphi[4];
  // bins of 100 GeV
  TH1F *h_HT[4], *h_Meff[4];
  // bins of 50 GeV
  TH1F *h_MET[4], *h_jet1et[4], *h_jet2et[4], *h_jetallet[4];
  TH1F *h_MHT[4], *h_MT[2];
  // bins of 25 and 1 GeV depdnding on plot
  TH1F *h_elecEt[4], *h_muonEt[4];
  //only plot these for pre cuts and post cuts
  TH1F *h_jet1phi[2], *h_jet2phi[2], *h_METphi[2];

  char histtitle[NUMHISTOS][128];
  string histname[NUMHISTOS] = {"01_jet1et","02_jet2et","03_jetallet",
				"11_MET","12_HT","13_MHT","14_Meff",
				"21_jet1metdphi","22_jet2metdphi","23_jet12dphi",
				"30_Njets","31_Ngoodjets","32_jet1eta","33_jet2eta",
				"41_jetFem","42_jet1emfrac","43_jet2emfrac",
				"51_Nelecs","52_Nmuons",
				"53_eleceta","54_muoneta",
				"55_elecet","56_muonet"};

  string histpre[4] = {"h_pre_cuts_","h_individual_cuts_","h_N1_cuts_","h_post_cuts_"};
                            
  double bins[4][NUMHISTOS] = {
    //pre cuts
  // j1et,  j2et,  jallet, met,   ht ,   mht,   meff
    {2500., 2500., 1500.,  5000., 5000., 2000., 5000., 
  // j1mdp,   j2mdp,   j12dp,   njets, ngood
     localpi, localpi, localpi, 20.,   20.,
  // jetetmultiplier, dummy values
     1.,   0.,   0.,   0.,   0.,
  // nelec, nmuon, dum6, dum7, elecet, muonet
     20.,   50.,   0.,   0.,   2000.,  5000.},
    
    //individual cuts
    {2500., 2500., 100.,  5000., 5000., 2000., 5000.,
  // j1mdp,   j2mdp,   j12dp,   njets, ngood
     localpi, localpi, localpi, 20.,   20.,
  // jetetmultiplier, dummy values
     10.,   0.,   0.,   0.,   0.,
  // nelec, nmuon, dum6, dum7, elecet, muonet
     15.,   15.,   0.,   0.,   20.,  20.},
    
    //N-1 cuts
    {2500., 2500., 1500.,  3000., 5000., 2000., 5000.,
  // j1mdp,   j2mdp,   j12dp,   njets, ngood
     localpi, localpi, localpi, 20.,   20.,
  // jetetmultiplier, dummy values
     1.,   0.,   0.,   0.,   0.,
  // nelec, nmuon, dum6, dum7, elecet, muonet
     15.,   30.,   0.,   0.,   1000.,  5000.},
    
    //post cuts
    {2500., 2500., 100.,  3000., 5000., 2000., 5000.,
  // j1mdp,   j2mdp,   j12dp,   njets, ngood
     localpi, localpi, localpi, 20.,   20.,
  // jetetmultiplier, dummy values
     10.,   0.,   0.,   0.,   0.,
  // nelec, nmuon, dum6, dum7, elecet, muonet
     15.,   15.,   0.,   0.,   20.,  20.}};

  double binsize = 0.;

  for (int tt = 0; tt < 4; tt++) {
    for (int hh = 0; hh < NUMHISTOS; hh++) {
      sprintf(histtitle[hh],"%s%s",histpre[tt].c_str(),histname[hh].c_str());
    }


    if (tt==0||tt==2) binsize = 25.;// bins of 25 GeV for pre/N-1
    else binsize = 5.;// bins of 5 GeV for individual/post
    h_jetallet[tt]    = new TH1F(histtitle[2],"",static_cast<int>(bins[tt][12]*bins[tt][2]/binsize),0.,bins[tt][2]);

    binsize = 50.; // bins of 50 GeV

    h_jet1et[tt]  = new TH1F(histtitle[0],"",static_cast<int>(bins[tt][0]/binsize),0.,bins[tt][0]);
    h_jet2et[tt]  = new TH1F(histtitle[1],"",static_cast<int>(bins[tt][1]/binsize),0.,bins[tt][1]);
    h_MET[tt]     = new TH1F(histtitle[3],"",static_cast<int>(bins[tt][3]/binsize),0.,bins[tt][3]);
    h_MHT[tt]     = new TH1F(histtitle[5],"",static_cast<int>(bins[tt][5]/binsize),0.,bins[tt][5]);

    binsize = 100.; // bins of 100 GeV
    h_HT[tt]      = new TH1F(histtitle[4],"",static_cast<int>(bins[tt][4]/binsize),0.,bins[tt][4]);
    h_Meff[tt]    = new TH1F(histtitle[6],"",static_cast<int>(bins[tt][6]/binsize),0.,bins[tt][6]);
    
    // fixed number of bins
    h_jet1metdphi[tt] = new TH1F(histtitle[7],"",50,0.,bins[tt][7]);
    h_jet2metdphi[tt] = new TH1F(histtitle[8],"",50,0.,bins[tt][8]);
    h_jet12dphi[tt]   = new TH1F(histtitle[9],"",50,0.,bins[tt][9]);

    h_Njets[tt][0]  = new TH1F(histtitle[10],"",static_cast<int>(bins[tt][10]),0.,bins[tt][10]);
    h_Njets[tt][1]  = new TH1F(histtitle[11],"",static_cast<int>(bins[tt][11]),0.,bins[tt][11]);
    h_jet1eta[tt]   = new TH1F(histtitle[12],"",100,-5.,5.);
    h_jet2eta[tt]   = new TH1F(histtitle[13],"",100,-5.,5.);

    h_jetFem[tt]     = new TH1F(histtitle[14],"",50,0.,1.);
    h_jet1emfrac[tt] = new TH1F(histtitle[15],"",50,0.,1.);
    h_jet2emfrac[tt] = new TH1F(histtitle[16],"",50,0.,1.);

    //lepton plots
    h_Nelec[tt]   = new TH1F(histtitle[17],"",static_cast<int>(bins[tt][17]),0.,bins[tt][17]);
    h_Nmuon[tt]   = new TH1F(histtitle[18],"",static_cast<int>(bins[tt][18]),0.,bins[tt][18]);
    h_elecEta[tt] = new TH1F(histtitle[19],"",100,-5.,5.);
    h_muonEta[tt] = new TH1F(histtitle[20],"",100,-5.,5.);
    
    if (tt==0||tt==2) binsize = 25.;// bins of 25 GeV for pre/N-1
    else binsize = 1.;// bins of 1 GeV for individual/post
    h_elecEt[tt]   = new TH1F(histtitle[21],"",static_cast<int>(bins[tt][21]/binsize),0.,bins[tt][21]);
    h_muonEt[tt]   = new TH1F(histtitle[22],"",static_cast<int>(bins[tt][22]/binsize),0.,bins[tt][22]);
  }

  //pre cuts
  h_jet1phi[0] = new TH1F("h_pre_cuts_9_jet1phi","",100,-localpi,localpi);
  h_jet2phi[0] = new TH1F("h_pre_cuts_9_jet2phi","",100,-localpi,localpi);
  h_METphi[0]  = new TH1F("h_pre_cuts_9_METphi", "",100,-localpi,localpi);
  h_MT[0]      = new TH1F("h_pre_cuts_9_MT",     "",1500/25,0,1500);
  //post cuts
  h_jet1phi[1] = new TH1F("h_post_cuts_9_jet1phi","",100,-localpi,localpi);
  h_jet2phi[1] = new TH1F("h_post_cuts_9_jet2phi","",100,-localpi,localpi);
  h_METphi[1]  = new TH1F("h_post_cuts_9_METphi", "",100,-localpi,localpi);
  h_MT[1]      = new TH1F("h_post_cuts_9_MT",     "",1500/25,0,1500);

  // individual cuts histos
  h_selections[0]   = new TH1F("h_selections","",19,0,19);
  // N-1 histos
  h_selections[1]   = new TH1F("h_N1_selections","",18,0,18);

  Long64_t nbytes = 0, nb = 0;

  double scale = luminosity_ * cross_section_ * efficiency_ / (double)nentries;

  cut_njet     = 2;  //require at least 2 jets
  cut_jet1et   = 100;//min et of jet 1
  cut_jet2et   = 100;//min et of jet 2
  cut_alljetet = 50; //max et of all other jets

  cut_jet1eta   = 2.5;//max eta of jet 1
  cut_jet2eta   = 2.5;//max eta of jet 2
  cut_alljeteta = 3.0;//max eta of all other jets

  cut_jet1phi  = 0.;
  cut_jet2phi  = 0.;
  
  cut_jet1emfrac[0] = 0.05;//min em fraction of jet 1
  cut_jet1emfrac[1] = 0.95;//max em fraction of jet 1
  cut_jet2emfrac[0] = 0.05;//min em fraction of jet 2
  cut_jet2emfrac[1] = 0.95;//max em fraction of jet 2
  
  cut_jetemfrac[0] = 0.05;//min event em fraction
  cut_jetemfrac[1] = 0.95;//min event em fraction
  cut_jet12dphi    = 0.3; //min dphi(jet1, jet2)
  
  cut_jet1metdphi = 0.6;//min dphi(jet1, met)
  cut_jet2metdphi = 0.4;//min dphi(jet2, met)
  
  cut_ht   = 250.;
  cut_mht  = 200.;
  cut_met  = 200.;
  cut_meff = 0.;

  cut_elecet = 15.0;
  cut_muonet = 15.0;

  // values to go into the selection table
  int  pscounter      = 0;
  int  fjcounter      = 0;
  int  htcounter      = 0;
  int  mhtcounter     = 0;
  int  dphicounter    = 0;
  int  metcounter     = 0;
  int  leptoncounter  = 0;
  //
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb; // nb = number of bytes read

    double METx   = MET_fullcorr_nocc[0];
    double METy   = MET_fullcorr_nocc[1];
    double jet12dphi = 0, jet1metdphi = 0, jet2metdphi = 0;
    double MET    = sqrt( METx*METx + METy*METy );

    jet12dphi   = Jetphi[0]-Jetphi[1];
    jet12dphi   = (jet12dphi<0)?-jet12dphi:jet12dphi;
    jet12dphi   = (jet12dphi>localpi)?2*localpi - jet12dphi:jet12dphi;

    jet1metdphi = Jetphi[0]-METphi_fullcorr_nocc;
    jet1metdphi = (jet1metdphi<0)?-jet1metdphi:jet1metdphi;
    jet1metdphi = (jet1metdphi>localpi)?2*localpi - jet1metdphi:jet1metdphi;

    jet2metdphi = Jetphi[1]-METphi_fullcorr_nocc;
    jet2metdphi = (jet2metdphi<0)?-jet2metdphi:jet2metdphi;
    jet2metdphi = (jet2metdphi>localpi)?2*localpi - jet2metdphi:jet2metdphi;

    //double METphi = atan2(METy,METx);

    //bool nPhotSelection[2]      = {false,false};
    //bool nElecSelection[2]      = {false,false};
    //bool nMuonSelection[2]      = {false,false};
    //bool nTauSelection[2]       = {false,false};

    bool preselection   = true;// jets pT>30GeV, |eta|<3, 0.05<emfrac<0.95
    bool dijetselection = true;// jet1 and jet2 pT>100GeV, |eta|<2.5
    bool finaljet       = true;// no jet pT>50GeV, combine with dijet selection
    bool htselection    = true;// HT of event > 250GeV
    bool mhtselection   = true;// MHT of event > 200GeV
    bool dphiselection  = true;// dphi(jet1, jet2, met) passes cuts
    bool metselection   = true;// MET of event > 200GeV//alternately 350GeV
    bool leptonveto     = true;// true if no leptons have pT>15GeV

    bool nJetSelection[2]       = {false,false};

    bool jet1EtSelection[2]     = {false,false};
    bool jet2EtSelection[2]     = {false,false};
    bool jet1EtaSelection[2]    = {false,false};
    bool jet2EtaSelection[2]    = {false,false};
    bool jet1EMFracSelection[2] = {false,false};
    bool jet2EMFracSelection[2] = {false,false};
    bool jet12dphiSelection[2]  = {false,false};

    //bool jetEMFracSelection[2]  = {false,false};
    bool excessiveJetVeto[2]    = {true,true};

    bool metSelection[2]         = {false,false};
    bool jet1metdphiSelection[2] = {false,false};
    bool jet2metdphiSelection[2] = {false,false};

    bool htSelection[2]   = {false,false};
    bool mhtSelection[2]  = {false,false};
    bool meffSelection[2] = {false,false};
    //bool mptSelection[2] = {false,false};

    bool electronVeto[2]    = {true,true};
    bool muonVeto[2]        = {true,true};
    
    nJetSelection[0]    = (Njets>=cut_njet)      ? true : false;
    

    if (Njets>0) jet1EtSelection[0]  = (Jetpt[0]>=cut_jet1et) ? true : false;
    if (Njets>1) jet2EtSelection[0]  = (Jetpt[1]>=cut_jet2et) ? true : false;

    if (Njets>0) jet1EtaSelection[0] = (fabs(Jeteta[0])<=cut_jet1eta)   ? true : false;
    if (Njets>1) jet2EtaSelection[0] = (fabs(Jeteta[1])<=cut_jet2eta)   ? true : false;

    if (Njets>0) jet1EMFracSelection[0] = ((JetFem[0]>=cut_jet1emfrac[0])&&(JetFem[0]<=cut_jet1emfrac[1])) ? true : false;
    if (Njets>1) jet2EMFracSelection[0] = ((JetFem[1]>=cut_jet2emfrac[0])&&(JetFem[1]<=cut_jet2emfrac[1])) ? true : false;
    
    if (Njets>1) jet12dphiSelection[0]   = (jet12dphi>=cut_jet12dphi)     ? true : false;

    if (Njets>0) jet1metdphiSelection[0] = (jet1metdphi>=cut_jet1metdphi) ? true : false;
    if (Njets>1) jet2metdphiSelection[0] = (jet2metdphi>=cut_jet2metdphi) ? true : false;

    metSelection[0]   = (MET>=cut_met)   ? true : false;
    htSelection[0]    = (Ht>=cut_ht)     ? true : false;
    mhtSelection[0]   = (MHt>=cut_mht)   ? true : false;

    UInt_t goodjetcount = 0;

    if (Njets>0) goodjetcount +=1;
    if (Njets>1) goodjetcount +=1;
    
    for (int ijet = 0; ijet < Njets; ijet++) {
      if (Jeteta[ijet]>cut_alljeteta || JetFem[ijet]<cut_jetemfrac[0] || JetFem[ijet]>cut_jetemfrac[1])
	preselection = false;}

    for (int ijet = 2; ijet < Njets; ijet++) {
      if (Jetpt[ijet]>cut_alljetet) 
	excessiveJetVeto[0] = false;
      else
	goodjetcount++;}

    for (int ielec = 0; ielec < Nelec; ielec++) {
      if (ElecEt[ielec]>cut_elecet) electronVeto[0] = false;
    }
    
    for (int imuon = 0; imuon < Nmuon; imuon++) {
      if (MuonEt[imuon]>cut_muonet) muonVeto[0] = false;
    }

    double Meff   = Ht + MHt;
    double MT     = Jetpt[0]+Jetpt[1];
    meffSelection[0]  = (Meff>=cut_meff) ? true : false;
      

    dijetselection = jet1EtSelection[0]&&jet1EMFracSelection[0]&&jet1EtaSelection[0]&&
                     jet2EtSelection[0]&&jet2EMFracSelection[0]&&jet2EtaSelection[0];
    finaljet       = dijetselection&&excessiveJetVeto[0];
    htselection    = htSelection[0];
    mhtselection   = mhtSelection[0];
    dphiselection  = jet1metdphiSelection[0]&&jet2metdphiSelection[0];
    metselection   = metSelection[0];
    leptonveto     = electronVeto[0]||muonVeto[0];

    if (preselection) {
      pscounter++;
      if (finaljet) {
	fjcounter++;
	if (dphiselection) {
	  dphicounter++;
	  if (htselection) {
	    htcounter++;
	    if (metselection) {
	      metcounter++;
	      if (mhtselection) {
		mhtcounter++;
		if (leptonveto) {
		  leptoncounter++;}}}}}}}
    

    nJetSelection[1] = jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];
    
    jet1EtSelection[1] = nJetSelection[0]&&
      jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];
    
    jet2EtSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];
    
    jet1EtaSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];
    
    jet2EtaSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];
    
    jet1EMFracSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet2EMFracSelection[0]&&jet1metdphiSelection[0]&&
      jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];
    
    jet2EMFracSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];
    
    jet1metdphiSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];
    
    jet2metdphiSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];

    metSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      mhtSelection[0]&&htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];

    htSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];

    mhtSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];

    meffSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];

    excessiveJetVeto[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      jet12dphiSelection[0];//&&
      //electronVeto[0]&&muonVeto[0];

    jet12dphiSelection[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0];//&&
      //electronVeto[0]&&muonVeto[0];

    electronVeto[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&jet12dphiSelection[0]&&muonVeto[0];

    muonVeto[1] = nJetSelection[0]&&
      jet1EtSelection[0]&&jet2EtSelection[0]&&
      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
      metSelection[0]&&mhtSelection[0]htSelection[0]&&meffSelection[0]&&
      excessiveJetVeto[0]&&electronVeto[0]&&jet12dphiSelection[0];


    ////////////////////////////////////
    //pre cut plots
    ////////////////////////////////////
    h_Njets[0][0]->Fill(Njets);
    
    if (Njets>0) {
      h_jet1et[0]->Fill(Jetpt[0]);
      h_jet1eta[0]->Fill(Jeteta[0]);
      h_jet1emfrac[0]->Fill(JetFem[0]);
      h_jet1phi[0]->Fill(Jetphi[0]);}
    if (Njets>1) {
      h_jet2et[0]->Fill(Jetpt[1]);
      h_jet2eta[0]->Fill(Jeteta[1]);
      h_jet2emfrac[0]->Fill(JetFem[1]);
      h_jet2phi[0]->Fill(Jetphi[1]);}
    if (Njets>1) h_MT[0]->Fill(MT);
    int ijet = 2;
    while (ijet < Njets) {
      h_jetFem[0]->Fill(JetFem[ijet]);
      h_jetallet[0]->Fill(Jetpt[ijet]);
      ijet++;}
    
    h_Njets[0][1]->Fill(goodjetcount);

    h_jet1metdphi[0]->Fill(jet1metdphi);
    h_jet2metdphi[0]->Fill(jet2metdphi);
    
    h_MET[0]->Fill(MET);
    h_HT[0]->Fill(Ht);
    h_MHT[0]->Fill(MHt);
    h_Meff[0]->Fill(Meff);

    if (Njets>1) h_jet12dphi[0]->Fill(sqrt((Jetphi[0]-Jetphi[1])*(Jetphi[0]-Jetphi[1])));
        
    h_METphi[0]->Fill(METphi_fullcorr_nocc);

    for (unsigned int ielec = 0; ielec < Nelec; ++ielec) h_elecEt[0]->Fill(ElecEt[ielec]);
    h_Nelec[0]->Fill(Nelec);
    for (unsigned int imuon = 0; imuon < Nmuon; ++imuon) h_muonEt[0]->Fill(MuonEt[imuon]);
    h_Nmuon[0]->Fill(Nmuon);

    ////////////////////////////////////
    //individual cut plots
    ////////////////////////////////////

    if(nJetSelection[0]) {
      h_Njets[1][0]->Fill(Njets);
      h_selections[0]->Fill(0.5);}
    
    if(jet1EtSelection[0]) {
      h_jet1et[1]->Fill(Jetpt[0]);
      h_selections[0]->Fill(1.5);}

    if(jet2EtSelection[0]) {
      h_jet2et[1]->Fill(Jetpt[1]);
      h_selections[0]->Fill(2.5);}

    if(metSelection[0]) {
      h_MET[1]->Fill(MET);
      h_selections[0]->Fill(9.5);}

    if(htSelection[0]) {
      h_HT[1]->Fill(Ht);
      h_selections[0]->Fill(10.5);}

    if(mhtSelection[0]) {
      h_MHT[1]->Fill(MHt);
      h_selections[0]->Fill(11.5);}

    if(meffSelection[0]) {
      h_Meff[1]->Fill(Meff);
      h_selections[0]->Fill(12.5);}

    if(jet1EtaSelection[0]) {
      h_jet1eta[1]->Fill(Jeteta[0]);
      h_selections[0]->Fill(3.5);}

    if(jet2EtaSelection[0]) {
      h_jet2eta[1]->Fill(Jeteta[1]);
      h_selections[0]->Fill(4.5);}

    if(jet1EMFracSelection[0]) {
      h_jet1emfrac[1]->Fill(JetFem[0]);
      h_selections[0]->Fill(5.5);}

    if(jet2EMFracSelection[0]) {
      h_jet2emfrac[1]->Fill(JetFem[1]);
      h_selections[0]->Fill(6.5);}
    

    ijet = 2;
    while (ijet < Njets) {
      if (Jetpt[ijet]>cut_minjetet) h_jetFem[1]->Fill(JetFem[ijet]);
      ijet++;}
    
    if(excessiveJetVeto[0]) {
      ijet = 2;
      while (ijet < Njets) {
	if (Jetpt[ijet]>cut_minjetet) h_jetallet[1]->Fill(Jetpt[ijet]);
	ijet++;}
      h_Njets[1][1]->Fill(goodjetcount);
      h_selections[0]->Fill(14.5);}

    if(electronVeto[0]) {
      //int eleccount = 0;
      for (unsigned int ielec = 0; ielec < Nelec; ++ielec) {
	//if(ElecEt[ielec]>cut_elecet) {
	h_elecEt[1]->Fill(ElecEt[ielec]);}
      //eleccount++;}
      h_Nelec[1]->Fill(Nelec);
      h_selections[1]->Fill(15.5);}
    //h_Nelec[1]->Fill(eleccount);}

    if(muonVeto[0]) {
      //int muoncount = 0;
      for (unsigned int imuon = 0; imuon < Nmuon; ++imuon) {
	//if(MuonEt[imuon]>cut_muonet) {
	h_muonEt[1]->Fill(MuonEt[imuon]);}
      //muoncount++;}
      h_Nmuon[1]->Fill(Nmuon);
      h_selections[0]->Fill(16.5);}
    //h_Nmuon[1]->Fill(muoncount);}


    if(jet1metdphiSelection[0]) {
      h_jet1metdphi[1]->Fill(jet1metdphi);
      h_selections[0]->Fill(7.5);}

    if(jet2metdphiSelection[0]) {
      h_jet2metdphi[1]->Fill(jet2metdphi);
      h_selections[0]->Fill(8.5);}
    
    if(jet12dphiSelection[0]) {
      h_jet12dphi[1]->Fill(sqrt((Jetphi[0]-Jetphi[1])*(Jetphi[0]-Jetphi[1])));
      h_selections[0]->Fill(13.5);}
    
    ///////////////////////////////////////////
    //N-1 plots
    //////////////////////////////////////////

    if(nJetSelection[1]) {
      h_Njets[2][0]->Fill(Njets);
      h_selections[1]->Fill(0.5);}
    
    if(jet1EtSelection[1]) {
      h_jet1et[2]->Fill(Jetpt[0]);
      h_selections[1]->Fill(1.5);}

    if(jet2EtSelection[1]) {
      h_jet2et[2]->Fill(Jetpt[1]);
      h_selections[1]->Fill(2.5);}

    if(jet1EtaSelection[1]) {
      h_jet1eta[2]->Fill(Jeteta[0]);
      h_selections[1]->Fill(3.5);}

    if(jet2EtaSelection[1]) {
      h_jet2eta[2]->Fill(Jeteta[1]);
      h_selections[1]->Fill(4.5);}

    if(jet1EMFracSelection[1]) {
      h_jet1emfrac[2]->Fill(JetFem[0]);
      h_selections[1]->Fill(5.5);}

    if(jet2EMFracSelection[1]) {
      h_jet2emfrac[2]->Fill(JetFem[1]);
      h_selections[1]->Fill(6.5);}
        
    if(jet1metdphiSelection[1]) {
      h_jet1metdphi[2]->Fill(jet1metdphi);
      h_selections[1]->Fill(7.5);}

    if(jet2metdphiSelection[1]) {
      h_jet2metdphi[2]->Fill(jet2metdphi);
      h_selections[1]->Fill(8.5);}
    
    if(metSelection[1]) {
      h_MET[2]->Fill(MET);
      h_selections[1]->Fill(9.5);}

    if(htSelection[1]) {
      h_HT[2]->Fill(Ht);
      h_selections[1]->Fill(10.5);}

    if(mhtSelection[1]) {
      h_MHT[2]->Fill(MHt);
      h_selections[1]->Fill(11.5);}

    if(meffSelection[1]) {
      h_Meff[2]->Fill(Meff);
      h_selections[1]->Fill(12.5);}

    if(jet12dphiSelection[1]) {
      h_jet12dphi[2]->Fill(sqrt((Jetphi[0]-Jetphi[1])*(Jetphi[0]-Jetphi[1])));
      h_selections[1]->Fill(13.5);}
    
    ijet = 2;
    while (ijet < Njets) {
      if (Jetpt[ijet]>cut_minjetet) h_jetFem[2]->Fill(JetFem[ijet]);
      ijet++;}

    if(excessiveJetVeto[1]) {
      ijet = 2;
      while (ijet < Njets) {
	if (Jetpt[ijet]>cut_minjetet) h_jetallet[2]->Fill(Jetpt[ijet]);
	ijet++;}
      h_Njets[2][1]->Fill(goodjetcount);
      //h_HT[2]->Fill(Ht);//what to do about recalculating HT when some jets are rejected?
      h_selections[1]->Fill(14.5);}
    
    if(electronVeto[1]) {
      for (unsigned int ielec = 0; ielec < Nelec; ++ielec) h_elecEt[2]->Fill(ElecEt[ielec]);
      h_Nelec[2]->Fill(Nelec);
      h_selections[1]->Fill(15.5);}
    
    if(muonVeto[1]) {
      for (unsigned int imuon = 0; imuon < Nmuon; ++imuon) h_muonEt[2]->Fill(MuonEt[imuon]);
      h_Nmuon[2]->Fill(Nmuon);
      h_selections[1]->Fill(16.5);}
    
    ////////////////////////////////////////    
    //Full selection
    ////////////////////////////////////////

    bool selections = excessiveJetVeto[0]&&nJetSelection[0]&&
                      jet1EtSelection[0]&&jet2EtSelection[0]&&
                      jet1EtaSelection[0]&&jet2EtaSelection[0]&&
                      jet1EMFracSelection[0]&&jet2EMFracSelection[0]&&
                      jet1metdphiSelection[0]&&jet2metdphiSelection[0]&&
                      jet12dphiSelection[0]&&metSelection[0]&&htSelection[0]&&
                      mhtSelection[0]&&meffSelection[0];
                      //&&electronVeto[0]&&muonVeto[0];

    if (selections) {
      h_selections[0]->Fill(17.5);
      h_Njets[3][0]->Fill(Njets);
      
      h_MT[1]->Fill(MT);

      h_jet1et[3]->Fill(Jetpt[0]);
      h_jet2et[3]->Fill(Jetpt[1]);

      h_MET[3]->Fill(MET);
      h_MHT[3]->Fill(MHt);
      h_HT[3]->Fill(Ht);
      h_Meff[3]->Fill(Meff);

      for (unsigned int ijet = 2; ijet < Njets; ++ijet) {
	if (Jetpt[ijet]>cut_minjetet) {
	  h_jetFem[3]->Fill(JetFem[ijet]);
	  h_jetallet[3]->Fill(Jetpt[ijet]);}}

      h_Njets[3][1]->Fill(goodjetcount);

      h_jet1eta[3]->Fill(Jeteta[0]);
      h_jet2eta[3]->Fill(Jeteta[1]);
      h_jet1emfrac[3]->Fill(JetFem[0]);
      h_jet2emfrac[3]->Fill(JetFem[1]);

      //total number of photons/leptons in event
      for (unsigned int ielec = 0; ielec < Nelec; ++ielec) h_elecEt[3]->Fill(ElecEt[ielec]);
      h_Nelec[3]->Fill(Nelec);

      for (unsigned int imuon = 0; imuon < Nmuon; ++imuon) h_muonEt[3]->Fill(MuonEt[imuon]);
      h_Nmuon[3]->Fill(Nmuon);

      ////number of good photons/leptons in event
      //eleccount = 0;
      //for (unsigned int ielec = 0; ielec < Nelec; ++ielec) if(ElecEt[ielec]>cut_elecet) {
      //	h_elecEt[3]->Fill(ElecEt[ielec]);
      //	eleccount++;}
      //h_Nelec[3]->Fill(eleccount);
      //
      //muoncount = 0;
      //for (unsigned int imuon = 0; imuon < Nmuon; ++imuon) if(MuonEt[imuon]>cut_muonet) {
      //	h_muonEt[3]->Fill(MuonEt[imuon]);
      //	muoncount++;}
      //h_Nmuon[3]->Fill(muoncount);
      
      h_jet12dphi[3]->Fill(sqrt((Jetphi[0]-Jetphi[1])*(Jetphi[0]-Jetphi[1])));      
      h_jet1metdphi[3]->Fill(jet1metdphi);
      h_jet2metdphi[3]->Fill(jet2metdphi);

      h_jet1phi[1]->Fill(Jetphi[0]);
      h_jet2phi[1]->Fill(Jetphi[1]);
      h_METphi[1]->Fill(METphi_fullcorr_nocc);
    }
  }

  // scale histograms to desired values
  char ytitle[128];
  sprintf(ytitle,"Events / %2.0f pb^{-1}",luminosity_);
  h_selections[0]->Scale(scale);
  h_selections[0]->GetYaxis()->SetTitle(ytitle);
  h_selections[1]->Scale(scale);
  h_selections[1]->GetYaxis()->SetTitle(ytitle);
  h_selections[0]->GetXaxis()->SetBinLabel(1,"nJets");
  h_selections[0]->GetXaxis()->SetBinLabel(2,"E_{T}^{Jet_{1}}");
  h_selections[0]->GetXaxis()->SetBinLabel(3,"E_{T}^{Jet_{2}}");
  h_selections[0]->GetXaxis()->SetBinLabel(4,"#eta_{Jet_{1}}");
  h_selections[0]->GetXaxis()->SetBinLabel(5,"#eta_{Jet_{2}}");
  h_selections[0]->GetXaxis()->SetBinLabel(6,"EM_{FRAC}^{Jet_{1}}");
  h_selections[0]->GetXaxis()->SetBinLabel(7,"EM_{FRAC}^{Jet_{2}}");
  h_selections[0]->GetXaxis()->SetBinLabel(8,"#Delta#phi_{Jet_{1}, #slashE_{T}}");
  h_selections[0]->GetXaxis()->SetBinLabel(9,"#Delta#phi_{Jet_{2}, #slashE_{T}}");
  h_selections[0]->GetXaxis()->SetBinLabel(10,"#slashE_{T}");
  h_selections[0]->GetXaxis()->SetBinLabel(11,"H_{T}");
  h_selections[0]->GetXaxis()->SetBinLabel(12,"#slashH_{T}");
  h_selections[0]->GetXaxis()->SetBinLabel(13,"M_{eff}");
  h_selections[0]->GetXaxis()->SetBinLabel(14,"#Delta#phi_{Jet_{1}, Jet_{2}}");
  h_selections[0]->GetXaxis()->SetBinLabel(15,"E^{MAX}_{T}^{jets}");
  h_selections[0]->GetXaxis()->SetBinLabel(16,"E^{e^{-}}_{T} < 15 GeV");
  h_selections[0]->GetXaxis()->SetBinLabel(17,"E^{#mu}_{T} < 15 GeV");
  h_selections[0]->GetXaxis()->SetBinLabel(18,"ALL");
  h_selections[0]->SetStats(kFALSE);

  h_selections[1]->GetXaxis()->SetBinLabel(1,"nJets");
  h_selections[1]->GetXaxis()->SetBinLabel(2,"E_{T}^{Jet_{1}}");
  h_selections[1]->GetXaxis()->SetBinLabel(3,"E_{T}^{Jet_{2}}");
  h_selections[1]->GetXaxis()->SetBinLabel(4,"#eta_{Jet_{1}}");
  h_selections[1]->GetXaxis()->SetBinLabel(5,"#eta_{Jet_{2}}");
  h_selections[1]->GetXaxis()->SetBinLabel(6,"EM_{FRAC}^{Jet_{1}}");
  h_selections[1]->GetXaxis()->SetBinLabel(7,"EM_{FRAC}^{Jet_{2}}");
  h_selections[1]->GetXaxis()->SetBinLabel(8,"#Delta#phi_{Jet_{1}, #slashE_{T}}");
  h_selections[1]->GetXaxis()->SetBinLabel(9,"#Delta#phi_{Jet_{2}, #slashE_{T}}");
  h_selections[1]->GetXaxis()->SetBinLabel(10,"#slashE_{T}");
  h_selections[1]->GetXaxis()->SetBinLabel(11,"H_{T}");
  h_selections[1]->GetXaxis()->SetBinLabel(12,"#slashH_{T}");
  h_selections[1]->GetXaxis()->SetBinLabel(13,"M_{eff}");
  h_selections[1]->GetXaxis()->SetBinLabel(14,"#Delta#phi_{Jet_{1}, Jet_{2}}");
  h_selections[1]->GetXaxis()->SetBinLabel(15,"E^{MAX}_{T}^{jets}");
  h_selections[1]->GetXaxis()->SetBinLabel(16,"E^{e^{-}}_{T} < 15 GeV");
  h_selections[1]->GetXaxis()->SetBinLabel(17,"E^{#mu}_{T} < 15 GeV");
  h_selections[1]->SetStats(kFALSE);

  for (int mine = 0; mine < 2; mine++) {
    sprintf(ytitle,"Events / %2.0f pb^{-1}",luminosity_);
    h_jet1phi[mine]->Scale(scale);
    h_jet1phi[mine]->GetYaxis()->SetTitle(ytitle);
    h_jet2phi[mine]->Scale(scale);
    h_jet2phi[mine]->GetYaxis()->SetTitle(ytitle);
    h_METphi[mine]->Scale(scale);
    h_METphi[mine]->GetYaxis()->SetTitle(ytitle);
    sprintf(ytitle,"Events /25 GeV / %2.0f pb^{-1}",luminosity_);
    h_MT[mine]->Scale(scale);
    h_MT[mine]->GetYaxis()->SetTitle(ytitle);}
  
  for (int mine = 0; mine < 4; mine++) {
    sprintf(ytitle,"Events / %2.0f pb^{-1}",luminosity_);
    for (int my = 0; my < 2; my++) {
      h_Njets[mine][my]->Scale(scale);
      h_Njets[mine][my]->GetYaxis()->SetTitle(ytitle);}

    h_Nelec[mine]->Scale(scale);
    h_Nelec[mine]->GetYaxis()->SetTitle(ytitle);
    h_Nmuon[mine]->Scale(scale);
    h_Nmuon[mine]->GetYaxis()->SetTitle(ytitle);

    h_jet1eta[mine]->Scale(scale);
    h_jet1eta[mine]->GetYaxis()->SetTitle(ytitle);
    h_jet2eta[mine]->Scale(scale);
    h_jet2eta[mine]->GetYaxis()->SetTitle(ytitle);
    h_jet1emfrac[mine]->Scale(scale);
    h_jet1emfrac[mine]->GetYaxis()->SetTitle(ytitle);
    h_jet2emfrac[mine]->Scale(scale);
    h_jet2emfrac[mine]->GetYaxis()->SetTitle(ytitle);
      
    h_jetFem[mine]->Scale(scale);
    h_jetFem[mine]->GetYaxis()->SetTitle(ytitle);
      
    h_jet12dphi[mine]->Scale(scale);
    h_jet12dphi[mine]->GetYaxis()->SetTitle(ytitle);
    h_jet1metdphi[mine]->Scale(scale);
    h_jet1metdphi[mine]->GetYaxis()->SetTitle(ytitle);
    h_jet2metdphi[mine]->Scale(scale);
    h_jet2metdphi[mine]->GetYaxis()->SetTitle(ytitle);

    sprintf(ytitle,"Events / 50 GeV / %2.0f pb^{-1}",luminosity_);

    h_jet1et[mine]->Scale(scale);
    h_jet1et[mine]->GetYaxis()->SetTitle(ytitle);
    h_jet2et[mine]->Scale(scale);
    h_jet2et[mine]->GetYaxis()->SetTitle(ytitle);      
    h_MET[mine]->Scale(scale);
    h_MET[mine]->GetYaxis()->SetTitle(ytitle);
    h_MHT[mine]->Scale(scale);
    h_MHT[mine]->GetYaxis()->SetTitle(ytitle);

    if (mine==1||mine==3) sprintf(ytitle,"Events / 5 GeV / %2.0f pb^{-1}",luminosity_);
    else                  sprintf(ytitle,"Events / 25 GeV / %2.0f pb^{-1}",luminosity_);

    h_jetallet[mine]->Scale(scale);
    h_jetallet[mine]->GetYaxis()->SetTitle(ytitle);      
    h_elecEt[mine]->Scale(scale);
    h_elecEt[mine]->GetYaxis()->SetTitle(ytitle);
    h_muonEt[mine]->Scale(scale);
    h_muonEt[mine]->GetYaxis()->SetTitle(ytitle);

    sprintf(ytitle,"Events / 100 GeV / %2.0f pb^{-1}",luminosity_);

    h_HT[mine]->Scale(scale);
    h_HT[mine]->GetYaxis()->SetTitle(ytitle);
    h_Meff[mine]->Scale(scale);
    h_Meff[mine]->GetYaxis()->SetTitle(ytitle);

  }

  //
  int Nevents = nentries;
  printf("\"Nevents\"   \"preselection\"   \"finaljet\"   \"dphiselection\"   \"htselection\"   \"metselection\"   \"mhtselection\"   \"leptonveto\"\n");
  printf("\"%7d\"   \"%12d\"   \"%8d\"   \"%13d\"   \"%11d\"   \"%12d\"   \"%11d\"   \"%10d\"\n",Nevents,pscounter,fjcounter,dphicounter,htcounter,metcounter,mhtcounter,leptoncounter);

  //---------------------------------------------------------
  //Write out root file with histograms
  //---------------------------------------------------------
  file->cd();
  file->Write();
   
}
