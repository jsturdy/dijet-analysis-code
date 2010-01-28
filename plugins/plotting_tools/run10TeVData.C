{
  const double NUMSAMPS = 21;

  gROOT->ProcessLine(".L allData.C++");

  double sigmaeff[NUMSAMPS][2] = {{110.0,1.0},{16.06,1.0},{1.94,1.0},                                //LM0  LM1  LM5
				  {30000.0,0.36},{3600.0,0.30},{6625.0,0.40},                        //WJets  ZJets  ZInvisibleJets
				  /*{292,0.24},*/{242.8,1.0}                                         //TTJets madgraph, TTJets pythia
				  {0.000008595,1.0},{0.001422,1.0},{0.1721,1.0},                     //PT3000  PT2200  PT1400
				  {11.94,1.0},/*{315.3,1.0},*/{3669,1.0},                            //PT800  PT470  PT300
				  {62510,1.0},{1936000,1.0},{109000000,1.0},{1458000000,1.0}         //PT170  PT80  PT30 PT15
				  {2.14,1.0},{0.402,10}, {0.106,1.0},                                //LQtoCMu M300, M400, M500
				  {11840,0.691},{1944,1.0}};                                         //Wmunu, Zmumu

  string infilenames[NUMSAMPS]  = {"/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/SUSY_LM0_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/SUSY_LM1_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/SUSY_LM5_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/WJets_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/ZJets_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/ZInvisibleJets_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/TTJets_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_3000_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_2200_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_1400_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_800_10TeV_PATtified.root",
				   //"/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_470_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_300_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_170_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_80_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_30_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_15_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/LQtoCMu_M300_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/LQtoCMu_M400_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/LQtoCMu_M500_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/Wmunu_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/Zmumu_10TeV_PATtified.root"};


  string outfilenames[NUMSAMPS]  = {"10TeV_SUSY_LM0.root","10TeV_SUSY_LM1.root","10TeV_SUSY_LM5.root",
				    "10TeV_WJets.root","10TeV_ZJets.root","10TeV_ZInvisibleJets.root","10TeV_TTJets.root",
				    "10TeV_QCD_pt_3000.root","10TeV_QCD_pt_2200.root","10TeV_QCD_pt_1400.root",
				    "10TeV_QCD_pt_800.root",/*"10TeV_QCD_pt_470.root",*/"10TeV_QCD_pt_300.root",
				    "10TeV_QCD_pt_170.root","10TeV_QCD_pt_80.root","10TeV_QCD_pt_30.root","10TeV_QCD_pt_15.root",
				    "10TeV_LQtoCMu_M300.root","10TeV_LQtoCMu_M400.root","10TeV_LQtoCMu_M500.root",
				    "10TeV_Wmunu.root","10TeV_Zmumu.root"};

  printf("\"sample\"   \"outfilename\"   \"cross section\"   \"efficiency\"   \"met cut\"   \"mht cut\"  \"metTag\"  \"jetTag\"\n");

  double  mymet = 200.;
  double  mymht = 0.;
  std::string mettag = "tcMET";   //CaloMET, pfMET, tcMET
  std::string jettag = "pfJets";  //CaloJets, pfJets
  char tmpfilename[128];

  for (int z = 15; z < NUMSAMPS; z++) {
    TChain *chainA = new TChain("jaredsusy/AllData");
    chainA->Add(infilenames[z].c_str());
    sprintf(tmpfilename,"%s_%s_%s",mettag.c_str(),jettag.c_str(),outfilenames[z].c_str());
    printf("\"%d\"   \"%s\"   \"%2.12f\"   \"%2.2f\"   \"%d\"   \"%d\"  \"%s\"  \"%s\"\n",z,tmpfilename,sigmaeff[z][0],sigmaeff[z][1],mymet,mymht,mettag.c_str(),jettag.c_str());
    //printf("\"%d\"   \"%s\"   \"%2.12f\"   \"%2.2f\"   \"%d\"   \"%d\"  \"%s\"  \"%s\"\n",z,outfilenames[z].c_str(),sigmaeff[z][0],sigmaeff[z][1],mymet,mymht,mettag.c_str(),jettag.c_str());

    TTree *treeA = chainA;
    
    //allData treeData(treeA, files,luminosity,cross section,efficiency,mymet,mymht,mettag,jettag);
    allData treeData(treeA, outfilenames[z],100,sigmaeff[z][0],sigmaeff[z][1],mymet,mymht,mettag,jettag);
    treeData.Loop();
  }
}

