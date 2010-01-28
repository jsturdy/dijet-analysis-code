{
  const double NUMSAMPS = 15;

  gROOT->ProcessLine(".L allData.C++");

  double sigmaeff[NUMSAMPS][2] = {{110.0,1.0},{16.06,1.0},//{1.94,1.0},                                     //LM0  LM1  LM5
				  {30000.0,0.36},{3600.0,0.30},//{2000.0,0.26},                            //WJets  ZJets  ZInvisibleJets
				  /*{292,0.24},*/{242.8,1}                                                   //TTJets madgraph, TTJets pythia
				  {0.000008595,1.0},{0.001422,1.0},{0.1721,1.0},                         //PT3000  PT2200  PT1400
				  {11.94,1.0},/*{315.3,1.0},*/{3669,1.0},                                    //PT800  PT470  PT300
				  {62510,1.0},{193600,1.0},{10090000,1.0},{104580000,1.0}};//,               //PT170  PT80  PT30 PT15
				  {5.880,1.0},{0.422,10}, {}};                                             //LQtoCMu M300, M400, M500


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
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/QCD_pt_15_10TeV_PATtified.root"
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/LQtoCMuM250_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/LQtoCMuM400_10TeV_PATtified.root",
				   "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/LQtoCMuM400_10TeV_PATtified.root"};

  string outfilenames[NUMSAMPS]  = {"10TeV_SUSY_LM0.root","10TeV_SUSY_LM1.root","10TeV_SUSY_LM5.root",
				    "10TeV_WJets.root","10TeV_ZJets.root","10TeV_ZInvisibleJets.root","10TeV_TTJets.root",
				    "10TeV_QCD_pt_3000.root","10TeV_QCD_pt_2200.root","10TeV_QCD_pt_1400.root",
				    "10TeV_QCD_pt_800.root",/*"10TeV_QCD_pt_470.root",*/"10TeV_QCD_pt_300.root",
				    "10TeV_QCD_pt_170.root","10TeV_QCD_pt_80.root","10TeV_QCD_pt_30.root","10TeV_QCD_pt_15.root",
				    "10TeV_LQtoCMu_M300.root","10TeV_LQtoCMu_M400.root","10TeV_LQtoCMu_M500.root"};

  printf("\"sample\"   \"outfilename\"   \"cross section\"   \"efficiency\"   \"met cut\"   \"mht cut\"\n");

  double  mymet = 350.;
  double  mymht = 200.;

  for (int z = 0; z < NUMSAMPS; z++) {
    TChain *chainA = new TChain("jaredsusy/AllData");
    chainA->Add(infilenames[z].c_str());
    printf("\"%d\"   \"%s\"   \"%2.12f\"   \"%2.2f\"   \"%d\"   \"%d\"\n",z,outfilenames[z].c_str(),sigmaeff[z][0],sigmaeff[z][1],mymet,mymht);
    TTree *treeA = chainA;
    
    allData treeData(treeA, outfilenames[z],100,sigmaeff[z][0],sigmaeff[z][1],mymet,mymht);
    treeData.Loop();
  }
}

