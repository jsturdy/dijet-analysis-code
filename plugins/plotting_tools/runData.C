{
  const double NUMSAMPS = 15;

  gROOT->ProcessLine(".L allData.C++");

  double sigmaeff[NUMSAMPS][2] = {{110.0,1.0},{16.06,1.0},//{1.94,1.0},                                     //LM0  LM1  LM5
				  {30000.0,0.36},{3600.0,0.30},//{2000.0,0.26},                            //WJets  ZJets  ZInvisibleJets
				  {292,0.24},//{242.8,1}                                                   //TTJets madgraph, TTJets pythia
				  {0.000008595,1.0},{0.001422,1.0},{0.1721,1.0},                         //PT3000  PT2200  PT1400
				  {11.94,1.0},{315.3,1.0},{3669,1.0},                                    //PT800  PT470  PT300
				  {62510,1.0},{193600,1.0},{10090000,1.0},{104580000,1.0}};//,               //PT170  PT80  PT30 PT15
				  //{5.880,1.0},{0.422,10}};                                             //LQtoCMu M250 and M400 full sim V11
                                  //{5.822,1.0},{0.432,1.0},                                               //LQtoCMu M250 and LQtoCMu M400 fast sim V12



  string infilenames[NUMSAMPS]  = {"/uscms_data/d2/sturdy07/SUSY/LM0/SUSY_LM0_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/SUSY/LM1/SUSY_LM1_AK5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/LM5/SUSY_LM5_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/WJets/WJets_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/ZJets/ZJets_AK5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/VectorBosons/ZInvisibleJets/ZInvisibleJets_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/TTJets/TTJets_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt3000/QCD_pt_3000_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt2200/QCD_pt_2200_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt1400/QCD_pt_1400_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt800/QCD_pt_800_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt470/QCD_pt_470_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt300/QCD_pt_300_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt170/QCD_pt_170_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt80/QCD_pt_80_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt30/QCD_pt_30_AK5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt15/QCD_pt_15_AK5MET_PATtified.root"};//,
				   //"/uscms_data/d2/sturdy07/Exotica/LQtoCMu/M250/LQtoCMuM250_FastSim_PATtified.root",
				   //"/uscms_data/d2/sturdy07/Exotica/LQtoCMu/M400/LQtoCMuM400_FastSim_PATtified.root",

  string outfilenames[NUMSAMPS]  = {"AK5_MET_350_SUSY_LM0.root","AK5_MET_350_SUSY_LM1.root",//"AK5_MET_350_SUSY_LM5.root",
				    "AK5_MET_350_WJets.root","AK5_MET_350_ZJets.root",/*"AK5_MET_350_ZInvisibleJets.root",*/"AK5_MET_350_TTJets.root",
				    "AK5_MET_350_QCD_pt_3000.root","AK5_MET_350_QCD_pt_2200.root","AK5_MET_350_QCD_pt_1400.root",
				    "AK5_MET_350_QCD_pt_800.root","AK5_MET_350_QCD_pt_470.root","AK5_MET_350_QCD_pt_300.root",
				    "AK5_MET_350_QCD_pt_170.root","AK5_MET_350_QCD_pt_80.root","AK5_MET_350_QCD_pt_30.root","AK5_MET_350_QCD_pt_15.root"};//,
				    //"AK5_MET_350_LQtoCMu_M250_FastSim.root","AK5_MET_350_LQtoCMu_M400_FastSim.root",

  printf("\"sample\"   \"outfilename\"   \"cross section\"   \"efficiency\"   \"met cut\"   \"mht cut\"\n");

  double  mymet = 350.;
  double  mymht = 0.;

  for (int z = 0; z < NUMSAMPS; z++) {
    TChain *chainA = new TChain("jaredsusy/AllData");
    chainA->Add(infilenames[z].c_str());
    printf("\"%d\"   \"%s\"   \"%2.12f\"   \"%2.2f\"   \"%d\"   \"%d\"\n",z,outfilenames[z].c_str(),sigmaeff[z][0],sigmaeff[z][1],mymet,mymht);
    TTree *treeA = chainA;
    
    allData treeData(treeA, outfilenames[z],100,sigmaeff[z][0],sigmaeff[z][1],mymet,mymht);
    treeData.Loop();
  }
}

