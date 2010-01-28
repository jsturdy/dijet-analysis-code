{
  const double NUMSAMPS = 39-12;

  gROOT->ProcessLine(".L allData.C++");

  double sigmaeff[NUMSAMPS][2] = {{110.0,1.0},{16.06,1.0},{2.4,1.0},                                     //LM0  LM1  LM2
				  {11.79,1.0},{6.70,1.0},{1.94,1.0},                                     //LM3  LM4  LM5
				  //{1.28,1.0},{2.90,1.0},{2.86,1.0},                                      //LM6  LM7  LM8
				  //{11.58,1.0},{6.55,1.0},{3.24,1.0},                                     //LM9  LM10 LM11
				  //{2.97,1.0},{0.843,1.0},{0.299,1.0},                                    //GM1b GM1c GM1d
				  //{0.0124,1.0},{0.00582,1.0},{0.00311,1.0},                              //GM1e GM1f GM1g
				  {40000.0,0.45},{3700.0,0.40},{2000.0,0.26},                            //WJets  ZJets  ZInvisibleJets
				  {317.0,0.33},                                                          //TTJets
				  {0.0000086008,1.0},{0.0014207778,1.0},{0.1720187189,1.0},              //PT3000  PT2200  PT1400
				  {11.9419745,1.0},{315.5131272,1.0},{3664.608301,1.0},                  //PT800  PT470  PT300
				  {62562.87713,1.0},{1934639.567,1.0},{109057220.4,1.0},{1457159248,1.0},//PT170  PT80  PT30 PT15
				  //{5.880,1.0},{0.422,10}};                                              //LQtoCMu M250 and M400 full sim V11
                                  {5.822,1.0},{0.432,1.0},                                               //LQtoCMu M250 and LQtoCMu M400 fast sim V12
				  {11850,0.738},{1233,0.701},                                            //Wenu  Zee  
				  {44.8,1.0},{17.4,1.0},{7.1,1.0}};                                      //WWinc WZinc ZZinc


  string infilenames[NUMSAMPS]  = {"/uscms_data/d2/sturdy07/SUSY/LM0/SUSY_LM0_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/SUSY/LM1/SUSY_LM1_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/SUSY/LM2/SUSY_LM2_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/SUSY/LM3/SUSY_LM3_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/SUSY/LM4/SUSY_LM4_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/SUSY/LM5/SUSY_LM5_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/LM6/SUSY_LM6_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/LM7/SUSY_LM7_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/LM8/SUSY_LM8_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/LM9/SUSY_LM9_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/LM10/SUSY_LM10_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/LM11/SUSY_LM11_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/GM1b/SUSY_GM1b_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/GM1c/SUSY_GM1c_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/GM1d/SUSY_GM1d_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/GM1e/SUSY_GM1e_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/GM1f/SUSY_GM1f_SC5MET_PATtified.root",
				   //"/uscms_data/d2/sturdy07/SUSY/GM1g/SUSY_GM1g_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/WJets/WJets_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/ZJets/ZJets_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/ZInvisibleJets/ZInvisibleJets_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/TTJets/TTJets_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt3000/QCD_pt_3000_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt2200/QCD_pt_2200_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt1400/QCD_pt_1400_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt800/QCD_pt_800_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt470/QCD_pt_470_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt300/QCD_pt_300_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt170/QCD_pt_170_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt80/QCD_pt_80_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt30/QCD_pt_30_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/QCD/Pythia/pt15/QCD_pt_15_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/Exotica/LQtoCMu/M250/LQtoCMuM250_FastSim_PATtified.root",
				   "/uscms_data/d2/sturdy07/Exotica/LQtoCMu/M400/LQtoCMuM400_FastSim_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/Wenu/Wenu_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/Zee/Zee_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/WWinc/WWinc_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/WZinc/WZinc_SC5MET_PATtified.root",
				   "/uscms_data/d2/sturdy07/VectorBosons/ZZinc/ZZinc_SC5MET_PATtified.root"};

  string outfilenames[NUMSAMPS]  = {"SC5_MET_200/SUSY_LM0.root","SC5_MET_200/SUSY_LM1.root","SC5_MET_200/SUSY_LM2.root",
				    "SC5_MET_200/SUSY_LM3.root","SC5_MET_200/SUSY_LM4.root","SC5_MET_200/SUSY_LM5.root",
				    //"SC5_MET_200/SUSY_LM6.root","SC5_MET_200/SUSY_LM7.root","SC5_MET_200/SUSY_LM8.root",
				    //"SC5_MET_200/SUSY_LM9.root","SC5_MET_200/SUSY_LM10.root","SC5_MET_200/SUSY_LM11.root",
				    //"SC5_MET_200/SUSY_GM1b.root","SC5_MET_200/SUSY_GM1c.root","SC5_MET_200/SUSY_GM1d.root",
				    //"SC5_MET_200/SUSY_GM1e.root","SC5_MET_200/SUSY_GM1f.root","SC5_MET_200/SUSY_GM1g.root",
				    "SC5_MET_200/WJets.root","SC5_MET_200/ZJets.root","SC5_MET_200/ZInvisibleJets.root","SC5_MET_200/TTJets.root",
				    "SC5_MET_200/QCD_pt_3000.root","SC5_MET_200/QCD_pt_2200.root","SC5_MET_200/QCD_pt_1400.root",
				    "SC5_MET_200/QCD_pt_800.root","SC5_MET_200/QCD_pt_470.root","SC5_MET_200/QCD_pt_300.root",
				    "SC5_MET_200/QCD_pt_170.root","SC5_MET_200/QCD_pt_80.root","SC5_MET_200/QCD_pt_30.root","SC5_MET_200/QCD_pt_15.root",
				    "SC5_MET_200/LQtoCMu_M250_FastSim.root","SC5_MET_200/LQtoCMu_M400_FastSim.root",
				    "SC5_MET_200/Wenu.root","SC5_MET_200/Zee.root","SC5_MET_200/WWinc.root",
				    "SC5_MET_200/WZinc.root","SC5_MET_200/ZZinc.root"};

  printf("\"sample\"   \"outfilename\"   \"cross section\"   \"efficiency\"   \"met cut\"   \"mht cut\"\n");

  double  mymet = 200.;
  double  mymht = 200.;

  for (int z = 0; z < NUMSAMPS-7; z++) {
    TChain *chainA = new TChain("jaredsusy/AllData");
    chainA->Add(infilenames[z].c_str());
    printf("\"%d\"   \"%s\"   \"%2.12f\"   \"%2.2f\"   \"%d\"   \"%d\"\n",z,outfilenames[z].c_str(),sigmaeff[z][0],sigmaeff[z][1],mymet,mymht);
    TTree *treeA = chainA;
    
    allData treeData(treeA, outfilenames[z],100,sigmaeff[z][0],sigmaeff[z][1],mymet,mymht);

//    cout<<"done!"<<endl;//<<"setting cut variables...  ";
//
//    cut_njet     = 2;//min number of jets
//    cut_jet1et   = 100;//min et of jet 1
//    cut_jet2et   = 100;//min et of jet 2
//    cut_alljetet = 50;//max et of all other jets
//    cut_jet1eta  = 2.5;//max eta of jet 1
//    cut_jet2eta  = 2.5;//max eta of jet 2
//    cut_jet1phi  = 0.;
//    cut_jet2phi  = 0.;
//  
//    cut_jet1emfrac[0] = 0.05;//min em fraction of jet 1
//    cut_jet1emfrac[1] = 0.95;//max em fraction of jet 1
//    cut_jet2emfrac[0] = 0.05;//min em fraction of jet 2
//    cut_jet2emfrac[1] = 0.95;//max em fraction of jet 2
//  
//    cut_jetemfrac[0] = 0.;//min event em fraction
//    cut_jetemfrac[1] = 1.;//min event em fraction
//    cut_jet12dphi = 0.0;//dphi(jet1, jet2)
//  
//    cut_jet1metdphi = 0.3;//min dphi(jet1, met)
//    cut_jet2metdphi = 0.3;//min dphi(jet2, met)
//    //cut_jet3metdphi = 0.3;//min dphi(jet3, met)
//  
//    cut_ht   = 250.0;
//    cut_mht  = 200.0;
//    cut_met  = 350.0;
//    cut_meff = 0.0;
//    cout<<"done!";

    treeData.Loop();
  }
}

