void analysisPlots10TeV(std::string mettype="CaloMET", std::string jettype="CaloJets", std::string energy="10TeV")
{
  static const double NUMHISTS  = 23;//23;
  static const double NUMHISTS2 = 6;
  static const double NUMVARS   = 4;
  static const double NUMFILES  = 6;
  static const double localpi  = acos(-1);

  TFile* file[NUMFILES];
  TString filenames[NUMFILES];

  Color_t linecolor[NUMFILES];
  UInt_t  linestyle[NUMFILES];
  UInt_t  linewidth;//[NUMFILES];

  TCanvas *mycanvas[NUMHISTS][NUMVARS];
  TCanvas *mycanvas2[NUMHISTS2][2];
 
  Double_t max[NUMHISTS][NUMVARS] = {0.0};
  Double_t max2[NUMHISTS2][2]     = {0.0};

  TString histnames[NUMHISTS][NUMVARS];
  TString histnames2[NUMHISTS2][2];
  TLegend *leg;
  TString histtitle[NUMHISTS];
  TString histtitle2[NUMHISTS2+1];

  TH1F *hist[NUMFILES][NUMHISTS][NUMVARS];
  TH1F *hist2[NUMFILES][NUMHISTS2][2];

  //filenames[0]  = "SUSY_LM0";
  filenames[0]  = "SUSY_LM1";
  //filenames[2]  = "SUSY_LM5";
  //filenames[3]  = "SM_Background";
  //filenames[4]  = "LQtoCMu_M300";
  //filenames[5]  = "LQtoCMu_M400";
  //filenames[6]  = "LQtoCMu_M500";

  filenames[1]  = "WJets";
  filenames[2]  = "ZJets";
  filenames[3]  = "ZInvisibleJets";
  filenames[4]  = "TTJets";  
  filenames[5]  = "QCD_Background";

  //filenames[6]  = "VectorBoson";

  //filenames[1]  = "Zmumu";
  //filenames[1]  = "Wmunu";

  histnames[0][0]  = "h_pre_cuts_01_jet1et";       histnames[0][1]  = "h_individual_cuts_01_jet1et";       histnames[0][2]  = "h_N1_cuts_01_jet1et";            histnames[0][3]  = "h_post_cuts_01_jet1et";
  histnames[1][0]  = "h_pre_cuts_02_jet2et";       histnames[1][1]  = "h_individual_cuts_02_jet2et";       histnames[1][2]  = "h_N1_cuts_02_jet2et";            histnames[1][3]  = "h_post_cuts_02_jet2et";
  histnames[2][0]  = "h_pre_cuts_03_jetallet";     histnames[2][1]  = "h_individual_cuts_03_jetallet";     histnames[2][2]  = "h_N1_cuts_03_jetallet";          histnames[2][3]  = "h_post_cuts_03_jetallet";
  histnames[3][0]  = "h_pre_cuts_11_MET";          histnames[3][1]  = "h_individual_cuts_11_MET";	   histnames[3][2]  = "h_N1_cuts_11_MET";	        histnames[3][3]  = "h_post_cuts_11_MET";
  histnames[4][0]  = "h_pre_cuts_12_HT";           histnames[4][1]  = "h_individual_cuts_12_HT";	   histnames[4][2]  = "h_N1_cuts_12_HT";	        histnames[4][3]  = "h_post_cuts_12_HT";
  histnames[5][0]  = "h_pre_cuts_13_MHT";          histnames[5][1]  = "h_individual_cuts_13_MHT";	   histnames[5][2]  = "h_N1_cuts_13_MHT";	        histnames[5][3]  = "h_post_cuts_13_MHT";
  histnames[6][0]  = "h_pre_cuts_14_Meff";         histnames[6][1]  = "h_individual_cuts_14_Meff";  	   histnames[6][2]  = "h_N1_cuts_14_Meff";	        histnames[6][3]  = "h_post_cuts_14_Meff";
  histnames[7][0]  = "h_pre_cuts_21_jet1metdphi";  histnames[7][1]  = "h_individual_cuts_21_jet1metdphi";  histnames[7][2]  = "h_N1_cuts_21_jet1metdphi";       histnames[7][3]  = "h_post_cuts_21_jet1metdphi";
  histnames[8][0]  = "h_pre_cuts_22_jet2metdphi";  histnames[8][1]  = "h_individual_cuts_22_jet2metdphi";  histnames[8][2]  = "h_N1_cuts_22_jet2metdphi";       histnames[8][3]  = "h_post_cuts_22_jet2metdphi";
  histnames[9][0]  = "h_pre_cuts_23_jet12dphi";    histnames[9][1]  = "h_individual_cuts_23_jet12dphi";    histnames[9][2]  = "h_N1_cuts_23_jet12dphi";         histnames[9][3]  = "h_post_cuts_23_jet12dphi";
  histnames[10][0] = "h_pre_cuts_30_Njets";        histnames[10][1] = "h_individual_cuts_30_Njets";	   histnames[10][2] = "h_N1_cuts_30_Njets";	        histnames[10][3] = "h_post_cuts_30_Njets";
  histnames[11][0] = "h_pre_cuts_31_Ngoodjets";    histnames[11][1] = "h_individual_cuts_31_Ngoodjets";    histnames[11][2] = "h_N1_cuts_31_Ngoodjets";         histnames[11][3] = "h_post_cuts_31_Ngoodjets";
  histnames[12][0] = "h_pre_cuts_32_jet1eta";      histnames[12][1] = "h_individual_cuts_32_jet1eta";      histnames[12][2] = "h_N1_cuts_32_jet1eta";           histnames[12][3] = "h_post_cuts_32_jet1eta";
  histnames[13][0] = "h_pre_cuts_33_jet2eta";      histnames[13][1] = "h_individual_cuts_33_jet2eta";      histnames[13][2] = "h_N1_cuts_33_jet2eta";           histnames[13][3] = "h_post_cuts_33_jet2eta";
  histnames[14][0] = "h_pre_cuts_41_jetFem";       histnames[14][1] = "h_individual_cuts_41_jetFem";       histnames[14][2] = "h_N1_cuts_41_jetFem";            histnames[14][3] = "h_post_cuts_41_jetFem";
  histnames[15][0] = "h_pre_cuts_42_jet1emfrac";   histnames[15][1] = "h_individual_cuts_42_jet1emfrac";   histnames[15][2] = "h_N1_cuts_42_jet1emfrac";        histnames[15][3] = "h_post_cuts_42_jet1emfrac";
  histnames[16][0] = "h_pre_cuts_43_jet2emfrac";   histnames[16][1] = "h_individual_cuts_43_jet2emfrac";   histnames[16][2] = "h_N1_cuts_43_jet2emfrac";        histnames[16][3] = "h_post_cuts_43_jet2emfrac";
  histnames[17][0] = "h_pre_cuts_51a_Nelecs";      histnames[17][1] = "h_individual_cuts_51a_Nelecs"; 	   histnames[17][2] = "h_N1_cuts_51a_Nelecs";           histnames[17][3] = "h_post_cuts_51a_Nelecs";
  histnames[18][0] = "h_pre_cuts_52a_Nmuons";      histnames[18][1] = "h_individual_cuts_52a_Nmuons"; 	   histnames[18][2] = "h_N1_cuts_52a_Nmuons";           histnames[18][3] = "h_post_cuts_52a_Nmuons";
  histnames[19][0] = "h_pre_cuts_51b_Ngoodelecs";  histnames[19][1] = "h_individual_cuts_51b_Ngoodelecs";  histnames[19][2] = "h_N1_cuts_51b_Ngoodelecs";        histnames[19][3] = "h_post_cuts_51b_Ngoodelecs";
  histnames[20][0] = "h_pre_cuts_52b_Ngoodmuons";  histnames[20][1] = "h_individual_cuts_52b_Ngoodmuons";  histnames[20][2] = "h_N1_cuts_52b_Ngoodmuons";        histnames[20][3] = "h_post_cuts_52b_Ngoodmuons";
  //histnames[19][0] = "h_pre_cuts_53_eleceta";      histnames[19][1] = "h_individual_cuts_53_eleceta";	   histnames[19][2] = "h_N1_cuts_53_eleceta";	   histnames[19][3] = "h_post_cuts_53_eleceta";
  //histnames[20][0] = "h_pre_cuts_54_muoneta";      histnames[20][1] = "h_individual_cuts_54_muoneta";	   histnames[20][2] = "h_N1_cuts_54_muoneta";	   histnames[20][3] = "h_post_cuts_54_muoneta";
  histnames[21][0] = "h_pre_cuts_55_elecet";       histnames[21][1] = "h_individual_cuts_55_elecet";	   histnames[21][2] = "h_N1_cuts_55_elecet";	        histnames[21][3] = "h_post_cuts_55_elecet";
  histnames[22][0] = "h_pre_cuts_56_muonet";       histnames[22][1] = "h_individual_cuts_56_muonet";	   histnames[22][2] = "h_N1_cuts_56_muonet";	        histnames[22][3] = "h_post_cuts_56_muonet";

  histnames2[0][0] = "h_pre_cuts_9_METphi";   histnames2[0][1] = "h_post_cuts_9_METphi";
  histnames2[1][0] = "h_pre_cuts_9_jet1phi";  histnames2[1][1] = "h_post_cuts_9_jet1phi";
  histnames2[2][0] = "h_pre_cuts_9_jet2phi";  histnames2[2][1] = "h_post_cuts_9_jet2phi";
  histnames2[3][0] = "h_pre_cuts_9_MT";       histnames2[3][1] = "h_post_cuts_9_MT";
  histnames2[4][0] = "h_pre_cuts_9_Minv";     histnames2[4][1] = "h_post_cuts_9_Minv";
  histnames2[5][0] = "h_selections";          histnames2[5][1] = "h_N1_selections";

  histtitle[0]  = "Jet 1 E_{T}";
  histtitle[1]  = "Jet 2 E_{T}";
  histtitle[2]  = "E^{Jet}_{T}";
  histtitle[3]  = "#slashE_{T}";
  histtitle[4]  = "H_{T}";
  histtitle[5]  = "#slashH_{T}";
  histtitle[6]  = "M_{eff}";
  histtitle[7]  = "#Delta#phi_{J_{1},#slashE_{T}}";
  histtitle[8]  = "#Delta#phi_{J_{2},#slashE_{T}}";
  histtitle[9]  = "#Delta#phi_{J_{1}J_{2}}";
  histtitle[10] = "N_{all jets}";
  histtitle[11] = "N_{good jets}";
  histtitle[12] = "#eta_{J_{1}}";
  histtitle[13] = "#eta_{J_{2}}";
  histtitle[14] = "Jet f_{EM}";
  histtitle[15] = "EM_{frac}^{J_{1}}";
  histtitle[16] = "EM_{frac}^{J_{2}}";
  histtitle[17] = "N_{e^{-}}";
  histtitle[18] = "N_{#mu}";
  histtitle[19] = "N^{good}_{e^{-}}";
  histtitle[20] = "N^{good}_{#mu}";
  //histtitle[19] = "#eta_{e^{-1}}";
  //histtitle[20] = "#eta_{#mu}";
  histtitle[21] = "E^{e^{-}}_{T}";
  histtitle[22] = "E^{#mu}_{T}";
  //histtitle[25] = "E^{#mu}_{T}";
  //histtitle[26] = "E^{#mu}_{T}";

  histtitle2[0] = "#slashE_{T}#phi";
  histtitle2[1] = "#phi_{J_{1}}";
  histtitle2[2] = "#phi_{J_{2}}";
  histtitle2[3] = "M_{T}";
  histtitle2[4] = "M_{INV}";
  histtitle2[5] = "selections";
  histtitle2[6] = "N-1 selections";

 
  linecolor[0] = kBlue;//LM0
  linecolor[1] = kRed+1;//LM1
  linecolor[2] = kGreen+1;//LM5
  linecolor[3] = kOrange+7;//SM
  linecolor[4] = kBlue+3;//LQ to CMu M300
  linecolor[5] = kYellow+3;//LQ to CMu M400
  //linecolor[6] = kBlue+8;//LQ to CMu M400

  linestyle[0] = 2;
  linestyle[1] = 3;
  linestyle[2] = 4;
  linestyle[3] = 8;
  linestyle[4] = 9;
  linestyle[5] = 10;
  //linestyle[6] = 1;

  linewidth = 5;

  cout<<"here is 1"<<endl;
  for (int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++) {
      mycanvas[z][tt]  = new TCanvas(histnames[z][tt], "", 600,600);}}
  
  cout<<"here is 2"<<endl;
  for (int ss = 0; ss < 2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++) {
      mycanvas2[y][ss] = new TCanvas(histnames2[y][ss], "", 600,600);}}

  TString metTag    = mettype;
  TString jetTag    = jettype;
  TString energyTag = energy;

  cout<<"here is 3"<<endl;
  for (int j = 0; j < NUMFILES; j++) {
    TString filepath = "./"+metTag+"_"+jetTag+"/"+energyTag+"/"+metTag+"_"+jetTag+"_"+energyTag+"_"+filenames[j]+".root";
    cout<<filepath<<endl;
    file[j] = new TFile(filepath);
    
    cout<<"here is 3.1"<<endl;
    for (int tt = 0; tt < NUMVARS; tt++) {
      cout<<"here is 3.1.1"<<endl;
      for (int z = 0; z < NUMHISTS; z++) {
	cout<<"tt and z "<<tt<<","<<z<<endl;
	hist[j][z][tt] = (TH1F*)gDirectory->Get(histnames[z][tt]);
	double histmax = hist[j][z][tt]->GetMaximum();
	max[z][tt] = (max[z][tt]>histmax)?max[z][tt]:histmax;}}
    
    cout<<"here is 3.2"<<endl;
    for (int ss = 0; ss < 2; ss++) {
      for (int y = 0; y < NUMHISTS2; y++) {
	hist2[j][y][ss] = (TH1F*)gDirectory->Get(histnames2[y][ss]);
	Double_t hist2max = hist2[j][y][ss]->GetMaximum();
	max2[y][ss] = (max2[y][ss]>hist2max)?max2[y][ss]:hist2max;}}}
  
  cout<<"here is 4"<<endl;
  for (int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++ ) max[z][tt] *= 1.1;}
  
  cout<<"here is 5"<<endl;
  for (int ss = 0; ss < 2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++ ) max2[y][ss] *= 1.1;}

  cout<<"here is 6"<<endl;
  for (int j = 0; j < NUMFILES; j++) {
    cout<<"here is 6.1"<<endl;
    for (int tt = 0; tt < NUMVARS; tt++) {
      for (int z = 0; z < NUMHISTS; z++) {
	//hist[j][z][tt]->SetFillStyle(3001);
	hist[j][z][tt]->SetLineColor(linecolor[j]);
	hist[j][z][tt]->SetLineStyle(linestyle[j]);
	hist[j][z][tt]->SetLineWidth(linewidth);
	//hist[j][z][tt]->SetLineStyle();
	//hist[j][z][tt]->SetFillColor(linecolor[j]);
	hist[j][z][tt]->SetStats(kFALSE);}}
    
    cout<<"here is 6.2"<<endl;
    for (int ss = 0; ss < 2; ss++) {
      for (int y = 0; y < NUMHISTS2; y++) {
	//hist2[j][y][ss]->SetFillStyle(3001);
	hist2[j][y][ss]->SetLineColor(linecolor[j]);
	hist2[j][y][ss]->SetLineStyle(linestyle[j]);
	hist2[j][y][ss]->SetLineWidth(linewidth);
	//hist2[j][y][ss]->SetLineStyle();
	//hist2[j][y][ss]->SetFillColor(linecolor[j]);
	hist2[j][y][ss]->SetStats(kFALSE);}}}

  leg = new TLegend(.70, .80, 1., 1.);
  cout<<"here is 7"<<endl;
  for (int j = 0; j < NUMFILES; j++) {
    leg->AddEntry(hist[j][0][0], metTag+"_"+jetTag+"_"+energyTag+"_"+filenames[j],"l");
    //leg->GetEntry()->SetTextSize(0);
    //leg->GetEntry()->SetEntrySeparation(0.075);
  }

  cout<<"here is 8"<<endl;
  for(int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++) {
      mycanvas[z][tt]->cd();
      hist[0][z][tt]->SetMaximum(max[z][tt]);
      hist[0][z][tt]->SetMinimum(0.0001);
      hist[0][z][tt]->GetXaxis()->SetTitle(histtitle[z]);
      hist[0][z][tt]->GetXaxis()->SetLabelSize(0.03);
      hist[0][z][tt]->GetYaxis()->SetTitleOffset(1.2);
      //hist[0][z][tt]->GetYaxis()->SetLabelOffset(-0.05);
      //if (z==3) {
      //  char ytitle[128];
      //  sprintf(ytitle,"Events / 25 GeV / %2.0f pb^{-1}",100);
      //  hist[0][z][tt]->GetYaxis()->SetTitle(ytitle);}
      //if (z==19||z==20) {
      //  if (tt==1||tt==3) {
      //    char ytitle[128];
      //    sprintf(ytitle,"Events / 1 GeV / %2.0f pb^{-1}",100);
      //    hist[0][z][tt]->GetYaxis()->SetTitle(ytitle);}}

      hist[0][z][tt]->GetYaxis()->SetLabelSize(0.03);
      gPad->SetLogy(kTRUE);
      hist[0][z][tt]->Draw();
      for (int j = 1; j < NUMFILES; j++) hist[j][z][tt]->Draw("same");
      leg->Draw();
      char outputimage[128];
      std::string outhistname = histnames[z][tt];
      sprintf(outputimage,"./%s_%s/%s/%s.ps",mettype.c_str(),jettype.c_str(),energy.c_str(),outhistname.c_str());
      mycanvas[z][tt]->SaveAs(outputimage);}}
      //mycanvas[z][tt]->SaveAs();}}

  cout<<"here is 9"<<endl;
  for(int ss = 0; ss <2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++) {
      mycanvas2[y][ss]->cd();
      hist2[0][y][ss]->SetMaximum(max2[y][ss]);
      hist2[0][y][ss]->SetMinimum(0.0001);
      hist2[0][y][ss]->GetXaxis()->SetTitle(histtitle2[y]);
      hist2[0][y][ss]->GetXaxis()->SetLabelSize(0.03);
      hist2[0][y][ss]->GetYaxis()->SetTitleOffset(1.2);
      //hist2[0][y][ss]->GetYaxis()->SetLabelOffset(-0.05);
      hist2[0][y][ss]->GetYaxis()->SetLabelSize(0.03);
      gPad->SetLogy(kTRUE);
      hist2[0][y][ss]->Draw();
      for (int j = 1; j < NUMFILES; j++) hist2[j][y][ss]->Draw("same");
      leg->Draw();
      char outputimage2[128];
      std::string outhistname2 = histnames2[y][ss];
      sprintf(outputimage2,"./%s_%s/%s/%s.ps",mettype.c_str(),jettype.c_str(),energy.c_str(),outhistname2.c_str());
      mycanvas2[y][ss]->SaveAs(outputimage2);}}  
      //mycanvas2[y][ss]->SaveAs();}}  
}
