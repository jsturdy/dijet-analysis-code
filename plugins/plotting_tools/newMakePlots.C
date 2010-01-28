void newMakePlots()
{
  static const double NUMHISTS  = 15;
  static const double NUMHISTS2 = 4;
  static const double NUMVARS   = 4;
  static const double NUMFILES  = 3;
  static const double localpi  = acos(-1);
  //bigcanvas1 = new TCanvas("bigcanvas1", "", 1440,900);
  //bigcanvas2 = new TCanvas("bigcanvas2", "", 1440,900);
  //bigcanvas3 = new TCanvas("bigcanvas3", "", 1440,900);
  TFile* file[NUMFILES];
  TString filenames[NUMFILES];

//  TFile *outfile = new TFile("mycanvasses.root","RECREATE");
//  outfile->cd();
  Color_t linecolor[NUMFILES];

  TCanvas *mycanvas[NUMHISTS][NUMVARS];
  TCanvas *mycanvas2[NUMHISTS2][2];
 
  double max[NUMHISTS][NUMVARS];
  double max2[NUMHISTS2][2];

  TString histnames[NUMHISTS][NUMVARS];
  TString histnames2[NUMHISTS2][2];
  TLegend *leg;
  //TLegend *leg[NUMHISTS];
  //TLegend *leg2[NUMHISTS2];
  TString histtitle[NUMHISTS];
  TString histtitle2[NUMHISTS2+1];

  TH1F *hist[NUMFILES][NUMHISTS][NUMVARS];
  TH1F *hist2[NUMFILES][NUMHISTS2][2];

  filenames[0] = "AK5_MET_350_SUSY_LM0";
  filenames[1] = "AK5_MET_350_SUSY_LM1";
//  filenames[2] = "SUSY_LM5";
  filenames[2] = "AK5_MET_SM_Background";
//  filenames[3] = "QCD_Background";
//  filenames[4] = "TTJets";
//  filenames[5] = "WJets";
//  filenames[6] = "ZJets";
//  filenames[7] = "ZInvisibleJets";

  histnames[0][0]  = "h_pre_cuts_Njets";          histnames[0][1]  = "h_individual_cuts_Njets_cut";	    histnames[0][2]  = "h_N1_cuts_Njets";	     histnames[0][3]  = "h_post_cuts_Njets";	     
  histnames[1][0]  = "h_pre_cuts_jet1et";         histnames[1][1]  = "h_individual_cuts_jet1et_cut";        histnames[1][2]  = "h_N1_cuts_jet1et";           histnames[1][3]  = "h_post_cuts_jet1et";     
  histnames[2][0]  = "h_pre_cuts_jet1eta";        histnames[2][1]  = "h_individual_cuts_jet1eta_cut";       histnames[2][2]  = "h_N1_cuts_jet1eta";          histnames[2][3]  = "h_post_cuts_jet1eta";     
  histnames[3][0]  = "h_pre_cuts_jet2et";         histnames[3][1]  = "h_individual_cuts_jet2et_cut";        histnames[3][2]  = "h_N1_cuts_jet2et";           histnames[3][3]  = "h_post_cuts_jet2et";     
  histnames[4][0]  = "h_pre_cuts_jet2eta";        histnames[4][1]  = "h_individual_cuts_jet2eta_cut";       histnames[4][2]  = "h_N1_cuts_jet2eta";          histnames[4][3]  = "h_post_cuts_jet2eta";     
  histnames[5][0]  = "h_pre_cuts_jetFem";         histnames[5][1]  = "h_individual_cuts_jetFem_cut";        histnames[5][2]  = "h_N1_cuts_jetFem";           histnames[5][3]  = "h_post_cuts_jetFem";     
  histnames[6][0]  = "h_pre_cuts_jet1emfrac";     histnames[6][1]  = "h_individual_cuts_jet1emfrac_cut";    histnames[6][2]  = "h_N1_cuts_jet1emfrac";       histnames[6][3]  = "h_post_cuts_jet1emfrac"; 
  histnames[7][0]  = "h_pre_cuts_jet2emfrac";     histnames[7][1]  = "h_individual_cuts_jet2emfrac_cut";    histnames[7][2]  = "h_N1_cuts_jet2emfrac";       histnames[7][3]  = "h_post_cuts_jet2emfrac"; 
  histnames[8][0]  = "h_pre_cuts_MET";            histnames[8][1]  = "h_individual_cuts_MET_cut";	    histnames[8][2]  = "h_N1_cuts_MET";	             histnames[8][3]  = "h_post_cuts_MET";	     
  histnames[9][0]  = "h_pre_cuts_HT";             histnames[9][1]  = "h_individual_cuts_HT_cut";	    histnames[9][2]  = "h_N1_cuts_HT";	             histnames[9][3]  = "h_post_cuts_HT";	     
  histnames[10][0] = "h_pre_cuts_MHT";            histnames[10][1] = "h_individual_cuts_MHT_cut";	    histnames[10][2] = "h_N1_cuts_MHT";	             histnames[10][3] = "h_post_cuts_MHT";	     
  histnames[11][0] = "h_pre_cuts_Meff";           histnames[11][1] = "h_individual_cuts_Meff_cut";	    histnames[11][2] = "h_N1_cuts_Meff";	     histnames[11][3] = "h_post_cuts_Meff";	     
  histnames[12][0] = "h_pre_cuts_jet12dphi";      histnames[12][1] = "h_individual_cuts_jet12dphi_cut";     histnames[12][2] = "h_N1_cuts_jet12dphi";        histnames[12][3] = "h_post_cuts_jet12dphi";  
  histnames[13][0] = "h_pre_cuts_jet1metdphi";    histnames[13][1] = "h_individual_cuts_jet1metdphi_cut";   histnames[13][2] = "h_N1_cuts_jet1metdphi";      histnames[13][3] = "h_post_cuts_jet1metdphi";
  histnames[14][0] = "h_pre_cuts_jet2metdphi";    histnames[14][1] = "h_individual_cuts_jet2metdphi_cut";   histnames[14][2] = "h_N1_cuts_jet2metdphi";      histnames[14][3] = "h_post_cuts_jet2metdphi";

  histnames2[0][0] = "h_pre_cuts_METphi";     histnames2[0][1] = "h_post_cuts_METphi_final";   
  histnames2[1][0] = "h_pre_cuts_jet1phi";    histnames2[1][1] = "h_post_cuts_jet1phi_final";  
  histnames2[2][0] = "h_pre_cuts_jet2phi";    histnames2[2][1] = "h_post_cuts_jet2phi_final";  
  histnames2[3][0] = "h_selections";          histnames2[3][1] = "h_N1_selections";  

  histtitle[0]  = "N_{jets}";
  histtitle[1]  = "Jet 1 E_{T}";
  histtitle[2]  = "#eta_{J_{1}}";
  histtitle[3]  = "Jet 2 E_{T}";
  histtitle[4]  = "#eta_{J_{2}}";
  histtitle[5]  = "Jet f_{EM}";
  histtitle[6]  = "EM_{frac}^{J_{1}}";
  histtitle[7]  = "EM_{frac}^{J_{2}}";
  histtitle[8]  = "#slashE_{T}";
  histtitle[9]  = "H_{T}";
  histtitle[10]  = "#slashH_{T}";
  histtitle[11]  = "M_{eff}";
  histtitle[12] = "#Delta#phi_{J_{1}J_{2}}";
  histtitle[13] = "#Delta#phi_{J_{1},#slashE_{T}}";
  histtitle[14] = "#Delta#phi_{J_{2},#slashE_{T}}";

  histtitle2[0] = "#slashE_{T}#phi";
  histtitle2[1] = "#phi_{J_{1}}";
  histtitle2[2] = "#phi_{J_{2}}";
  histtitle2[3] = "selections";
  histtitle2[4] = "N-1 selections";
 
  linecolor[0] = kBlue;//LM0
  linecolor[1] = kRed+1;//LM1
  linecolor[2] = kGreen+1;//LM5
//  linecolor[3] = kPink+7;//QCD
//  linecolor[4] = kBlue+3;//TTJets
//  linecolor[5] = kYellow+3;//WJets
//  linecolor[6] = kCyan+3;//ZJets
//  linecolor[7] = kOrange+7;//ZInvisibleJets
  
  for (int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++) {
      mycanvas[z][tt]  = new TCanvas(histnames[z][tt], "", 600,600);}}
  
  for (int ss = 0; ss < 2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++) {
      mycanvas2[y][ss] = new TCanvas(histnames2[y][ss], "", 600,600);}}
  
  for (int j = 0; j < NUMFILES; j++) {
    file[j] = new TFile(filenames[j]+".root");
    
    for (int tt = 0; tt < NUMVARS; tt++) {
      for (int z = 0; z < NUMHISTS; z++) {
	hist[j][z][tt] = (TH1F*)gDirectory->Get(histnames[z][tt]);
	max[z][tt] = (max[z][tt]>hist[j][z][tt]->GetMaximum())?max[z][tt]:hist[j][z][tt]->GetMaximum();}}
    
    for (int ss = 0; ss < 2; ss++) {
      for (int y = 0; y < NUMHISTS2; y++) {
	hist2[j][y][ss] = (TH1F*)gDirectory->Get(histnames2[y][ss]);
	max2[y][ss] = (max2[y][ss]>hist2[j][y][ss]->GetMaximum())?max2[y][ss]:hist2[j][y][ss]->GetMaximum();}}}
  
  for (int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++ ) max[z][tt] *= 1.1;}
  
  for (int ss = 0; ss < 2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++ ) max2[y][ss] *= 1.1;}

  for (int j = 0; j < NUMFILES; j++) {
    for (int tt = 0; tt < NUMVARS; tt++) {
      for (int z = 0; z < NUMHISTS; z++) {
	//hist[j][z][tt]->SetFillStyle(3001);
	hist[j][z][tt]->SetLineColor(linecolor[j]);
	//hist[j][z][tt]->SetLineStyle();
	//hist[j][z][tt]->SetFillColor(linecolor[j]);
	hist[j][z][tt]->SetStats(kFALSE);}}
    
    for (int ss = 0; ss < 2; ss++) {
      for (int y = 0; y < NUMHISTS2; y++) {
	//hist2[j][y][ss]->SetFillStyle(3001);
	hist2[j][y][ss]->SetLineColor(linecolor[j]);
	//hist2[j][y][ss]->SetLineStyle();
	//hist2[j][y][ss]->SetFillColor(linecolor[j]);
	hist2[j][y][ss]->SetStats(kFALSE);}}}

  for(int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++) {
      mycanvas[z][tt]->cd();
      hist[0][z][tt]->SetMaximum(max[z][tt]);
      hist[0][z][tt]->GetXaxis()->SetTitle(histtitle[z]);
      gPad->SetLogy(kTRUE);
      hist[0][z][tt]->Draw();
      for (int j = 1; j < NUMFILES; j++) { hist[j][z][tt]->Draw("same");}}}

  for(int ss = 0; ss <2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++) {
      mycanvas2[y][ss]->cd();
      hist2[0][y][ss]->SetMaximum(max2[y][ss]);
      hist2[0][y][ss]->GetXaxis()->SetTitle(histtitle2[y]);
      gPad->SetLogy(kTRUE);
      hist2[0][y][ss]->Draw();
      for (int j = 1; j < NUMFILES; j++) {
	hist2[j][y][ss]->Draw("same");}}}
  
  leg = new TLegend(.65, .5, .85, .8);
  for (int j = 0; j < NUMFILES; j++) {
    leg->AddEntry(hist[j][0][0], filenames[j],"f");}

  for(int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++) {
      mycanvas[z][tt]->cd();
      leg->Draw();
      mycanvas[z][tt]->SaveAs();}}

  for(int ss = 0; ss < 2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++) {
      mycanvas2[y][ss]->cd();
      leg->Draw();
      mycanvas2[y][ss]->SaveAs();}}

//  bigcanvas1->Divide(2,3);
//  bigcanvas1->cd(1);
//  mycanvas1->DrawClonePad();
//  bigcanvas1->cd(2);
//  mycanvas2->DrawClonePad();
//  bigcanvas1->cd(3);
//  mycanvas3->DrawClonePad();
//  bigcanvas1->cd(4);
//  mycanvas4->DrawClonePad();
//  bigcanvas1->cd(5);
//  mycanvas5->DrawClonePad();
//
//  bigcanvas2->Divide(2,3);
//  bigcanvas2->cd(1);
//  mycanvas6->DrawClonePad();
//  bigcanvas2->cd(2);
//  mycanvas7->DrawClonePad();
//  bigcanvas2->cd(3);
//  mycanvas8->DrawClonePad();
//  bigcanvas2->cd(4);
//  mycanvas9->DrawClonePad();
//  bigcanvas2->cd(5);
//  mycanvas0->DrawClonePad();
//
  //for (int tt = 0; tt < NUMVARS; tt++) for (z = 0; z < NUMHISTS; z++) delete hist[z][tt];
  //for (int ss = 0; ss < 2; ss++) for (y = 0; y < NUMHISTS2; y++) delete hist2[y][ss];
}
