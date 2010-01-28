void makeMyPlots()
{
  bigcanvas1 = new TCanvas("bigcanvas1", "", 1440,900);
  bigcanvas2 = new TCanvas("bigcanvas2", "", 1440,900);
  TCanvas *mycanvas[15];

//  mycanvas2 = new TCanvas("mycanvas2", "", 600,600); 
//  mycanvas3 = new TCanvas("mycanvas3", "", 600,600); 
//  mycanvas4 = new TCanvas("mycanvas4", "", 600,600); 
//  mycanvas5 = new TCanvas("mycanvas5", "", 600,600); 
//  mycanvas6 = new TCanvas("mycanvas6", "", 600,600); 
//  mycanvas7 = new TCanvas("mycanvas7", "", 600,600); 
//  mycanvas8 = new TCanvas("mycanvas8", "", 600,600); 
//  mycanvas9 = new TCanvas("mycanvas9", "", 600,600); 
//  mycanvas0 = new TCanvas("mycanvas0", "", 600,600); 
  mycanvas[] = new TCanvas("mycanvas1", "", 600,600); 
 
  double localpi  = acos(-1);
  double max[15];
  TFile* file[8];
  TString filenames[8];
  TString histnames[15];
  TH1F *hist[8][15];

  filenames[0] = "QCD_Background_cuts.root";
  filenames[1] = "TTJets_cuts.root";
  filenames[2] = "SUSY_LM5_cuts.root";
  filenames[3] = "SUSY_LM1_cuts.root";
  filenames[4] = "SUSY_LM0_cuts.root";
  filenames[5] = "WJets_cuts.root";
  filenames[6] = "ZJets_cuts.root";
  filenames[7] = "ZInvisibleJets_cuts.root";

  histnames[0]  = "h_Njets";
  histnames[1]  = "h_MET";
  histnames[2]  = "h_jet1et";
  histnames[3]  = "h_jet2et";
  histnames[4]  = "h_jet12dphi";
  histnames[5]  = "h_jet1metdphi";
  histnames[6]  = "h_jet2metdphi";
  histnames[7]  = "h_METphi";
  histnames[8]  = "h_jetFem";
  histnames[9]  = "h_HT";
  histnames[10] = "h_MHT";
  histnames[11] = "h_Meff";
  histnames[12] = "h_jet1emfrac";
  histnames[13] = "h_jet2emfrac";
  histnames[14] = "";
  // QCD cuts
  for (int j = 0; j < 8; j++) {
    file[j] = new TFile(filenames[j].c_str());
    for (int z = 0; z < 14; z++ {
      hist[j][z] = (TH1F*)gDirectory->Get(histnames[z].c_str());
      max[z] = (max[z]>hist[j][z]->GetMaximum())?max[z]:hist[j][z]->GetMaximum();
    }
  }
  file[0] = new TFile("QCD_Background_cuts.root"); // open 
  TH1F *hist11 = (TH1F*)gDirectory->Get("h_Njets");
  TH1F *hist12 = (TH1F*)gDirectory->Get("h_MET");
  TH1F *hist13 = (TH1F*)gDirectory->Get("h_jet1et");
  TH1F *hist14 = (TH1F*)gDirectory->Get("h_jet2et");
  TH1F *hist15 = (TH1F*)gDirectory->Get("h_jet12dphi");
  TH1F *hist16 = (TH1F*)gDirectory->Get("h_jet1metdphi");
  TH1F *hist17 = (TH1F*)gDirectory->Get("h_jet2metdphi");
  TH1F *hist18 = (TH1F*)gDirectory->Get("h_METphi");
  TH1F *hist19 = (TH1F*)gDirectory->Get("h_jetFem");
  TH1F *hist10 = (TH1F*)gDirectory->Get("h_MHT");

  max1 = hist11->GetMaximum();
  max2 = hist12->GetMaximum();
  max3 = hist13->GetMaximum();
  max4 = hist14->GetMaximum();
  max5 = hist15->GetMaximum();
  max6 = hist16->GetMaximum();
  max7 = hist17->GetMaximum();
  max8 = hist18->GetMaximum();
  max9 = hist19->GetMaximum();
  max0 = hist10->GetMaximum();

  // TT jets cuts
  TFile* f2 = new TFile("TTJets_cuts.root"); // open the file
  TH1F *hist21 = (TH1F*)gDirectory->Get("h_Njets");
  TH1F *hist22 = (TH1F*)gDirectory->Get("h_MET");
  TH1F *hist23 = (TH1F*)gDirectory->Get("h_jet1et");
  TH1F *hist24 = (TH1F*)gDirectory->Get("h_jet2et");
  TH1F *hist25 = (TH1F*)gDirectory->Get("h_jet12dphi");
  TH1F *hist26 = (TH1F*)gDirectory->Get("h_jet1metdphi");
  TH1F *hist27 = (TH1F*)gDirectory->Get("h_jet2metdphi");
  TH1F *hist28 = (TH1F*)gDirectory->Get("h_METphi");
  TH1F *hist29 = (TH1F*)gDirectory->Get("h_jetFem");
  TH1F *hist20 = (TH1F*)gDirectory->Get("h_MHT");

  max1 = (max1>hist21->GetMaximum())?max1:hist21->GetMaximum();
  max2 = (max2>hist22->GetMaximum())?max2:hist22->GetMaximum();
  max3 = (max3>hist23->GetMaximum())?max3:hist23->GetMaximum();
  max4 = (max4>hist24->GetMaximum())?max4:hist24->GetMaximum();
  max5 = (max5>hist25->GetMaximum())?max5:hist25->GetMaximum();
  max6 = (max6>hist26->GetMaximum())?max6:hist26->GetMaximum();
  max7 = (max7>hist27->GetMaximum())?max7:hist27->GetMaximum();
  max8 = (max8>hist28->GetMaximum())?max8:hist28->GetMaximum();
  max9 = (max9>hist29->GetMaximum())?max9:hist29->GetMaximum();
  max0 = (max0>hist20->GetMaximum())?max0:hist20->GetMaximum();

  // LM5 cuts
  TFile* f3 = new TFile("SUSY_LM5_cuts.root"); // open the file
  TH1F *hist31 = (TH1F*)gDirectory->Get("h_Njets");
  TH1F *hist32 = (TH1F*)gDirectory->Get("h_MET");
  TH1F *hist33 = (TH1F*)gDirectory->Get("h_jet1et");
  TH1F *hist34 = (TH1F*)gDirectory->Get("h_jet2et");
  TH1F *hist35 = (TH1F*)gDirectory->Get("h_jet12dphi");
  TH1F *hist36 = (TH1F*)gDirectory->Get("h_jet1metdphi");
  TH1F *hist37 = (TH1F*)gDirectory->Get("h_jet2metdphi");
  TH1F *hist38 = (TH1F*)gDirectory->Get("h_METphi");
  TH1F *hist39 = (TH1F*)gDirectory->Get("h_jetFem");
  TH1F *hist30 = (TH1F*)gDirectory->Get("h_MHT");

  max1 = (max1>hist31->GetMaximum())?max1:hist31->GetMaximum();
  max2 = (max2>hist32->GetMaximum())?max2:hist32->GetMaximum();
  max3 = (max3>hist33->GetMaximum())?max3:hist33->GetMaximum();
  max4 = (max4>hist34->GetMaximum())?max4:hist34->GetMaximum();
  max5 = (max5>hist35->GetMaximum())?max5:hist35->GetMaximum();
  max6 = (max6>hist36->GetMaximum())?max6:hist36->GetMaximum();
  max7 = (max7>hist37->GetMaximum())?max7:hist37->GetMaximum();
  max8 = (max8>hist38->GetMaximum())?max8:hist38->GetMaximum();
  max9 = (max9>hist39->GetMaximum())?max9:hist39->GetMaximum();
  max0 = (max0>hist30->GetMaximum())?max0:hist30->GetMaximum();

  // LM1 cuts
  TFile* f4 = new TFile("SUSY_LM1_cuts.root"); // open the file
  TH1F *hist41 = (TH1F*)gDirectory->Get("h_Njets");
  TH1F *hist42 = (TH1F*)gDirectory->Get("h_MET");
  TH1F *hist43 = (TH1F*)gDirectory->Get("h_jet1et");
  TH1F *hist44 = (TH1F*)gDirectory->Get("h_jet2et");
  TH1F *hist45 = (TH1F*)gDirectory->Get("h_jet12dphi");
  TH1F *hist46 = (TH1F*)gDirectory->Get("h_jet1metdphi");
  TH1F *hist47 = (TH1F*)gDirectory->Get("h_jet2metdphi");
  TH1F *hist48 = (TH1F*)gDirectory->Get("h_METphi");
  TH1F *hist49 = (TH1F*)gDirectory->Get("h_jetFem");
  TH1F *hist40 = (TH1F*)gDirectory->Get("h_MHT");

  max1 = (max1>hist41->GetMaximum())?max1:hist41->GetMaximum();
  max2 = (max2>hist42->GetMaximum())?max2:hist42->GetMaximum();
  max3 = (max3>hist43->GetMaximum())?max3:hist43->GetMaximum();
  max4 = (max4>hist44->GetMaximum())?max4:hist44->GetMaximum();
  max5 = (max5>hist45->GetMaximum())?max5:hist45->GetMaximum();
  max6 = (max6>hist46->GetMaximum())?max6:hist46->GetMaximum();
  max7 = (max7>hist47->GetMaximum())?max7:hist47->GetMaximum();
  max8 = (max8>hist48->GetMaximum())?max8:hist48->GetMaximum();
  max9 = (max9>hist49->GetMaximum())?max9:hist49->GetMaximum();
  max0 = (max0>hist40->GetMaximum())?max0:hist40->GetMaximum();

  // LM0 cuts
  TFile* f5 = new TFile("SUSY_LM0_cuts.root"); // open the file
  TH1F *hist51 = (TH1F*)gDirectory->Get("h_Njets");
  TH1F *hist52 = (TH1F*)gDirectory->Get("h_MET");
  TH1F *hist53 = (TH1F*)gDirectory->Get("h_jet1et");
  TH1F *hist54 = (TH1F*)gDirectory->Get("h_jet2et");
  TH1F *hist55 = (TH1F*)gDirectory->Get("h_jet12dphi");
  TH1F *hist56 = (TH1F*)gDirectory->Get("h_jet1metdphi");
  TH1F *hist57 = (TH1F*)gDirectory->Get("h_jet2metdphi");
  TH1F *hist58 = (TH1F*)gDirectory->Get("h_METphi");
  TH1F *hist59 = (TH1F*)gDirectory->Get("h_jetFem");
  TH1F *hist50 = (TH1F*)gDirectory->Get("h_MHT");

  max1 = (max1>hist51->GetMaximum())?max1:hist51->GetMaximum();
  max2 = (max2>hist52->GetMaximum())?max2:hist52->GetMaximum();
  max3 = (max3>hist53->GetMaximum())?max3:hist53->GetMaximum();
  max4 = (max4>hist54->GetMaximum())?max4:hist54->GetMaximum();
  max5 = (max5>hist55->GetMaximum())?max5:hist55->GetMaximum();
  max6 = (max6>hist56->GetMaximum())?max6:hist56->GetMaximum();
  max7 = (max7>hist57->GetMaximum())?max7:hist57->GetMaximum();
  max8 = (max8>hist58->GetMaximum())?max8:hist58->GetMaximum();
  max9 = (max9>hist59->GetMaximum())?max9:hist59->GetMaximum();
  max0 = (max0>hist50->GetMaximum())?max0:hist50->GetMaximum();


  //WJets
  TFile* f6 = new TFile("WJets_cuts.root"); // open the file
  TH1F *hist61 = (TH1F*)gDirectory->Get("h_Njets");
  TH1F *hist62 = (TH1F*)gDirectory->Get("h_MET");
  TH1F *hist63 = (TH1F*)gDirectory->Get("h_jet1et");
  TH1F *hist64 = (TH1F*)gDirectory->Get("h_jet2et");
  TH1F *hist65 = (TH1F*)gDirectory->Get("h_jet12dphi");
  TH1F *hist66 = (TH1F*)gDirectory->Get("h_jet1metdphi");
  TH1F *hist67 = (TH1F*)gDirectory->Get("h_jet2metdphi");
  TH1F *hist68 = (TH1F*)gDirectory->Get("h_METphi");
  TH1F *hist69 = (TH1F*)gDirectory->Get("h_jetFem");
  TH1F *hist60 = (TH1F*)gDirectory->Get("h_MHT");

  max1 = (max1>hist61->GetMaximum())?max1:hist61->GetMaximum();
  max2 = (max2>hist62->GetMaximum())?max2:hist62->GetMaximum();
  max3 = (max3>hist63->GetMaximum())?max3:hist63->GetMaximum();
  max4 = (max4>hist64->GetMaximum())?max4:hist64->GetMaximum();
  max5 = (max5>hist65->GetMaximum())?max5:hist65->GetMaximum();
  max6 = (max6>hist66->GetMaximum())?max6:hist66->GetMaximum();
  max7 = (max7>hist67->GetMaximum())?max7:hist67->GetMaximum();
  max8 = (max8>hist68->GetMaximum())?max8:hist68->GetMaximum();
  max9 = (max9>hist69->GetMaximum())?max9:hist69->GetMaximum();
  max0 = (max0>hist60->GetMaximum())?max0:hist60->GetMaximum();

  //ZJets
  TFile* f7 = new TFile("ZJets_cuts.root"); // open the file
  TH1F *hist71 = (TH1F*)gDirectory->Get("h_Njets");
  TH1F *hist72 = (TH1F*)gDirectory->Get("h_MET");
  TH1F *hist73 = (TH1F*)gDirectory->Get("h_jet1et");
  TH1F *hist74 = (TH1F*)gDirectory->Get("h_jet2et");
  TH1F *hist75 = (TH1F*)gDirectory->Get("h_jet12dphi");
  TH1F *hist76 = (TH1F*)gDirectory->Get("h_jet1metdphi");
  TH1F *hist77 = (TH1F*)gDirectory->Get("h_jet2metdphi");
  TH1F *hist78 = (TH1F*)gDirectory->Get("h_METphi");
  TH1F *hist79 = (TH1F*)gDirectory->Get("h_jetFem");
  TH1F *hist70 = (TH1F*)gDirectory->Get("h_MHT");

  max1 = (max1>hist71->GetMaximum())?max1:hist71->GetMaximum();
  max2 = (max2>hist72->GetMaximum())?max2:hist72->GetMaximum();
  max3 = (max3>hist73->GetMaximum())?max3:hist73->GetMaximum();
  max4 = (max4>hist74->GetMaximum())?max4:hist74->GetMaximum();
  max5 = (max5>hist75->GetMaximum())?max5:hist75->GetMaximum();
  max6 = (max6>hist76->GetMaximum())?max6:hist76->GetMaximum();
  max7 = (max7>hist77->GetMaximum())?max7:hist77->GetMaximum();
  max8 = (max8>hist78->GetMaximum())?max8:hist78->GetMaximum();
  max9 = (max9>hist79->GetMaximum())?max9:hist79->GetMaximum();
  max0 = (max0>hist70->GetMaximum())?max0:hist70->GetMaximum();

  //ZInvisibleJets
  TFile* f8 = new TFile("ZInvisibleJets_cuts.root"); // open the file
  TH1F *hist81 = (TH1F*)gDirectory->Get("h_Njets");
  TH1F *hist82 = (TH1F*)gDirectory->Get("h_MET");
  TH1F *hist83 = (TH1F*)gDirectory->Get("h_jet1et");
  TH1F *hist84 = (TH1F*)gDirectory->Get("h_jet2et");
  TH1F *hist85 = (TH1F*)gDirectory->Get("h_jet12dphi");
  TH1F *hist86 = (TH1F*)gDirectory->Get("h_jet1metdphi");
  TH1F *hist87 = (TH1F*)gDirectory->Get("h_jet2metdphi");
  TH1F *hist88 = (TH1F*)gDirectory->Get("h_METphi");
  TH1F *hist89 = (TH1F*)gDirectory->Get("h_jetFem");
  TH1F *hist80 = (TH1F*)gDirectory->Get("h_MHT");

  max1 = (max1>hist81->GetMaximum())?max1:hist81->GetMaximum();
  max2 = (max2>hist82->GetMaximum())?max2:hist82->GetMaximum();
  max3 = (max3>hist83->GetMaximum())?max3:hist83->GetMaximum();
  max4 = (max4>hist84->GetMaximum())?max4:hist84->GetMaximum();
  max5 = (max5>hist85->GetMaximum())?max5:hist85->GetMaximum();
  max6 = (max6>hist86->GetMaximum())?max6:hist86->GetMaximum();
  max7 = (max7>hist87->GetMaximum())?max7:hist87->GetMaximum();
  max8 = (max8>hist88->GetMaximum())?max8:hist88->GetMaximum();
  max9 = (max9>hist89->GetMaximum())?max9:hist89->GetMaximum();
  max0 = (max0>hist80->GetMaximum())?max0:hist80->GetMaximum();


  max1 *= 1.05;
  max2 *= 1.05;
  max3 *= 1.05;
  max4 *= 1.05;
  max5 *= 1.05;
  max6 *= 1.05;
  max7 *= 1.05;
  max8 *= 1.05;
  max9 *= 1.05;
  max0 *= 1.05;

  hist11->SetFillStyle(3001);    hist21->SetFillStyle(3001);    hist31->SetFillStyle(3001);     hist41->SetFillStyle(3001);
  hist11->SetLineColor(kRed);    hist21->SetLineColor(kBlue);   hist31->SetLineColor(kGreen);   hist41->SetLineColor(kBlack);
  hist11->SetStats(kFALSE);      hist21->SetStats(kFALSE);      hist31->SetStats(kFALSE);       hist41->SetStats(kFALSE);
  hist51->SetFillStyle(3001);    hist61->SetFillStyle(3001);    hist71->SetFillStyle(3001);     hist81->SetFillStyle(3001);
  hist51->SetLineColor(kRed);    hist61->SetLineColor(kBlue);   hist71->SetLineColor(kGreen);   hist81->SetLineColor(kBlack);
  hist51->SetStats(kFALSE);      hist61->SetStats(kFALSE);      hist71->SetStats(kFALSE);       hist81->SetStats(kFALSE);

  hist12->SetFillStyle(3001);    hist22->SetFillStyle(3001);    hist32->SetFillStyle(3001);     hist42->SetFillStyle(3001);
  hist12->SetLineColor(kRed);    hist22->SetLineColor(kBlue);   hist32->SetLineColor(kGreen);   hist42->SetLineColor(kBlack);
  hist12->SetStats(kFALSE);      hist22->SetStats(kFALSE);      hist32->SetStats(kFALSE);       hist42->SetStats(kFALSE);
  hist52->SetFillStyle(3001);    hist62->SetFillStyle(3001);    hist72->SetFillStyle(3001);     hist82->SetFillStyle(3001);
  hist52->SetLineColor(kRed);    hist62->SetLineColor(kBlue);   hist72->SetLineColor(kGreen);   hist82->SetLineColor(kBlack);
  hist52->SetStats(kFALSE);      hist62->SetStats(kFALSE);      hist72->SetStats(kFALSE);       hist82->SetStats(kFALSE);

  hist13->SetFillStyle(3001);    hist23->SetFillStyle(3001);    hist33->SetFillStyle(3001);     hist43->SetFillStyle(3001);
  hist13->SetLineColor(kRed);    hist23->SetLineColor(kBlue);   hist33->SetLineColor(kGreen);   hist43->SetLineColor(kBlack);
  hist13->SetStats(kFALSE);      hist23->SetStats(kFALSE);      hist33->SetStats(kFALSE);       hist43->SetStats(kFALSE);
  hist53->SetFillStyle(3001);    hist63->SetFillStyle(3001);    hist73->SetFillStyle(3001);     hist83->SetFillStyle(3001);
  hist53->SetLineColor(kRed);    hist63->SetLineColor(kBlue);   hist73->SetLineColor(kGreen);   hist83->SetLineColor(kBlack);
  hist53->SetStats(kFALSE);      hist63->SetStats(kFALSE);      hist73->SetStats(kFALSE);       hist83->SetStats(kFALSE);

  hist14->SetFillStyle(3001);    hist24->SetFillStyle(3001);    hist34->SetFillStyle(3001);     hist44->SetFillStyle(3001);
  hist14->SetLineColor(kRed);    hist24->SetLineColor(kBlue);   hist34->SetLineColor(kGreen);   hist44->SetLineColor(kBlack);
  hist14->SetStats(kFALSE);      hist24->SetStats(kFALSE);      hist34->SetStats(kFALSE);       hist44->SetStats(kFALSE);
  hist54->SetFillStyle(3001);    hist64->SetFillStyle(3001);    hist74->SetFillStyle(3001);     hist84->SetFillStyle(3001);
  hist54->SetLineColor(kRed);    hist64->SetLineColor(kBlue);   hist74->SetLineColor(kGreen);   hist84->SetLineColor(kBlack);
  hist54->SetStats(kFALSE);      hist64->SetStats(kFALSE);      hist74->SetStats(kFALSE);       hist84->SetStats(kFALSE);

  hist15->SetFillStyle(3001);    hist25->SetFillStyle(3001);    hist35->SetFillStyle(3001);     hist45->SetFillStyle(3001);
  hist15->SetLineColor(kRed);    hist25->SetLineColor(kBlue);   hist35->SetLineColor(kGreen);   hist45->SetLineColor(kBlack);
  hist15->SetStats(kFALSE);      hist25->SetStats(kFALSE);      hist35->SetStats(kFALSE);       hist45->SetStats(kFALSE);
  hist55->SetFillStyle(3001);    hist65->SetFillStyle(3001);    hist75->SetFillStyle(3001);     hist85->SetFillStyle(3001);
  hist55->SetLineColor(kRed);    hist65->SetLineColor(kBlue);   hist75->SetLineColor(kGreen);   hist85->SetLineColor(kBlack);
  hist55->SetStats(kFALSE);      hist65->SetStats(kFALSE);      hist75->SetStats(kFALSE);       hist85->SetStats(kFALSE);

  hist16->SetFillStyle(3001);    hist26->SetFillStyle(3001);    hist36->SetFillStyle(3001);     hist46->SetFillStyle(3001);
  hist16->SetLineColor(kRed);    hist26->SetLineColor(kBlue);   hist36->SetLineColor(kGreen);   hist46->SetLineColor(kBlack);
  hist16->SetStats(kFALSE);      hist26->SetStats(kFALSE);      hist36->SetStats(kFALSE);       hist46->SetStats(kFALSE);
  hist56->SetFillStyle(3001);    hist66->SetFillStyle(3001);    hist76->SetFillStyle(3001);     hist86->SetFillStyle(3001);
  hist56->SetLineColor(kRed);    hist66->SetLineColor(kBlue);   hist76->SetLineColor(kGreen);   hist86->SetLineColor(kBlack);
  hist56->SetStats(kFALSE);      hist66->SetStats(kFALSE);      hist76->SetStats(kFALSE);       hist86->SetStats(kFALSE);

  hist17->SetFillStyle(3001);    hist27->SetFillStyle(3001);    hist37->SetFillStyle(3001);     hist47->SetFillStyle(3001);
  hist17->SetLineColor(kRed);    hist27->SetLineColor(kBlue);   hist37->SetLineColor(kGreen);   hist47->SetLineColor(kBlack);
  hist17->SetStats(kFALSE);      hist27->SetStats(kFALSE);      hist37->SetStats(kFALSE);       hist47->SetStats(kFALSE);
  hist57->SetFillStyle(3001);    hist67->SetFillStyle(3001);    hist77->SetFillStyle(3001);     hist87->SetFillStyle(3001);
  hist57->SetLineColor(kRed);    hist67->SetLineColor(kBlue);   hist77->SetLineColor(kGreen);   hist87->SetLineColor(kBlack);
  hist57->SetStats(kFALSE);      hist67->SetStats(kFALSE);      hist77->SetStats(kFALSE);       hist87->SetStats(kFALSE);

  hist18->SetFillStyle(3001);    hist28->SetFillStyle(3001);    hist38->SetFillStyle(3001);     hist48->SetFillStyle(3001);
  hist18->SetLineColor(kRed);    hist28->SetLineColor(kBlue);   hist38->SetLineColor(kGreen);   hist48->SetLineColor(kBlack);
  hist18->SetStats(kFALSE);      hist28->SetStats(kFALSE);      hist38->SetStats(kFALSE);       hist48->SetStats(kFALSE);
  hist58->SetFillStyle(3001);    hist68->SetFillStyle(3001);    hist78->SetFillStyle(3001);     hist88->SetFillStyle(3001);
  hist58->SetLineColor(kRed);    hist68->SetLineColor(kBlue);   hist78->SetLineColor(kGreen);   hist88->SetLineColor(kBlack);
  hist58->SetStats(kFALSE);      hist68->SetStats(kFALSE);      hist78->SetStats(kFALSE);       hist88->SetStats(kFALSE);

  hist19->SetFillStyle(3001);    hist29->SetFillStyle(3001);    hist39->SetFillStyle(3001);     hist49->SetFillStyle(3001);
  hist19->SetLineColor(kRed);    hist29->SetLineColor(kBlue);   hist39->SetLineColor(kGreen);   hist49->SetLineColor(kBlack);
  hist19->SetStats(kFALSE);      hist29->SetStats(kFALSE);      hist39->SetStats(kFALSE);       hist49->SetStats(kFALSE);
  hist59->SetFillStyle(3001);    hist69->SetFillStyle(3001);    hist79->SetFillStyle(3001);     hist89->SetFillStyle(3001);
  hist59->SetLineColor(kRed);    hist69->SetLineColor(kBlue);   hist79->SetLineColor(kGreen);   hist89->SetLineColor(kBlack);
  hist59->SetStats(kFALSE);      hist69->SetStats(kFALSE);      hist79->SetStats(kFALSE);       hist89->SetStats(kFALSE);

  hist10->SetFillStyle(3001);    hist20->SetFillStyle(3001);    hist30->SetFillStyle(3001);     hist40->SetFillStyle(3001);
  hist10->SetLineColor(kRed);    hist20->SetLineColor(kBlue);   hist30->SetLineColor(kGreen);   hist40->SetLineColor(kBlack);
  hist10->SetStats(kFALSE);      hist20->SetStats(kFALSE);      hist30->SetStats(kFALSE);       hist40->SetStats(kFALSE);
  hist50->SetFillStyle(3001);    hist60->SetFillStyle(3001);    hist70->SetFillStyle(3001);     hist80->SetFillStyle(3001);
  hist50->SetLineColor(kRed);    hist60->SetLineColor(kBlue);   hist70->SetLineColor(kGreen);   hist80->SetLineColor(kBlack);
  hist50->SetStats(kFALSE);      hist60->SetStats(kFALSE);      hist70->SetStats(kFALSE);       hist80->SetStats(kFALSE);

  
  mycanvas1->cd();
  hist11->SetMaximum(max1);
  hist11->GetXaxis()->SetTitle("N_{Jets}");
  gPad->SetLogy(kTRUE);
  hist11->Draw();
  hist21->Draw("same");
  hist31->Draw("same");
  hist41->Draw("same");
  hist51->Draw("same");
  hist61->Draw("same");
  hist71->Draw("same");
  hist81->Draw("same");

  mycanvas2->cd();
  hist12->SetMaximum(max2);
  hist12->GetXaxis()->SetTitle("#slashE_{T}");
  gPad->SetLogy(kTRUE);
  hist12->Draw();
  hist22->Draw("same");
  hist32->Draw("same");
  hist42->Draw("same");
  hist52->Draw("same");
  hist62->Draw("same");
  hist72->Draw("same");
  hist82->Draw("same");

  mycanvas3->cd();
  hist13->SetMaximum(max3);
  hist13->GetXaxis()->SetTitle("Jet 1 E_{T}");
  gPad->SetLogy(kTRUE);
  hist13->Draw();
  hist23->Draw("same");
  hist33->Draw("same");
  hist43->Draw("same");
  hist53->Draw("same");
  hist63->Draw("same");
  hist73->Draw("same");
  hist83->Draw("same");

  mycanvas4->cd();
  hist14->SetMaximum(max4);
  hist14->GetXaxis()->SetTitle("Jet 2 E_{T}");
  gPad->SetLogy(kTRUE);
  hist14->Draw();
  hist24->Draw("same");
  hist34->Draw("same");
  hist44->Draw("same");
  hist54->Draw("same");
  hist64->Draw("same");
  hist74->Draw("same");
  hist84->Draw("same");

  mycanvas5->cd();
  hist15->SetMaximum(max5);
  hist15->GetXaxis()->SetTitle("#Delta#phi_{J_{1}J_{2}}");
  gPad->SetLogy(kTRUE);
  hist15->Draw();
  hist25->Draw("same");
  hist35->Draw("same");
  hist45->Draw("same");
  hist55->Draw("same");
  hist65->Draw("same");
  hist75->Draw("same");
  hist85->Draw("same");

  mycanvas6->cd();
  hist16->SetMaximum(max6);
  hist16->GetXaxis()->SetTitle("#Delta#phi_{J_{1},#slashE_{T}}");
  gPad->SetLogy(kTRUE);
  hist16->Draw();
  hist26->Draw("same");
  hist36->Draw("same");
  hist46->Draw("same");
  hist56->Draw("same");
  hist66->Draw("same");
  hist76->Draw("same");
  hist86->Draw("same");

  mycanvas7->cd();
  hist17->SetMaximum(max7);
  hist17->GetXaxis()->SetTitle("#Delta#phi_{J_{2},#slashE_{T}}");
  gPad->SetLogy(kTRUE);
  hist17->Draw();
  hist27->Draw("same");
  hist37->Draw("same");
  hist47->Draw("same");
  hist57->Draw("same");
  hist67->Draw("same");
  hist77->Draw("same");
  hist87->Draw("same");

  mycanvas8->cd();
  hist18->SetMaximum(max8);
  hist18->GetXaxis()->SetTitle("#slashE_{T}#phi");
  gPad->SetLogy(kTRUE);
  hist18->Draw();
  hist28->Draw("same");
  hist38->Draw("same");
  hist48->Draw("same");
  hist58->Draw("same");
  hist68->Draw("same");
  hist78->Draw("same");
  hist88->Draw("same");

  mycanvas9->cd();
  hist19->SetMaximum(max9);
  hist19->GetXaxis()->SetTitle("Jet f_{EM}");
  gPad->SetLogy(kTRUE);
  hist19->Draw();
  hist29->Draw("same");
  hist39->Draw("same");
  hist49->Draw("same");
  hist59->Draw("same");
  hist69->Draw("same");
  hist79->Draw("same");
  hist89->Draw("same");

  mycanvas0->cd();
  hist10->SetMaximum(max9);
  hist12->GetXaxis()->SetTitle("#slashH_{T}");
  gPad->SetLogy(kTRUE);
  hist10->Draw();
  hist20->Draw("same");
  hist30->Draw("same");
  hist40->Draw("same");
  hist50->Draw("same");
  hist60->Draw("same");
  hist70->Draw("same");
  hist80->Draw("same");


  TLegend *leg1 = new TLegend(.65, .5, .85, .8);
  leg1->AddEntry(hist11, "QCD Background", "f");
  leg1->AddEntry(hist21, "TT Jets", "f");
  leg1->AddEntry(hist31, "LM5", "f");
  leg1->AddEntry(hist41, "LM1", "f");
  leg1->AddEntry(hist51, "LM0", "f");
  leg1->AddEntry(hist61, "W Jets", "f");
  leg1->AddEntry(hist71, "Z Jets", "f");
  leg1->AddEntry(hist81, "ZInvisible Jets", "f");

  TLegend *leg2 = new TLegend(.65, .5, .85, .8);
  leg2->AddEntry(hist12, "QCD Background", "f");
  leg2->AddEntry(hist22, "TT Jets", "f");
  leg2->AddEntry(hist32, "LM5", "f");
  leg2->AddEntry(hist42, "LM1", "f");
  leg2->AddEntry(hist52, "LM0", "f");
  leg2->AddEntry(hist62, "W Jets", "f");
  leg2->AddEntry(hist72, "Z Jets", "f");
  leg2->AddEntry(hist82, "ZInvisible Jets", "f");

  TLegend *leg3 = new TLegend(.65, .5, .85, .8);
  leg3->AddEntry(hist13, "QCD Background", "f");
  leg3->AddEntry(hist23, "TT Jets", "f");
  leg3->AddEntry(hist33, "LM5", "f");
  leg3->AddEntry(hist43, "LM1", "f");
  leg3->AddEntry(hist53, "LM0", "f");
  leg3->AddEntry(hist63, "W Jets", "f");
  leg3->AddEntry(hist73, "Z Jets", "f");
  leg3->AddEntry(hist83, "ZInvisible Jets", "f");

  TLegend *leg4 = new TLegend(.65, .5, .85, .8);
  leg4->AddEntry(hist14, "QCD Background", "f");
  leg4->AddEntry(hist24, "TT Jets", "f");
  leg4->AddEntry(hist34, "LM5", "f");
  leg4->AddEntry(hist44, "LM1", "f");
  leg4->AddEntry(hist54, "LM0", "f");
  leg4->AddEntry(hist64, "W Jets", "f");
  leg4->AddEntry(hist74, "Z Jets", "f");
  leg4->AddEntry(hist84, "ZInvisible Jets", "f");

  TLegend *leg5 = new TLegend(.65, .5, .85, .8);
  leg5->AddEntry(hist15, "QCD Background", "f");
  leg5->AddEntry(hist25, "TT Jets", "f");
  leg5->AddEntry(hist35, "LM5", "f");
  leg5->AddEntry(hist45, "LM1", "f");
  leg5->AddEntry(hist55, "LM0", "f");
  leg5->AddEntry(hist65, "W Jets", "f");
  leg5->AddEntry(hist75, "Z Jets", "f");
  leg5->AddEntry(hist85, "ZInvisible Jets", "f");

  TLegend *leg6 = new TLegend(.65, .5, .85, .8);
  leg6->AddEntry(hist16, "QCD Background", "f");
  leg6->AddEntry(hist26, "TT Jets", "f");
  leg6->AddEntry(hist36, "LM5", "f");
  leg6->AddEntry(hist46, "LM1", "f");
  leg6->AddEntry(hist56, "LM0", "f");
  leg6->AddEntry(hist66, "W Jets", "f");
  leg6->AddEntry(hist76, "Z Jets", "f");
  leg6->AddEntry(hist86, "ZInvisible Jets", "f");

  TLegend *leg7 = new TLegend(.65, .5, .85, .8);
  leg7->AddEntry(hist17, "QCD Background", "f");
  leg7->AddEntry(hist27, "TT Jets", "f");
  leg7->AddEntry(hist37, "LM5", "f");
  leg7->AddEntry(hist47, "LM1", "f");
  leg7->AddEntry(hist57, "LM0", "f");
  leg7->AddEntry(hist67, "W Jets", "f");
  leg7->AddEntry(hist77, "Z Jets", "f");
  leg7->AddEntry(hist87, "ZInvisible Jets", "f");

  TLegend *leg8 = new TLegend(.65, .5, .85, .8);
  leg8->AddEntry(hist18, "QCD Background", "f");
  leg8->AddEntry(hist28, "TT Jets", "f");
  leg8->AddEntry(hist38, "LM5", "f");
  leg8->AddEntry(hist48, "LM1", "f");
  leg8->AddEntry(hist58, "LM0", "f");
  leg8->AddEntry(hist68, "W Jets", "f");
  leg8->AddEntry(hist78, "Z Jets", "f");
  leg8->AddEntry(hist88, "ZInvisible Jets", "f");

  TLegend *leg9 = new TLegend(.65, .5, .85, .8);
  leg9->AddEntry(hist19, "QCD Background", "f");
  leg9->AddEntry(hist29, "TT Jets", "f");
  leg9->AddEntry(hist39, "LM5", "f");
  leg9->AddEntry(hist49, "LM1", "f");
  leg9->AddEntry(hist59, "LM0", "f");
  leg9->AddEntry(hist69, "W Jets", "f");
  leg9->AddEntry(hist79, "Z Jets", "f");
  leg9->AddEntry(hist89, "ZInvisible Jets", "f");

  TLegend *leg0 = new TLegend(.65, .5, .85, .8);
  leg0->AddEntry(hist10, "QCD Background", "f");
  leg0->AddEntry(hist20, "TT Jets", "f");
  leg0->AddEntry(hist30, "LM5", "f");
  leg0->AddEntry(hist40, "LM1", "f");
  leg0->AddEntry(hist50, "LM0", "f");
  leg0->AddEntry(hist60, "W Jets", "f");
  leg0->AddEntry(hist70, "Z Jets", "f");
  leg0->AddEntry(hist80, "ZInvisible Jets", "f");

  mycanvas1->cd();
  leg1->Draw();
  mycanvas2->cd();
  leg2->Draw();
  mycanvas3->cd();
  leg3->Draw();
  mycanvas4->cd();
  leg4->Draw();
  mycanvas5->cd();
  leg5->Draw();
  mycanvas6->cd();
  leg6->Draw();
  mycanvas7->cd();
  leg7->Draw();
  mycanvas8->cd();
  leg8->Draw();
  mycanvas9->cd();
  leg9->Draw();
  mycanvas0->cd();
  leg0->Draw();

  bigcanvas1->Divide(2,3);
  bigcanvas1->cd(1);
  mycanvas1->DrawClonePad();
  bigcanvas1->cd(2);
  mycanvas2->DrawClonePad();
  bigcanvas1->cd(3);
  mycanvas3->DrawClonePad();
  bigcanvas1->cd(4);
  mycanvas4->DrawClonePad();
  bigcanvas1->cd(5);
  mycanvas5->DrawClonePad();

  bigcanvas2->Divide(2,3);
  bigcanvas2->cd(1);
  mycanvas6->DrawClonePad();
  bigcanvas2->cd(2);
  mycanvas7->DrawClonePad();
  bigcanvas2->cd(3);
  mycanvas8->DrawClonePad();
  bigcanvas2->cd(4);
  mycanvas9->DrawClonePad();
  bigcanvas2->cd(5);
  mycanvas0->DrawClonePad();

  delete mycanvas1;
  delete mycanvas2;
  delete mycanvas3;
  delete mycanvas4;
  delete mycanvas5;
  delete mycanvas6;
  delete mycanvas7;
  delete mycanvas8;
  delete mycanvas9;
  delete mycanvas0;
}
