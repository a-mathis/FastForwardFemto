
std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7};
std::vector<int> fColors     = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2};
std::vector<int> fMarkers    = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette=false, bool title=false)
{
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);

  const int NRGBs = 6;
  Double_t stops[NRGBs];
  for(int i=0; i<NRGBs; ++i) stops[i] = float(i)/(NRGBs-1);

  Double_t red[NRGBs]   = { 1.,  29./255., 25./255., 27./255., 32./255., 24./255.};
  Double_t green[NRGBs] = { 1., 221./255., 160./255., 113./255., 74./255., 37./255.};
  Double_t blue[NRGBs] = {  1., 221./255., 184./255., 154./255., 129./255., 98./255.};
  TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHisto(TH1 *histo, int marker, int color)
{
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* Calculate_CF(TH1F* histRE_relK,TH1F* histME_relK, TString CFname,Double_t normleft,Double_t normright)
{
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  Double_t norm_relK = histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) / histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright));

  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  Hist_CF->Divide(histRE_relK,histME_relK,1,norm_relK);

  return Hist_CF;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* add_CF(TH1F* hist_CF1, TH1F* hist_CF2,TString HistName)
{
  //Calculate CFs with error weighting
  TH1F* hist_CF_sum = (TH1F*)hist_CF1->Clone(HistName.Data());

  Int_t NBins = hist_CF_sum->GetNbinsX();

  for(Int_t i=0;i<NBins;i++)
  {
    double CF1_val = hist_CF1->GetBinContent(i+1);
    double CF1_err = hist_CF1->GetBinError(i+1);
    double CF2_val = hist_CF2->GetBinContent(i+1);
    double CF2_err = hist_CF2->GetBinError(i+1);

    //average for bin i:
    if(CF1_val != 0. && CF2_val != 0.)
    {
      double CF1_err_weight = 1./TMath::Power(CF1_err,2.);
      double CF2_err_weight = 1./TMath::Power(CF2_err,2.);

      double CF_sum_average = (CF1_err_weight*CF1_val + CF2_err_weight*CF2_val)/(CF1_err_weight+CF2_err_weight);
      double CF_sum_err = 1./TMath::Sqrt(CF1_err_weight+CF2_err_weight);

      hist_CF_sum->SetBinContent(i+1,CF_sum_average);
      hist_CF_sum->SetBinError(i+1,CF_sum_err);
    }
    else if(CF1_val == 0. && CF2_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF2_val);
      hist_CF_sum->SetBinError(i+1,CF2_err);
    }
    else if(CF2_val == 0 && CF1_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF1_val);
      hist_CF_sum->SetBinError(i+1,CF1_err);
    }

  }
  return hist_CF_sum;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *DrawSystematicError(TH1F* histexp,TH1F *histerr,double errorwidth)
{
  //Input are the experimental histogram with statistical errors only and a histogram containing the errors

  const int histbins = histexp->GetNbinsX();

  TGraphErrors *ge_SysError_C2 = new TGraphErrors();

  for(int i=0;i<histbins;i++)
  {
    if(histexp->GetBinCenter(i+1) > 0.2) continue;
    ge_SysError_C2->SetPoint(i, histexp->GetBinCenter(i+1), histexp->GetBinContent(i+1));
    ge_SysError_C2->SetPointError(i, errorwidth, histerr->GetBinContent(i+1));
  }

  return ge_SysError_C2;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *FemtoModelFitBands(TGraph *grMedian1, TGraph *grLower, TGraph *grUpper)
{
  TGraphErrors *grFemtoModel = new TGraphErrors();
  grFemtoModel->SetName(grMedian1->GetName());
  double x, yM1, yLo, yUp;
  int count = 0;
  for(int i=0; i< grMedian1->GetN(); ++i) {
    grMedian1->GetPoint(i, x, yM1);
    grLower->GetPoint(i, x, yLo);
    grUpper->GetPoint(i, x, yUp);
    std::vector<float> yAll;
    yAll.push_back(yM1);
    yAll.push_back(yLo);
    yAll.push_back(yUp);
    std::sort(yAll.begin(), yAll.end());
    grFemtoModel->SetPoint(count, x/1000.f, (yAll[2]+yAll[0])/2.f);
    grFemtoModel->SetPointError(count++, 0, (yAll[2]+yAll[0])/2.f - yAll[0]);
  }
  return grFemtoModel;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraph *convertInGev(TGraph *gr) {
  TGraph *grOut = new TGraph();
  grOut->SetName(gr->GetName());
  double x,y;
  for(int i=0; i< gr->GetN(); ++i) {
    gr->GetPoint(i, x, y);
    grOut->SetPoint(i, x/1000.f, y);
  }
  return grOut;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *convertHistoInGev(TH1F *gr) {
  TGraphErrors *grOut = new TGraphErrors();
  grOut->SetName(gr->GetName());
  for(int i=0; i< gr->GetNbinsX(); ++i) {
    grOut->SetPoint(i, gr->GetBinCenter(i)/1000.f, gr->GetBinContent(i));
    grOut->SetPointError(i, 0, gr->GetBinError(i));
  }
  return grOut;
}
void PlotMCCF(const char *filename, const char *prefix) {
  SetStyle(false,true);
  TH1F *SEALL[6][6];
  TH1F *SEALL_CommonAncest[6][6];
  TH1F *SEALL_NonCommonAncest[6][6];
  TH1F *MEALL[6][6];
  TH1F *CFALL[6][6];
  TH1F *SEACCEPT[6][6];
  TH1F *SEACCEPT_CommonAncest[6][6];
  TH1F *SEACCEPT_NonCommonAncest[6][6];
  TH1F *MEACCEPT[6][6];
  TH1F *CFACCEPT[6][6];
  TFile *file=TFile::Open(filename);
  TString QAName="";
  QAName+=prefix;
  TDirectoryFile *dirQA=(TDirectoryFile*)(file->FindObjectAny(QAName.Data()));
  if (dirQA) {
    TList *tmp=dirQA->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *QA;
    dirQA->GetObject(name,QA);
    if (QA) {
      TList *CollAll=(TList*)QA->FindObject("CollAll");
      TList *CollAllResults;
      if (CollAll) {
        CollAllResults=(TList*)CollAll->FindObject("Results");
      }
      TList *CollAccept=(TList*)QA->FindObject("CollAccept");
      TList *CollAcceptResults;
      if (CollAccept) {
        CollAcceptResults=(TList*)CollAccept->FindObject("Results");
      }
      for (int iPart1=0;iPart1<6;iPart1++) {
        for (int iPart2=iPart1;iPart2<6;iPart2++) {
          TString ListName=Form("Particle%i_Particle%i",iPart1,iPart2);
          TList *ListALL=(TList*)CollAllResults->FindObject(ListName.Data());
          TString SEName=Form("SEDist_Particle%i_Particle%i",iPart1,iPart2);
          SEALL[iPart1][iPart2]=(TH1F*)ListALL->FindObject(SEName.Data());
          if (!SEALL[iPart1][iPart2]) {
            std::cout << SEName.Data() << " not Found\n";
            std::cout << SEName.Data() << " not Found\n";
          }
          TString NewSEALLName=Form("%s%s",SEName.Data(),"ALL");
          SEALL[iPart1][iPart2]->SetNameTitle(NewSEALLName.Data(),NewSEALLName.Data());

          TString SE_Common=Form("SECommonAncestDist_Particle%i_Particle%i",iPart1,iPart2);
          SEALL_CommonAncest[iPart1][iPart2]=(TH1F*)ListALL->FindObject(SE_Common.Data());
          if (!SEALL_CommonAncest[iPart1][iPart2]) {
            std::cout << SE_Common.Data() << " not Found\n";
            std::cout << SE_Common.Data() << " not Found\n";
          }
          TString NewSEALLCommonName=Form("%s%s",SE_Common.Data(),"ALL");
          SEALL_CommonAncest[iPart1][iPart2]->SetNameTitle(NewSEALLCommonName.Data(),NewSEALLCommonName.Data());

          TString SE_NonCommon=Form("SENonCommonAncestDist_Particle%i_Particle%i",iPart1,iPart2);
          SEALL_NonCommonAncest[iPart1][iPart2]=(TH1F*)ListALL->FindObject(SE_NonCommon.Data());
          if (!SEALL_NonCommonAncest[iPart1][iPart2]) {
            std::cout << SE_NonCommon.Data() << " not Found\n";
            std::cout << SE_NonCommon.Data() << " not Found\n";
          }
          TString NewSEALLNonCommonName=Form("%s%s",SE_NonCommon.Data(),"ALL");
          SEALL_NonCommonAncest[iPart1][iPart2]->SetNameTitle(NewSEALLNonCommonName.Data(),NewSEALLNonCommonName.Data());

          TString MEName=Form("MEDist_Particle%i_Particle%i",iPart1,iPart2);
          MEALL[iPart1][iPart2]=(TH1F*)ListALL->FindObject(MEName.Data());
          if (!MEALL[iPart1][iPart2]) {
            std::cout << SEName.Data() << " not Found\n";
            std::cout << SEName.Data() << " not Found\n";
          }
          TList *ListACCEPT=(TList*)CollAcceptResults->FindObject(ListName.Data());

          SEACCEPT[iPart1][iPart2]=(TH1F*)ListACCEPT->FindObject(SEName.Data());
          if (!SEACCEPT[iPart1][iPart2]) {
            std::cout << SEName.Data() << " not Found\n";
            std::cout << SEName.Data() << " not Found\n";
          }
          TString NewSEACCEPTName=Form("%s%s",SEName.Data(),"ACCEPT");
          SEACCEPT[iPart1][iPart2]->SetNameTitle(NewSEACCEPTName.Data(),NewSEACCEPTName.Data());

          SEACCEPT_CommonAncest[iPart1][iPart2]=(TH1F*)ListACCEPT->FindObject(SE_Common.Data());
          if (!SEACCEPT_CommonAncest[iPart1][iPart2]) {
            std::cout << SE_Common.Data() << " not Found\n";
            std::cout << SE_Common.Data() << " not Found\n";
          }
          TString NewSEACCEPTCommonName=Form("%s%s",SE_Common.Data(),"ACCEPT");
          SEACCEPT_CommonAncest[iPart1][iPart2]->SetNameTitle(NewSEACCEPTCommonName.Data(),NewSEACCEPTCommonName.Data());

          SEACCEPT_NonCommonAncest[iPart1][iPart2]=(TH1F*)ListACCEPT->FindObject(SE_NonCommon.Data());
          if (!SEACCEPT_NonCommonAncest[iPart1][iPart2]) {
            std::cout << SE_NonCommon.Data() << " not Found\n";
            std::cout << SE_NonCommon.Data() << " not Found\n";
          }
          TString NewSEACCEPTNonCommonName=Form("%s%s",SE_NonCommon.Data(),"ACCEPT");
          SEACCEPT_NonCommonAncest[iPart1][iPart2]->SetNameTitle(NewSEACCEPTNonCommonName.Data(),NewSEACCEPTNonCommonName.Data());

          MEACCEPT[iPart1][iPart2]=(TH1F*)ListACCEPT->FindObject(MEName.Data());
          if (!MEACCEPT[iPart1][iPart2]) {
            std::cout << SEName.Data() << " not Found\n";
            std::cout << SEName.Data() << " not Found\n";
          }
          TString NewMEACCEPTName=Form("%s%s",MEName.Data(),"ACCEPT");
          MEACCEPT[iPart1][iPart2]->SetNameTitle(NewMEACCEPTName.Data(),NewMEACCEPTName.Data());

        }
      }
    } else {
      std::cout << "No QA List \n";
    }
  } else {
    std::cout << "No Dir \n";
  }
  const float normleft = 0.2;
  const float normright = 0.4;
  TH1F *ALL_hist_CF_Lp_ALAp_exp[3];
  TH1F *ALLCA_hist_CF_Lp_ALAp_exp[3];
  TH1F *ALLNCA_hist_CF_Lp_ALAp_exp[3];

  TH1F *ALL_hist_CF_LL_ALAL_exp[3];
  TH1F *ALLCA_hist_CF_LL_ALAL_exp[3];
  TH1F *ALLNCA_hist_CF_LL_ALAL_exp[3];

  TH1F *ALL_hist_CF_pp_ApAp_exp[3];
  TH1F *ALLCA_hist_CF_pp_ApAp_exp[3];
  TH1F *ALLNCA_hist_CF_pp_ApAp_exp[3];

  TH1F *ALL_hist_CF_pXi_ApAXi_exp[3];
  TH1F *ALLCA_hist_CF_pXi_ApAXi_exp[3];
  TH1F *ALLNCA_hist_CF_pXi_ApAXi_exp[3];

  TH1F *ACCEPT_hist_CF_Lp_ALAp_exp[3];
  TH1F *ACCEPTCA_hist_CF_Lp_ALAp_exp[3];
  TH1F *ACCEPTNCA_hist_CF_Lp_ALAp_exp[3];

  TH1F *ACCEPT_hist_CF_LL_ALAL_exp[3];
  TH1F *ACCEPTCA_hist_CF_LL_ALAL_exp[3];
  TH1F *ACCEPTNCA_hist_CF_LL_ALAL_exp[3];

  TH1F *ACCEPT_hist_CF_pp_ApAp_exp[3];
  TH1F *ACCEPTCA_hist_CF_pp_ApAp_exp[3];
  TH1F *ACCEPTNCA_hist_CF_pp_ApAp_exp[3];

  TH1F *ACCEPT_hist_CF_pXi_ApAXi_exp[3];
  TH1F *ACCEPTCA_hist_CF_pXi_ApAXi_exp[3];
  TH1F *ACCEPTNCA_hist_CF_pXi_ApAXi_exp[3];

  auto *c1= new TCanvas("Can_CF","Can_CF",0,0,2000,1100);
  c1->Divide(2,2);

  ALL_hist_CF_Lp_ALAp_exp[0] = Calculate_CF(SEALL[0][2],MEALL[0][2],"CF_PL_ALL",normleft,normright);
  ALLCA_hist_CF_Lp_ALAp_exp[0] = Calculate_CF(SEALL_CommonAncest[0][2],MEALL[0][2],"CF_PL_ALLCommonAncest",normleft,normright);
  ALLNCA_hist_CF_Lp_ALAp_exp[0] = Calculate_CF(SEALL_NonCommonAncest[0][2],MEALL[0][2],"CF_PL_ALLNonCommonAncest",normleft,normright);

  ALL_hist_CF_Lp_ALAp_exp[1] = Calculate_CF(SEALL[1][3],MEALL[1][3],"CF_APAL_ALL",normleft,normright);
  ALLCA_hist_CF_Lp_ALAp_exp[1] = Calculate_CF(SEALL_CommonAncest[1][3],MEALL[1][3],"CF_APAL_ALLCommonAncest",normleft,normright);
  ALLNCA_hist_CF_Lp_ALAp_exp[1] = Calculate_CF(SEALL_NonCommonAncest[1][3],MEALL[1][3],"CF_APAL_ALLNonCommonAncest",normleft,normright);

  ALL_hist_CF_Lp_ALAp_exp[2] = add_CF(ALL_hist_CF_Lp_ALAp_exp[0],ALL_hist_CF_Lp_ALAp_exp[1],"CF_PLAPAL_ALL");
  ALLCA_hist_CF_Lp_ALAp_exp[2] = add_CF(ALLCA_hist_CF_Lp_ALAp_exp[0],ALLCA_hist_CF_Lp_ALAp_exp[1],"CF_PLAPAL_ALLCommonAncest");
  ALLNCA_hist_CF_Lp_ALAp_exp[2] = add_CF(ALLNCA_hist_CF_Lp_ALAp_exp[0],ALLNCA_hist_CF_Lp_ALAp_exp[1],"CF_PLAPAL_ALLNonCommonAncest");
  c1->cd(2);
  ALL_hist_CF_Lp_ALAp_exp[2]->SetLineColor(2);
  ALL_hist_CF_Lp_ALAp_exp[2]->SetTitle("Proton-Lambda ALL");
  ALL_hist_CF_Lp_ALAp_exp[2]->DrawCopy();
  ALLCA_hist_CF_Lp_ALAp_exp[2]->SetLineColor(3);
  ALLCA_hist_CF_Lp_ALAp_exp[2]->DrawCopy("SAME");
  ALLNCA_hist_CF_Lp_ALAp_exp[2]->SetLineColor(4);
  ALLNCA_hist_CF_Lp_ALAp_exp[2]->DrawCopy("SAME");

  ALL_hist_CF_LL_ALAL_exp[0] = Calculate_CF(SEALL[2][2],MEALL[2][2],"CF_LL_ALL",normleft,normright);
  ALLCA_hist_CF_LL_ALAL_exp[0] = Calculate_CF(SEALL_CommonAncest[2][2],MEALL[2][2],"CF_LL_ALL",normleft,normright);
  ALLNCA_hist_CF_LL_ALAL_exp[0] = Calculate_CF(SEALL_NonCommonAncest[2][2],MEALL[2][2],"CF_LL_ALL",normleft,normright);

  ALL_hist_CF_LL_ALAL_exp[1] = Calculate_CF(SEALL[3][3],MEALL[3][3],"CF_ALAL_ALL",normleft,normright);
  ALLCA_hist_CF_LL_ALAL_exp[1] = Calculate_CF(SEALL_CommonAncest[3][3],MEALL[3][3],"CF_ALAL_ALL",normleft,normright);
  ALLNCA_hist_CF_LL_ALAL_exp[1] = Calculate_CF(SEALL_NonCommonAncest[3][3],MEALL[3][3],"CF_ALAL_ALL",normleft,normright);

  ALL_hist_CF_LL_ALAL_exp[2] = add_CF(ALL_hist_CF_LL_ALAL_exp[0],ALL_hist_CF_LL_ALAL_exp[1],"CF_LLALAL_ALL");
  ALLCA_hist_CF_LL_ALAL_exp[2] = add_CF(ALLCA_hist_CF_LL_ALAL_exp[0],ALLCA_hist_CF_LL_ALAL_exp[1],"CF_LLALAL_ALL");
  ALLNCA_hist_CF_LL_ALAL_exp[2] = add_CF(ALLNCA_hist_CF_LL_ALAL_exp[0],ALLNCA_hist_CF_LL_ALAL_exp[1],"CF_LLALAL_ALL");
  c1->cd(3);
  ALL_hist_CF_LL_ALAL_exp[2]->SetLineColor(2);
  ALL_hist_CF_LL_ALAL_exp[2]->SetTitle("Lambda-Lambda ALL");
  ALL_hist_CF_LL_ALAL_exp[2]->DrawCopy();
  ALLCA_hist_CF_LL_ALAL_exp[2]->SetLineColor(3);
  ALLCA_hist_CF_LL_ALAL_exp[2]->DrawCopy("SAME");
  ALLNCA_hist_CF_LL_ALAL_exp[2]->SetLineColor(4);
  ALLNCA_hist_CF_LL_ALAL_exp[2]->DrawCopy("SAME");


  ALL_hist_CF_pp_ApAp_exp[0] = Calculate_CF(SEALL[0][0],MEALL[0][0],"CF_PP_ALL",normleft,normright);
  ALLCA_hist_CF_pp_ApAp_exp[0] = Calculate_CF(SEALL_CommonAncest[0][0],MEALL[0][0],"CF_PP_ALL",normleft,normright);
  ALLNCA_hist_CF_pp_ApAp_exp[0] = Calculate_CF(SEALL_NonCommonAncest[0][0],MEALL[0][0],"CF_PP_ALL",normleft,normright);

  ALL_hist_CF_pp_ApAp_exp[1] = Calculate_CF(SEALL[1][1],MEALL[1][1],"CF_APAP_ALL",normleft,normright);
  ALLCA_hist_CF_pp_ApAp_exp[1] = Calculate_CF(SEALL_CommonAncest[1][1],MEALL[1][1],"CF_APAP_ALL",normleft,normright);
  ALLNCA_hist_CF_pp_ApAp_exp[1] = Calculate_CF(SEALL_NonCommonAncest[1][1],MEALL[1][1],"CF_APAP_ALL",normleft,normright);

  ALL_hist_CF_pp_ApAp_exp[2] = add_CF(ALL_hist_CF_pp_ApAp_exp[0],ALL_hist_CF_pp_ApAp_exp[1],"CF_PPAPAP_ALL");
  ALLCA_hist_CF_pp_ApAp_exp[2] = add_CF(ALLCA_hist_CF_pp_ApAp_exp[0],ALLCA_hist_CF_pp_ApAp_exp[1],"CF_PPAPAP_ALL");
  ALLNCA_hist_CF_pp_ApAp_exp[2] = add_CF(ALLNCA_hist_CF_pp_ApAp_exp[0],ALLNCA_hist_CF_pp_ApAp_exp[1],"CF_PPAPAP_ALL");
  c1->cd(1);
  ALL_hist_CF_pp_ApAp_exp[2]->SetLineColor(2);
  ALL_hist_CF_pp_ApAp_exp[2]->SetTitle("Proton-Proton ALL");
  ALL_hist_CF_pp_ApAp_exp[2]->DrawCopy();
  ALLCA_hist_CF_pp_ApAp_exp[2]->SetLineColor(3);
  ALLCA_hist_CF_pp_ApAp_exp[2]->DrawCopy("SAME");
  ALLNCA_hist_CF_pp_ApAp_exp[2]->SetLineColor(4);
  ALLNCA_hist_CF_pp_ApAp_exp[2]->DrawCopy("SAME");

  ALL_hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(SEALL[0][4],MEALL[0][4],"CF_PXI_ALL",normleft,normright);
  ALLCA_hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(SEALL_CommonAncest[0][4],MEALL[0][4],"CF_PXI_ALL",normleft,normright);
  ALLNCA_hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(SEALL_NonCommonAncest[0][4],MEALL[0][4],"CF_PXI_ALL",normleft,normright);

  ALL_hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(SEALL[1][5],MEALL[1][5],"CF_APAXI_ALL",normleft,normright);
  ALLCA_hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(SEALL_CommonAncest[1][5],MEALL[1][5],"CF_APAXI_ALL",normleft,normright);
  ALLNCA_hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(SEALL_NonCommonAncest[1][5],MEALL[1][5],"CF_APAXI_ALL",normleft,normright);

  ALL_hist_CF_pXi_ApAXi_exp[2] = add_CF(ALL_hist_CF_pXi_ApAXi_exp[0],ALL_hist_CF_pXi_ApAXi_exp[1],"CF_PXI_APAXI_ALL");
  ALLCA_hist_CF_pXi_ApAXi_exp[2] = add_CF(ALLCA_hist_CF_pXi_ApAXi_exp[0],ALLCA_hist_CF_pXi_ApAXi_exp[1],"CF_PXI_APAXI_ALL");
  ALLNCA_hist_CF_pXi_ApAXi_exp[2] = add_CF(ALLNCA_hist_CF_pXi_ApAXi_exp[0],ALLNCA_hist_CF_pXi_ApAXi_exp[1],"CF_PXI_APAXI_ALL");
  c1->cd(4);
  ALL_hist_CF_pXi_ApAXi_exp[2]->SetLineColor(2);
  ALL_hist_CF_pXi_ApAXi_exp[2]->SetTitle("Proton-Xi ALL");
  ALL_hist_CF_pXi_ApAXi_exp[2]->DrawCopy();
  ALLCA_hist_CF_pXi_ApAXi_exp[2]->SetLineColor(3);
  ALLCA_hist_CF_pXi_ApAXi_exp[2]->DrawCopy("SAME");
  ALLNCA_hist_CF_pXi_ApAXi_exp[2]->SetLineColor(4);
  ALLNCA_hist_CF_pXi_ApAXi_exp[2]->DrawCopy("SAME");

  auto *c2= new TCanvas("Can_ACEEPT","Can_ACCEPT",0,0,2000,1100);
  c2->Divide(2,2);

  ACCEPT_hist_CF_Lp_ALAp_exp[0] = Calculate_CF(SEACCEPT[0][2],MEACCEPT[0][2],"CF_PL_ACCEPT",normleft,normright);
  ACCEPTCA_hist_CF_Lp_ALAp_exp[0] = Calculate_CF(SEACCEPT_CommonAncest[0][2],MEACCEPT[0][2],"CF_PL_ACCEPTCommonAncest",normleft,normright);
  ACCEPTNCA_hist_CF_Lp_ALAp_exp[0] = Calculate_CF(SEACCEPT_NonCommonAncest[0][2],MEACCEPT[0][2],"CF_PL_ACCEPTNonCommonAncest",normleft,normright);

  ACCEPT_hist_CF_Lp_ALAp_exp[1] = Calculate_CF(SEACCEPT[1][3],MEACCEPT[1][3],"CF_APAL_ACCEPT",normleft,normright);
  ACCEPTCA_hist_CF_Lp_ALAp_exp[1] = Calculate_CF(SEACCEPT_CommonAncest[1][3],MEACCEPT[1][3],"CF_APAL_ACCEPTCommonAncest",normleft,normright);
  ACCEPTNCA_hist_CF_Lp_ALAp_exp[1] = Calculate_CF(SEACCEPT_NonCommonAncest[1][3],MEACCEPT[1][3],"CF_APAL_ACCEPTNonCommonAncest",normleft,normright);

  ACCEPT_hist_CF_Lp_ALAp_exp[2] = add_CF(ACCEPT_hist_CF_Lp_ALAp_exp[0],ACCEPT_hist_CF_Lp_ALAp_exp[1],"CF_PLAPAL_ACCEPT");
  ACCEPTCA_hist_CF_Lp_ALAp_exp[2] = add_CF(ACCEPTCA_hist_CF_Lp_ALAp_exp[0],ACCEPTCA_hist_CF_Lp_ALAp_exp[1],"CF_PLAPAL_ACCEPTCommonAncest");
  ACCEPTNCA_hist_CF_Lp_ALAp_exp[2] = add_CF(ACCEPTNCA_hist_CF_Lp_ALAp_exp[0],ACCEPTNCA_hist_CF_Lp_ALAp_exp[1],"CF_PLAPAL_ACCEPTNonCommonAncest");
  c2->cd(2);
  ACCEPT_hist_CF_Lp_ALAp_exp[2]->SetLineColor(2);
  ACCEPT_hist_CF_Lp_ALAp_exp[2]->SetTitle("Proton-Lambda ACCEPT");
  ACCEPT_hist_CF_Lp_ALAp_exp[2]->DrawCopy();
  ACCEPTCA_hist_CF_Lp_ALAp_exp[2]->SetLineColor(3);
  ACCEPTCA_hist_CF_Lp_ALAp_exp[2]->DrawCopy("SAME");
  ACCEPTNCA_hist_CF_Lp_ALAp_exp[2]->SetLineColor(4);
  ACCEPTNCA_hist_CF_Lp_ALAp_exp[2]->DrawCopy("SAME");

  ACCEPT_hist_CF_LL_ALAL_exp[0] = Calculate_CF(SEACCEPT[2][2],MEACCEPT[2][2],"CF_LL_ACCEPT",normleft,normright);
  ACCEPTCA_hist_CF_LL_ALAL_exp[0] = Calculate_CF(SEACCEPT_CommonAncest[2][2],MEACCEPT[2][2],"CF_LL_ACCEPT",normleft,normright);
  ACCEPTNCA_hist_CF_LL_ALAL_exp[0] = Calculate_CF(SEACCEPT_NonCommonAncest[2][2],MEACCEPT[2][2],"CF_LL_ACCEPT",normleft,normright);

  ACCEPT_hist_CF_LL_ALAL_exp[1] = Calculate_CF(SEACCEPT[3][3],MEACCEPT[3][3],"CF_ALAL_ACCEPT",normleft,normright);
  ACCEPTCA_hist_CF_LL_ALAL_exp[1] = Calculate_CF(SEACCEPT_CommonAncest[3][3],MEACCEPT[3][3],"CF_ALAL_ACCEPT",normleft,normright);
  ACCEPTNCA_hist_CF_LL_ALAL_exp[1] = Calculate_CF(SEACCEPT_NonCommonAncest[3][3],MEACCEPT[3][3],"CF_ALAL_ACCEPT",normleft,normright);

  ACCEPT_hist_CF_LL_ALAL_exp[2] = add_CF(ACCEPT_hist_CF_LL_ALAL_exp[0],ACCEPT_hist_CF_LL_ALAL_exp[1],"CF_LLALAL_ACCEPT");
  ACCEPTCA_hist_CF_LL_ALAL_exp[2] = add_CF(ACCEPTCA_hist_CF_LL_ALAL_exp[0],ACCEPTCA_hist_CF_LL_ALAL_exp[1],"CF_LLALAL_ACCEPT");
  ACCEPTNCA_hist_CF_LL_ALAL_exp[2] = add_CF(ACCEPTNCA_hist_CF_LL_ALAL_exp[0],ACCEPTNCA_hist_CF_LL_ALAL_exp[1],"CF_LLALAL_ACCEPT");
  c2->cd(3);
  ACCEPT_hist_CF_LL_ALAL_exp[2]->SetLineColor(2);
  ACCEPT_hist_CF_LL_ALAL_exp[2]->SetTitle("Lambda-Lambda ACCEPT");
  ACCEPT_hist_CF_LL_ALAL_exp[2]->DrawCopy();
  ACCEPTCA_hist_CF_LL_ALAL_exp[2]->SetLineColor(3);
  ACCEPTCA_hist_CF_LL_ALAL_exp[2]->DrawCopy("SAME");
  ACCEPTNCA_hist_CF_LL_ALAL_exp[2]->SetLineColor(4);
  ACCEPTNCA_hist_CF_LL_ALAL_exp[2]->DrawCopy("SAME");


  ACCEPT_hist_CF_pp_ApAp_exp[0] = Calculate_CF(SEACCEPT[0][0],MEACCEPT[0][0],"CF_PP_ACCEPT",normleft,normright);
  ACCEPTCA_hist_CF_pp_ApAp_exp[0] = Calculate_CF(SEACCEPT_CommonAncest[0][0],MEACCEPT[0][0],"CF_PP_ACCEPT",normleft,normright);
  ACCEPTNCA_hist_CF_pp_ApAp_exp[0] = Calculate_CF(SEACCEPT_NonCommonAncest[0][0],MEACCEPT[0][0],"CF_PP_ACCEPT",normleft,normright);

  ACCEPT_hist_CF_pp_ApAp_exp[1] = Calculate_CF(SEACCEPT[1][1],MEACCEPT[1][1],"CF_APAP_ACCEPT",normleft,normright);
  ACCEPTCA_hist_CF_pp_ApAp_exp[1] = Calculate_CF(SEACCEPT_CommonAncest[1][1],MEACCEPT[1][1],"CF_APAP_ACCEPT",normleft,normright);
  ACCEPTNCA_hist_CF_pp_ApAp_exp[1] = Calculate_CF(SEACCEPT_NonCommonAncest[1][1],MEACCEPT[1][1],"CF_APAP_ACCEPT",normleft,normright);

  ACCEPT_hist_CF_pp_ApAp_exp[2] = add_CF(ACCEPT_hist_CF_pp_ApAp_exp[0],ACCEPT_hist_CF_pp_ApAp_exp[1],"CF_PPAPAP_ACCEPT");
  ACCEPTCA_hist_CF_pp_ApAp_exp[2] = add_CF(ACCEPTCA_hist_CF_pp_ApAp_exp[0],ACCEPTCA_hist_CF_pp_ApAp_exp[1],"CF_PPAPAP_ACCEPT");
  ACCEPTNCA_hist_CF_pp_ApAp_exp[2] = add_CF(ACCEPTNCA_hist_CF_pp_ApAp_exp[0],ACCEPTNCA_hist_CF_pp_ApAp_exp[1],"CF_PPAPAP_ACCEPT");
  c2->cd(1);
  ACCEPT_hist_CF_pp_ApAp_exp[2]->SetLineColor(2);
  ACCEPT_hist_CF_pp_ApAp_exp[2]->SetTitle("Proton-Proton ACCEPT");
  ACCEPT_hist_CF_pp_ApAp_exp[2]->DrawCopy();
  ACCEPTCA_hist_CF_pp_ApAp_exp[2]->SetLineColor(3);
  ACCEPTCA_hist_CF_pp_ApAp_exp[2]->DrawCopy("SAME");
  ACCEPTNCA_hist_CF_pp_ApAp_exp[2]->SetLineColor(4);
  ACCEPTNCA_hist_CF_pp_ApAp_exp[2]->DrawCopy("SAME");

  ACCEPT_hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(SEACCEPT[0][4],MEACCEPT[0][4],"CF_PXI_ACCEPT",normleft,normright);
  ACCEPTCA_hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(SEACCEPT_CommonAncest[0][4],MEACCEPT[0][4],"CF_PXI_ACCEPT",normleft,normright);
  ACCEPTNCA_hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(SEACCEPT_NonCommonAncest[0][4],MEACCEPT[0][4],"CF_PXI_ACCEPT",normleft,normright);

  ACCEPT_hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(SEACCEPT[1][5],MEACCEPT[1][5],"CF_APAXI_ACCEPT",normleft,normright);
  ACCEPTCA_hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(SEACCEPT_CommonAncest[1][5],MEACCEPT[1][5],"CF_APAXI_ACCEPT",normleft,normright);
  ACCEPTNCA_hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(SEACCEPT_NonCommonAncest[1][5],MEACCEPT[1][5],"CF_APAXI_ACCEPT",normleft,normright);

  ACCEPT_hist_CF_pXi_ApAXi_exp[2] = add_CF(ACCEPT_hist_CF_pXi_ApAXi_exp[0],ACCEPT_hist_CF_pXi_ApAXi_exp[1],"CF_PXI_APAXI_ACCEPT");
  ACCEPTCA_hist_CF_pXi_ApAXi_exp[2] = add_CF(ACCEPTCA_hist_CF_pXi_ApAXi_exp[0],ACCEPTCA_hist_CF_pXi_ApAXi_exp[1],"CF_PXI_APAXI_ACCEPT");
  ACCEPTNCA_hist_CF_pXi_ApAXi_exp[2] = add_CF(ACCEPTNCA_hist_CF_pXi_ApAXi_exp[0],ACCEPTNCA_hist_CF_pXi_ApAXi_exp[1],"CF_PXI_APAXI_ACCEPT");
  c2->cd(4);
  ACCEPT_hist_CF_pXi_ApAXi_exp[2]->SetLineColor(2);
  ACCEPT_hist_CF_pXi_ApAXi_exp[2]->SetTitle("Proton-Xi ACCEPT");
  ACCEPT_hist_CF_pXi_ApAXi_exp[2]->DrawCopy();
  ACCEPTCA_hist_CF_pXi_ApAXi_exp[2]->SetLineColor(3);
  ACCEPTCA_hist_CF_pXi_ApAXi_exp[2]->DrawCopy("SAME");
  ACCEPTNCA_hist_CF_pXi_ApAXi_exp[2]->SetLineColor(4);
  ACCEPTNCA_hist_CF_pXi_ApAXi_exp[2]->DrawCopy("SAME");

}
