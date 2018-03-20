
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


void CompOldNew(const char *fileNameNew,const char *fileNameOld,const char *CForSEorME, const char *prefix) {
  SetStyle(true,true);
  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;
  TFile *file=TFile::Open(fileNameNew);

  TString RelKNames[6]={"Proton","AntiProton","Lambda","AntiLambda","Xi","AntiXi"};
  TString ResultsName="";
  ResultsName+=prefix;
  ResultsName+="Results";

  TH1F *SEDist[6][6];
  TH1F *MEDist[6][6];
  TH1F *CFDist[6][6];

  TDirectoryFile *dirResults=(TDirectoryFile*)(file->FindObjectAny(ResultsName.Data()));
  if (dirResults) {
    TList *tmp=dirResults->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *Results;
    dirResults->GetObject(name,Results);
    if (!Results) {
      std::cout << "No Results \n";
    } else {
      for (int iPart1=0;iPart1<6;++iPart1) {
        for (int iPart2=iPart1;iPart2<6;++iPart2) {
          TString folderName=Form("Particle%i_Particle%i",iPart1,iPart2);
          TList* tmpFolder=(TList*)Results->FindObject(folderName.Data());
          if (!tmpFolder) {
            std::cout << folderName.Data() <<" not Found\n";
          } else {
            TString SEName=Form("SEDist_Particle%i_Particle%i",iPart1,iPart2);
            SEDist[iPart1][iPart2]=(TH1F*)tmpFolder->FindObject(SEName.Data());
            SEDist[iPart1][iPart2]->SetDirectory(0);
            if (!SEDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }

            TString MEName=Form("MEDist_Particle%i_Particle%i",iPart1,iPart2);
            MEDist[iPart1][iPart2]=(TH1F*)tmpFolder->FindObject(MEName.Data());
            MEDist[iPart1][iPart2]->SetDirectory(0);
            if (!MEDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            TString CFName = Form("CF_Part%i_Part%i",iPart1,iPart2);
            CFDist[iPart1][iPart2]=Calculate_CF(SEDist[iPart1][iPart2],MEDist[iPart1][iPart2],CFName.Data(),normleft,normright);
          }
        }
      }
    }
  }

  TH1F *CFSumPP=add_CF(CFDist[0][0],CFDist[1][1],"ProtonProtonCF");
  TH1F *CFSumPL=add_CF(CFDist[0][2],CFDist[1][3],"ProtonLambdaCF");
  TH1F *CFSumLL=add_CF(CFDist[2][2],CFDist[3][3],"LambdaLambdaCF");
  TH1F *CFSumPXi=add_CF(CFDist[0][4],CFDist[1][5],"ProtonXiCF");

  TFile *fileOld=TFile::Open(fileNameOld);
  TList *listTP=0;
  fileOld->GetObject("/PWGCF_PLFemto_0/TP_dir_0",listTP);
  if(!listTP) fileOld->GetObject("/PWGCF_PLFemto_0/TPdir_0",listTP);
  TH1F* histRE_relK_Lp = (TH1F*)listTP->FindObject("ProtonLambda_relK");
  if(!histRE_relK_Lp) histRE_relK_Lp = (TH1F*)listTP->FindObject("fProtonLambdaRelK");
  TH1F* histME_relK_Lp = (TH1F*)listTP->FindObject("ProtonLambda_relK_ME");
  if(!histME_relK_Lp) histME_relK_Lp = (TH1F*)listTP->FindObject("fProtonLambdaRelKME");
  TH1F* histRE_relK_ALAp = (TH1F*)listTP->FindObject("AntiProtonAntiLambda_relK");
  if(!histRE_relK_ALAp) histRE_relK_ALAp = (TH1F*)listTP->FindObject("fAntiProtonAntiLambdaRelK");
  TH1F* histME_relK_ALAp = (TH1F*)listTP->FindObject("AntiProtonAntiLambda_relK_ME");
  if(!histME_relK_ALAp) histME_relK_ALAp = (TH1F*)listTP->FindObject("fAntiProtonAntiLambdaRelKME");
  TH1F* histRE_relK_LL = (TH1F*)listTP->FindObject("LambdaLambda_relK");
  if(!histRE_relK_LL) histRE_relK_LL = (TH1F*)listTP->FindObject("fLambdaLambdaRelK");
  TH1F* histME_relK_LL = (TH1F*)listTP->FindObject("LambdaLambda_relK_ME");
  if(!histME_relK_LL) histME_relK_LL = (TH1F*)listTP->FindObject("fLambdaLambdaRelKME");
  TH1F* histRE_relK_ALAL = (TH1F*)listTP->FindObject("AntiLambdaAntiLambda_relK");
  if(!histRE_relK_ALAL) histRE_relK_ALAL = (TH1F*)listTP->FindObject("fAntiLambdaAntiLambdaRelK");
  TH1F* histME_relK_ALAL = (TH1F*)listTP->FindObject("AntiLambdaAntiLambda_relK_ME");
  if(!histME_relK_ALAL) histME_relK_ALAL = (TH1F*)listTP->FindObject("fAntiLambdaAntiLambdaRelKME");
  TH1F* histRE_relK_pp = (TH1F*)listTP->FindObject("ProtonProton_relK_FB");
  if(!histRE_relK_pp) histRE_relK_pp = (TH1F*)listTP->FindObject("fProtonProtonRelK");
  TH1F* histME_relK_pp = (TH1F*)listTP->FindObject("ProtonProton_relK_FB_ME");
  if(!histME_relK_pp) histME_relK_pp = (TH1F*)listTP->FindObject("fProtonProtonRelKME");
  TH1F* histRE_relK_ApAp = (TH1F*)listTP->FindObject("AntiProtonAntiProton_relK_FB");
  if(!histRE_relK_ApAp) histRE_relK_ApAp = (TH1F*)listTP->FindObject("fAntiProtonAntiProtonRelK");
  TH1F* histME_relK_ApAp = (TH1F*)listTP->FindObject("AntiProtonAntiProton_relK_FB_ME");
  if(!histME_relK_ApAp) histME_relK_ApAp = (TH1F*)listTP->FindObject("fAntiProtonAntiProtonRelKME");
  TH1F* histRE_relK_Xip  =    (TH1F*)listTP->FindObject("fProtonXiRelK");
  TH1F* histME_relK_Xip   =   (TH1F*)listTP->FindObject("fProtonXiRelKME");
  TH1F* histRE_relK_AXiAp    = (TH1F*)listTP->FindObject("fAntiProtonXiRelK");
  TH1F* histME_relK_AXiAp    = (TH1F*)listTP->FindObject("fAntiProtonXiRelKME");
  TH1F *hist_CF_Lp_ALAp_exp[3];
  TH1F *hist_CF_LL_ALAL_exp[3];
  TH1F *hist_CF_pp_ApAp_exp[3];
  TH1F *hist_CF_pXi_ApAXi_exp[3];

  hist_CF_Lp_ALAp_exp[0] = Calculate_CF(histRE_relK_Lp,histME_relK_Lp,"hist_CF_Lp_exp",normleft,normright);
  hist_CF_Lp_ALAp_exp[1] = Calculate_CF(histRE_relK_ALAp,histME_relK_ALAp,"hist_CF_ALAp_exp",normleft,normright);
  hist_CF_Lp_ALAp_exp[2] = add_CF(hist_CF_Lp_ALAp_exp[0],hist_CF_Lp_ALAp_exp[1],"hist_CF_Lp_ALAp_exp_sum");
  hist_CF_LL_ALAL_exp[0] = Calculate_CF(histRE_relK_LL,histME_relK_LL,"hist_CF_LL_exp",normleft,normright);
  hist_CF_LL_ALAL_exp[1] = Calculate_CF(histRE_relK_ALAL,histME_relK_ALAL,"hist_CF_LL_exp",normleft,normright);
  hist_CF_LL_ALAL_exp[2] = add_CF(hist_CF_LL_ALAL_exp[0],hist_CF_LL_ALAL_exp[1],"hist_CF_LL_ALAL_exp_sum");
  hist_CF_pp_ApAp_exp[0] = Calculate_CF(histRE_relK_pp,histME_relK_pp,"hist_CF_pp",normleft,normright);
  hist_CF_pp_ApAp_exp[1] = Calculate_CF(histRE_relK_ApAp,histME_relK_ApAp,"hist_CF_ApAp",normleft,normright);
  hist_CF_pp_ApAp_exp[2] = add_CF(hist_CF_pp_ApAp_exp[0],hist_CF_pp_ApAp_exp[1],"hist_CF_pp_ApAp_exp_sum");
  hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(histRE_relK_Xip,histME_relK_Xip,"hist_CF_pXi",normleft,normright);
  hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(histRE_relK_AXiAp,histME_relK_AXiAp,"hist_CF_ApAXi",normleft,normright);
  hist_CF_pXi_ApAXi_exp[2] = add_CF(hist_CF_pXi_ApAXi_exp[0],hist_CF_pXi_ApAXi_exp[1],"hist_CF_pXi_ApAXi_exp_sum");

  TCanvas *Can_CF = new TCanvas("Can_CF","Can_CF",0,0,2000,2000);
  TString modus=Form("%s",CForSEorME);
  if (modus=="CF") {
    Can_CF->Divide(4,3);
    Can_CF->cd(1);
    CFDist[0][0]->SetTitle("Proton-Proton CF");
    CFDist[0][0]->DrawCopy();
    hist_CF_pp_ApAp_exp[0]->SetLineColor(2);
    hist_CF_pp_ApAp_exp[0]->DrawCopy("SAME");

    Can_CF->cd(2);
    CFDist[0][2]->SetTitle("Proton-Lambda CF");
    CFDist[0][2]->DrawCopy();
    hist_CF_Lp_ALAp_exp[0]->SetLineColor(2);
    hist_CF_Lp_ALAp_exp[0]->DrawCopy("SAME");

    Can_CF->cd(3);
    CFDist[2][2]->SetTitle("Lambda-Lambda CF");
    CFDist[2][2]->DrawCopy();
    hist_CF_LL_ALAL_exp[0]->SetLineColor(2);
    hist_CF_LL_ALAL_exp[0]->DrawCopy("SAME");

    Can_CF->cd(4);
    CFDist[0][4]->SetTitle("Proton-Xi CF");
    CFDist[0][4]->DrawCopy();
    hist_CF_pXi_ApAXi_exp[0]->SetLineColor(2);
    hist_CF_pXi_ApAXi_exp[0]->DrawCopy("SAME");

    Can_CF->cd(5);
    CFDist[1][1]->SetTitle("AntiProton-AntiProton CF");
    CFDist[1][1]->DrawCopy();
    hist_CF_pp_ApAp_exp[1]->SetLineColor(2);
    hist_CF_pp_ApAp_exp[1]->DrawCopy("SAME");

    Can_CF->cd(6);
    CFDist[1][3]->SetTitle("AntiProton-AntiLambda CF");
    CFDist[1][3]->DrawCopy();
    hist_CF_Lp_ALAp_exp[1]->SetLineColor(2);
    hist_CF_Lp_ALAp_exp[1]->DrawCopy("SAME");

    Can_CF->cd(7);
    CFDist[3][3]->SetTitle("AntiLambda-AntiLambda CF");
    CFDist[3][3]->DrawCopy();
    hist_CF_LL_ALAL_exp[1]->SetLineColor(2);
    hist_CF_LL_ALAL_exp[1]->DrawCopy("SAME");

    Can_CF->cd(8);
    CFDist[1][5]->SetTitle("AntiProton-AntiXi CF");
    CFDist[1][5]->DrawCopy();
    hist_CF_pXi_ApAXi_exp[1]->SetLineColor(2);
    hist_CF_pXi_ApAXi_exp[1]->DrawCopy("SAME");

    Can_CF->cd(9);
    CFSumPP->SetTitle("Proton-Proton Sum CF");
    CFSumPP->DrawCopy();
    hist_CF_pp_ApAp_exp[2]->SetLineColor(2);
    hist_CF_pp_ApAp_exp[2]->DrawCopy("SAME");

    Can_CF->cd(10);
    CFSumPL->SetTitle("Proton-Lambda Sum CF");
    CFSumPL->DrawCopy();
    hist_CF_Lp_ALAp_exp[2]->SetLineColor(2);
    hist_CF_Lp_ALAp_exp[2]->DrawCopy("SAME");

    Can_CF->cd(11);
    CFSumLL->SetTitle("Lambda-Lambda Sum CF");
    CFSumLL->DrawCopy();
    hist_CF_LL_ALAL_exp[2]->SetLineColor(2);
    hist_CF_LL_ALAL_exp[2]->DrawCopy("SAME");

    Can_CF->cd(12);
    CFSumPXi->SetTitle("Proton-Xi Sum CF");
    CFSumPXi->DrawCopy();
    hist_CF_pXi_ApAXi_exp[2]->SetLineColor(2);
    hist_CF_pXi_ApAXi_exp[2]->DrawCopy("SAME");

  } else if (modus=="SE") {
    Can_CF->Divide(4,4);
    Can_CF->cd(1);
    SEDist[0][0]->SetTitle("Proton-Proton SE");
    SEDist[0][0]->DrawCopy();
    histRE_relK_pp->SetLineColor(2);
    histRE_relK_pp->DrawCopy("SAME");

    Can_CF->cd(2);
    SEDist[0][2]->SetTitle("Proton-Lambda SE");
    SEDist[0][2]->DrawCopy();
    histRE_relK_Lp->SetLineColor(2);
    histRE_relK_Lp->DrawCopy("SAME");

    Can_CF->cd(3);
    SEDist[2][2]->SetTitle("Lambda-Lambda SE");
    SEDist[2][2]->DrawCopy();
    histRE_relK_LL->SetLineColor(2);
    histRE_relK_LL->DrawCopy("SAME");

    Can_CF->cd(4);
    SEDist[0][4]->SetTitle("Proton-Xi SE");
    SEDist[0][4]->GetYaxis()->SetRangeUser(0,7000);
    SEDist[0][4]->DrawCopy();
    histRE_relK_Xip->SetLineColor(2);
    histRE_relK_Xip->DrawCopy("SAME");

    Can_CF->cd(5);
    TH1F* RatioSEPP = (TH1F*)SEDist[0][0]->Clone("RatioSEPP");
    RatioSEPP->Divide(histRE_relK_pp);
    RatioSEPP->SetTitle("Ratio PP SE");
    RatioSEPP->DrawCopy();

    Can_CF->cd(6);
    TH1F* RatioSEPL = (TH1F*)SEDist[0][2]->Clone("RatioSEPL");
    RatioSEPL->Divide(histRE_relK_Lp);
    RatioSEPL->SetTitle("Ratio PL SE");
    RatioSEPL->DrawCopy();

    Can_CF->cd(7);
    TH1F* RatioSELL = (TH1F*)SEDist[2][2]->Clone("RatioSELL");
    RatioSELL->Divide(histRE_relK_LL);
    RatioSELL->SetTitle("Ratio LL SE");
    RatioSELL->DrawCopy();

    Can_CF->cd(8);
    TH1F* RatioSEPXi = (TH1F*)SEDist[0][4]->Clone("RatioSEPXi");
    RatioSEPXi->Divide(histRE_relK_Xip);
    RatioSEPXi->SetTitle("Ratio PXi SE");
    RatioSEPXi->DrawCopy();

    Can_CF->cd(9);
    SEDist[1][1]->SetTitle("AntiProton-AntiProton SE");
    SEDist[1][1]->DrawCopy();
    histRE_relK_ApAp->SetLineColor(2);
    histRE_relK_ApAp->DrawCopy("SAME");

    Can_CF->cd(10);
    SEDist[1][3]->SetTitle("AntiProton-AntiLambda SE");
    SEDist[1][3]->DrawCopy();
    histRE_relK_ALAp->SetLineColor(2);
    histRE_relK_ALAp->DrawCopy("SAME");

    Can_CF->cd(11);
    SEDist[3][3]->SetTitle("AntiLambda-AntiLambda SE");
    SEDist[3][3]->DrawCopy();
    histRE_relK_ALAL->SetLineColor(2);
    histRE_relK_ALAL->DrawCopy("SAME");

    Can_CF->cd(12);
    SEDist[1][5]->SetTitle("AntiProton-AntiXi SE");
    SEDist[1][5]->GetYaxis()->SetRangeUser(0,5000);
    SEDist[1][5]->DrawCopy();
    histRE_relK_AXiAp->SetLineColor(2);
    histRE_relK_AXiAp->DrawCopy("SAME");

    Can_CF->cd(13);
    TH1F* RatioSEAPAP = (TH1F*)SEDist[1][1]->Clone("RatioSEAPAP");
    RatioSEAPAP->Divide(histRE_relK_ApAp);
    RatioSEAPAP->SetTitle("Ratio APAP SE");
    RatioSEAPAP->DrawCopy();

    Can_CF->cd(14);
    TH1F* RatioSEAPAL = (TH1F*)SEDist[1][3]->Clone("RatioSEAPAL");
    RatioSEAPAL->Divide(histRE_relK_ALAp);
    RatioSEAPAL->SetTitle("Ratio APAL SE");
    RatioSEAPAL->DrawCopy();

    Can_CF->cd(15);
    TH1F* RatioSEALAL = (TH1F*)SEDist[3][3]->Clone("RatioSEALAL");
    RatioSEALAL->Divide(histRE_relK_ALAL);
    RatioSEALAL->SetTitle("Ratio ALAL SE");
    RatioSEALAL->DrawCopy();

    Can_CF->cd(16);
    TH1F* RatioSEAPAXi = (TH1F*)SEDist[1][5]->Clone("RatioSEAPAXi");
    RatioSEAPAXi->Divide(histRE_relK_AXiAp);
    RatioSEAPAXi->SetTitle("Ratio APAXi SE");
    RatioSEAPAXi->DrawCopy();

  } else if (modus=="ME") {
    Can_CF->Divide(4,4);
    Can_CF->cd(1);
    MEDist[0][0]->SetTitle("Proton-Proton ME");
    MEDist[0][0]->Scale(1/(MEDist[0][0]->Integral(MEDist[0][0]->FindBin(normleft),MEDist[0][0]->FindBin(normright))));
    MEDist[0][0]->DrawCopy();
    histME_relK_pp->Scale(1/(histME_relK_pp->Integral(histME_relK_pp->FindBin(normleft),histME_relK_pp->FindBin(normright))));
    histME_relK_pp->SetLineColor(2);
    histME_relK_pp->DrawCopy("SAME");

    Can_CF->cd(2);
    MEDist[0][2]->SetTitle("Proton-Lambdas ME");
    MEDist[0][2]->Scale(1/(MEDist[0][2]->Integral(MEDist[0][2]->FindBin(normleft),MEDist[0][2]->FindBin(normright))));
    MEDist[0][2]->DrawCopy();
    histME_relK_Lp->Scale(1/(histME_relK_Lp->Integral(histME_relK_Lp->FindBin(normleft),histME_relK_Lp->FindBin(normright))));
    histME_relK_Lp->SetLineColor(2);
    histME_relK_Lp->DrawCopy("SAME");

    Can_CF->cd(3);
    MEDist[2][2]->SetTitle("Lambda-Lambda ME");
    MEDist[2][2]->Scale(1/(MEDist[2][2]->Integral(MEDist[2][2]->FindBin(normleft),MEDist[2][2]->FindBin(normright))));
    MEDist[2][2]->DrawCopy();
    histME_relK_LL->Scale(1/(histME_relK_LL->Integral(histME_relK_LL->FindBin(normleft),histME_relK_LL->FindBin(normright))));
    histME_relK_LL->SetLineColor(2);
    histME_relK_LL->DrawCopy("SAME");

    Can_CF->cd(4);
    MEDist[0][4]->SetTitle("Proton-Xi ME");
    MEDist[0][4]->Scale(1/(MEDist[0][4]->Integral(MEDist[0][4]->FindBin(normleft),MEDist[0][4]->FindBin(normright))));
    MEDist[0][4]->DrawCopy();
    histME_relK_Xip->Scale(1/(histME_relK_Xip->Integral(histME_relK_Xip->FindBin(normleft),histME_relK_Xip->FindBin(normright))));
    histME_relK_Xip->SetLineColor(2);
    histME_relK_Xip->DrawCopy("SAME");

    Can_CF->cd(5);
    TH1F* RatioMEPP = (TH1F*)MEDist[0][0]->Clone("RatioMEPP");
    RatioMEPP->Divide(histME_relK_pp);
    RatioMEPP->SetTitle("Ratio PP ME");
    RatioMEPP->DrawCopy();

    Can_CF->cd(6);
    TH1F* RatioMEPL = (TH1F*)MEDist[0][2]->Clone("RatioMEPL");
    RatioMEPL->Divide(histME_relK_Lp);
    RatioMEPL->SetTitle("Ratio PL ME");
    RatioMEPL->DrawCopy();

    Can_CF->cd(7);
    TH1F* RatioMELL = (TH1F*)MEDist[2][2]->Clone("RatioMELL");
    RatioMELL->Divide(histME_relK_LL);
    RatioMELL->SetTitle("Ratio LL ME");
    RatioMELL->DrawCopy();

    Can_CF->cd(8);
    TH1F* RatioMEPXi = (TH1F*)MEDist[0][4]->Clone("RatioMEPXi");
    RatioMEPXi->Divide(histME_relK_Xip);
    RatioMEPXi->SetTitle("Ratio PXi ME");
    RatioMEPXi->DrawCopy();

    Can_CF->cd(9);
    MEDist[1][1]->SetTitle("AntiProton-AntiProton ME");
    MEDist[1][1]->Scale(1/(MEDist[1][1]->Integral(MEDist[1][1]->FindBin(normleft),MEDist[1][1]->FindBin(normright))));
    MEDist[1][1]->DrawCopy();
    histME_relK_ApAp->Scale(1/(histME_relK_ApAp->Integral(histME_relK_ApAp->FindBin(normleft),histME_relK_ApAp->FindBin(normright))));
    histME_relK_ApAp->SetLineColor(2);
    histME_relK_ApAp->DrawCopy("SAME");

    Can_CF->cd(10);
    MEDist[1][3]->SetTitle("AntiProton-AntiLambda ME");
    MEDist[1][3]->Scale(1/(MEDist[1][3]->Integral(MEDist[1][3]->FindBin(normleft),MEDist[1][3]->FindBin(normright))));
    MEDist[1][3]->DrawCopy();
    histME_relK_ALAp->Scale(1/(histME_relK_ALAp->Integral(histME_relK_ALAp->FindBin(normleft),histME_relK_ALAp->FindBin(normright))));
    histME_relK_ALAp->SetLineColor(2);
    histME_relK_ALAp->DrawCopy("SAME");

    Can_CF->cd(11);
    MEDist[3][3]->SetTitle("AntiLambda-AntiLambda ME");
    MEDist[3][3]->Scale(1/(MEDist[3][3]->Integral(MEDist[3][3]->FindBin(normleft),MEDist[3][3]->FindBin(normright))));
    MEDist[3][3]->DrawCopy();
    histME_relK_ALAL->Scale(1/(histME_relK_ALAL->Integral(histME_relK_ALAL->FindBin(normleft),histME_relK_ALAL->FindBin(normright))));
    histME_relK_ALAL->SetLineColor(2);
    histME_relK_ALAL->DrawCopy("SAME");

    Can_CF->cd(12);
    MEDist[1][5]->SetTitle("AntiProton-AntiXi ME");
    MEDist[1][5]->Scale(1/(MEDist[1][5]->Integral(MEDist[1][5]->FindBin(normleft),MEDist[1][5]->FindBin(normright))));
    MEDist[1][5]->DrawCopy();
    histME_relK_AXiAp->Scale(1/(histME_relK_AXiAp->Integral(histME_relK_AXiAp->FindBin(normleft),histME_relK_AXiAp->FindBin(normright))));
    histME_relK_AXiAp->SetLineColor(2);
    histME_relK_AXiAp->DrawCopy("SAME");

    Can_CF->cd(13);
    TH1F* RatioMEAPAP = (TH1F*)MEDist[1][1]->Clone("RatioMEAPAP");
    RatioMEAPAP->Divide(histME_relK_ApAp);
    RatioMEAPAP->SetTitle("Ratio APAP ME");
    RatioMEAPAP->DrawCopy();

    Can_CF->cd(14);
    TH1F* RatioMEAPAL = (TH1F*)MEDist[1][3]->Clone("RatioMEAPAL");
    RatioMEAPAL->Divide(histME_relK_ALAp);
    RatioMEAPAL->SetTitle("Ratio APAL ME");
    RatioMEAPAL->DrawCopy();

    Can_CF->cd(15);
    TH1F* RatioMEALAL = (TH1F*)MEDist[3][3]->Clone("RatioMEALAL");
    RatioMEALAL->Divide(histME_relK_ALAL);
    RatioMEALAL->SetTitle("Ratio ALAL ME");
    RatioMEALAL->DrawCopy();

    Can_CF->cd(16);
    TH1F* RatioMEAPAXi = (TH1F*)MEDist[1][5]->Clone("RatioMEAPAXi");
    RatioMEAPAXi->Divide(histME_relK_AXiAp);
    RatioMEAPAXi->SetTitle("Ratio APAXi ME");
    RatioMEAPAXi->DrawCopy();

  }

}
