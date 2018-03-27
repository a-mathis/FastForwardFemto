
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
TH1F* Calculate_CF(TH1F* histRE_relK,TH1F* histME_relK, TString CFname,Double_t normleft,Double_t normright, const char* folder, float spinningDepth)
{
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  if(strcmp(folder, "") == 0) {
    Double_t norm_relK = histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) / histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright));
    Hist_CF->Divide(histRE_relK,histME_relK,1,norm_relK);
  }
  else {
    histME_relK->Scale(1.f/spinningDepth);
    Hist_CF->Divide(histRE_relK,histME_relK);
  }

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
void pPb_ppCFinMult(const char* input,const char *prefix,const char *suffix,bool doCent) {
  TList *TPdir=new TList();
  TPdir->SetName("TPdir_0");
  TList *TPdir2=new TList();
  TPdir2->SetName("NPairs");
  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;
  const float kTSplit[6][3]={{0,1.,4.5},  {0,1.,4.5}, {0,1.,4.5},  {0,1.,4.5}, {0,1.,4.5}, {0,1.,4.5}};
  const float mTSplit[21][3]=
  {
      {0,1.15,4.5},//pp
      {0,1.15,4.5},//pbar p
      {0,1.2,4.5},//p lambda
      {0,1.2,4.5},//p bar lambda
      {0,1.35,4.5},//p Xi
      {0,1.35,4.5},//p bar Xi
      {0,1.15,4.5},//bar p bar p
      {0,1.2,4.5},//bar p lambda
      {0,1.2,4.5},//bar p bar lambda
      {0,1.35,4.5},//bar p Xi
      {0,1.35,4.5},//bar p bar Xi
      {0,1.3,4.5},//lambda lambda
      {0,1.3,4.5},//lambda bar lambda
      {0,1.4,4.5},//lambda Xi
      {0,1.4,4.5},//lambda bar Xi
      {0,1.4,4.5},//bar Lambda bar lambda
      {0,1.4,4.5},//bar Lambda Xi
      {0,1.4,4.5},//bar Lambda bar Xi
      {0,1.4,4.5},//Xi Xi
      {0,1.4,4.5},//Xi bar Xi
      {0,1.4,4.5},//bar Xi bar Xi
  };

  gStyle->SetCanvasPreferGL(1);
  TFile *file=TFile::Open(input);
  TH2F* SEMultDist[6][6];
  TH2F* MEMultDist[6][6];
  TH1F* SEDist[6][6][3];
  TH1F* MEDist[6][6][3];

  TH2F* SEkTDist[6][6];
  TH2F* MEkTDist[6][6];
  TH1F* SEkTDistBin[6][6][2];
  TH1F* MEkTDistBin[6][6][2];

  TH2F* SEkTCentDist[6][6][3];
  TH2F* MEkTCentDist[6][6][3];
  TH1F* SEkTCentDistBin[6][6][3][2];
  TH1F* MEkTCentDistBin[6][6][3][2];

  TH2F* SEmTDist[6][6];
  TH2F* MEmTDist[6][6];
  TH1F* SEmTDistBin[6][6][2];
  TH1F* MEmTDistBin[6][6][2];

  TString ResultsName="";
  ResultsName+=prefix;
  ResultsName+="Results";
  ResultsName+=suffix;
  TString PartName[6]={"Proton","AntiProton","Lambda","AntiLambda","Xi","AntiXi"};
  int Centralities[4]={0,20,40,90};
  int MultBins[4]={1,8,12,21};
  int iPartCounter=0;
  TDirectoryFile *dirResults=(TDirectoryFile*)(file->FindObjectAny(ResultsName.Data()));
  if (dirResults) {
    TList *tmp=dirResults->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *Results;
    dirResults->GetObject(name,Results);
    if (!Results) {
      std::cout << "No Results \n";
    } else {
      for (int iPart1=0;iPart1<6;++iPart1){
        for (int iPart2=iPart1;iPart2<6;++iPart2){
          TString folderName=Form("Particle%i_Particle%i",iPart1,iPart2);
          TList* tmpFolder=(TList*)Results->FindObject(folderName.Data());
          if (!tmpFolder) {
            std::cout << folderName.Data() <<" not Found\n";
          } else {
            TString SEName=Form("SEMultDist_Particle%i_Particle%i",iPart1,iPart2);
            SEMultDist[iPart1][iPart2]=(TH2F*)tmpFolder->FindObject(SEName.Data());
            SEMultDist[iPart1][iPart2]->SetDirectory(0);
            if (!SEMultDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            TString MEName=Form("MEMultDist_Particle%i_Particle%i",iPart1,iPart2);
            MEMultDist[iPart1][iPart2]=(TH2F*)tmpFolder->FindObject(MEName.Data());
            MEMultDist[iPart1][iPart2]->SetDirectory(0);
            if (!MEMultDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            for (int iMult=0;iMult<3;++iMult) {
              TString projSEName=Form("SEDist%s%s_Mult_%i",PartName[iPart1].Data(),PartName[iPart2].Data(),iMult);
              SEDist[iPart1][iPart2][iMult]=(TH1F*)SEMultDist[iPart1][iPart2]->ProjectionX(projSEName.Data(),MultBins[iMult],MultBins[iMult+1]-1);

              TString projMEName=Form("MEDist%s%s_Mult_%i",PartName[iPart1].Data(),PartName[iPart2].Data(),iMult);
              MEDist[iPart1][iPart2][iMult]=(TH1F*)MEMultDist[iPart1][iPart2]->ProjectionX(projMEName.Data(),MultBins[iMult],MultBins[iMult+1]-1);
            }

            TString SEkTName=Form("SEkTDist_Particle%i_Particle%i",iPart1,iPart2);
            SEkTDist[iPart1][iPart2]=(TH2F*)tmpFolder->FindObject(SEkTName.Data());
            SEkTDist[iPart1][iPart2]->SetDirectory(0);
            if (!SEkTDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            TString MEkTName=Form("MEkTDist_Particle%i_Particle%i",iPart1,iPart2);
            MEkTDist[iPart1][iPart2]=(TH2F*)tmpFolder->FindObject(MEkTName.Data());
            MEkTDist[iPart1][iPart2]->SetDirectory(0);
            if (!MEkTDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            if (doCent) {
              for (int iCent=0;iCent<3;++iCent) {
                TString SEkTCentName=Form("SEkTCentDist_Cent_%i_Particle%i_Particle%i",Centralities[iCent+1],iPart1,iPart2);
                SEkTCentDist[iPart1][iPart2][iCent]=(TH2F*)tmpFolder->FindObject(SEkTCentName.Data());
                SEkTCentDist[iPart1][iPart2][iCent]->SetDirectory(0);
                if (!SEkTCentDist[iPart1][iPart2][iCent]) {
                  std::cout << SEName.Data() << " not Found\n";
                  std::cout << SEName.Data() << " not Found\n";
                }
                TString MEkTCentName=Form("MEkTCentDistCent%i_Particle%i_Particle%i",Centralities[iCent+1],iPart1,iPart2);
                MEkTCentDist[iPart1][iPart2][iCent]=(TH2F*)tmpFolder->FindObject(MEkTCentName.Data());
                MEkTCentDist[iPart1][iPart2][iCent]->SetDirectory(0);
                if (!MEkTCentDist[iPart1][iPart2][iCent]) {
                  std::cout << SEName.Data() << " not Found\n";
                  std::cout << SEName.Data() << " not Found\n";
                }
                for (int iBins=0;iBins<2;++iBins) {
                  TString projkTCentSEName=Form("SEkTCentDist%s%s_kTCent_%i_bin_%i",PartName[iPart1].Data(),PartName[iPart2].Data(),iCent,iBins);
                  SEkTCentDistBin[iPart1][iPart2][iCent][iBins]=(TH1F*)SEkTCentDist[iPart1][iPart2][iCent]->ProjectionX(
                      projkTCentSEName.Data(),SEkTCentDist[iPart1][iPart2][iCent]->GetXaxis()->FindBin(kTSplit[iPart2][iBins]),
                      SEkTCentDist[iPart1][iPart2][iCent]->FindBin(kTSplit[iPart2][iBins+1]));

                  TString projkTCentMEName=Form("MEkTCentDist%s%s_kTCent_%i_bin_%i",PartName[iPart1].Data(),PartName[iPart2].Data(),iCent,iBins);
                  MEkTCentDistBin[iPart1][iPart2][iCent][iBins]=(TH1F*)MEkTCentDist[iPart1][iPart2][iCent]->ProjectionX(
                      projkTCentMEName.Data(),MEkTCentDist[iPart1][iPart2][iCent]->GetXaxis()->FindBin(kTSplit[iPart2][iBins]),
                      MEkTCentDist[iPart1][iPart2][iCent]->GetXaxis()->FindBin(kTSplit[iPart2][iBins+1]));
                }
              }
            }
            TString SEmTName=Form("SEmTDist_Particle%i_Particle%i",iPart1,iPart2);
            SEmTDist[iPart1][iPart2]=(TH2F*)tmpFolder->FindObject(SEmTName.Data());
            SEmTDist[iPart1][iPart2]->SetDirectory(0);
            if (!SEmTDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            TString MEmTName=Form("MEmTDist_Particle%i_Particle%i",iPart1,iPart2);
            MEmTDist[iPart1][iPart2]=(TH2F*)tmpFolder->FindObject(MEmTName.Data());
            MEmTDist[iPart1][iPart2]->SetDirectory(0);
            if (!MEmTDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            for (int iBins=0;iBins<2;++iBins) {
              TString projkTSEName=Form("SEkTDist%s%s_kTbin_%i",PartName[iPart1].Data(),PartName[iPart2].Data(),iBins);
              SEkTDistBin[iPart1][iPart2][iBins]=(TH1F*)SEkTDist[iPart1][iPart2]->ProjectionX(
                  projkTSEName.Data(),SEkTDist[iPart1][iPart2]->GetXaxis()->FindBin(kTSplit[iPart2][iBins]),
                  SEkTDist[iPart1][iPart2]->GetXaxis()->FindBin(kTSplit[iPart2][iBins+1]));

              TString projkTMEName=Form("MEkTDist%s%s_kTbin_%i",PartName[iPart1].Data(),PartName[iPart2].Data(),iBins);
              MEkTDistBin[iPart1][iPart2][iBins]=(TH1F*)MEkTDist[iPart1][iPart2]->ProjectionX(
                  projkTMEName.Data(),MEkTDist[iPart1][iPart2]->GetXaxis()->FindBin(kTSplit[iPart2][iBins]),
                  MEkTDist[iPart1][iPart2]->GetXaxis()->FindBin(kTSplit[iPart2][iBins+1]));
              TString projmTSEName=Form("SEmTDist%s%s_mTbin_%i",PartName[iPart1].Data(),PartName[iPart2].Data(),iBins);
              SEmTDistBin[iPart1][iPart2][iBins]=(TH1F*)SEmTDist[iPart1][iPart2]->ProjectionX(
                  projmTSEName.Data(),SEmTDist[iPart1][iPart2]->GetXaxis()->FindBin(mTSplit[iPartCounter][iBins]),
                  SEmTDist[iPart1][iPart2]->GetXaxis()->FindBin(mTSplit[iPartCounter][iBins+1]));

              TString projmTMEName=Form("MEmTDist%s%s_mTbin_%i",PartName[iPart1].Data(),PartName[iPart2].Data(),iBins);
              MEmTDistBin[iPart1][iPart2][iBins]=(TH1F*)MEmTDist[iPart1][iPart2]->ProjectionX(
                  projmTMEName.Data(),MEmTDist[iPart1][iPart2]->GetXaxis()->FindBin(mTSplit[iPartCounter][iBins]),
                  MEmTDist[iPart1][iPart2]->GetXaxis()->FindBin(mTSplit[iPartCounter][iBins+1]));
            }
          }
          iPartCounter++;
        }
      }
    }
  } else {
    std::cout << "No directory!\n";
  }
  TH1F *hist_CF_pp_ApAp_exp[3][3];
  TH1F *hist_CF_pL_ApAL_exp[3][3];
  TH1F *hist_CF_LL_ALAL_exp[3][3];
  TH1F *hist_CF_pXi_ApAXi_exp[3][3];
  TH1F *hist_Pairs[4][3];
  hist_Pairs[0][0]=new TH1F("ppApApPairs02","ppApApPairs02",3,0,3);
  hist_Pairs[0][1]=new TH1F("ppApApPairs04","ppApApPairs02",3,0,3);
  hist_Pairs[0][2]=new TH1F("ppApApPairsInt","ppApApPairsInt",3,0,3);

  hist_Pairs[1][0]=new TH1F("pLApALPairs02","ppApApPairs02",3,0,3);
  hist_Pairs[1][1]=new TH1F("pLApALPairs04","ppApApPairs02",3,0,3);
  hist_Pairs[1][2]=new TH1F("pLApALPairsInt","ppApApPairsInt",3,0,3);

  hist_Pairs[2][0]=new TH1F("LLALALPairs02","ppApApPairs02",3,0,3);
  hist_Pairs[2][1]=new TH1F("LLALALPairs04","ppApApPairs02",3,0,3);
  hist_Pairs[2][2]=new TH1F("LLALALPairsInt","ppApApPairsInt",3,0,3);

  hist_Pairs[3][0]=new TH1F("pXiApAXiPairs02","ppApApPairs02",3,0,3);
  hist_Pairs[3][1]=new TH1F("pXiApAXiPairs04","ppApApPairs02",3,0,3);
  hist_Pairs[3][2]=new TH1F("pXiApAXiPairsInt","ppApApPairsInt",3,0,3);

  TCanvas *cSumPP=new TCanvas("CanSumPP","CanSumPP",0,0,1500,1100);
  TCanvas *cSumPL=new TCanvas("CanSumPL","CanSumPL",0,0,1500,1100);
  TCanvas *cSumLL=new TCanvas("CanSumLL","CanSumLL",0,0,1500,1100);
  TCanvas *cSumPXi=new TCanvas("CanSumPXi","CanSumPXi",0,0,1500,1100);
  auto* leg= new TLegend(0.45, 0.7, 0.95, 0.95);
  for (int iMult=0;iMult<3;++iMult) {
    TString ppMultName = Form("kThist_CF_pp_%i",iMult);
    hist_CF_pp_ApAp_exp[0][iMult] = Calculate_CF(SEDist[0][0][iMult],MEDist[0][0][iMult],ppMultName.Data(),normleft,normright,suffix,10);

    hist_Pairs[0][0]->SetBinContent(iMult+1,SEDist[0][0][iMult]->Integral(
        SEDist[0][0][iMult]->FindBin(0),SEDist[0][0][iMult]->FindBin(0.2))+
        SEDist[1][1][iMult]->Integral(SEDist[1][1][iMult]->FindBin(0),SEDist[1][1][iMult]->FindBin(0.2)));

    hist_Pairs[0][1]->SetBinContent(iMult+1,SEDist[0][0][iMult]->Integral(
        SEDist[0][0][iMult]->FindBin(0),SEDist[0][0][iMult]->FindBin(0.4))+
        SEDist[1][1][iMult]->Integral(SEDist[1][1][iMult]->FindBin(0),SEDist[1][1][iMult]->FindBin(0.4)));

    hist_Pairs[0][2]->SetBinContent(iMult+1,SEDist[0][0][iMult]->Integral(
        SEDist[0][0][iMult]->FindBin(0),SEDist[0][0][iMult]->FindBin(3))+
        SEDist[1][1][iMult]->Integral(SEDist[1][1][iMult]->FindBin(0),SEDist[1][1][iMult]->FindBin(3)));

    TString ApApMultName = Form("kThist_CF_ApAp_%i",iMult);
    hist_CF_pp_ApAp_exp[1][iMult] = Calculate_CF(SEDist[1][1][iMult],MEDist[1][1][iMult],ApApMultName.Data(),normleft,normright,suffix,10);
    TString SumppApApMultName = Form("kThist_CF_pp_ApAp_exp_sum_%i",iMult);
    hist_CF_pp_ApAp_exp[2][iMult] = add_CF(hist_CF_pp_ApAp_exp[0][iMult],hist_CF_pp_ApAp_exp[1][iMult],SumppApApMultName.Data());
    hist_CF_pp_ApAp_exp[2][iMult]->SetStats(0);
    hist_CF_pp_ApAp_exp[2][iMult]->GetXaxis()->SetTitleOffset(.9);
    hist_CF_pp_ApAp_exp[2][iMult]->GetXaxis()->SetRangeUser(0,0.125);
    hist_CF_pp_ApAp_exp[2][iMult]->GetYaxis()->SetTitleOffset(0.9);
    hist_CF_pp_ApAp_exp[2][iMult]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");
    hist_CF_pp_ApAp_exp[2][iMult]->SetMarkerColor(fColors[iMult]);
    hist_CF_pp_ApAp_exp[2][iMult]->SetLineColor(fColors[iMult]);
    leg->AddEntry(hist_CF_pp_ApAp_exp[2][iMult],Form("Mult. Bin %i",iMult),"lp");

    TPdir->Add(hist_CF_pp_ApAp_exp[2][iMult]);

    TString pLMultName = Form("kThist_CF_pL_%i",iMult);
    hist_CF_pL_ApAL_exp[0][iMult] = Calculate_CF(SEDist[0][2][iMult],MEDist[0][2][iMult],pLMultName.Data(),normleft,normright,suffix,10);

    hist_Pairs[1][0]->SetBinContent(iMult+1,SEDist[0][2][iMult]->Integral(
        SEDist[0][2][iMult]->FindBin(0),SEDist[0][2][iMult]->FindBin(0.2))+
        SEDist[1][3][iMult]->Integral(SEDist[1][3][iMult]->FindBin(0),SEDist[1][3][iMult]->FindBin(0.2)));

    hist_Pairs[1][1]->SetBinContent(iMult+1,SEDist[0][2][iMult]->Integral(
        SEDist[0][2][iMult]->FindBin(0),SEDist[0][2][iMult]->FindBin(0.4))+
        SEDist[1][3][iMult]->Integral(SEDist[1][3][iMult]->FindBin(0),SEDist[1][3][iMult]->FindBin(0.2)));

    hist_Pairs[1][2]->SetBinContent(iMult+1,SEDist[0][2][iMult]->Integral(
        SEDist[0][2][iMult]->FindBin(0),SEDist[0][2][iMult]->FindBin(3))+
        SEDist[1][3][iMult]->Integral(SEDist[1][3][iMult]->FindBin(0),SEDist[1][3][iMult]->FindBin(3)));

    TString ApALMultName = Form("kThist_CF_ApAL_%i",iMult);
    hist_CF_pL_ApAL_exp[1][iMult] = Calculate_CF(SEDist[1][3][iMult],MEDist[1][3][iMult],ApALMultName.Data(),normleft,normright,suffix,10);
    TString SumpLApALMultName = Form("kThist_CF_pL_ApAL_exp_sum_%i",iMult);
    hist_CF_pL_ApAL_exp[2][iMult] = add_CF(hist_CF_pL_ApAL_exp[0][iMult],hist_CF_pL_ApAL_exp[1][iMult],SumpLApALMultName.Data());
    hist_CF_pL_ApAL_exp[2][iMult]->SetStats(0);
    hist_CF_pL_ApAL_exp[2][iMult]->GetXaxis()->SetTitleOffset(.9);
    hist_CF_pL_ApAL_exp[2][iMult]->GetXaxis()->SetRangeUser(0,0.125);
    hist_CF_pL_ApAL_exp[2][iMult]->GetYaxis()->SetTitleOffset(0.9);
    hist_CF_pL_ApAL_exp[2][iMult]->SetMarkerColor(fColors[iMult]);
    hist_CF_pL_ApAL_exp[2][iMult]->SetLineColor(fColors[iMult]);
    hist_CF_pL_ApAL_exp[2][iMult]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

    TPdir->Add(hist_CF_pL_ApAL_exp[2][iMult]);

    TString LLMultName = Form("kThist_CF_LL_%i",iMult);
    hist_CF_LL_ALAL_exp[0][iMult] = Calculate_CF(SEDist[2][2][iMult],MEDist[2][2][iMult],LLMultName.Data(),normleft,normright,suffix,10);
    hist_Pairs[2][0]->SetBinContent(iMult+1,SEDist[2][2][iMult]->Integral(
        SEDist[2][2][iMult]->FindBin(0),SEDist[2][2][iMult]->FindBin(0.2))+
        SEDist[3][3][iMult]->Integral(SEDist[3][3][iMult]->FindBin(0),SEDist[3][3][iMult]->FindBin(0.2)));

    hist_Pairs[2][1]->SetBinContent(iMult+1,SEDist[2][2][iMult]->Integral(
        SEDist[2][2][iMult]->FindBin(0),SEDist[2][2][iMult]->FindBin(0.4))+
        SEDist[3][3][iMult]->Integral(SEDist[3][3][iMult]->FindBin(0),SEDist[3][3][iMult]->FindBin(0.4)));

    hist_Pairs[2][2]->SetBinContent(iMult+1,SEDist[2][2][iMult]->Integral(
        SEDist[2][2][iMult]->FindBin(0),SEDist[2][2][iMult]->FindBin(3))+
        SEDist[3][3][iMult]->Integral(SEDist[3][3][iMult]->FindBin(0),SEDist[3][3][iMult]->FindBin(3)));

    TString ALALMultName = Form("kThist_CF_ALAL_%i",iMult);
    hist_CF_LL_ALAL_exp[1][iMult] = Calculate_CF(SEDist[3][3][iMult],MEDist[3][3][iMult],ALALMultName.Data(),normleft,normright,suffix,10);
    TString SumLLALALMultName = Form("kThist_CF_LL_ALAL_exp_sum_%i",iMult);
    hist_CF_LL_ALAL_exp[2][iMult] = add_CF(hist_CF_LL_ALAL_exp[0][iMult],hist_CF_LL_ALAL_exp[1][iMult],SumLLALALMultName.Data());
    hist_CF_LL_ALAL_exp[2][iMult]->SetStats(0);
    hist_CF_LL_ALAL_exp[2][iMult]->GetXaxis()->SetTitleOffset(.9);
    hist_CF_LL_ALAL_exp[2][iMult]->GetXaxis()->SetRangeUser(0,0.125);
    hist_CF_LL_ALAL_exp[2][iMult]->GetYaxis()->SetTitleOffset(0.9);
    hist_CF_LL_ALAL_exp[2][iMult]->SetMarkerColor(fColors[iMult]);
    hist_CF_LL_ALAL_exp[2][iMult]->SetLineColor(fColors[iMult]);
    hist_CF_LL_ALAL_exp[2][iMult]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

    TPdir->Add(hist_CF_LL_ALAL_exp[2][iMult]);

    TString pXiMultName = Form("kThist_CF_pXi_%i",iMult);
    hist_CF_pXi_ApAXi_exp[0][iMult] = Calculate_CF(SEDist[0][4][iMult],MEDist[0][4][iMult],pXiMultName.Data(),normleft,normright,suffix,10);
    hist_Pairs[3][0]->SetBinContent(iMult+1,SEDist[0][4][iMult]->Integral(
        SEDist[0][4][iMult]->FindBin(0),SEDist[0][4][iMult]->FindBin(0.2))+
        SEDist[1][5][iMult]->Integral(SEDist[1][5][iMult]->FindBin(0),SEDist[1][5][iMult]->FindBin(0.2)));

    hist_Pairs[3][1]->SetBinContent(iMult+1,SEDist[0][2][iMult]->Integral(
        SEDist[0][4][iMult]->FindBin(0),SEDist[0][4][iMult]->FindBin(0.4))+
        SEDist[1][5][iMult]->Integral(SEDist[1][5][iMult]->FindBin(0),SEDist[1][5][iMult]->FindBin(0.4)));

    hist_Pairs[3][2]->SetBinContent(iMult+1,SEDist[0][2][iMult]->Integral(
        SEDist[0][4][iMult]->FindBin(0),SEDist[0][4][iMult]->FindBin(3))+
        SEDist[1][5][iMult]->Integral(SEDist[1][5][iMult]->FindBin(0),SEDist[1][5][iMult]->FindBin(3)));

    TString ApAXiMultName = Form("kThist_CF_ApAXi_%i",iMult);
    hist_CF_pXi_ApAXi_exp[1][iMult] = Calculate_CF(SEDist[1][5][iMult],MEDist[1][5][iMult],ApAXiMultName.Data(),normleft,normright,suffix,10);
    TString SumpXiApAXiMultName = Form("kThist_CF_pXi_ApAXi_exp_sum_%i",iMult);
    hist_CF_pXi_ApAXi_exp[2][iMult] = add_CF(hist_CF_pXi_ApAXi_exp[0][iMult],hist_CF_pXi_ApAXi_exp[1][iMult],SumpXiApAXiMultName.Data());
    hist_CF_pXi_ApAXi_exp[2][iMult]->SetStats(0);
    hist_CF_pXi_ApAXi_exp[2][iMult]->GetXaxis()->SetTitleOffset(.9);
    hist_CF_pXi_ApAXi_exp[2][iMult]->GetXaxis()->SetRangeUser(0,0.125);
    hist_CF_pXi_ApAXi_exp[2][iMult]->GetYaxis()->SetTitleOffset(0.9);
    hist_CF_pXi_ApAXi_exp[2][iMult]->SetMarkerColor(fColors[iMult]);
    hist_CF_pXi_ApAXi_exp[2][iMult]->SetLineColor(fColors[iMult]);
    hist_CF_pXi_ApAXi_exp[2][iMult]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

    TPdir->Add(hist_CF_pXi_ApAXi_exp[2][iMult]);

    if(iMult==0){
      cSumPP->cd();
      hist_CF_pp_ApAp_exp[2][iMult]->DrawCopy();
      cSumPL->cd();
      hist_CF_pL_ApAL_exp[2][iMult]->DrawCopy();
      cSumLL->cd();
      hist_CF_LL_ALAL_exp[2][iMult]->DrawCopy();
      cSumPXi->cd();
      hist_CF_pXi_ApAXi_exp[2][iMult]->DrawCopy();
    } else {
      cSumPP->cd();
      hist_CF_pp_ApAp_exp[2][iMult]->DrawCopy("SAME");
      cSumPL->cd();
      hist_CF_pL_ApAL_exp[2][iMult]->DrawCopy("SAME");
      cSumLL->cd();
      hist_CF_LL_ALAL_exp[2][iMult]->DrawCopy("SAME");
      cSumPXi->cd();
      hist_CF_pXi_ApAXi_exp[2][iMult]->DrawCopy("SAME");
    }
  }
  for (int i=0;i<4;++i) {
    for (int j=0;j<3;++j) {
      TPdir2->Add(hist_Pairs[i][j]);
    }
  }

  cSumPP->cd();
  leg->Draw("SAME");
  cSumPL->cd();
  leg->Draw("SAME");
  cSumLL->cd();
  leg->Draw("SAME");
  cSumPXi->cd();
  leg->Draw("SAME");

  TH1F *kThist_CF_pp_ApAp_exp[3][2];
  TH1F *kThist_CF_pL_ApAL_exp[3][2];
  TH1F *kThist_CF_LL_ALAL_exp[3][2];
  TH1F *kThist_CF_pXi_ApAXi_exp[3][2];
  TH1F *kTCenthist_CF_pp_ApAp_exp[3][2][3];
  TH1F *kTCenthist_CF_pL_ApAL_exp[3][2][3];
  TH1F *kTCenthist_CF_LL_ALAL_exp[3][2][3];
  TH1F *kTCenthist_CF_pXi_ApAXi_exp[3][2][3];
  TCanvas *cSumkTPP=new TCanvas(  "CanSumkTPP", "CanSumkTPP",0,0,1500,1100);
  TCanvas *cSumkTPL=new TCanvas(  "CanSumkTPL", "CanSumkTPL",0,0,1500,1100);
  TCanvas *cSumkTLL=new TCanvas(  "CanSumkTLL", "CanSumkTLL",0,0,1500,1100);
  TCanvas *cSumkTPXi=new TCanvas( "CanSumkTPXi","CanSumkTPXi",0,0,1500,1100);
  auto* legkT= new TLegend(0.65, 0.7, 0.95, 0.95);

  for (int iBins=0;iBins<2;++iBins) {
    TString ppMultName = Form("kThist_CF_pp_%i",iBins);
    kThist_CF_pp_ApAp_exp[0][iBins] = Calculate_CF(SEkTDistBin[0][0][iBins],MEkTDistBin[0][0][iBins],ppMultName.Data(),normleft,normright,suffix,10);
    TString ApApMultName = Form("kThist_CF_ApAp_%i",iBins);
    kThist_CF_pp_ApAp_exp[1][iBins] = Calculate_CF(SEkTDistBin[1][1][iBins],MEkTDistBin[1][1][iBins],ApApMultName.Data(),normleft,normright,suffix,10);
    TString SumppApApMultName = Form("kThist_CF_pp_ApAp_exp_sum_%i",iBins);
    kThist_CF_pp_ApAp_exp[2][iBins] = add_CF(kThist_CF_pp_ApAp_exp[0][iBins],kThist_CF_pp_ApAp_exp[1][iBins],SumppApApMultName.Data());
    kThist_CF_pp_ApAp_exp[2][iBins]->SetStats(0);
    kThist_CF_pp_ApAp_exp[2][iBins]->GetXaxis()->SetTitleOffset(.9);
    kThist_CF_pp_ApAp_exp[2][iBins]->GetXaxis()->SetRangeUser(0,0.125);
    kThist_CF_pp_ApAp_exp[2][iBins]->GetYaxis()->SetTitleOffset(0.9);
    kThist_CF_pp_ApAp_exp[2][iBins]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");
    kThist_CF_pp_ApAp_exp[2][iBins]->SetMarkerColor(fColors[iBins]);
    kThist_CF_pp_ApAp_exp[2][iBins]->SetLineColor(fColors[iBins]);
    legkT->AddEntry(kThist_CF_pp_ApAp_exp[2][iBins],Form("%.1f < kT < %.1f",kTSplit[0][iBins],kTSplit[0][iBins+1]),"lp");

    TPdir->Add(kThist_CF_pp_ApAp_exp[2][iBins]);

    TString pLMultName = Form("kThist_CF_pL_%i",iBins);
    kThist_CF_pL_ApAL_exp[0][iBins] = Calculate_CF(SEkTDistBin[0][2][iBins],MEkTDistBin[0][2][iBins],pLMultName.Data(),normleft,normright,suffix,10);
    TString ApALMultName = Form("kThist_CF_ApAL_%i",iBins);
    kThist_CF_pL_ApAL_exp[1][iBins] = Calculate_CF(SEkTDistBin[1][3][iBins],MEkTDistBin[1][3][iBins],ApALMultName.Data(),normleft,normright,suffix,10);
    TString SumpLApALMultName = Form("kThist_CF_pL_ApAL_exp_sum_%i",iBins);
    kThist_CF_pL_ApAL_exp[2][iBins] = add_CF(kThist_CF_pL_ApAL_exp[0][iBins],kThist_CF_pL_ApAL_exp[1][iBins],SumpLApALMultName.Data());
    kThist_CF_pL_ApAL_exp[2][iBins]->SetStats(0);
    kThist_CF_pL_ApAL_exp[2][iBins]->GetXaxis()->SetTitleOffset(.9);
    kThist_CF_pL_ApAL_exp[2][iBins]->GetXaxis()->SetRangeUser(0,0.125);
    kThist_CF_pL_ApAL_exp[2][iBins]->GetYaxis()->SetTitleOffset(0.9);
    kThist_CF_pL_ApAL_exp[2][iBins]->SetMarkerColor(fColors[iBins]);
    kThist_CF_pL_ApAL_exp[2][iBins]->SetLineColor(fColors[iBins]);
    kThist_CF_pL_ApAL_exp[2][iBins]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

    TPdir->Add(kThist_CF_pL_ApAL_exp[2][iBins]);

    TString LLMultName = Form("kThist_CF_LL_%i",iBins);
    kThist_CF_LL_ALAL_exp[0][iBins] = Calculate_CF(SEkTDistBin[2][2][iBins],MEkTDistBin[2][2][iBins],LLMultName.Data(),normleft,normright,suffix,10);
    TString ALALMultName = Form("kThist_CF_ALAL_%i",iBins);
    kThist_CF_LL_ALAL_exp[1][iBins] = Calculate_CF(SEkTDistBin[3][3][iBins],MEkTDistBin[3][3][iBins],ALALMultName.Data(),normleft,normright,suffix,10);
    TString SumLLALALMultName = Form("kThist_CF_LL_ALAL_exp_sum_%i",iBins);
    kThist_CF_LL_ALAL_exp[2][iBins] = add_CF(kThist_CF_LL_ALAL_exp[0][iBins],kThist_CF_LL_ALAL_exp[1][iBins],SumLLALALMultName.Data());
    kThist_CF_LL_ALAL_exp[2][iBins]->SetStats(0);
    kThist_CF_LL_ALAL_exp[2][iBins]->GetXaxis()->SetTitleOffset(.9);
    kThist_CF_LL_ALAL_exp[2][iBins]->GetXaxis()->SetRangeUser(0,0.125);
    kThist_CF_LL_ALAL_exp[2][iBins]->GetYaxis()->SetTitleOffset(0.9);
    kThist_CF_LL_ALAL_exp[2][iBins]->SetMarkerColor(fColors[iBins]);
    kThist_CF_LL_ALAL_exp[2][iBins]->SetLineColor(fColors[iBins]);
    kThist_CF_LL_ALAL_exp[2][iBins]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

    TPdir->Add(kThist_CF_LL_ALAL_exp[2][iBins]);

    TString pXiBinsName = Form("kThist_CF_pXi_%i",iBins);
    kThist_CF_pXi_ApAXi_exp[0][iBins] = Calculate_CF(SEkTDistBin[0][4][iBins],MEkTDistBin[0][4][iBins],pXiBinsName.Data(),normleft,normright,suffix,10);
    TString ApAXiBinsName = Form("kThist_CF_ApAXi_%i",iBins);
    kThist_CF_pXi_ApAXi_exp[1][iBins] = Calculate_CF(SEkTDistBin[1][5][iBins],MEkTDistBin[1][5][iBins],ApAXiBinsName.Data(),normleft,normright,suffix,10);
    TString SumpXiApAXiBinsName = Form("kThist_CF_pXi_ApAXi_exp_sum_%i",iBins);
    kThist_CF_pXi_ApAXi_exp[2][iBins] = add_CF(kThist_CF_pXi_ApAXi_exp[0][iBins],kThist_CF_pXi_ApAXi_exp[1][iBins],SumpXiApAXiBinsName.Data());
    kThist_CF_pXi_ApAXi_exp[2][iBins]->SetStats(0);
    kThist_CF_pXi_ApAXi_exp[2][iBins]->GetXaxis()->SetTitleOffset(.9);
    kThist_CF_pXi_ApAXi_exp[2][iBins]->GetXaxis()->SetRangeUser(0,0.125);
    kThist_CF_pXi_ApAXi_exp[2][iBins]->GetYaxis()->SetTitleOffset(0.9);
    kThist_CF_pXi_ApAXi_exp[2][iBins]->SetMarkerColor(fColors[iBins]);
    kThist_CF_pXi_ApAXi_exp[2][iBins]->SetLineColor(fColors[iBins]);
    kThist_CF_pXi_ApAXi_exp[2][iBins]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

    TPdir->Add(kThist_CF_pXi_ApAXi_exp[2][iBins]);

    for (int iCent=0;iCent<3;++iCent) {
      TString ppMultCentName = Form("kTCent%ihist_CF_pp_%i",iCent,iBins);
      kTCenthist_CF_pp_ApAp_exp[0][iBins][iCent] = Calculate_CF(SEkTCentDistBin[0][0][iCent][iBins],MEkTCentDistBin[0][0][iCent][iBins],ppMultCentName.Data(),normleft,normright,suffix,10);
      TString ApApMultCentName = Form("kTCent%ihist_CF_ApAp_%i",iCent,iBins);
      kTCenthist_CF_pp_ApAp_exp[1][iBins][iCent] = Calculate_CF(SEkTCentDistBin[1][1][iCent][iBins],MEkTCentDistBin[1][1][iCent][iBins],ApApMultCentName.Data(),normleft,normright,suffix,10);
      TString SumppApApMultCentName = Form("kTCent%ihist_CF_pp_ApAp_exp_sum_%i",iCent,iBins);
      kTCenthist_CF_pp_ApAp_exp[2][iBins][iCent] = add_CF(kTCenthist_CF_pp_ApAp_exp[0][iBins][iCent],kTCenthist_CF_pp_ApAp_exp[1][iBins][iCent],SumppApApMultCentName.Data());
      kTCenthist_CF_pp_ApAp_exp[2][iBins][iCent]->SetStats(0);
      kTCenthist_CF_pp_ApAp_exp[2][iBins][iCent]->GetXaxis()->SetTitleOffset(.9);
      kTCenthist_CF_pp_ApAp_exp[2][iBins][iCent]->GetXaxis()->SetRangeUser(0,0.125);
      kTCenthist_CF_pp_ApAp_exp[2][iBins][iCent]->GetYaxis()->SetTitleOffset(0.9);
      kTCenthist_CF_pp_ApAp_exp[2][iBins][iCent]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");
      kTCenthist_CF_pp_ApAp_exp[2][iBins][iCent]->SetMarkerColor(fColors[iBins]);
      kTCenthist_CF_pp_ApAp_exp[2][iBins][iCent]->SetLineColor(fColors[iBins]);

      TPdir->Add(kTCenthist_CF_pp_ApAp_exp[2][iBins][iCent]);

      TString pLMultCentName = Form("kTCent%ihist_CF_pL_%i",iCent,iBins);
      kTCenthist_CF_pL_ApAL_exp[0][iBins][iCent] = Calculate_CF(SEkTCentDistBin[0][2][iBins][iCent],MEkTCentDistBin[0][2][iBins][iCent],pLMultCentName.Data(),normleft,normright,suffix,10);
      TString ApALMultCentName = Form("kTCent%ihist_CF_ApAL_%i",iCent,iBins);
      kTCenthist_CF_pL_ApAL_exp[1][iBins][iCent] = Calculate_CF(SEkTCentDistBin[1][3][iBins][iCent],MEkTCentDistBin[1][3][iBins][iCent],ApALMultCentName.Data(),normleft,normright,suffix,10);
      TString SumpLApALMultCentName = Form("kTCent%ihist_CF_pL_ApAL_exp_sum_%i",iCent,iBins);
      kTCenthist_CF_pL_ApAL_exp[2][iBins][iCent] = add_CF(kTCenthist_CF_pL_ApAL_exp[0][iBins][iCent],kTCenthist_CF_pL_ApAL_exp[1][iBins][iCent],SumpLApALMultCentName.Data());
      kTCenthist_CF_pL_ApAL_exp[2][iBins][iCent]->SetStats(0);
      kTCenthist_CF_pL_ApAL_exp[2][iBins][iCent]->GetXaxis()->SetTitleOffset(.9);
      kTCenthist_CF_pL_ApAL_exp[2][iBins][iCent]->GetXaxis()->SetRangeUser(0,0.125);
      kTCenthist_CF_pL_ApAL_exp[2][iBins][iCent]->GetYaxis()->SetTitleOffset(0.9);
      kTCenthist_CF_pL_ApAL_exp[2][iBins][iCent]->SetMarkerColor(fColors[iBins]);
      kTCenthist_CF_pL_ApAL_exp[2][iBins][iCent]->SetLineColor(fColors[iBins]);
      kTCenthist_CF_pL_ApAL_exp[2][iBins][iCent]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

      TPdir->Add(kTCenthist_CF_pL_ApAL_exp[2][iBins][iCent]);

      TString LLMultCentName = Form("kTCent%ihist_CF_LL_%i",iCent,iBins);
      kTCenthist_CF_LL_ALAL_exp[0][iBins][iCent] = Calculate_CF(SEkTCentDistBin[2][2][iCent][iBins],MEkTCentDistBin[2][2][iCent][iBins],LLMultCentName.Data(),normleft,normright,suffix,10);
      TString ALALMultCentName = Form("kTCent%ihist_CF_ALAL_%i",iCent,iBins);
      kTCenthist_CF_LL_ALAL_exp[1][iBins][iCent] = Calculate_CF(SEkTCentDistBin[3][3][iCent][iBins],MEkTCentDistBin[3][3][iCent][iBins],ALALMultCentName.Data(),normleft,normright,suffix,10);
      TString SumLLALALMultCentName = Form("kTCent%ihist_CF_LL_ALAL_exp_sum_%i",iCent,iBins);
      kTCenthist_CF_LL_ALAL_exp[2][iBins][iCent] = add_CF(kTCenthist_CF_LL_ALAL_exp[0][iBins][iCent],kTCenthist_CF_LL_ALAL_exp[1][iBins][iCent],SumLLALALMultCentName.Data());
      kTCenthist_CF_LL_ALAL_exp[2][iBins][iCent]->SetStats(0);
      kTCenthist_CF_LL_ALAL_exp[2][iBins][iCent]->GetXaxis()->SetTitleOffset(.9);
      kTCenthist_CF_LL_ALAL_exp[2][iBins][iCent]->GetXaxis()->SetRangeUser(0,0.125);
      kTCenthist_CF_LL_ALAL_exp[2][iBins][iCent]->GetYaxis()->SetTitleOffset(0.9);
      kTCenthist_CF_LL_ALAL_exp[2][iBins][iCent]->SetMarkerColor(fColors[iBins]);
      kTCenthist_CF_LL_ALAL_exp[2][iBins][iCent]->SetLineColor(fColors[iBins]);
      kTCenthist_CF_LL_ALAL_exp[2][iBins][iCent]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

      TPdir->Add(kTCenthist_CF_LL_ALAL_exp[2][iBins][iCent]);

      TString pXiBinsCentName = Form("kTCent%ihist_CF_pXi_%i",iCent,iBins);
      kTCenthist_CF_pXi_ApAXi_exp[0][iBins][iCent] = Calculate_CF(SEkTCentDistBin[0][4][iBins][iCent],MEkTCentDistBin[0][4][iBins][iCent],pXiBinsCentName.Data(),normleft,normright,suffix,10);
      TString ApAXiBinsCentName = Form("kTCent%ihist_CF_ApAXi_%i",iCent,iBins);
      kTCenthist_CF_pXi_ApAXi_exp[1][iBins][iCent] = Calculate_CF(SEkTCentDistBin[1][5][iBins][iCent],MEkTCentDistBin[1][5][iBins][iCent],ApAXiBinsCentName.Data(),normleft,normright,suffix,10);
      TString SumpXiApAXiBinsCentName = Form("kTCent%ihist_CF_pXi_ApAXi_exp_sum_%i",iCent,iBins);
      kTCenthist_CF_pXi_ApAXi_exp[2][iBins][iCent] = add_CF(kTCenthist_CF_pXi_ApAXi_exp[0][iBins][iCent],kTCenthist_CF_pXi_ApAXi_exp[1][iBins][iCent],SumpXiApAXiBinsCentName.Data());
      kTCenthist_CF_pXi_ApAXi_exp[2][iBins][iCent]->SetStats(0);
      kTCenthist_CF_pXi_ApAXi_exp[2][iBins][iCent]->GetXaxis()->SetTitleOffset(.9);
      kTCenthist_CF_pXi_ApAXi_exp[2][iBins][iCent]->GetXaxis()->SetRangeUser(0,0.125);
      kTCenthist_CF_pXi_ApAXi_exp[2][iBins][iCent]->GetYaxis()->SetTitleOffset(0.9);
      kTCenthist_CF_pXi_ApAXi_exp[2][iBins][iCent]->SetMarkerColor(fColors[iBins]);
      kTCenthist_CF_pXi_ApAXi_exp[2][iBins][iCent]->SetLineColor(fColors[iBins]);
      kTCenthist_CF_pXi_ApAXi_exp[2][iBins][iCent]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

      TPdir->Add(kTCenthist_CF_pXi_ApAXi_exp[2][iBins][iCent]);

    }

    if(iBins==0){
      cSumkTPP->cd();
      kThist_CF_pp_ApAp_exp[2][iBins]->DrawCopy();
      cSumkTPL->cd();
      kThist_CF_pL_ApAL_exp[2][iBins]->DrawCopy();
      cSumkTLL->cd();
      kThist_CF_LL_ALAL_exp[2][iBins]->DrawCopy();
      cSumkTPXi->cd();
      hist_CF_pXi_ApAXi_exp[2][iBins]->DrawCopy();
    } else {
      cSumkTPP->cd();
      kThist_CF_pp_ApAp_exp[2][iBins]->DrawCopy("SAME");
      cSumkTPL->cd();
      kThist_CF_pL_ApAL_exp[2][iBins]->DrawCopy("SAME");
      cSumkTLL->cd();
      kThist_CF_LL_ALAL_exp[2][iBins]->DrawCopy("SAME");
      cSumkTPXi->cd();
      kThist_CF_pXi_ApAXi_exp[2][iBins]->DrawCopy("SAME");
    }
  }
  cSumkTPP->cd();
  legkT->Draw("SAME");
  cSumkTPL->cd();
  legkT->Draw("SAME");
  cSumkTLL->cd();
  legkT->Draw("SAME");
  cSumkTPXi->cd();
  legkT->Draw("SAME");

  TH1F *mThist_CF_pp_ApAp_exp[3][2];
  TH1F *mThist_CF_pL_ApAL_exp[3][2];
  TH1F *mThist_CF_LL_ALAL_exp[3][2];
  TH1F *mThist_CF_pXi_ApAXi_exp[3][2];
  TCanvas *cSummTPP=new TCanvas(  "CanSummTPP", "CanSummTPP",0,0,1500,1100);
  TCanvas *cSummTPL=new TCanvas(  "CanSummTPL", "CanSummTPL",0,0,1500,1100);
  TCanvas *cSummTLL=new TCanvas(  "CanSummTLL", "CanSummTLL",0,0,1500,1100);
  TCanvas *cSummTPXi=new TCanvas( "CanSummTPXi","CanSummTPXi",0,0,1500,1100);

  for (int iBins=0;iBins<2;++iBins) {
    TString ppMultName = Form("mThist_CF_pp_%i",iBins);
    mThist_CF_pp_ApAp_exp[0][iBins] = Calculate_CF(SEmTDistBin[0][0][iBins],MEmTDistBin[0][0][iBins],ppMultName.Data(),normleft,normright,suffix,10);
    TString ApApMultName = Form("mThist_CF_ApAp_%i",iBins);
    mThist_CF_pp_ApAp_exp[1][iBins] = Calculate_CF(SEmTDistBin[1][1][iBins],MEmTDistBin[1][1][iBins],ApApMultName.Data(),normleft,normright,suffix,10);
    TString SumppApApMultName = Form("mThist_CF_pp_ApAp_exp_sum_%i",iBins);
    mThist_CF_pp_ApAp_exp[2][iBins] = add_CF(mThist_CF_pp_ApAp_exp[0][iBins],mThist_CF_pp_ApAp_exp[1][iBins],SumppApApMultName.Data());
    mThist_CF_pp_ApAp_exp[2][iBins]->SetStats(0);
    mThist_CF_pp_ApAp_exp[2][iBins]->GetXaxis()->SetTitleOffset(.9);
    mThist_CF_pp_ApAp_exp[2][iBins]->GetXaxis()->SetRangeUser(0,0.125);
    mThist_CF_pp_ApAp_exp[2][iBins]->GetYaxis()->SetTitleOffset(0.9);
    mThist_CF_pp_ApAp_exp[2][iBins]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");
    mThist_CF_pp_ApAp_exp[2][iBins]->SetMarkerColor(fColors[iBins]);
    mThist_CF_pp_ApAp_exp[2][iBins]->SetLineColor(fColors[iBins]);

    TPdir->Add(mThist_CF_pp_ApAp_exp[2][iBins]);

    TString pLMultName = Form("mThist_CF_pL_%i",iBins);
    mThist_CF_pL_ApAL_exp[0][iBins] = Calculate_CF(SEmTDistBin[0][2][iBins],MEmTDistBin[0][2][iBins],pLMultName.Data(),normleft,normright,suffix,10);
    TString ApALMultName = Form("mThist_CF_ApAL_%i",iBins);
    mThist_CF_pL_ApAL_exp[1][iBins] = Calculate_CF(SEmTDistBin[1][3][iBins],MEmTDistBin[1][3][iBins],ApALMultName.Data(),normleft,normright,suffix,10);
    TString SumpLApALMultName = Form("mThist_CF_pL_ApAL_exp_sum_%i",iBins);
    mThist_CF_pL_ApAL_exp[2][iBins] = add_CF(mThist_CF_pL_ApAL_exp[0][iBins],mThist_CF_pL_ApAL_exp[1][iBins],SumpLApALMultName.Data());
    mThist_CF_pL_ApAL_exp[2][iBins]->SetStats(0);
    mThist_CF_pL_ApAL_exp[2][iBins]->GetXaxis()->SetTitleOffset(.9);
    mThist_CF_pL_ApAL_exp[2][iBins]->GetXaxis()->SetRangeUser(0,0.125);
    mThist_CF_pL_ApAL_exp[2][iBins]->GetYaxis()->SetTitleOffset(0.9);
    mThist_CF_pL_ApAL_exp[2][iBins]->SetMarkerColor(fColors[iBins]);
    mThist_CF_pL_ApAL_exp[2][iBins]->SetLineColor(fColors[iBins]);
    mThist_CF_pL_ApAL_exp[2][iBins]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

    TPdir->Add(mThist_CF_pL_ApAL_exp[2][iBins]);

    TString LLMultName = Form("mThist_CF_LL_%i",iBins);
    mThist_CF_LL_ALAL_exp[0][iBins] = Calculate_CF(SEmTDistBin[2][2][iBins],MEmTDistBin[2][2][iBins],LLMultName.Data(),normleft,normright,suffix,10);
    TString ALALMultName = Form("mThist_CF_ALAL_%i",iBins);
    mThist_CF_LL_ALAL_exp[1][iBins] = Calculate_CF(SEmTDistBin[3][3][iBins],MEmTDistBin[3][3][iBins],ALALMultName.Data(),normleft,normright,suffix,10);
    TString SumLLALALMultName = Form("mThist_CF_LL_ALAL_exp_sum_%i",iBins);
    mThist_CF_LL_ALAL_exp[2][iBins] = add_CF(mThist_CF_LL_ALAL_exp[0][iBins],mThist_CF_LL_ALAL_exp[1][iBins],SumLLALALMultName.Data());
    mThist_CF_LL_ALAL_exp[2][iBins]->SetStats(0);
    mThist_CF_LL_ALAL_exp[2][iBins]->GetXaxis()->SetTitleOffset(.9);
    mThist_CF_LL_ALAL_exp[2][iBins]->GetXaxis()->SetRangeUser(0,0.125);
    mThist_CF_LL_ALAL_exp[2][iBins]->GetYaxis()->SetTitleOffset(0.9);
    mThist_CF_LL_ALAL_exp[2][iBins]->SetMarkerColor(fColors[iBins]);
    mThist_CF_LL_ALAL_exp[2][iBins]->SetLineColor(fColors[iBins]);
    mThist_CF_LL_ALAL_exp[2][iBins]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

    TPdir->Add(mThist_CF_LL_ALAL_exp[2][iBins]);

    TString pXiBinsName = Form("mThist_CF_pXi_%i",iBins);
    mThist_CF_pXi_ApAXi_exp[0][iBins] = Calculate_CF(SEmTDistBin[0][4][iBins],MEmTDistBin[0][4][iBins],pXiBinsName.Data(),normleft,normright,suffix,10);
    TString ApAXiBinsName = Form("mThist_CF_ApAXi_%i",iBins);
    mThist_CF_pXi_ApAXi_exp[1][iBins] = Calculate_CF(SEmTDistBin[1][5][iBins],MEmTDistBin[1][5][iBins],ApAXiBinsName.Data(),normleft,normright,suffix,10);
    TString SumpXiApAXiBinsName = Form("mThist_CF_pXi_ApAXi_exp_sum_%i",iBins);
    mThist_CF_pXi_ApAXi_exp[2][iBins] = add_CF(mThist_CF_pXi_ApAXi_exp[0][iBins],mThist_CF_pXi_ApAXi_exp[1][iBins],SumpXiApAXiBinsName.Data());
    mThist_CF_pXi_ApAXi_exp[2][iBins]->SetStats(0);
    mThist_CF_pXi_ApAXi_exp[2][iBins]->GetXaxis()->SetTitleOffset(.9);
    mThist_CF_pXi_ApAXi_exp[2][iBins]->GetXaxis()->SetRangeUser(0,0.125);
    mThist_CF_pXi_ApAXi_exp[2][iBins]->GetYaxis()->SetTitleOffset(0.9);
    mThist_CF_pXi_ApAXi_exp[2][iBins]->SetMarkerColor(fColors[iBins]);
    mThist_CF_pXi_ApAXi_exp[2][iBins]->SetLineColor(fColors[iBins]);
    mThist_CF_pXi_ApAXi_exp[2][iBins]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");

    TPdir->Add(mThist_CF_pXi_ApAXi_exp[2][iBins]);

    if(iBins==0){
      cSummTPP->cd();
      mThist_CF_pp_ApAp_exp[2][iBins]->DrawCopy();
      cSummTPL->cd();
      mThist_CF_pL_ApAL_exp[2][iBins]->DrawCopy();
      cSummTLL->cd();
      mThist_CF_LL_ALAL_exp[2][iBins]->DrawCopy();
      cSummTPXi->cd();
      hist_CF_pXi_ApAXi_exp[2][iBins]->DrawCopy();
    } else {
      cSummTPP->cd();
      mThist_CF_pp_ApAp_exp[2][iBins]->DrawCopy("SAME");
      cSummTPL->cd();
      mThist_CF_pL_ApAL_exp[2][iBins]->DrawCopy("SAME");
      cSummTLL->cd();
      mThist_CF_LL_ALAL_exp[2][iBins]->DrawCopy("SAME");
      cSummTPXi->cd();
      mThist_CF_pXi_ApAXi_exp[2][iBins]->DrawCopy("SAME");
    }
  }
  cSummTPP->cd();
  cSummTPL->cd();
  cSummTLL->cd();
  cSummTPXi->cd();

  TFile *output=new TFile("ppMuti.root","RECREATE");
  TDirectoryFile* dirOutput=new TDirectoryFile("PWGCF_PLFemto_0","PWGCF_PLFemto_0");
  dirOutput->Add(TPdir);
  dirOutput->Add(TPdir2);
  dirOutput->Write("PWGCF_PLFemto_0",TObject::kSingleKey);
}
