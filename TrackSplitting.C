
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

void TrackSplitting(const char *filename,const char *prefix) {
  SetStyle(false,true);
  const float normleft = 0.2;
  const float normright = 0.4;
  TFile *file=TFile::Open(filename);

  TString RelKNames[6]={"Proton","AntiProton","Lambda","AntiLambda","Xi","AntiXi"};
  TString DaugNames[3]={"#pi","p","#pi"};
  TString ResultsName="";
  ResultsName+=prefix;
  ResultsName+="ResultQA";

  TH2F *SEDist[6][6][3][9];
  TH2F *SESumDist[4][3][9];
  TH2F *MEDist[6][6][3][9];
  TH2F *MESumDist[4][3][9];
  TH2F *CFDist[6][6][3][9];
  TH2F *CFSumDist[4][3][9];
  TCanvas *c1[5][3];
  TDirectoryFile *dirQAResults=(TDirectoryFile*)(file->FindObjectAny(ResultsName.Data()));
  if (dirQAResults) {
    TList *tmp=dirQAResults->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *Results;
    dirQAResults->GetObject(name,Results);
    if (!Results) {
      std::cout << "No QAResult\n";
    } else {
      for (int iPart1=0;iPart1<6;++iPart1) {
        for (int iPart2=iPart1;iPart2<6;++iPart2) {
          TString folderName=Form("QA_Particle%i_Particle%i",iPart1,iPart2);
          TList* tmpFolder=(TList*)Results->FindObject(folderName.Data());
          if (!tmpFolder) {
            std::cout << folderName.Data() <<" not Found\n";
          } else {
            for (int iDaug=0;iDaug<3;++iDaug) {\
              if (iPart1==0&&iPart2==0&&iDaug<1) {
                c1[iPart2][iDaug]=new TCanvas("TracksPP","TracksPP",2000,2000);
                c1[iPart2][iDaug]->Divide(3,3);
              } else if (iPart1==0&&iPart2==2&&iDaug<2) {
                TString CanName=Form("TracksPL_Daug%i",iDaug);
                c1[iPart2][iDaug]=new TCanvas(CanName.Data(),CanName.Data(),2000,2000);
                c1[iPart2][iDaug]->Divide(3,3);
              } else if (iPart1==0&&iPart2==4) {
                TString CanName=Form("TracksPXi_Daug%i",iDaug);
                c1[iPart2][iDaug]=new TCanvas(CanName.Data(),CanName.Data(),2000,2000);
                c1[iPart2][iDaug]->Divide(3,3);
              }
            for (int iRad=0;iRad<9;++iRad) {
              TString SEName=Form("SERad_%i_Particle%i_Particle%i_Daug%i",iRad,iPart1,iPart2,iDaug);
              SEDist[iPart1][iPart2][iDaug][iRad]=(TH2F*)tmpFolder->FindObject(SEName.Data());
              SEDist[iPart1][iPart2][iDaug][iRad]->SetDirectory(0);
              SEDist[iPart1][iPart2][iDaug][iRad]->Scale(1./SEDist[iPart1][iPart2][iDaug][iRad]->GetEntries());
              if (!SEDist[iPart1][iPart2]) {
                std::cout << SEName.Data() << " not Found\n";
                std::cout << SEName.Data() << " not Found\n";
              }

              TString MEName=Form("MERad_%i_Particle%i_Particle%i_Daug%i",iRad,iPart1,iPart2,iDaug);
              MEDist[iPart1][iPart2][iDaug][iRad]=(TH2F*)tmpFolder->FindObject(MEName.Data());
              MEDist[iPart1][iPart2][iDaug][iRad]->SetDirectory(0);
              MEDist[iPart1][iPart2][iDaug][iRad]->Scale(1./MEDist[iPart1][iPart2][iDaug][iRad]->GetEntries());
              if (!MEDist[iPart1][iPart2]) {
                std::cout << MEName.Data() << " not Found\n";
                std::cout << MEName.Data() << " not Found\n";
              }
              TString CFName="";
              if ((iPart1==0&&iPart2==0)||(iPart1==1&&iPart2==1)) {
                CFName+= Form("%s-%s d#eta d#Phi at Radius %i",RelKNames[iPart1].Data(),RelKNames[iPart2].Data(),iRad);
              } else {
                CFName+= Form("%s-%s %s Daugter d#eta d#Phi at Radius %i",RelKNames[iPart1].Data(),RelKNames[iPart2].Data(),DaugNames[iDaug].Data(),iRad);
              }
              CFDist[iPart1][iPart2][iDaug][iRad]=(TH2F*)SEDist[iPart1][iPart2][iDaug][iRad]->Clone(CFName.Data());
              CFDist[iPart1][iPart2][iDaug][iRad]->SetTitle(CFName.Data());
              CFDist[iPart1][iPart2][iDaug][iRad]->Divide(MEDist[iPart1][iPart2][iDaug][iRad]);
              CFDist[iPart1][iPart2][iDaug][iRad]->GetXaxis()->SetRangeUser(0,1.3);
              if (iPart1==0&&iPart2==0&&iDaug<1) {
                c1[iPart2][iDaug]->cd(iRad+1);
                CFDist[iPart1][iPart2][iDaug][iRad]->DrawCopy("COLZ");
              } else if (iPart1==0&&iPart2==2&&iDaug<2) {
                TString CanName=Form("TracksPL_Daug%i",iDaug);
                c1[iPart2][iDaug]->cd(iRad+1);
                CFDist[iPart1][iPart2][iDaug][iRad]->DrawCopy("COLZ");
              } else if (iPart1==0&&iPart2==4) {
                TString CanName=Form("TracksPXi_Daug%i",iDaug);
                c1[iPart2][iDaug]->cd(iRad+1);
                CFDist[iPart1][iPart2][iDaug][iRad]->DrawCopy("COLZ");
              }
            }
            }
          }
        }
      }
    }
  } else {
    std::cout << "No directory \n";
  }
}
