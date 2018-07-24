#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include <vector>

std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-3, kYellow-7};
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
  histo->SetMarkerSize(1.5);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleGraph(TGraph *histo, int marker, int color)
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
  histo->SetMarkerSize(1.5);
  histo->SetLineWidth(2);
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
  Double_t norm_relK = 0;
  if(strcmp(folder, "") == 0) {
    double IntegralSE = histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright));
    double IntegralME = histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright));
     norm_relK =  IntegralSE/IntegralME;

    std::cout << IntegralSE << '\t' << IntegralME << '\t' << norm_relK << std::endl;
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


TH1F* ConvertToOtherUnit(TH1F* HistCF,int Scale,TString name) {
  int nBins = HistCF->GetNbinsX();
  double kMin=HistCF->GetXaxis()->GetXmin();
  double kMax=HistCF->GetXaxis()->GetXmax();
  TH1F* HistScaled = new TH1F(name.Data(),name.Data(),nBins,kMin*Scale,kMax*Scale);
  for (int iBin = 1; iBin <= nBins; ++iBin ) {
    HistScaled->SetBinContent(iBin,HistCF->GetBinContent(iBin));
    HistScaled->SetBinError(iBin,HistCF->GetBinError(iBin));
  }
  return HistScaled;
}

void ShiftForEmptyAndRebin(TH1F* InputSE, TH1F* InputME, int rebin,
                           float normleft, float normright,
                           TString outname, TH1F* output[3]) {
  TString nameSE = Form("%s_Shifted",InputSE->GetName());
  TString nameME = Form("%s_Shifted",InputME->GetName());
  int nBins = InputSE->GetXaxis()->GetNbins();
  int nBinsRebin = (int)nBins/rebin;
//  std::cout << nBinsRebin << std::endl;
  float kMin = InputSE->GetXaxis()->GetXmin();
  float kMax = InputSE->GetXaxis()->GetXmax();
  float binWidth = InputSE->GetXaxis()->GetBinWidth(1);
  float binWidthnew = binWidth*rebin;
//  TH1F* Shifted[3];
  int binNew=1;
  int effCounter=0;
  int nEntriesSE = 0;
  int nEntriesME = 0;
  bool burnin = true;
  bool construct = true;
  for (int iBins = 1; iBins <= nBins; ++iBins) {
    nEntriesSE+=InputSE->GetBinContent(iBins);
    if (burnin && (nEntriesSE == 0)) {
      kMin+=binWidth;
//      std::cout << iBins << '\t' << rebin << std::endl;
      if (iBins%rebin==0) {
        nBinsRebin--;
      }
      continue;
    } else {
      if (construct) {
//        std::cout << kMin << '\t' << binWidthnew << '\t' << nBinsRebin << std::endl;
        kMax = kMin+binWidthnew*nBinsRebin;
        std::cout << kMax << std::endl;
        output[0]=new TH1F(nameSE.Data(),nameSE.Data(),nBinsRebin,kMin,kMax);
        output[0]->Sumw2();
        output[1]=new TH1F(nameME.Data(),nameME.Data(),nBinsRebin,kMin,kMax);
        output[1]->Sumw2();
        burnin = false;
        construct = false;
      }
    }
    effCounter++;
    nEntriesME+=InputME->GetBinContent(iBins);

    if (effCounter%rebin==0) {
      output[0]->SetBinContent(binNew,nEntriesSE);
      output[0]->SetBinError(binNew,TMath::Sqrt(nEntriesSE));
      output[1]->SetBinContent(binNew,nEntriesME);
      output[1]->SetBinError(binNew,TMath::Sqrt(nEntriesME));

      nEntriesSE=0;
      nEntriesME=0;

      effCounter=0;
      binNew++;
    }
  }
  output[2] = Calculate_CF(output[0],output[1],outname,normleft,normright, "", 0);
//  output=Shifted;
}
//
void VariableBinningCF(TH1F* InputSE, TH1F* InputME, double relErrMin,
                        float normleft, float normright,
                        TString outname, TH1F* output[3]) {
  TString nameSE = Form("%s_IndepBin",InputSE->GetName());
  TString nameME = Form("%s_IndepBin",InputME->GetName());
  if (relErrMin > 1 || relErrMin <= 0) {
    std::cout << "are you sure you want relErrMin too be larger than 1? Should be a value larger 0 and less than 1. \n";
  }
  int threshold = 1/(float)::pow(relErrMin,2);
  std::vector<float> xValues;
  std::vector<float> ySECounts;
  std::vector<float> yMECounts;
  bool burnyIn=true;
  int binContentSE = 0;
  int binContentME = 0;
  double xMin = 0;
  for (int iBin = 1; iBin <= InputSE->GetNbinsX(); ++iBin) {
    binContentSE += InputSE->GetBinContent(iBin);
//    make sure we omit the first few empty bins
    if (burnyIn && (binContentSE == 0)) {
      continue;
    } else {
      if (burnyIn) {
        xMin = InputSE->GetXaxis()->GetBinLowEdge(iBin);
      }
      burnyIn=false;
    }
//    as soon as we are done with this we can start looking for the number of bins we need to get the relative error of our dreams
    binContentME += InputME->GetBinContent(iBin);
    if (threshold < binContentSE || iBin == InputSE->GetNbinsX()) {
      // if we have enough entries to minimize our rel err create a new bin
      ySECounts.push_back(binContentSE);
      yMECounts.push_back(binContentME);
      xValues.push_back(xMin);
//      std::cout << xMin << std::endl;
      xMin = InputSE->GetXaxis()->GetBinLowEdge(iBin+1);
      binContentSE=0;
      binContentME=0;
    }
  }
  xValues.push_back(InputSE->GetXaxis()->GetBinLowEdge(InputSE->GetNbinsX()));
  const int size = xValues.size();
  float arrXval[size];
  std::copy(xValues.begin(), xValues.end(), arrXval);
//  for (int iBin=1;iBin<=size;++iBin) {
//    std::cout << arrXval[iBin-1] << '\t';
//    if (iBin %20 == 0) {
//      std::cout << std::endl;
//    }
//  }
  output[0] = new TH1F(nameSE.Data(),nameSE.Data(),size-1,arrXval);
  output[0]->Sumw2();
  output[1] = new TH1F(nameME.Data(),nameME.Data(),size-1,arrXval);
  output[1]->Sumw2();
  for (int iBin=1;iBin<size;++iBin) {
    output[0]->SetBinContent(iBin,ySECounts.at(iBin-1));
    output[0]->SetBinError(iBin,TMath::Sqrt(ySECounts.at(iBin-1)));
    output[1]->SetBinContent(iBin,yMECounts.at(iBin-1));
    output[1]->SetBinError(iBin,TMath::Sqrt(yMECounts.at(iBin-1)));
  }
  output[2] = Calculate_CF(output[0],output[1],outname,normleft,normright, "", 0);
}

void Rebinned(TH1F* InputppSE, TH1F* InputppME, TH1F* InputApApSE, TH1F* InputApApME,
              int rebin=1,TString partPair="")
{
  const float normleft = 200;
  const float normright = 400;
  TH1F* ppSE=(TH1F*)InputppSE->Clone("ppSE");
  TH1F* ppME=(TH1F*)InputppME->Clone("ppME");
  TH1F* ApApSE=(TH1F*)InputApApSE->Clone("ApApSE");
  TH1F* ApApME=(TH1F*)InputApApME->Clone("ApApME");

  TH1F* ppSESum = (TH1F*)InputppSE->Clone("SumppSE");
  ppSESum->Add(ApApSE);

  TH1F* ppMESum = (TH1F*)InputppME->Clone("SumppME");
  ppMESum->Add(ApApME);

  TH1F* ShiftedCF_Summed[3];

  ShiftForEmptyAndRebin(
      ppSESum,ppMESum,rebin,normleft,normright,"hSummedCkNormShifted",ShiftedCF_Summed);

  ppSESum->Rebin(rebin);
  ppMESum->Rebin(rebin);

  TH1F* ppApApCFSum = Calculate_CF(
      ppSESum,ppMESum,"hCKPartSumNorm",normleft,normright,"",0);

  TH1F *ppApApCF[3];
  TH1F* ShiftedCF_pp[3];

  ShiftForEmptyAndRebin(
      ppSE,ppME,rebin,normleft,normright,"hCkPartNormShifted",ShiftedCF_pp);

  ppSE->Rebin(rebin);
  ppME->Rebin(rebin);
  ppApApCF[0] = Calculate_CF(
      ppSE,ppME,"hCkPartNorm",normleft,normright, "", 0);

  TH1F *ShiftedCF_ApAp[3];
  ShiftForEmptyAndRebin(
      ApApSE,ApApME,rebin,normleft,normright,"hCkAntiPartNormShifted",ShiftedCF_ApAp);

  ApApSE->Rebin(rebin);
  ApApME->Rebin(rebin);
  ppApApCF[1] = Calculate_CF(
      ApApSE,ApApME,"hCkAntiPartNorm",normleft,normright, "", 0);
  ppApApCF[2] = add_CF(
      ppApApCF[0],ppApApCF[1],"hCkTotNormWeight");

  TH1F *ShiftedCF_Sum = add_CF(
      ShiftedCF_pp[2],ShiftedCF_ApAp[2],"hCkTotNormWeight_Shifted");



  TFile *output;

  output=TFile::Open(Form("CFOutput_%s_Rebin_%d.root",partPair.Data(),rebin),"RECREATE");
  ppSE->Write();
  ppME->Write();
  ApApSE->Write();
  ApApME->Write();
  ppApApCF[0]->Write();
  ppApApCF[1]->Write();
  ppApApCF[2]->Write();
  ShiftedCF_pp[0]->Write();
  ShiftedCF_pp[1]->Write();
  ShiftedCF_pp[2]->Write();
  ShiftedCF_ApAp[0]->Write();
  ShiftedCF_ApAp[1]->Write();
  ShiftedCF_ApAp[2]->Write();
  ShiftedCF_Sum->Write();
  ppSESum->Write();
  ppMESum->Write();
  ppApApCFSum->Write();
  ShiftedCF_Summed[0]->Write();
  ShiftedCF_Summed[1]->Write();
  ShiftedCF_Summed[2]->Write();

  output->Close();

  delete ppSE;
  delete ppME;
  delete ApApSE;
  delete ApApME;
  delete ppApApCF[0];
  delete ppApApCF[1];
  delete ppApApCF[2];
  delete ShiftedCF_pp[0];
  delete ShiftedCF_pp[1];
  delete ShiftedCF_pp[2];
  delete ShiftedCF_ApAp[0];
  delete ShiftedCF_ApAp[1];
  delete ShiftedCF_ApAp[2];
  delete ShiftedCF_Sum;
  delete output;

  return;
}

void BlindBinning(TH1F* InputppSE, TH1F* InputppME, TH1F* InputApApSE, TH1F* InputApApME,
                  float confidence=0.3, TString partPair="")
{
  const float normleft = 200;
  const float normright = 400;

  TH1F* ppSE=(TH1F*)InputppSE->Clone("ppSE");
  TH1F* ppME=(TH1F*)InputppME->Clone("ppME");
  TH1F* ApApSE=(TH1F*)InputApApSE->Clone("ApApSE");
  TH1F* ApApME=(TH1F*)InputApApME->Clone("ApApME");

  TH1F* IndepBinnedCF_pp[3];
  VariableBinningCF(
      ppSE,ppME,confidence,normleft,normright,"hCkPartNormIndepBinning",IndepBinnedCF_pp);

  TH1F *IndepBinnedCF_ApAp[3];
  VariableBinningCF(
      ApApSE,ApApME,confidence,normleft,normright,"hCkAntiPartNormIndepBinning",IndepBinnedCF_ApAp);

  TH1F *IndepBinningCF_Sum =
      add_CF(IndepBinnedCF_pp[2],IndepBinnedCF_ApAp[2],"hCkTotNormWeightIndepBinning");
  TFile *output;
  output=TFile::Open(Form("CFOutput_%s_RelFac_%.2f.root",partPair.Data(),confidence),"RECREATE");

  ppSE->Write();
  ppME->Write();
  ApApSE->Write();
  ApApME->Write();
  IndepBinnedCF_pp[0]->Write();
  IndepBinnedCF_pp[1]->Write();
  IndepBinnedCF_pp[2]->Write();
  IndepBinnedCF_ApAp[0]->Write();
  IndepBinnedCF_ApAp[1]->Write();
  IndepBinnedCF_ApAp[2]->Write();
  IndepBinningCF_Sum->Write();
  output->Close();

  delete ppSE;
  delete ppME;
  delete ApApSE;
  delete ApApME;
  delete IndepBinnedCF_pp[0];
  delete IndepBinnedCF_pp[1];
  delete IndepBinnedCF_pp[2];
  delete IndepBinnedCF_ApAp[0];
  delete IndepBinnedCF_ApAp[1];
  delete IndepBinnedCF_ApAp[2];
  delete IndepBinningCF_Sum;
  delete output;

  return;
}


void ConvertToCats(const char *fileExp, TString partPair) {
  TFile* _file0=TFile::Open(fileExp,"READ");
  const float normleft = 200;
  const float normright = 400;

 TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny("HMResults"));
  TList *Results;
  dirResults->GetObject("HMResults",Results);
  TList* tmpFolder;

  TH1F *ppSE;
  TH1F *ppME;

  TH1F *ApApSE;
  TH1F *ApApME;


  TString FolderNameParticle;
  TString FolderNameAntiParticle;
  if (partPair.Contains("pp")) {
    FolderNameParticle="Particle0_Particle0";
    FolderNameAntiParticle="Particle1_Particle1";
  } else if (partPair.Contains("pL")) {
    FolderNameParticle="Particle0_Particle2";
    FolderNameAntiParticle="Particle1_Particle3";
  } else if (partPair.Contains("LL")) {
    FolderNameParticle="Particle2_Particle2";
    FolderNameAntiParticle="Particle3_Particle3";
  } else if (partPair.Contains("pXi")) {
    FolderNameParticle="Particle0_Particle4";
    FolderNameAntiParticle="Particle1_Particle5";
  }

  tmpFolder=(TList*)Results->FindObject(FolderNameParticle.Data());
  ppSE = ConvertToOtherUnit(
      ((TH1F*)tmpFolder->FindObject(Form("SEDist_%s",FolderNameParticle.Data()))),
      1000,"hPartSe");
  if (!ppSE ) {
    std::cout << "No pp SE Hist \n";
  }

  ppME = ConvertToOtherUnit(
      ((TH1F*)tmpFolder->FindObject(Form("MEDist_%s",FolderNameParticle.Data()))),
      1000,"hPartMe");
  if (!ppME ) {
    std::cout << "No pp ME Hist \n";
  }

  tmpFolder=(TList*)Results->FindObject(FolderNameAntiParticle.Data());
  ApApSE = ConvertToOtherUnit(
      ((TH1F*)tmpFolder->FindObject(Form("SEDist_%s",FolderNameAntiParticle.Data()))),
      1000,"hAntiPartSe");
  if (!ApApSE ) {
    std::cout << "No ApAp SE Hist \n";
  }
  ApApME = ConvertToOtherUnit(
      ((TH1F*)tmpFolder->FindObject(Form("MEDist_%s",FolderNameAntiParticle.Data()))),
      1000,"hAntiPartMe");
  if (!ApApME ) {
    std::cout << "No ApAp ME Hist \n";
  }

  for (int iRebin = 1; iRebin < 6; ++iRebin) {
//    if (iRebin == 4) continue;
    Rebinned(ppSE,ppME,ApApSE,ApApME,iRebin,partPair);
  }
  float confLevel = 0.1;
  for (int iCon = 0; iCon < 5; ++iCon ) {
    BlindBinning(ppSE,ppME,ApApSE,ApApME,confLevel,partPair);
    confLevel+=0.05;
  }
  return; 
}
