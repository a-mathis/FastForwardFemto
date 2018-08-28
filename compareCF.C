#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

std::vector<int> fFillColors = {kGray + 1,  kRed - 10,    kBlue - 9,
                                kGreen - 8, kMagenta - 9, kOrange - 9,
                                kCyan - 3,  kYellow - 7};
std::vector<int> fColors = {kBlack,       kRed + 1,    kBlue + 2, kGreen + 3,
                            kMagenta + 1, kOrange - 1, kCyan + 2, kYellow + 2};
std::vector<int> fMarkers = {
    kFullCircle, kFullSquare, kOpenCircle,  kOpenSquare, kOpenDiamond,
    kOpenCross,  kFullCross,  kFullDiamond, kFullStar,   kOpenStar};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette = false, bool title = false) {
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if (graypalette)
    gStyle->SetPalette(8, 0);
  else
    gStyle->SetPalette(1);
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
  gStyle->SetLabelSize(0.045, "xyz");
  gStyle->SetLabelOffset(0.01, "y");
  gStyle->SetLabelOffset(0.01, "x");
  gStyle->SetLabelColor(kBlack, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);

  const int NRGBs = 6;
  Double_t stops[NRGBs];
  for (int i = 0; i < NRGBs; ++i) stops[i] = float(i) / (NRGBs - 1);

  Double_t red[NRGBs] = {1.,         29. / 255., 25. / 255.,
                         27. / 255., 32. / 255., 24. / 255.};
  Double_t green[NRGBs] = {1.,          221. / 255., 160. / 255.,
                           113. / 255., 74. / 255.,  37. / 255.};
  Double_t blue[NRGBs] = {1.,          221. / 255., 184. / 255.,
                          154. / 255., 129. / 255., 98. / 255.};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHisto(TH1* histo, int marker, int color) {
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
void SetStyleGraph(TGraph* histo, int marker, int color) {
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
TH1F* Calculate_CF(TH1F* histRE_relK, TH1F* histME_relK, TString CFname,
                   Double_t normleft, Double_t normright, const char* folder,
                   float spinningDepth) {
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  if (strcmp(folder, "") == 0) {
    Double_t norm_relK =
        histRE_relK->Integral(histRE_relK->FindBin(normleft),
                              histRE_relK->FindBin(normright)) /
        histME_relK->Integral(histME_relK->FindBin(normleft),
                              histME_relK->FindBin(normright));
    Hist_CF->Divide(histRE_relK, histME_relK, 1, norm_relK);
  } else {
    histME_relK->Scale(1.f / spinningDepth);
    Hist_CF->Divide(histRE_relK, histME_relK);
  }

  return Hist_CF;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *RatioUncertainties(TH1F *histOld, TH1F *histNew) {
  TH1F *hist = (TH1F*)histNew->Clone(Form("%s_Clone", histNew->GetName()));
  hist->GetYaxis()->SetTitle("Ratio uncertainties old/new");
  for(int i = 0; i<hist->GetNbinsX(); ++i) {
    hist->SetBinContent(i+1, histOld->GetBinError(i+1) / histNew->GetBinError(i+1));
  }
  return hist;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* add_CF(TH1F* hist_CF1, TH1F* hist_CF2, TString HistName) {
  // Calculate CFs with error weighting
  TH1F* hist_CF_sum = (TH1F*)hist_CF1->Clone(HistName.Data());

  Int_t NBins = hist_CF_sum->GetNbinsX();

  for (Int_t i = 0; i < NBins; i++) {
    double CF1_val = hist_CF1->GetBinContent(i + 1);
    double CF1_err = hist_CF1->GetBinError(i + 1);
    double CF2_val = hist_CF2->GetBinContent(i + 1);
    double CF2_err = hist_CF2->GetBinError(i + 1);

    // average for bin i:
    if (CF1_val != 0. && CF2_val != 0.) {
      double CF1_err_weight = 1. / TMath::Power(CF1_err, 2.);
      double CF2_err_weight = 1. / TMath::Power(CF2_err, 2.);

      double CF_sum_average =
          (CF1_err_weight * CF1_val + CF2_err_weight * CF2_val) /
          (CF1_err_weight + CF2_err_weight);
      double CF_sum_err = 1. / TMath::Sqrt(CF1_err_weight + CF2_err_weight);

      hist_CF_sum->SetBinContent(i + 1, CF_sum_average);
      hist_CF_sum->SetBinError(i + 1, CF_sum_err);
    } else if (CF1_val == 0. && CF2_val != 0.) {
      hist_CF_sum->SetBinContent(i + 1, CF2_val);
      hist_CF_sum->SetBinError(i + 1, CF2_err);
    } else if (CF2_val == 0 && CF1_val != 0.) {
      hist_CF_sum->SetBinContent(i + 1, CF1_val);
      hist_CF_sum->SetBinError(i + 1, CF1_err);
    }
  }
  return hist_CF_sum;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors* DrawSystematicError(TH1F* histexp, TH1F* histerr,
                                  double errorwidth) {
  // Input are the experimental histogram with statistical errors only and a
  // histogram containing the errors

  const int histbins = histexp->GetNbinsX();

  TGraphErrors* ge_SysError_C2 = new TGraphErrors();

  for (int i = 0; i < histbins; i++) {
    if (histexp->GetBinCenter(i + 1) > 0.2) continue;
    ge_SysError_C2->SetPoint(i, histexp->GetBinCenter(i + 1),
                             histexp->GetBinContent(i + 1));
    ge_SysError_C2->SetPointError(i, errorwidth, histerr->GetBinContent(i + 1));
  }

  return ge_SysError_C2;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors* FemtoModelFitBands(TGraph* grMedian1, TGraph* grLower,
                                 TGraph* grUpper) {
  TGraphErrors* grFemtoModel = new TGraphErrors();
  grFemtoModel->SetName(grMedian1->GetName());
  double x, yM1, yLo, yUp;
  int count = 0;
  for (int i = 0; i < grMedian1->GetN(); ++i) {
    grMedian1->GetPoint(i, x, yM1);
    grLower->GetPoint(i, x, yLo);
    grUpper->GetPoint(i, x, yUp);
    std::vector<float> yAll;
    yAll.push_back(yM1);
    yAll.push_back(yLo);
    yAll.push_back(yUp);
    std::sort(yAll.begin(), yAll.end());
    grFemtoModel->SetPoint(count, x / 1000.f, (yAll[2] + yAll[0]) / 2.f);
    grFemtoModel->SetPointError(count++, 0,
                                (yAll[2] + yAll[0]) / 2.f - yAll[0]);
  }
  return grFemtoModel;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraph* convertInGev(TGraph* gr) {
  TGraph* grOut = new TGraph();
  grOut->SetName(gr->GetName());
  double x, y;
  for (int i = 0; i < gr->GetN(); ++i) {
    gr->GetPoint(i, x, y);
    grOut->SetPoint(i, x / 1000.f, y);
  }
  return grOut;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors* convertHistoInGev(TH1F* gr) {
  TGraphErrors* grOut = new TGraphErrors();
  grOut->SetName(gr->GetName());
  for (int i = 0; i < gr->GetNbinsX(); ++i) {
    grOut->SetPoint(i, gr->GetBinCenter(i) / 1000.f, gr->GetBinContent(i));
    grOut->SetPointError(i, 0, gr->GetBinError(i));
  }
  return grOut;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void compareCF(const char* expfile = "~/Results/LHC17p_fast/AnalysisResults.root",
            const char* oldfile = "", const char* CATSfile = "") {
  gStyle->SetCanvasPreferGL(1);
  const float right = 0.025;
  const float top = 0.025;

  const float spinningDepth = 10;

  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;

  // for Data
  const float rebinDataNew = 5.f;
  const float rebinDataOld = 1.f;

  const char* prefix = "MB";
  const char* prefixold = prefix;
  const char* addon = "";
  const char* addonold = addon;

  SetStyle();

  // new DATA
  TFile* _file0 = TFile::Open(expfile);
  TDirectoryFile* dirResults = (TDirectoryFile*)(_file0->FindObjectAny(
      Form("%sResults%s", prefix, addon)));
  TList* Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon), Results);
  TList* tmpFolder = (TList*)Results->FindObject("Particle0_Particle0");
  TH1F* histRE_relK_pp =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle0");
  TH1F* histME_relK_pp =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle0");
  tmpFolder = (TList*)Results->FindObject("Particle1_Particle1");
  TH1F* histRE_relK_ApAp =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle1");
  TH1F* histME_relK_ApAp =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle1");
  tmpFolder = (TList*)Results->FindObject("Particle0_Particle2");
  TH1F* histRE_relK_Lp =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle2");
  histRE_relK_Lp->Rebin(rebinDataNew);
  TH1F* histME_relK_Lp =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle2");
  histME_relK_Lp->Rebin(rebinDataNew);
  tmpFolder = (TList*)Results->FindObject("Particle1_Particle3");
  TH1F* histRE_relK_ALAp =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle3");
  histRE_relK_ALAp->Rebin(rebinDataNew);
  TH1F* histME_relK_ALAp =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle3");
  histME_relK_ALAp->Rebin(rebinDataNew);
  tmpFolder = (TList*)Results->FindObject("Particle2_Particle2");
  TH1F* histRE_relK_LL =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle2_Particle2");
  histRE_relK_LL->Rebin(rebinDataNew);
  TH1F* histME_relK_LL =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle2_Particle2");
  histME_relK_LL->Rebin(rebinDataNew);
  tmpFolder = (TList*)Results->FindObject("Particle3_Particle3");
  TH1F* histRE_relK_ALAL =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle3_Particle3");
  histRE_relK_ALAL->Rebin(rebinDataNew);
  TH1F* histME_relK_ALAL =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle3_Particle3");
  histME_relK_ALAL->Rebin(rebinDataNew);
  tmpFolder = (TList*)Results->FindObject("Particle0_Particle4");
  TH1F* histRE_relK_Xip =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle4");
  histRE_relK_Xip->Rebin(rebinDataNew);
  TH1F* histME_relK_Xip =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle4");
  histME_relK_Xip->Rebin(rebinDataNew);
  tmpFolder = (TList*)Results->FindObject("Particle1_Particle5");
  TH1F* histRE_relK_AXiAp =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle5");
  histRE_relK_AXiAp->Rebin(rebinDataNew);
  TH1F* histME_relK_AXiAp =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle5");
  histME_relK_AXiAp->Rebin(rebinDataNew);
  TH1F* hist_CF_Lp_ALAp_exp[3];
  TH1F* hist_CF_LL_ALAL_exp[3];
  TH1F* hist_CF_pp_ApAp_exp[3];
  TH1F* hist_CF_pXi_ApAXi_exp[3];

  hist_CF_Lp_ALAp_exp[0] =
      Calculate_CF(histRE_relK_Lp, histME_relK_Lp, "hist_CF_Lp_exp", normleft,
                   normright, addon, spinningDepth);
  hist_CF_Lp_ALAp_exp[1] =
      Calculate_CF(histRE_relK_ALAp, histME_relK_ALAp, "hist_CF_ALAp_exp",
                   normleft, normright, addon, spinningDepth);
  hist_CF_Lp_ALAp_exp[2] =
      add_CF(hist_CF_Lp_ALAp_exp[0], hist_CF_Lp_ALAp_exp[1],
             "hist_CF_Lp_ALAp_exp_sum");
  hist_CF_LL_ALAL_exp[0] =
      Calculate_CF(histRE_relK_LL, histME_relK_LL, "hist_CF_LL_exp", normleft,
                   normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[1] =
      Calculate_CF(histRE_relK_ALAL, histME_relK_ALAL, "hist_CF_LL_exp",
                   normleft, normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[2] =
      add_CF(hist_CF_LL_ALAL_exp[0], hist_CF_LL_ALAL_exp[1],
             "hist_CF_LL_ALAL_exp_sum");
  hist_CF_pp_ApAp_exp[0] =
      Calculate_CF(histRE_relK_pp, histME_relK_pp, "hist_CF_pp", normleft,
                   normright, addon, spinningDepth);
  hist_CF_pp_ApAp_exp[1] =
      Calculate_CF(histRE_relK_ApAp, histME_relK_ApAp, "hist_CF_ApAp", normleft,
                   normright, addon, spinningDepth);
  hist_CF_pp_ApAp_exp[2] =
      add_CF(hist_CF_pp_ApAp_exp[0], hist_CF_pp_ApAp_exp[1],
             "hist_CF_pp_ApAp_exp_sum");
  hist_CF_pXi_ApAXi_exp[0] =
      Calculate_CF(histRE_relK_Xip, histME_relK_Xip, "hist_CF_pXi", normleft,
                   normright, addon, spinningDepth);
  hist_CF_pXi_ApAXi_exp[1] =
      Calculate_CF(histRE_relK_AXiAp, histME_relK_AXiAp, "hist_CF_ApAXi",
                   normleft, normright, addon, spinningDepth);
  hist_CF_pXi_ApAXi_exp[2] =
      add_CF(hist_CF_pXi_ApAXi_exp[0], hist_CF_pXi_ApAXi_exp[1],
             "hist_CF_pXi_ApAXi_exp_sum");

  SetStyleHisto(hist_CF_pp_ApAp_exp[0], 1, 1);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[0], 1, 1);
  SetStyleHisto(hist_CF_LL_ALAL_exp[0], 1, 1);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[0], 1, 1);
  SetStyleHisto(hist_CF_pp_ApAp_exp[1], 0, 2);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[1], 0, 2);
  SetStyleHisto(hist_CF_LL_ALAL_exp[1], 0, 2);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[1], 0, 2);
  SetStyleHisto(hist_CF_pp_ApAp_exp[2], 0, 0);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[2], 0, 0);
  SetStyleHisto(hist_CF_LL_ALAL_exp[2], 0, 0);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[2], 0, 0);

  // old DATA
  TFile* _file0old = TFile::Open(oldfile);
  TH1F* hist_CF_Lp_ALAp_old[3];
  TH1F* hist_CF_LL_ALAL_old[3];
  TH1F* hist_CF_pp_ApAp_old[3];
  TH1F* hist_CF_pXi_ApAXi_old[3];
  TDirectoryFile* diroldResults = (TDirectoryFile*)(_file0old->FindObjectAny(
      Form("%sResults%s", prefixold, addonold)));
  TList* oldResults;
  diroldResults->GetObject(Form("%sResults%s", prefixold, addonold),
                           oldResults);
  tmpFolder = (TList*)oldResults->FindObject("Particle0_Particle0");
  TH1F* histRE_relK_ppold =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle0");
  TH1F* histME_relK_ppold =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle0");
  tmpFolder = (TList*)oldResults->FindObject("Particle1_Particle1");
  TH1F* histRE_relK_ApApold =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle1");
  TH1F* histME_relK_ApApold =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle1");
  tmpFolder = (TList*)oldResults->FindObject("Particle0_Particle2");
  TH1F* histRE_relK_Lpold =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle2");
  histRE_relK_Lpold->Rebin(rebinDataOld);
  TH1F* histME_relK_Lpold =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle2");
  histME_relK_Lpold->Rebin(rebinDataOld);
  tmpFolder = (TList*)oldResults->FindObject("Particle1_Particle3");
  TH1F* histRE_relK_ALApold =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle3");
  histRE_relK_ALApold->Rebin(rebinDataOld);
  TH1F* histME_relK_ALApold =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle3");
  histME_relK_ALApold->Rebin(rebinDataOld);
  tmpFolder = (TList*)oldResults->FindObject("Particle2_Particle2");
  TH1F* histRE_relK_LLold =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle2_Particle2");
  histRE_relK_LLold->Rebin(rebinDataOld);
  TH1F* histME_relK_LLold =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle2_Particle2");
  histME_relK_LLold->Rebin(rebinDataOld);
  tmpFolder = (TList*)oldResults->FindObject("Particle3_Particle3");
  TH1F* histRE_relK_ALALold =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle3_Particle3");
  histRE_relK_ALALold->Rebin(rebinDataOld);
  TH1F* histME_relK_ALALold =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle3_Particle3");
  histME_relK_ALALold->Rebin(rebinDataOld);
  tmpFolder = (TList*)oldResults->FindObject("Particle0_Particle4");
  TH1F* histRE_relK_Xipold =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle4");
  histRE_relK_Xipold->Rebin(rebinDataOld);
  TH1F* histME_relK_Xipold =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle4");
  histME_relK_Xipold->Rebin(rebinDataOld);
  tmpFolder = (TList*)oldResults->FindObject("Particle1_Particle5");
  TH1F* histRE_relK_AXiApold =
      (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle5");
  histRE_relK_AXiApold->Rebin(rebinDataOld);
  TH1F* histME_relK_AXiApold =
      (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle5");
  histME_relK_AXiApold->Rebin(rebinDataOld);

  hist_CF_Lp_ALAp_old[0] =
      Calculate_CF(histRE_relK_Lpold, histME_relK_Lpold, "hist_CF_Lp_old",
                   normleft, normright, addonold, spinningDepth);
  hist_CF_Lp_ALAp_old[1] =
      Calculate_CF(histRE_relK_ALApold, histME_relK_ALApold, "hist_CF_ALAp_old",
                   normleft, normright, addonold, spinningDepth);
  hist_CF_Lp_ALAp_old[2] =
      add_CF(hist_CF_Lp_ALAp_old[0], hist_CF_Lp_ALAp_old[1],
             "hist_CF_Lp_ALAp_old_sum");
  hist_CF_LL_ALAL_old[0] =
      Calculate_CF(histRE_relK_LLold, histME_relK_LLold, "hist_CF_LL_old",
                   normleft, normright, addonold, spinningDepth);
  hist_CF_LL_ALAL_old[1] =
      Calculate_CF(histRE_relK_ALALold, histME_relK_ALALold, "hist_CF_LL_old",
                   normleft, normright, addonold, spinningDepth);
  hist_CF_LL_ALAL_old[2] =
      add_CF(hist_CF_LL_ALAL_old[0], hist_CF_LL_ALAL_old[1],
             "hist_CF_LL_ALAL_old_sum");
  hist_CF_pp_ApAp_old[0] =
      Calculate_CF(histRE_relK_ppold, histME_relK_ppold, "hist_CF_ppold",
                   normleft, normright, addonold, spinningDepth);
  hist_CF_pp_ApAp_old[1] =
      Calculate_CF(histRE_relK_ApApold, histME_relK_ApApold, "hist_CF_ApApold",
                   normleft, normright, addonold, spinningDepth);
  hist_CF_pp_ApAp_old[2] =
      add_CF(hist_CF_pp_ApAp_old[0], hist_CF_pp_ApAp_old[1],
             "hist_CF_pp_ApAp_old_sum");
  hist_CF_pXi_ApAXi_old[0] =
      Calculate_CF(histRE_relK_Xipold, histME_relK_Xipold, "hist_CF_pXiold",
                   normleft, normright, addonold, spinningDepth);
  hist_CF_pXi_ApAXi_old[1] = Calculate_CF(
      histRE_relK_AXiApold, histME_relK_AXiApold, "hist_CF_ApAXi_old", normleft,
      normright, addonold, spinningDepth);
  hist_CF_pXi_ApAXi_old[2] =
      add_CF(hist_CF_pXi_ApAXi_old[0], hist_CF_pXi_ApAXi_old[1],
             "hist_CF_pXi_ApAXi_old_sum");
  SetStyleHisto(hist_CF_pp_ApAp_old[2], 1, 2);
  SetStyleHisto(hist_CF_Lp_ALAp_old[2], 1, 2);
  SetStyleHisto(hist_CF_LL_ALAL_old[2], 1, 2);
  SetStyleHisto(hist_CF_pXi_ApAXi_old[2], 1, 2);

  auto hist_ratio_uncert_pp = RatioUncertainties(hist_CF_pp_ApAp_old[2], hist_CF_pp_ApAp_exp[2]);
  auto hist_ratio_uncert_pL = RatioUncertainties(hist_CF_Lp_ALAp_old[2], hist_CF_Lp_ALAp_exp[2]);
  auto hist_ratio_uncert_LL = RatioUncertainties(hist_CF_LL_ALAL_old[2], hist_CF_LL_ALAL_exp[2]);
  auto hist_ratio_uncert_pXi = RatioUncertainties(hist_CF_pXi_ApAXi_old[2], hist_CF_pXi_ApAXi_exp[2]);
  SetStyleHisto(hist_ratio_uncert_pp, 0, 0);
  SetStyleHisto(hist_ratio_uncert_pL, 0, 0);
  SetStyleHisto(hist_ratio_uncert_LL, 0, 0);
  SetStyleHisto(hist_ratio_uncert_pXi, 0, 0);

  // PLOTTING
  TCanvas* Can_CF = new TCanvas("Can_CF", "Can_CF", 0, 0, 2000, 1100);
  Can_CF->Divide(4, 2);
  Can_CF->cd(1);
  Can_CF->cd(1)->SetRightMargin(right);
  Can_CF->cd(1)->SetTopMargin(top);
  hist_CF_pp_ApAp_exp[2]->Draw("pe");
  hist_CF_pp_ApAp_old[2]->Draw("pe same");
  hist_CF_pp_ApAp_exp[2]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_pp_ApAp_exp[2]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_pp_ApAp_exp[2]->GetYaxis()->SetRangeUser(0, 4);
  auto* leg = new TLegend(0.38, 0.7, 0.85, 0.85);
  leg->AddEntry(hist_CF_pp_ApAp_exp[2], "New", "pe");
  leg->AddEntry(hist_CF_pp_ApAp_old[2], "ALICE preliminary", "pe");
  leg->Draw("same");
  Can_CF->cd(2);
  Can_CF->cd(2)->SetRightMargin(right);
  Can_CF->cd(2)->SetTopMargin(top);
  hist_CF_Lp_ALAp_exp[2]->Draw("pe");
  hist_CF_Lp_ALAp_old[2]->Draw("pe same");
  hist_CF_Lp_ALAp_exp[2]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_Lp_ALAp_exp[2]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_Lp_ALAp_exp[2]->GetXaxis()->SetNdivisions(505);
  hist_CF_Lp_ALAp_exp[2]->GetYaxis()->SetRangeUser(0.8, 2.5);
  Can_CF->cd(3);
  Can_CF->cd(3)->SetRightMargin(right);
  Can_CF->cd(3)->SetTopMargin(top);
  hist_CF_LL_ALAL_exp[2]->Draw("pe");
  hist_CF_LL_ALAL_old[2]->Draw("pe same");
  hist_CF_LL_ALAL_exp[2]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_LL_ALAL_exp[2]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_LL_ALAL_exp[2]->GetXaxis()->SetNdivisions(505);
  hist_CF_LL_ALAL_exp[2]->GetYaxis()->SetRangeUser(0.25, 3);
  Can_CF->cd(4);
  Can_CF->cd(4)->SetRightMargin(right);
  Can_CF->cd(4)->SetTopMargin(top);
  hist_CF_pXi_ApAXi_exp[2]->Draw("pe");
  hist_CF_pXi_ApAXi_old[2]->Draw("pe same");
  hist_CF_pXi_ApAXi_exp[2]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_pXi_ApAXi_exp[2]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_pXi_ApAXi_exp[2]->GetXaxis()->SetNdivisions(505);
  hist_CF_pXi_ApAXi_exp[2]->GetYaxis()->SetRangeUser(0.25, 3);
  Can_CF->cd(5);
  Can_CF->cd(5)->SetRightMargin(right);
  Can_CF->cd(5)->SetTopMargin(top);
  TH1F* ppRatio = (TH1F*)hist_CF_pp_ApAp_exp[2]->Clone();
  auto* line = new TLine(0, 1, 0.4, 1);
  line->SetLineColor(kGray + 3);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  ppRatio->Divide(hist_CF_pp_ApAp_old[2]);
  ppRatio->GetYaxis()->SetRangeUser(0, 2);
  ppRatio->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)_{new}/#it{C}(k*)_{old}");
  ppRatio->SetMarkerStyle(fMarkers[2]);
  ppRatio->Draw("pe");
  line->Draw("same");
  Can_CF->cd(6);
  Can_CF->cd(6)->SetRightMargin(right);
  Can_CF->cd(6)->SetTopMargin(top);
  TH1F* pLRatio = (TH1F*)hist_CF_Lp_ALAp_exp[2]->Clone();
  pLRatio->Divide(hist_CF_Lp_ALAp_old[2]);
  pLRatio->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)_{new}/#it{C}(k*)_{old}");
  pLRatio->GetYaxis()->SetRangeUser(0, 2);
  pLRatio->SetMarkerStyle(fMarkers[2]);
  pLRatio->Draw("pe");
  line->Draw("same");
  Can_CF->cd(7);
  Can_CF->cd(7)->SetRightMargin(right);
  Can_CF->cd(7)->SetTopMargin(top);
  TH1F* LLRatio = (TH1F*)hist_CF_LL_ALAL_exp[2]->Clone();
  LLRatio->Divide(hist_CF_LL_ALAL_old[2]);
  LLRatio->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)_{new}/#it{C}(k*)_{old}");
  LLRatio->GetYaxis()->SetRangeUser(0, 2);
  LLRatio->SetMarkerStyle(fMarkers[2]);
  LLRatio->Draw("pe");
  line->Draw("same");
  Can_CF->cd(8);
  Can_CF->cd(8)->SetRightMargin(right);
  Can_CF->cd(8)->SetTopMargin(top);
  TH1F* pXiRatio = (TH1F*)hist_CF_pXi_ApAXi_exp[2]->Clone();
  pXiRatio->Divide(hist_CF_pXi_ApAXi_old[2]);
  pXiRatio->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)_{new}/#it{C}(k*)_{old}");
  pXiRatio->GetYaxis()->SetRangeUser(0, 2);
  pXiRatio->SetMarkerStyle(fMarkers[2]);
  pXiRatio->Draw("pe");
  line->Draw("same");
  Can_CF->Print("ANplot/CF_old-new.pdf");

  TCanvas* Can_Uncertainties = new TCanvas("Can_Uncertainties", "Can_Uncertainties", 0, 0, 2000, 1100);
  Can_Uncertainties->Divide(4,1);
  Can_Uncertainties->cd(1);
  Can_Uncertainties->cd(1)->SetRightMargin(right);
  Can_Uncertainties->cd(1)->SetTopMargin(top);
  hist_ratio_uncert_pp->Draw("pe");
  hist_ratio_uncert_pp->GetXaxis()->SetRangeUser(0, 0.4);
  Can_Uncertainties->cd(2);
  Can_Uncertainties->cd(2)->SetRightMargin(right);
  Can_Uncertainties->cd(2)->SetTopMargin(top);
  hist_ratio_uncert_pL->Draw("pe");
  hist_ratio_uncert_pL->GetXaxis()->SetRangeUser(0, 0.4);
  Can_Uncertainties->cd(3);
  Can_Uncertainties->cd(3)->SetRightMargin(right);
  Can_Uncertainties->cd(3)->SetTopMargin(top);
  hist_ratio_uncert_LL->Draw("pe");
  hist_ratio_uncert_LL->GetXaxis()->SetRangeUser(0, 0.4);
  Can_Uncertainties->cd(4);
  Can_Uncertainties->cd(4)->SetRightMargin(right);
  Can_Uncertainties->cd(4)->SetTopMargin(top);
  hist_ratio_uncert_pXi->Draw("pe");
  hist_ratio_uncert_pXi->GetXaxis()->SetRangeUser(0, 0.4);

}
