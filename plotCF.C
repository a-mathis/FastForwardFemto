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
  histo->SetMarkerSize(1.5);
  histo->SetLineWidth(2.5);
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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void plotCF(const char *expfile = "~/Results/LHC17p_fast/AnalysisResults.root", const char *simfile = "", const char *CATSfile = "")
{
  gStyle->SetCanvasPreferGL(1);
  const float right = 0.025;
  const float top = 0.025;

  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;
  const float spinningDepth = 10.f;

  // for pythia comparison
  const float rebin = 2;

  const int energy = 5; // TeV
  const char *system = "p-Pb";
//  const int energy = 13; // TeV
//  const char *system = "pp";

  bool EPOS = false;

  const char *prefix = "MB";
  const char *prefixSim = prefix;
  const char *addon = "";
  const char *addonSim = addon;

  TString data = "Data";
  TString sim = "DPMJET";
//  TString sim = "Pythia 8";

  float r = 1.437;
  float rErr = 0.011;
  float rSystErrUp = 0.013;
  float rSystErrDown = 0.006;
//  float r = 1.188;
//  float rErr = 0.009;
//  float rSystErrUp = 0.016;
//  float rSystErrDown = 0.009;

  float ppBL0 = 0.924;
  float ppBL1 = 0.386;
  float pLBL0 = 0.898;
  float pLBL1 = 0.355;
  float LLBL0 = 0.894;
  float LLBL1 = 0.320;
  float pXiBL0 = 0.919;
  float pXiBL1 = 0.453;

  //  EPOS
  if(EPOS) {
    r = 1.476;
    rErr = 0.015;
    rSystErrUp = 0.007;
    rSystErrDown = 0.014;

    ppBL0 = 0.931;
    ppBL1 = 0.399;
    pLBL0 = 0.938;
    pLBL1 = 0.328;
    LLBL0 = 0.893;
    LLBL1 = 0.321;
    pXiBL0 = 0.921;
    pXiBL1 = 0.443;
  }
  SetStyle();

  // EXP DATA
  TFile* _file0=TFile::Open(expfile);
  TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon),Results);
  TList* tmpFolder=(TList*)Results->FindObject("Particle0_Particle0");
  TH1F* histRE_relK_pp = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle0");
  TH1F* histME_relK_pp = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle0");
  tmpFolder=(TList*)Results->FindObject("Particle1_Particle1");
  TH1F* histRE_relK_ApAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle1");
  TH1F* histME_relK_ApAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle1");
  tmpFolder=(TList*)Results->FindObject("Particle0_Particle2");
  TH1F* histRE_relK_Lp = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle2");
  TH1F* histME_relK_Lp = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle2");
  tmpFolder=(TList*)Results->FindObject("Particle1_Particle3");
  TH1F* histRE_relK_ALAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle3");
  TH1F* histME_relK_ALAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle3");
  tmpFolder=(TList*)Results->FindObject("Particle2_Particle2");
  TH1F* histRE_relK_LL = (TH1F*)tmpFolder->FindObject("SEDist_Particle2_Particle2");
  TH1F* histME_relK_LL = (TH1F*)tmpFolder->FindObject("MEDist_Particle2_Particle2");
  tmpFolder=(TList*)Results->FindObject("Particle3_Particle3");
  TH1F* histRE_relK_ALAL = (TH1F*)tmpFolder->FindObject("SEDist_Particle3_Particle3");
  TH1F* histME_relK_ALAL = (TH1F*)tmpFolder->FindObject("MEDist_Particle3_Particle3");
  tmpFolder=(TList*)Results->FindObject("Particle0_Particle4");
  TH1F* histRE_relK_Xip = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle4");
  TH1F* histME_relK_Xip = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle4");
  tmpFolder=(TList*)Results->FindObject("Particle1_Particle5");
  TH1F* histRE_relK_AXiAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle5");
  TH1F* histME_relK_AXiAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle5");
  TH1F *hist_CF_Lp_ALAp_exp[3];
  TH1F *hist_CF_LL_ALAL_exp[3];
  TH1F *hist_CF_pp_ApAp_exp[3];
  TH1F *hist_CF_pXi_ApAXi_exp[3];

  hist_CF_Lp_ALAp_exp[0] = Calculate_CF(histRE_relK_Lp,histME_relK_Lp,"hist_CF_Lp_exp",normleft,normright, addon, spinningDepth);
  hist_CF_Lp_ALAp_exp[1] = Calculate_CF(histRE_relK_ALAp,histME_relK_ALAp,"hist_CF_ALAp_exp",normleft,normright, addon, spinningDepth);
  hist_CF_Lp_ALAp_exp[2] = add_CF(hist_CF_Lp_ALAp_exp[0],hist_CF_Lp_ALAp_exp[1],"hist_CF_Lp_ALAp_exp_sum");
  hist_CF_LL_ALAL_exp[0] = Calculate_CF(histRE_relK_LL,histME_relK_LL,"hist_CF_LL_exp",normleft,normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[1] = Calculate_CF(histRE_relK_ALAL,histME_relK_ALAL,"hist_CF_LL_exp",normleft,normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[2] = add_CF(hist_CF_LL_ALAL_exp[0],hist_CF_LL_ALAL_exp[1],"hist_CF_LL_ALAL_exp_sum");
  hist_CF_pp_ApAp_exp[0] = Calculate_CF(histRE_relK_pp,histME_relK_pp,"hist_CF_pp",normleft,normright, addon, spinningDepth);
  hist_CF_pp_ApAp_exp[1] = Calculate_CF(histRE_relK_ApAp,histME_relK_ApAp,"hist_CF_ApAp",normleft,normright, addon, spinningDepth);
  hist_CF_pp_ApAp_exp[2] = add_CF(hist_CF_pp_ApAp_exp[0],hist_CF_pp_ApAp_exp[1],"hist_CF_pp_ApAp_exp_sum");
  hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(histRE_relK_Xip,histME_relK_Xip,"hist_CF_pXi",normleft,normright, addon, spinningDepth);
  hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(histRE_relK_AXiAp,histME_relK_AXiAp,"hist_CF_ApAXi",normleft,normright, addon, spinningDepth);
  hist_CF_pXi_ApAXi_exp[2] = add_CF(hist_CF_pXi_ApAXi_exp[0],hist_CF_pXi_ApAXi_exp[1],"hist_CF_pXi_ApAXi_exp_sum");

  SetStyleHisto(hist_CF_pp_ApAp_exp[0], 1,1);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[0], 1,1);
  SetStyleHisto(hist_CF_LL_ALAL_exp[0], 1,1);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[0], 1,1);
  SetStyleHisto(hist_CF_pp_ApAp_exp[1], 0,2);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[1], 0,2);
  SetStyleHisto(hist_CF_LL_ALAL_exp[1], 0,2);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[1], 0,2);
  SetStyleHisto(hist_CF_pp_ApAp_exp[2], 0,0);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[2], 0,0);
  SetStyleHisto(hist_CF_LL_ALAL_exp[2], 0,0);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[2], 0,0);

  // SYSTEMATIC UNCERTAINTIES
  TString input_sys_pp = "C2totalsysPP.root";
  TString input_sys_pL = "C2totalsysPL.root";
  TString input_sys_LL = "C2totalsysLL.root";
  TString input_sys_pXi = "C2totalsysPXi.root";
  TFile* file_sys_pp = new TFile(input_sys_pp.Data());
  TFile* file_sys_pL = new TFile(input_sys_pL.Data());
  TFile* file_sys_pXi = new TFile(input_sys_pXi.Data());
  TFile* file_sys_LL = new TFile(input_sys_LL.Data());
  TH1F* hist_sys_pp = (TH1F*)file_sys_pp->Get("C2totalsysPP");
  TH1F* hist_sys_pL = (TH1F*)file_sys_pL->Get("C2totalsysPL");
  TH1F* hist_sys_LL = (TH1F*)file_sys_LL->Get("C2totalsysLL");
  TH1F* hist_sys_pXi = (TH1F*)file_sys_pXi->Get("C2totalsysPXi");

  // Sim DATA
  TFile* _file0sim=TFile::Open(simfile);
  TH1F *hist_CF_Lp_ALAp_sim[3];
  TH1F *hist_CF_LL_ALAL_sim[3];
  TH1F *hist_CF_pp_ApAp_sim[3];
  TH1F *hist_CF_pXi_ApAXi_sim[3];
  if(_file0sim) {
    TDirectoryFile *dirSimResults=(TDirectoryFile*)(_file0sim->FindObjectAny(Form("%sResults%s", prefixSim, addonSim)));
    TList *SimResults;
    dirSimResults->GetObject(Form("%sResults%s", prefixSim, addonSim),SimResults);
    tmpFolder=(TList*)SimResults->FindObject("Particle0_Particle0");
    TH1F* histRE_relK_ppsim = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle0");
    TH1F* histME_relK_ppsim = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle0");
    tmpFolder=(TList*)SimResults->FindObject("Particle1_Particle1");
    TH1F* histRE_relK_ApApsim = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle1");
    TH1F* histME_relK_ApApsim = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle1");
    tmpFolder=(TList*)SimResults->FindObject("Particle0_Particle2");
    TH1F* histRE_relK_Lpsim = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle2");
    TH1F* histME_relK_Lpsim = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle2");
    tmpFolder=(TList*)SimResults->FindObject("Particle1_Particle3");
    TH1F* histRE_relK_ALApsim = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle3");
    TH1F* histME_relK_ALApsim = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle3");
    tmpFolder=(TList*)SimResults->FindObject("Particle2_Particle2");
    TH1F* histRE_relK_LLsim = (TH1F*)tmpFolder->FindObject("SEDist_Particle2_Particle2");
    TH1F* histME_relK_LLsim = (TH1F*)tmpFolder->FindObject("MEDist_Particle2_Particle2");
    tmpFolder=(TList*)SimResults->FindObject("Particle3_Particle3");
    TH1F* histRE_relK_ALALsim = (TH1F*)tmpFolder->FindObject("SEDist_Particle3_Particle3");
    TH1F* histME_relK_ALALsim = (TH1F*)tmpFolder->FindObject("MEDist_Particle3_Particle3");
    tmpFolder=(TList*)SimResults->FindObject("Particle0_Particle4");
    TH1F* histRE_relK_Xipsim = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle4");
    TH1F* histME_relK_Xipsim = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle4");
    tmpFolder=(TList*)SimResults->FindObject("Particle1_Particle5");
    TH1F* histRE_relK_AXiApsim = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle5");
    TH1F* histME_relK_AXiApsim = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle5");

    hist_CF_Lp_ALAp_sim[0] = Calculate_CF(histRE_relK_Lpsim,histME_relK_Lpsim,"hist_CF_Lp_sim",normleft,normright, addonSim, spinningDepth);
    hist_CF_Lp_ALAp_sim[1] = Calculate_CF(histRE_relK_ALApsim,histME_relK_ALApsim,"hist_CF_ALAp_sim",normleft,normright, addonSim, spinningDepth);
    hist_CF_Lp_ALAp_sim[2] = add_CF(hist_CF_Lp_ALAp_sim[0],hist_CF_Lp_ALAp_sim[1],"hist_CF_Lp_ALAp_sim_sum");
    hist_CF_LL_ALAL_sim[0] = Calculate_CF(histRE_relK_LLsim,histME_relK_LLsim,"hist_CF_LL_sim",normleft,normright, addonSim, spinningDepth);
    hist_CF_LL_ALAL_sim[1] = Calculate_CF(histRE_relK_ALALsim,histME_relK_ALALsim,"hist_CF_LL_sim",normleft,normright, addonSim, spinningDepth);
    hist_CF_LL_ALAL_sim[2] = add_CF(hist_CF_LL_ALAL_sim[0],hist_CF_LL_ALAL_sim[1],"hist_CF_LL_ALAL_sim_sum");
    hist_CF_pp_ApAp_sim[0] = Calculate_CF(histRE_relK_ppsim,histME_relK_ppsim,"hist_CF_ppsim",normleft,normright, addonSim, spinningDepth);
    hist_CF_pp_ApAp_sim[1] = Calculate_CF(histRE_relK_ApApsim,histME_relK_ApApsim,"hist_CF_ApApsim",normleft,normright, addonSim, spinningDepth);
    hist_CF_pp_ApAp_sim[2] = add_CF(hist_CF_pp_ApAp_sim[0],hist_CF_pp_ApAp_sim[1],"hist_CF_pp_ApAp_sim_sum");
    hist_CF_pXi_ApAXi_sim[0] = Calculate_CF(histRE_relK_Xipsim,histME_relK_Xipsim,"hist_CF_pXisim",normleft,normright, addonSim, spinningDepth);
    hist_CF_pXi_ApAXi_sim[1] = Calculate_CF(histRE_relK_AXiApsim,histME_relK_AXiApsim,"hist_CF_ApAXisim",normleft,normright, addonSim, spinningDepth);
    hist_CF_pXi_ApAXi_sim[2] = add_CF(hist_CF_pXi_ApAXi_sim[0],hist_CF_pXi_ApAXi_sim[1],"hist_CF_pXi_ApAXi_sim_sum");
    SetStyleHisto(hist_CF_pp_ApAp_sim[2], 1,2);
    SetStyleHisto(hist_CF_Lp_ALAp_sim[2], 1,2);
    SetStyleHisto(hist_CF_LL_ALAL_sim[2], 1,2);
    SetStyleHisto(hist_CF_pXi_ApAXi_sim[2], 1,2);
  }

  // UNCERTAINTIES ON THE FIT
  TFile *systFit = TFile::Open(CATSfile);
  TGraph *grppDefault = nullptr;
  TGraph *grppLow = nullptr;
  TGraph *grppUp = nullptr;
  TGraph *grpLDefaultNLO = nullptr;
  TGraph *grpLLowNLO = nullptr;
  TGraph *grpLUpNLO = nullptr;
  TGraph *grpLDefaultLO = nullptr;
  TGraph *grpLLowLO = nullptr;
  TGraph *grpLUpLO = nullptr;
  TGraph *grLLDefault = nullptr;
  TGraph *grLLLow = nullptr;
  TGraph *grLLUp = nullptr;
  TGraph *grpXiDefault = nullptr;
  TGraph *grpXiLower = nullptr;
  TGraph *grpXiUpper = nullptr;
  TGraph *grpXiDefaultCoulomb = nullptr;
  TGraph *grpXiLowerCoulomb = nullptr;
  TGraph *grpXiUpperCoulomb = nullptr;
  TGraphErrors *grFemtopp = nullptr;
  TGraphErrors *grFemtopLNLO = nullptr;
  TGraphErrors *grFemtopLLO = nullptr;
  TGraphErrors *grFemtoLL = nullptr;
  TGraphErrors *grFemtopXi = nullptr;
  TGraphErrors *grFemtopXiCoulomb = nullptr;
  if(systFit) {
    grppDefault = (TGraph*)systFit->Get("ppGraphDefault");
    grppLow = (TGraph*)systFit->Get("ppGraphLowerLim");
    grppUp = (TGraph*)systFit->Get("ppGraphUpperLim");

    const char* nlostring = (EPOS) ? "" : "_NLO";
    const char* lostring = (EPOS) ? "_LO" : "";
    const char* prefix = (EPOS) ? "" : "Copy_";
    const char* prefixLO = (!EPOS) ? "" : "Copy_";
    grpLDefaultNLO = (TGraph*)systFit->Get(Form("pLamGraphDefault%s", nlostring));
    grpLLowNLO = (TGraph*)systFit->Get(Form("%spLamGraphLowerLim%s", prefix, nlostring));
    grpLUpNLO = (TGraph*)systFit->Get(Form("pLamGraphUpperLim%s", nlostring));
    grpLDefaultLO = (TGraph*)systFit->Get(Form("pLamGraphDefault%s", lostring));
    grpLLowLO = (TGraph*)systFit->Get(Form("%spLamGraphLowerLim%s", prefixLO, lostring));
    grpLUpLO = (TGraph*)systFit->Get(Form("pLamGraphUpperLim%s", lostring));

    grLLDefault = (TGraph*)systFit->Get("LamLamGraphDefault");
    grLLLow = (TGraph*)systFit->Get("LamLamGraphLowerLim");
    grLLUp = (TGraph*)systFit->Get("LamLamGraphUpperLim");
    grpXiDefault = (TGraph*)systFit->Get("pXimGraphDefault");
    grpXiLower = (TGraph*)systFit->Get("pXimGraphLowerLim");
    grpXiUpper = (TGraph*)systFit->Get("pXimGraphUpperLim");

    grpXiDefaultCoulomb = (TGraph*)systFit->Get("pXimGraphDefault_COULOMB");
    grpXiLowerCoulomb = (TGraph*)systFit->Get("pXimGraphLowerLim_COULOMB");
    grpXiUpperCoulomb = (TGraph*)systFit->Get("pXimGraphUpperLim_COULOMB");

    grFemtopp = FemtoModelFitBands(grppDefault, grppLow, grppUp);
    grFemtopp->SetFillColor(fColors[2]);
    grFemtopp->SetLineColor(fColors[2]);
    grFemtopp->SetLineWidth(3);
    grFemtopLNLO = FemtoModelFitBands(grpLDefaultNLO, grpLLowNLO, grpLUpNLO);
    grFemtopLNLO->SetFillColor(fColors[1]);
    grFemtopLNLO->SetLineColor(fColors[1]);
    grFemtopLNLO->SetLineWidth(3);
    grFemtopLLO = FemtoModelFitBands(grpLDefaultLO, grpLLowLO, grpLUpLO);
    grFemtopLLO->SetFillColor(fColors[3]);
    grFemtopLLO->SetLineColor(fColors[3]);
    grFemtopLLO->SetLineWidth(3);
    grFemtoLL = FemtoModelFitBands(grLLDefault, grLLLow, grLLUp);
    grFemtoLL->SetFillColor(fColors[5]);
    grFemtoLL->SetLineColor(fColors[5]);
    grFemtoLL->SetLineWidth(3);
    grFemtopXi = FemtoModelFitBands(grpXiDefault, grpXiLower, grpXiUpper);
    grFemtopXi->SetFillColor(fColors[6]);
    grFemtopXi->SetLineColor(fColors[6]);
    grFemtopXi->SetLineWidth(3);
    grFemtopXiCoulomb = FemtoModelFitBands(grpXiDefaultCoulomb, grpXiLowerCoulomb, grpXiUpperCoulomb);
    grFemtopXiCoulomb->SetFillColor(fColors[7]);
    grFemtopXiCoulomb->SetLineColor(fColors[7]);
    grFemtopXiCoulomb->SetLineWidth(3);
  }

  TGraph *grFakeSyst = new TGraph();
  grFakeSyst->SetFillColor(fFillColors[0]);
  grFakeSyst->SetLineColor(fFillColors[0]);

  TGraph *grFakePP = new TGraph();
  grFakePP->SetLineColor(fColors[7]);
  grFakePP->SetLineWidth(4);
  TGraph *grFakePLnlo = new TGraph();
  grFakePLnlo->SetLineColor(fColors[1]);
  grFakePLnlo->SetLineWidth(4);
  TGraph *grFakePLlo = new TGraph();
  grFakePLlo->SetLineColor(fColors[3]);
  grFakePLlo->SetLineWidth(4);
  TGraph *grFakeLL = new TGraph();
  grFakeLL->SetLineColor(fColors[5]);
  grFakeLL->SetLineWidth(4);
  TGraph *grFakeLLStar = new TGraph();
  grFakeLLStar->SetLineColor(fColors[6]);
  grFakeLLStar->SetLineWidth(4);
  TGraph *grFakepXi= new TGraph();
  grFakepXi->SetLineColor(fColors[6]);
  grFakepXi->SetLineWidth(4);
  TGraph *grFakeXiCoulomb= new TGraph();
  grFakeXiCoulomb->SetLineColor(fColors[7]);
  grFakeXiCoulomb->SetLineWidth(4);

  // BASELINE
  TF1 *baselinePP = new TF1("baselinePP", "pol1", 0, 1);
  TF1 *baselinePL = new TF1("baselinePL", "pol1", 0, 1);
  TF1 *baselineLL = new TF1("baselineLL", "pol1", 0, 1);
  TF1 *baselinePXI = new TF1("baselinePXI", "pol1", 0, 1);
  baselinePP->SetParameter(0, ppBL0);
  baselinePP->SetParameter(1, ppBL1);
  baselinePL->SetParameter(0, pLBL0);
  baselinePL->SetParameter(1, pLBL1);
  baselineLL->SetParameter(0, LLBL0);
  baselineLL->SetParameter(1, LLBL1);
  baselinePXI->SetParameter(0, pXiBL0);
  baselinePXI->SetParameter(1, pXiBL1);
  baselinePP->SetLineStyle(2);
  baselinePL->SetLineStyle(2);
  baselineLL->SetLineStyle(2);
  baselinePXI->SetLineStyle(2);
  baselinePP->SetLineColor(fFillColors[6]);
  baselinePL->SetLineColor(fFillColors[6]);
  baselineLL->SetLineColor(fFillColors[6]);
  baselinePXI->SetLineColor(fFillColors[6]);

  // PLOTTING
  TCanvas *Can_CF = new TCanvas("Can_CF","Can_CF",0,0,2000,1100);
  Can_CF->Divide(4,2);
  Can_CF->cd(1);
  Can_CF->cd(1)->SetRightMargin(right);
  Can_CF->cd(1)->SetTopMargin(top);
  hist_CF_pp_ApAp_exp[0]->Draw("pe");
  hist_CF_pp_ApAp_exp[1]->Draw("pe same");
  hist_CF_pp_ApAp_exp[0]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_pp_ApAp_exp[0]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_pp_ApAp_exp[0]->GetYaxis()->SetRangeUser(0, 4);
  auto* leg= new TLegend(0.25, 0.7, 0.95, 0.95);
  leg->AddEntry(hist_CF_pp_ApAp_exp[0], "Particle-particle CF", "pe");
  leg->AddEntry(hist_CF_pp_ApAp_exp[1], "Antiparticle-antiparticle CF", "pe");
  leg->Draw("same");
  Can_CF->cd(2);
  Can_CF->cd(2)->SetRightMargin(right);
  Can_CF->cd(2)->SetTopMargin(top);
  hist_CF_Lp_ALAp_exp[0]->Draw("pe");
  hist_CF_Lp_ALAp_exp[1]->Draw("pe same");
  hist_CF_Lp_ALAp_exp[0]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_Lp_ALAp_exp[0]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_Lp_ALAp_exp[0]->GetXaxis()->SetNdivisions(505);
  hist_CF_Lp_ALAp_exp[0]->GetYaxis()->SetRangeUser(0.8, 2.5);
  Can_CF->cd(3);
  Can_CF->cd(3)->SetRightMargin(right);
  Can_CF->cd(3)->SetTopMargin(top);
  hist_CF_LL_ALAL_exp[0]->Draw("pe");
  hist_CF_LL_ALAL_exp[1]->Draw("pe same");
  hist_CF_LL_ALAL_exp[0]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_LL_ALAL_exp[0]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_LL_ALAL_exp[0]->GetXaxis()->SetNdivisions(505);
  hist_CF_LL_ALAL_exp[0]->GetYaxis()->SetRangeUser(0.25, 3);
  Can_CF->cd(4);
  Can_CF->cd(4)->SetRightMargin(right);
  Can_CF->cd(4)->SetTopMargin(top);
  hist_CF_pXi_ApAXi_exp[0]->Draw("pe");
  hist_CF_pXi_ApAXi_exp[1]->Draw("pe same");
  hist_CF_pXi_ApAXi_exp[0]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_pXi_ApAXi_exp[0]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_pXi_ApAXi_exp[0]->GetXaxis()->SetNdivisions(505);
  hist_CF_pXi_ApAXi_exp[0]->GetYaxis()->SetRangeUser(0.25, 3);
  Can_CF->cd(5);
  Can_CF->cd(5)->SetRightMargin(right);
  Can_CF->cd(5)->SetTopMargin(top);
  TH1F *ppRatio = (TH1F*)hist_CF_pp_ApAp_exp[0]->Clone();
  auto* line = new TLine(0,1,0.4,1);
  line->SetLineColor(kGray+3);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  ppRatio->Divide(hist_CF_pp_ApAp_exp[1]);
  ppRatio->GetYaxis()->SetRangeUser(0,2);
  ppRatio->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)_{pp}/#it{C}(k*)_{#bar{pp}}");
  ppRatio->SetMarkerStyle(fMarkers[2]);
  ppRatio->Draw("pe");
  line->Draw("same");
  Can_CF->cd(6);
  Can_CF->cd(6)->SetRightMargin(right);
  Can_CF->cd(6)->SetTopMargin(top);
  TH1F *pLRatio = (TH1F*)hist_CF_Lp_ALAp_exp[0]->Clone();
  pLRatio->Divide(hist_CF_Lp_ALAp_exp[1]);
  pLRatio->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)_{p#Lambda}/#it{C}(k*)_{#bar{p}#bar{#Lambda}}");
  pLRatio->GetYaxis()->SetRangeUser(0,2);
  pLRatio->SetMarkerStyle(fMarkers[2]);
  pLRatio->Draw("pe");
  line->Draw("same");
  Can_CF->cd(7);
  Can_CF->cd(7)->SetRightMargin(right);
  Can_CF->cd(7)->SetTopMargin(top);
  TH1F *LLRatio = (TH1F*)hist_CF_LL_ALAL_exp[0]->Clone();
  LLRatio->Divide(hist_CF_LL_ALAL_exp[1]);
  LLRatio->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)_{#Lambda#Lambda}/#it{C}(k*)_{#bar{#Lambda}#bar{#Lambda}}");
  LLRatio->GetYaxis()->SetRangeUser(0,2);
  LLRatio->SetMarkerStyle(fMarkers[2]);
  LLRatio->Draw("pe");
  line->Draw("same");
  Can_CF->cd(8);
  Can_CF->cd(8)->SetRightMargin(right);
  Can_CF->cd(8)->SetTopMargin(top);
  TH1F *pXiRatio = (TH1F*)hist_CF_pXi_ApAXi_exp[0]->Clone();
  pXiRatio->Divide(hist_CF_pXi_ApAXi_exp[1]);
  pXiRatio->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)_{p#Xi^{-}}/#it{C}(k*)_{#bar{p}#Xi^{+}}");
  pXiRatio->GetYaxis()->SetRangeUser(0,2);
  pXiRatio->SetMarkerStyle(fMarkers[2]);
  pXiRatio->Draw("pe");
  line->Draw("same");
  Can_CF->Print("ANplot/CF_pp-apap.pdf");

  TCanvas *Can_CF_fitting = new TCanvas("Can_CF_fitting","Can_CF_fitting",0,0,1100,1000);
  Can_CF_fitting->Divide(2,2);
  Can_CF_fitting->cd(1);
  Can_CF_fitting->cd(1)->SetRightMargin(right);
  Can_CF_fitting->cd(1)->SetTopMargin(top);
  TGraphErrors *Tgraph_syserror_pp_ApAp = DrawSystematicError(hist_CF_pp_ApAp_exp[2], hist_sys_pp, 0.002);
  Tgraph_syserror_pp_ApAp->SetLineColor(kWhite);
  Tgraph_syserror_pp_ApAp->Draw("Ap");
  baselinePP->Draw("same");
  Tgraph_syserror_pp_ApAp->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pp_ApAp->GetXaxis()->SetRangeUser(0, 0.125);
  Tgraph_syserror_pp_ApAp->GetYaxis()->SetRangeUser(0.5, 3.5);
  if(grFemtopp) grFemtopp->Draw("L3 same");
  Tgraph_syserror_pp_ApAp->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_pp_ApAp->Draw("2 same");
  hist_CF_pp_ApAp_exp[2]->Draw("pe same");
  TLegend *legpp = new TLegend(0.48,0.595,0.85,0.81);
  legpp->SetBorderSize(0);
  legpp->SetTextFont(42);
  legpp->SetTextSize(gStyle->GetTextSize()*0.75);
  legpp->AddEntry(hist_CF_pp_ApAp_exp[2], "pp #oplus #bar{p}#bar{p} pairs", "pe");
  legpp->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legpp->AddEntry(grFakePP,"Femtoscopic fit","l");
  legpp->Draw("same");
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText.SetNDC(kTRUE);
  if(strcmp(system, "pp") == 0) BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}} = %i TeV", system, energy));
  else BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}_{NN}} = %i TeV", system, energy));
  TLatex text;
  text.SetTextSize(gStyle->GetTextSize()*0.75);
  text.SetNDC();
  text.SetTextColor(1);
  if(!EPOS) text.DrawLatex(0.5, 0.825, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.5, 0.825, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));

  Can_CF_fitting->cd(2);
  Can_CF_fitting->cd(2)->SetRightMargin(right);
  Can_CF_fitting->cd(2)->SetTopMargin(top);
  TGraphErrors *Tgraph_syserror_pL_ApAL = DrawSystematicError(hist_CF_Lp_ALAp_exp[2], hist_sys_pL, 0.01);
  Tgraph_syserror_pL_ApAL->SetLineColor(kWhite);
  Tgraph_syserror_pL_ApAL->Draw("ap");
  baselinePL->Draw("same");
  Tgraph_syserror_pL_ApAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_pL_ApAL->GetYaxis()->SetRangeUser(0.8, 2.);
  if(grFemtopLNLO) grFemtopLNLO->Draw("l3 same");
  if(!EPOS && grFemtopLLO) grFemtopLLO->Draw("l3 same");
  Tgraph_syserror_pL_ApAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_pL_ApAL->Draw("2 same");
  hist_CF_Lp_ALAp_exp[2]->Draw("pe same");
  TLegend *legLp;
  if(EPOS) legLp = new TLegend(0.385,0.595,0.75,0.81);
  else legLp     = new TLegend(0.385,0.545,0.75,0.81);
  legLp->SetBorderSize(0);
  legLp->SetTextFont(42);
  legLp->SetTextSize(gStyle->GetTextSize()*0.75);
  legLp->AddEntry(hist_CF_Lp_ALAp_exp[2], "p#Lambda #oplus #bar{p}#bar{#Lambda} pairs", "pe");
  legLp->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legLp->AddEntry(grFakePLnlo,"Femtoscopic fit (#chiEFT NLO)","l");
  if(!EPOS) legLp->AddEntry(grFakePLlo,"Femtoscopic fit (#chiEFT LO)","l");
  legLp->Draw("same");
  if(strcmp(system, "pp") == 0) BeamText.DrawLatex(0.4, 0.875, Form("ALICE %s #sqrt{#it{s}} = %i TeV", system, energy));
  else BeamText.DrawLatex(0.4, 0.875, Form("ALICE %s #sqrt{#it{s}_{NN}} = %i TeV", system, energy));
  if(!EPOS) text.DrawLatex(0.4, 0.825, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.4, 0.825, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  TLatex ref;
  ref.SetTextSize(gStyle->GetTextSize()*0.6);
  ref.SetNDC(kTRUE);
  if(EPOS) ref.DrawLatex(0.48, 0.565, "Nucl. Phys. A915 (2013) 24.");
  else ref.DrawLatex(0.48, 0.515, "Nucl. Phys. A915 (2013) 24.");

  Can_CF_fitting->cd(3);
  Can_CF_fitting->cd(3)->SetRightMargin(right);
  Can_CF_fitting->cd(3)->SetTopMargin(top);
  TGraphErrors *Tgraph_syserror_LL_ALAL = DrawSystematicError(hist_CF_LL_ALAL_exp[2], hist_sys_LL, 0.01);
  Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);
  Tgraph_syserror_LL_ALAL->Draw("ap");
  baselineLL->Draw("same");
  Tgraph_syserror_LL_ALAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_LL_ALAL->GetYaxis()->SetRangeUser(0.35, 2.);
  if(grFemtoLL) grFemtoLL->Draw("l3 same");
  Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_LL_ALAL->Draw("2 same");
  hist_CF_LL_ALAL_exp[2]->Draw("pe same");
  TLegend *legLL = new TLegend(0.48,0.595,0.85,0.81);
  legLL->SetBorderSize(0);
  legLL->SetTextFont(42);
  legLL->SetTextSize(gStyle->GetTextSize()*0.75);
  legLL->AddEntry(hist_CF_pp_ApAp_exp[2], "#Lambda#Lambda #oplus #bar{#Lambda}#bar{#Lambda} pairs", "pe");
  legLL->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legLL->AddEntry(grFakeLL,"Femtoscopic fit","l");
  legLL->Draw("same");
  if(strcmp(system, "pp") == 0) BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}} = %i TeV", system, energy));
  else BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}_{NN}} = %i TeV", system, energy));
  if(!EPOS) text.DrawLatex(0.5, 0.825, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.5, 0.825, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));

  Can_CF_fitting->cd(4);
  Can_CF_fitting->cd(4)->SetRightMargin(right);
  Can_CF_fitting->cd(4)->SetTopMargin(top);
  TGraphErrors *Tgraph_syserror_pXi_ApAXi = DrawSystematicError(hist_CF_pXi_ApAXi_exp[2], hist_sys_pXi, 0.01);
  Tgraph_syserror_pXi_ApAXi->SetLineColor(kWhite);
  Tgraph_syserror_pXi_ApAXi->Draw("ap");
  baselinePXI->Draw("same");
  Tgraph_syserror_pXi_ApAXi->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pXi_ApAXi->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_pXi_ApAXi->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_pXi_ApAXi->GetYaxis()->SetRangeUser(0.75, 3.);
  Tgraph_syserror_pXi_ApAXi->SetFillColorAlpha(kBlack, 0.4);
  if(grFemtopXi) grFemtopXi->Draw("l3 same");
  if(grFemtopXiCoulomb) grFemtopXiCoulomb->Draw("l3 same");
  Tgraph_syserror_pXi_ApAXi->Draw("2 same");
  hist_CF_pXi_ApAXi_exp[2]->Draw("pe same");
  TLegend *legpXi = new TLegend(0.385,0.545,0.75,0.81);
  legpXi->SetBorderSize(0);
  legpXi->SetTextFont(42);
  legpXi->SetTextSize(gStyle->GetTextSize()*0.75);
  legpXi->AddEntry(hist_CF_pp_ApAp_exp[2], "p#Xi^{-} #oplus #bar{p}#Xi^{+} pairs", "pe");
  legpXi->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legpXi->AddEntry(grFakeXiCoulomb, "Femtoscopic fit (Coulomb)","l");
  legpXi->AddEntry(grFakepXi,"Femtoscopic fit (HAL QCD)","l");
  legpXi->Draw("same");
  ref.DrawLatex(0.48, 0.515, "Nucl. Phys. A967 (2017) 856.");
  ref.DrawLatex(0.48, 0.475, "PoS LATTICE2016 (2017) 116.");
  if(strcmp(system, "pp") == 0) BeamText.DrawLatex(0.4, 0.875, Form("ALICE %s #sqrt{#it{s}} = %i TeV", system, energy));
  else BeamText.DrawLatex(0.4, 0.875, Form("ALICE %s #sqrt{#it{s}_{NN}} = %i TeV", system, energy));
  if(!EPOS) text.DrawLatex(0.4, 0.825, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.4, 0.825, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) Can_CF_fitting->Print("ANplot/CF_EPOS_prelim.pdf");
  else Can_CF_fitting->Print("ANplot/CF_Gauss_prelim.pdf");


  // Compare to Pythia
  if(_file0sim) {
    auto* grDummy1 = new TH1F("hDummy1", "; k* (GeV/#it{c}); #it{C}(k*)", 100, 0,2);
    auto* grDummy2 = new TH1F("hDummy2", "; k* (GeV/#it{c}); #it{C}(k*)", 100, 0,2);
    auto* grDummy3 = new TH1F("hDummy3", "; k* (GeV/#it{c}); #it{C}(k*)", 100, 0,2);
    auto* grDummy4 = new TH1F("hDummy4", "; k* (GeV/#it{c}); #it{C}(k*)", 100, 0,2);

    TCanvas *Can_CF_Pythia= new TCanvas("Can_CF_Pythia","Can_CF_Pythia",0,0,1000,1100);
    Can_CF_Pythia->Divide(2,2);
    Can_CF_Pythia->cd(1);
    Can_CF_Pythia->cd(1)->SetRightMargin(right);
    Can_CF_Pythia->cd(1)->SetTopMargin(top);
    grDummy1->Draw();
    grDummy1->GetYaxis()->SetRangeUser(0, 4);
    hist_CF_pp_ApAp_sim[2]->Rebin(rebin);
    hist_CF_pp_ApAp_sim[2]->Scale(1/rebin);
    hist_CF_pp_ApAp_sim[2]->Draw("pe same");
    auto *hist_pp_exp = (TH1F*)hist_CF_pp_ApAp_exp[2]->Clone("pp");
    hist_pp_exp->Rebin(rebin);
    hist_pp_exp->Scale(1/rebin);
    hist_pp_exp->Draw("pe same");
    auto *legpppy = new TLegend(0.65,0.695,0.85,0.91);
    legpppy->SetBorderSize(0);
    legpppy->SetTextFont(42);
    legpppy->SetTextSize(gStyle->GetTextSize()*0.75);
    legpppy->SetHeader("pp #oplus #bar{p}#bar{p}");
    legpppy->AddEntry(hist_pp_exp, data, "pe");
    legpppy->AddEntry(hist_CF_pp_ApAp_sim[2], sim, "pe");
    legpppy->Draw("same");

    Can_CF_Pythia->cd(2);
    Can_CF_Pythia->cd(2)->SetRightMargin(right);
    Can_CF_Pythia->cd(2)->SetTopMargin(top);
    grDummy2->Draw();
    grDummy2->GetYaxis()->SetRangeUser(0.25, 3);
    hist_CF_Lp_ALAp_sim[2]->Rebin(rebin);
    hist_CF_Lp_ALAp_sim[2]->Scale(1/rebin);
    hist_CF_Lp_ALAp_sim[2]->Draw("pe same");
    auto *hist_pL_exp = (TH1F*)hist_CF_Lp_ALAp_exp[2]->Clone("pL");
    hist_pL_exp->Rebin(rebin);
    hist_pL_exp->Scale(1/rebin);
    hist_pL_exp->Draw("pe same");
    auto *legLppy = new TLegend(0.65,0.695,0.85,0.91);
    legLppy->SetBorderSize(0);
    legLppy->SetTextFont(42);
    legLppy->SetTextSize(gStyle->GetTextSize()*0.75);
    legLppy->SetHeader("p#Lambda #oplus #bar{p}#bar{#Lambda}");
    legLppy->AddEntry(hist_pL_exp, data, "pe");
    legLppy->AddEntry(hist_CF_Lp_ALAp_sim[2], sim, "pe");
    legLppy->Draw("same");


    Can_CF_Pythia->cd(3);
    Can_CF_Pythia->cd(3)->SetRightMargin(right);
    Can_CF_Pythia->cd(3)->SetTopMargin(top);
    grDummy3->Draw();
    grDummy3->GetYaxis()->SetRangeUser(0.25, 3);
    hist_CF_LL_ALAL_sim[2]->Rebin(rebin);
    hist_CF_LL_ALAL_sim[2]->Scale(1/rebin);
    hist_CF_LL_ALAL_sim[2]->Draw("pe same");
    auto *hist_LL_exp = (TH1F*)hist_CF_LL_ALAL_exp[2]->Clone("LL");
    hist_LL_exp->Rebin(rebin);
    hist_LL_exp->Scale(1/rebin);
    hist_LL_exp->Draw("pe same");
    auto *legLLpy = new TLegend(0.65,0.695,0.85,0.91);
    legLLpy->SetBorderSize(0);
    legLLpy->SetTextFont(42);
    legLLpy->SetTextSize(gStyle->GetTextSize()*0.75);
    legLLpy->SetHeader("#Lambda#Lambda #oplus #bar{#Lambda}#bar{#Lambda}");
    legLLpy->AddEntry(hist_LL_exp, data, "pe");
    legLLpy->AddEntry(hist_CF_LL_ALAL_sim[2], sim, "pe");
    legLLpy->Draw("same");

    Can_CF_Pythia->cd(4);
    Can_CF_Pythia->cd(4)->SetRightMargin(right);
    Can_CF_Pythia->cd(4)->SetTopMargin(top);
    grDummy4->GetYaxis()->SetRangeUser(0.25, 3);
    grDummy4->Draw();
    hist_CF_pXi_ApAXi_sim[2]->Rebin(rebin);
    hist_CF_pXi_ApAXi_sim[2]->Scale(1/rebin);
    hist_CF_pXi_ApAXi_sim[2]->Draw("pe same");
    auto *hist_pXi_exp = (TH1F*)hist_CF_pXi_ApAXi_exp[2]->Clone("pXi");
    hist_pXi_exp->Rebin(rebin);
    hist_pXi_exp->Scale(1/rebin);
    hist_pXi_exp->Draw("pe same");
    auto *legpXipy = new TLegend(0.65,0.695,0.85,0.91);
    legpXipy->SetBorderSize(0);
    legpXipy->SetTextFont(42);
    legpXipy->SetTextSize(gStyle->GetTextSize()*0.75);
    legpXipy->SetHeader("p#Xi^{-} #oplus #bar{p}#Xi^{+}");
    legpXipy->AddEntry(hist_pXi_exp, data, "pe");
    legpXipy->AddEntry(hist_CF_pXi_ApAXi_sim[2], sim, "pe");
    legpXipy->Draw("same");

    Can_CF_Pythia->Print("ANplot/CF_pythia.pdf");
  }


  // PRELIMINARY PLOTS

  TCanvas *Can_CF_pp = new TCanvas("pp","pp", 0,0,650,550);
  Can_CF_pp->SetRightMargin(right);
  Can_CF_pp->SetTopMargin(top);
  Tgraph_syserror_pp_ApAp->SetLineColor(kWhite);
  Tgraph_syserror_pp_ApAp->Draw("Ap");
  baselinePP->Draw("same");
  Tgraph_syserror_pp_ApAp->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pp_ApAp->GetXaxis()->SetRangeUser(0, 0.125);
  Tgraph_syserror_pp_ApAp->GetYaxis()->SetRangeUser(0.5, 3.5);
  if(grFemtopp) grFemtopp->Draw("L3 same");
  Tgraph_syserror_pp_ApAp->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_pp_ApAp->Draw("2 same");
  hist_CF_pp_ApAp_exp[2]->Draw("pe same");
  TLegend *legpp2 = new TLegend(0.53,0.545,0.85,0.76);
  legpp2->SetBorderSize(0);
  legpp2->SetTextFont(42);
  legpp2->SetTextSize(gStyle->GetTextSize()*0.75);
  legpp2->AddEntry(hist_CF_pp_ApAp_exp[2], "pp #oplus #bar{p}#bar{p} pairs", "pe");
  legpp2->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legpp2->AddEntry(grFakePP,"Femtoscopic fit","l");
  legpp2->Draw("same");
  BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.55, 0.875, "ALICE Preliminary");
  if(strcmp(system, "pp") == 0) BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}} = %i TeV", system, energy));
  else BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %i TeV", system, energy));
  text.SetTextSize(gStyle->GetTextSize()*0.75);
  text.SetNDC();
  text.SetTextColor(1);
  if(!EPOS) text.DrawLatex(0.55, 0.775, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.55, 0.775, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) Can_CF_pp->Print("ANplot/CF_pp_EPOS_prelim.pdf");
  else Can_CF_pp->Print("ANplot/CF_pp_Gauss_prelim.pdf");

  TCanvas *Can_CF_pL = new TCanvas("pL","pL", 0,0,650,550);
  Can_CF_pL->SetRightMargin(right);
  Can_CF_pL->SetTopMargin(top);
  Tgraph_syserror_pL_ApAL->SetLineColor(kWhite);
  Tgraph_syserror_pL_ApAL->Draw("ap");
  baselinePL->Draw("same");
  Tgraph_syserror_pL_ApAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_pL_ApAL->GetYaxis()->SetRangeUser(0.8, 2.);
  if(grFemtopLNLO) grFemtopLNLO->Draw("l3 same");
  if(!EPOS && grFemtopLLO) grFemtopLLO->Draw("l3 same");
  Tgraph_syserror_pL_ApAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_pL_ApAL->Draw("2 same");
  hist_CF_Lp_ALAp_exp[2]->Draw("pe same");
  TLegend *legLp2;
  if(EPOS) legLp2 = new TLegend(0.53,0.545,0.85,0.76);
  else legLp2 = new TLegend(0.53,0.495,0.85,0.76);
  legLp2->SetBorderSize(0);
  legLp2->SetTextFont(42);
  legLp2->SetTextSize(gStyle->GetTextSize()*0.75);
  legLp2->AddEntry(hist_CF_Lp_ALAp_exp[2], "p#Lambda #oplus #bar{p}#bar{#Lambda} pairs", "pe");
  legLp2->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legLp2->AddEntry(grFakePLnlo,"Femtoscopic fit (#chiEFT NLO)","l");
  if(!EPOS) legLp2->AddEntry(grFakePLlo,"Femtoscopic fit (#chiEFT LO)","l");
  legLp2->Draw("same");
  BeamText.DrawLatex(0.55, 0.875, "ALICE Preliminary");
  if(strcmp(system, "pp") == 0) BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}} = %i TeV", system, energy));
  else BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %i TeV", system, energy));
  if(!EPOS) text.DrawLatex(0.55, 0.775, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.55, 0.775, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) ref.DrawLatex(0.61, 0.515, "Nucl. Phys. A915 (2013) 24.");
  else ref.DrawLatex(0.61, 0.465, "Nucl. Phys. A915 (2013) 24.");
  if(EPOS) Can_CF_pL->Print("ANplot/CF_pL_EPOS_prelim.pdf");
  else Can_CF_pL->Print("ANplot/CF_pL_Gauss_prelim.pdf");


  TCanvas *Can_CF_LL = new TCanvas("LL","LL", 0,0,650,550);
  Can_CF_LL->SetRightMargin(right);
  Can_CF_LL->SetTopMargin(top);
  Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);
  Tgraph_syserror_LL_ALAL->Draw("ap");
  baselineLL->Draw("same");
  Tgraph_syserror_LL_ALAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_LL_ALAL->GetYaxis()->SetRangeUser(0.35, 2.);
  if(grFemtoLL) grFemtoLL->Draw("l3 same");
  Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_LL_ALAL->Draw("2 same");
  hist_CF_LL_ALAL_exp[2]->Draw("pe same");
  TLegend *legLL2 = new TLegend(0.53,0.545,0.85,0.76);
  legLL2->SetBorderSize(0);
  legLL2->SetTextFont(42);
  legLL2->SetTextSize(gStyle->GetTextSize()*0.75);
  legLL2->AddEntry(hist_CF_pp_ApAp_exp[2], "#Lambda#Lambda #oplus #bar{#Lambda}#bar{#Lambda} pairs", "pe");
  legLL2->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legLL2->AddEntry(grFakeLL,"Femtoscopic fit","l");
  legLL2->Draw("same");
  BeamText.DrawLatex(0.55, 0.875, "ALICE Preliminary");
  if(strcmp(system, "pp") == 0) BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}} = %i TeV", system, energy));
  else BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %i TeV", system, energy));
  if(!EPOS) text.DrawLatex(0.55, 0.775, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.55, 0.775, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) Can_CF_LL->Print("ANplot/CF_LL_EPOS_prelim.pdf");
  else Can_CF_LL->Print("ANplot/CF_LL_Gauss_prelim.pdf");

  TCanvas *Can_CF_pXi = new TCanvas("pXi","pXi", 0,0,650,550);
  Can_CF_pXi->SetRightMargin(right);
  Can_CF_pXi->SetTopMargin(top);
  Tgraph_syserror_pXi_ApAXi->SetLineColor(kWhite);
  Tgraph_syserror_pXi_ApAXi->Draw("ap");
  baselinePXI->Draw("same");
  Tgraph_syserror_pXi_ApAXi->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pXi_ApAXi->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_pXi_ApAXi->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_pXi_ApAXi->GetYaxis()->SetRangeUser(0.75, 5.);
  Tgraph_syserror_pXi_ApAXi->SetFillColorAlpha(kBlack, 0.4);
  if(grFemtopXi) grFemtopXi->Draw("l3 same");
  if(grFemtopXiCoulomb) grFemtopXiCoulomb->Draw("l3 same");
  Tgraph_syserror_pXi_ApAXi->Draw("2 same");
  hist_CF_pXi_ApAXi_exp[2]->Draw("pe same");
  TLegend *legpXi2 = new TLegend(0.53,0.495,0.85,0.76);
  legpXi2->SetBorderSize(0);
  legpXi2->SetTextFont(42);
  legpXi2->SetTextSize(gStyle->GetTextSize()*0.75);
  legpXi2->AddEntry(hist_CF_pp_ApAp_exp[2], "p#Xi^{-} #oplus #bar{p}#Xi^{+} pairs", "pe");
  legpXi2->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legpXi2->AddEntry(grFakeXiCoulomb, "Femtoscopic fit (Coulomb)","l");
  legpXi2->AddEntry(grFakepXi,"Femtoscopic fit (HAL QCD)","l");
  legpXi2->Draw("same");
  ref.DrawLatex(0.61, 0.465, "Nucl. Phys. A967 (2017) 856.");
  ref.DrawLatex(0.61, 0.415, "PoS LATTICE2016 (2017) 116.");
  BeamText.DrawLatex(0.55, 0.875, "ALICE Preliminary");
  if(strcmp(system, "pp") == 0) BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}} = %i TeV", system, energy));
  else BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %i TeV", system, energy));
  if(!EPOS) text.DrawLatex(0.55, 0.775, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.55, 0.775, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) Can_CF_pXi->Print("ANplot/CF_pXi_EPOS_prelim.pdf");
  else Can_CF_pXi->Print("ANplot/CF_pXi_Gauss_prelim.pdf");

  auto *outfile = new TFile("CFout.root", "RECREATE");
  hist_CF_pp_ApAp_exp[2]->Write("pp");
  hist_CF_Lp_ALAp_exp[2]->Write("pL");
  hist_CF_LL_ALAL_exp[2]->Write("LL");
  hist_CF_pXi_ApAXi_exp[2]->Write("pXi");  
}
