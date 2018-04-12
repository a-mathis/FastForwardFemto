#include "TPad.h"
Color_t fFillColors[8] = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9,
    kOrange-9, kCyan-8, kYellow-7};
Color_t fColors[8]  = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1,
    kOrange-1, kCyan+2, kYellow+2};
Style_t fMarkers[10] = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare,
    kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a and c are the weighting entities, i.e. wMean = (weightA*A + weightB*d) / (a+weightB)
float weightedMean(float weightA, float A, float weightB, float B)
{
  return (weightA*A + weightB*B)/(weightA+weightB);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a and weightB are the weighting entities, i.e. wMean = (weightA*A + weightB*B) / (a+weightB)
float weightedMeanError(float weightA, float A, float weightB, float B,
                        float weightAErr, float AErr, float weightBErr, float BErr)
{
  return std::sqrt(weightAErr*weightAErr*std::pow(((weightB*(A -
      B))/((weightA + weightB)*(weightA + weightB))), 2) + AErr*AErr*
                   std::pow(weightA/(weightA + weightB), 2)  +
                   weightBErr*weightBErr*std::pow((weightA* (-A + B))/((weightA +
                       weightB)*(weightA + weightB)) , 2) + BErr*BErr*
                       std::pow(weightB/(weightA + weightB), 2));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *getSignalHisto(TF1 *function, TH1F *histo, float
                     rangeLow, float rangeHigh, const char *name)
{
  const int firstBin = histo->FindBin(rangeLow);
  const int lastBin  = histo->FindBin(rangeHigh);
  TH1F *result = new TH1F(Form("result_%.2f_%.2f_%s", rangeLow, rangeHigh,
                               name), "", histo->GetNbinsX(), histo->GetXaxis()->GetXmin(),
                          histo->GetXaxis()->GetXmax());
  for(int i = firstBin; i<lastBin; ++i) {
    float weight = histo->GetBinContent(i) -
        function->Eval(histo->GetBinCenter(i));
    result->Fill(histo->GetBinCenter(i), weight);
    result->SetBinError(i, histo->GetBinError(i));
  }
  result->SetFillColor(kGray+1);
  result->SetLineColor(kGray+1);
  return result;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette, bool title)
{
  gStyle->Reset("Plain");
  gStyle->SetCanvasPreferGL(1);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  auto ci = 1179;
  auto color = new TColor(ci, 1, 0, 0, " ", 0.);
  gStyle->SetStatColor(ci);
  gStyle->SetTitleColor(ci);
  gStyle->SetCanvasColor(ci);
  gStyle->SetPadColor(ci);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //  gStyle->SetPadBottomMargin(0.15);
  //  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.16);
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
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
  histo->GetYaxis()->SetTitleOffset(1.2);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// second order polynomial + double gaus to Lambda peak
void FitLambda(TH1F* histo, float &signal, float &signalErr, float &background, float &backgroundErr, float lowerBound, float upperBound)
{
  histo->Sumw2();
  // Fit Background with second order polynomial, excluding Mlambda +/- 10 MeV
  TF1 *fBackground = new TF1("fBackground", [&](double *x, double *p) { if (x[0] > 1.1075 && x[0] < 1.1235) {TF1::RejectPoint(); return (double)0; } return p[0] + p[1]*x[0] + p[2]*x[0]*x[0]; }, 1.095, 1.15, 3);
  TFitResultPtr backgroundR = histo->Fit("fBackground", "SRQ0", "", 1.095, 1.15);

  // parse then to proper TF1
  TF1 *fBackground2 = new TF1("fBackground2","pol2", 0, 1.5);
  fBackground2->SetParameter(0, fBackground->GetParameter(0));
  fBackground2->SetParameter(1, fBackground->GetParameter(1));
  fBackground2->SetParameter(2, fBackground->GetParameter(2));

  // remove background from signal
  TH1F *signalOnly = getSignalHisto(fBackground2, histo, 1.0, 1.3, Form("%s_signal_only", histo->GetName()));
  signalOnly->Sumw2();
  signalOnly->Draw("same");

  // fit signal only
  TF1 *fSignalSingleGauss = new TF1("fSignalSingleGauss", "gaus", 1.095, 1.15);
//  fSignalSingleGauss->SetParameter(1, 1.115);
  signalOnly->Fit("fSignalSingleGauss", "SRQ0", "", 1.1075, 1.1235);

  TF1 *fSignalGauss = new TF1("fSignalGauss", "gaus(0) + gaus(3)", 1.1, 1.3);
  fSignalGauss->SetParameter(0, 0.05 * histo->GetMaximum());
  fSignalGauss->SetParameter(1, fSignalSingleGauss->GetParameter(1));
  fSignalGauss->SetParLimits(1, 1.115-0.01, 1.115+0.01);
  fSignalGauss->SetParameter(2, 5.f*fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParameter(3, 0.95 * histo->GetMaximum());
  fSignalGauss->SetParameter(4, fSignalSingleGauss->GetParameter(1));
  fSignalGauss->SetParLimits(4, 1.115-0.01, 1.115+0.01);
  fSignalGauss->SetParameter(5, fSignalSingleGauss->GetParameter(2));
  TFitResultPtr r = signalOnly->Fit("fSignalGauss", "SRQ0", "", 1.1075, 1.1235);

  // Extract signal as integral
  signal = fSignalGauss->Integral(lowerBound, upperBound) /double(histo->GetBinWidth(1));
  signalErr = fSignalGauss->IntegralError(lowerBound, upperBound, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()) /double(histo->GetBinWidth(1));

  TF1 *fLambda = new TF1("fLambda", "fBackground2 + fSignalGauss", 1.1, 1.13);
  fLambda->SetNpx(1000);
  fLambda->SetParameter(3, 0.75 * histo->GetMaximum());
  fLambda->SetParameter(4, fSignalGauss->GetParameter((1)));
  fLambda->SetParameter(5, fSignalGauss->GetParameter((2)));
  fLambda->SetParameter(6, 0.2 * histo->GetMaximum());
  fLambda->SetParameter(7, fSignalGauss->GetParameter((4)));
  fLambda->SetParameter(8, fSignalGauss->GetParameter((5)));
  fLambda->SetLineColor(fColors[1]);
  histo->Fit("fLambda", "SRQ", "", 1.095, 1.15);

  TF1 *fLambda_background = new TF1("fLambda_background", "pol2(0)", 1.05, 1.25);
  fLambda_background->SetParameter(0, fLambda->GetParameter(0));
  fLambda_background->SetParameter(1, fLambda->GetParameter(1));
  fLambda_background->SetParameter(2, fLambda->GetParameter(2));
  fLambda_background->SetLineStyle(3);
  fLambda_background->SetLineColor(fColors[1]);

  background = fLambda_background->Integral(lowerBound, upperBound) /double(histo->GetBinWidth(1));
  backgroundErr = fLambda_background->IntegralError(lowerBound, upperBound, backgroundR->GetParams(), backgroundR->GetCovarianceMatrix().GetMatrixArray()) /double(histo->GetBinWidth(1));

  histo->GetListOfFunctions()->Add(fLambda_background);

  delete signalOnly;
  delete fSignalGauss;
  delete fSignalSingleGauss;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// second order polynomial + double gaus to Xi peak
void FitXi(TH1F* histo, float &signal, float
               &signalErr, float &background, float &backgroundErr, float lowerBound,
               float upperBound)
{
  // Fit Background with second order polynomial, excluding Mlambda +/- 10 MeV
  TF1 *fBackground = new TF1("fBackground", [&](double *x, double *p) {
    if (x[0] > 1.31 && x[0] < 1.335) {TF1::RejectPoint(); return
        (double)0; } return p[0] + p[1]*x[0] + p[2]*x[0]*x[0]; }, 1.29, 1.35, 3);
  TFitResultPtr backgroundR = histo->Fit("fBackground", "SRQ0", "",1.295,1.349);

  // parse then to proper TF1
  TF1 *fBackground2 = new TF1("fBackground2","pol2", 0, 1.4);
  fBackground2->SetParameter(0, fBackground->GetParameter(0));
  fBackground2->SetParameter(1, fBackground->GetParameter(1));
  fBackground2->SetParameter(2, fBackground->GetParameter(2));

  // remove background from signal
  TH1F *signalOnly = getSignalHisto(fBackground2, histo, 1.3, 1.34,
                                    Form("%s_signal_only", histo->GetName()));

  // fit signal only
  TF1 *fSignalSingleGauss = new TF1("fSignalSingleGauss", "gaus(0)",1.31,1.34);
  //  signalOnly->DrawCopy();
  signalOnly->Fit("fSignalSingleGauss", "Q");
  TF1 *fSignalGauss = new TF1("fSignalGauss", "gaus(0) + gaus(3)", 1.3,
                              1.4);
  fSignalGauss->SetParameter(0, 0.75 * histo->GetMaximum());
  fSignalGauss->SetParameter(1, fSignalSingleGauss->GetParameter(1));
  fSignalGauss->SetParameter(2, 2.f*fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParLimits(2, 0.5*fSignalSingleGauss->GetParameter(2),1e2*2.f*fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParameter(3, 0.2 * histo->GetMaximum());
  fSignalGauss->SetParameter(4, fSignalSingleGauss->GetParameter(1));
  fSignalGauss->SetParLimits(4,
                             fSignalSingleGauss->GetParameter(1)-fSignalSingleGauss->GetParameter(2),
                             fSignalSingleGauss->GetParameter(1)+fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParameter(5, 0.5*fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParLimits(5, 0.5*fSignalSingleGauss->GetParameter(2),1e2*2.f*fSignalSingleGauss->GetParameter(2));
  TFitResultPtr r = signalOnly->Fit("fSignalGauss", "SRQ0", "", 1.29,
                                    1.38);

  // Extract signal as integral
  signal = fSignalGauss->Integral(lowerBound, upperBound)
                /double(histo->GetBinWidth(1));
  signalErr = fSignalGauss->IntegralError(lowerBound, upperBound,
                                          r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())
                /double(histo->GetBinWidth(1));

  TF1 *fLambda = new TF1("fLambda", "fBackground2 + fSignalGauss", 1.25,
                         1.4);
  fLambda->SetNpx(1000);
  fLambda->SetParameter(3, 0.75 * histo->GetMaximum());
  fLambda->SetParameter(4, fSignalGauss->GetParameter((1)));
  fLambda->SetParameter(5, fSignalGauss->GetParameter((2)));
  fLambda->SetParameter(6, 0.2 * histo->GetMaximum());
  fLambda->SetParameter(7, fSignalGauss->GetParameter((4)));
  fLambda->SetParameter(8, fSignalGauss->GetParameter((5)));
  fLambda->SetLineColor(fColors[1]);
  histo->Fit("fLambda", "SRQ", "", 1.28, 1.4);

  TF1 *fLambda_background = new TF1("fLambda_background", "pol2(0)",
                                    1.28, 1.45);
  fLambda_background->SetParameter(0, fLambda->GetParameter(0));
  fLambda_background->SetParameter(1, fLambda->GetParameter(1));
  fLambda_background->SetParameter(2, fLambda->GetParameter(2));
  fLambda_background->SetLineStyle(3);
  fLambda_background->SetLineColor(fColors[1]);

  background = fLambda_background->Integral(lowerBound, upperBound)
                /double(histo->GetBinWidth(1));
  backgroundErr = fLambda_background->IntegralError(lowerBound,
                                                    upperBound, backgroundR->GetParams(),
                                                    backgroundR->GetCovarianceMatrix().GetMatrixArray())
                /double(histo->GetBinWidth(1));

  histo->GetListOfFunctions()->Add(fLambda_background);

  delete signalOnly;
  delete fSignalGauss;
  delete fSignalSingleGauss;
}

void PurityByPeriod(const char *fileName, const char *prefix) {
  SetStyle(false,true);
  TH2F *V0PurityPerRun;
  TH2F *AntiV0PurityPerRun;
  TH2F *CascPurityPerRun;
  TH2F *AntiCascPurityPerRun;

  TFile *file=TFile::Open(fileName);

  TString v0CutsName="";
  v0CutsName+=prefix;
  v0CutsName+="v0Cuts";

  TDirectoryFile *dirv0Cuts=(TDirectoryFile*)(file->FindObjectAny(v0CutsName.Data()));
  if (dirv0Cuts) {
    TList *tmp=dirv0Cuts->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *v0CutList;
    dirv0Cuts->GetObject(name,v0CutList);
    if (v0CutList) {
      TList *v0Cuts = (TList*)v0CutList->FindObject("v0Cuts");
      if (v0Cuts) {
        V0PurityPerRun=(TH2F*)v0Cuts->FindObject("InvMassPerRunnumber");
        V0PurityPerRun->SetNameTitle("InvMassPerRunnumberV0","InvMassPerRunnumberV0");
      }
    }
  }

  TString Antiv0CutsName="";
  Antiv0CutsName+=prefix;
  Antiv0CutsName+="Antiv0Cuts";

  TDirectoryFile *dirAntiv0Cuts=(TDirectoryFile*)(file->FindObjectAny(Antiv0CutsName.Data()));
  if (dirAntiv0Cuts) {
    TList *tmp=dirAntiv0Cuts->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *Antiv0CutList;
    dirAntiv0Cuts->GetObject(name,Antiv0CutList);
    if (Antiv0CutList) {
      TList *Antiv0Cuts=(TList*)Antiv0CutList->FindObject("v0Cuts");
      if (Antiv0Cuts) {
        AntiV0PurityPerRun=(TH2F*)Antiv0Cuts->FindObject("InvMassPerRunnumber");
        AntiV0PurityPerRun->SetNameTitle("InvMassPerRunnumberAV0","InvMassPerRunnumberAV0");
      }
    }
  }

  TString CascCutsName="";
  CascCutsName+=prefix;
  CascCutsName+="CascadeCuts";

  TDirectoryFile *dirCascCuts=(TDirectoryFile*)(file->FindObjectAny(CascCutsName.Data()));
  if (dirCascCuts) {
    TList *tmp=dirCascCuts->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *CascCutList;
    dirCascCuts->GetObject(name,CascCutList);
    if (CascCutList) {
      TList *CascCuts = (TList*)CascCutList->FindObject("Cascade");
      if (CascCuts) {
        CascPurityPerRun=(TH2F*)CascCuts->FindObject("InvMassPerRunnumber");
        CascPurityPerRun->SetNameTitle("InvMassPerRunnumberCasc","InvMassPerRunnumberCasc");
      }
    }
  }

  TString AntiCascCutsName="";
  AntiCascCutsName+=prefix;
  AntiCascCutsName+="AntiCascadeCuts";

  TDirectoryFile *dirAntiCascCuts=(TDirectoryFile*)(file->FindObjectAny(AntiCascCutsName.Data()));
  if (dirAntiCascCuts) {
    TList *tmp=dirAntiCascCuts->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *AntiCascCutList;
    dirAntiCascCuts->GetObject(name,AntiCascCutList);
    if (AntiCascCutList) {
      TList *AntiCascCuts = (TList*)AntiCascCutList->FindObject("Cascade");
      if (AntiCascCuts) {
        AntiCascPurityPerRun=(TH2F*)AntiCascCuts->FindObject("InvMassPerRunnumber");
        AntiCascPurityPerRun->SetNameTitle("InvMassPerRunnumberAntiCasc","InvMassPerRunnumberAntiCasc");
      }
    }
  }

  const int nOutput= V0PurityPerRun->GetXaxis()->GetNbins();
  TH2F *V0purity;
  TH2F *AntiV0purity;
  TH2F *Cascpurity;
  TH2F *AntiCascpurity;

  TH1F* V0iRun[nOutput];
  TH1F* AntiV0iRun[nOutput];
  TH1F* CasciRun[nOutput];
  TH1F* AntiCasciRun[nOutput];

  TH1F* FinalPurityV0;
  TH1F* FinalPurityAV0;
  TH1F* FinalPurityCasc;
  TH1F* FinalPurityACasc;

  std::vector<int> iRuns;
  std::vector<float> _valV0Purity;
  std::vector<float> _valAntiV0Purity;
  std::vector<float> _valCascPurity;
  std::vector<float> _valAntiCascPurity;

  int histCounterV0=0;
  TCanvas *cLambda = new TCanvas("cLambda","cLambda",0,0,1000,550);
  cLambda->Divide(2,2);

  int histCounterAV0=0;
  TCanvas *cALambda = new TCanvas("cALambda","cALambda",0,0,1000,550);
  cALambda->Divide(2,2);

  int histCounterCasc=0;
  TCanvas *cXi= new TCanvas("cXi","cXi",0,0,1000,550);
  cXi->Divide(2,2);

  int histCounterACasc=0;
  TCanvas *cAXi = new TCanvas("cAXi","cAXi",0,0,1000,550);
  cAXi->Divide(2,2);

  for (int iPeriod=1;iPeriod<nOutput;++iPeriod) {
    TString RunNumber=Form("%.0f",V0PurityPerRun->GetXaxis()->GetBinLowEdge(iPeriod));
    V0iRun[iPeriod]=(TH1F*)V0PurityPerRun->ProjectionY(Form("V0%s",RunNumber.Data()),iPeriod,iPeriod);
    if (V0iRun[iPeriod]->GetEffectiveEntries() > 0) {
      std::cout << RunNumber << '\n';
      iRuns.push_back((int)V0PurityPerRun->GetXaxis()->GetBinLowEdge(iPeriod));
      histCounterV0++;
      float signal=0;
      float signalErr=0;
      float background=0;
      float backgroundErr=0;
      float lowerBound=1.322-0.005;
      float upperBound=1.322+0.005;
      const float marginLambda = 0.004;
      const float massLambda = 1.115;
      V0iRun[iPeriod]->GetXaxis()->SetRangeUser(1.05,1.18);
      V0iRun[iPeriod]->GetXaxis()->SetTitle("IM_{#pi p} (GeV/#it{c}^{2})");
      V0iRun[iPeriod]->SetTitle(Form("%s",RunNumber.Data()));
      cLambda->cd(histCounterV0);
//      V0iRun[iPeriod]->DrawCopy();
      FitLambda(V0iRun[iPeriod], signal, signalErr, background, backgroundErr, massLambda-marginLambda,massLambda+marginLambda);
      float purityV0=signal/(signal+background)*100.f;
      std::cout << purityV0 << '\n';
      _valV0Purity.push_back(purityV0);
    }

    AntiV0iRun[iPeriod]=(TH1F*)AntiV0PurityPerRun->ProjectionY(Form("AV0%s",RunNumber.Data()),iPeriod,iPeriod);
    if (AntiV0iRun[iPeriod]->GetEffectiveEntries() > 0) {
      histCounterAV0++;
      float signal=0;
      float signalErr=0;
      float background=0;
      float backgroundErr=0;
      float lowerBound=1.322-0.005;
      float upperBound=1.322+0.005;
      const float marginLambda = 0.004;
      const float massLambda = 1.115;
      AntiV0iRun[iPeriod]->GetXaxis()->SetRangeUser(1.05,1.18);
      AntiV0iRun[iPeriod]->GetXaxis()->SetTitle("IM_{#pi p} (GeV/#it{c}^{2})");
      AntiV0iRun[iPeriod]->SetTitle(Form("%s",RunNumber.Data()));
      cALambda->cd(histCounterAV0);
//      AntiV0iRun[iPeriod]->DrawCopy();
      FitLambda(AntiV0iRun[iPeriod], signal, signalErr, background, backgroundErr, massLambda-marginLambda,massLambda+marginLambda);
      float purityAntiV0=signal/(signal+background)*100.f;
      _valAntiV0Purity.push_back(purityAntiV0);
    }
    CasciRun[iPeriod]=(TH1F*)CascPurityPerRun->ProjectionY(Form("Casc%s",RunNumber.Data()),iPeriod,iPeriod);
    if (CasciRun[iPeriod]->GetEffectiveEntries() > 0) {
      histCounterCasc++;
      float signal=0;
      float signalErr=0;
      float background=0;
      float backgroundErr=0;
      float lowerBound=1.322-0.005;
      float upperBound=1.322+0.005;
      const float marginLambda = 0.004;
      const float massLambda = 1.115;
      CasciRun[iPeriod]->GetXaxis()->SetRangeUser(1.05,1.18);
      CasciRun[iPeriod]->GetXaxis()->SetTitle("IM_{#pi#Lambda} (GeV/#it{c}^{2})");
      CasciRun[iPeriod]->SetTitle(Form("%s",RunNumber.Data()));
      cXi->cd(histCounterCasc);
      CasciRun[iPeriod]->DrawCopy();
      FitXi(CasciRun[iPeriod], signal, signalErr, background, backgroundErr, lowerBound,upperBound);
      float purityCasc=signal/(signal+background)*100.f;
      _valCascPurity.push_back(purityCasc);
    }
    AntiCasciRun[iPeriod]=(TH1F*)AntiCascPurityPerRun->ProjectionY(Form("ACasc%s",RunNumber.Data()),iPeriod,iPeriod);
    if (AntiCasciRun[iPeriod]->GetEffectiveEntries() > 0) {
      histCounterACasc++;
      float signal=0;
      float signalErr=0;
      float background=0;
      float backgroundErr=0;
      float lowerBound=1.322-0.005;
      float upperBound=1.322+0.005;
      const float marginLambda = 0.004;
      const float massLambda = 1.115;
      AntiCasciRun[iPeriod]->GetXaxis()->SetRangeUser(1.05,1.18);
      AntiCasciRun[iPeriod]->GetXaxis()->SetTitle("IM_{#pi#Lambda} (GeV/#it{c}^{2})");
      AntiCasciRun[iPeriod]->SetTitle(Form("%s",RunNumber.Data()));
      cAXi->cd(histCounterACasc);
      AntiCasciRun[iPeriod]->DrawCopy();
      FitXi(AntiCasciRun[iPeriod], signal, signalErr, background, backgroundErr, lowerBound,upperBound);
      float purityAntiCasc=signal/(signal+background)*100.f;
      _valAntiCascPurity.push_back(purityAntiCasc);
    }
  }
  const int nRuns=iRuns.size();
  float runNumbArray[nRuns];
  std::copy(iRuns.begin(),iRuns.end(),runNumbArray);

  TCanvas *cPurity= new TCanvas("cPurity","cPurity",0,0,1000,550);
  cPurity->Divide(2,2);

  FinalPurityV0=new TH1F("Purity Lambda","Purity Lambda",nRuns,0,nRuns);
  FinalPurityV0->GetXaxis()->SetTitle("Run Number");
  FinalPurityV0->GetYaxis()->SetTitle("Purity");
  FinalPurityV0->GetYaxis()->SetRangeUser(80,105);
  FinalPurityAV0=new TH1F("Purity Anti-Lambda","Purity Anti-Lambda",nRuns,0,nRuns);
  FinalPurityAV0->GetXaxis()->SetTitle("Run Number");
  FinalPurityAV0->GetYaxis()->SetTitle("Purity");
//  FinalPurityAV0->GetYaxis()->SetRangeUser(0.8,1.05);
  FinalPurityCasc=new TH1F("Purity Xi","Purity Xi",nRuns,0,nRuns);
  FinalPurityCasc->GetXaxis()->SetTitle("Run Number");
  FinalPurityCasc->GetYaxis()->SetTitle("Purity");
  FinalPurityCasc->GetYaxis()->SetRangeUser(0.8,1.05);
  FinalPurityACasc=new TH1F("Purity Anti-Xi","Purity Anti-Xi",nRuns,0,nRuns);
  FinalPurityACasc->GetXaxis()->SetTitle("Run Number");
  FinalPurityACasc->GetYaxis()->SetTitle("Purity");
  FinalPurityACasc->GetYaxis()->SetRangeUser(0.8,1.05);
  for (int i=1;i<nRuns+1;++i) {
    FinalPurityV0->SetBinContent(i,_valV0Purity.at(i-1));
    FinalPurityV0->GetXaxis()->SetBinLabel(i,Form("%i",iRuns.at(i-1)));
    FinalPurityAV0->SetBinContent(i,_valAntiV0Purity.at(i-1));
    FinalPurityAV0->GetXaxis()->SetBinLabel(i,Form("%i",iRuns.at(i-1)));
//    FinalPurityCasc->SetBinContent(i+1,_valCascPurity.at(i));
//    FinalPurityACasc->SetBinContent(i+1,_valAntiCascPurity.at(i));
  }



  cPurity->cd(1);
  FinalPurityV0->DrawCopy();
  cPurity->cd(2);
  FinalPurityAV0->DrawCopy();
//  cPurity->cd(3);
//  FinalPurityCasc->DrawCopy();
//  cPurity->cd(4);
//  FinalPurityACasc->DrawCopy();
}
