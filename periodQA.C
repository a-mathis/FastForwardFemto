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
void SetStyle(bool graypalette=false, bool title=false)
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
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.055);
  histo->GetXaxis()->SetLabelOffset(0.02);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleSize(0.055);
  histo->GetYaxis()->SetLabelOffset(0.02);
  histo->GetYaxis()->SetTitleOffset(1.2);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// second order polynomial + double gaus to Lambda peak
void FitLambda(TH1F* histo, float &signal, float &signalErr, float &background, float &backgroundErr, float lowerBound, float upperBound)
{
//  histo->Sumw2();
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
//  signalOnly->Sumw2();
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
void SetStyleGraph(TGraph *graph, int marker, int color) {
  graph->GetXaxis()->SetLabelSize(0.05);
  graph->GetXaxis()->SetTitleSize(0.055);
  graph->GetXaxis()->SetLabelOffset(0.02);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetXaxis()->SetLabelFont(42);
  graph->GetYaxis()->SetLabelSize(0.05);
  graph->GetYaxis()->SetTitleSize(0.055);
  graph->GetYaxis()->SetLabelOffset(0.02);
  graph->GetYaxis()->SetTitleOffset(1.2);
  graph->SetMarkerStyle(fMarkers[marker]);
  graph->SetMarkerColor(fColors[color]);
  graph->SetLineColor(fColors[color]);
  graph->GetXaxis()->SetTitleOffset(1.15);
  graph->GetYaxis()->SetTitleOffset(1.25);
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
   TFitResultPtr r = signalOnly->Fit("fSignalGauss", "SRQ0", "", 1.29,1.38);

   // Extract signal as integral
   signal = fSignalGauss->Integral(lowerBound, upperBound)
                 /double(histo->GetBinWidth(1));
   signalErr = fSignalGauss->IntegralError(lowerBound, upperBound,
                                           r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())
                 /double(histo->GetBinWidth(1));

   TF1 *fLambda = new TF1("fLambda", "fBackground2 + fSignalGauss", 1.25,
                          1.4);
   fLambda->SetNpx(1000);
   fLambda->SetParameter(0, fBackground->GetParameter(0));
   fLambda->SetParameter(1, fBackground->GetParameter(1));
   fLambda->SetParameter(2, fBackground->GetParameter(2));
   fLambda->SetParameter(3, 0.75 * histo->GetMaximum());
   fLambda->SetParameter(4, fSignalGauss->GetParameter((1)));
   fLambda->SetParameter(5, fSignalGauss->GetParameter((2)));
   fLambda->SetParameter(6, 0.2 * histo->GetMaximum());
   fLambda->SetParameter(7, fSignalGauss->GetParameter((4)));
   fLambda->SetParameter(8, fSignalGauss->GetParameter((5)));
   fLambda->SetLineColor(fColors[2]);
   histo->Fit("fLambda", "SRQ", "", 1.28, 1.4);

   TF1 *fLambda_background = new TF1("fLambda_background", "pol2(0)",
                                     1.28, 1.45);
   fLambda_background->SetParameter(0, fBackground->GetParameter(0));
   fLambda_background->SetParameter(1, fBackground->GetParameter(1));
   fLambda_background->SetParameter(2, fBackground->GetParameter(2));
   fLambda_background->SetLineStyle(3);
   fLambda_background->SetLineColor(fColors[2]);

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

void periodQA(const char *path = ".", const char *prefix = "MB", const char *addon = "") {
  SetStyle();

  std::vector<const char *> periods = {
      {"LHC16d", "LHC16e", "LHC16g", "LHC16h", "LHC16i", "LHC16j", "LHC16k",
       "LHC16l", "LHC16o", "LHC16p", "LHC17c", "LHC17e", "LHC17f", "LHC17h",
       "LHC17i", "LHC17j", "LHC17k", "LHC17l", "LHC17m", "LHC17o", "LHC17r",
       "LHC18b", "LHC18d", "LHC18e", "LHC18f", "LHC18g", "LHC18h", "LHC18i",
       "LHC18j"}};
  auto histPeriod =
      new TH1F("histPeriod", ";;", periods.size(), 0, periods.size());
  TH1F *hInvMassSigma[periods.size()];
  TH1F *hInvMassSigmaCummul[periods.size()];

  auto grNevts = new TGraphErrors();
  auto grNProtons = new TGraphErrors();
  auto grNAntiProtons = new TGraphErrors();
  auto grNLambda = new TGraphErrors();
  auto grNAntiLambda = new TGraphErrors();
  auto grPurityLambda = new TGraphErrors();
  auto grPurityAntiLambda = new TGraphErrors();
  auto grMeanLambda = new TGraphErrors();
  auto grMeanAntiLambda = new TGraphErrors();
  auto grSigmaLambda = new TGraphErrors();
  auto grSigmaAntiLambda = new TGraphErrors();
  auto grNXi = new TGraphErrors();
  auto grNAntiXi = new TGraphErrors();
  auto grPurityXi = new TGraphErrors();
  auto grPurityAntiXi = new TGraphErrors();

  int counter = 0;
  float signal, signalErr, background, backgroundErr;
  const float lambdaLower = 1.115-0.004;
  const float lambdaUpper = 1.115+0.004;
  const float xiLower = 1.322-0.005;
  const float xiUpper = 1.322+0.005;

  for (int p = 0; p < static_cast<int>(periods.size()); ++p) {
    std::cout << "Processing period " << periods[p] << "\n";
    histPeriod->GetXaxis()->SetBinLabel(p + 1, Form("%s", periods[p]));
    auto _file0 = TFile::Open(Form("%s/AnalysisResults_%s.root", path, periods[p]));
    if (!_file0) continue;

    auto dirEvtCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sEvtCuts%s", prefix, addon)));
    TList *EvtCuts;
    dirEvtCuts->GetObject(Form("%sEvtCuts%s", prefix, addon), EvtCuts);
    auto hEvtCounter = (TH1F*)EvtCuts->FindObject("EventCounter");
//    const float nEvt = hEvtCounter->GetBinContent(hEvtCounter->GetNbinsX());
    const float nEvt = hEvtCounter->GetBinContent(2);

    auto dirTrackCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sTrackCuts%s", prefix, addon)));
    TList *TrackCuts;
    dirTrackCuts->GetObject(Form("%sTrackCuts%s", prefix, addon), TrackCuts);
    auto tmpFolder=(TList*)TrackCuts->FindObject("after");
    auto hNProton = (TH1F*)tmpFolder->FindObject("pTDist_after");
    const float nProton = hNProton->GetEntries();

    auto dirAntiTrackCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sAntiTrackCuts%s", prefix, addon)));
    TList *AntiTrackCuts;
    dirAntiTrackCuts->GetObject(Form("%sAntiTrackCuts%s", prefix, addon), AntiTrackCuts);
    tmpFolder=(TList*)AntiTrackCuts->FindObject("after");
    auto hNAntiProton = (TH1F*)tmpFolder->FindObject("pTDist_after");
    const float nAntiProton = hNAntiProton->GetEntries();

    auto dirv0Cuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sv0Cuts%s", prefix, addon)));
    TList *v0Cuts;
    dirv0Cuts->GetObject(Form("%sv0Cuts%s", prefix, addon), v0Cuts);
    tmpFolder=(TList*)v0Cuts->FindObject("v0Cuts");
    auto hLambdaMass = (TH1F*)tmpFolder->FindObject("InvMasswithCuts");
    tmpFolder=(TList*)tmpFolder->FindObject("after");
    auto hNLambda = (TH1F*)tmpFolder->FindObject("pTDist_after");
    const float nLambda = hNLambda->GetEntries();
    FitLambda(hLambdaMass, signal, signalErr, background, backgroundErr, lambdaLower, lambdaUpper);
    auto signalFit = (TF1*)hLambdaMass->GetListOfFunctions()->FindObject("fLambda");
    float mean2 = signalFit->GetParameter(4);
    float sigma2 = signalFit->GetParameter(5);
    float mean2err = signalFit->GetParError(4);
    float sigma2err = signalFit->GetParError(5);
    const float meanLambda = mean2;
    const float meanLambdaErr = mean2err;
    const float sigmaLambda = std::abs(sigma2);
    const float sigmaLambdaErr = std::abs(sigma2err);

    const float purityLambda = signal/(signal+background);
    const float purityLambdaErr = std::sqrt( signalErr*signalErr*background*background/std::pow(signal+background, 4) + backgroundErr*backgroundErr*signal*signal/std::pow(signal+background, 4) );

    auto dirAntiv0Cuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sAntiv0Cuts%s", prefix, addon)));
    TList *Antiv0Cuts;
    dirAntiv0Cuts->GetObject(Form("%sAntiv0Cuts%s", prefix, addon), Antiv0Cuts);
    tmpFolder=(TList*)Antiv0Cuts->FindObject("v0Cuts");
    auto hAntiLambdaMass = (TH1F*)tmpFolder->FindObject("InvMasswithCuts");
    tmpFolder=(TList*)tmpFolder->FindObject("after");
    auto hNAntiLambda = (TH1F*)tmpFolder->FindObject("pTDist_after");
    const float nAntiLambda = hNAntiLambda->GetEntries();
    FitLambda(hAntiLambdaMass, signal, signalErr, background, backgroundErr, lambdaLower, lambdaUpper);
    const float purityAntiLambda = signal/(signal+background);
    const float purityAntiLambdaErr = std::sqrt( signalErr*signalErr*background*background/std::pow(signal+background, 4) + backgroundErr*backgroundErr*signal*signal/std::pow(signal+background, 4) );
    signalFit = (TF1*)hAntiLambdaMass->GetListOfFunctions()->FindObject("fLambda");
    mean2 = signalFit->GetParameter(4);
    sigma2 = signalFit->GetParameter(5);
    mean2err = signalFit->GetParError(4);
    sigma2err = signalFit->GetParError(5);
    const float meanAntiLambda = mean2;
    const float meanAntiLambdaErr = mean2err;
    const float sigmaAntiLambda = std::abs(sigma2);
    const float sigmaAntiLambdaErr = std::abs(sigma2err);

    auto dirCascadeCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sCascadeCuts%s", prefix, addon)));
    TList *CascadeCuts;
    dirCascadeCuts->GetObject(Form("%sCascadeCuts%s", prefix, addon), CascadeCuts);
    tmpFolder=(TList*)CascadeCuts->FindObject("Cascade");
    auto hXiMass = (TH1F*)((TH2F*)tmpFolder->FindObject("InvMassXi"))->ProjectionY();
    tmpFolder=(TList*)tmpFolder->FindObject("after");
    auto hNXi = (TH1F*)tmpFolder->FindObject("XiPt_after");
    const float nXi = hNXi->GetEntries();
    FitXi(hXiMass, signal, signalErr, background, backgroundErr, xiLower, xiUpper);
    const float purityXi = signal/(signal+background);
    const float purityXiErr = std::sqrt( signalErr*signalErr*background*background/std::pow(signal+background, 4) + backgroundErr*backgroundErr*signal*signal/std::pow(signal+background, 4) );


    auto dirAntiCascadeCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sAntiCascadeCuts%s", prefix, addon)));
    TList *AntiCascadeCuts;
    dirAntiCascadeCuts->GetObject(Form("%sAntiCascadeCuts%s", prefix, addon), AntiCascadeCuts);
    tmpFolder=(TList*)AntiCascadeCuts->FindObject("Cascade");
    auto hAntiXiMass = (TH1F*)((TH2F*)tmpFolder->FindObject("InvMassXi"))->ProjectionY();
    tmpFolder=(TList*)tmpFolder->FindObject("after");
    auto hNAntiXi = (TH1F*)tmpFolder->FindObject("XiPt_after");
    const float nAntiXi = hNAntiXi->GetEntries();
    FitXi(hAntiXiMass, signal, signalErr, background, backgroundErr, xiLower, xiUpper);
    const float purityAntiXi = signal/(signal+background);
    const float purityAntiXiErr = std::sqrt( signalErr*signalErr*background*background/std::pow(signal+background, 4) + backgroundErr*backgroundErr*signal*signal/std::pow(signal+background, 4) );


    grNevts->SetPoint(counter, p, nEvt);
    grNProtons->SetPoint(counter, p, nProton/nEvt);
    grNAntiProtons->SetPoint(counter, p, nAntiProton/nEvt);
    grNLambda->SetPoint(counter, p, nLambda/nEvt);
    grPurityLambda->SetPoint(counter, p, purityLambda * 100.f);
    grPurityLambda->SetPointError(counter, 0, purityLambdaErr * 100.f);
    grMeanLambda->SetPoint(counter, p, meanLambda);
    grMeanLambda->SetPointError(counter, 0, meanLambdaErr);
    grSigmaLambda->SetPoint(counter, p, sigmaLambda);
    grSigmaLambda->SetPointError(counter, 0, sigmaLambdaErr);
    grNAntiLambda->SetPoint(counter, p, nAntiLambda/nEvt);
    grPurityAntiLambda->SetPoint(counter, p, purityAntiLambda * 100.f);
    grPurityAntiLambda->SetPointError(counter, 0, purityAntiLambdaErr * 100.f);
    grMeanAntiLambda->SetPoint(counter, p, meanAntiLambda);
    grMeanAntiLambda->SetPointError(counter, 0, meanAntiLambdaErr);
    grSigmaAntiLambda->SetPoint(counter, p, sigmaAntiLambda);
    grSigmaAntiLambda->SetPointError(counter, 0, sigmaAntiLambdaErr);
    grNXi->SetPoint(counter, p, nXi/nEvt);
    grPurityXi->SetPoint(counter, p, purityXi * 100.f);
    grPurityXi->SetPointError(counter, 0, purityXiErr * 100.f);
    grNAntiXi->SetPoint(counter, p, nAntiXi/nEvt);
    grPurityAntiXi->SetPoint(counter, p, purityAntiXi * 100.f);
    grPurityAntiXi->SetPointError(counter, 0, purityAntiXiErr * 100.f);

    counter++;
  }

  SetStyleHisto(histPeriod, 0, 1);

  auto nEvt = new TCanvas();
  nEvt->SetLogy();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; Number of events");
  SetStyleGraph(grNevts, 0, 1);
  grNevts->Draw("pez");
  histPeriod->SetMaximum(grNevts->GetHistogram()->GetMaximum());
  histPeriod->SetMinimum(grNevts->GetHistogram()->GetMinimum());
  nEvt->Print("PeriodQA/nEvts.pdf");

  auto nProton = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #(p)/evt.");
  SetStyleGraph(grNProtons, 0, 1);
  grNProtons->Draw("pez");
  histPeriod->SetMaximum(0.15);
  histPeriod->SetMinimum(0.05);
  nProton->Print("PeriodQA/nProton.pdf");

  auto nAntiProton = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #(#bar{p})/evt.");
  SetStyleGraph(grNAntiProtons, 0, 1);
  grNAntiProtons->Draw("pez");
  histPeriod->SetMaximum(0.15);
  histPeriod->SetMinimum(0.05);
  nAntiProton->Print("PeriodQA/nAntiProton.pdf");

  auto nLambda = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #(#Lambda)/evt.");
  SetStyleGraph(grNLambda, 0, 1);
  grNLambda->Draw("pez");
  histPeriod->SetMaximum(0.05);
  histPeriod->SetMinimum(0.0);
  nLambda->Print("PeriodQA/nLambda.pdf");

  auto purityLambda = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; Purity (#Lambda) (%)");
  SetStyleGraph(grPurityLambda, 0, 1);
  grPurityLambda->Draw("pez");
  histPeriod->SetMaximum(105);
  histPeriod->SetMinimum(80);
  purityLambda->Print("PeriodQA/PurityLambda.pdf");

  auto meanLambda = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #mu (#Lambda) (GeV/#it{c}^{2})");
  SetStyleGraph(grMeanLambda, 0, 1);
  grMeanLambda->Draw("pez");
  histPeriod->SetMaximum(1.13);
  histPeriod->SetMinimum(1.1);
  meanLambda->Print("PeriodQA/MeanLambda.pdf");

  auto sigmaLambda = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #sigma (#Lambda) (GeV/#it{c}^{2})");
  SetStyleGraph(grSigmaLambda, 0, 1);
  grSigmaLambda->Draw("pez");
  histPeriod->SetMaximum(0.005);
  histPeriod->SetMinimum(0);
  sigmaLambda->Print("PeriodQA/SigmaLambda.pdf");

  auto nAntiLambda = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #(#bar{#Lambda})/evt.");
  SetStyleGraph(grNAntiLambda, 0, 1);
  grNAntiLambda->Draw("pez");
  histPeriod->SetMaximum(0.05);
  histPeriod->SetMinimum(0.0);
  nAntiLambda->Print("PeriodQA/nAntiLambda.pdf");

  auto purityAntiLambda = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; Purity (#bar{#Lambda}) (%)");
  SetStyleGraph(grPurityAntiLambda, 0, 1);
  grPurityAntiLambda->Draw("pez");
  histPeriod->SetMaximum(105);
  histPeriod->SetMinimum(80);
  purityAntiLambda->Print("PeriodQA/PurityAntiLambda.pdf");

  auto meanAntiLambda = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #mu (#bar{#Lambda}) (GeV/#it{c}^{2})");
  SetStyleGraph(grMeanAntiLambda, 0, 1);
  grMeanAntiLambda->Draw("pez");
  histPeriod->SetMaximum(1.13);
  histPeriod->SetMinimum(1.1);
  meanAntiLambda->Print("PeriodQA/MeanAntiLambda.pdf");

  auto sigmaAntiLambda = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #sigma (#bar{#Lambda}) (GeV/#it{c}^{2})");
  SetStyleGraph(grSigmaAntiLambda, 0, 1);
  grSigmaAntiLambda->Draw("pez");
  histPeriod->SetMaximum(0.005);
  histPeriod->SetMinimum(0);
  sigmaAntiLambda->Print("PeriodQA/SigmaAntiLambda.pdf");

  auto nXi = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #(#Xi)/evt.");
  SetStyleGraph(grNXi, 0, 1);
  grNXi->Draw("pez");
  histPeriod->SetMaximum(0.002);
  histPeriod->SetMinimum(0.0);
  nXi->Print("PeriodQA/nXi.pdf");

  auto purityXi= new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; Purity (#Xi) (%)");
  SetStyleGraph(grPurityXi, 0, 1);
  grPurityXi->Draw("pez");
  histPeriod->SetMaximum(105);
  histPeriod->SetMinimum(80);
  purityXi->Print("PeriodQA/PurityXi.pdf");

  auto nAntiXi = new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; #(#bar{#Xi})/evt.");
  SetStyleGraph(grNAntiXi, 0, 1);
  grNAntiXi->Draw("pez");
  histPeriod->SetMaximum(0.002);
  histPeriod->SetMinimum(0.0);
  nAntiXi->Print("PeriodQA/nAntiXi.pdf");

  auto purityAntiXi= new TCanvas();
  histPeriod->Draw();
  histPeriod->SetTitle("; ; Purity (#bar{#Xi}) (%)");
  SetStyleGraph(grPurityAntiXi, 0, 1);
  grPurityAntiXi->Draw("pez");
  histPeriod->SetMaximum(105);
  histPeriod->SetMinimum(80);
  purityAntiXi->Print("PeriodQA/PurityAntiXi.pdf");
}
