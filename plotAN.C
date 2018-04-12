#include "TPad.h"
std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7};
std::vector<int> fColors     = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2};
std::vector<int> fMarkers    = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette=false, bool title=true)
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
  gStyle->SetPalette(57);
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
  histo->SetLineWidth(2);
  histo->Sumw2();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHistoMany(TH1 *histo, int marker, int color)
{
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.055);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleSize(0.055);
  histo->GetYaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHisto(TH2 *histo)
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
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleGraph(TGraph *graph, int marker, int color)
{
  graph->SetMarkerStyle(fMarkers[marker]);
  graph->SetMarkerColor(fColors[color]);
  graph->SetLineColor(fColors[color]);
  graph->GetYaxis()->SetTitleOffset(1.25);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *getSignalHisto(TF1 *function, TH1F *histo, float rangeLow, float rangeHigh, const char *name)
{
  const int firstBin = histo->FindBin(rangeLow);
  const int lastBin  = histo->FindBin(rangeHigh);
  TH1F *result = new TH1F(Form("result_%f_%f_%s", rangeLow, rangeHigh, name), "", histo->GetNbinsX(), histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
  for(int i = firstBin; i<lastBin; ++i) {
    float weight = histo->GetBinContent(i) - function->Eval(histo->GetBinCenter(i));
    result->Fill(histo->GetBinCenter(i), weight);
    result->SetBinError(i, histo->GetBinError(i));
  }
  result->SetFillColor(fFillColors[0]);
  result->SetLineColor(fFillColors[0]);
  result->Sumw2();
  return result;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a and c are the weighting entities, i.e. wMean = (weightA*A + weightB*d) / (a+weightB)
float weightedMean(float weightA, float A, float weightB, float B)
{
  return (weightA*A + weightB*B)/(weightA+weightB);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a and weightB are the weighting entities, i.e. wMean = (weightA*A + weightB*B) / (a+weightB)
float weightedMeanError(float weightA, float A, float weightB, float B, float weightAErr, float AErr, float weightBErr, float BErr)
{
  return std::sqrt(weightAErr*weightAErr*std::pow(((weightB*(A - B))/((weightA + weightB)*(weightA + weightB))), 2) + AErr*AErr* std::pow(weightA/(weightA + weightB), 2)  + weightBErr*weightBErr*std::pow((weightA* (-A + B))/((weightA + weightB)*(weightA + weightB)) , 2) + BErr*BErr* std::pow(weightB/(weightA + weightB), 2));
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void plotAN(const char* file="~/Results/LHC17p_fast/AnalysisResults.root")
{
  SetStyle();

  float nTriggers = 0.f;
  float nLambda = 0.f;
  float nAntiLambda = 0.f;
  float purLambda = 0.f;
  float purAntiLambda = 0.f;
  float nEvt = 0.f;
  float nProton = 0.f;
  float nAntiProton = 0.f;

  TFile* _file0=TFile::Open(file);
  TList *listSP=0;
  TList *listPID=0;
  TList *listTP=0;
  TList *listAliEvent=0;
  _file0->GetObject("/PWGCF_PLFemto_0/SPdir_0",listSP);
  _file0->GetObject("/PWGCF_PLFemto_0/PIDdir_0",listPID);
  _file0->GetObject("/PWGCF_PLFemto_0/TP_dir_0",listTP);
  _file0->GetObject("/PWGCF_PLFemto_0/AliEventCuts_0",listAliEvent);

  // Event variables
  auto* histCutStats = (TH1F*)listAliEvent->FindObject("fCutStats");
  auto* histZvertex = (TH1F*)listAliEvent->FindObject("Vtz_selected");
  auto* histNtracklets = (TH1F*)listSP->FindObject("fNTracklets");

  // PID variables
  auto* protonTPC = (TH2F*)listPID->FindObject("fProtonNSigmaTPC");
  auto* protonCombined = (TH2F*)listPID->FindObject("fProtonNSigmaCombined");

  // Single particle plots
  auto* protonDCAxy = (TH1F*)listSP->FindObject("fProtonDCAxy");
  auto* protonDCAxyCutZ = (TH1F*)listSP->FindObject("fProtonDCAxyCutz");
  auto* protonPt = (TH1F*)listSP->FindObject("fProtonPt");
  nProton = protonPt->GetEntries();
  auto* antiprotonPt = (TH1F*)listSP->FindObject("fAntiProtonPt");
  nAntiProton = antiprotonPt->GetEntries();
  auto* protonPhi = (TH1F*)listSP->FindObject("fProtonPhi");
  auto* protonEta = (TH1F*)listSP->FindObject("fProtonEta");
  auto* lambdaPt = (TH1F*)listSP->FindObject("fLambdaPt");
  auto* lambdaPhi = (TH1F*)listSP->FindObject("fLambdaPhi");
  auto* lambdaEta = (TH1F*)listSP->FindObject("fLambdaEta");
  auto* lambdaDaugherDCA = (TH1F*)listSP->FindObject("fLambdaDCADaughterTracks");
  auto* lambdaPosDaughterDCA = (TH1F*)listSP->FindObject("fLambdaDCAPosdaughPrimVertex");
  auto* lambdaNegDaughterDCA = (TH1F*)listSP->FindObject("fLambdaDCANegdaughPrimVertex");
  auto* lambdaTransRadius = (TH1F*)listSP->FindObject("fLambdaTransverseRadius");
  auto* lambdaInvMass = (TH1F*)listSP->FindObject("fInvMassLambdawCuts");
  auto* antilambdaInvMass = (TH1F*)listSP->FindObject("fInvMassAntiLambdawCuts");
  auto* lambdaK0InvMass = (TH1F*)listSP->FindObject("fInvMassMissIDK0s");
  auto* lambdaK0InvMassCut = (TH1F*)listSP->FindObject("fInvMassMissIDK0swCuts");

  // Particle pair plots
  auto* V0V0shared = (TH1F*)listSP->FindObject("fNV0TrackSharing");
  auto* aV0aV0shared = (TH1F*)listSP->FindObject("fNAntiV0TrackSharing");
  auto* V0pshared = (TH1F*)listSP->FindObject("fNV0protonSharedTracks");
  auto* aV0apshared = (TH1F*)listSP->FindObject("fNAntiV0AntiprotonSharedTracks");
  auto* Xipshared = (TH1F*)listSP->FindObject("fNXiSharedTracks");
  auto* aXiapshared = (TH1F*)listSP->FindObject("fNAntiXiAntiprotonSharedTracks");

  auto *cCutStats = new TCanvas();
  SetStyleHisto(histCutStats, 0, 1);
  histCutStats->Draw();
  nTriggers = histCutStats->GetBinContent(1);
  nEvt = histCutStats->GetBinContent(16);
  cCutStats->Print("ANplot/CutStats.pdf");

  auto *cEventProps= new TCanvas();
  cEventProps->Divide(2,1);
  cEventProps->cd(1);
  cEventProps->cd(1)->SetLogy();
  SetStyleHisto(histNtracklets, 0, 1);
  histNtracklets->SetTitle("; Number of SPD tracklets in |#eta|<0.8; Entries");
  histNtracklets->Draw("hist");
  cEventProps->cd(2);
  SetStyleHisto(histZvertex, 0, 1);
  histZvertex->Draw("hist");
  cEventProps->Print("ANplot/eventProperties.pdf");

  auto *cPIDproton= new TCanvas();
  SetStyleHisto(protonTPC);
  protonTPC->SetTitle("; #it{p} (GeV/#it{c}); |n#sigma_{TPC}|");
  SetStyleHisto(protonCombined);
  protonCombined->SetTitle("; #it{p} (GeV/#it{c}); n#sigma_{comb} = #sqrt{n#sigma_{TPC}^{2} + n#sigma_{TOF}^{2}}");
  cPIDproton->Divide(2,1);
  cPIDproton->cd(1);
  protonTPC->Draw("colz");
  auto *cutLineA = new TLine(0, 3, 0.75,3);
  cutLineA->SetLineColor(kGray);
  cutLineA->SetLineWidth(2);
  cutLineA->SetLineStyle(2);
  cutLineA->Draw("same");
  cPIDproton->cd(2);
  protonCombined->Draw("colz");
  auto *cutLineB = new TLine(0.75, 3, 4.05,3);
  cutLineB->SetLineColor(kGray);
  cutLineB->SetLineWidth(2);
  cutLineB->SetLineStyle(2);
  cutLineB->Draw("same");
  cPIDproton->Print("ANplot/ProtonPID.pdf");

  auto *cQAproton= new TCanvas("cQAproton", "cQAproton", 1250,1000);
  SetStyleHisto(protonPt, 0, 1);
  protonPt->SetTitle("; #it{p}_{T} (GeV/#it{c}); Entries");
  SetStyleHisto(protonPhi, 0, 1);
  protonPhi->GetYaxis()->SetRangeUser(0,1.1*protonPhi->GetMaximum());
  protonPhi->SetTitle("; #phi (rad); Entries");
  SetStyleHisto(protonEta, 0, 1);
  protonEta->GetXaxis()->SetRangeUser(-1.0,1.0);
  protonEta->SetTitle("; #eta; Entries");
  SetStyleHisto(protonDCAxy, 0, 1);
  protonDCAxy->SetTitle("; DCA_{xy} (cm); Entries");
  SetStyleHisto(protonDCAxyCutZ, 0, 2);
  cQAproton->Divide(2,2);
  cQAproton->cd(1);
  cQAproton->cd(1)->SetLogy();
  protonPt->Draw("hist");
  cQAproton->cd(2);
  protonPhi->Draw("hist");
  cQAproton->cd(3);
  protonEta->Draw("hist");
  cQAproton->cd(4);
  cQAproton->cd(4)->SetLogy();
  protonDCAxy->Draw("hist");
  protonDCAxyCutZ->Draw("hist same");
  auto *cutLineC = new TLine(-0.1, 0, -0.1, protonDCAxy->GetMaximum());
  cutLineC->SetLineColor(kGray+3);
  cutLineC->SetLineWidth(2);
  cutLineC->SetLineStyle(2);
  cutLineC->Draw("same");
  auto *cutLineD = new TLine(0.1, 0, 0.1, protonDCAxy->GetMaximum());
  cutLineD->SetLineColor(kGray+3);
  cutLineD->SetLineWidth(2);
  cutLineD->SetLineStyle(2);
  cutLineD->Draw("same");
  cQAproton->Print("ANplot/ProtonQA.pdf");

  auto *cQAlambda = new TCanvas("cQAlambda", "cQAlambda", 1250,1000);
  SetStyleHisto(lambdaPt, 0, 1);
  lambdaPt->SetTitle("; #it{p}_{T} (GeV/#it{c}); Entries");
  SetStyleHisto(lambdaPhi, 0, 1);
  lambdaPhi->GetYaxis()->SetRangeUser(0,1.1*lambdaPhi->GetMaximum());
  lambdaPhi->SetTitle("; #phi (rad); Entries");
  SetStyleHisto(lambdaEta, 0, 1);
  lambdaEta->GetXaxis()->SetRangeUser(-1.0,1.0);
  lambdaEta->SetTitle("; #eta; Entries");
  cQAlambda->Divide(2,2);
  cQAlambda->cd(1);
  cQAlambda->cd(1)->SetLogy();
  lambdaPt->Draw("hist");
  cQAlambda->cd(2);
  lambdaPhi->Draw("hist");
  cQAlambda->cd(3);
  lambdaEta->Draw("hist");
  cQAlambda->Print("ANplot/LambdaQA.pdf");

  auto *cQAlambdaDaughter = new TCanvas("cQAlambdaDaughter", "cQAlambdaDaughter", 1250,1000);
  SetStyleHisto(lambdaDaugherDCA, 0, 1);
  lambdaDaugherDCA->SetTitle("; DCA(p, #pi) (cm); Entries");
  SetStyleHisto(lambdaPosDaughterDCA, 0, 1);
  lambdaPosDaughterDCA->SetTitle("; DCA(p, PV) (cm); Entries");
  SetStyleHisto(lambdaNegDaughterDCA, 0, 1);
  lambdaNegDaughterDCA->SetTitle("; DCA(#pi, PV) (cm); Entries");
  SetStyleHisto(lambdaTransRadius, 0,1);
  lambdaTransRadius->SetTitle("; #it{r}_{xy} (cm); Entries");
  cQAlambdaDaughter->Divide(2,2);
  cQAlambdaDaughter->cd(1);
  cQAlambdaDaughter->cd(1)->SetLogy();
  lambdaDaugherDCA->Draw("hist");
  cQAlambdaDaughter->cd(2);
  cQAlambdaDaughter->cd(2)->SetLogy();
  lambdaPosDaughterDCA->Draw("hist");
  cQAlambdaDaughter->cd(3);
  cQAlambdaDaughter->cd(3)->SetLogy();
  lambdaNegDaughterDCA->Draw("hist");
  cQAlambdaDaughter->cd(4);
  lambdaTransRadius->Draw("hist");
  cQAlambdaDaughter->Print("ANplot/LambdaDaughterQA.pdf");
  auto* cLambdaMass = new TCanvas();
  SetStyleHisto(lambdaInvMass, 0, 1);
  SetStyleHisto(antilambdaInvMass, 0, 1);
  lambdaInvMass->SetTitle("; M_{p#pi^{-}} (GeV/#it{c}^{2}); Entries");
  antilambdaInvMass->SetTitle("; M_{#bar{p}#pi^{+}} (GeV/#it{c}^{2}); Entries");
  SetStyleHisto(lambdaK0InvMass, 0, 1);
  lambdaK0InvMass->SetTitle("; M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); Entries");
  SetStyleHisto(lambdaK0InvMassCut, 0, 2);
  cLambdaMass->Divide(2,1);
  cLambdaMass->cd(1);
  lambdaInvMass->Draw("hist");
  lambdaInvMass->GetXaxis()->SetRangeUser(1.1, 1.13);
  // antilambdaInvMass->Draw("hist same");
  cLambdaMass->cd(2);
  lambdaK0InvMass->Draw("hist");
  lambdaK0InvMass->GetXaxis()->SetRangeUser(0.4, 0.6);
  lambdaK0InvMassCut->Draw("hist same");
  cLambdaMass->Print("ANplot/LambdaInvMass.pdf");

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // LAMBDA

  std::vector<float> ptBinLambdaRange = {0.3, 0.8, 1.3, 1.8, 2.3, 2.8, 3.3, 3.8, 4.3};
  const int nPtBinsLambda = ptBinLambdaRange.size();
  const float marginLambda = 0.004;
  const float massLambda = 1.115;

  auto* cLambda = new TCanvas();
  auto* hRecoLambdaM = (TH1F*)lambdaInvMass->Clone();
  hRecoLambdaM->Draw("PE");
  float lambdaSignalAll, lambdaSignalAllErr, lambdaBackgroundAll, lambdaBackgroundAllErr;
  TList *funListLambda;
  FitLambda(hRecoLambdaM, lambdaSignalAll, lambdaSignalAllErr, lambdaBackgroundAll, lambdaBackgroundAllErr, massLambda-marginLambda, massLambda+marginLambda);
  std::cout << "Lambda \n";
  std::cout << "Signal " << lambdaSignalAll << " Background " << lambdaBackgroundAll << " S/B " << lambdaSignalAll/lambdaBackgroundAll << " Purity " << lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f << "\n";
  nLambda = lambdaSignalAll;
  purLambda = lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f;
  hRecoLambdaM->GetXaxis()->SetRangeUser(1.1, 1.13);
  funListLambda = hRecoLambdaM->GetListOfFunctions();
  TF1 *fLambdaTotal = (TF1*)funListLambda->FindObject("fLambda");

  float amp1 = fLambdaTotal->GetParameter(3);
  float amp2 = fLambdaTotal->GetParameter(6);
  float mean1 = fLambdaTotal->GetParameter(4);
  float mean2 = fLambdaTotal->GetParameter(7);
  float width1 = fLambdaTotal->GetParameter(5);
  float width2 = fLambdaTotal->GetParameter(8);

  float meanMass = weightedMean(amp1, mean1, amp2, mean2);
  float meanWidth = weightedMean(amp1, width1, amp2, width2);
  TLatex LambdaLabel;
  LambdaLabel.SetNDC(kTRUE);
  LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
  LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,
                        Form("#splitline{#splitline{#splitline{#Lambda: %.0f}{m_{#Lambda} = %.1f (MeV/#it{c}^{2})}}{#sigma_{#Lambda} = %.1f (MeV/#it{c}^{2})}}{Purity = %.1f %%}", lambdaSignalAll, meanMass*1000.f, meanWidth*1000.f, lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f));
  TLatex BeamTextLambda;
  BeamTextLambda.SetNDC(kTRUE);
  BeamTextLambda.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.2, "pp #sqrt{#it{s}} = 13 TeV");
  cLambda->Print("ANplot/InvMassLambda.pdf");


  auto* cAntiLambda = new TCanvas();
  auto* hRecoAntiLambdaM = (TH1F*)antilambdaInvMass->Clone();
  hRecoAntiLambdaM->Draw("PE");
  FitLambda(hRecoAntiLambdaM, lambdaSignalAll, lambdaSignalAllErr, lambdaBackgroundAll, lambdaBackgroundAllErr, massLambda-marginLambda, massLambda+marginLambda);
  std::cout << "Lambda \n";
  std::cout << "Signal " << lambdaSignalAll << " Background " << lambdaBackgroundAll << " S/B " << lambdaSignalAll/lambdaBackgroundAll << " Purity " << lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f << "\n";
  nAntiLambda = lambdaSignalAll;
  purAntiLambda = lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f;
  hRecoAntiLambdaM->GetXaxis()->SetRangeUser(1.1, 1.13);
  funListLambda = hRecoAntiLambdaM->GetListOfFunctions();
  fLambdaTotal = (TF1*)funListLambda->FindObject("fLambda");

  amp1 = fLambdaTotal->GetParameter(3);
  amp2 = fLambdaTotal->GetParameter(6);
  mean1 = fLambdaTotal->GetParameter(4);
  mean2 = fLambdaTotal->GetParameter(7);
  width1 = fLambdaTotal->GetParameter(5);
  width2 = fLambdaTotal->GetParameter(8);

  meanMass = weightedMean(amp1, mean1, amp2, mean2);
  meanWidth = weightedMean(amp1, width1, amp2, width2);
  LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,
                        Form("#splitline{#splitline{#splitline{#bar{#Lambda}: %.0f}{m_{#bar{#Lambda}} = %.1f (MeV/#it{c}^{2})}}{#sigma_{#bar{#Lambda}} = %.1f (MeV/#it{c}^{2})}}{Purity = %.1f %%}", lambdaSignalAll, meanMass*1000.f, meanWidth*1000.f, lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f));
  BeamTextLambda.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.2, "pp #sqrt{#it{s}} = 13 TeV");
  cAntiLambda->Print("ANplot/InvMassAntiLambda.pdf");

  TH1F *hRecoLambdaPt[nPtBinsLambda-1];
  TGraphErrors *grLambdaMass = new TGraphErrors();
  TGraphErrors *grLambdaWidth = new TGraphErrors();
  TGraphErrors *grLambdaPurity = new TGraphErrors();
  TGraphErrors *grLambdaYield = new TGraphErrors();

  TCanvas *cLambdaPt = new TCanvas("cLambdaPt", "", 1500, 1000);
  cLambdaPt->Divide(4,2);
  TList *funListLambdaPt[nPtBinsLambda];
  float lambdaSignal[nPtBinsLambda];
  float lambdaSignalErr[nPtBinsLambda];
  float lambdaBackground[nPtBinsLambda];
  float lambdaBackgroundErr[nPtBinsLambda];
  float nLambdaPt[nPtBinsLambda];
  float nLambdaPtErr[nPtBinsLambda];

  for(int i=0; i<nPtBinsLambda-1; ++i) {
    cLambdaPt->cd(i+1);
    hRecoLambdaPt[i] = (TH1F*)listSP->FindObject(Form("fInvMassLambdawCutsPtBin%i", i));
    hRecoLambdaPt[i]->SetTitle(Form("%.2f < #it{p}_{T} < %.2f GeV/#it{c}; M_{p#pi^{-}} (GeV/#it{c}^{2}); Entries", ptBinLambdaRange[i], ptBinLambdaRange[i+1]));
    SetStyleHistoMany(hRecoLambdaPt[i], 0, 0);
    hRecoLambdaPt[i]->Draw();
    FitLambda(hRecoLambdaPt[i], lambdaSignal[i], lambdaSignalErr[i], lambdaBackground[i], lambdaBackgroundErr[i], massLambda-marginLambda, massLambda+marginLambda);
    hRecoLambdaPt[i]->GetXaxis()->SetRangeUser(1.1, 1.13);
    funListLambdaPt[i] = hRecoLambdaPt[i]->GetListOfFunctions();
    TF1 *fLambdaPtTotal = (TF1*)funListLambdaPt[i]->FindObject("fLambda");

    float amp1Err = fLambdaPtTotal->GetParError(3);
    float amp2Err = fLambdaPtTotal->GetParError(6);
    float mean1Err = fLambdaPtTotal->GetParError(4);
    float mean2Err = fLambdaPtTotal->GetParError(7);
    float width1Err = fLambdaPtTotal->GetParError(5);
    float width2Err = fLambdaPtTotal->GetParError(8);

    float amp1 = fLambdaPtTotal->GetParameter(3);
    float amp2 = fLambdaPtTotal->GetParameter(6);
    float mean1 = fLambdaPtTotal->GetParameter(4);
    float mean2 = fLambdaPtTotal->GetParameter(7);
    float width1 = fLambdaPtTotal->GetParameter(5);
    float width2 = fLambdaPtTotal->GetParameter(8);

    float meanMass = weightedMean(amp1, mean1, amp2, mean2);
    float meanMassErr = weightedMeanError(amp1, mean1, amp2, mean2, amp1Err, mean1Err, amp2Err, mean2Err);
    float meanWidth = weightedMean(amp1, width1, amp2, width2);
    float meanWidthErr = weightedMeanError(amp1, width1, amp2, width2, amp1Err, width1Err, amp2Err, width2Err);

    lambdaSignal[i] += lambdaBackground[i];
    nLambdaPtErr[i] = lambdaSignalErr[i];
    lambdaSignalErr[i] = std::sqrt(lambdaSignalErr[i]*lambdaSignalErr[i] + lambdaBackgroundErr[i]*lambdaBackgroundErr[i]);
    nLambdaPt[i] = lambdaSignal[i] - lambdaBackground[i];

    grLambdaMass->SetPoint(i, (ptBinLambdaRange[i+1]-ptBinLambdaRange[i])/2.f + ptBinLambdaRange[i],  meanMass);
    grLambdaMass->SetPointError(i, (ptBinLambdaRange[i+1]-ptBinLambdaRange[i])/2.f, meanMassErr);
    grLambdaWidth->SetPoint(i, (ptBinLambdaRange[i+1]-ptBinLambdaRange[i])/2.f + ptBinLambdaRange[i],  meanWidth);
    grLambdaWidth->SetPointError(i, (ptBinLambdaRange[i+1]-ptBinLambdaRange[i])/2.f, meanWidthErr);
    grLambdaPurity->SetPoint(i, (ptBinLambdaRange[i+1]-ptBinLambdaRange[i])/2.f + ptBinLambdaRange[i], nLambdaPt[i]/lambdaSignal[i] *100.f);
    grLambdaPurity->SetPointError(i, (ptBinLambdaRange[i+1]-ptBinLambdaRange[i])/2.f, 100.f * std::sqrt(nLambdaPtErr[i]*nLambdaPtErr[i]/(lambdaSignal[i]*lambdaSignal[i]) + nLambdaPt[i]*nLambdaPt[i]*lambdaSignalErr[i]*lambdaSignalErr[i]/std::pow(lambdaSignal[i], 4)));
    grLambdaYield->SetPoint(i, (ptBinLambdaRange[i+1]-ptBinLambdaRange[i])/2.f + ptBinLambdaRange[i], nLambdaPt[i]);
    grLambdaYield->SetPointError(i, (ptBinLambdaRange[i+1]-ptBinLambdaRange[i])/2.f, nLambdaPtErr[i]);
  }

  auto *cLambdaProperties = new TCanvas();
  cLambdaProperties->Divide(2,2);
  cLambdaProperties->cd(1);
  SetStyleGraph(grLambdaMass , 0, 0);
  grLambdaMass->SetTitle("; #it{p}_{T} (GeV/#it{c}); M_{#Lambda} (GeV/#it{c^{2}})");
  grLambdaMass->Draw("APEZ");

  cLambdaProperties->cd(2);
  SetStyleGraph(grLambdaWidth , 0, 0);
  grLambdaWidth->SetTitle("; #it{p}_{T} (GeV/#it{c}); #sigma(M_{#Lambda}) (GeV/#it{c^{2}})");
  grLambdaWidth->Draw("APEZ");

  cLambdaProperties->cd(3);
  SetStyleGraph(grLambdaPurity , 0, 0);
  grLambdaPurity->SetTitle("; #it{p}_{T} (GeV/#it{c}); Purity (%)");
  grLambdaPurity->Draw("APEZ");
  TF1 *fLambdaPurity = new TF1("fLambdaPurity", "pol0", 0.3, 4.3);
  fLambdaPurity->SetLineColor(fColors[0]);
  fLambdaPurity->SetLineStyle(2);
  //  fLambdaPurity->SetParameter(0, computeWeightedAverage(hRecoLambdaPtPurity, grLambdaPurity, 0.3, 4.3));
  fLambdaPurity->Draw("same");
  TLatex purity;
  purity.SetTextFont(43);
  purity.SetTextSize(21);
  purity.SetNDC(kTRUE);
  purity.DrawLatex(0.65, 0.25, Form("Purity (%.1f %%)", fLambdaPurity->GetParameter(0)));


  auto* cShared = new TCanvas("cShared", "cShared", 1250,1000);
  SetStyleHisto(V0V0shared, 0, 1);
  V0V0shared->SetTitle("; # V0 pairs with shared tracks / event; Entries");
  SetStyleHisto(aV0aV0shared, 0, 1);
  aV0aV0shared->SetTitle("; # #bar{V0} pairs with shared tracks / event; Entries");
  SetStyleHisto(V0pshared, 0, 1);
  V0pshared->SetTitle("; # p-V0 pairs with shared tracks / event; Entries");
  SetStyleHisto(aV0apshared, 0,1);
  aV0apshared->SetTitle("; # #bar{p}-#bar{V0} pairs with shared tracks / event; Entries");
  cShared->Divide(2,2);
  cShared->cd(1);
  cShared->cd(1)->SetLogy();
  V0V0shared->Draw("hist");
  cShared->cd(2);
  cShared->cd(2)->SetLogy();
  aV0aV0shared->Draw("hist");
  cShared->cd(3);
  cShared->cd(3)->SetLogy();
  V0pshared->Draw("hist");
  cShared->cd(4);
  cShared->cd(4)->SetLogy();
  aV0apshared->Draw("hist");
  cShared->Print("ANplot/SharedTracks.pdf");
  if (Xipshared) {
    auto *cSharedpXi=new TCanvas("cSharedXi", "cSharedXi", 1250,1000);
    cSharedpXi->Divide(2,1);
    SetStyleHisto(Xipshared, 0, 1);
    Xipshared->SetTitle("; # p-#Xi pairs with shared tracks; Entries");
    SetStyleHisto(aXiapshared, 0, 1);
    aXiapshared->SetTitle("; # #bar{p}-#bar{#Xi} pairs with shared tracks; Entries");
    cSharedpXi->cd(1);
    cSharedpXi->cd(1)->SetLogy();
    Xipshared->Draw("hist");
    cSharedpXi->cd(2);
    cSharedpXi->cd(2)->SetLogy();
    aXiapshared->Draw("hist");
    cSharedpXi->Print("ANplot/SharedTracksXi.pdf");
  }
  // OUTPUT
  std::cout << "==========================\n";
  std::cout << "kINT7 triggers \t" << Form("%.1f", nTriggers/(1E6)) << " x E6 \n";
  std::cout << "nEvts \t\t" << Form("%.1f", nEvt/(1E6)) << " x E6 \n";
  std::cout << "nLambda \t" << nLambda << "\tPurity " << purLambda << "\tFraction " << nLambda/nEvt << "\tfrac err " << std::sqrt(nLambda/(nEvt*nEvt) + nLambda*nLambda/(nEvt*nEvt*nEvt)) << "\n";
  std::cout << "nAntiLambda \t" << nAntiLambda << "\tPurity " << purAntiLambda << "\tFraction " << nAntiLambda/nEvt << "\tfrac err " << std::sqrt(nAntiLambda/(nEvt*nEvt) + nAntiLambda*nAntiLambda/(nEvt*nEvt*nEvt)) << "\n";
  std::cout << "nProton \t" << nProton << "\tFraction " << nProton/nEvt << "\tfrac err " << std::sqrt(nProton/(nEvt*nEvt) + nProton*nProton/(nEvt*nEvt*nEvt)) << "\n";
  std::cout << "nAntiProton \t" << nAntiProton << "\tFraction " << nAntiProton/nEvt << "\tfrac err " << std::sqrt(nAntiProton/(nEvt*nEvt) + nAntiProton*nAntiProton/(nEvt*nEvt*nEvt)) << "\n";
}
