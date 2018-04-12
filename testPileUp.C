
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
//  histRE_relK->Sumw2();
//  histME_relK->Sumw2();
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
//  result->Sumw2();
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
void testPileUp(const char *expfile = "~/Downloads/AnalysisResults.root", const char* prefix ="MB", const char* folder = "0") {

  SetStyle();

  TFile* _file0=TFile::Open(expfile);

  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;
  const float right = 0.025;
  const float top = 0.025;

  TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sResults%s", prefix, folder)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, folder), Results);
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

  const float marginLambda = 0.004;
  const float massLambda = 1.115;

  TDirectoryFile *dirv0Cuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sv0Cuts%s", prefix, folder)));
  TList *tmp=dirv0Cuts->GetListOfKeys();
  TString name=tmp->At(0)->GetName();
  TList *v0CutList;
  dirv0Cuts->GetObject(name,v0CutList);
  TList *v0Cuts = (TList*)v0CutList->FindObject("MinimalBooking");
  if(!v0Cuts) v0Cuts = (TList*)v0CutList->FindObject("v0Cuts");
  TH1F *InvMassPtLambda= (TH1F*)((TH2F*)v0Cuts->FindObject("InvMassPt"))->ProjectionY();
  SetStyleHisto(InvMassPtLambda, 0,1);
  float lambdaSignalAll, lambdaSignalAllErr, lambdaBackgroundAll, lambdaBackgroundAllErr;
  FitLambda(InvMassPtLambda, lambdaSignalAll, lambdaSignalAllErr, lambdaBackgroundAll, lambdaBackgroundAllErr, massLambda-marginLambda, massLambda+marginLambda);
  const int LambdaSignal = lambdaSignalAll;
  const int LambdaBackground = lambdaBackgroundAll;
  const float LambdaPurity = lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f;

  std::cout << "Lambda \n";
  std::cout << "Signal " << lambdaSignalAll << " Background " << lambdaBackgroundAll << " S/B " << lambdaSignalAll/lambdaBackgroundAll << " Purity " << lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f << "\n";

  TDirectoryFile *dirAv0Cuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sAntiv0Cuts%s", prefix, folder)));
  tmp=dirAv0Cuts->GetListOfKeys();
  name=tmp->At(0)->GetName();
  TList *Av0CutList;
  dirAv0Cuts->GetObject(name,Av0CutList);
  TList *Av0Cuts = (TList*)Av0CutList->FindObject("MinimalBooking");
  if(!Av0Cuts) Av0Cuts = (TList*)Av0CutList->FindObject("v0Cuts");
  TH1F *InvMassPtAntiLambda= (TH1F*)((TH2F*)Av0Cuts->FindObject("InvMassPt"))->ProjectionY();
  SetStyleHisto(InvMassPtAntiLambda, 0,1);
  FitLambda(InvMassPtAntiLambda, lambdaSignalAll, lambdaSignalAllErr, lambdaBackgroundAll, lambdaBackgroundAllErr, massLambda-marginLambda, massLambda+marginLambda);
  std::cout << "Anti-Lambda \n";
  std::cout << "Signal " << lambdaSignalAll << " Background " << lambdaBackgroundAll << " S/B " << lambdaSignalAll/lambdaBackgroundAll << " Purity " << lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f << "\n";
  const int AntiLambdaSignal = lambdaSignalAll;
  const int AntiLambdaBackground = lambdaBackgroundAll;
  const float AntiLambdaPurity = lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f;

  TH2F* xiMass2D;
  TH1F* xiMass;
  TH2F* AntiXiMass2D;
  TH1F* AntiXiMass;
  TDirectoryFile *dirExp=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sCascadeCuts%s", prefix, folder)));
  tmp=dirExp->GetListOfKeys();
  name=tmp->At(0)->GetName();
  TList *output;
  dirExp->GetObject(name,output);
  TList *Cascade=(TList*)output->FindObject("MinimalBooking");
  if(!Cascade) Cascade=(TList*)output->FindObject("Cascade");
  xiMass2D=(TH2F*)Cascade->FindObject("InvMassXiPt");
  xiMass=(TH1F*)xiMass2D->ProjectionY("InvMassXi");
  SetStyleHisto(xiMass,0,1);
  xiMass->GetXaxis()->SetRangeUser(1.285,1.345);
  xiMass->GetXaxis()->SetTitle("M_{#pi#Lambda} (GeV/#it{c}^{2})");

  float signalXi, signalErrXi, backgroundXi, backgroundErrXi;
  float lowerBound=1.322-0.005;
  float upperBound=1.322+0.005;
  FitXi(xiMass, signalXi, signalErrXi, backgroundXi, backgroundErrXi, lowerBound,upperBound);
  std::cout << "Xi\n";
  std::cout << "Signal " << signalXi << " Background " << backgroundXi << " S/B " << signalXi/backgroundXi << " Purity " << signalXi/(signalXi+backgroundXi)*100.f << "\n";
  const int XiSignal = signalXi;
  const int XiBackground = backgroundXi;
  const float XiPurity = signalXi/(signalXi+backgroundXi)*100.f;

  dirExp=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sAntiCascadeCuts%s", prefix, folder)));
  tmp=dirExp->GetListOfKeys();
  name=tmp->At(0)->GetName();
  dirExp->GetObject(name,output);
  TList *AntiCascade=(TList*)output->FindObject("MinimalBooking");
  if(!AntiCascade) AntiCascade=(TList*)output->FindObject("Cascade");
  AntiXiMass2D=(TH2F*)AntiCascade->FindObject("InvMassXiPt");
  AntiXiMass=(TH1F*)AntiXiMass2D->ProjectionY("InvMassXi");
  SetStyleHisto(AntiXiMass,0,1);
  AntiXiMass->GetXaxis()->SetRangeUser(1.285,1.345);
  AntiXiMass->GetXaxis()->SetTitle("M_{#pi#Lambda} (GeV/#it{c}^{2})");

  float signalAntiXi, signalErrAntiXi, backgroundAntiXi, backgroundErrAntiXi;
  FitXi(AntiXiMass, signalAntiXi, signalErrAntiXi, backgroundAntiXi, backgroundErrAntiXi, lowerBound,upperBound);
  std::cout << "Anti-Xi\n";
  std::cout << "Signal " << signalAntiXi << " Background " << backgroundAntiXi << " S/B " << signalAntiXi/backgroundAntiXi << " Purity " << signalAntiXi/(signalAntiXi+backgroundAntiXi)*100.f << "\n";
  const int AntiXiSignal = signalAntiXi;
  const int AntiXiBackground = backgroundAntiXi;
  const float AntiXiPurity = signalAntiXi/(signalAntiXi+backgroundAntiXi)*100.f;

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
  hist_CF_pXi_ApAXi_exp[0]->GetYaxis()->SetRangeUser(0.5, 7.5);
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
  Can_CF->Print(Form("PileUp/CF_pp-apap_%s.pdf", folder));

  TCanvas *Can_CF_comb = new TCanvas("Can_CF_comb","Can_CF_comb",0,0,1000,550);
  Can_CF_comb->Divide(4,1);
  Can_CF_comb->cd(1);
  Can_CF_comb->cd(1)->SetRightMargin(right);
  Can_CF_comb->cd(1)->SetTopMargin(top);
  hist_CF_pp_ApAp_exp[2]->Draw("pe");
  hist_CF_pp_ApAp_exp[2]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_pp_ApAp_exp[2]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_pp_ApAp_exp[2]->GetYaxis()->SetRangeUser(0, 4);
  Can_CF_comb->cd(2);
  Can_CF_comb->cd(2)->SetRightMargin(right);
  Can_CF_comb->cd(2)->SetTopMargin(top);
  hist_CF_Lp_ALAp_exp[2]->Draw("pe");
  hist_CF_Lp_ALAp_exp[2]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_Lp_ALAp_exp[2]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_Lp_ALAp_exp[2]->GetXaxis()->SetNdivisions(505);
  hist_CF_Lp_ALAp_exp[2]->GetYaxis()->SetRangeUser(0.8, 2.5);
  Can_CF_comb->cd(3);
  Can_CF_comb->cd(3)->SetRightMargin(right);
  Can_CF_comb->cd(3)->SetTopMargin(top);
  hist_CF_LL_ALAL_exp[2]->Draw("pe");
  hist_CF_LL_ALAL_exp[2]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_LL_ALAL_exp[2]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_LL_ALAL_exp[2]->GetXaxis()->SetNdivisions(505);
  hist_CF_LL_ALAL_exp[2]->GetYaxis()->SetRangeUser(0.25, 3);
  Can_CF_comb->cd(4);
  Can_CF_comb->cd(4)->SetRightMargin(right);
  Can_CF_comb->cd(4)->SetTopMargin(top);
  hist_CF_pXi_ApAXi_exp[2]->Draw("pe");
  hist_CF_pXi_ApAXi_exp[2]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  hist_CF_pXi_ApAXi_exp[2]->GetXaxis()->SetRangeUser(0, 0.4);
  hist_CF_pXi_ApAXi_exp[2]->GetXaxis()->SetNdivisions(505);
  hist_CF_pXi_ApAXi_exp[2]->GetYaxis()->SetRangeUser(0.5, 7.5);
  Can_CF_comb->Print(Form("PileUp/CF_%s.pdf", folder));

  TCanvas *cLambda = new TCanvas("cLambda","cLambda",0,0,1000,550);
  cLambda->Divide(2,1);
  cLambda->cd(1);
  InvMassPtLambda->Draw();
  InvMassPtLambda->GetXaxis()->SetRangeUser(1.1, 1.13);
  TLatex LambdaLabel;
  LambdaLabel.SetNDC(kTRUE);
  LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
  LambdaLabel.DrawLatex(0.2, 0.8, Form("Purity: %.1f %%", LambdaPurity));
  LambdaLabel.DrawLatex(0.2, 0.75, Form("Signal: %.2f #times 10^{6}", float(LambdaSignal)/1000000.f));
  LambdaLabel.DrawLatex(0.2, 0.7, Form("Background: %.2f #times 10^{6}",float(LambdaBackground)/1000000.f));
  cLambda->cd(2);
  InvMassPtAntiLambda->Draw();
  InvMassPtAntiLambda->GetXaxis()->SetRangeUser(1.1, 1.13);
  LambdaLabel.DrawLatex(0.2, 0.8, Form("Purity: %.1f %%", AntiLambdaPurity));
  LambdaLabel.DrawLatex(0.2, 0.75, Form("Signal: %.2f #times 10^{6}", float(AntiLambdaSignal)/1000000.f));
  LambdaLabel.DrawLatex(0.2, 0.7, Form("Background: %.2f #times 10^{6}",float(AntiLambdaBackground)/1000000.f));
  cLambda->Print(Form("PileUp/Lambda_%s.pdf", folder));


  TCanvas *Xi= new TCanvas("cXi","cXi",0,0,1000,550);
  Xi->Divide(2,1);
  Xi->cd(1);
  xiMass->Draw();
  LambdaLabel.SetNDC(kTRUE);
  LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
  LambdaLabel.DrawLatex(0.2, 0.8, Form("Purity: %.1f %%", XiPurity));
  LambdaLabel.DrawLatex(0.2, 0.75, Form("Signal: %.2f #times 10^{3}", float(XiSignal)/1000.f));
  LambdaLabel.DrawLatex(0.2, 0.7, Form("Background: %.2f #times 10^{3}",float(XiBackground)/1000.f));
  Xi->cd(2);
  AntiXiMass->Draw();
  LambdaLabel.SetNDC(kTRUE);
  LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
  LambdaLabel.DrawLatex(0.2, 0.8, Form("Purity: %.1f %%", AntiXiPurity));
  LambdaLabel.DrawLatex(0.2, 0.75, Form("Signal: %.2f #times 10^{3}", float(AntiXiSignal)/1000.f));
  LambdaLabel.DrawLatex(0.2, 0.7, Form("Background: %.2f #times 10^{3}",float(AntiXiBackground)/1000.f));
  Xi->Print(Form("PileUp/Xi_%s.pdf", folder));
}
