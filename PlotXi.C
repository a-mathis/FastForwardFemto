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
  histo->GetXaxis()->SetLabelSize(0.06);
  histo->GetXaxis()->SetTitleSize(0.07);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.0);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.06);
  histo->GetYaxis()->SetTitleSize(0.07);
  histo->GetYaxis()->SetLabelOffset(0.02);
  histo->GetYaxis()->SetTitleOffset(1.0);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// second order polynomial + double gaus to Lambda peak
void FitTheLambda(TH1F* histo, float &signal, float
               &signalErr, float &background, float &backgroundErr, float lowerBound,
               float upperBound, TCanvas *c1)
{
  c1->cd();
  //A solid fit.
  double xValue1 = 1.108;
  double yValue1 = histo->GetBinContent(histo->FindBin(xValue1));
  double xValue2 = 1.124;
  double yValue2 = histo->GetBinContent(histo->FindBin(xValue2));

  double a1 = (yValue2-yValue1)/(xValue2-xValue1);
  double a0 = yValue1 - a1*xValue1;
  // parse then to proper TF1
  TF1 *theBackground = new TF1("theBackground","pol1", 0, 1.4);
  theBackground->SetParameter(0, a0);
  theBackground->SetParameter(1, a1);

  // remove background from signal
  TH1F *signalOnly = getSignalHisto(theBackground, histo, xValue1, xValue2,
                                    Form("%s_Thesignal_only", histo->GetName()));
  // fit signal only
  TF1 *fTheSignalSingleGauss = new TF1("fTheSignalSingleGauss", "gaus(0)",0.8*xValue1,1.2*xValue2);
  //  signalOnly->DrawCopy();
  signalOnly->Fit("fTheSignalSingleGauss");
  TF1 *fTheSignalGauss = new TF1("fTheSignalGauss", "gaus(0) + gaus(3)", 0.8*xValue1,1.2*xValue2);
  fTheSignalGauss->SetParameter(0, 0.75 * histo->GetMaximum());
  fTheSignalGauss->SetParameter(1, fTheSignalSingleGauss->GetParameter(1));
  fTheSignalGauss->SetParameter(2, 2.f*fTheSignalSingleGauss->GetParameter(2));
  fTheSignalGauss->SetParLimits(2, 0.5*fTheSignalSingleGauss->GetParameter(2),1e2*2.f*fTheSignalSingleGauss->GetParameter(2));
  fTheSignalGauss->SetParameter(3, 0.2 * histo->GetMaximum());
  fTheSignalGauss->SetParameter(4, fTheSignalSingleGauss->GetParameter(1));
  fTheSignalGauss->SetParLimits(4,
                             fTheSignalSingleGauss->GetParameter(1)-fTheSignalSingleGauss->GetParameter(2),
                             fTheSignalSingleGauss->GetParameter(1)+fTheSignalSingleGauss->GetParameter(2));
  fTheSignalGauss->SetParameter(5, 0.5*fTheSignalSingleGauss->GetParameter(2));
  fTheSignalGauss->SetParLimits(5, 0.5*fTheSignalSingleGauss->GetParameter(2),1e2*2.f*fTheSignalSingleGauss->GetParameter(2));
  TFitResultPtr r = signalOnly->Fit("fTheSignalGauss", "SRQ0", "",xValue1,xValue2);

  // Extract signal as integral
  signal = fTheSignalGauss->Integral(lowerBound, upperBound)
                /double(histo->GetBinWidth(1));
//  signalErr = fTheSignalGauss->IntegralError(lowerBound, upperBound,
//                                          r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())
//                /double(histo->GetBinWidth(1));
//  c1->cd();
//  histo->GetXaxis()->SetRangeUser(0.8*xValue1,1.2*xValue2);
//  histo->Draw();
  TF1 *fTheLambda = new TF1("fTheLambda", "theBackground + fTheSignalGauss", 0.8*xValue1,1.2*xValue2);
  fTheLambda->SetNpx(1000);
  fTheLambda->SetParameter(2, 0.75 * histo->GetMaximum());
  fTheLambda->SetParameter(3, fTheSignalGauss->GetParameter((1)));
  fTheLambda->SetParameter(4, fTheSignalGauss->GetParameter((2)));
  fTheLambda->SetParameter(5, 0.2 * histo->GetMaximum());
  fTheLambda->SetParameter(6, fTheSignalGauss->GetParameter((4)));
  fTheLambda->SetParameter(7, fTheSignalGauss->GetParameter((5)));
  fTheLambda->SetLineColor(fColors[2]);

  histo->Fit("fTheLambda", "SRQ", "", xValue1,xValue2);

  TF1 *fTheLambda_background = new TF1("fTheLambda_background", "pol1(0)",0.8*xValue1,1.2*xValue2);
  fTheLambda_background->SetParameter(0, fTheLambda->GetParameter(0));
  fTheLambda_background->SetParameter(1, fTheLambda->GetParameter(1));
  fTheLambda_background->SetLineStyle(3);
  fTheLambda_background->SetLineColor(fColors[2]);

  background = fTheLambda_background->Integral(lowerBound, upperBound)
                /double(histo->GetBinWidth(1));
}
// second order polynomial + double gaus to Lambda peak
void FitLambda(TH1F* histo, float &signal, float
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
  signalOnly->Fit("fSignalSingleGauss");
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
  fLambda_background->SetParameter(0, fLambda->GetParameter(0));
  fLambda_background->SetParameter(1, fLambda->GetParameter(1));
  fLambda_background->SetParameter(2, fLambda->GetParameter(2));
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


void PlotXi(TString fileName, const char *prefix) {
  SetStyle(false,false);
  TH1F* xiPt[13];
  TGaxis::SetMaxDigits(2);
  TH2F* xiMass2D;
  TH1F* xiMass;
  TH2F* LambdaMass2D;
  TH1F* LambdaMass;
  TFile *file=TFile::Open(fileName.Data());
  //  file->ls();
  TString dirName="";
  dirName+=prefix;
  dirName+="AntiCascadeCuts";
  TDirectoryFile *dirExp=(TDirectoryFile*)(file->FindObjectAny(dirName.Data()));
  TCanvas *Xi= new TCanvas("Xi","Xi",0,0,1500,1100);
  TCanvas *Lambda= new TCanvas("Lambda","Lambda",0,0,1500,1100);
  TCanvas *PtXi= new TCanvas("PtXi","PtXi",0,0,1500,1100);
  TH1F *Purity=0;
  PtXi->Divide(4,3);
  if (dirExp) {
    TList *tmp=dirExp->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *output;
    dirExp->GetObject(name,output);
    if (!output) {
      std::cout << "No Output \n";
    }
    //output->ls();
    TList *Cascade=(TList*)output->FindObject("Cascade");
    if(!Cascade) {
      std::cout << "No Cascade \n";
    }
    xiMass2D=(TH2F*)Cascade->FindObject("InvMassXi");
    xiMass=(TH1F*)xiMass2D->ProjectionY("InvMassXioneBin");
    LambdaMass2D=(TH2F*)Cascade->FindObject("InvMassv0Pt");
    LambdaMass = (TH1F*)LambdaMass2D->ProjectionY("InvMassLamoneBin");

    Purity=new TH1F("Purity","Purity",14,xiMass2D->GetXaxis()->GetBinLowEdge(1),xiMass2D->GetXaxis()->GetBinUpEdge(13));
    std::cout << xiMass2D->GetXaxis()->GetBinLowEdge(1) << '\t' << xiMass2D->GetXaxis()->GetBinUpEdge(13);
    for (int iBin = 1; iBin<13; ++iBin) {
      TPad *padPt=(TPad*)PtXi->cd(iBin);
      padPt->SetRightMargin(0.0);
      padPt->SetTopMargin(0.065);
      padPt->SetBottomMargin(0.15);
      padPt->SetLeftMargin(0.09);
      padPt->Draw();
      TString XiNamePt=Form("XiPt_%i",iBin);
      xiPt[iBin]=(TH1F*)xiMass2D->ProjectionY(XiNamePt.Data(),iBin+1,iBin+1);
      SetStyleHisto(xiPt[iBin],2,1);
      xiPt[iBin]->GetXaxis()->SetRangeUser(1.285,1.345);
      xiPt[iBin]->GetXaxis()->SetTitle("M_{#pi#Lambda} (GeV/#it{c}^{2})");
      xiPt[iBin]->DrawCopy();
      float signal=0;
      float signalErr=0;
      float background=0;
      float backgroundErr=0;
      float lowerBound=1.322-0.005;
      float upperBound=1.322+0.005;
      FitLambda(xiPt[iBin], signal, signalErr, background, backgroundErr, lowerBound,upperBound);

      TList *funListLambda = xiPt[iBin]->GetListOfFunctions();
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
      LambdaLabel.SetTextSize(gStyle->GetTextSize()*1.2);
      LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.45,
                            Form("#splitline{#splitline{#splitline{#splitline{%.2f < p_{T} < %.2f}{#Xi^{-}: %.0f}}{m_{#Xi} = %.1fMeV/#it{c}^{2}}}{#sigma_{#Xi}= %.1fMeV/#it{c}^{2}}}{Purity = %.1f %%}",
                                 xiMass2D->GetXaxis()->GetBinCenter(iBin+1),xiMass2D->GetXaxis()->GetBinCenter(iBin+2),
                                 signal, meanMass*1000.f, meanWidth*1000.f,
                                 signal/(signal+background)*100.f));//#splitline{#splitline{
      Purity->SetBinContent(iBin+1,signal/(signal+background)*100.f);
    }
    SetStyleHisto(xiMass,20,1);
    xiMass->GetXaxis()->SetRangeUser(1.285,1.345);
    xiMass->GetXaxis()->SetTitle("M_{#pi#Lambda} (GeV/#it{c}^{2})");
    if (!xiMass) {
      std::cout << "Missing Hist" << std::endl;
    }
  } else {
    std::cout << "dirExp missing \n";
  }

  float signal=0;
  float signalErr=0;
  float background=0;
  float backgroundErr=0;
  float lowerBound=1.322-0.005;
  float upperBound=1.322+0.005;

//  xiMass->DrawCopy();
//  SetStyle(false,false);
  Xi->cd();
  FitLambda(xiMass, signal, signalErr, background, backgroundErr, lowerBound,upperBound);

  TList *funListLambda = xiMass->GetListOfFunctions();
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
                        Form("#splitline{#splitline{#splitline{#Xi^{-}: %.0f}{m_{#Xi} = %.1fMeV/#it{c}^{2}}}{#sigma_{#Xi}= %.1fMeV/#it{c}^{2}}}{Purity = %.1f %%}",
                             signal, meanMass*1000.f, meanWidth*1000.f,
                             signal/(signal+background)*100.f));
  TLatex BeamTextLambda;
  BeamTextLambda.SetNDC(kTRUE);
  BeamTextLambda.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.2,
                           "p-Pb (2016) #sqrt{#it{s_{NN}}} = 5.02 TeV");
  //  TH1F *dummy=new TH1F("dummy","dummy",9,0.6,7.121051);//,10,xiMass2D->GetXaxis()->GetBinLowEdge(1),xiMass2D->GetXaxis()->GetBinUpEdge(21));
  //  dummy->GetXaxis()->SetTitle("P_{T} (GeV/#it{c})");
  //  dummy->GetYaxis()->SetTitle("Purity");
  //  dummy->GetYaxis()->SetRangeUser(80,100);
  //  dummy->SetTitle(";P;");
  TCanvas *cPurity=new TCanvas("cPurity","cPurity");
  cPurity->cd();
  SetStyleHisto(Purity,21,2);
  //  Purity->SetLineColor(fColors[2]);
  //  Purity->SetMarkerColor(fColors[2]);
  //  Purity->SetMarkerStyle(20);
  //  Purity->SetMarkerSize(2);
  Purity->GetYaxis()->SetRangeUser(80,100);
  Purity->SetTitle(";P_{T} (GeV/#it{c});Purity (%)");
  Purity->SetStats(0);
  Purity->GetYaxis()->SetTitleOffset(1.);
  Purity->GetXaxis()->SetTitleOffset(.9);
  //  dummy->DrawCopy();
//  TGaxis::SetMaxDigits(4);
  Purity->DrawCopy("p");
  std::cout << " now the very good lambda fit. D.J.Trump \n";
  LambdaMass->GetXaxis()->SetRangeUser(1.105,1.13);
  lowerBound=1.116-0.005;
  upperBound=1.116+0.005;
  LambdaMass->GetXaxis()->SetTitle("M_{#pip} (GeV/#it{c}^{2})");

  FitTheLambda(LambdaMass, signal, signalErr, background, backgroundErr, lowerBound,upperBound,Lambda);

  LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,
                        Form("#splitline{#splitline{#Lambda: %.0f}{Purity = %.1f %%}}{S/B=%.1f}",
                             signal, signal/(signal+background)*100.f,signal/background));
}
