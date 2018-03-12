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
  fSignalGauss->SetParameter(3, 0.2 * histo->GetMaximum());
  fSignalGauss->SetParameter(4, fSignalSingleGauss->GetParameter(1));
  fSignalGauss->SetParLimits(4,
                             fSignalSingleGauss->GetParameter(1)-fSignalSingleGauss->GetParameter(2),
                             fSignalSingleGauss->GetParameter(1)+fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParameter(5, 0.5*fSignalSingleGauss->GetParameter(2));
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


void PlotXi(TString fileName) {
  SetStyle(false,false);
  TH1F* xiPt[9];

  TH2F* xiMass2D;
  TH1F* xiMass;
  TFile *file=TFile::Open(fileName.Data());
  //  file->ls();
  TDirectoryFile *dirExp=(TDirectoryFile*)(file->FindObjectAny("CascadeCuts"));
  TCanvas *Xi= new TCanvas("Xi","Xi",0,0,1500,1100);
  TCanvas *PtXi= new TCanvas("PtXi","PtXi",0,0,1500,1100);
  PtXi->Divide(3,3);
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
    xiMass2D=(TH2F*)Cascade->FindObject("InvMassXiPt");
    xiMass=(TH1F*)xiMass2D->ProjectionY("InvMassXi");
    for (int iBin = 0; iBin<9; ++iBin) {
      PtXi->cd(1+iBin);
      TString XiNamePt=Form("XiPt_%i",iBin);
      std::cout << 2*iBin+1 << '\t' << 2*iBin+2 << std::endl;
      xiPt[iBin]=(TH1F*)xiMass2D->ProjectionY(XiNamePt.Data(),2*iBin+1,2*iBin+2);

      SetStyleHisto(xiPt[iBin],20,1);
      std::cout << Form("%2.2f < p_{T} < %2.2f",xiMass2D->GetXaxis()->GetBinLowEdge(2*iBin+1),xiMass2D->GetXaxis()->GetBinUpEdge(2*iBin+2)) << std::endl;
      xiPt[iBin]->GetXaxis()->SetRangeUser(1.285,1.345);
      xiPt[iBin]->GetXaxis()->SetTitle("IM_{#pi#Lambda} (GeV/#it{c}^{2})");
      float signal=0;
      float signalErr=0;
      float background=0;
      float backgroundErr=0;
      float lowerBound=1.322-0.005;
      float upperBound=1.322+0.005;
      if (iBin > 1) {
        //  xiMass->DrawCopy();
        //  SetStyle(false,false);
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
        LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
        LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,
                              Form("#splitline{#splitline{#splitline{#splitline{%.2f < p_{T} < %.2f}{#Xi^{-}: %.0f}}{m_{#Xi} = %.1fMeV/#it{c}^{2}}}{#sigma_{#Xi}= %.1fMeV/#it{c}^{2}}}{Purity = %.1f %%}",
                                   xiMass2D->GetXaxis()->GetBinLowEdge(2*iBin+1),xiMass2D->GetXaxis()->GetBinUpEdge(2*iBin+2),
                                   signal, meanMass*1000.f, meanWidth*1000.f,
                                   signal/(signal+background)*100.f));//#splitline{#splitline{
      } else {
        xiPt[iBin]->DrawCopy();
        TLatex LambdaLabel;
        LambdaLabel.SetNDC(kTRUE);
        LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
        LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,
                              Form("%.2f < p_{T} < %.2f}",xiMass2D->GetXaxis()->GetBinLowEdge(2*iBin+1),xiMass2D->GetXaxis()->GetBinUpEdge(2*iBin+2)));
      }
    }
    SetStyleHisto(xiMass,20,1);
    xiMass->GetXaxis()->SetRangeUser(1.285,1.345);
    xiMass->GetXaxis()->SetTitle("IM_{#pi#Lambda} (GeV/#it{c}^{2})");
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
                           "p-p (2017) #sqrt{#it{s}} = 13 TeV");
}
