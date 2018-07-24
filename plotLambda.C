#include "TH2F.h"
#include "TSystem.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7};
std::vector<int> fColors     = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2};
std::vector<int> fMarkers    = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

/* Data from https://drupal.star.bnl.gov/STAR/files/starpublications/221/data.html
Data for Figure 4
1/a0 (fm -1)                          1/|a0| (fm -1)                               reff
Lm-Lm(STAR)  0.91± 0.31(stat)±0.07(sys)±0.56(sys) 0.91± 0.31(stat)±0.07(sys)±0.56(sys) 8.52±2.56(stat)±0.74(sys)±2.09 (sys)
Lm-Lm(Nagara)                    -1.74                                        1.74                                              6.45
Lm-Lm(Nagara)                    -1.30                                        1.30                                              6.59
*/

/* Data from https://journals-aps-org.eaccess.ub.tum.de/prc/pdf/10.1103/PhysRevC.91.024916
Data for Table 1 and 2
model   a0 [fm] rEff [fm]
ND46    4.621  1.300
ND48    14.394  1.633
ND50    -10.629 2.042
ND52    −3.483  2.592
ND54    −1.893  3.389
ND56    −1.179  4.656
ND58    −0.764  6.863

NF42    3.659   0.975
NF44    23.956  1.258
NF46    −3.960  1.721
NF48    −1.511  2.549
NF50    −0.772  4.271
NF52    −0.406  8.828

NSC89-1020  −0.250  7.200
NSC89-920   −2.100  1.900
NSC89-820   −1.110  3.200

NSC97a  −0.329  12.370
NSC97b  −0.397  10.360
NSC97c  −0.476  9.130
NSC97d  −0.401  1.150
NSC97e  −0.501  9.840
NSC97f  −0.350  16.330

Ehime   −4.21   2.41

fss2    −0.81   3.99

ESC08   −0.97   3.86

HKMYY   −0.575  6.45

FG      −0.77   6.59
*/
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
  Int_t palette[4];
  palette[0] = kWhite;
  palette[1] = 33;
  palette[2] = 38;
  palette[3] = kGray;
  gStyle->SetPalette(4,palette);
//
//  const int NRGBs = 5;
//  Double_t stops[NRGBs];
//  for(int i=0; i<NRGBs; ++i) stops[i] = float(i)/(NRGBs-1);
//  Double_t red[NRGBs]   = { 1.,  29./255., 25./255., 27./255., 32./255.};
//  Double_t green[NRGBs] = { 1., 221./255., 160./255., 113./255., 74./255.};
//  Double_t blue[NRGBs] = {  1., 221./255., 184./255., 154./255., 129./255.};
//  TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
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
void plotLambda(const int flag = 0) {

  SetStyle();
  gStyle->SetCanvasPreferGL(1);

  TH1D* hist_emptyhist = new TH1D("hist_emptyhist","; 1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm)",50,-2,5);

  // Data
  TFile *file2 = TFile::Open("CombinedMap.root");
  TH2F *sigma, *ledniSucks;
  TGraphAsymmErrors *ledniSucks2;
  bool GetInvLedniSucks = true;
  TFile *file3;
  if (GetInvLedniSucks) file3 = TFile::Open("Friendship.root");
  // Combined fit
  if(flag == 0) {
    sigma = (TH2F*)file2->Get("hGlobNsigma");
    ledniSucks = (TH2F*)file2->Get("hGlobMinCk");
    if (GetInvLedniSucks) {
      ledniSucks2 = (TGraphAsymmErrors*)file3->Get("invLedniSucks");
    }
  }
  // pp 13 TeV
  if(flag == 1) {
    sigma = (TH2F*)file2->Get("hNsigma_1");
    ledniSucks = (TH2F*)file2->Get("hMinCk_1");
  }

  // p-Pb 5.02 TeV
  if(flag == 2) {
    sigma = (TH2F*)file2->Get("hNsigma_2");
    ledniSucks = (TH2F*)file2->Get("hMinCk_2");
  }

  SetStyleHisto(sigma);
  sigma->SetMaximum(5);
  sigma->SetMinimum(1);
  sigma->SetTitle(";;;#it{n_{#sigma}}");
  sigma->GetZaxis()->SetLabelFont(hist_emptyhist->GetXaxis()->GetLabelFont());
  sigma->GetZaxis()->SetTitleFont(hist_emptyhist->GetXaxis()->GetTitleFont());
  sigma->GetZaxis()->SetLabelSize(hist_emptyhist->GetXaxis()->GetLabelSize());
  sigma->GetZaxis()->SetTitleSize(hist_emptyhist->GetXaxis()->GetTitleSize());
  TH2F *sigma2 = (TH2F*)sigma->Clone();
  const int nContours = 4;
  double contours[nContours];
  contours[0] = 0;
  contours[1] = 3;
  contours[2] = 5;
  contours[3] = 999;
  sigma->SetContour(nContours,contours);
  sigma2->SetContour(nContours,contours);


  const int nContoursLedni = 1;
  double contoursLedni[nContoursLedni];
  contoursLedni[0] = 0.00000001;
  ledniSucks->SetContour(nContoursLedni, contoursLedni);
  ledniSucks->SetLineColorAlpha(kGray,0.5);
  ledniSucks->SetLineWidth(2);

  if (GetInvLedniSucks) {
    const int nContoursLedni2 = 4;
    double contoursLedni2[nContoursLedni2];
    contoursLedni2[0] = 0;
    contoursLedni2[1] = 3;
    contoursLedni2[2] = 5;
    contoursLedni2[3] = 14;
//    ledniSucks2->SetContour(nContoursLedni2, contoursLedni2);

  }

  // PRL C02 (2015) 022301.
  TGraphErrors *grStar = new TGraphErrors();
  grStar->SetPoint(0, -0.91, 8.52);
  grStar->SetPointError(0, 0.31, 2.56);
  SetStyleGraph(grStar, 8, 0);
  grStar->SetMarkerSize(2);
  TGraphAsymmErrors *grStar_syst = new TGraphAsymmErrors();
  grStar_syst->SetPoint(0, -0.91, 8.52);
  grStar_syst->SetPointError(0, 0.56, 0.07, 0.74, 2.09);
  grStar_syst->SetFillColorAlpha(kBlack, 0.4);

//  TGraph *grNagara = new TGraph();
//  grNagara->SetPoint(0, -1.74, 6.45);
//  grNagara->SetPoint(1, -1.3, 6.59);
//  SetStyleGraph(grNagara, 8, 1);
//  grNagara->SetMarkerSize(2);

  // ND
  // M. M. Nagels, T. A. Rijken, and J. J. de Swart, Phys. Rev. D 15, 2547 (1977).
  // Phys. Rev. D 15, 2547 (1977).
  float NDoneOvera0[] = {- 1./4.621, - 1./14.394, 1./10.629, 1./3.483, 1./1.893, 1./1.179, 1./0.764};
  float NDrEff[] = {1.300, 1.633, 2.042, 2.592, 3.389, 4.656, 6.863};
  TGraph *grND = new TGraph(7, NDoneOvera0, NDrEff);
  SetStyleGraph(grND, 2, 2);
  grND->SetLineStyle(3);

  // NF
  // M. M. Nagels, T. A. Rijken, and J. J. de Swart, Phys. Rev. D 20, 1633 (1979).
  // Phys. Rev. D 20, 1633 (1979).
  float NFoneOvera0[] = {- 1./3.659, - 1./23.956, 1./3.960, 1./1.511, 1./0.772, 1./0.406};
  float NFrEff[] = {0.975, 1.258, 1.721, 2.549, 4.271, 8.828};
  TGraph *grNF = new TGraph(6, NFoneOvera0, NFrEff);
  SetStyleGraph(grNF, 3, 3);
  grNF->SetLineStyle(3);

  // NSC89
  // P. M. M. Maessen, T. A. Rijken, and J. J. de Swart, Phys. Rev. C 40, 2226 (1989).
  // Phys. Rev. C 40, 2226 (1989).
  float NSC89oneOvera0[] = {1./0.25, 1./2.1, 1./1.1};
  float NSC89rEff[] = {7.2, 1.9, 3.2};
  TGraph *grNSC89 = new TGraph(3, NSC89oneOvera0, NSC89rEff);
  SetStyleGraph(grNSC89, 4, 4);
  grNSC89->SetLineStyle(2);
  grNSC89->SetMarkerSize(1.4);

  // NSC97
  // T. A. Rijken, V. G. J. Stoks, and Y. Yamamoto, Phys. Rev. C 59, 21 (1999).
  // Phys. Rev. C 59, 21 (1999).
  float NSC97oneOvera0[] = {1./0.329, 1./0.397, 1./0.476, 1./0.401, 1./0.501, 1/0.350};
  float NSC97rEff[] = {12.370, 10.360, 9.130, 1.150, 9.840, 16.330};
  TGraph *grNSC97 = new TGraph(6, NSC97oneOvera0, NSC97rEff);
  SetStyleGraph(grNSC97, 4, 5);
  grNSC97->SetLineStyle(2);
  grNSC97->SetMarkerSize(1.4);

  // Ehime
  // T. Ueda, K. Tominaga, M. Yamaguchi, N. Kijima, D. Okamoto, K. Miyagawa, and T. Yamada, Prog. Theor. Phys. 99, 891 (1998);
  // K. Tominaga, T. Ueda, M. Yamaguchi, N. Kijima, D. Okamoto, K. Miyagawa, and T. Yamada, Nucl. Phys. A 642, 483 (1998).
  // Prog. Theor. Phys. 99, 891 (1998); Nucl. Phys. A 642, 483 (1998).
  TGraph *grEhime = new TGraph();
  grEhime->SetPoint(0, 1./4.21, 2.51);
  SetStyleGraph(grEhime, 0, 2);

  // fss2
  // Y. Fujiwara, Y. Suzuki, and C. Nakamoto, Prog. Part. Nucl. Phys. 58, 439 (2007);
  // Y. Fujiwara, M. Kohno, C. Nakamoto, and Y. Suzuki, Phys. Rev. C 64, 054001 (2001).
  // Part. Nucl. Phys. 58, 439 (2007); Phys. Rev. C 64, 054001 (2001).
  TGraph *grfss2 = new TGraph();
  grfss2->SetPoint(0, 1./0.81, 3.99);
  SetStyleGraph(grfss2, 1, 3);

  // ESC8
  // T. A. Rijken, M. M. Nagels, and Y. Yamamoto, Prog. Theor. Phys. Suppl. 185, 14 (2010).
  // Prog. Theor. Phys. Suppl. 185, 14 (2010).
  TGraph *grESC8 = new TGraph();
  grESC8->SetPoint(0, 1./0.97, 3.86);
  SetStyleGraph(grESC8, 7, 4);
  grESC8->SetMarkerSize(1.4);

  // HKMYY - NAGARA
  // E. Hiyama, M. Kamimura, T. Motoba, T. Yamada, and Y. Yamamoto, Phys.Rev. C 66, 024007 (2002).
  // Phys.Rev. C 66, 024007 (2002).
  TGraph *grHKMYY = new TGraph();
  grHKMYY->SetPoint(0, 1./0.575, 6.45);
  SetStyleGraph(grHKMYY, 8, 5);
  grHKMYY->SetMarkerSize(2);

  // FG - NAGARA
  // I.N. Filikhin and A. Gal, Nucl. Phys. A 707, 491 (2002).
  // Nucl. Phys. A 707, 491 (2002).
  TGraph *grFG = new TGraph();
  grFG->SetPoint(0, 1./0.77, 6.59);
  SetStyleGraph(grFG, 8, 6);
  grFG->SetMarkerSize(2);

  const float splitFraction = 0.625;
  TCanvas *c = new TCanvas("c", "LambdaLambda", 700./splitFraction, 500);
  TPad *pad1 = new TPad("lambda", "", 0.0,0.,splitFraction,1);
  TPad *pad2 = new TPad("lambda", "", splitFraction,0.,1,1);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  pad1->SetTopMargin(0.025);
  pad1->SetRightMargin(0.01);
  hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
  hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);
  hist_emptyhist->Draw();
  sigma->Draw("cont3 same");
  sigma->SetLineColor(kGray+3);
  sigma->SetLineWidth(2);
  grStar_syst->Draw("2same");
  grStar->Draw("PZ same");
//  grNagara->Draw("P same");
  grND->Draw("PL same");
  grNF->Draw("PL same");
  grNSC89->Draw("PL same");
  grNSC97->Draw("PL same");
  grEhime->Draw("PL same");
  grfss2->Draw("PL same");
  grESC8->Draw("PL same");
  grHKMYY->Draw("PL same");
  grFG->Draw("PL same");

  TGraphErrors *grStarFake = new TGraphErrors();
  grStarFake->SetMarkerStyle(grStar->GetMarkerStyle());
  grStarFake->SetMarkerColor(grStar->GetMarkerColor());
  grStarFake->SetMarkerSize(grStar->GetMarkerSize());
  grStarFake->SetFillColor(fFillColors[0]);
  grStarFake->SetLineColor(fFillColors[0]);

  TGraphErrors *grEntry1Fake = new TGraphErrors();
  grEntry1Fake->SetMarkerStyle(kOpenSquare);
  grEntry1Fake->SetMarkerColor(33);
  grEntry1Fake->SetMarkerSize(1.1);
  grEntry1Fake->SetFillColor(33);
  grEntry1Fake->SetLineColor(kBlack);

  TGraphErrors *grEntry2Fake = new TGraphErrors();
  grEntry2Fake->SetMarkerStyle(grStar->GetMarkerStyle());
  grEntry2Fake->SetMarkerColor(38);
  grEntry2Fake->SetMarkerSize(1.1);
  grEntry2Fake->SetFillColor(38);
  grEntry2Fake->SetLineColor(kBlack);

  TGraphErrors *grEntry3Fake = new TGraphErrors();
  grEntry3Fake->SetMarkerStyle(grStar->GetMarkerStyle());
  grEntry3Fake->SetMarkerColor(kWhite);
  grEntry3Fake->SetMarkerSize(grStar->GetMarkerSize());
  grEntry3Fake->SetFillColor(kWhite);
  grEntry3Fake->SetLineColor(kBlack);

  pad2->cd();
  pad2->SetTopMargin(0.025);
  TLegend *leg = new TLegend(0.01, 0.18, 0.4, 0.95);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize()*0.85);
  leg->AddEntry(grStarFake, "STAR", "PF");
  leg->AddEntry(grHKMYY, "HKMYY", "P");
  leg->AddEntry(grFG, "FG", "P");
  leg->AddEntry(grND, "ND", "LP");
  leg->AddEntry(grNF, "NF", "LP");
  leg->AddEntry(grNSC89, "NSC89", "LP");
  leg->AddEntry(grNSC97, "NSC97", "LP");
  leg->AddEntry(grEhime, "Ehime", "P");
  leg->AddEntry(grESC8, "ESC08", "P");
  leg->AddEntry(grfss2, "fss2", "P");
  leg->Draw("same");

  TLatex ref;
  ref.SetTextSize(gStyle->GetTextSize()*0.68);
  ref.SetNDC(kTRUE);

  const float xAll = 0.35;
  ref.DrawLatex(xAll, 0.895, "PRL C02 (2015) 022301."); // STAR
  ref.DrawLatex(xAll, 0.81828, "Phys.Rev. C 66, 024007 (2002)."); // HKMYY
  ref.DrawLatex(xAll, 0.74162, "Nucl. Phys. A 707, 491 (2002).");  //  FG
  ref.DrawLatex(xAll, 0.66496, "Phys. Rev. D 15, 2547 (1977).");  //  ND
  ref.DrawLatex(xAll, 0.58833, "Phys. Rev. D 20, 1633 (1979).");  //  NF
  ref.DrawLatex(xAll, 0.51164, "Phys. Rev. C 40, 2226 (1989).");  //  NSC89
  ref.DrawLatex(xAll, 0.43599, "Phys. Rev. C 59, 21 (1999).");  //  NSC97
  ref.DrawLatex(xAll, 0.35833, "Nucl. Phys. A 642, 483 (1998).");  //  Ehime Prog. Theor. Phys. 99, 891 (1998);
  ref.DrawLatex(xAll, 0.28166, "Phys. Rev. C 64, 054001 (2001).");  //  fss2 Part. Nucl. Phys. 58, 439 (2007);
  ref.DrawLatex(xAll, 0.205, "Prog. Theor. Phys. Suppl. 185, 14 (2010).");  //  ESC08

  grStar_syst->SetFillColorAlpha(kGray+1, 0.5);
  const float splitFraction2 = 0.8;
  TCanvas *c2 = new TCanvas("c2", "LambdaLambda", 700./splitFraction2, 500);
  TPad *pad12 = new TPad("lambda", "", 0.0,0.,splitFraction2,1);
  TPad *pad22 = new TPad("lambda", "", splitFraction2,0.,1,1);
  pad12->Draw();
  pad22->Draw();
  pad12->cd();
  pad12->SetTopMargin(0.025);
  pad12->SetRightMargin(0.01);
  hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
  hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);
  hist_emptyhist->Draw();
  sigma2->DrawCopy("col same");
  sigma->Draw("cont3 same");
  sigma->SetLineColor(kWhite);
  sigma->SetLineWidth(1);
  sigma->SetLineStyle(2);

  if (GetInvLedniSucks) {
    ledniSucks2->SetFillStyle(3344);
    ledniSucks2->SetLineStyle(3344);
    ledniSucks2->SetFillColorAlpha(kGray,0.7 );
    ledniSucks2->SetLineColorAlpha(kGray,0.7 );


    ledniSucks2->Draw("F SAME");
  } else {
    ledniSucks->Draw("cont3 same");
  }


  grStar_syst->Draw("2same");
  grStar->Draw("PZ same");
//  grNagara->Draw("P same");
  grND->Draw("PL same");
  grNF->Draw("PL same");
  grNSC89->Draw("PL same");
  grNSC97->Draw("PL same");
  grEhime->Draw("PL same");
  grfss2->Draw("PL same");
  grESC8->Draw("PL same");
  grHKMYY->Draw("PL same");
  grFG->Draw("PL same");

  pad12->RedrawAxis();

  TLatex sigmaLabel;
  sigmaLabel.SetTextSize(gStyle->GetTextSize()*0.9);
  sigmaLabel.SetTextFont(62);
  sigmaLabel.SetTextColor(kBlack);
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText.SetNDC(kTRUE);

  if(flag == 0) {
//    sigmaLabel.SetTextAngle(32);
//    sigmaLabel.DrawLatex(1.35, 7.4, "1#kern[0.3]{#sigma}");
    sigmaLabel.SetTextAngle(-50);
    sigmaLabel.DrawLatex(0.7, 1.5, "3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-55);
    sigmaLabel.DrawLatex(0.4, 1.5, "5#kern[0.2]{#sigma}");

    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.5, 1.2, "3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.35, 2.4, "5#kern[0.2]{#sigma}");

    BeamText.SetTextSize(gStyle->GetTextSize()*0.75);
    BeamText.DrawLatex(0.75, 0.365, "ALICE Preliminary");
    BeamText.DrawLatex(0.75, 0.32, "pp #sqrt{#it{s}} = 7 TeV");
    BeamText.DrawLatex(0.75, 0.275, "pp #sqrt{#it{s}} = 13 TeV");
    BeamText.DrawLatex(0.75, 0.23, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    BeamText.DrawLatex(0.75, 0.185, "#Lambda#Lambda #oplus #bar{#Lambda}#bar{#Lambda} pairs");
  }
  if(flag == 1) {
    sigmaLabel.SetTextAngle(12.5);
    sigmaLabel.DrawLatex(1.35, 7.5, "1#kern[0.3]{#sigma}");
    sigmaLabel.SetTextAngle(-50);
    sigmaLabel.DrawLatex(0.7, 1.5, "3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-65);
    sigmaLabel.DrawLatex(0.35, 1.5, "5#kern[0.2]{#sigma}");

    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.2, 1.5, "3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.05, 2.7, "5#kern[0.2]{#sigma}");

    BeamText.DrawLatex(0.6, 0.305, "ALICE Preliminary");
    BeamText.DrawLatex(0.6, 0.255, "pp #sqrt{#it{s}} = 13 TeV");
    BeamText.DrawLatex(0.6, 0.205, "#Lambda#Lambda #oplus #bar{#Lambda}#bar{#Lambda} pairs");
  }
  if(flag == 2) {
    sigmaLabel.SetTextAngle(320);
    sigmaLabel.DrawLatex(-1.74,1.16,"3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.74,3.41,"5#kern[0.2]{#sigma}");

    sigmaLabel.SetTextAngle(300);
    sigmaLabel.DrawLatex(0.41,1.25,"3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-65);
    sigmaLabel.DrawLatex(0.15,1.25, "5#kern[0.2]{#sigma}");

    BeamText.DrawLatex(0.6, 0.305, "ALICE Preliminary");
    BeamText.SetTextSize(gStyle->GetTextSize()*0.65);
    BeamText.DrawLatex(0.6, 0.255, "p-Pb #sqrt{#it{s_{NN}}} = 5.02 TeV");
    BeamText.DrawLatex(0.6, 0.205, "#Lambda#Lambda #oplus #bar{#Lambda}#bar{#Lambda} pairs");
  }

  pad22->cd();
  pad22->SetTopMargin(0.025);
  TLegend *leg2 = new TLegend(0.01, 0.01, 0.99, 0.89);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(gStyle->GetTextSize()*2.25);
//  leg2->AddEntry(grEntry3Fake, "Exclusion", "");
//  leg2->AddEntry(grEntry3Fake, "#sigma < 3", "PF");
  leg2->AddEntry(grEntry1Fake, "3 < #sigma < 5", "PF");
  leg2->AddEntry(grEntry2Fake, "#sigma >5 ", "PF");
  leg2->AddEntry(grEntry3Fake, " ", "");
  leg2->AddEntry(grStarFake, "STAR", "PF");
  leg2->AddEntry(grHKMYY, "HKMYY", "P");
  leg2->AddEntry(grFG, "FG", "P");
  leg2->AddEntry(grND, "ND", "LP");
  leg2->AddEntry(grNF, "NF", "LP");
  leg2->AddEntry(grNSC89, "NSC89", "LP");
  leg2->AddEntry(grNSC97, "NSC97", "LP");
  leg2->AddEntry(grEhime, "Ehime", "P");
  leg2->AddEntry(grESC8, "ESC08", "P");
  leg2->AddEntry(grfss2, "fss2", "P");
  leg2->Draw("same");
  if(flag == 0) c2->SaveAs("LambdaLambda_global.pdf");
  if(flag == 1) c2->SaveAs("LambdaLambda_13TeV.pdf");
  if(flag == 2) c2->SaveAs("LambdaLambda_pPb502TeV.pdf");


}
