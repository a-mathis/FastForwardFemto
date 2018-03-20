#include <iostream>
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TFile.h"
#include "TFractionFitter.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TText.h"

std::vector<int> fillColors = {kGray + 1,  kRed - 10,    kBlue - 9,
                               kGreen - 8, kMagenta - 9, kOrange - 9,
                               kCyan - 8,  kYellow - 7};
std::vector<int> colors = {kBlack,       kRed + 1,    kBlue + 2, kGreen + 3,
                           kMagenta + 1, kOrange - 1, kCyan + 2, kYellow + 2};
std::vector<int> markers = {
    kFullCircle, kFullSquare, kOpenCircle,  kOpenSquare, kOpenDiamond,
    kOpenCross,  kFullCross,  kFullDiamond, kFullStar,   kOpenStar};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette = false, bool title = true) {
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
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
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
  histo->SetMarkerStyle(markers[marker]);
  histo->SetMarkerColor(colors[color]);
  histo->SetLineColor(colors[color]);
  histo->SetLineWidth(2);
  histo->Sumw2();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawHist(TH1F* hist1, TString HistTitle, TString XTitle, TString YTitle,
              Int_t option, Int_t MarkerColor, Int_t MarkerSize,
              Bool_t ShowTitle) {
  hist1->SetStats(0);
  if (!ShowTitle) {
    hist1->SetTitle(0);
  } else {
    hist1->SetTitle(HistTitle.Data());
  }
  hist1->GetXaxis()->SetTitle(XTitle.Data());
  hist1->GetYaxis()->SetTitle(YTitle.Data());
  hist1->SetMarkerColor(colors[0]);
  hist1->SetLineColor(colors[0]);
  hist1->SetMarkerStyle(markers[0]);
  if (option == 0)
    hist1->DrawCopy("PE");
  else if (option == 1)
    hist1->DrawCopy("histo");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawHist(TH2F* hist1, TString HistTitle, TString XTitle, TString YTitle,
              Int_t option) {
  hist1->SetStats(0);
  // hist1->SetTitle(0);
  hist1->SetTitle(HistTitle.Data());
  hist1->GetXaxis()->SetTitle(XTitle.Data());
  hist1->GetYaxis()->SetTitle(YTitle.Data());

  if (option == 0)
    hist1->DrawCopy("colz");
  else if (option == 1)
    hist1->DrawCopy("Text0");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawHist(TH1F* hist1, TH1F* hist2, TString HistTitle, TString XTitle,
              TString YTitle, Int_t option, Int_t MarkerColor, Int_t MarkerSize,
              Int_t MarkerColor2, Int_t MarkerSize2, Bool_t drawTitle) {
  hist1->SetStats(0);
  if (drawTitle) {
    hist1->SetTitle(HistTitle.Data());
  } else {
    hist1->SetTitle(0);
  }
  hist1->GetXaxis()->SetTitle(XTitle.Data());
  hist1->GetYaxis()->SetTitle(YTitle.Data());

  hist1->SetMarkerSize(MarkerSize);
  hist1->SetMarkerColor(colors[MarkerColor - 1]);
  hist1->SetLineColor(colors[MarkerColor - 1]);
  hist1->SetMarkerStyle(markers[0]);

  hist2->SetMarkerSize(MarkerSize2);
  hist2->SetMarkerColor(colors[MarkerColor2 - 1]);
  hist2->SetLineColor(colors[MarkerColor2 - 1]);
  hist2->SetMarkerStyle(markers[1]);

  if (option == 0) {
    hist1->DrawCopy("PE");
    hist2->DrawCopy("PE same");
  } else if (option == 1) {
    hist1->DrawCopy("histo");
    // hist2->SetLineColor(2);
    hist2->DrawCopy("histo same");
  }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawHist(TH1F* hist1, TH1F* hist2, TH1F* hist3, TString HistTitle,
              TString XTitle, TString YTitle, Int_t option, Int_t MarkerColor,
              Int_t MarkerSize, Int_t MarkerColor2, Int_t MarkerSize2,
              Int_t MarkerColor3, Int_t MarkerSize3, Bool_t DrawTitle) {
  hist1->SetStats(0);
  if (!DrawTitle) {
    hist1->SetTitle(0);
  } else {
    hist1->SetTitle(HistTitle.Data());
  }
  hist1->GetXaxis()->SetTitle(XTitle.Data());
  hist1->GetYaxis()->SetTitle(YTitle.Data());
  hist1->SetMarkerSize(MarkerSize);
  hist1->SetMarkerColor(colors[MarkerColor - 1]);
  hist1->SetLineColor(colors[MarkerColor - 1]);
  hist1->SetMarkerStyle(markers[0]);

  hist2->SetMarkerSize(MarkerSize2);
  hist2->SetMarkerColor(colors[MarkerColor2 - 1]);
  hist2->SetLineColor(colors[MarkerColor2 - 1]);
  hist2->SetMarkerStyle(markers[1]);

  hist3->SetMarkerSize(MarkerSize3);
  hist3->SetMarkerColor(colors[MarkerColor3 - 1]);
  hist3->SetLineColor(colors[MarkerColor3 - 1]);
  hist3->SetMarkerStyle(markers[2]);

  if (option == 0) {
    hist1->DrawCopy("PE");
    hist2->DrawCopy("PE same");
    hist3->DrawCopy("PE same");
  } else if (option == 1) {
    hist1->DrawCopy("histo");
    hist2->DrawCopy("histo same");
    hist3->DrawCopy("histo same");
  }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawHist(TH2F* hist1, TString XTitle, TString YTitle, Int_t option,
              Float_t rangeleftx, Float_t rangerightx, Float_t rangelefty,
              Float_t rangerighty, Bool_t ShowTitle) {
  hist1->SetStats(0);
  if (!ShowTitle) hist1->SetTitle(0);
  hist1->GetXaxis()->SetTitle(XTitle.Data());
  hist1->GetYaxis()->SetTitle(YTitle.Data());
  hist1->GetXaxis()->SetRangeUser(rangeleftx, rangerightx);
  hist1->GetYaxis()->SetRangeUser(rangelefty, rangerighty);

  if (option == 0) hist1->DrawCopy("colz");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawHist(TH1F* hist1, TH1F* hist2, TH1F* hist3, TH1F* hist4,
              TString HistTitle, TString XTitle, TString YTitle, Int_t option,
              Int_t MarkerColor, Int_t MarkerSize, Int_t MarkerColor2,
              Int_t MarkerSize2, Int_t MarkerColor3, Int_t MarkerSize3,
              Int_t MarkerColor4, Int_t MarkerSize4, Bool_t showTitle, bool singlepanel=false) {
  hist1->SetStats(0);
  if (!showTitle) {
    hist1->SetTitle(0);
  } else {
    hist1->SetTitle(HistTitle.Data());
  }
  hist1->GetXaxis()->SetTitle(XTitle.Data());
  hist1->GetYaxis()->SetTitle(YTitle.Data());

  if(singlepanel) {
    hist1->GetXaxis()->SetLabelSize(0.045);
    hist1->GetXaxis()->SetTitleSize(0.05);
    hist1->GetXaxis()->SetLabelOffset(0.01);
    hist1->GetXaxis()->SetTitleOffset(1.2);
    hist1->GetXaxis()->SetLabelFont(42);
    hist1->GetYaxis()->SetLabelSize(0.045);
    hist1->GetYaxis()->SetTitleSize(0.05);
    hist1->GetYaxis()->SetLabelOffset(0.01);
    hist1->GetYaxis()->SetTitleOffset(1.25);
  }
  else {
    hist1->GetXaxis()->SetLabelSize(0.065);
    hist1->GetXaxis()->SetTitleSize(0.075);
    hist1->GetXaxis()->SetLabelFont(42);
    hist1->GetYaxis()->SetLabelSize(0.065);
    hist1->GetYaxis()->SetTitleSize(0.075);
    hist1->GetYaxis()->SetLabelFont(42);
    hist1->GetYaxis()->SetTitleOffset(1.25);
  }
  hist1->SetLineColor(colors[MarkerColor - 1]);
  hist2->SetLineColor(colors[MarkerColor2 - 1]);
  hist3->SetLineColor(colors[MarkerColor3 - 1]);
  hist4->SetLineColor(colors[MarkerColor4 - 1]);

  if (option == 0) {
    hist1->SetMarkerSize(MarkerSize);
    hist1->SetMarkerColor(colors[MarkerColor - 1]);
    hist1->SetMarkerStyle(markers[0]);

    hist2->SetMarkerSize(MarkerSize2);
    hist2->SetMarkerColor(colors[MarkerColor2 - 1]);
    hist2->SetMarkerStyle(markers[1]);

    hist3->SetMarkerSize(MarkerSize3);
    hist3->SetMarkerColor(colors[MarkerColor3 - 1]);
    hist3->SetMarkerStyle(markers[2]);

    hist4->SetMarkerSize(MarkerSize4);
    hist4->SetMarkerColor(colors[MarkerColor4 - 1]);
    hist4->SetMarkerStyle(markers[4]);

    hist1->DrawCopy("PE");
    hist2->DrawCopy("PE same");
    hist3->DrawCopy("PE same");
    hist4->DrawCopy("PE same");
  } else if (option == 1) {
    hist1->DrawCopy("histo");
    hist2->DrawCopy("histo same");
    hist3->DrawCopy("histo same");
    hist4->DrawCopy("histo same");
  } else if (option == 2) {
    hist1->SetMarkerSize(MarkerSize);
    hist1->SetMarkerColor(colors[MarkerColor - 1]);
    hist1->SetLineColor(colors[MarkerColor - 1]);
    hist1->SetMarkerStyle(markers[0]);

    hist2->SetFillStyle(3004);
    hist2->SetFillColor(colors[MarkerColor2 - 1]);
    hist2->SetLineColor(colors[MarkerColor2 - 1]);

    hist3->SetFillStyle(3005);
    hist3->SetFillColor(colors[MarkerColor3 - 1]);
    hist3->SetLineColor(colors[MarkerColor3 - 1]);

    hist4->SetFillStyle(3013);
    hist4->SetFillColor(colors[MarkerColor4 - 1]);
    hist4->SetLineColor(colors[MarkerColor4 - 1]);

    hist1->DrawCopy("PE");
    hist2->DrawCopy("histo same");
    hist3->DrawCopy("histo same");
    hist4->DrawCopy("histo same");
  }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawHist(TH1F* hist1, TH1F* hist2, TH1F* hist3, TH1F* hist4, TH1F* hist5,
              TString Title, TString XTitle, TString YTitle, Int_t option,
              Int_t MarkerColor, Int_t MarkerSize, Int_t MarkerColor2,
              Int_t MarkerSize2, Int_t MarkerColor3, Int_t MarkerSize3,
              Int_t MarkerColor4, Int_t MarkerSize4, Int_t MarkerColor5,
              Int_t MarkerSize5, Bool_t showTitle, bool singlepanel=false) {
  hist1->SetStats(0);
  if (!showTitle) {
    hist1->SetTitle(0);
  } else {
    hist1->SetTitle(Title.Data());
  }
  hist1->GetXaxis()->SetTitle(XTitle.Data());
  hist1->GetYaxis()->SetTitle(YTitle.Data());

  hist1->SetLineColor(colors[MarkerColor - 1]);
  hist2->SetLineColor(colors[MarkerColor2 - 1]);
  hist3->SetLineColor(colors[MarkerColor3 - 1]);
  hist4->SetLineColor(colors[MarkerColor4 - 1]);
  hist5->SetLineColor(colors[MarkerColor5 - 1]);

  if(singlepanel) {
    hist1->GetXaxis()->SetLabelSize(0.045);
    hist1->GetXaxis()->SetTitleSize(0.05);
    hist1->GetXaxis()->SetLabelOffset(0.01);
    hist1->GetXaxis()->SetTitleOffset(1.2);
    hist1->GetXaxis()->SetLabelFont(42);
    hist1->GetYaxis()->SetLabelSize(0.045);
    hist1->GetYaxis()->SetTitleSize(0.05);
    hist1->GetYaxis()->SetLabelOffset(0.01);
    hist1->GetYaxis()->SetTitleOffset(1.25);
  }
  else {
    hist1->GetXaxis()->SetNdivisions(505);
    hist1->GetXaxis()->SetLabelSize(0.065);
    hist1->GetXaxis()->SetTitleSize(0.075);
    hist1->GetXaxis()->SetLabelFont(42);
    hist1->GetYaxis()->SetLabelSize(0.065);
    hist1->GetYaxis()->SetTitleSize(0.075);
    hist1->GetYaxis()->SetLabelFont(42);
    hist1->GetYaxis()->SetTitleOffset(1.25);
  }

  if (option == 0) {
    hist1->SetMarkerSize(MarkerSize);
    hist1->SetMarkerColor(colors[MarkerColor - 1]);
    hist1->SetMarkerStyle(markers[0]);

    hist2->SetMarkerSize(MarkerSize2);
    hist2->SetMarkerColor(colors[MarkerColor2 - 1]);
    hist2->SetMarkerStyle(markers[1]);

    hist3->SetMarkerSize(MarkerSize3);
    hist3->SetMarkerColor(colors[MarkerColor3 - 1]);
    hist3->SetMarkerStyle(markers[2]);

    hist4->SetMarkerSize(MarkerSize4);
    hist4->SetMarkerColor(colors[MarkerColor4 - 1]);
    hist4->SetMarkerStyle(markers[3]);

    hist5->SetMarkerSize(MarkerSize5);
    hist5->SetMarkerColor(colors[MarkerColor5]);
    hist5->SetMarkerStyle(markers[4]);

    hist1->DrawCopy("PE");
    hist2->DrawCopy("PE same");
    hist3->DrawCopy("PE same");
    hist4->DrawCopy("PE same");
    hist5->DrawCopy("PE same");
  } else if (option == 1) {
    hist1->DrawCopy("histo");
    hist2->DrawCopy("histo same");
    hist3->DrawCopy("histo same");
    hist4->DrawCopy("histo same");
    hist5->DrawCopy("histo same");
  } else if (option == 2) {
    hist1->SetMarkerSize(MarkerSize);
    hist1->SetMarkerColor(colors[MarkerColor - 1]);
    hist1->SetMarkerStyle(markers[0]);

    hist2->SetFillStyle(3004);
    hist2->SetFillColor(colors[MarkerColor2 - 1]);
    hist2->SetLineColor(colors[MarkerColor2 - 1]);

    hist3->SetFillStyle(3005);
    hist3->SetFillColor(colors[MarkerColor3 - 1]);
    hist3->SetLineColor(colors[MarkerColor3 - 1]);

    hist4->SetFillStyle(3013);
    hist4->SetFillColor(colors[MarkerColor4 - 1]);
    hist4->SetLineColor(colors[MarkerColor4 - 1]);

    hist5->SetFillStyle(3006);
    hist5->SetFillColor(colors[MarkerColor5]);
    hist5->SetLineColor(colors[MarkerColor5]);

    hist1->DrawCopy("PE");
    hist2->DrawCopy("histo same");
    hist3->DrawCopy("histo same");
    hist4->DrawCopy("histo same");
    hist5->DrawCopy("histo same");
  }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawProf(TProfile* prof, TString XTitle, TString YTitle, Int_t option,
              Int_t MarkerColor, Int_t MarkerSize) {
  prof->SetStats(0);
  prof->SetTitle(0);
  prof->GetXaxis()->SetTitle(XTitle.Data());
  prof->GetYaxis()->SetTitle(YTitle.Data());
  prof->SetMarkerSize(MarkerSize);
  prof->SetMarkerColor(colors[MarkerColor - 1]);
  prof->SetLineColor(colors[MarkerColor - 1]);
  prof->SetMarkerStyle(20);
  if (option == 0)
    prof->DrawCopy("PE");
  else if (option == 1)
    prof->DrawCopy("histo");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void printval(TString val_name, TString val_dim, double val, double x,
              double y) {
  TString val_print = val_name;
  val_print += " = ";

  char buffer[200];

  sprintf(buffer, "%2.2f", val);
  val_print += buffer;
  val_print += " ";
  val_print += val_dim;

  TLatex* text = new TLatex(x, y, val_print);  // = new TLatex(x,y,val_print);
  text->SetTextFont(42);
  text->SetTextSize(0.06);
  text->SetNDC();
  text->SetTextColor(1);
  text->Draw("same");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void printval(TString val_name, double val, double valerr, double x, double y) {
  // val[0]: mean
  // val[1]: error

  TString val_print = val_name;
  val_print += " = ";

  char buffer[200];

  sprintf(buffer, "%2.5f", val);
  val_print += buffer;
  val_print += " #pm ";
  sprintf(buffer, "%2.5f", valerr);
  val_print += buffer;
  val_print += " ";

  TLatex* text = new TLatex(x, y, val_print);  // = new TLatex(x,y,val_print);
  text->SetTextFont(42);
  text->SetTextSize(0.06);
  text->SetNDC();
  text->SetTextColor(1);
  text->Draw("same");
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t GaussPolyFit(double* x, double* par) {
  double xval = x[0];
  double par0 = par[0];
  double par1 = par[1];
  double par2 = par[2];
  double par3 = par[3];
  double par4 = par[4];
  double par5 = par[5];
  double par6 = par[6];
  double par7 = par[7];
  double par8 = par[8];
  double par9 = par[9];
  double par10 = par[10];
  double par11 = par[11];

  double sig = par0 * (par1 * TMath::Gaus(xval, par2, par3) +
                       (1 - par1) * TMath::Gaus(xval, par4, par5));
  double bkg = par6 + par7 * xval + par8 * xval * xval +
               par9 * xval * xval * xval + par10 * pow(xval, 4.) +
               par11 * pow(xval, 5.);

  return sig + bkg;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t GaussFit(double* x, double* par) {
  double xval = x[0];
  double par0 = par[0];
  double par1 = par[1];
  double par2 = par[2];
  double par3 = par[3];
  double par4 = par[4];
  double par5 = par[5];

  double sig = par0 * (par1 * TMath::Gaus(xval, par2, par3) +
                       (1 - par1) * TMath::Gaus(xval, par4, par5));

  return sig;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t Gaussfit(Double_t* x, Double_t* par) {
  Double_t xval = x[0];
  Double_t par0 = par[0];
  Double_t par1 = par[1];
  Double_t par2 = par[2];

  Double_t y = par0 * TMath::Gaus(xval, par1, par2);

  return y;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t TwoGaussfit(Double_t* x, Double_t* par) {
  Double_t xval = x[0];
  Double_t par0 = par[0];
  Double_t par1 = par[1];
  Double_t par2 = par[2];
  Double_t par3 = par[3];
  Double_t par4 = par[4];
  Double_t par5 = par[5];

  Double_t y = par0 * TMath::Gaus(xval, par1, par2) +
               par3 * TMath::Gaus(xval, par4, par5);

  return y;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t CrystalBallbkg(double* x, double* par) {
  double xcur = x[0];
  xcur *= -1;
  double alpha = par[0];
  double n = par[1];
  double mu = par[2];
  double sigma = par[3];
  double N = par[4];

  double a = par[5];
  double b = par[6];
  double c = par[7];
  double d = par[8];
  double e = par[9];

  double poly = a + b * xcur + c * xcur * xcur + d * xcur * xcur * xcur +
                e * xcur * xcur * xcur * xcur;
  double A = 0.;
  double B = 0.;
  if (alpha < 0.) {
    A = pow((n / (-1. * alpha)), n) * TMath::Exp((-1.) * alpha * alpha / 2.);
    B = n / (-1. * alpha) + alpha;
  } else {
    A = pow((n / alpha), n) * TMath::Exp((-1.) * alpha * alpha / 2.);
    B = n / alpha - alpha;
  }

  double f;
  if ((xcur - mu) / sigma > (-1.) * alpha)
    f = N *
        TMath::Exp((-1.) * (xcur - mu) * (xcur - mu) / (2. * sigma * sigma));
  else
    f = N * A * pow((B - (xcur - mu) / sigma), (-1 * n));

  return f + poly;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1D* hist_projection(TH2F* hist_2D, Double_t projection_low,
                      Double_t projection_up, TString Projection,
                      TString HistName) {
  TH1D* hist_projected = 0;

  if (Projection == "X") {
    Int_t bin_low = hist_2D->GetYaxis()->FindBin(projection_low);
    Int_t bin_up = hist_2D->GetYaxis()->FindBin(projection_up);
    hist_projected =
        hist_2D->ProjectionX(HistName.Data(), bin_low, bin_up,
                             "e");  //"e" for recalculation of the error
  } else if (Projection == "Y") {
    Int_t bin_low = hist_2D->GetXaxis()->FindBin(projection_low);
    Int_t bin_up = hist_2D->GetXaxis()->FindBin(projection_up);
    hist_projected =
        hist_2D->ProjectionY(HistName.Data(), bin_low, bin_up, "e");
  } else {
    std::cout << "specifiy your axes projection" << std::endl;
  }

  return hist_projected;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void print_ranges(Float_t start, Float_t range, int bin, TString observable,
                  TString unit, Double_t x, Double_t y) {
  char buffer[1000];
  TString total = "";

  Float_t bound_low = start + range * bin;
  Float_t bound_up = start + range * (bin + 1);
  sprintf(buffer, "%.2f", bound_low);
  total += buffer;
  total += " < ";
  total += observable;
  total += " < ";
  sprintf(buffer, "%.2f", bound_up);
  total += buffer;
  total += " " + unit;
  TLatex* t = new TLatex();
  //-1.,80.5,total.Data());
  // t->SetTextAlign(22);
  // t->SetTextColor(kRed+2);
  t->SetTextFont(43);
  t->SetTextSize(20);
  // t->SetTextAngle(45);
  t->DrawLatex(x, y, total.Data());
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TString print_ranges(Float_t start, Float_t range, int bin, TString observable,
                     TString unit) {
  char buffer[1000];
  TString total = "";

  Float_t bound_low = start + range * bin;
  Float_t bound_up = start + range * (bin + 1);
  sprintf(buffer, "%.2f", bound_low);
  total += buffer;
  total += " < ";
  total += observable;
  total += " < ";
  sprintf(buffer, "%.2f", bound_up);
  total += buffer;
  total += " " + unit;

  return total;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TString Create_path(TString basedir, TString subdir, TString filename,
                    TString MCorData) {
  TString printpath = basedir;
  printpath += "/" + subdir + "/";
  TString printname = filename + "_" + MCorData;
  printname.ReplaceAll(".root", "");
  printname += ".pdf";
  printpath += printname;

  return printpath;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double getDCAxyvalue(double pt) {
  double dcaxy = 0.0182 + 0.035 / pow(pt, 1.01);

  return dcaxy;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int templateFitter(const char* file, const char* file2) {
  SetStyle();

  std::cout << "file sim: " << file << std::endl;
  std::cout << "file exp: " << file2 << std::endl;

  Bool_t doDCAxy_templateFit = kTRUE;
  Bool_t doCPA_templateFit = kTRUE;

  // load simulation data
  TFile* _file0 = TFile::Open(file, "READ");
  _file0->ls();
  if (!_file0) return 0;

  TList* listSP = 0;
  _file0->GetObject("/PWGCF_PLFemto_0/SPdir_0", listSP);

  TList* listPID = 0;
  _file0->GetObject("/PWGCF_PLFemto_0/PIDdir_0", listPID);

  if (!listSP) {
    std::cout << "list object sim not found" << std::endl;
    return 0;
  }

  if (!listPID) {
    std::cout << "list object sim not found" << std::endl;
    return 0;
  }

  // Proton configurations
  const Int_t DCA_Ptbins = 20;
  const Float_t pt_proton_threshold_low = 0.5;  // GeV/#it{c}
  const Float_t pt_proton_threshold_up = 4.05;  // GeV/#it{c}
  const Float_t pt_proton_switch = 0.75;        // GeV/#it{c}
  const Float_t Delta_pt_proton =
      (pt_proton_threshold_up - pt_proton_threshold_low) / (Float_t)DCA_Ptbins;
  TH2F* hist_proton_DCAxy_DCAz_all[DCA_Ptbins];
  TH2F* hist_proton_DCAxy_DCAz_primary[DCA_Ptbins];
  TH2F* hist_proton_DCAxy_DCAz_secondary[DCA_Ptbins];
  TH2F* hist_proton_DCAxy_DCAz_secondaryLambda[DCA_Ptbins];
  TH2F* hist_proton_DCAxy_DCAz_secondarySigma[DCA_Ptbins];
  TH2F* hist_proton_DCAxy_DCAz_material[DCA_Ptbins];

  // V0 configurations
  const Float_t V0_CPA_norm_low = 0.99;
  const Float_t V0_CPA_norm_up = 1.;
  const Float_t V0_InvMass_width = 0.004;
  const Float_t V0_InvMass_lowfit = 1.115683 - V0_InvMass_width;
  const Float_t V0_InvMass_upfit = 1.115683 + V0_InvMass_width;
  const Int_t V0_InvMass_Ptbins = 8;
  const Float_t V0_InvMass_threshold_low = 0.3;
  const Float_t V0_InvMass_threshold_up = 4.3;
  const Float_t Delta_pt_V0 =
      (V0_InvMass_threshold_up - V0_InvMass_threshold_low) /
      (Float_t)V0_InvMass_Ptbins;

  for (int ptbin = 0; ptbin < DCA_Ptbins; ptbin++) {
    TString HistName = "fProtonDCAxyDCAzMCCase0PtBin";
    HistName += ptbin;
    hist_proton_DCAxy_DCAz_all[ptbin] =
        (TH2F*)listSP->FindObject(HistName.Data());

    HistName = "fProtonDCAxyDCAzMCCase1PtBin";
    HistName += ptbin;
    hist_proton_DCAxy_DCAz_primary[ptbin] =
        (TH2F*)listSP->FindObject(HistName.Data());

    HistName = "fProtonDCAxyDCAzMCCase2PtBin";
    HistName += ptbin;
    hist_proton_DCAxy_DCAz_secondary[ptbin] =
        (TH2F*)listSP->FindObject(HistName.Data());

    HistName = "fProtonDCAxyDCAzMCCase3PtBin";
    HistName += ptbin;
    hist_proton_DCAxy_DCAz_material[ptbin] =
        (TH2F*)listSP->FindObject(HistName.Data());

    HistName = "fProtonDCAxyDCAzMCPtBinLambdaPtBin";
    HistName += ptbin;
    hist_proton_DCAxy_DCAz_secondaryLambda[ptbin] =
        (TH2F*)listSP->FindObject(HistName.Data());

    HistName = "fProtonDCAxyDCAzMCPtBinSigmaPtBin";
    HistName += ptbin;
    hist_proton_DCAxy_DCAz_secondarySigma[ptbin] =
        (TH2F*)listSP->FindObject(HistName.Data());
  }

  const Int_t CPA_Ptbins = 8;
  TH1F* hist_lambda_CPA[CPA_Ptbins];
  TH1F* hist_lambda_CPA_primary[CPA_Ptbins];
  TH1F* hist_lambda_CPA_primary_normalized_ExpEvts[CPA_Ptbins];
  TH1F* hist_lambda_CPA_secondary[CPA_Ptbins];
  TH1F* hist_lambda_CPA_secondary_normalized_ExpEvts[CPA_Ptbins];
  TH1F* hist_lambda_CPA_material[CPA_Ptbins];
  TH1F* hist_lambda_CPA_material_normalized_ExpEvts[CPA_Ptbins];
  TH1F* hist_lambda_CPA_bkg[CPA_Ptbins];
  TH1F* hist_lambda_CPA_bkg_normalized_ExpEvts[CPA_Ptbins];

  for (Int_t ptbin = 0; ptbin < CPA_Ptbins; ptbin++) {
    TString HistName = "fLambdaCPAPtBin";
    HistName += ptbin;
    hist_lambda_CPA[ptbin] = (TH1F*)listSP->FindObject(HistName.Data());

    HistName = "fLambdaCPAPrimaryPtBin";
    HistName += ptbin;
    hist_lambda_CPA_primary[ptbin] = (TH1F*)listSP->FindObject(HistName.Data());

    HistName = "fLambdaCPASecondaryPtBin";
    HistName += ptbin;
    hist_lambda_CPA_secondary[ptbin] =
        (TH1F*)listSP->FindObject(HistName.Data());

    HistName = "fLambdaCPAMaterialPtBin";
    HistName += ptbin;
    hist_lambda_CPA_material[ptbin] =
        (TH1F*)listSP->FindObject(HistName.Data());

    HistName = "fLambdaCPABkgPtBin";
    HistName += ptbin;
    hist_lambda_CPA_bkg[ptbin] = (TH1F*)listSP->FindObject(HistName.Data());
  }

  const Int_t Purity_Ptbins = 33;

  // load experimental data:
  TFile* _file2 = TFile::Open(file2, "READ");
  if (!_file2) return 0;
  _file2->ls();

  TList* listPID_exp = 0;
  _file2->GetObject("/PWGCF_PLFemto_0/PIDdir_0", listPID_exp);

  TList* listSP_exp = 0;
  _file2->GetObject("/PWGCF_PLFemto_0/SPdir_0", listSP_exp);

  if (!listPID_exp) {
    std::cout << "list object exp not found" << std::endl;
    return 0;
  }

  if (!listSP_exp) {
    std::cout << "list object exp not found" << std::endl;
    return 0;
  }

  // DCAxy
  TH2F* hist_proton_DCAxy_DCAz_all_exp[DCA_Ptbins];

  for (int ptbin = 0; ptbin < DCA_Ptbins; ptbin++) {
    TString HistName = "fProtonDCAxyDCAzPt";
    HistName += ptbin;
    hist_proton_DCAxy_DCAz_all_exp[ptbin] =
        (TH2F*)listSP_exp->FindObject(HistName.Data());
    hist_proton_DCAxy_DCAz_all_exp[ptbin]->Sumw2();
  }

  // CPA
  TH1F* hist_lambda_CPA_exp[CPA_Ptbins];
  TH1F* hist_lambda_InvMass[CPA_Ptbins];

  for (Int_t ptbin = 0; ptbin < CPA_Ptbins; ptbin++) {
    TString HistName = "fLambdaCPAPtBin";
    HistName += ptbin;
    hist_lambda_CPA_exp[ptbin] = (TH1F*)listSP_exp->FindObject(HistName.Data());

    if (!hist_lambda_CPA_exp[ptbin]) std::cout << "histo missing" << std::endl;

    if (!hist_lambda_CPA_primary[ptbin])
      std::cout << "histo missing" << std::endl;

    HistName = "hist_Lambda_CPA_primary_normalized_ExpEvts_";
    HistName += ptbin;
    hist_lambda_CPA_primary_normalized_ExpEvts[ptbin] =
        (TH1F*)hist_lambda_CPA_primary[ptbin]->Clone(HistName.Data());
    Float_t scale_exp_sim =
        hist_lambda_CPA_exp[ptbin]->Integral(
            hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_low),
            hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_up)) /
        hist_lambda_CPA_primary[ptbin]->Integral(
            hist_lambda_CPA_primary[ptbin]->FindBin(V0_CPA_norm_low),
            hist_lambda_CPA_primary[ptbin]->FindBin(V0_CPA_norm_up));
    std::cout << "scale: " << scale_exp_sim << std::endl;
    hist_lambda_CPA_primary_normalized_ExpEvts[ptbin]->Scale(scale_exp_sim);

    HistName = "hist_Lambda_CPA_secondary_normalized_ExpEvts_";
    HistName += ptbin;
    hist_lambda_CPA_secondary_normalized_ExpEvts[ptbin] =
        (TH1F*)hist_lambda_CPA_secondary[ptbin]->Clone(HistName.Data());
    scale_exp_sim =
        hist_lambda_CPA_exp[ptbin]->Integral(
            hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_low),
            hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_up)) /
        hist_lambda_CPA_secondary[ptbin]->Integral(
            hist_lambda_CPA_secondary[ptbin]->FindBin(V0_CPA_norm_low),
            hist_lambda_CPA_secondary[ptbin]->FindBin(V0_CPA_norm_up));
    hist_lambda_CPA_secondary_normalized_ExpEvts[ptbin]->Scale(scale_exp_sim);

    HistName = "hist_Lambda_CPA_material_normalized_ExpEvts_";
    HistName += ptbin;
    hist_lambda_CPA_material_normalized_ExpEvts[ptbin] =
        (TH1F*)hist_lambda_CPA_material[ptbin]->Clone(HistName.Data());
    scale_exp_sim =
        hist_lambda_CPA_exp[ptbin]->Integral(
            hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_low),
            hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_up)) /
        hist_lambda_CPA_material[ptbin]->Integral(
            hist_lambda_CPA_material[ptbin]->FindBin(V0_CPA_norm_low),
            hist_lambda_CPA_material[ptbin]->FindBin(V0_CPA_norm_up));
    hist_lambda_CPA_material_normalized_ExpEvts[ptbin]->Scale(scale_exp_sim);

    HistName = "hist_Lambda_CPA_bkg_normalized_ExpEvts_";
    HistName += ptbin;
    hist_lambda_CPA_bkg_normalized_ExpEvts[ptbin] =
        (TH1F*)hist_lambda_CPA_bkg[ptbin]->Clone(HistName.Data());
    scale_exp_sim = hist_lambda_CPA_exp[ptbin]->Integral(
                        hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_low),
                        hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_up)) /
                    hist_lambda_CPA_bkg[ptbin]->Integral(
                        hist_lambda_CPA_bkg[ptbin]->FindBin(V0_CPA_norm_low),
                        hist_lambda_CPA_bkg[ptbin]->FindBin(V0_CPA_norm_up));
    hist_lambda_CPA_bkg_normalized_ExpEvts[ptbin]->Scale(scale_exp_sim);

    HistName = "fInvMassLambdawCutsPtBin";
    HistName += ptbin;

    hist_lambda_InvMass[ptbin] = (TH1F*)listSP_exp->FindObject(HistName.Data());
    hist_lambda_InvMass[ptbin]->Sumw2();
    hist_lambda_InvMass[ptbin]->Scale(
        1. / hist_lambda_InvMass[ptbin]->GetBinWidth(1));  // Scale to bin width
  }

  // single particle plots:
  TH1F* hist_proton_pt_spectrum = (TH1F*)listSP_exp->FindObject(
      "fProtonPt");  // pt spectrum of the produced protons from exp
  hist_proton_pt_spectrum->Scale(1. / hist_proton_pt_spectrum->GetBinWidth(1));
  TH1F* hist_Lambda_pt_spectrum = (TH1F*)listSP_exp->FindObject(
      "fLambdaPt");  // pt spectrum of the produced Lambdas from exp
  hist_Lambda_pt_spectrum->Scale(1. / hist_Lambda_pt_spectrum->GetBinWidth(1));

  // Purity
  TH1F* hist_proton_purity = new TH1F(
      "hist_proton_purity", "Purity of protons estimated from Monte Carlo",
      Purity_Ptbins, pt_proton_threshold_low, pt_proton_threshold_up);
  TH1F* hist_proton_all = new TH1F(
      "hist_proton_purity2", "Purity of protons estimated from Monte Carlo",
      Purity_Ptbins, pt_proton_threshold_low, pt_proton_threshold_up);

  Double_t purity_proton_ptweighted = 0.;
  Double_t purity_proton_weights = 0.;

  TString HistName = "fProtonsCorrectlyIdentified";
  hist_proton_purity = (TH1F*)listPID->FindObject(HistName.Data());
  hist_proton_purity->Sumw2();

  HistName = "fProtonsTotallyIdentified";
  hist_proton_all = (TH1F*)listPID->FindObject(HistName.Data());
  hist_proton_all->Sumw2();
  hist_proton_purity->Divide(hist_proton_all);

  for (int ptbin = 0; ptbin < hist_proton_purity->GetEntries(); ptbin++) {
    Float_t purity = hist_proton_purity->GetBinContent(ptbin + 1);

    Float_t ptval = hist_proton_purity->GetBinCenter(ptbin + 1);

    purity_proton_ptweighted += purity *
                                hist_proton_pt_spectrum->GetBinContent(
                                    hist_proton_pt_spectrum->FindBin(ptval));
    purity_proton_weights += hist_proton_pt_spectrum->GetBinContent(
        hist_proton_pt_spectrum->FindBin(ptval));
  }

  purity_proton_ptweighted /= purity_proton_weights;

  TLine* horline = new TLine();
  horline->SetLineWidth(2);
  horline->SetLineStyle(2);

  // start plotting:
  TCanvas* Can_proton_DCAxy_DCAz_ptbins =
      new TCanvas("Can_proton_DCAxy_DCAz_ptbins",
                  "Can_proton_DCAxy_DCAz_ptbins", 0, 0, 1300, 900);
  Can_proton_DCAxy_DCAz_ptbins->Divide(4, 5);
  for (int i = 0; i < DCA_Ptbins; i++) {
    Can_proton_DCAxy_DCAz_ptbins->cd(i + 1);
    gPad->SetLogz();
    DrawHist(hist_proton_DCAxy_DCAz_all[i], "", "DCA_{xy} (cm)", "DCA_{z} (cm)",
             0);
  }

  TCanvas* Can_proton_DCAxy_DCAz_ptbins_exp =
      new TCanvas("Can_proton_DCAxy_DCAz_ptbins_exp",
                  "Can_proton_DCAxy_DCAz_ptbins_exp", 0, 0, 1300, 900);
  Can_proton_DCAxy_DCAz_ptbins_exp->Divide(4, 5);
  for (int i = 0; i < DCA_Ptbins; i++) {
    Can_proton_DCAxy_DCAz_ptbins_exp->cd(i + 1);
    gPad->SetLogz();
    DrawHist(hist_proton_DCAxy_DCAz_all_exp[i], "", "DCA_{xy} (cm)",
             "DCA_{z} (cm)", 0);
  }

  // DCAxy template fits:
  TH1D* hist_proton_DCAxy_projection_all_exp[DCA_Ptbins];
  TH1D* hist_proton_DCAxy_projection_all_normalized_Entries
      [DCA_Ptbins];  // normalized to same entries as exp
  TH1D* hist_proton_DCAxy_projection_primary_normalized_ExpEvts
      [DCA_Ptbins];  // normalized to same entries as exp
  TH1D* hist_proton_DCAxy_projection_secondary_normalized_ExpEvts
      [DCA_Ptbins];  // normalized to same entries as exp
  TH1D* hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts
      [DCA_Ptbins];  // normalized to same entries as exp
  TH1D* hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts
      [DCA_Ptbins];  // normalized to same entries as exp
  TH1D* hist_proton_DCAxy_projection_material_normalized_ExpEvts
      [DCA_Ptbins];  // normalized to same entries as exp
  Float_t NEvts_exp_ptbin[DCA_Ptbins];

  const Double_t DCAxy_bestval = 0.1;  // cm
  const Double_t DCAz_bestval = 0.2;   // cm

  // const Double_t DCAxy_bestval = 2.4; //cm
  // const Double_t DCAz_bestval = 3.2; //cm

  Double_t protons_tot_modelPred[DCA_Ptbins];        // model prediction before
                                                     // template fit
  Double_t protons_primary_modelPred[DCA_Ptbins];    // model prediction before
                                                     // template fit
  Double_t protons_secondary_modelPred[DCA_Ptbins];  // model prediction before
                                                     // template fit
  Double_t protons_material_modelPred[DCA_Ptbins];   // model prediction before
                                                     // template fit

  TH1D* hist_proton_fractions_primary_modelPred =
      new TH1D("hist_proton_fractions_primary_modelPred",
               "Fraction of primary protons_modelPred", DCA_Ptbins, 0.5,
               4.05);  // model prediction
  TH1D* hist_proton_fractions_secondary_modelPred =
      new TH1D("hist_proton_fractions_secondary_modelPred",
               "Fraction of secondary protons_modelPred", DCA_Ptbins, 0.5,
               4.05);  // model prediction
  TH1D* hist_proton_fractions_material_modelPred =
      new TH1D("hist_proton_fractions_material_modelPred",
               "Fraction of material protons_modelPred", DCA_Ptbins, 0.5,
               4.05);  // model prediction

  TCanvas* Can_proton_DCAxy_projections_ptbins_expSimComp = new TCanvas(
      "Can_proton_DCAxy_projections_ptbins_expSimComp",
      "Can_proton_DCAxy_projections_ptbins_expSimComp", 0, 0, 1300, 900);
  Can_proton_DCAxy_projections_ptbins_expSimComp->Divide(4, 5);
  TCanvas* Can_proton_DCAxy_projections_ptbins_Simdecomp = new TCanvas(
      "Can_proton_DCAxy_projections_ptbins_Simdecomp",
      "Can_proton_DCAxy_projections_ptbins_Simdecomp", 0, 0, 1300, 900);
  Can_proton_DCAxy_projections_ptbins_Simdecomp->Divide(4, 5);
  TCanvas* Can_proton_DCAxy_projections_ptbins_secondaryComp = new TCanvas(
      "Can_proton_DCAxy_projections_ptbins_secondaryComp",
      "Can_proton_DCAxy_projections_ptbins_secondaryComp", 0, 0, 1300, 900);
  Can_proton_DCAxy_projections_ptbins_secondaryComp->Divide(4, 5);

  for (int i = 0; i < DCA_Ptbins; i++) {
    TString Projection_Name = "hist_proton_DCAxy_projection_all_exp_ptbin_";
    Projection_Name += i;
    hist_proton_DCAxy_projection_all_exp[i] =
        hist_projection(hist_proton_DCAxy_DCAz_all_exp[i], -DCAz_bestval,
                        DCAz_bestval, "X", Projection_Name.Data());
    NEvts_exp_ptbin[i] = hist_proton_DCAxy_projection_all_exp[i]->Integral();

    Projection_Name = "hist_proton_DCAxy_projection_all_ptbin_";
    Projection_Name += i;
    TH1D* hist_proton_DCAxy_projection_all =
        hist_projection(hist_proton_DCAxy_DCAz_all[i], -DCAz_bestval,
                        DCAz_bestval, "X", Projection_Name.Data());
    Projection_Name += "_normalized";
    // Scale exp to sim value:
    Double_t scale_exp_sim =
        hist_proton_DCAxy_projection_all_exp[i]->Integral() /
        hist_proton_DCAxy_projection_all->Integral();
    hist_proton_DCAxy_projection_all_normalized_Entries[i] =
        (TH1D*)hist_proton_DCAxy_projection_all->Clone(Projection_Name.Data());
    hist_proton_DCAxy_projection_all_normalized_Entries[i]->Scale(
        scale_exp_sim);  // exp sim to same entries as sim

    protons_tot_modelPred[i] = hist_proton_DCAxy_projection_all->Integral(
        hist_proton_DCAxy_projection_all->FindBin(-DCAxy_bestval),
        hist_proton_DCAxy_projection_all->FindBin(DCAxy_bestval));

    Projection_Name = "hist_proton_DCAxy_projection_primary_ptbin_";
    Projection_Name += i;
    TH1D* hist_proton_DCAxy_projection_primary =
        hist_projection(hist_proton_DCAxy_DCAz_primary[i], -DCAz_bestval,
                        DCAz_bestval, "X", Projection_Name.Data());
    Projection_Name += "_normalized";
    scale_exp_sim = hist_proton_DCAxy_projection_all_exp[i]->Integral() /
                    hist_proton_DCAxy_projection_primary->Integral();
    hist_proton_DCAxy_projection_primary_normalized_ExpEvts[i] =
        (TH1D*)hist_proton_DCAxy_projection_primary->Clone(
            Projection_Name.Data());
    hist_proton_DCAxy_projection_primary_normalized_ExpEvts[i]->Scale(
        scale_exp_sim);

    protons_primary_modelPred[i] =
        hist_proton_DCAxy_projection_primary->Integral(
            hist_proton_DCAxy_projection_primary->FindBin(-DCAxy_bestval),
            hist_proton_DCAxy_projection_primary->FindBin(DCAxy_bestval));
    hist_proton_fractions_primary_modelPred->SetBinContent(
        i + 1, protons_primary_modelPred[i] / protons_tot_modelPred[i]);

    Projection_Name = "hist_proton_DCAxy_projection_secondary_ptbin_";
    Projection_Name += i;
    TH1D* hist_proton_DCAxy_projection_secondary =
        hist_projection(hist_proton_DCAxy_DCAz_secondary[i], -DCAz_bestval,
                        DCAz_bestval, "X", Projection_Name.Data());
    Projection_Name += "_normalized";
    scale_exp_sim = hist_proton_DCAxy_projection_all_exp[i]->Integral() /
                    hist_proton_DCAxy_projection_secondary->Integral();
    hist_proton_DCAxy_projection_secondary_normalized_ExpEvts[i] =
        (TH1D*)hist_proton_DCAxy_projection_secondary->Clone(
            Projection_Name.Data());
    hist_proton_DCAxy_projection_secondary_normalized_ExpEvts[i]->Scale(
        scale_exp_sim);

    protons_secondary_modelPred[i] =
        hist_proton_DCAxy_projection_secondary->Integral(
            hist_proton_DCAxy_projection_secondary->FindBin(-DCAxy_bestval),
            hist_proton_DCAxy_projection_secondary->FindBin(DCAxy_bestval));
    hist_proton_fractions_secondary_modelPred->SetBinContent(
        i + 1, protons_secondary_modelPred[i] / protons_tot_modelPred[i]);

    // project Sigma+ and Lambda secondary contributions seperately with the
    // same scale factor calculated before for comparison:

    Projection_Name = "hist_proton_DCAxy_projection_secondaryLambda_ptbin_";
    Projection_Name += i;
    TH1D* hist_proton_DCAxy_projection_secondaryLambda = hist_projection(
        hist_proton_DCAxy_DCAz_secondaryLambda[i], -DCAz_bestval, DCAz_bestval,
        "X", Projection_Name.Data());
    Projection_Name += "_normalized";
    hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts[i] =
        (TH1D*)hist_proton_DCAxy_projection_secondaryLambda->Clone(
            Projection_Name.Data());
    hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts[i]->Scale(
        scale_exp_sim);

    Projection_Name = "hist_proton_DCAxy_projection_secondarySigma_ptbin_";
    Projection_Name += i;
    TH1D* hist_proton_DCAxy_projection_secondarySigma =
        hist_projection(hist_proton_DCAxy_DCAz_secondarySigma[i], -DCAz_bestval,
                        DCAz_bestval, "X", Projection_Name.Data());
    Projection_Name += "_normalized";
    hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts[i] =
        (TH1D*)hist_proton_DCAxy_projection_secondarySigma->Clone(
            Projection_Name.Data());
    hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts[i]->Scale(
        scale_exp_sim);

    Projection_Name = "hist_proton_DCAxy_projection_material_ptbin_";
    Projection_Name += i;
    TH1D* hist_proton_DCAxy_projection_material =
        hist_projection(hist_proton_DCAxy_DCAz_material[i], -DCAz_bestval,
                        DCAz_bestval, "X", Projection_Name.Data());
    Projection_Name += "_normalized";
    scale_exp_sim = hist_proton_DCAxy_projection_all_exp[i]->Integral() /
                    hist_proton_DCAxy_projection_material->Integral();
    hist_proton_DCAxy_projection_material_normalized_ExpEvts[i] =
        (TH1D*)hist_proton_DCAxy_projection_material->Clone(
            Projection_Name.Data());
    hist_proton_DCAxy_projection_material_normalized_ExpEvts[i]->Scale(
        scale_exp_sim);

    protons_material_modelPred[i] =
        hist_proton_DCAxy_projection_material->Integral(
            hist_proton_DCAxy_projection_material->FindBin(-DCAxy_bestval),
            hist_proton_DCAxy_projection_material->FindBin(DCAxy_bestval));
    hist_proton_fractions_material_modelPred->SetBinContent(
        i + 1, protons_material_modelPred[i] / protons_tot_modelPred[i]);

    // Plotting:

    Can_proton_DCAxy_projections_ptbins_expSimComp->cd(i + 1);
    gPad->SetLogy();
    hist_proton_DCAxy_projection_all_exp[i]->GetXaxis()->SetRangeUser(-3, 3);
    DrawHist((TH1F*)hist_proton_DCAxy_projection_all_exp[i],
             (TH1F*)hist_proton_DCAxy_projection_all_normalized_Entries[i], "",
             "DCA_{xy} (cm)", "Counts", 1, 1, 1, 4, 1, kFALSE);

    Can_proton_DCAxy_projections_ptbins_Simdecomp->cd(i + 1);
    gPad->SetLogy();
    hist_proton_DCAxy_projection_all->GetXaxis()->SetRangeUser(-3, 3);
    DrawHist((TH1F*)hist_proton_DCAxy_projection_all,
             (TH1F*)hist_proton_DCAxy_projection_primary,
             (TH1F*)hist_proton_DCAxy_projection_secondary,
             (TH1F*)hist_proton_DCAxy_projection_material, "", "DCA_{xy} (cm)",
             "Counts", 1, 1, 1, 4, 1, 2, 1, 3, 1,
             kFALSE);  // The sum of all contributions gives the total
                       // contribution: checked

    Can_proton_DCAxy_projections_ptbins_secondaryComp->cd(i + 1);
    gPad->SetLogy();

    TString CloneName = "histSum_";
    CloneName += i;

    TH1D* histSum =
        (TH1D*)
            hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts[i]
                ->Clone(CloneName.Data());
    histSum->Add(
        hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts[i]);

    DrawHist(
        (TH1F*)hist_proton_DCAxy_projection_secondary_normalized_ExpEvts[i],
        (TH1F*)
            hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts[i],
        (TH1F*)
            hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts[i],
        "", "DCA_{xy} (cm)", "Counts", 1, 1, 1, 4, 1, 2, 1, kFALSE);

    histSum->DrawCopy("same");
  }

  TCanvas* Can_proton_DCAxy_representation =
      new TCanvas("Can_proton_DCAxy_representation",
                  "Can_proton_DCAxy_representation", 0, 0, 1300, 900);
  Can_proton_DCAxy_representation->cd();
  gPad->SetLogy();
  DrawHist((TH1F*)hist_proton_DCAxy_projection_all_exp[2],
           (TH1F*)hist_proton_DCAxy_projection_primary_normalized_ExpEvts[2],
           (TH1F*)hist_proton_DCAxy_projection_secondary_normalized_ExpEvts[2],
           (TH1F*)hist_proton_DCAxy_projection_material_normalized_ExpEvts[2],
           "", "DCA_{xy} (cm)", "Counts",  2, 1, 1, 2, 1, 3, 1, 4, 1,
           kFALSE);  // The sum of all contributions gives the total
                     // contribution: checked

  // test MC template fit:
  TH1D* hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight
      [DCA_Ptbins];
  TH1D* hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight
      [DCA_Ptbins];
  TH1D* hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts_weight
      [DCA_Ptbins];
  TH1D* hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts_weight
      [DCA_Ptbins];
  TH1D* hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight
      [DCA_Ptbins];
  TH1D* hist_proton_DCAxy_MC_templatesum[DCA_Ptbins];
  TH1D* hist_proton_secondary_DCAxy_MC_templatesum[DCA_Ptbins];
  TH1D* hist_proton_DCAxy_MC_ratio_data_MC[DCA_Ptbins];
  TH1D* hist_proton_fractions_primary =
      new TH1D("hist_proton_fractions_primary", "Fraction of primary protons",
               DCA_Ptbins, 0.5, 4.05);
  TH1D* hist_proton_fractions_secondary =
      new TH1D("hist_proton_fractions_secondary",
               "Fraction of secondary protons", DCA_Ptbins, 0.5, 4.05);
  TH1D* hist_proton_fractions_material =
      new TH1D("hist_proton_fractions_material", "Fraction of material protons",
               DCA_Ptbins, 0.5, 4.05);
  TH1D* hist_proton_fractions_primary_wDCAxyCut = new TH1D(
      "hist_proton_fractions_primary_wDCAxyCut",
      "Fraction of primary protons with DCAxyCut", DCA_Ptbins, 0.5, 4.05);
  TH1D* hist_proton_fractions_secondary_wDCAxyCut = new TH1D(
      "hist_proton_fractions_secondary_wDCAxyCut",
      "Fraction of secondary protons with DCAxyCut", DCA_Ptbins, 0.5, 4.05);
  TH1D* hist_proton_fractions_secondaryLambda_wDCAxyCut =
      new TH1D("hist_proton_fractions_secondaryLambda_wDCAxyCut",
               "Fraction of secondary protons with DCAxyCut coming from Lambda",
               DCA_Ptbins, 0.5, 4.05);
  TH1D* hist_proton_fractions_secondarySigma_wDCAxyCut =
      new TH1D("hist_proton_fractions_secondarySigma_wDCAxyCut",
               "Fraction of secondary protons with DCAxyCut coming from Sigma",
               DCA_Ptbins, 0.5, 4.05);
  TH1D* hist_proton_fractions_secondarySum_wDCAxyCut =
      new TH1D("hist_proton_fractions_secondarySum_wDCAxyCut",
               "Fraction of secondary protons with DCAxyCut coming from Sigma "
               "and Lambda",
               DCA_Ptbins, 0.5, 4.05);
  TH1D* hist_proton_fractions_material_wDCAxyCut = new TH1D(
      "hist_proton_fractions_material_wDCAxyCut",
      "Fraction of material protons with DCAxyCut", DCA_Ptbins, 0.5, 4.05);
  TH1D* hist_protonproton_lambdapar = new TH1D(
      "hist_protonproton_lambdapar", "Lambda parameter for proton-proton",
      DCA_Ptbins, pt_proton_threshold_low, pt_proton_threshold_up);

  TCanvas* Can_proton_DCAxy_templateFit =
      new TCanvas("Can_proton_DCAxy_templateFit",
                  "Can_proton_DCAxy_templateFit", 0, 0, 1300, 900);
  Can_proton_DCAxy_templateFit->Divide(4, 5);

  TCanvas* Can_proton_DCAxy_templateFit_SecondaryDecomp = new TCanvas(
      "Can_proton_DCAxy_templateFit_SecondaryDecomp",
      "Can_proton_DCAxy_templateFit_SecondaryDecomp", 0, 0, 1300, 900);
  Can_proton_DCAxy_templateFit_SecondaryDecomp->Divide(4, 5);

  TCanvas* Can_proton_ratio_data_MC = new TCanvas(
      "Can_proton_ratio_data_MC", "Can_proton_ratio_data_MC", 0, 0, 1300, 900);
  Can_proton_ratio_data_MC->Divide(4, 5);

  Double_t protons_tot[DCA_Ptbins];
  Double_t protons_primary[DCA_Ptbins];
  Double_t protons_secondary[DCA_Ptbins];
  Double_t protons_secondaryLambda[DCA_Ptbins];
  Double_t protons_secondarySigma[DCA_Ptbins];
  Double_t protons_material[DCA_Ptbins];

  // Parameters to obtain the average weighted by the pt production spectrum of
  // protons
  Double_t ptweight_proton_pt = 0;
  Double_t fraction_proton_primary_ptweighted = 0.;
  Double_t fraction_Lambda_feeddown_ptweighted = 0.;
  Double_t fraction_proton_material_ptweighted = 0.;

  TLegend* leg_DCA = new TLegend(0.1, 0.7, 0.48, 0.9);
  leg_DCA->SetTextSize(0.055);

  for (int ptbin = 0; ptbin < DCA_Ptbins; ptbin++) {
    std::cout << "*******************************" << std::endl;
    std::cout << "*******************************" << std::endl;
    std::cout << "ptbin: " << ptbin << std::endl;
    std::cout << "*******************************" << std::endl;
    std::cout << "*******************************" << std::endl;

    TObjArray* mcDCA_array =
        new TObjArray(3);  // MC histograms are put in this array
    mcDCA_array->Add(
        hist_proton_DCAxy_projection_primary_normalized_ExpEvts[ptbin]);
    mcDCA_array->Add(
        hist_proton_DCAxy_projection_secondary_normalized_ExpEvts[ptbin]);
    mcDCA_array->Add(
        hist_proton_DCAxy_projection_material_normalized_ExpEvts[ptbin]);

    TFractionFitter* fit =
        new TFractionFitter(hist_proton_DCAxy_projection_all_exp[ptbin],
                            mcDCA_array);  // initialise
    // fit->Constrain(0,0.,1.);               // constrain fraction 0 to be
    // between 0 and 1
    // fit->Constrain(1,0.,1.);               // constrain fraction 1 to be
    // between 0 and 1
    // fit->Constrain(2,0.,1.);               // constrain fraction 2 to be
    // between 0 and 1
    if(ptbin == 13)  fit->Constrain(0,0.,1.);
    if(ptbin == 15)  fit->Constrain(0,0.,1.);
    if(ptbin == 18)  fit->Constrain(0,0.,10.);
    fit->SetRangeX(hist_proton_DCAxy_projection_all_exp[ptbin]->FindBin(-2.4),
                   hist_proton_DCAxy_projection_all_exp[ptbin]->FindBin(
                       2.4));  // use only the first 15 bins in the fit

    Int_t status = -1;

    if (doDCAxy_templateFit) status = fit->Fit();  // perform the fit

    if (status == 0) {
      // obtain weights:
      Double_t par0, errpar0, par1, errpar1, par2, errpar2;
      fit->GetResult(0, par0, errpar0);
      fit->GetResult(1, par1, errpar1);
      fit->GetResult(2, par2, errpar2);

      Double_t weight0 = par0;  // times Nevts?
      Double_t weight1 = par1;
      Double_t weight2 = par2;

      hist_proton_fractions_primary->SetBinContent(ptbin + 1, par0);
      hist_proton_fractions_primary->SetBinError(ptbin + 1, errpar0);
      hist_proton_fractions_secondary->SetBinContent(ptbin + 1, par1);
      hist_proton_fractions_secondary->SetBinError(ptbin + 1, errpar1);
      hist_proton_fractions_material->SetBinContent(ptbin + 1, par2);
      hist_proton_fractions_material->SetBinError(ptbin + 1, errpar2);

      if (errpar0 > 0.1) errpar0 = 0.f;
      if (errpar1 > 0.1) errpar1 = 0.f;
      if (errpar2 > 0.1) errpar2 = 0.f;

      TString CloneName =
          "hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight_";
      CloneName += ptbin;
      hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight[ptbin] =
          (TH1D*)hist_proton_DCAxy_projection_primary_normalized_ExpEvts[ptbin]
              ->Clone(CloneName.Data());
      hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight[ptbin]
          ->Scale(weight0);

      CloneName =
          "hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight_";
      CloneName += ptbin;
      hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight[ptbin] =
          (TH1D*)
              hist_proton_DCAxy_projection_secondary_normalized_ExpEvts[ptbin]
                  ->Clone(CloneName.Data());
      hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight[ptbin]
          ->Scale(weight1);

      // Scale also the individual Sigma+/Lambda contributions:

      CloneName =
          "hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts_"
          "weight_";
      CloneName += ptbin;
      hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts_weight
          [ptbin] =
              (TH1D*)
                  hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts
                      [ptbin]
                          ->Clone(CloneName.Data());
      hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts_weight
          [ptbin]
              ->Scale(weight1);

      CloneName =
          "hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts_"
          "weight_";
      CloneName += ptbin;
      hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts_weight
          [ptbin] =
              (TH1D*)
                  hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts
                      [ptbin]
                          ->Clone(CloneName.Data());
      hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts_weight
          [ptbin]
              ->Scale(weight1);

      CloneName = "hist_proton_secondary_DCAxy_MC_templatesum_";
      CloneName += ptbin;

      hist_proton_secondary_DCAxy_MC_templatesum[ptbin] =
          (TH1D*)
              hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts_weight
                  [ptbin]
                      ->Clone(CloneName.Data());
      ;
      hist_proton_secondary_DCAxy_MC_templatesum[ptbin]->Add(
          hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts_weight
              [ptbin]);
      hist_proton_secondary_DCAxy_MC_templatesum[ptbin]->SetLineColor(2);

      CloneName =
          "hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight_";
      CloneName += ptbin;
      hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight[ptbin] =
          (TH1D*)hist_proton_DCAxy_projection_material_normalized_ExpEvts[ptbin]
              ->Clone(CloneName.Data());
      hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight[ptbin]
          ->Scale(weight2);

      CloneName = "hist_proton_DCAxy_MC_templatesum_";
      CloneName += ptbin;
      hist_proton_DCAxy_MC_templatesum[ptbin] =
          (TH1D*)hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight
              [ptbin]
                  ->Clone(CloneName.Data());
      hist_proton_DCAxy_MC_templatesum[ptbin]->Add(
          hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight
              [ptbin]);
      hist_proton_DCAxy_MC_templatesum[ptbin]->Add(
          hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight
              [ptbin]);

      CloneName = "hist_proton_DCAxy_MC_ratio_data_MC_";
      CloneName += ptbin;
      hist_proton_DCAxy_MC_ratio_data_MC[ptbin] =
          (TH1D*)hist_proton_DCAxy_MC_templatesum[ptbin]->Clone(
              CloneName.Data());
      hist_proton_DCAxy_MC_ratio_data_MC[ptbin]->Divide(
          (TH1F*)hist_proton_DCAxy_projection_all_exp[ptbin]);

      // Estimate fractions with DCAxy cut:
      protons_tot[ptbin] =
          hist_proton_DCAxy_projection_all_exp[ptbin]->Integral(
              hist_proton_DCAxy_projection_all_exp[ptbin]->FindBin(
                  -DCAxy_bestval),
              hist_proton_DCAxy_projection_all_exp[ptbin]->FindBin(
                  DCAxy_bestval));
      protons_primary[ptbin] =
          hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight[ptbin]
              ->Integral(
                  hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight
                      [ptbin]
                          ->FindBin(-DCAxy_bestval),
                  hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight
                      [ptbin]
                          ->FindBin(DCAxy_bestval));
      protons_secondary[ptbin] =
          hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight[ptbin]
              ->Integral(
                  hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight
                      [ptbin]
                          ->FindBin(-DCAxy_bestval),
                  hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight
                      [ptbin]
                          ->FindBin(DCAxy_bestval));
      protons_secondaryLambda[ptbin] =
          hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts_weight
              [ptbin]
                  ->Integral(
                      hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts_weight
                          [ptbin]
                              ->FindBin(-DCAxy_bestval),
                      hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts_weight
                          [ptbin]
                              ->FindBin(DCAxy_bestval));
      protons_secondarySigma[ptbin] =
          hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts_weight
              [ptbin]
                  ->Integral(
                      hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts_weight
                          [ptbin]
                              ->FindBin(-DCAxy_bestval),
                      hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts_weight
                          [ptbin]
                              ->FindBin(DCAxy_bestval));
      protons_material[ptbin] =
          hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight[ptbin]
              ->Integral(
                  hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight
                      [ptbin]
                          ->FindBin(-DCAxy_bestval),
                  hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight
                      [ptbin]
                          ->FindBin(DCAxy_bestval));
      hist_proton_fractions_primary_wDCAxyCut->SetBinContent(
          ptbin + 1, protons_primary[ptbin] / protons_tot[ptbin]);
      hist_proton_fractions_primary_wDCAxyCut->SetBinError(ptbin + 1, errpar0);
      hist_proton_fractions_secondary_wDCAxyCut->SetBinContent(
          ptbin + 1, protons_secondary[ptbin] / protons_tot[ptbin]);
      hist_proton_fractions_secondary_wDCAxyCut->SetBinError(ptbin + 1,
                                                             errpar1);
      hist_proton_fractions_secondaryLambda_wDCAxyCut->SetBinContent(
          ptbin + 1, protons_secondaryLambda[ptbin] / protons_secondary[ptbin]);
      hist_proton_fractions_secondaryLambda_wDCAxyCut->SetBinError(ptbin + 1,
                                                                   errpar1);
      hist_proton_fractions_secondarySigma_wDCAxyCut->SetBinContent(
          ptbin + 1, protons_secondarySigma[ptbin] / protons_secondary[ptbin]);
      hist_proton_fractions_secondarySigma_wDCAxyCut->SetBinError(ptbin + 1,
                                                                  errpar1);
      hist_proton_fractions_secondarySum_wDCAxyCut->SetBinContent(
          ptbin + 1,
          (protons_secondaryLambda[ptbin] + protons_secondarySigma[ptbin]) /
              protons_secondary[ptbin]);
      hist_proton_fractions_material_wDCAxyCut->SetBinContent(
          ptbin + 1, protons_material[ptbin] / protons_tot[ptbin]);
      hist_proton_fractions_material_wDCAxyCut->SetBinError(ptbin + 1, errpar2);

      Double_t fraction_primary_protons_squared =
          protons_primary[ptbin] / protons_tot[ptbin];
      fraction_primary_protons_squared *= fraction_primary_protons_squared;
      hist_protonproton_lambdapar->SetBinContent(
          ptbin + 1, fraction_primary_protons_squared);
      hist_protonproton_lambdapar->SetBinError(ptbin + 1, 0);

      // Calculate the average parameters needed for the analysis:
      //***********************************************************
      Double_t ptval = hist_proton_fractions_primary->GetBinCenter(ptbin + 1);

      fraction_proton_primary_ptweighted +=
          protons_primary[ptbin] / protons_tot[ptbin] *
          hist_proton_pt_spectrum->GetBinContent(
              hist_proton_pt_spectrum->FindBin(ptval));
      ptweight_proton_pt += hist_proton_pt_spectrum->GetBinContent(
          hist_proton_pt_spectrum->FindBin(ptval));

      fraction_Lambda_feeddown_ptweighted +=
          protons_secondary[ptbin] / protons_tot[ptbin] *
          hist_proton_pt_spectrum->GetBinContent(
              hist_proton_pt_spectrum->FindBin(ptval));
      fraction_proton_material_ptweighted +=
          protons_material[ptbin] / protons_tot[ptbin] *
          hist_proton_pt_spectrum->GetBinContent(
              hist_proton_pt_spectrum->FindBin(ptval));
      //***********************************************************

      Can_proton_DCAxy_templateFit->cd(ptbin + 1);
      gPad->SetLogy();
      TString HistTitle = print_ranges(pt_proton_threshold_low, Delta_pt_proton,
                                       ptbin, "#it{p}_{T}", "GeV/#it{c}");
      hist_proton_DCAxy_projection_all_exp[ptbin]->SetTitle(HistTitle);
      hist_proton_DCAxy_projection_all_exp[ptbin]->SetTitleSize(1.5);
      DrawHist(
          (TH1F*)hist_proton_DCAxy_projection_all_exp[ptbin],
          (TH1F*)hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight
              [ptbin],
          (TH1F*)
              hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight
                  [ptbin],
          (TH1F*)hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight
              [ptbin],
          HistTitle, "DCA_{xy} (cm)", "Counts", 2, 1, 1, 2, 1, 3, 1, 4, 1,
          kTRUE);  // The sum of all contributions gives the total
                    // contribution: checked
      hist_proton_DCAxy_MC_templatesum[ptbin]->SetLineWidth(2);
      hist_proton_DCAxy_MC_templatesum[ptbin]->SetLineColor(1);
      hist_proton_DCAxy_MC_templatesum[ptbin]->DrawCopy("same");
      // print_ranges(pt_proton_threshold_low,Delta_pt_proton,ptbin,"#it{p}_{T}","GeV/#it{c}",-2.,90.);
      if (ptbin == 0) {
        leg_DCA->AddEntry(
            hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight
                [ptbin],
            "primary proton", "f");
        leg_DCA->AddEntry(
            hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight
                [ptbin],
            "secondary proton", "f");
        leg_DCA->AddEntry(
            hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight
                [ptbin],
            "material proton", "f");
      }
      leg_DCA->Draw();


      Can_proton_DCAxy_templateFit_SecondaryDecomp->cd(ptbin + 1);
      gPad->SetLogy();

      DrawHist(
          (TH1F*)
              hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight
                  [ptbin],
          (TH1F*)
              hist_proton_DCAxy_projection_secondaryLambda_normalized_ExpEvts_weight
                  [ptbin],
          (TH1F*)
              hist_proton_DCAxy_projection_secondarySigma_normalized_ExpEvts_weight
                  [ptbin],
          "", "DCA_{xy} (cm)", "Counts", 1, 1, 1, 4, 1, 2, 1, kFALSE);

      hist_proton_secondary_DCAxy_MC_templatesum[ptbin]->SetLineColor(1);
      hist_proton_secondary_DCAxy_MC_templatesum[ptbin]->SetMarkerColor(1);
      hist_proton_secondary_DCAxy_MC_templatesum[ptbin]->DrawCopy("histo same");

      Can_proton_ratio_data_MC->cd(ptbin + 1);
      DrawHist((TH1F*)hist_proton_DCAxy_MC_ratio_data_MC[ptbin], "",
               "DCA_{xy} (cm)", "Counts", 1, 1, 1, kFALSE);
    }
  }
  // Can_proton_DCAxy_templateFit->Print(Create_path(alinotedir,"MC_template","DCAxy_template",""));

  Can_proton_DCAxy_templateFit->Print("ANplot/protonTemplateFit.pdf");

  fraction_proton_primary_ptweighted /=
      ptweight_proton_pt;  // average fraction of primary protons
  fraction_Lambda_feeddown_ptweighted /=
      ptweight_proton_pt;  // average fraction of Lambdas feeding to protons
  fraction_proton_material_ptweighted /=
      ptweight_proton_pt;  // average fraction of material protons

  TFile* outFile = new TFile("templateOutput.root", "RECREATE");

  const Int_t PIDBinf_for_repo = 1;

  TLegend* leg_proton_summary = new TLegend(0.3, 0.4, 0.78, 0.6);
  leg_proton_summary->SetTextSize(0.055);

  if (doDCAxy_templateFit) {
    TCanvas* Can_proton_DCAxy_representation_fit =
        new TCanvas("Can_proton_DCAxy_representation_fit",
                    "Can_proton_DCAxy_representation_fit");
    Can_proton_DCAxy_representation_fit->cd();
    gPad->SetLogy();
    DrawHist(
        (TH1F*)hist_proton_DCAxy_projection_all_exp[PIDBinf_for_repo],
        (TH1F*)hist_proton_DCAxy_projection_primary_normalized_ExpEvts_weight
            [PIDBinf_for_repo],
        (TH1F*)hist_proton_DCAxy_projection_secondary_normalized_ExpEvts_weight
            [PIDBinf_for_repo],
        (TH1F*)hist_proton_DCAxy_projection_material_normalized_ExpEvts_weight
            [PIDBinf_for_repo],
        "", "DCA_{xy} (cm)", "Counts",  1, 1, 1, 2, 1, 3, 1, 4, 1, kFALSE, true);
    hist_proton_DCAxy_MC_templatesum[PIDBinf_for_repo]->SetLineWidth(2);
    hist_proton_DCAxy_MC_templatesum[PIDBinf_for_repo]->SetLineColor(colors[1]);
    hist_proton_DCAxy_MC_templatesum[PIDBinf_for_repo]->DrawCopy("histo same");
    print_ranges(pt_proton_threshold_low, Delta_pt_proton, PIDBinf_for_repo,
                 "#it{p}_{T}", "GeV/#it{c}", -2.65, 2E6);
    Can_proton_DCAxy_representation_fit->Print("ANplot/protonTemplateFit_single.pdf");

    TF1* fraction_fit_primary = new TF1("fraction_fit_primary", "[0]", 0.5, 4.);
    fraction_fit_primary->SetParameter(0, 0.92);

    TCanvas* Can_proton_fractions = new TCanvas(
        "Can_proton_fractions", "Can_proton_fractions", 0, 0, 1300, 900);
    Can_proton_fractions->Divide(2, 2);
    Can_proton_fractions->cd(1);
    gPad->SetLogy();
    SetStyleHisto(hist_proton_pt_spectrum, 0,0);
    SetStyleHisto(hist_proton_purity, 0,0);
    hist_proton_purity->SetTitle("; #it{p}_{T} (GeV/#it{c}); Purity");
    hist_proton_pt_spectrum->GetXaxis()->SetRangeUser(0,
                                                      pt_proton_threshold_up);
    DrawHist((TH1F*)hist_proton_pt_spectrum, "", "#it{p}_{T} (GeV/#it{c})",
             "dN/d#it{p}_{T} 1/(GeV/#it{c})", 0, 1, 1, kFALSE);
    Can_proton_fractions->cd(2);
    hist_proton_fractions_primary_wDCAxyCut->GetYaxis()->SetRangeUser(-0.1,
                                                                      1.2);
    DrawHist((TH1F*)hist_proton_fractions_primary_wDCAxyCut,
             (TH1F*)hist_proton_fractions_secondary_wDCAxyCut,
             (TH1F*)hist_proton_fractions_material_wDCAxyCut, "",
             "#it{p}_{T} (GeV/#it{c})", "fraction", 0, 1, 1, 4, 1, 2, 1,
             kFALSE);
    horline->DrawLine(
        pt_proton_threshold_low, fraction_proton_primary_ptweighted,
        pt_proton_threshold_up, fraction_proton_primary_ptweighted);
    horline->SetLineColor(colors[0]);

    horline->DrawLine(
        pt_proton_threshold_low, fraction_Lambda_feeddown_ptweighted,
        pt_proton_threshold_up, fraction_Lambda_feeddown_ptweighted);
    horline->SetLineColor(colors[1]);

    horline->DrawLine(
        pt_proton_threshold_low, fraction_proton_material_ptweighted,
        pt_proton_threshold_up, fraction_proton_material_ptweighted);
    horline->SetLineColor(colors[0]);

    leg_proton_summary->AddEntry(hist_proton_fractions_primary_wDCAxyCut,
                                 "primary proton", "p");
    leg_proton_summary->AddEntry(hist_proton_fractions_secondary_wDCAxyCut,
                                 "secondary proton", "p");
    leg_proton_summary->AddEntry(hist_proton_fractions_material_wDCAxyCut,
                                 "material proton", "p");
    leg_proton_summary->Draw();

    hist_proton_fractions_primary_wDCAxyCut->Write();
    hist_proton_fractions_secondary_wDCAxyCut->Write();
    hist_proton_fractions_material_wDCAxyCut->Write();

    Can_proton_fractions->cd(3);
    hist_proton_purity->Draw("pe");
    hist_proton_purity->SetMarkerStyle(20);
    hist_proton_purity->SetLineColor(kBlack);
    //      DrawHist((TH1F*)hist_proton_purity,"","#it{p}_{T}
    //      (GeV/#it{c})","Purity(Proton)",0,1,1, kFALSE);
    horline->DrawLine(pt_proton_switch, purity_proton_ptweighted,
                      pt_proton_threshold_up, purity_proton_ptweighted);
    horline->SetLineColor(colors[2]);
    // Can_proton_fractions->Print(Create_path(alinotedir,"MC_template","Proton_summary",""));
    hist_proton_purity->Write();
    Can_proton_fractions->Print("ANplot/protonFraction.pdf");

    TLegend* SigmaLambdaPercentage = new TLegend(0.2, 0.6, 0.65, 0.8);
    SigmaLambdaPercentage->SetTextSize(0.055);


    SigmaLambdaPercentage->AddEntry(
        hist_proton_fractions_secondaryLambda_wDCAxyCut, "#Lambda percentage",
        "p");
    SigmaLambdaPercentage->AddEntry(
        hist_proton_fractions_secondarySigma_wDCAxyCut, "#Sigma percentage",
        "p");

    TCanvas* Can_LambdaSigma_fractions =
        new TCanvas("Can_LambdaSigma_fractions", "Can_LambdaSigma_fractions");
    Can_LambdaSigma_fractions->cd();
    DrawHist((TH1F*)hist_proton_fractions_secondaryLambda_wDCAxyCut,
             (TH1F*)hist_proton_fractions_secondarySigma_wDCAxyCut, "",
             "#it{p}_{T} (GeV/#it{c})", "Percentage to proton feeddown", 0, 1,
             1, 4, 1, kFALSE);
    hist_proton_fractions_secondaryLambda_wDCAxyCut->GetYaxis()->SetRangeUser(0,1);
    hist_proton_fractions_secondarySigma_wDCAxyCut->GetYaxis()->SetRangeUser(0,1);
    SigmaLambdaPercentage->Draw("same");
    Can_LambdaSigma_fractions->Print("ANplot/SigmaLambdaPercentage.pdf");
  }

  //***************************************************************************
  // Cosine pointing angle of Lambdas:
  //***************************************************************************

  TCanvas* Can_Lambda_simexpComp = new TCanvas(
      "Can_Lambda_simexpComp", "Can_Lambda_simexpComp", 0, 0, 1300, 900);
  Can_Lambda_simexpComp->Divide(4, 2);

  TCanvas* Can_Lambda_simDecomp = new TCanvas(
      "Can_Lambda_simDecomp", "Can_Lambda_simDecomp", 0, 0, 1300, 900);
  Can_Lambda_simDecomp->Divide(4, 2);

  for (Int_t ptbin = 0; ptbin < CPA_Ptbins; ptbin++) {
    Float_t sim_exp_scale =
        hist_lambda_CPA_exp[ptbin]->Integral(
            hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_low),
            hist_lambda_CPA_exp[ptbin]->FindBin(V0_CPA_norm_up)) /
        hist_lambda_CPA[ptbin]->Integral(
            hist_lambda_CPA[ptbin]->FindBin(V0_CPA_norm_low),
            hist_lambda_CPA[ptbin]->FindBin(V0_CPA_norm_up));
    hist_lambda_CPA[ptbin]->Scale(sim_exp_scale);
    Can_Lambda_simexpComp->cd(ptbin + 1);
    gPad->SetLogy();
    DrawHist(hist_lambda_CPA_exp[ptbin], hist_lambda_CPA[ptbin], "", "CPA",
             "Counts", 0, 1, 1, 4, 1, kFALSE);

    Can_Lambda_simDecomp->cd(ptbin + 1);
    gPad->SetLogy();
    DrawHist(hist_lambda_CPA_primary[ptbin], hist_lambda_CPA_secondary[ptbin],
             hist_lambda_CPA_material[ptbin], hist_lambda_CPA_bkg[ptbin], "",
             "#it{p}_{T} (GeV/#it{c})", "fraction", 0, 1, 1, 4, 1, 2, 1, 7, 1,
             kFALSE);
    hist_lambda_CPA_exp[ptbin]->DrawCopy("same");
  }

  TH1D* hist_lambda_CPA_primary_normalized_ExpEvts_weight[CPA_Ptbins];
  TH1D* hist_lambda_CPA_secondary_normalized_ExpEvts_weight[CPA_Ptbins];
  TH1D* hist_lambda_CPA_material_normalized_ExpEvts_weight[CPA_Ptbins];
  TH1D* hist_lambda_CPA_bkg_normalized_ExpEvts_weight[CPA_Ptbins];
  TH1D* hist_lambda_CPA_templatesum[CPA_Ptbins];

  TH1F* hist_lambda_fraction_primary =
      new TH1F("hist_lambda_fraction_primary",
               "Fraction of primary lambdas in the sample", V0_InvMass_Ptbins,
               V0_InvMass_threshold_low, V0_InvMass_threshold_up);
  TH1F* hist_lambda_fraction_secondary =
      new TH1F("hist_lambda_fraction_secondary",
               "Fraction of secondary lambdas in the sample", V0_InvMass_Ptbins,
               V0_InvMass_threshold_low, V0_InvMass_threshold_up);
  TH1F* hist_lambda_fraction_material =
      new TH1F("hist_lambda_fraction_material",
               "Fraction of material lambdas in the sample", V0_InvMass_Ptbins,
               V0_InvMass_threshold_low, V0_InvMass_threshold_up);
  TH1F* hist_lambda_fraction_bkg = new TH1F(
      "hist_lambda_fraction_bkg", "Fraction of fake lambdas in the sample",
      V0_InvMass_Ptbins, V0_InvMass_threshold_low, V0_InvMass_threshold_up);

  TCanvas* Can_lambda_CPA_templateFit =
      new TCanvas("Can_lambda_CPA_templateFit", "Can_lambda_CPA_templateFit", 0,
                  0, 1300, 900);
  Can_lambda_CPA_templateFit->Divide(4, 2);

  // Parameters to obtain the average weighted by the pt production spectrum of
  // protons

  Double_t fraction_Lambda_primary_ptweighted = 0.;
  Double_t fraction_Lambda_secondary_ptweighted = 0.;
  Double_t fraction_Lambda_material_ptweighted = 0.;
  Double_t fraction_Lambda_background_ptweighted = 0.;

  Double_t ptweight_Lambda_pt = 0;

  TLegend* leg_CPA = new TLegend(0.2, 0.7, 0.58, 0.9);
  leg_CPA->SetTextSize(0.055);
  // leg->SetHeader("The Legend Title","C"); // option "C" allows to center the
  // header

  std::cout << "*************************" << std::endl;
  std::cout << "*************************" << std::endl;
  std::cout << "start CPA template fit" << std::endl;
  std::cout << "*************************" << std::endl;
  std::cout << "*************************" << std::endl;
  if (doCPA_templateFit) {
    for (int ptbin = 0; ptbin < CPA_Ptbins; ptbin++) {
      std::cout << "*************************" << std::endl;
      std::cout << "*************************" << std::endl;
      std::cout << "ptBin" << ptbin << std::endl;
      std::cout << "*************************" << std::endl;
      std::cout << "*************************" << std::endl;

      TObjArray* mcDCA_array_CPA =
          new TObjArray(4);  // MC histograms are put in this array
      mcDCA_array_CPA->Add(hist_lambda_CPA_primary_normalized_ExpEvts[ptbin]);
      mcDCA_array_CPA->Add(hist_lambda_CPA_secondary_normalized_ExpEvts[ptbin]);
      mcDCA_array_CPA->Add(hist_lambda_CPA_material_normalized_ExpEvts[ptbin]);
      mcDCA_array_CPA->Add(hist_lambda_CPA_bkg_normalized_ExpEvts[ptbin]);

      TFractionFitter* fit_CPA = new TFractionFitter(
          hist_lambda_CPA_exp[ptbin], mcDCA_array_CPA);  // initialise
      if(ptbin == 1) fit_CPA->Constrain(0,0.,1.);
      if(ptbin == 2) fit_CPA->Constrain(0,0.,10.);
      if(ptbin == 3) fit_CPA->Constrain(0,0.,10.);

      Int_t status = fit_CPA->Fit();  // perform the fit

      if (status == 0) {
        Double_t par0, errpar0, par1, errpar1, par2, errpar2, par3, errpar3;

        fit_CPA->GetResult(0, par0, errpar0);
        fit_CPA->GetResult(1, par1, errpar1);
        fit_CPA->GetResult(2, par2, errpar2);
        fit_CPA->GetResult(3, par3, errpar3);

        hist_lambda_CPA_primary_normalized_ExpEvts_weight[ptbin] =
            (TH1D*)hist_lambda_CPA_primary_normalized_ExpEvts[ptbin]->Clone(
                "hist_lambda_CPA_primary_normalized_ExpEvts_weight");
        hist_lambda_CPA_primary_normalized_ExpEvts_weight[ptbin]->Scale(par0);

        hist_lambda_CPA_secondary_normalized_ExpEvts_weight[ptbin] =
            (TH1D*)hist_lambda_CPA_secondary_normalized_ExpEvts[ptbin]->Clone(
                "hist_lambda_CPA_secondary_normalized_ExpEvts_weight");
        hist_lambda_CPA_secondary_normalized_ExpEvts_weight[ptbin]->Scale(par1);

        hist_lambda_CPA_material_normalized_ExpEvts_weight[ptbin] =
            (TH1D*)hist_lambda_CPA_material_normalized_ExpEvts[ptbin]->Clone(
                "hist_lambda_CPA_material_normalized_ExpEvts_weight");
        hist_lambda_CPA_material_normalized_ExpEvts_weight[ptbin]->Scale(par2);

        hist_lambda_CPA_bkg_normalized_ExpEvts_weight[ptbin] =
            (TH1D*)hist_lambda_CPA_bkg_normalized_ExpEvts[ptbin]->Clone(
                "hist_lambda_CPA_bkg_normalized_ExpEvts_weight");
        hist_lambda_CPA_bkg_normalized_ExpEvts_weight[ptbin]->Scale(par3);

        hist_lambda_CPA_templatesum[ptbin] =
            (TH1D*)hist_lambda_CPA_primary_normalized_ExpEvts_weight[ptbin]
                ->Clone("hist_lambda_CPA_templatesum");
        hist_lambda_CPA_templatesum[ptbin]->Add(
            hist_lambda_CPA_secondary_normalized_ExpEvts_weight[ptbin]);
        hist_lambda_CPA_templatesum[ptbin]->Add(
            hist_lambda_CPA_material_normalized_ExpEvts_weight[ptbin]);
        hist_lambda_CPA_templatesum[ptbin]->Add(
            hist_lambda_CPA_bkg_normalized_ExpEvts_weight[ptbin]);

        Can_lambda_CPA_templateFit->cd(ptbin + 1);
        gPad->SetLogy();

        TString HistTitle = print_ranges(V0_InvMass_threshold_low, Delta_pt_V0,
                                         ptbin, "#it{p}_{T}", "GeV/#it{c}");
        hist_lambda_CPA_exp[ptbin]->SetTitle(HistTitle);
        hist_lambda_CPA_exp[ptbin]->SetTitleSize(1.5);
        DrawHist(
            (TH1F*)hist_lambda_CPA_exp[ptbin],
            (TH1F*)hist_lambda_CPA_primary_normalized_ExpEvts_weight[ptbin],
            (TH1F*)hist_lambda_CPA_secondary_normalized_ExpEvts_weight[ptbin],
            (TH1F*)hist_lambda_CPA_material_normalized_ExpEvts_weight[ptbin],
            (TH1F*)hist_lambda_CPA_bkg_normalized_ExpEvts_weight[ptbin], HistTitle,
            "cos #alpha", "Counts", 2, 1, 1, 2, 1, 3, 1, 4, 1, 5, 1, kTRUE);
        hist_lambda_CPA_templatesum[ptbin]->SetLineWidth(2);
        hist_lambda_CPA_templatesum[ptbin]->SetLineColor(1);
        hist_lambda_CPA_templatesum[ptbin]->DrawCopy("same");

        if (ptbin == 3) {
          leg_CPA->AddEntry(
              hist_lambda_CPA_primary_normalized_ExpEvts_weight[ptbin],
              "primary #Lambda", "f");
          leg_CPA->AddEntry(
              hist_lambda_CPA_secondary_normalized_ExpEvts_weight[ptbin],
              "secondary #Lambda", "f");
          leg_CPA->AddEntry(
              hist_lambda_CPA_material_normalized_ExpEvts_weight[ptbin],
              "material #Lambda", "f");
          leg_CPA->AddEntry(
              hist_lambda_CPA_bkg_normalized_ExpEvts_weight[ptbin],
              "fake #Lambda", "f");
        }
        leg_CPA->Draw("same");

        hist_lambda_fraction_primary->SetBinContent(ptbin + 1, par0);
        hist_lambda_fraction_primary->SetBinError(ptbin + 1, errpar0);

        hist_lambda_fraction_secondary->SetBinContent(ptbin + 1, par1);
        hist_lambda_fraction_secondary->SetBinError(ptbin + 1, errpar1);

        hist_lambda_fraction_material->SetBinContent(ptbin + 1, par2);
        hist_lambda_fraction_material->SetBinError(ptbin + 1, errpar2);

        hist_lambda_fraction_bkg->SetBinContent(ptbin + 1, par3);
        hist_lambda_fraction_bkg->SetBinError(ptbin + 1, errpar3);

        // average with production pt
        //*************************************************
        Double_t ptval = hist_lambda_fraction_primary->GetBinCenter(ptbin + 1);

        fraction_Lambda_primary_ptweighted +=
            par0 *
            hist_Lambda_pt_spectrum->GetBinContent(
                hist_Lambda_pt_spectrum->FindBin(ptval));
        fraction_Lambda_secondary_ptweighted +=
            par1 *
            hist_Lambda_pt_spectrum->GetBinContent(
                hist_Lambda_pt_spectrum->FindBin(ptval));
        fraction_Lambda_material_ptweighted +=
            par2 *
            hist_Lambda_pt_spectrum->GetBinContent(
                hist_Lambda_pt_spectrum->FindBin(ptval));
        fraction_Lambda_background_ptweighted +=
            par3 *
            hist_Lambda_pt_spectrum->GetBinContent(
                hist_Lambda_pt_spectrum->FindBin(ptval));

        ptweight_Lambda_pt += hist_Lambda_pt_spectrum->GetBinContent(
            hist_Lambda_pt_spectrum->FindBin(ptval));

        std::cout << "par0: " << par0 << std::endl;

        //*************************************************
      }
    }
    Can_lambda_CPA_templateFit->Print("ANplot/lambdaTemplateFit.pdf");

    // Can_lambda_CPA_templateFit->Print(Create_path(alinotedir,"MC_template","CPA_template",""));

    fraction_Lambda_primary_ptweighted /= ptweight_Lambda_pt;
    fraction_Lambda_secondary_ptweighted /= ptweight_Lambda_pt;
    fraction_Lambda_material_ptweighted /= ptweight_Lambda_pt;
    fraction_Lambda_background_ptweighted /= ptweight_Lambda_pt;

    const Int_t PtbinCPA_repo = 2;

    TCanvas* Can_proton_CPA_representation_fit =
        new TCanvas("Can_proton_CPA_representation_fit",
                    "Can_proton_CPA_representation_fit");
    gPad->SetLogy();
    DrawHist(
        (TH1F*)hist_lambda_CPA_exp[PtbinCPA_repo],
        (TH1F*)hist_lambda_CPA_primary_normalized_ExpEvts_weight[PtbinCPA_repo],
        (TH1F*)
            hist_lambda_CPA_secondary_normalized_ExpEvts_weight[PtbinCPA_repo],
        (TH1F*)
            hist_lambda_CPA_material_normalized_ExpEvts_weight[PtbinCPA_repo],
        (TH1F*)hist_lambda_CPA_bkg_normalized_ExpEvts_weight[PtbinCPA_repo], "",
        "cos #alpha", "Counts",  2, 1, 1, 2, 1, 3, 1, 4, 1, 5, 1, kFALSE, true);
    hist_lambda_CPA_templatesum[PtbinCPA_repo]->DrawCopy("histo same");
    hist_lambda_CPA_templatesum[PtbinCPA_repo]->SetLineWidth(2);
    hist_lambda_CPA_templatesum[PtbinCPA_repo]->SetLineColor(colors[1]);
    hist_lambda_CPA_templatesum[PtbinCPA_repo]->DrawCopy("histo same");
    print_ranges(V0_InvMass_threshold_low, Delta_pt_V0, PtbinCPA_repo,
                 "#it{p}_{T}", "GeV/#it{c}", 0.992, 5E5);
    Can_proton_CPA_representation_fit->Print("ANplot/lambdaTemplateFit_single.pdf");

  }

  //***************************************************************************
  // Purity of Lambda Sample:
  //***************************************************************************

  Double_t purity_Lambda_ptweighted = 0.;
  Double_t purity_Lambda_weights = 0.;

  TH1F* hist_lambda_purity = new TH1F(
      "hist_lambda_purity", "purity of lambda sample", V0_InvMass_Ptbins,
      V0_InvMass_threshold_low, V0_InvMass_threshold_up);
  TH1F* hist_lambda_signal_over_bkg = new TH1F(
      "hist_lambda_signal_over_bkg", "signal over background of lambda sample",
      V0_InvMass_Ptbins, V0_InvMass_threshold_low, V0_InvMass_threshold_up);
  TH1F* hist_lambda_signal = new TH1F(
      "hist_lambda_signal", "signal of lambda sample as a function of pt",
      V0_InvMass_Ptbins, V0_InvMass_threshold_low, V0_InvMass_threshold_up);

  TF1* InvMassFit_spectrum =
      new TF1("InvMassFit_spectrum", GaussPolyFit, 1., 2., 12);

  InvMassFit_spectrum->SetParameter(
      0, hist_lambda_InvMass[0]->GetBinContent(
             hist_lambda_InvMass[0]->GetMaximumBin()) *
             1.3);
  InvMassFit_spectrum->SetParLimits(
      0, 0, hist_lambda_InvMass[0]->GetBinContent(
                hist_lambda_InvMass[0]->GetMaximumBin()) *
                3.7);
  InvMassFit_spectrum->SetParameter(1, 0.5);
  InvMassFit_spectrum->SetParLimits(1, 0, 1.);
  InvMassFit_spectrum->SetParameter(2, 1.116);
  InvMassFit_spectrum->SetParLimits(2, 1.116 - 0.006, 1.116 + 0.006);
  InvMassFit_spectrum->SetParameter(3, 3);  // 3 sigma
  InvMassFit_spectrum->SetParLimits(3, 0, 6);
  InvMassFit_spectrum->SetParameter(4, 1.116);
  InvMassFit_spectrum->SetParLimits(4, 1.116 - 0.006, 1.116 + 0.006);
  InvMassFit_spectrum->SetParameter(5, 5);  // 5
  InvMassFit_spectrum->SetParLimits(5, 0, 10);
  InvMassFit_spectrum->SetParameter(6, 10);
  InvMassFit_spectrum->SetParameter(7, 10);
  InvMassFit_spectrum->SetParameter(8, 0.1);
  InvMassFit_spectrum->SetParameter(9, 0);
  InvMassFit_spectrum->SetParameter(10, 0);
  InvMassFit_spectrum->SetParameter(11, 0);

  InvMassFit_spectrum->SetLineWidth(2);
  InvMassFit_spectrum->SetLineColor(colors[2]);

  // Double_t InvMassFit_parameter[12];

  TCanvas* Can_lambda_InvMass =
      new TCanvas("Can_lambda_InvMass", "Can_lambda_InvMass", 0, 0, 1300, 900);
  Can_lambda_InvMass->Divide(4, 2);
  double totalsignal = 0.;
  for (Int_t ptbin = 0; ptbin < CPA_Ptbins; ptbin++) {
    Can_lambda_InvMass->cd(ptbin + 1);
    hist_lambda_InvMass[ptbin]->GetXaxis()->SetRangeUser(1.09, 1.14);
    TString HistTitle = print_ranges(V0_InvMass_threshold_low, Delta_pt_V0,
                                     ptbin, "#it{p}_{T}", "GeV/#it{c}");
    hist_lambda_InvMass[ptbin]->SetTitle(HistTitle);
    hist_lambda_InvMass[ptbin]->SetTitleSize(1.5);
    hist_lambda_InvMass[ptbin]->GetXaxis()->SetRangeUser(1.1, 1.13);
    DrawHist((TH1F*)hist_lambda_InvMass[ptbin], "",
             "IM_{p#pi} [GeV/#it{c}^{2}]", "Counts/(1 MeV/c^{2})", 0, 1, 1,
             kTRUE);
    InvMassFit_spectrum->SetRange(1.09, 1.14);
    hist_lambda_InvMass[ptbin]->Fit(InvMassFit_spectrum, "N", "", 1.09, 1.14);
    InvMassFit_spectrum->SetLineColor(colors[2]);
    InvMassFit_spectrum->SetLineStyle(1);
    InvMassFit_spectrum->DrawCopy("same");

    // Determine signal + background:
    Double_t integral_InvMass_sigbkg =
        InvMassFit_spectrum->Integral(V0_InvMass_lowfit, V0_InvMass_upfit);
    Double_t integral_InvMass_sigbkgErr =
        InvMassFit_spectrum->IntegralError(V0_InvMass_lowfit, V0_InvMass_upfit);
    std::cout << integral_InvMass_sigbkg << " " << integral_InvMass_sigbkgErr
              << "\n";
    InvMassFit_spectrum->SetParameter(0, 0.);  // Set the amplitude == 0, to
                                               // obtain only the background
                                               // contribution
    InvMassFit_spectrum->SetLineStyle(2);
    InvMassFit_spectrum->DrawCopy("same");

    Double_t integral_InvMass_bkg =
        InvMassFit_spectrum->Integral(V0_InvMass_lowfit, V0_InvMass_upfit);
    Double_t integral_InvMass_bkgErr =
        InvMassFit_spectrum->IntegralError(V0_InvMass_lowfit, V0_InvMass_upfit);
    Double_t integral_InvMass_sig =
        integral_InvMass_sigbkg - integral_InvMass_bkg;

    Double_t purity = integral_InvMass_sig / integral_InvMass_sigbkg;
    Double_t purityErr =
        std::sqrt(integral_InvMass_sigbkgErr * integral_InvMass_sigbkgErr *
                      integral_InvMass_bkg * integral_InvMass_bkg /
                      std::pow(integral_InvMass_sigbkg, 4) +
                  integral_InvMass_bkgErr * integral_InvMass_bkgErr /
                      (integral_InvMass_sigbkg * integral_InvMass_sigbkg));

    Double_t sigma_1 = InvMassFit_spectrum->GetParameter(3);
    Double_t sigma_2 = InvMassFit_spectrum->GetParameter(5);

    totalsignal += integral_InvMass_sig;
    std::cout << "chi2: " << InvMassFit_spectrum->GetChisquare() << std::endl;
    std::cout << "NDF: " << InvMassFit_spectrum->GetNDF() << std::endl;
    std::cout << "signal: " << integral_InvMass_sig << std::endl;
    std::cout << "sigma_1: " << sigma_1 << std::endl;
    std::cout << "sigma_2: " << sigma_2 << std::endl;
    std::cout << "........................." << std::endl;

    TString HistName = "hist_lambda_InvMass_signal_";
    HistName += ptbin;
    TH1F* hist_lambda_InvMass_signal =
        (TH1F*)hist_lambda_InvMass[ptbin]->Clone(HistName.Data());
    hist_lambda_InvMass_signal->Add(InvMassFit_spectrum,
                                    -1);  // subract background contribution

    hist_lambda_purity->SetBinContent(ptbin + 1, purity);
    hist_lambda_purity->SetBinError(ptbin + 1, purityErr);

    hist_lambda_signal_over_bkg->SetBinContent(
        ptbin + 1, integral_InvMass_sig / integral_InvMass_bkg);
    hist_lambda_signal_over_bkg->SetBinError(ptbin + 1, 0);

    hist_lambda_signal->SetBinContent(ptbin + 1, integral_InvMass_sig);

    // average parameters with production p.d.f
    //***************************************************
    Double_t ptval = hist_lambda_signal_over_bkg->GetBinCenter(ptbin + 1);

    purity_Lambda_ptweighted +=
        purity *
        hist_lambda_signal->GetBinContent(hist_lambda_signal->FindBin(ptval));
    purity_Lambda_weights +=
        hist_lambda_signal->GetBinContent(hist_lambda_signal->FindBin(ptval));

    printval("P_{#Lambda}", "", purity, 6, 90000);
    //***************************************************
  }

  // Can_lambda_InvMass->Print(Create_path(alinotedir,"MC_template","Lambda_InvMass",""));

  purity_Lambda_ptweighted /= purity_Lambda_weights;

  // hist_lambda_signal_over_bkg->GetYaxis()->SetRangeUser(0,45);
  // DrawHist(hist_lambda_signal_over_bkg,"#it{p}_{T}
  // (GeV/#it{c})","S/B",0,1,1);

  TLegend* leg_Lambda_summary = new TLegend(0.3, 0.4, 0.78, 0.6);
  leg_Lambda_summary->SetTextSize(0.055);

  TCanvas* Can_lambda_fractions = new TCanvas(
      "Can_lambda_fractions", "Can_lambda_fractions", 0, 0, 1300, 900);
  Can_lambda_fractions->Divide(2, 2);
  Can_lambda_fractions->cd(1);
  SetStyleHisto(hist_Lambda_pt_spectrum, 1,1);
  SetStyleHisto(hist_lambda_purity, 1,1);
  gPad->SetLogy();
  DrawHist(hist_Lambda_pt_spectrum, "", "#it{p}_{T} (GeV/#it{c})",
           "dN/d#it{p}_{T} 1/(GeV/#it{c})", 0, 1, 1, kFALSE);
  Can_lambda_fractions->cd(3);
  hist_lambda_purity->GetYaxis()->SetRangeUser(0.8, 1.1);
  hist_lambda_purity->Write();
  DrawHist(hist_lambda_purity, "", "#it{p}_{T} (GeV/#it{c})",
           "Purity(#Lambda) = S/S+B", 0, 1, 1, kFALSE);
  horline->DrawLine(V0_InvMass_threshold_low, purity_Lambda_ptweighted,
                    V0_InvMass_threshold_up, purity_Lambda_ptweighted);
  horline->SetLineColor(colors[0]);
  Can_lambda_fractions->cd(2);
  hist_lambda_fraction_primary->GetYaxis()->SetRangeUser(0,1);
  DrawHist(hist_lambda_fraction_primary, hist_lambda_fraction_secondary,
           hist_lambda_fraction_material, hist_lambda_fraction_bkg, "",
           "#it{p}_{T} (GeV/#it{c})", "Fraction", 0, 1, 1, 2, 1, 3, 1, 4, 1,
           kFALSE, true);
  horline->DrawLine(V0_InvMass_threshold_low,
                    fraction_Lambda_primary_ptweighted, V0_InvMass_threshold_up,
                    fraction_Lambda_primary_ptweighted);
  horline->SetLineColor(colors[0]);
  horline->DrawLine(
      V0_InvMass_threshold_low, fraction_Lambda_secondary_ptweighted,
      V0_InvMass_threshold_up, fraction_Lambda_secondary_ptweighted);
  horline->SetLineColor(colors[1]);
  horline->DrawLine(
      V0_InvMass_threshold_low, fraction_Lambda_background_ptweighted,
      V0_InvMass_threshold_up, fraction_Lambda_background_ptweighted);
  horline->SetLineColor(colors[0]);

  leg_Lambda_summary->AddEntry(hist_lambda_fraction_primary, "primary #Lambda",
                               "p");
  leg_Lambda_summary->AddEntry(hist_lambda_fraction_secondary,
                               "secondary #Lambda", "p");
  leg_Lambda_summary->AddEntry(hist_lambda_fraction_material,
                               "material #Lambda", "p");
  leg_Lambda_summary->AddEntry(hist_lambda_fraction_bkg, "fake #Lambda", "p");
  leg_Lambda_summary->Draw();
  Can_lambda_fractions->Print("ANplot/lambdaFraction.pdf");

  hist_lambda_fraction_primary->Write();
  hist_lambda_fraction_secondary->Write();
  hist_lambda_fraction_material->Write();
  hist_lambda_fraction_bkg->Write();

  TCanvas* Can_lambda_signal =
      new TCanvas("Can_lambda_signal", "Can_lambda_signal");
  Can_lambda_signal->cd();
  hist_lambda_signal->Draw();
  hist_lambda_signal->SetTitle("; #it{p}_{T} (GeV/#it{c}); Counts");

  std::cout << "Parameters" << std::endl;
  std::cout << "---------- For proton-proton CF ----------" << std::endl;
  std::cout << "fraction of primary protons: "
            << fraction_proton_primary_ptweighted << std::endl;
  std::cout << "fraction of lambdas feeding to protons: "
            << fraction_Lambda_feeddown_ptweighted << std::endl;
  std::cout << "fraction of protons produced at detector material: "
            << fraction_proton_material_ptweighted << std::endl;
  std::cout << "total fraction: "
            << fraction_proton_primary_ptweighted +
                   fraction_Lambda_feeddown_ptweighted +
                   fraction_proton_material_ptweighted
            << std::endl;
  std::cout << "purity of proton sample: " << purity_proton_ptweighted
            << std::endl;
  std::cout << "---------- For Lambda-proton CF ----------" << std::endl;
  std::cout << "total number of lambdas: " << totalsignal << std::endl;
  std::cout << "fraction of primary Lambdas: "
            << fraction_Lambda_primary_ptweighted << std::endl;
  std::cout << "fraction of secondary Lambdas: "
            << fraction_Lambda_secondary_ptweighted << std::endl;
  std::cout << "fraction of material Lambdas: "
            << fraction_Lambda_material_ptweighted << std::endl;
  std::cout << "fraction of background Lambdas: "
            << fraction_Lambda_background_ptweighted << std::endl;
  std::cout << "total fraction: "
            << fraction_Lambda_primary_ptweighted +
                   fraction_Lambda_secondary_ptweighted +
                   fraction_Lambda_material_ptweighted +
                   fraction_Lambda_background_ptweighted
            << std::endl;
  std::cout << "purity of Lambda sample: " << purity_Lambda_ptweighted
            << "\n\n\n\n";

  /// TODO: Fraction Lambda/Sigma (0.7, 0.3)
  ///
  const float fractionLambda = 0.7;
  const float fractionSigma = 0.3;

  const bool Compact = true;

  // for the proton:
  // 0 = primary
  // 1 = from Lambda
  // 2 = other feeddown (flat)
  // 3 = missidentified
  const unsigned NumChannels_p = 4;
  const double ProtonPurity = purity_proton_ptweighted;
  const double Purities_p[NumChannels_p] = {ProtonPurity, ProtonPurity,
                                            ProtonPurity, 1. - ProtonPurity};
  const double Fraction_p[NumChannels_p] = {
      fraction_proton_primary_ptweighted /
          (fraction_proton_primary_ptweighted +
           fraction_Lambda_feeddown_ptweighted),
      fraction_Lambda_feeddown_ptweighted /
          (fraction_proton_primary_ptweighted +
           fraction_Lambda_feeddown_ptweighted) *
          fractionLambda,
      fraction_Lambda_feeddown_ptweighted /
          (fraction_proton_primary_ptweighted +
           fraction_Lambda_feeddown_ptweighted) *
          fractionSigma,
      1};

  // for the Lambda:
  // 0 = primary
  // 1 = from Sigma0
  // 2 = from Xim
  // 3 = from Xi0
  // 4 = missidentified
  const unsigned NumChannels_L = 5;
  const double LambdaPurity = purity_Lambda_ptweighted;
  const double Purities_L[NumChannels_L] = {LambdaPurity, LambdaPurity,
                                            LambdaPurity, LambdaPurity,
                                            1. - LambdaPurity};
  const double Fraction_L[NumChannels_L] = {
      fraction_Lambda_primary_ptweighted /
          (fraction_Lambda_primary_ptweighted +
           fraction_Lambda_secondary_ptweighted) *
          3. / 4.,
      fraction_Lambda_primary_ptweighted /
          (fraction_Lambda_primary_ptweighted +
           fraction_Lambda_secondary_ptweighted) *
          1. / 4.,
      fraction_Lambda_secondary_ptweighted /
          (fraction_Lambda_primary_ptweighted +
           fraction_Lambda_secondary_ptweighted) *
          0.5,
      fraction_Lambda_secondary_ptweighted /
          (fraction_Lambda_primary_ptweighted +
           fraction_Lambda_secondary_ptweighted) *
          0.5,
      1};

  // for the Xi:
  // 0 = primary
  // 1 = from Xi-(1530)
  // 2 = from Xi0(1530)
  // 3 = from Omega
  // 4 = missidentified
  const unsigned NumChannels_Xim = 5;
  const double XimPurity = 0.9;
  // ratio Xi-(1530) to Xi0 (m=minus)
  const double Xim1530_to_Xim = 0.32 * (1. / 3.);
  // ratio Xi0(1530) to Xi0 (n=neutral)
  const double Xin1530_to_Xim = 0.32 * (2. / 3.);
  const double Omegam_to_Xim = 0.1;   // ratio Î©/Îž
  const double OmegamXim_BR = 0.086;  // branching ratio Î©->Îž
  const double Purities_Xim[NumChannels_Xim] = {XimPurity, XimPurity, XimPurity,
                                                XimPurity, 1. - XimPurity};
  const double Fraction_Xim[NumChannels_Xim] = {
      1. - Xim1530_to_Xim - Xin1530_to_Xim - Omegam_to_Xim * OmegamXim_BR,
      Xim1530_to_Xim, Xin1530_to_Xim, Omegam_to_Xim * OmegamXim_BR, 1};


  std::cout << "pp purity : " << purity_proton_ptweighted << "\n";
  for(int i=0; i<NumChannels_p; ++i) {
    std::cout << "pp fraction " << i << " " << Fraction_p[i] << "\n";
  }

  std::cout << "pL purity : " << purity_Lambda_ptweighted << "\n";
  for(int i=0; i<NumChannels_L; ++i) {
    std::cout << "pL fraction " << i << " " << Fraction_L[i] << "\n";
  }

  std::cout << "\n\n\n";

  printf("Legend:\n");
  printf(" A : primary particle\n");
  printf(" A_B : feed-down from B to A\n");
  printf(" A' : missidentified particle\n");
  printf("\n");

  TString* ChannelName_p = new TString[NumChannels_p];
  ChannelName_p[0] = "p";
  ChannelName_p[1] = "p_Lambda";
  ChannelName_p[2] = "p_Sigma+";
  ChannelName_p[3] = "p'";

  TString* ChannelName_L = new TString[NumChannels_L];
  ChannelName_L[0] = "Lambda";
  ChannelName_L[1] = "Lambda_Sigma0";
  ChannelName_L[2] = "Lambda_Xi-";
  ChannelName_L[3] = "Lambda_Xi0";
  ChannelName_L[4] = "Lambda'";

  TString* ChannelName_Xim = new TString[NumChannels_Xim];
  ChannelName_Xim[0] = "Xi-";
  ChannelName_Xim[1] = "Xi-_Xi-(1530)";
  ChannelName_Xim[2] = "Xi-_Xi0(1530)";
  ChannelName_Xim[3] = "Xi-_Omega-";
  ChannelName_Xim[4] = "Xi-'";

  double LambdaPar;
  double TotProb = 0;
  printf("p_p:\n");
  for (unsigned uCh_p = 0; uCh_p < NumChannels_p; uCh_p++) {
    for (unsigned uCh_p2 = 0; uCh_p2 < NumChannels_p; uCh_p2++) {
      if (Compact && uCh_p > uCh_p2) continue;
      LambdaPar = Purities_p[uCh_p] * Fraction_p[uCh_p] * Purities_p[uCh_p2] *
                  Fraction_p[uCh_p2];
      if (Compact && uCh_p != uCh_p2) LambdaPar *= 2;
      TotProb += LambdaPar;
      printf(" lambda(%s,%s) = %.2f%\n", ChannelName_p[uCh_p].Data(),
             ChannelName_p[uCh_p2].Data(), LambdaPar * 100);
    }
  }
  // printf(" TotProb = %.3f\n", TotProb*100);
  printf("\np_Lambda:\n");
  TotProb = 0;
  for (unsigned uCh_p = 0; uCh_p < NumChannels_p; uCh_p++) {
    for (unsigned uCh_L = 0; uCh_L < NumChannels_L; uCh_L++) {
      LambdaPar = Purities_p[uCh_p] * Fraction_p[uCh_p] * Purities_L[uCh_L] *
                  Fraction_L[uCh_L];
      TotProb += LambdaPar;
      printf(" lambda(%s,%s) = %.2f%\n", ChannelName_p[uCh_p].Data(),
             ChannelName_L[uCh_L].Data(), LambdaPar * 100);
    }
  }

  printf("\nLambda_Lambda:\n");
  TotProb = 0;
  for (unsigned uCh_L = 0; uCh_L < NumChannels_L; uCh_L++) {
    for (unsigned uCh_L2 = 0; uCh_L2 < NumChannels_L; uCh_L2++) {
      if (Compact && uCh_L > uCh_L2) continue;
      LambdaPar = Purities_L[uCh_L] * Fraction_L[uCh_L] * Purities_L[uCh_L2] *
                  Fraction_L[uCh_L2];
      if (Compact && uCh_L != uCh_L2) LambdaPar *= 2;
      TotProb += LambdaPar;
      printf(" lambda(%s,%s) = %.2f%\n", ChannelName_L[uCh_L].Data(),
             ChannelName_L[uCh_L2].Data(), LambdaPar * 100);
    }
  }
  // printf(" TotProb = %.3f\n", TotProb*100);

  printf("\np_Xi-:\n");
  TotProb = 0;
  for (unsigned uCh_p = 0; uCh_p < NumChannels_p; uCh_p++) {
    for (unsigned uCh_Xim = 0; uCh_Xim < NumChannels_Xim; uCh_Xim++) {
      LambdaPar = Purities_p[uCh_p] * Fraction_p[uCh_p] *
                  Purities_Xim[uCh_Xim] * Fraction_Xim[uCh_Xim];
      TotProb += LambdaPar;
      printf(" lambda(%s,%s) = %.2f%\n", ChannelName_p[uCh_p].Data(),
             ChannelName_Xim[uCh_Xim].Data(), LambdaPar * 100);
    }
  }

  return 0;
}
