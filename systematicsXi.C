# include "TH1F.h"
# include <iostream>
# include "TAxis.h"
# include "TCanvas.h"
# include "TLine.h"
# include "TLatex.h"
# include "TProfile.h"
# include "TProfile2D.h"
# include "TGraph.h"
# include "TF1.h"
# include "TList.h"
# include "TMath.h"
# include "TFile.h"
# include "TStyle.h"
# include <stdio.h>
# include <unistd.h>
# include "TSystem.h"
# include "TGraphErrors.h"


std::vector<int> fillColors = {kGray + 1,  kRed - 10,    kBlue - 9,
                               kGreen - 8, kMagenta - 9, kOrange - 9,
                               kCyan - 8,  kYellow - 7};
std::vector<int> colors = {kBlack,       kRed + 1,    kBlue + 2, kGreen + 3,
                           kMagenta + 1, kOrange - 1, kCyan + 2, kYellow + 2};
std::vector<int> markers = {
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

//DrawHist(hCFPPAPAP_Barlow[protonVariation],"k* (GeV/#it{c})","n#sigma",0,1,1.);
void DrawHist(TH1F* hist1,TString XTitle, TString YTitle,Int_t option,Int_t MarkerColor,Int_t MarkerSize,Bool_t ShowTitle=false)
{
  hist1->SetStats(0);
  if(!ShowTitle)hist1->SetTitle(0);
  hist1->GetXaxis()->SetTitle(XTitle.Data());
  hist1->GetYaxis()->SetTitle(YTitle.Data());
  hist1->SetMarkerColor(colors[MarkerColor]);
  hist1->SetLineColor(colors[MarkerColor]);
  hist1->SetMarkerStyle(20);
  if(option == 0) hist1->Draw("PE");
  else if(option == 1) hist1->Draw("histo");
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
              Int_t MarkerColor4, Int_t MarkerSize4, Bool_t showTitle) {
  hist1->SetStats(0);
  if (!showTitle) {
    hist1->SetTitle(0);
  } else {
    hist1->SetTitle(HistTitle.Data());
  }
  hist1->GetXaxis()->SetTitle(XTitle.Data());
  hist1->GetYaxis()->SetTitle(YTitle.Data());

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
    hist1->SetMarkerStyle(markers[0]);

    hist2->SetFillStyle(3004);
    hist2->SetFillColor(fillColors[MarkerColor2 - 1]);

    hist3->SetFillStyle(3005);
    hist3->SetFillColor(fillColors[MarkerColor3 - 1]);

    hist4->SetFillStyle(3013);
    hist4->SetFillColor(fillColors[MarkerColor4 - 1]);

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
              Int_t MarkerSize5, Bool_t showTitle) {
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
    hist5->SetMarkerColor(colors[MarkerColor5 - 1]);
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
    hist2->SetFillColor(fillColors[MarkerColor2 - 1]);

    hist3->SetFillStyle(3005);
    hist3->SetFillColor(fillColors[MarkerColor3 - 1]);

    hist4->SetFillStyle(3013);
    hist4->SetFillColor(fillColors[MarkerColor4 - 1]);

    hist5->SetFillStyle(3006);
    hist5->SetFillColor(fillColors[MarkerColor5 - 1]);

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
void HistoStyleMany(TH1F *histo) {
  histo->GetXaxis()->SetLabelSize(0.075);
  histo->GetXaxis()->SetTitleSize(0.08);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.075);
  histo->GetYaxis()->SetTitleSize(0.08);
  histo->GetYaxis()->SetLabelFont(42);
  histo->SetTitleFont(42);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void StyleGraph(TGraph *gr) {
  gr->GetYaxis()->SetTitleSize(0.07);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetTitleSize(0.07);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->SetMarkerStyle(20);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* add_CF(TH1F* hist_CF1, TH1F* hist_CF2,TString HistName)
{
  //Calculate CFs with error weighting
  TH1F* hist_CF_sum = (TH1F*)hist_CF1->Clone(HistName.Data());

  Int_t NBins = hist_CF_sum->GetNbinsX();
  Int_t NBins2 = hist_CF1->GetNbinsX();

  if(NBins != NBins2) std::cout << "Number of bins is different" << std::endl;

  for(Int_t i=0;i<NBins;i++)
  {
    Float_t CF1_val = hist_CF1->GetBinContent(i+1);
    Float_t CF1_err = hist_CF1->GetBinError(i+1);
    Float_t CF2_val = hist_CF2->GetBinContent(i+1);
    Float_t CF2_err = hist_CF2->GetBinError(i+1);

    //average for bin i:
    if(CF1_val != 0. && CF2_val != 0.)
    {
      Float_t CF1_err_weight = 1./TMath::Power(CF1_err,2.);
      Float_t CF2_err_weight = 1./TMath::Power(CF2_err,2.);

      Float_t CF_sum_average = (CF1_err_weight*CF1_val + CF2_err_weight*CF2_val)/(CF1_err_weight+CF2_err_weight);
      Float_t CF_sum_err = 1./TMath::Sqrt(CF1_err_weight+CF2_err_weight);

      hist_CF_sum->SetBinContent(i+1,CF_sum_average);
      hist_CF_sum->SetBinError(i+1,CF_sum_err);
    }
    else if(CF1_val == 0. && CF2_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF2_val);
      hist_CF_sum->SetBinError(i+1,CF2_err);
    }
    else if(CF2_val == 0. && CF1_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF1_val);
      hist_CF_sum->SetBinError(i+1,CF1_err);
    }

  }
  
  return hist_CF_sum;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t Calculate_CFnormalization(TH1F* histRE_relK, TH1F* histME_relK,Double_t normleft,Double_t normright)
{
  Double_t norm_relK = histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) / histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright));
  return norm_relK;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* Calculate_CF(TH1F* histRE_relK,TH1F* histME_relK, TString CFname,Double_t normleft,Double_t normright)
{
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  Double_t norm_relK = Calculate_CFnormalization(histRE_relK,histME_relK,normleft,normright);

  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  Hist_CF->Divide(histRE_relK,histME_relK,1,norm_relK);
  
  return Hist_CF;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int systematicsXi(const char* file,const char *output="systematics_output")
{
  gStyle->SetTitleFontSize(0.1);

  const Double_t fitRangeRightPPAPAP = 0.5;
  const Double_t fitRangeRightPLAPAL = 0.5;
  const Double_t fitRangeRightLLALAL = 0.5;
  
  const int rebinPP = 10;
  const int rebinLP = 10;
  const int rebinLL = 10;
  const float thresholdBarlow = 0.;
  const char *fit = "pol2";
  bool useFit = false;

  TH1F *hSystPP = new TH1F("hSystPP", "", 750, 0, 3);
  TH1F *hSystPL = new TH1F("hSystPL", "", 150, 0, 3);
  TH1F *hSystLL = new TH1F("hSystLL", "", 150, 0, 3);

  TH1F *hSystErrorPP = new TH1F("hSystErrorPP", "", 750, 0, 3);
  TH1F *hSystErrorPL = new TH1F("hSystErrorPL", "", 150, 0, 3);
  TH1F *hSystErrorLL = new TH1F("hSystErrorLL", "", 150, 0, 3);


  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);
  
  TFile* _file0=TFile::Open(file,"READ");
  if(!_file0) return 0;

  const Int_t NTasks = 22;//Number of different total variations (including default)
  const Int_t NTasksProtonVariations = 9; //Number of different proton variations
  const Int_t NTasksV0Variations = 10 + NTasksProtonVariations; //Number of different V0 variations (but proton selection does also affect p-V0 correlation function
  
  TList *listTP[NTasks] = {0};
  
  TString outputfile = "";
  TString TPDir = "";
  Int_t taskNumberAddon[NTasks]={0,1,2,3,4,5,6,7,8,9,19,20,21,22,23,24,25,26,27,28,29,30};
  for(int taskNumber = 0;taskNumber < NTasks; taskNumber++)
  {
    if(taskNumber == 0) outputfile += "/PWGCF_PLFemto_0";
    //proton vars
    else if(taskNumber == 1) outputfile = "/PWGCF_PLFemto_1";
    else if(taskNumber == 2) outputfile = "/PWGCF_PLFemto_2";
    else if(taskNumber == 3) outputfile = "/PWGCF_PLFemto_3";
    else if(taskNumber == 4) outputfile = "/PWGCF_PLFemto_4";
    else if(taskNumber == 5) outputfile = "/PWGCF_PLFemto_5";
    else if(taskNumber == 6) outputfile = "/PWGCF_PLFemto_6";
    else if(taskNumber == 7) outputfile = "/PWGCF_PLFemto_7";
    else if(taskNumber == 8) outputfile = "/PWGCF_PLFemto_8";
    else if(taskNumber == 9) outputfile = "/PWGCF_PLFemto_9";

    //Xi vars
    else if(taskNumber == 10) outputfile = "/PWGCF_PLFemto_19";
    else if(taskNumber == 11) outputfile = "/PWGCF_PLFemto_20";
    else if(taskNumber == 12) outputfile = "/PWGCF_PLFemto_21";
    else if(taskNumber == 13) outputfile = "/PWGCF_PLFemto_22";
    else if(taskNumber == 14) outputfile = "/PWGCF_PLFemto_23";
    else if(taskNumber == 15) outputfile = "/PWGCF_PLFemto_24";
    else if(taskNumber == 16) outputfile = "/PWGCF_PLFemto_25";
    else if(taskNumber == 17) outputfile = "/PWGCF_PLFemto_26";
    else if(taskNumber == 18) outputfile = "/PWGCF_PLFemto_27";
    else if(taskNumber == 19) outputfile = "/PWGCF_PLFemto_28";
    else if(taskNumber == 20) outputfile = "/PWGCF_PLFemto_29";
    else if(taskNumber == 21) outputfile = "/PWGCF_PLFemto_30";

    outputfile += "/TPdir_";
    outputfile += taskNumberAddon[taskNumber];

    _file0->GetObject(outputfile.Data(),listTP[taskNumber]);
    if(!listTP[taskNumber]) {std::cout << "list missing: " << outputfile.Data() <<  std::endl; return 0;}
  }
  std::cout << "lists loaded" << std::endl;
  TString HistName = "";
  Double_t normLeft = 0.2;
  Double_t normRight = 0.4;

  // ORIGINAL CF

  TH1F *hREPPorig = (TH1F*)listTP[0]->FindObject("fProtonProtonRelK");
  TH1F *hMEPPorig = (TH1F*)listTP[0]->FindObject("fProtonProtonRelKME");
  HistName = "hCFPP";
  TH1F *hCFPPorig = Calculate_CF(hREPPorig,hMEPPorig,HistName.Data(),normLeft,normRight);

  TH1F *hREPLorig = (TH1F*)listTP[0]->FindObject("fProtonXiRelK");
  TH1F *hMEPLorig = (TH1F*)listTP[0]->FindObject("fProtonXiRelKME");
  HistName = "hCFPXi";
  TH1F *hCFPLorig = Calculate_CF(hREPLorig,hMEPLorig,HistName.Data(),normLeft,normRight);

  TH1F *hRELLorig = (TH1F*)listTP[0]->FindObject("fLambdaLambdaRelK");
  TH1F *hMELLorig = (TH1F*)listTP[0]->FindObject("fLambdaLambdaRelKME");
  HistName = "hCFLL";
  TH1F *hCFLLorig = Calculate_CF(hRELLorig,hMELLorig,HistName.Data(),normLeft,normRight);

  //Second AntiBaryon-AntiBaryon Pairs:
  TH1F *hREAPAPorig = (TH1F*)listTP[0]->FindObject("fAntiProtonAntiProtonRelK");
  TH1F *hMEAPAPorig = (TH1F*)listTP[0]->FindObject("fAntiProtonAntiProtonRelKME");
  HistName = "hCFAPAP";
  TH1F *hCFAPAPorig = Calculate_CF(hREAPAPorig,hMEAPAPorig,HistName.Data(),normLeft,normRight);

  TH1F *hREAPALorig = (TH1F*)listTP[0]->FindObject("fAntiProtonAntiXiRelK");
  TH1F *hMEAPALorig = (TH1F*)listTP[0]->FindObject("fAntiProtonAntiXiRelKME");
  HistName = "hCFAPAXi";
  TH1F *hCFAPALorig = Calculate_CF(hREAPALorig,hMEAPALorig,HistName.Data(),normLeft,normRight);

  TH1F *hREALALorig = (TH1F*)listTP[0]->FindObject("fAntiLambdaAntiLambdaRelK");
  TH1F *hMEALALorig = (TH1F*)listTP[0]->FindObject("fAntiLambdaAntiLambdaRelKME");
  HistName = "hCFALAL";
  TH1F *hCFALALorig = Calculate_CF(hREALALorig,hMEALALorig,HistName.Data(),normLeft,normRight);

  //Last but not least add them:
  HistName = "hCFPPAPAP";
  TH1F *hCFPPAPAPorig = add_CF(hCFPPorig,hCFAPAPorig,HistName.Data());

  HistName = "hCFPXiAPAXi";
  TH1F *hCFPLAPALorig = add_CF(hCFPLorig,hCFAPALorig,HistName.Data());

  HistName = "hCFLLALAL";
  TH1F *hCFLLALALorig = add_CF(hCFLLorig,hCFALALorig,HistName.Data());


  //Lets define the histos:
  //Baryon-Baryon pairs:
  TH1F* hREPP[NTasks];
  TH1F* hMEPP[NTasks];
  TH1F* hCFPP[NTasks];

  TH1F* hREPL[NTasks];
  TH1F* hMEPL[NTasks];
  TH1F* hCFPL[NTasks];

  TH1F* hRELL[NTasks];
  TH1F* hMELL[NTasks];
  TH1F *hCFLL[NTasks];

  //AntiBaryon-AntiBaryon pairs:
  TH1F* hREAPAP[NTasks];
  TH1F* hMEAPAP[NTasks];
  TH1F* hCFAPAP[NTasks];

  TH1F* hREAPAL[NTasks];
  TH1F* hMEAPAL[NTasks];
  TH1F* hCFAPAL[NTasks];

  TH1F* hREALAL[NTasks];
  TH1F* hMEALAL[NTasks];
  TH1F *hCFALAL[NTasks];

  //Sum of Baryon-Baryon and AntiBaryon-AntiBaryon pairs:
  TH1F* hCFPPAPAP[NTasks];
  TH1F* hCFPLAPAL[NTasks];
  TH1F* hCFLLALAL[NTasks];

  //Lets iterate a second time over the lists to get the histos:
  for(int taskNumber = 0;taskNumber < NTasks; taskNumber++)
  {
    //First Baryon-Baryon Pairs:
    hREPP[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fProtonProtonRelK");
    hMEPP[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fProtonProtonRelKME");
    hREPP[taskNumber]->Rebin(rebinPP);
    hMEPP[taskNumber]->Rebin(rebinPP);
    HistName = "hCFPP_task";
    HistName += taskNumber;
    hCFPP[taskNumber] = Calculate_CF(hREPP[taskNumber],hMEPP[taskNumber],HistName.Data(),normLeft,normRight);
    if(!hCFPP[taskNumber]) std::cout << "PP with task number " << taskNumber << " does not exist" << std::endl;

    hREPL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fProtonXiRelK");
    hMEPL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fProtonXiRelKME");
    hREPL[taskNumber]->Rebin(rebinLP);
    hMEPL[taskNumber]->Rebin(rebinLP);
    HistName = "hCFPL_task";
    HistName += taskNumber;
    hCFPL[taskNumber] = Calculate_CF(hREPL[taskNumber],hMEPL[taskNumber],HistName.Data(),normLeft,normRight);
    if(!hCFPL[taskNumber]) std::cout << "PL with task number " << taskNumber << " does not exist" << std::endl;

    hRELL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fLambdaLambdaRelK");
    hMELL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fLambdaLambdaRelKME");
    hRELL[taskNumber]->Rebin(rebinLL);
    hMELL[taskNumber]->Rebin(rebinLL);
    HistName = "hCFLL_task";
    HistName += taskNumber;
    hCFLL[taskNumber] = Calculate_CF(hRELL[taskNumber],hMELL[taskNumber],HistName.Data(),normLeft,normRight);
    if(!hCFLL[taskNumber]) std::cout << "LL with task number " << taskNumber << " does not exist" << std::endl;

    //Second AntiBaryon-AntiBaryon Pairs:
    hREAPAP[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fAntiProtonAntiProtonRelK");
    hMEAPAP[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fAntiProtonAntiProtonRelKME");
    hREAPAP[taskNumber]->Rebin(rebinPP);
    hMEAPAP[taskNumber]->Rebin(rebinPP);
    HistName = "hCFAPAP_task";
    HistName += taskNumber;
    hCFAPAP[taskNumber] = Calculate_CF(hREAPAP[taskNumber],hMEAPAP[taskNumber],HistName.Data(),normLeft,normRight);
    if(!hCFAPAP[taskNumber]) std::cout << "APAP with task number " << taskNumber << " does not exist" << std::endl;

    hREAPAL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fAntiProtonAntiXiRelK");
    hMEAPAL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fAntiProtonAntiXiRelKME");
    hREAPAL[taskNumber]->Rebin(rebinLP);
    hMEAPAL[taskNumber]->Rebin(rebinLP);
    HistName = "hCFAPAL_task";
    HistName += taskNumber;
    hCFAPAL[taskNumber] = Calculate_CF(hREAPAL[taskNumber],hMEAPAL[taskNumber],HistName.Data(),normLeft,normRight);
    if(!hCFAPAL[taskNumber]) std::cout << "APAL with task number " << taskNumber << " does not exist" << std::endl;

    hREALAL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fAntiLambdaAntiLambdaRelK");
    hMEALAL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fAntiLambdaAntiLambdaRelKME");
    hREALAL[taskNumber]->Rebin(rebinLL);
    hMEALAL[taskNumber]->Rebin(rebinLL);
    HistName = "hCFALAL_task";
    HistName += taskNumber;
    hCFALAL[taskNumber] = Calculate_CF(hREALAL[taskNumber],hMEALAL[taskNumber],HistName.Data(),normLeft,normRight);
    if(!hCFALAL[taskNumber]) std::cout << "ALAL with task number " << taskNumber << " does not exist" << std::endl;

    //Last but not least add them:
    HistName = "hCFPPAPAP_task";
    HistName += taskNumber;
    hCFPPAPAP[taskNumber] = add_CF(hCFPP[taskNumber],hCFAPAP[taskNumber],HistName.Data());
    if(!hCFPPAPAP[taskNumber]) std::cout << "Sum of PP and APAP with task number " << taskNumber << "does not exist" << std::endl;

    HistName = "hCFPLAPAL_task";
    HistName += taskNumber;
    hCFPLAPAL[taskNumber] = add_CF(hCFPL[taskNumber],hCFAPAL[taskNumber],HistName.Data());
    if(!hCFPLAPAL[taskNumber]) std::cout << "Sum of PL and APAL with task number " << taskNumber << "does not exist" << std::endl;

    HistName = "hCFLLALAL_task";
    HistName += taskNumber;
    hCFLLALAL[taskNumber] = add_CF(hCFLL[taskNumber],hCFALAL[taskNumber],HistName.Data());
    if(!hCFLLALAL[taskNumber]) std::cout << "Sum of LL and ALAL with task number " << taskNumber << "does not exist" << std::endl;
  }

  //Also the sum, in which we are finally interested
  TCanvas *CanCFsummarySum = new TCanvas("CanCFsummarySum","CanCFsummarySum",0,0,1500,500);
  CanCFsummarySum->Divide(3,1);

  for(int taskNumber = 0;taskNumber < NTasks; taskNumber++)
  {
    if(taskNumber < NTasksProtonVariations) {
      CanCFsummarySum->cd(1);
      hCFPPAPAP[taskNumber]->GetXaxis()->SetRangeUser(0,0.2);
      hCFPPAPAP[taskNumber]->SetLineColor(taskNumber+1);
      hCFPPAPAP[taskNumber]->SetMarkerColor(taskNumber+1);
      if(taskNumber == 0) { hCFPPAPAP[taskNumber]->SetMarkerStyle(23);  hCFPPAPAP[taskNumber]->SetMarkerSize(1.3); hCFPPAPAP[taskNumber]->SetTitle("pp #oplus #bar{p}#bar{p}"); hCFPPAPAP[taskNumber]->Draw(); }
      else hCFPPAPAP[taskNumber]->Draw("same");
    }

    CanCFsummarySum->cd(2);
    hCFPLAPAL[taskNumber]->GetXaxis()->SetRangeUser(0,0.5);
    hCFPLAPAL[taskNumber]->SetLineColor(taskNumber+1);
    hCFPLAPAL[taskNumber]->SetMarkerColor(taskNumber+1);
    if(taskNumber == 0) {hCFPLAPAL[taskNumber]->SetMarkerStyle(23); hCFPLAPAL[taskNumber]->SetMarkerSize(1.3);hCFPLAPAL[taskNumber]->SetTitle("p#Lambda #oplus #bar{p}#bar{#Lambda}"); hCFPLAPAL[taskNumber]->Draw();}
    else hCFPLAPAL[taskNumber]->Draw("same");

    CanCFsummarySum->cd(3);
    hCFLLALAL[taskNumber]->GetXaxis()->SetRangeUser(0,0.5);
    hCFLLALAL[taskNumber]->SetLineColor(taskNumber+1);
    hCFLLALAL[taskNumber]->SetMarkerColor(taskNumber+1);
    if(taskNumber == 0) {hCFLLALAL[taskNumber]->SetMarkerStyle(23); hCFLLALAL[taskNumber]->SetMarkerSize(1.3);hCFLLALAL[taskNumber]->SetTitle("#Lambda#Lambda #oplus #bar{#Lambda}#bar{#Lambda}"); hCFLLALAL[taskNumber]->Draw();}
    else hCFLLALAL[taskNumber]->Draw("same");
  }
   if(rebinPP == 1 && rebinLP == 1 && rebinLL == 1) CanCFsummarySum->Print("CF.pdf");
   else CanCFsummarySum->Print(Form("CF_%i_%i_%i.pdf", rebinPP, rebinLP, rebinLL));


  //Lets build all ratios:
  
  TH1F *hCFRatioPPAPAP[NTasksProtonVariations];
  TH1F *hCFRatioPLAPAL[NTasksV0Variations];
  TH1F *hCFRatioLLALAL[NTasksV0Variations];
  
  //All ratios of default to variations C_default / C_variations for pp
  for(int protonVariation = 0; protonVariation < NTasksProtonVariations; protonVariation++)
  {
    //First clone the default correlation function:
    HistName = "hCFRatioPPAPAP_task_0_";
    HistName += protonVariation+1;

    hCFRatioPPAPAP[protonVariation] = (TH1F*)hCFPPAPAP[0]->Clone(HistName.Data());
    hCFRatioPPAPAP[protonVariation]->Divide(hCFPPAPAP[protonVariation+1]);
  }

  //All ratios of default to variations C_default / C_variations for pL
  for(int V0Variation = 0; V0Variation < NTasksV0Variations; V0Variation++)
  {
    //First clone the default correlation function:
    HistName = "hCFRatioPLAPAL_task_0_";
    HistName += V0Variation+1;

    hCFRatioPLAPAL[V0Variation] = (TH1F*)hCFPLAPAL[0]->Clone(HistName.Data());
    hCFRatioPLAPAL[V0Variation]->Divide(hCFPLAPAL[V0Variation+1]);
  }

  for(int V0Variation = 0; V0Variation < NTasksV0Variations; V0Variation++)
  {
    //First clone the default correlation function:
    HistName = "hCFRatioLLALAL_task_0_";
    HistName += V0Variation+1;

    hCFRatioLLALAL[V0Variation] = (TH1F*)hCFLLALAL[0]->Clone(HistName.Data());
    hCFRatioLLALAL[V0Variation]->Divide(hCFLLALAL[V0Variation+1]);
  }

  
  
  //Evaluate the systematic error for pp pairs. Since we vary the values up and down we take the larger systematic uncertainty per bin
  //To not confuse oneself we use for every bin an own array to see what happens

  std::vector<const char*> protonCut(8);
  protonCut[0] = "#it{p}_{T} down";
  protonCut[1] = "#it{p}_{T} up";
  protonCut[2] = "#eta up";
  protonCut[3] = "#eta down";
  protonCut[4] = "n#sigma up";
  protonCut[5] = "n#sigma down";
  protonCut[6] = "FilterBit";
  protonCut[7] = "TPC cluster down";
  protonCut[8] = "TPC cluster up";

  const int NkstarBinsPPAPAP = hCFRatioPPAPAP[0]->GetXaxis()->FindBin(fitRangeRightPPAPAP);
  double sysErrorPPAPAP_ptValueDown[NkstarBinsPPAPAP];
  double sysErrorPPAPAP_ptValueUp[NkstarBinsPPAPAP];
  double sysErrorPPAPAP_EtaValueUp[NkstarBinsPPAPAP];
  double sysErrorPPAPAP_EtaValueDown[NkstarBinsPPAPAP];
  double sysErrorPPAPAP_NSigmaValueUp[NkstarBinsPPAPAP];
  double sysErrorPPAPAP_NSigmaValueDown[NkstarBinsPPAPAP];
  double sysErrorPPAPAP_FilterBit[NkstarBinsPPAPAP];
  double sysErrorPPAPAP_TPCcluster[NkstarBinsPPAPAP];
  double sysErrorTotalPPAPAP[NkstarBinsPPAPAP];
  TH1F *hCFPPAPAP_errBudget[NTasksProtonVariations];
  TF1* fFitRatioPPAPAP[NTasksProtonVariations];

  for(int kstarBin = 0; kstarBin < NkstarBinsPPAPAP; kstarBin++) {
    sysErrorPPAPAP_ptValueDown[kstarBin] = 0;
    sysErrorPPAPAP_ptValueUp[kstarBin] = 0;
    sysErrorPPAPAP_EtaValueUp[kstarBin] = 0;
    sysErrorPPAPAP_EtaValueDown[kstarBin] = 0;
    sysErrorPPAPAP_NSigmaValueUp[kstarBin] = 0;
    sysErrorPPAPAP_NSigmaValueDown[kstarBin] = 0;
    sysErrorPPAPAP_FilterBit[kstarBin] = 0;
    sysErrorPPAPAP_TPCcluster[kstarBin] = 0;
    sysErrorTotalPPAPAP[kstarBin] = 0;
  }


  // NEW : BARLOW CHECK

  TH1F *hCFPPAPAP_Barlow[NTasksProtonVariations];
  //Plot all variations for proton pairs:
  TCanvas *CanCFPPAPAPBarlow = new TCanvas("CanCFPPAPAPBarlow","CanCFPPAPAPBarlow",0,0,1500,1000);
  CanCFPPAPAPBarlow->Divide(3,3);

  for(int protonVariation = 0; protonVariation < NTasksProtonVariations; protonVariation++) {
//    for(int protonVariation = 0; protonVariation < NTasksProtonVariations; protonVariation++) {
    if(protonVariation > NTasksProtonVariations-1) break;
    CanCFPPAPAPBarlow->cd(protonVariation+1);
    HistName = "hCFPPAPAP_Barlow_";
    HistName += protonVariation+1;
    hCFPPAPAP_Barlow[protonVariation] = new TH1F(HistName.Data(),"contains relative error",NkstarBinsPPAPAP,0,hCFRatioPPAPAP[protonVariation]->GetBinCenter(NkstarBinsPPAPAP));

    //evaluate the systematic error from the fit:
    for(int kstarBin = 0; kstarBin < NkstarBinsPPAPAP; kstarBin++)
    {
      double nominalVal = hCFPPAPAP[0]->GetBinContent(kstarBin+1);
      double variationVal= hCFPPAPAP[protonVariation+1]->GetBinContent(kstarBin+1);
      double statErrVariation = std::sqrt(hCFPPAPAP[protonVariation+1]->GetBinError(kstarBin+1)*hCFPPAPAP[protonVariation+1]->GetBinError(kstarBin+1) + hCFPPAPAP[0]->GetBinError(kstarBin+1)*hCFPPAPAP[0]->GetBinError(kstarBin+1));
      if(statErrVariation < 1E-35) continue;
      hCFPPAPAP_Barlow[protonVariation]->SetBinContent(kstarBin+1, std::abs(variationVal - nominalVal)/statErrVariation);
      hCFPPAPAP_Barlow[protonVariation]->SetBinError(kstarBin+1, 0);
    }

    DrawHist(hCFPPAPAP_Barlow[protonVariation],"k* (GeV/#it{c})","n#sigma",0,1,1.);
    HistoStyleMany(hCFPPAPAP_Barlow[protonVariation]);
    hCFPPAPAP_Barlow[protonVariation]->SetTitle(protonCut[protonVariation]);
    hCFPPAPAP_Barlow[protonVariation]->GetXaxis()->SetRangeUser(0, 0.12);
  }
   if(rebinPP == 1) CanCFPPAPAPBarlow->Print("PP_barlow.pdf");
  else CanCFPPAPAPBarlow->Print(Form("PP_barlow_%i.pdf", rebinPP));


  //As systematic error we define the difference between the default and the variation correlation function:
  // Delta = |C_default - C_variation| -> Delta = C_default |1 - C_variation/C_default|
  
  //Plot all variations for proton pairs:
  TCanvas *CanCFPPAPAPRatios = new TCanvas("CanCFPPAPAPRatios","CanCFPPAPAPRatios",0,0,1500,1000);
  CanCFPPAPAPRatios->Divide(3,3);
  for(int protonVariation = 0; protonVariation < NTasksProtonVariations; protonVariation++)
  {
    fFitRatioPPAPAP[protonVariation] = new TF1(Form("fFitRatioPPAPAP_%i", protonVariation),fit,0., 0.12);
    CanCFPPAPAPRatios->cd(protonVariation+1);
    hCFRatioPPAPAP[protonVariation]->GetXaxis()->SetRangeUser(0.,0.12);
    DrawHist(hCFRatioPPAPAP[protonVariation],"k* (GeV/#it{c})","#it{C}(k*)_{default}/C(k*)_{variation}",0,1,1.);
    hCFRatioPPAPAP[protonVariation]->SetTitle(protonCut[protonVariation]);
    if(useFit) hCFRatioPPAPAP[protonVariation]->Fit(fFitRatioPPAPAP[protonVariation],"NRQ", "", 0, fitRangeRightPPAPAP);
    fFitRatioPPAPAP[protonVariation]->Draw("same");

    HistName = "hCFPPAPAP_errBudget_";
    HistName += protonVariation;
    hCFPPAPAP_errBudget[protonVariation] = new TH1F(HistName.Data(),"contains relative error",NkstarBinsPPAPAP,0,hCFRatioPPAPAP[protonVariation]->GetBinCenter(NkstarBinsPPAPAP));

    if(hCFPPAPAP_Barlow[protonVariation]->GetMaximum() < thresholdBarlow) continue;

    //evaluate the systematic error from the fit:
    for(int kstarBin = 0; kstarBin < NkstarBinsPPAPAP; kstarBin++)
    {
      double kstarValue = hCFRatioPPAPAP[protonVariation]->GetBinCenter(kstarBin+1);
      double SysError;
      if(useFit) SysError = fFitRatioPPAPAP[protonVariation]->Eval(kstarValue);
      else SysError = hCFRatioPPAPAP[protonVariation]->GetBinContent(kstarBin+1);

      if(SysError-1.f < 1E-35) continue;

      //We have to take the inverse in the calculations of syserror since the ratio is defined as C_default / C_variation
      if(protonVariation == 0)//Pt down
      {
        sysErrorPPAPAP_ptValueDown[kstarBin] = 0.;
        sysErrorPPAPAP_ptValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPPAPAP[0]->GetBinContent(kstarBin+1);
        hCFPPAPAP_errBudget[protonVariation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(protonVariation == 1)//Pt up
      {
        sysErrorPPAPAP_ptValueUp[kstarBin] = 0.;
        sysErrorPPAPAP_ptValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPPAPAP[0]->GetBinContent(kstarBin+1);
        hCFPPAPAP_errBudget[protonVariation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(protonVariation == 2)//Eta up
      {
        sysErrorPPAPAP_EtaValueUp[kstarBin] = 0.;
        sysErrorPPAPAP_EtaValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPPAPAP[0]->GetBinContent(kstarBin+1);
        hCFPPAPAP_errBudget[protonVariation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(protonVariation == 3)//Eta down
      {
        sysErrorPPAPAP_EtaValueDown[kstarBin] = 0.;
        sysErrorPPAPAP_EtaValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPPAPAP[0]->GetBinContent(kstarBin+1);
        hCFPPAPAP_errBudget[protonVariation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(protonVariation == 4)//Nsigma up
      {
        sysErrorPPAPAP_NSigmaValueUp[kstarBin] = 0.;
        sysErrorPPAPAP_NSigmaValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPPAPAP[0]->GetBinContent(kstarBin+1);
        hCFPPAPAP_errBudget[protonVariation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(protonVariation == 5)//Nsigma down
      {
        sysErrorPPAPAP_NSigmaValueDown[kstarBin] = 0.;
        sysErrorPPAPAP_NSigmaValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPPAPAP[0]->GetBinContent(kstarBin+1);
        hCFPPAPAP_errBudget[protonVariation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(protonVariation == 6)//Filterbit
      {
        sysErrorPPAPAP_FilterBit[kstarBin] = 0.;
        sysErrorPPAPAP_FilterBit[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPPAPAP[0]->GetBinContent(kstarBin+1);
        hCFPPAPAP_errBudget[protonVariation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(protonVariation == 7)//TPC cluster
      {
        sysErrorPPAPAP_TPCcluster[kstarBin] = 0.;
        sysErrorPPAPAP_TPCcluster[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPPAPAP[0]->GetBinContent(kstarBin+1);
        hCFPPAPAP_errBudget[protonVariation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    }
  }
   if(rebinPP == 1) CanCFPPAPAPRatios->Print("PP.pdf");
  else CanCFPPAPAPRatios->Print(Form("PP_%i.pdf", rebinPP));
   std::cout << "Made it up to here \n";
  
  //Evaluate the systematic error for pL pairs. Since we vary the values up and down we take the larger systematic uncertainty per bin
  //To not confuse oneself we use for every bin an own array to see what happens

  std::vector<const char*> v0Cut(20
                                 );
  v0Cut[0] = "#it{p}_{T} down";
  v0Cut[1] = "#it{p}_{T} up";
  v0Cut[2] = "#eta up";
  v0Cut[3] = "#eta down";
  v0Cut[4] = "n#sigma up";
  v0Cut[5] = "n#sigma down";
  v0Cut[6] = "FilterBit";
  v0Cut[7] = "TPC cluster down";
  v0Cut[8] = "TPC cluster up";
  v0Cut[9] = "#Xi d_{daug} down";
  v0Cut[10] = "#Xi d_{bach,PV} up";
  v0Cut[11] = "#Xi CPA up";
  v0Cut[12] = "#Xi Min. Transverse Radius Up";
  v0Cut[13] = "#Lambda d_{daug} down";
  v0Cut[14] = "#Lambda CPA up";
  v0Cut[15] = "#Lambda Min. Transverse Radius Up";
  v0Cut[16] = "#Xi d_{#Lambda,PV} up";
  v0Cut[17] = "#Lambda d_{#Lambda,PV} up";
  v0Cut[18] = "V0 #eta down";
  v0Cut[19] = "V0 #eta up";
  v0Cut[20] = "V0 n#sigma down";

  const int NkstarBinsPLAPAL = hCFRatioPLAPAL[0]->GetXaxis()->FindBin(fitRangeRightPLAPAL);
  //systematic error due to proton selection:
  double sysErrorPLAPAL_ptValueDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_ptValueUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_EtaValueUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_EtaValueDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_NSigmaValueUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_NSigmaValueDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_FilterBit[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_TPCcluster[NkstarBinsPLAPAL];
  //systematic error due to V0 selection:
  double sysErrorPLAPAL_V0ptValueDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0ptValueUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0CPAUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0NsigmaDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0TPCCluster80[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0EtaValueUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0EtaValueDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0trackDistanceDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0trackDistanctToPVUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0trackDistanctToPVUp10[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0trackDistanctToPVUp11[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_V0trackDistanctToPVUp12[NkstarBinsPLAPAL];
  double sysErrorTotalPLAPAL[NkstarBinsPLAPAL];
  TH1F *hCFPLAPAL_errBudget[NTasksV0Variations];
  TF1* fFitRatioPLAPAL[NTasksV0Variations];

  for(int kstarBin = 0; kstarBin < NkstarBinsPLAPAL; kstarBin++) {
    sysErrorPLAPAL_ptValueDown[kstarBin] = 0;
    sysErrorPLAPAL_ptValueUp[kstarBin] = 0;
    sysErrorPLAPAL_EtaValueUp[kstarBin] = 0;
    sysErrorPLAPAL_EtaValueDown[kstarBin] = 0;
    sysErrorPLAPAL_NSigmaValueUp[kstarBin] = 0;
    sysErrorPLAPAL_NSigmaValueDown[kstarBin] = 0;
    sysErrorPLAPAL_FilterBit[kstarBin] = 0;
    sysErrorPLAPAL_TPCcluster[kstarBin] = 0;
    //systematic error due to V0 selection:
    sysErrorPLAPAL_V0ptValueDown[kstarBin] = 0;
    sysErrorPLAPAL_V0ptValueUp[kstarBin] = 0;
    sysErrorPLAPAL_V0CPAUp[kstarBin] = 0;
    sysErrorPLAPAL_V0NsigmaDown[kstarBin] = 0;
    sysErrorPLAPAL_V0TPCCluster80[kstarBin] = 0;
    sysErrorPLAPAL_V0EtaValueUp[kstarBin] = 0;
    sysErrorPLAPAL_V0EtaValueDown[kstarBin] = 0;
    sysErrorPLAPAL_V0trackDistanceDown[kstarBin] = 0;
    sysErrorPLAPAL_V0trackDistanctToPVUp[kstarBin] = 0;
    sysErrorPLAPAL_V0trackDistanctToPVUp10[kstarBin] = 0;
    sysErrorPLAPAL_V0trackDistanctToPVUp11[kstarBin] = 0;
    sysErrorPLAPAL_V0trackDistanctToPVUp12[kstarBin] = 0;
    sysErrorTotalPLAPAL[kstarBin] = 0;
  }
  std::cout << "In Front of the BARLOW Check \n";
  // NEW : BARLOW CHECK

  TH1F *hCFPLAPAL_Barlow[NTasksV0Variations];
  //Plot all variations for proton pairs:
  TCanvas *CanCFPLAPALBarlow = new TCanvas("CanCFPLAPALBarlow","CanCFPLAPALBarlow",0,0,1500,1000);
  CanCFPLAPALBarlow->Divide(5,4);
  std::cout << "NTasksV0Variations:   " << NTasksV0Variations << std::endl;
  for(int V0Variation = 0; V0Variation < NTasksV0Variations; V0Variation++)
  {
    CanCFPLAPALBarlow->cd(V0Variation+1);

    HistName = "hCFPLAPAL_Barlow_";
    HistName += V0Variation;
    hCFPLAPAL_Barlow[V0Variation] = new TH1F(HistName.Data(),"contains relative error",NkstarBinsPLAPAL,0,hCFRatioPLAPAL[V0Variation]->GetBinCenter(NkstarBinsPLAPAL));

    //evaluate the systematic error from the fit:
    for(int kstarBin = 0; kstarBin < NkstarBinsPLAPAL; kstarBin++) {
      double nominalVal = hCFPLAPAL[0]->GetBinContent(kstarBin+1);
      double variationVal= hCFPLAPAL[V0Variation+1]->GetBinContent(kstarBin+1);
      double statErrVariation = std::sqrt(hCFPLAPAL[V0Variation+1]->GetBinError(kstarBin+1)*hCFPLAPAL[V0Variation+1]->GetBinError(kstarBin+1) + hCFPLAPAL[0]->GetBinError(kstarBin+1)*hCFPLAPAL[0]->GetBinError(kstarBin+1));
      hCFPLAPAL_Barlow[V0Variation]->SetBinContent(kstarBin+1, std::abs(variationVal - nominalVal)/statErrVariation);
      hCFPLAPAL_Barlow[V0Variation]->SetBinError(kstarBin+1, 0);
    }

    DrawHist(hCFPLAPAL_Barlow[V0Variation],"k* (GeV/#it{c})","n#sigma",0,1,1.);
    HistoStyleMany(hCFPLAPAL_Barlow[V0Variation]);
    hCFPLAPAL_Barlow[V0Variation]->GetXaxis()->SetRangeUser(0, 0.25);
    hCFPLAPAL_Barlow[V0Variation]->SetTitle(v0Cut[V0Variation]);
  }
   if(rebinLP == 1) CanCFPLAPALBarlow->Print("PL_barlow.pdf");
  else CanCFPLAPALBarlow->Print(Form("PL_barlow_%i.pdf", rebinLP));

  //Plot all variations for proton-V0 pairs:
  TCanvas *CanCFPLAPALRatios = new TCanvas("CanCFPLAPALRatios","CanCFPLAPALRatios",0,0,1500,1000);
  CanCFPLAPALRatios->Divide(5,4);
  for(int V0Variation = 0; V0Variation < NTasksV0Variations; V0Variation++)
  {
    CanCFPLAPALRatios->cd(V0Variation+1);
    fFitRatioPLAPAL[V0Variation] = new TF1(Form("fFitRatioPLAPAL_%i", V0Variation) ,fit,0., 0.2);
    hCFRatioPLAPAL[V0Variation]->GetXaxis()->SetRangeUser(0.,0.2);
    DrawHist(hCFRatioPLAPAL[V0Variation],"k* (GeV/#it{c})","#it{C}(k*)_{default}/C(k*)_{variation}",0,1,1.);
    hCFRatioPLAPAL[V0Variation]->SetTitle(v0Cut[V0Variation]);
    if(useFit) hCFRatioPLAPAL[V0Variation]->Fit(fFitRatioPLAPAL[V0Variation],"NRQ", "", 0.f, fitRangeRightPLAPAL);
    fFitRatioPLAPAL[V0Variation]->Draw("same");

    HistName = "hCFPLAPAL_errBudget_";
    HistName += V0Variation;
    hCFPLAPAL_errBudget[V0Variation] = new TH1F(HistName.Data(),"contains relative error",NkstarBinsPLAPAL,0,hCFRatioPLAPAL[V0Variation]->GetBinCenter(NkstarBinsPLAPAL));

    if(hCFPLAPAL_Barlow[V0Variation]->GetMaximum() < thresholdBarlow) continue;

    //evaluate the systematic error from the fit:
    for(int kstarBin = 0; kstarBin < NkstarBinsPLAPAL; kstarBin++)
    {
      double kstarValue = hCFRatioPLAPAL[V0Variation]->GetBinCenter(kstarBin+1);
      double SysError;
      if(useFit) SysError = fFitRatioPLAPAL[V0Variation]->Eval(kstarValue);
      else SysError = hCFRatioPLAPAL[V0Variation]->GetBinContent(kstarBin+1);

      //We have to take the inverse in the calculations of syserror since the ratio is defined as C_default / C_variation

      if(V0Variation == 0)//proton Pt down
      {
        sysErrorPLAPAL_ptValueDown[kstarBin] = 0.;
        sysErrorPLAPAL_ptValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 1)//proton Pt up
      {
        sysErrorPLAPAL_ptValueUp[kstarBin] = 0.;
        sysErrorPLAPAL_ptValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 2)//proton Eta up
      {
        sysErrorPLAPAL_EtaValueUp[kstarBin] = 0.;
        sysErrorPLAPAL_EtaValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 3)//proton Eta down
      {
        sysErrorPLAPAL_EtaValueDown[kstarBin] = 0.;
        sysErrorPLAPAL_EtaValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 4)//proton Nsigma up
      {
        sysErrorPLAPAL_NSigmaValueUp[kstarBin] = 0.;
        sysErrorPLAPAL_NSigmaValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 5)//proton Nsigma down
      {
        sysErrorPLAPAL_NSigmaValueDown[kstarBin] = 0.;
        sysErrorPLAPAL_NSigmaValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 6)//proton Filterbit
      {
        sysErrorPLAPAL_FilterBit[kstarBin] = 0.;
        sysErrorPLAPAL_FilterBit[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 7)//proton TPC cluster
      {
        sysErrorPLAPAL_TPCcluster[kstarBin] = 0.;
        sysErrorPLAPAL_TPCcluster[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 8)//V0 Pt down
      {
        sysErrorPLAPAL_V0ptValueDown[kstarBin] = 0.;
        sysErrorPLAPAL_V0ptValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 9)//V0 Pt up
      {
        sysErrorPLAPAL_V0ptValueUp[kstarBin] = 0.;
        sysErrorPLAPAL_V0ptValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 10)//V0 CPA up
      {
        sysErrorPLAPAL_V0CPAUp[kstarBin] = 0.;
        sysErrorPLAPAL_V0CPAUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 11)//V0 NSigma down
      {
        sysErrorPLAPAL_V0NsigmaDown[kstarBin] = 0.;
        sysErrorPLAPAL_V0NsigmaDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 12)//V0 TPC cluster up
      {
        sysErrorPLAPAL_V0TPCCluster80[kstarBin] = 0.;
        sysErrorPLAPAL_V0TPCCluster80[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 13)//V0 eta up
      {
        sysErrorPLAPAL_V0EtaValueUp[kstarBin] = 0.;
        sysErrorPLAPAL_V0EtaValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 14)//V0 eta down
      {
        sysErrorPLAPAL_V0EtaValueDown[kstarBin] = 0.;
        sysErrorPLAPAL_V0EtaValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 15)//V0 track distance down
      {
        sysErrorPLAPAL_V0trackDistanceDown[kstarBin] = 0.;
        sysErrorPLAPAL_V0trackDistanceDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 16)//V0 impact parameter up
      {
        sysErrorPLAPAL_V0trackDistanctToPVUp[kstarBin] = 0.;
        sysErrorPLAPAL_V0trackDistanctToPVUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 17)//V0 impact parameter up
      {
        sysErrorPLAPAL_V0trackDistanctToPVUp10[kstarBin] = 0.;
        sysErrorPLAPAL_V0trackDistanctToPVUp10[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 18)//V0 impact parameter up
      {
        sysErrorPLAPAL_V0trackDistanctToPVUp11[kstarBin] = 0.;
        sysErrorPLAPAL_V0trackDistanctToPVUp11[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 19)//V0 impact parameter up
      {
        sysErrorPLAPAL_V0trackDistanctToPVUp12[kstarBin] = 0.;
        sysErrorPLAPAL_V0trackDistanctToPVUp12[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    }

  }
   if(rebinLP == 1) CanCFPLAPALRatios->Print("PL_.pdf");
  else CanCFPLAPALRatios->Print(Form("PL_%i.pdf", rebinLP));
   std::cout << "Done with pXi\n";
  // LAMBDA - LAMBDA


  const int NkstarBinsLLALAL = hCFRatioLLALAL[0]->GetXaxis()->FindBin(fitRangeRightLLALAL);
  double sysErrorLLALAL_V0ptValueDown[NkstarBinsLLALAL];
  double sysErrorLLALAL_V0ptValueUp[NkstarBinsLLALAL];
  double sysErrorLLALAL_V0CPAUp[NkstarBinsLLALAL];
  double sysErrorLLALAL_V0NsigmaDown[NkstarBinsLLALAL];
  double sysErrorLLALAL_V0TPCCluster80[NkstarBinsLLALAL];
  double sysErrorLLALAL_V0EtaValueUp[NkstarBinsLLALAL];
  double sysErrorLLALAL_V0EtaValueDown[NkstarBinsLLALAL];
  double sysErrorLLALAL_V0trackDistanceDown[NkstarBinsLLALAL];
  double sysErrorLLALAL_V0trackDistanctToPVUp[NkstarBinsLLALAL];
  double sysErrorTotalLLALAL[NkstarBinsLLALAL];
  TH1F *hCFLLALAL_errBudget[NTasksV0Variations];
  TF1* fFitRatioLLALAL[NTasksV0Variations];

  for(int kstarBin = 0; kstarBin < NkstarBinsLLALAL; kstarBin++) {
    sysErrorLLALAL_V0ptValueDown[kstarBin] = 0;
    sysErrorLLALAL_V0ptValueUp[kstarBin] = 0;
    sysErrorLLALAL_V0CPAUp[kstarBin] = 0;
    sysErrorLLALAL_V0NsigmaDown[kstarBin] = 0;
    sysErrorLLALAL_V0TPCCluster80[kstarBin] = 0;
    sysErrorLLALAL_V0EtaValueUp[kstarBin] = 0;
    sysErrorLLALAL_V0EtaValueDown[kstarBin] = 0;
    sysErrorLLALAL_V0trackDistanceDown[kstarBin] = 0;
    sysErrorLLALAL_V0trackDistanctToPVUp[kstarBin] = 0;
    sysErrorTotalLLALAL[kstarBin] = 0;
  }

  TH1F *hCFLLALAL_Barlow[NTasksV0Variations];
  //Plot all variations for proton pairs:
  TCanvas *CanCFLLALALBarlow = new TCanvas("CanCFLLALALBarlow","CanCFLLALALBarlow",0,0,1500,1000);
  CanCFLLALALBarlow->Divide(3,3);

  for(int V0Variation = NTasksProtonVariations; V0Variation < NTasksV0Variations; V0Variation++)
  {
    CanCFLLALALBarlow->cd(V0Variation-NTasksProtonVariations+1);

    HistName = "hCFLLALAL_Barlow_";
    HistName += V0Variation;
    hCFLLALAL_Barlow[V0Variation] = new TH1F(HistName.Data(),"contains relative error",NkstarBinsLLALAL,0,hCFRatioLLALAL[V0Variation]->GetBinCenter(NkstarBinsLLALAL));

    //evaluate the systematic error from the fit:
    for(int kstarBin = 0; kstarBin < NkstarBinsLLALAL; kstarBin++)
    {
      double nominalVal = hCFLLALAL[0]->GetBinContent(kstarBin+1);
      double variationVal= hCFLLALAL[V0Variation+1]->GetBinContent(kstarBin+1);
      double statErrVariation = std::sqrt(hCFLLALAL[V0Variation+1]->GetBinError(kstarBin+1)*hCFLLALAL[V0Variation+1]->GetBinError(kstarBin+1) + hCFLLALAL[0]->GetBinError(kstarBin+1)*hCFLLALAL[0]->GetBinError(kstarBin+1));
      hCFLLALAL_Barlow[V0Variation]->SetBinContent(kstarBin+1, std::abs(variationVal - nominalVal)/statErrVariation);
      hCFLLALAL_Barlow[V0Variation]->SetBinError(kstarBin+1, 0);
    }

    DrawHist(hCFLLALAL_Barlow[V0Variation],"k* (GeV/#it{c})","n#sigma",0,1,1.);
    HistoStyleMany(hCFLLALAL_Barlow[V0Variation]);
    hCFLLALAL_Barlow[V0Variation]->GetXaxis()->SetRangeUser(0, 0.25);
    hCFLLALAL_Barlow[V0Variation]->SetTitle(v0Cut[V0Variation]);
  }
   if(rebinLL == 1) CanCFLLALALBarlow->Print("LL_barlow.pdf");
  else CanCFLLALALBarlow->Print(Form("LL_barlow_%i.pdf", rebinLL));


  //Plot all variations for proton-V0 pairs:
  TCanvas *CanCFLLALALRatios = new TCanvas("CanCFLLALALRatios","CanCFLLALALRatios",0,0,1500,1000);
  CanCFLLALALRatios->Divide(3,3);
  for(int V0Variation = NTasksProtonVariations; V0Variation < NTasksV0Variations; V0Variation++)
  {
    fFitRatioLLALAL[V0Variation] = new TF1(Form("fFitRatioLLALAL_%i", V0Variation),fit,0., 0.2);
    CanCFLLALALRatios->cd(V0Variation-NTasksProtonVariations+1);
    hCFRatioLLALAL[V0Variation]->GetXaxis()->SetRangeUser(0.,0.2);
    DrawHist(hCFRatioLLALAL[V0Variation],"k* (GeV/#it{c})","#it{C}(k*)_{default}/C(k*)_{variation}",0,1,1.);
    hCFRatioLLALAL[V0Variation]->SetTitle(v0Cut[V0Variation]);

    if(useFit) hCFRatioLLALAL[V0Variation]->Fit(fFitRatioLLALAL[V0Variation],"NRQ", "", 0.f, fitRangeRightLLALAL);
    fFitRatioLLALAL[V0Variation]->Draw("same");

    HistName = "hCFLLALAL_errBudget_";
    HistName += V0Variation;
    hCFLLALAL_errBudget[V0Variation] = new TH1F(HistName.Data(),"contains relative error",NkstarBinsLLALAL,0,hCFRatioLLALAL[V0Variation]->GetBinCenter(NkstarBinsLLALAL));

    if(hCFLLALAL_Barlow[V0Variation]->GetMaximum() < thresholdBarlow ) continue;

    //evaluate the systematic error from the fit:
    for(int kstarBin = 0; kstarBin < NkstarBinsLLALAL; kstarBin++)
    {
      double kstarValue = hCFRatioLLALAL[V0Variation]->GetBinCenter(kstarBin+1);
      double SysError;
      if(useFit) SysError = fFitRatioLLALAL[V0Variation]->Eval(kstarValue);
      else SysError = hCFRatioLLALAL[V0Variation]->GetBinContent(kstarBin+1);

      //We have to take the inverse in the calculations of syserror since the ratio is defined as C_default / C_variation

    if(V0Variation == NTasksProtonVariations+0)//V0 Pt down
      {
        sysErrorLLALAL_V0ptValueDown[kstarBin] = 0.;
        sysErrorLLALAL_V0ptValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFLLALAL[0]->GetBinContent(kstarBin+1);
        hCFLLALAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    else if(V0Variation == NTasksProtonVariations+1)//V0 Pt up
      {
        sysErrorLLALAL_V0ptValueUp[kstarBin] = 0.;
        sysErrorLLALAL_V0ptValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFLLALAL[0]->GetBinContent(kstarBin+1);
        hCFLLALAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    else if(V0Variation == NTasksProtonVariations+2)//V0 CPA up
      {
        sysErrorLLALAL_V0CPAUp[kstarBin] = 0.;
        sysErrorLLALAL_V0CPAUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFLLALAL[0]->GetBinContent(kstarBin+1);
        hCFLLALAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    else if(V0Variation == NTasksProtonVariations+3)//V0 NSigma down
      {
        sysErrorLLALAL_V0NsigmaDown[kstarBin] = 0.;
        sysErrorLLALAL_V0NsigmaDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFLLALAL[0]->GetBinContent(kstarBin+1);
        hCFLLALAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    else if(V0Variation == NTasksProtonVariations+4)//V0 TPC cluster up
      {
        sysErrorLLALAL_V0TPCCluster80[kstarBin] = 0.;
        sysErrorLLALAL_V0TPCCluster80[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFLLALAL[0]->GetBinContent(kstarBin+1);
        hCFLLALAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    else if(V0Variation == NTasksProtonVariations+5)//V0 eta up
      {
        sysErrorLLALAL_V0EtaValueUp[kstarBin] = 0.;
        sysErrorLLALAL_V0EtaValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFLLALAL[0]->GetBinContent(kstarBin+1);
        hCFLLALAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    else if(V0Variation == NTasksProtonVariations+6)//V0 eta down
      {
        sysErrorLLALAL_V0EtaValueDown[kstarBin] = 0.;
        sysErrorLLALAL_V0EtaValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFLLALAL[0]->GetBinContent(kstarBin+1);
        hCFLLALAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    else if(V0Variation == NTasksProtonVariations+7)//V0 track distance down
      {
        sysErrorLLALAL_V0trackDistanceDown[kstarBin] = 0.;
        sysErrorLLALAL_V0trackDistanceDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFLLALAL[0]->GetBinContent(kstarBin+1);
        hCFLLALAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    else if(V0Variation == NTasksProtonVariations+8)//V0 impact parameter up
      {
        sysErrorLLALAL_V0trackDistanctToPVUp[kstarBin] = 0.;
        sysErrorLLALAL_V0trackDistanctToPVUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFLLALAL[0]->GetBinContent(kstarBin+1);
        hCFLLALAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
    }
  }
   if(rebinLL == 1) CanCFLLALALRatios->Print("LL.pdf");
  else CanCFLLALALRatios->Print(Form("LL_%i.pdf", rebinLL));

  double xval_CFPPAPAP[NkstarBinsPPAPAP];
  double yval_CFPPAPAP[NkstarBinsPPAPAP];
  double xerr_CFPPAPAP[NkstarBinsPPAPAP];
  
  TH1F *sysErrorC2PPAPAP = (TH1F*)hCFAPAP[0]->Clone("C2totalsysPP");//Clone one pp histogram
  
  //Lets combine them together to a total systematic uncertainty for pp:
  for(int kstarBin = 0; kstarBin < NkstarBinsPPAPAP; kstarBin++)
  {
    //calculate the average of two variations if applicable:

    double sysErrorAveragePt = 0.5*(sysErrorPPAPAP_ptValueUp[kstarBin] + sysErrorPPAPAP_ptValueDown[kstarBin]);
    double sysErrorAverageEta = 0.5*(sysErrorPPAPAP_EtaValueUp[kstarBin] + sysErrorPPAPAP_EtaValueDown[kstarBin]);
    double sysErrorAverageNSigma = 0.5*(sysErrorPPAPAP_NSigmaValueUp[kstarBin] + sysErrorPPAPAP_NSigmaValueDown[kstarBin]);

    //Add them up quadratically:

    sysErrorTotalPPAPAP[kstarBin] = pow(sysErrorAveragePt,2.) + pow(sysErrorAverageEta,2.) + pow(sysErrorAverageNSigma,2.) + pow(sysErrorPPAPAP_FilterBit[kstarBin],2.) + pow(sysErrorPPAPAP_TPCcluster[kstarBin],2.);
    sysErrorTotalPPAPAP[kstarBin] = TMath::Sqrt(sysErrorTotalPPAPAP[kstarBin]);

    //Put this information into arrays for plotting:
    xval_CFPPAPAP[kstarBin] = hCFPPAPAP[0]->GetBinCenter(kstarBin+1);
    yval_CFPPAPAP[kstarBin] = hCFPPAPAP[0]->GetBinContent(kstarBin+1);
    //xerr_CFPPAPAP[kstarBin] = 0.0015;
    xerr_CFPPAPAP[kstarBin] = hCFPPAPAP[0]->GetBinWidth(1)/2.;

    hSystPP->Fill(hCFPPAPAP[0]->GetBinCenter(kstarBin+1), sysErrorTotalPPAPAP[kstarBin]);

  }
  
  //Set them in a histogram which can be later used as input to the fitting code:
  for(int kstarBin = 0;kstarBin< sysErrorC2PPAPAP->GetNbinsX(); kstarBin++)
  {
    if(kstarBin < NkstarBinsPPAPAP)
    {
      sysErrorC2PPAPAP->SetBinContent(kstarBin+1,sysErrorTotalPPAPAP[kstarBin]);
      sysErrorC2PPAPAP->SetBinError(kstarBin+1,0.);
    }
    else
    {
      sysErrorC2PPAPAP->SetBinContent(kstarBin+1,0.);
      sysErrorC2PPAPAP->SetBinError(kstarBin+1,0.);
    }
  }

  double xval_CFPLAPAL[NkstarBinsPLAPAL];
  double yval_CFPLAPAL[NkstarBinsPLAPAL];
  double xerr_CFPLAPAL[NkstarBinsPLAPAL];
  
  TH1F *sysErrorC2PLAPAL = (TH1F*)hCFAPAL[0]->Clone("C2totalsysPL");

  //Lets combine them together to a total systematic uncertainty for pp:
  for(int kstarBin = 0; kstarBin < NkstarBinsPLAPAL; kstarBin++)
  {
    //calculate the average of two variations if applicable:
    //systematic error from proton selection averaged:
    double sysErrorAveragePt = 0.5*(sysErrorPLAPAL_ptValueUp[kstarBin] + sysErrorPLAPAL_ptValueDown[kstarBin]);
    double sysErrorAverageEta = 0.5*(sysErrorPLAPAL_EtaValueUp[kstarBin] + sysErrorPLAPAL_EtaValueDown[kstarBin]);
    double sysErrorAverageNSigma = 0.5*(sysErrorPLAPAL_NSigmaValueUp[kstarBin] + sysErrorPLAPAL_NSigmaValueDown[kstarBin]);
    //systematic error from V0 selection averaged:
    double sysErrorAverageV0Pt = 0.5*(sysErrorPLAPAL_V0ptValueUp[kstarBin] + sysErrorPLAPAL_V0ptValueDown[kstarBin]);
    double sysErrorAverageV0Eta = 0.5*(sysErrorPLAPAL_V0EtaValueUp[kstarBin] + sysErrorPLAPAL_V0EtaValueDown[kstarBin]);

    //Add them up quadratically:

    sysErrorTotalPLAPAL[kstarBin] = pow(sysErrorAverageV0Pt,2.) + pow(sysErrorAveragePt,2.) + pow(sysErrorAverageV0Eta,2.) + pow(sysErrorAverageEta,2.) + pow(sysErrorAverageNSigma,2.) + pow(sysErrorPLAPAL_FilterBit[kstarBin],2.) + pow(sysErrorPLAPAL_TPCcluster[kstarBin],2.)  + pow(sysErrorPLAPAL_V0CPAUp[kstarBin],2.) + pow(sysErrorPLAPAL_V0NsigmaDown[kstarBin],2.) + pow(sysErrorPLAPAL_V0TPCCluster80[kstarBin],2.) + pow(sysErrorPLAPAL_V0trackDistanceDown[kstarBin],2.) + pow(sysErrorPLAPAL_V0trackDistanctToPVUp[kstarBin],2.);

    sysErrorTotalPLAPAL[kstarBin] = TMath::Sqrt(sysErrorTotalPLAPAL[kstarBin]);

    //Put this information into arrays for plotting:
    xval_CFPLAPAL[kstarBin] = hCFPLAPAL[0]->GetBinCenter(kstarBin+1);
    yval_CFPLAPAL[kstarBin] = hCFPLAPAL[0]->GetBinContent(kstarBin+1);
    //xerr_CFPLAPAL[kstarBin] = 0.005;
    xerr_CFPLAPAL[kstarBin] = hCFPLAPAL[0]->GetBinWidth(1)/2.;

    sysErrorC2PLAPAL->SetBinContent(kstarBin+1,sysErrorTotalPLAPAL[kstarBin]);
    sysErrorC2PLAPAL->SetBinError(kstarBin+1,0.);

    hSystPL->Fill(hCFPLAPAL[0]->GetBinCenter(kstarBin+1), sysErrorTotalPLAPAL[kstarBin]);
  }

  for(int kstarBin = 0;kstarBin< sysErrorC2PLAPAL->GetNbinsX(); kstarBin++)
  {
    if(kstarBin < NkstarBinsPLAPAL)
    {
      sysErrorC2PLAPAL->SetBinContent(kstarBin+1,sysErrorTotalPLAPAL[kstarBin]);
      sysErrorC2PLAPAL->SetBinError(kstarBin+1,0.);
    }
    else
    {
      sysErrorC2PLAPAL->SetBinContent(kstarBin+1,0.);
      sysErrorC2PLAPAL->SetBinError(kstarBin+1,0.);
    }
  }

  double xval_CFLLALAL[NkstarBinsLLALAL];
  double yval_CFLLALAL[NkstarBinsLLALAL];
  double xerr_CFLLALAL[NkstarBinsLLALAL];
  
  TH1F *sysErrorC2LLALAL = (TH1F*)hCFAPAL[0]->Clone("C2totalsysPL");

  //Lets combine them together to a total systematic uncertainty for pp:
  for(int kstarBin = 0; kstarBin < NkstarBinsLLALAL; kstarBin++)
  {
    //calculate the average of two variations if applicable:
    //systematic error from V0 selection averaged:
    double sysErrorAverageV0Pt = 0.5*(sysErrorLLALAL_V0ptValueUp[kstarBin] + sysErrorLLALAL_V0ptValueDown[kstarBin]);
    double sysErrorAverageV0Eta = 0.5*(sysErrorLLALAL_V0EtaValueUp[kstarBin] + sysErrorLLALAL_V0EtaValueDown[kstarBin]);

    //Add them up quadratically:

    sysErrorTotalLLALAL[kstarBin] = pow(sysErrorAverageV0Pt,2.) + pow(sysErrorAverageV0Eta,2.) + pow(sysErrorLLALAL_V0CPAUp[kstarBin],2.) + pow(sysErrorLLALAL_V0NsigmaDown[kstarBin],2.) + pow(sysErrorLLALAL_V0TPCCluster80[kstarBin],2.) + pow(sysErrorLLALAL_V0trackDistanceDown[kstarBin],2.) + pow(sysErrorLLALAL_V0trackDistanctToPVUp[kstarBin],2.);

    sysErrorTotalLLALAL[kstarBin] = TMath::Sqrt(sysErrorTotalLLALAL[kstarBin]);

    //Put this information into arrays for plotting:
    xval_CFLLALAL[kstarBin] = hCFLLALAL[0]->GetBinCenter(kstarBin+1);
    yval_CFLLALAL[kstarBin] = hCFLLALAL[0]->GetBinContent(kstarBin+1);
    //xerr_CFLLALAL[kstarBin] = 0.005;
    xerr_CFLLALAL[kstarBin] = hCFLLALAL[0]->GetBinWidth(1)/2.;

    sysErrorC2LLALAL->SetBinContent(kstarBin+1,sysErrorTotalLLALAL[kstarBin]);
    sysErrorC2LLALAL->SetBinError(kstarBin+1,0.);

    hSystLL->Fill(hCFLLALAL[0]->GetBinCenter(kstarBin+1), sysErrorTotalLLALAL[kstarBin]);
  }

  for(int kstarBin = 0;kstarBin< sysErrorC2LLALAL->GetNbinsX(); kstarBin++)
  {
    if(kstarBin < NkstarBinsLLALAL)
    {
      sysErrorC2LLALAL->SetBinContent(kstarBin+1,sysErrorTotalLLALAL[kstarBin]);
      sysErrorC2LLALAL->SetBinError(kstarBin+1,0.);
    }
    else
    {
      sysErrorC2LLALAL->SetBinContent(kstarBin+1,0.);
      sysErrorC2LLALAL->SetBinError(kstarBin+1,0.);
    }
  }



  //********************************************************
  //********************************************************
  //Here starts the plotting stuff:
  //********************************************************
  //********************************************************
  
  TCanvas *CanSysErrorBudgetPPAPAP = new TCanvas("CanSysErrorBudgetPPAPAP","CanSysErrorBudgetPPAPAP",0,0,1500,1000);
  CanSysErrorBudgetPPAPAP->Divide(3,3);
  for(int protonVariation = 0; protonVariation < NTasksProtonVariations; protonVariation++)
  {
    CanSysErrorBudgetPPAPAP->cd(protonVariation+1);
    DrawHist(hCFPPAPAP_errBudget[protonVariation],"k* (GeV/#it{c})","|1 - #it{C}(k*)_{variation}/#it{C}(k*)_{default}|",1,1,1.);
    hCFPPAPAP_errBudget[protonVariation]->SetTitle(protonCut[protonVariation]);
  }
   if(rebinPP == 1) CanSysErrorBudgetPPAPAP->Print("PP_Budget.pdf");
  else CanSysErrorBudgetPPAPAP->Print(Form("PP_Budget_%i.pdf", rebinPP));

  TCanvas *CanSysErrorBudgetPLAPAL = new TCanvas("CanSysErrorBudgetPLAPAL","CanSysErrorBudgetPLAPAL",0,0,1500,1000);
  CanSysErrorBudgetPLAPAL->Divide(5,4);
  for(int V0Variation = 0; V0Variation < NTasksV0Variations; V0Variation++)
  {
    CanSysErrorBudgetPLAPAL->cd(V0Variation+1);
    DrawHist(hCFPLAPAL_errBudget[V0Variation],"k* (GeV/#it{c})","|1 - #it{C}(k*)_{variation}/#it{C}(k*)_{default}|",1,1,1.);
    hCFPLAPAL_errBudget[V0Variation]->SetTitle(v0Cut[V0Variation]);
  }
   if(rebinLP == 1) CanSysErrorBudgetPLAPAL->Print("PL_Budget.pdf");
  else CanSysErrorBudgetPLAPAL->Print(Form("PL_Budget_%i.pdf", rebinLP));

  TCanvas *CanSysErrorBudgetLLALAL = new TCanvas("CanSysErrorBudgetLLALAL","CanSysErrorBudgetLLALAL",0,0,1500,1000);
  CanSysErrorBudgetLLALAL->Divide(3,3);
  for(int V0Variation = NTasksProtonVariations; V0Variation < NTasksV0Variations; V0Variation++)
  {
    CanSysErrorBudgetLLALAL->cd(V0Variation-NTasksProtonVariations+1);
    DrawHist(hCFLLALAL_errBudget[V0Variation],"k* (GeV/#it{c})","|1 - #it{C}(k*)_{variation}/#it{C}(k*)_{default}|",1,1,1.);
    hCFLLALAL_errBudget[V0Variation]->SetTitle(v0Cut[V0Variation]);
  }
   if(rebinLL == 1) CanSysErrorBudgetLLALAL->Print("LL_Budget.pdf");
  else CanSysErrorBudgetLLALAL->Print(Form("LL_Budget_%i.pdf", rebinLL));
  
  TGraphErrors *graphError_CPPAPAP = new TGraphErrors(NkstarBinsPPAPAP,xval_CFPPAPAP,yval_CFPPAPAP,xerr_CFPPAPAP,sysErrorTotalPPAPAP);
  graphError_CPPAPAP->SetFillStyle(0);

  TGraphErrors *graphError_CPLAPAL = new TGraphErrors(NkstarBinsPLAPAL,xval_CFPLAPAL,yval_CFPLAPAL,xerr_CFPLAPAL,sysErrorTotalPLAPAL);
  graphError_CPLAPAL->SetFillStyle(0);

  TGraphErrors *graphError_CLLALAL = new TGraphErrors(NkstarBinsLLALAL,xval_CFLLALAL,yval_CFLLALAL,xerr_CFLLALAL,sysErrorTotalLLALAL);
  graphError_CLLALAL->SetFillStyle(0);
  
  TGraph *graphError_pp = new TGraph(NkstarBinsPPAPAP,xval_CFPPAPAP,sysErrorTotalPPAPAP);
  graphError_pp->SetName("systError_pp");
  TGraph *graphError_pL = new TGraph(NkstarBinsPLAPAL,xval_CFPLAPAL,sysErrorTotalPLAPAL);
  graphError_pL->SetName("systError_pL");
  TGraph *graphError_LL = new TGraph(NkstarBinsLLALAL,xval_CFLLALAL,sysErrorTotalLLALAL);
  graphError_LL->SetName("systError_LL");
  
  TFile *outfile = new TFile(Form("%sXi.root", output), "RECREATE");
  graphError_pp->Write();
  graphError_pL->Write();
  graphError_LL->Write();
//  hSystPP->Write("C2totalsysPP");
  hSystPL->Write("C2totalsysPXi");
//  hSystLL->Write("C2totalsysLL");
  CanCFPLAPALRatios->Write("CanCFPXiAPAXiRatios");
//  CanCFLLALALRatios->Write();
//  CanCFPPAPAPRatios->Write();
  
  // ratio total error to total value
  Double_t ratioPP[NkstarBinsPPAPAP];
  Double_t ratioPL[NkstarBinsPLAPAL];
  Double_t ratioLL[NkstarBinsPLAPAL];

  for(int kstarBin = 0; kstarBin < NkstarBinsPPAPAP; kstarBin++)
  {
    ratioPP[kstarBin] = sysErrorTotalPPAPAP[kstarBin]/yval_CFPPAPAP[kstarBin];
  }

  for(int kstarBin = 0; kstarBin < NkstarBinsPLAPAL; kstarBin++)
  {
    ratioPL[kstarBin] = sysErrorTotalPLAPAL[kstarBin]/yval_CFPLAPAL[kstarBin];
    ratioLL[kstarBin] = sysErrorTotalLLALAL[kstarBin]/yval_CFLLALAL[kstarBin];
  }

  TGraph *graphError_Ratio_pp = new TGraph(NkstarBinsPPAPAP,xval_CFPPAPAP,ratioPP);
  graphError_Ratio_pp->SetTitle("; k* (GeV/#it{c}); Ratio syst. error/data");
  TGraph *graphError_Ratio_pL = new TGraph(NkstarBinsPLAPAL,xval_CFPLAPAL,ratioPL);
  graphError_Ratio_pL->SetTitle("; k* (GeV/#it{c}); Ratio syst. error/data");
  TGraph *graphError_Ratio_LL = new TGraph(NkstarBinsLLALAL,xval_CFLLALAL,ratioLL);
  graphError_Ratio_LL->SetTitle("; k* (GeV/#it{c}); Ratio syst. error/data");

  TF1 *fRatioPP = new TF1("fRatioPP", "pol2", 0, 10);
  TF1 *fRatioPL = new TF1("fRatioPL", "pol2", 0, 10);
  TF1 *fRatioLL = new TF1("fRatioLL", "pol2", 0, 10);

  gStyle->SetTitleFontSize(0.05);
  TCanvas *CanSysErrorPlotting = new TCanvas("CanSysErrorPlotting","CanSysErrorPlotting",0,0,1500,1000);
  CanSysErrorPlotting->Divide(3,2);
  CanSysErrorPlotting->cd(1);
  DrawHist(hCFPPAPAP[0],"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.);
  hCFPPAPAP[0]->SetTitle("p-p #oplus #bar{p}#bar{p}");
  hCFPPAPAP[0]->GetXaxis()->SetRangeUser(0, 0.15);
  //hCFPPAPAP[0]->Draw();
  graphError_CPPAPAP->Draw("2same");
  CanSysErrorPlotting->cd(2);
  DrawHist(hCFPLAPAL[0],"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.);
  hCFPLAPAL[0]->SetTitle("p-#Lambda #oplus #bar{p}#bar{#Lambda}");
  //hCFPLAPAL[0]->Draw();
  graphError_CPLAPAL->Draw("2same");
  CanSysErrorPlotting->cd(3);
  DrawHist(hCFLLALAL[0],"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.);
  hCFLLALAL[0]->SetTitle("#Lambda-#Lambda #oplus #bar{#Lambda}#bar{#Lambda}");
  graphError_CLLALAL->Draw("2same");
  CanSysErrorPlotting->cd(4);
  StyleGraph(graphError_Ratio_pp);
  graphError_Ratio_pp->Fit("fRatioPP", "R", "", 0.01, 0.2);
  graphError_Ratio_pp->Draw("ape");
  graphError_Ratio_pp->SetMarkerStyle(20);
  graphError_Ratio_pp->GetXaxis()->SetRangeUser(0,0.16);
  graphError_Ratio_pp->GetYaxis()->SetRangeUser(0,0.05);
  graphError_Ratio_pp->GetYaxis()->SetTitleOffset(1.1);
  CanSysErrorPlotting->cd(5);
  StyleGraph(graphError_Ratio_pL);
  graphError_Ratio_pL->Fit("fRatioPL", "R", "", 0, 1);
  graphError_Ratio_pL->Draw("ape");
  graphError_Ratio_pL->SetMarkerStyle(20);
  graphError_Ratio_pL->GetXaxis()->SetRangeUser(0,0.5);
  graphError_Ratio_pL->GetYaxis()->SetRangeUser(0,0.02);
  graphError_Ratio_pL->GetYaxis()->SetTitleOffset(1.1);
  CanSysErrorPlotting->cd(6);
  StyleGraph(graphError_Ratio_LL);
  graphError_Ratio_LL->Fit("fRatioLL", "R", "", 0, 0.5);
  graphError_Ratio_LL->Draw("ape");
  graphError_Ratio_LL->SetMarkerStyle(20);
  graphError_Ratio_LL->GetXaxis()->SetRangeUser(0,0.5);
  graphError_Ratio_LL->GetYaxis()->SetRangeUser(0,0.05);
  graphError_Ratio_LL->GetYaxis()->SetTitleOffset(1.1);
  if(rebinPP == 1 && rebinLP == 1 && rebinLL == 1) CanSysErrorPlotting->Print("CF_sys.pdf");
  else CanSysErrorPlotting->Print(Form("CF_sys_%i_%i_%i.pdf", rebinPP, rebinLP, rebinLL));

  // NOW COMPUTE THE ERRORS

  TGraphErrors *grFinalErrorPP = new TGraphErrors();
  const float xerrPP = hCFPPAPAPorig->GetBinWidth(1)/2.;
  for(int kstarBin = 0; kstarBin < hCFPPAPAPorig->GetXaxis()->FindBin(fitRangeRightPPAPAP); kstarBin++)
  {
    const float x = hCFPPAPAPorig->GetBinCenter(kstarBin+1);
    const float y = hCFPPAPAPorig->GetBinContent(kstarBin+1);
    grFinalErrorPP->SetPoint(kstarBin, x, y);
    grFinalErrorPP->SetPointError(kstarBin, xerrPP,  y * fRatioPP->Eval(x));
    hSystErrorPP->Fill(x, y * fRatioPP->Eval(x));
  }
  grFinalErrorPP->SetFillColor(kGray+1);
  grFinalErrorPP->SetMarkerStyle(20);
  grFinalErrorPP->SetLineWidth(3);

  TGraphErrors *grFinalErrorPL = new TGraphErrors();
  const float xerrPL = hCFPLAPALorig->GetBinWidth(1)/2.;
  for(int kstarBin = 0; kstarBin < hCFPLAPALorig->GetXaxis()->FindBin(fitRangeRightPPAPAP); kstarBin++)
  {
    const float x = hCFPLAPALorig->GetBinCenter(kstarBin+1);
    const float y = hCFPLAPALorig->GetBinContent(kstarBin+1);
    grFinalErrorPL->SetPoint(kstarBin, x, y);
    grFinalErrorPL->SetPointError(kstarBin, xerrPL,  y * fRatioPL->Eval(x));
    hSystErrorPL->Fill(x, y * fRatioPL->Eval(x));
  }
  grFinalErrorPL->SetFillColor(kGray+1);
  grFinalErrorPL->SetMarkerStyle(20);
  grFinalErrorPL->SetLineWidth(3);

  TGraphErrors *grFinalErrorLL = new TGraphErrors();
  const float xerrLL = hCFLLALALorig->GetBinWidth(1)/2.;
  for(int kstarBin = 0; kstarBin < hCFLLALALorig->GetXaxis()->FindBin(fitRangeRightPPAPAP); kstarBin++)
  {
    const float x = hCFLLALALorig->GetBinCenter(kstarBin+1);
    const float y = hCFLLALALorig->GetBinContent(kstarBin+1);
    grFinalErrorLL->SetPoint(kstarBin, x, y);
    grFinalErrorLL->SetPointError(kstarBin, xerrLL,  y * fRatioLL->Eval(x));
    hSystErrorLL->Fill(x, y * fRatioLL->Eval(x));
  }
  grFinalErrorLL->SetFillColor(kGray+1);
  grFinalErrorLL->SetMarkerStyle(20);
  grFinalErrorLL->SetLineWidth(3);


  TCanvas *CanPlotting = new TCanvas("CanPlotting","CanPlotting",0,0,1500,500);
  CanPlotting->Divide(3,1);
  CanPlotting->cd(1);
  DrawHist(hCFPPAPAPorig,"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.);
  hCFPPAPAPorig->SetTitle("p-p #oplus #bar{p}#bar{p}");
  hCFPPAPAPorig->GetXaxis()->SetRangeUser(0,0.12);
  hCFPPAPAPorig->GetYaxis()->SetRangeUser(0,4);
  grFinalErrorPP->Draw("2 same");
  hCFPPAPAPorig->Draw("same");
  CanPlotting->cd(2);
  DrawHist(hCFPLAPALorig,"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.);
  hCFPLAPALorig->SetTitle("p-#Lambda #oplus #bar{p}#bar{#Lambda}");
  hCFPLAPALorig->GetXaxis()->SetRangeUser(0,0.25);
  hCFPLAPALorig->GetYaxis()->SetRangeUser(0.5,2.1);
  grFinalErrorPL->Draw("2 same");
  hCFPLAPALorig->Draw("same");
  CanPlotting->cd(3);
  DrawHist(hCFLLALALorig,"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.);
  hCFLLALALorig->SetTitle("#Lambda-#Lambda #oplus #bar{#Lambda}#bar{#Lambda}");
  hCFLLALALorig->GetXaxis()->SetRangeUser(0,0.25);
  hCFLLALALorig->GetYaxis()->SetRangeUser(0.5,2.5);
  grFinalErrorLL->Draw("2 same");
  hCFLLALALorig->Draw("same");
  CanPlotting->Print(Form("CF_systnew_%i_%i_%i.pdf", rebinPP, rebinLP, rebinLL));

//  TFile *ppFile = new TFile("C2totalsysPP.root", "RECREATE");
//  hSystErrorPP->Write("C2totalsysPP");

  TFile *plFile = new TFile("C2totalsysPXi.root", "RECREATE");
  hSystErrorPL->Write("C2totalsysPXi");

//  TFile *llFile = new TFile("C2totalsysLL.root", "RECREATE");
//  hSystErrorLL->Write("C2totalsysLL");

  return 0;
}
