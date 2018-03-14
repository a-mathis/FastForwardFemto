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
  TH1F *hSystErrorPL = new TH1F("hSystErrorPXi", "", 150, 0, 3);
  TH1F *hSystErrorLL = new TH1F("hSystErrorLL", "", 150, 0, 3);


  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);
  
  TFile* _file0=TFile::Open(file,"READ");
  if(!_file0) return 0;

  const Int_t NTasks = 22;//Number of different total variations (including default)
  const Int_t NTasksProtonVariations = 9; //Number of different proton variations
  const Int_t NTasksV0Variations = 12 + NTasksProtonVariations; //Number of different V0 variations (but proton selection does also affect p-V0 correlation function
  
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
    else if(taskNumber == 10) outputfile = "/PWGCF_PLFemto_19";//DCA Xi Daughters
    else if(taskNumber == 11) outputfile = "/PWGCF_PLFemto_20";//Min Dist Bach to PV
    else if(taskNumber == 12) outputfile = "/PWGCF_PLFemto_21";//CPA Xi
    else if(taskNumber == 13) outputfile = "/PWGCF_PLFemto_22";//Transverse Radius Xi
    else if(taskNumber == 14) outputfile = "/PWGCF_PLFemto_23";//DCA v0 Daughters
    else if(taskNumber == 15) outputfile = "/PWGCF_PLFemto_24";//CPA v0
    else if(taskNumber == 16) outputfile = "/PWGCF_PLFemto_25";//Transverse Radius v0
    else if(taskNumber == 17) outputfile = "/PWGCF_PLFemto_26";//Min Dist v0 to PV
    else if(taskNumber == 18) outputfile = "/PWGCF_PLFemto_27";//Min Dist v0 Daughters to PV
    else if(taskNumber == 19) outputfile = "/PWGCF_PLFemto_28";//Xi Daughters Eta down
    else if(taskNumber == 20) outputfile = "/PWGCF_PLFemto_29";//Xi Daughrers Eta up
    else if(taskNumber == 21) outputfile = "/PWGCF_PLFemto_30";//nSigma Daughters down

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

    hREPL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fProtonXiRelK");
    hMEPL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fProtonXiRelKME");
    hREPL[taskNumber]->Rebin(rebinLP);
    hMEPL[taskNumber]->Rebin(rebinLP);
    HistName = "hCFPL_task";
    HistName += taskNumber;
    hCFPL[taskNumber] = Calculate_CF(hREPL[taskNumber],hMEPL[taskNumber],HistName.Data(),normLeft,normRight);
    if(!hCFPL[taskNumber]) std::cout << "PL with task number " << taskNumber << " does not exist" << std::endl;

    hREAPAL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fAntiProtonAntiXiRelK");
    hMEAPAL[taskNumber] = (TH1F*)listTP[taskNumber]->FindObject("fAntiProtonAntiXiRelKME");
    hREAPAL[taskNumber]->Rebin(rebinLP);
    hMEAPAL[taskNumber]->Rebin(rebinLP);
    HistName = "hCFAPAL_task";
    HistName += taskNumber;
    hCFAPAL[taskNumber] = Calculate_CF(hREAPAL[taskNumber],hMEAPAL[taskNumber],HistName.Data(),normLeft,normRight);
    if(!hCFAPAL[taskNumber]) std::cout << "APAL with task number " << taskNumber << " does not exist" << std::endl;

    HistName = "hCFPLAPAL_task";
    HistName += taskNumber;
    hCFPLAPAL[taskNumber] = add_CF(hCFPL[taskNumber],hCFAPAL[taskNumber],HistName.Data());
    if(!hCFPLAPAL[taskNumber]) std::cout << "Sum of PL and APAL with task number " << taskNumber << "does not exist" << std::endl;
  }

  //Also the sum, in which we are finally interested
  TCanvas *CanCFsummarySum = new TCanvas("CanCFsummarySum","CanCFsummarySum",0,0,1500,500);

  for(int taskNumber = 0;taskNumber < NTasks; taskNumber++)
  {

    hCFPLAPAL[taskNumber]->GetXaxis()->SetRangeUser(0,0.5);
    hCFPLAPAL[taskNumber]->SetLineColor(taskNumber+1);
    hCFPLAPAL[taskNumber]->SetMarkerColor(taskNumber+1);
    if(taskNumber == 0) {hCFPLAPAL[taskNumber]->SetMarkerStyle(23); hCFPLAPAL[taskNumber]->SetMarkerSize(1.3);hCFPLAPAL[taskNumber]->SetTitle("p#Lambda #oplus #bar{p}#bar{#Lambda}"); hCFPLAPAL[taskNumber]->Draw();}
    else hCFPLAPAL[taskNumber]->Draw("same");

  }
   if(rebinPP == 1 && rebinLP == 1 && rebinLL == 1) CanCFsummarySum->Print("CF.pdf");
   else CanCFsummarySum->Print(Form("CF_%i_%i_%i.pdf", rebinPP, rebinLP, rebinLL));


  //Lets build all ratios:
  
  TH1F *hCFRatioPLAPAL[NTasksV0Variations];
  

  //All ratios of default to variations C_default / C_variations for pL
  for(int V0Variation = 0; V0Variation < NTasksV0Variations; V0Variation++)
  {
    //First clone the default correlation function:
    HistName = "hCFRatioPLAPAL_task_0_";
    HistName += V0Variation+1;

    hCFRatioPLAPAL[V0Variation] = (TH1F*)hCFPLAPAL[0]->Clone(HistName.Data());
    hCFRatioPLAPAL[V0Variation]->Divide(hCFPLAPAL[V0Variation+1]);
  }


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
  v0Cut[9] = "#Xi d_{daug} down";                     //DCA Xi Daughters
  v0Cut[10] = "#Xi d_{bach,PV} up";                   //Min Dist Bach to PV
  v0Cut[11] = "#Xi CPA up";                           //CPA Xi
  v0Cut[12] = "#Xi Min. Transverse Radius Up";        //Transverse Radius Xi
  v0Cut[13] = "#Lambda d_{daug} down";                //DCA v0 Daughters
  v0Cut[14] = "#Lambda CPA up";                       //CPA v0
  v0Cut[15] = "#Lambda Min. Transverse Radius Up";    //Transverse Radius v0
  v0Cut[16] = "#Xi d_{#Lambda,PV} up";                //Min Dist v0 to PV
  v0Cut[17] = "#Lambda d_{#Lambda,PV} up";            //Min Dist v0 Daughters to PV
  v0Cut[18] = "V0 #eta down";                         //Xi Daughters Eta down
  v0Cut[19] = "V0 #eta up";                           //Xi Daughrers Eta up
  v0Cut[20] = "V0 n#sigma down";                      //nSigma Daughters down

  const int NkstarBinsPLAPAL = hCFRatioPLAPAL[0]->GetXaxis()->FindBin(fitRangeRightPLAPAL);
  //systematic error due to proton selection:
  double sysErrorPLAPAL_ptValueDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_ptValueUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_EtaValueUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_EtaValueDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_NSigmaValueUp[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_NSigmaValueDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_FilterBit[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_TPCclusterDown[NkstarBinsPLAPAL];
  double sysErrorPLAPAL_TPCclusterUp[NkstarBinsPLAPAL];
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
    sysErrorPLAPAL_TPCclusterDown[kstarBin] = 0;
    sysErrorPLAPAL_TPCclusterUp[kstarBin] = 0;
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
  CanCFPLAPALRatios->Divide(5,5);
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
        sysErrorPLAPAL_TPCclusterDown[kstarBin] = 0.;
        sysErrorPLAPAL_TPCclusterDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 8)//proton TPC cluster
      {
        sysErrorPLAPAL_TPCclusterUp[kstarBin] = 0.;
        sysErrorPLAPAL_TPCclusterUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 9)//V0 Pt down
      {
        sysErrorPLAPAL_V0ptValueDown[kstarBin] = 0.;
        sysErrorPLAPAL_V0ptValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 10)//V0 Pt up
      {
        sysErrorPLAPAL_V0ptValueUp[kstarBin] = 0.;
        sysErrorPLAPAL_V0ptValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 11)//V0 CPA up
      {
        sysErrorPLAPAL_V0CPAUp[kstarBin] = 0.;
        sysErrorPLAPAL_V0CPAUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 12)//V0 NSigma down
      {
        sysErrorPLAPAL_V0NsigmaDown[kstarBin] = 0.;
        sysErrorPLAPAL_V0NsigmaDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 13)//V0 TPC cluster up
      {
        sysErrorPLAPAL_V0TPCCluster80[kstarBin] = 0.;
        sysErrorPLAPAL_V0TPCCluster80[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 14)//V0 eta up
      {
        sysErrorPLAPAL_V0EtaValueUp[kstarBin] = 0.;
        sysErrorPLAPAL_V0EtaValueUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 15)//V0 eta down
      {
        sysErrorPLAPAL_V0EtaValueDown[kstarBin] = 0.;
        sysErrorPLAPAL_V0EtaValueDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 16)//V0 track distance down
      {
        sysErrorPLAPAL_V0trackDistanceDown[kstarBin] = 0.;
        sysErrorPLAPAL_V0trackDistanceDown[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 17)//V0 impact parameter up
      {
        sysErrorPLAPAL_V0trackDistanctToPVUp[kstarBin] = 0.;
        sysErrorPLAPAL_V0trackDistanctToPVUp[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 18)//V0 impact parameter up
      {
        sysErrorPLAPAL_V0trackDistanctToPVUp10[kstarBin] = 0.;
        sysErrorPLAPAL_V0trackDistanctToPVUp10[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 19)//V0 impact parameter up
      {
        sysErrorPLAPAL_V0trackDistanctToPVUp11[kstarBin] = 0.;
        sysErrorPLAPAL_V0trackDistanctToPVUp11[kstarBin] = fabs(1. - pow(SysError,-1)) * hCFPLAPAL[0]->GetBinContent(kstarBin+1);
        hCFPLAPAL_errBudget[V0Variation]->SetBinContent(kstarBin+1,fabs(1. - pow(SysError,-1)));
      }
      else if(V0Variation == 20)//V0 impact parameter up
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

  double xval_CFPLAPAL[NkstarBinsPLAPAL];
  double yval_CFPLAPAL[NkstarBinsPLAPAL];
  double xerr_CFPLAPAL[NkstarBinsPLAPAL];
  
  TH1F *sysErrorC2PLAPAL = (TH1F*)hCFAPAL[0]->Clone("C2totalsysPXi");

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

    sysErrorTotalPLAPAL[kstarBin] = pow(sysErrorAverageV0Pt,2.) + pow(sysErrorAveragePt,2.) +
        pow(sysErrorAverageV0Eta,2.) + pow(sysErrorAverageEta,2.) + pow(sysErrorAverageNSigma,2.) +
        pow(sysErrorPLAPAL_FilterBit[kstarBin],2.) + pow(sysErrorPLAPAL_TPCclusterDown[kstarBin],2.) + pow(sysErrorPLAPAL_TPCclusterUp[kstarBin],2.)  + pow(sysErrorPLAPAL_V0CPAUp[kstarBin],2.) + pow(sysErrorPLAPAL_V0NsigmaDown[kstarBin],2.) + pow(sysErrorPLAPAL_V0TPCCluster80[kstarBin],2.) + pow(sysErrorPLAPAL_V0trackDistanceDown[kstarBin],2.) + pow(sysErrorPLAPAL_V0trackDistanctToPVUp[kstarBin],2.);

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

  TCanvas *CanSysErrorBudgetPLAPAL = new TCanvas("CanSysErrorBudgetPLAPAL","CanSysErrorBudgetPLAPAL",0,0,1500,1000);
  CanSysErrorBudgetPLAPAL->Divide(5,5);
  for(int V0Variation = 0; V0Variation < NTasksV0Variations; V0Variation++)
  {
    CanSysErrorBudgetPLAPAL->cd(V0Variation+1);
    DrawHist(hCFPLAPAL_errBudget[V0Variation],"k* (GeV/#it{c})","|1 - #it{C}(k*)_{variation}/#it{C}(k*)_{default}|",1,1,1.);
    hCFPLAPAL_errBudget[V0Variation]->SetTitle(v0Cut[V0Variation]);
    TF1 *f1 = new TF1("f1","pol0",0,1);
    hCFPLAPAL_errBudget[V0Variation]->Fit("f1","R","",0,0.2);
    f1->Draw("SAME");
    TLatex LambdaLabel;
    LambdaLabel.SetNDC(kTRUE);
    LambdaLabel.SetTextSize(gStyle->GetTextSize()*1.2);
    LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,Form("Par Fit: %.2f",f1->GetParameter(0)*100));
  }
   if(rebinLP == 1) CanSysErrorBudgetPLAPAL->Print("PL_Budget.pdf");
  else CanSysErrorBudgetPLAPAL->Print(Form("PL_Budget_%i.pdf", rebinLP));

  TGraphErrors *graphError_CPLAPAL = new TGraphErrors(NkstarBinsPLAPAL,xval_CFPLAPAL,yval_CFPLAPAL,xerr_CFPLAPAL,sysErrorTotalPLAPAL);
  graphError_CPLAPAL->SetFillStyle(0);

  TGraph *graphError_pL = new TGraph(NkstarBinsPLAPAL,xval_CFPLAPAL,sysErrorTotalPLAPAL);
  graphError_pL->SetName("systError_pL");
  
  TFile *outfile = new TFile(Form("%sXi.root", output), "RECREATE");
  graphError_pL->Write();
  hSystPL->Write("C2totalsysPXi");
  CanCFPLAPALRatios->Write("CanCFPXiAPAXiRatios");
  
  // ratio total error to total value
  Double_t ratioPL[NkstarBinsPLAPAL];


  for(int kstarBin = 0; kstarBin < NkstarBinsPLAPAL; kstarBin++)
  {
    ratioPL[kstarBin] = sysErrorTotalPLAPAL[kstarBin]/yval_CFPLAPAL[kstarBin];
  }

  TGraph *graphError_Ratio_pL = new TGraph(NkstarBinsPLAPAL,xval_CFPLAPAL,ratioPL);
  graphError_Ratio_pL->SetTitle("; k* (GeV/#it{c}); Ratio syst. error/data");

  TF1 *fRatioPL = new TF1("fRatioPL", "pol2", 0, 10);

  gStyle->SetTitleFontSize(0.05);
  TCanvas *CanSysErrorPlotting = new TCanvas("CanSysErrorPlotting","CanSysErrorPlotting",0,0,1500,1000);
  CanSysErrorPlotting->Divide(2,1);
  CanSysErrorPlotting->cd(1);
  DrawHist(hCFPLAPAL[0],"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.);
  hCFPLAPAL[0]->SetTitle("p-#Lambda #oplus #bar{p}#bar{#Lambda}");
  graphError_CPLAPAL->Draw("2same");
  CanSysErrorPlotting->cd(2);
  StyleGraph(graphError_Ratio_pL);
  graphError_Ratio_pL->Fit("fRatioPL", "R", "", 0, 1);
  graphError_Ratio_pL->Draw("ape");
  graphError_Ratio_pL->SetMarkerStyle(20);
  graphError_Ratio_pL->GetXaxis()->SetRangeUser(0,0.5);
  graphError_Ratio_pL->GetYaxis()->SetRangeUser(0,0.02);
  graphError_Ratio_pL->GetYaxis()->SetTitleOffset(1.1);
  if(rebinPP == 1 && rebinLP == 1 && rebinLL == 1) CanSysErrorPlotting->Print("CF_sys.pdf");
  else CanSysErrorPlotting->Print(Form("CF_sys_%i_%i_%i.pdf", rebinPP, rebinLP, rebinLL));

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

  TCanvas *CanPlotting = new TCanvas("CanPlotting","CanPlotting",0,0,1500,500);
  DrawHist(hCFPLAPALorig,"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.);
  hCFPLAPALorig->SetTitle("p-#Lambda #oplus #bar{p}#bar{#Lambda}");
  hCFPLAPALorig->GetXaxis()->SetRangeUser(0,0.25);
  hCFPLAPALorig->GetYaxis()->SetRangeUser(0.5,2.1);
  grFinalErrorPL->Draw("2 same");
  hCFPLAPALorig->Draw("same");
  CanPlotting->Print(Form("CF_systnew_%i_%i_%i.pdf", rebinPP, rebinLP, rebinLL));


  TFile *plFile = new TFile("C2totalsysPXi.root", "RECREATE");
  hSystErrorPL->Write("C2totalsysPXi");


  return 0;
}
