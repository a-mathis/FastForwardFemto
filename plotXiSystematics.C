#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"

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
  histo->GetXaxis()->SetTitleOffset(1.05);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.05);
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

void plotXiSystematics(const char *expfile) {
  // EXP DATA

  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;

  TFile* _file0=TFile::Open(expfile);
  //  TCanvas *XiCanvas=new TCanvas("XiCan","XiCan",0,0,1500,1100);
  //  XiCanvas->Divide(1,2);
  //  XiCanvas->cd(1);
  //  hist_CF_pXi_ApAXi_exp[0]->DrawCopy();
  //  hist_CF_pXi_ApAXi_exp[1]->DrawCopy("SAME");
  //  XiCanvas->cd(2);
  //  hist_CF_pXi_ApAXi_exp[2]->DrawCopy();
  TFile *errFile=new TFile("RelErr.root","RECREATE");
//  _file0->ls();
  int iVar=0;
  TH1F* RelError[22];
  TCanvas *relErr = new TCanvas("relErr","relErr",0,0,1500,1100);
  TCanvas *XiCorr= new TCanvas("XiCorr","XiCorr",0,0,1500,1100);
  XiCorr->Divide(2,1);
  for (int i = 0; i<=30;++i) {
    if (i==0||i>17) {
      TString inputDirName=Form("PWGCF_PLFemto_%i",i);
      TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny(inputDirName.Data()));
      if (!dirResults) {
        std::cout << "No Directory: " << inputDirName.Data() << std::endl;
      } else {
        TList *tmp=dirResults->GetListOfKeys();
        TString name=tmp->At(0)->GetName();
        TList *listTP;
        dirResults->GetObject(name,listTP);
        if (listTP) {

          TH1F* histRE_relK_Xip = (TH1F*)listTP->FindObject("fProtonXiRelK");
          TH1F* histME_relK_Xip = (TH1F*)listTP->FindObject("fProtonXiRelKME");
          TH1F* histRE_relK_AXiAp = (TH1F*)listTP->FindObject("fAntiProtonAntiXiRelK");
          TH1F* histME_relK_AXiAp = (TH1F*)listTP->FindObject("fAntiProtonAntiXiRelKME");

          TH1F *hist_CF_pXi_ApAXi_exp[3];
          hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(histRE_relK_Xip,histME_relK_Xip,"hist_CF_pp",normleft,normright);
          hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(histRE_relK_AXiAp,histME_relK_AXiAp,"hist_CF_ApAp",normleft,normright);
          hist_CF_pXi_ApAXi_exp[2] = add_CF(hist_CF_pXi_ApAXi_exp[0],hist_CF_pXi_ApAXi_exp[1],"hist_CF_pp_ApAp_exp_sum");
          hist_CF_pXi_ApAXi_exp[0]->SetStats(0);
          hist_CF_pXi_ApAXi_exp[0]->SetTitle("Comparison p#Xi and #bar{p#Xi}");
          hist_CF_pXi_ApAXi_exp[0]->GetXaxis()->SetTitle("k* (GeV/#it{c})");
          hist_CF_pXi_ApAXi_exp[0]->GetXaxis()->SetRangeUser(0.,0.4);
          hist_CF_pXi_ApAXi_exp[0]->GetYaxis()->SetTitle("#it{C}(k*)");
          hist_CF_pXi_ApAXi_exp[1]->SetStats(0);
          hist_CF_pXi_ApAXi_exp[1]->SetTitle("");
          hist_CF_pXi_ApAXi_exp[1]->GetXaxis()->SetTitle("k* (GeV/#it{c})");
          hist_CF_pXi_ApAXi_exp[1]->GetXaxis()->SetRangeUser(0.,0.4);
          hist_CF_pXi_ApAXi_exp[1]->GetYaxis()->SetTitle("#it{C}(k*)");
          hist_CF_pXi_ApAXi_exp[2]->SetTitle("Combined p#Xi and #bar{p#Xi}");
          hist_CF_pXi_ApAXi_exp[2]->SetStats(0);
          hist_CF_pXi_ApAXi_exp[2]->GetXaxis()->SetTitle("k* (GeV/#it{c})");
          hist_CF_pXi_ApAXi_exp[2]->GetXaxis()->SetRangeUser(0.,0.4);
          hist_CF_pXi_ApAXi_exp[2]->GetYaxis()->SetTitle("#it{C}(k*)");
          SetStyleHisto(hist_CF_pXi_ApAXi_exp[0], 1,1);
          SetStyleHisto(hist_CF_pXi_ApAXi_exp[1], 0,2);
          SetStyleHisto(hist_CF_pXi_ApAXi_exp[2], 0,0);
          if (i==0){
            XiCorr->cd(1);
            hist_CF_pXi_ApAXi_exp[0]->DrawCopy();
            hist_CF_pXi_ApAXi_exp[1]->DrawCopy("SAME");
            XiCorr->cd(2);
            hist_CF_pXi_ApAXi_exp[2]->DrawCopy();
          }


          TString RelErrName=Form("RelErr_%i",i);
          int maxBins=hist_CF_pXi_ApAXi_exp[2]->GetXaxis()->GetNbins();
          int nBins=maxBins*30;
          RelError[iVar]=new TH1F(RelErrName.Data(),RelErrName.Data(),nBins,0,3);
          for (int iBins=1;iBins<maxBins;++iBins) {
            float RelErrVal=hist_CF_pXi_ApAXi_exp[2]->GetBinError(iBins)/hist_CF_pXi_ApAXi_exp[2]->GetBinContent(iBins);
            RelError[iVar]->SetBinContent(30*iBins+iVar,RelErrVal);
          }
          RelError[iVar]->SetStats(0);
          RelError[iVar]->SetTitle("");
          RelError[iVar]->GetXaxis()->SetTitle("k* (GeV/#it{c})");
          RelError[iVar]->GetYaxis()->SetTitle("Rel. Err #sigma/mean");
          SetStyleHisto(RelError[iVar],1,iVar%8);
          RelError[iVar]->GetYaxis()->SetRangeUser(0.,1.2);
          RelError[iVar]->GetXaxis()->SetRangeUser(0.,0.25);
          relErr->cd();
          if(iVar==0)RelError[iVar]->DrawCopy();
          else RelError[iVar]->DrawCopy("SAME");
          RelError[iVar]->Write(RelErrName);
          if (i==0) {
            iVar+=8;
          }
          iVar++;
        }
      }
    }
  }
  relErr->Write("ErrorCanvas");
}
