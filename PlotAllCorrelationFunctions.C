
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
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
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
void PlotAllCorrelationFunctions(const char *filename) {
  SetStyle(false,false);
  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;

  TH1F* SEDist[6][6];
  TH1F* MEDist[6][6];
  TH1F* CF[6][6];
  TFile *file=TFile::Open(filename);
  TString RelKNames[6]={"Proton","AntiProton","Lambda","AntiLambda","Xi","AntiXi"};
  TDirectoryFile *dirResults=(TDirectoryFile*)(file->FindObjectAny("Results"));
  TCanvas *c1[6];
  if (dirResults) {
    TList *tmp=dirResults->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *Results;
    dirResults->GetObject(name,Results);
    if (!Results) {
      std::cout << "No Results \n";
    } else {
      for (int iPart1=0;iPart1<6;++iPart1) {
        c1[iPart1]= new TCanvas(Form("CF_%s",RelKNames[iPart1].Data()),Form("CF_%s",RelKNames[iPart1].Data()),0,0,1000,1100);
//        c1[iPart1]->Divide(6-iPart1,1);
        if (iPart1 == 0 || iPart1 == 1) {
          c1[iPart1]->Divide(3,2);
        } else if (iPart1 == 2) {
          c1[iPart1]->Divide(2,2);
        } else if (iPart1 == 3 ) {
          c1[iPart1]->Divide(3,1);
        } else if (iPart1 == 4 ) {
          c1[iPart1]->Divide(2,1);
        }
        for (int iPart2=iPart1;iPart2<6;++iPart2) {
          c1[iPart1]->cd(iPart2-iPart1+1);
          TString folderName=Form("Particle%i_Particle%i",iPart1,iPart2);
          TList* tmpFolder=(TList*)Results->FindObject(folderName.Data());
          if (!tmpFolder) {
            std::cout << folderName.Data() <<" not Found\n";
          } else {
            TString SEName=Form("SEDist_Particle%i_Particle%i",iPart1,iPart2);
            SEDist[iPart1][iPart2]=(TH1F*)tmpFolder->FindObject(SEName.Data());
            SEDist[iPart1][iPart2]->SetDirectory(0);
            if (!SEDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            TString NewSEName=Form("f%s%sRelK",RelKNames[iPart1].Data(),RelKNames[iPart2].Data());
            SEDist[iPart1][iPart2]->SetNameTitle(NewSEName.Data(),NewSEName.Data());

            TString MEName=Form("MEDist_Particle%i_Particle%i",iPart1,iPart2);
            MEDist[iPart1][iPart2]=(TH1F*)tmpFolder->FindObject(MEName.Data());
            MEDist[iPart1][iPart2]->SetDirectory(0);
            if (!MEDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            TString NewMEName=Form("f%s%sRelKME",RelKNames[iPart1].Data(),RelKNames[iPart2].Data());
            MEDist[iPart1][iPart2]->SetNameTitle(NewMEName.Data(),NewMEName.Data());
            CF[iPart1][iPart2]=
                Calculate_CF(SEDist[iPart1][iPart2],MEDist[iPart1][iPart2],
                             Form("CF_%s_%s",RelKNames[iPart1].Data(),
                                  RelKNames[iPart2].Data()),normleft,normright);
            SetStyleHisto(CF[iPart1][iPart2],8,fColors[0]);
            TString setNameStuff=Form("%s-%s; k* (GeV/#it{c}); #it{C}(k*)",RelKNames[iPart1].Data(),RelKNames[iPart2].Data());
            CF[iPart1][iPart2]->SetTitle(setNameStuff.Data());
            CF[iPart1][iPart2]->GetXaxis()->SetRangeUser(0,0.4);
            CF[iPart1][iPart2]->DrawCopy();
            TLatex BeamText;
            BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
            BeamText.SetNDC(kTRUE);
            BeamText.DrawLatex(0.5, 0.7, Form("%s-%s",RelKNames[iPart1].Data(),RelKNames[iPart2].Data()));
          }
        }
      }
    }
  }
}
