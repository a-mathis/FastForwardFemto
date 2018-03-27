
std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7};
std::vector<int> fColors     = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2};
std::vector<int> fMarkers    = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette=false, bool title=false)
{
  const int NCont = 25;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetMarkerSize(1);
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
  gStyle->SetPalette(kRainBow);
  //  const int NRGBs = 6;
  //  Double_t stops[NRGBs];
  //  for(int i=0; i<NRGBs; ++i) stops[i] = float(i)/(NRGBs-1);
  //
  //  Double_t red[NRGBs]   = { 1.,  29./255., 25./255., 27./255., 32./255., 24./255.};
  //  Double_t green[NRGBs] = { 1., 221./255., 160./255., 113./255., 74./255., 37./255.};
  //  Double_t blue[NRGBs] = {  1., 221./255., 184./255., 154./255., 129./255., 98./255.};
  //  TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *DrawSystematicError(TH1F* histexp,TH1F *histerr,double errorwidth)
{
  //Input are the experimental histogram with statistical errors only and a histogram containing the errors

  const int histbins = histexp->GetNbinsX();

  TGraphErrors *ge_SysError_C2 = new TGraphErrors();

  for(int i=0;i<histbins;i++)
  {
    if(histexp->GetBinCenter(i+1) > 0.2) continue;
    ge_SysError_C2->SetPoint(i, histexp->GetBinCenter(i+1), histexp->GetBinContent(i+1));
    ge_SysError_C2->SetPointError(i, errorwidth, histerr->GetBinContent(i+1));
  }

  return ge_SysError_C2;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *FemtoModelFitBands(TGraph *grMedian1, TGraph *grLower, TGraph *grUpper)
{
  TGraphErrors *grFemtoModel = new TGraphErrors();
  grFemtoModel->SetName(grMedian1->GetName());
  double x, yM1, yLo, yUp;
  int count = 0;
  for(int i=0; i< grMedian1->GetN(); ++i) {
    grMedian1->GetPoint(i, x, yM1);
    grLower->GetPoint(i, x, yLo);
    grUpper->GetPoint(i, x, yUp);
    std::vector<float> yAll;
    yAll.push_back(yM1);
    yAll.push_back(yLo);
    yAll.push_back(yUp);
    std::sort(yAll.begin(), yAll.end());
    grFemtoModel->SetPoint(count, x/1000.f, (yAll[2]+yAll[0])/2.f);
    grFemtoModel->SetPointError(count++, 0, (yAll[2]+yAll[0])/2.f - yAll[0]);
  }
  return grFemtoModel;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraph *convertInGev(TGraph *gr) {
  TGraph *grOut = new TGraph();
  grOut->SetName(gr->GetName());
  double x,y;
  for(int i=0; i< gr->GetN(); ++i) {
    gr->GetPoint(i, x, y);
    grOut->SetPoint(i, x/1000.f, y);
  }
  return grOut;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *convertHistoInGev(TH1F *gr) {
  TGraphErrors *grOut = new TGraphErrors();
  grOut->SetName(gr->GetName());
  for(int i=0; i< gr->GetNbinsX(); ++i) {
    grOut->SetPoint(i, gr->GetBinCenter(i)/1000.f, gr->GetBinContent(i));
    grOut->SetPointError(i, 0, gr->GetBinError(i));
  }
  return grOut;
}

void dEtadPhiDist(const char *filename,const char *prefix) {
  SetStyle(false,true);
  const float normleft = 0.2;
  const float normright = 0.4;
  const float lowerBoundPhi=-TMath::Pi()/3;
  const float upperBoundPhi=5*TMath::Pi()/3;
  TFile *file=TFile::Open(filename);

  TString RelKNames[6]={"Proton","AntiProton","Lambda","AntiLambda","Xi","AntiXi"};
  TString ResultsName="";
  ResultsName+=prefix;
  ResultsName+="Results";

  TH2F *SEDistTmp[6][6];
  TH2F *SEDist[6][6];
  TH2F *SESumDist[4];
  TH2F *MEDistTmp[6][6];
  TH2F *MEDist[6][6];
  TH2F *MESumDist[4];
  TH2F *CFDist[6][6];
  TH1F *CFSumProjDist[4];
  TH2F *CFSumDist[4];
  int rebin=3;
  const char *draw="surf2";
  TDirectoryFile *dirResults=(TDirectoryFile*)(file->FindObjectAny(ResultsName.Data()));
  if (dirResults) {
    TList *tmp=dirResults->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *Results;
    dirResults->GetObject(name,Results);
    if (!Results) {
      std::cout << "No ResultsSample \n";
    } else {
      for (int iPart1=0;iPart1<6;++iPart1) {
        for (int iPart2=iPart1;iPart2<6;++iPart2) {
          TString folderName=Form("Particle%i_Particle%i",iPart1,iPart2);
          TList* tmpFolder=(TList*)Results->FindObject(folderName.Data());
          if (!tmpFolder) {
            std::cout << folderName.Data() <<" not Found\n";
          } else {
            TString SEName=Form("SEdPhidEtaDist_Particle%i_Particle%i",iPart1,iPart2);
            SEDist[iPart1][iPart2]=(TH2F*)tmpFolder->FindObject(SEName.Data());
            TString MEName=Form("MEdPhidEtaDist_Particle%i_Particle%i",iPart1,iPart2);
            MEDist[iPart1][iPart2]=(TH2F*)tmpFolder->FindObject(MEName.Data());
            float scaleSE=1./SEDist[iPart1][iPart2]->GetEntries();
            float scaleME=1./MEDist[iPart1][iPart2]->GetEntries();
            for (int iXBin=1;iXBin<SEDist[iPart1][iPart2]->GetNbinsX();++iXBin) {
              float dEta=SEDist[iPart1][iPart2]->GetXaxis()->GetBinCenter(iXBin);
              for (int iYBin=1;iYBin<SEDist[iPart1][iPart2]->GetNbinsY();++iYBin) {
                float dPhi=SEDist[iPart1][iPart2]->GetYaxis()->GetBinCenter(iYBin);
                if (dPhi<lowerBoundPhi) {
                  int SECounts=SEDist[iPart1][iPart2]->GetBinContent(iXBin,iYBin);
                  SEDist[iPart1][iPart2]->SetBinContent(iXBin,iYBin,0);
                  SEDist[iPart1][iPart2]->AddBinContent(SEDist[iPart1][iPart2]->FindBin(dEta,dPhi+2*TMath::Pi()),SECounts);

                  int MECounts=MEDist[iPart1][iPart2]->GetBinContent(iXBin,iYBin);
                  MEDist[iPart1][iPart2]->SetBinContent(iXBin,iYBin,0);
                  MEDist[iPart1][iPart2]->AddBinContent(MEDist[iPart1][iPart2]->FindBin(dEta,dPhi+2*TMath::Pi()),MECounts);
                } else if (dPhi > upperBoundPhi) {
                  int SECounts=SEDist[iPart1][iPart2]->GetBinContent(iXBin,iYBin);
                  SEDist[iPart1][iPart2]->SetBinContent(iXBin,iYBin,0);
                  SEDist[iPart1][iPart2]->AddBinContent(SEDist[iPart1][iPart2]->FindBin(dEta,dPhi-2*TMath::Pi()),SECounts);

                  int MECounts=MEDist[iPart1][iPart2]->GetBinContent(iXBin,iYBin);
                  MEDist[iPart1][iPart2]->SetBinContent(iXBin,iYBin,0);
                  MEDist[iPart1][iPart2]->AddBinContent(MEDist[iPart1][iPart2]->FindBin(dEta,dPhi-2*TMath::Pi()),MECounts);
                } else {
                  continue;
                }
              }
            }

            SEDist[iPart1][iPart2]->SetDirectory(0);
            SEDist[iPart1][iPart2]->Rebin2D(rebin,rebin);
            SEDist[iPart1][iPart2]->Scale(scaleSE);
            if (!SEDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }

            MEDist[iPart1][iPart2]->Rebin2D(rebin,rebin);
            MEDist[iPart1][iPart2]->SetDirectory(0);
            MEDist[iPart1][iPart2]->Scale(scaleME);
            if (!MEDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
          }
        }
      }
    }
  }
  auto *c1 = new TCanvas("dEtadPhi","dEtadPhi",2000,2000);

  c1->Divide(2,2);
  SESumDist[0]=(TH2F*)SEDist[0][0]->Clone("ProtonProtonSE");
  SESumDist[0]->Add(SEDist[1][1]);
  MESumDist[0]=(TH2F*)MEDist[0][0]->Clone("ProtonProtonME");
  MESumDist[0]->Add(MEDist[1][1]);
  CFSumDist[0]=(TH2F*)SESumDist[0]->Clone("Proton Proton");
  CFSumDist[0]->Divide(MESumDist[0]);
  CFSumDist[0]->SetTitle("p-p #oplus #bar{p}-#bar{p}");
  CFSumDist[0]->GetXaxis()->SetRangeUser(-1.4,1.4);
  CFSumDist[0]->GetYaxis()->SetRangeUser(.95*lowerBoundPhi,upperBoundPhi);
  c1->cd(1);
  CFSumDist[0]->DrawCopy(draw);

  SESumDist[1]=(TH2F*)SEDist[0][2]->Clone("ProtonLambdaSE");
  SESumDist[1]->Add(SEDist[1][3]);
  MESumDist[1]=(TH2F*)MEDist[0][2]->Clone("ProtonLambdaME");
  MESumDist[1]->Add(MEDist[1][3]);
  CFSumDist[1]=(TH2F*)SESumDist[1]->Clone("Proton Lambda");
  CFSumDist[1]->Divide(MESumDist[1]);
  CFSumDist[1]->SetTitle("p-#Lambda #oplus #bar{p}-#bar{#Lambda}");
  CFSumDist[1]->GetXaxis()->SetRangeUser(-1.4,1.4);
  CFSumDist[1]->GetYaxis()->SetRangeUser(.95*lowerBoundPhi,upperBoundPhi);
  c1->cd(2);
  CFSumDist[1]->DrawCopy(draw);


  SESumDist[2]=(TH2F*)SEDist[2][2]->Clone("LambdaLambdaSE");
  SESumDist[2]->Add(SEDist[3][3]);
  MESumDist[2]=(TH2F*)MEDist[2][2]->Clone("LambdaLambdaME");
  MESumDist[2]->Add(MEDist[3][3]);
  CFSumDist[2]=(TH2F*)SESumDist[2]->Clone("Lambda Lambda");
  CFSumDist[2]->Divide(MESumDist[2]);
  CFSumDist[2]->SetTitle("#Lambda-#Lambda #oplus #bar{#Lambda}-#bar{#Lambda}");
  CFSumDist[2]->GetXaxis()->SetRangeUser(-1.4,1.4);
  CFSumDist[2]->GetYaxis()->SetRangeUser(.95*lowerBoundPhi,upperBoundPhi);
  c1->cd(3);
  CFSumDist[2]->DrawCopy(draw);


  SESumDist[3]=(TH2F*)SEDist[0][4]->Clone("ProtonXiSE");
  SESumDist[3]->Add(SEDist[1][5]);
  MESumDist[3]=(TH2F*)MEDist[0][4]->Clone("ProtonXiME");
  MESumDist[3]->Add(MEDist[1][5]);
  CFSumDist[3]=(TH2F*)SESumDist[3]->Clone("Proton Xi");
  CFSumDist[3]->Divide(MESumDist[3]);
  CFSumDist[3]->SetTitle("p-#Xi #oplus #bar{p}-#bar{#Xi}");
  CFSumDist[3]->GetXaxis()->SetRangeUser(-1.4,1.4);
  CFSumDist[3]->GetYaxis()->SetRangeUser(.95*lowerBoundPhi,upperBoundPhi);
  c1->cd(4);
  CFSumDist[3]->DrawCopy(draw);

  auto *c2 = new TCanvas("dPhi","dPhi",1000,1000);
  c2->Divide(2,2);
  CFSumProjDist[0]=(TH1F*)(SESumDist[0]->ProjectionY("phiProtonProton",SESumDist[0]->GetXaxis()->FindBin(-1.2),SESumDist[0]->GetXaxis()->FindBin(1.2)));
  CFSumProjDist[0]->Divide((TH1F*)(MESumDist[0]->ProjectionY("phiMEProtonProton",MESumDist[0]->GetXaxis()->FindBin(-1.2),MESumDist[0]->GetXaxis()->FindBin(1.2))));
  CFSumProjDist[0]->SetTitle("p-p #oplus #bar{p}-#bar{p}");
  CFSumProjDist[0]->GetYaxis()->SetTitle("#it{C}(#Delta#phi)");
  CFSumProjDist[0]->GetXaxis()->SetRangeUser(.95*lowerBoundPhi,upperBoundPhi);
  SetStyleHisto(CFSumProjDist[0],0,2);
  c2->cd(1);
  CFSumProjDist[0]->DrawCopy();

  CFSumProjDist[1]=(TH1F*)(SESumDist[1]->ProjectionY("phiProtonLambda",SESumDist[1]->GetXaxis()->FindBin(-1.2),SESumDist[1]->GetXaxis()->FindBin(1.2)));
  CFSumProjDist[1]->Divide((TH1F*)(MESumDist[1]->ProjectionY("phiMEProtonLambda",MESumDist[1]->GetXaxis()->FindBin(-1.2),MESumDist[1]->GetXaxis()->FindBin(1.2))));
  CFSumProjDist[1]->SetTitle("p-#Lambda #oplus #bar{p}-#bar{#Lambda}");
  CFSumProjDist[1]->GetYaxis()->SetTitle("#it{C}(#Delta#phi)");
  CFSumProjDist[1]->GetXaxis()->SetRangeUser(.95*lowerBoundPhi,upperBoundPhi);
  SetStyleHisto(CFSumProjDist[1],0,2);
  c2->cd(2);
  CFSumProjDist[1]->DrawCopy();

  CFSumProjDist[2]=(TH1F*)(SESumDist[2]->ProjectionY("phiLambdaLambda",SESumDist[2]->GetXaxis()->FindBin(-1.2),SESumDist[2]->GetXaxis()->FindBin(1.2)));
  CFSumProjDist[2]->Divide((TH1F*)(MESumDist[2]->ProjectionY("phiMELambdaLambda",MESumDist[2]->GetXaxis()->FindBin(-1.2),MESumDist[2]->GetXaxis()->FindBin(1.2))));
  CFSumProjDist[2]->SetTitle("#Lambda-#Lambda #oplus #bar{#Lambda}-#bar{#Lambda}");
  CFSumProjDist[2]->GetYaxis()->SetTitle("#it{C}(#Delta#phi)");
  CFSumProjDist[2]->GetXaxis()->SetRangeUser(.95*lowerBoundPhi,upperBoundPhi);
  SetStyleHisto(CFSumProjDist[2],0,2);
  c2->cd(3);
  CFSumProjDist[2]->DrawCopy();

  CFSumProjDist[3]=(TH1F*)(SESumDist[3]->ProjectionY("phiProtonXi",SESumDist[3]->GetXaxis()->FindBin(-1.2),SESumDist[3]->GetXaxis()->FindBin(1.2)));
  CFSumProjDist[3]->Divide((TH1F*)(MESumDist[3]->ProjectionY("phiMEProtonXi",MESumDist[3]->GetXaxis()->FindBin(-1.2),MESumDist[3]->GetXaxis()->FindBin(1.2))));
  CFSumProjDist[3]->SetTitle("p-#Xi #oplus #bar{p}-#bar{#Xi}");
  CFSumProjDist[3]->GetYaxis()->SetTitle("#it{C}(#Delta#phi)");
  CFSumProjDist[3]->GetXaxis()->SetRangeUser(.95*lowerBoundPhi,upperBoundPhi);
  SetStyleHisto(CFSumProjDist[3],0,2);
  c2->cd(4);
  CFSumProjDist[3]->DrawCopy();

}
