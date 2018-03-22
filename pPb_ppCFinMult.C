
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
void pPb_ppCFinMult(const char* input) {
  TList *TPdir=new TList();
  TPdir->SetName("TPdir_0");
  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;
  gStyle->SetCanvasPreferGL(1);
  TFile *file=TFile::Open(input);
  TH2F* SEMultDist[2];
  TH2F* MEMultDist[2];
  TH1F* SEDist[2][3];
  TH1F* MEDist[2][3];
  TDirectoryFile *dirResults=(TDirectoryFile*)(file->FindObjectAny("Results"));
  TString PartName[2]={"Proton","AntiProton"};
  int MultBins[4]={1,8,12,21};
  if (dirResults) {
    TList *tmp=dirResults->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *Results;
    dirResults->GetObject(name,Results);
    if (!Results) {
      std::cout << "No Results \n";
    } else {
      for (int iPart=0;iPart<2;++iPart){
        TString folderName=Form("Particle%i_Particle%i",iPart,iPart);
        TList* tmpFolder=(TList*)Results->FindObject(folderName.Data());
        if (!tmpFolder) {
          std::cout << folderName.Data() <<" not Found\n";
        } else {
          TString SEName=Form("SEMultDist_Particle%i_Particle%i",iPart,iPart);
          SEMultDist[iPart]=(TH2F*)tmpFolder->FindObject(SEName.Data());
          SEMultDist[iPart]->SetDirectory(0);
          if (!SEMultDist[iPart]) {
            std::cout << SEName.Data() << " not Found\n";
            std::cout << SEName.Data() << " not Found\n";
          }
          TString MEName=Form("MEMultDist_Particle%i_Particle%i",iPart,iPart);
          MEMultDist[iPart]=(TH2F*)tmpFolder->FindObject(MEName.Data());
          MEMultDist[iPart]->SetDirectory(0);
          if (!MEMultDist[iPart]) {
            std::cout << SEName.Data() << " not Found\n";
            std::cout << SEName.Data() << " not Found\n";
          }
          for (int iMult=0;iMult<3;++iMult) {
            TString projSEName=Form("f%s%sRelK_%i",PartName[iPart].Data(),PartName[iPart].Data(),iMult);
            SEDist[iPart][iMult]=(TH1F*)SEMultDist[iPart]->ProjectionX(projSEName.Data(),MultBins[iMult],MultBins[iMult+1]-1);
            TPdir->Add(SEDist[iPart][iMult]);
            TString projMEName=Form("f%s%sRelKME_%i",PartName[iPart].Data(),PartName[iPart].Data(),iMult);
            MEDist[iPart][iMult]=(TH1F*)MEMultDist[iPart]->ProjectionX(projMEName.Data(),MultBins[iMult],MultBins[iMult+1]-1);
            TPdir->Add(MEDist[iPart][iMult]);
          }
        }
      }
    }
  }
  TH1F *hist_CF_pp_ApAp_exp[3][3];
  TCanvas *cApAp_pp=new TCanvas("cPAP","cPAP",0,0,1500,1100);
  cApAp_pp->Divide(2,1);
  TCanvas *cSum=new TCanvas("CanSum","CanSum",0,0,1500,1100);
  auto* leg= new TLegend(0.45, 0.7, 0.95, 0.95);
  for (int iMult=0;iMult<3;++iMult) {
    TString ppMultName = Form("hist_CF_pp_%i",iMult);
    hist_CF_pp_ApAp_exp[0][iMult] = Calculate_CF(SEDist[0][iMult],MEDist[0][iMult],ppMultName.Data(),normleft,normright);
    SetStyleHisto(hist_CF_pp_ApAp_exp[0][iMult],21,iMult+2);
    cApAp_pp->cd(1);
    if(iMult==0)hist_CF_pp_ApAp_exp[0][iMult]->DrawCopy();
    else hist_CF_pp_ApAp_exp[0][iMult]->DrawCopy("SAME");
    TString ApApMultName = Form("hist_CF_ApAp_%i",iMult);
    hist_CF_pp_ApAp_exp[1][iMult] = Calculate_CF(SEDist[0][iMult],MEDist[0][iMult],ApApMultName.Data(),normleft,normright);
    SetStyleHisto(hist_CF_pp_ApAp_exp[1][iMult],21,iMult+2);
    cApAp_pp->cd(2);
    if(iMult==0)hist_CF_pp_ApAp_exp[1][iMult]->DrawCopy();
    else hist_CF_pp_ApAp_exp[1][iMult]->DrawCopy("SAME");
    TString SumppApApMultName = Form("hist_CF_pp_ApAp_exp_sum_%i",iMult);
    hist_CF_pp_ApAp_exp[2][iMult] = add_CF(hist_CF_pp_ApAp_exp[0][iMult],hist_CF_pp_ApAp_exp[1][iMult],SumppApApMultName.Data());
    SetStyleHisto(hist_CF_pp_ApAp_exp[2][iMult],20,iMult+2);
    hist_CF_pp_ApAp_exp[2][iMult]->SetStats(0);
    hist_CF_pp_ApAp_exp[2][iMult]->GetXaxis()->SetTitleOffset(.9);
    hist_CF_pp_ApAp_exp[2][iMult]->GetXaxis()->SetRangeUser(0,0.125);
    hist_CF_pp_ApAp_exp[2][iMult]->GetYaxis()->SetTitleOffset(0.9);
    hist_CF_pp_ApAp_exp[2][iMult]->SetTitle(";k* (GeV/#it{c}); #it{C}(k*)");
    cSum->cd();
    leg->AddEntry(hist_CF_pp_ApAp_exp[2][iMult],Form("Mult. Bin %i",iMult),"lp");
    if(iMult==0)hist_CF_pp_ApAp_exp[2][iMult]->DrawCopy();
    else hist_CF_pp_ApAp_exp[2][iMult]->DrawCopy("SAME");
  }
  cSum->cd();
  leg->Draw("SAME");
  TFile *output=new TFile("ppMuti.root","RECREATE");
  TDirectoryFile* dirOutput=new TDirectoryFile("PWGCF_PLFemto_0","PWGCF_PLFemto_0");
  dirOutput->Add(TPdir);
  dirOutput->Write("PWGCF_PLFemto_0",TObject::kSingleKey);
}
