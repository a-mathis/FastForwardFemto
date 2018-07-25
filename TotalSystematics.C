# include "TCanvas.h"
std::vector<int> fillColors = {kGray + 1,  kRed - 10,    kBlue - 9,
    kGreen - 8, kMagenta - 9, kOrange - 9,
    kCyan - 8,  kYellow - 7};
std::vector<int> colors = {kBlack,       kRed + 1,    kBlue + 2, kGreen + 3,
    kMagenta + 1, kOrange - 1, kCyan + 2, kYellow + 2};
std::vector<int> markers = {
    kFullCircle, kFullSquare, kOpenCircle,  kOpenSquare, kOpenDiamond,
    kOpenCross,  kFullCross,  kFullDiamond, kFullStar,   kOpenStar};

static double normleft = 0.2;
static double normright = 0.4;
static const char *prefix = "MB";
static int nSysFiles=2;
static bool useFit = false;
static const char* fit = "pol2";
static double fitRange = 1.;
//Set here the names of the variations in the order they are added to the
//output file, same for the ranges. Later this is used to loop over them, so do
//this properly! This was created in conjunction with the AddTaskFemtoDreamSysVar,
//so look there for inspiration. up generall refers to 'stricter', down to 'looser'
static int BinningPPDefault = 1; //Default Binnign for the CF
static int BinningPPSyst = 10; //Binning for estimation of Uncertainties
static int BinningPLDefault = 5; //Default Binnign for the CF
static int BinningPLSyst = 50; //Binning for estimation of Uncertainties
static int BinningLLDefault = 5; //Default Binnign for the CF
static int BinningLLSyst = 50; //Binning for estimation of Uncertainties
static int BinningPXiDefault = 5; //Default Binnign for the CF
static int BinningPXiSyst = 50; //Binning for estimation of Uncertainties

static int PartNumberProton = 0;
static int PartNumberAntiProton = 1;
static int ProtonVarTaskMin = 1;
static int ProtonVarTaskMax = 9;//-> delta 8
static std::vector<const char*> ProtonVarName = {
    "Trk #it{p}_{T} down","Trk #it{p}_{T} up",
    "Trk #eta up","Trk #eta down",
    "Trk n#sigma up","Trk n#sigma down",
    "Trk FilterBit",
    "Trk TPC cluster down","Trk TPC cluster up"
};

static int PartNumberLambda= 2;
static int PartNumberAntiLambda= 3;
static int LambdaVarTaskMin = 10;
static int LambdaVarTaskMax = 18; //delta -> 9
static std::vector<const char*> LambdaVarName = {
    "V0 #it{p}_{T} down",
    "V0 #it{p}_{T} up",
    "V0 CPA up",
    "V0 Daug n#sigma down",
    "V0 Daug TPC ncls 80",
    "V0 Daug #eta up",
    "V0 Daug #eta down",
    "V0 Daug d_{track} down",
    "V0 Daug d_{Track, PV} up"
};

static int PartNumberXi = 4;
static int PartNumberAntiXi= 5;
static int XiVarTaskMin = 19;
static int XiVarTaskMax = 41;//delta -> 23
static std::vector<const char*> XiVarName = {
    "Casc Daug d_{track} up","Casc Daug d_{track} down",
    "Casc Bach d_{Track, PV} down","Casc Bach d_{Track, PV} up",
    "Casc CPA up 1", "Casc CPA up 2",
    "Casc Transverse Radius down", "Casc Transverse Radius up",
    "Casc V0 Daug d_{track} up 1", "Casc V0 Daug d_{track} up 2",
    "Casc V0 CPA down", "Casc V0 CPA up",
    "Casc V0 Transverse Radius down", "Casc V0 Transverse Radius up",
    "Casc V0 d_{V0, PV} up", "Casc V0 d_{V0, PV} down",
    "Casc V0 Daug d_{Track, PV} down", "Casc V0 Daug d_{Track, PV} up",
    "Casc Track #eta up", "Casc Track #eta down",
    "Casc Track n#sigma up","Casc Track n#sigma down",
    "Casc #it{p}_{T} down"
};
static int Rebin = 1;//Dont touch this one, will be set around in the functions, change binning above!
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
void StyleGraph(TGraph *gr) {
  gr->GetYaxis()->SetTitleSize(0.07);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetTitleSize(0.07);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->SetMarkerStyle(20);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

      Float_t CF_sum_average =
          (CF1_err_weight*CF1_val + CF2_err_weight*CF2_val)/
          (CF1_err_weight+CF2_err_weight);
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


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Calculate_CFnormalization(TH1F* histRE_relK, TH1F* histME_relK)
{
  Double_t norm_relK =
      histRE_relK->Integral(
          histRE_relK->FindBin(normleft),
          histRE_relK->FindBin(normright)) /
          histME_relK->Integral(
              histME_relK->FindBin(normleft),
              histME_relK->FindBin(normright));
  return norm_relK;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* Calculate_CF(TH1F* histRE_relK,TH1F* histME_relK, TString CFname)
{
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  Double_t norm_relK =
      Calculate_CFnormalization(histRE_relK,histME_relK);

  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  Hist_CF->Divide(histRE_relK,histME_relK,1,norm_relK);

  return Hist_CF;
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TH1F* ReadHistograms (TString *FileName,int nFiles,TString Variation,
                      int Part1, int Part2)
{
  TFile* _file0=nullptr;
  TDirectoryFile *dirResults= nullptr;
  for (int iFiles=0;iFiles<nFiles;++iFiles) {
    _file0=TFile::Open(FileName[iFiles]);
    if (!_file0) {
      std::cout << FileName[iFiles] <<
          " does not exist, so there is no file, have a nice day \n";
      return nullptr;
    }
    if (_file0->FindObjectAny(Form("%sResults%s", prefix,Variation.Data()))) {
      dirResults=
          (TDirectoryFile*)(_file0->FindObjectAny(Form(
              "%sResults%s", prefix,Variation.Data())));
      break;
    }
  }
  if (!dirResults) {
    std::cout << Form("%sResults%s", prefix,Variation.Data()) <<
        " does not exist in any file, so there is no directory, "
        "have a nice day \n";
    return nullptr;
  }
  TList *Results=nullptr;
  dirResults->GetObject(
      Form("%sResults%s", prefix, Variation.Data()),Results);
  if (!Results) {
    std::cout << Form("%sResults%s", prefix,Variation.Data()) <<
        " does not exist, so there is no List, have a nice day \n";
    return nullptr;
  }
  TList* ParticlePair=nullptr;
  ParticlePair=
      (TList*)Results->FindObject(Form("Particle%d_Particle%d",Part1,Part2));
  if (!ParticlePair) {
    std::cout << Form("Particle%d_Particle%d",Part1,Part2) <<
        " does not exist, so there is no Hist List, have a nice day \n";
    return nullptr;
  }
  TH1F* histSE=nullptr;
  histSE=(TH1F*)ParticlePair->FindObject(
      Form("SEDist_Particle%d_Particle%d",Part1,Part2));
  if (!histSE) {
    std::cout << Form("SEDist_Particle%d_Particle%d",Part1,Part2) <<
        " does not exist, so there is no Histogram, have a nice day \n";
    return nullptr;
  }
  histSE->Rebin(Rebin);

  TH1F* histME=nullptr;
  histME=(TH1F*)ParticlePair->FindObject(
      Form("MEDist_Particle%d_Particle%d",Part1,Part2));
  if (!histME) {
    std::cout << Form("MEDist_Particle%d_Particle%d",Part1,Part2) <<
        " does not exist, so there is no Histogram, have a nice day \n";
    return nullptr;
  }
  histME->Rebin(Rebin);
  return Calculate_CF(
      histSE,histME,Form("CF_Part%d_Part%d",Part1,Part2));

}

TH1F *GetSummedCF(TString *FileName,int nFiles,TString Variation,
                  int Part1, int Part2,int APart1, int APart2)
{
  TH1F *CFParticlePair=
      ReadHistograms(FileName,nFiles,Variation,Part1,Part2);
  if (!CFParticlePair) {
    std::cout << Form("CFPair_Part%d_Part%d",Part1,Part2) <<
        " does not exist, so there is no Histogram, have a nice day \n";
    return nullptr;
  }
  TH1F *CFAntiParticlePair=
      ReadHistograms(FileName,nFiles,Variation,APart1,APart2);
  if (!CFAntiParticlePair) {
    std::cout << Form("CFPair_Part%d_Part%d",APart1,APart2) <<
        " does not exist, so there is no Histogram, have a nice day \n";
    return nullptr;
  }
  return add_CF(
      CFParticlePair,CFAntiParticlePair,
      Form("CF_Pair%d%d_APair%d%d",Part1,Part2,APart1,APart2));
}


void Barlow(TH1F* DefaultCF, TH1F* Variations,TH1F *Barlows,const char* varNames) {
  int nBins=DefaultCF->GetXaxis()->FindBin(fitRange);

  Barlows->SetTitle(varNames);
  Barlows->GetXaxis()->SetRangeUser(0.,0.2);
  for (int iKstar=1;iKstar<nBins;++iKstar) {
    double nominalVal = DefaultCF->GetBinContent(iKstar);
    double nominalErr =  DefaultCF->GetBinError(iKstar);
    double variationVal= Variations->GetBinContent(iKstar);
    double variationErr = Variations->GetBinError(iKstar);
    double statErrVariation = std::sqrt(
        variationErr*variationErr + nominalErr*nominalErr);
    Barlows->SetBinContent(iKstar, std::abs(variationVal - nominalVal)/statErrVariation);
    Barlows->SetBinError(iKstar, 0);
  }
  DrawHist(Barlows,"k* (GeV/#it{c})","n#sigma",0,1,1.,true);
  TLine Line = TLine();
  Line.SetLineColor(colors[4]);
  Line.SetLineWidth(2);
  Line.DrawLine(0,2.0,0.2,2.0);
  return;
}

void RatioWithDefaultAndBarlow(TH1F* DefaultCF, TH1F** Variations,TH1F** Ratios,
                               std::vector<const char*> varNames,
                               int iSplitX,int iSplitY, TString Particle)
{
  const int nVariations=varNames.size();
  int nBins=DefaultCF->GetXaxis()->FindBin(fitRange);
  double kMax = DefaultCF->GetXaxis()->GetBinUpEdge(DefaultCF->FindBin(fitRange));

  TString canNameRatio=DefaultCF->GetName();
  canNameRatio+="RatioCan";
  TCanvas *c1 = new TCanvas(canNameRatio.Data(),canNameRatio.Data(),2000,1000);
  c1->Divide(iSplitX,iSplitY);

  TH1F* Barlows[nVariations];
  TString canNameBarlow=DefaultCF->GetName();
  canNameBarlow+="BarlowCan";
  TCanvas *c2 = new TCanvas(canNameBarlow.Data(),canNameBarlow.Data(),2000,1000);
  c2->Divide(iSplitX,iSplitY);

  if (!DefaultCF) {
    std::cout << "Default CF missing! \n";
    return;
  }
  int iVar=0;
  for (auto &it:varNames) {
    c1->cd(iVar+1);
    TString HistNameRatio=Form("Ratio_%s",it);
    Ratios[iVar]=(TH1F*)DefaultCF->Clone(HistNameRatio.Data());
    if (!Variations[iVar]) {
      std::cout << "Variation " << varNames[iVar] << " missing!\n";
      return;
    }
    Ratios[iVar]->Divide(Variations[iVar]);
    Ratios[iVar]->SetTitle(varNames[iVar]);
    Ratios[iVar]->GetXaxis()->SetRangeUser(0.,0.2);
    DrawHist(Ratios[iVar],"k* (GeV/#it{c})","#it{C}(k*)_{default}/C(k*)_{variation}",0,1,1.,true);
    TF1 *fitRatio = new TF1(Form("fFitRatio_%s_%i",DefaultCF->GetName(),iVar) ,fit,0., 1.0);
    if(useFit) {
      Ratios[iVar]->Fit(fitRatio,"NRQ", "", 0, fitRange);
      Ratios[iVar]->GetListOfFunctions()->Add(fitRatio);
      fitRatio->Draw("SAME");
    }
    c2->cd(iVar+1);
    TString HistNameBarlow=Form("Barlow_%s",it);
    Barlows[iVar]=new TH1F(HistNameBarlow.Data(),HistNameBarlow.Data(),nBins,0,kMax);
    Barlow(DefaultCF,Variations[iVar],Barlows[iVar],it);

    iVar++;
  }
  c1->SaveAs(Form("%s_ratios.pdf",Particle.Data()));
  c2->SaveAs(Form("%s_barlow.pdf",Particle.Data()));
  return;
}

void SystematicAndErrBudget(TH1F* DefaultCF, TH1F* VariationRatio,
                            TH1F** ErrBudget, TH1F** AbsErr, const char* VarName,
                            int iVar)
{
  int nBins=DefaultCF->GetXaxis()->FindBin(fitRange);
  double kMax = DefaultCF->GetXaxis()->GetBinUpEdge(DefaultCF->FindBin(fitRange));

  TString ErrBudgetName=Form("ErrBudget%s_%i",DefaultCF->GetName(),iVar);
  ErrBudget[iVar]=new TH1F(ErrBudgetName.Data(),ErrBudgetName.Data(),nBins,0,kMax);

  TString AbsErrName=Form("AbsErr%s_%d",DefaultCF->GetName(),iVar);
  AbsErr[iVar]=new TH1F(AbsErrName.Data(),AbsErrName.Data(),nBins,0,kMax);

  for (int iKstar=1;iKstar<nBins;++iKstar) {
    double RatioVal=0;
    if (useFit) {
      TList *functList = VariationRatio->GetListOfFunctions();
      TF1 *ratioFit = nullptr;
      ratioFit= (TF1*)functList->FindObject(Form("fFitRatio_%s_%i",DefaultCF->GetName(),iVar));
      if (!ratioFit) {
        std::cout << "No Ratio Fit found for " <<
            Form("fFitRatio_%s_%i",DefaultCF->GetName(),iVar) << std::endl;
        return;
      }
      double kStarValue=DefaultCF->GetBinCenter(iKstar);
      RatioVal=ratioFit->Eval(kStarValue);
    } else {
      RatioVal = VariationRatio->GetBinContent(iKstar);
    }
    AbsErr[iVar]->SetBinContent(iKstar,fabs(1. - pow(RatioVal,-1)) * DefaultCF->GetBinContent(iKstar));
    ErrBudget[iVar]->SetBinContent(iKstar,fabs(1. - pow(RatioVal,-1)));
  }
}

void PlotSystematicsPerVariation(TH1F* DefaultCF, TH1F** VariationRatio,
                                 TH1F** ErrBudget, TH1F** AbsErr,
                                 std::vector<const char*> VarName,
                                 int iSplitX,int iSplitY,
                                 TString Particle)
{
  TString CanBudgetName=Form("cErrBudget%s",DefaultCF->GetName());
  TCanvas* cErrBudget = new TCanvas(CanBudgetName.Data(),CanBudgetName.Data(),2000,1000);
  cErrBudget->Divide(iSplitX,iSplitY);
  TString CanAbsName=Form("cAbsErr%s",DefaultCF->GetName());
  TCanvas* cAbsErr= new TCanvas(CanAbsName.Data(),CanAbsName.Data(),2000,1000);
  cAbsErr->Divide(iSplitX,iSplitY);
  const int Variations=VarName.size();
  for (int iVar=0;iVar<Variations;++iVar) {
    AbsErr[iVar]= nullptr;
    ErrBudget[iVar]=nullptr;
    SystematicAndErrBudget(DefaultCF,VariationRatio[iVar],ErrBudget,AbsErr,ProtonVarName[iVar],iVar);
    if (!ErrBudget[iVar]) {
      std::cout << ProtonVarName[iVar] << "is Missing its Err Budget in the Proton Proton Systematic \n";
      return;
    }
    if (!AbsErr[iVar]) {
      std::cout << ProtonVarName[iVar] << "is Missing its Abs Err in the Proton Proton Systematic \n";
      return;
    }
    cErrBudget->cd(iVar+1);
    ErrBudget[iVar]->SetTitle(VarName[iVar]);
    if(Particle!="PP") {
      ErrBudget[iVar]->GetXaxis()->SetRangeUser(0,0.6);
    } else {
      ErrBudget[iVar]->GetXaxis()->SetRangeUser(0,0.2);
    }
    DrawHist(ErrBudget[iVar],"k* (GeV/#it{c})","|1 - #it{C}(k*)_{variation}/#it{C}(k*)_{default}|",1,2,1.,true);
    TF1 *f1 = new TF1("f1","pol0",0,1);
    ErrBudget[iVar]->Fit("f1","WR","",0,0.2); //does not take the errors of the ratio into account
    f1->Draw("SAME");
    TLatex LambdaLabel;
    LambdaLabel.SetNDC(kTRUE);
    LambdaLabel.SetTextSize(gStyle->GetTextSize()*1.2);
    LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,Form("Par Fit: %.2f",f1->GetParameter(0)*100));

    cAbsErr->cd(iVar+1);
    AbsErr[iVar]->SetTitle(VarName[iVar]);
    if(Particle!="PP") {
      AbsErr[iVar]->GetXaxis()->SetRangeUser(0,0.6);
    } else {
      AbsErr[iVar]->GetXaxis()->SetRangeUser(0,0.2);
    }
    TF1 *f2 = new TF1("f2","pol0",0,1);
    DrawHist(AbsErr[iVar],"k* (GeV/#it{c})","|#it{C}(k*)_{default} - #it{C}(k*)_{variation}|",1,2,1.,true);
    AbsErr[iVar]->Fit("f2","WR","",0,0.04); //does not take the errors of the ratio into account
    f2->Draw("SAME");
    LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,Form("Par Fit: %.2f",f2->GetParameter(0)));
  }
  cAbsErr->SaveAs(Form("AbsoluteErr%s.pdf",Particle.Data()));
  cErrBudget->SaveAs(Form("Budget%s.pdf",Particle.Data()));
}

void PlottingSysErrors(TH1F* DefaultCF, TH1F* OrgCF, TGraphErrors* SysErr,
                       TGraph* RelErr,TString Particle)
{
  TString CanName = Form("Can%s_PlottingSysErr",DefaultCF->GetName());
  TCanvas* c1 = new TCanvas(CanName.Data(),CanName.Data(),2000,1000);
  c1->Divide(2,1);
  c1->cd(1);
  c1->cd(1)->SetLeftMargin(0.1);
  c1->cd(1)->SetRightMargin(0.01);
  DrawHist(DefaultCF,"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.,true);
  SysErr->SetFillStyle(0);
  SysErr->Draw("2same");
  c1->cd(2);
  c1->cd(2)->SetLeftMargin(0.2);

  TF1 *fRatio = new TF1(Form("Ratio_%s",DefaultCF->GetName()), "pol2", 0, fitRange*1.5);
//  StyleGraph(RelErr);
  RelErr->SetMarkerStyle(20);
  RelErr->SetTitle("Relative Error");
  RelErr->GetXaxis()->SetTitle("k* (GeV/#it{c})");
  RelErr->GetYaxis()->SetTitle("#sigma_{Rel,Sum}");
  if (Particle!="PP") {
    fRatio->SetParameter(2,100);
    fRatio->SetParLimits(2,0,1e6);
    RelErr->Fit(Form("Ratio_%s",DefaultCF->GetName()), "R", "", 0.01, fitRange);
    RelErr->GetXaxis()->SetRangeUser(0,0.5);
  }  else {
    RelErr->Fit(Form("Ratio_%s",DefaultCF->GetName()), "R", "", 0.01, fitRange);
    RelErr->GetXaxis()->SetRangeUser(0,0.2);
  }
  RelErr->SetMarkerStyle(20);

  RelErr->Draw("ape");

  TH1F* outputParam =
      new TH1F(
          Form("Para%s",Particle.Data()),
          Form("Para%s",Particle.Data()),3,0,3);
  outputParam->SetBinContent(1,fRatio->GetParameter(0));
  outputParam->SetBinContent(2,fRatio->GetParameter(1));
  outputParam->SetBinContent(3,fRatio->GetParameter(2));

  TGraphErrors *grFinalError= new TGraphErrors();
  TString HistName=Form("hSystError%s",Particle.Data());
  double maxX = OrgCF->GetXaxis()->GetBinUpEdge(OrgCF->GetNbinsX());
  TH1F *hSystError = new TH1F(HistName.Data(),"",OrgCF->GetNbinsX(), 0,maxX);

  const float xerr = OrgCF->GetBinWidth(1)/2.;
  for(int kstarBin = 0; kstarBin < OrgCF->GetXaxis()->FindBin(fitRange); kstarBin++)
  {
    const float x = OrgCF->GetBinCenter(kstarBin+1);
    const float y = OrgCF->GetBinContent(kstarBin+1);
    grFinalError->SetPoint(kstarBin, x, y);
    grFinalError->SetPointError(kstarBin, xerr,  y * fRatio->Eval(x));
    hSystError->Fill(x, y * fRatio->Eval(x));
  }
  grFinalError->SetFillColor(kGray+1);
  grFinalError->SetMarkerStyle(20);
  grFinalError->SetLineWidth(3);

  TString OutName = Form("C2total%s.root",Particle.Data());
  TFile *OutFile = new TFile(OutName.Data(), "RECREATE");
  hSystError->Write(Form("C2totalsys%s",Particle.Data()));
  outputParam->Write(Form("SysParam%s",Particle.Data()));

  TString CanName2 = Form("CanPlotting%s",Particle.Data());
  TCanvas *CanPlotting = new TCanvas(CanName2.Data(),CanName2.Data(),0,0,500,500);
  DrawHist(OrgCF,"k* (GeV/#it{c})","#it{C}(k*)",0,1,1.);
  if (Particle == "PP") {
    OrgCF->SetTitle("p-p #oplus #bar{p}#bar{p}");
    OrgCF->GetXaxis()->SetRangeUser(0,0.12);
  } else if (Particle == "PL") {
    OrgCF->SetTitle("p-#Lambda #oplus #bar{p}#bar{#Lambda}");
    OrgCF->GetXaxis()->SetRangeUser(0,0.2);
  } else if (Particle == "LL") {
    OrgCF->SetTitle("#Lambda-#Lambda #oplus #bar{#Lambda}#bar{#Lambda}");
    OrgCF->GetXaxis()->SetRangeUser(0,0.2);
  } else if (Particle == "PXi") {
    OrgCF->SetTitle("p-#Xi #oplus #bar{p}#bar{#Xi}");
    OrgCF->GetXaxis()->SetRangeUser(0,0.2);
  } else {
    std::cout << "You are using strange Particles! \n";
  }
  //  OrgCF->GetYaxis()->SetRangeUser(0,4);
  grFinalError->Draw("2 same");
  OrgCF->Draw("same");
  c1->SaveAs(Form("%s_CFSysErrRebinned.pdf",Particle.Data()));
  CanPlotting->SaveAs(Form("%s_FinalCF.pdf",Particle.Data()));
  return ;
}


void ProtonProtonSystematic(TString FileDefault, TString *FileSystematics)
{
  //Get the Original Binned CF
  Rebin=BinningPPDefault;
  TH1F* ProtonOriginal =
      GetSummedCF(&FileDefault,1,"",
                  PartNumberProton,PartNumberProton,
                  PartNumberAntiProton,PartNumberAntiProton);
  //Get the Default CF
  Rebin=BinningPPSyst;
  TH1F* ProtonDefault =
      GetSummedCF(&FileDefault,1,"",
                  PartNumberProton,PartNumberProton,
                  PartNumberAntiProton,PartNumberAntiProton);
  const int Variations=ProtonVarTaskMax-ProtonVarTaskMin+1;
  if (!(Variations==ProtonVarName.size())) {
    std::cout<<
        "Number of Proton Variations are: "<<Variations<<
        " Number of Names given are: "<<ProtonVarName.size()<<
        " Please fix! \n";
    return;
  }
  TH1F* ProtonVariations[Variations];
  for (int iVar=0;iVar<Variations;++iVar) {
    //Get the summed CF for the variations
    ProtonVariations[iVar]=
        GetSummedCF(FileSystematics,nSysFiles,Form("%d",iVar+1),
                    PartNumberProton,PartNumberProton,
                    PartNumberAntiProton,PartNumberAntiProton);
  }
  TH1F* Ratios[Variations];
  RatioWithDefaultAndBarlow(ProtonDefault,ProtonVariations,Ratios,ProtonVarName,3,3,"PP");

  TH1F* ErrBudget[Variations];
  TH1F* AbsErr[Variations];
  PlotSystematicsPerVariation(ProtonDefault,Ratios,ErrBudget,AbsErr,ProtonVarName,3,3,"PP");
  //now that we have relative errors & aboslute values, we have to think about which
  //ones to combine since they are up/down variations:
  //Track Cut -- Variations (as stored in the arrays)
  //    pT          0,1
  //    eta         2,3
  //    nSigma      4,5
  //    TPC Cls     7,8
  //After that loop over all k* bins of interest to combine them
  int nBins=ProtonDefault->GetXaxis()->FindBin(fitRange);
  const int EffVariations=5;
  double addMeUp[EffVariations];
  TGraphErrors *graphSysErr = new TGraphErrors();
  TGraph *graphRelErr = new TGraph(nBins);
  for (int iKstar=1;iKstar<nBins;++iKstar) {
    double SysErrTotal = 0;
    addMeUp[0] = (AbsErr[0]->GetBinContent(iKstar)+AbsErr[1]->GetBinContent(iKstar))/2.;//pT
    addMeUp[1] = (AbsErr[2]->GetBinContent(iKstar)+AbsErr[3]->GetBinContent(iKstar))/2.;//Eta
    addMeUp[2] = (AbsErr[4]->GetBinContent(iKstar)+AbsErr[5]->GetBinContent(iKstar))/2.;//nSigma
    addMeUp[3] = AbsErr[6]->GetBinContent(iKstar); // Filter Bit 96
    addMeUp[4] = (AbsErr[7]->GetBinContent(iKstar)+AbsErr[8]->GetBinContent(iKstar))/2.;//TPC Cls
    for (int iAdd=0;iAdd<EffVariations;++iAdd) {
      if (iKstar == 1) {
        std::cout << iAdd << '\t'<< addMeUp[iAdd] << '\n';
      }
      SysErrTotal+=(addMeUp[iAdd]*addMeUp[iAdd]);
    }
    SysErrTotal=TMath::Sqrt(SysErrTotal);
    double xValue=ProtonDefault->GetBinCenter(iKstar);
    graphSysErr->SetPoint(iKstar,xValue,ProtonDefault->GetBinContent(iKstar));
    graphSysErr->SetPointError(iKstar, ProtonDefault->GetBinWidth(1)/2.,SysErrTotal);
    double relErr = SysErrTotal/ProtonDefault->GetBinContent(iKstar);
    graphRelErr->SetPoint(iKstar,xValue,relErr);
  }
  ProtonDefault->SetTitle("p-p #oplus #bar{p}#bar{p}");
  ProtonDefault->GetXaxis()->SetRangeUser(0, 0.13);
  ProtonDefault->GetYaxis()->SetRangeUser(0.85,2.2);
  ProtonDefault->GetYaxis()->SetTitleOffset(1.2);
  PlottingSysErrors(ProtonDefault,ProtonOriginal,graphSysErr,graphRelErr,"PP");
}


void ProtonLambdaSystematic(TString FileDefault, TString *FileSystematics)
{
  //Get the Original Binned CF
  Rebin=BinningPLDefault;
  TH1F* ProtonLambdaOriginal =
      GetSummedCF(&FileDefault,1,"",
                  PartNumberProton,PartNumberLambda,
                  PartNumberAntiProton,PartNumberAntiLambda);
  //Get the Default CF
  Rebin=BinningPLSyst;
  TH1F* ProtonLambdaDefault =
      GetSummedCF(&FileDefault,1,"",
                  PartNumberProton,PartNumberLambda,
                  PartNumberAntiProton,PartNumberAntiLambda);
  const int Variations=ProtonVarTaskMax-ProtonVarTaskMin+1+LambdaVarTaskMax-LambdaVarTaskMin+1;
  std::vector<const char *>ProtonLambdaVarName = ProtonVarName;
  ProtonLambdaVarName.insert(std::end(ProtonLambdaVarName),std::begin(LambdaVarName),std::end(LambdaVarName));
  if (!(Variations==ProtonLambdaVarName.size())) {
    std::cout<<
        "Number of Proton Lambda Variations are: "<<Variations<<
        " Number of Names given are: "<<ProtonLambdaVarName.size()<<
        " Please fix! \n";
    return;
  }
  TH1F* ProtonLambdaVariations[Variations];
  for (int iVar=0;iVar<Variations;++iVar) {
    //Get the summed CF for the variations
    ProtonLambdaVariations[iVar]=
        GetSummedCF(FileSystematics,nSysFiles,Form("%d",iVar+1),
                    PartNumberProton,PartNumberLambda,
                    PartNumberAntiProton,PartNumberAntiLambda);
  }
  TH1F* Ratios[Variations];
  RatioWithDefaultAndBarlow(ProtonLambdaDefault,ProtonLambdaVariations,Ratios,ProtonLambdaVarName,5,4,"PL");

  TH1F* ErrBudget[Variations];
  TH1F* AbsErr[Variations];
  PlotSystematicsPerVariation(ProtonLambdaDefault,Ratios,ErrBudget,AbsErr,ProtonLambdaVarName,5,4,"PL");
  //now that we have relative errors & aboslute values, we have to think about which
  //ones to combine since they are up/down variations:
  //    Cut         -- Variations (as stored in the arrays)
  //    pT                    0,1
  //    eta                   2,3
  //    nSigma                4,5
  //    TPC Cls               7,8
  //    V0 Pt                 9,10
  //    V0 Daugh Eta          14,15
  //After that loop over all k* bins of interest to combine them
  int nBins=ProtonLambdaDefault->GetXaxis()->FindBin(fitRange);
  const int EffVariations=12;
  double addMeUp[EffVariations];
  TGraphErrors *graphSysErr = new TGraphErrors();
  TGraph *graphRelErr = new TGraph(nBins);
  for (int iKstar=1;iKstar<nBins;++iKstar) {
    double SysErrTotal = 0;
    addMeUp[0] = (AbsErr[0]->GetBinContent(iKstar)+AbsErr[1]->GetBinContent(iKstar))/2.;//pT
    addMeUp[1] = (AbsErr[2]->GetBinContent(iKstar)+AbsErr[3]->GetBinContent(iKstar))/2.;//Eta
    addMeUp[2] = (AbsErr[4]->GetBinContent(iKstar)+AbsErr[5]->GetBinContent(iKstar))/2.;//nSigma
    addMeUp[3] = AbsErr[6]->GetBinContent(iKstar); // Filter Bit 96
    addMeUp[4] = (AbsErr[7]->GetBinContent(iKstar)+AbsErr[8]->GetBinContent(iKstar))/2.;//TPC Cls
    addMeUp[5] = (AbsErr[9]->GetBinContent(iKstar)+AbsErr[10]->GetBinContent(iKstar))/2.;//V0 Pt
    addMeUp[6] = AbsErr[11]->GetBinContent(iKstar); // CPA
    addMeUp[7] = AbsErr[12]->GetBinContent(iKstar); // Daughter PID
    addMeUp[8] = AbsErr[13]->GetBinContent(iKstar); // Daughter TPC Cls
    addMeUp[9] = (AbsErr[14]->GetBinContent(iKstar)+AbsErr[15]->GetBinContent(iKstar))/2.;//V0 Daugh Eta
    addMeUp[10] = AbsErr[15]->GetBinContent(iKstar); // Daughter DCA to V0 Vtx
    addMeUp[11] = AbsErr[16]->GetBinContent(iKstar); // Daughter DCA to Prim Vtx
    for (int iAdd=0;iAdd<EffVariations;++iAdd) {
      SysErrTotal+=(addMeUp[iAdd]*addMeUp[iAdd]);
    }
    SysErrTotal=TMath::Sqrt(SysErrTotal);
    double xValue=ProtonLambdaDefault->GetBinCenter(iKstar);
    graphSysErr->SetPoint(iKstar,xValue,ProtonLambdaDefault->GetBinContent(iKstar));
    graphSysErr->SetPointError(iKstar, ProtonLambdaDefault->GetBinWidth(1)/2.,SysErrTotal);
    double relErr = SysErrTotal/ProtonLambdaDefault->GetBinContent(iKstar);
    graphRelErr->SetPoint(iKstar,xValue,relErr);
  }
  ProtonLambdaDefault->SetTitle("p-#Lambda #oplus #bar{p}#bar{#Lambda}");
  ProtonLambdaDefault->GetXaxis()->SetRangeUser(0, 0.5);
  ProtonLambdaDefault->GetYaxis()->SetRangeUser(0.96, 1.02);
  ProtonLambdaDefault->GetYaxis()->SetTitleOffset(1.3);
  graphRelErr->GetYaxis()->SetRangeUser(0,0.04);
  PlottingSysErrors(ProtonLambdaDefault,ProtonLambdaOriginal,graphSysErr,graphRelErr,"PL");
}

void LambdaLambdaSystematic(TString FileDefault, TString *FileSystematics)
{
  //Get the Original Binned CF
  Rebin=BinningLLDefault;
  TH1F* LambdaOriginal =
      GetSummedCF(&FileDefault,1,"",
                  PartNumberLambda,PartNumberLambda,
                  PartNumberAntiLambda,PartNumberAntiLambda);
  //Get the Default CF
  Rebin=BinningLLSyst;
  TH1F* LambdaDefault =
      GetSummedCF(&FileDefault,1,"",
                  PartNumberLambda,PartNumberLambda,
                  PartNumberAntiLambda,PartNumberAntiLambda);
  const int Variations=LambdaVarTaskMax-LambdaVarTaskMin+1;
  //  std::vector<const char *>ProtonLambdaVarName = ProtonVarName;
  //  ProtonLambdaVarName.insert(std::end(ProtonLambdaVarName),std::begin(LambdaVarName),std::end(LambdaVarName));
  if (!(Variations==LambdaVarName.size())) {
    std::cout<<
        "Number of Lambda Variations are: "<<Variations<<
        " Number of Names given are: "<<LambdaVarName.size()<<
        " Please fix! \n";
    return;
  }
  TH1F* LambdaVariations[Variations];
  for (int iVar=LambdaVarTaskMin;iVar<LambdaVarTaskMax+1;++iVar) {
    //Get the summed CF for the variations
    LambdaVariations[iVar-LambdaVarTaskMin]=nullptr;
    LambdaVariations[iVar-LambdaVarTaskMin]=
        GetSummedCF(FileSystematics,nSysFiles,Form("%d",iVar),
                    PartNumberLambda,PartNumberLambda,
                    PartNumberAntiLambda,PartNumberAntiLambda);
    if (!LambdaVariations[iVar-LambdaVarTaskMin]) {
      std::cout<<
          " Variation:" << iVar-LambdaVarTaskMin << " Missing! \n";
      return;
    }
  }
  TH1F* Ratios[Variations];
  RatioWithDefaultAndBarlow(LambdaDefault,LambdaVariations,Ratios,LambdaVarName,3,3,"LL");

  TH1F* ErrBudget[Variations];
  TH1F* AbsErr[Variations];
  PlotSystematicsPerVariation(LambdaDefault,Ratios,ErrBudget,AbsErr,LambdaVarName,3,3,"LL");
  //now that we have relative errors & aboslute values, we have to think about which
  //ones to combine since they are up/down variations:
  //    Cut         -- Variations (as stored in the arrays)
  //    V0 Pt                 0,1
  //    V0 Daugh Eta          5,6
  //After that loop over all k* bins of interest to combine them
  int nBins=LambdaDefault->GetXaxis()->FindBin(fitRange);
  const int EffVariations=7;
  double addMeUp[EffVariations];
  TGraphErrors *graphSysErr = new TGraphErrors();
  TGraph *graphRelErr = new TGraph(nBins);
  for (int iKstar=1;iKstar<nBins;++iKstar) {
    double SysErrTotal = 0;
    addMeUp[0] = (AbsErr[0]->GetBinContent(iKstar)+AbsErr[1]->GetBinContent(iKstar))/2.;//V0 pT
    addMeUp[1] = AbsErr[2]->GetBinContent(iKstar); // CPA
    addMeUp[2] = AbsErr[3]->GetBinContent(iKstar); // Daughter PID
    addMeUp[3] = AbsErr[4]->GetBinContent(iKstar); // Daughter TPC Cls
    addMeUp[4] = (AbsErr[5]->GetBinContent(iKstar)+AbsErr[6]->GetBinContent(iKstar))/2.;//V0 Daugh Eta
    addMeUp[5] = AbsErr[7]->GetBinContent(iKstar); // Daughter DCA to V0 Vtx
    addMeUp[6] = AbsErr[8]->GetBinContent(iKstar); // Daughter DCA to Prim Vtx
    for (int iAdd=0;iAdd<EffVariations;++iAdd) {
      SysErrTotal+=(addMeUp[iAdd]*addMeUp[iAdd]);
    }
    SysErrTotal=TMath::Sqrt(SysErrTotal);
    double xValue=LambdaDefault->GetBinCenter(iKstar);
    graphSysErr->SetPoint(iKstar,xValue,LambdaDefault->GetBinContent(iKstar));
    graphSysErr->SetPointError(iKstar, LambdaDefault->GetBinWidth(1)/2.,SysErrTotal);
    double relErr = SysErrTotal/LambdaDefault->GetBinContent(iKstar);
    graphRelErr->SetPoint(iKstar,xValue,relErr);
  }
  LambdaDefault->SetTitle("#Lambda-#Lambda #oplus #bar{#Lambda}#bar{#Lambda}");
  LambdaDefault->GetXaxis()->SetRangeUser(0, 0.5);
  LambdaDefault->GetYaxis()->SetTitleOffset(1.3);
  LambdaDefault->GetYaxis()->SetRangeUser(0.93,1.03);
  graphRelErr->GetYaxis()->SetRangeUser(0,0.03);
  PlottingSysErrors(LambdaDefault,LambdaOriginal,graphSysErr,graphRelErr,"LL");
}

void ProtonXiSystematic(TString FileDefault, TString *FileSystematics)
{
  //Get the Original Binned CF
  Rebin=BinningPXiDefault;
  TH1F* ProtonXiOriginal =
      GetSummedCF(&FileDefault,1,"",
                  PartNumberProton,PartNumberXi,
                  PartNumberAntiProton,PartNumberAntiXi);
  //Get the Default CF
  Rebin=BinningPXiSyst;
  TH1F* ProtonXiDefault =
      GetSummedCF(&FileDefault,1,"",
                  PartNumberProton,PartNumberXi,
                  PartNumberAntiProton,PartNumberAntiXi);
  const int Variations=ProtonVarTaskMax-ProtonVarTaskMin+1+XiVarTaskMax-XiVarTaskMin+1;
  std::vector<const char *>ProtonXiVarName = ProtonVarName;
  ProtonXiVarName.insert(std::end(ProtonXiVarName),std::begin(XiVarName),std::end(XiVarName));
  if (!(Variations==ProtonXiVarName.size())) {
    std::cout<<
        "Number of Proton Xi Variations are: "<<Variations<<
        " Number of Names given are: "<<ProtonXiVarName.size()<<
        " Please fix! \n";
    return;
  }
  int iHist=0;
  TH1F* ProtonXiVariations[Variations];
  std::cout << Variations << std::endl;
  for (int iVar=ProtonVarTaskMin;iVar<=XiVarTaskMax;++iVar) {
    if ((iVar > ProtonVarTaskMax+1) && (iVar < XiVarTaskMin) ) {
      continue;
    }
    //Get the summed CF for the variations
    ProtonXiVariations[iHist]=
        GetSummedCF(FileSystematics,nSysFiles,Form("%d",iVar),
                    PartNumberProton,PartNumberXi,
                    PartNumberAntiProton,PartNumberAntiXi);
    iHist++;
  }
  TH1F* Ratios[Variations];
  RatioWithDefaultAndBarlow(ProtonXiDefault,ProtonXiVariations,Ratios,ProtonXiVarName,8,4,"PXi");

  TH1F* ErrBudget[Variations];
  TH1F* AbsErr[Variations];
  PlotSystematicsPerVariation(ProtonXiDefault,Ratios,ErrBudget,AbsErr,ProtonXiVarName,8,4,"PXi");
  //now that we have relative errors & aboslute values, we have to think about which
  //ones to combine since they are up/down variations:
  //    Cut         -- Variations (as stored in the arrays)
  //    pT                    0,1
  //    eta                   2,3
  //    nSigma                4,5
  //    TPC Cls               7,8
  //    Daugh DCA Vtx         9,10
  //    Bach Min Dist PV      11,12
  //    Xi CPA                13,14
  //    Xi TranRadius         15,16
  //    v0 Daugh DCA Vtx      17,18
  //    v0 CPA                19,20
  //    v0 TranRadius         21,22
  //    v0 Min Dist PV        23,24
  //    v0 Daugh Min Dist PV  25,26
  //    Xi Daugh Eta          27,28
  //    Xi Daugh PID          29,30
  //    Xi Pt                 31

  //After that loop over all k* bins of interest to combine them
  int nBins=ProtonXiDefault->GetXaxis()->FindBin(fitRange);
  const int EffVariations=17;
  double addMeUp[EffVariations];
  TGraphErrors *graphSysErr = new TGraphErrors();
  TGraph *graphRelErr = new TGraph(nBins);
  for (int iKstar=1;iKstar<nBins;++iKstar) {
    double SysErrTotal = 0;
    addMeUp[0] = (AbsErr[0]->GetBinContent(iKstar)+AbsErr[1]->GetBinContent(iKstar))/2.;//pT
    addMeUp[1] = (AbsErr[2]->GetBinContent(iKstar)+AbsErr[3]->GetBinContent(iKstar))/2.;//Eta
    addMeUp[2] = (AbsErr[4]->GetBinContent(iKstar)+AbsErr[5]->GetBinContent(iKstar))/2.;//nSigma
    addMeUp[3] = AbsErr[6]->GetBinContent(iKstar); // Filter Bit 96
    addMeUp[4] = (AbsErr[7]->GetBinContent(iKstar)+AbsErr[8]->GetBinContent(iKstar))/2.;//TPC Cls
    addMeUp[5] = (AbsErr[9]->GetBinContent(iKstar)+AbsErr[10]->GetBinContent(iKstar))/2.;//Daughter DCA to Casc Vtx
    addMeUp[6] = (AbsErr[11]->GetBinContent(iKstar)+AbsErr[12]->GetBinContent(iKstar))/2.;//Bach DCA to PV
    addMeUp[7] = (AbsErr[13]->GetBinContent(iKstar)+AbsErr[14]->GetBinContent(iKstar))/2.;//Xi CPA
    addMeUp[8] = (AbsErr[15]->GetBinContent(iKstar)+AbsErr[16]->GetBinContent(iKstar))/2.;//Xi Transverse Radius
    addMeUp[9] = (AbsErr[17]->GetBinContent(iKstar)+AbsErr[18]->GetBinContent(iKstar))/2.;//Daugh DCA to V0 Vtx
    addMeUp[10] = (AbsErr[19]->GetBinContent(iKstar)+AbsErr[20]->GetBinContent(iKstar))/2.;//V0 CPA
    addMeUp[11] = (AbsErr[21]->GetBinContent(iKstar)+AbsErr[22]->GetBinContent(iKstar))/2.;//V0 Transverse Radius
    addMeUp[12] = (AbsErr[23]->GetBinContent(iKstar)+AbsErr[24]->GetBinContent(iKstar))/2.;//V0 DCA to PV
    addMeUp[13] = (AbsErr[25]->GetBinContent(iKstar)+AbsErr[26]->GetBinContent(iKstar))/2.;//Daugh DCA to PV
    addMeUp[14] = (AbsErr[27]->GetBinContent(iKstar)+AbsErr[28]->GetBinContent(iKstar))/2.;//Xi Tracks Eta
    addMeUp[15] = (AbsErr[29]->GetBinContent(iKstar)+AbsErr[30]->GetBinContent(iKstar))/2.;//Xi Tracks PID
    addMeUp[16] = AbsErr[31]->GetBinContent(iKstar); // Xi Pt
    for (int iAdd=0;iAdd<EffVariations;++iAdd) {
      SysErrTotal+=(addMeUp[iAdd]*addMeUp[iAdd]);
    }
    SysErrTotal=TMath::Sqrt(SysErrTotal);
    double xValue=ProtonXiDefault->GetBinCenter(iKstar);
    graphSysErr->SetPoint(iKstar,xValue,ProtonXiDefault->GetBinContent(iKstar));
    graphSysErr->SetPointError(iKstar, ProtonXiDefault->GetBinWidth(1)/2.,SysErrTotal);
    double relErr = SysErrTotal/ProtonXiDefault->GetBinContent(iKstar);
    graphRelErr->SetPoint(iKstar,xValue,relErr);
  }
  ProtonXiDefault->SetTitle("p-#Xi #oplus #bar{p}#bar{#Xi}");
  ProtonXiDefault->GetXaxis()->SetRangeUser(0, 0.5);
  ProtonXiDefault->GetYaxis()->SetRangeUser(0.98, 1.12);
  ProtonXiDefault->GetYaxis()->SetTitleOffset(1.2);
  graphRelErr->GetYaxis()->SetRangeUser(0,0.06);
  PlottingSysErrors(ProtonXiDefault,ProtonXiOriginal,graphSysErr,graphRelErr,"PXi");
}

void TotalSystematics()
{
  TString FileDefault = "DataFullest/AnalysisResults.root";
  TString FileVariations[2] = {"Systematics/AnalysisResults_1.root","Systematics/AnalysisResults_2.root"};
  //  TString FileDefault = "AnalysisResults.root";
  //  TString FileVariations[1] = {"AnalysisResults.root"};
  ProtonProtonSystematic(FileDefault,FileVariations);
  ProtonLambdaSystematic(FileDefault,FileVariations);
  LambdaLambdaSystematic(FileDefault,FileVariations);
  ProtonXiSystematic(FileDefault,FileVariations);
  return;
}
