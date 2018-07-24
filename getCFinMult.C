Color_t fFillColors[8] = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9,
    kOrange-9, kCyan-8, kYellow-7};
Color_t fColors[8]  = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1,
    kOrange-1, kCyan+2, kYellow+2};
Style_t fMarkers[10] = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare,
    kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHisto(TH1 *histo, int marker, int color)
{
  histo->GetXaxis()->SetLabelSize(0.06);
  histo->GetXaxis()->SetTitleSize(0.07);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.0);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.06);
  histo->GetYaxis()->SetTitleSize(0.07);
  histo->GetYaxis()->SetLabelOffset(0.02);
  histo->GetYaxis()->SetTitleOffset(1.0);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* Calculate_CF(TH1F* histRE_relK,TH1F* histME_relK, TString CFname,Double_t normleft,Double_t normright, const char* folder, float spinningDepth)
{
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  if(strcmp(folder, "") == 0) {
    Double_t norm_relK = histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) / histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright));
    std::cout << histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) << '\t' <<
        histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright)) << '\t' <<
        norm_relK << std::endl;
    Hist_CF->Divide(histRE_relK,histME_relK,1,norm_relK);
  }
  else {
    histME_relK->Scale(1.f/spinningDepth);
    Hist_CF->Divide(histRE_relK,histME_relK);
  }

  return Hist_CF;
}

void getCFinMult(const char* filename) {
  const float normleft = 0.2;
  const float normright = 0.4;

  TFile *input=TFile::Open(filename);
  TString FolderName = "MBResults";
  TH2F* CentVsRefMult=nullptr;
  TDirectoryFile *dir=(TDirectoryFile*)(input->FindObjectAny(FolderName.Data()));
  TString centOrMult = "Cent";
  std::vector<int> binLimits;
  if (centOrMult == "Mult") {
    binLimits.push_back(-1);
    binLimits.push_back(11);
    binLimits.push_back(15);
    binLimits.push_back(27);
  } else if (centOrMult == "Cent") {
    binLimits.push_back(-1);
    binLimits.push_back(7);
    binLimits.push_back(25);
    binLimits.push_back(100);
  } else {
    std::cout << "No defined estimator \n";
    return;
  }
  const int nMultBins = 3;
  TH2F* SEIntPart;
  TH2F* MEIntPart;
  TH1F* SEPartBin[nMultBins];
  TH1F* MEPartBin[nMultBins];
  TH1F* CFPartBin[nMultBins];

  if (dir) {
    TList *tmp=dir->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *OutList;
    dir->GetObject(name,OutList);
    if (!OutList) {
      std::cout << "No QAResult\n";
    } else {
      TList* CFPartList = (TList*)OutList->FindObject("Particle0_Particle0");
      SEIntPart =
          (TH2F*)CFPartList->FindObject(Form("SE%sDist_Particle0_Particle0",
                                             centOrMult.Data()));
      if (!SEIntPart) {
        std::cout << "No SEIntPart found \n";
        return;
      }

      MEIntPart =
          (TH2F*)CFPartList->FindObject(Form("ME%sDist_Particle0_Particle0",
                                             centOrMult.Data()));
      if (!MEIntPart) {
        std::cout << "No MEIntPart found \n";
        return;
      }
      for (int iMult = 0; iMult < nMultBins; ++iMult) {
        int firstBin=SEIntPart->GetYaxis()->FindBin(binLimits.at(iMult)+1);
        int lastBin=SEIntPart->GetYaxis()->FindBin(binLimits.at(iMult+1));
        std::cout << firstBin << '\t' << lastBin << std::endl;
        TString HistNameSE = Form("SE_Mult_%d",iMult);
        TString HistNameME = Form("ME_Mult_%d",iMult);
        TString HistNameCF = Form("CF_Mult_%d",iMult);
        SEPartBin[iMult] = (TH1F*)SEIntPart->ProjectionX(HistNameSE.Data(),firstBin,lastBin);
        MEPartBin[iMult] = (TH1F*)MEIntPart->ProjectionX(HistNameME.Data(),firstBin,lastBin);

        CFPartBin[iMult]=Calculate_CF(SEPartBin[iMult],MEPartBin[iMult],HistNameCF.Data(),normleft,normright,"",0);
      }
    }
  } else {
    std::cout << "Dir not found \n" ;
    return;
  }

  TCanvas *MultDist= new TCanvas("MultDist","MultDist",0,0,1500,1100);
  SetStyleHisto(SEPartBin[0],1,1);
  SEPartBin[0]->DrawCopy();
  SetStyleHisto(SEPartBin[1],1,2);
  SEPartBin[1]->DrawCopy("Same");
  SetStyleHisto(SEPartBin[2],1,3);
  SEPartBin[2]->DrawCopy("Same");

  TCanvas *MultDist2= new TCanvas("MultDist2","MultDist2",0,0,1500,1100);
  SetStyleHisto(CFPartBin[0],1,1);
  CFPartBin[0]->GetXaxis()->SetRangeUser(0,0.2);
  CFPartBin[0]->DrawCopy();
  SetStyleHisto(CFPartBin[1],1,2);
  CFPartBin[1]->DrawCopy("Same");
  SetStyleHisto(CFPartBin[2],1,3);
  CFPartBin[2]->DrawCopy("Same");



  TFile *output = TFile::Open(Form("%sBinnedProtonProtonCF.root",centOrMult.Data()),"RECREATE");

  SEPartBin[0]->Write();
  MEPartBin[0]->Write();
  CFPartBin[0]->Write();
  SEPartBin[1]->Write();
  MEPartBin[1]->Write();
  CFPartBin[1]->Write();
  SEPartBin[2]->Write();
  MEPartBin[2]->Write();
  CFPartBin[2]->Write();

}
