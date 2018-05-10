void Friendship() {
  // Data
  TFile *file2 = TFile::Open("CombinedMap_ForIlya.root");
  TH2F *sigma, *ledniSucks;
  ledniSucks = (TH2F*)file2->Get("hGlobMinCk");
  int nBinsX=ledniSucks->GetXaxis()->GetNbins();
  float xMin=ledniSucks->GetXaxis()->GetBinLowEdge(1);
  float xMax=ledniSucks->GetXaxis()->GetBinUpEdge(nBinsX);
  int nBinsY=ledniSucks->GetYaxis()->GetNbins();
  float yMin=ledniSucks->GetYaxis()->GetBinLowEdge(1);
  float yMax=ledniSucks->GetYaxis()->GetBinUpEdge(nBinsY);
  std::cout << nBinsX << '\t' << xMin << '\t' << xMax << '\t' <<
      nBinsY << '\t' << yMin << '\t' << yMax << std::endl;
  TH2F* invLedniSucks=new TH2F("invLedniSucks","invLedniSucks",nBinsX,xMin,xMax
                               ,nBinsY,yMin,yMax);
  for (int iXbin=0;iXbin<=nBinsX;++iXbin) {
    for (int iYbin=0;iYbin<=nBinsY;++iYbin) {
      if (ledniSucks->GetBinContent(iXbin,iYbin)==0) {
        invLedniSucks->SetBinContent(iXbin,iYbin,15);
      } else {
        invLedniSucks->SetBinContent(iXbin,iYbin,0);
      }
    }
  }

  TFile *output=new TFile("Friendship.root","RECREATE");
  invLedniSucks->Write();
}
