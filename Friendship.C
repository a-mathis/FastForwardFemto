void Friendship(int iMode=0) {
  // Data
  TFile *file2 = TFile::Open("CombinedMap.root");
  TH2F *sigma, *ledniSucks;
  if (iMode==1)  ledniSucks = (TH2F*)file2->Get("hMinCk_0");
  if (iMode==2)  ledniSucks = (TH2F*)file2->Get("hMinCk_1");
  if (iMode==3)  ledniSucks = (TH2F*)file2->Get("hMinCk_2");
  if (iMode==4)  ledniSucks = (TH2F*)file2->Get("hGlobMinCk");
  int nBinsX=ledniSucks->GetXaxis()->GetNbins();
  float xMin=ledniSucks->GetXaxis()->GetBinLowEdge(1);
  float xMax=ledniSucks->GetXaxis()->GetBinUpEdge(nBinsX);
  int nBinsY=ledniSucks->GetYaxis()->GetNbins();
  float yMin=ledniSucks->GetYaxis()->GetBinLowEdge(1);
  float yMax=ledniSucks->GetYaxis()->GetBinUpEdge(nBinsY);
  TH2F* invLedniSucks=new TH2F("invLedniSucks","invLedniSucks",nBinsX,xMin,xMax
                               ,nBinsY,yMin,yMax);
  TGraphAsymmErrors *lednickyReallySucks = new TGraphAsymmErrors();
  int iPoint = 0;
  for (int iXbin=1;iXbin<=nBinsX;++iXbin) {
    float yBinMin = 0;
    float yBinMax = yMax;
    for (int iYbin=1;iYbin<=nBinsY;++iYbin) {
      if (ledniSucks->GetBinContent(iXbin,iYbin)==0) {
        invLedniSucks->SetBinContent(iXbin,iYbin,15);
        if (yBinMin == 0) {
          yBinMin = ledniSucks->GetYaxis()->GetBinLowEdge(iYbin);
        }
      } else {
        invLedniSucks->SetBinContent(iXbin,iYbin,0);
      }
    }
    if (yBinMin == 0) {
      break;
    }
    float median = yBinMin+yBinMax;
    median/=2;
    float errLow = yBinMax - median;
    float errUp = median - yBinMin;
    float xValue = invLedniSucks->GetXaxis()->GetBinLowEdge(iXbin);
    lednickyReallySucks->SetPoint(iPoint,xValue,median);
    lednickyReallySucks->SetPointError(iPoint++,0,0,errLow,errUp);
  }
  TFile *output=new TFile(Form("Friendship_%d.root",iMode),"RECREATE");
  invLedniSucks->Write();
  lednickyReallySucks->Write("myfriend");
}
