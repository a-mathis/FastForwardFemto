using namespace std;
void fit_forL_sum(){



//16MeV
//-----

 TFile *_file0 = TFile::Open("SYSTEMATICS_CutVarAdd_Global_Radius_Normal_16MeV.root");
 TGraph* pXimGraphDefault = (TGraph*) _file0->Get("pXimGraphDefault");

 //check the points of the model
 Double_t xval,yval;
 //for(int i=0;i<pXimGraphDefault->GetN();i++){
 // pXimGraphDefault->GetPoint(i,xval,yval);
 // cout<<" p"<<i<<"    y="<<yval<<"   x="<<xval<<endl;
 //}
 //p0    y=2.12357   x=16
 //p1    y=1.56054   x=32
 //p2    y=1.22629   x=48
 //p3    y=1.15273   x=64
 //p4    y=1.10786   x=80
 //p5    y=1.08411   x=96
 //p6    y=1.07307   x=112
 //p7    y=1.06325   x=128
 //p8    y=1.05872   x=144
 //p9    y=1.05501   x=160
 //p10    y=1.05246   x=176
 //p11    y=1.05098   x=192
 //p12    y=1.0497   x=208

 //histo model
 Double_t binwidth=16.;
 const Int_t nbk = 14;
 Double_t kbins[nbk+1];
 kbins[0]=0.;
 kbins[1]=8.;
 for(int i=2;i<=nbk;i++){
  kbins[i]=kbins[i-1]+binwidth;
//  cout<<" b"<<i<<"    val="<<kbins[i]<<endl;
 }
 TH1F* model = new TH1F("model","model",nbk,kbins);
 //check the points of the model
 for(int i=0;i<13;i++){
  pXimGraphDefault->GetPoint(i,xval,yval);
  model->SetBinContent(i+2,yval);
 }
 model->Draw("");
 pXimGraphDefault->Draw("same");


//fit the model with an exponential + pol1
TF1* f0 = new TF1("f0","pol1(0)+expo(2)",8.,220.);
//   1  p0           1.33166e+00   8.03323e-02   2.09270e-05   1.69976e-02
//   2  p1          -1.73434e-03   4.98923e-04   1.34600e-07   2.21926e+00
//   3  p2           1.48398e+00   7.57108e-01   4.03957e-04   7.00589e-04
//   4  p3          -7.09006e-02   3.06780e-02   1.43620e-05   2.56536e-02
f0->SetParameters(1.33166e+00,-1.73434e-03,1.48398e+00,-7.09006e-02);
f0->SetLineStyle(2);
f0->SetLineColor(3);
 model->Fit(f0,"R");


//data, add se and me for pp,ApAp, for both collision systems
//pPb
TFile *_file1 = TFile::Open("pPb_502TeV_CFOutput_pXi_Rebin_2.root");
//TH1F*  se_pp = (TH1F*) _file1->Get("ppSE_Shifted")->Clone("se_pp");
//TH1F*  se_ApAp = (TH1F*) _file1->Get("ApApSE_Shifted")->Clone("se_ApAp");
//TH1F*  me_pp = (TH1F*) _file1->Get("ppME_Shifted")->Clone("me_pp");
//TH1F*  me_ApAp = (TH1F*) _file1->Get("ApApME_Shifted")->Clone("me_ApAp");
    TH1F*  se_pp = (TH1F*) _file1->Get("ppSE")->Clone("se_pp");
    TH1F*  se_ApAp = (TH1F*) _file1->Get("ApApSE")->Clone("se_ApAp");
    TH1F*  me_pp = (TH1F*) _file1->Get("ppME")->Clone("me_pp");
    TH1F*  me_ApAp = (TH1F*) _file1->Get("ApApME")->Clone("me_ApAp");
se_pp->Add(se_ApAp);
me_pp->Add(me_ApAp);

//pp
TFile *_file2 = TFile::Open("pp_13TeV_CFOutput_pXi_Rebin_2.root");
//TH1F*  sepp_pp = (TH1F*) _file2->Get("ppSE_Shifted")->Clone("sepp_pp");
//TH1F*  sepp_ApAp = (TH1F*) _file2->Get("ApApSE_Shifted")->Clone("sepp_ApAp");
//TH1F*  mepp_pp = (TH1F*) _file2->Get("ppME_Shifted")->Clone("mepp_pp");
//TH1F*  mepp_ApAp = (TH1F*) _file2->Get("ApApME_Shifted")->Clone("mepp_ApAp");
   TH1F*  sepp_pp = (TH1F*) _file2->Get("ppSE")->Clone("sepp_pp");
   TH1F*  sepp_ApAp = (TH1F*) _file2->Get("ApApSE")->Clone("sepp_ApAp");
   TH1F*  mepp_pp = (TH1F*) _file2->Get("ppME")->Clone("mepp_pp");
   TH1F*  mepp_ApAp = (TH1F*) _file2->Get("ApApME")->Clone("mepp_ApAp");
se_pp->Add(sepp_pp);
se_pp->Add(sepp_ApAp);
me_pp->Add(mepp_pp);
me_pp->Add(mepp_ApAp);
 
//normalize so integral of k* in 200-400 is as in ME of gen
Double_t n1=200.0001;
Double_t n2=399.9999;
se_pp->Scale( me_pp->Integral(me_pp->FindBin(n1),me_pp->FindBin(n2))  /  se_pp->Integral(se_pp->FindBin(n1),se_pp->FindBin(n2))   );
TH1F* data = (TH1F*) se_pp->Clone("data");
data->Divide(me_pp);
data->Draw("");
data->SetLineColor(4);data->SetLineWidth(2);data->SetMarkerStyle(20);data->SetMarkerColor(4);
data->Draw("same");

//FIT 1: WITH SIGNAL (f1)
//now try to fit the data with the same exponential and changing only the pol1
cout<<endl<<" FIT WITH SIGNAL F1 "<<endl;
TF1* f1 = new TF1("f1","pol1(0)+expo(2)+gaus(4)",8.,220.);
f1->SetParameter(0,f0->GetParameter(0));
f1->SetParameter(1,f0->GetParameter(1));
f1->SetParameter(2,f0->GetParameter(2));
f1->SetParameter(3,f0->GetParameter(3));
f1->SetParameter(4,.01);
//f1->SetParameter(5,75.);
f1->SetParameter(5,80.);
f1->SetParameter(6,5.);

//   5  p4           1.61071e-01   9.00714e-02   1.03581e-04   3.31302e-03
//   6  p5           7.57188e+01   4.26082e+00   5.65461e-03   4.40653e-05
//   7  p6           6.15229e+00   3.33497e+00   3.78218e-03  -1.18676e-05


f1->SetParameter(4,1.61071e-01);
f1->SetParameter(5,7.57188e+01);
f1->SetParameter(6,6.15229e+00);


data->GetXaxis()->SetRangeUser(0.,200.);
 data->Fit(f1,"MR");
f0->Draw("same");

//FIT 2: no signal, exluding signal region (f2)
//bins 8, 9, 10, error infinite:
TH1F* dataFAKE = (TH1F*) data->Clone("dataFAKE");
dataFAKE->SetBinError(8+1,100000.);
dataFAKE->SetBinError(9+1,100000.);
dataFAKE->SetBinError(10+1,100000.);
TF1* f2 = new TF1("f2","pol1(0)+expo(2)",8.,220.);
f2->SetParameter(0,f0->GetParameter(0));
f2->SetParameter(1,f0->GetParameter(1));
f2->SetParameter(2,f0->GetParameter(2));
f2->SetParameter(3,f0->GetParameter(3));
f2->SetLineColor(1);
dataFAKE->SetLineColor(1);dataFAKE->SetLineWidth(1);
dataFAKE->Fit(f2,"MR");
data->Draw("");
dataFAKE->Draw("same");
data->Draw("same");
f0->Draw("same");
f1->Draw("same");

//FIT 3: no signal, whole data (f3)
cout<<endl<<" FIT whole data W/O SIGNAL F3 "<<endl;
TH1F* dataAGAIN = (TH1F*) data->Clone("dataAGAIN");
TF1* f3 = new TF1("f3","pol1(0)+expo(2)",8.,220.);
f3->SetParameter(0,f0->GetParameter(0));
f3->SetParameter(1,f0->GetParameter(1));
f3->SetParameter(2,f0->GetParameter(2));
f3->SetParameter(3,f0->GetParameter(3));
f3->SetLineColor(1);
f3->SetLineStyle(2);


gStyle->SetOptStat(0);
dataAGAIN->Fit(f3,"MR");
data->Draw("");
data->SetTitle("");
data->GetXaxis()->SetTitle("k* (MeV)");
data->GetYaxis()->SetTitle("C(k*)");
data->GetYaxis()->SetRangeUser(0.7,3.3);
//dataFAKE->Draw("same");
dataAGAIN->Draw("same");
data->Draw("same");
f0->Draw("same");
f2->Draw("same");
f1->Draw("same");



TLegend *leg2 = new TLegend( 0.45, 0.45, 0.85, 0.85);
leg2->AddEntry(data,"ALICE p-Pb 5TeV + pp 13TeV");
leg2->AddEntry(f1,"Expo+Pol1 + Gaus ");
leg2->AddEntry(f3,"Expo+Pol1 ");
leg2->AddEntry(f2,"Expo+Pol1 (exclude 64-88 MeV)");
leg2->AddEntry(f0,"Femtoscopic fit HAL QCD (p-Pb)");
leg2->SetLineColor(0);
leg2->Draw();




//evaluate simple nsigma
//data
//f3 (whole fit)
//f2 (excluding signal)
//f0 (CATS)
//compute single bin (bins 8,9,10)

int bin1=8+1;
int bin2=10+1;
Double_t data_val=0.;
Double_t data_e=0.;
Double_t model_val0=0.;
Double_t model_val2=0.;
Double_t model_val3=0.;
for (int i=bin1;i<=bin2;i++){
 data_val+=data->GetBinContent(i);
 data_e+=pow(data->GetBinError(i),2);
 model_val0+=f0->Eval(data->GetBinCenter(i));
 model_val2+=f2->Eval(data->GetBinCenter(i));
 model_val3+=f3->Eval(data->GetBinCenter(i));
}
data_e=sqrt(data_e);
//model_val0/=3.;
//model_val2/=3.;
//model_val3/=3.;

cout<<" data="<<data_val<<"+-"<<data_e<<endl;
cout<<" data="<<data_val<<"+-"<<data_e<<endl;

cout<<" MODEL 0 (cats) nsigma="<<fabs(data_val-model_val0)/data_e<<endl;
cout<<" MODEL 2 (excluding signal) nsigma="<<fabs(data_val-model_val2)/data_e<<endl;
cout<<" MODEL 3 (fit whole range) nsigma="<<fabs(data_val-model_val3)/data_e<<endl;

}
