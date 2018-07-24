void theUnholySigma(const char* inputCF, const char * inputSys, const char* FitResult) {
  TFile* ExpCF =  TFile::Open(inputCF);
  TFile* ExpSys =  TFile::Open(inputSys);
  TFile* CatsFile = TFile::Open(FitResult);

  const char* cfName = "hCkTotNormWeight_Shifted";
  const char* sysName = "";

  return;
}
