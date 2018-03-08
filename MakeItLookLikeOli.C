void MakeItLookLikeOli(const char *fileName) {
  TList *TPdir=new TList();
//  TPdir->SetOwner();
  TPdir->SetName("TPdir_0");
  TString RelKNames[6]={"Proton","AntiProton","Lambda","AntiLambda","Xi","AntiXi"};
  TH1F* SEDist[6][6];
  TH1F* MEDist[6][6];
  TFile *file=TFile::Open(fileName);
  TDirectoryFile *dirResults=(TDirectoryFile*)(file->FindObjectAny("Results"));
  if (dirResults) {
    TList *tmp=dirResults->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *Results;
    dirResults->GetObject(name,Results);
    if (!Results) {
      std::cout << "No Results \n";
    } else {
      for (int iPart1=0;iPart1<6;++iPart1) {
        for (int iPart2=iPart1;iPart2<6;++iPart2) {
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
            TPdir->Add(SEDist[iPart1][iPart2]);

            TString MEName=Form("MEDist_Particle%i_Particle%i",iPart1,iPart2);
            MEDist[iPart1][iPart2]=(TH1F*)tmpFolder->FindObject(MEName.Data());
            MEDist[iPart1][iPart2]->SetDirectory(0);
            if (!MEDist[iPart1][iPart2]) {
              std::cout << SEName.Data() << " not Found\n";
              std::cout << SEName.Data() << " not Found\n";
            }
            TString NewMEName=Form("f%s%sRelKME",RelKNames[iPart1].Data(),RelKNames[iPart2].Data());
            MEDist[iPart1][iPart2]->SetNameTitle(NewMEName.Data(),NewMEName.Data());
            TPdir->Add(MEDist[iPart1][iPart2]);
          }
        }
      }
    }
  }
  TFile *output=new TFile("Renamed.root","RECREATE");
  TDirectoryFile* dirOutput=new TDirectoryFile("PWGCF_PLFemto_0","PWGCF_PLFemto_0");
//  TList* dirOutput=new TList();
//  dirOutput->SetName("PWGCF_PLFemto_0");
//  dirOutput->SetOwner();
  dirOutput->Add(TPdir);
  dirOutput->Write("PWGCF_PLFemto_0",TObject::kSingleKey);
//  TPdir->Write("TPdir_0");
}
