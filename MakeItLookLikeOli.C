void MakeItLookLikeOli(const char *fileName) {
  TList *TPdir=new TList();
  TPdir->SetName("TPdir_0");
  TList *PIDdir=new TList();
  PIDdir->SetName("PIDdir_0");
  TList *AliEventCuts;
  TList *SPdir=new TList();
  SPdir->SetName("SPdir_0");

  TString RelKNames[6]={"Proton","AntiProton","Lambda","AntiLambda","Xi","AntiXi"};
  TH1F* SEDist[6][6];
  TH1F* MEDist[6][6];
  TFile *file=TFile::Open(fileName);

  TDirectoryFile *dirQA=(TDirectoryFile*)(file->FindObjectAny("QA"));
  if (dirQA) {
    TList *tmp=dirQA->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *QA;
    dirQA->GetObject(name,QA);
    if (QA) {
      AliEventCuts=(TList*)(QA->FindObject("AliEventCuts"));
      AliEventCuts->SetName("AliEventCuts_0");
      if (!AliEventCuts) {
        std::cout << "No Event Cuts \n";
      }
      TString v0TSharedName[4]={"fNV0protonSharedTracks","fNAntiV0AntiprotonSharedTracks","fNXiSharedTracks","fNAntiXiAntiprotonSharedTracks"};
      TString v0v0SharedName[4]={"fNV0TrackSharing","fNAntiV0TrackSharing","fNXiTrackSharing","fNAntiXiTrackSharing"};

      TList *PairCleaner = (TList*)QA->FindObject("PairCleaner");
      if (PairCleaner) {
        for (int i=0; i<4;++i) {
          TH1F* DaugTrack=(TH1F*)PairCleaner->FindObject(Form("DaugthersSharedTracks_%i",i));
          DaugTrack->SetName(v0TSharedName[i].Data());
          SPdir->Add(DaugTrack);
          TH1F* DaugDaug=(TH1F*)PairCleaner->FindObject(Form("DaugthersSharedDaughters_%i",i));
          DaugDaug->SetName(v0v0SharedName[i].Data());
          SPdir->Add(DaugDaug);
        }
      } else {
        std::cout << "No Pair Cleaner \n";
      }
    } else {
      std::cout << "No QA List \n";
    }
  } else {
    std::cout << "No Dir QA!\n";
  }
  TDirectoryFile *dirEvent=(TDirectoryFile*)(file->FindObjectAny("EvtCuts"));
  if (dirEvent) {
    TList *tmp=dirEvent->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *EvtCuts;
    dirEvent->GetObject(name,EvtCuts);
    if (EvtCuts) {
      TList *after = (TList*)EvtCuts->FindObject("after");
      if (after) {
        TH1F* MultSPD=(TH1F*)after->FindObject("MultiplicitySPD_after");
        MultSPD->SetName("fNTracklets");
        SPdir->Add(MultSPD);
      } else {
        std::cout << "No Evt Cuts After \n";
      }
    } else {
      std::cout << "No Evt Cuts \n";
    }
  }
  TDirectoryFile *dirProtonCuts=(TDirectoryFile*)(file->FindObjectAny("TrackCuts"));
  if (dirProtonCuts) {
    TList *tmp=dirProtonCuts->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *TrackCuts;
    dirProtonCuts->GetObject(name,TrackCuts);
    if (TrackCuts) {
      TH2F *dcaXYProton=(TH2F*)TrackCuts->FindObject("DCAXYPtBinningTot");
      TH1F *projDCAXY=(TH1F*)dcaXYProton->ProjectionY("fProtonDCAxy");
      SPdir->Add(projDCAXY);
      //we need to trick this into a TH2F* for every pt bin
      for (int i =0; i<20;++i) {
        TString DCAXYOliHistName=Form("fProtonDCAxyDCAzPt%i",i);
        TH2F *outputHist=new TH2F(DCAXYOliHistName.Data(),DCAXYOliHistName.Data(),500,-5,5,500,-5,5);
        for (int iBin=1;iBin<=500;++iBin) {
          outputHist->SetBinContent(iBin,outputHist->GetYaxis()->FindBin(0.),dcaXYProton->GetBinContent(i+1,iBin));
          SPdir->Add(outputHist);
        }
      }
      TList *after=(TList*)TrackCuts->FindObject("after");
      if (after) {
        TH2F *AfterDCAXY=(TH2F*)after->FindObject("DCAXY_after");
        TH1F *DCAXY1D=(TH1F*)AfterDCAXY->ProjectionY("fProtonDCAxyCutz");
        SPdir->Add(DCAXY1D);

        TH1F *AfterPt=(TH1F*)after->FindObject("pTDist_after");
        AfterPt->SetName("fProtonPt");
        SPdir->Add(AfterPt);

        TH1F *AfterPhi=(TH1F*)after->FindObject("phiDist_after");
        AfterPhi->SetName("fProtonPhi");
        SPdir->Add(AfterPhi);

        TH1F *AfterEta=(TH1F*)after->FindObject("EtaDist_after");
        AfterEta->SetName("fProtonEta");
        SPdir->Add(AfterEta);

        TH1F *AfterSigmaTPC=(TH1F*)after->FindObject("NSigTPC_after");
        AfterSigmaTPC->SetName("fProtonNSigmaTPC");
        PIDdir->Add(AfterSigmaTPC);

        TH1F *AfterSigmaTOF=(TH1F*)after->FindObject("NSigTOF_after");
        AfterSigmaTOF->SetName("fProtonNSigmaCombined");
        PIDdir->Add(AfterSigmaTOF);

      } else {
        std::cout << "No After Track Cuts \n";
      }
    } else {
      std::cout << "No Track Cuts \n";
    }
  } else {
    std::cout << "No Proton Cuts \n";
  }

  TDirectoryFile *dirv0Cuts=(TDirectoryFile*)(file->FindObjectAny("v0Cuts"));
  if (dirv0Cuts) {
    TList *tmp=dirv0Cuts->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *v0CutList;
    dirv0Cuts->GetObject(name,v0CutList);
    if (v0CutList) {
      TList *v0Cuts = (TList*)v0CutList->FindObject("v0Cuts");
      if (v0Cuts) {
        TH1F *AfterInvMassLambda=(TH1F*)v0Cuts->FindObject("InvMasswithCuts");
        AfterInvMassLambda->SetName("fInvMassLambdawCuts");
        SPdir->Add(AfterInvMassLambda);

        TH1F *InvMassKaon=(TH1F*)v0Cuts->FindObject("InvMassKaon");
        InvMassKaon->SetName("fInvMassMissIDK0s");
        SPdir->Add(InvMassKaon);

        TH1F *fakeRejKaon= new TH1F("fInvMassMissIDK0swCuts","fInvMassMissIDK0swCuts",
                                    InvMassKaon->GetNbinsX(),
                                    InvMassKaon->GetXaxis()->GetBinCenter(1),
                                    InvMassKaon->GetXaxis()->GetBinCenter(InvMassKaon->GetNbinsX()));

        SPdir->Add(fakeRejKaon);

        TH2F* CPALamDataTot=(TH2F*)v0Cuts->FindObject("CPAPtBinsTot");
        for (int i=0;i<CPALamDataTot->GetXaxis()->GetNbins();++i) {
          SPdir->Add(CPALamDataTot->ProjectionY(Form("fLambdaCPAPtBin%i",i),i+1,i+2));
        }
        TList *after=(TList*)v0Cuts->FindObject("after");
        if (after) {
          TH1F *AfterPt=(TH1F*)after->FindObject("pTDist_after");
          AfterPt->SetName("fLambdaPt");
          SPdir->Add(AfterPt);

          TH1F *AfterPhi=(TH1F*)after->FindObject("PhiDist_after");
          AfterPhi->SetName("fLambdaPhi");
          SPdir->Add(AfterPhi);

          TH1F *AfterEta=(TH1F*)after->FindObject("EtaDist_after");
          AfterEta->SetName("fLambdaEta");
          SPdir->Add(AfterEta);

          TH1F *AfterDCADaughterTracks=(TH1F*)after->FindObject("DCADauToVtx_after");
          AfterDCADaughterTracks->SetName("fLambdaDCADaughterTracks");
          SPdir->Add(AfterDCADaughterTracks);

          TH1F *AfterPDaugPrimVtx=(TH1F*)after->FindObject("DCADauPToPV_after");
          AfterPDaugPrimVtx->SetName("fLambdaDCAPosdaughPrimVertex");
          SPdir->Add(AfterPDaugPrimVtx);

          TH1F *AfterNDaugPrimVtx=(TH1F*)after->FindObject("DCADauNToPV_after");
          AfterNDaugPrimVtx->SetName("fLambdaDCANegdaughPrimVertex");
          SPdir->Add(AfterNDaugPrimVtx);

          TH1F *AfterTransRadius=(TH1F*)after->FindObject("TransverseRadius_after");
          AfterTransRadius->SetName("fLambdaTransverseRadius");
          SPdir->Add(AfterTransRadius);

        } else {
          std::cout << "No After v0 Cuts \n";
        }
      } else {
        std::cout << "No v0 Cuts \n";
      }
    } else {
      std::cout << "No v0 Cut List \n";
    }
  }

  TDirectoryFile *dirAntiv0Cuts=(TDirectoryFile*)(file->FindObjectAny("Antiv0Cuts"));
  if (dirAntiv0Cuts) {
    TList *tmp=dirAntiv0Cuts->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *Antiv0CutList;
    dirAntiv0Cuts->GetObject(name,Antiv0CutList);
    if (Antiv0CutList) {
      TList *Antiv0Cuts=(TList*)Antiv0CutList->FindObject("v0Cuts");
      if (Antiv0Cuts) {
        TH1F *AfterInvMassAntiLambda=(TH1F*)Antiv0Cuts->FindObject("InvMasswithCuts");
        AfterInvMassAntiLambda->SetName("fInvMassAntiLambdawCuts");
        SPdir->Add(AfterInvMassAntiLambda);
      }
    }
  }

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
  dirOutput->Add(TPdir);
  dirOutput->Add(PIDdir);
  dirOutput->Add(AliEventCuts);
  dirOutput->Add(SPdir);
  dirOutput->Write("PWGCF_PLFemto_0",TObject::kSingleKey);


  //for the systematics:

  for (int i=1;i<=30;++i) {
    TString dirName=Form("PWGCF_PLFemto_%i",i);
    TDirectoryFile* sysdirOutput=new TDirectoryFile(dirName.Data(),dirName.Data());

    TList *sysTPList=new TList();
    sysTPList->SetName(Form("TPdir_%i",i));
    TString resultsName=Form("Results_%i",i);
    TDirectoryFile *sysdirResults=(TDirectoryFile*)(file->FindObjectAny(resultsName.Data()));
    if (sysdirResults) {
      TList *tmp=sysdirResults->GetListOfKeys();
      TString name=tmp->At(0)->GetName();
      TList *Results;
      sysdirResults->GetObject(name,Results);
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
              sysTPList->Add(SEDist[iPart1][iPart2]);

              TString MEName=Form("MEDist_Particle%i_Particle%i",iPart1,iPart2);
              MEDist[iPart1][iPart2]=(TH1F*)tmpFolder->FindObject(MEName.Data());
              MEDist[iPart1][iPart2]->SetDirectory(0);
              if (!MEDist[iPart1][iPart2]) {
                std::cout << SEName.Data() << " not Found\n";
                std::cout << SEName.Data() << " not Found\n";
              }
              TString NewMEName=Form("f%s%sRelKME",RelKNames[iPart1].Data(),RelKNames[iPart2].Data());
              MEDist[iPart1][iPart2]->SetNameTitle(NewMEName.Data(),NewMEName.Data());
              sysTPList->Add(MEDist[iPart1][iPart2]);
            }
          }
        }
      }
    }
    sysdirOutput->Add(sysTPList);
    sysdirOutput->Write("PWGCF_PLFemto_0",TObject::kSingleKey);
  }

}
