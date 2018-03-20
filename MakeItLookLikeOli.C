void MakeItLookLikeOli(const char *fileName,const char *prefix,bool isMC) {
  TList *TPdir=new TList();
  TPdir->SetName("TPdir_0");
  TList *PIDdir=new TList();
  PIDdir->SetName("PIDdir_0");
  TList *AliEventCuts;
  TList *SPdir=new TList();
  SPdir->SetName("SPdir_0");
  TList *outResults=new TList();
  outResults->SetName("Results");
  TString RelKNames[6]={"Proton","AntiProton","Lambda","AntiLambda","Xi","AntiXi"};
  TH1F* SEDist[6][6];
  TH1F* MEDist[6][6];
  TFile *file=TFile::Open(fileName);
  TString QAName="";
  QAName+=prefix;
  QAName+="QA";

  TDirectoryFile *dirQA=(TDirectoryFile*)(file->FindObjectAny(QAName.Data()));
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
  TString EvtCutsName="";
  EvtCutsName+=prefix;
  EvtCutsName+="EvtCuts";

  TDirectoryFile *dirEvent=(TDirectoryFile*)(file->FindObjectAny(EvtCutsName.Data()));
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
  TString TrackCutsName="";
  TrackCutsName+=prefix;
  TrackCutsName+="TrackCuts";

  TDirectoryFile *dirProtonCuts=(TDirectoryFile*)(file->FindObjectAny(TrackCutsName.Data()));
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
        TString DCAXYOliHistName;
        DCAXYOliHistName=Form("fProtonDCAxyDCAzPt%i",i);
        TH2F *outputHist=new TH2F(DCAXYOliHistName.Data(),DCAXYOliHistName.Data(),500,-5,5,500,-5,5);
        for (int iBin=1;iBin<=500;++iBin) {
          outputHist->SetBinContent(iBin,outputHist->GetYaxis()->FindBin(0.),dcaXYProton->GetBinContent(i+1,iBin));
        }
        SPdir->Add(outputHist);
        if (isMC) {
          DCAXYOliHistName=Form("fProtonDCAxyDCAzMCCase0PtBin%i",i);
          outputHist->SetName(DCAXYOliHistName.Data());
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
  TString AntiTrackCutsName="";
  AntiTrackCutsName+=prefix;
  AntiTrackCutsName+="AntiTrackCuts";

  TDirectoryFile *dirAntiProtonCuts=(TDirectoryFile*)(file->FindObjectAny(AntiTrackCutsName.Data()));
  if (dirAntiProtonCuts) {
    TList *tmp=dirAntiProtonCuts->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *AntiTrackCuts;
    dirAntiProtonCuts->GetObject(name,AntiTrackCuts);
    if (AntiTrackCuts) {
      TH2F *dcaXYProton=(TH2F*)AntiTrackCuts->FindObject("DCAXYPtBinningTot");
      TH1F *projDCAXY=(TH1F*)dcaXYProton->ProjectionY("fAntiProtonDCAxy");
      SPdir->Add(projDCAXY);
      //we need to trick this into a TH2F* for every pt bin
      for (int i =0; i<20;++i) {
        TString DCAXYOliHistName=Form("fAntiProtonDCAxyDCAzPt%i",i);
        TH2F *outputHist=new TH2F(DCAXYOliHistName.Data(),DCAXYOliHistName.Data(),500,-5,5,500,-5,5);
        for (int iBin=1;iBin<=500;++iBin) {
          outputHist->SetBinContent(iBin,outputHist->GetYaxis()->FindBin(0.),dcaXYProton->GetBinContent(i+1,iBin));
          SPdir->Add(outputHist);
        }
      }
      TList *after=(TList*)AntiTrackCuts->FindObject("after");
      if (after) {
        TH2F *AfterDCAXY=(TH2F*)after->FindObject("DCAXY_after");
        TH1F *DCAXY1D=(TH1F*)AfterDCAXY->ProjectionY("fAntiProtonDCAxyCutz");
        SPdir->Add(DCAXY1D);

        TH1F *AfterPt=(TH1F*)after->FindObject("pTDist_after");
        AfterPt->SetName("fAntiProtonPt");
        SPdir->Add(AfterPt);

        TH1F *AfterPhi=(TH1F*)after->FindObject("phiDist_after");
        AfterPhi->SetName("fAntiProtonPhi");
        SPdir->Add(AfterPhi);

        TH1F *AfterEta=(TH1F*)after->FindObject("EtaDist_after");
        AfterEta->SetName("fAntiProtonEta");
        SPdir->Add(AfterEta);

        TH1F *AfterSigmaTPC=(TH1F*)after->FindObject("NSigTPC_after");
        AfterSigmaTPC->SetName("fAntiProtonNSigmaTPC");
        PIDdir->Add(AfterSigmaTPC);

        TH1F *AfterSigmaTOF=(TH1F*)after->FindObject("NSigTOF_after");
        AfterSigmaTOF->SetName("fAntiProtonNSigmaCombined");
        PIDdir->Add(AfterSigmaTOF);

      } else {
        std::cout << "No After Anti  Track Cuts \n";
      }
    } else {
      std::cout << "No Anti Track Cuts \n";
    }
  } else {
    std::cout << "No AntiProton Cuts \n";
  }
  TString v0CutsName="";
  v0CutsName+=prefix;
  v0CutsName+="v0Cuts";

  TDirectoryFile *dirv0Cuts=(TDirectoryFile*)(file->FindObjectAny(v0CutsName.Data()));
  if (dirv0Cuts) {
    TList *tmp=dirv0Cuts->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *v0CutList;
    dirv0Cuts->GetObject(name,v0CutList);
    if (v0CutList) {
      TList *v0Cuts = (TList*)v0CutList->FindObject("v0Cuts");
      if (v0Cuts) {
        TH2F *InvMassPt=(TH2F*)v0Cuts->FindObject("InvMassPt");
        for (int i =0; i<8;++i) {
          TString namePtLambda=Form("fInvMassLambdawCutsPtBin%i",i);
          TH1F* proj=(TH1F*)InvMassPt->ProjectionY(namePtLambda.Data(),i+1,i+2);
          SPdir->Add(proj);
        }
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
  TString Antiv0CutsName="";
  Antiv0CutsName+=prefix;
  Antiv0CutsName+="EvtCuts";

  TDirectoryFile *dirAntiv0Cuts=(TDirectoryFile*)(file->FindObjectAny(Antiv0CutsName.Data()));
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
  int iCounter=0;
  TList *OuttmpList[2];
  OuttmpList[0]=new TList();
  OuttmpList[0]->SetName(Form("Particle%i_Particle%i",0,4));
  outResults->Add(OuttmpList[0]);
  OuttmpList[1]=new TList();
  OuttmpList[1]->SetName(Form("Particle%i_Particle%i",1,5));
  outResults->Add(OuttmpList[1]);

  TString ResultsName="";
  ResultsName+=prefix;
  ResultsName+="Results";


  TDirectoryFile *dirResults=(TDirectoryFile*)(file->FindObjectAny(ResultsName.Data()));
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
            if ((iPart1==0&&iPart2==4)||(iPart1==1&&iPart2==5)) {
              OuttmpList[iCounter]->Add(SEDist[iPart1][iPart2]->Clone(SEName.Data()));
              OuttmpList[iCounter]->Add(MEDist[iPart1][iPart2]->Clone(MEName.Data()));
              iCounter++;
            }
          }
        }
      }
    }
  }


  //Some MC Stuff

  if (isMC) {
    TString MCTrackName="";
    MCTrackName+=prefix;
    MCTrackName+="TrkCutsMC";


    TDirectoryFile *dirMCTrkCuts=(TDirectoryFile*)(file->FindObjectAny(MCTrackName.Data()));
    if (dirMCTrkCuts) {
      TList *tmp=dirMCTrkCuts->GetListOfKeys();
      TString name=tmp->At(0)->GetName();
      TList *mcTrkCuts;
      dirMCTrkCuts->GetObject(name,mcTrkCuts);
      if (!mcTrkCuts) {
        std::cout << "No MCTrkCuts \n";
      } else {
        TH1F *IDed=(TH1F*)mcTrkCuts->FindObject("IdentPartPt");
        IDed->SetName("fProtonsTotallyIdentified");
        PIDdir->Add(IDed);

        TH1F *CorrIDed=(TH1F*)mcTrkCuts->FindObject("CorrParPt");
        CorrIDed->SetName("fProtonsCorrectlyIdentified");
        PIDdir->Add(CorrIDed);

        TList *DCAPtBinning=(TList*)mcTrkCuts->FindObject("DCAPtBinning");
        if (DCAPtBinning) {
          TH2F *dcaXYProtonPri=(TH2F*)DCAPtBinning->FindObject("DCAPtBinningPri");
          TH2F *dcaXYProtonMat=(TH2F*)DCAPtBinning->FindObject("DCAPtBinningMat");
          TH2F *dcaXYProtonSec=(TH2F*)DCAPtBinning->FindObject("DCAPtBinningSec");
          TH2F *dcaXYProtonSecLam=(TH2F*)DCAPtBinning->FindObject("DCAPtBinningSecLam");
          TH2F *dcaXYProtonSecSig=(TH2F*)DCAPtBinning->FindObject("DCAPtBinningSecSig");
          //we need to trick this into a TH2F* for every pt bin
          for (int i =0; i<20;++i) {
            TString DCAXYOliHistNameCasePri=Form("fProtonDCAxyDCAzMCCase1PtBin%i",i);
            TString DCAXYOliHistNameCaseSec=Form("fProtonDCAxyDCAzMCCase2PtBin%i",i);
            TString DCAXYOliHistNameCaseMat=Form("fProtonDCAxyDCAzMCCase3PtBin%i",i);
            TString DCAXYOliHistNameCaseSecLam=Form("fProtonDCAxyDCAzMCPtBinLambdaPtBin%i",i);
            TString DCAXYOliHistNameCaseSecSig=Form("fProtonDCAxyDCAzMCPtBinSigmaPtBin%i",i);

            TH2F *outputHistPri=new TH2F(DCAXYOliHistNameCasePri.Data(),DCAXYOliHistNameCasePri.Data(),500,-5,5,500,-5,5);
            TH2F *outputHistSec=new TH2F(DCAXYOliHistNameCaseSec.Data(),DCAXYOliHistNameCaseSec.Data(),500,-5,5,500,-5,5);
            TH2F *outputHistMat=new TH2F(DCAXYOliHistNameCaseMat.Data(),DCAXYOliHistNameCaseMat.Data(),500,-5,5,500,-5,5);
            TH2F *outputHistSecLam=new TH2F(DCAXYOliHistNameCaseSecLam.Data(),DCAXYOliHistNameCaseSecLam.Data(),500,-5,5,500,-5,5);
            TH2F *outputHistSecSig=new TH2F(DCAXYOliHistNameCaseSecSig.Data(),DCAXYOliHistNameCaseSecSig.Data(),500,-5,5,500,-5,5);
            for (int iBin=1;iBin<=500;++iBin) {
              outputHistPri->SetBinContent(iBin,outputHistPri->GetYaxis()->FindBin(0.),dcaXYProtonPri->GetBinContent(i+1,iBin));
              outputHistSec->SetBinContent(iBin,dcaXYProtonSec->GetYaxis()->FindBin(0.),dcaXYProtonSec->GetBinContent(i+1,iBin));
              outputHistMat->SetBinContent(iBin,outputHistMat->GetYaxis()->FindBin(0.),dcaXYProtonMat->GetBinContent(i+1,iBin));
              outputHistSecLam->SetBinContent(iBin,outputHistSecLam->GetYaxis()->FindBin(0.),dcaXYProtonSecLam->GetBinContent(i+1,iBin));
              outputHistSecSig->SetBinContent(iBin,outputHistSecSig->GetYaxis()->FindBin(0.),dcaXYProtonSecSig->GetBinContent(i+1,iBin));
            }
            SPdir->Add(outputHistPri);
            SPdir->Add(outputHistSec);
            SPdir->Add(outputHistMat);
            SPdir->Add(outputHistSecLam);
            SPdir->Add(outputHistSecSig);
          }
        } else {
          std::cout << "No DCAPtBinning \n";
        }
      }
    }

    TString MCv0Name="";
    MCv0Name+=prefix;
    MCv0Name+="v0CutsMC";


    TDirectoryFile *dirMCv0Cuts=(TDirectoryFile*)(file->FindObjectAny(MCv0Name.Data()));
    if (dirMCv0Cuts) {
      TList *tmp=dirMCv0Cuts->GetListOfKeys();
      TString name=tmp->At(0)->GetName();
      TList *mcv0Cuts;
      dirMCv0Cuts->GetObject(name,mcv0Cuts);
      if (!mcv0Cuts) {
        std::cout << "No MCv0Cuts \n";
      } else {
        TList *v0MC=(TList*)mcv0Cuts->FindObject("v0MonteCarlo");
        if (v0MC) {
          TList *CPAPtBinning=(TList*)v0MC->FindObject("CPAPtBinning");
          if (CPAPtBinning) {
            TH2F *CPALamPri=(TH2F*)CPAPtBinning->FindObject("CPAPtBinningPri");
            TH2F *CPALamMat=(TH2F*)CPAPtBinning->FindObject("CPAPtBinningMat");
            TH2F *CPALamSec=(TH2F*)CPAPtBinning->FindObject("CPAPtBinningSec");
            TH2F *CPALamCont=(TH2F*)CPAPtBinning->FindObject("CPAPtBinningCont");

            for (int iBins=0;iBins<CPALamPri->GetXaxis()->GetNbins();++iBins) {
              TString PriName=Form("fLambdaCPAPrimaryPtBin%i",iBins);
              TH1F *CPALamPriPt=(TH1F*)CPALamPri->ProjectionY(PriName.Data(),iBins+1,iBins+2);
              SPdir->Add(CPALamPriPt);
              TString SecName=Form("fLambdaCPASecondaryPtBin%i",iBins);
              TH1F *CPALamSecPt=(TH1F*)CPALamMat->ProjectionY(SecName.Data(),iBins+1,iBins+2);
              SPdir->Add(CPALamSecPt);
              TString MatName=Form("fLambdaCPAMaterialPtBin%i",iBins);
              TH1F *CPALamMatPt=(TH1F*)CPALamSec->ProjectionY(MatName.Data(),iBins+1,iBins+2);
              SPdir->Add(CPALamMatPt);
              TString BkgName=Form("fLambdaCPABkgPtBin%i",iBins);
              TH1F *CPALamBkgPt=(TH1F*)CPALamCont->ProjectionY(BkgName.Data(),iBins+1,iBins+2);
              SPdir->Add(CPALamBkgPt);
            }

          } else {
            std::cout << "No CPA PtBinning \n";
          }
        } else {
          std::cout << "No v0MC \n";
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

  TDirectoryFile* dirOutputMine=new TDirectoryFile("Results","Results");
  dirOutputMine->Add(outResults);
  dirOutputMine->Write("Results",TObject::kSingleKey);

  //for the systematics:
  if (!isMC) {
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
}
