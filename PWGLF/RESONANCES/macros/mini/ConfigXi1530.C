/***************************************************************************
 ** Bong-Hwi Lim
 ** Configuration script for Xi1530 analysis based on ConfigPhiPP13Tev_PID.C
 ****************************************************************************/

Bool_t ConfigXi1530
(
 AliRsnMiniAnalysisTask *task,
 Bool_t                 isMC,
 Bool_t                 isPP,
 const char             *suffix,
 AliRsnCutSet           *cutsPair,
 Int_t                  aodFilterBit=5,
 Int_t                  customQualityCutsID=AliRsnCutSetDaughterParticle::kDisableCustom,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate=AliRsnCutSetDaughterParticle::kQualityStd2011,
 Float_t                nsigmaPr=3.,
 Bool_t                 enableMonitor=kTRUE,
 Bool_t                 IsMcTrueOnly=kFALSE,
 TString                monitorOpt="",
 Bool_t                 useMixLS=0,
 Bool_t                 checkReflex=0,
 AliRsnMiniValue::EType yaxisVar=AliRsnMiniValue::kPt,
 TString                polarizationOpt="", /* J - Jackson,T - Transversity */
 UInt_t                 triggerMask=AliVEvent::kINT7
 )
{
    // manage suffix
    if(strlen(suffix)>0) suffix=Form("_%s",suffix);
    
    // set daughter cuts
    AliRsnCutSetDaughterParticle* cutSetPrimPion;
    AliRsnCutSetDaughterParticle* cutSetXiPion;
    AliRsnCutSetDaughterParticle* cutSetV0Pion;
    AliRsnCutSetDaughterParticle* cutSetV0Proton;
    
    Int_t MultBins=aodFilterBit/100;
    aodFilterBit=aodFilterBit%100;
    
    Bool_t misIDpion=false;
    if(customQualityCutsID==1000){
        customQualityCutsID=1;
        misIDpion=true;
    }
    
    Float_t nsigmaPrTPC=fmod(nsigmaPr,1000.);
    Float_t nsigmaPrTOF=(nsigmaPr-fmod(nsigmaPr,1000.))/1000.;
    if(nsigmaPrTOF<1.e-10) nsigmaPrTOF=-1.;
    
    AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
    if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,aodFilterBit)){
        //Set custom quality cuts for systematic checks
        cutSetPrimPion=new AliRsnCutSetDaughterParticle(Form("cutPrimPion_bit%i",aodFilterBit),
                                                        trkQualityCut,
                                                        AliRsnCutSetDaughterParticle::kQualityStd2011,
                                                        AliPID::kPion,
                                                        -1.);
        cutSetXiPion=new AliRsnCutSetDaughterParticle(Form("cutXiPion_bit%i",aodFilterBit),
                                                        trkQualityCut,
                                                        AliRsnCutSetDaughterParticle::kQualityStd2011,
                                                        AliPID::kPion,
                                                        -1.);
        cutSetV0Pion=new AliRsnCutSetDaughterParticle(Form("cutV0Pion_bit%i",aodFilterBit),
                                                        trkQualityCut,
                                                        AliRsnCutSetDaughterParticle::kQualityStd2011,
                                                        AliPID::kPion,
                                                        -1.);
        cutSetV0Proton=new AliRsnCutSetDaughterParticle(Form("cutV0Proton%i_%2.1fsigma",cutPrCandidate, nsigmaPr),
                                                        trkQualityCut,
                                                        cutPrCandidate,
                                                        AliPID::kProton,
                                                        nsigmaPrTPC,
                                                        nsigmaPrTOF);
    }else{
        //use default quality cuts std 2011 with crossed rows TPC
        Bool_t useCrossedRows = 1;
        cutSetPrimPion=new AliRsnCutSetDaughterParticle(Form("cutPrimPion_bit%i",aodFilterBit),
                                                        AliRsnCutSetDaughterParticle::kQualityStd2011,
                                                        AliPID::kPion,
                                                        -1.,
                                                        aodFilterBit,
                                                        useCrossedRows);
        cutSetXiPion=new AliRsnCutSetDaughterParticle(Form("cutXiPion_bit%i",aodFilterBit),
                                                        AliRsnCutSetDaughterParticle::kQualityStd2011,
                                                        AliPID::kPion,
                                                        -1.,
                                                        aodFilterBit,
                                                        useCrossedRows);
        cutSetV0Pion=new AliRsnCutSetDaughterParticle(Form("cutC0Pion_bit%i",aodFilterBit),
                                                        AliRsnCutSetDaughterParticle::kQualityStd2011,
                                                        AliPID::kPion,
                                                        -1.,
                                                        aodFilterBit,
                                                        useCrossedRows);
        cutSetV0Proton=new AliRsnCutSetDaughterParticle(Form("cutV0Proton%i_%2.1fsigma",cutPrCandidate,nsigmaPr),
                                                        cutPrCandidate,
                                                        AliPID::kPrton,
                                                        nsigmaPr,
                                                        aodFilterBit,
                                                        useCrossedRows);
    }
    
    Int_t iCutPrimPion=task->AddTrackCuts(cutSetPrimPion);
    Int_t iCutXiPion=task->AddTrackCuts(cutSetXiPion);
    Int_t iCutV0Pion=task->AddTrackCuts(cutSetV0Pion);
    Int_t iCutV0Prton=task->AddTrackCuts(cutSetV0Proton);
    
    if(enableMonitor){
        Printf("======== Cut monitoring enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC, cutSetPrimPion->GetMonitorOutput(), monitorOpt.Data());
        AddMonitorOutput(isMC, cutSetXiPion->GetMonitorOutput(), monitorOpt.Data());
        AddMonitorOutput(isMC, cutSetV0Pion->GetMonitorOutput(), monitorOpt.Data());
        AddMonitorOutput(isMC, cutSetV0Proton->GetMonitorOutput()), monitorOpt.Data();
    }
    
    // -- Values ------------------------------------------------------------------------------------
    /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,kFALSE);
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP,kFALSE);
    /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP,kFALSE);
    /* cos(theta) J     */ Int_t ctjID  = task->CreateValue(AliRsnMiniValue::kCosThetaJackson,kFALSE);
    /* cos(theta) J (MC)*/ Int_t ctjmID  = task->CreateValue(AliRsnMiniValue::kCosThetaJackson,kTRUE);
    /* cos(theta) T     */ Int_t cttID  = task->CreateValue(AliRsnMiniValue::kCosThetaTransversity,kFALSE);
    /* cos(theta) T (MC)*/ Int_t cttmID  = task->CreateValue(AliRsnMiniValue::kCosThetaTransversity,kTRUE);
    
    Double_t multbins[200];
    int j,nmult=0;
    if(isMC){
        for(j=0;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
        for(j=1;j<50;j++){multbins[nmult]=0.01*j; nmult++;}
        for(j=5;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
        for(j=1;j<10;j++){multbins[nmult]=j; nmult++;}
        for(j=2;j<=20;j++){multbins[nmult]=5.*j; nmult++;}
    }else if(triggerMask==AliVEvent::kHighMultV0){
        for(j=0;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
        for(j=1;j<50;j++){multbins[nmult]=0.01*j; nmult++;}
        for(j=5;j<=10;j++){multbins[nmult]=0.1*j; nmult++;}
    }else{
        for(j=0;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
        for(j=1;j<10;j++){multbins[nmult]=j; nmult++;}
        for(j=2;j<=20;j++){multbins[nmult]=5.*j; nmult++;}
    }
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    // [0] = unlike (+,-)
    // [1] = unlike (-,+)
    // [2] = mixing (+,-)
    // [3] = mixing (+,-)
    // [4] = likesign ++
    // [5] = likesign --
    
    Bool_t  use    [13]={!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,isMC,isMC,isMC,isMC,isMC,useMixLS,useMixLS};
    Int_t   useIM  [13]={      1      ,      1      ,      1      ,      1      ,      1      ,      1      ,  1 ,  1 ,  2 ,  2 ,  0 ,   1    ,   1    };
    TString name   [13]={  "Unlike"   ,  "Unlike"   ,   "Mixing"  ,   "Mixing"  ,   "LikePP"  ,   "LikeMM"  ,"Trues","TruesFine","TruesMM","TruesFineMM","Res","MixingPP","MixingMM"};
    TString comp   [13]={   "PAIR"    ,   "PAIR"    ,    "MIX"    ,    "MIX"    ,    "PAIR"   ,    "PAIR"   ,"TRUE","TRUE","TRUE","TRUE","TRUE","MIX","MIX"};
    TString output [13]={   "HIST"    ,   "HIST"    ,    "HIST"   ,    "HIST"   ,    "HIST"   ,    "HIST"   ,"HIST","HIST","HIST","HIST","HIST","HIST","HIST"};
    Int_t   pdgCode[13]={    3324     ,    3324     ,     3324    ,     3324    ,     3324    ,     3324    ,3324,3324,3324,3324,3324,3324,3324};
    Char_t  charge1[13]={    '+'      ,    '-'      ,     '+'     ,     '-'     ,     '+'     ,     '-'     ,'+','+','+','+','+','+','-'};
    Char_t  charge2[13]={    '-'      ,    '+'      ,     '-'     ,     '+'     ,     '+'     ,     '-'     ,'-','-','-','-','-','+','-'};
    
    for(Int_t i=0;i<12;i++){
        if(!use[i]) continue;
        AliRsnMiniOutput* out=task->CreateOutput(Form("Xi1530_%s%s",name[i].Data(),suffix),output[i].Data(),comp[i].Data());
        out->SetCutID(0,iCutPrimPion);
        out->SetCutID(1,iCutK); // NEED TO CHECK
        out->SetDaughter(0,AliRsnDaughter::kPion);
        out->SetDaughter(1,AliRsnDaughter::kXi);
        out->SetCharge(0,charge1[i]);
        out->SetCharge(1,charge2[i]);
        out->SetMotherPDG(pdgCode[i]);
        out->SetMotherMass(1.53178); // Xi1530 PDG Mass
        out->SetPairCuts(cutsPair);
        
        //axis X: invmass (or resolution)
        if(useIM[i]==1) out->AddAxis(imID,1000,1.5,2.5);
        else if(useIM[i]==2) out->AddAxis(mmID,1000,1.5,2.5);
        else out->AddAxis(diffID,200,-0.02,0.02);
        
        //axis Y: transverse momentum of pair as default - else chosen value
        if(yaxisVar==AliRsnMiniValue::kFirstDaughterPt) out->AddAxis(fdpt,150,0.,15.);
        else if(yaxisVar==AliRsnMiniValue::kSecondDaughterPt) out->AddAxis(sdpt,150,0.,15.);
        else if(yaxisVar==AliRsnMiniValue::kFirstDaughterP) out->AddAxis(fdp,150,0.,15.);
        else if(yaxisVar==AliRsnMiniValue::kSecondDaughterP)  out->AddAxis(sdp,150,0.,15.);
        else if(isMC && (i==5 || i==7)) out->AddAxis(ptID,300,0.,3.);//fine binning for efficiency weighting
        else out->AddAxis(ptID,200,0.,20.);//default use mother pt
        
        // axis Z: centrality-multiplicity
        if(!isPP || MultBins) out->AddAxis(centID,nmult,multbins);//out->AddAxis(centID,100,0.,100.);
        else out->AddAxis(centID,nmult,multbins);//out->AddAxis(centID,161,-0.5,160.5);
        // axis W: pseudorapidity
        // out->AddAxis(etaID, 20, -1.0, 1.0);
        // axis J: rapidity
        // out->AddAxis(yID, 10, -0.5, 0.5);
        
        if (polarizationOpt.Contains("J")) out->AddAxis(ctjID,21,-1.,1);
        if (polarizationOpt.Contains("T")) out->AddAxis(cttID,21,-1.,1);
    }
    
    if(isMC){
        //get mothers for phi PDG = 3324
        AliRsnMiniOutput* outm=task->CreateOutput(Form("Xi1530_Mother%s", suffix),"HIST","MOTHER");
        outm->SetDaughter(0,AliRsnDaughter::kPion);
        outm->SetDaughter(1,AliRsnDaughter::kXi);
        outm->SetMotherPDG(3324);
        outm->SetMotherMass(1.53178); // Xi1530 PDG Mass
        outm->SetPairCuts(cutsPair);
        outm->AddAxis(imID,500,1.5,2.0);
        outm->AddAxis(ptID,200,0.,20.);
        outm->AddAxis(centID,nmult,multbins);
        //if(!isPP || MultBins) outm->AddAxis(centID,100,0.,100.);
        //else outm->AddAxis(centID,161,-0.5,160.5);
        if (polarizationOpt.Contains("J")) outm->AddAxis(ctjmID,21,-1.,1.);
        if (polarizationOpt.Contains("T")) outm->AddAxis(cttmID,21,-1.,1.);
        
        //get phase space of the decay from mothers
        AliRsnMiniOutput* outps=task->CreateOutput(Form("Xi1530_phaseSpace%s", suffix),"HIST","TRUE");
        outps->SetDaughter(0,AliRsnDaughter::kPion);
        outps->SetDaughter(1,AliRsnDaughter::kXi);
        outps->SetCutID(0,iCutPrimPion);
        outps->SetCutID(1,iCutK); // NEED TO CHECK
        outps->SetMotherPDG(3324);
        outps->SetMotherMass(1.53178); // Xi1530 PDG Mass
        outps->SetPairCuts(cutsPair);
        outps->AddAxis(fdpt,100,0.,10.);
        outps->AddAxis(sdpt,100,0.,10.);
        outps->AddAxis(ptID,200,0.,20.);
        
        AliRsnMiniOutput* outpsf=task->CreateOutput(Form("Xi1530_phaseSpaceFine%s", suffix),"HIST","TRUE");
        outpsf->SetDaughter(0,AliRsnDaughter::kPion);
        outpsf->SetDaughter(1,AliRsnDaughter::kXi);
        outpsf->SetCutID(0,iCutPrimPion);
        outpsf->SetCutID(1,iCutK); // NEED TO CHECK
        outpsf->SetMotherPDG(3324);
        outpsf->SetMotherMass(1.53178); // Xi1530 PDG Mass
        outpsf->SetPairCuts(cutsPair);
        outpsf->AddAxis(fdpt,30,0.,3.);
        outpsf->AddAxis(sdpt,30,0.,3.);
        outpsf->AddAxis(ptID,300,0.,3.);
        
        //get reflections
        if(checkReflex){
            AliRsnMiniOutput* outreflex=task->CreateOutput(Form("Xi1530_reflex%s", suffix),"SPARSE","TRUE");
            outreflex->SetDaughter(0,AliRsnDaughter::kPion);
            outreflex->SetDaughter(1,AliRsnDaughter::kXi);
            outreflex->SetCutID(0,iCutPrimPion);
            outreflex->SetCutID(1,iCutK);
            outreflex->SetMotherPDG(3324);
            outreflex->SetMotherMass(1.53178); // Xi1530 PDG Mass
            outreflex->SetPairCuts(cutsPair);
            outreflex->AddAxis(imID,500,1.5,2.0);
            outreflex->AddAxis(ptID,200,0.,20.);
            if(!isPP) outreflex->AddAxis(centID,100,0.,100.);
            else outreflex->AddAxis(centID,400,0.5,400.5);
            if (polarizationOpt.Contains("J")) outreflex->AddAxis(ctjID,21,-1.,1.);
            if (polarizationOpt.Contains("T")) outreflex->AddAxis(cttID,21,-1.,1.);
        }//end reflections
    }//end MC
    
    //-------------------------------------------------------
    // misidentified daughters in simulation
    
    TString mName;
    Double_t mMass;
    Int_t mPDG0,mPDG1,mPDG2;
    
    if(isMC){
        for(Int_t i=0;i<10;i++){
            if(!i){
                mName.Form("pi0_ee");
                mMass=0.134977;
                mPDG0=111;
                mPDG1=mPDG2=11;
            }else if(i==1){
                mName.Form("Kx_pipi");
                mMass=0.493677;
                mPDG0=321;
                mPDG1=mPDG2=211;
            }else if(i==2){
                mName.Form("K0S_pipi");
                mMass=0.497614;
                mPDG0=310;
                mPDG1=mPDG2=211;
            }else if(i==3){
                mName.Form("K0L_pipi");
                mMass=0.497614;
                mPDG0=130;
                mPDG1=mPDG2=211;
            }else if(i==4){
                mName.Form("eta_pipi");
                mMass=0.547862;
                mPDG0=221;
                mPDG1=mPDG2=211;
            }else if(i==5){
                mName.Form("etaprime_pipi");
                mMass=0.95778;
                mPDG0=331;
                mPDG1=mPDG2=211;
            }else if(i==6){
                mName.Form("Lambdap_ppi");
                mMass=1.115683;
                mPDG0=3122;
                mPDG1=2212;
                mPDG2=211;
            }else if(i==7){
                mName.Form("Lambdaa_ppi");
                mMass=1.115683;
                mPDG0=-3122;
                mPDG1=211;
                mPDG2=2212;
            }else if(i==8){
                mName.Form("Lambda1520p_pK");
                mMass=1.5915;
                mPDG0=3124;
                mPDG1=2212;
                mPDG2=321;
            }else if(i==9){
                mName.Form("Lambda1520a_pK");
                mMass=1.5915;
                mPDG0=-3124;
                mPDG1=321;
                mPDG2=2212;
            }
            
            AliRsnMiniOutput* out=task->CreateOutput(Form("%s_true",mName.Data()),"HIST","TRUE");
            out->SetCutID(0,iCutPrimPion);
            out->SetCutID(1,iCutK); //NEED TO FIX
            out->SetDaughter(0,AliRsnDaughter::kPion);
            out->SetDaughter(1,AliRsnDaughter::kXi);
            out->SetCharge(0,'+');
            out->SetCharge(1,'-');
            out->SetMotherPDG(mPDG0);
            out->SetMotherMass(mMass);
            out->SetPairCuts(cutsPair);
            
            if(mPDG1==11) out->SetDaughterTrue(0,AliRsnDaughter::kElectron);
            else if(mPDG1==211) out->SetDaughterTrue(0,AliRsnDaughter::kPion);
            else if(mPDG1==321) out->SetDaughterTrue(0,AliRsnDaughter::kKaon);
            else if(mPDG1==2212) out->SetDaughterTrue(0,AliRsnDaughter::kProton);
            else if(mPDG1==3312) out->SetDaughterTrue(0,AliRsnDaughter::kXi);
            
            if(mPDG2==11) out->SetDaughterTrue(1,AliRsnDaughter::kElectron);
            else if(mPDG2==211) out->SetDaughterTrue(1,AliRsnDaughter::kPion);
            else if(mPDG2==321) out->SetDaughterTrue(1,AliRsnDaughter::kKaon);
            else if(mPDG2==2212) out->SetDaughterTrue(1,AliRsnDaughter::kProton);
            else if(mPDG2==3312) out->SetDaughterTrue(1,AliRsnDaughter::kXi);
            
            out->AddAxis(imID,215,0.985,1.2);
            out->AddAxis(ptID,100,0.,20.);
        }
    }
    
    return kTRUE;
}

//-------------------------------------------------------
Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID = 0, Int_t customFilterBit = 0)
{
    //Sets configuration for track quality object different from std quality cuts.
    //Returns kTRUE if track quality cut object is successfully defined,
    //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
    //object to be configured does not exist.
    
    if ((!trkQualityCut)){
        Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
        return kFALSE;
    }
    
    if(customQualityCutsID>=1 && customQualityCutsID<100 && customQualityCutsID!=2){
        trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
        Printf(Form("::::: SetCustomQualityCut:: using standard 2011 track quality cuts"));
        
        if(!customFilterBit){//ESD
            if(customQualityCutsID==3){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.0150+0.0500/pt^1.1");}
            else if(customQualityCutsID==4){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.006+0.0200/pt^1.1");}
            else if(customQualityCutsID==5){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(5.);}
            else if(customQualityCutsID==6){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);}
            else if(customQualityCutsID==7){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(5.);}
            else if(customQualityCutsID==8){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(2.3);}
            else if(customQualityCutsID==9){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(60);}
            else if(customQualityCutsID==10){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);}
            else if(customQualityCutsID==11){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
            else if(customQualityCutsID==12){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}
            else if(customQualityCutsID==13){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(49.);}
            else if(customQualityCutsID==14){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(4.);}
            else if(customQualityCutsID==15){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(49.);}
            else if(customQualityCutsID==16){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(25.);}
            else if(customQualityCutsID==17){trkQualityCut->GetESDtrackCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);}
            else if(customQualityCutsID==56){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.);}
            else if(customQualityCutsID==58){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3.);}
            else if(customQualityCutsID==60){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(80);}
            else if(customQualityCutsID==64){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(25.);}
        }else{//AOD
            trkQualityCut->SetCheckOnlyFilterBit(kFALSE);
            if(customQualityCutsID==4){trkQualityCut->SetDCARPtFormula("0.006+0.0200/pt^1.1");}
            else if(customQualityCutsID==6){trkQualityCut->SetDCAZmax(0.2);}
            else if(customQualityCutsID==8){trkQualityCut->SetTrackMaxChi2(2.3);}
            else if(customQualityCutsID==10){trkQualityCut->SetMinNCrossedRowsTPC(100,kTRUE);}
            else if(customQualityCutsID==12){trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.9,kTRUE);}
            else if(customQualityCutsID==56){trkQualityCut->SetDCAZmax(1.);}
            else if(customQualityCutsID==58){trkQualityCut->SetTrackMaxChi2(3.5);}
            else if(customQualityCutsID==60){trkQualityCut->SetMinNCrossedRowsTPC(80,kTRUE);}
        }
        
        trkQualityCut->Print();
        return kTRUE;
    }else if(customQualityCutsID==2 || (customQualityCutsID>=100 && customQualityCutsID<200)){
        trkQualityCut->SetDefaultsTPCOnly(kTRUE);
        Printf(Form("::::: SetCustomQualityCut:: using TPC-only track quality cuts"));
        
        if(customQualityCutsID==103){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXY(3.);}
        else if(customQualityCutsID==104){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXY(1.);}
        else if(customQualityCutsID==105){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(4.);}
        else if(customQualityCutsID==106){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.);}
        else if(customQualityCutsID==107){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(7.);}
        else if(customQualityCutsID==108){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(2.5);}
        else if(customQualityCutsID==109){trkQualityCut->GetESDtrackCuts()->SetMinNClustersTPC(30);}
        else if(customQualityCutsID==110){trkQualityCut->GetESDtrackCuts()->SetMinNClustersTPC(85);}
        
        trkQualityCut->Print();
        return kTRUE;
    }else{
        Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
        return kFALSE;
    }
    
    //for pA 2013
    //trkQualityCut->SetDefaults2011();//with filter bit=10
    //reset filter bit to very loose cuts
    trkQualityCut->SetAODTestFilterBit(customFilterBit);
    //apply all other cuts "by hand"
    trkQualityCut->SetCheckOnlyFilterBit(kFALSE);
    trkQualityCut->SetMinNCrossedRowsTPC(70, kTRUE);
    trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.8, kTRUE);
    trkQualityCut->SetMaxChi2TPCConstrainedGlobal(36);//used for ESD only - for AOD does not correspond to any cut
    trkQualityCut->SetTPCmaxChi2(4.0); //already in filter bit 0
    trkQualityCut->SetRejectKinkDaughters(kTRUE); //already in filter bit 0
    trkQualityCut->SetSPDminNClusters(AliESDtrackCuts::kAny);
    trkQualityCut->SetITSmaxChi2(36);
    trkQualityCut->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);//already in defaults 2011
    trkQualityCut->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);//already in defaults 2011
    trkQualityCut->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);//already in defaults 2011
    
    if (customQualityCutsID==AliRsnCutSetDaughterParticle::kFilterBitCustom) {
        trkQualityCut->SetCheckOnlyFilterBit(kTRUE);
    }
    
    if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdLooserDCAXY){
        trkQualityCut->SetDCARmax(2.4);
    } else {
        trkQualityCut->SetDCARPtFormula("0.0105+0.0350/pt^1.1");
    }
    
    if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdLooserDCAZ){
        trkQualityCut->SetDCAZmax(3.2);
    } else {
        trkQualityCut->SetDCAZmax(2.0);
    }
    
    if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdCrossedRows60){
        trkQualityCut->SetMinNCrossedRowsTPC(60, kTRUE);
    }
    
    if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdCrossedRows80){
        trkQualityCut->SetMinNCrossedRowsTPC(80, kTRUE);
    }
    
    if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdRowsToCls075){
        trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.75, kTRUE);
    }
    
    if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdRowsToCls085){
        trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.85, kTRUE);
    }
    
    if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdCls70){
        trkQualityCut->SetAODTestFilterBit(10);
        trkQualityCut->SetTPCminNClusters(70);
    }
    
    if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdChi2TPCCls35){
        trkQualityCut->SetTPCmaxChi2(3.5);
    }
    
    trkQualityCut->SetPtRange(0.15, 20.0);
    trkQualityCut->SetEtaRange(-0.8, 0.8);
    
    Printf(Form("::::: SetCustomQualityCut:: using custom track quality cuts #%i",customQualityCutsID));
    trkQualityCut->Print();
    return kTRUE;
}
