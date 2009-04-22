#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, TString * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

void analysisClass::Loop()
{
  std::cout << "analysisClass::Loop() begins" <<std::endl;   
    
  if (fChain == 0) return;
   
  //////////book histos here

  TH1F *h_ele_E = new TH1F ("ele_E","ele_E",400,0,4000); 
  h_ele_E->Sumw2();
  TH1F *h_ele_Pt  = new TH1F ("ele_Pt","ele_Pt",400,0,4000); 
  h_ele_Pt->Sumw2();
  TH1F *h_ele_Phi  = new TH1F ("ele_Phi","ele_Phi",71,-3.55,3.55); 
  h_ele_Phi->Sumw2();
  TH1F *h_ele_Eta  = new TH1F ("ele_Eta","ele_Eta",201,-10.05,10.05); 
  h_ele_Eta->Sumw2();
  TH1F *h_ele_CaloEnergy  = new TH1F ("ele_CaloEnergy","ele_CaloEnergy",400,0,4000); 
  h_ele_CaloEnergy->Sumw2();

  TH1F *h_N_ele_matched = new TH1F ("N_ele_matched","N_ele_matched",5,-0.5,4.5); 
  h_N_ele_matched->Sumw2();
  TH1F *h_ele_E_matched = new TH1F ("ele_E_matched","ele_E_matched",400,0,4000); 
  h_ele_E_matched->Sumw2();
  TH1F *h_ele_Pt_matched  = new TH1F ("ele_Pt_matched","ele_Pt_matched",400,0,4000); 
  h_ele_Pt_matched->Sumw2();
  TH1F *h_ele_Phi_matched  = new TH1F ("ele_Phi_matched","ele_Phi_matched",71,-3.55,3.55); 
  h_ele_Phi_matched->Sumw2();
  TH1F *h_ele_Eta_matched  = new TH1F ("ele_Eta_matched","ele_Eta_matched",201,-10.05,10.05); 
  h_ele_Eta_matched->Sumw2();
  TH1F *h_ele_CaloEnergy_matched  = new TH1F ("ele_CaloEnergy_matched","ele_CaloEnergy_matched",400,0,4000); 
  h_ele_CaloEnergy_matched->Sumw2();

  TH1F *h_ele_E_ID_ISO = new TH1F ("ele_E_ID_ISO","ele_E_ID_ISO",400,0,4000); 
  h_ele_E_ID_ISO->Sumw2();
  TH1F *h_ele_Pt_ID_ISO  = new TH1F ("ele_Pt_ID_ISO","ele_Pt_ID_ISO",400,0,4000); 
  h_ele_Pt_ID_ISO->Sumw2();
  TH1F *h_ele_Phi_ID_ISO  = new TH1F ("ele_Phi_ID_ISO","ele_Phi_ID_ISO",71,-3.55,3.55); 
  h_ele_Phi_ID_ISO->Sumw2();
  TH1F *h_ele_Eta_ID_ISO  = new TH1F ("ele_Eta_ID_ISO","ele_Eta_ID_ISO",201,-10.05,10.05); 
  h_ele_Eta_ID_ISO->Sumw2();
  TH1F *h_ele_CaloEnergy_ID_ISO  = new TH1F ("ele_CaloEnergy_ID_ISO","ele_CaloEnergy_ID_ISO",400,0,4000); 
  h_ele_CaloEnergy_ID_ISO->Sumw2();

  TH1F *h_N_ele_Gen = new TH1F ("N_ele_Gen","N_ele_Gen",5,-0.5,4.5); 
  h_N_ele_Gen->Sumw2();
  TH1F *h_ele_Pt_Gen  = new TH1F ("ele_Pt_Gen","ele_Pt_Gen",400,0,4000); 
  h_ele_Pt_Gen->Sumw2();
  TH1F *h_ele_Pt_Gen_etaCut  = new TH1F ("ele_Pt_Gen_etaCut","ele_Pt_Gen_etaCut",400,0,4000); 
  h_ele_Pt_Gen_etaCut->Sumw2();
  TH1F *h_ele_Pt_Gen_matched  = new TH1F ("ele_Pt_Gen_matched","ele_Pt_Gen_matched",400,0,4000); 
  h_ele_Pt_Gen_matched->Sumw2();
  TH1F *h_ele_Pt_Gen_matched_ID  = new TH1F ("ele_Pt_Gen_matched_ID","ele_Pt_Gen_matched_ID",400,0,4000); 
  h_ele_Pt_Gen_matched_ID->Sumw2();
  TH1F *h_ele_Pt_Gen_matched_ID_ISO = new TH1F ("ele_Pt_Gen_matched_ID_ISO","ele_Pt_Gen_matched_ID_ISO",400,0,4000); 
  h_ele_Pt_Gen_matched_ID_ISO->Sumw2();
  TH1F *h_ele_Eta_Gen  = new TH1F ("ele_Eta_Gen","ele_Eta_Gen",201,-10.05,10.05); 
  h_ele_Eta_Gen->Sumw2();
  TH1F *h_ele_Eta_Gen_lowPt  = new TH1F ("ele_Eta_Gen_lowPt","ele_Eta_Gen_lowPt",201,-10.05,10.05); 
  h_ele_Eta_Gen_lowPt->Sumw2();
  TH1F *h_ele_Eta_Gen_etaCut  = new TH1F ("ele_Eta_Gen_etaCut","ele_Eta_Gen_etaCut",201,-10.05,10.05); 
  h_ele_Eta_Gen_etaCut->Sumw2();
  TH1F *h_ele_Eta_Gen_matched  = new TH1F ("ele_Eta_Gen_matched","ele_Eta_Gen_matched",201,-10.05,10.05); 
  h_ele_Eta_Gen_matched->Sumw2();
  TH1F *h_ele_Eta_Gen_matched_lowPt  = new TH1F ("ele_Eta_Gen_matched_lowPt","ele_Eta_Gen_matched_lowPt",201,-10.05,10.05); 
  h_ele_Eta_Gen_matched_lowPt->Sumw2();
  TH1F *h_ele_Eta_Gen_matched_ID  = new TH1F ("ele_Eta_Gen_matched_ID","ele_Eta_Gen_matched_ID",201,-10.05,10.05); 
  h_ele_Eta_Gen_matched_ID->Sumw2();
  TH1F *h_ele_Eta_Gen_matched_ID_ISO  = new TH1F ("ele_Eta_Gen_matched_ID_ISO","ele_Eta_Gen_matched_ID_ISO",201,-10.05,10.05); 
  h_ele_Eta_Gen_matched_ID_ISO->Sumw2();
  TH2F *h_ele_Eta_Pt_Gen = new TH2F ("ele_Eta_Pt_Gen","ele_Eta_Pt_Gen",400,0,4000,201,-10.05,10.05);


  //extra stdies
  TH1F *h_ele_Pt_Gen_matched_ID_ISO_1 = new TH1F ("ele_Pt_Gen_matched_ID_ISO_1","ele_Pt_Gen_matched_ID_ISO_1",400,0,4000); 
  h_ele_Pt_Gen_matched_ID_ISO_1->Sumw2();
  TH1F *h_ele_Pt_Gen_matched_ID_ISO_2 = new TH1F ("ele_Pt_Gen_matched_ID_ISO_2","ele_Pt_Gen_matched_ID_ISO_2",400,0,4000); 
  h_ele_Pt_Gen_matched_ID_ISO_2->Sumw2();
  TH1F *h_ele_Pt_Gen_matched_ID_ISO_3 = new TH1F ("ele_Pt_Gen_matched_ID_ISO_3","ele_Pt_Gen_matched_ID_ISO_3",400,0,4000); 
  h_ele_Pt_Gen_matched_ID_ISO_3->Sumw2();
  TH1F *h_ele_Pt_Gen_matched_ID_ISO_4 = new TH1F ("ele_Pt_Gen_matched_ID_ISO_4","ele_Pt_Gen_matched_ID_ISO_4",400,0,4000); 
  h_ele_Pt_Gen_matched_ID_ISO_4->Sumw2();


  TH1F *h_DeltaR_Gen_Reco = new TH1F("DeltaR_Gen_Reco","DeltaR_Gen_Reco",500,0,10.); 
  h_DeltaR_Gen_Reco->Sumw2();
  TH1F *h_DeltaR_Gen_2ndReco = new TH1F("DeltaR_Gen_2ndReco","DeltaR_Gen_2ndReco",500,0,10.); 
  h_DeltaR_Gen_2ndReco->Sumw2();

  TH1F *h_N_ele_Pre = new TH1F ("N_ele_Pre","N_ele_Pre",11,-0.5,10.5); 
  h_N_ele_Pre->Sumw2();
  TH1F *h_N_ele_Post_ID = new TH1F ("N_ele_Post_ID","N_ele_Post_ID",11,-0.5,10.5); 
  h_N_ele_Post_ID->Sumw2();
  TH1F *h_N_ele_Post_ID_ISO = new TH1F ("N_ele_Post_ID_ISO","N_ele_Post_ID_ISO",11,-0.5,10.5); 
  h_N_ele_Post_ID_ISO->Sumw2();

  TH1F *h_Energy_Res = new TH1F ("Energy_Res","Energy_Res",150,0,1.5); 
  h_Energy_Res->Sumw2();
  TH1F *h_Energy_Res_barrel = new TH1F ("Energy_Res_barrel","Energy_Res_barrel",150,0,1.5); 
  h_Energy_Res_barrel->Sumw2();
  TH1F *h_Energy_Res_endcap = new TH1F ("Energy_Res_endcap","Energy_Res_endcap",150,0,1.5); 
  h_Energy_Res_endcap->Sumw2();

  //////matched
  TH1F *h_eleHoE_barrel_matched = new TH1F("eleHoE_barrel_matched","eleHoE_barrel_matched",100,0,0.1); 
  h_eleHoE_barrel_matched->Sumw2();
  TH1F *h_eleSigmaEE_barrel_matched = new TH1F("eleSigmaEE_barrel_matched","eleSigmaEE_barrel_matched",400,0,0.2); 
  h_eleSigmaEE_barrel_matched->Sumw2();
  TH1F *h_eleDeltaPhiTrkSC_barrel_matched = new TH1F("eleDeltaPhiTrkSC_barrel_matched","eleDeltaPhiTrkSC_barrel_matched",1000,-0.2,0.2); 
  h_eleDeltaPhiTrkSC_barrel_matched->Sumw2();
  TH1F *h_eleDeltaEtaTrkSC_barrel_matched = new TH1F("eleDeltaEtaTrkSC_barrel_matched","eleDeltaEtaTrkSC_barrel_matched",1000,-0.2,0.2); 
  h_eleDeltaEtaTrkSC_barrel_matched->Sumw2();
  TH1F *h_eleClassif_barrel_matched = new TH1F("eleClassif_barrel_matched","eleClassif_barrel_matched",201,-0.5,200.5); 
  h_eleClassif_barrel_matched->Sumw2();

  TH1F *h_eleNumTrkIso_barrel_matched = new TH1F("eleNumTrkIso_barrel_matched","eleNumTrkIso_barrel_matched",101,-0.5,100.5); 
  h_eleNumTrkIso_barrel_matched->Sumw2();
  TH1F *h_eleTrkIso_barrel_matched = new TH1F("eleTrkIso_barrel_matched","eleTrkIso_barrel_matched",1006,-5.5,1000.5); 
  h_eleTrkIso_barrel_matched->Sumw2();
  TH1F *h_eleEcalRecHitIso_barrel_matched = new TH1F("eleEcalRecHitIso_barrel_matched","eleEcalRecHitIso_barrel_matched",1006,-5.5,1000.5); 
  h_eleEcalRecHitIso_barrel_matched->Sumw2();
  TH1F *h_eleHcalRecHitIso_barrel_matched = new TH1F("eleHcalRecHitIso_barrel_matched","eleHcalRecHitIso_barrel_matched",1006,-5.5,1000.5); 
  h_eleHcalRecHitIso_barrel_matched->Sumw2();

  //extra
  TH2F *h_eleEcalRecHitIso_barrel_matched_vs_pt 
    = new TH2F("eleEcalRecHitIso_barrel_matched_vs_pt","eleEcalRecHitIso_barrel_matched_vs_pt",1006,-5.5,1000.5, 1000, 0, 1000);
  TH2F *h_eleHcalRecHitIso_barrel_matched_vs_pt 
    = new TH2F("eleHcalRecHitIso_barrel_matched_vs_pt","eleHcalRecHitIso_barrel_matched_vs_pt",1006,-5.5,1000.5, 1000, 0, 1000);
  

  TH1F *h_eleHoE_endcap_matched = new TH1F("eleHoE_endcap_matched","eleHoE_endcap_matched",100,0,0.1); 
  h_eleHoE_endcap_matched->Sumw2();
  TH1F *h_eleSigmaEE_endcap_matched = new TH1F("eleSigmaEE_endcap_matched","eleSigmaEE_endcap_matched",400,0,0.2); 
  h_eleSigmaEE_endcap_matched->Sumw2();
  TH1F *h_eleDeltaPhiTrkSC_endcap_matched = new TH1F("eleDeltaPhiTrkSC_endcap_matched","eleDeltaPhiTrkSC_endcap_matched",1000,-0.2,0.2); 
  h_eleDeltaPhiTrkSC_endcap_matched->Sumw2();
  TH1F *h_eleDeltaEtaTrkSC_endcap_matched = new TH1F("eleDeltaEtaTrkSC_endcap_matched","eleDeltaEtaTrkSC_endcap_matched",1000,-0.2,0.2); 
  h_eleDeltaEtaTrkSC_endcap_matched->Sumw2();
  TH1F *h_eleClassif_endcap_matched = new TH1F("eleClassif_endcap_matched","eleClassif_endcap_matched",201,-0.5,200.5); 
  h_eleClassif_endcap_matched->Sumw2();

  TH1F *h_eleNumTrkIso_endcap_matched = new TH1F("eleNumTrkIso_endcap_matched","eleNumTrkIso_endcap_matched",101,-0.5,100.5); 
  h_eleNumTrkIso_endcap_matched->Sumw2();
  TH1F *h_eleTrkIso_endcap_matched = new TH1F("eleTrkIso_endcap_matched","eleTrkIso_endcap_matched",1006,-5.5,1000.5); 
  h_eleTrkIso_endcap_matched->Sumw2();
  TH1F *h_eleEcalRecHitIso_endcap_matched = new TH1F("eleEcalRecHitIso_endcap_matched","eleEcalRecHitIso_endcap_matched",1006,-5.5,1000.5); 
  h_eleEcalRecHitIso_endcap_matched->Sumw2();
  TH1F *h_eleHcalRecHitIso_endcap_matched = new TH1F("eleHcalRecHitIso_endcap_matched","eleHcalRecHitIso_endcap_matched",1006,-5.5,1000.5); 
  h_eleHcalRecHitIso_endcap_matched->Sumw2();

  //extra
  TH2F *h_eleEcalRecHitIso_endcap_matched_vs_pt 
    = new TH2F("eleEcalRecHitIso_endcap_matched_vs_pt","eleEcalRecHitIso_endcap_matched_vs_pt",1006,-5.5,1000.5, 1000, 0, 1000);
  TH2F *h_eleHcalRecHitIso_endcap_matched_vs_pt 
    = new TH2F("eleHcalRecHitIso_endcap_matched_vs_pt","eleHcalRecHitIso_endcap_matched_vs_pt",1006,-5.5,1000.5, 1000, 0, 1000);


  ////////unmatched
  TH1F *h_eleHoE_barrel_unmatched = new TH1F("eleHoE_barrel_unmatched","eleHoE_barrel_unmatched",100,0,0.1); 
  h_eleHoE_barrel_unmatched->Sumw2();
  TH1F *h_eleSigmaEE_barrel_unmatched = new TH1F("eleSigmaEE_barrel_unmatched","eleSigmaEE_barrel_unmatched",400,0,0.2); 
  h_eleSigmaEE_barrel_unmatched->Sumw2();
  TH1F *h_eleDeltaPhiTrkSC_barrel_unmatched = new TH1F("eleDeltaPhiTrkSC_barrel_unmatched","eleDeltaPhiTrkSC_barrel_unmatched",1000,-0.2,0.2); 
  h_eleDeltaPhiTrkSC_barrel_unmatched->Sumw2();
  TH1F *h_eleDeltaEtaTrkSC_barrel_unmatched = new TH1F("eleDeltaEtaTrkSC_barrel_unmatched","eleDeltaEtaTrkSC_barrel_unmatched",1000,-0.2,0.2); 
  h_eleDeltaEtaTrkSC_barrel_unmatched->Sumw2();
  TH1F *h_eleClassif_barrel_unmatched = new TH1F("eleClassif_barrel_unmatched","eleClassif_barrel_unmatched",201,-0.5,200.5); 
  h_eleClassif_barrel_unmatched->Sumw2();

  TH1F *h_eleNumTrkIso_barrel_unmatched = new TH1F("eleNumTrkIso_barrel_unmatched","eleNumTrkIso_barrel_unmatched",101,-0.5,100.5); 
  h_eleNumTrkIso_barrel_unmatched->Sumw2();
  TH1F *h_eleTrkIso_barrel_unmatched = new TH1F("eleTrkIso_barrel_unmatched","eleTrkIso_barrel_unmatched",1006,-5.5,1000.5); 
  h_eleTrkIso_barrel_unmatched->Sumw2();
  TH1F *h_eleEcalRecHitIso_barrel_unmatched = new TH1F("eleEcalRecHitIso_barrel_unmatched","eleEcalRecHitIso_barrel_unmatched",1006,-5.5,1000.5); 
  h_eleEcalRecHitIso_barrel_unmatched->Sumw2();
  TH1F *h_eleHcalRecHitIso_barrel_unmatched = new TH1F("eleHcalRecHitIso_barrel_unmatched","eleHcalRecHitIso_barrel_unmatched",1006,-5.5,1000.5); 
  h_eleHcalRecHitIso_barrel_unmatched->Sumw2();

  TH1F *h_eleHoE_endcap_unmatched = new TH1F("eleHoE_endcap_unmatched","eleHoE_endcap_unmatched",100,0,0.1); 
  h_eleHoE_endcap_unmatched->Sumw2();
  TH1F *h_eleSigmaEE_endcap_unmatched = new TH1F("eleSigmaEE_endcap_unmatched","eleSigmaEE_endcap_unmatched",400,0,0.2); 
  h_eleSigmaEE_endcap_unmatched->Sumw2();
  TH1F *h_eleDeltaPhiTrkSC_endcap_unmatched = new TH1F("eleDeltaPhiTrkSC_endcap_unmatched","eleDeltaPhiTrkSC_endcap_unmatched",1000,-0.2,0.2); 
  h_eleDeltaPhiTrkSC_endcap_unmatched->Sumw2();
  TH1F *h_eleDeltaEtaTrkSC_endcap_unmatched = new TH1F("eleDeltaEtaTrkSC_endcap_unmatched","eleDeltaEtaTrkSC_endcap_unmatched",1000,-0.2,0.2); 
  h_eleDeltaEtaTrkSC_endcap_unmatched->Sumw2();
  TH1F *h_eleClassif_endcap_unmatched = new TH1F("eleClassif_endcap_unmatched","eleClassif_endcap_unmatched",201,-0.5,200.5); 
  h_eleClassif_endcap_unmatched->Sumw2();

  TH1F *h_eleNumTrkIso_endcap_unmatched = new TH1F("eleNumTrkIso_endcap_unmatched","eleNumTrkIso_endcap_unmatched",101,-0.5,100.5); 
  h_eleNumTrkIso_endcap_unmatched->Sumw2();
  TH1F *h_eleTrkIso_endcap_unmatched = new TH1F("eleTrkIso_endcap_unmatched","eleTrkIso_endcap_unmatched",1006,-5.5,1000.5); 
  h_eleTrkIso_endcap_unmatched->Sumw2();
  TH1F *h_eleEcalRecHitIso_endcap_unmatched = new TH1F("eleEcalRecHitIso_endcap_unmatched","eleEcalRecHitIso_endcap_unmatched",1006,-5.5,1000.5); 
  h_eleEcalRecHitIso_endcap_unmatched->Sumw2();
  TH1F *h_eleHcalRecHitIso_endcap_unmatched = new TH1F("eleHcalRecHitIso_endcap_unmatched","eleHcalRecHitIso_endcap_unmatched",1006,-5.5,1000.5); 
  h_eleHcalRecHitIso_endcap_unmatched->Sumw2();

  TH1F *h_eleEff_Pt = new TH1F ("eleEff_Pt","eleEff_Pt",400,0,4000); 
  h_eleEff_Pt->Sumw2();
  TH1F *h_eleEff_Eta = new TH1F ("eleEff_Eta","eleEff_Eta",201,-10.05,10.05); 
  h_eleEff_Eta->Sumw2();
  TH1F *h_eleEff_Pt_ID = new TH1F ("eleEff_Pt_ID","eleEff_Pt_ID",400,0,4000); 
  h_eleEff_Pt_ID->Sumw2();
  TH1F *h_eleEff_Eta_ID = new TH1F ("eleEff_Eta_ID","eleEff_Eta_ID",201,-10.05,10.05); 
  h_eleEff_Eta_ID->Sumw2();
  TH1F *h_eleEff_Pt_ID_ISO = new TH1F ("eleEff_Pt_ID_ISO","eleEff_Pt_ID_ISO",400,0,4000); 
  h_eleEff_Pt_ID_ISO->Sumw2();
  TH1F *h_eleEff_Eta_ID_ISO = new TH1F ("eleEff_Eta_ID_ISO","eleEff_Eta_ID_ISO",201,-10.05,10.05); 
  h_eleEff_Eta_ID_ISO->Sumw2();

  //extra
  TH1F *h_eleEff_Pt_ID_ISO_1 = new TH1F ("eleEff_Pt_ID_ISO_1","eleEff_Pt_ID_ISO_1",400,0,4000); 
  h_eleEff_Pt_ID_ISO_1->Sumw2();
  TH1F *h_eleEff_Pt_ID_ISO_2 = new TH1F ("eleEff_Pt_ID_ISO_2","eleEff_Pt_ID_ISO_2",400,0,4000); 
  h_eleEff_Pt_ID_ISO_2->Sumw2();
  TH1F *h_eleEff_Pt_ID_ISO_3 = new TH1F ("eleEff_Pt_ID_ISO_3","eleEff_Pt_ID_ISO_3",400,0,4000); 
  h_eleEff_Pt_ID_ISO_3->Sumw2();
  TH1F *h_eleEff_Pt_ID_ISO_4 = new TH1F ("eleEff_Pt_ID_ISO_4","eleEff_Pt_ID_ISO_4",400,0,4000); 
  h_eleEff_Pt_ID_ISO_4->Sumw2();

  //getPreCutValue1("");
   
  /////////initialize variables
  int electron_PID=int(getPreCutValue1("electronPID"));
  int MotherPID=int(getPreCutValue1("motherPID"));
  float ConeSizeMCmatch_cut=getPreCutValue1("coneSizeMCmatchCut");

  /////////////Barrel Cuts ////////////////
  float eleHoE_barrel_cut=getPreCutValue1("ID_HoE_bar"); //0.05 = HEEP
  float eleSigmaEE_barrel_cut=getPreCutValue1("ID_sigEtaEta_bar"); //0.011 = HEEP
  float eleDeltaPhiTrkSC_barrel_cut=getPreCutValue1("ID_deltaPhi_bar"); //0.09 = HEEP
  float eleDeltaEtaTrkSC_barrel_cut=getPreCutValue1("ID_deltaEta_bar"); //0.005=HEEP

  int eleNumTrkIso_barrel_cut=int(getPreCutValue1("ISO_NumTrack_bar")); //4=HEEP
  float eleTrkIso_barrel_cut=getPreCutValue1("ISO_TrackIso_bar"); //7.5=HEEP
  float eleEcalRecHitIso_barrel_cut=getPreCutValue1("ISO_EcalIso_bar"); // 6+0.01*Et = HEEP
  float eleEcalRecHitIso_barrel_cut2=getPreCutValue2("ISO_EcalIso_bar"); // 6+0.01*Et = HEEP
  float eleHcalRecHitIso_barrel_cut=getPreCutValue1("ISO_HcalIso_bar"); // 4+0.005*Et = HEEP
  float eleHcalRecHitIso_barrel_cut2=getPreCutValue2("ISO_HcalIso_bar"); // 4+0.005*Et = HEEP

  float eleClassification_barrel_cut=getPreCutValue1("eleClass_bar"); // <40 = HEEP
  float eleEta_barrel_cut=getPreCutValue1("eleEta_bar");  //>1.422 HEEP

  /////////////End Cap Cuts ////////////////
  float eleHoE_endcap_cut=getPreCutValue1("ID_HoE_end"); //0.1 = HEEP
  float eleSigmaEE_endcap_cut=getPreCutValue1("ID_sigEtaEta_end"); //0.0275 = HEEP
  float eleDeltaPhiTrkSC_endcap_cut=getPreCutValue1("ID_deltaPhi_end"); //0.09 = HEEP
  float eleDeltaEtaTrkSC_endcap_cut=getPreCutValue1("ID_deltaEta_end"); //0.007=HEEP

  int eleNumTrkIso_endcap_cut=int(getPreCutValue1("ISO_NumTrack_end"));//4=HEEP
  float eleTrkIso_endcap_cut=getPreCutValue1("ISO_TrackIso_end"); //15=HEEP
  float eleEcalRecHitIso_endcap_cut=getPreCutValue1("ISO_EcalIso_end"); // 6+0.01*Et= HEEP 
  float eleEcalRecHitIso_endcap_cut2=getPreCutValue2("ISO_EcalIso_end"); // 6+0.01*Et= HEEP 
  float eleHcalRecHitIso_endcap_cut=getPreCutValue1("ISO_HcalIso_end"); // 4+0.005*Et = HEEP
  float eleHcalRecHitIso_endcap_cut2=getPreCutValue2("ISO_HcalIso_end"); // 4+0.005*Et = HEEP

  float eleClassification_endcap_cut=getPreCutValue1("eleClass_end"); // >=100 = HEEP
  float eleEta_endcap_cut=getPreCutValue1("eleEta_end");  //>1.560 HEEP
  float eleEta_endcap_cut2=getPreCutValue2("eleEta_end"); //<2.5 HEEP

  //getPreCutValue1("");

  float elePt_cut=getPreCutValue1("elePtCut");
  float eleEta_cut=getPreCutValue1("eleEtaCut");
  float genPartPt_cut=getPreCutValue1("genPartPtCut");

  //-----

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
    // if (Cut(ientry) < 0) continue;

    ////////////////////// User's code starts here ///////////////////////

    ///Stuff to be done every event

    int nEleMatched = 0;
    int nGenEle = 0;
    int nEle_pre = 0;
    int nEle_post_ID = 0;
    int nEle_post_ID_ISO = 0;

    vector<int> v_idx_reco_matched;

    //////////Finding gen particles from "mother"
    for(int igen=0;igen<GenParticleCount;igen++)
      {
 
	//select gen particles from "mother" decay
	if(abs(GenParticlePdgId[igen])==electron_PID
	   && abs(GenParticlePdgId[GenParticleMotherIndex[igen]])==MotherPID)
	  {
	    nGenEle++;
	    h_ele_Pt_Gen->Fill(GenParticlePt[igen]);
	    h_ele_Eta_Gen->Fill(GenParticleEta[igen]);
	    h_ele_Eta_Pt_Gen->Fill(GenParticlePt[igen],GenParticleEta[igen]);
	    
	    //if (from_LQ_idx[0]==99) from_LQ_idx[0]=igen;
	    //else from_LQ_idx[1]=igen;

	    if (fabs(GenParticleEta[igen])<eleEta_cut){
	      h_ele_Pt_Gen_etaCut->Fill(GenParticlePt[igen]);
	      h_ele_Eta_Gen_etaCut->Fill(GenParticleEta[igen]);
	    }

	    if (GenParticlePt[igen]<genPartPt_cut)
	      h_ele_Eta_Gen_lowPt->Fill(GenParticleEta[igen]);

	    //skip events with gen particle with low pT or with high eta
	    if (GenParticlePt[igen]<genPartPt_cut || fabs(GenParticleEta[igen])>eleEta_cut)
	      continue;

	    TVector3 elegen;
	    elegen.SetPtEtaPhi(GenParticlePt[igen],
			       GenParticleEta[igen],
			       GenParticlePhi[igen]);	

	    ///////Calculating DeltaR and matching
	    float minDeltaR = 999;
	    float min2ndDeltaR = 999;
	    int idx_minDeltaR = -1;

	    for (int iele=0; iele<eleCount; iele++) {
	      if (elePt[iele]<elePt_cut) continue;
	      
	      TVector3 ele;
	      ele.SetPtEtaPhi(elePt[iele],
			      eleEta[iele],
			      elePhi[iele]);	
	      
	      float DeltaR_ele_elegen = elegen.DeltaR(ele);
	      
	      if (DeltaR_ele_elegen < minDeltaR) {
		min2ndDeltaR=minDeltaR;
		minDeltaR=DeltaR_ele_elegen;
		idx_minDeltaR = iele;
	      }
	      else if (DeltaR_ele_elegen < min2ndDeltaR){
		min2ndDeltaR=DeltaR_ele_elegen;
	      }

	    } // end loop over reco particles
	    
	    h_DeltaR_Gen_Reco->Fill(minDeltaR);
	    h_DeltaR_Gen_2ndReco->Fill(min2ndDeltaR);

	    //gen particle matched with reco candidate
	    if ( minDeltaR < ConeSizeMCmatch_cut )
	      {
		nEleMatched++;
		// the index of the matched reco particle is idx_minDeltaR;

		//to be used later
		v_idx_reco_matched.push_back(idx_minDeltaR);

		h_ele_Pt_Gen_matched->Fill(GenParticlePt[igen]);
		h_ele_Eta_Gen_matched->Fill(GenParticleEta[igen]);

		if (GenParticlePt[igen]<genPartPt_cut)
		  h_ele_Eta_Gen_matched_lowPt->Fill(GenParticleEta[igen]);
		h_ele_E_matched->Fill(eleEnergy[idx_minDeltaR]);
		h_ele_Pt_matched->Fill(elePt[idx_minDeltaR]);
		h_ele_Phi_matched->Fill(elePhi[idx_minDeltaR]);
		h_ele_Eta_matched->Fill(eleEta[idx_minDeltaR]);
		h_ele_CaloEnergy_matched->Fill(eleCaloEnergy[idx_minDeltaR]);
		h_Energy_Res->Fill(eleEnergy[idx_minDeltaR]/GenParticleE[igen]);


		//--

		bool InBarrel = false;
		bool InEndCap = false;
		bool PassID = false;
		bool PassISO = false;

		//check for barrel or endcap
		float electronEta = sqrt(eleEta[idx_minDeltaR]*eleEta[idx_minDeltaR]);
		if ( electronEta<getPreCutValue1("eleEta_bar") )
		  InBarrel = true;
		if ( (electronEta>getPreCutValue1("eleEta_end") ) && (electronEta<getPreCutValue2("eleEta_end")))
		  InEndCap = true;
		
		if (InBarrel){
		  h_eleHoE_barrel_matched->Fill(eleHoE[idx_minDeltaR]);
		  h_eleSigmaEE_barrel_matched->Fill(eleSigmaEE[idx_minDeltaR]);
		  h_eleDeltaPhiTrkSC_barrel_matched->Fill(eleDeltaPhiTrkSC[idx_minDeltaR]);
		  h_eleDeltaEtaTrkSC_barrel_matched->Fill(eleDeltaEtaTrkSC[idx_minDeltaR]);
		  h_eleClassif_barrel_matched->Fill(eleClassif[idx_minDeltaR]);
		  h_eleNumTrkIso_barrel_matched->Fill(eleNumTrkIso[idx_minDeltaR]);
		  h_eleTrkIso_barrel_matched->Fill(eleTrkIso[idx_minDeltaR]);
		  h_eleEcalRecHitIso_barrel_matched->Fill(eleEcalRecHitIso[idx_minDeltaR]);
		  h_eleHcalRecHitIso_barrel_matched->Fill(eleHcalRecHitIso[idx_minDeltaR]);
		  h_Energy_Res_barrel->Fill(eleEnergy[idx_minDeltaR]/GenParticleE[igen]);
		}
		
		if (InEndCap){
		  h_eleHoE_endcap_matched->Fill(eleHoE[idx_minDeltaR]);
		  h_eleSigmaEE_endcap_matched->Fill(eleSigmaEE[idx_minDeltaR]);
		  h_eleDeltaPhiTrkSC_endcap_matched->Fill(eleDeltaPhiTrkSC[idx_minDeltaR]);
		  h_eleDeltaEtaTrkSC_endcap_matched->Fill(eleDeltaEtaTrkSC[idx_minDeltaR]);
		  h_eleClassif_endcap_matched->Fill(eleClassif[idx_minDeltaR]);
		  h_eleNumTrkIso_endcap_matched->Fill(eleNumTrkIso[idx_minDeltaR]);
		  h_eleTrkIso_endcap_matched->Fill(eleTrkIso[idx_minDeltaR]);
		  h_eleEcalRecHitIso_endcap_matched->Fill(eleEcalRecHitIso[idx_minDeltaR]);
		  h_eleHcalRecHitIso_endcap_matched->Fill(eleHcalRecHitIso[idx_minDeltaR]);
		  h_Energy_Res_endcap->Fill(eleEnergy[idx_minDeltaR]/GenParticleE[igen]);
		}
		

		///check for ID 
		if (InBarrel)
		  if ((eleHoE[idx_minDeltaR]<eleHoE_barrel_cut)&&(eleSigmaEE[idx_minDeltaR]<eleSigmaEE_barrel_cut)
		      &&(eleDeltaPhiTrkSC[idx_minDeltaR]*eleDeltaPhiTrkSC[idx_minDeltaR])<(eleDeltaPhiTrkSC_barrel_cut*eleDeltaPhiTrkSC_barrel_cut)
		      &&(eleDeltaEtaTrkSC[idx_minDeltaR]*eleDeltaEtaTrkSC[idx_minDeltaR])<(eleDeltaEtaTrkSC_barrel_cut*eleDeltaEtaTrkSC_barrel_cut)
		      &&(eleClassif[idx_minDeltaR]<eleClassification_barrel_cut)
		      )
		    PassID=true;
		
		if (InEndCap)
		  if ((eleHoE[idx_minDeltaR]<eleHoE_endcap_cut)&&(eleSigmaEE[idx_minDeltaR]<eleSigmaEE_endcap_cut)
		      &&(eleDeltaPhiTrkSC[idx_minDeltaR]*eleDeltaPhiTrkSC[idx_minDeltaR])<(eleDeltaPhiTrkSC_endcap_cut*eleDeltaPhiTrkSC_endcap_cut)
		      &&(eleDeltaEtaTrkSC[idx_minDeltaR]*eleDeltaEtaTrkSC[idx_minDeltaR])<(eleDeltaEtaTrkSC_endcap_cut*eleDeltaEtaTrkSC_endcap_cut)
		      &&(eleClassif[idx_minDeltaR]>=eleClassification_endcap_cut)
		      )
		    PassID=true;
		
		if ( PassID ){
		  nEle_post_ID++;
		  h_ele_Pt_Gen_matched_ID->Fill(GenParticlePt[igen]);
		  h_ele_Eta_Gen_matched_ID->Fill(GenParticleEta[igen]);
		}
		
		///check for ISO
		if (InBarrel)
		  if ((eleNumTrkIso[idx_minDeltaR]<eleNumTrkIso_barrel_cut)&&(eleTrkIso[idx_minDeltaR]<eleTrkIso_barrel_cut)
		      &&(eleEcalRecHitIso[idx_minDeltaR]<(eleEcalRecHitIso_barrel_cut+(eleEcalRecHitIso_barrel_cut2*elePt[idx_minDeltaR])))
		      &&(eleHcalRecHitIso[idx_minDeltaR]<(eleHcalRecHitIso_barrel_cut+(eleHcalRecHitIso_barrel_cut2*elePt[idx_minDeltaR])))
		      )
		    PassISO=true;
		
		if (InEndCap)
		  if ((eleNumTrkIso[idx_minDeltaR]<eleNumTrkIso_endcap_cut)&&(eleTrkIso[idx_minDeltaR]<eleTrkIso_endcap_cut)
		      &&(eleEcalRecHitIso[idx_minDeltaR]<(eleEcalRecHitIso_endcap_cut+(eleEcalRecHitIso_endcap_cut2*elePt[idx_minDeltaR])))
		      &&(eleHcalRecHitIso[idx_minDeltaR]<(eleHcalRecHitIso_endcap_cut+(eleHcalRecHitIso_endcap_cut2*elePt[idx_minDeltaR])))
		      )
		    PassISO=true;


		if ( (PassID)&&(PassISO) ){
		  nEle_post_ID_ISO++;
		  h_ele_Pt_Gen_matched_ID_ISO->Fill(GenParticlePt[igen]);
		  h_ele_Eta_Gen_matched_ID_ISO->Fill(GenParticleEta[igen]);
		}



		//extra studies

		bool passISO_NumTrkIso=false;
		bool passISO_TrkIso=false;
		bool passISO_EcalIso=false;
		bool passISO_HcalIso=false;

		if (InBarrel){
		  if( eleNumTrkIso[idx_minDeltaR]<eleNumTrkIso_barrel_cut )
		    passISO_NumTrkIso=true;

		  if( eleTrkIso[idx_minDeltaR]<eleTrkIso_barrel_cut )
		    passISO_TrkIso=true;

		  if( eleEcalRecHitIso[idx_minDeltaR]<(eleEcalRecHitIso_barrel_cut+(eleEcalRecHitIso_barrel_cut2*elePt[idx_minDeltaR])) )
		    passISO_EcalIso=true;

		  if( eleHcalRecHitIso[idx_minDeltaR]<(eleHcalRecHitIso_endcap_cut+(eleHcalRecHitIso_endcap_cut2*elePt[idx_minDeltaR])) )
		    passISO_HcalIso=true;
		}

		if (InEndCap){
		  if( eleNumTrkIso[idx_minDeltaR]<eleNumTrkIso_endcap_cut )
		    passISO_NumTrkIso=true;

		  if( eleTrkIso[idx_minDeltaR]<eleTrkIso_endcap_cut )
		    passISO_TrkIso=true;

		  if( eleEcalRecHitIso[idx_minDeltaR]<(eleEcalRecHitIso_endcap_cut+(eleEcalRecHitIso_endcap_cut2*elePt[idx_minDeltaR])) )
		    passISO_EcalIso=true;

		  if( eleHcalRecHitIso[idx_minDeltaR]<(eleHcalRecHitIso_endcap_cut+(eleHcalRecHitIso_endcap_cut2*elePt[idx_minDeltaR])) )
		    passISO_HcalIso=true;
		}


		//1
		if(PassID && passISO_EcalIso)
		  h_ele_Pt_Gen_matched_ID_ISO_1->Fill(GenParticlePt[igen]);

		//2
		if(PassID && passISO_EcalIso && passISO_HcalIso)
		  h_ele_Pt_Gen_matched_ID_ISO_2->Fill(GenParticlePt[igen]);

		//3
		if(PassID && passISO_EcalIso && passISO_HcalIso && passISO_NumTrkIso)
		  h_ele_Pt_Gen_matched_ID_ISO_3->Fill(GenParticlePt[igen]);

		//4
		if(PassID && passISO_EcalIso && passISO_HcalIso && passISO_NumTrkIso && passISO_TrkIso )
		  h_ele_Pt_Gen_matched_ID_ISO_4->Fill(GenParticlePt[igen]);


		if (PassID && InBarrel)
		  {
		    //extra
		    h_eleEcalRecHitIso_barrel_matched_vs_pt->Fill(eleEcalRecHitIso[idx_minDeltaR],elePt[idx_minDeltaR]); 
		    h_eleHcalRecHitIso_barrel_matched_vs_pt->Fill(eleHcalRecHitIso[idx_minDeltaR],elePt[idx_minDeltaR]); 
		  }
		
		if (PassID && InEndCap)
		  {
		    //extra
		    h_eleEcalRecHitIso_endcap_matched_vs_pt->Fill(eleEcalRecHitIso[idx_minDeltaR],elePt[idx_minDeltaR]); 
		    h_eleHcalRecHitIso_endcap_matched_vs_pt->Fill(eleHcalRecHitIso[idx_minDeltaR],elePt[idx_minDeltaR]); 
		  }

		
	      }// end matched reco candidate


	  }// end selecting electrons from mother decay
      }// end loop over gen particles


    //loop over reco particles
    for (int iele=0; iele<eleCount; iele++) {

      if (elePt[iele]<elePt_cut) continue;

      nEle_pre++;

      bool EleIsMatched = false;
      bool InBarrel = false;
      bool InEndCap = false;
      bool PassID = false;
      bool PassISO = false;
      
      //check for barrel or endcap
      float electronEta = fabs(eleEta[iele]);
      if ( electronEta<getPreCutValue1("eleEta_bar") )
	InBarrel = true;
      if ( (electronEta>getPreCutValue1("eleEta_end")) && (electronEta<getPreCutValue2("eleEta_end")) )
	InEndCap = true;
      
      ///check for ID 
      if (InBarrel)
	if ((eleHoE[iele]<eleHoE_barrel_cut)&&(eleSigmaEE[iele]<eleSigmaEE_barrel_cut)
	    &&(eleDeltaPhiTrkSC[iele]*eleDeltaPhiTrkSC[iele])<(eleDeltaPhiTrkSC_barrel_cut*eleDeltaPhiTrkSC_barrel_cut)
	    &&(eleDeltaEtaTrkSC[iele]*eleDeltaEtaTrkSC[iele])<(eleDeltaEtaTrkSC_barrel_cut*eleDeltaEtaTrkSC_barrel_cut)
	    &&(eleClassif[iele]<eleClassification_barrel_cut)
	    )
	  PassID=true;
      
      if (InEndCap)
	if ((eleHoE[iele]<eleHoE_endcap_cut)&&(eleSigmaEE[iele]<eleSigmaEE_endcap_cut)
	    &&(eleDeltaPhiTrkSC[iele]*eleDeltaPhiTrkSC[iele])<(eleDeltaPhiTrkSC_endcap_cut*eleDeltaPhiTrkSC_endcap_cut)
	    &&(eleDeltaEtaTrkSC[iele]*eleDeltaEtaTrkSC[iele])<(eleDeltaEtaTrkSC_endcap_cut*eleDeltaEtaTrkSC_endcap_cut)
	    &&(eleClassif[iele]>=eleClassification_endcap_cut)
	    )
	  PassID=true;
      
      ///check for ISO
      if (InBarrel)
	if ((eleNumTrkIso[iele]<eleNumTrkIso_barrel_cut)&&(eleTrkIso[iele]<eleTrkIso_barrel_cut)
	    &&(eleEcalRecHitIso[iele]<(eleEcalRecHitIso_barrel_cut+(eleEcalRecHitIso_barrel_cut2*elePt[iele])))
	    &&(eleHcalRecHitIso[iele]<(eleHcalRecHitIso_barrel_cut+(eleHcalRecHitIso_barrel_cut2*elePt[iele])))
	    )
	  PassISO=true;
      
      if (InEndCap)
	if ((eleNumTrkIso[iele]<eleNumTrkIso_endcap_cut)&&(eleTrkIso[iele]<eleTrkIso_endcap_cut)
	    &&(eleEcalRecHitIso[iele]<(eleEcalRecHitIso_endcap_cut+(eleEcalRecHitIso_endcap_cut2*elePt[iele])))
	    &&(eleHcalRecHitIso[iele]<(eleHcalRecHitIso_endcap_cut+(eleHcalRecHitIso_endcap_cut2*elePt[iele])))
	    )
	  PassISO=true;

      h_ele_E->Fill(eleEnergy[iele]);
      h_ele_Pt->Fill(elePt[iele]);
      h_ele_Phi->Fill(elePhi[iele]);
      h_ele_Eta->Fill(eleEta[iele]);
      h_ele_CaloEnergy->Fill(eleCaloEnergy[iele]);
      
      if ( (PassID) && (PassISO) ){
	h_ele_E_ID_ISO->Fill(eleEnergy[iele]);
	h_ele_Pt_ID_ISO->Fill(elePt[iele]);
	h_ele_Phi_ID_ISO->Fill(elePhi[iele]);
	h_ele_Eta_ID_ISO->Fill(eleEta[iele]);
	h_ele_CaloEnergy_ID_ISO->Fill(eleCaloEnergy[iele]);
      }


      //see if the reco particle is unmatched
      for (int m=0; m<v_idx_reco_matched.size(); m++)
	{
	  if( iele == v_idx_reco_matched[m] )
	    {
	      EleIsMatched=true;
	      break;
	    }
	}

      //NOT matched reco candidate
      if(EleIsMatched == false)
	{
	  //cout << "Found unmatched ele" << endl;
	  if (InBarrel){
	    h_eleHoE_barrel_unmatched->Fill(eleHoE[iele]);
	    h_eleSigmaEE_barrel_unmatched->Fill(eleSigmaEE[iele]);
	    h_eleDeltaPhiTrkSC_barrel_unmatched->Fill(eleDeltaPhiTrkSC[iele]);
	    h_eleDeltaEtaTrkSC_barrel_unmatched->Fill(eleDeltaEtaTrkSC[iele]);
	    h_eleClassif_barrel_unmatched->Fill(eleClassif[iele]);
	    h_eleNumTrkIso_barrel_unmatched->Fill(eleNumTrkIso[iele]);
	    h_eleTrkIso_barrel_unmatched->Fill(eleTrkIso[iele]);
	    h_eleEcalRecHitIso_barrel_unmatched->Fill(eleEcalRecHitIso[iele]);
	    h_eleHcalRecHitIso_barrel_unmatched->Fill(eleHcalRecHitIso[iele]);
	  }
	  
	  if (InEndCap){
	    h_eleHoE_endcap_unmatched->Fill(eleHoE[iele]);
	    h_eleSigmaEE_endcap_unmatched->Fill(eleSigmaEE[iele]);
	    h_eleDeltaPhiTrkSC_endcap_unmatched->Fill(eleDeltaPhiTrkSC[iele]);
	    h_eleDeltaEtaTrkSC_endcap_unmatched->Fill(eleDeltaEtaTrkSC[iele]);
	    h_eleClassif_endcap_unmatched->Fill(eleClassif[iele]);
	    h_eleNumTrkIso_endcap_unmatched->Fill(eleNumTrkIso[iele]);
	    h_eleTrkIso_endcap_unmatched->Fill(eleTrkIso[iele]);
	    h_eleEcalRecHitIso_endcap_unmatched->Fill(eleEcalRecHitIso[iele]);
	    h_eleHcalRecHitIso_endcap_unmatched->Fill(eleHcalRecHitIso[iele]);
	  }
	} // end else unmatched
      
    }// end ele loop
 
    h_N_ele_Pre->Fill(nEle_pre);
    h_N_ele_Post_ID->Fill(nEle_post_ID);
    h_N_ele_Post_ID_ISO->Fill(nEle_post_ID_ISO);
 
    h_N_ele_matched->Fill(nEleMatched);
    h_N_ele_Gen->Fill(nGenEle);

    ////////////////////// User's code ends here ///////////////////////

  } // End loop over events

  //////////write histos 
  h_ele_E->Write();
  h_ele_Pt->Write();
  h_ele_Phi->Write();
  h_ele_Eta->Write();
  h_ele_CaloEnergy->Write();

  h_N_ele_matched->Write();
  h_ele_E_matched->Write();
  h_ele_Pt_matched->Write();
  h_ele_Phi_matched->Write();
  h_ele_Eta_matched->Write();
  h_ele_CaloEnergy_matched->Write();

  h_ele_E_ID_ISO->Write();
  h_ele_Pt_ID_ISO->Write();
  h_ele_Phi_ID_ISO->Write();
  h_ele_Eta_ID_ISO->Write();
  h_ele_CaloEnergy_ID_ISO->Write();

  h_N_ele_Gen->Write();
  h_ele_Pt_Gen->Write();
  h_ele_Pt_Gen_etaCut->Write();
  h_ele_Pt_Gen_matched->Write();
  h_ele_Pt_Gen_matched_ID->Write();
  h_ele_Pt_Gen_matched_ID_ISO->Write();
  h_ele_Eta_Gen->Write();
  h_ele_Eta_Gen_lowPt->Write();
  h_ele_Eta_Gen_etaCut->Write();
  h_ele_Eta_Gen_matched->Write();
  h_ele_Eta_Gen_matched_lowPt->Write();
  h_ele_Eta_Gen_matched_ID->Write();
  h_ele_Eta_Gen_matched_ID_ISO->Write();
  h_ele_Eta_Pt_Gen->Write();

  h_DeltaR_Gen_Reco->Write();
  h_DeltaR_Gen_2ndReco->Write();

  h_N_ele_Pre->Write();
  h_N_ele_Post_ID->Write();
  h_N_ele_Post_ID_ISO->Write();

  h_Energy_Res->Write();
  h_Energy_Res_barrel->Write();
  h_Energy_Res_endcap->Write();

  ////matched
  h_eleHoE_barrel_matched->Write();
  h_eleSigmaEE_barrel_matched->Write();
  h_eleDeltaPhiTrkSC_barrel_matched->Write();
  h_eleDeltaEtaTrkSC_barrel_matched->Write();
  h_eleClassif_barrel_matched->Write();
  h_eleNumTrkIso_barrel_matched->Write();
  h_eleTrkIso_barrel_matched->Write();
  h_eleEcalRecHitIso_barrel_matched->Write();
  h_eleHcalRecHitIso_barrel_matched->Write();

  h_eleHoE_endcap_matched->Write();
  h_eleSigmaEE_endcap_matched->Write();
  h_eleDeltaPhiTrkSC_endcap_matched->Write();
  h_eleDeltaEtaTrkSC_endcap_matched->Write();
  h_eleClassif_endcap_matched->Write();
  h_eleNumTrkIso_endcap_matched->Write();
  h_eleTrkIso_endcap_matched->Write();
  h_eleEcalRecHitIso_endcap_matched->Write();
  h_eleHcalRecHitIso_endcap_matched->Write();

  ///unmathced
  h_eleHoE_barrel_unmatched->Write();
  h_eleSigmaEE_barrel_unmatched->Write();
  h_eleDeltaPhiTrkSC_barrel_unmatched->Write();
  h_eleDeltaEtaTrkSC_barrel_unmatched->Write();
  h_eleClassif_barrel_unmatched->Write();
  h_eleNumTrkIso_barrel_unmatched->Write();
  h_eleTrkIso_barrel_unmatched->Write();
  h_eleEcalRecHitIso_barrel_unmatched->Write();
  h_eleHcalRecHitIso_barrel_unmatched->Write();

  h_eleHoE_endcap_unmatched->Write();
  h_eleSigmaEE_endcap_unmatched->Write();
  h_eleDeltaPhiTrkSC_endcap_unmatched->Write();
  h_eleDeltaEtaTrkSC_endcap_unmatched->Write();
  h_eleClassif_endcap_unmatched->Write();
  h_eleNumTrkIso_endcap_unmatched->Write();
  h_eleTrkIso_endcap_unmatched->Write();
  h_eleEcalRecHitIso_endcap_unmatched->Write();
  h_eleHcalRecHitIso_endcap_unmatched->Write();

  h_eleEff_Pt->Divide(h_ele_Pt_Gen_matched,h_ele_Pt_Gen,1,1);
  h_eleEff_Pt->Write();
  h_eleEff_Eta->Divide(h_ele_Eta_Gen_matched,h_ele_Eta_Gen,1,1);
  h_eleEff_Eta->Write();

  h_eleEff_Pt_ID->Divide(h_ele_Pt_Gen_matched_ID,h_ele_Pt_Gen,1,1);
  h_eleEff_Pt_ID->Write();
  h_eleEff_Eta_ID->Divide(h_ele_Eta_Gen_matched_ID,h_ele_Eta_Gen,1,1);
  h_eleEff_Eta_ID->Write();

  h_eleEff_Pt_ID_ISO->Divide(h_ele_Pt_Gen_matched_ID_ISO,h_ele_Pt_Gen,1,1);
  h_eleEff_Pt_ID_ISO->Write();
  h_eleEff_Eta_ID_ISO->Divide(h_ele_Eta_Gen_matched_ID_ISO,h_ele_Eta_Gen,1,1);
  h_eleEff_Eta_ID_ISO->Write();

  //extra
  h_ele_Pt_Gen_matched_ID_ISO_1->Write(); 
  h_ele_Pt_Gen_matched_ID_ISO_2->Write(); 
  h_ele_Pt_Gen_matched_ID_ISO_3->Write(); 
  h_ele_Pt_Gen_matched_ID_ISO_4->Write(); 

  h_eleEff_Pt_ID_ISO_1->Divide(h_ele_Pt_Gen_matched_ID_ISO_1,h_ele_Pt_Gen,1,1);
  h_eleEff_Pt_ID_ISO_1->Write();
  h_eleEff_Pt_ID_ISO_2->Divide(h_ele_Pt_Gen_matched_ID_ISO_2,h_ele_Pt_Gen,1,1);
  h_eleEff_Pt_ID_ISO_2->Write();
  h_eleEff_Pt_ID_ISO_3->Divide(h_ele_Pt_Gen_matched_ID_ISO_3,h_ele_Pt_Gen,1,1);
  h_eleEff_Pt_ID_ISO_3->Write();
  h_eleEff_Pt_ID_ISO_4->Divide(h_ele_Pt_Gen_matched_ID_ISO_4,h_ele_Pt_Gen,1,1);
  h_eleEff_Pt_ID_ISO_4->Write();
    
  h_eleEcalRecHitIso_barrel_matched_vs_pt->Write();
  h_eleHcalRecHitIso_barrel_matched_vs_pt->Write();
  h_eleEcalRecHitIso_endcap_matched_vs_pt->Write();
  h_eleHcalRecHitIso_endcap_matched_vs_pt->Write();
  
  std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
