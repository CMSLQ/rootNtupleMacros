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

  TH1F *h_muon_E = new TH1F ("muon_E","muon_E",100,0,1000);
  TH1F *h_muon_Pt  = new TH1F ("muon_Pt","muon_Pt",100,0,1000);
  TH1F *h_muon_Phi  = new TH1F ("muon_Phi","muon_Phi",71,-3.55,3.55);
  TH1F *h_muon_Eta  = new TH1F ("muon_Eta","muon_Eta",121,-6.05,6.05);

  TH1F *h_N_muon_matched = new TH1F ("N_muon_matched","N_muon_matched",5,-0.5,4.5);
  TH1F *h_muon_E_matched = new TH1F ("muon_E_matched","muon_E_matched",100,0,1000);
  TH1F *h_muon_Pt_matched  = new TH1F ("muon_Pt_matched","muon_Pt_matched",100,0,1000);
  TH1F *h_muon_Phi_matched  = new TH1F ("muon_Phi_matched","muon_Phi_matched",71,-3.55,3.55);
  TH1F *h_muon_Eta_matched  = new TH1F ("muon_Eta_matched","muon_Eta_matched",121,-6.05,6.05);

  TH1F *h_muon_E_ID_ISO = new TH1F ("muon_E_ID_ISO","muon_E_ID_ISO",100,0,1000);
  TH1F *h_muon_Pt_ID_ISO  = new TH1F ("muon_Pt_ID_ISO","muon_Pt_ID_ISO",100,0,1000);
  TH1F *h_muon_Phi_ID_ISO  = new TH1F ("muon_Phi_ID_ISO","muon_Phi_ID_ISO",71,-3.55,3.55);
  TH1F *h_muon_Eta_ID_ISO  = new TH1F ("muon_Eta_ID_ISO","muon_Eta_ID_ISO",121,-6.05,6.05);

  TH1F *h_N_muon_Gen = new TH1F ("N_muon_Gen","N_muon_Gen",5,-0.5,4.5);
  TH1F *h_muon_Pt_Gen  = new TH1F ("muon_Pt_Gen","muon_Pt_Gen",100,0,1000);
  TH1F *h_muon_Pt_Gen_etaCut  = new TH1F ("muon_Pt_Gen_etaCut","muon_Pt_Gen_etaCut",100,0,1000);
  TH1F *h_muon_Pt_Gen_matched  = new TH1F ("muon_Pt_Gen_matched","muon_Pt_Gen_matched",100,0,1000);
  TH1F *h_muon_Pt_Gen_matched_ID  = new TH1F ("muon_Pt_Gen_matched_ID","muon_Pt_Gen_matched_ID",100,0,1000);
  TH1F *h_muon_Pt_Gen_matched_ID_ISO = new TH1F ("muon_Pt_Gen_matched_ID_ISO","muon_Pt_Gen_matched_ID_ISO",100,0,1000);
  TH1F *h_muon_Eta_Gen  = new TH1F ("muon_Eta_Gen","muon_Eta_Gen",121,-6.05,6.05);
  TH1F *h_muon_Eta_Gen_lowPt  = new TH1F ("muon_Eta_Gen_lowPt","muon_Eta_Gen_lowPt",121,-6.05,6.05);
  TH1F *h_muon_Eta_Gen_etaCut  = new TH1F ("muon_Eta_Gen_etaCut","muon_Eta_Gen_etaCut",121,-6.05,6.05);
  TH1F *h_muon_Eta_Gen_matched  = new TH1F ("muon_Eta_Gen_matched","muon_Eta_Gen_matched",121,-6.05,6.05);
  TH1F *h_muon_Eta_Gen_matched_lowPt  = new TH1F ("muon_Eta_Gen_matched_lowPt","muon_Eta_Gen_matched_lowPt",121,-6.05,6.05);
  TH1F *h_muon_Eta_Gen_matched_ID  = new TH1F ("muon_Eta_Gen_matched_ID","muon_Eta_Gen_matched_ID",121,-6.05,6.05);
  TH1F *h_muon_Eta_Gen_matched_ID_ISO  = new TH1F ("muon_Eta_Gen_matched_ID_ISO","muon_Eta_Gen_matched_ID_ISO",121,-6.05,6.05);
  TH2F *h_muon_Eta_Pt_Gen = new TH2F ("muon_Eta_Pt_Gen","muon_Eta_Pt_Gen",100,0,1000,121,-6.05,6.05);

  TH1F *h_DeltaR_Gen_Reco = new TH1F("DeltaR_Gen_Reco","DeltaR_Gen_Reco",100,0,0.2);
  TH1F *h_DeltaR_Gen_2ndReco = new TH1F("DeltaR_Gen_2ndReco","DeltaR_Gen_2ndReco",200,0,1.0);

  TH1F *h_N_muon_Pre = new TH1F ("N_muon_Pre","N_muon_Pre",11,-0.5,10.5);
  TH1F *h_N_muon_Post_ID = new TH1F ("N_muon_Post_ID","N_muon_Post_ID",11,-0.5,10.5);
  TH1F *h_N_muon_Post_ID_ISO = new TH1F ("N_muon_Post_ID_ISO","N_muon_Post_ID_ISO",11,-0.5,10.5);

  TH1F *h_Energy_Res = new TH1F ("Energy_Res","Energy_Res",150,0,1.5);
  TH1F *h_Energy_Res_barrel = new TH1F ("Energy_Res_barrel","Energy_Res_barrel",150,0,1.5);
  TH1F *h_Energy_Res_endcap = new TH1F ("Energy_Res_endcap","Energy_Res_endcap",150,0,1.5);

  //////matched
  TH1F *h_muonGlobalChi2_barrel_matched = new TH1F("muonGlobalChi2_barrel_matched","muonGlobalChi2_barrel_matched",100,0,50);
  TH1F *h_muonTrkD0_barrel_matched = new TH1F("h_muonTrkD0_barrel_matched","h_muonTrkD0_barrel_matched",100,0,2);
  TH1F *h_muonTrkHits_barrel_matched = new TH1F("h_muonTrkHits_barrel_matched","h_muonTrkHits_barrel_matched",101,-0.5,100.5);
  TH1F *h_muonTrkIso_barrel_matched = new TH1F("muonTrkIso_barrel_matched","muonTrkIso_barrel_matched",200,-0.02,0.02);

  TH1F *h_muonGlobalChi2_endcap_matched = new TH1F("muonGlobalChi2_endcap_matched","muonGlobalChi2_endcap_matched",100,0,50);
  TH1F *h_muonTrkD0_endcap_matched = new TH1F("h_muonTrkD0_endcap_matched","h_muonTrkD0_endcap_matched",100,0,2);
  TH1F *h_muonTrkHits_endcap_matched = new TH1F("h_muonTrkHits_endcap_matched","h_muonTrkHits_endcap_matched",101,-0.5,100.5);
  TH1F *h_muonTrkIso_endcap_matched = new TH1F("muonTrkIso_endcap_matched","muonTrkIso_endcap_matched",200,-0.02,0.02);

  //////unmatched
  TH1F *h_muonGlobalChi2_barrel_unmatched = new TH1F("muonGlobalChi2_barrel_unmatched","muonGlobalChi2_barrel_unmatched",100,0,50);
  TH1F *h_muonTrkD0_barrel_unmatched = new TH1F("h_muonTrkD0_barrel_unmatched","h_muonTrkD0_barrel_unmatched",100,0,2);
  TH1F *h_muonTrkHits_barrel_unmatched = new TH1F("h_muonTrkHits_barrel_unmatched","h_muonTrkHits_barrel_unmatched",101,-0.5,100.5);
  TH1F *h_muonTrkIso_barrel_unmatched = new TH1F("muonTrkIso_barrel_unmatched","muonTrkIso_barrel_unmatched",200,-0.02,0.02);

  TH1F *h_muonGlobalChi2_endcap_unmatched = new TH1F("muonGlobalChi2_endcap_unmatched","muonGlobalChi2_endcap_unmatched",100,0,50);
  TH1F *h_muonTrkD0_endcap_unmatched = new TH1F("h_muonTrkD0_endcap_unmatched","h_muonTrkD0_endcap_unmatched",100,0,2);
  TH1F *h_muonTrkHits_endcap_unmatched = new TH1F("h_muonTrkHits_endcap_unmatched","h_muonTrkHits_endcap_unmatched",101,-0.5,100.5);
  TH1F *h_muonTrkIso_endcap_unmatched = new TH1F("muonTrkIso_endcap_unmatched","muonTrkIso_endcap_unmatched",200,-0.02,0.02);

  TH1F *h_muonEff_Pt = new TH1F ("muonEff_Pt","muonEff_Pt",100,0,1000);
  TH1F *h_muonEff_Eta = new TH1F ("muonEff_Eta","muonEff_Eta",121,-6.05,6.05);
  TH1F *h_muonEff_Pt_ID = new TH1F ("muonEff_Pt_ID","muonEff_Pt_ID",100,0,1000);
  TH1F *h_muonEff_Eta_ID = new TH1F ("muonEff_Eta_ID","muonEff_Eta_ID",121,-6.05,6.05);
  TH1F *h_muonEff_Pt_ID_ISO = new TH1F ("muonEff_Pt_ID_ISO","muonEff_Pt_ID_ISO",100,0,1000);
  TH1F *h_muonEff_Eta_ID_ISO = new TH1F ("muonEff_Eta_ID_ISO","muonEff_Eta_ID_ISO",121,-6.05,6.05);


  // -  The muon ID is the official one for muons, the same that the TeV
  // muon group is using and has studied in their note. They have studied
  // the other cuts as well, but found them not so efficient and didn't use
  // them.
  // 1.  if (!(muonGlobalChi2[0]<5.0 && muonGlobalChi2[1]<5.0)  || !(
  // fabs(muonTrkD0[0])< 0.25 && fabs(muonTrkD0[1])< 0.25) ||
  // !(muonTrkHits[0]>=7 && muonTrkHits[1]>=7) ) continue;
  // 2. muonTrkIso[0]<10.0 && muonTrkIso[1]<10.0

  /////////initialize variables
  int muon_PID=int(getPreCutValue1("muonPID"));
  int LQ_PID=int(getPreCutValue1("leptoquarkPID"));
  float ConeSizeMCmatch_cut=getPreCutValue1("coneSizeMCmatchCut");

  /////////////Barrel Cuts ////////////////
  float muonGlobalChi2_barrel_cut=getPreCutValue1("ID_muonGlobalChi2_bar"); //<5.0 Oana
  float muonTrkD0_barrel_cut=getPreCutValue1("ID_muonTrkD0_bar"); //<0.25 Oana
  float muonTrkHits_barrel_cut=getPreCutValue1("ID_muonTrkHits_bar"); //>=7 Oana
  float muonTrkIso_barrel_cut=getPreCutValue1("ISO_muonTrkIso_bar"); // <10.0 Oana
  float muonEta_barrel_cut=getPreCutValue1("muonEta_bar");  //|eta|<1.2 ? waiting for Oana

  /////////////End Cap Cuts ////////////////
  float muonGlobalChi2_endcap_cut=getPreCutValue1("ID_muonGlobalChi2_end"); //<5.0 Oana
  float muonTrkD0_endcap_cut=getPreCutValue1("ID_muonTrkD0_end"); //<0.25 Oana
  float muonTrkHits_endcap_cut=getPreCutValue1("ID_muonTrkHits_end"); //>=7 Oana
  float muonTrkIso_endcap_cut=getPreCutValue1("ISO_muonTrkIso_end"); // <10.0 Oana
  float muonEta_endcap_cut=getPreCutValue1("muonEta_end");  //|eta|<1.2 ? waiting for Oana
  float muonEta_endcap_cut2=getPreCutValue2("muonEta_end");  //|eta|>1.2 && |eta|<2.4 ? waiting for Oana

  float muonPt_cut=getPreCutValue1("muonPtCut");
  float muonEta_cut=getPreCutValue1("muonEtaCut");
  float genPartPt_cut=getPreCutValue1("genPartPtCut");

  //-----

  int n_muon_doubleMatched = 0; //counter

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

    int nMuonMatched = 0;
    int nGenMuon = 0;
    bool FirstLoop = true;
    int from_LQ_idx[2]={99,99};
    int matched_reco_idx[2]={99,99};
    int nMuon_pre = 0;
    int nMuon_post_ID = 0;
    int nMuon_post_ID_ISO = 0;

    //////////Finding muonc from LQ
    for(int igen=0;igen<GenParticleCount;igen++)
      {
 
	//muons from LQ decay
	if(abs(GenParticlePdgId[igen])==muon_PID
	   && abs(GenParticlePdgId[GenParticleMotherIndex[igen]])==LQ_PID)
	  {
	    nGenMuon++;
	    h_muon_Pt_Gen->Fill(GenParticlePt[igen]);
	    h_muon_Eta_Gen->Fill(GenParticleEta[igen]);
	    h_muon_Eta_Pt_Gen->Fill(GenParticlePt[igen],GenParticleEta[igen]);
	    if (from_LQ_idx[0]==99) from_LQ_idx[0]=igen;
	    else from_LQ_idx[1]=igen;
	    if (fabs(GenParticleEta[igen])<muonEta_cut){
	      h_muon_Pt_Gen_etaCut->Fill(GenParticlePt[igen]);
	      h_muon_Eta_Gen_etaCut->Fill(GenParticleEta[igen]);
	    }
	    if (GenParticlePt[igen]<genPartPt_cut)
	      h_muon_Eta_Gen_lowPt->Fill(GenParticleEta[igen]);
	  }
      }
    
    TVector3 muon1gen, muon2gen;
    muon1gen.SetPtEtaPhi(GenParticlePt[from_LQ_idx[0]],
			 GenParticleEta[from_LQ_idx[0]],
			 GenParticlePhi[from_LQ_idx[0]]);	
    
    muon2gen.SetPtEtaPhi(GenParticlePt[from_LQ_idx[1]],
			 GenParticleEta[from_LQ_idx[1]],
			 GenParticlePhi[from_LQ_idx[1]]);	
    
    ///////Calculating DeltaR and matching
    float minDeltaR1 = 99;
    float minDeltaR2 = 99;
    float min2ndDeltaR1 = 99;
    float min2ndDeltaR2 = 99;
    int RecoIndex1 = 99;
    int RecoIndex2 = 99;
    
    for (int imuon=0; imuon<muonCount; imuon++) {
      if (muonPt[imuon]<muonPt_cut) continue;
      TVector3 muon;
      muon.SetPtEtaPhi(muonPt[imuon],
		       muonEta[imuon],
		       muonPhi[imuon]);	
      
      float DeltaR_muon_muon1gen = muon1gen.DeltaR(muon);
      float DeltaR_muon_muon2gen = muon2gen.DeltaR(muon);

      if (DeltaR_muon_muon1gen < minDeltaR1) {
	min2ndDeltaR1=minDeltaR1;
	minDeltaR1=DeltaR_muon_muon1gen;
	RecoIndex1 = imuon;
      }
      else if (DeltaR_muon_muon1gen < min2ndDeltaR1){
	min2ndDeltaR1=DeltaR_muon_muon1gen;
      }
      if (DeltaR_muon_muon2gen < minDeltaR2) {
	min2ndDeltaR2=minDeltaR2;
	minDeltaR2=DeltaR_muon_muon2gen;
	RecoIndex2 = imuon;
      }
      else if (DeltaR_muon_muon2gen < min2ndDeltaR2){
	min2ndDeltaR2=DeltaR_muon_muon2gen;
      }
    }
    h_DeltaR_Gen_Reco->Fill(minDeltaR1);
    h_DeltaR_Gen_Reco->Fill(minDeltaR2);
    h_DeltaR_Gen_2ndReco->Fill(min2ndDeltaR1);
    h_DeltaR_Gen_2ndReco->Fill(min2ndDeltaR2);
    if (minDeltaR1< ConeSizeMCmatch_cut) matched_reco_idx[0]=RecoIndex1;
    if (minDeltaR2< ConeSizeMCmatch_cut) matched_reco_idx[1]=RecoIndex2;
    //      cout << "GenMuonc: " << from_LQ_idx[0] << " matched to: " << matched_reco_idx[0] << endl;
    //      cout << "GenMuonc: " << from_LQ_idx[1] << " matched to: " << matched_reco_idx[1] << endl;
    //      cout << "////////////////////" << endl;


    /////////Fill histograms appropriately
    for (int imuon=0; imuon<muonCount; imuon++) {
      if (muonPt[imuon]<muonPt_cut) continue;
      nMuon_pre++;
      bool MuonIsMatched = false;
      bool InBarrel = false;
      bool InEndCap = false;
      bool PassID = false;
      bool PassISO = false;
      int igen_matched = -1;

      h_muon_E->Fill(muonEnergy[imuon]);
      h_muon_Pt->Fill(muonPt[imuon]);
      h_muon_Phi->Fill(muonPhi[imuon]);
      h_muon_Eta->Fill(muonEta[imuon]);

      ///////////do matching 
      TVector3 muon;
      muon.SetPtEtaPhi(muonPt[imuon],
		      muonEta[imuon],
		      muonPhi[imuon]);	
 
      float DeltaR_muon_muon1gen = muon1gen.DeltaR(muon);
      float DeltaR_muon_muon2gen = muon2gen.DeltaR(muon);
 
      if((imuon==matched_reco_idx[0])||(imuon==matched_reco_idx[1])) MuonIsMatched = true;
      if(imuon==matched_reco_idx[0])igen_matched=from_LQ_idx[0] ;
      if(imuon==matched_reco_idx[1])igen_matched=from_LQ_idx[1] ;
      if ((DeltaR_muon_muon1gen<ConeSizeMCmatch_cut)&&(DeltaR_muon_muon2gen<ConeSizeMCmatch_cut)) n_muon_doubleMatched++;
 
      if (MuonIsMatched == true){
	h_muon_Pt_Gen_matched->Fill(GenParticlePt[igen_matched]);
	h_muon_Eta_Gen_matched->Fill(GenParticleEta[igen_matched]);
	if (GenParticlePt[igen_matched]<genPartPt_cut)
	  h_muon_Eta_Gen_matched_lowPt->Fill(GenParticleEta[igen_matched]);
	h_muon_E_matched->Fill(muonEnergy[imuon]);
	h_muon_Pt_matched->Fill(muonPt[imuon]);
	h_muon_Phi_matched->Fill(muonPhi[imuon]);
	h_muon_Eta_matched->Fill(muonEta[imuon]);
	nMuonMatched++;
	h_Energy_Res->Fill(muonEnergy[imuon]/GenParticleE[igen_matched]);
	//cout << "GenMuon " << igen_matched << " eta: " << GenParticleEta[igen_matched] << " is match to RecoMuon " << imuon << " eta: " << muonEta[imuon] << endl;
      }

      //check for barrel or endcap
      float muon_Eta = fabs(muonEta[imuon]);
      if ( muon_Eta < muonEta_barrel_cut )
	InBarrel = true;
      if ( (muon_Eta > muonEta_endcap_cut ) && (muon_Eta < muonEta_endcap_cut2) )
	InEndCap = true;
 
      ////matched
      if (MuonIsMatched){
	if (InBarrel){
	  h_muonGlobalChi2_barrel_matched->Fill(muonGlobalChi2[imuon]);
	  h_muonTrkD0_barrel_matched->Fill(muonTrkD0[imuon]);
	  h_muonTrkHits_barrel_matched->Fill(muonTrkHits[imuon]);
	  h_muonTrkIso_barrel_matched->Fill(muonTrkIso[imuon]);
	  h_Energy_Res_barrel->Fill(muonEnergy[imuon]/GenParticleE[igen_matched]);
	}
  
	if (InEndCap){
	  h_muonGlobalChi2_endcap_matched->Fill(muonGlobalChi2[imuon]);
	  h_muonTrkD0_endcap_matched->Fill(muonTrkD0[imuon]);
	  h_muonTrkHits_endcap_matched->Fill(muonTrkHits[imuon]);
	  h_muonTrkIso_endcap_matched->Fill(muonTrkIso[imuon]);
	  h_Energy_Res_endcap->Fill(muonEnergy[imuon]/GenParticleE[igen_matched]);
	}
      }
       
      else { ////unmatched

	//cout << "Found unmatched muon" << endl;
	if (InBarrel){
	  h_muonGlobalChi2_barrel_unmatched->Fill(muonGlobalChi2[imuon]);
	  h_muonTrkD0_barrel_unmatched->Fill(muonTrkD0[imuon]);
	  h_muonTrkHits_barrel_unmatched->Fill(muonTrkHits[imuon]);
	  h_muonTrkIso_barrel_unmatched->Fill(muonTrkIso[imuon]);
	}

	if (InEndCap){
	  h_muonGlobalChi2_endcap_unmatched->Fill(muonGlobalChi2[imuon]);
	  h_muonTrkD0_endcap_unmatched->Fill(muonTrkD0[imuon]);
	  h_muonTrkHits_endcap_unmatched->Fill(muonTrkHits[imuon]);
	  h_muonTrkIso_endcap_unmatched->Fill(muonTrkIso[imuon]);
	}

      } // end else unmatched


      ///check for ID 
      if (InBarrel)
	if ((muonGlobalChi2[imuon]<muonGlobalChi2_barrel_cut)
	    && (muonTrkD0[imuon]<muonTrkD0_barrel_cut)
	    && (muonTrkHits[imuon]>=muonTrkHits_barrel_cut)
	    && (muonTrkHits[imuon]>=muonTrkHits_barrel_cut)
	    )
	  PassID=true;
 
      if (InEndCap)
	if ((muonGlobalChi2[imuon]<muonGlobalChi2_endcap_cut)
	    && (muonTrkD0[imuon]<muonTrkD0_endcap_cut)
	    && (muonTrkHits[imuon]>=muonTrkHits_endcap_cut)
	    && (muonTrkHits[imuon]>=muonTrkHits_endcap_cut)
	    )
	  PassID=true;
 
      if (PassID) nMuon_post_ID++;
 
      if ((MuonIsMatched)&&(PassID)){
	h_muon_Pt_Gen_matched_ID->Fill(GenParticlePt[igen_matched]);
	h_muon_Eta_Gen_matched_ID->Fill(GenParticleEta[igen_matched]);
      }

      ///check for ISO
      if (InBarrel)
	if ((muonTrkIso[imuon]<muonTrkIso_barrel_cut)
	    )
	  PassISO=true;
      
      if (InEndCap)
	if ((muonTrkIso[imuon]<muonTrkIso_endcap_cut)
	    )
	  PassISO=true;
      
      if ((MuonIsMatched)&&(PassID)&&(PassISO)){
	h_muon_Pt_Gen_matched_ID_ISO->Fill(GenParticlePt[igen_matched]);
	h_muon_Eta_Gen_matched_ID_ISO->Fill(GenParticleEta[igen_matched]);
      }
      
      if ((PassISO)&&(PassID)) {
	nMuon_post_ID_ISO++;
	if (FirstLoop){
	  h_muon_E_ID_ISO->Fill(muonEnergy[imuon]);
	  h_muon_Pt_ID_ISO->Fill(muonPt[imuon]);
	  h_muon_Phi_ID_ISO->Fill(muonPhi[imuon]);
	  h_muon_Eta_ID_ISO->Fill(muonEta[imuon]);
	}
      }
    }// end muon loop
 
    h_N_muon_Pre->Fill(nMuon_pre);
    h_N_muon_Post_ID->Fill(nMuon_post_ID);
    h_N_muon_Post_ID_ISO->Fill(nMuon_post_ID_ISO);
 
    h_N_muon_matched->Fill(nMuonMatched);
    h_N_muon_Gen->Fill(nGenMuon);

    ////////////////////// User's code ends here ///////////////////////

  } // End loop over events

  //////////write histos 
  h_muon_E->Write();
  h_muon_Pt->Write();
  h_muon_Phi->Write();
  h_muon_Eta->Write();

  h_N_muon_matched->Write();
  h_muon_E_matched->Write();
  h_muon_Pt_matched->Write();
  h_muon_Phi_matched->Write();
  h_muon_Eta_matched->Write();

  h_muon_E_ID_ISO->Write();
  h_muon_Pt_ID_ISO->Write();
  h_muon_Phi_ID_ISO->Write();
  h_muon_Eta_ID_ISO->Write();

  h_N_muon_Gen->Write();
  h_muon_Pt_Gen->Write();
  h_muon_Pt_Gen_etaCut->Write();
  h_muon_Pt_Gen_matched->Write();
  h_muon_Pt_Gen_matched_ID->Write();
  h_muon_Pt_Gen_matched_ID_ISO->Write();
  h_muon_Eta_Gen->Write();
  h_muon_Eta_Gen_lowPt->Write();
  h_muon_Eta_Gen_etaCut->Write();
  h_muon_Eta_Gen_matched->Write();
  h_muon_Eta_Gen_matched_lowPt->Write();
  h_muon_Eta_Gen_matched_ID->Write();
  h_muon_Eta_Gen_matched_ID_ISO->Write();
  h_muon_Eta_Pt_Gen->Write();

  h_DeltaR_Gen_Reco->Write();
  h_DeltaR_Gen_2ndReco->Write();

  h_N_muon_Pre->Write();
  h_N_muon_Post_ID->Write();
  h_N_muon_Post_ID_ISO->Write();

  h_Energy_Res->Write();
  h_Energy_Res_barrel->Write();
  h_Energy_Res_endcap->Write();

  //////matched
  h_muonGlobalChi2_barrel_matched->Write();
  h_muonTrkD0_barrel_matched->Write();
  h_muonTrkHits_barrel_matched->Write();
  h_muonTrkIso_barrel_matched->Write();

  h_muonGlobalChi2_endcap_matched->Write();
  h_muonTrkD0_endcap_matched->Write();
  h_muonTrkHits_endcap_matched->Write();
  h_muonTrkIso_endcap_matched->Write();

  //////unmatched
  h_muonGlobalChi2_barrel_unmatched->Write();
  h_muonTrkD0_barrel_unmatched->Write();
  h_muonTrkHits_barrel_unmatched->Write();
  h_muonTrkIso_barrel_unmatched->Write();

  h_muonGlobalChi2_endcap_unmatched->Write();
  h_muonTrkD0_endcap_unmatched->Write();
  h_muonTrkHits_endcap_unmatched->Write();
  h_muonTrkIso_endcap_unmatched->Write();

  h_muonEff_Pt->Divide(h_muon_Pt_Gen_matched,h_muon_Pt_Gen,1,1);
  h_muonEff_Pt->Write();
  h_muonEff_Eta->Divide(h_muon_Eta_Gen_matched,h_muon_Eta_Gen,1,1);
  h_muonEff_Eta->Write();

  h_muonEff_Pt_ID->Divide(h_muon_Pt_Gen_matched_ID,h_muon_Pt_Gen,1,1);
  h_muonEff_Pt_ID->Write();
  h_muonEff_Eta_ID->Divide(h_muon_Eta_Gen_matched_ID,h_muon_Eta_Gen,1,1);
  h_muonEff_Eta_ID->Write();

  h_muonEff_Pt_ID_ISO->Divide(h_muon_Pt_Gen_matched_ID_ISO,h_muon_Pt_Gen,1,1);
  h_muonEff_Pt_ID_ISO->Write();
  h_muonEff_Eta_ID_ISO->Divide(h_muon_Eta_Gen_matched_ID_ISO,h_muon_Eta_Gen,1,1);
  h_muonEff_Eta_ID_ISO->Write();

  cout << "Number of events where recoMuon matches BOTH genMuon: " << n_muon_doubleMatched << endl;

  std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
