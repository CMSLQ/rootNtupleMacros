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

   TH1F *h_N_muon_Gen = new TH1F ("N_muon_Gen","N_muon_Gen",5,-0.5,4.5);
   TH1F *h_muon_Pt_Gen  = new TH1F ("muon_Pt_Gen","muon_Pt_Gen",100,0,1000);
   TH1F *h_muon_Pt_Gen_etaCut  = new TH1F ("muon_Pt_Gen_etaCut","muon_Pt_Gen_etaCut",100,0,1000);
   TH1F *h_muon_Pt_Gen_matched  = new TH1F ("muon_Pt_Gen_matched","muon_Pt_Gen_matched",100,0,1000);
   TH1F *h_muon_Eta_Gen  = new TH1F ("muon_Eta_Gen","muon_Eta_Gen",121,-6.05,6.05);
   TH1F *h_muon_Eta_Gen_lowPt  = new TH1F ("muon_Eta_Gen_lowPt","muon_Eta_Gen_lowPt",121,-6.05,6.05);
   TH1F *h_muon_Eta_Gen_etaCut  = new TH1F ("muon_Eta_Gen_etaCut","muon_Eta_Gen_etaCut",121,-6.05,6.05);
   TH1F *h_muon_Eta_Gen_matched  = new TH1F ("muon_Eta_Gen_matched","muon_Eta_Gen_matched",121,-6.05,6.05);
   TH1F *h_muon_Eta_Gen_matched_lowPt  = new TH1F ("muon_Eta_Gen_matched_lowPt","muon_Eta_Gen_matched_lowPt",121,-6.05,6.05);
   TH2F *h_muon_Eta_Pt_Gen = new TH2F ("muon_Eta_Pt_Gen","muon_Eta_Pt_Gen",100,0,1000,121,-6.05,6.05);

   TH1F *h_DeltaR_Gen_Reco = new TH1F("DeltaR_Gen_Reco","DeltaR_Gen_Reco",100,0,0.2);
   TH1F *h_DeltaR_Gen_2ndReco = new TH1F("DeltaR_Gen_2ndReco","DeltaR_Gen_2ndReco",200,0,1.0);

   TH1F *h_Energy_Res = new TH1F ("Energy_Res","Energy_Res",150,0,1.5);
   TH1F *h_Energy_Res_barrel = new TH1F ("Energy_Res_barrel","Energy_Res_barrel",150,0,1.5);
   TH1F *h_Energy_Res_endcap = new TH1F ("Energy_Res_endcap","Energy_Res_endcap",150,0,1.5);

   TH1F *h_muonEff_Pt = new TH1F ("muonEff_Pt","muonEff_Pt",100,0,1000);
   TH1F *h_muonEff_Eta = new TH1F ("muonEff_Eta","muonEff_Eta",121,-6.05,6.05);


   /////////initialize variables
   //int electron_PID=11;
   int muon_PID=13;
   int LQ_PID=42;
   //int LQ_PID=23;
   float ConeSizeMCmatch_cut=0.07;

   float muonPt_cut=10.;
   float muonEta_cut=2.6;
   int n_muon_doubleMatched = 0;

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<100     ;jentry++) {
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

     //////////Finding elec from LQ
     for(int igen=0;igen<GenParticleCount;igen++)
       {
 
	 //select muons from LQ decay
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
	     if (GenParticlePt[igen]<40)
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
//      cout << "GenElec: " << from_LQ_idx[0] << " matched to: " << matched_reco_idx[0] << endl;
//      cout << "GenElec: " << from_LQ_idx[1] << " matched to: " << matched_reco_idx[1] << endl;
//      cout << "////////////////////" << endl;


     /////////Fill histograms appropriately
     for (int imuon=0; imuon<muonCount; imuon++) {
       if (muonPt[imuon]<muonPt_cut) continue;
       bool MuonIsMatched = false;
       bool InBarrel = false;
       bool InEndCap = false;
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
	 if (GenParticlePt[igen_matched]<40)
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
       float MuonEta = sqrt(muonEta[imuon]*muonEta[imuon]);
       if (MuonEta<1.422)InBarrel = true;
       if ((MuonEta>1.560)&&(MuonEta<2.5))InEndCap = true;
 
       ////matched
       if (MuonIsMatched){
         if (InBarrel){
	   h_Energy_Res_barrel->Fill(muonEnergy[imuon]/GenParticleE[igen_matched]);
         }
  
         if (InEndCap){
	   h_Energy_Res_endcap->Fill(muonEnergy[imuon]/GenParticleE[igen_matched]);
         }
       }
       
  }// end muon loop
 
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

   h_N_muon_Gen->Write();
   h_muon_Pt_Gen->Write();
   h_muon_Pt_Gen_etaCut->Write();
   h_muon_Pt_Gen_matched->Write();
   h_muon_Eta_Gen->Write();
   h_muon_Eta_Gen_lowPt->Write();
   h_muon_Eta_Gen_etaCut->Write();
   h_muon_Eta_Gen_matched->Write();
   h_muon_Eta_Gen_matched_lowPt->Write();
   h_muon_Eta_Pt_Gen->Write();

   h_DeltaR_Gen_Reco->Write();
   h_DeltaR_Gen_2ndReco->Write();

   h_Energy_Res->Write();
   h_Energy_Res_barrel->Write();
   h_Energy_Res_endcap->Write();

   h_muonEff_Pt->Divide(h_muon_Pt_Gen_matched,h_muon_Pt_Gen,1,1);
   h_muonEff_Pt->Write();
   h_muonEff_Eta->Divide(h_muon_Eta_Gen_matched,h_muon_Eta_Gen,1,1);
   h_muonEff_Eta->Write();

   cout << "Number of events where recoMuon matches BOTH genMuon: " << n_muon_doubleMatched << endl;

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
