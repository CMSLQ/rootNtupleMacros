#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
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

   //Combinations
   TH1F *h_Meq_gen = new TH1F ("Meq_gen","Meq_gen",200,0,800);  h_Meq_gen->Sumw2();
   TH1F *h_Meq_inside_gen = new TH1F ("Meq_inside_gen","Meq_inside_gen",200,0,800);  h_Meq_inside_gen->Sumw2();
   TH1F *h_Meq_above_gen = new TH1F ("Meq_above_gen","Meq_above_gen",200,0,800);  h_Meq_above_gen->Sumw2();

   TH2F *h_pTee_pTqq = new TH2F ("pTee_pTqq","pTee_pTqq",100,0,400,200,0,800);
   TH2F *h_pTee_pTqq_inside = new TH2F ("pTee_pTqq_inside","pTee_pTqq_inside",100,0,400,200,0,800);
   TH2F *h_pTee_pTqq_above = new TH2F ("pTee_pTqq_above","pTee_pTqq_above",100,0,400,200,0,800);

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<10;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     // Set the evaluation of the cuts to false and clear the variable values and filled status
     resetCuts();

     ///Gen Level quatities
     int pdgId_Mom = 23; //Z Boson
     TLorentzVector genEle1, genEle2, Z_vec;
     int N_Z_Ele=0;
     int idx_q_1=99, idx_q_2=99;
     int idx_e_1=99, idx_e_2=99;
     double STee=-9999;
     double Z_mass=-99, pTZ=-999;
     for (int igen=0; igen<GenParticleCount; igen++){
       if (GenParticlePt[igen]==0) continue;
       if (GenParticlePdgId[igen]==pdgId_Mom){
	 fillVariableWithValue("pTZ_gen",GenParticlePt[igen]);
	 Z_mass = sqrt((GenParticleE[igen]*GenParticleE[igen])-(GenParticleP[igen]*GenParticleP[igen]));
	 fillVariableWithValue("M_Z",Z_mass);
	 if ((Z_mass>80)&&(Z_mass<100))fillVariableWithValue("M_Z_inside",Z_mass);
	 if (Z_mass>100)fillVariableWithValue("M_Z_above",Z_mass);
       }
       if ((abs(GenParticlePdgId[igen])==11)&&
	   (GenParticlePdgId[GenParticleMotherIndex[igen]]==pdgId_Mom)){  // elec from Z
	 if (N_Z_Ele==0){
	   genEle1.SetPtEtaPhiM(GenParticlePt[igen],
				GenParticleEta[igen],
				GenParticlePhi[igen],0);
	   idx_e_1=igen;
	 }
	 else {
	   genEle2.SetPtEtaPhiM(GenParticlePt[igen],
				GenParticleEta[igen],
				GenParticlePhi[igen],0);
	   idx_e_2=igen;
	 }

	 N_Z_Ele++;
       }
     }
     if (N_Z_Ele>0){
       fillVariableWithValue("pTele1_gen",genEle1.Pt());
       fillVariableWithValue("eta_ele1_gen",genEle1.Eta());
       if ((Z_mass>80)&&(Z_mass<100)){
	 fillVariableWithValue("pTele1_gen_inside",genEle1.Pt());
	 fillVariableWithValue("eta_ele1_gen_inside",genEle1.Eta());
       }
       if (Z_mass>100){
	 fillVariableWithValue("pTele1_gen_above",genEle1.Pt());
	 fillVariableWithValue("eta_ele1_gen_above",genEle1.Eta());
       }

     }
     if (N_Z_Ele>1){
       Z_vec = genEle1 + genEle2;
       pTZ=Z_vec.Pt();
       fillVariableWithValue("pTZ_gen",pTZ);
       fillVariableWithValue("pTele2_gen",genEle2.Pt());
       fillVariableWithValue("eta_ele2_gen",genEle2.Eta());
       STee = genEle1.Pt()+genEle2.Pt();
       fillVariableWithValue("STee_gen",STee);
       if ((Z_mass>80)&&(Z_mass<100)){
	 fillVariableWithValue("pTZ_gen_inside",Z_vec.Pt());
	 fillVariableWithValue("pTele2_gen_inside",genEle2.Pt());
	 fillVariableWithValue("eta_ele2_gen_inside",genEle2.Eta());
	 fillVariableWithValue("STee_gen_inside",STee);
       }
       if (Z_mass>100){
	 fillVariableWithValue("pTZ_gen_above",Z_vec.Pt());
	 fillVariableWithValue("pTele2_gen_above",genEle2.Pt());
	 fillVariableWithValue("eta_ele2_gen_above",genEle2.Eta());
	 fillVariableWithValue("STee_gen_above",STee);
       }
     }

     //////Get GenJets and combine with e_1 and e_2
     int idx_genJet1=99, idx_genJet2=99, N_genJet=0;
     float deltaR_minCut = getPreCutValue1("jet_ele_DeltaRcut");
     for (int igenJet=0; igenJet<genJetCount; igenJet++){
       if (N_genJet>=2) break;
       if (genJetPt[igenJet]==0) continue;
	 float minDeltaR=9999;
	 TVector3 gen_jet_vec;
	 gen_jet_vec.SetPtEtaPhi(genJetPt[igenJet],genJetEta[igenJet],genJetPhi[igenJet]);
	 for (int i=0; i<GenParticleCount; i++){
	   if (abs(GenParticlePdgId[i])==11){
	     TVector3 gen_ele_vec;
	     gen_ele_vec.SetPtEtaPhi(GenParticlePt[i],GenParticleEta[i],GenParticlePhi[i]);
	     double distance = gen_jet_vec.DeltaR(gen_ele_vec);
	     if (distance<minDeltaR) minDeltaR=distance;
	   }//if electron
	 }// for all gen particles
	 if (minDeltaR > deltaR_minCut){
	   if (N_genJet==0) idx_genJet1=igenJet;
	   else idx_genJet2=igenJet;
	   N_genJet++;
	 } // if minDeltaR larger
     }//for all gen Jets

     TLorentzVector genJet1, genJet2;
     double M11_gen = -1, M12_gen = -1, M21_gen = -1, M22_gen = -1;
     double diff_11_22_gen, diff_12_21_gen;
     double  STqq=-99, pTqq=-999;
     if (N_genJet>0){
       genJet1.SetPtEtaPhiM(genJetPt[idx_genJet1],genJetEta[idx_genJet1],genJetPhi[idx_genJet1],0);
       fillVariableWithValue("pTq1_gen",genJet1.Pt());
       fillVariableWithValue("eta_q1_gen",genJet1.Eta());
       if ((Z_mass>80)&&(Z_mass<100)){
	 fillVariableWithValue("pTq1_gen_inside",genJet1.Pt());
	 fillVariableWithValue("eta_q1_gen_inside",genJet1.Eta());
       }
       if (Z_mass>100){
	 fillVariableWithValue("pTq1_gen_above",genJet1.Pt());
	 fillVariableWithValue("eta_q1_gen_above",genJet1.Eta());
       }
     }
     if (N_genJet>1){
       genJet2.SetPtEtaPhiM(genJetPt[idx_genJet2],genJetEta[idx_genJet2],genJetPhi[idx_genJet2],0);
       pTqq=(genJet1+genJet2).Pt();
       fillVariableWithValue("pTqq_gen",pTqq);
       fillVariableWithValue("pTq2_gen",genJet2.Pt());
       fillVariableWithValue("eta_q2_gen",genJet2.Eta());
       fillVariableWithValue("STqq_gen",genJet1.Pt()+genJet2.Pt() );
       if ((Z_mass>80)&&(Z_mass<100)){
	 fillVariableWithValue("pTqq_gen_inside",(genJet1+genJet2).Pt());
	 fillVariableWithValue("pTq2_gen_inside",genJet2.Pt());
	 fillVariableWithValue("eta_q2_gen_inside",genJet2.Eta());
	 fillVariableWithValue("STqq_gen_inside",genJet1.Pt()+genJet2.Pt() );
       }
       if (Z_mass>100){
	 fillVariableWithValue("pTqq_gen_above",(genJet1+genJet2).Pt());
	 fillVariableWithValue("pTq2_gen_above",genJet2.Pt());
	 fillVariableWithValue("eta_q2_gen_above",genJet2.Eta());
	 fillVariableWithValue("STqq_gen_above",genJet1.Pt()+genJet2.Pt() );
       }
     }
     if ((N_genJet>=2)&&(N_Z_Ele>=2)){
       h_pTee_pTqq->Fill(pTZ,pTqq);
       if ((Z_mass>80)&&(Z_mass<100))h_pTee_pTqq_inside->Fill(pTZ,pTqq);
       if (Z_mass>100) h_pTee_pTqq_above->Fill(pTZ,pTqq);
       fillVariableWithValue("sT",STee+STqq);
       fillVariableWithValue("sT_low",STee+STqq);
       if ((Z_mass>80)&&(Z_mass<100)) {
	 fillVariableWithValue("sT_inside",STee+STqq);
	 fillVariableWithValue("sT_low_inside",STee+STqq);
       }
       if (Z_mass>100){
	 fillVariableWithValue("sT_above",STee+STqq);
	 fillVariableWithValue("sT_low_above",STee+STqq);
       }
       TLorentzVector jet1ele1_gen, jet1ele2_gen, jet2ele1_gen, jet2ele2_gen;
       jet1ele1_gen= genEle1+genJet1;
       jet1ele2_gen= genEle2+genJet1;
       jet2ele1_gen= genEle1+genJet2;
       jet2ele2_gen= genEle2+genJet2;
       M11_gen = jet1ele1_gen.M();
       M12_gen = jet1ele2_gen.M();
       M21_gen = jet2ele1_gen.M();
       M22_gen = jet2ele2_gen.M();
       diff_11_22_gen=fabs(M11_gen-M22_gen);
       diff_12_21_gen=fabs(M12_gen-M21_gen);
	 if (diff_11_22_gen<diff_12_21_gen) 
	   {
	     h_Meq_gen->Fill(M11_gen);
	     h_Meq_gen->Fill(M22_gen);
	     if ((Z_mass>80)&&(Z_mass<100)){
	       h_Meq_inside_gen->Fill(M11_gen);
	       h_Meq_inside_gen->Fill(M22_gen);
	     }
	     if (Z_mass>100){
	       h_Meq_above_gen->Fill(M11_gen);
	       h_Meq_above_gen->Fill(M22_gen);
	     }
	   }
	 else 
	   {
	     h_Meq_gen->Fill(M21_gen);
	     h_Meq_gen->Fill(M12_gen);
	     if ((Z_mass>80)&&(Z_mass<100)){
	       h_Meq_inside_gen->Fill(M21_gen);
	       h_Meq_inside_gen->Fill(M12_gen);
	     }
	     if (Z_mass>100){
	       h_Meq_above_gen->Fill(M21_gen);
	       h_Meq_above_gen->Fill(M12_gen);
	     }
	   }

     }



     //--------------------------------------------------

     // Evaluate cuts (but do not apply them)
     evaluateCuts();

     // Fill histograms and do analysis based on cut evaluation

     if( passedCut("1") )
       {
       }
     
     if( (passedCut("1"))&&(passedCut("2")) )
       {
        }
     
     if( (passedCut("1"))&&(passedCut("3")) )
       {
       }
     
     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 
 h_Meq_gen->Write();            
 h_Meq_inside_gen->Write();     
 h_Meq_above_gen->Write();      
 			        
 h_pTee_pTqq->Write();                 
 h_pTee_pTqq_inside->Write();                 
 h_pTee_pTqq_above->Write();                 

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
