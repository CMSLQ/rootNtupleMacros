#define analysisClass_cxx
#include "analysisClass.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "vector"

//############################################################################################
// Plots to be done :
// - take the energy of the second electron after the corrections 
//  LQ --> e+q  ; e --> e+g
//############################################################################################

int LQ_pdgID = 42;
int Ele_pdgID = 11;
int Muon_pdgID = 13;
int QuarkD_pdgID = 1;
int QuarkU_pdgID = 2;
int QuarkS_pdgID = 3;
int QuarkC_pdgID = 4;
int QuarkB_pdgID = 5;
int QuarkT_pdgID = 6;

float cut__deltaRmin_genJet_ele = 0.5;
float cut__deltaRmin_genJet_quark = 0.5;
float cut__genJetPtForMatch = 0;

analysisClass::analysisClass(string * inputList, string * treeName, TString * outputFileName)
  :baseClass(inputList, treeName, outputFileName)
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

   //## Book histograms
   TH1F *h_pthat = new TH1F("h_pthat","pthat",200,0,2000);
   TH1F *h_processID = new TH1F("h_processID","processID",501,-0.5,500.5);
   TH1F *h_GenParticlePdgIdAll = new TH1F("h_GenParticlePdgIdAll","GenParticlePdgIdAll",10001,-5000.5,5000.5);
   TH1F *h_NumLQPerEvent = new TH1F("h_NumLQPerEvent","NumLQPerEvent",4,-0.5,3.5);
   TH1F *h_NumEleFromLQPerEvent = new TH1F("h_NumEleFromLQPerEvent","NumEleFromLQPerEvent",4,-0.5,3.5);
   TH1F *h_NumQuarkFromLQPerEvent = new TH1F("h_NumQuarkFromLQPerEvent","NumQuarkFromLQPerEvent",4,-0.5,3.5);
   TH1F *h_MassLQeqMinusMassLQ = new TH1F("h_MassLQeqMinusMassLQ","MassLQeqMinusMassLQ",2100,-10.5,10.5);

   TH1F *h_pTLQ = new TH1F("h_pTLQ","pTLQ",200,0,2000);
   TH1F *h_EnergyLQ = new TH1F("h_EnergyLQ","EnergyLQ",200,0,2000);
   TH1F *h_EtaLQ = new TH1F("h_EtaLQ","EtaLQ",100,-6,6);
   TH1F *h_PhiLQ = new TH1F("h_PhiLQ","PhiLQ",100,-TMath::Pi(),TMath::Pi());
   TH2F *h_EnergyVsMassLQ = new TH2F("h_EnergyVsMassLQ","EnergyVsMassLQ",200,0,2000,100,0,1000);

   TH1F *h_pTEleFromLQ = new TH1F("h_pTEleFromLQ","pTEleFromLQ",200,0,2000);
   TH1F *h_EnergyEleFromLQ = new TH1F("h_EnergyEleFromLQ","EnergyEleFromLQ",200,0,2000);
   TH1F *h_EtaEleFromLQ = new TH1F("h_EtaEleFromLQ","EtaEleFromLQ",100,-6,6);
   TH1F *h_PhiEleFromLQ = new TH1F("h_PhiEleFromLQ","PhiEleFromLQ",100,-TMath::Pi(),TMath::Pi());

   TH1F *h_pTQuarkFromLQ = new TH1F("h_pTQuarkFromLQ","pTQuarkFromLQ",200,0,2000);
   TH1F *h_EnergyQuarkFromLQ = new TH1F("h_EnergyQuarkFromLQ","EnergyQuarkFromLQ",200,0,2000);
   TH1F *h_EtaQuarkFromLQ = new TH1F("h_EtaQuarkFromLQ","EtaQuarkFromLQ",100,-6,6);
   TH1F *h_PhiQuarkFromLQ = new TH1F("h_PhiQuarkFromLQ","PhiQuarkFromLQ",100,-TMath::Pi(),TMath::Pi());

   TH1F *h_deltaRmin_genJet_ele = new TH1F("h_deltaRmin_genJet_ele","deltaRmin_genJet_ele",1000,0,10);
   TH1F *h_deltaRmin_genJet_quark = new TH1F("h_deltaRmin_genJet_quark","deltaRmin_genJet_quark",1000,0,10);
   TH1F *h_deltaRMin_genJet_genJetNearestEle = new TH1F("h_deltaRMin_genJet_genJetNearestEle",
							"h_deltaRMin_genJet_genJetNearestEle"
							,1000,0,10);
   TH1F *h_deltaRMin_genJet_genJetNearestQuark = new TH1F("h_deltaRMin_genJet_genJetNearestQuark",
							  "h_deltaRMin_genJet_genJetNearestQuark"
							  ,1000,0,10);
   
   TH1F *h_EgenJetOverEgenEle_ForDeltaRmin = new TH1F("h_EgenJetOverEgenEle_ForDeltaRmin","EgenJetOverEgenEle_ForDeltaRmin",1000,0,2);
   TH1F *h_EgenJetOverEgenQuark_ForDeltaRmin = new TH1F("h_EgenJetOverEgenQuark_ForDeltaRmin","EgenJetOverEgenQuark_ForDeltaRmin",1000,0,2);

   TH1F *h_EgenJetOverEgenEle_MatchCands = new TH1F("h_EgenJetOverEgenEle_MatchCands","h_EgenJetOverEgenEle_MatchCands",1000,0,2);
   TH1F *h_EgenJetOverEgenQuark_MatchCands = new TH1F("h_EgenJetOverEgenQuark_MatchCands","h_EgenJetOverEgenQuark_MatchCands",1000,0,2);

   TH1F *h_NumGenJetMatchEleFromLQ = new TH1F("h_NumGenJetMatchEleFromLQ","h_NumGenJetMatchEleFromLQ",21,-0.5,20.5);
   TH1F *h_NumGenJetMatchQuarkFromLQ = new TH1F("h_NumGenJetMatchQuarkFromLQ","h_NumGenJetMatchQuarkFromLQ",21,-0.5,20.5);
   TH1F *h_NumGenJetNoMatchFromLQ = new TH1F("h_NumGenJetNoMatchFromLQ","h_NumGenJetNoMatchFromLQ",21,-0.5,20.5);

   TH1F *h_PxTwoLQSystem = new TH1F("h_PxTwoLQSystem","PxTwoLQSystem",1001,-500.5,500.5);
   TH1F *h_PyTwoLQSystem = new TH1F("h_PyTwoLQSystem","PyTwoLQSystem",1001,-500.5,500.5);
   TH1F *h_PzTwoLQSystem = new TH1F("h_PzTwoLQSystem","PzTwoLQSystem",1001,-500.5,500.5);

   TH1F *h_PxEventRecoil = new TH1F("h_PxEventRecoil","PxEventRecoil",1001,-500.5,500.5);
   TH1F *h_PyEventRecoil = new TH1F("h_PyEventRecoil","PyEventRecoil",1001,-500.5,500.5);
   TH1F *h_PzEventRecoil = new TH1F("h_PzEventRecoil","PzEventRecoil",1001,-500.5,500.5);

   TH2F *h_PxBalance = new TH2F("h_PxBalance","PxBalance",101,-500.5,500.5,101,-500.5,500.5);
   TH2F *h_PyBalance = new TH2F("h_PyBalance","PyBalance",101,-500.5,500.5,101,-500.5,500.5);
   TH2F *h_PzBalance = new TH2F("h_PzBalance","PzBalance",101,-500.5,500.5,101,-500.5,500.5);

   TH2F *h_PxBalance1stGenJetNoMatch = new TH2F("h_PxBalance1stGenJetNoMatch","PxBalance1stGenJetNoMatch"
						,101,-500.5,500.5,101,-500.5,500.5);
   TH2F *h_PyBalance1stGenJetNoMatch = new TH2F("h_PyBalance1stGenJetNoMatch","PyBalance1stGenJetNoMatch"
						,101,-500.5,500.5,101,-500.5,500.5);
   TH2F *h_PzBalance1stGenJetNoMatch = new TH2F("h_PzBalance1stGenJetNoMatch","PzBalance1stGenJetNoMatch"
						,101,-500.5,500.5,101,-500.5,500.5);

   TH2F *h_PxBalance2ndGenJetNoMatch = new TH2F("h_PxBalance2ndGenJetNoMatch","PxBalance2ndGenJetNoMatch"
						,101,-500.5,500.5,101,-500.5,500.5);
   TH2F *h_PyBalance2ndGenJetNoMatch = new TH2F("h_PyBalance2ndGenJetNoMatch","PyBalance2ndGenJetNoMatch"
						,101,-500.5,500.5,101,-500.5,500.5);
   TH2F *h_PzBalance2ndGenJetNoMatch = new TH2F("h_PzBalance2ndGenJetNoMatch","PzBalance2ndGenJetNoMatch"
						,101,-500.5,500.5,101,-500.5,500.5);

   TH1F *h_DeltaPt_genJetMinPtQuarkMatchFromLQ_genJet1stPtNoMatch = new TH1F("h_DeltaPt_genJetMinPtQuarkMatchFromLQ_genJet1stPtNoMatch"
									,"DeltaPt_genJetMinPtQuarkMatchFromLQ_genJet1stPtNoMatch"
									,400,-2000,2000);


   //====================================================


   //## Initialize variables

   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<1000;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%100 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////


     //## index of gen particles for LQ / eleFromLQ / quarkFromLQ

     vector<int> v_idx_LQ;
     vector<int> v_idx_elefromLQ;
     vector<int> v_idx_quarkfromLQ;


     //## loop over gen particles

     for(int genPart=0; genPart<GenParticleCount; genPart++)
       {

	 // print the ptyhia tree
	 // 	 cout << "index: " << genPart <<  " , id: " << GenParticlePdgId[genPart] 
	 // 	      << " , mother_index: " << GenParticleMotherIndex[genPart] 
	 // 	      << " , mother_id: " << GenParticlePdgId[GenParticleMotherIndex[genPart]] 
	 // 	      << " , energy: " << GenParticleE[genPart] << endl;

	 bool IsLQ=false;
	 bool IsEleFromLQ=false;
	 bool IsQuarkFromLQ=false;

	 h_GenParticlePdgIdAll->Fill(GenParticlePdgId[genPart]);

	 if( abs( GenParticlePdgId[genPart] ) == LQ_pdgID )
	   {
	     IsLQ=true;
	     v_idx_LQ.push_back(genPart);
	   }
	 
	 if(abs( GenParticlePdgId[genPart] ) == Ele_pdgID && 
	    abs( GenParticlePdgId[GenParticleMotherIndex[genPart]] ) == LQ_pdgID  
	    )
	   {
	     IsEleFromLQ=true;
	     v_idx_elefromLQ.push_back(genPart);
	   }

	 if(abs( GenParticlePdgId[genPart] ) >= QuarkD_pdgID && 
	    abs( GenParticlePdgId[genPart] ) <= QuarkT_pdgID &&
	    abs( GenParticlePdgId[GenParticleMotherIndex[genPart]] ) == LQ_pdgID  
	    )
	   {
	     IsQuarkFromLQ=true;
	     v_idx_quarkfromLQ.push_back(genPart);
	   }


	 //## Fill histograms per gen particle

	 if(IsLQ)
	   {
	     h_pTLQ->Fill(GenParticlePt[genPart]);
	     h_EnergyLQ->Fill(GenParticleE[genPart]);
	     h_EtaLQ->Fill(GenParticleEta[genPart]);
	     h_PhiLQ->Fill(GenParticlePhi[genPart]);
	   }

	 if(IsEleFromLQ)
	   {
	     h_pTEleFromLQ->Fill(GenParticlePt[genPart]);
	     h_EnergyEleFromLQ->Fill(GenParticleE[genPart]);
	     h_EtaEleFromLQ->Fill(GenParticleEta[genPart]);
	     h_PhiEleFromLQ->Fill(GenParticlePhi[genPart]);
	   }

	 if(IsQuarkFromLQ)
	   {
	     h_pTQuarkFromLQ->Fill(GenParticlePt[genPart]);
	     h_EnergyQuarkFromLQ->Fill(GenParticleE[genPart]);
	     h_EtaQuarkFromLQ->Fill(GenParticleEta[genPart]);
	     h_PhiQuarkFromLQ->Fill(GenParticlePhi[genPart]);
	   }

       }


     //## TLorentz vectors for LQ / eleFromLQ / quarkFromLQ
     
     vector<TLorentzVector> v_TLV_LQ;
     vector<TLorentzVector> v_TLV_eleFromLQ;
     vector<TLorentzVector> v_TLV_quarkFromLQ;

     for(int particle=0; particle < v_idx_elefromLQ.size(); particle++)
       {

	 TLorentzVector v_TLV_LQTemp;
	 TLorentzVector TLV_eleTemp;
	 TLorentzVector TLV_quarkTemp;

	 v_TLV_LQTemp.SetPtEtaPhiE(GenParticlePt[v_idx_LQ[particle]],
				 GenParticleEta[v_idx_LQ[particle]],
				 GenParticlePhi[v_idx_LQ[particle]],
				 GenParticleE[v_idx_LQ[particle]]);  

	 TLV_eleTemp.SetPtEtaPhiE(GenParticlePt[v_idx_elefromLQ[particle]],
				  GenParticleEta[v_idx_elefromLQ[particle]],
				  GenParticlePhi[v_idx_elefromLQ[particle]],
				  GenParticleE[v_idx_elefromLQ[particle]]);  

	 TLV_quarkTemp.SetPtEtaPhiE(GenParticlePt[v_idx_quarkfromLQ[particle]],
				  GenParticleEta[v_idx_quarkfromLQ[particle]],
				  GenParticlePhi[v_idx_quarkfromLQ[particle]],
				  GenParticleE[v_idx_quarkfromLQ[particle]]);  

	 v_TLV_LQ.push_back(v_TLV_LQTemp);
	 v_TLV_eleFromLQ.push_back(TLV_eleTemp);
	 v_TLV_quarkFromLQ.push_back(TLV_quarkTemp);

       }


     //## LQ mass

     float massLQ0_eq = (v_TLV_eleFromLQ[0]+v_TLV_quarkFromLQ[0]).M();
     float massLQ1_eq = (v_TLV_eleFromLQ[1]+v_TLV_quarkFromLQ[1]).M();
     float massLQ0 = v_TLV_LQ[0].M();
     float massLQ1 = v_TLV_LQ[1].M();

     
     //## index of gen jets matched with gen particles eleFromLQ / quarkFromLQ
     
     vector<int> v_idx_genJetMatchEleFromLQ;
     vector<int> v_idx_genJetMatchQuarkFromLQ;
     vector<int> v_idx_genJetNoMatchFromLQ;

     //## match genJets with eleFromLQ 
 
     for(int ele=0; ele<v_TLV_eleFromLQ.size(); ele++)
       {
	 
	 float deltaRmin_genJet_ele=99999.;
	 int idx_deltaRmin_genJet_ele=-1;

	 for(int genjet=0; genjet<genJetCount; genjet++)
	   {

	     if(genJetPt[genjet]<cut__genJetPtForMatch)
	       continue;

	     TLorentzVector TLV_genJetTemp;
	     TLV_genJetTemp.SetPtEtaPhiE(genJetPt[genjet],
					 genJetEta[genjet],
					 genJetPhi[genjet],
					 genJetEnergy[genjet]);  

	     float deltaR_genJet_ele = TLV_genJetTemp.DeltaR(v_TLV_eleFromLQ[ele]);
	     if(deltaR_genJet_ele < deltaRmin_genJet_ele)
	       {
		 deltaRmin_genJet_ele=deltaR_genJet_ele;
		 idx_deltaRmin_genJet_ele=genjet;
	       }
	   }

	 float EgenJetOverEgenEle = genJetEnergy[idx_deltaRmin_genJet_ele] / v_TLV_eleFromLQ[ele].E();

	 h_EgenJetOverEgenEle_ForDeltaRmin->Fill(EgenJetOverEgenEle);
	 h_deltaRmin_genJet_ele->Fill(deltaRmin_genJet_ele);

	 if(deltaRmin_genJet_ele < cut__deltaRmin_genJet_ele)
	   {
	     v_idx_genJetMatchEleFromLQ.push_back(idx_deltaRmin_genJet_ele);
	   }

	 //DeltaR between the 1st and 2nd jet more distant to the gen electron

	 float deltaRMin_genJet_genJetNearest=99999.;
	 int idx_deltaRMin_genJet_genJetNearest=-1;

	 for(int genjet=0; genjet<genJetCount; genjet++)
	   {

	     if(genjet==idx_deltaRmin_genJet_ele)
	       continue;

 	     if(genJetPt[genjet]<cut__genJetPtForMatch)
 	       continue;
	     
	     TLorentzVector TLV_genJetTemp;
	     TLV_genJetTemp.SetPtEtaPhiE(genJetPt[genjet],
					 genJetEta[genjet],
					 genJetPhi[genjet],
					 genJetEnergy[genjet]);  

	     TLorentzVector TLV_genJetNearest;
	     TLV_genJetNearest.SetPtEtaPhiE(genJetPt[idx_deltaRmin_genJet_ele],
					    genJetEta[idx_deltaRmin_genJet_ele],
					    genJetPhi[idx_deltaRmin_genJet_ele],
					    genJetEnergy[idx_deltaRmin_genJet_ele]);  

	     float deltaR_genJet_genJetNearest = TLV_genJetTemp.DeltaR(TLV_genJetNearest);
	     if(deltaR_genJet_genJetNearest < deltaRMin_genJet_genJetNearest)
	       {
		 deltaRMin_genJet_genJetNearest = deltaR_genJet_genJetNearest;
		 idx_deltaRMin_genJet_genJetNearest = genjet;
	       }

	   }

	 h_deltaRMin_genJet_genJetNearestEle->Fill(deltaRMin_genJet_genJetNearest);


       }


     //## match genJets with quarkFromLQ

     for(int quark=0; quark<v_TLV_quarkFromLQ.size(); quark++)
       {
	 
	 float deltaRmin_genJet_quark=99999.;
	 int idx_deltaRmin_genJet_quark=-1;

	 for(int genjet=0; genjet<genJetCount; genjet++)
	   {

	     if(genJetPt[genjet]<cut__genJetPtForMatch)
	       continue;

	     TLorentzVector TLV_genJetTemp;
	     TLV_genJetTemp.SetPtEtaPhiE(genJetPt[genjet],
					 genJetEta[genjet],
					 genJetPhi[genjet],
					 genJetEnergy[genjet]);  

	     float deltaR_genJet_quark = TLV_genJetTemp.DeltaR(v_TLV_quarkFromLQ[quark]);
	     if(deltaR_genJet_quark<deltaRmin_genJet_quark)
	       {
		 deltaRmin_genJet_quark=deltaR_genJet_quark;
		 idx_deltaRmin_genJet_quark=genjet;
	       }
	   }

	 float EgenJetOverEgenQuark = genJetEnergy[idx_deltaRmin_genJet_quark] / v_TLV_quarkFromLQ[quark].E();

	 h_EgenJetOverEgenQuark_ForDeltaRmin->Fill(EgenJetOverEgenQuark);
	 h_deltaRmin_genJet_quark->Fill(deltaRmin_genJet_quark);

	 if(deltaRmin_genJet_quark < cut__deltaRmin_genJet_quark)
	   {
	     v_idx_genJetMatchQuarkFromLQ.push_back(idx_deltaRmin_genJet_quark);
	   }



	 //DeltaR between the 1st and 2nd jet more distant to the gen quark

	 float deltaRMin_genJet_genJetNearest=99999.;
	 int idx_deltaRMin_genJet_genJetNearest=-1;

	 for(int genjet=0; genjet<genJetCount; genjet++)
	   {

	     if(genjet==idx_deltaRmin_genJet_quark)
	       continue;

 	     if(genJetPt[genjet]<cut__genJetPtForMatch)
 	       continue;
	     
	     TLorentzVector TLV_genJetTemp;
	     TLV_genJetTemp.SetPtEtaPhiE(genJetPt[genjet],
					 genJetEta[genjet],
					 genJetPhi[genjet],
					 genJetEnergy[genjet]);  

	     TLorentzVector TLV_genJetNearest;
	     TLV_genJetNearest.SetPtEtaPhiE(genJetPt[idx_deltaRmin_genJet_quark],
					    genJetEta[idx_deltaRmin_genJet_quark],
					    genJetPhi[idx_deltaRmin_genJet_quark],
					    genJetEnergy[idx_deltaRmin_genJet_quark]);  

	     float deltaR_genJet_genJetNearest = TLV_genJetTemp.DeltaR(TLV_genJetNearest);
	     if(deltaR_genJet_genJetNearest < deltaRMin_genJet_genJetNearest)
	       {
		 deltaRMin_genJet_genJetNearest = deltaR_genJet_genJetNearest;
		 idx_deltaRMin_genJet_genJetNearest = genjet;
	       }

	   }

	 h_deltaRMin_genJet_genJetNearestQuark->Fill(deltaRMin_genJet_genJetNearest);

	 //

       }


     //## genJets not matched with eleFromLQ or quarkFromLQ

     vector<int> v_idx_genJetMatchEleQuarkFromLQ;
     for(int eleMatch=0; eleMatch<v_idx_genJetMatchEleFromLQ.size() ; eleMatch++)
       v_idx_genJetMatchEleQuarkFromLQ.push_back(v_idx_genJetMatchEleFromLQ[eleMatch]);
     for(int quarkMatch=0; quarkMatch<v_idx_genJetMatchQuarkFromLQ.size() ; quarkMatch++)
       v_idx_genJetMatchEleQuarkFromLQ.push_back(v_idx_genJetMatchQuarkFromLQ[quarkMatch]);

     for(int genjet=0; genjet<genJetCount; genjet++)
       {
	 bool IsMatch=false;
	 for(int anyMatch=0; anyMatch<v_idx_genJetMatchEleQuarkFromLQ.size() ; anyMatch++)
	   {
	     if( genjet == v_idx_genJetMatchEleQuarkFromLQ[anyMatch])
	       {
		 IsMatch=true;
		 break;
	       }
	   }
	 if(IsMatch==false)
	   v_idx_genJetNoMatchFromLQ.push_back(genjet);
       }

     //            //check
     //            for(int genjet=0; genjet<genJetCount; genjet++)
     //              cout << "all gen jets : " << genjet << endl; 
     //            for(int anyMatch=0; anyMatch<v_idx_genJetMatchEleQuarkFromLQ.size() ; anyMatch++)
     //              cout << "any match: " << v_idx_genJetMatchEleQuarkFromLQ[anyMatch] << endl; 
     //            for(int eleMatch=0; eleMatch<v_idx_genJetMatchEleFromLQ.size() ; eleMatch++)
     //              cout << "ele match: " << v_idx_genJetMatchEleFromLQ[eleMatch] << endl; 
     //            for(int quarkMatch=0; quarkMatch<v_idx_genJetMatchQuarkFromLQ.size() ; quarkMatch++)
     //              cout << "quark match: " << v_idx_genJetMatchQuarkFromLQ[quarkMatch] << endl; 
     //            for(int noMatch=0; noMatch<v_idx_genJetNoMatchFromLQ.size() ; noMatch++)
     //              cout << "gen jets no match: " << v_idx_genJetNoMatchFromLQ[noMatch] << endl; 


     //## EgenJetOverEgen(Ele/Quark) with genJets matched with gen Particles

     if( v_TLV_eleFromLQ.size() == v_idx_genJetMatchEleFromLQ.size() )
       {
	 for(int eleMatch=0; eleMatch<v_TLV_eleFromLQ.size() ; eleMatch++)
	   {
	     float EgenJetOverEgenEle_MatchCands 
	       = genJetEnergy[ v_idx_genJetMatchEleFromLQ[eleMatch] ] / v_TLV_eleFromLQ[eleMatch].E();
	     h_EgenJetOverEgenEle_MatchCands->Fill(EgenJetOverEgenEle_MatchCands);
	   }
       }

     if( v_TLV_quarkFromLQ.size() == v_idx_genJetMatchQuarkFromLQ.size() )
       {
	 for(int quarkMatch=0; quarkMatch<v_TLV_quarkFromLQ.size() ; quarkMatch++)
	   {
	     float EgenJetOverEgenQuark_MatchCands 
	       = genJetEnergy[ v_idx_genJetMatchQuarkFromLQ[quarkMatch] ] / v_TLV_quarkFromLQ[quarkMatch].E();
	     h_EgenJetOverEgenQuark_MatchCands->Fill(EgenJetOverEgenQuark_MatchCands);
	   }
       }


     //## calculate TLV event recoil using gen jets 
     //## and TLV for not matched genjets

     TLorentzVector TLV_EventRecoil;
     vector<TLorentzVector> v_TLV_genJetNoMatchFromLQ;

     for(int noMatch=0; noMatch<v_idx_genJetNoMatchFromLQ.size(); noMatch++)
       {
	 TLorentzVector TLV_genJetTemp;
	 TLV_genJetTemp.SetPtEtaPhiE(genJetPt[v_idx_genJetNoMatchFromLQ[noMatch]],
				     genJetEta[v_idx_genJetNoMatchFromLQ[noMatch]],
				     genJetPhi[v_idx_genJetNoMatchFromLQ[noMatch]],
				     genJetEnergy[v_idx_genJetNoMatchFromLQ[noMatch]]);  
	 
	 v_TLV_genJetNoMatchFromLQ.push_back(TLV_genJetTemp);

	 TLV_EventRecoil+=TLV_genJetTemp;

       }

     //## TLV for genjets matched with ele and quarks
     vector<TLorentzVector> v_TLV_genJetMatchEleFromLQ;
     vector<TLorentzVector> v_TLV_genJetMatchQuarkFromLQ;

     for(int ele=0; ele<v_idx_genJetMatchEleFromLQ.size(); ele++)
       {
	 TLorentzVector TLV_genJetTemp;
	 TLV_genJetTemp.SetPtEtaPhiE(genJetPt[v_idx_genJetMatchEleFromLQ[ele]],
				     genJetEta[v_idx_genJetMatchEleFromLQ[ele]],
				     genJetPhi[v_idx_genJetMatchEleFromLQ[ele]],
				     genJetEnergy[v_idx_genJetMatchEleFromLQ[ele]]);  
	 
	 v_TLV_genJetMatchEleFromLQ.push_back(TLV_genJetTemp);
       }

     for(int quark=0; quark<v_idx_genJetMatchEleFromLQ.size(); quark++)
       {
	 TLorentzVector TLV_genJetTemp;
	 TLV_genJetTemp.SetPtEtaPhiE(genJetPt[v_idx_genJetMatchEleFromLQ[quark]],
				     genJetEta[v_idx_genJetMatchEleFromLQ[quark]],
				     genJetPhi[v_idx_genJetMatchEleFromLQ[quark]],
				     genJetEnergy[v_idx_genJetMatchEleFromLQ[quark]]);  
	 
	 v_TLV_genJetMatchQuarkFromLQ.push_back(TLV_genJetTemp);
       }
     

     //## calculate LQ system
     TLorentzVector TLV_LQsystem;
     for(int LQ=0; LQ<v_TLV_LQ.size(); LQ++)
       TLV_LQsystem+=v_TLV_LQ[LQ];

     //## Compare Pt of mathced genJets from LQ quarks and un-matched ones
     float DeltaPt_genJetMinPtQuarkMatchFromLQ_genJet1stPtNoMatch;

     if(v_TLV_genJetNoMatchFromLQ.size() >=1 )
       {
	 float minPt=999;
	 for(int quark=0; quark < v_TLV_genJetMatchQuarkFromLQ.size(); quark++)
	   {
	     if(v_TLV_genJetMatchQuarkFromLQ[quark].Pt()<minPt)
	       minPt=v_TLV_genJetMatchQuarkFromLQ[quark].Pt();
	   }

	 DeltaPt_genJetMinPtQuarkMatchFromLQ_genJet1stPtNoMatch 
	   = minPt - v_TLV_genJetNoMatchFromLQ[0].Pt(); 
       }


     //## Fill histograms per event

     h_pthat->Fill(pthat);
     h_processID->Fill(processID);
     h_NumLQPerEvent->Fill(v_idx_LQ.size());
     h_NumEleFromLQPerEvent->Fill(v_idx_elefromLQ.size());
     h_NumQuarkFromLQPerEvent->Fill(v_idx_quarkfromLQ.size());

     h_NumGenJetMatchEleFromLQ->Fill(v_idx_genJetMatchEleFromLQ.size());
     h_NumGenJetMatchQuarkFromLQ->Fill(v_idx_genJetMatchQuarkFromLQ.size());
     h_NumGenJetNoMatchFromLQ->Fill(v_idx_genJetNoMatchFromLQ.size());

     h_EnergyVsMassLQ->Fill( v_TLV_LQ[0].E() , massLQ0 );
     h_EnergyVsMassLQ->Fill( v_TLV_LQ[1].E() , massLQ1 );
     h_MassLQeqMinusMassLQ->Fill( massLQ0_eq - massLQ0 );
     h_MassLQeqMinusMassLQ->Fill( massLQ1_eq - massLQ1 );


     if(v_TLV_quarkFromLQ.size() == v_idx_genJetMatchQuarkFromLQ.size() 
	&& v_TLV_eleFromLQ.size() == v_idx_genJetMatchEleFromLQ.size()
	&& v_TLV_quarkFromLQ.size() == v_TLV_eleFromLQ.size() 
	&& v_idx_genJetNoMatchFromLQ.size() >=2 )
       {

	 h_PxTwoLQSystem->Fill( TLV_LQsystem.Px() );
	 h_PyTwoLQSystem->Fill( TLV_LQsystem.Py() );
	 h_PzTwoLQSystem->Fill( TLV_LQsystem.Pz() );
	 
	 h_PxEventRecoil->Fill( TLV_EventRecoil.Px() );
	 h_PyEventRecoil->Fill( TLV_EventRecoil.Py() );
	 h_PzEventRecoil->Fill( TLV_EventRecoil.Pz() );
	 
	 h_PxBalance->Fill( TLV_LQsystem.Px() , TLV_EventRecoil.Px() );
	 h_PyBalance->Fill( TLV_LQsystem.Py() , TLV_EventRecoil.Py() );
	 h_PzBalance->Fill( TLV_LQsystem.Pz() , TLV_EventRecoil.Pz() );

	 h_PxBalance1stGenJetNoMatch->Fill( TLV_LQsystem.Px() , v_TLV_genJetNoMatchFromLQ[0].Px() );
	 h_PyBalance1stGenJetNoMatch->Fill( TLV_LQsystem.Py() , v_TLV_genJetNoMatchFromLQ[0].Py() );
	 h_PzBalance1stGenJetNoMatch->Fill( TLV_LQsystem.Pz() , v_TLV_genJetNoMatchFromLQ[0].Pz() );

	 h_PxBalance2ndGenJetNoMatch->Fill( TLV_LQsystem.Px() , v_TLV_genJetNoMatchFromLQ[1].Px() );
	 h_PyBalance2ndGenJetNoMatch->Fill( TLV_LQsystem.Py() , v_TLV_genJetNoMatchFromLQ[1].Py() );
	 h_PzBalance2ndGenJetNoMatch->Fill( TLV_LQsystem.Pz() , v_TLV_genJetNoMatchFromLQ[1].Pz() );

       }

     if(v_TLV_quarkFromLQ.size() == v_idx_genJetMatchQuarkFromLQ.size() 
	&& v_TLV_eleFromLQ.size() == v_idx_genJetMatchEleFromLQ.size()
	&& v_TLV_quarkFromLQ.size() == v_TLV_eleFromLQ.size() 
	&& v_idx_genJetNoMatchFromLQ.size() >=1 )
       {

	 h_DeltaPt_genJetMinPtQuarkMatchFromLQ_genJet1stPtNoMatch->Fill(DeltaPt_genJetMinPtQuarkMatchFromLQ_genJet1stPtNoMatch);

	 //printout
	 if(DeltaPt_genJetMinPtQuarkMatchFromLQ_genJet1stPtNoMatch<-200)
	   {
	     for(int genjet=0; genjet<genJetCount; genjet++)
	       {

		 bool IsMatchedWithEle=0;
		 bool IsMatchedWithQuark=0;

		 if(genjet == v_idx_genJetMatchEleFromLQ[0] || genjet == v_idx_genJetMatchEleFromLQ[1])
		   IsMatchedWithEle=1;

		 if(genjet == v_idx_genJetMatchQuarkFromLQ[0] || genjet == v_idx_genJetMatchQuarkFromLQ[1])
		   IsMatchedWithQuark=1;

 		 cout << "idx: " << genjet << " "  
 		      << "genJetPt: " << genJetPt[genjet] << " "  
 		      << "genJetEta: " << genJetEta[genjet] << " "  
 		      << "genJetPhi: " << genJetPhi[genjet] << " " 
 		      << "MatchedWithEleFromLQ?: " << IsMatchedWithEle << " "  
 		      << "MatchedWithQuarkFromLQ?: " << IsMatchedWithQuark << " "  
 		      << endl;

	       }

 	     cout << "genEle1: " << " "  
 		  << "genEle1Pt: " << v_TLV_eleFromLQ[0].Pt() << " "  
 		  << "genEle1Eta: " << v_TLV_eleFromLQ[0].Eta() << " "  
 		  << "genEle1Phi: " << v_TLV_eleFromLQ[0].Phi() << " " 
 		  << endl;

 	     cout << "genEle2: " << " "  
 		  << "genEle2Pt: " << v_TLV_eleFromLQ[1].Pt() << " "  
 		  << "genEle2Eta: " << v_TLV_eleFromLQ[1].Eta() << " "  
 		  << "genEle2Phi: " << v_TLV_eleFromLQ[1].Phi() << " " 
 		  << endl;

 	     cout << "genQuark1: " << " "  
 		  << "genQuark1Pt: " << v_TLV_quarkFromLQ[0].Pt() << " "  
 		  << "genQuark1Eta: " << v_TLV_quarkFromLQ[0].Eta() << " "  
 		  << "genQuark1Phi: " << v_TLV_quarkFromLQ[0].Phi() << " " 
 		  << endl;

 	     cout << "genQuark2: " << " "  
 		  << "genQuark2Pt: " << v_TLV_quarkFromLQ[1].Pt() << " "  
 		  << "genQuark2Eta: " << v_TLV_quarkFromLQ[1].Eta() << " "  
 		  << "genQuark2Phi: " << v_TLV_quarkFromLQ[1].Phi() << " " 
 		  << endl;

	     cout << endl;

	   }

       }

     ////////////////////// User's code ends here ///////////////////////
   }
   std::cout << "analysisClass::Loop() ends" <<std::endl;   


   //## Write histograms

   h_pthat->Write();
   h_processID->Write();
   h_GenParticlePdgIdAll->Write();
   h_NumLQPerEvent->Write();
   h_NumEleFromLQPerEvent->Write();
   h_NumQuarkFromLQPerEvent->Write();
   h_MassLQeqMinusMassLQ->Write();

   h_NumGenJetMatchEleFromLQ->Write();
   h_NumGenJetMatchQuarkFromLQ->Write();
   h_NumGenJetNoMatchFromLQ->Write();

   h_PxTwoLQSystem->Write();
   h_PyTwoLQSystem->Write();
   h_PzTwoLQSystem->Write();
   h_PxEventRecoil->Write();
   h_PyEventRecoil->Write();
   h_PzEventRecoil->Write();
   h_PxBalance->Write();
   h_PyBalance->Write();
   h_PzBalance->Write();
   h_PxBalance1stGenJetNoMatch->Write();
   h_PyBalance1stGenJetNoMatch->Write();
   h_PzBalance1stGenJetNoMatch->Write();
   h_PxBalance2ndGenJetNoMatch->Write();
   h_PyBalance2ndGenJetNoMatch->Write();
   h_PzBalance2ndGenJetNoMatch->Write();

   h_pTLQ->Write();
   h_EnergyLQ->Write();
   h_EtaLQ->Write();
   h_PhiLQ->Write();
   h_EnergyVsMassLQ->Write();

   h_pTEleFromLQ->Write();
   h_EnergyEleFromLQ->Write();
   h_EtaEleFromLQ->Write();
   h_PhiEleFromLQ->Write();
   
   h_pTQuarkFromLQ->Write();
   h_EnergyQuarkFromLQ->Write();
   h_EtaQuarkFromLQ->Write();
   h_PhiQuarkFromLQ->Write();

   h_deltaRmin_genJet_ele->Write();
   h_deltaRmin_genJet_quark->Write();

   h_EgenJetOverEgenEle_ForDeltaRmin->Write();
   h_EgenJetOverEgenQuark_ForDeltaRmin->Write();

   h_EgenJetOverEgenEle_MatchCands->Write();
   h_EgenJetOverEgenQuark_MatchCands->Write();

   h_deltaRMin_genJet_genJetNearestEle->Write();
   h_deltaRMin_genJet_genJetNearestQuark->Write();

   h_DeltaPt_genJetMinPtQuarkMatchFromLQ_genJet1stPtNoMatch->Write();

}
