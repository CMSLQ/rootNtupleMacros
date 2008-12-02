#define analysisClass_cxx
#include "analysisClass.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "vector"

int LQ_pdgID=42;
int Ele_pdgID=11;
int Muon_pdgID=13;
int QuarkD_pdgID=1;
int QuarkU_pdgID=2;
int QuarkS_pdgID=3;
int QuarkC_pdgID=4;
int QuarkB_pdgID=5;
int QuarkT_pdgID=6;

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
   TH1F *h_PxTwoLQSystem = new TH1F("h_PxTwoLQSystem","PxTwoLQSystem",1001,-500.5,500.5);
   TH1F *h_PyTwoLQSystem = new TH1F("h_PyTwoLQSystem","PyTwoLQSystem",1001,-500.5,500.5);
   TH1F *h_PzTwoLQSystem = new TH1F("h_PzTwoLQSystem","PzTwoLQSystem",1001,-500.5,500.5);


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

     vector<int> v_idx_LQ;
     vector<int> v_idx_elefromLQ;
     vector<int> v_idx_quarkfromLQ;

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

     //## TLorentz vectors for LQ/ele/quark
     
     vector<TLorentzVector> TLV_LQ;
     vector<TLorentzVector> TLV_eleFromLQ;
     vector<TLorentzVector> TLV_quarkFromLQ;

     for(int particle=0; particle < v_idx_elefromLQ.size(); particle++)
       {

	 TLorentzVector TLV_LQTemp;
	 TLorentzVector TLV_eleTemp;
	 TLorentzVector TLV_quarkTemp;

	 TLV_LQTemp.SetPtEtaPhiE(GenParticlePt[v_idx_LQ[particle]],
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

	 TLV_LQ.push_back(TLV_LQTemp);
	 TLV_eleFromLQ.push_back(TLV_eleTemp);
	 TLV_quarkFromLQ.push_back(TLV_quarkTemp);

       }

     float massLQ0_eq = (TLV_eleFromLQ[0]+TLV_quarkFromLQ[0]).M();
     float massLQ1_eq = (TLV_eleFromLQ[1]+TLV_quarkFromLQ[1]).M();
     float massLQ0 = TLV_LQ[0].M();
     float massLQ1 = TLV_LQ[1].M();


     //## loop over genJets

//      for(int genjet=0; genjet<genJetCount; genjet++)
//        {
	 
// 	 TLorentzVector TLV_genJetTemp;
	 
// 	 TLV_genJetTemp.SetPtEtaPhiE(genJetPt[genjet],
// 				     genJetEta[genjet],
// 				     genJetPhi[genjet],
// 				     genJetE[genjet]);  
	 
// 	 for(int qrk=0; qrk<TLV_eleFromLQ.size(); qrk++)
// 	   {
	     
// 	     float deltaR_genJet_ele = TLV_genJetTemp.DeltaR(TLV_eleFromLQ[qrk]);
// 	     float deltaR_genJet_quark = TLV_genJetTemp.DeltaR(TLV_quarkFromLQ[qrk]);

// 	     if(deltaR_genJet_ele<0.02)
// 	       genJetMatchEle=true;

// 	     if(deltaR_genJet_quark<0.02)
// 	       genJetMatchQuark=true;
	     
// 	     break;
// 	   }
	 
//        }

     //## Fill histograms per event

     h_pthat->Fill(pthat);
     h_processID->Fill(processID);
     h_NumLQPerEvent->Fill(v_idx_LQ.size());
     h_NumEleFromLQPerEvent->Fill(v_idx_elefromLQ.size());
     h_NumQuarkFromLQPerEvent->Fill(v_idx_quarkfromLQ.size());

     h_EnergyVsMassLQ->Fill( TLV_LQ[0].E() , massLQ0 );
     h_EnergyVsMassLQ->Fill( TLV_LQ[1].E() , massLQ1 );
     h_MassLQeqMinusMassLQ->Fill( massLQ0_eq - massLQ0 );
     h_MassLQeqMinusMassLQ->Fill( massLQ1_eq - massLQ1 );

     h_PxTwoLQSystem->Fill( ( TLV_LQ[0] + TLV_LQ[1] ).Px() );
     h_PyTwoLQSystem->Fill( ( TLV_LQ[0] + TLV_LQ[1] ).Py() );
     h_PzTwoLQSystem->Fill( ( TLV_LQ[0] + TLV_LQ[1] ).Pz() );

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
   h_PxTwoLQSystem->Write();
   h_PyTwoLQSystem->Write();
   h_PzTwoLQSystem->Write();

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

}
