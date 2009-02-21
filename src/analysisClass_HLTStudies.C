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

   TH1F *h_ele_E_bit1 = new TH1F ("ele_E_bit1","ele_E_bit1",100,0,1000);
   TH1F *h_ele_Pt_bit1  = new TH1F ("ele_Pt_bit1","ele_Pt_bit1",100,0,1000);
   TH1F *h_ele_Phi_bit1  = new TH1F ("ele_Phi_bit1","ele_Phi_bit1",71,-3.55,3.55);
   TH1F *h_ele_Eta_bit1  = new TH1F ("ele_Eta_bit1","ele_Eta_bit1",121,-6.05,6.05);
   TH1F *h_ele_CaloEnergy_bit1  = new TH1F ("ele_CaloEnergy_bit1","ele_CaloEnergy_bit1",100,0,1000);

   TH1F *h_ele_E_bit2 = new TH1F ("ele_E_bit2","ele_E_bit2",100,0,1000);
   TH1F *h_ele_Pt_bit2  = new TH1F ("ele_Pt_bit2","ele_Pt_bit2",100,0,1000);
   TH1F *h_ele_Phi_bit2  = new TH1F ("ele_Phi_bit2","ele_Phi_bit2",71,-3.55,3.55);
   TH1F *h_ele_Eta_bit2  = new TH1F ("ele_Eta_bit2","ele_Eta_bit2",121,-6.05,6.05);
   TH1F *h_ele_CaloEnergy_bit2  = new TH1F ("ele_CaloEnergy_bit2","ele_CaloEnergy_bit2",100,0,1000);

   TH1F *h_ele_E_bit3 = new TH1F ("ele_E_bit3","ele_E_bit3",100,0,1000);
   TH1F *h_ele_Pt_bit3  = new TH1F ("ele_Pt_bit3","ele_Pt_bit3",100,0,1000);
   TH1F *h_ele_Phi_bit3  = new TH1F ("ele_Phi_bit3","ele_Phi_bit3",71,-3.55,3.55);
   TH1F *h_ele_Eta_bit3  = new TH1F ("ele_Eta_bit3","ele_Eta_bit3",121,-6.05,6.05);
   TH1F *h_ele_CaloEnergy_bit3  = new TH1F ("ele_CaloEnergy_bit3","ele_CaloEnergy_bit3",100,0,1000);

   TH1F *h_ele_E_bit4 = new TH1F ("ele_E_bit4","ele_E_bit4",100,0,1000);
   TH1F *h_ele_Pt_bit4  = new TH1F ("ele_Pt_bit4","ele_Pt_bit4",100,0,1000);
   TH1F *h_ele_Phi_bit4  = new TH1F ("ele_Phi_bit4","ele_Phi_bit4",71,-3.55,3.55);
   TH1F *h_ele_Eta_bit4  = new TH1F ("ele_Eta_bit4","ele_Eta_bit4",121,-6.05,6.05);
   TH1F *h_ele_CaloEnergy_bit4  = new TH1F ("ele_CaloEnergy_bit4","ele_CaloEnergy_bit4",100,0,1000);

   /////////initialize variables
   double EventsPassed_bit1=0;
   double EventsPassed_bit2=0;
   double EventsPassed_bit3=0;
   double EventsPassed_bit4=0;
   double n_Events=0;

   int electron_PID=11;
   int LQ_PID=42;
   //int LQ_PID=23; //PdgID for Z
   float elePt_cut=10.;
   float eleEta_cut=2.6;

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

     //HLT
     int TrigBit1=getPreCutValue1("TriggerBits1");
     int TrigBit2=getPreCutValue1("TriggerBits2");
     int TrigBit3=getPreCutValue1("TriggerBits3");
     int TrigBit4=getPreCutValue1("TriggerBits4");

     /////To print out list of trigger names:
//      int results_index=0;
//      string tmp="";
//      for (int itrig=0;itrig<hltNamesLen;itrig++){
//        if (HLTNames[itrig]==':') {
// 	 cout << tmp << "   " << HLTResults[results_index] << endl;
// 	 tmp.clear(); //reset temporary string of HLT name
// 	 results_index++; //move to next HLT result
//        }
//        else tmp.push_back(HLTNames[itrig]) ;  //build up HLT name until ":" is reached, indicating end of name
//      }

	 if (HLTResults[TrigBit1]) {
	   EventsPassed_bit1++;
	   for (int iele=0; iele<eleCount; iele++) {
	     if (elePt[iele]<elePt_cut) continue;
	     h_ele_E_bit1->Fill(eleEnergy[iele]);
	     h_ele_Pt_bit1->Fill(elePt[iele]);
	     h_ele_Phi_bit1->Fill(elePhi[iele]);
	     h_ele_Eta_bit1->Fill(eleEta[iele]);
	     h_ele_CaloEnergy_bit1->Fill(eleCaloEnergy[iele]);
	   } //for iele	   
	 }


	 if (HLTResults[TrigBit2]) {
	   EventsPassed_bit2++;
	   for (int iele=0; iele<eleCount; iele++) {
	     if (elePt[iele]<elePt_cut) continue;
	     h_ele_E_bit2->Fill(eleEnergy[iele]);
	     h_ele_Pt_bit2->Fill(elePt[iele]);
	     h_ele_Phi_bit2->Fill(elePhi[iele]);
	     h_ele_Eta_bit2->Fill(eleEta[iele]);
	     h_ele_CaloEnergy_bit2->Fill(eleCaloEnergy[iele]);
	   } //for iele	   
	 }

	 if (HLTResults[TrigBit3]) {
	   EventsPassed_bit3++;
	   for (int iele=0; iele<eleCount; iele++) {
	     if (elePt[iele]<elePt_cut) continue;
	     h_ele_E_bit3->Fill(eleEnergy[iele]);
	     h_ele_Pt_bit3->Fill(elePt[iele]);
	     h_ele_Phi_bit3->Fill(elePhi[iele]);
	     h_ele_Eta_bit3->Fill(eleEta[iele]);
	     h_ele_CaloEnergy_bit3->Fill(eleCaloEnergy[iele]);
	   } //for iele	   
	 }

	 if (HLTResults[TrigBit4]) {
	   EventsPassed_bit4++;
	   for (int iele=0; iele<eleCount; iele++) {
	     if (elePt[iele]<elePt_cut) continue;
	     h_ele_E_bit4->Fill(eleEnergy[iele]);
	     h_ele_Pt_bit4->Fill(elePt[iele]);
	     h_ele_Phi_bit4->Fill(elePhi[iele]);
	     h_ele_Eta_bit4->Fill(eleEta[iele]);
	     h_ele_CaloEnergy_bit4->Fill(eleCaloEnergy[iele]);
	   } //for iele	   
	 }

     resetCuts();
     if (HLTResults[TrigBit1]) fillVariableWithValue("HLT_bit1",1) ;
     if (HLTResults[TrigBit2]) fillVariableWithValue("HLT_bit2",1) ;
     if (HLTResults[TrigBit3]) fillVariableWithValue("HLT_bit3",1) ;
     if (HLTResults[TrigBit4]) fillVariableWithValue("HLT_bit4",1) ;
     evaluateCuts();

    n_Events++;

   } // End loop over events

   //////////write histos 

   h_ele_E_bit1->Write();
   h_ele_Pt_bit1->Write();
   h_ele_Phi_bit1->Write();
   h_ele_Eta_bit1->Write();
   h_ele_CaloEnergy_bit1->Write();

   h_ele_E_bit2->Write();
   h_ele_Pt_bit2->Write();
   h_ele_Phi_bit2->Write();
   h_ele_Eta_bit2->Write();
   h_ele_CaloEnergy_bit2->Write();

   h_ele_E_bit3->Write();
   h_ele_Pt_bit3->Write();
   h_ele_Phi_bit3->Write();
   h_ele_Eta_bit3->Write();
   h_ele_CaloEnergy_bit3->Write();

   h_ele_E_bit4->Write();
   h_ele_Pt_bit4->Write();
   h_ele_Phi_bit4->Write();
   h_ele_Eta_bit4->Write();
   h_ele_CaloEnergy_bit4->Write();

   cout << "Events Passed Trigger EM_bit1: " << EventsPassed_bit1 << "  Fraction: " << EventsPassed_bit1/n_Events << endl;
   cout << "Events Passed Trigger EM_bit2: " << EventsPassed_bit2 << "  Fraction: " << EventsPassed_bit2/n_Events <<endl;
   cout << "Events Passed Trigger EM_bit3: " << EventsPassed_bit3 << "  Fraction: " << EventsPassed_bit3/n_Events << endl;
   cout << "Events Passed Trigger EM_bit4: " << EventsPassed_bit4 << "  Fraction: " << EventsPassed_bit4/n_Events <<endl;

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
