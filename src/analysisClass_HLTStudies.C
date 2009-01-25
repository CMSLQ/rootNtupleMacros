#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <vector>
#include <string.h>

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
    
   if (fChain == 0) return;
   
   //////////book histos here
   TH1F *h_ele_E_80 = new TH1F ("ele_E_80","ele_E_80",100,0,1000);
   TH1F *h_ele_Pt_80  = new TH1F ("ele_Pt_80","ele_Pt_80",100,0,1000);
   TH1F *h_ele_Phi_80  = new TH1F ("ele_Phi_80","ele_Phi_80",71,-3.55,3.55);
   TH1F *h_ele_Eta_80  = new TH1F ("ele_Eta_80","ele_Eta_80",121,-6.05,6.05);
   TH1F *h_ele_CaloEnergy_80  = new TH1F ("ele_CaloEnergy_80","ele_CaloEnergy_80",100,0,1000);

   TH1F *h_ele_E_200 = new TH1F ("ele_E_200","ele_E_200",100,0,1000);
   TH1F *h_ele_Pt_200  = new TH1F ("ele_Pt_200","ele_Pt_200",100,0,1000);
   TH1F *h_ele_Phi_200  = new TH1F ("ele_Phi_200","ele_Phi_200",71,-3.55,3.55);
   TH1F *h_ele_Eta_200  = new TH1F ("ele_Eta_200","ele_Eta_200",121,-6.05,6.05);
   TH1F *h_ele_CaloEnergy_200  = new TH1F ("ele_CaloEnergy_200","ele_CaloEnergy_200",100,0,1000);

   TH1F *h_ele_E_Photon = new TH1F ("ele_E_Photon","ele_E_Photon",100,0,1000);
   TH1F *h_ele_Pt_Photon  = new TH1F ("ele_Pt_Photon","ele_Pt_Photon",100,0,1000);
   TH1F *h_ele_Phi_Photon  = new TH1F ("ele_Phi_Photon","ele_Phi_Photon",71,-3.55,3.55);
   TH1F *h_ele_Eta_Photon  = new TH1F ("ele_Eta_Photon","ele_Eta_Photon",121,-6.05,6.05);
   TH1F *h_ele_CaloEnergy_Photon  = new TH1F ("ele_CaloEnergy_Photon","ele_CaloEnergy_Photon",100,0,1000);


   /////////initialize variables
   double EventsPassed_80=0;
   double EventsPassed_200=0;
   double EventsPassed_80_OR_200=0;
   double EventsPassed_Photon=0;
   double n_Events=0;

   int electron_PID=11;
   int LQ_PID=42;
   //int LQ_PID=23; //PdgID for Z
   float elePt_cut=30.;
   float eleEta_cut=2.6;
 
   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<1;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     string Trig1="HLT_EM80";
     string Trig2="HLT_EM200";
     string Trig3="HLT_Photon25"; 

     string tmp="";
     int results_index=0;
     for (int itrig=0;itrig<hltNamesLen;itrig++){
     //for (int itrig=0;itrig<100;itrig++){
       //cout << HLTNames[itrig];
       if (HLTNames[itrig]==':') {
	 //cout << tmp << "   " << HLTResults[results_index] << endl;
	 bool IsEM80 =  !strncmp(Trig1.c_str(),tmp.c_str(),8); // returns zero if first 8 characters in 2 strings are same
	 bool IsEM200 = !strncmp(Trig2.c_str(),tmp.c_str(),9);
	 bool IsPhoton = !strncmp(Trig3.c_str(),tmp.c_str(),11);
	 //if ((IsEM80)||(IsEM200)) cout << "Trigger:  " << tmp << "  HLTResults " << HLTResults[results_index] << endl;
	 //if (IsPhoton) cout << "Trigger:  " << tmp << "  HLTResults " << HLTResults[results_index] << endl;

	 if ((IsEM80)&&(HLTResults[results_index]==true)) {
	   EventsPassed_80++;
	   for (int iele=0; iele<eleCount; iele++) {
	     if (elePt[iele]<elePt_cut) continue;
	     h_ele_E_80->Fill(eleEnergy[iele]);
	     h_ele_Pt_80->Fill(elePt[iele]);
	     h_ele_Phi_80->Fill(elePhi[iele]);
	     h_ele_Eta_80->Fill(eleEta[iele]);
	     h_ele_CaloEnergy_80->Fill(eleCaloEnergy[iele]);
	   } //for iele	   
	 }


	 if ((IsEM200)&&(HLTResults[results_index]==true)) {
	   EventsPassed_200++;
	   for (int iele=0; iele<eleCount; iele++) {
	     if (elePt[iele]<elePt_cut) continue;
	     h_ele_E_200->Fill(eleEnergy[iele]);
	     h_ele_Pt_200->Fill(elePt[iele]);
	     h_ele_Phi_200->Fill(elePhi[iele]);
	     h_ele_Eta_200->Fill(eleEta[iele]);
	     h_ele_CaloEnergy_200->Fill(eleCaloEnergy[iele]);
	   } //for iele	   
	 }

	 if ((IsPhoton)&&(HLTResults[results_index]==true)) {
	   EventsPassed_Photon++;
	   for (int iele=0; iele<eleCount; iele++) {
	     if (elePt[iele]<elePt_cut) continue;
	     h_ele_E_Photon->Fill(eleEnergy[iele]);
	     h_ele_Pt_Photon->Fill(elePt[iele]);
	     h_ele_Phi_Photon->Fill(elePhi[iele]);
	     h_ele_Eta_Photon->Fill(eleEta[iele]);
	     h_ele_CaloEnergy_Photon->Fill(eleCaloEnergy[iele]);
	   } //for iele	   
	 }

	 if (((IsEM80)&&(HLTResults[results_index]==true))||((IsEM80)&&(HLTResults[results_index]==true))) {
	   EventsPassed_80_OR_200++;
	 }

	 tmp.clear(); //reset temporary string of HLT name
	 results_index++; //move to next HLT result
       }
       else tmp.push_back(HLTNames[itrig]) ;  //build up HLT name until ":" is reached, indicating end of name
     } //end of loop through trigger names

     n_Events++;
     ////////////////////// User's code ends here ///////////////////////
   }   

   //////////write histos 
   h_ele_E_80->Write();
   h_ele_Pt_80->Write();
   h_ele_Phi_80->Write();
   h_ele_Eta_80->Write();
   h_ele_CaloEnergy_80->Write();

   h_ele_E_200->Write();
   h_ele_Pt_200->Write();
   h_ele_Phi_200->Write();
   h_ele_Eta_200->Write();
   h_ele_CaloEnergy_200->Write();

   h_ele_E_Photon->Write();
   h_ele_Pt_Photon->Write();
   h_ele_Phi_Photon->Write();
   h_ele_Eta_Photon->Write();
   h_ele_CaloEnergy_Photon->Write();


   cout << "Events Passed Trigger EM_80: " << EventsPassed_80 << "  Fraction: " << EventsPassed_80/n_Events << endl;
   cout << "Events Passed Trigger EM_200: " << EventsPassed_200 << "  Fraction: " << EventsPassed_200/n_Events <<endl;
   cout << "Events Passed Trigger EM_80 OR EM_200: " << EventsPassed_80_OR_200 << "  Fraction: " << EventsPassed_80_OR_200/n_Events << endl;
   cout << "Events Passed Trigger Photon: " << EventsPassed_Photon << "  Fraction: " << EventsPassed_Photon/n_Events <<endl;

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
