#define analysisClass_cxx
#include "analysisClass.h"

#include "TH1F.h"

int LQ_pdgID=42;

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
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%100 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     int NumLQPerEvent=0;


     for(int genPart=0; genPart<GenParticleCount; genPart++)
       {

	 bool IsLQ=false;

	 h_GenParticlePdgIdAll->Fill(GenParticlePdgId[genPart]);

	 if(abs(GenParticlePdgId[genPart])==LQ_pdgID)
	   {
	     IsLQ=true;
	   }

	 if(IsLQ)
	   {
	     NumLQPerEvent++;
	   }
       }


     //## Fill histograms per event

     h_pthat->Fill(pthat);
     h_processID->Fill(processID);
     h_NumLQPerEvent->Fill(NumLQPerEvent);

     ////////////////////// User's code ends here ///////////////////////
   }
   std::cout << "analysisClass::Loop() ends" <<std::endl;   



   //## Write histograms

   h_pthat->Write();
   h_processID->Write();
   h_GenParticlePdgIdAll->Write();
   h_NumLQPerEvent->Write();

}
