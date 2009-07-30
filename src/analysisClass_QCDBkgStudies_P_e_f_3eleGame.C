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

   // Histos to combine for N_events

   TH1F *h_Pt1stEle = new TH1F ("Pt1stEle","Pt1stEle",200,0,2000);  h_Pt1stEle->Sumw2();
   TH1F *h_Pt2ndEle = new TH1F ("Pt2ndEle","Pt2ndEle",200,0,2000);  h_Pt2ndEle->Sumw2();
   TH1F *h_Pt1stJet = new TH1F ("Pt1stJet","Pt1stJet",200,0,2000);  h_Pt1stJet->Sumw2();
   TH1F *h_Pt2ndJet = new TH1F ("Pt2ndJet","Pt2ndJet",200,0,2000);  h_Pt2ndJet->Sumw2();

   TH1F *h_Eta1stEle = new TH1F ("Eta1stEle","Eta1stEle",200,-3,3);  h_Eta1stEle->Sumw2();
   TH1F *h_Eta2ndEle = new TH1F ("Eta2ndEle","Eta2ndEle",200,-3,3);  h_Eta2ndEle->Sumw2();
   TH1F *h_Eta1stJet = new TH1F ("Eta1stJet","Eta1stJet",200,-5,5);  h_Eta1stJet->Sumw2();
   TH1F *h_Eta2ndJet = new TH1F ("Eta2ndJet","Eta2ndJet",200,-5,5);  h_Eta2ndJet->Sumw2();

   TH1F *h_Pt1stEle_12 = new TH1F ("Pt1stEle_12","Pt1stEle_12",200,0,2000);  h_Pt1stEle_12->Sumw2();
   TH1F *h_Pt1stEle_13 = new TH1F ("Pt1stEle_13","Pt1stEle_13",200,0,2000);  h_Pt1stEle_13->Sumw2();
   TH1F *h_Pt1stEle_23 = new TH1F ("Pt1stEle_23","Pt1stEle_23",200,0,2000);  h_Pt1stEle_23->Sumw2();

   TH1F *h_Pt1stJet_12 = new TH1F ("Pt1stJet_12","Pt1stJet_12",200,0,2000);  h_Pt1stJet_12->Sumw2();
   TH1F *h_Pt1stJet_13 = new TH1F ("Pt1stJet_13","Pt1stJet_13",200,0,2000);  h_Pt1stJet_13->Sumw2();
   TH1F *h_Pt1stJet_23 = new TH1F ("Pt1stJet_23","Pt1stJet_23",200,0,2000);  h_Pt1stJet_23->Sumw2();

   TH1F *h_Pt1stEle_barrel = new TH1F ("Pt1stEle_barrel","Pt1stEle_barrel",200,0,2000);  h_Pt1stEle_barrel->Sumw2();
   TH1F *h_Pt1stEle_12_barrel = new TH1F ("Pt1stEle_12_barrel","Pt1stEle_12_barrel",200,0,2000);  h_Pt1stEle_12_barrel->Sumw2();
   TH1F *h_Pt1stEle_13_barrel = new TH1F ("Pt1stEle_13_barrel","Pt1stEle_13_barrel",200,0,2000);  h_Pt1stEle_13_barrel->Sumw2();
   TH1F *h_Pt1stEle_23_barrel = new TH1F ("Pt1stEle_23_barrel","Pt1stEle_23_barrel",200,0,2000);  h_Pt1stEle_23_barrel->Sumw2();

   TH1F *h_Pt1stEle_endcap = new TH1F ("Pt1stEle_endcap","Pt1stEle_endcap",200,0,2000);  h_Pt1stEle_endcap->Sumw2();
   TH1F *h_Pt1stEle_12_endcap = new TH1F ("Pt1stEle_12_endcap","Pt1stEle_12_endcap",200,0,2000);  h_Pt1stEle_12_endcap->Sumw2();
   TH1F *h_Pt1stEle_13_endcap = new TH1F ("Pt1stEle_13_endcap","Pt1stEle_13_endcap",200,0,2000);  h_Pt1stEle_13_endcap->Sumw2();
   TH1F *h_Pt1stEle_23_endcap = new TH1F ("Pt1stEle_23_endcap","Pt1stEle_23_endcap",200,0,2000);  h_Pt1stEle_23_endcap->Sumw2();

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<100;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done every event

     //## HLT
     int PassTrig=0;
     int TrigBit1_1=int(getPreCutValue1("TriggerBits1"));
     int TrigBit1_2=int(getPreCutValue2("TriggerBits1"));
     int TrigBit1_3=int(getPreCutValue3("TriggerBits1"));
     int TrigBit1_4=int(getPreCutValue4("TriggerBits1"));

     if ( (HLTResults[TrigBit1_1]) || (HLTResults[TrigBit1_2]) )
	 PassTrig=1;

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

     //cout << "Electrons" << endl;

     //## Electrons
     vector<int> v_idx_ele_noCut;
     // here pre-cut on pT is applied
     vector<int> v_idx_ele;
     vector<int> v_idx_ele_ID;
     vector<int> v_idx_ele_ISO;
    vector<int> v_idx_ele_ID_ISO;
     vector<int> v_idx_ele_final;
     vector<int> v_idx_ele_final_12;
     vector<int> v_idx_ele_final_13;
     vector<int> v_idx_ele_final_23;
     int nElePtCut=0;

     for(int iele=0;iele<eleCount;iele++)
       {

	 //no cut on reco electrons
	 v_idx_ele_noCut.push_back(iele);

	 //pT pre-cut on reco electrons
	 if ( elePt[iele] < getPreCutValue1("ele_PtPreCut") ) continue;

	 v_idx_ele.push_back(iele);

	 //barrel-endcap definition
	 bool in_Barrel=false;
         bool in_Endcap=false;
	 if( fabs(eleEta[iele]) < getPreCutValue1("eleEta_bar") )  in_Barrel=true;
	 if( ( fabs(eleEta[iele]) > getPreCutValue1("eleEta_end") )
	     && 
	     ( fabs(eleEta[iele]) < getPreCutValue2("eleEta_end") ) 
	     )  in_Endcap=true;

	 //ID and ISO booleans
	 //%%%%%%% ALWAYS TRUE HERE %%%%%%%%%
	 bool pass_ISO=true;
	 bool pass_ID=true;
	 //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	 //ID cuts
	 bool pass_HoE=false;
	 if( (in_Barrel)&&((eleHoE[iele]) < getPreCutValue1("ID_HoE_bar")) )  pass_HoE=true;
	 if( (in_Endcap)&&((eleHoE[iele]) < getPreCutValue1("ID_HoE_end")) )  pass_HoE=true;
	 bool pass_sigEtaEta=false;
	 if( (in_Barrel)&&((eleSigmaEE[iele]) < getPreCutValue1("ID_sigEtaEta_bar")) )  pass_sigEtaEta=true;
	 if( (in_Endcap)&&((eleSigmaEE[iele]) < getPreCutValue1("ID_sigEtaEta_end")) )  pass_sigEtaEta=true;
	 bool pass_deltaEta=false;
	 if( (in_Barrel)&&(fabs(eleDeltaEtaTrkSC[iele]) < getPreCutValue1("ID_deltaEta_bar")) )  pass_deltaEta=true;
	 if( (in_Endcap)&&(fabs(eleDeltaEtaTrkSC[iele]) < getPreCutValue1("ID_deltaEta_end")) )  pass_deltaEta=true;
	 bool pass_deltaPhi=false;
	 if( (in_Barrel)&&(fabs(eleDeltaPhiTrkSC[iele]) < getPreCutValue1("ID_deltaPhi_bar")) )  pass_deltaPhi=true;
	 if( (in_Endcap)&&(fabs(eleDeltaPhiTrkSC[iele]) < getPreCutValue1("ID_deltaPhi_end")) )  pass_deltaPhi=true;

	 if ( (pass_HoE) && (pass_sigEtaEta) && (pass_deltaEta) && (pass_deltaPhi) ) 
	   {
	     v_idx_ele_ID.push_back(iele);
	     pass_ID=true;
	   }

	 //ISO cuts
	 bool pass_NumTrack=false;
	 if( (in_Barrel)&&((eleNumTrkIso[iele]) < getPreCutValue1("ISO_NumTrack_bar")) )  pass_NumTrack=true;
	 if( (in_Endcap)&&((eleNumTrkIso[iele]) < getPreCutValue1("ISO_NumTrack_end")) )  pass_NumTrack=true;
	 bool pass_TrackIso=false;
	 if( (in_Barrel)&&((eleTrkIso[iele]) < getPreCutValue1("ISO_TrackIso_bar")) )  pass_TrackIso=true;
	 if( (in_Endcap)&&((eleTrkIso[iele]) < getPreCutValue1("ISO_TrackIso_end")) )  pass_TrackIso=true;
	 bool pass_EcalIso=false;
	 float ecal_iso_bar = getPreCutValue1("ISO_EcalIso_bar")+(getPreCutValue2("ISO_EcalIso_bar")*elePt[iele]);
	 float ecal_iso_end = getPreCutValue1("ISO_EcalIso_end")+(getPreCutValue2("ISO_EcalIso_bar")*elePt[iele]);
	 if( (in_Barrel)&&((eleEcalRecHitIso[iele]) < (ecal_iso_bar)) )  pass_EcalIso=true;
	 if( (in_Endcap)&&((eleEcalRecHitIso[iele]) < (ecal_iso_end)) )  pass_EcalIso=true;
	 bool pass_HcalIso=false;
	 float hcal_iso_bar = getPreCutValue1("ISO_HcalIso_bar")+(getPreCutValue2("ISO_HcalIso_bar")*elePt[iele]);
	 float hcal_iso_end = getPreCutValue1("ISO_HcalIso_end")+(getPreCutValue2("ISO_HcalIso_bar")*elePt[iele]);
	 if( (in_Barrel)&&((eleHcalRecHitIso[iele]) < (hcal_iso_bar)) )  pass_HcalIso=true;
	 if( (in_Endcap)&&((eleHcalRecHitIso[iele]) < (hcal_iso_end)) )  pass_HcalIso=true;

	 if ( (pass_NumTrack) && ( pass_TrackIso) && (pass_EcalIso) && (pass_HcalIso) ) 
	   {
	     v_idx_ele_ISO.push_back(iele);
	     pass_ISO=true;
	   }

	 if ( (pass_NumTrack) && ( pass_TrackIso) && (pass_EcalIso) && (pass_HcalIso)
	      && (pass_HoE) && (pass_sigEtaEta) && (pass_deltaEta) && (pass_deltaPhi) ) 
	   {
	     v_idx_ele_ID_ISO.push_back(iele);
	   }

  
	 //%%%%%%%%% ID AND ISO   N O T   ACTUALLY APPLIED HERE %%%%%%%%%%%%
	 if ( (pass_ID) && (pass_ISO) ) v_idx_ele_final.push_back(iele);
	 
       } //loop over electrons     

     if (v_idx_ele_final.size()>2){
       v_idx_ele_final_12.push_back(v_idx_ele_final[0]);
       v_idx_ele_final_12.push_back(v_idx_ele_final[1]);
       v_idx_ele_final_13.push_back(v_idx_ele_final[0]);
       v_idx_ele_final_13.push_back(v_idx_ele_final[2]);
       v_idx_ele_final_23.push_back(v_idx_ele_final[1]);
       v_idx_ele_final_23.push_back(v_idx_ele_final[2]);
     }

     //cout << "Jets" << endl;

     //## Jets
     vector<int> v_idx_jet_noCut;
     vector<int> v_idx_jet;
     vector<int> v_idx_jet_final;
     vector<int> v_idx_jet_final_12;
     vector<int> v_idx_jet_final_23;
     vector<int> v_idx_jet_final_13;
     float deltaR_minCut = getPreCutValue1("jet_ele_DeltaRcut");

     for(int ijet=0;ijet<caloJetIC5Count;ijet++)
       {
	 //no cut on reco jets (no disambiguation)
	 v_idx_jet_noCut.push_back(ijet);

	 //pT pre-cut on reco jets (no disambiguation)
	 if ( caloJetIC5Pt[ijet] < getPreCutValue1("jet_PtPreCut") ) continue;

	 v_idx_jet.push_back(ijet);

	 //Disambiguation of electrons from jets
	 float minDeltaR=9999;
	 TVector3 jet_vec;
	 jet_vec.SetPtEtaPhi(caloJetIC5Pt[ijet],caloJetIC5Eta[ijet],caloJetIC5Phi[ijet]);

	 if (v_idx_ele_final.size()<3){
	   for (int i=0; i < v_idx_ele_final.size(); i++){
	     TVector3 ele_vec;
	     ele_vec.SetPtEtaPhi(elePt[v_idx_ele_final[i]],eleEta[v_idx_ele_final[i]],elePhi[v_idx_ele_final[i]]);
	     double distance = jet_vec.DeltaR(ele_vec);
	     if (distance<minDeltaR) minDeltaR=distance;
	   }
	 
	   if ( minDeltaR > deltaR_minCut )  v_idx_jet_final.push_back(ijet);
	 }

	 minDeltaR=9999;
	 if (v_idx_ele_final.size()>2){
	   for (int i=0; i < v_idx_ele_final_12.size(); i++){
	     TVector3 ele_vec;
	     ele_vec.SetPtEtaPhi(elePt[v_idx_ele_final_12[i]],eleEta[v_idx_ele_final_12[i]],elePhi[v_idx_ele_final_12[i]]);
	     double distance = jet_vec.DeltaR(ele_vec);
	     if (distance<minDeltaR) minDeltaR=distance;
	   }
	 
	   if ( minDeltaR > deltaR_minCut )  v_idx_jet_final_12.push_back(ijet);
	 }

	 minDeltaR=9999;
	 if (v_idx_ele_final.size()>2){
	   for (int i=0; i < v_idx_ele_final_13.size(); i++){
	     TVector3 ele_vec;
	     ele_vec.SetPtEtaPhi(elePt[v_idx_ele_final_13[i]],eleEta[v_idx_ele_final_13[i]],elePhi[v_idx_ele_final_13[i]]);
	     double distance = jet_vec.DeltaR(ele_vec);
	     if (distance<minDeltaR) minDeltaR=distance;
	   }
	 
	   if ( minDeltaR > deltaR_minCut )  v_idx_jet_final_13.push_back(ijet);
	 }

	 minDeltaR=9999;
	 if (v_idx_ele_final.size()>2){
	   for (int i=0; i < v_idx_ele_final_23.size(); i++){
	     TVector3 ele_vec;
	     ele_vec.SetPtEtaPhi(elePt[v_idx_ele_final_23[i]],eleEta[v_idx_ele_final_23[i]],elePhi[v_idx_ele_final_23[i]]);
	     double distance = jet_vec.DeltaR(ele_vec);
	     if (distance<minDeltaR) minDeltaR=distance;
	   }
	 
	   if ( minDeltaR > deltaR_minCut )  v_idx_jet_final_23.push_back(ijet);
	 }
       }

      

     // Set the evaluation of the cuts to false and clear the variable values and filled status
     resetCuts();
     
     // Set the value of the variableNames listed in the cutFile to their current value

     //cout << "HLT" << endl;

     //## HLT
     fillVariableWithValue( "HLT", PassTrig ) ;

     //cout << "nEle" << endl;

     //## nEle
     fillVariableWithValue( "nEle_noCut", v_idx_ele_noCut.size() ) ;
     fillVariableWithValue( "nEle_PtPreCut", v_idx_ele.size() ) ;

     //cout << "nJet" << endl;

     //## nJet
     fillVariableWithValue( "nJet_noCut", v_idx_jet_noCut.size() ) ;
     fillVariableWithValue( "nJet_PtPreCut", v_idx_jet.size() ) ;

     fillVariableWithValue( "nJet_PtPreCut_DIS", v_idx_jet_final.size() ) ;
     fillVariableWithValue( "nJet_PtPreCut_DIS_12", v_idx_jet_final_12.size() ) ;
     fillVariableWithValue( "nJet_PtPreCut_DIS_13", v_idx_jet_final_13.size() ) ;
     fillVariableWithValue( "nJet_PtPreCut_DIS_23", v_idx_jet_final_23.size() ) ;

     //cout << "1st Ele" << endl;

     //## 1st ele
     if( v_idx_ele_final.size() >= 1 ) 
       {
	 fillVariableWithValue( "Pt1stEle", elePt[v_idx_ele_final[0]] );
	 fillVariableWithValue( "Eta1stEle", eleEta[v_idx_ele_final[0]] );
       }


     //cout << "2nd Ele" << endl;

     //## 2nd ele
     if( v_idx_ele_final.size() >= 2 ) 
       {
	 fillVariableWithValue( "Pt2ndEle", elePt[v_idx_ele_final[1]] );
	 fillVariableWithValue( "Eta2ndEle", eleEta[v_idx_ele_final[1]] );
       }

     //cout << "1st Jet" << endl;

     //## 1st jet
     if( v_idx_jet_final.size() >= 1 ) 
       {
	 fillVariableWithValue( "Pt1stJet_DIS", caloJetIC5Pt[v_idx_jet_final[0]] );
	 fillVariableWithValue( "Eta1stJet_DIS", caloJetIC5Eta[v_idx_jet_final[0]] );
       }

     //cout << "2nd Jet" << endl;

     //## 2nd jet
     if( v_idx_jet_final.size() >= 2 ) 
       {
	 fillVariableWithValue( "Pt2ndJet_DIS", caloJetIC5Pt[v_idx_jet_final[1]] );
	 fillVariableWithValue( "Eta2ndJet_DIS", caloJetIC5Eta[v_idx_jet_final[1]] );
       }

     //## 1st jet
     if( v_idx_jet_final_12.size() >= 1 ) 
       {
	 fillVariableWithValue( "Pt1stJet_DIS_12", caloJetIC5Pt[v_idx_jet_final_12[0]] );
	 fillVariableWithValue( "Eta1stJet_DIS_12", caloJetIC5Eta[v_idx_jet_final_12[0]] );
       }

     //## 2nd jet
     if( v_idx_jet_final_12.size() >= 2 ) 
       {
	 fillVariableWithValue( "Pt2ndJet_DIS_12", caloJetIC5Pt[v_idx_jet_final_12[1]] );
	 fillVariableWithValue( "Eta2ndJet_DIS_12", caloJetIC5Eta[v_idx_jet_final_12[1]] );
       }

     //## 1st jet
     if( v_idx_jet_final_13.size() >= 1 ) 
       {
	 fillVariableWithValue( "Pt1stJet_DIS_13", caloJetIC5Pt[v_idx_jet_final_13[0]] );
	 fillVariableWithValue( "Eta1stJet_DIS_13", caloJetIC5Eta[v_idx_jet_final_13[0]] );
       }

     //## 2nd jet
     if( v_idx_jet_final_13.size() >= 2 ) 
       {
	 fillVariableWithValue( "Pt2ndJet_DIS_13", caloJetIC5Pt[v_idx_jet_final_13[1]] );
	 fillVariableWithValue( "Eta2ndJet_DIS_13", caloJetIC5Eta[v_idx_jet_final_13[1]] );
       }

     //## 1st jet
     if( v_idx_jet_final_23.size() >= 1 ) 
       {
	 fillVariableWithValue( "Pt1stJet_DIS_23", caloJetIC5Pt[v_idx_jet_final_23[0]] );
	 fillVariableWithValue( "Eta1stJet_DIS_23", caloJetIC5Eta[v_idx_jet_final_23[0]] );
       }

     //## 2nd jet
     if( v_idx_jet_final_23.size() >= 2 ) 
       {
	 fillVariableWithValue( "Pt2ndJet_DIS_23", caloJetIC5Pt[v_idx_jet_final_23[1]] );
	 fillVariableWithValue( "Eta2ndJet_DIS_23", caloJetIC5Eta[v_idx_jet_final_23[1]] );
       }


     //## define "2ele" and "2jets" booleans
     bool TwoEles=false;
     bool TwoJets=false;
     bool TwoJets_12=false;
     bool TwoJets_13=false;
     bool TwoJets_23=false;
     if( v_idx_ele_final.size() >= 2 ) TwoEles = true;
     if( v_idx_jet_final.size() >= 2 ) TwoJets = true;
     if( v_idx_jet_final_12.size() >= 2 ) TwoJets_12 = true;
     if( v_idx_jet_final_13.size() >= 2 ) TwoJets_13 = true;
     if( v_idx_jet_final_23.size() >= 2 ) TwoJets_23 = true;

     //cout << "Mee" << endl;

     //## Mee
     if( TwoEles ) 
       {
	 TLorentzVector v_ee, ele1, ele2;
	 ele1.SetPtEtaPhiM(elePt[v_idx_ele_final[0]],
			   eleEta[v_idx_ele_final[0]],
			   elePhi[v_idx_ele_final[0]],0);
	 ele2.SetPtEtaPhiM(elePt[v_idx_ele_final[1]],
			   eleEta[v_idx_ele_final[1]],
			   elePhi[v_idx_ele_final[1]],0);
	 v_ee = ele1 + ele2;
	 fillVariableWithValue( "invMass_ee", v_ee.M() ) ;
       }

     if( v_idx_ele_final_12.size()>1 ) 
       {
	 TLorentzVector v_ee, ele1, ele2;
	 ele1.SetPtEtaPhiM(elePt[v_idx_ele_final_12[0]],
			   eleEta[v_idx_ele_final_12[0]],
			   elePhi[v_idx_ele_final_12[0]],0);
	 ele2.SetPtEtaPhiM(elePt[v_idx_ele_final_12[1]],
			   eleEta[v_idx_ele_final_12[1]],
			   elePhi[v_idx_ele_final_12[1]],0);
	 v_ee = ele1 + ele2;
	 fillVariableWithValue( "invMass_ee_12", v_ee.M() ) ;
       }

     if( v_idx_ele_final_13.size()>1 ) 
       {
	 TLorentzVector v_ee, ele1, ele2;
	 ele1.SetPtEtaPhiM(elePt[v_idx_ele_final_13[0]],
			   eleEta[v_idx_ele_final_13[0]],
			   elePhi[v_idx_ele_final_13[0]],0);
	 ele2.SetPtEtaPhiM(elePt[v_idx_ele_final_13[1]],
			   eleEta[v_idx_ele_final_13[1]],
			   elePhi[v_idx_ele_final_13[1]],0);
	 v_ee = ele1 + ele2;
	 fillVariableWithValue( "invMass_ee_13", v_ee.M() ) ;
       }

     if( v_idx_ele_final_23.size()>1 ) 
       {
	 TLorentzVector v_ee, ele1, ele2;
	 ele1.SetPtEtaPhiM(elePt[v_idx_ele_final_23[0]],
			   eleEta[v_idx_ele_final_23[0]],
			   elePhi[v_idx_ele_final_23[0]],0);
	 ele2.SetPtEtaPhiM(elePt[v_idx_ele_final_23[1]],
			   eleEta[v_idx_ele_final_23[1]],
			   elePhi[v_idx_ele_final_23[1]],0);
	 v_ee = ele1 + ele2;
	 fillVariableWithValue( "invMass_ee_23", v_ee.M() ) ;
       }

     //cout << "ST" << endl;

     //## ST
     if ( (TwoEles) && (TwoJets) ) 
       {
	 double calc_sT = 
	   elePt[v_idx_ele_final[0]]
	   + elePt[v_idx_ele_final[1]]
	   + caloJetIC5Pt[v_idx_jet_final[0]]
	   + caloJetIC5Pt[v_idx_jet_final[1]];
	 fillVariableWithValue("sT", calc_sT);
       }

     if ( (TwoEles) && (TwoJets_12) ) 
       {
	 double calc_sT = 
	   elePt[v_idx_ele_final_12[0]]
	   + elePt[v_idx_ele_final_12[1]]
	   + caloJetIC5Pt[v_idx_jet_final_12[0]]
	   + caloJetIC5Pt[v_idx_jet_final_12[1]];
	 fillVariableWithValue("sT_12", calc_sT);
       }

     if ( (TwoEles) && (TwoJets_13) ) 
       {
	 double calc_sT = 
	   elePt[v_idx_ele_final_13[0]]
	   + elePt[v_idx_ele_final_13[1]]
	   + caloJetIC5Pt[v_idx_jet_final_13[0]]
	   + caloJetIC5Pt[v_idx_jet_final_13[1]];
	 fillVariableWithValue("sT_13", calc_sT);
       }

     if ( (TwoEles) && (TwoJets_23) ) 
       {
	 double calc_sT = 
	   elePt[v_idx_ele_final_23[0]]
	   + elePt[v_idx_ele_final_23[1]]
	   + caloJetIC5Pt[v_idx_jet_final_23[0]]
	   + caloJetIC5Pt[v_idx_jet_final_23[1]];
	 fillVariableWithValue("sT_23", calc_sT);
       }


     //--------------------------------------------------

     // Evaluate cuts (but do not apply them)
     evaluateCuts();

     // Fill histograms and do analysis based on cut evaluation
     //double probBarrel = getPreCutValue1("probBarrel");
     double probEndcap = 0.10605;  // average value for endcap electrons 

     vector< pair<int,double> > probBarrel; //first=lower bin edge, second=fake rate
     probBarrel.push_back(make_pair(0,0.05));
     probBarrel.push_back(make_pair(50,0.052));
     probBarrel.push_back(make_pair(100,0.065));
     probBarrel.push_back(make_pair(150,0.092));
     probBarrel.push_back(make_pair(200,0.105));
     probBarrel.push_back(make_pair(250,0.12));
     probBarrel.push_back(make_pair(300,0.125));
     probBarrel.push_back(make_pair(350,0.135));
     probBarrel.push_back(make_pair(400,0.155));
     probBarrel.push_back(make_pair(450,0.120));
     probBarrel.push_back(make_pair(550,0.085));
     probBarrel.push_back(make_pair(700,0.0));

     double p1=0;
     double p2=0;


     if( (passedCut("0")) && (passedCut("1")) && (v_idx_ele_final.size()<3) )
       {
	 if( fabs(eleEta[v_idx_ele_final[0]]) < getPreCutValue1("eleEta_bar") ) 
	   { for (int i=0;i<probBarrel.size();i++)
	       if (elePt[v_idx_ele_final[0]]> probBarrel[i].first) p1=probBarrel[i].second;
	   }
	 if( ( fabs(eleEta[v_idx_ele_final[0]]) > getPreCutValue1("eleEta_end") )
	     && 
	     ( fabs(eleEta[v_idx_ele_final[0]]) < getPreCutValue2("eleEta_end") ) 
	     )  p1=probEndcap;
	 if( fabs(eleEta[v_idx_ele_final[1]]) < getPreCutValue1("eleEta_bar") )
	   { for (int i=0;i<probBarrel.size();i++)
	       if (elePt[v_idx_ele_final[1]]> probBarrel[i].first) p2=probBarrel[i].second;
	   }
	 if( ( fabs(eleEta[v_idx_ele_final[1]]) > getPreCutValue1("eleEta_end") )
	     && 
	     ( fabs(eleEta[v_idx_ele_final[1]]) < getPreCutValue2("eleEta_end") ) 
	     )  p2=probEndcap;

	 h_Pt1stEle->Fill(elePt[v_idx_ele_final[0]],p1*p2);
	 h_Pt1stJet->Fill(caloJetIC5Pt[v_idx_jet_final[0]],p1*p2);
	 h_Pt2ndEle->Fill(elePt[v_idx_ele_final[1]],p1*p2);
	 h_Pt2ndJet->Fill(caloJetIC5Pt[v_idx_jet_final[1]],p1*p2);

	 h_Eta1stEle->Fill(eleEta[v_idx_ele_final[0]],p1*p2);
	 h_Eta1stJet->Fill(caloJetIC5Eta[v_idx_jet_final[0]],p1*p2);
	 h_Eta2ndEle->Fill(eleEta[v_idx_ele_final[1]],p1*p2);
	 h_Eta2ndJet->Fill(caloJetIC5Eta[v_idx_jet_final[1]],p1*p2);

       }
     
     p1=0;
     p2=0;
     if( (passedCut("0")) && (passedCut("2")) )
       {
	 h_Pt1stEle_12->Fill(elePt[v_idx_ele_final_12[0]]);
	 h_Pt1stJet_12->Fill(caloJetIC5Pt[v_idx_jet_final_12[0]]);

	 if( fabs(eleEta[v_idx_ele_final_12[0]]) < getPreCutValue1("eleEta_bar") )
	   { for (int i=0;i<probBarrel.size();i++)
	       if (elePt[v_idx_ele_final_12[0]]> probBarrel[i].first) p1=probBarrel[i].second;
	   }
	 if( ( fabs(eleEta[v_idx_ele_final_12[0]]) > getPreCutValue1("eleEta_end") )
	     && 
	     ( fabs(eleEta[v_idx_ele_final_12[0]]) < getPreCutValue2("eleEta_end") ) 
	     )  p1=probBarrel[0].second;
	 if( fabs(eleEta[v_idx_ele_final_12[1]]) < getPreCutValue1("eleEta_bar") )
	   { for (int i=0;i<probBarrel.size();i++)
	       if (elePt[v_idx_ele_final_12[1]]> probBarrel[i].first) p2=probBarrel[i].second;
	   }
	 if( ( fabs(eleEta[v_idx_ele_final_12[1]]) > getPreCutValue1("eleEta_end") )
	     && 
	     ( fabs(eleEta[v_idx_ele_final_12[1]]) < getPreCutValue2("eleEta_end") ) 
	     )  p2=probBarrel[0].second;

	 h_Pt1stEle->Fill(elePt[v_idx_ele_final_12[0]],p1*p2);
	 h_Pt1stJet->Fill(caloJetIC5Pt[v_idx_jet_final_12[0]],p1*p2);
	 h_Pt2ndEle->Fill(elePt[v_idx_ele_final_12[1]],p1*p2);
	 h_Pt2ndJet->Fill(caloJetIC5Pt[v_idx_jet_final_12[1]],p1*p2);

	 h_Eta1stEle->Fill(eleEta[v_idx_ele_final_12[0]],p1*p2);
	 h_Eta1stJet->Fill(caloJetIC5Eta[v_idx_jet_final_12[0]],p1*p2);
	 h_Eta2ndEle->Fill(eleEta[v_idx_ele_final_12[1]],p1*p2);
	 h_Eta2ndJet->Fill(caloJetIC5Eta[v_idx_jet_final_12[1]],p1*p2);

       }
     
     p1=0;
     p2=0;
     if( (passedCut("0")) && (passedCut("3")) )
       {
	 h_Pt1stEle_13->Fill(elePt[v_idx_ele_final_13[0]]);
	 h_Pt1stJet_13->Fill(caloJetIC5Pt[v_idx_jet_final_13[0]]);

	 if( fabs(eleEta[v_idx_ele_final_13[0]]) < getPreCutValue1("eleEta_bar") ) 
	   { for (int i=0;i<probBarrel.size();i++)
	       if (elePt[v_idx_ele_final_13[0]]> probBarrel[i].first) p1=probBarrel[i].second;
	   }
	 if( ( fabs(eleEta[v_idx_ele_final_13[0]]) > getPreCutValue1("eleEta_end") )
	     && 
	     ( fabs(eleEta[v_idx_ele_final_13[0]]) < getPreCutValue2("eleEta_end") ) 
	     )  p1=probEndcap;
	 if( fabs(eleEta[v_idx_ele_final_13[1]]) < getPreCutValue1("eleEta_bar") )
	   { for (int i=0;i<probBarrel.size();i++)
	       if (elePt[v_idx_ele_final_13[1]]> probBarrel[i].first) p2=probBarrel[i].second;
	   }
	 if( ( fabs(eleEta[v_idx_ele_final_13[1]]) > getPreCutValue1("eleEta_end") )
	     && 
	     ( fabs(eleEta[v_idx_ele_final_13[1]]) < getPreCutValue2("eleEta_end") ) 
	     )  p2=probEndcap;

	 h_Pt1stEle->Fill(elePt[v_idx_ele_final_13[0]],p1*p2);
	 h_Pt1stJet->Fill(caloJetIC5Pt[v_idx_jet_final_13[0]],p1*p2);
	 h_Pt2ndEle->Fill(elePt[v_idx_ele_final_13[1]],p1*p2);
	 h_Pt2ndJet->Fill(caloJetIC5Pt[v_idx_jet_final_13[1]],p1*p2);

	 h_Eta1stEle->Fill(eleEta[v_idx_ele_final_13[0]],p1*p2);
	 h_Eta1stJet->Fill(caloJetIC5Eta[v_idx_jet_final_13[0]],p1*p2);
	 h_Eta2ndEle->Fill(eleEta[v_idx_ele_final_13[1]],p1*p2);
	 h_Eta2ndJet->Fill(caloJetIC5Eta[v_idx_jet_final_13[1]],p1*p2);

       }
     
     p1=0;
     p2=0;
     if( (passedCut("0")) && (passedCut("4")) )
       {
	 h_Pt1stEle_23->Fill(elePt[v_idx_ele_final_23[0]]);
	 h_Pt1stJet_23->Fill(caloJetIC5Pt[v_idx_jet_final_23[0]]);
 
	 if( fabs(eleEta[v_idx_ele_final_23[0]]) < getPreCutValue1("eleEta_bar") )
	   { for (int i=0;i<probBarrel.size();i++)
	       if (elePt[v_idx_ele_final_23[0]]> probBarrel[i].first) p1=probBarrel[i].second;
	   }
	 if( ( fabs(eleEta[v_idx_ele_final_23[0]]) > getPreCutValue1("eleEta_end") )
	     && 
	     ( fabs(eleEta[v_idx_ele_final_23[0]]) < getPreCutValue2("eleEta_end") ) 
	     )  p1=probEndcap;
	 if( fabs(eleEta[v_idx_ele_final_23[1]]) < getPreCutValue1("eleEta_bar") )
	   { for (int i=0;i<probBarrel.size();i++)
	       if (elePt[v_idx_ele_final_23[1]]> probBarrel[i].first) p2=probBarrel[i].second;
	   }
	 if( ( fabs(eleEta[v_idx_ele_final_23[1]]) > getPreCutValue1("eleEta_end") )
	     && 
	     ( fabs(eleEta[v_idx_ele_final_23[1]]) < getPreCutValue2("eleEta_end") ) 
	     )  p2=probEndcap;

	 h_Pt1stEle->Fill(elePt[v_idx_ele_final_23[0]],p1*p2);
	 h_Pt1stJet->Fill(caloJetIC5Pt[v_idx_jet_final_23[0]],p1*p2);
	 h_Pt2ndEle->Fill(elePt[v_idx_ele_final_23[1]],p1*p2);
	 h_Pt2ndJet->Fill(caloJetIC5Pt[v_idx_jet_final_23[1]],p1*p2);

	 h_Eta1stEle->Fill(eleEta[v_idx_ele_final_23[0]],p1*p2);
	 h_Eta1stJet->Fill(caloJetIC5Eta[v_idx_jet_final_23[0]],p1*p2);
	 h_Eta2ndEle->Fill(eleEta[v_idx_ele_final_23[1]],p1*p2);
	 h_Eta2ndJet->Fill(caloJetIC5Eta[v_idx_jet_final_23[1]],p1*p2);

      }
     


     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 

   h_Pt1stEle->Write();
   h_Pt2ndEle->Write();
   h_Pt1stJet->Write();
   h_Pt2ndJet->Write();

   h_Eta1stEle->Write();
   h_Eta2ndEle->Write();
   h_Eta1stJet->Write();
   h_Eta2ndJet->Write();

   h_Pt1stEle_12->Write();
   h_Pt1stEle_13->Write();
   h_Pt1stEle_23->Write();

   h_Pt1stJet_12->Write();
   h_Pt1stJet_13->Write();
   h_Pt1stJet_23->Write();

   h_Pt1stEle_barrel->Write();
   h_Pt1stEle_12_barrel->Write();
   h_Pt1stEle_13_barrel->Write();
   h_Pt1stEle_23_barrel->Write();

   h_Pt1stEle_endcap->Write();
   h_Pt1stEle_12_endcap->Write();
   h_Pt1stEle_13_endcap->Write();
   h_Pt1stEle_23_endcap->Write();


   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
