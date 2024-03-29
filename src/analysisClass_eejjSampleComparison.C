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

   //Gen Particles
   TH1F *h_GenPart_eleEta = new TH1F ("GenPart_eleEta","GenPart_eleEta",101,-5.05,5.05); h_GenPart_eleEta->Sumw2();
   TH1F *h_GenPart_elePt = new TH1F ("GenPart_elePt","GenPart_elePt",100,0,1000); h_GenPart_elePt->Sumw2();
   TH1F *h_GenPart_quarkEta = new TH1F ("GenPart_quarkEta","GenPart_quarkEta",101,-5.05,5.05); h_GenPart_quarkEta->Sumw2();
   TH1F *h_GenPart_quarkPt = new TH1F ("GenPart_quarkPt","GenPart_quarkPt",100,0,1000); h_GenPart_quarkPt->Sumw2();
   TH1F *h_GenPart_LQEta = new TH1F ("GenPart_LQEta","GenPart_LQEta",101,-5.05,5.05); h_GenPart_LQEta->Sumw2();
   TH1F *h_GenPart_LQPt = new TH1F ("GenPart_LQPt","GenPart_LQPt",100,0,1000); h_GenPart_LQPt->Sumw2();
   TH1F *h_GenJet_pT_1st = new TH1F ("GenJet_pT_1st","GenJet_pT_1st",100,0,1000); h_GenJet_pT_1st->Sumw2();
   TH1F *h_GenJet_eta_1st = new TH1F ("GenJet_eta_1st","GenJet_eta_1st",101,-5.05,5.05); h_GenJet_eta_1st->Sumw2();
   TH1F *h_GenJet_pT_2nd = new TH1F ("GenJet_pT_2nd","GenJet_pT_2nd",100,0,1000); h_GenJet_pT_2nd->Sumw2();
   TH1F *h_GenJet_eta_2nd = new TH1F ("GenJet_eta_2nd","GenJet_eta_2nd",101,-5.05,5.05); h_GenJet_eta_2nd->Sumw2();

   //electrons
   TH1F *h_elec_sigEtaEta_barrel = new TH1F ("elec_sigEtaEta_barrel","Sigma Eta Eta (ID var) Barrel",60,0,0.06); h_elec_sigEtaEta_barrel->Sumw2();
   TH1F *h_elec_sigEtaEta_endcap = new TH1F ("elec_sigEtaEta_endcap","Sigma Eta Eta (ID var) Endcap",60,0,0.06); h_elec_sigEtaEta_endcap->Sumw2();
   TH1F *h_elec_HoE_barrel = new TH1F ("elec_HoE_barrel","H over E (ID var) Barrel",200,0,0.2); h_elec_HoE_barrel->Sumw2();
   TH1F *h_elec_HoE_endcap = new TH1F ("elec_HoE_endcap","H over E (ID var) Endcap",200,0,0.2); h_elec_HoE_endcap->Sumw2();
   TH1F *h_nEle_ptCut = new TH1F ("nEle_ptCut","Number elecs passing pt cut",6,-0.5,5.5); h_nEle_ptCut->Sumw2();
   TH1F *h_nEleID_30GeV = new TH1F ("nEleID_30GeV","Number of IDed elecs",11,-0.5,10.5); h_nEleID_30GeV->Sumw2();
   TH1F *h_nEleISO_30GeV = new TH1F ("nEleISO_30GeV","Number of Isolated elecs",11,-0.5,10.5); h_nEleISO_30GeV->Sumw2();
   TH1F *h_nEleFinal = new TH1F ("nEleFinal","Number of elecs passing ID+ISO",11,-0.5,10.5);  h_nEleFinal->Sumw2();
   TH1F *h_elec_pT_1st = new TH1F ("elec_pT_1st","elec_pT_1st",100,0,1000); h_elec_pT_1st->Sumw2();
   TH1F *h_elec_pT_2nd = new TH1F ("elec_pT_2nd","elec_pT_2nd",100,0,1000);  h_elec_pT_2nd->Sumw2();
   TH1F *h_elec_pT_1st_ID_ISO = new TH1F ("elec_pT_1st_ID_ISO","pT of highest Reco Electron",100,0,1000); h_elec_pT_1st_ID_ISO->Sumw2();
   TH1F *h_elec_pT_2nd_ID_ISO = new TH1F ("elec_pT_2nd_ID_ISO","pT of second highest Reco Elec",100,0,1000);  h_elec_pT_2nd_ID_ISO->Sumw2();
   TH1F *h_elec_Eta = new TH1F ("elec_Eta","elec_Eta",61,-3.05, 3.05); h_elec_Eta->Sumw2();
   TH1F *h_elec_Eta_1st = new TH1F ("elec_Eta_1st","elec_Eta_1st",61,-3.05, 3.05); h_elec_Eta_1st->Sumw2();
   TH1F *h_elec_Eta_2nd = new TH1F ("elec_Eta_2nd","elec_Eta_2nd",61,-3.05, 3.05); h_elec_Eta->Sumw2();

   // Jets
   TH1F *h_nJet = new TH1F ("nJet","nJet",21,-0.5,20.5);  h_nJet->Sumw2();
   TH1F *h_nJet_ptCut = new TH1F ("nJet_ptCut","nJet_ptCut",11,-0.5,10.5);  h_nJet_ptCut->Sumw2();
   TH1F *h_jet_pT_1st = new TH1F ("jet_pT_1st","jet_pT_1st",100,0,1000); h_jet_pT_1st->Sumw2();
   TH1F *h_jet_pT_2nd = new TH1F ("jet_pT_2nd","jet_pT_2nd",100,0,1000);  h_jet_pT_2nd->Sumw2();
   TH1F *h_jet_Eta_1st = new TH1F ("jet_Eta_1st","jet_Eta_1st",61,-3.05, 3.05);  h_jet_Eta_1st->Sumw2();
   TH1F *h_jet_Eta_2nd = new TH1F ("jet_Eta_2nd","jet_Eta_2nd",61,-3.05, 3.05);  h_jet_Eta_2nd->Sumw2();
   TH1F *h_jet_Eta = new TH1F ("jet_Eta","jet_Eta",61,-3.05, 3.05);  h_jet_Eta->Sumw2();

   //Combinations
   TH1F *h_sT = new TH1F ("sT","sT",100,0,1000);  h_sT->Sumw2();
   TH1F *h_Mej = new TH1F ("Mej","Mej",100,0,1000);  h_Mej->Sumw2();



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

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done every event

     //HLT
     bool PassTrig=false;
     int TrigBit1=getPreCutValue1("FirstTrigBit");
     int TrigBit2=getPreCutValue1("SecondTrigBit");
     if ((HLTResults[TrigBit1])||(HLTResults[TrigBit2]))
	 PassTrig=true;

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


     //GenPart Match
     int electron_PID=11;
     int quark_PID=2;
     int LQ_PID=42;
     //int LQ_PID=23;
     float ConeSizeMCmatch_cut=0.07;
     int nEleMatched = 0;
     int nGenEle = 0;
     int from_LQ_idx[2]={99,99};
     int matched_reco_idx[2]={99,99};

     //////////Finding elec from LQ
     for(int igen=0;igen<GenParticleCount;igen++)
       {
	 //select electrons from LQ decay
	 if(abs(GenParticlePdgId[igen])==electron_PID
	    && abs(GenParticlePdgId[GenParticleMotherIndex[igen]])==LQ_PID)
	   {
	     nGenEle++;
	     h_GenPart_eleEta->Fill(GenParticleEta[igen]);
	     h_GenPart_elePt->Fill(GenParticlePt[igen]);
	     if (from_LQ_idx[0]==99) from_LQ_idx[0]=igen;
	     else from_LQ_idx[1]=igen;
	   }
	 if(abs(GenParticlePdgId[igen])==quark_PID
	    && abs(GenParticlePdgId[GenParticleMotherIndex[igen]])==LQ_PID)
	   {
	     h_GenPart_quarkEta->Fill(GenParticleEta[igen]);
	     h_GenPart_quarkPt->Fill(GenParticlePt[igen]);
	   }
	 if(abs(GenParticlePdgId[igen])==LQ_PID)
	   {
	     h_GenPart_LQEta->Fill(GenParticleEta[igen]);
	     h_GenPart_LQPt->Fill(GenParticlePt[igen]);
	   }
       }

     TVector3 ele1gen, ele2gen;
     ele1gen.SetPtEtaPhi(GenParticlePt[from_LQ_idx[0]],
 			   GenParticleEta[from_LQ_idx[0]],
 			   GenParticlePhi[from_LQ_idx[0]]);	
 	
     ele2gen.SetPtEtaPhi(GenParticlePt[from_LQ_idx[1]],
 			   GenParticleEta[from_LQ_idx[1]],
 			   GenParticlePhi[from_LQ_idx[1]]);	
 	
     ///////Calculating DeltaR and matching
     float minDeltaR1 = 99;
     float minDeltaR2 = 99;
     float min2ndDeltaR1 = 99;
     float min2ndDeltaR2 = 99;
     int RecoIndex1 = 99;
     int RecoIndex2 = 99;

     for (int iele=0; iele<eleCount; iele++) {
       //       if (elePt[iele]<elePt_cut) continue;
       TVector3 ele;
       ele.SetPtEtaPhi(elePt[iele],
 	       eleEta[iele],
 	       elePhi[iele]);	
 
       float DeltaR_ele_ele1gen = ele1gen.DeltaR(ele);
       float DeltaR_ele_ele2gen = ele2gen.DeltaR(ele);

       if (DeltaR_ele_ele1gen < minDeltaR1) {
	 min2ndDeltaR1=minDeltaR1;
	 minDeltaR1=DeltaR_ele_ele1gen;
	 RecoIndex1 = iele;
       }
       else if (DeltaR_ele_ele1gen < min2ndDeltaR1){
	 min2ndDeltaR1=DeltaR_ele_ele1gen;
       }
       if (DeltaR_ele_ele2gen < minDeltaR2) {
	 min2ndDeltaR2=minDeltaR2;
	 minDeltaR2=DeltaR_ele_ele2gen;
	 RecoIndex2 = iele;
       }
       else if (DeltaR_ele_ele2gen < min2ndDeltaR2){
	 min2ndDeltaR2=DeltaR_ele_ele2gen;
       }
     }
     if (minDeltaR1< ConeSizeMCmatch_cut) {
       matched_reco_idx[0]=RecoIndex1;
       nEleMatched++;
     }
     if (minDeltaR2< ConeSizeMCmatch_cut) {
       matched_reco_idx[1]=RecoIndex2;
       nEleMatched++;
     }

     //GenJets
     if(genJetCount>0) {
       h_GenJet_pT_1st->Fill(genJetPt[0]);
       h_GenJet_eta_1st->Fill(genJetEta[0]);
     }
     if(genJetCount>1) {
       h_GenJet_pT_2nd->Fill(genJetPt[1]);
       h_GenJet_eta_2nd->Fill(genJetEta[1]);
     }

     // Electrons
     vector<int> v_idx_ele_FR;
     vector<int> v_idx_ele_HoE;
     vector<int> v_idx_ele_sigEta;
     vector<int> v_idx_ele_deltaEta;
     vector<int> v_idx_ele_deltaPhi;
     vector<int> v_idx_ele_ID;
     vector<int> v_idx_ele_ISO;
     vector<int> v_idx_ele_final;
     int nElePtCut=0;

     for(int iele=0;iele<eleCount;iele++)
       {
	 if (elePt[iele]>10) h_elec_Eta->Fill(eleEta[iele]);

	 // ECAL barrel fiducial region
	 bool pass_ECAL_FR=false;
	 if( fabs(eleEta[iele]) < getPreCutValue1("eleFidRegion") )  pass_ECAL_FR=true;
	 if (pass_ECAL_FR) v_idx_ele_FR.push_back(iele);

	 if (elePt[iele]< getPreCutValue1("ele_pt_cut")) continue;
	 nElePtCut++;

	 bool pass_ISO=false;
	 bool pass_ID=false;

	 bool in_Barrel=false;
         bool in_Endcap=false;
	 if( fabs(eleEta[iele]) < 1.442 )  in_Barrel=true;
	 if( (fabs(eleEta[iele]) > 1.560)&&(fabs(eleEta[iele] < 2.5)) )  in_Endcap=true;

	 // ID cuts
	 bool pass_HoE=false;
	 if( (in_Barrel)&&((eleHoE[iele]) < getPreCutValue1("ID_HoE_bar")) )  pass_HoE=true;
	 if( (in_Endcap)&&((eleHoE[iele]) < getPreCutValue1("ID_HoE_end")) )  pass_HoE=true;
	 if (pass_HoE) v_idx_ele_HoE.push_back(iele);
	 if (in_Barrel) h_elec_HoE_barrel->Fill(eleHoE[iele]);
	 if (in_Endcap) h_elec_HoE_endcap->Fill(eleHoE[iele]);
	 bool pass_sigEtaEta=false;
	 if( (in_Barrel)&&((eleSigmaEE[iele]) < getPreCutValue1("ID_sigEtaEta_bar")) )  pass_sigEtaEta=true;
	 if( (in_Endcap)&&((eleSigmaEE[iele]) < getPreCutValue1("ID_sigEtaEta_end")) )  pass_sigEtaEta=true;
	 if (pass_sigEtaEta) v_idx_ele_sigEta.push_back(iele);
	 if (in_Barrel) h_elec_sigEtaEta_barrel->Fill(eleSigmaEE[iele]);
	 if (in_Endcap) h_elec_sigEtaEta_endcap->Fill(eleSigmaEE[iele]);
	 bool pass_deltaEta=false;
	 if( (in_Barrel)&&(fabs(eleDeltaEtaTrkSC[iele]) < getPreCutValue1("ID_deltaEta_bar")) )  pass_deltaEta=true;
	 if( (in_Endcap)&&(fabs(eleDeltaEtaTrkSC[iele]) < getPreCutValue1("ID_deltaEta_end")) )  pass_deltaEta=true;
	 if (pass_deltaEta) v_idx_ele_deltaEta.push_back(iele);
	 bool pass_deltaPhi=false;
	 if( (in_Barrel)&&(fabs(eleDeltaPhiTrkSC[iele]) < getPreCutValue1("ID_deltaPhi_bar")) )  pass_deltaPhi=true;
	 if( (in_Endcap)&&(fabs(eleDeltaPhiTrkSC[iele]) < getPreCutValue1("ID_deltaPhi_end")) )  pass_deltaPhi=true;
	 if (pass_deltaPhi) v_idx_ele_deltaPhi.push_back(iele);

	 if ((pass_ECAL_FR)&&(pass_HoE)&&(pass_sigEtaEta)&&(pass_deltaEta)&&(pass_deltaPhi)) {
	   v_idx_ele_ID.push_back(iele);
	   pass_ID=true;
	 }

	 // ISO cuts
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

	 ///////For OLD RootTuples comment out Ecal and Hcal Iso above and add these lines in instead:
// 	 bool pass_EcalIso=false;
// 	 float ecal_iso_bar = getPreCutValue1("ISO_EcalIso_bar")+(getPreCutValue2("ISO_EcalIso_bar")*elePt[iele]);
// 	 float ecal_iso_end = getPreCutValue1("ISO_EcalIso_end")+(getPreCutValue2("ISO_EcalIso_bar")*elePt[iele]);
// 	 if( (in_Barrel)&&((eleEcalIso[iele]) < (ecal_iso_bar)) )  pass_EcalIso=true;
// 	 if( (in_Endcap)&&((eleEcalIso[iele]) < (ecal_iso_end)) )  pass_EcalIso=true;
//  	 bool pass_HcalIso=true;

	 if ((pass_ECAL_FR)&&(pass_NumTrack)&&( pass_TrackIso)&&(pass_EcalIso)&&(pass_HcalIso)) {
	   v_idx_ele_ISO.push_back(iele);
	   pass_ISO=true;
	 }
	 
	 if ((pass_ID)&&(pass_ISO)) v_idx_ele_final.push_back(iele);

       } //loop over electrons     
     h_nEle_ptCut->Fill(nElePtCut);

     // Jets
     vector<int> v_idx_jet_final;
     float deltaR_minCut = 0.5;
     int nJet=0;
     for(int ijet=0;ijet<caloJetIC5Count;ijet++)
       {

	 // HCAL barrel fiducial region
	 bool pass_HCAL_FR=false;
	 if( fabs(caloJetIC5Eta[ijet]) < getPreCutValue1("jetFidRegion") ) pass_HCAL_FR=true ;
	 if (caloJetIC5Pt[ijet]<1) continue;

	 ///Disambiguation of electrons from jets
	 float minDeltaR=99;
	 TVector3 jet_vec;
	 jet_vec.SetPtEtaPhi(caloJetIC5Pt[ijet],caloJetIC5Eta[ijet],caloJetIC5Phi[ijet]);
	 for (int i=0; i<v_idx_ele_final.size();i++){
	   TVector3 ele_vec;
	   ele_vec.SetPtEtaPhi(elePt[v_idx_ele_final[i]],eleEta[v_idx_ele_final[i]],elePhi[v_idx_ele_final[i]]);
	   double distance=jet_vec.DeltaR(ele_vec);
	   if (distance<minDeltaR) minDeltaR=distance;
 	 }
	 if (minDeltaR>deltaR_minCut) {
	   h_jet_Eta->Fill(caloJetIC5Eta[ijet]);
	   nJet++;
	 }
	 if (caloJetIC5Pt[ijet]< getPreCutValue1("jet_pt_cut")) continue;
	 if ((pass_HCAL_FR)&&(minDeltaR>deltaR_minCut)) v_idx_jet_final.push_back(ijet);
       }     
     h_nJet_ptCut->Fill(v_idx_jet_final.size());

      // Set the evaluation of the cuts to false and clear the variable values and filled status
     resetCuts();
     
     // Set the value of the variableNames listed in the cutFile to their current value
     if (PassTrig) fillVariableWithValue("HLT",1) ;
     fillVariableWithValue("nEle_FR", v_idx_ele_FR.size()) ;
     fillVariableWithValue("nEle_HoE", v_idx_ele_HoE.size()) ;
     fillVariableWithValue("nEle_sigEtaEta", v_idx_ele_sigEta.size()) ;
     fillVariableWithValue("nEle_deltaEta", v_idx_ele_deltaEta.size()) ;
     fillVariableWithValue("nEle_deltaPhi", v_idx_ele_deltaPhi.size()) ;
     fillVariableWithValue("nEleID_30GeV", v_idx_ele_ID.size()) ;
     fillVariableWithValue("nEleISO_30GeV", v_idx_ele_ISO.size()) ;

     bool TwoEles=false;
     bool TwoJets=false;
     if( v_idx_ele_final.size() >= 1 ) 
       {
	 fillVariableWithValue( "elec_pT_1st", elePt[v_idx_ele_final[0]] );
       }
     if( v_idx_ele_final.size() >= 2 ) 
       {
	 // Calculate Mee
	 TLorentzVector v_ee, ele1, ele2;
	 ele1.SetPtEtaPhiM(elePt[v_idx_ele_final[0]],eleEta[v_idx_ele_final[0]],elePhi[v_idx_ele_final[0]],0);
	 ele2.SetPtEtaPhiM(elePt[v_idx_ele_final[1]],eleEta[v_idx_ele_final[1]],elePhi[v_idx_ele_final[1]],0);
	 v_ee = ele1 + ele2;
	 fillVariableWithValue( "invMass_ee", v_ee.M() ) ;
	 TwoEles=true;
       }
     if( v_idx_jet_final.size() >= 2 )  TwoJets=true;
     if ((TwoEles)&&(TwoJets)){
       double calc_sT=elePt[v_idx_ele_final[0]]+elePt[v_idx_ele_final[1]]+caloJetIC5Pt[v_idx_jet_final[0]]+caloJetIC5Pt[v_idx_jet_final[1]];
       fillVariableWithValue("sT",calc_sT);
       fillVariableWithValue("2e30_2j50", 1) ;
     }
     ///Calculate ej mass and take closer pair (old way)
     double M11, M12, M21, M22=-100;
     double diff_11_22, diff_12_21;
     if ((TwoEles)&&(TwoJets)){
 	 TLorentzVector jet1, jet2, ele1, ele2;
	 ele1.SetPtEtaPhiM(elePt[v_idx_ele_final[0]],eleEta[v_idx_ele_final[0]],elePhi[v_idx_ele_final[0]],0);
	 ele2.SetPtEtaPhiM(elePt[v_idx_ele_final[1]],eleEta[v_idx_ele_final[1]],elePhi[v_idx_ele_final[1]],0);
	 jet1.SetPtEtaPhiM(caloJetIC5Pt[v_idx_jet_final[0]],caloJetIC5Eta[v_idx_jet_final[0]],caloJetIC5Phi[v_idx_jet_final[0]],0);
	 jet2.SetPtEtaPhiM(caloJetIC5Pt[v_idx_jet_final[1]],caloJetIC5Eta[v_idx_jet_final[1]],caloJetIC5Phi[v_idx_jet_final[1]],0);
	 TLorentzVector jet1ele1, jet1ele2, jet2ele1, jet2ele2;
	 jet1ele1 = jet1 + ele1;
	 jet1ele2 = jet1 + ele2;
	 jet2ele1 = jet2 + ele1;
	 jet2ele2 = jet2 + ele2;
	 M11=jet1ele1.M();
	 M12=jet1ele2.M();
	 M21=jet2ele1.M();
	 M22=jet2ele2.M();
	 //cout << "M11:  " << M11 << "  M12:  " << M12 << "  M21:  " << M21 << "  M22:  " << M22 << endl;
	 diff_11_22=fabs(M11-M22);
	 diff_12_21=fabs(M12-M21);
	 if (diff_11_22<diff_12_21) {
	   fillVariableWithValue("Mej",M11);
	 }
	 else {
	   fillVariableWithValue("Mej",M12);
	 }
     }
     // Evaluate cuts (but do not apply them)
     evaluateCuts();

     if( passedCut("all") ){
	 if (diff_11_22<diff_12_21) {
	   h_Mej->Fill(M11);
	   h_Mej->Fill(M22);
	 }
	 else {
	   h_Mej->Fill(M12);
	   h_Mej->Fill(M21);
	 }

     }
     // Fill histograms and do analysis based on cut evaluation
     if (eleCount>0) {
       h_elec_pT_1st->Fill(elePt[0]);
       h_elec_Eta_1st->Fill(eleEta[0]);
     }
     if (eleCount>1) {
       h_elec_pT_2nd->Fill(elePt[1]);
       h_elec_Eta_2nd->Fill(eleEta[1]);
     }
     if (v_idx_ele_final.size()>0) h_elec_pT_1st_ID_ISO->Fill(elePt[v_idx_ele_final[0]]);
     if (v_idx_ele_final.size()>1)h_elec_pT_2nd_ID_ISO->Fill(elePt[v_idx_ele_final[1]]);
     if (v_idx_jet_final.size()>0) {
       h_jet_pT_1st->Fill(caloJetIC5Pt[v_idx_jet_final[0]]);
       h_jet_Eta_1st->Fill(caloJetIC5Eta[v_idx_jet_final[0]]);
     }
     if (v_idx_jet_final.size()>1){
       h_jet_pT_2nd->Fill(caloJetIC5Pt[v_idx_jet_final[1]]);
       h_jet_Eta_2nd->Fill(caloJetIC5Eta[v_idx_jet_final[1]]);
     }
     h_nEleID_30GeV->Fill(v_idx_ele_ID.size());
     h_nEleISO_30GeV->Fill(v_idx_ele_ISO.size());
     h_nEleFinal->Fill(v_idx_ele_final.size());
     h_nJet->Fill(nJet);
     
     // reject events that did not pass level 0 cuts
     if( !passedCut("0") ) continue;
     // ......
     
     // reject events that did not pass level 1 cuts
     if( !passedCut("1") ) continue;
     // ......

     // reject events that did not pass the full cut list
     if( !passedCut("all") ) continue;
     // ......



     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 
   h_GenPart_eleEta->Write();
   h_GenPart_elePt->Write();
   h_GenPart_quarkEta->Write();
   h_GenPart_quarkPt->Write();
   h_GenPart_LQEta->Write();
   h_GenPart_LQPt->Write();
   h_GenJet_pT_1st->Write();
   h_GenJet_eta_1st->Write();
   h_GenJet_pT_2nd->Write();
   h_GenJet_eta_2nd->Write();

   h_elec_Eta->Write();
   h_elec_sigEtaEta_barrel->Write();
   h_elec_sigEtaEta_endcap->Write();
   h_elec_HoE_barrel->Write();
   h_elec_HoE_endcap->Write();
   h_nEle_ptCut->Write();
   h_nEleID_30GeV->Write();
   h_elec_pT_1st->Write();
   h_elec_pT_2nd->Write();
   h_elec_Eta_1st->Write();
   h_elec_Eta_2nd->Write();
   h_elec_pT_1st_ID_ISO->Write();
   h_elec_pT_2nd_ID_ISO->Write();

   h_nJet->Write();
   h_nJet_ptCut->Write();
   h_jet_pT_1st->Write();
   h_jet_pT_2nd->Write();
   h_jet_Eta->Write();
   h_jet_Eta_1st->Write();
   h_jet_Eta_2nd->Write();

   h_Mej->Write();

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
