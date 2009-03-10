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
   TH1F *h_Mej = new TH1F ("Mej","Mej",150,0,300);  h_Mej->Sumw2();
   TH1F *h_Mej_inside = new TH1F ("Mej_inside","Mej_inside",200,0,400);  h_Mej_inside->Sumw2();
   TH1F *h_Mej_above = new TH1F ("Mej_above","Mej_above",200,0,400);  h_Mej_above->Sumw2();
   TH1F *h_Mee = new TH1F ("Mee","Mee",150,0,300);  h_Mee->Sumw2();
   TH2F *h_Mee_Mej = new TH2F ("Mee_Mej","Mee_Mej",45,50,140,200,0,400);

   /////////initialize variables

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
     vector<int> v_idx_ele_final;
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
	 bool pass_ISO=false;
	 bool pass_ID=false;

	 //ID cuts
	 bool pass_HoE=false;
	 if( (in_Barrel)&&((eleHoE[iele]) < getPreCutValue1("ID_HoE_bar")) )  pass_HoE=true;
	 if( (in_Endcap)&&((eleHoE[iele]) < getPreCutValue1("ID_HoE_end")) )  pass_HoE=true;
	 bool pass_sigEtaEta=false;
	 if( (in_Barrel)&&((eleSigmaEE[iele]) < getPreCutValue1("ID_sigEtaEta_bar")) )  pass_sigEtaEta=true;
	 if( (in_Endcap)&&((eleSigmaEE[iele]) < getPreCutValue1("ID_sigEtaEta_end")) )  pass_sigEtaEta=true;
	 bool pass_deltaEta=false;
	 if( (in_Barrel)&&((eleDeltaEtaTrkSC[iele]) < getPreCutValue1("ID_deltaEta_bar")) )  pass_deltaEta=true;
	 if( (in_Endcap)&&((eleDeltaEtaTrkSC[iele]) < getPreCutValue1("ID_deltaEta_end")) )  pass_deltaEta=true;
	 bool pass_deltaPhi=false;
	 if( (in_Barrel)&&((eleDeltaPhiTrkSC[iele]) < getPreCutValue1("ID_deltaPhi_bar")) )  pass_deltaPhi=true;
	 if( (in_Endcap)&&((eleDeltaPhiTrkSC[iele]) < getPreCutValue1("ID_deltaPhi_end")) )  pass_deltaPhi=true;

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
	 
	 if ( (pass_ID) && (pass_ISO) ) v_idx_ele_final.push_back(iele);

       } //loop over electrons     

     //cout << "Jets" << endl;

     //## Jets
     vector<int> v_idx_jet_noCut;
     vector<int> v_idx_jet;
     vector<int> v_idx_jet_final;
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
	 for (int i=0; i < v_idx_ele_final.size(); i++){
	   TVector3 ele_vec;
	   ele_vec.SetPtEtaPhi(elePt[v_idx_ele_final[i]],eleEta[v_idx_ele_final[i]],elePhi[v_idx_ele_final[i]]);
	   double distance = jet_vec.DeltaR(ele_vec);
	   if (distance<minDeltaR) minDeltaR=distance;
 	 }
	 
	 if ( minDeltaR > deltaR_minCut )  v_idx_jet_final.push_back(ijet);

       }     

     // Set the evaluation of the cuts to false and clear the variable values and filled status
     resetCuts();
     
     // Set the value of the variableNames listed in the cutFile to their current value

     fillVariableWithValue( "HLT", PassTrig ) ;

     if( v_idx_ele_final.size() >= 1 ) 
       {
	 fillVariableWithValue( "Pt1stEleIDISO", elePt[v_idx_ele_final[0]] );
	 fillVariableWithValue( "Eta1stEleIDISO", eleEta[v_idx_ele_final[0]] );
       }
     if( v_idx_ele_final.size() >= 2 ) 
       {
	 fillVariableWithValue( "Pt2ndEleIDISO", elePt[v_idx_ele_final[1]] );
	 fillVariableWithValue( "Eta2ndEleIDISO", eleEta[v_idx_ele_final[1]] );
       }
     if( v_idx_jet_final.size() >= 1 ) 
       {
	 fillVariableWithValue( "Pt1stJet_DIS", caloJetIC5Pt[v_idx_jet_final[0]] );
	 fillVariableWithValue( "Eta1stJet_DIS", caloJetIC5Eta[v_idx_jet_final[0]] );
       }
     if( v_idx_jet_final.size() >= 2 ) 
       {
	 fillVariableWithValue( "Pt2ndJet_DIS", caloJetIC5Pt[v_idx_jet_final[1]] );
	 fillVariableWithValue( "Eta2ndJet_DIS", caloJetIC5Eta[v_idx_jet_final[1]] );
       }


     //## define "2ele" and "2jets" booleans
     bool TwoEles=false;
     bool TwoJets=false;
     if( v_idx_ele_final.size() >= 2 ) TwoEles = true;
     if( v_idx_jet_final.size() >= 2 ) TwoJets = true;

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

     //## Mej "old" algorithm
     double M11, M12, M21, M22, Mee = -999;
     double diff_11_22, diff_12_21;
     if ( (TwoEles) && (TwoJets) ) 
       {
 	 TLorentzVector jet1, jet2, ele1, ele2;
	 ele1.SetPtEtaPhiM(elePt[v_idx_ele_final[0]],
			   eleEta[v_idx_ele_final[0]],
			   elePhi[v_idx_ele_final[0]],0);
	 ele2.SetPtEtaPhiM(elePt[v_idx_ele_final[1]],
			   eleEta[v_idx_ele_final[1]],
			   elePhi[v_idx_ele_final[1]],0);
	 jet1.SetPtEtaPhiM(caloJetIC5Pt[v_idx_jet_final[0]],
			   caloJetIC5Eta[v_idx_jet_final[0]],
			   caloJetIC5Phi[v_idx_jet_final[0]],0);
	 jet2.SetPtEtaPhiM(caloJetIC5Pt[v_idx_jet_final[1]],
			   caloJetIC5Eta[v_idx_jet_final[1]],
			   caloJetIC5Phi[v_idx_jet_final[1]],0);
	 TLorentzVector jet1ele1, jet1ele2, jet2ele1, jet2ele2, ele1ele2;
	 jet1ele1 = jet1 + ele1;
	 jet1ele2 = jet1 + ele2;
	 jet2ele1 = jet2 + ele1;
	 jet2ele2 = jet2 + ele2;
	 ele1ele2 = ele1 + ele2;
	 M11 = jet1ele1.M();
	 M12 = jet1ele2.M();
	 M21 = jet2ele1.M();
	 M22 = jet2ele2.M();
	 Mee = ele1ele2.M();
	 fillVariableWithValue("Mee_control",Mee);
	 fillVariableWithValue("Mee",Mee);
	 //cout << "M11:  " << M11 << "  M12:  " << M12 << "  M21:  " << M21 << "  M22:  " << M22 << endl;
	 diff_11_22=fabs(M11-M22);
	 diff_12_21=fabs(M12-M21);
	 if (diff_11_22<diff_12_21) 
	   {
	     if(M11>M22)
	       {
		 fillVariableWithValue("MejMax", M11);	       
		 fillVariableWithValue("MejMin", M22);
	       }
	     else
	       {
		 fillVariableWithValue("MejMax", M22);	       
		 fillVariableWithValue("MejMin", M11);
	       }
	   }
	 else 
	   {
	     if(M12>M21)
	       {
		 fillVariableWithValue("MejMax", M12);	       
		 fillVariableWithValue("MejMin", M21);
	       }
	     else
	       {
		 fillVariableWithValue("MejMax", M21);	       
		 fillVariableWithValue("MejMin", M12);
	       }
	   }

       }

     //--------------------------------------------------

     // Evaluate cuts (but do not apply them)
     evaluateCuts();

     // Fill histograms and do analysis based on cut evaluation

     if( passedCut("1") )
       {
	 if (diff_11_22<diff_12_21) 
	   {
	     h_Mej->Fill(M11);
	     h_Mej->Fill(M22);
	     h_Mee_Mej->Fill(Mee,M11);
	     h_Mee_Mej->Fill(Mee,M22);
	   }
	 else 
	   {
	     h_Mej->Fill(M12);
	     h_Mej->Fill(M21);
	     h_Mee_Mej->Fill(Mee,M12);
	     h_Mee_Mej->Fill(Mee,M21);

	   }
	 h_Mee->Fill(Mee);
       }
     
     if( (passedCut("1"))&&(passedCut("2")) )
       {
	 if (diff_11_22<diff_12_21) 
	   {
	     h_Mej_inside->Fill(M11);
	     h_Mej_inside->Fill(M22);
	   }
	 else 
	   {
	     h_Mej_inside->Fill(M12);
	     h_Mej_inside->Fill(M21);
	   }
       }
     
     if( (passedCut("1"))&&(passedCut("3")) )
       {
	 if (diff_11_22<diff_12_21) 
	   {
	     h_Mej_above->Fill(M11);
	     h_Mej_above->Fill(M22);
	   }
	 else 
	   {
	     h_Mej_above->Fill(M12);
	     h_Mej_above->Fill(M21);
	   }
       }
     
     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 
   h_Mee->Write();
   h_Mej->Write();
   h_Mej_inside->Write();
   h_Mej_above->Write();
   h_Mee_Mej->Write();

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
