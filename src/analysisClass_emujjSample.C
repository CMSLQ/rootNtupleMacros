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

   //Combinations
   TH1F *h_Mlj = new TH1F ("Mlj","Mlj",200,0,2000);  h_Mlj->Sumw2();

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<50000;jentry++) {
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
     vector<int> v_idx_ele_final;

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

     //## Muons
     vector<int> v_idx_muon_noCut;
     // here pre-cut on pT is applied
     vector<int> v_idx_muon;
     vector<int> v_idx_muon_ID;
     vector<int> v_idx_muon_ISO;
     vector<int> v_idx_muon_ID_ISO;
     vector<int> v_idx_muon_DIS;
     vector<int> v_idx_muon_final;
     float deltaR_minCut_2 = getPreCutValue1("muon_jet_DeltaRcut");     
     for(int imuon=0;imuon<muonCount;imuon++)
       {
	 //no cut on reco muons
	 v_idx_muon_noCut.push_back(imuon);
	 
	 //pT pre-cut on reco muons
	 if ( muonPt[imuon] < getPreCutValue1("muon_PtPreCut") ) continue;
	 
	 v_idx_muon.push_back(imuon);
	 
	 //barrel-endcap definition
	 bool in_Barrel=false;
         bool in_Endcap=false;
	 if( fabs(muonEta[imuon]) < getPreCutValue1("muonEta_bar") )  in_Barrel=true;
	 if( ( fabs(muonEta[imuon]) > getPreCutValue1("muonEta_end") )
	     && 
	     ( fabs(muonEta[imuon]) < getPreCutValue2("muonEta_end") ) 
	     )  in_Endcap=true;
	 
	 //ID and ISO booleans
	 bool pass_ISO=false;
	 bool pass_ID=false;
	 bool pass_DIS=false;
	 
	 //ID cuts
	 bool pass_TrkD0=false;
	 if( (in_Barrel)&&((muonTrkD0[imuon]) < getPreCutValue1("ID_muonTrkD0_bar") ) )  pass_TrkD0=true;
	 if( (in_Endcap)&&((muonTrkD0[imuon]) < getPreCutValue1("ID_muonTrkD0_end") ) )  pass_TrkD0=true;
	 
	 bool pass_TrkHits=false;
	 if( (in_Barrel)&&((muonTrkHits[imuon]) >= getPreCutValue1("ID_muonTrkHits_bar") ) )  pass_TrkHits=true;
	 if( (in_Endcap)&&((muonTrkHits[imuon]) >= getPreCutValue1("ID_muonTrkHits_end") ) )  pass_TrkHits=true;
	 
	 bool pass_GlobalChi2=false;
	 if( (in_Barrel)&&((muonGlobalChi2[imuon]) >= getPreCutValue1("ID_muonGlobalChi2_bar") ) )  pass_GlobalChi2=true;
	 if( (in_Endcap)&&((muonGlobalChi2[imuon]) >= getPreCutValue1("ID_muonGlobalChi2_end") ) )  pass_GlobalChi2=true;
	 
	 if ( (pass_TrkD0) && (pass_TrkHits) /*&& (pass_GlobalChi2) */ ) 
	   {
	     v_idx_muon_ID.push_back(imuon);
	     pass_ID=true;
	   }

	 //ISO cuts
	 bool pass_TrkIso=false;
	 if( (in_Barrel)&&((muonTrkIso[imuon]) < getPreCutValue1("ISO_muonTrkIso_bar") ) )  pass_TrkIso=true;
	 if( (in_Endcap)&&((muonTrkIso[imuon]) < getPreCutValue1("ISO_muonTrkIso_end") ) )  pass_TrkIso=true;

	 //ignore it for the moment
	 // 	 if ( pass_TrkIso ) 
	 // 	   {
	 // 	     v_idx_muon_ISO.push_back(imuon);
	 // 	     pass_ISO=true;
	 // 	   }
	 pass_ISO=true;
	 v_idx_muon_ISO.push_back(imuon);

	 if ( (pass_ID) && (pass_ISO) ) v_idx_muon_ID_ISO.push_back(imuon);
	
	 //Disambiguation of jets from muons
	 float minDeltaR=9999;
	 TVector3 muon_vec;
	 muon_vec.SetPtEtaPhi(muonPt[imuon],muonEta[imuon],muonPhi[imuon]);
	 for (int i=0; i < v_idx_jet_final.size(); i++){
	   TVector3 jet_vec;
	   jet_vec.SetPtEtaPhi(caloJetIC5Pt[v_idx_jet_final[i]],caloJetIC5Eta[v_idx_jet_final[i]],caloJetIC5Phi[v_idx_jet_final[i]]);
	   double distance = muon_vec.DeltaR(jet_vec);
	   if (distance<minDeltaR) minDeltaR=distance;
	 }
	 
	 if ( minDeltaR > deltaR_minCut_2 )  
	   pass_DIS = true;

	 if ( (pass_DIS) ) v_idx_muon_DIS.push_back(imuon);
	 
	 if ( (pass_ID) && (pass_ISO) && (pass_DIS) ) v_idx_muon_final.push_back(imuon);

       }// end loop over muons


     //-----------

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
     fillVariableWithValue( "nEle_PtPreCut_ID", v_idx_ele_ID.size() ) ;
     fillVariableWithValue( "nEle_PtPreCut_ISO", v_idx_ele_ISO.size() ) ;
     fillVariableWithValue( "nEle_PtPreCut_IDISO", v_idx_ele_final.size() ) ;

     //cout << "nJet" << endl;

     //## nJet
     fillVariableWithValue( "nJet_noCut", v_idx_jet_noCut.size() ) ;
     fillVariableWithValue( "nJet_PtPreCut", v_idx_jet.size() ) ;
     fillVariableWithValue( "nJet_PtPreCut_DIS", v_idx_jet_final.size() ) ;

     //cout << "nMuon" << endl;

     //## nMuon
     fillVariableWithValue( "nMuon_noCut", v_idx_muon_noCut.size() ) ;
     fillVariableWithValue( "nMuon_PtPreCut", v_idx_muon.size() ) ;
     fillVariableWithValue( "nMuon_PtPreCut_ID", v_idx_muon_ID.size() ) ;
     fillVariableWithValue( "nMuon_PtPreCut_ISO", v_idx_muon_ISO.size() ) ;
     fillVariableWithValue( "nMuon_PtPreCut_DIS", v_idx_muon_DIS.size() ) ;
     fillVariableWithValue( "nMuon_PtPreCut_IDISO_DIS", v_idx_muon_final.size() ) ;

     //## nEMu
     if(v_idx_muon_noCut.size()>=1 && v_idx_ele_noCut.size()>=1 )
       fillVariableWithValue( "nEMu_noCut", v_idx_muon_noCut.size()+v_idx_ele_noCut.size() ) ;

     if(v_idx_muon.size()>=1 && v_idx_ele.size()>=1 )
       fillVariableWithValue( "nEMu_PtPreCut", v_idx_muon.size()+v_idx_ele.size() ) ;

     if(v_idx_muon_ID.size()>=1 && v_idx_ele_ID.size()>=1 )
       fillVariableWithValue( "nEMu_PtPreCut_ID", v_idx_muon_ID.size()+v_idx_ele_ID.size() ) ;

     if(v_idx_muon_ISO.size()>=1 && v_idx_ele_ISO.size()>=1 )
       fillVariableWithValue( "nEMu_PtPreCut_ISO", v_idx_muon_ISO.size()+v_idx_ele_ISO.size() ) ;

     if(v_idx_muon_DIS.size()>=1 && v_idx_ele.size()>=1 )
       fillVariableWithValue( "nEMu_PtPreCut_DIS", v_idx_muon_DIS.size()+v_idx_ele.size() ) ;
     
     if(v_idx_muon_final.size()>=1 && v_idx_ele_final.size()>=1 )
       fillVariableWithValue( "nEMu_PtPreCut_IDISO_DIS", v_idx_muon_final.size()+v_idx_ele_final.size() ) ;

     //cout << "1st and 2nd Emu" << endl;

     //## 1st and 2nd Emu
     if( v_idx_ele.size() >=1 && v_idx_muon.size()>=1 ) 
       {
	 if( muonPt[v_idx_muon[0]] > elePt[v_idx_ele[0]] )
	   {
	     fillVariableWithValue( "Pt1stEMu", muonPt[v_idx_muon[0]] );
	     fillVariableWithValue( "Pt2ndEMu", elePt[v_idx_ele[0]] );
	   }
	 else
	   {
	     fillVariableWithValue( "Pt1stEMu", elePt[v_idx_ele[0]] );
	     fillVariableWithValue( "Pt2ndEMu", muonPt[v_idx_muon[0]] );
	   }
       }
     if( v_idx_ele_ID.size() >=1 && v_idx_muon_ID.size()>=1 ) 
       {
	 if( muonPt[v_idx_muon_ID[0]] > elePt[v_idx_ele_ID[0]] )
	   {
	     fillVariableWithValue( "Pt1stEMuID", muonPt[v_idx_muon_ID[0]] );
	     fillVariableWithValue( "Pt2ndEMuID", elePt[v_idx_ele_ID[0]] );
	   }
	 else
	   {
	     fillVariableWithValue( "Pt1stEMuID", elePt[v_idx_ele_ID[0]] );
	     fillVariableWithValue( "Pt2ndEMuID", muonPt[v_idx_muon_ID[0]] );
	   }
       }
     if( v_idx_ele_ISO.size() >=1 && v_idx_muon_ISO.size()>=1 ) 
       {
	 if( muonPt[v_idx_muon_ISO[0]] > elePt[v_idx_ele_ISO[0]] )
	   {
	     fillVariableWithValue( "Pt1stEMuISO", muonPt[v_idx_muon_ISO[0]] );
	     fillVariableWithValue( "Pt2ndEMuISO", elePt[v_idx_ele_ISO[0]] );
	   }
	 else
	   {
	     fillVariableWithValue( "Pt1stEMuISO", elePt[v_idx_ele_ISO[0]] );
	     fillVariableWithValue( "Pt2ndEMuISO", muonPt[v_idx_muon_ISO[0]] );
	   }
       }
     if( v_idx_ele.size() >=1 && v_idx_muon_DIS.size()>=1 ) 
       {
	 if( muonPt[v_idx_muon_DIS[0]] > elePt[v_idx_ele[0]] )
	   {
	     fillVariableWithValue( "Pt1stEMuDIS", muonPt[v_idx_muon_DIS[0]] );
	     fillVariableWithValue( "Pt2ndEMuDIS", elePt[v_idx_ele[0]] );
	   }
	 else
	   {
	     fillVariableWithValue( "Pt1stEMuDIS", elePt[v_idx_ele[0]] );
	     fillVariableWithValue( "Pt2ndEMuDIS", muonPt[v_idx_muon_DIS[0]] );
	   }
       }
     if( v_idx_ele_final.size() >=1 && v_idx_muon_final.size()>=1 ) 
       {
	 if( muonPt[v_idx_muon_final[0]] > elePt[v_idx_ele_final[0]] )
	   {
	     fillVariableWithValue( "Pt1stEMuIDISO_DIS", muonPt[v_idx_muon_final[0]] );
	     fillVariableWithValue( "Pt2ndEMuIDISO_DIS", elePt[v_idx_ele_final[0]] );
	     fillVariableWithValue( "Eta1stEMuIDISO_DIS", muonEta[v_idx_muon_final[0]] );
	     fillVariableWithValue( "Eta2ndEMuIDISO_DIS", eleEta[v_idx_ele_final[0]] );
	   }
	 else
	   {
	     fillVariableWithValue( "Pt1stEMuIDISO_DIS", elePt[v_idx_ele_final[0]] );
	     fillVariableWithValue( "Pt2ndEMuIDISO_DIS", muonPt[v_idx_muon_final[0]] );
	     fillVariableWithValue( "Eta1stEMuIDISO_DIS", eleEta[v_idx_ele_final[0]] );
	     fillVariableWithValue( "Eta2ndEMuIDISO_DIS", muonEta[v_idx_muon_final[0]] );
	   }
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


     //## define "2ele" and "2jets" booleans
     bool OneEleOneMu=false;
     bool TwoJets=false;
     if( v_idx_ele_final.size() >= 1 && v_idx_muon_final.size() >= 1 ) OneEleOneMu = true;
     if( v_idx_jet_final.size() >= 2 ) TwoJets = true;

     //cout << "Mee" << endl;

     //## Mee
     if( OneEleOneMu ) 
       {
	 TLorentzVector v_emu, ele, muon;
	 ele.SetPtEtaPhiM(elePt[v_idx_ele_final[0]],
			  eleEta[v_idx_ele_final[0]],
			  elePhi[v_idx_ele_final[0]],0);
	 muon.SetPtEtaPhiM(muonPt[v_idx_ele_final[0]],
			   muonEta[v_idx_ele_final[0]],
			   muonPhi[v_idx_ele_final[0]],0);
	 v_emu = ele + muon;
	 fillVariableWithValue( "invMass_emu", v_emu.M() ) ;
       }

     //cout << "ST" << endl;

     //## ST
     if ( (OneEleOneMu) && (TwoJets) ) 
       {
	 double calc_sT = 
	   elePt[v_idx_ele_final[0]]
	   + muonPt[v_idx_muon_final[0]]
	   + caloJetIC5Pt[v_idx_jet_final[0]]
	   + caloJetIC5Pt[v_idx_jet_final[1]];
	 fillVariableWithValue("sT", calc_sT);
       }

     //cout << "Mlj" << endl;

     //## Mlj "old" algorithm
     double M11, M12, M21, M22 = -999;
     double diff_11_22, diff_12_21;
     if ( (OneEleOneMu) && (TwoJets) ) 
       {
 	 TLorentzVector jet1, jet2, ele, muon;
	 ele.SetPtEtaPhiM(elePt[v_idx_ele_final[0]],
			   eleEta[v_idx_ele_final[0]],
			   elePhi[v_idx_ele_final[0]],0);
	 muon.SetPtEtaPhiM(muonPt[v_idx_muon_final[0]],
			   muonEta[v_idx_muon_final[0]],
			   muonPhi[v_idx_muon_final[0]],0);
	 jet1.SetPtEtaPhiM(caloJetIC5Pt[v_idx_jet_final[0]],
			   caloJetIC5Eta[v_idx_jet_final[0]],
			   caloJetIC5Phi[v_idx_jet_final[0]],0);
	 jet2.SetPtEtaPhiM(caloJetIC5Pt[v_idx_jet_final[1]],
			   caloJetIC5Eta[v_idx_jet_final[1]],
			   caloJetIC5Phi[v_idx_jet_final[1]],0);
	 TLorentzVector jet1ele, jet1muon, jet2ele, jet2muon;
	 jet1ele = jet1 + ele;
	 jet1muon = jet1 + muon;
	 jet2ele = jet2 + ele;
	 jet2muon = jet2 + muon;
	 M11 = jet1ele.M();
	 M12 = jet1muon.M();
	 M21 = jet2ele.M();
	 M22 = jet2muon.M();
	 //cout << "M11:  " << M11 << "  M12:  " << M12 << "  M21:  " << M21 << "  M22:  " << M22 << endl;
	 diff_11_22=fabs(M11-M22);
	 diff_12_21=fabs(M12-M21);
	 if (diff_11_22<diff_12_21) 
	   {
	     if(M11>M22)
	       {
		 fillVariableWithValue("MljMax", M11);	       
		 fillVariableWithValue("MljMin", M22);
	       }
	     else
	       {
		 fillVariableWithValue("MljMax", M22);	       
		 fillVariableWithValue("MljMin", M11);
	       }
	   }
	 else 
	   {
	     if(M12>M21)
	       {
		 fillVariableWithValue("MljMax", M12);	       
		 fillVariableWithValue("MljMin", M21);
	       }
	     else
	       {
		 fillVariableWithValue("MljMax", M21);	       
		 fillVariableWithValue("MljMin", M12);
	       }
	   }

       }

     //--------------------------------------------------

     // Evaluate cuts (but do not apply them)
     evaluateCuts();

     // Fill histograms and do analysis based on cut evaluation

     if( passedCut("all") )
       {
	 if (diff_11_22<diff_12_21) 
	   {
	     h_Mlj->Fill(M11);
	     h_Mlj->Fill(M22);
	   }
	 else 
	   {
	     h_Mlj->Fill(M12);
	     h_Mlj->Fill(M21);
	   }
       }
     
     // reject events that did not pass level 0 cuts
     // if( !passedCut("0") ) continue;
     // ......
     
     // reject events that did not pass level 1 cuts
     // if( !passedCut("1") ) continue;
     // ......

     // reject events that did not pass the full cut list
     // if( !passedCut("all") ) continue;
     // ......



     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 
   h_Mlj->Write();

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
