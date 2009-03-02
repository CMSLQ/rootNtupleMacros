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

   TH1F *h_Nevent = new TH1F ("Nevent","N events passing",1000000,0,1000000);
   h_Nevent->Sumw2();

   /////////initialize variables

   ////initialize optimization multi-dimensional array
   int opt[10][10][10][10][10][10];

     for (int i=0;i<10;i++){
	 for (int j=0;j<10;j++){
	     for (int k=0;k<10;k++){
		 for (int l=0;l<10;l++){
		   for (int m=0;m<10;m++){
		     for (int n=0;n<10;n++){
		     opt[i][j][k][l][m][n]=0;
		     }
		   }
		 }
	     }
	 }
     }

     /////arrays of thresholds to be tested 

     double ele1pTvec[10] = {0,0,0,0,0,0,0,0,0,0};
     double ele1pTlow = getPreCutValue1("OptEle1Pt");
     double ele1pThigh = getPreCutValue2("OptEle1Pt");
     double ele1pTsteps = getPreCutValue3("OptEle1Pt");
     for (int step=0;step<10;step++){
      double tmp=ele1pTlow + (step*ele1pTsteps);
      if (tmp>=ele1pThigh)tmp=ele1pThigh;
      ele1pTvec[step]=tmp;
     }

     double ele2pTvec[10] = {0,0,0,0,0,0,0,0,0,0};
     double ele2pTlow = getPreCutValue1("OptEle2Pt");
     double ele2pThigh = getPreCutValue2("OptEle2Pt");
     double ele2pTsteps = getPreCutValue3("OptEle2Pt");
     for (int step=0;step<10;step++){
      double tmp=ele2pTlow + (step*ele2pTsteps);
      if (tmp>=ele2pThigh)tmp=ele2pThigh;
      ele2pTvec[step]=tmp;
     }

     double jet1pTvec[10] = {0,0,0,0,0,0,0,0,0,0};
     double jet1pTlow = getPreCutValue1("OptJet1Pt");
     double jet1pThigh = getPreCutValue2("OptJet1Pt");
     double jet1pTsteps = getPreCutValue3("OptJet1Pt");
     for (int step=0;step<10;step++){
      double tmp=jet1pTlow + (step*jet1pTsteps);
      if (tmp>=jet1pThigh)tmp=jet1pThigh;
      jet1pTvec[step]=tmp;
     }

     double jet2pTvec[10] = {0,0,0,0,0,0,0,0,0,0};
     double jet2pTlow = getPreCutValue1("OptJet2Pt");
     double jet2pThigh = getPreCutValue2("OptJet2Pt");
     double jet2pTsteps = getPreCutValue3("OptJet2Pt");
     for (int step=0;step<10;step++){
      double tmp=jet2pTlow + (step*jet2pTsteps);
      if (tmp>=jet2pThigh)tmp=jet2pThigh;
      jet2pTvec[step]=tmp;
     }

     double MeeLowvec[10] = {0,0,0,0,0,0,0,0,0,0};
     double MeeLowlow = getPreCutValue1("OptMeeLow");
     double MeeLowhigh = getPreCutValue2("OptMeeLow");
     double MeeLowsteps = getPreCutValue3("OptMeeLow");
     for (int step=0;step<10;step++){
      double tmp=MeeLowlow + (step*MeeLowsteps);
      if (tmp>=MeeLowhigh)tmp=MeeLowhigh;
      MeeLowvec[step]=tmp;
     }

     double sTvec[10] = {0,0,0,0,0,0,0,0,0,0};
     double sTlow = getPreCutValue1("OptSt");
     double sThigh = getPreCutValue2("OptSt");
     double sTsteps = getPreCutValue3("OptSt");
     for (int step=0;step<10;step++){
      double tmp=sTlow + (step*sTsteps);
      if (tmp>=sThigh)tmp=sThigh;
      sTvec[step]=tmp;
     }

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

     //////////////HLT
     bool PassTrig=false;
     int TrigBit1=getPreCutValue1("FirstTrigBit");
     int TrigBit2=getPreCutValue1("SecondTrigBit");
     if ((HLTResults[TrigBit1])||(HLTResults[TrigBit2]))
	 PassTrig=true;
     
     //////////////Identify Good Electrons (pass ID and ISO)
     vector<int> v_idx_ele_ID_ISO;

     for(int iele=0;iele<eleCount;iele++)
       {
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

	 bool pass_sigEtaEta=false;
	 if( (in_Barrel)&&((eleSigmaEE[iele]) < getPreCutValue1("ID_sigEtaEta_bar")) )  pass_sigEtaEta=true;
	 if( (in_Endcap)&&((eleSigmaEE[iele]) < getPreCutValue1("ID_sigEtaEta_end")) )  pass_sigEtaEta=true;

	 bool pass_deltaEta=false;
	 if( (in_Barrel)&&((eleDeltaEtaTrkSC[iele]) < getPreCutValue1("ID_deltaEta_bar")) )  pass_deltaEta=true;
	 if( (in_Endcap)&&((eleDeltaEtaTrkSC[iele]) < getPreCutValue1("ID_deltaEta_end")) )  pass_deltaEta=true;

	 bool pass_deltaPhi=false;
	 if( (in_Barrel)&&((eleDeltaPhiTrkSC[iele]) < getPreCutValue1("ID_deltaPhi_bar")) )  pass_deltaPhi=true;
	 if( (in_Endcap)&&((eleDeltaPhiTrkSC[iele]) < getPreCutValue1("ID_deltaPhi_end")) )  pass_deltaPhi=true;

	 if ((pass_HoE)&&(pass_sigEtaEta)&&(pass_deltaEta)&&(pass_deltaPhi)) pass_ID=true;

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

	 if ((pass_NumTrack)&&( pass_TrackIso)&&(pass_EcalIso)&&(pass_HcalIso))  pass_ISO=true;
	 
	 if ((pass_ID)&&(pass_ISO)) v_idx_ele_ID_ISO.push_back(iele);

       } //loop over electrons     


     /////////////Clean Jet collection of electrons

     // Jets
     vector<int> v_idx_jet_final;
     float deltaR_minCut = 0.5;

     for(int ijet=0;ijet<caloJetIC5Count;ijet++)
       {
	 if (caloJetIC5Pt[ijet]==0) continue;

	 ///Disambiguation of electrons from jets
	 float minDeltaR=99;
	 TVector3 jet_vec;
	 jet_vec.SetPtEtaPhi(caloJetIC5Pt[ijet],caloJetIC5Eta[ijet],caloJetIC5Phi[ijet]);
	 for (int i=0; i<v_idx_ele_ID_ISO.size();i++){
	   TVector3 ele_vec;
	   ele_vec.SetPtEtaPhi(elePt[v_idx_ele_ID_ISO[i]],eleEta[v_idx_ele_ID_ISO[i]],elePhi[v_idx_ele_ID_ISO[i]]);
	   double distance=jet_vec.DeltaR(ele_vec);
	   if (distance<minDeltaR) minDeltaR=distance;
 	 }
	 if (minDeltaR>deltaR_minCut) v_idx_jet_final.push_back(ijet);
       }     


     /////////////Decide here which jets and electrons you want to pair
     /////////////For now take leading 2 jets and electrons and choose pair with smallest mass difference

     ///Make sure we have at least 2 good electrons and 2 good jets
     
     bool TwoEles=false;
     bool TwoJets=false;
     if( v_idx_jet_final.size() >= 2 )  TwoJets=true;
     if( v_idx_ele_ID_ISO.size() >= 2 )  TwoEles=true;

     ///Calculate ej mass and take closer pair (old way)
     double M11, M12, M21, M22=-100;
     double diff_11_22, diff_12_21;
     int pair_1_ele_idx=99, pair_1_jet_idx=99, pair_2_ele_idx=99, pair_2_jet_idx=99;
     if ((TwoEles)&&(TwoJets)){
       //cout << elePt[v_idx_ele_ID_ISO[0]] << "\t" << elePt[v_idx_ele_ID_ISO[1]] << "\t" << caloJetIC5Pt[v_idx_jet_final[0]] << "\t" << caloJetIC5Pt[v_idx_jet_final[1]] << endl;
 	 TLorentzVector jet1, jet2, ele1, ele2;
	 ele1.SetPtEtaPhiM(elePt[v_idx_ele_ID_ISO[0]],eleEta[v_idx_ele_ID_ISO[0]],elePhi[v_idx_ele_ID_ISO[0]],0);
	 ele2.SetPtEtaPhiM(elePt[v_idx_ele_ID_ISO[1]],eleEta[v_idx_ele_ID_ISO[1]],elePhi[v_idx_ele_ID_ISO[1]],0);
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
	   pair_1_ele_idx=v_idx_ele_ID_ISO[0];
	   pair_1_jet_idx=v_idx_jet_final[0];
	   pair_2_ele_idx=v_idx_ele_ID_ISO[1];
	   pair_2_jet_idx=v_idx_jet_final[1];
	 }
	 else {
	   pair_1_ele_idx=v_idx_ele_ID_ISO[0];
	   pair_1_jet_idx=v_idx_jet_final[1];
	   pair_2_ele_idx=v_idx_ele_ID_ISO[1];
	   pair_2_jet_idx=v_idx_jet_final[0];
	 }
     }  // end if 2 ele and 2 jets


     //////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////  OPTIMIZATION CODE  /////////////////////////////////////////
     /////////////////////////////////////////////////////////////////////////////////////////

     ////get Jet Energy Factor from cut file (for systematic error estimate)
     float PtScale = getPreCutValue1("EnergyFactor");

     //// indecies of the electrons and jets, paired according to selection algorithm

     int First_Pair_ele_idx = 99;
     int First_Pair_jet_idx = 99;
     int Sec_Pair_ele_idx = 99;
     int Sec_Pair_jet_idx = 99;

     //////assign these indicies according to what was determined above...

     First_Pair_ele_idx = pair_1_ele_idx;
     First_Pair_jet_idx = pair_1_jet_idx;
     Sec_Pair_ele_idx = pair_2_ele_idx;
     Sec_Pair_jet_idx = pair_2_jet_idx;

//      for (int i=0;i<10;i++) 
//        cout << ele1pTvec[i] << "\t" ;
//      cout << endl;
//      for (int i=0;i<10;i++) 
//        cout << ele2pTvec[i] << "\t" ;
//      cout << endl;
//      for (int i=0;i<10;i++) 
//        cout << jet1pTvec[i] << "\t" ;
//      cout << endl;
//      for (int i=0;i<10;i++) 
//        cout << jet2pTvec[i] << "\t" ;
//      cout << endl;
//      for (int i=0;i<10;i++) 
//        cout << MeeLowvec[i] << "\t" ;
//      cout << endl;
//      for (int i=0;i<10;i++) 
//        cout << sTvec[i] << "\t" ;
//      cout << endl;


   if (TwoEles&&TwoJets){
     ////get object values just once
     double ele1pT = elePt[First_Pair_ele_idx];
     double ele2pT = elePt[Sec_Pair_ele_idx];
     double jet1pT = PtScale * caloJetIC5Pt[First_Pair_jet_idx];
     double jet2pT = PtScale * caloJetIC5Pt[Sec_Pair_jet_idx];
     //     cout << jet1pT << "\t" << jet2pT << endl;

     TLorentzVector v_ee, ele1, ele2;
     ele1.SetPtEtaPhiM(elePt[First_Pair_ele_idx],
  		   eleEta[First_Pair_ele_idx],
  		   elePhi[First_Pair_ele_idx],0);
     ele2.SetPtEtaPhiM(elePt[Sec_Pair_ele_idx],
  		   eleEta[Sec_Pair_ele_idx],
  		   elePhi[Sec_Pair_ele_idx],0);
     v_ee = ele1 + ele2;
     double Mee = v_ee.M();

     double sT = elePt[First_Pair_ele_idx] + elePt[Sec_Pair_ele_idx] + caloJetIC5Pt[First_Pair_jet_idx] + caloJetIC5Pt[Sec_Pair_jet_idx];
     

     ///increment array value if this event passes
     for (int i=0;i<10;i++){
       if (ele1pT>ele1pTvec[i]){ //check that the first electron is above the threshold
	 for (int j=0;j<10;j++){
	   if (ele2pT>ele2pTvec[j]){
	     for (int k=0;k<10;k++){
	       if (jet1pT>jet1pTvec[k]){
		 for (int l=0;l<10;l++){
		   if (jet2pT>jet2pTvec[l]){
		     for (int m=0;m<10;m++){
		       if (Mee>MeeLowvec[m]){
			 for (int n=0;n<10;n++){
			   if (sT>sTvec[n])
			   ++opt[i][j][k][l][m][n];
			 }
		       }
		     }
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }

   } // end if TwoEles and TwoJets
    //////////////////////////Fill cut histograms///////////////
   resetCuts();
   
   if( v_idx_ele_ID_ISO.size() >= 1 )   fillVariableWithValue( "pT1stEle", elePt[First_Pair_ele_idx] );
   if( v_idx_ele_ID_ISO.size() >= 2 )   fillVariableWithValue( "pT2ndEle", elePt[Sec_Pair_ele_idx] );
   if( v_idx_jet_final.size() >= 1 )  fillVariableWithValue( "pT1stJet", PtScale * caloJetIC5Pt[First_Pair_jet_idx] );
   if( v_idx_jet_final.size() >= 2 )  fillVariableWithValue( "pT2ndJet", PtScale * caloJetIC5Pt[Sec_Pair_jet_idx] );

   // Evaluate cuts (but do not apply them)
   evaluateCuts();


    ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

     ////Fill opt histogram
     for (int i=0;i<10;i++){
	 for (int j=0;j<10;j++){
	     for (int k=0;k<10;k++){
		 for (int l=0;l<10;l++){
		     for (int m=0;m<10;m++){
			 for (int n=0;n<10;n++){
			   int bin = (100000*i)+(10000*j)+(1000*k)+(100*l)+(10*m)+n+1; 
			   h_Nevent->SetBinContent(bin,opt[i][j][k][l][m][n]);
			 }
		     }
		 }
	     }
	 }
     }

   //////////write histos 
   h_Nevent->Write();


   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
