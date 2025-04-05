#define ssb_analysis_cxx

#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "TMath.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

#include "./../interface/ssb_analysis.hpp"
#include "./../CommonTools.hpp"

using namespace std;

void ssb_analysis::SetConfig()
{
 
   // Lep Infor ( Di Muon or Di Electron )
   lepisotype = SSBConfReader->GetText( "Iso_type" );
   lepId      = SSBConfReader->GetText( "Lep_ID"   );
   n_elep_Id  = SSBConfReader->Size( "Lep_ID_cut"  );
  
   // Muon Infor for "Only MuEl channel"
   Muonisotype = SSBConfReader->GetText( "MuonIso_type" );
   Muonid      = SSBConfReader->GetText( "Muon_ID"      );

   // Electron Infor "Only for MuEl channel"
   Elecisotype = SSBConfReader->GetText( "ElecIso_type" );
   Elecid      = SSBConfReader->GetText( "Elec_ID"      );
   n_elec_Id   = SSBConfReader->Size(    "Elec_ID_cut"  );
   // Jet Infor
   jetId      = SSBConfReader->GetText("Jet_ID");
   jetbtag    = SSBConfReader->GetText("Jet_btag");


   // Kinematic variables for Object //(not MuEl)
   lep_pt     = SSBConfReader->GetNumber( "Lep_pt"  );
   lep_eta    = SSBConfReader->GetNumber( "Lep_eta" );
   lepisocut  = SSBConfReader->GetNumber( "Iso_cut" );
   jet_pt     = SSBConfReader->GetNumber( "Jet_pt"  );
   jet_eta    = SSBConfReader->GetNumber( "Jet_eta" );
   met_cut    = SSBConfReader->GetNumber( "MET_cut" );

   // Kinematic variables for Object (MuEl)
   muon_pt     = SSBConfReader->GetNumber( "Muon_pt"     );
   muon_eta    = SSBConfReader->GetNumber( "Muon_eta"    );
   muonisocut  = SSBConfReader->GetNumber( "MuonIso_cut" );
   elec_pt     = SSBConfReader->GetNumber( "Elec_pt"     );
   elec_eta    = SSBConfReader->GetNumber( "Elec_eta"    );
   elecisocut  = SSBConfReader->GetNumber( "ElecIso_cut" );
   

   // variables for Jet Cleaning 
   muiso_type = SSBConfReader->GetText( "Muiso_for_Jetcleanning" );
   muId       = SSBConfReader->GetText( "Muid_for_Jetcleanning"  );
   eliso_type = SSBConfReader->GetText( "Eliso_for_Jetcleanning" );
   elId       = SSBConfReader->GetText( "Elid_for_Jetcleanning"  );

   // Kinematic variables for Jet Cleaning 
   mu_pt_jetcl   = SSBConfReader->GetNumber( "Mupt_for_Jetcleanning"     );
   mu_eta_jetcl  = SSBConfReader->GetNumber( "Mueta_for_Jetcleanning"    );
   mu_iso_jetcl  = SSBConfReader->GetNumber( "Muisocut_for_Jetcleanning" );
   el_pt_jetcl   = SSBConfReader->GetNumber( "Elpt_for_Jetcleanning"     );
   el_eta_jetcl  = SSBConfReader->GetNumber( "Eleta_for_Jetcleanning"    );
   el_iso_jetcl  = SSBConfReader->GetNumber( "Elisocut_for_Jetcleanning" );
   n_el_id_jetcl = SSBConfReader->Size(      "Elidcut_for_Jetcleanning"  );

   if (n_el_id_jetcl == 1)
   {
      el_id_jetcl_1 = SSBConfReader->GetNumber("Elidcut_for_Jetcleanning");
      el_id_jetcl_2 = -999;
   }
   else if (n_el_id_jetcl == 2 )
   {
      el_id_jetcl_1 = SSBConfReader->GetNumber("Elidcut_for_Jetcleanning", 1);
      el_id_jetcl_2 = SSBConfReader->GetNumber("Elidcut_for_Jetcleanning", 2);
   }
   else {cout << " Electron Jetcleaning error" <<endl; }


   /////////////////////////////
   // For KinematicSolver !!! //
   /////////////////////////////

//   double nupair[] = { 30.7137, 56.2880, 23.0744, 59.1015,24.9145 }; // for 8TeV
   double nupair[] = { 1/4.53731e+06,  5.81649e+01,  2.31911e+01,  1.19038e+04, 
                       -9.78134e+01 , -2.76894e+01, -2.01114e-02, -5.07004e+00,
                         5.76608e+01,  2.34419e+01,  1.18119e+04,  5.26711e+03,
                         1.64080e+03,  5.86524e+00, -8.07173e+01, -5.44618e+00 };
   int nupair_size = sizeof(nupair)/sizeof(nupair[0]);
   
   for (int i = 0; i < nupair_size; ++i ){ nupars_.push_back(nupair[i]); }

   topmass_begin = 50.0;
   topmass_end   = 500.0;
   topmass_step_ = 1.0;


}
void ssb_analysis::GetVariables()
{

   ///////////////////////
   /// Case of Di-Muon ///
   /////////////////////// 
   if ( TString(Decaymode).Contains( "dimuon" ) ) 
   { 
      // for each decay modes, this part is getting lepton type from config file
      // 각 decay mode에서, 이 부분은 lepton 타입들을 config 파일에서 가져오는 부분임

      //////////////////////
      /// Muon condition ///
      //////////////////////

      /// Muon Iso type
      if ( TString(lepisotype).Contains( "PFIsodbeta03" ) )       { v_lepton_iso = Muon_PFIsodBeta03; }
      else if ( TString(lepisotype).Contains( "PFIsodbeta04" ) )  { v_lepton_iso = Muon_PFIsodBeta04; }
      else { cout << "Muon Condition Iso Error" << endl; } 

      /// Muon ID
      if      (TString(lepId).Contains( "Loose"  ) ) { v_lepton_Id = Muon_isLoose;  }
      else if (TString(lepId).Contains( "Tight"  ) ) { v_lepton_Id = Muon_isTight;  }
      else { cout << "Muon ID Error" << endl; } 

   }
   ///////////////////////////
   /// Case of Di-Electron ///
   ///////////////////////////
   else if (TString(Decaymode).Contains( "dielec" ) )
   {
      /// Electron iso type
      if      ( TString( lepisotype ).Contains( "PFIsoRho03"   ) ) { v_lepton_iso = Elec_PFIsoRho03;   }
      else if ( TString( lepisotype ).Contains( "PFIsoRho04"   ) ) { v_lepton_iso = Elec_PFIsoRho04;   }
      else { cout << "Electron Iso type Error at Case of Di-Electron" << endl; }
   
      // Electron ID

      /// simple cut based electron ID
      if      (TString(lepId).Contains( "Loose"  ) ) { v_lepton_Id = Elec_SCB_Loose;  }
      else if (TString(lepId).Contains( "Medium" ) ) { v_lepton_Id = Elec_SCB_Medium; }
      else if (TString(lepId).Contains( "Tight"  ) ) { v_lepton_Id = Elec_SCB_Tight;  }
      else if (TString(lepId).Contains( "Veto"   ) ) { v_lepton_Id = Elec_SCB_Veto;   }
      else { cout << "Electron ID Error" << endl; }

   } 
   ////////////////////
   /// Case of MuEl ///
   ////////////////////
   else if (TString(Decaymode).Contains( "muel" ) )
   {
      /// Muon Iso type
      if      ( TString(Muonisotype).Contains( "PFIsodbeta03" ) )  { v_mlep_iso = Muon_PFIsodBeta03; }
      else if ( TString(Muonisotype).Contains( "PFIsodbeta04" ) )  { v_mlep_iso = Muon_PFIsodBeta04; }
      else { cout << "Muon Condition Iso Error at Case of MuEl" << endl; } 

      /// Muon ID
      if      ( TString(Muonid).Contains( "Loose"  ) ) { v_muon_Id = Muon_isLoose;  }
      else if ( TString(Muonid).Contains( "Tight"  ) ) { v_muon_Id = Muon_isTight;  }
      else { cout << "Muon Condition ID Error" << endl; } 

      /// Electron iso type
      if      ( TString( Elecisotype ).Contains( "PFIsoRho03" ) )   { v_elec_iso = Elec_PFIsoRho03;   }
      else if ( TString( Elecisotype ).Contains( "PFIsoRho04" ) )   { v_elec_iso = Elec_PFIsoRho04;   }
      else { cout << "Electron Iso type Error at Case of MuEl" << endl; }
   
      // Electron ID SCB
     
      if      (TString(Elecid).Contains( "Loose"  ) ) { v_electron_Id = Elec_SCB_Loose;  }
      else if (TString(Elecid).Contains( "Medium" ) ) { v_electron_Id = Elec_SCB_Medium; }
      else if (TString(Elecid).Contains( "Tight"  ) ) { v_electron_Id = Elec_SCB_Tight;  }
      else if (TString(Elecid).Contains( "Veto"   ) ) { v_electron_Id = Elec_SCB_Veto;   }
      else { cout << "Electron ID Error" << endl; } 

   }
   else { cout << "Error Decaymode !! " << endl; }

   /////////////////////////////
   /// *** Jet condition *** ///
   /////////////////////////////

   if      ( TString(jetId).Contains( "PFLoose" ) ) { v_jet_Id = Jet_PFId; jet_id = 1; }
   else if ( TString(jetId).Contains( "PFTight" ) ) { v_jet_Id = Jet_PFId; jet_id = 2; }
   else {cout << "Jet condition error" << endl;}

   if      ( TString(jetbtag).Contains( "CSVL"  ) )    { bdisccut = 0.244; }
   else if ( TString(jetbtag).Contains( "CSVM"  ) )    { bdisccut = 0.679; }
   else if ( TString(jetbtag).Contains( "CSVT"  ) )    { bdisccut = 0.898; }
   else if ( TString(jetbtag).Contains( "CISVL" ) )    { bdisccut = 0.423; }
   else if ( TString(jetbtag).Contains( "CISVM" ) )    { bdisccut = 0.814; }
   else if ( TString(jetbtag).Contains( "CISVT" ) )    { bdisccut = 0.941; }
   else if ( TString(jetbtag).Contains( "pfCSVV2L" ) ) { bdisccut = 0.5426; }
   else if ( TString(jetbtag).Contains( "pfCSVV2M" ) ) { bdisccut = 0.8484; }
   else if ( TString(jetbtag).Contains( "pfCSVV2T" ) ) { bdisccut = 0.9535; }
   else { cout << "bscriminator error !!" << endl; }
   //////////////////////////
   // *** JetCleanning ***///
   //////////////////////////
   // This process is removal untrustful jets to reduce noise, choose trustful jet, and improve reconstruction
   // 신뢰할수없는 jet을 제거하는 과정이고 노이즈 제거, 신뢰가능한 jet 선택, 그리고 reconstruction 향상등의 목적이 있다
   // iso lepton이 jet으로 라벨링(?) 경우(즉 뮤온과 젯이 같은 객체일때) 그 이벤트는 버린다

   ///////////////////////////
   // Muon for Jet Cleaning //
   ///////////////////////////

   /// Muon iso type
   if      ( TString( muiso_type ).Contains( "PFIsodbeta03" ) ) { v_mu_iso_jcl = Muon_PFIsodBeta03; }
   else if ( TString( muiso_type ).Contains( "PFIsodbeta04" ) ) { v_mu_iso_jcl = Muon_PFIsodBeta04; }
   else { cout << "Muon Iso type error for Jet cleanning" << endl;}

   /// Muon ID
   if      ( TString( muId ).Contains( "Loose" ) ) { v_mu_Id_jcl = Muon_isLoose; }
   else if ( TString( muId ).Contains( "Tight" ) ) { v_mu_Id_jcl = Muon_isTight; }
   else { cout << "Muon ID error for Jet cleanning" << endl;}
 
   ///////////////////////////////
   // Electron for Jet Cleaning //
   ///////////////////////////////

   /// Electron iso type
   if      ( TString( eliso_type ).Contains( "PFIsoRho03"   ) ) { v_el_iso_jcl = Elec_PFIsoRho03;   }
   else if ( TString( eliso_type ).Contains( "PFIsoRho04"   ) ) { v_el_iso_jcl = Elec_PFIsoRho04;   }
   else { cout << "Electron Iso type Error for Jet Cleaning" << endl;}

   // Electron ID SCB 

   if      (TString(elId).Contains( "Loose"  ) ) { v_el_Id_jcl = Elec_SCB_Loose;  }
   else if (TString(elId).Contains( "Medium" ) ) { v_el_Id_jcl = Elec_SCB_Medium; }
   else if (TString(elId).Contains( "Tight"  ) ) { v_el_Id_jcl = Elec_SCB_Tight;  }
   else if (TString(elId).Contains( "Veto"   ) ) { v_el_Id_jcl = Elec_SCB_Veto;   }
   else { cout << "Electron ID Error" << endl; }

}

void ssb_analysis::TLVInitial()
{
   Lep1      = new TLorentzVector(); // Leading Lepton 
   Lep2      = new TLorentzVector(); // Second Leading Lepton

   TMuon     = new TLorentzVector(); // Leading Lepton 
   TElectron = new TLorentzVector(); // Second Leading Lepton

   Jet1      = new TLorentzVector(); // Leading Jet 
   Jet2      = new TLorentzVector(); // Second Leading Jet

   Jet1Up      = new TLorentzVector(); // Leading Jet 
   Jet2Up      = new TLorentzVector(); // Second Leading Jet
   Jet1Dn      = new TLorentzVector(); // Leading Jet 
   Jet2Dn      = new TLorentzVector(); // Second Leading Jet
   Jet1JERUp      = new TLorentzVector(); // Leading Jet 
   Jet2JERUp      = new TLorentzVector(); // Second Leading Jet
   Jet1JERDn      = new TLorentzVector(); // Leading Jet 
   Jet2JERDn      = new TLorentzVector(); // Second Leading Jet

   bJet1     = new TLorentzVector(); // Leading Jet 
   bJet2     = new TLorentzVector(); // Second Leading Jet
   bJet1Up     = new TLorentzVector(); // Leading Jet 
   bJet2Up     = new TLorentzVector(); // Second Leading Jet
   bJet1Dn     = new TLorentzVector(); // Leading Jet 
   bJet2Dn     = new TLorentzVector(); // Second Leading Jet

   bJet1JERUp     = new TLorentzVector(); // Leading Jet 
   bJet2JERUp     = new TLorentzVector(); // Second Leading Jet
   bJet1JERDn     = new TLorentzVector(); // Leading Jet 
   bJet2JERDn     = new TLorentzVector(); // Second Leading Jet

   Met       = new TLorentzVector(); // MET
   MetJESUp  = new TLorentzVector(); // MET
   MetJESDn  = new TLorentzVector(); // MET
   MetJERUp  = new TLorentzVector(); // MET
   MetJERDn  = new TLorentzVector(); // MET
  
   Nu1       = new TLorentzVector(); // Nutrino 1 ( For Lep1 )
   Nu2       = new TLorentzVector(); // Nutrino 2 ( For Lep2 )

   W1        = new TLorentzVector(); // W boson 1 ( For Lep1 )  
   W2        = new TLorentzVector(); // W boson 2 ( For Lep2 )

   Top       = new TLorentzVector(); // Top
   AnTop     = new TLorentzVector(); // Anti-top

   Top1      = new TLorentzVector(); // Top Leading Pt
   Top2      = new TLorentzVector(); // Top Second Leading Pt 

   bJet      = new TLorentzVector(); // bJet 
   AnbJet    = new TLorentzVector(); // Anti-bJet

   Lep       = new TLorentzVector(); // Lepton 
   AnLep     = new TLorentzVector(); // Anti-Lepton

   Nu       = new TLorentzVector(); // Lepton 
   AnNu     = new TLorentzVector(); // Anti-Lepton

   GenLep       = new TLorentzVector(); // Lepton 
   GenAnLep     = new TLorentzVector(); // Anti-Lepton

   // Lepton selection
   Mu_lep_sel  = new TLorentzVector();
   Ele_lep_sel = new TLorentzVector();

   // Veto muon and electron
   vetoMu  = new TLorentzVector();
   vetoEle = new TLorentzVector();

   // Jet
   Mu_jet_clean           = new TLorentzVector();     
   Ele_jet_clean          = new TLorentzVector();
   selected_Ele_jet_clean = new TLorentzVector();
   selected_Mu_jet_clean  = new TLorentzVector();
   TJet                   = new TLorentzVector();

   // new algorithm ( )
   highPtJet = new TLorentzVector();
   leading_bTagJet = new TLorentzVector();
   pair_bJet = new TLorentzVector();
   pair_AnbJet = new TLorentzVector();
   pair_Lep = new TLorentzVector();
   pair_AnLep = new TLorentzVector();
   // O3 version 2
   v2_pair_Jet = new TLorentzVector();
   v2_pair_AnJet = new TLorentzVector();
   v2_pair_Lep = new TLorentzVector();
   v2_pair_AnLep = new TLorentzVector();

   Genb = new TLorentzVector();
   GenAnb = new TLorentzVector();
   GenMu = new TLorentzVector();
   GenAnMu = new TLorentzVector();
}


void ssb_analysis::Loop( char *logfile )
{

   GetTotalEvent();
   //////////
   if (fChain == 0) return;
   //////////

   //////////
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   //////////

   ///My variables
   Long64_t __tot_evt = 0;

   /// Check Total Event
   cout << "Ntuple Total Event Check !! " << NtupletotalEvent << endl;

   MCSF();
   cout << "MC_SF Check !! " << mc_sf_ << endl;

   Np_eventw_2 = 0.0;
   Nm_eventw_2 = 0.0;
   Num_cpo3_p  = 0; 
   Num_cpo3_m  = 0;

   // before dupleevet is __tot_evt
   Long64_t af_dupleevt = 0;

   int count_emu = 0;
  
   ////////////////////////
   /// start event loop ///
   ////////////////////////

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
      {
         printf("ERROR: Could not load tree!!!\n");
         break;
      }

      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
      if (jentry % 10000 == 0) printf("Event %lld\n", jentry); //%lld supports Long64_t

      __tot_evt++;
      ////////////////////////////////////////
      /// start Main Loop Function and Cut ///
      ////////////////////////////////////////

      // Make TL to use SSBTree //
      MakeVecforTL();

      // initailizing TLorentzVector
      TLVInitial();

      evt_weight_ = 1;

      GetVariables();
      // Apply MC Scale Factor //
      MCSFApply();
      // Apply GenWeight //
      GenWeightApply();
      // Apply Fact & Reno //
      FactRenoApply();
      // Apply PileUpReWeight
      PDFWeightApply();
      FragmentApply();
      DecayTableApply();
      PileUpReWeightApply();
      L1PreFireApply();
      int Info_RunNumber_ =Info_RunNumber;
      int Info_Luminosity_ =Info_Luminosity;
      int Info_EventNumber_ =Info_EventNumber;
      TopPtReweightApply(); /// Apply Top pT Reweight ///
      ////////////////////////////////////////
      /// ** Di-Lepton Channel Analysis ** ///
      ////////////////////////////////////////
    
      if ( TString(Decaymode).Contains( "dielec" ) ||
           TString(Decaymode).Contains( "dimuon" ) || 
           TString(Decaymode).Contains( "muel" )     )
      {
      
         /////////////////////////////////////////////
         /// Finding Di-Lep. Channel at Gen.Level. ///
         /////////////////////////////////////////////

         if ( ChannelIndex() == false) {continue;}
         //cout << "L1_PreFire_Central : " << L1_PreFire_Central << endl;
         //cout << "L1_PreFire_Up : " << L1_PreFire_Up << endl;
         //cout << "L1_PreFire_Down : " << L1_PreFire_Down << endl;
         // Num. primary vertex counter
         NumPVCount();
         // Lepton Define //
         LeptonSelector();
         LeptonOrder();
         // Jet Define //
         JetSelector();
         // Met Define //
         METDefiner();
         BDsicApply();
         BJetDefiner(); 

         /// To check up the vertex distribution before pre-selection ///
         FillHisto(h_Num_PV_BeforePreSel, num_pv, evt_weight_);

         if (METFilterAPP() == true){
            FillHisto(h_Num_PV_AfterMetFilter, num_pv, evt_weight_);
         }
         if (Trigger() == true){
            FillHisto(h_Num_PV_AfterTrigger, num_pv, evt_weight_);
         }
         ////////////////////
         /// Event Filter ///
         ////////////////////
         if ( METFilterAPP() == false ) {continue;}
         if ( Filter_PV->at(0) == false ) {continue;}

         //////////////////////////////////
         /// trigger requirement step 0 ///
         //////////////////////////////////
         if ( Trigger() == false ) {continue;}
         // Apply Trigger SF //
         if ( FileName_.Contains("Data_Single") )
         {
            if ( !RMDuplEvt(Info_RunNumber,Info_Luminosity,Info_EventNumber) ) { continue; }
            // if문에서 RMDupleEvt이 true가 아니라면 현재 이벤트는 중복이벤트이기 때문에 스킵하게 됨
         }
         af_dupleevt++;

         FillHisto( h_EventWeight[0] , evt_weight_  );
         FillHisto( h_cf_NLeptons[0], Muon_Count, evt_weight_);
         FillHisto( h_cf_NJets[0], v_jet_idx.size(), evt_weight_ );
         FillHisto( h_cf_NPV[0]    , num_pv, evt_weight_ );
         FillHisto( h_Num_PV[0]    , num_pv, evt_weight_ );
      
         if ( NumIsoLeptons() == false ){ continue;}

         TriggerSFApply();
         ///////////////////////
         /// 3rd Lepton Veto ///
         ///////////////////////
         if ( ThirdLeptonVeto() == false ){continue;}
         if (LeptonsPtAddtional() == false ) {continue;}
         //if (evt_weight_ != 1.){cout << "evt_weight_ : " << evt_weight_ << endl;}
         // TLorentzVector define

         ////////////////////////////////
         /// DiLeptonMassCut() step 1 ///
         ////////////////////////////////
         if ( DiLeptonMassCut() == false ) {continue;} //((*Lep1)+(*Lep2)).M() > 20 ==> true
         LeptonSFApply();

         FillHisto( h_EventWeight[1] , evt_weight_  );
         FillHisto( h_cf_NPV[1]    , num_pv        , evt_weight_ );
         FillHisto( h_cf_NLeptons[1], v_lepton_idx.size(), evt_weight_ );
         FillHisto( h_cf_Lep1pt[1] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep1eta[1], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep1phi[1], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_cf_Lep2pt[1] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep2eta[1], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep2phi[1], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_cf_dilep_inv_mass[1], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_cf_NJets[1]  , v_jet_idx.size(), evt_weight_ );
         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_cf_Jet1pt[1] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet1eta[1], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet1phi[1], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_cf_Jet2pt[1] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet2eta[1], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet2phi[1], (*Jet2).Phi() , evt_weight_ );
            }
         }

/*         cout << " Info_EventNumber : " << Info_EventNumber << " step 1 evt_weight_ : " <<  v_SystEvt[1] 
         << " Lep1 Pt : " << Lep1->Pt() 
         << " Lep2 Pt : " << Lep2->Pt() 
         << endl;*/

         FillHisto( h_cf_metpt[1] ,Met->Pt()  , evt_weight_ );
         FillHisto( h_cf_metphi[1],Met->Phi() , evt_weight_ );
         FillHisto( h_DiLepMass[1], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_Num_PV[1], num_pv, evt_weight_ );
         FillHisto( h_Lep1pt[1] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_Lep1eta[1], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_Lep1phi[1], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_Lep2pt[1] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_Lep2eta[1], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_Lep2phi[1], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_METpt[1]   , Met->Pt()  , evt_weight_ );
         FillHisto( h_METphi[1]  , Met->Phi()  , evt_weight_ );
         if(v_jet_idx.size() > 0 ) //after jet requirement
         {
            FillHisto( h_Jet1pt[1] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_Jet1eta[1], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_Jet1phi[1], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_Jet2pt[1] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_Jet2eta[1], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_Jet2phi[1], (*Jet2).Phi() , evt_weight_ );
            }
         }
         if (isAllSyst == true){
            for ( int i = 0; i < v_SystFullName.size(); ++i )
            {  
               FillHisto( h_cf_sys_NLeptons[i][1], v_lepton_idx.size(), v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1pt[i][1] , (*Lep1).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1eta[i][1], (*Lep1).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1phi[i][1], (*Lep1).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2pt[i][1] , (*Lep2).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2eta[i][1], (*Lep2).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2phi[i][1], (*Lep2).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_NPV[i][1]    , num_pv        , v_SystEvt[i] );
               FillHisto( h_cf_sys_NJets[i][1]  , v_jet_idx.size(), v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1pt[i][1] , (*Jet1).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1eta[i][1], (*Jet1).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1phi[i][1], (*Jet1).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2pt[i][1] , (*Jet2).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2eta[i][1], (*Jet2).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2phi[i][1], (*Jet2).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_metpt[i][1] ,Met->Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_metphi[i][1],Met->Phi() , v_SystEvt[i] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[i][1], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
               FillHisto( h_cf_sys_metpt[i][1] , Met->Pt() , v_SystEvt[i] );
               FillHisto( h_cf_sys_metphi[i][1], Met->Phi(), v_SystEvt[i] );
               
               FillHisto( h_sys_Num_PV[i][1], num_pv, v_SystEvt[i] );
               FillHisto( h_sys_DiLepMass[i][1], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
               FillHisto( h_sys_Lep1pt[i][1] , (*Lep1).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Lep1eta[i][1], (*Lep1).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Lep1phi[i][1], (*Lep1).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Lep2pt[i][1] , (*Lep2).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Lep2eta[i][1], (*Lep2).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Lep2phi[i][1], (*Lep2).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Jet1pt[i][1] , (*Jet1).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Jet1eta[i][1], (*Jet1).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Jet1phi[i][1], (*Jet1).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Jet2pt[i][1] , (*Jet2).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Jet2eta[i][1], (*Jet2).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Jet2phi[i][1], (*Jet2).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Num_Jets[i][1], v_jet_idx.size(), v_SystEvt[i] );
               FillHisto( h_sys_METpt[i][1]   , Met->Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_METphi[i][1]  , Met->Phi()  , v_SystEvt[i] );
            }     
         }
         /////////////////////
         // ZVetocut step 2 //
         /////////////////////
         if ( ZVetoCut() == false ) {continue;} // range(76~206 GeV) veto

/*         cout << " Info_EventNumber : " << Info_EventNumber << " step 2 evt_weight_ : " <<  v_SystEvt[1] 
         << " Lep1 Pt : " << Lep1->Pt() 
         << " Lep2 Pt : " << Lep2->Pt() 
         << endl;*/

         FillHisto( h_EventWeight[2], evt_weight_  );
         FillHisto( h_cf_NLeptons[2], v_lepton_idx.size(), evt_weight_ );
         FillHisto( h_cf_Lep1pt[2] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep1eta[2], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep1phi[2], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_cf_Lep2pt[2] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep2eta[2], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep2phi[2], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_cf_NPV[2]    , num_pv        , evt_weight_ );
         FillHisto( h_cf_NJets[2]  , v_jet_idx.size(), evt_weight_ );
         FillHisto( h_cf_dilep_inv_mass[2], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_cf_metpt[2] , Met->Pt() , evt_weight_ );
         FillHisto( h_cf_metphi[2], Met->Phi(), evt_weight_ );
         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_cf_Jet1pt[2] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet1eta[2], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet1phi[2], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_cf_Jet2pt[2] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet2eta[2], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet2phi[2], (*Jet2).Phi() , evt_weight_ );
            }
         }

         FillHisto( h_Num_PV[2], num_pv, evt_weight_ );
         FillHisto( h_DiLepMass[2], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_Lep1pt[2] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_Lep1eta[2], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_Lep1phi[2], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_Lep2pt[2] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_Lep2eta[2], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_Lep2phi[2], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_Num_Jets[2]  , v_jet_idx.size(), evt_weight_ );
         FillHisto( h_METpt[2]   , Met->Pt()  , evt_weight_ );
         FillHisto( h_METphi[2]  , Met->Phi()  , evt_weight_ );
         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_Jet1pt[2] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_Jet1eta[2], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_Jet1phi[2], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_Jet2pt[2] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_Jet2eta[2], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_Jet2phi[2], (*Jet2).Phi() , evt_weight_ );
            }
         }
         //FillHisto( h_cf_met[2], , evt_weight_ );
         //FillHisto( h_cf_Nbjets[2], , evt_weight_ );
         if (isAllSyst == true){
            for ( int i = 0; i < v_SystFullName.size(); ++i )
            {  
               FillHisto( h_cf_sys_NLeptons[i][2], v_lepton_idx.size(), v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1pt[i][2] , (*Lep1).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1eta[i][2], (*Lep1).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1phi[i][2], (*Lep1).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2pt[i][2] , (*Lep2).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2eta[i][2], (*Lep2).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2phi[i][2], (*Lep2).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_NPV[i][2]    , num_pv        , v_SystEvt[i] );
               FillHisto( h_cf_sys_NJets[i][2]  , v_jet_idx.size(), v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1pt[i][2] , (*Jet1).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1eta[i][2], (*Jet1).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1phi[i][2], (*Jet1).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2pt[i][2] , (*Jet2).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2eta[i][2], (*Jet2).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2phi[i][2], (*Jet2).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_metpt[i][2] ,Met->Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_metphi[i][2],Met->Phi() , v_SystEvt[i] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[i][2], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
               FillHisto( h_cf_sys_metpt[i][2] , Met->Pt() , v_SystEvt[i] );
               FillHisto( h_cf_sys_metphi[i][2], Met->Phi(), v_SystEvt[i] );
               
               FillHisto( h_sys_Num_PV[i][2], num_pv, v_SystEvt[i] );
               FillHisto( h_sys_DiLepMass[i][2], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
               FillHisto( h_sys_Lep1pt[i][2] , (*Lep1).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Lep1eta[i][2], (*Lep1).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Lep1phi[i][2], (*Lep1).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Lep2pt[i][2] , (*Lep2).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Lep2eta[i][2], (*Lep2).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Lep2phi[i][2], (*Lep2).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Jet1pt[i][2] , (*Jet1).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Jet1eta[i][2], (*Jet1).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Jet1phi[i][2], (*Jet1).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Jet2pt[i][2] , (*Jet2).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Jet2eta[i][2], (*Jet2).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Jet2phi[i][2], (*Jet2).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Num_Jets[i][2], v_jet_idx.size(), v_SystEvt[i] );
               FillHisto( h_sys_METpt[i][2]   , Met->Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_METphi[i][2]  , Met->Phi()  , v_SystEvt[i] );
            }
         }
         /////////////////////////      
         // Num. Jet cut step 3 //
         /////////////////////////
         JetDefiner();
         //if ( NumJetCut(v_jet_idx) == false ) {continue;}
         if ( NumJetCut(v_jet_idx) == true )// if number of jet is over 2, true
         {
/*            cout  << " Info_EventNumber : " << Info_EventNumber << " step 3 evt_weight_ : " << v_SystEvt[1] 
                  << " Lep1 Pt : " << (*Lep1).Pt()
                  << " Lep2 Pt : " << (*Lep2).Pt()
                  << " Jet1 Pt : " << (*Jet1).Pt()
                  << " Jet2 Pt : " << (*Jet2).Pt()
                  << " MET : " << Met->Pt()
                  << endl;*/
            // Define Leading Jet and Second Leading Jet 
            FillHisto( h_EventWeight[3] , evt_weight_  );
            FillHisto( h_cf_NLeptons[3], v_lepton_idx.size(), evt_weight_ );
            FillHisto( h_cf_Lep1pt[3] , (*Lep1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Lep1eta[3], (*Lep1).Eta() , evt_weight_ );
            FillHisto( h_cf_Lep1phi[3], (*Lep1).Phi() , evt_weight_ );
            FillHisto( h_cf_Lep2pt[3] , (*Lep2).Pt()  , evt_weight_ );
            FillHisto( h_cf_Lep2eta[3], (*Lep2).Eta() , evt_weight_ );
            FillHisto( h_cf_Lep2phi[3], (*Lep2).Phi() , evt_weight_ );
            FillHisto( h_cf_NPV[3]    , num_pv        , evt_weight_ );
            FillHisto( h_cf_NJets[3]  , v_jet_idx.size(), evt_weight_ );
            FillHisto( h_cf_Jet1pt[3] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet1eta[3], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet1phi[3], (*Jet1).Phi() , evt_weight_ );
            FillHisto( h_cf_Jet2pt[3] , (*Jet2).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet2eta[3], (*Jet2).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet2phi[3], (*Jet2).Phi() , evt_weight_ );
            FillHisto( h_cf_metpt[3] ,Met->Pt()  , evt_weight_ );
            FillHisto( h_cf_metphi[3],Met->Phi() , evt_weight_ );
            
            FillHisto( h_cf_dilep_inv_mass[3], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
            FillHisto( h_cf_metpt[3] , Met->Pt() , evt_weight_ );
            FillHisto( h_cf_metphi[3], Met->Phi(), evt_weight_ );
            
            FillHisto( h_Num_PV[3], num_pv, evt_weight_ );
            FillHisto( h_DiLepMass[3], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
            FillHisto( h_Lep1pt[3] , (*Lep1).Pt()  , evt_weight_ );
            FillHisto( h_Lep1eta[3], (*Lep1).Eta() , evt_weight_ );
            FillHisto( h_Lep1phi[3], (*Lep1).Phi() , evt_weight_ );
            FillHisto( h_Lep2pt[3] , (*Lep2).Pt()  , evt_weight_ );
            FillHisto( h_Lep2eta[3], (*Lep2).Eta() , evt_weight_ );
            FillHisto( h_Lep2phi[3], (*Lep2).Phi() , evt_weight_ );
            FillHisto( h_Jet1pt[3] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_Jet1eta[3], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_Jet1phi[3], (*Jet1).Phi() , evt_weight_ );
            FillHisto( h_Jet2pt[3] , (*Jet2).Pt()  , evt_weight_ );
            FillHisto( h_Jet2eta[3], (*Jet2).Eta() , evt_weight_ );
            FillHisto( h_Jet2phi[3], (*Jet2).Phi() , evt_weight_ );
            FillHisto( h_Num_Jets[3], v_jet_idx.size(), evt_weight_ );
            FillHisto( h_METpt[3]   , Met->Pt()  , evt_weight_ );
            FillHisto( h_METphi[3]  , Met->Phi()  , evt_weight_ );

            if (isAllSyst == true){
               for ( int i = 0; i < v_SystFullName.size(); ++i )
               {  
                  if (TString(v_SystFullName[i]).Contains("Jet")  ){continue;}
                  FillHisto( h_cf_sys_NLeptons[i][3], v_lepton_idx.size(), v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep1pt[i][3] , (*Lep1).Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep1eta[i][3], (*Lep1).Eta() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep1phi[i][3], (*Lep1).Phi() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep2pt[i][3] , (*Lep2).Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep2eta[i][3], (*Lep2).Eta() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep2phi[i][3], (*Lep2).Phi() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_NPV[i][3]    , num_pv        , v_SystEvt[i] );
                  FillHisto( h_cf_sys_NJets[i][3]  , v_jet_idx.size(), v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet1pt[i][3] , (*Jet1).Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet1eta[i][3], (*Jet1).Eta() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet1phi[i][3], (*Jet1).Phi() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet2pt[i][3] , (*Jet2).Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet2eta[i][3], (*Jet2).Eta() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet2phi[i][3], (*Jet2).Phi() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_metpt[i][3] ,Met->Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_metphi[i][3],Met->Phi() , v_SystEvt[i] );
                  
                  FillHisto( h_cf_sys_dilep_inv_mass[i][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                  FillHisto( h_cf_sys_metpt[i][3] , Met->Pt() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_metphi[i][3], Met->Phi(), v_SystEvt[i] );
                  
                  FillHisto( h_sys_Num_PV[i][3], num_pv, v_SystEvt[i] );
                  FillHisto( h_sys_DiLepMass[i][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                  FillHisto( h_sys_Lep1pt[i][3] , (*Lep1).Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_Lep1eta[i][3], (*Lep1).Eta() , v_SystEvt[i] );
                  FillHisto( h_sys_Lep1phi[i][3], (*Lep1).Phi() , v_SystEvt[i] );
                  FillHisto( h_sys_Lep2pt[i][3] , (*Lep2).Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_Lep2eta[i][3], (*Lep2).Eta() , v_SystEvt[i] );
                  FillHisto( h_sys_Lep2phi[i][3], (*Lep2).Phi() , v_SystEvt[i] );
                  FillHisto( h_sys_Jet1pt[i][3] , (*Jet1).Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_Jet1eta[i][3], (*Jet1).Eta() , v_SystEvt[i] );
                  FillHisto( h_sys_Jet1phi[i][3], (*Jet1).Phi() , v_SystEvt[i] );
                  FillHisto( h_sys_Jet2pt[i][3] , (*Jet2).Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_Jet2eta[i][3], (*Jet2).Eta() , v_SystEvt[i] );
                  FillHisto( h_sys_Jet2phi[i][3], (*Jet2).Phi() , v_SystEvt[i] );
                  FillHisto( h_sys_Num_Jets[i][3], v_jet_idx.size(), v_SystEvt[i] );
                  FillHisto( h_sys_METpt[i][3]   , Met->Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_METphi[i][3]  , Met->Phi()  , v_SystEvt[i] );
               }
            }
            //////////////////// 
            // MET cut step 4 //
            ////////////////////
            //if ( METCut(Met) == false ) {continue;}
            if ( METCut(Met) == true ) { //MET > 40 GeV
/*            cout  << " Info_EventNumber : " << Info_EventNumber << " step 4 evt_weight_ : " << v_SystEvt[1] 
                  << " Lep1 Pt : " << (*Lep1).Pt()
                  << " Lep2 Pt : " << (*Lep2).Pt()
                  << " Jet1 Pt : " << (*Jet1).Pt()
                  << " Jet2 Pt : " << (*Jet2).Pt()
                  << " MET : " << Met->Pt()
                  << endl; */
               FillHisto( h_EventWeight[4] , evt_weight_  );
               FillHisto( h_cf_NLeptons[4], v_lepton_idx.size(), evt_weight_ );
               
               FillHisto( h_cf_Lep1pt[4] , (*Lep1).Pt()  , evt_weight_ );
               FillHisto( h_cf_Lep1eta[4], (*Lep1).Eta() , evt_weight_ );
               FillHisto( h_cf_Lep1phi[4], (*Lep1).Phi() , evt_weight_ );
               FillHisto( h_cf_Lep2pt[4] , (*Lep2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Lep2eta[4], (*Lep2).Eta() , evt_weight_ );
               FillHisto( h_cf_Lep2phi[4], (*Lep2).Phi() , evt_weight_ );
               FillHisto( h_cf_NPV[4]    , num_pv        , evt_weight_ );
               FillHisto( h_cf_NJets[4]  , v_jet_idx.size(), evt_weight_ );
               FillHisto( h_cf_Jet1pt[4] , (*Jet1).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet1eta[4], (*Jet1).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet1phi[4], (*Jet1).Phi() , evt_weight_ );
               FillHisto( h_cf_Jet2pt[4] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet2eta[4], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet2phi[4], (*Jet2).Phi() , evt_weight_ );
               
               FillHisto( h_cf_dilep_inv_mass[4], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
               FillHisto( h_cf_metpt[4] , Met->Pt() , evt_weight_ );
               FillHisto( h_cf_metphi[4], Met->Phi(), evt_weight_ );
               
               FillHisto( h_Lep1pt[4]  , Lep1->Pt() , evt_weight_ );
               FillHisto( h_Lep2pt[4]  , Lep2->Pt() , evt_weight_ );
               FillHisto( h_Lep1eta[4] , Lep1->Eta(), evt_weight_ );
               FillHisto( h_Lep2eta[4] , Lep2->Eta(), evt_weight_ );
               FillHisto( h_Lep1phi[4] , Lep1->Phi(), evt_weight_ );
               FillHisto( h_Lep2phi[4] , Lep2->Phi(), evt_weight_ );
               
               FillHisto( h_Jet1pt[4]  , Jet1->Pt() , evt_weight_ );
               FillHisto( h_Jet2pt[4]  , Jet2->Pt() , evt_weight_ );
               FillHisto( h_Jet1eta[4] , Jet1->Eta(), evt_weight_ );
               FillHisto( h_Jet2eta[4] , Jet2->Eta(), evt_weight_ );
               FillHisto( h_Jet1phi[4] , Jet1->Phi(), evt_weight_ );
               FillHisto( h_Jet2phi[4] , Jet2->Phi(), evt_weight_ );
               FillHisto( h_METpt[4]   , Met->Pt()  , evt_weight_ );
               FillHisto( h_METphi[4]  , Met->Phi() , evt_weight_ );
               
               FillHisto( h_DiLepMass[4], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
               FillHisto( h_Num_PV[4], num_pv, evt_weight_ );
               FillHisto( h_Num_Jets[4], v_jet_idx.size(), evt_weight_ );
               if (isAllSyst == true)
               {
                  for ( int i = 0; i < v_SystFullName.size(); ++i )
                  {  
                     if (TString(v_SystFullName[i]).Contains("Jet")){continue;}
                     FillHisto( h_cf_sys_NLeptons[i][4], v_lepton_idx.size(), v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep1pt[i][4] , (*Lep1).Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep1eta[i][4], (*Lep1).Eta() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep1phi[i][4], (*Lep1).Phi() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep2pt[i][4] , (*Lep2).Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep2eta[i][4], (*Lep2).Eta() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep2phi[i][4], (*Lep2).Phi() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_NPV[i][4]    , num_pv        , v_SystEvt[i] );
                     FillHisto( h_cf_sys_NJets[i][4]  , v_jet_idx.size(), v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet1pt[i][4] , (*Jet1).Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet1eta[i][4], (*Jet1).Eta() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet1phi[i][4], (*Jet1).Phi() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet2pt[i][4] , (*Jet2).Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet2eta[i][4], (*Jet2).Eta() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet2phi[i][4], (*Jet2).Phi() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_metpt[i][4] ,Met->Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_metphi[i][4],Met->Phi() , v_SystEvt[i] );
                     
                     FillHisto( h_cf_sys_dilep_inv_mass[i][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                     FillHisto( h_cf_sys_metpt[i][4] , Met->Pt() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_metphi[i][4], Met->Phi(), v_SystEvt[i] );
                     
                     FillHisto( h_sys_Num_PV[i][4], num_pv, v_SystEvt[i] );
                     FillHisto( h_sys_DiLepMass[i][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                     FillHisto( h_sys_Lep1pt[i][4] , (*Lep1).Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_Lep1eta[i][4], (*Lep1).Eta() , v_SystEvt[i] );
                     FillHisto( h_sys_Lep1phi[i][4], (*Lep1).Phi() , v_SystEvt[i] );
                     FillHisto( h_sys_Lep2pt[i][4] , (*Lep2).Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_Lep2eta[i][4], (*Lep2).Eta() , v_SystEvt[i] );
                     FillHisto( h_sys_Lep2phi[i][4], (*Lep2).Phi() , v_SystEvt[i] );
                     FillHisto( h_sys_Jet1pt[i][4] , (*Jet1).Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_Jet1eta[i][4], (*Jet1).Eta() , v_SystEvt[i] );
                     FillHisto( h_sys_Jet1phi[i][4], (*Jet1).Phi() , v_SystEvt[i] );
                     FillHisto( h_sys_Jet2pt[i][4] , (*Jet2).Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_Jet2eta[i][4], (*Jet2).Eta() , v_SystEvt[i] );
                     FillHisto( h_sys_Jet2phi[i][4], (*Jet2).Phi() , v_SystEvt[i] );
                     FillHisto( h_sys_Num_Jets[i][4], v_jet_idx.size(), v_SystEvt[i] );
                     FillHisto( h_sys_METpt[i][4]   , Met->Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_METphi[i][4]  , Met->Phi()  , v_SystEvt[i] );
                  }
               }
               double AllJetpt = 0;
               for ( int ijet = 0; ijet < v_jet_idx.size(); ++ijet )
               { 
//                  TLorentzVector *htJet = (TLorentzVector*)Jet->At( v_jet_idx[ijet] );
                  TLorentzVector *htJet = v_jet_TL[ijet];
                  AllJetpt += htJet->Pt(); 
               }
               
               FillHisto( h_HT[4], AllJetpt, evt_weight_);
               /////////////////////////////////////////
               /// One or more b-Tagging Requirement ///
               /////////////////////////////////////////
               if (JetPtPhiDil) { ApplyJetPtPhiDilution();}
               
               BTaggigSFApply();
               //if ( BJetCut(v_bjet_idx) == false ) {continue;} // one or more b-tagging 
               if ( BJetCut(v_bjet_idx) == true ) // one or more b-tagging (Ref(Line: 5477):if (v_bjets.size() >= 1) { bjetcut = true; })
               { 
/*            cout  << " Info_EventNumber : " << Info_EventNumber << " step 5 evt_weight_ : " << v_SystEvt[1] 
                  << " Lep1 Pt : " << (*Lep1).Pt()
                  << " Lep2 Pt : " << (*Lep2).Pt()
                  << " Jet1 Pt : " << (*Jet1).Pt()
                  << " Jet2 Pt : " << (*Jet2).Pt()
                  << " bJet1 Pt : " << (*bJet1).Pt()
                  << " bJet2 Pt : " << (*bJet2).Pt()
                  << " MET : " << Met->Pt()
                  << endl; */

                  FillHisto(h_bTagWeight,evt_weight_/evt_weight_beforeBtag_); 
                  FillHisto( h_EventWeight[5] , evt_weight_  );
                  FillHisto( h_cf_NLeptons[5], v_lepton_idx.size(), evt_weight_ );
                  
                  FillHisto( h_cf_Lep1pt[5] , (*Lep1).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Lep1eta[5], (*Lep1).Eta() , evt_weight_ );
                  FillHisto( h_cf_Lep1phi[5], (*Lep1).Phi() , evt_weight_ );
                  FillHisto( h_cf_Lep2pt[5] , (*Lep2).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Lep2eta[5], (*Lep2).Eta() , evt_weight_ );
                  FillHisto( h_cf_Lep2phi[5], (*Lep2).Phi() , evt_weight_ );
                  FillHisto( h_cf_NPV[5]    , num_pv        , evt_weight_ );
                  FillHisto( h_cf_NJets[5]  , v_jet_idx.size(), evt_weight_ );
                  FillHisto( h_cf_Jet1pt[5] , (*Jet1).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Jet1eta[5], (*Jet1).Eta() , evt_weight_ );
                  FillHisto( h_cf_Jet1phi[5], (*Jet1).Phi() , evt_weight_ );
                  FillHisto( h_cf_Jet2pt[5] , (*Jet2).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Jet2eta[5], (*Jet2).Eta() , evt_weight_ );
                  FillHisto( h_cf_Jet2phi[5], (*Jet2).Phi() , evt_weight_ );
                  
                  FillHisto( h_cf_dilep_inv_mass[5], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                  FillHisto( h_cf_metpt[5] , Met->Pt() , evt_weight_ );
                  FillHisto( h_cf_metphi[5], Met->Phi(), evt_weight_ );
                  
                  
                  

                  
                  FillHisto( h_Lep1pt[5]  , Lep1->Pt() , evt_weight_ );
                  FillHisto( h_Lep2pt[5]  , Lep2->Pt() , evt_weight_ );
                  FillHisto( h_Lep1eta[5] , Lep1->Eta(), evt_weight_ );
                  FillHisto( h_Lep2eta[5] , Lep2->Eta(), evt_weight_ );
                  FillHisto( h_Lep1phi[5] , Lep1->Phi(), evt_weight_ );
                  FillHisto( h_Lep2phi[5] , Lep2->Phi(), evt_weight_ );
                  
                  if (TString(Decaymode).Contains("muel"))
                  { 
                     FillHisto( h_Muonpt[5]  , TMuon->Pt()      , evt_weight_ );
                     FillHisto( h_Elecpt[5]  , TElectron->Pt()  , evt_weight_ );
                     FillHisto( h_Muoneta[5] , TMuon->Eta()     , evt_weight_ );
                     FillHisto( h_Eleceta[5] , TElectron->Eta() , evt_weight_ );
                     FillHisto( h_Muonphi[5] , TMuon->Phi()     , evt_weight_ );
                     FillHisto( h_Elecphi[5] , TElectron->Phi() , evt_weight_ );
                  }
                  
                  FillHisto( h_Jet1pt[5]  , Jet1->Pt() , evt_weight_ );
                  FillHisto( h_Jet2pt[5]  , Jet2->Pt() , evt_weight_ );
                  FillHisto( h_Jet1eta[5] , Jet1->Eta(), evt_weight_ );
                  FillHisto( h_Jet2eta[5] , Jet2->Eta(), evt_weight_ );
                  FillHisto( h_Jet1phi[5] , Jet1->Phi(), evt_weight_ );
                  FillHisto( h_Jet2phi[5] , Jet2->Phi(), evt_weight_ );
                  FillHisto( h_METpt[5]   , Met->Pt()  , evt_weight_ );
                  FillHisto( h_METphi[5]  , Met->Phi() , evt_weight_ );
                  
                  
                  FillHisto( h_DiLepMass[5], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                  
                  FillHisto( h_Num_PV[5], num_pv, evt_weight_ );
                  FillHisto( h_Num_Jets[5], v_jet_idx.size(), evt_weight_ );
                  FillHisto( h_Num_bJets[5], nbtagged, evt_weight_ );
                  
                  FillHisto( h_HT[5], AllJetpt, evt_weight_);
//////////////////////////////////////////////////////////////////// new algorithm start ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  //       << New algorithm version 1 >>                                                                           //
                  //       - use 1 highest pT b-tagging jet among b-tagging jets, leading jet                                      //
                  //       - calcuate dR with b-tagging jet & lepton, Anti-lepton                                                  //
                  //       - if dR btw b-tagging jet & Lepton, b-tagging jet becomes anti-b-tagging jet vise versa                 //
                  //       - the rest is, for example, Anti-lepton, it will be paired with leading jet                             //
                  //       - however, in pT of leading jet and b-tagging jet is the same case, subleading jet will be paired       //
                  //       - get O3 with different dR requirements (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)                                  //
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  ///////// find highest pT b tagging jet & pair b tagging jet(highest pT) and lepton /////////
                  // 1) find highest pT b tagging jet
                  // (v_jet_TL: 전체 jet, v_bjet_TL: b-tagging jet)
                  double maxPt_bJet = 0;
                  int bTagJet_total = 0;
                  int bTagJet_is_leadingJet = 0;
                  int bTagJet_idx = -1;
                  //cout << "=============== b tagged jet & jet idx check ===============" << endl;
                  for (int ibjet = 0; ibjet < v_bjet_TL.size(); ++ibjet)
                  {
                     if (v_bjet_TL[ibjet]->Pt() > maxPt_bJet)
                     {
                        maxPt_bJet = v_bjet_TL[ibjet]->Pt();
                        leading_bTagJet = v_bjet_TL[ibjet];
                        bTagJet_idx = ibjet;
                     }
                     //cout << "bjet idx: " << ibjet << " pT: " << v_bjet_TL[ibjet]->Pt() << endl;
                  }
                  for(int ijet = 0; ijet < v_jet_TL.size(); ijet++)
                  {
                     //cout << "jet idx: " << ijet << " jet pT: " << v_jet_TL[ijet]->Pt() << endl;
                  }
                  //cout << "leading_bTagJet_idx: " << bTagJet_idx << ", leading_bTagJet: " << leading_bTagJet->Pt() << endl;
                  //cout << "leading jet pT: " << Jet1->Pt() << ", subleading jet pT: " << Jet2->Pt() << endl;
                  // fraction for leading b tagging jet of whole b tagging jet
                  if (v_bjet_TL.size() > 0) 
                  {
                     bTagJet_total = v_bjet_TL.size();
                     // leading b-tagging jet has the same pT with leading jet(Jet1)
                     if (leading_bTagJet->Pt() == Jet1->Pt()) {bTagJet_is_leadingJet = 1;}
                     else {bTagJet_is_leadingJet = 0;}
                  }
                  FillHisto(h_frac_bTagJet_total            , bTagJet_total            ,evt_weight_);
                  FillHisto(h_frac_LeadingbTagJet_is_Leading, bTagJet_is_leadingJet    ,evt_weight_);
                  // highest pT b tagging jet is leading b tagging jet
                  // Jet2 is subleading jet
                  FillHisto(h_leading_bTagJetpt    ,leading_bTagJet->Pt()     ,evt_weight_);
                  FillHisto(h_leading_bTagJeteta   ,leading_bTagJet->Eta()    ,evt_weight_);
                  FillHisto(h_leading_bTagJetphi   ,leading_bTagJet->Phi()    ,evt_weight_);
                  FillHisto(h_leading_jetpt        ,Jet1->Pt()                ,evt_weight_);
                  FillHisto(h_leading_jeteta       ,Jet1->Eta()               ,evt_weight_);
                  FillHisto(h_leading_jetphi       ,Jet1->Phi()               ,evt_weight_);
                  FillHisto(h_subleading_jetpt     ,Jet2->Pt()                ,evt_weight_);
                  FillHisto(h_subleading_jeteta    ,Jet2->Eta()               ,evt_weight_);
                  FillHisto(h_subleading_jetphi    ,Jet2->Phi()               ,evt_weight_);
                  // 2) paring between b tagging jet & lepton, anti lepton
                  double dR_bTagJet_Lep = 999;
                  double dR_bTagJet_AnLep = 999;
                  double dR_pair_bJet_AnLep = 999;
                  double dR_pair_AnbJet_Lep = 999;
                  dR_bTagJet_Lep = leading_bTagJet->DeltaR(*((TLorentzVector*)Lep)); //Lep has minus charge
                  dR_bTagJet_AnLep = leading_bTagJet->DeltaR(*((TLorentzVector*)AnLep)); //AnLep has plus charge
                  //cout << "=============== dR check ================" << endl;
                  //cout << "dR_bTagJet_Lep: " << dR_bTagJet_Lep << ", dR_bTagJet_AnLep: " << dR_bTagJet_AnLep << endl;
                  if(dR_bTagJet_Lep < dR_bTagJet_AnLep) //bTagJet is close to Lep, it means bTagJet has minus charge, that is AnbJet
                  {
                     pair_AnbJet = leading_bTagJet;
                     pair_Lep = Lep;
                     pair_AnLep = AnLep;
                     // if bTagJet's pT and jet_leading's pT are the same, select subleading jet as a AnbJet paired with Lep2
                     if(leading_bTagJet->Pt() != Jet1->Pt()) {pair_bJet = Jet1;}
                     else {pair_bJet = Jet2;}

                     dR_pair_bJet_AnLep = pair_bJet->DeltaR(*((TLorentzVector*)pair_AnLep));
                     dR_pair_AnbJet_Lep = pair_AnbJet->DeltaR(*((TLorentzVector*)pair_Lep));
                     //cout << "dR_bTagJet_Lep < dR_bTagJet_AnLep: leading b tagging jet is close to Lep" << endl;
                     //cout << "leading_bTagJet becomes pair_AnbJet" << endl;
                     //cout << "dR: " << dR_bTagJet_Lep << ", pair_AnbJet: " << pair_AnbJet->Pt() << ", pair_bJet: " << pair_bJet->Pt() << ", pair_Lep: " << pair_Lep->Pt() << ", pair_AnLep: " << pair_AnLep->Pt() << endl;
                  }
                  else //bTagJet is close to AnLep
                  {
                     pair_bJet = leading_bTagJet;
                     pair_Lep = Lep;
                     pair_AnLep = AnLep;
                     if(leading_bTagJet->Pt() != Jet1->Pt()) {pair_AnbJet = Jet1;}
                     else {pair_AnbJet = Jet2;}

                     dR_pair_bJet_AnLep = pair_bJet->DeltaR(*((TLorentzVector*)pair_AnLep));
                     dR_pair_AnbJet_Lep = pair_AnbJet->DeltaR(*((TLorentzVector*)pair_Lep));
                     //cout << "dR_bTagJet_Lep > dR_bTagJet_AnLep: leading b tagging jet is close to AnLep" << endl;
                     //cout << "leading_bTagJet becomes pair_bJet" << endl;
                     //cout << "dR: " << dR_bTagJet_AnLep << ", pair_AnbJet: " << pair_AnbJet->Pt() << ", pair_bJet: " << pair_bJet->Pt() << ", pair_Lep: " << pair_Lep->Pt() << ", pair_AnLep: " << pair_AnLep->Pt() << endl;
                  }
                  // delta R distibution: before pairing
                  FillHisto(h_dR_bTagJet_Lep    ,dR_bTagJet_Lep                     ,evt_weight_);
                  FillHisto(h_dR_bTagJet_AnLep  ,dR_bTagJet_AnLep                   ,evt_weight_);
                  FillHisto(h_dR_bJet_AnbJet    ,dR_bTagJet_Lep  ,dR_bTagJet_AnLep  ,evt_weight_);
                  // delta R distribution: After pairing
                  FillHisto(h_dR_pair_bJet_AnLep   ,dR_pair_bJet_AnLep                       ,evt_weight_);   //dR distribution b quark & lepton
                  FillHisto(h_dR_pair_AnbJet_Lep   ,dR_pair_AnbJet_Lep                       ,evt_weight_);   //dR distribution b-bar quark & anti-lepton
                  FillHisto(h_dR_pair_bJet_AnbJet  ,dR_pair_AnbJet_Lep ,dR_pair_bJet_AnLep   ,evt_weight_);

                  vector<double> v_dR_value = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 10.0}; //0, 1, 2, 3, 4, 5, 6
                  // 10.0 is the same as no dR requirement
                  if(pair_bJet && pair_AnbJet && pair_Lep && pair_AnLep)
                  {
                     for(int idr=0; idr < v_dR_value.size(); idr++)
                     {
                        if(dR_pair_bJet_AnLep < v_dR_value[idr] && dR_pair_AnbJet_Lep < v_dR_value[idr]) // bJet-Lep, AnbJet-AnLep간의 거리가 v_dR_value[idr](0.5~3.0) 이하일때
                        {
                           FillHisto(h_dRNum_pair_bJet_AnLep[idr]  ,dR_pair_bJet_AnLep  ,evt_weight_);
                           FillHisto(h_dRNum_pair_AnbJet_Lep[idr]  ,dR_pair_AnbJet_Lep  ,evt_weight_);
                           FillHisto(h_dRNum_pair_bJet_AnbJet[idr] ,dR_pair_AnbJet_Lep  ,dR_pair_bJet_AnLep  ,evt_weight_);
                           FillHisto(h_Pair_bJetpt[idr]     ,pair_bJet->Pt()    ,evt_weight_);
                           FillHisto(h_Pair_bJeteta[idr]    ,pair_bJet->Eta()   ,evt_weight_);
                           FillHisto(h_Pair_bJetphi[idr]    ,pair_bJet->Phi()   ,evt_weight_);
                           FillHisto(h_Pair_AnbJetpt[idr]   ,pair_AnbJet->Pt()  ,evt_weight_);
                           FillHisto(h_Pair_AnbJeteta[idr]  ,pair_AnbJet->Eta() ,evt_weight_);
                           FillHisto(h_Pair_AnbJetphi[idr]  ,pair_AnbJet->Phi() ,evt_weight_);
                           FillHisto(h_pair_Leppt[idr]      ,pair_Lep->Pt()     ,evt_weight_);
                           FillHisto(h_pair_Lepeta[idr]     ,pair_Lep->Eta()    ,evt_weight_);
                           FillHisto(h_pair_Lepphi[idr]     ,pair_Lep->Phi()    ,evt_weight_);
                           FillHisto(h_pair_AnLeppt[idr]    ,pair_AnLep->Pt()   ,evt_weight_);
                           FillHisto(h_pair_AnLepeta[idr]   ,pair_AnLep->Eta()  ,evt_weight_);
                           FillHisto(h_pair_AnLepphi[idr]   ,pair_AnLep->Phi()  ,evt_weight_);                                       
                           //invariant mass distribution
                           // 1) (An)bJet - (An)Mu
                           FillHisto(h_pair_b_Anb_Mass[idr]   ,( (*pair_bJet)+(*pair_AnbJet) ).M()  ,evt_weight_ );
                           FillHisto(h_pair_Mu_AnMu_Mass[idr] ,( (*pair_Lep)+(*pair_AnLep) ).M()    ,evt_weight_ );
                           // 2) bJet - AnbJet , Mu - AnMu
                           FillHisto(h_pair_b_AnMu_Mass[idr]  ,( (*pair_bJet)+(*pair_AnLep) ).M()   ,evt_weight_ );
                           FillHisto(h_pair_Anb_Mu_Mass[idr]  ,( (*pair_AnbJet)+(*pair_Lep) ).M()   ,evt_weight_ );
                           //cout << "Reco O3 v1 dR: " << v_dR_value[idr] << endl;
                           //cout << "bJet: " << (*pair_bJet).M() << ", AnbJet: " << (*pair_AnbJet).M() << ", Lep: " << (*pair_Lep).M() << ", AnLep: " << (*pair_AnLep).M() << endl;
                           //cout << "bJet-AnLep: " << ( (*pair_bJet)+(*pair_AnLep) ).M() << ", AnbJet-Lep: " << ( (*pair_AnbJet)+(*pair_Lep) ).M() << endl;
                           // calculate O3
                           double n_cpO3 = ssbcpviol->getO3Vari( pair_bJet , pair_AnbJet, pair_AnLep, pair_Lep );
                           if(n_cpO3 == 0.0 ){cout << "--00000000000000--" << n_cpO3 << endl;}
                           // Calculating O3 error.
                           if(n_cpO3 > 0.0) {FillHisto(h_CPO3_Pos[idr]  ,n_cpO3  ,evt_weight_);}
                           else             {FillHisto(h_CPO3_Neg[idr]  ,n_cpO3  ,evt_weight_);}
                           FillHisto(h_CPO3_bfReco[idr]        ,n_cpO3  ,evt_weight_);
                           FillHisto(h_CPO3_bfReco_ReRange[idr],n_cpO3  ,evt_weight_);
                        }
                     }
                  }
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  //       << New algorithm version 2 >>                                                                           //
                  //       - use all jets                                                                                          //
                  //       - calcuate dR with jet & (anti-)lepton then, save the values                                            //
                  //       - find jet idx that has minimum dR with (anti-)lepton                                                   //
                  //       - if they have different index, pair them                                                               //
                  //       - if indices are the same, follow substeps                                                              //
                  //          - sub 1) dR w/ anti-lepton is smaller than dR w/ lepton                                              //
                  //             - pair the jet with anti-lepton as pair_bJet                                                      //
                  //             - find second minimum dR value and pair them with lepton                                          //
                  //          - sub 2) dR w/ lepton is smaller than dR w/ anti-lepton                                              //
                  //             - pair the jet with lepton as pair_AnbJet                                                         //
                  //             - find second minimum dR value and pair them with anti-lepton                                     //
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  // 1) All jets; before pairing   //
                  ///////////////////////////////////
                  for(int ijet = 0; ijet < v_jet_TL.size(); ijet++)
                  {
                     FillHisto(h_AllJet_idx  ,ijet                     ,evt_weight_);
                     FillHisto(h_AllJet_pt   ,v_jet_TL[ijet]->Pt()     ,evt_weight_);
                     FillHisto(h_AllJet_eta  ,v_jet_TL[ijet]->Eta()    ,evt_weight_);
                     FillHisto(h_AllJet_phi  ,v_jet_TL[ijet]->Phi()    ,evt_weight_);
                  }
                  //-----------------------------------------------------------//
                  // 2) finding Jet1 & Jet2 which have minimum dR with AnLep & Lep
                  //-----------------------------------------------------------//
                  // dR 계산 및 매칭 로직 구현
                  // 제트와 Lep, AnLep 간의 dR 값을 제트 인덱스 순서대로 저장
                  dR_jet_Lep.clear();
                  dR_jet_AnLep.clear();
                  //cout << " =============== dR with Lep & AnLep =============== " << endl;
                  for(int ijet = 0; ijet < v_jet_TL.size(); ijet++)
                  {
                     double dR_Lep = v_jet_TL[ijet]->DeltaR(*((TLorentzVector*)Lep));
                     double dR_AnLep = v_jet_TL[ijet]->DeltaR(*((TLorentzVector*)AnLep));
                     dR_jet_Lep.push_back(dR_Lep);
                     dR_jet_AnLep.push_back(dR_AnLep);
                     // all of dR distribution
                     FillHisto(h_v2_dR_Lep          ,dR_Lep            ,evt_weight_);
                     FillHisto(h_v2_dR_AnLep        ,dR_AnLep          ,evt_weight_);
                     FillHisto(h_v2_dR_Lep_AnLep    ,dR_Lep            ,dR_AnLep         ,evt_weight_);
                     //cout << "dR_jet_Lep[" << ijet << "]: " << dR_Lep << ", dR_jet_AnLep[" << ijet << "]: " << dR_AnLep << endl;
                  }
                  // 최소 dR 값과 해당 제트 인덱스 찾기
                  auto min_iter_jet_Lep = std::min_element(dR_jet_Lep.begin(), dR_jet_Lep.end());
                  int min_index_jet_Lep = std::distance(dR_jet_Lep.begin(), min_iter_jet_Lep);
                  double min_dR_jet_Lep = *min_iter_jet_Lep;
                  auto min_iter_jet_AnLep = std::min_element(dR_jet_AnLep.begin(), dR_jet_AnLep.end());
                  int min_index_jet_AnLep = std::distance(dR_jet_AnLep.begin(), min_iter_jet_AnLep);
                  double min_dR_jet_AnLep = *min_iter_jet_AnLep;
                  // 매칭 로직
                  if(min_index_jet_Lep != min_index_jet_AnLep)
                  {
                     // 경우 1: 서로 다른 제트가 최소 dR을 가짐
                     v2_pair_AnJet = v_jet_TL[min_index_jet_Lep]; // Lep과 매칭된 제트를 AnJet으로
                     v2_pair_Jet = v_jet_TL[min_index_jet_AnLep]; // AnLep과 매칭된 제트를 Jet으로
                     v2_pair_Lep = Lep;
                     v2_pair_AnLep = AnLep;
                  }
                  else
                  {
                     // 경우 2: 같은 제트가 Lep와 AnLep에 대한 최소 dR을 가짐
                     if(min_dR_jet_Lep < min_dR_jet_AnLep)
                     {
                        // 서브 경우 (1): Lep와의 dR이 더 작음
                        v2_pair_AnJet = v_jet_TL[min_index_jet_Lep];
                        v2_pair_Lep = Lep;
                        v2_pair_AnLep = AnLep;
                        // 두 번째로 작은 dR_AnLep 찾기
                        double second_min_dR_jet_AnLep = std::numeric_limits<double>::max();
                        int second_min_index_jet_AnLep = -1;
                        for(int ijet = 0; ijet < v_jet_TL.size(); ijet++)
                        {
                           if(ijet == min_index_jet_AnLep) continue; // 이미 매칭된 제트 제외
                           if(dR_jet_AnLep[ijet] < second_min_dR_jet_AnLep)
                           {
                              second_min_dR_jet_AnLep = dR_jet_AnLep[ijet];
                              second_min_index_jet_AnLep = ijet;
                           }
                        }
                        if(second_min_index_jet_AnLep != -1)
                        {
                           v2_pair_Jet = v_jet_TL[second_min_index_jet_AnLep];
                           min_dR_jet_AnLep = second_min_dR_jet_AnLep; //update min_dR_jet_AnLep
                           min_index_jet_AnLep = second_min_index_jet_AnLep; //update min_index_jet_AnLep
                        }
                        else
                        {
                           // 두 번째 매칭 제트가 없는 경우 처리
                           v2_pair_Jet = nullptr;
                           cout << "can't find second min dR" << endl;
                        }
                     }
                     else
                     {
                        // 서브 경우 (2): AnLep과의 dR이 더 작음
                        v2_pair_Jet = v_jet_TL[min_index_jet_AnLep];
                        v2_pair_AnLep = AnLep;
                        v2_pair_Lep = Lep;
                        // 두 번째로 작은 dR_Lep 찾기
                        double second_min_dR_jet_Lep = std::numeric_limits<double>::max();
                        int second_min_index_jet_Lep = -1;
                        for(int ijet = 0; ijet < v_jet_TL.size(); ijet++)
                        {
                           if(ijet == min_index_jet_Lep) continue; // 이미 매칭된 제트 제외
                           if(dR_jet_Lep[ijet] < second_min_dR_jet_Lep)
                           {
                              second_min_dR_jet_Lep = dR_jet_Lep[ijet];
                              second_min_index_jet_Lep = ijet;
                           }
                        }
                        if(second_min_index_jet_Lep != -1)
                        {
                           v2_pair_AnJet = v_jet_TL[second_min_index_jet_Lep];
                           min_dR_jet_Lep = second_min_dR_jet_Lep; //update min_dR_jet_Lep
                           min_index_jet_Lep = second_min_index_jet_Lep; //update min_index_jet_Lep
                        }
                        else
                        {
                           // 두 번째 매칭 제트가 없는 경우 처리
                           v2_pair_AnJet = nullptr;
                           cout << "can't find second min dR" << endl;
                        }
                     }
                  }
                  //cout << "v2_pair_Jet_idx: " << min_index_jet_AnLep << ", dR: " << min_dR_jet_AnLep << ", pt: " << v2_pair_Jet->Pt() << ", v2_pair_AnJet_idx: " << min_index_jet_Lep << ", dR: " << min_dR_jet_Lep << ", pt: " << v2_pair_AnJet->Pt() << endl;
                  // minimum dR distribution
                  FillHisto(h_v2_dR_Min_Lep        ,min_dR_jet_Lep         ,evt_weight_);
                  FillHisto(h_v2_dR_Min_AnLep      ,min_dR_jet_AnLep       ,evt_weight_);
                  FillHisto(h_v2_dR_Min_Lep_AnLep  ,min_dR_jet_Lep         ,min_dR_jet_AnLep     ,evt_weight_);
                  ////////////////////
                  // after pairing  //
                  ////////////////////------------------------------------------------------------------------------//
                  // 3) Jet1 & Jet2 which have minimum dR with AnLep & Lep are matched with pair_bJet & pair_AnbJet
                  //------------------------------------------------------------------------------------------------//
                  // Jet & AnLep
                  if(min_index_jet_AnLep != -1 && v2_pair_Jet != nullptr && v2_pair_AnLep)
                  {
                     // check dR values final_dR_jet_AnLep is the same as min_dR_jet_AnLep and index is the same as min_index_jet_AnLep
                     double final_dR_jet_AnLep = v2_pair_Jet->DeltaR((*((TLorentzVector*)AnLep)));
                     //if(min_dR_jet_AnLep != final_dR_jet_AnLep)
                     //{cout << "Jet & AnLep: after pairing: min_index_jet_AnLep = " << min_index_jet_AnLep << ", min_dR_jet_AnLep = " << min_dR_jet_AnLep << ", final_dR_jet_AnLep = " << final_dR_jet_AnLep << ", Jet1 pt = " << v2_pair_Jet->Pt() << endl;}
                     FillHisto(h_v2_pair_Jet_dR_AnLep        ,final_dR_jet_AnLep    ,evt_weight_);
                     FillHisto(h_v2_pair_Jet_idx             ,min_index_jet_AnLep    ,evt_weight_);

                     FillHisto(h_v2_pair_Jetpt              ,v2_pair_Jet->Pt()      ,evt_weight_);
                     FillHisto(h_v2_pair_Jeteta             ,v2_pair_Jet->Eta()     ,evt_weight_);
                     FillHisto(h_v2_pair_Jetphi             ,v2_pair_Jet->Phi()     ,evt_weight_);
                     FillHisto(h_v2_pair_AnLeppt            ,v2_pair_AnLep->Pt()    ,evt_weight_);
                     FillHisto(h_v2_pair_AnLepeta           ,v2_pair_AnLep->Eta()   ,evt_weight_);
                     FillHisto(h_v2_pair_AnLepphi           ,v2_pair_AnLep->Phi()   ,evt_weight_);
                  }
                  // AnJet & Lep
                  if(min_index_jet_Lep != -1 && v2_pair_AnJet != nullptr && v2_pair_Lep)
                  {
                     v2_pair_AnJet = v_jet_TL[min_index_jet_Lep];
                     v2_pair_Lep = Lep;
                     // check dR values final_dR_jet_Lep is the same as min_dR_jet_Lep and index is the same as min_index_jet_Lep
                     double final_dR_jet_Lep = v2_pair_AnJet->DeltaR((*((TLorentzVector*)Lep)));
                     //if(min_dR_jet_Lep != final_dR_jet_Lep)
                     //{cout << "AnJet & Lep: after pairing: min_index_jet_Lep = " << min_index_jet_Lep << ", min_dR_jet_Lep = " << min_dR_jet_Lep << ", final_dR_jet_Lep = " << final_dR_jet_Lep << ", Jet2 pt = " << v2_pair_AnJet->Pt() << endl;}
                     FillHisto(h_v2_pair_AnJet_dR_Lep       ,final_dR_jet_Lep      ,evt_weight_);
                     FillHisto(h_v2_pair_AnJet_idx          ,min_index_jet_Lep      ,evt_weight_);

                     FillHisto(h_v2_pair_AnJetpt            ,v2_pair_AnJet->Pt()    ,evt_weight_);
                     FillHisto(h_v2_pair_AnJeteta           ,v2_pair_AnJet->Eta()   ,evt_weight_);
                     FillHisto(h_v2_pair_AnJetphi           ,v2_pair_AnJet->Phi()   ,evt_weight_);
                     FillHisto(h_v2_pair_Leppt              ,v2_pair_Lep->Pt()      ,evt_weight_);
                     FillHisto(h_v2_pair_Lepeta             ,v2_pair_Lep->Eta()     ,evt_weight_);
                     FillHisto(h_v2_pair_Lepphi             ,v2_pair_Lep->Phi()     ,evt_weight_);
                  }
                  // index &invariant mass & O3 calculation
                  if( (min_index_jet_AnLep != -1 && v2_pair_Jet != nullptr && v2_pair_AnLep) ||
                      (min_index_jet_Lep != -1 && v2_pair_AnJet != nullptr && v2_pair_Lep) )
                  {
                     // distribution for index of Jet and AnJet
                     FillHisto(h_v2_pair_Jet_AnJet_idx         ,min_index_jet_Lep     ,min_index_jet_AnLep  ,evt_weight_);
                     FillHisto(h_v2_pair_comp_Jet_pairbJet     ,v2_pair_Jet->Pt()     ,pair_bJet->Pt()      ,evt_weight_); //v2_pair_Jet & pair_bJet are matched AnMu
                     FillHisto(h_v2_pair_comp_AnJet_pairAnbJet ,v2_pair_AnJet->Pt()   ,pair_AnbJet->Pt()    ,evt_weight_); //v2_pair_AnJet & pair_AnbJet are matched Mu
                     // invariant mass
                     FillHisto(h_v2_pair_Jet_AnMu_Mass         ,( (*v2_pair_Jet)+(*v2_pair_AnLep) ).M()     ,evt_weight_ );
                     FillHisto(h_v2_pair_AnJet_Mu_Mass         ,( (*v2_pair_AnJet)+(*v2_pair_Lep) ).M()     ,evt_weight_ );
                     FillHisto(h_v2_pair_Jet_AnbJet_Mass   ,( (*v2_pair_Jet)+(*v2_pair_AnJet) ).M()     ,( (*pair_bJet)+(*pair_AnbJet) ).M()   ,evt_weight_ );
                     //cout << "Reco O3 v2" << endl;
                     //cout << "Jet1: " << (*v2_pair_Jet).M() << ", Jet2: " << (*v2_pair_AnJet).M() << ", Lep: " << (*v2_pair_Lep).M() << ", AnLep: " << (*v2_pair_AnLep).M() << endl;
                     //cout << "Jet1-AnLep: " << ( (*v2_pair_Jet)+(*v2_pair_AnLep) ).M() << ", Jet2-Lep: " << ( (*v2_pair_AnJet)+(*v2_pair_Lep) ).M() << endl;
                     //cout << "Jet1-Jet2: " << ( (*v2_pair_Jet)+(*v2_pair_AnJet) ).M() << ", pair_bJet-pair_AnbJet: " << ( (*pair_bJet)+(*pair_AnbJet) ).M() << endl;
                     // calcualte O3
                     double n_cpO3_v2 = ssbcpviol->getO3Vari( v2_pair_Jet , v2_pair_AnJet, v2_pair_AnLep, v2_pair_Lep );
                     if(n_cpO3_v2 == 0.0 ){cout << "O3 version2 --00000000000000--" << n_cpO3_v2 << endl;}
                     // Calculating O3 error.
                     if(n_cpO3_v2 > 0.0)  {FillHisto(h_v2_CPO3_Pos  ,n_cpO3_v2  ,evt_weight_);}
                     else                 {FillHisto(h_v2_CPO3_Neg  ,n_cpO3_v2  ,evt_weight_);}
                     FillHisto(h_v2_CPO3_bfReco          ,n_cpO3_v2  ,evt_weight_);
                     FillHisto(h_v2_CPO3_bfReco_ReRange  ,n_cpO3_v2  ,evt_weight_);
                  }
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  //       << Generator level >>                                                                                   //
                  //       - find b quarks and muon came from mother particle (top)                                                //
                  //       - because they(b & muon) has the same mother particle, regarded as they are already paired              //
                  //       - with the same way, get O3 when dR is under 3                                                          //
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  //SetGenLepAnLep();
                  int mom1_bidx = -99;
                  int mom2_bidx = -99;
                  int mu_mom1_widx = -99;
                  int mu_mom2_widx = -99;
                  int mu_mom1_topidx = -99;
                  int mu_mom2_topidx = -99;
                  double Gen_dR_Genb_GenAnLep = 999.0;
                  double Gen_dR_GenAnb_GenLep = 999.0;
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  //                               start finding mother particle of b quark and muon                               //
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  //    muon : 13         //    anti-muon : -13   //    electorn : 11  //    positrion : -11                       //
                  //    W- boson : -24    //    W+ boson : 24                                                                      //
                  //    top : 6           //    anti-top : -6                                                                      //
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  //cout << "\n event number: " << jentry << endl;
                  for(int ipdg=0; ipdg<GenPar_pdgId->size(); ipdg++) // find mother particle at generator level
                  {
                     //cout << "GenPar_pdgId loop: " << ipdg << endl;
                     if(GenPar_pdgId->at(ipdg) == 13) // gen lv muon
                     {
                        mu_mom1_widx = GenPar_Mom1_Idx->at(ipdg);
                        for(int iw=0; iw<GenPar_Idx->size(); iw++) //find W-
                        {
                           if(mu_mom1_widx == GenPar_Idx->at(iw))
                           {
                              mu_mom1_topidx = GenPar_Mom1_Idx->at(iw);
                              for(int itop=0; itop<GenPar_Idx->size(); itop++) //find anti-top
                              {
                                 if(mu_mom1_topidx == GenPar_Idx->at(itop))
                                 {
                                    GenMu = (TLorentzVector*)GenPar->At(ipdg);
                                    //cout << "muon?: " << GenPar_pdgId->at(ipdg) << ", W- boson?: " << GenPar_pdgId->at(iw) << ", anti-top?: " << GenPar_pdgId->at(itop) << endl;
                                    FillHisto(h_GenMuon_pt     ,GenMu->Pt()     ,evt_weight_);
                                    FillHisto(h_GenMuon_eta    ,GenMu->Eta()    ,evt_weight_);
                                    FillHisto(h_GenMuon_phi    ,GenMu->Phi()    ,evt_weight_);
                                 }//end top
                              }//end top loop
                           }//end w boson
                        }//endl w boson loop
                     }//end muon
                     if(GenPar_pdgId->at(ipdg) == -13) // gen lv anti-muon
                     {
                        mu_mom2_widx = GenPar_Mom1_Idx->at(ipdg);
                        for(int iw=0; iw<GenPar_Idx->size(); iw++) //find W+
                        {
                           if(mu_mom2_widx == GenPar_Idx->at(iw))
                           {
                              mu_mom2_topidx = GenPar_Mom1_Idx->at(iw);
                              for(int itop=0; itop<GenPar_Idx->size(); itop++) //find anti-top
                              {
                                 if(mu_mom2_topidx == GenPar_Idx->at(itop))
                                 {
                                    GenAnMu = (TLorentzVector*)GenPar->At(ipdg);
                                    //cout << "anti muon?: " << GenPar_pdgId->at(ipdg) << ", W+ boson?: " << GenPar_pdgId->at(iw) << ", top?: " << GenPar_pdgId->at(itop) << endl;
                                    FillHisto(h_GenAnMuon_pt   ,GenAnMu->Pt()   ,evt_weight_);
                                    FillHisto(h_GenAnMuon_eta  ,GenAnMu->Eta()  ,evt_weight_);
                                    FillHisto(h_GenAnMuon_phi  ,GenAnMu->Phi()  ,evt_weight_);
                                 }//end top
                              }//end top loop
                           }//end w boson
                        }//endl w boson loop
                     }//end muon
                     if(GenPar_pdgId->at(ipdg) == 5) // gen lv b quark
                     {
                        mom1_bidx = GenPar_Mom1_Idx->at(ipdg); // obtain the index of mother particle of pdgid 5 that is b jet
                        for(int itop=0; itop<GenPar_Idx->size(); itop++) // compare mom_idx with GenPar_Idx
                        {
                           if(mom1_bidx == GenPar_Idx->at(itop)) // if value of index of both is same, check pdgid whether it is top or not using iteration number
                           {
                              //cout << "=================== b =================" << endl;
                              //cout << "GenPar_pdgId: " << GenPar_pdgId->at(ipdg) << ", ipdg: " << ipdg << ", itop: " << itop << ", mom1 idx: " << mom1_bidx << ", Genpar_Idx: " << GenPar_Idx->at(itop) << endl;
                              //cout << "mother1 particle(GenPar_pdgId): " << GenPar_pdgId->at(itop) << endl;
                              //cout << "b?: " << GenPar_pdgId->at(ipdg) << ", top?: " << GenPar_pdgId->at(itop) << endl;
                              Genb = (TLorentzVector*)GenPar->At(ipdg); // at the moment, b quark came from top is named Gen lv bjet
                              FillHisto(h_Genb_pt     ,Genb->Pt()    ,evt_weight_);
                              FillHisto(h_Genb_eta    ,Genb->Eta()   ,evt_weight_);
                              FillHisto(h_Genb_phi    ,Genb->Phi()   ,evt_weight_);
                           }
                        }
                     }
                     if(GenPar_pdgId->at(ipdg) == -5) // gen lv anti bjet
                     {
                        mom2_bidx = GenPar_Mom1_Idx->at(ipdg);
                        for(int itop=0; itop<GenPar_Idx->size(); itop++)
                        {
                           if(mom2_bidx == GenPar_Idx->at(itop))
                           {
                              //cout << "=================== anti b =================" << endl;
                              //cout << "GenPar_pdgId: " << GenPar_pdgId->at(ipdg) << ", ipdg: " << ipdg << ", itop: " << itop << ", mom2 idx: " << mom2_bidx << ", Genpar_Idx: " << GenPar_Idx->at(itop) << endl;
                              //cout << "mother2 particle(GenPar_pdgId): " << GenPar_pdgId->at(itop) << endl;
                              //cout << "anti-b?: " << GenPar_pdgId->at(ipdg) << ", anti-top?: " << GenPar_pdgId->at(itop) << endl;
                              GenAnb = (TLorentzVector*)GenPar->At(ipdg);
                              FillHisto(h_GenAnb_pt   ,GenAnb->Pt()   ,evt_weight_);
                              FillHisto(h_GenAnb_eta  ,GenAnb->Eta()  ,evt_weight_);
                              FillHisto(h_GenAnb_phi  ,GenAnb->Phi()  ,evt_weight_);
                           }
                        }
                     }
                  }//end finding mother particle
                  if(GenMu && GenAnMu && Genb && GenAnb)
                  {
                     Gen_dR_Genb_GenAnLep = Genb->DeltaR(*((TLorentzVector*)GenAnMu));
                     Gen_dR_GenAnb_GenLep = GenAnb->DeltaR(*((TLorentzVector*)GenMu));
                     FillHisto(h_Gen_dR_Genb_GenAnLep    ,Gen_dR_Genb_GenAnLep   ,evt_weight_);
                     FillHisto(h_Gen_dR_GenAnb_GenLep    ,Gen_dR_GenAnb_GenLep   ,evt_weight_);
                     FillHisto(h_Gen_dR_Genb_GenAnb      ,Gen_dR_GenAnb_GenLep  ,Gen_dR_Genb_GenAnLep    ,evt_weight_);
                     if(Gen_dR_Genb_GenAnLep < 3 && Gen_dR_GenAnb_GenLep < 3)
                     {
                        FillHisto(h_Gen_dR3_pair_Genb_GenAnLep ,Gen_dR_Genb_GenAnLep   ,evt_weight_);
                        FillHisto(h_Gen_dR3_pair_GenAnb_GenLep ,Gen_dR_GenAnb_GenLep   ,evt_weight_);
                        FillHisto(h_Gen_dR3_af_Genb_GenAnb     ,Gen_dR_GenAnb_GenLep ,Gen_dR_Genb_GenAnLep  ,evt_weight_);
                        FillHisto(h_Pair_Genbpt                ,Genb->Pt()           ,evt_weight_);
                        FillHisto(h_Pair_Genbeta               ,Genb->Eta()          ,evt_weight_);
                        FillHisto(h_Pair_Genbphi               ,Genb->Phi()          ,evt_weight_);
                        FillHisto(h_Pair_GenAnbpt              ,GenAnb->Pt()         ,evt_weight_);
                        FillHisto(h_Pair_GenAnbeta             ,GenAnb->Eta()        ,evt_weight_);
                        FillHisto(h_Pair_GenAnbphi             ,GenAnb->Phi()        ,evt_weight_);
                        FillHisto(h_pair_GenMupt               ,GenMu->Pt()          ,evt_weight_);
                        FillHisto(h_pair_GenMueta              ,GenMu->Eta()         ,evt_weight_);
                        FillHisto(h_pair_GenMuphi              ,GenMu->Phi()         ,evt_weight_);
                        FillHisto(h_pair_GenAnMupt             ,GenAnMu->Pt()        ,evt_weight_);
                        FillHisto(h_pair_GenAnMueta            ,GenAnMu->Eta()       ,evt_weight_);
                        FillHisto(h_pair_GenAnMuphi            ,GenAnMu->Phi()       ,evt_weight_);
                        
                        //invariant mass distribution
                        // 1) (An)bJet - (An)Mu
                        FillHisto(h_pair_Gen_b_Anb_Mass  , ( (*Genb)+(*GenAnb) ).M()   , evt_weight_ );
                        FillHisto(h_pair_Gen_Mu_AnMu_Mass, ( (*GenMu)+(*GenAnMu) ).M() , evt_weight_ );
                        // 2) bJet - AnbJet , Mu - AnMu
                        FillHisto(h_pair_Gen_b_AnMu_Mass , ( (*Genb)+(*GenAnMu) ).M()  , evt_weight_ );
                        FillHisto(h_pair_Gen_Anb_Mu_Mass , ( (*GenAnb)+(*GenMu) ).M()  , evt_weight_ );
                        //cout << "invariant mass: b-b: " << ( (*Genb)+(*GenAnb) ).M() << ", mu-mu: " << ( (*GenMu)+(*GenAnMu) ).M() << ", b-anmu: " << ( (*Genb)+(*GenAnMu) ).M() << ", anb-mu: " << ( (*GenAnb)+(*GenMu) ).M() << endl;
                        // 3) b & anti-b, Mu & anti-mu mass
                        FillHisto(h_pair_Gen_b_Mass       ,(*Genb).M()       ,evt_weight_ );
                        FillHisto(h_pair_Gen_Anb_Mass     ,(*GenAnb).M()     ,evt_weight_ );
                        FillHisto(h_pair_Gen_Mu_Mass      ,(*GenMu).M()      ,evt_weight_ );
                        FillHisto(h_pair_Gen_AnMu_Mass    ,(*GenAnMu).M()    ,evt_weight_ );
                        //cout << "Gen" << endl;
                        //cout << "b: " << (*Genb).M() << ", anti-b: " << (*GenAnb).M() << ", mu: " << (*GenMu).M() << ", anti-mu: " << (*GenAnMu).M() << endl;
                        //cout << "b-anti-mu: " << ( (*Genb)+(*GenAnMu) ).M() << ", anti-b-mu: " << ( (*GenAnb)+(*GenMu) ).M() << endl;
                        // calculate O3
                        double n_GencpO3 = ssbcpviol->getO3Vari( Genb , GenAnb, GenAnMu, GenMu );
                        if (n_GencpO3 == 0.0 ) {cout << "--GEN: 00000000000000--: " << n_GencpO3 << endl;}
                        // Calculating O3 error.
                        if (n_GencpO3 > 0.0) {FillHisto(h_GenCPO3_Pos   ,n_GencpO3    ,evt_weight_);}
                        else                 {FillHisto(h_GenCPO3_Neg   ,n_GencpO3    ,evt_weight_);}
                        FillHisto(h_GenCPO3_bfReco          ,n_GencpO3 ,evt_weight_);
                        FillHisto(h_GenCPO3_bfReco_ReRange  ,n_GencpO3 ,evt_weight_);
                     }
                  }
//////////////////////////////////////////////////////////////////// end new algorithm ////////////////////////////////////////////////////////////////////////////////////////////////////////////

                  if (isAllSyst == true) {
                     for ( int i = 0; i < v_SystFullName.size(); ++i )
                     {  
                        if (TString(v_SystFullName[i]).Contains("Jet")){continue;}
                        FillHisto( h_cf_sys_NLeptons[i][5], v_lepton_idx.size(), v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep1pt[i][5] , (*Lep1).Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep1eta[i][5], (*Lep1).Eta() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep1phi[i][5], (*Lep1).Phi() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep2pt[i][5] , (*Lep2).Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep2eta[i][5], (*Lep2).Eta() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep2phi[i][5], (*Lep2).Phi() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_NPV[i][5]    , num_pv        , v_SystEvt[i] );
                        FillHisto( h_cf_sys_NJets[i][5]  , v_jet_idx.size(), v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet1pt[i][5] , (*Jet1).Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet1eta[i][5], (*Jet1).Eta() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet1phi[i][5], (*Jet1).Phi() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet2pt[i][5] , (*Jet2).Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet2eta[i][5], (*Jet2).Eta() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet2phi[i][5], (*Jet2).Phi() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_metpt[i][5] ,Met->Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_metphi[i][5],Met->Phi() , v_SystEvt[i] );
                        
                        FillHisto( h_cf_sys_dilep_inv_mass[i][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                        FillHisto( h_cf_sys_metpt[i][5] , Met->Pt() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_metphi[i][5], Met->Phi(), v_SystEvt[i] );
                        
                        FillHisto( h_sys_Num_PV[i][5], num_pv, v_SystEvt[i] );
                        FillHisto( h_sys_DiLepMass[i][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                        FillHisto( h_sys_Lep1pt[i][5] , (*Lep1).Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_Lep1eta[i][5], (*Lep1).Eta() , v_SystEvt[i] );
                        FillHisto( h_sys_Lep1phi[i][5], (*Lep1).Phi() , v_SystEvt[i] );
                        FillHisto( h_sys_Lep2pt[i][5] , (*Lep2).Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_Lep2eta[i][5], (*Lep2).Eta() , v_SystEvt[i] );
                        FillHisto( h_sys_Lep2phi[i][5], (*Lep2).Phi() , v_SystEvt[i] );
                        FillHisto( h_sys_Jet1pt[i][5] , (*Jet1).Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_Jet1eta[i][5], (*Jet1).Eta() , v_SystEvt[i] );
                        FillHisto( h_sys_Jet1phi[i][5], (*Jet1).Phi() , v_SystEvt[i] );
                        FillHisto( h_sys_Jet2pt[i][5] , (*Jet2).Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_Jet2eta[i][5], (*Jet2).Eta() , v_SystEvt[i] );
                        FillHisto( h_sys_Jet2phi[i][5], (*Jet2).Phi() , v_SystEvt[i] );
                        FillHisto( h_sys_Num_Jets[i][5], v_jet_idx.size(), v_SystEvt[i] );
                        FillHisto( h_sys_METpt[i][5]   , Met->Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_METphi[i][5]  , Met->Phi()  , v_SystEvt[i] );
                     }                 
                  }
                  TVector3 lepvec( Lep1->Px() + Lep2->Px(), Lep1->Py() + Lep2->Py(), Lep1->Pz() + Lep2->Pz() );
                  TMatrixD matrix(3,3);
                  
                  double  vec[3];
                  double  lep2 = lepvec.Mag2();
                  double  _apla = 0;
                  double  _sphe = 0;
                  double  _plan = 0;
                  
                  for ( int i=0; i<3; i++)
                  {
                     for ( int j=0; j<3; j++ )
                     {
                        matrix(i,j) = lepvec(i)*lepvec(j);
                        double norm = lep2;
                        for ( int k=0; k<v_jet_idx.size(); k++ )
                        {
                           //TLorentzVector *topoJet = (TLorentzVector*)Jet->At( v_jet_idx[k] );
                           TLorentzVector *topoJet =  v_jet_TL[k];
                           vec[0] = topoJet->Px();
                           vec[1] = topoJet->Py();
                           vec[2] = topoJet->Pz();
                           matrix(i,j) += vec[i]*vec[j];
                           norm        += vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                        }
                        if (norm > 0) matrix(i,j) /= norm;
                     }
                  }
                  if ( v_jet_idx.size() > 0 )
                  {
                     TMatrixDEigen evmatrix(matrix);
                     TVectorD eigenv=evmatrix.GetEigenValuesRe();
                     matrix.EigenVectors(eigenv);
                  
                     _apla  = 1.5*eigenv(2);
                     _sphe  = 1.5*(eigenv(2) + eigenv(1));
                     _plan  = eigenv(1) - eigenv(2);
                  }
                  
                  FillHisto(h_Topo_Apla[0] , _apla, evt_weight_);
                  FillHisto(h_Topo_Sphe[0] , _sphe, evt_weight_);
                  FillHisto(h_Topo_Plan[0] , _plan, evt_weight_);
                  
                  ///////////////////////////////////////
                  /// 2 or More B-Tagging Requriement ///
                  /////////////////////////////////////// 
                  
                  if ( DoubleBtag(v_bjet_idx) == true ) 
                  {
/*            cout  << " Info_EventNumber : " << Info_EventNumber << " step 6 evt_weight_ : " << v_SystEvt[1] 
                  << " Lep1 Pt : " << (*Lep1).Pt()
                  << " Lep2 Pt : " << (*Lep2).Pt()
                  << " Jet1 Pt : " << (*Jet1).Pt()
                  << " Jet2 Pt : " << (*Jet2).Pt()
                  << " bJet1 Pt : " << (*bJet1).Pt()
                  << " bJet2 Pt : " << (*bJet2).Pt()
                  << " MET : " << Met->Pt()
                  << endl;*/
                     FillHisto( h_EventWeight[6], evt_weight_  );
                     FillHisto( h_Lep1pt[6] , Lep1->Pt() , evt_weight_ );
                     FillHisto( h_Lep2pt[6] , Lep2->Pt() , evt_weight_ );
                     FillHisto( h_Lep1eta[6], Lep1->Eta(), evt_weight_ );
                     FillHisto( h_Lep2eta[6], Lep2->Eta(), evt_weight_ );
                     FillHisto( h_Lep1phi[6], Lep1->Phi(), evt_weight_ );
                     FillHisto( h_Lep2phi[6], Lep2->Phi(), evt_weight_ );
                     
                     if (TString(Decaymode).Contains("muel"))
                     { 
                        FillHisto( h_Muonpt[6]  , TMuon->Pt()      , evt_weight_ );
                        FillHisto( h_Elecpt[6]  , TElectron->Pt()  , evt_weight_ );
                        FillHisto( h_Muoneta[6] , TMuon->Eta()     , evt_weight_ );
                        FillHisto( h_Eleceta[6] , TElectron->Eta() , evt_weight_ );
                        FillHisto( h_Muonphi[6] , TMuon->Phi()     , evt_weight_ );
                        FillHisto( h_Elecphi[6] , TElectron->Phi() , evt_weight_ );
                     }
                     FillHisto( h_Jet1pt[6] , Jet1->Pt() , evt_weight_ );
                     FillHisto( h_Jet2pt[6] , Jet2->Pt() , evt_weight_ );
                     FillHisto( h_Jet1eta[6], Jet1->Eta(), evt_weight_ );
                     FillHisto( h_Jet2eta[6], Jet2->Eta(), evt_weight_ );
                     FillHisto( h_Jet1phi[6], Jet1->Phi(), evt_weight_ );
                     FillHisto( h_Jet2phi[6], Jet2->Phi(), evt_weight_ );
                     FillHisto( h_METpt[6]  , Met->Pt()  , evt_weight_ );
                     FillHisto( h_METphi[6] , Met->Phi() , evt_weight_ );
                     FillHisto( h_HT[6] , AllJetpt, evt_weight_);
                   
                     FillHisto( h_DiLepMass[6], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                     FillHisto( h_Num_PV[6]   , num_pv          , evt_weight_ );
                     FillHisto( h_Num_Jets[6] , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_Num_bJets[6], nbtagged        , evt_weight_ );
                   
                     FillHisto( h_cf_NLeptons[6], v_lepton_idx.size(), evt_weight_ );
                     FillHisto( h_cf_Lep1pt[6] , (*Lep1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep1eta[6], (*Lep1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep1phi[6], (*Lep1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Lep2pt[6] , (*Lep2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep2eta[6], (*Lep2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep2phi[6], (*Lep2).Phi() , evt_weight_ );
                     FillHisto( h_cf_NPV[6]    , num_pv        , evt_weight_ );
                     FillHisto( h_cf_NJets[6]  , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_cf_Jet1pt[6] , (*Jet1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet1eta[6], (*Jet1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet1phi[6], (*Jet1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Jet2pt[6] , (*Jet2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet2eta[6], (*Jet2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet2phi[6], (*Jet2).Phi() , evt_weight_ );
                   
                     FillHisto( h_cf_dilep_inv_mass[6], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                     FillHisto( h_cf_metpt[6] , Met->Pt() , evt_weight_ );
                     FillHisto( h_cf_metphi[6], Met->Phi(), evt_weight_ );
                   
                   
                     FillHisto( h_Topo_Apla[1] , _apla, evt_weight_);
                     FillHisto( h_Topo_Sphe[1] , _sphe, evt_weight_);
                     FillHisto( h_Topo_Plan[1] , _plan, evt_weight_);
                     if (isAllSyst == true) {
                        for ( int i = 0; i < v_SystFullName.size(); ++i )
                        {  
                           if (TString(v_SystFullName[i]).Contains("Jet")){continue;}
//                           if ( TString(v_SystFullName[i]).Contains("PileUpUp") ) { cout << "evt_weight_ : " << v_SystEvt[i] << " Jet1 Pt : " << Jet1->Pt() << endl;}
                           FillHisto( h_cf_sys_NLeptons[i][6], v_lepton_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1pt[i][6] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1eta[i][6], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1phi[i][6], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2pt[i][6] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2eta[i][6], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2phi[i][6], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NPV[i][6]    , num_pv        , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NJets[i][6]  , v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1pt[i][6] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1eta[i][6], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1phi[i][6], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2pt[i][6] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2eta[i][6], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2phi[i][6], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][6] ,Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][6],Met->Phi() , v_SystEvt[i] );
                           
                           FillHisto( h_cf_sys_dilep_inv_mass[i][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][6] , Met->Pt() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][6], Met->Phi(), v_SystEvt[i] );
                           
                           FillHisto( h_sys_Num_PV[i][6], num_pv, v_SystEvt[i] );
                           FillHisto( h_sys_DiLepMass[i][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_sys_Lep1pt[i][6] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1eta[i][6], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1phi[i][6], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2pt[i][6] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2eta[i][6], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2phi[i][6], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1pt[i][6] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1eta[i][6], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1phi[i][6], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2pt[i][6] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2eta[i][6], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2phi[i][6], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Num_Jets[i][6], v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_sys_METpt[i][6]   , Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_METphi[i][6]  , Met->Phi()  , v_SystEvt[i] );
                        }
                     }
                  } // One or More b-tagging // Step 5 
                   
                  if (nbtagged ==2)// exactly 2b-tagging ..
                  {
//                     ApplyBTagWeight(-2);
/*            cout  << " Info_EventNumber : " << Info_EventNumber << " step 7 evt_weight_ : " << v_SystEvt[1]
                  << " Lep1 Pt : " << (*Lep1).Pt()
                  << " Lep2 Pt : " << (*Lep2).Pt()
                  << " Jet1 Pt : " << (*Jet1).Pt()
                  << " Jet2 Pt : " << (*Jet2).Pt()
                  << " bJet1 Pt : " << (*bJet1).Pt()
                  << " bJet2 Pt : " << (*bJet2).Pt()
                  << " MET : " << Met->Pt()
                  << endl;*/
                     FillHisto( h_EventWeight[7], evt_weight_  );
                     FillHisto( h_cf_NLeptons[7], v_lepton_idx.size(), evt_weight_ );
                  
                     FillHisto( h_cf_Lep1pt[7] , (*Lep1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep1eta[7], (*Lep1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep1phi[7], (*Lep1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Lep2pt[7] , (*Lep2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep2eta[7], (*Lep2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep2phi[7], (*Lep2).Phi() , evt_weight_ );
                     FillHisto( h_cf_NPV[7]    , num_pv        , evt_weight_ );
                     FillHisto( h_cf_NJets[7]  , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_cf_Jet1pt[7] , (*Jet1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet1eta[7], (*Jet1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet1phi[7], (*Jet1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Jet2pt[7] , (*Jet2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet2eta[7], (*Jet2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet2phi[7], (*Jet2).Phi() , evt_weight_ );
                  
                     FillHisto( h_cf_dilep_inv_mass[7], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                     FillHisto( h_cf_metpt[7] , Met->Pt() , evt_weight_ );
                     FillHisto( h_cf_metphi[7], Met->Phi(), evt_weight_ );
                  
                     FillHisto( h_Lep1pt[7] , Lep1->Pt() , evt_weight_ );
                     FillHisto( h_Lep2pt[7] , Lep2->Pt() , evt_weight_ );
                     FillHisto( h_Lep1eta[7], Lep1->Eta(), evt_weight_ );
                     FillHisto( h_Lep2eta[7], Lep2->Eta(), evt_weight_ );
                     FillHisto( h_Lep1phi[7], Lep1->Phi(), evt_weight_ );
                     FillHisto( h_Lep2phi[7], Lep2->Phi(), evt_weight_ );
                    
                     if (TString(Decaymode).Contains("muel"))
                     { 
                        FillHisto( h_Muonpt[7]  , TMuon->Pt()      , evt_weight_ );
                        FillHisto( h_Elecpt[7]  , TElectron->Pt()  , evt_weight_ );
                        FillHisto( h_Muoneta[7] , TMuon->Eta()     , evt_weight_ );
                        FillHisto( h_Eleceta[7] , TElectron->Eta() , evt_weight_ );
                        FillHisto( h_Muonphi[7] , TMuon->Phi()     , evt_weight_ );
                        FillHisto( h_Elecphi[7] , TElectron->Phi() , evt_weight_ );
                     }
                     FillHisto( h_Jet1pt[7] , Jet1->Pt() , evt_weight_ );
                     FillHisto( h_Jet2pt[7] , Jet2->Pt() , evt_weight_ );
                     FillHisto( h_Jet1eta[7], Jet1->Eta(), evt_weight_ );
                     FillHisto( h_Jet2eta[7], Jet2->Eta(), evt_weight_ );
                     FillHisto( h_Jet1phi[7], Jet1->Phi(), evt_weight_ );
                     FillHisto( h_Jet2phi[7], Jet2->Phi(), evt_weight_ );
                     FillHisto( h_METpt[7]  , Met->Pt()  , evt_weight_ );
                     FillHisto( h_METphi[7] , Met->Phi() , evt_weight_ );
                     FillHisto( h_HT[7] , AllJetpt, evt_weight_);
                  
                     FillHisto( h_DiLepMass[7], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                  
                     FillHisto( h_Num_PV[7]   , num_pv          , evt_weight_ );
                     FillHisto( h_Num_Jets[7] , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_Num_bJets[7], nbtagged        , evt_weight_ );
                  
                     FillHisto( h_Topo_Apla[2] , _apla, evt_weight_);
                     FillHisto( h_Topo_Sphe[2] , _sphe, evt_weight_);
                     FillHisto( h_Topo_Plan[2] , _plan, evt_weight_);
                     if (isAllSyst == true) {
                        for ( int i = 0; i < v_SystFullName.size(); ++i )
                        {  
                           if (TString(v_SystFullName[i]).Contains("Jet")){continue;}
                           FillHisto( h_cf_sys_NLeptons[i][7], v_lepton_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1pt[i][7] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1eta[i][7], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1phi[i][7], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2pt[i][7] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2eta[i][7], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2phi[i][7], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NPV[i][7]    , num_pv        , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NJets[i][7]  , v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1pt[i][7] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1eta[i][7], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1phi[i][7], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2pt[i][7] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2eta[i][7], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2phi[i][7], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][7] ,Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][7],Met->Phi() , v_SystEvt[i] );
                           
                           FillHisto( h_cf_sys_dilep_inv_mass[i][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][7] , Met->Pt() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][7], Met->Phi(), v_SystEvt[i] );
                           
                           FillHisto( h_sys_Num_PV[i][7], num_pv, v_SystEvt[i] );
                           FillHisto( h_sys_DiLepMass[i][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_sys_Lep1pt[i][7] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1eta[i][7], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1phi[i][7], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2pt[i][7] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2eta[i][7], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2phi[i][7], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1pt[i][7] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1eta[i][7], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1phi[i][7], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2pt[i][7] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2eta[i][7], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2phi[i][7], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Num_Jets[i][7], v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_sys_METpt[i][7]   , Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_METphi[i][7]  , Met->Phi()  , v_SystEvt[i] );
                        }
                     }
                  }
                  ////////////////////////////////////////
                  // Reconstruction of Top with KinSol. //
                  ////////////////////////////////////////
                  //KinSol();   
                  SetUpKINObs(); 
                  //SetGenLepAnLep();
                  //if ( ksolweight_ != -1 )
                  if (isKinSol)
                  {
                  /*cout  << " First Info_EventNumber : " << Info_EventNumber << " step 8 evt_weight_ : " << evt_weight_ 
                  << " Lep1 Pt : " << (*Lep1).Pt()
                  << " Lep2 Pt : " << (*Lep2).Pt()
                  << " Jet1 Pt : " << (*Jet1).Pt()
                  << " Jet2 Pt : " << (*Jet2).Pt()
                  << " bJet Pt : " << (*bJet).Pt()
                  << " AnbJet Pt : " << (*AnbJet).Pt()
                  << " MET : " << Met->Pt()
                  << " Top : " << Top->Pt()
                  << " AnTop : " << AnTop->Pt()
                  << " cp01 : " << ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep )
                  << " cp03 : " << ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep )
                  << endl;*/
                     BandBbarJetDiff();
/*            cout  << " Second Info_EventNumber : " << Info_EventNumber << " step 8 evt_weight_ : " << evt_weight_ 
                  << " Lep1 Pt : " << (*Lep1).Pt()
                  << " Lep2 Pt : " << (*Lep2).Pt()
                  << " Jet1 Pt : " << (*Jet1).Pt()
                  << " Jet2 Pt : " << (*Jet2).Pt()
                  << " bJet Pt : " << (*bJet).Pt()
                  << " AnbJet Pt : " << (*AnbJet).Pt()
                  << " MET : " << Met->Pt()
                  << " Top : " << Top->Pt()
                  << " AnTop : " << AnTop->Pt()
                  << " cp03 : " << ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep )
                  << endl;*/
                     FillHisto( h_EventWeight[8], evt_weight_  );
                     FillHisto( h_cf_NLeptons[8], v_lepton_idx.size(), evt_weight_ );
                  
                     FillHisto( h_cf_Lep1pt[8] , (*Lep1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep1eta[8], (*Lep1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep1phi[8], (*Lep1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Lep2pt[8] , (*Lep2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep2eta[8], (*Lep2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep2phi[8], (*Lep2).Phi() , evt_weight_ );
                     FillHisto( h_cf_NPV[8]    , num_pv        , evt_weight_ );
                     FillHisto( h_cf_NJets[8]  , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_cf_Jet1pt[8] , (*Jet1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet1eta[8], (*Jet1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet1phi[8], (*Jet1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Jet2pt[8] , (*Jet2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet2eta[8], (*Jet2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet2phi[8], (*Jet2).Phi() , evt_weight_ );
                  
                     FillHisto( h_cf_dilep_inv_mass[8], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                     FillHisto( h_cf_metpt[8] , Met->Pt() , evt_weight_ );
                     FillHisto( h_cf_metphi[8], Met->Phi(), evt_weight_ );
                  
                     FillHisto( h_Lep1pt[8] , Lep1->Pt() , evt_weight_ );
                     FillHisto( h_Lep2pt[8] , Lep2->Pt() , evt_weight_ );
                     FillHisto( h_Lep1eta[8], Lep1->Eta(), evt_weight_ );
                     FillHisto( h_Lep2eta[8], Lep2->Eta(), evt_weight_ );
                     FillHisto( h_Lep1phi[8], Lep1->Phi(), evt_weight_ );
                     FillHisto( h_Lep2phi[8], Lep2->Phi(), evt_weight_ );
                    
                     if (TString(Decaymode).Contains("muel"))
                     { 
                        FillHisto( h_Muonpt[8]  , TMuon->Pt()      , evt_weight_ );
                        FillHisto( h_Elecpt[8]  , TElectron->Pt()  , evt_weight_ );
                        FillHisto( h_Muoneta[8] , TMuon->Eta()     , evt_weight_ );
                        FillHisto( h_Eleceta[8] , TElectron->Eta() , evt_weight_ );
                        FillHisto( h_Muonphi[8] , TMuon->Phi()     , evt_weight_ );
                        FillHisto( h_Elecphi[8] , TElectron->Phi() , evt_weight_ );
                     }
                     FillHisto( h_Jet1pt[8] , Jet1->Pt() , evt_weight_ );
                     FillHisto( h_Jet2pt[8] , Jet2->Pt() , evt_weight_ );
                     FillHisto( h_Jet1eta[8], Jet1->Eta(), evt_weight_ );
                     FillHisto( h_Jet2eta[8], Jet2->Eta(), evt_weight_ );
                     FillHisto( h_Jet1phi[8], Jet1->Phi(), evt_weight_ );
                     FillHisto( h_Jet2phi[8], Jet2->Phi(), evt_weight_ );
                     FillHisto( h_METpt[8]  , Met->Pt()  , evt_weight_ );
                     FillHisto( h_METphi[8] , Met->Phi() , evt_weight_ );
                     FillHisto( h_HT[8]     , AllJetpt   , evt_weight_);
                  
                     FillHisto( h_DiLepMass[8], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                  
                     FillHisto( h_Num_PV[8]   , num_pv          , evt_weight_ );
                     FillHisto( h_Num_Jets[8] , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_Num_bJets[8], nbtagged        , evt_weight_ );
                  
                     FillHisto( h_Topo_Apla[3] , _apla, evt_weight_);
                     FillHisto( h_Topo_Sphe[3] , _sphe, evt_weight_);
                     FillHisto( h_Topo_Plan[3] , _plan, evt_weight_);
                     if ( Top->Pt() > AnTop->Pt() ) { (*Top1) = (*Top); (*Top2) = (*AnTop); }
                     else { (*Top1) = (*AnTop); (*Top2) = (*Top); }

                     if (isAllSyst == true){
                        for ( int i = 0; i < v_SystFullName.size(); ++i )
                        {  
                           if (TString(v_SystFullName[i]).Contains("Jet")){continue;}
                           FillHisto( h_cf_sys_NLeptons[i][8], v_lepton_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1pt[i][8] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1eta[i][8], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1phi[i][8], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2pt[i][8] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2eta[i][8], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2phi[i][8], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NPV[i][8]    , num_pv        , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NJets[i][8]  , v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1pt[i][8] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1eta[i][8], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1phi[i][8], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2pt[i][8] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2eta[i][8], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2phi[i][8], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][8] ,Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][8],Met->Phi() , v_SystEvt[i] );
                           
                           FillHisto( h_cf_sys_dilep_inv_mass[i][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][8] , Met->Pt() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][8], Met->Phi(), v_SystEvt[i] );
                           
                           FillHisto( h_sys_Num_PV[i][8], num_pv, v_SystEvt[i] );
                           FillHisto( h_sys_DiLepMass[i][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_sys_Lep1pt[i][8] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1eta[i][8], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1phi[i][8], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2pt[i][8] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2eta[i][8], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2phi[i][8], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1pt[i][8] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1eta[i][8], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1phi[i][8], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2pt[i][8] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2eta[i][8], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2phi[i][8], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Num_Jets[i][8], v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_sys_METpt[i][8]   , Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_METphi[i][8]  , Met->Phi()  , v_SystEvt[i] );

                           FillHisto( h_sys_Top1Mass_[i]    , Top1->M()        , v_SystEvt[i] );
                           FillHisto( h_sys_Top1pt_[i]      , Top1->Pt()       , v_SystEvt[i] );
                           FillHisto( h_sys_Top1phi_[i]     , Top1->Phi()      , v_SystEvt[i] );
                           FillHisto( h_sys_Top1Rapidity_[i], Top1->Rapidity() , v_SystEvt[i] );
                           FillHisto( h_sys_Top1Energy_[i]  , Top1->Energy()   , v_SystEvt[i] );
                           FillHisto( h_sys_Top2Mass_[i]    , Top2->M()        , v_SystEvt[i] );
                           FillHisto( h_sys_Top2pt_[i]      , Top2->Pt()       , v_SystEvt[i] );
                           FillHisto( h_sys_Top2phi_[i]     , Top2->Phi()      , v_SystEvt[i] );
                           FillHisto( h_sys_Top2Rapidity_[i], Top2->Rapidity() , v_SystEvt[i] );
                           FillHisto( h_sys_Top2Energy_[i]  , Top2->Energy()   , v_SystEvt[i] );

                           FillHisto( h_sys_TopMass_[i]      , Top->M()         , v_SystEvt[i] );
                           FillHisto( h_sys_Toppt_[i]        , Top->Pt()        , v_SystEvt[i] );
                           FillHisto( h_sys_Topphi_[i]       , Top->Phi()       , v_SystEvt[i] );
                           FillHisto( h_sys_TopRapidity_[i]  , Top->Rapidity()  , v_SystEvt[i] );
                           FillHisto( h_sys_TopEnergy_[i]    , Top->Energy()    , v_SystEvt[i] );
                           FillHisto( h_sys_AnTopMass_[i]    , AnTop->M()       , v_SystEvt[i] );
                           FillHisto( h_sys_AnToppt_[i]      , AnTop->Pt()      , v_SystEvt[i] );
                           FillHisto( h_sys_AnTopphi_[i]     , AnTop->Phi()     , v_SystEvt[i] );
                           FillHisto( h_sys_AnTopRapidity_[i], AnTop->Rapidity(), v_SystEvt[i] );
                           FillHisto( h_sys_AnTopEnergy_[i]  , AnTop->Energy()  , v_SystEvt[i] );

                        }
                     }
                  
                     FillHisto( h_TopMass      , Top->M()         , evt_weight_ );
                     FillHisto( h_Toppt        , Top->Pt()        , evt_weight_ );
                     FillHisto( h_Topphi       , Top->Phi()       , evt_weight_ );
                     FillHisto( h_TopRapidity  , Top->Rapidity()  , evt_weight_ );
                     FillHisto( h_TopEnergy    , Top->Energy()    , evt_weight_ );
                     FillHisto( h_AnTopMass    , AnTop->M()       , evt_weight_ );
                     FillHisto( h_AnToppt      , AnTop->Pt()      , evt_weight_ );
                     FillHisto( h_AnTopphi     , AnTop->Phi()     , evt_weight_ );
                     FillHisto( h_AnTopRapidity, AnTop->Rapidity(), evt_weight_ );
                     FillHisto( h_AnTopEnergy  , AnTop->Energy()  , evt_weight_ );
                  
                  
                     FillHisto( h_W1Mass , W1->M()  , evt_weight_ );
                     FillHisto( h_W2Mass , W2->M()  , evt_weight_ );
                  
                     FillHisto( h_W1Mt , W1->Mt()  , evt_weight_ );
                     FillHisto( h_W2Mt , W2->Mt()  , evt_weight_ );
                  
                     FillHisto( h_bJet1Energy , bJet1->Energy()  , evt_weight_ );
                     FillHisto( h_bJet2Energy , bJet2->Energy()  , evt_weight_ );
                  
                     FillHisto( h_bJetEnergy   , bJet->Energy()   , evt_weight_ );
                     FillHisto( h_AnbJetEnergy , AnbJet->Energy() , evt_weight_ );
                     FillHisto( h_LepEnergy    , Lep->Energy()    , evt_weight_ );
                     FillHisto( h_AnLepEnergy  , AnLep->Energy()  , evt_weight_ );
                     FillHisto( h_NuEnergy     , Nu->Energy()     , evt_weight_ );
                     FillHisto( h_AnNuEnergy   , AnNu->Energy()   , evt_weight_ );
                     // kinematics of objects after top reco
                     FillHisto( h_bJetpt       , bJet->Pt()       , evt_weight_ );
                     FillHisto( h_bJeteta      , bJet->Eta()      , evt_weight_ );
                     FillHisto( h_bJetphi      , bJet->Phi()      , evt_weight_ );
                     FillHisto( h_AnbJetpt     , AnbJet->Pt()     , evt_weight_ );
                     FillHisto( h_AnbJeteta    , AnbJet->Eta()    , evt_weight_ );
                     FillHisto( h_AnbJetphi    , AnbJet->Phi()    , evt_weight_ );
                     FillHisto( h_Leppt        , Lep->Pt()        , evt_weight_ );
                     FillHisto( h_Lepeta       , Lep->Eta()       , evt_weight_ );
                     FillHisto( h_Lepphi       , Lep->Phi()       , evt_weight_ );
                     FillHisto( h_AnLeppt      , AnLep->Pt()      , evt_weight_ );
                     FillHisto( h_AnLepeta     , AnLep->Eta()     , evt_weight_ );
                     FillHisto( h_AnLepphi     , AnLep->Phi()     , evt_weight_ );
                  
                     /// Kin Solver Purity ... ///
                     FillHisto( h2_TopMassVsLepBMass    , Top->M() ,((*Lep)+(*AnbJet)).M()   , evt_weight_ );
                     FillHisto( h2_AnTopMassVsLepBMass  , AnTop->M() ,((*AnLep)+(*bJet)).M()   , evt_weight_ );
                     FillHisto( h2_AnLepBMassVsLepBMass , ((*Lep)+(*AnbJet)).M() ,((*AnLep)+(*bJet)).M()   , evt_weight_ );
                  
                     FillHisto( h_LepbJetMass   , ((*Lep)+(*AnbJet)).M() , evt_weight_ );
                     FillHisto( h_AnLepbJetMass , ((*AnLep)+(*bJet)).M() , evt_weight_ );
                  
                     LepAnLepMisCharge();
                     
                  
                     // Calculating O3 variable
                     double cpO3 = ssbcpviol->getO3Vari( bJet , AnbJet, AnLep, Lep );
                     double cpO3_jprup = ssbcpviol->getO3VariJPRUp( bJet ,AnbJet , AnLep ,Lep );
                     double cpO3_jprdown = ssbcpviol->getO3VariJPRDown( bJet ,AnbJet , AnLep ,Lep );
                  
                     double cpOb = ssbcpviol->getObVari( bJet , AnbJet, AnLep, Lep );
                     double cpO5 = ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep );
                     if (cpO3 == 0.0 ){cout << "--00000000000000--" << cpO3 << endl;}
                     // Calculating O3 error.
                     if (cpO3 > 0.0) {Num_cpo3_p++;
      //                  cout << "evt_weight_ ^2 " << evt_weight_*evt_weight_ << endl;
                        FillHisto( h_CPO3_evtweight,0.5,evt_weight_);
                        FillHisto( h_CPO3_evtweight,6.5,evt_weight_*evt_weight_);
                        FillHisto( h_CPO3_Plus,cpO3,evt_weight_);
                     }
                     else           { Num_cpo3_m++; Nm_eventw_2 += (evt_weight_*evt_weight_);
                        FillHisto( h_CPO3_evtweight,2.5,evt_weight_);
                        FillHisto( h_CPO3_evtweight,8.5,evt_weight_*evt_weight_);
                        laecO3 += 50;
                        FillHisto( h_CPO3_Minus,cpO3,evt_weight_);
                     }
                  
                     if (cpOb > 0.0) {
                        FillHisto( h_CPOb_evtweight,0.5,evt_weight_);
                        FillHisto( h_CPOb_evtweight,6.5,evt_weight_*evt_weight_);
                     }
                     else           { 
                        FillHisto( h_CPOb_evtweight,2.5,evt_weight_);
                        FillHisto( h_CPOb_evtweight,8.5,evt_weight_*evt_weight_);
                        laecOb += 50;
                     }
                  
                     if (cpO5 > 0.0) {
                        FillHisto( h_CPO5_evtweight,0.5,evt_weight_);
                        FillHisto( h_CPO5_evtweight,6.5,evt_weight_*evt_weight_);
                     }
                     else           { 
                        FillHisto( h_CPO5_evtweight,2.5,evt_weight_);
                        FillHisto( h_CPO5_evtweight,8.5,evt_weight_*evt_weight_);
                        laecO5 += 50;
                     }
                  
                     /// Jet pT Variation Study For Dilution
                     if (cpO3 > 0.0 && cpO3_jprup > 0.0 )FillHisto( h_CPO3_JPRUp_Plus_Plus ,  cpO3_jprup, evt_weight_ );
                     if (cpO3 > 0.0 && cpO3_jprup < 0.0 )FillHisto( h_CPO3_JPRUp_Plus_Minus,  cpO3_jprup, evt_weight_ );
                     if (cpO3 < 0.0 && cpO3_jprup < 0.0 )FillHisto( h_CPO3_JPRUp_Minus_Minus, cpO3_jprup, evt_weight_ );
                     if (cpO3 < 0.0 && cpO3_jprup > 0.0 )FillHisto( h_CPO3_JPRUp_Minus_Plus,  cpO3_jprup, evt_weight_ );
                  
                     if (cpO3 > 0.0 && cpO3_jprdown > 0.0 )FillHisto( h_CPO3_JPRDown_Plus_Plus ,  cpO3_jprdown, evt_weight_ );
                     if (cpO3 > 0.0 && cpO3_jprdown < 0.0 )FillHisto( h_CPO3_JPRDown_Plus_Minus,  cpO3_jprdown, evt_weight_ );
                     if (cpO3 < 0.0 && cpO3_jprdown < 0.0 )FillHisto( h_CPO3_JPRDown_Minus_Minus, cpO3_jprdown, evt_weight_ );
                     if (cpO3 < 0.0 && cpO3_jprdown > 0.0 )FillHisto( h_CPO3_JPRDown_Minus_Plus,  cpO3_jprdown, evt_weight_ );
                  
                  
//                     cout << "laec ? " << laec << endl;
                     
                     // Calculating Bjorken variables
                     Bjorken( Lep , AnLep , bJet, AnbJet , Nu , AnNu ); 
                  
                     FillHisto( h_CPO3_reco        , cpO3         , evt_weight_ ); // for the test, not use now
                     FillHisto( h_CPO3_reco_JPRUp  , cpO3_jprup   , evt_weight_ );
                     FillHisto( h_CPO3_reco_JPRDown, cpO3_jprdown , evt_weight_ );
                     FillHisto( h_CPOb_reco, ssbcpviol->getObVari( bJet , AnbJet, AnLep, Lep ) , evt_weight_ );
                     FillHisto( h_CPO5_reco, ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) , evt_weight_ );
                     FillHisto( h_BjorkenX1, x1_bj, evt_weight_ );
                     FillHisto( h_BjorkenX2, x2_bj, evt_weight_ );
                     FillHisto( h_BjorkenX3, x3_bj, evt_weight_ );
                     
                     // Get CP-Violation Variables //
                     v_recocp_O.clear();
                     v_recocp_O.push_back( ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO2Vari( Top, AnTop, bJet, AnbJet ) );
                     v_recocp_O.push_back( ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO4Vari( AnbJet, bJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO6Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO7Vari( Top , AnTop, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO8Vari( Top, AnTop, bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO9Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO10Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO11Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO12Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO13Vari( bJet , AnbJet, AnLep, Lep )  );

                     NewLepAnLepMisCharge();
                     for (int i = 0; i < v_recocp_O.size(); ++ i)
                     {
                        FillHisto( h_Reco_CPO_[i], v_recocp_O[i] , evt_weight_ );
                        FillHisto( h_Reco_CPO_ReRange_[i], v_recocp_O[i] , evt_weight_ );
                        FillHisto( h_LepAnLepEngCheck_[i], v_MisCharge[i] , evt_weight_ );
                     }
                     if (isAllSyst == true){
                        for (int i = 0; i <  v_SystEvt.size(); ++ i)
                        {
                           if (TString(v_SystFullName[i]).Contains("Jet")){continue;}
                           for (int j = 0; j < v_recocp_O.size(); ++ j )
                           {
                              FillHisto( h_sys_Reco_CPO_[i][j], v_recocp_O[j] , v_SystEvt[i]   );
                              FillHisto( h_sys_Reco_CPO_ReRange_[i][j], v_recocp_O[j] , v_SystEvt[i]  );
                           }
                        }
                     } 
                     FillHisto( h_Top1Mass    , Top1->M()        , evt_weight_ );
                     FillHisto( h_Top1pt      , Top1->Pt()       , evt_weight_ );
                     FillHisto( h_Top1phi     , Top1->Phi()      , evt_weight_ );
                     FillHisto( h_Top1Rapidity, Top1->Rapidity() , evt_weight_ );
                     FillHisto( h_Top1Energy  , Top1->Energy()   , evt_weight_ );
                  
                     FillHisto( h_Top2Mass    , Top2->M()        , evt_weight_ );
                     FillHisto( h_Top2pt      , Top2->Pt()       , evt_weight_ );
                     FillHisto( h_Top2phi     , Top2->Phi()      , evt_weight_ );
                     FillHisto( h_Top2Rapidity, Top2->Rapidity() , evt_weight_ );
                     FillHisto( h_Top2Energy  , Top2->Energy()   , evt_weight_ );
                     // BJet - BJet //
                     FillHisto( h_Diff_BBar_Pt    , fabs( bJet->Pt() - AnbJet->Pt() )         , evt_weight_ );
                     FillHisto( h_Diff_BBar_Energy, fabs( bJet->Energy() - AnbJet->Energy() ) , evt_weight_ );
                     FillHisto( h_Diff_BBar_P     , fabs( bJet->P() - AnbJet->P() )           , evt_weight_ );
                  
                     FillHisto( h_Diff_LepAnLep_Pt    , fabs( Lep->Pt() - AnLep->Pt() )         , evt_weight_ );
                     FillHisto( h_Diff_LepAnLep_Energy, fabs( Lep->Energy() - AnLep->Energy() ) , evt_weight_ );
                     FillHisto( h_Diff_LepAnLep_P     , fabs( Lep->P() - AnLep->P() )           , evt_weight_ );
                  
                     FillHisto( h_LepAnLepEngCheckO3    , laecO3 , evt_weight_ );
                     FillHisto( h_LepAnLepEngCheckOb    , laecOb , evt_weight_ );
                     FillHisto( h_LepAnLepEngCheckO5    , laecO5 , evt_weight_ );
                  
                     ////////////////////////////////
                     /// JetAngular Study Area ...///
                     ////////////////////////////////
                     double recocpO3_etavariup   = ssbcpviol->getO3Vari( ssbcpviol->JetAngEta(bJet,"up"  ),ssbcpviol->JetAngEta(AnbJet,"up"  ), AnLep, Lep );
                     double recocpO3_etavaridown = ssbcpviol->getO3Vari( ssbcpviol->JetAngEta(bJet,"down"),ssbcpviol->JetAngEta(AnbJet,"down"), AnLep, Lep );
                     double recocpO3_phivariup   = ssbcpviol->getO3Vari( ssbcpviol->JetAngPhi(bJet,"up"  ),ssbcpviol->JetAngPhi(AnbJet,"up"  ), AnLep, Lep );
                     double recocpO3_phivaridown = ssbcpviol->getO3Vari( ssbcpviol->JetAngPhi(bJet,"down"),ssbcpviol->JetAngPhi(AnbJet,"down"), AnLep, Lep );
                  
                     FillHisto( h_CPO3_EtaVariUp,  recocpO3_etavariup,  evt_weight_);
                     FillHisto( h_CPO3_EtaVariDown,recocpO3_etavaridown,evt_weight_);
                     FillHisto( h_CPO3_PhiVariUp,  recocpO3_phivariup,  evt_weight_);
                     FillHisto( h_CPO3_PhiVariDown,recocpO3_phivaridown,evt_weight_);
                     if (cpO3 > 0.0) {
                        FillHisto( h_CPO3reco_evtweight,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight,8.5,evt_weight_*evt_weight_);
                     }
                  
                     //////////////////////////////////////////////////////////// 
                     /// CPO3 with Eta Variation (Up) at Reconstruction Level ///
                     //////////////////////////////////////////////////////////// 
                     if ( recocpO3_etavariup > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_EtaVariUp,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_EtaVariUp,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_EtaVariUp,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_EtaVariUp,8.5,evt_weight_*evt_weight_);
                     }
                  
                     if ( cpO3 > 0.0  && recocpO3_etavariup > 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariUp_Plus_Plus,recocpO3_etavariup,evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_etavariup < 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariUp_Plus_Minus,recocpO3_etavariup,evt_weight_);
                     }
                  
                     if ( cpO3 < 0.0  && recocpO3_etavariup > 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariUp_Minus_Plus,recocpO3_etavariup,evt_weight_);
                     }
                     if ( cpO3 < 0.0  && recocpO3_etavariup < 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariUp_Minus_Minus,recocpO3_etavariup,evt_weight_);
                     }
                  
                     ////////////////////////////////////////////////////////////// 
                     /// CPO3 with Eta Variation (Down) at Reconstruction Level ///
                     //////////////////////////////////////////////////////////////
                     if (recocpO3_etavaridown > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_EtaVariDown,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_EtaVariDown,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_EtaVariDown,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_EtaVariDown,8.5,evt_weight_*evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_etavaridown > 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariDown_Plus_Plus,recocpO3_etavaridown,evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_etavaridown < 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariDown_Plus_Minus,recocpO3_etavaridown,evt_weight_);
                     }
                  
                     if ( cpO3 < 0.0  && recocpO3_etavaridown > 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariDown_Minus_Plus,recocpO3_etavaridown,evt_weight_);
                     }
                     if ( cpO3 < 0.0  && recocpO3_etavaridown < 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariDown_Minus_Minus,recocpO3_etavaridown,evt_weight_);
                     }
                  
                     //////////////////////////////////////////////////////////// 
                     /// CPO3 with Phi Variation (Up) at Reconstruction Level ///
                     //////////////////////////////////////////////////////////// 
                     if (recocpO3_phivariup > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,8.5,evt_weight_*evt_weight_);
                     }
                     if ( recocpO3_phivariup > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,8.5,evt_weight_*evt_weight_);
                     }
                  
                     if ( cpO3 > 0.0  && recocpO3_phivariup > 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariUp_Plus_Plus,recocpO3_phivariup,evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_phivariup < 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariUp_Plus_Minus,recocpO3_phivariup,evt_weight_);
                     }
                  
                     if ( cpO3 < 0.0  && recocpO3_phivariup > 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariUp_Minus_Plus,recocpO3_phivariup,evt_weight_);
                     }
                     if ( cpO3 < 0.0  && recocpO3_phivariup < 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariUp_Minus_Minus,recocpO3_phivariup,evt_weight_);
                     }
                     ////////////////////////////////////////////////////////////// 
                     /// CPO3 with Phi Variation (Down) at Reconstruction Level ///
                     //////////////////////////////////////////////////////////////
                  
                     if (recocpO3_phivaridown > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_PhiVariDown,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariDown,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_PhiVariDown,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariDown,8.5,evt_weight_*evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_phivaridown > 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariDown_Plus_Plus,recocpO3_phivaridown,evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_phivaridown < 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariDown_Plus_Minus,recocpO3_phivaridown,evt_weight_);
                     }
                  
                     if ( cpO3 < 0.0  && recocpO3_phivaridown > 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariDown_Minus_Plus,recocpO3_phivaridown,evt_weight_);
                     }
                     if ( cpO3 < 0.0  && recocpO3_phivaridown < 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariDown_Minus_Minus,recocpO3_phivaridown,evt_weight_);
                     }

                  } // KinSolver // for central !!
               } // bJetcut 1- or More //
            }  // MET Cut 
         } // Num Jet Cut 

         ///////////////////////////////
         /// JES Up Systematic Study ///
         ///////////////////////////////
         if ( NumJetCut(v_jetup_idx) == true && isAllSyst == true && idx_jecup != -2 )
         {

            FillHisto( h_cf_sys_NLeptons[idx_jecup][3], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep1pt[idx_jecup][3] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep1eta[idx_jecup][3], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep1phi[idx_jecup][3], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep2pt[idx_jecup][3] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep2eta[idx_jecup][3], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep2phi[idx_jecup][3], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_NPV[idx_jecup][3]    , num_pv        , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_NJets[idx_jecup][3]  , v_jetup_idx.size(), v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet1pt[idx_jecup][3] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet1eta[idx_jecup][3], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet1phi[idx_jecup][3], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet2pt[idx_jecup][3] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet2eta[idx_jecup][3], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet2phi[idx_jecup][3], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_metpt[idx_jecup][3] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_metphi[idx_jecup][3],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
            
            FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_metpt[idx_jecup][3] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_metphi[idx_jecup][3], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
            
            FillHisto( h_sys_Num_PV[idx_jecup][3], num_pv, v_SystEvt[idx_jecup] );
            FillHisto( h_sys_DiLepMass[idx_jecup][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep1pt[idx_jecup][3] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep1eta[idx_jecup][3], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep1phi[idx_jecup][3], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep2pt[idx_jecup][3] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep2eta[idx_jecup][3], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep2phi[idx_jecup][3], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet1pt[idx_jecup][3] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet1eta[idx_jecup][3], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet1phi[idx_jecup][3], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet2pt[idx_jecup][3] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet2eta[idx_jecup][3], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet2phi[idx_jecup][3], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Num_Jets[idx_jecup][3], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
            FillHisto( h_sys_METpt[idx_jecup][3]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_METphi[idx_jecup][3]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );
            if ( METCut(MetJESUp) == true ) 
            {
               FillHisto( h_cf_sys_NLeptons[idx_jecup][4], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep1pt[idx_jecup][4] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep1eta[idx_jecup][4], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep1phi[idx_jecup][4], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep2pt[idx_jecup][4] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep2eta[idx_jecup][4], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep2phi[idx_jecup][4], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_NPV[idx_jecup][4]    , num_pv        , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_NJets[idx_jecup][4]  , v_jetup_idx.size(), v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet1pt[idx_jecup][4] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet1eta[idx_jecup][4], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet1phi[idx_jecup][4], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet2pt[idx_jecup][4] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet2eta[idx_jecup][4], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet2phi[idx_jecup][4], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_metpt[idx_jecup][4] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_metphi[idx_jecup][4],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_metpt[idx_jecup][4] , Met->Pt() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_metphi[idx_jecup][4], Met->Phi(), v_SystEvt[idx_jecup] );
               
               FillHisto( h_sys_Num_PV[idx_jecup][4], num_pv, v_SystEvt[idx_jecup] );
               FillHisto( h_sys_DiLepMass[idx_jecup][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep1pt[idx_jecup][4] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep1eta[idx_jecup][4], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep1phi[idx_jecup][4], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep2pt[idx_jecup][4] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep2eta[idx_jecup][4], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep2phi[idx_jecup][4], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet1pt[idx_jecup][4] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet1eta[idx_jecup][4], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet1phi[idx_jecup][4], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet2pt[idx_jecup][4] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet2eta[idx_jecup][4], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet2phi[idx_jecup][4], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Num_Jets[idx_jecup][4], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
               FillHisto( h_sys_METpt[idx_jecup][4]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_METphi[idx_jecup][4]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );
               
               /// Apply BTagging Scale Factor for JES Up
               //BTaggigSFApplyJESR("JetEnUp");
               BTaggigSFApplyJESR(v_jetup_idx, v_jetup_TL, idx_jecup);
               if ( BJetCut(v_bjetup_idx) == true ) // one or more b-tagging 
               {


                  FillHisto( h_cf_sys_NLeptons[idx_jecup][5], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep1pt[idx_jecup][5] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep1eta[idx_jecup][5], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep1phi[idx_jecup][5], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep2pt[idx_jecup][5] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep2eta[idx_jecup][5], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep2phi[idx_jecup][5], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_NPV[idx_jecup][5]    , num_pv        , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_NJets[idx_jecup][5]  , v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet1pt[idx_jecup][5] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet1eta[idx_jecup][5], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet1phi[idx_jecup][5], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet2pt[idx_jecup][5] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet2eta[idx_jecup][5], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet2phi[idx_jecup][5], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_metpt[idx_jecup][5] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_metphi[idx_jecup][5],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
                  
                  FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_metpt[idx_jecup][5] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_metphi[idx_jecup][5], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
                  
                  FillHisto( h_sys_Num_PV[idx_jecup][5], num_pv, v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_DiLepMass[idx_jecup][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep1pt[idx_jecup][5] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep1eta[idx_jecup][5], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep1phi[idx_jecup][5], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep2pt[idx_jecup][5] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep2eta[idx_jecup][5], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep2phi[idx_jecup][5], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet1pt[idx_jecup][5] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet1eta[idx_jecup][5], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet1phi[idx_jecup][5], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet2pt[idx_jecup][5] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet2eta[idx_jecup][5], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet2phi[idx_jecup][5], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Num_Jets[idx_jecup][5], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_METpt[idx_jecup][5]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_METphi[idx_jecup][5]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );

                  ///////////////////////////////////////
                  /// 2 or More B-Tagging Requriement ///
                  /////////////////////////////////////// 
                  if ( DoubleBtag(v_bjetup_idx) == true )
                  {

                     FillHisto( h_cf_sys_NLeptons[idx_jecup][6], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecup][6] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecup][6], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecup][6], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecup][6] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecup][6], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecup][6], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NPV[idx_jecup][6]    , num_pv        , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NJets[idx_jecup][6]  , v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecup][6] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecup][6], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecup][6], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecup][6] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecup][6], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecup][6], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][6] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][6],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
                     
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][6] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][6], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
                     
                     FillHisto( h_sys_Num_PV[idx_jecup][6], num_pv, v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_DiLepMass[idx_jecup][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1pt[idx_jecup][6] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1eta[idx_jecup][6], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1phi[idx_jecup][6], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2pt[idx_jecup][6] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2eta[idx_jecup][6], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2phi[idx_jecup][6], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1pt[idx_jecup][6] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1eta[idx_jecup][6], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1phi[idx_jecup][6], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2pt[idx_jecup][6] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2eta[idx_jecup][6], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2phi[idx_jecup][6], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_Jets[idx_jecup][6], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METpt[idx_jecup][6]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METphi[idx_jecup][6]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );

                  }//// JES Up  2 or more b-tagging  
                  ////
                  
                  if (v_bjetup_idx.size() ==2)// exactly 2b-tagging ..
                  {


                     FillHisto( h_cf_sys_NLeptons[idx_jecup][7], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecup][7] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecup][7], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecup][7], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecup][7] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecup][7], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecup][7], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NPV[idx_jecup][7]    , num_pv        , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NJets[idx_jecup][7]  , v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecup][7] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecup][7], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecup][7], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecup][7] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecup][7], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecup][7], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][7] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][7],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][7] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][7], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_PV[idx_jecup][7], num_pv, v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_DiLepMass[idx_jecup][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1pt[idx_jecup][7] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1eta[idx_jecup][7], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1phi[idx_jecup][7], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2pt[idx_jecup][7] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2eta[idx_jecup][7], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2phi[idx_jecup][7], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1pt[idx_jecup][7] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1eta[idx_jecup][7], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1phi[idx_jecup][7], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2pt[idx_jecup][7] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2eta[idx_jecup][7], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2phi[idx_jecup][7], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_Jets[idx_jecup][7], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METpt[idx_jecup][7]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METphi[idx_jecup][7]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );
                  }// JES Up  exactly 2b-tagging ..
                  /// Top Reconstruction for JES Up 
                  //KinSolSys(Lep,AnLep,bJet1Up,bJet2Up,MetJESUp,"JetEnUp");
                  SetUpKINObsSyst(v_jetup_idx, v_jetup_TL ,MetJESUp);
                  //if ( ksolweight_ != -1 )
                  if (isKinSol)
                  {
                     if ( Top->Pt() > AnTop->Pt() ) { (*Top1) = (*Top); (*Top2) = (*AnTop); }
                     else { (*Top1) = (*AnTop); (*Top2) = (*Top); }
                     FillHisto( h_cf_sys_NLeptons[idx_jecup][8], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecup][8] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecup][8], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecup][8], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecup][8] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecup][8], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecup][8], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NPV[idx_jecup][8]    , num_pv        , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NJets[idx_jecup][8]  , v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecup][8] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecup][8], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecup][8], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecup][8] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecup][8], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecup][8], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][8] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][8],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][8] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][8], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_PV[idx_jecup][8], num_pv, v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_DiLepMass[idx_jecup][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1pt[idx_jecup][8] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1eta[idx_jecup][8], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1phi[idx_jecup][8], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2pt[idx_jecup][8] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2eta[idx_jecup][8], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2phi[idx_jecup][8], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1pt[idx_jecup][8] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1eta[idx_jecup][8], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1phi[idx_jecup][8], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2pt[idx_jecup][8] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2eta[idx_jecup][8], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2phi[idx_jecup][8], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_Jets[idx_jecup][8], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METpt[idx_jecup][8]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METphi[idx_jecup][8]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );

                     FillHisto( h_sys_Top1Mass_[idx_jecup]    , Top1->M()        , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top1pt_[idx_jecup]      , Top1->Pt()       , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top1phi_[idx_jecup]     , Top1->Phi()      , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top1Rapidity_[idx_jecup], Top1->Rapidity() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top1Energy_[idx_jecup]  , Top1->Energy()   , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2Mass_[idx_jecup]    , Top2->M()        , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2pt_[idx_jecup]      , Top2->Pt()       , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2phi_[idx_jecup]     , Top2->Phi()      , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2Rapidity_[idx_jecup], Top2->Rapidity() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2Energy_[idx_jecup]  , Top2->Energy()   , v_SystEvt[idx_jecup] );

                     FillHisto( h_sys_TopMass_[idx_jecup]      , Top->M()         , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Toppt_[idx_jecup]        , Top->Pt()        , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Topphi_[idx_jecup]       , Top->Phi()       , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_TopRapidity_[idx_jecup]  , Top->Rapidity()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_TopEnergy_[idx_jecup]    , Top->Energy()    , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnTopMass_[idx_jecup]    , AnTop->M()       , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnToppt_[idx_jecup]      , AnTop->Pt()      , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnTopphi_[idx_jecup]     , AnTop->Phi()     , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnTopRapidity_[idx_jecup], AnTop->Rapidity(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnTopEnergy_[idx_jecup]  , AnTop->Energy()  , v_SystEvt[idx_jecup] );



                     v_recocp_O.clear();
                     v_recocp_O.push_back( ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO2Vari( Top, AnTop, bJet, AnbJet ) );
                     v_recocp_O.push_back( ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO4Vari( AnbJet, bJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO6Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO7Vari( Top , AnTop, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO8Vari( Top, AnTop, bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO9Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO10Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO11Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO12Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO13Vari( bJet , AnbJet, AnLep, Lep )  );
                     
                     for (int j = 0; j < v_recocp_O.size(); ++ j )
                     {
                        FillHisto( h_sys_Reco_CPO_[idx_jecup][j], v_recocp_O[j] , v_SystEvt[idx_jecup]   );
                        FillHisto( h_sys_Reco_CPO_ReRange_[idx_jecup][j], v_recocp_O[j] , v_SystEvt[idx_jecup]  );
                     }
                  } // Kinmatic Solver //
               }// JES Up  one or more b-tagging 
            }// JES Up MET Cut   
         } // JES Up Num Jet Cut
         ////////////////
         /// JES Down ///
         ////////////////
         if ( NumJetCut(v_jetdn_idx) == true && isAllSyst == true && idx_jecdn != -2) 
         {


            FillHisto( h_cf_sys_NLeptons[idx_jecdn][3], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep1pt[idx_jecdn][3] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep1eta[idx_jecdn][3], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep1phi[idx_jecdn][3], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep2pt[idx_jecdn][3] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep2eta[idx_jecdn][3], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep2phi[idx_jecdn][3], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_NPV[idx_jecdn][3]    , num_pv        , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_NJets[idx_jecdn][3]  , v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet1pt[idx_jecdn][3] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet1eta[idx_jecdn][3], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet1phi[idx_jecdn][3], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet2pt[idx_jecdn][3] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet2eta[idx_jecdn][3], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet2phi[idx_jecdn][3], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_metpt[idx_jecdn][3] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_metphi[idx_jecdn][3],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
            
            FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_metpt[idx_jecdn][3] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_metphi[idx_jecdn][3], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
            
            FillHisto( h_sys_Num_PV[idx_jecdn][3], num_pv, v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_DiLepMass[idx_jecdn][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep1pt[idx_jecdn][3] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep1eta[idx_jecdn][3], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep1phi[idx_jecdn][3], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep2pt[idx_jecdn][3] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep2eta[idx_jecdn][3], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep2phi[idx_jecdn][3], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet1pt[idx_jecdn][3] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet1eta[idx_jecdn][3], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet1phi[idx_jecdn][3], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet2pt[idx_jecdn][3] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet2eta[idx_jecdn][3], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet2phi[idx_jecdn][3], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Num_Jets[idx_jecdn][3], v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_METpt[idx_jecdn][3]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_METphi[idx_jecdn][3]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );

            if ( METCut(MetJESDn) == true )
            {

               FillHisto( h_cf_sys_NLeptons[idx_jecdn][4], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep1pt[idx_jecdn][4] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep1eta[idx_jecdn][4], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep1phi[idx_jecdn][4], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep2pt[idx_jecdn][4] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep2eta[idx_jecdn][4], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep2phi[idx_jecdn][4], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_NPV[idx_jecdn][4]    , num_pv        , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_NJets[idx_jecdn][4]  , v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet1pt[idx_jecdn][4] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet1eta[idx_jecdn][4], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet1phi[idx_jecdn][4], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet2pt[idx_jecdn][4] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet2eta[idx_jecdn][4], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet2phi[idx_jecdn][4], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_metpt[idx_jecdn][4] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_metphi[idx_jecdn][4],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_metpt[idx_jecdn][4] , Met->Pt() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_metphi[idx_jecdn][4], Met->Phi(), v_SystEvt[idx_jecdn] );
               
               FillHisto( h_sys_Num_PV[idx_jecdn][4], num_pv, v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_DiLepMass[idx_jecdn][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep1pt[idx_jecdn][4] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep1eta[idx_jecdn][4], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep1phi[idx_jecdn][4], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep2pt[idx_jecdn][4] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep2eta[idx_jecdn][4], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep2phi[idx_jecdn][4], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet1pt[idx_jecdn][4] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet1eta[idx_jecdn][4], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet1phi[idx_jecdn][4], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet2pt[idx_jecdn][4] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet2eta[idx_jecdn][4], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet2phi[idx_jecdn][4], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Num_Jets[idx_jecdn][4], v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_METpt[idx_jecdn][4]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_METphi[idx_jecdn][4]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );
               /// Apply BTagging Scale Factor for JES Dn
               //BTaggigSFApplyJESR("JetEnDown");
               BTaggigSFApplyJESR(v_jetdn_idx, v_jetdn_TL, idx_jecdn);
               if ( BJetCut(v_bjetdn_idx) == true ) // one or more b-tagging 
               {

                  FillHisto( h_cf_sys_NLeptons[idx_jecdn][5], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep1pt[idx_jecdn][5] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep1eta[idx_jecdn][5], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep1phi[idx_jecdn][5], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep2pt[idx_jecdn][5] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep2eta[idx_jecdn][5], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep2phi[idx_jecdn][5], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_NPV[idx_jecdn][5]    , num_pv        , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_NJets[idx_jecdn][5]  , v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet1pt[idx_jecdn][5] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet1eta[idx_jecdn][5], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet1phi[idx_jecdn][5], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet2pt[idx_jecdn][5] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet2eta[idx_jecdn][5], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet2phi[idx_jecdn][5], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_metpt[idx_jecdn][5] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_metphi[idx_jecdn][5],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
                  
                  FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_metpt[idx_jecdn][5] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_metphi[idx_jecdn][5], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
                  
                  FillHisto( h_sys_Num_PV[idx_jecdn][5], num_pv, v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_DiLepMass[idx_jecdn][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep1pt[idx_jecdn][5] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep1eta[idx_jecdn][5], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep1phi[idx_jecdn][5], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep2pt[idx_jecdn][5] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep2eta[idx_jecdn][5], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep2phi[idx_jecdn][5], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet1pt[idx_jecdn][5] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet1eta[idx_jecdn][5], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet1phi[idx_jecdn][5], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet2pt[idx_jecdn][5] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet2eta[idx_jecdn][5], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet2phi[idx_jecdn][5], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Num_Jets[idx_jecdn][5], v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_METpt[idx_jecdn][5]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_METphi[idx_jecdn][5]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );
                  ///////////////////////////////////////
                  /// 2 or More B-Tagging Requriement ///
                  /////////////////////////////////////// 
                  if ( DoubleBtag(v_bjetdn_idx) == true )
                  {
                     FillHisto( h_cf_sys_NLeptons[idx_jecdn][6], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecdn][6] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecdn][6], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecdn][6], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecdn][6] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecdn][6], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecdn][6], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NPV[idx_jecdn][6]    , num_pv        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NJets[idx_jecdn][6]  , v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecdn][6] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecdn][6], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecdn][6], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecdn][6] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecdn][6], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecdn][6], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][6] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][6],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
                     
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][6] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][6], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
                     
                     FillHisto( h_sys_Num_PV[idx_jecdn][6], num_pv, v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_DiLepMass[idx_jecdn][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1pt[idx_jecdn][6] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1eta[idx_jecdn][6], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1phi[idx_jecdn][6], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2pt[idx_jecdn][6] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2eta[idx_jecdn][6], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2phi[idx_jecdn][6], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1pt[idx_jecdn][6] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1eta[idx_jecdn][6], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1phi[idx_jecdn][6], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2pt[idx_jecdn][6] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2eta[idx_jecdn][6], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2phi[idx_jecdn][6], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_Jets[idx_jecdn][6], v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METpt[idx_jecdn][6]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METphi[idx_jecdn][6]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );

                  }//// JES Down  2 or more b-tagging  
                  if (v_bjetdn_idx.size() ==2)// exactly 2b-tagging ..
                  {

                     FillHisto( h_cf_sys_NLeptons[idx_jecdn][7], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecdn][7] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecdn][7], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecdn][7], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecdn][7] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecdn][7], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecdn][7], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NPV[idx_jecdn][7]    , num_pv        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NJets[idx_jecdn][7]  , v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecdn][7] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecdn][7], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecdn][7], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecdn][7] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecdn][7], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecdn][7], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][7] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][7],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][7] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][7], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_PV[idx_jecdn][7], num_pv, v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_DiLepMass[idx_jecdn][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1pt[idx_jecdn][7] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1eta[idx_jecdn][7], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1phi[idx_jecdn][7], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2pt[idx_jecdn][7] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2eta[idx_jecdn][7], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2phi[idx_jecdn][7], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1pt[idx_jecdn][7] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1eta[idx_jecdn][7], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1phi[idx_jecdn][7], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2pt[idx_jecdn][7] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2eta[idx_jecdn][7], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2phi[idx_jecdn][7], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_Jets[idx_jecdn][7], v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METpt[idx_jecdn][7]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METphi[idx_jecdn][7]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );
                  }// JES Dn  exactly 2b-tagging ..
               
                  /// Top Reconstruction for JES Down 
                  //KinSolSys(Lep,AnLep,bJet1Dn,bJet2Dn,MetJESDn,"JetEnDown");
                  SetUpKINObsSyst(v_jetdn_idx, v_jetdn_TL ,MetJESDn);
                  //if ( ksolweight_ != -1 )
                  if ( isKinSol )
                  {
                     if ( Top->Pt() > AnTop->Pt() ) { (*Top1) = (*Top); (*Top2) = (*AnTop); }
                     else { (*Top1) = (*AnTop); (*Top2) = (*Top); }


                     FillHisto( h_cf_sys_NLeptons[idx_jecdn][8], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecdn][8] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecdn][8], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecdn][8], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecdn][8] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecdn][8], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecdn][8], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NPV[idx_jecdn][8]    , num_pv        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NJets[idx_jecdn][8]  , v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecdn][8] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecdn][8], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecdn][8], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecdn][8] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecdn][8], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecdn][8], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][8] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][8],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][8] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][8], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_PV[idx_jecdn][8], num_pv, v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_DiLepMass[idx_jecdn][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1pt[idx_jecdn][8] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1eta[idx_jecdn][8], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1phi[idx_jecdn][8], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2pt[idx_jecdn][8] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2eta[idx_jecdn][8], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2phi[idx_jecdn][8], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1pt[idx_jecdn][8] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1eta[idx_jecdn][8], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1phi[idx_jecdn][8], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2pt[idx_jecdn][8] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2eta[idx_jecdn][8], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2phi[idx_jecdn][8], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_Jets[idx_jecdn][8], v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METpt[idx_jecdn][8]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METphi[idx_jecdn][8]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );

                     FillHisto( h_sys_Top1Mass_[idx_jecdn]    , Top1->M()        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top1pt_[idx_jecdn]      , Top1->Pt()       , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top1phi_[idx_jecdn]     , Top1->Phi()      , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top1Rapidity_[idx_jecdn], Top1->Rapidity() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top1Energy_[idx_jecdn]  , Top1->Energy()   , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2Mass_[idx_jecdn]    , Top2->M()        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2pt_[idx_jecdn]      , Top2->Pt()       , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2phi_[idx_jecdn]     , Top2->Phi()      , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2Rapidity_[idx_jecdn], Top2->Rapidity() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2Energy_[idx_jecdn]  , Top2->Energy()   , v_SystEvt[idx_jecdn] );

                     FillHisto( h_sys_TopMass_[idx_jecdn]      , Top->M()         , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Toppt_[idx_jecdn]        , Top->Pt()        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Topphi_[idx_jecdn]       , Top->Phi()       , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_TopRapidity_[idx_jecdn]  , Top->Rapidity()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_TopEnergy_[idx_jecdn]    , Top->Energy()    , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnTopMass_[idx_jecdn]    , AnTop->M()       , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnToppt_[idx_jecdn]      , AnTop->Pt()      , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnTopphi_[idx_jecdn]     , AnTop->Phi()     , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnTopRapidity_[idx_jecdn], AnTop->Rapidity(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnTopEnergy_[idx_jecdn]  , AnTop->Energy()  , v_SystEvt[idx_jecdn] );

                     v_recocp_O.clear();
                     v_recocp_O.push_back( ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO2Vari( Top, AnTop, bJet, AnbJet ) );
                     v_recocp_O.push_back( ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO4Vari( AnbJet, bJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO6Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO7Vari( Top , AnTop, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO8Vari( Top, AnTop, bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO9Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO10Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO11Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO12Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO13Vari( bJet , AnbJet, AnLep, Lep )  );
                     for (int j = 0; j < v_recocp_O.size(); ++ j )
                     {
                        FillHisto( h_sys_Reco_CPO_[idx_jecdn][j], v_recocp_O[j] , v_SystEvt[idx_jecdn]   );
                        FillHisto( h_sys_Reco_CPO_ReRange_[idx_jecdn][j], v_recocp_O[j] , v_SystEvt[idx_jecdn]  );
                     }
                  } // Kinmatic Solver //
               }// Num bJet Cut // one or more b-tagging //
            } // JES Down MET Cut // 
         } // JES Down Num Jet Cut 

         ///////////////////////////////
         /// JER Up Systematic Study ///
         ///////////////////////////////
         if ( NumJetCut(v_jetresup_idx) == true && isAllSyst == true && idx_jerup != -2 ) 
         {


            FillHisto( h_cf_sys_NLeptons[idx_jerup][3], v_lepton_idx.size(), v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Lep1pt[idx_jerup][3] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Lep1eta[idx_jerup][3], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Lep1phi[idx_jerup][3], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Lep2pt[idx_jerup][3] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Lep2eta[idx_jerup][3], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Lep2phi[idx_jerup][3], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_NPV[idx_jerup][3]    , num_pv        , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_NJets[idx_jerup][3]  , v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Jet1pt[idx_jerup][3] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Jet1eta[idx_jerup][3], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Jet1phi[idx_jerup][3], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Jet2pt[idx_jerup][3] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Jet2eta[idx_jerup][3], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_Jet2phi[idx_jerup][3], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_metpt[idx_jerup][3] ,MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_metphi[idx_jerup][3],MetJERUp->Phi() , v_SystEvt[idx_jerup] );
            
            FillHisto( h_cf_sys_dilep_inv_mass[idx_jerup][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_metpt[idx_jerup][3] , MetJERUp->Pt() , v_SystEvt[idx_jerup] );
            FillHisto( h_cf_sys_metphi[idx_jerup][3], MetJERUp->Phi(), v_SystEvt[idx_jerup] );
            
            FillHisto( h_sys_Num_PV[idx_jerup][3], num_pv, v_SystEvt[idx_jerup] );
            FillHisto( h_sys_DiLepMass[idx_jerup][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Lep1pt[idx_jerup][3] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Lep1eta[idx_jerup][3], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Lep1phi[idx_jerup][3], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Lep2pt[idx_jerup][3] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Lep2eta[idx_jerup][3], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Lep2phi[idx_jerup][3], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Jet1pt[idx_jerup][3] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Jet1eta[idx_jerup][3], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Jet1phi[idx_jerup][3], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Jet2pt[idx_jerup][3] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Jet2eta[idx_jerup][3], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Jet2phi[idx_jerup][3], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_Num_Jets[idx_jerup][3], v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
            FillHisto( h_sys_METpt[idx_jerup][3]   , MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
            FillHisto( h_sys_METphi[idx_jerup][3]  , MetJERUp->Phi()  , v_SystEvt[idx_jerup] );
            if ( METCut(MetJERUp) == true ) 
            {


               FillHisto( h_cf_sys_NLeptons[idx_jerup][4], v_lepton_idx.size(), v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Lep1pt[idx_jerup][4] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Lep1eta[idx_jerup][4], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Lep1phi[idx_jerup][4], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Lep2pt[idx_jerup][4] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Lep2eta[idx_jerup][4], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Lep2phi[idx_jerup][4], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_NPV[idx_jerup][4]    , num_pv        , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_NJets[idx_jerup][4]  , v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Jet1pt[idx_jerup][4] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Jet1eta[idx_jerup][4], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Jet1phi[idx_jerup][4], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Jet2pt[idx_jerup][4] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Jet2eta[idx_jerup][4], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_Jet2phi[idx_jerup][4], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_metpt[idx_jerup][4] ,MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_metphi[idx_jerup][4],MetJERUp->Phi() , v_SystEvt[idx_jerup] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[idx_jerup][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_metpt[idx_jerup][4] , Met->Pt() , v_SystEvt[idx_jerup] );
               FillHisto( h_cf_sys_metphi[idx_jerup][4], Met->Phi(), v_SystEvt[idx_jerup] );
               
               FillHisto( h_sys_Num_PV[idx_jerup][4], num_pv, v_SystEvt[idx_jerup] );
               FillHisto( h_sys_DiLepMass[idx_jerup][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Lep1pt[idx_jerup][4] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Lep1eta[idx_jerup][4], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Lep1phi[idx_jerup][4], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Lep2pt[idx_jerup][4] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Lep2eta[idx_jerup][4], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Lep2phi[idx_jerup][4], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Jet1pt[idx_jerup][4] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Jet1eta[idx_jerup][4], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Jet1phi[idx_jerup][4], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Jet2pt[idx_jerup][4] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Jet2eta[idx_jerup][4], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Jet2phi[idx_jerup][4], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_Num_Jets[idx_jerup][4], v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
               FillHisto( h_sys_METpt[idx_jerup][4]   , MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
               FillHisto( h_sys_METphi[idx_jerup][4]  , MetJERUp->Phi()  , v_SystEvt[idx_jerup] );
               
               /// Apply BTagging Scale Factor for JER Up
               //BTaggigSFApplyJESR("JetEnUp");
               BTaggigSFApplyJESR(v_jetresup_idx, v_jetresup_TL, idx_jerup);
               if ( BJetCut(v_bjetresup_idx) == true ) // one or more b-tagging 
               {


                  FillHisto( h_cf_sys_NLeptons[idx_jerup][5], v_lepton_idx.size(), v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Lep1pt[idx_jerup][5] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Lep1eta[idx_jerup][5], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Lep1phi[idx_jerup][5], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Lep2pt[idx_jerup][5] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Lep2eta[idx_jerup][5], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Lep2phi[idx_jerup][5], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_NPV[idx_jerup][5]    , num_pv        , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_NJets[idx_jerup][5]  , v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Jet1pt[idx_jerup][5] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Jet1eta[idx_jerup][5], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Jet1phi[idx_jerup][5], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Jet2pt[idx_jerup][5] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Jet2eta[idx_jerup][5], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_Jet2phi[idx_jerup][5], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_metpt[idx_jerup][5] ,MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_metphi[idx_jerup][5],MetJERUp->Phi() , v_SystEvt[idx_jerup] );
                  
                  FillHisto( h_cf_sys_dilep_inv_mass[idx_jerup][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_metpt[idx_jerup][5] , MetJERUp->Pt() , v_SystEvt[idx_jerup] );
                  FillHisto( h_cf_sys_metphi[idx_jerup][5], MetJERUp->Phi(), v_SystEvt[idx_jerup] );
                  
                  FillHisto( h_sys_Num_PV[idx_jerup][5], num_pv, v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_DiLepMass[idx_jerup][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Lep1pt[idx_jerup][5] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Lep1eta[idx_jerup][5], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Lep1phi[idx_jerup][5], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Lep2pt[idx_jerup][5] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Lep2eta[idx_jerup][5], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Lep2phi[idx_jerup][5], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Jet1pt[idx_jerup][5] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Jet1eta[idx_jerup][5], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Jet1phi[idx_jerup][5], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Jet2pt[idx_jerup][5] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Jet2eta[idx_jerup][5], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Jet2phi[idx_jerup][5], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_Num_Jets[idx_jerup][5], v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_METpt[idx_jerup][5]   , MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
                  FillHisto( h_sys_METphi[idx_jerup][5]  , MetJERUp->Phi()  , v_SystEvt[idx_jerup] );

                  ///////////////////////////////////////
                  /// 2 or More B-Tagging Requriement ///
                  /////////////////////////////////////// 
                  if ( DoubleBtag(v_bjetresup_idx) == true )
                  {

                     FillHisto( h_cf_sys_NLeptons[idx_jerup][6], v_lepton_idx.size(), v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jerup][6] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jerup][6], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jerup][6], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jerup][6] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jerup][6], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jerup][6], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_NPV[idx_jerup][6]    , num_pv        , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_NJets[idx_jerup][6]  , v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jerup][6] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jerup][6], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jerup][6], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jerup][6] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jerup][6], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jerup][6], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metpt[idx_jerup][6] ,MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metphi[idx_jerup][6],MetJERUp->Phi() , v_SystEvt[idx_jerup] );
                     
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jerup][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metpt[idx_jerup][6] , MetJERUp->Pt() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metphi[idx_jerup][6], MetJERUp->Phi(), v_SystEvt[idx_jerup] );
                     
                     FillHisto( h_sys_Num_PV[idx_jerup][6], num_pv, v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_DiLepMass[idx_jerup][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep1pt[idx_jerup][6] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep1eta[idx_jerup][6], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep1phi[idx_jerup][6], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep2pt[idx_jerup][6] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep2eta[idx_jerup][6], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep2phi[idx_jerup][6], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet1pt[idx_jerup][6] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet1eta[idx_jerup][6], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet1phi[idx_jerup][6], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet2pt[idx_jerup][6] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet2eta[idx_jerup][6], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet2phi[idx_jerup][6], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Num_Jets[idx_jerup][6], v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_METpt[idx_jerup][6]   , MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_METphi[idx_jerup][6]  , MetJERUp->Phi()  , v_SystEvt[idx_jerup] );

                  }//// JER Up  2 or more b-tagging  
                  ////
                  
                  if (v_bjetresup_idx.size() ==2)// exactly 2b-tagging ..
                  {
                     FillHisto( h_cf_sys_NLeptons[idx_jerup][7], v_lepton_idx.size(), v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jerup][7] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jerup][7], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jerup][7], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jerup][7] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jerup][7], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jerup][7], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_NPV[idx_jerup][7]    , num_pv        , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_NJets[idx_jerup][7]  , v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jerup][7] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jerup][7], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jerup][7], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jerup][7] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jerup][7], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jerup][7], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metpt[idx_jerup][7] ,MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metphi[idx_jerup][7],MetJERUp->Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jerup][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metpt[idx_jerup][7] , MetJERUp->Pt() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metphi[idx_jerup][7], MetJERUp->Phi(), v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Num_PV[idx_jerup][7], num_pv, v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_DiLepMass[idx_jerup][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep1pt[idx_jerup][7] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep1eta[idx_jerup][7], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep1phi[idx_jerup][7], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep2pt[idx_jerup][7] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep2eta[idx_jerup][7], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep2phi[idx_jerup][7], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet1pt[idx_jerup][7] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet1eta[idx_jerup][7], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet1phi[idx_jerup][7], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet2pt[idx_jerup][7] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet2eta[idx_jerup][7], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet2phi[idx_jerup][7], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Num_Jets[idx_jerup][7], v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_METpt[idx_jerup][7]   , MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_METphi[idx_jerup][7]  , MetJERUp->Phi()  , v_SystEvt[idx_jerup] );
                  }// JER Up  exactly 2b-tagging ..
                  /// Top Reconstruction for JER Up 
                  //KinSolSys(Lep,AnLep,bJet1JERUp,bJet2JERUp,MetJERUp,"JetResUp");
                  SetUpKINObsSyst(v_jetresup_idx, v_jetresup_TL ,MetJERUp);
                  //if ( ksolweight_ != -1 )
                  if ( isKinSol )
                  {
                     if ( Top->Pt() > AnTop->Pt() ) { (*Top1) = (*Top); (*Top2) = (*AnTop); }
                     else { (*Top1) = (*AnTop); (*Top2) = (*Top); }

                     FillHisto( h_cf_sys_NLeptons[idx_jerup][8], v_lepton_idx.size(), v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jerup][8] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jerup][8], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jerup][8], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jerup][8] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jerup][8], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jerup][8], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_NPV[idx_jerup][8]    , num_pv        , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_NJets[idx_jerup][8]  , v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jerup][8] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jerup][8], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jerup][8], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jerup][8] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jerup][8], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jerup][8], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metpt[idx_jerup][8] ,MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metphi[idx_jerup][8],MetJERUp->Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jerup][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metpt[idx_jerup][8] , MetJERUp->Pt() , v_SystEvt[idx_jerup] );
                     FillHisto( h_cf_sys_metphi[idx_jerup][8], MetJERUp->Phi(), v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Num_PV[idx_jerup][8], num_pv, v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_DiLepMass[idx_jerup][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep1pt[idx_jerup][8] , (*Lep1).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep1eta[idx_jerup][8], (*Lep1).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep1phi[idx_jerup][8], (*Lep1).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep2pt[idx_jerup][8] , (*Lep2).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep2eta[idx_jerup][8], (*Lep2).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Lep2phi[idx_jerup][8], (*Lep2).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet1pt[idx_jerup][8] , (*Jet1JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet1eta[idx_jerup][8], (*Jet1JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet1phi[idx_jerup][8], (*Jet1JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet2pt[idx_jerup][8] , (*Jet2JERUp).Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet2eta[idx_jerup][8], (*Jet2JERUp).Eta() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Jet2phi[idx_jerup][8], (*Jet2JERUp).Phi() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Num_Jets[idx_jerup][8], v_jetresup_idx.size(), v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_METpt[idx_jerup][8]   , MetJERUp->Pt()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_METphi[idx_jerup][8]  , MetJERUp->Phi()  , v_SystEvt[idx_jerup] );

                     FillHisto( h_sys_Top1Mass_[idx_jerup]    , Top1->M()        , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Top1pt_[idx_jerup]      , Top1->Pt()       , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Top1phi_[idx_jerup]     , Top1->Phi()      , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Top1Rapidity_[idx_jerup], Top1->Rapidity() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Top1Energy_[idx_jerup]  , Top1->Energy()   , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Top2Mass_[idx_jerup]    , Top2->M()        , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Top2pt_[idx_jerup]      , Top2->Pt()       , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Top2phi_[idx_jerup]     , Top2->Phi()      , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Top2Rapidity_[idx_jerup], Top2->Rapidity() , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Top2Energy_[idx_jerup]  , Top2->Energy()   , v_SystEvt[idx_jerup] );

                     FillHisto( h_sys_TopMass_[idx_jerup]      , Top->M()         , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Toppt_[idx_jerup]        , Top->Pt()        , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_Topphi_[idx_jerup]       , Top->Phi()       , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_TopRapidity_[idx_jerup]  , Top->Rapidity()  , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_TopEnergy_[idx_jerup]    , Top->Energy()    , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_AnTopMass_[idx_jerup]    , AnTop->M()       , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_AnToppt_[idx_jerup]      , AnTop->Pt()      , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_AnTopphi_[idx_jerup]     , AnTop->Phi()     , v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_AnTopRapidity_[idx_jerup], AnTop->Rapidity(), v_SystEvt[idx_jerup] );
                     FillHisto( h_sys_AnTopEnergy_[idx_jerup]  , AnTop->Energy()  , v_SystEvt[idx_jerup] );

                     v_recocp_O.clear();
                     v_recocp_O.push_back( ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO2Vari( Top, AnTop, bJet, AnbJet ) );
                     v_recocp_O.push_back( ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO4Vari( AnbJet, bJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO6Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO7Vari( Top , AnTop, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO8Vari( Top, AnTop, bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO9Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO10Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO11Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO12Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO13Vari( bJet , AnbJet, AnLep, Lep )  );
                     
                     for (int j = 0; j < v_recocp_O.size(); ++ j )
                     {
                        FillHisto( h_sys_Reco_CPO_[idx_jerup][j], v_recocp_O[j] , v_SystEvt[idx_jerup]   );
                        FillHisto( h_sys_Reco_CPO_ReRange_[idx_jerup][j], v_recocp_O[j] , v_SystEvt[idx_jerup]  );
                     }
                  } // Kinmatic Solver //
               }// JER Up  one or more b-tagging 
            }// JER Up MET Cut   
         } // JER Up Num Jet Cut
         if ( NumJetCut(v_jetresdn_idx) == true && isAllSyst == true && idx_jerdn != -2) 
         {
            FillHisto( h_cf_sys_NLeptons[idx_jerdn][3], v_lepton_idx.size(), v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Lep1pt[idx_jerdn][3] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Lep1eta[idx_jerdn][3], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Lep1phi[idx_jerdn][3], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Lep2pt[idx_jerdn][3] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Lep2eta[idx_jerdn][3], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Lep2phi[idx_jerdn][3], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_NPV[idx_jerdn][3]    , num_pv        , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_NJets[idx_jerdn][3]  , v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Jet1pt[idx_jerdn][3] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Jet1eta[idx_jerdn][3], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Jet1phi[idx_jerdn][3], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Jet2pt[idx_jerdn][3] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Jet2eta[idx_jerdn][3], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_Jet2phi[idx_jerdn][3], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_metpt[idx_jerdn][3] ,MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_metphi[idx_jerdn][3],MetJERDn->Phi() , v_SystEvt[idx_jerdn] );
            
            FillHisto( h_cf_sys_dilep_inv_mass[idx_jerdn][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_metpt[idx_jerdn][3] , MetJERDn->Pt() , v_SystEvt[idx_jerdn] );
            FillHisto( h_cf_sys_metphi[idx_jerdn][3], MetJERDn->Phi(), v_SystEvt[idx_jerdn] );
            
            FillHisto( h_sys_Num_PV[idx_jerdn][3], num_pv, v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_DiLepMass[idx_jerdn][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Lep1pt[idx_jerdn][3] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Lep1eta[idx_jerdn][3], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Lep1phi[idx_jerdn][3], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Lep2pt[idx_jerdn][3] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Lep2eta[idx_jerdn][3], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Lep2phi[idx_jerdn][3], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Jet1pt[idx_jerdn][3] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Jet1eta[idx_jerdn][3], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Jet1phi[idx_jerdn][3], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Jet2pt[idx_jerdn][3] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Jet2eta[idx_jerdn][3], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Jet2phi[idx_jerdn][3], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_Num_Jets[idx_jerdn][3], v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_METpt[idx_jerdn][3]   , MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
            FillHisto( h_sys_METphi[idx_jerdn][3]  , MetJERDn->Phi()  , v_SystEvt[idx_jerdn] );

            if ( METCut(MetJERDn) == true )
            {
               FillHisto( h_cf_sys_NLeptons[idx_jerdn][4], v_lepton_idx.size(), v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Lep1pt[idx_jerdn][4] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Lep1eta[idx_jerdn][4], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Lep1phi[idx_jerdn][4], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Lep2pt[idx_jerdn][4] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Lep2eta[idx_jerdn][4], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Lep2phi[idx_jerdn][4], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_NPV[idx_jerdn][4]    , num_pv        , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_NJets[idx_jerdn][4]  , v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Jet1pt[idx_jerdn][4] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Jet1eta[idx_jerdn][4], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Jet1phi[idx_jerdn][4], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Jet2pt[idx_jerdn][4] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Jet2eta[idx_jerdn][4], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_Jet2phi[idx_jerdn][4], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_metpt[idx_jerdn][4] ,MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_metphi[idx_jerdn][4],MetJERDn->Phi() , v_SystEvt[idx_jerdn] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[idx_jerdn][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_metpt[idx_jerdn][4] , Met->Pt() , v_SystEvt[idx_jerdn] );
               FillHisto( h_cf_sys_metphi[idx_jerdn][4], Met->Phi(), v_SystEvt[idx_jerdn] );
               
               FillHisto( h_sys_Num_PV[idx_jerdn][4], num_pv, v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_DiLepMass[idx_jerdn][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Lep1pt[idx_jerdn][4] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Lep1eta[idx_jerdn][4], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Lep1phi[idx_jerdn][4], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Lep2pt[idx_jerdn][4] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Lep2eta[idx_jerdn][4], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Lep2phi[idx_jerdn][4], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Jet1pt[idx_jerdn][4] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Jet1eta[idx_jerdn][4], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Jet1phi[idx_jerdn][4], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Jet2pt[idx_jerdn][4] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Jet2eta[idx_jerdn][4], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Jet2phi[idx_jerdn][4], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_Num_Jets[idx_jerdn][4], v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_METpt[idx_jerdn][4]   , MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
               FillHisto( h_sys_METphi[idx_jerdn][4]  , MetJERDn->Phi()  , v_SystEvt[idx_jerdn] );
               /// Apply BTagging Scale Factor for JER Dn
               BTaggigSFApplyJESR(v_jetresdn_idx, v_jetresdn_TL, idx_jerdn);
               if ( BJetCut(v_bjetresdn_idx) == true ) // one or more b-tagging 
               {
                  FillHisto( h_cf_sys_NLeptons[idx_jerdn][5], v_lepton_idx.size(), v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Lep1pt[idx_jerdn][5] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Lep1eta[idx_jerdn][5], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Lep1phi[idx_jerdn][5], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Lep2pt[idx_jerdn][5] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Lep2eta[idx_jerdn][5], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Lep2phi[idx_jerdn][5], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_NPV[idx_jerdn][5]    , num_pv        , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_NJets[idx_jerdn][5]  , v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Jet1pt[idx_jerdn][5] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Jet1eta[idx_jerdn][5], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Jet1phi[idx_jerdn][5], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Jet2pt[idx_jerdn][5] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Jet2eta[idx_jerdn][5], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_Jet2phi[idx_jerdn][5], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_metpt[idx_jerdn][5] ,MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_metphi[idx_jerdn][5],MetJERDn->Phi() , v_SystEvt[idx_jerdn] );
                  
                  FillHisto( h_cf_sys_dilep_inv_mass[idx_jerdn][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_metpt[idx_jerdn][5] , MetJERDn->Pt() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_cf_sys_metphi[idx_jerdn][5], MetJERDn->Phi(), v_SystEvt[idx_jerdn] );
                  
                  FillHisto( h_sys_Num_PV[idx_jerdn][5], num_pv, v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_DiLepMass[idx_jerdn][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Lep1pt[idx_jerdn][5] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Lep1eta[idx_jerdn][5], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Lep1phi[idx_jerdn][5], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Lep2pt[idx_jerdn][5] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Lep2eta[idx_jerdn][5], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Lep2phi[idx_jerdn][5], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Jet1pt[idx_jerdn][5] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Jet1eta[idx_jerdn][5], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Jet1phi[idx_jerdn][5], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Jet2pt[idx_jerdn][5] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Jet2eta[idx_jerdn][5], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Jet2phi[idx_jerdn][5], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_Num_Jets[idx_jerdn][5], v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_METpt[idx_jerdn][5]   , MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
                  FillHisto( h_sys_METphi[idx_jerdn][5]  , MetJERDn->Phi()  , v_SystEvt[idx_jerdn] );
                  ///////////////////////////////////////
                  /// 2 or More B-Tagging Requriement ///
                  /////////////////////////////////////// 
                  if ( DoubleBtag(v_bjetresdn_idx) == true )
                  {
                     FillHisto( h_cf_sys_NLeptons[idx_jerdn][6], v_lepton_idx.size(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jerdn][6] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jerdn][6], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jerdn][6], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jerdn][6] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jerdn][6], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jerdn][6], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_NPV[idx_jerdn][6]    , num_pv        , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_NJets[idx_jerdn][6]  , v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jerdn][6] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jerdn][6], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jerdn][6], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jerdn][6] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jerdn][6], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jerdn][6], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metpt[idx_jerdn][6] ,MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metphi[idx_jerdn][6],MetJERDn->Phi() , v_SystEvt[idx_jerdn] );
                     
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jerdn][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metpt[idx_jerdn][6] , MetJERDn->Pt() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metphi[idx_jerdn][6], MetJERDn->Phi(), v_SystEvt[idx_jerdn] );
                     
                     FillHisto( h_sys_Num_PV[idx_jerdn][6], num_pv, v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_DiLepMass[idx_jerdn][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep1pt[idx_jerdn][6] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep1eta[idx_jerdn][6], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep1phi[idx_jerdn][6], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep2pt[idx_jerdn][6] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep2eta[idx_jerdn][6], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep2phi[idx_jerdn][6], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet1pt[idx_jerdn][6] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet1eta[idx_jerdn][6], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet1phi[idx_jerdn][6], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet2pt[idx_jerdn][6] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet2eta[idx_jerdn][6], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet2phi[idx_jerdn][6], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Num_Jets[idx_jerdn][6], v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_METpt[idx_jerdn][6]   , MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_METphi[idx_jerdn][6]  , MetJERDn->Phi()  , v_SystEvt[idx_jerdn] );

                  }//// JER Down  2 or more b-tagging  
                  if (v_bjetresdn_idx.size() ==2)// exactly 2b-tagging ..
                  {
                     FillHisto( h_cf_sys_NLeptons[idx_jerdn][7], v_lepton_idx.size(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jerdn][7] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jerdn][7], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jerdn][7], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jerdn][7] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jerdn][7], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jerdn][7], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_NPV[idx_jerdn][7]    , num_pv        , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_NJets[idx_jerdn][7]  , v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jerdn][7] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jerdn][7], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jerdn][7], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jerdn][7] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jerdn][7], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jerdn][7], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metpt[idx_jerdn][7] ,MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metphi[idx_jerdn][7],MetJERDn->Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jerdn][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metpt[idx_jerdn][7] , MetJERDn->Pt() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metphi[idx_jerdn][7], MetJERDn->Phi(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Num_PV[idx_jerdn][7], num_pv, v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_DiLepMass[idx_jerdn][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep1pt[idx_jerdn][7] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep1eta[idx_jerdn][7], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep1phi[idx_jerdn][7], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep2pt[idx_jerdn][7] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep2eta[idx_jerdn][7], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep2phi[idx_jerdn][7], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet1pt[idx_jerdn][7] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet1eta[idx_jerdn][7], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet1phi[idx_jerdn][7], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet2pt[idx_jerdn][7] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet2eta[idx_jerdn][7], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet2phi[idx_jerdn][7], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Num_Jets[idx_jerdn][7], v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_METpt[idx_jerdn][7]   , MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_METphi[idx_jerdn][7]  , MetJERDn->Phi()  , v_SystEvt[idx_jerdn] );
                  }// JER Dn  exactly 2b-tagging ..
               
                  /// Top Reconstruction for JER Down 
                  //KinSolSys(Lep,AnLep,bJet1JERDn,bJet2JERDn,MetJERDn,"JetResDown");
                  SetUpKINObsSyst(v_jetresdn_idx, v_jetresdn_TL ,MetJERDn);
                  //if ( ksolweight_ != -1 )
                  if ( isKinSol )
                  {
                     if ( Top->Pt() > AnTop->Pt() ) { (*Top1) = (*Top); (*Top2) = (*AnTop); }
                     else { (*Top1) = (*AnTop); (*Top2) = (*Top); }
                     FillHisto( h_cf_sys_NLeptons[idx_jerdn][8], v_lepton_idx.size(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jerdn][8] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jerdn][8], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jerdn][8], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jerdn][8] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jerdn][8], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jerdn][8], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_NPV[idx_jerdn][8]    , num_pv        , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_NJets[idx_jerdn][8]  , v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jerdn][8] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jerdn][8], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jerdn][8], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jerdn][8] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jerdn][8], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jerdn][8], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metpt[idx_jerdn][8] ,MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metphi[idx_jerdn][8],MetJERDn->Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jerdn][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metpt[idx_jerdn][8] , MetJERDn->Pt() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_cf_sys_metphi[idx_jerdn][8], MetJERDn->Phi(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Num_PV[idx_jerdn][8], num_pv, v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_DiLepMass[idx_jerdn][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep1pt[idx_jerdn][8] , (*Lep1).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep1eta[idx_jerdn][8], (*Lep1).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep1phi[idx_jerdn][8], (*Lep1).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep2pt[idx_jerdn][8] , (*Lep2).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep2eta[idx_jerdn][8], (*Lep2).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Lep2phi[idx_jerdn][8], (*Lep2).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet1pt[idx_jerdn][8] , (*Jet1JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet1eta[idx_jerdn][8], (*Jet1JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet1phi[idx_jerdn][8], (*Jet1JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet2pt[idx_jerdn][8] , (*Jet2JERDn).Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet2eta[idx_jerdn][8], (*Jet2JERDn).Eta() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Jet2phi[idx_jerdn][8], (*Jet2JERDn).Phi() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Num_Jets[idx_jerdn][8], v_jetresdn_idx.size(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_METpt[idx_jerdn][8]   , MetJERDn->Pt()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_METphi[idx_jerdn][8]  , MetJERDn->Phi()  , v_SystEvt[idx_jerdn] );

                     FillHisto( h_sys_Top1Mass_[idx_jerdn]    , Top1->M()        , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Top1pt_[idx_jerdn]      , Top1->Pt()       , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Top1phi_[idx_jerdn]     , Top1->Phi()      , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Top1Rapidity_[idx_jerdn], Top1->Rapidity() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Top1Energy_[idx_jerdn]  , Top1->Energy()   , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Top2Mass_[idx_jerdn]    , Top2->M()        , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Top2pt_[idx_jerdn]      , Top2->Pt()       , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Top2phi_[idx_jerdn]     , Top2->Phi()      , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Top2Rapidity_[idx_jerdn], Top2->Rapidity() , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Top2Energy_[idx_jerdn]  , Top2->Energy()   , v_SystEvt[idx_jerdn] );

                     FillHisto( h_sys_TopMass_[idx_jerdn]      , Top->M()         , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Toppt_[idx_jerdn]        , Top->Pt()        , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_Topphi_[idx_jerdn]       , Top->Phi()       , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_TopRapidity_[idx_jerdn]  , Top->Rapidity()  , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_TopEnergy_[idx_jerdn]    , Top->Energy()    , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_AnTopMass_[idx_jerdn]    , AnTop->M()       , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_AnToppt_[idx_jerdn]      , AnTop->Pt()      , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_AnTopphi_[idx_jerdn]     , AnTop->Phi()     , v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_AnTopRapidity_[idx_jerdn], AnTop->Rapidity(), v_SystEvt[idx_jerdn] );
                     FillHisto( h_sys_AnTopEnergy_[idx_jerdn]  , AnTop->Energy()  , v_SystEvt[idx_jerdn] );

                     v_recocp_O.clear();
                     v_recocp_O.push_back( ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO2Vari( Top, AnTop, bJet, AnbJet ) );
                     v_recocp_O.push_back( ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO4Vari( AnbJet, bJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO6Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO7Vari( Top , AnTop, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO8Vari( Top, AnTop, bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO9Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO10Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO11Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO12Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO13Vari( bJet , AnbJet, AnLep, Lep )  );
                     for (int j = 0; j < v_recocp_O.size(); ++ j )
                     {
                        FillHisto( h_sys_Reco_CPO_[idx_jerdn][j], v_recocp_O[j] , v_SystEvt[idx_jerdn]   );
                        FillHisto( h_sys_Reco_CPO_ReRange_[idx_jerdn][j], v_recocp_O[j] , v_SystEvt[idx_jerdn]  );
                     }
                  } // Kinmatic Solver //
               }// Num bJet Cut // one or more b-tagging //
            } // JER Down MET Cut // 
         } // JER Down Num Jet Cut 


         //////////////////////
         /// Fill Histogram ///
         //////////////////////
         //ClearVectors();*/
     
      }//Di-Lepton Analysis//

      ///////////////////////////////////
      /// ** Lepton + Jet Analysis ** ///
      ///////////////////////////////////
      else if ( TString(Decaymode).Contains( "muonJet" )  ) 
      {
         cout << " -- We don't have Muon + Jet code ---" << endl;   
      }
      else {cout << "Decaymode Error : " << Decaymode << endl; }
   }//event loop
  
   printf("Total processed number of events: %lld\n", __tot_evt);
   cout << "Total Event after duplicated event filtering: " << af_dupleevt << endl;
}//end Loop function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ssb_analysis::GetNtupleTotalEvent( unsigned int totevent )// Not Used Function
{
   NtupletotalEvent = totevent;
}

void ssb_analysis::Start( int genLoopon )
{
   if      ( genLoopon == 0 ){ fout = new TFile(Form("output/%s",outfile),"RECREATE");}
   else if      ( genLoopon == 1 ){ fout = new TFile(Form("output/%s",outfile),"UPDATE");}
   //if      ( genLoopon == 0 ){ fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/SSB_CPviolation/output/%s",outfile),"RECREATE");}
   //else if ( genLoopon == 1 ){ fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/SSB_CPviolation/output/%s",outfile),"UPDATE"  );}
   else {cout << "genLoopon error" << endl;}
   fout->cd("");

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();

   /////////////////////////////////
   /// For All-in-One Systematic ///
   /////////////////////////////////
   if (isAllSyst == true)
   {
      for (int i = 0; i < v_SystFullName.size(); ++i )
      {
         gSystem->mkdir(Form("output/%s",v_SystFullName[i].Data()));
         a_fout[i] = new TFile(Form("output/%s/%s",v_SystFullName[i].Data(),v_outName[i].Data()),"RECREATE");
         // a_fout[i] = new TFile(Form("%s",v_outName[i].Data()),"RECREATE");
         a_fout[i]->cd();
         DeclareHistosSyst(i);
      }
   }
   ReadDupleList(); 
}

void ssb_analysis::DeclareHistos()
{

   /// Test For Systematic All-in-One Code ///

   for (int i =0 ; i < 10 ; i++)
   {
      h_cf_NLeptons[i]       = new TH1D(Form("_h_cf_NLeptons_%d_"      , i),Form("Num. Lepton %s"                ,cutflowName[i].Data()), 20  , 0    , 20  ); h_cf_NLeptons[i]->Sumw2();         
      h_cf_Lep1pt[i]         = new TH1D(Form("_h_cf_Lep1pt_%d_"        , i),Form("Leading Lepton pT %s"          ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Lep1pt[i]->Sumw2();           
      h_cf_Lep1phi[i]        = new TH1D(Form("_h_cf_Lep1phi_%d_"       , i),Form("Leading Lepton phi %s"        ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Lep1phi[i]->Sumw2();          
      h_cf_Lep1eta[i]        = new TH1D(Form("_h_cf_Lep1eta_%d_"       , i),Form("Leading Lepton eta %s"        ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Lep1eta[i]->Sumw2();          
      h_cf_Lep2pt[i]         = new TH1D(Form("_h_cf_Lep2pt_%d_"        , i),Form("Second Leading Lepton pT %s"   ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Lep2pt[i]->Sumw2();           
      h_cf_Lep2phi[i]        = new TH1D(Form("_h_cf_Lep2phi_%d_"       , i),Form("Second Leading Lepton phi %s"  ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Lep2phi[i]->Sumw2();      
      h_cf_Lep2eta[i]        = new TH1D(Form("_h_cf_Lep2eta_%d_"       , i),Form("Second Leading Lepton eta %s"  ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Lep2eta[i]->Sumw2();       
      h_cf_dilep_inv_mass[i] = new TH1D(Form("_h_cf_dilep_inv_mass_%d_", i),Form("Invariant Mass of Dilepton %s" ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_dilep_inv_mass[i]->Sumw2();
      h_cf_Jet1pt[i]         = new TH1D(Form("_h_cf_Jet1pt_%d_"        , i),Form("Leading Jet pT %s"             ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Jet1pt[i]->Sumw2();       
      h_cf_Jet1phi[i]        = new TH1D(Form("_h_cf_Jet1phi_%d_"       , i),Form("Leading Jet phi %s"            ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Jet1phi[i]->Sumw2();       
      h_cf_Jet1eta[i]        = new TH1D(Form("_h_cf_Jet1eta_%d_"       , i),Form("Leading Jet eta %s"            ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Jet1eta[i]->Sumw2();       
      h_cf_Jet2pt[i]         = new TH1D(Form("_h_cf_Jet2pt_%d_"        , i),Form("Second Leading Jet pT %s"      ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Jet2pt[i]->Sumw2();        
      h_cf_Jet2phi[i]        = new TH1D(Form("_h_cf_Jet2phi_%d_"       , i),Form("Second Leading Jet phi %s"     ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Jet2phi[i]->Sumw2();       
      h_cf_Jet2eta[i]        = new TH1D(Form("_h_cf_Jet2eta_%d_"       , i),Form("Second Leading Jet eta %s"     ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Jet2eta[i]->Sumw2();      
      h_cf_NJets[i]          = new TH1D(Form("_h_cf_NJets_%d_"         , i),Form("Num. of Jets %s"               ,cutflowName[i].Data()), 20  , 0    , 20  ); h_cf_NJets[i]->Sumw2();        
      h_cf_metpt[i]          = new TH1D(Form("_h_cf_metpt_%d_"         , i),Form("MET %s"                        ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_metpt[i]->Sumw2();        
      h_cf_metphi[i]         = new TH1D(Form("_h_cf_metphi_%d_"        , i),Form("MET Phi %s"                    ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_metphi[i]->Sumw2();       
      h_cf_Nbjets[i]         = new TH1D(Form("_h_cf_Nbjets_%d_"        , i),Form("Num. of b-jets %s"             ,cutflowName[i].Data()), 20  , 0    , 20  ); h_cf_Nbjets[i]->Sumw2();        
      h_cf_NPV[i]            = new TH1D(Form("_h_cf_NPV_%d_"           , i),Form("Num. of Primary Vertex %s"     ,cutflowName[i].Data()), 100 , 0    , 100 ); h_cf_NPV[i]->Sumw2();          


      h_EventWeight[i] = new TH1D(Form("h_EventWeight_%d",i), Form("Event Weight %s",cutflowName[i].Data()), 1000, -5, 5); h_EventWeight[i]->Sumw2();

      h_Lep1pt[i]  = new TH1D(Form("h_Lep1pt_%d" ,i), Form("Leading Lepton pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep1pt[i]->Sumw2(); 
      h_Lep2pt[i]  = new TH1D(Form("h_Lep2pt_%d" ,i), Form("Second Leading Lepton pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep2pt[i]->Sumw2();
      h_Lep1eta[i] = new TH1D(Form("h_Lep1eta_%d",i), Form("Leading Lepton Eta    %s"    ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep1eta[i]->Sumw2();
      h_Lep2eta[i] = new TH1D(Form("h_Lep2eta_%d",i), Form("Second Leading Lepton Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep2eta[i]->Sumw2();
      h_Lep1phi[i] = new TH1D(Form("h_Lep1phi_%d",i), Form("Leading Lepton Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Lep1phi[i]->Sumw2();
      h_Lep2phi[i] = new TH1D(Form("h_Lep2phi_%d",i), Form("Second Leading Lepton Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Lep2phi[i]->Sumw2();
   
      h_Muonpt[i]  = new TH1D(Form("h_Muonpt_%d",i ), Form("Muon pT %s"         ,cutflowName[i].Data()), 250, 0.0, 250); h_Muonpt[i]->Sumw2();
      h_Elecpt[i]  = new TH1D(Form("h_Elecpt_%d",i ), Form("Electron pT %s"     ,cutflowName[i].Data()), 250, 0.0, 250); h_Elecpt[i]->Sumw2();
      h_Muoneta[i] = new TH1D(Form("h_Muoneta_%d",i), Form("Muon Eta %s"        ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Muoneta[i]->Sumw2();
      h_Eleceta[i] = new TH1D(Form("h_Eleceta_%d",i), Form("Electron Eta %s"    ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Eleceta[i]->Sumw2();
      h_Muonphi[i] = new TH1D(Form("h_Muonphi_%d",i), Form("Muon Phi %s"        ,cutflowName[i].Data()), 24, -1*pi, pi); h_Muonphi[i]->Sumw2();
      h_Elecphi[i] = new TH1D(Form("h_Elecphi_%d",i), Form("Electron Phi %s"    ,cutflowName[i].Data()), 24, -1*pi, pi); h_Elecphi[i]->Sumw2();
   
      h_Jet1pt[i]  = new TH1D(Form("h_Jet1pt_%d", i), Form("Leading Jet pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet1pt[i]->Sumw2();
      h_Jet2pt[i]  = new TH1D(Form("h_Jet2pt_%d", i), Form("Second Leading Jet pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet2pt[i]->Sumw2();
      h_Jet1eta[i] = new TH1D(Form("h_Jet1eta_%d",i), Form("Leading Jet Eta %s"       ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet1eta[i]->Sumw2();
      h_Jet2eta[i] = new TH1D(Form("h_Jet2eta_%d",i), Form("Second Leading Jet Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet2eta[i]->Sumw2();
      h_Jet1phi[i] = new TH1D(Form("h_Jet1phi_%d",i), Form("Leading Jet Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Jet1phi[i]->Sumw2();
      h_Jet2phi[i] = new TH1D(Form("h_Jet2phi_%d",i), Form("Second Leading Jet Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Jet2phi[i]->Sumw2();
    
     
      h_METpt[i]  = new TH1D(Form("h_METpt_%d",i), Form("MET pT %s" ,cutflowName[i].Data()), 200, 0.0, 200); h_METpt[i]->Sumw2();
      h_METphi[i] = new TH1D(Form("h_METphi_%d",i), Form("MET Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_METphi[i]->Sumw2();

      h_DiLepMass[i] = new TH1D(Form("h_DiLepMass_%d",i),Form("Di-Lepton Invariant Mass %s",cutflowName[i].Data()), 300, 0.0, 300); h_DiLepMass[i]->Sumw2();
      h_Num_PV[i]    = new TH1D(Form("h_Num_PV_%d",i),     Form("Num of Primary Vertex after %s",cutflowName[i].Data()), 100, 0.0, 100); h_Num_PV[i]->Sumw2();
      h_Num_Jets[i]  = new TH1D(Form("h_Num_Jets_%d",i), Form("Num. of Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_Jets[i]->Sumw2();
      h_Num_bJets[i] = new TH1D(Form("h_Num_bJets_%d",i),Form("Num. of b Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_bJets[i]->Sumw2();
      
      if(i>3){h_HT[i]      = new TH1D(Form("h_HT_%d",i), Form("Sum . pT of All Jet %s",cutflowName[i].Data()), 1000, 0.0, 1000); h_HT[i]->Sumw2();}
   }
   for (int i =0; i < 5; ++i)
   {
      h_Topo_Apla[i] = new TH1D(Form("h_Topo_Apla_%d",i), Form("Topologycal Vari. Aplanarity after %s",cutflowName[i+5].Data()), 300, -1, 2); h_Topo_Apla[i]->Sumw2();
      h_Topo_Sphe[i] = new TH1D(Form("h_Topo_Sphe_%d",i), Form("Topologycal Vari. Sphericity after %s",cutflowName[i+5].Data()), 300, -1, 2); h_Topo_Sphe[i]->Sumw2();
      h_Topo_Plan[i] = new TH1D(Form("h_Topo_Plan_%d",i), Form("Topologycal Vari. Planarity after %s",cutflowName[i+5].Data()), 300, -1, 2); h_Topo_Plan[i]->Sumw2();
      h_bTagEff[i]  = new TH1D(Form("h_bTagEff_%d",i), Form("BTagging Weight %s",cutflowName[i+5].Data()), 1000, -5, 5); h_bTagEff[i]->Sumw2();
   }
   h_PileUp      = new TH1D(Form("h_PileUp"), Form("hPileUp"), 1000, 0.0, 100); h_PileUp->Sumw2();// For Closure test // 
   h_PileUp_Up   = new TH1D(Form("h_PileUp_Up"), Form("hPileUp_Up"), 1000, 0.0, 100); h_PileUp_Up->Sumw2();// For Closure test // 
   h_PileUp_Down = new TH1D(Form("h_PileUp_Down"), Form("hPileUp_Down"), 1000, 0.0, 100); h_PileUp_Down->Sumw2();// For Closure test // 

   h_NumEl  = new TH1D(Form("h_NumEl" ), Form("hNumEl" ), 10, 0.0, 10); h_NumEl->Sumw2(); 
   h_NumMu  = new TH1D(Form("h_NumMu" ), Form("hNumMu" ), 10, 0.0, 10); h_NumMu->Sumw2(); 

   h_Num_PV_BeforePreSel = new TH1D(Form("h_Num_PV_BeforePreSel"), Form("Num. of Primary Vertex before Pre-Selection"), 100, 0.0, 100); h_Num_PV_BeforePreSel->Sumw2();// For Closure test // 
   h_Num_PV_AfterMetFilter = new TH1D(Form("h_Num_PV_AfterMetFilter"), Form("Num. of Primary Vertex After MetFilter"), 100, 0.0, 100); h_Num_PV_AfterMetFilter->Sumw2();// For Closure test // 
   h_Num_PV_AfterTrigger = new TH1D(Form("h_Num_PV_AfterTrigger"), Form("Num. of Primary Vertex After Trigger"), 100, 0.0, 100); h_Num_PV_AfterTrigger->Sumw2();// For Closure test // 

   h_Top1Mass     = new TH1D(Form("h_Top1Mass"   ), Form("Top1 Mass"   ), 1000, 0.0, 1000); h_Top1Mass->Sumw2(); 
   h_Top1pt       = new TH1D(Form("h_Top1pt"   ), Form("Top1 pt"   ), 1000, 0.0, 1000); h_Top1pt->Sumw2(); 
   h_Top1Rapidity = new TH1D(Form("h_Top1Rapidity"   ), Form("Top1 Rapidity"   ), 100, -5, 5); h_Top1Rapidity->Sumw2(); 
   h_Top1phi      = new TH1D(Form("h_Top1phi"   ), Form("Top1 phi"   ), 24, -1*pi, pi); h_Top1phi->Sumw2(); 
   h_Top1Energy   = new TH1D(Form("h_Top1Energy"   ), Form("Top1 Energy"   ), 1000, 0.0, 1000); h_Top1Energy->Sumw2(); 
   h_Top2Mass     = new TH1D(Form("h_Top2Mass" ), Form("Top2 Mass" ), 1000, 0.0, 1000); h_Top2Mass->Sumw2();
   h_Top2pt       = new TH1D(Form("h_Top2pt"   ), Form("Top2 pt"   ), 1000, 0.0, 1000); h_Top2pt->Sumw2(); 
   h_Top2Rapidity = new TH1D(Form("h_Top2Rapidity"   ), Form("Top2 Rapidity"   ), 100, -5, 5); h_Top2Rapidity->Sumw2(); 
   h_Top2phi      = new TH1D(Form("h_Top2phi"   ), Form("Top2 phi"   ), 24, -1*pi, pi); h_Top2phi->Sumw2(); 
   h_Top2Energy = new TH1D(Form("h_Top2Energy"   ), Form("Top2 Energy"   ), 1000, 0.0, 1000); h_Top2Energy->Sumw2(); 

   h_TopMass       = new TH1D(Form("h_TopMass"   ), Form("Top Mass"   ), 1000, 0.0, 1000); h_TopMass->Sumw2(); 
   h_Toppt         = new TH1D(Form("h_Toppt"   ), Form("Top pt"   ), 1000, 0.0, 1000); h_Toppt->Sumw2(); 
   h_TopRapidity   = new TH1D(Form("h_TopRapidity"   ), Form("Top Rapidity"   ), 100, -5, 5); h_TopRapidity->Sumw2(); 
   h_Topphi        = new TH1D(Form("h_Topphi"   ), Form("Top phi"   ), 24, -1*pi, pi); h_Topphi->Sumw2(); 
   h_TopEnergy     = new TH1D(Form("h_TopEnergy"   ), Form("Top Energy"   ), 1000, 0.0, 1000); h_TopEnergy->Sumw2(); 
   h_AnTopMass     = new TH1D(Form("h_AnTopMass" ), Form("AnTop Mass" ), 1000, 0.0, 1000); h_AnTopMass->Sumw2();
   h_AnToppt       = new TH1D(Form("h_AnToppt"   ), Form("AnTop pt"   ), 1000, 0.0, 1000); h_AnToppt->Sumw2(); 
   h_AnTopRapidity = new TH1D(Form("h_AnTopRapidity"   ), Form("AnTop Rapidity"   ), 100, -5, 5); h_AnTopRapidity->Sumw2(); 
   h_AnTopphi      = new TH1D(Form("h_AnTopphi"   ), Form("AnTop phi"   ), 24, -1*pi, pi); h_AnTopphi->Sumw2(); 
   h_AnTopEnergy   = new TH1D(Form("h_AnTopEnergy"   ), Form("AnTop Energy"   ), 1000, 0.0, 1000); h_AnTopEnergy->Sumw2(); 


   h_W1Mass     = new TH1D(Form("h_W1Mass"  ), Form("W1 Mass" ), 300, 0.0, 300); h_W1Mass->Sumw2();
   h_W2Mass     = new TH1D(Form("h_W2Mass"  ), Form("W2 Mass" ), 300, 0.0, 300); h_W2Mass->Sumw2();

   h_W1Mt     = new TH1D(Form("h_W1Mt"  ), Form("W1 Transverse Mass" ), 300, 0.0, 300); h_W1Mt->Sumw2();
   h_W2Mt     = new TH1D(Form("h_W2Mt"  ), Form("W2 Transverse Mass" ), 300, 0.0, 300); h_W2Mt->Sumw2();

   h_bJet1Energy = new TH1D(Form("h_bJet1Energy" ), Form("Leading bJet Energy"   ), 500, 0.0, 500); h_bJet1Energy->Sumw2(); 
   h_bJet2Energy = new TH1D(Form("h_bJet2Energy" ), Form("Second Leading bJet Energy" ), 500, 0.0, 500); h_bJet2Energy->Sumw2();

   h_bJetEnergy   = new TH1D(Form("h_bJetEnergy" ), Form("bJet Energy"   ), 1000, 0.0, 1000); h_bJetEnergy->Sumw2(); 
   h_AnbJetEnergy = new TH1D(Form("h_AnbJetEnergy" ), Form("b-barJet Energy" ), 1000, 0.0, 1000); h_AnbJetEnergy->Sumw2(); 

   h_Lep1Energy = new TH1D(Form("h_Lep1Energy" ), Form("Leading Lepton Energy"   ), 400, 0.0, 400); h_Lep1Energy->Sumw2(); 
   h_Lep2Energy = new TH1D(Form("h_Lep2Energy" ), Form("Second Leading Lepton Energy" ), 400, 0.0, 400); h_Lep2Energy->Sumw2();

   h_LepEnergy   = new TH1D(Form("h_LepEnergy" ), Form("Lepton Energy"   ), 400, 0.0, 400); h_LepEnergy->Sumw2(); 
   h_AnLepEnergy = new TH1D(Form("h_AnLepEnergy" ), Form("Anti-Lepton Energy" ), 400, 0.0, 400); h_AnLepEnergy->Sumw2();

   h_Nu1Energy = new TH1D(Form("h_Nu1Energy" ), Form("Leading Nuetrino Energy"   ), 400, 0.0, 400); h_Nu1Energy->Sumw2(); 
   h_Nu2Energy = new TH1D(Form("h_Nu2Energy" ), Form("Second Leading Nuetrino Energy" ), 400, 0.0, 400); h_Nu2Energy->Sumw2(); 

   h_NuEnergy   = new TH1D(Form("h_NuEnergy" ), Form("Nuetrino Energy"   ), 400, 0.0, 400); h_NuEnergy->Sumw2(); 
   h_AnNuEnergy = new TH1D(Form("h_AnNuEnergy" ), Form("anti-Nuetrino Energy" ), 400, 0.0, 400); h_AnNuEnergy->Sumw2();

   h_GenNuEnergy   = new TH1D(Form("h_GenNuEnergy" ), Form("Nuetrino Energy At Genrator Level"   ), 2000, 0.0, 2000); h_GenNuEnergy->Sumw2();
   h_GenAnNuEnergy = new TH1D(Form("h_GenAnNuEnergy" ), Form("anti-Nuetrino Energy At Genrator Level" ), 2000, 0.0, 2000); h_GenAnNuEnergy->Sumw2();
   // kinematics of objects after top reco
   h_bJetpt    = new TH1D(Form("_h_bJetpt_")     ,Form("bJet pT")           ,1000,    0.0,    1000);   h_bJetpt->Sumw2(); 
   h_bJeteta   = new TH1D(Form("_h_bJeteta_")    ,Form("bJet Eta")          ,50,      -2.5,   2.5);    h_bJeteta->Sumw2(); 
   h_bJetphi   = new TH1D(Form("_h_bJetphi_")    ,Form("bJet Phi")          ,24,      -1*pi,  pi);     h_bJetphi->Sumw2(); 
   h_AnbJetpt  = new TH1D(Form("_h_AnbJetpt_")   ,Form("b-barJet pT")       ,1000,    0.0,    1000);   h_AnbJetpt->Sumw2();
   h_AnbJeteta = new TH1D(Form("_h_AnbJeteta_")  ,Form("b-barJet Eta")      ,50,      -2.5,   2.5);    h_AnbJeteta->Sumw2();
   h_AnbJetphi = new TH1D(Form("_h_AnbJetphi_")  ,Form("b-barJet Phi")      ,24,      -1*pi,  pi);     h_AnbJetphi->Sumw2();
   h_Leppt     = new TH1D(Form("_h_Leppt_")      ,Form("Lepton pT")         ,1000,    0.0,    1000);   h_Leppt->Sumw2();
   h_Lepeta    = new TH1D(Form("_h_Lepeta_")     ,Form("Lepton Eta")        ,50,      -2.5,   2.5);    h_Lepeta->Sumw2();
   h_Lepphi    = new TH1D(Form("_h_Lepphi_")     ,Form("Lepton Phi")        ,24,      -1*pi,  pi);     h_Lepphi->Sumw2();
   h_AnLeppt   = new TH1D(Form("_h_AnLeppt_")    ,Form("anti-Lepton pT")    ,1000,    0.0,    1000);   h_AnLeppt->Sumw2();
   h_AnLepeta  = new TH1D(Form("_h_AnLepeta_")   ,Form("anti-Lepton Eta")   ,50,      -2.5,   2.5);    h_AnLepeta->Sumw2();
   h_AnLepphi  = new TH1D(Form("_h_AnLepphi_")   ,Form("anti-Lepton Phi")   ,24,      -1*pi,  pi);     h_AnLepphi->Sumw2();

   h_CPO3_reco         = new TH1D(Form("h_CPO3_reco"   ), Form("CPO3_reco"   ), 200, -10, 10); h_CPO3_reco->Sumw2();
   h_CPO3_reco_JPRUp   = new TH1D(Form("h_CPO3_reco_JPRUp"     ), Form("CPO3_reco_JPRUp"     ), 200, -10, 10); h_CPO3_reco_JPRUp->Sumw2();
   h_CPO3_reco_JPRDown = new TH1D(Form("h_CPO3_reco_JPRDown"   ), Form("CPO3_reco_JPRDown"   ), 200, -10, 10); h_CPO3_reco_JPRDown->Sumw2();
   h_CPO3_reco_TopRes  = new TH1D(Form("h_CPO3_reco_TopRes"   ), Form("CPO3_reco in the Top Mass Window"   ), 200, -10, 10); h_CPO3_reco_TopRes->Sumw2();
   h_CPOb_reco         = new TH1D(Form("h_CPOb_reco"   ), Form("CPOb_reco"   ), 200, -10, 10); h_CPOb_reco->Sumw2();
   h_CPO5_reco         = new TH1D(Form("h_CPO5_reco"   ), Form("CPO5_reco"   ), 200, -10, 10); h_CPO5_reco->Sumw2();

   h_CPO3_evtweight  = new TH1D(Form("h_CPO3_evtweight"   ), Form("CPO3_evtweight"   ), 10, 0, 10); h_CPO3_evtweight->Sumw2();
   h_CPOb_evtweight  = new TH1D(Form("h_CPOb_evtweight"   ), Form("CPOb_evtweight"   ), 10, 0, 10); h_CPOb_evtweight->Sumw2();
   h_CPO5_evtweight  = new TH1D(Form("h_CPO5_evtweight"   ), Form("CPO5_evtweight"   ), 10, 0, 10); h_CPO5_evtweight->Sumw2();

   h_BjorkenX1 = new TH1D(Form("h_BjorkenX1" ), Form("Bjorken X1" ), 120, -0.1  , 1.1); h_BjorkenX1->Sumw2();  
   h_BjorkenX2 = new TH1D(Form("h_BjorkenX2" ), Form("Bjorken X2" ), 120, -0.1, 1.1); h_BjorkenX2->Sumw2(); 
   h_BjorkenX3 = new TH1D(Form("h_BjorkenX3" ), Form("Bjorken X3" ), 2000, 0.0  , 2000); h_BjorkenX3->Sumw2();

   ///////////////////////////////////////////////////
   /// New Algorithm O3 without top reconstruction ///
   ///////////////////////////////////////////////////
   // check the number of hieghest pT jet and the ratio that leading b-tagging jet is the same as the highest pT jet
   h_frac_bTagJet_total                = new TH1D("_h_frac_bTagJet_total_"                , "fraction of leading jet"         ,10, 0, 10); h_frac_bTagJet_total->Sumw2();
   h_frac_LeadingbTagJet_is_Leading    = new TH1D("_h_frac_LeadingbTagJet_is_Leading_"    , "fraction of b-tagging jet"       ,2, 0, 2); h_frac_LeadingbTagJet_is_Leading->Sumw2();
   // find hieght pT b tagging jet
   h_leading_bTagJetpt     = new TH1D(Form("_h_leading_bTagJetpt_")     ,Form("highest pT of b tagging jet pT")   ,1000   ,0        ,1000   ); h_leading_bTagJetpt->Sumw2();
   h_leading_bTagJeteta    = new TH1D(Form("_h_leading_bTagJeteta_")    ,Form("highest pT of b tagging jet eta")  ,50     ,-2.5     ,2.5    ); h_leading_bTagJeteta->Sumw2();
   h_leading_bTagJetphi    = new TH1D(Form("_h_leading_bTagJetphi_")    ,Form("highest pT of b tagging jet phi")  ,24     ,-1*pi    ,1*pi   ); h_leading_bTagJetphi->Sumw2();
   h_leading_jetpt         = new TH1D(Form("_h_leading_jetpt_")         ,Form("leading jet pT (af cut)")          ,1000   ,0        ,1000   ); h_leading_jetpt->Sumw2();
   h_leading_jeteta        = new TH1D(Form("_h_leading_jeteta_")        ,Form("leading jet eta (af cut)")         ,50     ,-2.5     ,2.5    ); h_leading_jeteta->Sumw2();
   h_leading_jetphi        = new TH1D(Form("_h_leading_jetphi_")        ,Form("leading jet phi (af cut)")         ,24     ,-1*pi    ,1*pi   ); h_leading_jetphi->Sumw2();
   h_subleading_jetpt      = new TH1D(Form("_h_subleading_jetpt_")      ,Form("subleading jet pT (af cut)")       ,1000   ,0        ,1000   ); h_subleading_jetpt->Sumw2();
   h_subleading_jeteta     = new TH1D(Form("_h_subleading_jeteta_")     ,Form("subleading jet eta (af cut)")      ,50     ,-2.5     ,2.5    ); h_subleading_jeteta->Sumw2();
   h_subleading_jetphi     = new TH1D(Form("_h_subleading_jetphi_")     ,Form("subleading jet phi (af cut)")      ,24     ,-1*pi    ,1*pi   ); h_subleading_jetphi->Sumw2();

   // before pairing: dR value between leading_bTagJet & Lep1,Lep2 ||  2D: // dR value for (pair_bJet & Lep) & (pair_AnbJet & AnLep)
   h_dR_bTagJet_Lep        = new TH1D("_h_dR_bTagJet_Lep_"        ,"dR value btw leading_bTagJet & Lep"    , 100, 0, 10); h_dR_bTagJet_Lep->Sumw2();
   h_dR_bTagJet_AnLep      = new TH1D("_h_dR_bTagJet_AnLep_"      ,"dR value btw leading_bTagJet & AnLep"  , 100, 0, 10); h_dR_bTagJet_AnLep->Sumw2();
   h_dR_bJet_AnbJet        = new TH2D("_h_dR_bJet_AnbJet_"        ,"delta R value scatter: (pair_bJet&pair_Lep) & (pair_AnbJet&AnLep)", 100, 0, 10, 100, 0, 10); h_dR_bJet_AnbJet->Sumw2();
   // after pairing: dR value between pair_bJet & Lep1,Lep2
   h_dR_pair_bJet_AnLep    = new TH1D("_h_dR_pair_bJet_AnLep_"    ,"delta R btw pair_bJet & Lep"           , 100, 0, 10); h_dR_pair_bJet_AnLep->Sumw2();
   h_dR_pair_AnbJet_Lep    = new TH1D("_h_dR_pair_AnbJet_Lep_"    ,"delta R btw pair_bJet & AnLep"         , 100, 0, 10); h_dR_pair_AnbJet_Lep->Sumw2();
   h_dR_pair_bJet_AnbJet   = new TH2D("_h_dR_pair_bJet_AnbJet_"   ,"delta R value scatter: (pair_bJet&pair_Lep) & (pair_AnbJet&AnLep)", 100, 0, 10, 100, 0, 10); h_dR_pair_bJet_AnbJet->Sumw2();
   // after pairing about each dR values : 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 10.0
   for(int i=0; i<7; i++)
   {
      // dR distribution
      h_dRNum_pair_bJet_AnLep[i]   = new TH1D(Form("_h_dRNum_pair_bJet_AnLep_%d_",i)    ,Form("delta R btw pair_bJet & Lep: dR < %s",cutdRName[i].Data())         , 100, 0, 10); h_dRNum_pair_bJet_AnLep[i]->Sumw2();
      h_dRNum_pair_AnbJet_Lep[i]   = new TH1D(Form("_h_dRNum_pair_AnbJet_Lep_%d_",i)    ,Form("delta R btw pair_bJet & AnLep: dR < %s",cutdRName[i].Data())       , 100, 0, 10); h_dRNum_pair_AnbJet_Lep[i]->Sumw2();
      h_dRNum_pair_bJet_AnbJet[i]  = new TH2D(Form("_h_dRNum_pair_bJet_AnbJet_%d_",i)   ,Form("delta R value scatter: (pair_bJet&pair_Lep) & (pair_AnbJet&AnLep): dR < %s",cutdRName[i].Data()), 100, 0, 10, 100, 0, 10); h_dRNum_pair_bJet_AnbJet[i]->Sumw2();
      // kinematics of bJet, AnbJet, Lep, AnLep 
      h_Pair_bJetpt[i]     = new TH1D(Form("_h_Pair_bJetpt_%d_",i)     ,Form("b tagging jet(bJet) pT after pairing: dR < %s",cutdRName[i].Data())     ,1000   ,0        ,1000   ); h_Pair_bJetpt[i]->Sumw2();
      h_Pair_bJeteta[i]    = new TH1D(Form("_h_Pair_bJeteta_%d_",i)    ,Form("b tagging jet(bJet) eta after pairing: dR < %s",cutdRName[i].Data())    ,50     ,-2.5     ,2.5    ); h_Pair_bJeteta[i]->Sumw2();
      h_Pair_bJetphi[i]    = new TH1D(Form("_h_Pair_bJetphi_%d_",i)    ,Form("b tagging jet(bJet) phi after pairing: dR < %s",cutdRName[i].Data())    ,24     ,-1*pi    ,1*pi   ); h_Pair_bJetphi[i]->Sumw2();
      h_Pair_AnbJetpt[i]   = new TH1D(Form("_h_Pair_AnbJetpt_%d_",i)   ,Form("b tagging jet pT (AnbJet after pairing): dR < %s",cutdRName[i].Data())  ,1000   ,0        ,1000   ); h_Pair_AnbJetpt[i]->Sumw2();
      h_Pair_AnbJeteta[i]  = new TH1D(Form("_h_Pair_AnbJeteta_%d_",i)  ,Form("b tagging jet eta (AnbJet after pairing): dR < %s",cutdRName[i].Data()) ,50     ,-2.5     ,2.5    ); h_Pair_AnbJeteta[i]->Sumw2();
      h_Pair_AnbJetphi[i]  = new TH1D(Form("_h_Pair_AnbJetphi_%d_",i)  ,Form("b tagging jet phi (AnbJet after pairing): dR < %s",cutdRName[i].Data()) ,24     ,-1*pi    ,1*pi   ); h_Pair_AnbJetphi[i]->Sumw2();
      h_pair_Leppt[i]      = new TH1D(Form("_h_pair_Leppt_%d_",i)      ,Form("Lep pT after pairing: dR < %s",cutdRName[i].Data())                     ,1000   ,0        ,1000   ); h_pair_Leppt[i]->Sumw2();
      h_pair_Lepeta[i]     = new TH1D(Form("_h_pair_Lepeta_%d_",i)     ,Form("Lep eta after pairing: dR < %s",cutdRName[i].Data())                    ,50     ,-2.5     ,2.5    ); h_pair_Lepeta[i]->Sumw2();
      h_pair_Lepphi[i]     = new TH1D(Form("_h_pair_Lepphi_%d_",i)     ,Form("Lep phi after pairing: dR < %s",cutdRName[i].Data())                    ,24     ,-1*pi    ,1*pi   ); h_pair_Lepphi[i]->Sumw2();
      h_pair_AnLeppt[i]    = new TH1D(Form("_h_pair_AnLeppt_%d_",i)    ,Form("AnLep pT after pairing: dR < %s",cutdRName[i].Data())                   ,1000   ,0        ,1000   ); h_pair_AnLeppt[i]->Sumw2();
      h_pair_AnLepeta[i]   = new TH1D(Form("_h_pair_AnLepeta_%d_",i)   ,Form("AnLep eta after pairing: dR < %s",cutdRName[i].Data())                  ,50     ,-2.5     ,2.5    ); h_pair_AnLepeta[i]->Sumw2();
      h_pair_AnLepphi[i]   = new TH1D(Form("_h_pair_AnLepphi_%d_",i)   ,Form("AnLep phi after pairing: dR < %s",cutdRName[i].Data())                  ,24     ,-1*pi    ,1*pi   ); h_pair_AnLepphi[i]->Sumw2();
      // invariant mass
      h_pair_b_Anb_Mass[i]    = new TH1D(Form("_h_pair_b_Anb_Mass_%d_",i)    ,Form("b & Anti-b Invariant Mass after pairing: dR < %s",cutdRName[i].Data())   ,1000   ,0.0      ,1000   ); h_pair_b_Anb_Mass[i]->Sumw2();
      h_pair_Mu_AnMu_Mass[i]  = new TH1D(Form("_h_pair_Mu_AnMu_Mass_%d_",i)  ,Form("Mu & Anti-Mu Invariant Mass after pairing: dR < %s",cutdRName[i].Data()) ,1000   ,0.0      ,1000   ); h_pair_Mu_AnMu_Mass[i]->Sumw2();
      h_pair_b_AnMu_Mass[i]   = new TH1D(Form("_h_pair_b_AnMu_Mass_%d_",i)   ,Form("b & Anti-Mu Invariant Mass after pairing: dR < %s",cutdRName[i].Data())  ,1000   ,0.0      ,1000   ); h_pair_b_AnMu_Mass[i]->Sumw2();
      h_pair_Anb_Mu_Mass[i]   = new TH1D(Form("_h_pair_Anb_Mu_Mass_%d_",i)   ,Form("Anti-b & Mu Invariant Mass after pairing: dR < %s",cutdRName[i].Data())  ,1000   ,0.0      ,1000   ); h_pair_Anb_Mu_Mass[i]->Sumw2();
      // O3 before top reconstruction
      h_CPO3_bfReco[i]           = new TH1D(Form("_h_CPO3_bfReco_%d_",i)           , Form(""), 200   , -10    , 10);   h_CPO3_bfReco[i]->Sumw2();
      h_CPO3_bfReco_ReRange[i]   = new TH1D(Form("_h_CPO3_bfReco_ReRange_%d_",i)   , Form(""), 40    , -2     , 2);    h_CPO3_bfReco_ReRange[i]->Sumw2();
      h_CPO3_Pos[i]              = new TH1D(Form("_h_CPO3_Pos_%d_",i)              , Form(""), 200   , -10    , 10);   h_CPO3_Pos[i]->Sumw2();
      h_CPO3_Neg[i]              = new TH1D(Form("_h_CPO3_Neg_%d_",i)              , Form(""), 200   , -10    , 10);   h_CPO3_Neg[i]->Sumw2();
   }
   ///// generator level /////
   h_GenMuon_pt      = new TH1D(Form("_h_GenMuon_pt_")      ,Form("Gen muon pT")          ,1000   ,0        ,1000   ); h_GenMuon_pt->Sumw2();
   h_GenMuon_eta     = new TH1D(Form("_h_GenMuon_eta_")     ,Form("Gen muon eta")         ,50     ,-2.5     ,2.5    ); h_GenMuon_eta->Sumw2();
   h_GenMuon_phi     = new TH1D(Form("_h_GenMuon_phi_")     ,Form("Gen muon phi")         ,24     ,-1*pi    ,1*pi   ); h_GenMuon_phi->Sumw2();
   h_GenAnMuon_pt    = new TH1D(Form("_h_GenAnMuon_pt_")    ,Form("Gen anti-muon pT")     ,1000   ,0        ,1000   ); h_GenAnMuon_pt->Sumw2();
   h_GenAnMuon_eta   = new TH1D(Form("_h_GenAnMuon_eta_")   ,Form("Gen anti-muon eta")    ,50     ,-2.5     ,2.5    ); h_GenAnMuon_eta->Sumw2();
   h_GenAnMuon_phi   = new TH1D(Form("_h_GenAnMuon_phi_")   ,Form("Gen anti-muon phi")    ,24     ,-1*pi    ,1*pi   ); h_GenAnMuon_phi->Sumw2();
   h_Genb_pt         = new TH1D(Form("_h_Genb_pt_")         ,Form("Gen b quark pT")       ,1000   ,0        ,1000   ); h_Genb_pt->Sumw2();
   h_Genb_eta        = new TH1D(Form("_h_Genb_eta_")        ,Form("Gen b quark eta")      ,50     ,-2.5     ,2.5    ); h_Genb_eta->Sumw2();
   h_Genb_phi        = new TH1D(Form("_h_Genb_phi_")        ,Form("Gen b quark phi")      ,24     ,-1*pi    ,1*pi   ); h_Genb_phi->Sumw2();
   h_GenAnb_pt       = new TH1D(Form("_h_GenAnb_pt_")       ,Form("Gen anti-b quark pT")  ,1000   ,0        ,1000   ); h_GenAnb_pt->Sumw2();
   h_GenAnb_eta      = new TH1D(Form("_h_GenAnb_eta_")      ,Form("Gen anti-b quark eta") ,50     ,-2.5     ,2.5    ); h_GenAnb_eta->Sumw2();
   h_GenAnb_phi      = new TH1D(Form("_h_GenAnb_phi_")      ,Form("Gen anti-b quark phi") ,24     ,-1*pi    ,1*pi   ); h_GenAnb_phi->Sumw2();
   // dR distribution
   h_Gen_dR_Genb_GenAnLep   = new TH1D("_h_Gen_dR_Genb_GenAnLep_"   , "delta R btw Gen b & anti-muon"         , 100, 0, 10); h_Gen_dR_Genb_GenAnLep->Sumw2();
   h_Gen_dR_GenAnb_GenLep   = new TH1D("_h_Gen_dR_GenAnb_GenLep_"   , "delta R btw Gen anti-b & muon"         , 100, 0, 10); h_Gen_dR_GenAnb_GenLep->Sumw2();
   h_Gen_dR_Genb_GenAnb  = new TH2D("_h_Gen_dR_Genb_GenAnb_"  , "delta R value scatter: (b & anti-muon) & (anti-b & muon)", 100, 0, 10, 100, 0, 10); h_Gen_dR_Genb_GenAnb->Sumw2();
   // dR is under 3 case
   h_Gen_dR3_pair_Genb_GenAnLep   = new TH1D("_h_Gen_dR3_pair_Genb_GenAnLep_"   , "delta R btw Gen b & anti-muon"                   , 100, 0, 10); h_Gen_dR3_pair_Genb_GenAnLep->Sumw2();
   h_Gen_dR3_pair_GenAnb_GenLep   = new TH1D("_h_Gen_dR3_pair_GenAnb_GenLep_"   , "delta R btw Gen anti-b & muon"                   , 100, 0, 10); h_Gen_dR3_pair_GenAnb_GenLep->Sumw2();
   h_Gen_dR3_af_Genb_GenAnb  = new TH2D("_h_Gen_dR3_af_Genb_GenAnb_"  , "delta R value scatter: (b & anti-muon) & (anti-b & muon)"  , 100, 0, 10, 100, 0, 10); h_Gen_dR3_af_Genb_GenAnb->Sumw2();
   // when dR is under 3
   h_Pair_Genbpt       = new TH1D(Form("_h_Pair_Genbpt_")       ,Form("Gen b quark pT after pairing")           ,1000   ,0        ,1000   ); h_Pair_Genbpt->Sumw2();
   h_Pair_Genbeta      = new TH1D(Form("_h_Pair_Genbeta_")      ,Form("Gen b quark eta after pairing")          ,50     ,-2.5     ,2.5    ); h_Pair_Genbeta->Sumw2();
   h_Pair_Genbphi      = new TH1D(Form("_h_Pair_Genbphi_")      ,Form("Gen b quark phi after pairing")          ,24     ,-1*pi    ,1*pi   ); h_Pair_Genbphi->Sumw2();
   h_Pair_GenAnbpt     = new TH1D(Form("_h_Pair_GenAnbpt_")     ,Form("Gen anti-b quark pT after pairing")      ,1000   ,0        ,1000   ); h_Pair_GenAnbpt->Sumw2();
   h_Pair_GenAnbeta    = new TH1D(Form("_h_Pair_GenAnbeta_")    ,Form("Gen anti-b quark eta after pairing")     ,50     ,-2.5     ,2.5    ); h_Pair_GenAnbeta->Sumw2();
   h_Pair_GenAnbphi    = new TH1D(Form("_h_Pair_GenAnbphi_")    ,Form("Gen anti-b quark phi after pairing")     ,24     ,-1*pi    ,1*pi   ); h_Pair_GenAnbphi->Sumw2();
   h_pair_GenMupt      = new TH1D(Form("_h_pair_GenMupt_")      ,Form("Gen muon pT after pairing")              ,1000   ,0        ,1000   ); h_pair_GenMupt->Sumw2();
   h_pair_GenMueta     = new TH1D(Form("_h_pair_GenMueta_")     ,Form("Gen muon eta after pairing")             ,50     ,-2.5     ,2.5    ); h_pair_GenMueta->Sumw2();
   h_pair_GenMuphi     = new TH1D(Form("_h_pair_GenMuphi_")     ,Form("Gen muon phi after pairing")             ,24     ,-1*pi    ,1*pi   ); h_pair_GenMuphi->Sumw2();
   h_pair_GenAnMupt    = new TH1D(Form("_h_pair_GenAnMupt_")    ,Form("Gen anti-muon quark pT after pairing")   ,1000   ,0        ,1000   ); h_pair_GenAnMupt->Sumw2();
   h_pair_GenAnMueta   = new TH1D(Form("_h_pair_GenAnMueta_")   ,Form("Gen anti-muon quark eta after pairing")  ,50     ,-2.5     ,2.5    ); h_pair_GenAnMueta->Sumw2();
   h_pair_GenAnMuphi   = new TH1D(Form("_h_pair_GenAnMuphi_")   ,Form("Gen anti-muon quark phi after pairing")  ,24     ,-1*pi    ,1*pi   ); h_pair_GenAnMuphi->Sumw2();
   // invariant mass
   h_pair_Gen_b_Anb_Mass    = new TH1D(Form("_h_pair_Gen_b_Anb_Mass_")    ,Form("Gen: b & Anti-b Invariant Mass after pairing")   , 1000, 0.0, 1000); h_pair_Gen_b_Anb_Mass->Sumw2();
   h_pair_Gen_Mu_AnMu_Mass  = new TH1D(Form("_h_pair_Gen_Mu_AnMu_Mass_")  ,Form("Gen: Mu & Anti-Mu Invariant Mass after pairing") , 1000, 0.0, 1000); h_pair_Gen_Mu_AnMu_Mass->Sumw2();
   h_pair_Gen_b_AnMu_Mass   = new TH1D(Form("_h_pair_Gen_b_AnMu_Mass_")   ,Form("Gen: b & Anti-Mu Invariant Mass after pairing")  , 1000, 0.0, 1000); h_pair_Gen_b_AnMu_Mass->Sumw2();
   h_pair_Gen_Anb_Mu_Mass   = new TH1D(Form("_h_pair_Gen_Anb_Mu_Mass_")   ,Form("Gen: Anti-b & Mu Invariant Mass after pairing")  , 1000, 0.0, 1000); h_pair_Gen_Anb_Mu_Mass->Sumw2();
   h_pair_Gen_b_Mass        = new TH1D(Form("_h_pair_Gen_b_Mass_")        ,Form("Gen: b Invariant Mass after pairing")            , 1000, 0.0, 1000); h_pair_Gen_b_Mass->Sumw2();
   h_pair_Gen_Anb_Mass      = new TH1D(Form("_h_pair_Gen_Anb_Mass_")      ,Form("Gen: Anti-b Invariant Mass after pairing")       , 1000, 0.0, 1000); h_pair_Gen_Anb_Mass->Sumw2();
   h_pair_Gen_Mu_Mass       = new TH1D(Form("_h_pair_Gen_Mu_Mass_")       ,Form("Gen: Mu Invariant Mass after pairing")           , 1000, 0.0, 1000); h_pair_Gen_Mu_Mass->Sumw2();
   h_pair_Gen_AnMu_Mass     = new TH1D(Form("_h_pair_Gen_AnMu_Mass_")     ,Form("Gen: Anti-Mu Invariant Mass after pairing")      , 1000, 0.0, 1000); h_pair_Gen_AnMu_Mass->Sumw2();
   // O3 before top reconstruction
   h_GenCPO3_bfReco           = new TH1D(Form("_h_GenCPO3_bfReco_")           , Form(""), 200   , -10    , 10);   h_GenCPO3_bfReco->Sumw2();
   h_GenCPO3_bfReco_ReRange   = new TH1D(Form("_h_GenCPO3_bfReco_ReRange_")   , Form(""), 40    , -2     , 2);    h_GenCPO3_bfReco_ReRange->Sumw2();
   h_GenCPO3_Pos              = new TH1D(Form("_h_GenCPO3_Pos_")              , Form(""), 200   , -10    , 10);   h_GenCPO3_Pos->Sumw2();
   h_GenCPO3_Neg              = new TH1D(Form("_h_GenCPO3_Neg_")              , Form(""), 200   , -10    , 10);   h_GenCPO3_Neg->Sumw2();

   ///// O3 version 2 /////
   ////////////////////////
   // before pairing     //
   ////////////////////////
   h_AllJet_idx            = new TH1D(Form("_h_AllJet_idx_")     ,Form("Dist for idx of Jet2 after pairing")       ,10      ,0        ,10);   h_AllJet_idx->Sumw2();
   h_AllJet_pt             = new TH1D(Form("_h_AllJet_pt_")      ,Form("Jet1 pT when dR is minimum after pairing") ,1000    ,0        ,1000); h_AllJet_pt->Sumw2();
   h_AllJet_eta            = new TH1D(Form("_h_AllJet_eta_")     ,Form("Jet1 eta when dR is minimum")              ,50      ,-2.5     ,2.5 ); h_AllJet_eta->Sumw2();
   h_AllJet_phi            = new TH1D(Form("_h_AllJet_phi_")     ,Form("Jet1 phi when dR is minimum")              ,24      ,-1*pi    ,1*pi); h_AllJet_phi->Sumw2();
   ////////////////////////////////
   // calculate dR with all jets //
   ////////////////////////////////
   // all of dR distribution
   h_v2_dR_AnLep           = new TH1D(Form("_h_v2_dR_AnLep_")    ,Form("dR value btw all Jets & AnLep")                          ,100  ,0 ,10);                h_v2_dR_AnLep->Sumw2();
   h_v2_dR_Lep             = new TH1D(Form("_h_v2_dR_Lep_")      ,Form("dR value btw all Jets & Lep")                            ,100  ,0 ,10);                h_v2_dR_Lep->Sumw2();
   h_v2_dR_Lep_AnLep       = new TH2D(Form("_h_v2_dR_Lep_AnLep_"),Form("delta R value scatter: (Jet1 & AnLep) & (Jet2 & Lep)")   ,100  ,0 ,10   ,100  ,0 ,10); h_v2_dR_Lep_AnLep->Sumw2();
   // find minimum dR
   h_v2_dR_Min_Lep         = new TH1D(Form("_h_v2_dR_Min_AnJet_Lep_")   ,Form("min delta R btw AnJet & Lep")                              ,100  ,0 ,10);                h_v2_dR_Min_Lep->Sumw2();
   h_v2_dR_Min_AnLep       = new TH1D(Form("_h_v2_dR_Min_Jet_AnLep_")   ,Form("min delta R btw Jet & AnLep")                              ,100  ,0 ,10);                h_v2_dR_Min_AnLep->Sumw2();
   h_v2_dR_Min_Lep_AnLep   = new TH2D(Form("_h_v2_dR_Min_Jet_AnJet_")   ,Form("min delta R value scatter: (Jet & AnLep) & (AnJet & Lep)") ,100  ,0 ,10   ,100  ,0 ,10); h_v2_dR_Min_Lep_AnLep->Sumw2();
   ///////////////////////
   // after pairing     //
   ///////////////////////
   // Jet & AnLep //
   h_v2_pair_Jet_dR_AnLep       = new TH1D(Form("_h_v2_pair_dR_Jet_AnLep_")      ,Form("dR value btw Jet & AnLep")                  ,100     ,0        ,10);   h_v2_pair_Jet_dR_AnLep->Sumw2();
   h_v2_pair_Jet_idx            = new TH1D(Form("_h_v2_pair_Jet_idx_")           ,Form("Dist for idx of Jet after pairing")         ,10      ,0        ,10);   h_v2_pair_Jet_idx->Sumw2();
   h_v2_pair_Jetpt              = new TH1D(Form("_h_v2_pair_Jetpt_")             ,Form("Jet pT when dR is minimum after pairing")   ,1000    ,0        ,1000); h_v2_pair_Jetpt->Sumw2();
   h_v2_pair_Jeteta             = new TH1D(Form("_h_v2_pair_Jeteta_")            ,Form("Jet eta when dR is minimum")                ,50      ,-2.5     ,2.5 ); h_v2_pair_Jeteta->Sumw2();
   h_v2_pair_Jetphi             = new TH1D(Form("_h_v2_pair_Jetphi_")            ,Form("Jet phi when dR is minimum")                ,24      ,-1*pi    ,1*pi); h_v2_pair_Jetphi->Sumw2();
   h_v2_pair_AnLeppt            = new TH1D(Form("_h_v2_pair_AnLeppt_")           ,Form("AnLep pT when dR is minimum")               ,1000    ,0        ,1000); h_v2_pair_AnLeppt->Sumw2();
   h_v2_pair_AnLepeta           = new TH1D(Form("_h_v2_pair_AnLepeta_")          ,Form("AnLep eta when dR is minimum")              ,50      ,-2.5     ,2.5 ); h_v2_pair_AnLepeta->Sumw2();
   h_v2_pair_AnLepphi           = new TH1D(Form("_h_v2_pair_AnLepphi_")          ,Form("AnLep phi when dR is minimum")              ,24      ,-1*pi    ,1*pi); h_v2_pair_AnLepphi->Sumw2();
   // AnJet & Lep //
   h_v2_pair_AnJet_dR_Lep       = new TH1D(Form("_h_v2_pair_dR_AnJet_Lep_")      ,Form("dR value btw AnJet & Lep")                  ,100     ,0        ,10);   h_v2_pair_AnJet_dR_Lep->Sumw2();
   h_v2_pair_AnJet_idx          = new TH1D(Form("_h_v2_pair_AnJet_idx_")         ,Form("Dist for idx of AnJet after pairing")       ,10      ,0        ,10);   h_v2_pair_AnJet_idx->Sumw2();
   h_v2_pair_AnJetpt            = new TH1D(Form("_h_v2_pair_AnJetpt_")           ,Form("AnJet pT when dR is minimum")               ,1000    ,0        ,1000); h_v2_pair_AnJetpt->Sumw2();
   h_v2_pair_AnJeteta           = new TH1D(Form("_h_v2_pair_AnJeteta_")          ,Form("AnJet eta when dR is minimum")              ,50      ,-2.5     ,2.5 ); h_v2_pair_AnJeteta->Sumw2();
   h_v2_pair_AnJetphi           = new TH1D(Form("_h_v2_pair_AnJetphi_")          ,Form("AnJet phi when dR is minimum")              ,24      ,-1*pi    ,1*pi); h_v2_pair_AnJetphi->Sumw2();
   h_v2_pair_Leppt              = new TH1D(Form("_h_v2_pair_Leppt_")             ,Form("Lep pT when dR is minimum")                 ,1000    ,0        ,1000); h_v2_pair_Leppt->Sumw2();
   h_v2_pair_Lepeta             = new TH1D(Form("_h_v2_pair_Lepeta_")            ,Form("Lep eta when dR is minimum")                ,50      ,-2.5     ,2.5 ); h_v2_pair_Lepeta->Sumw2();
   h_v2_pair_Lepphi             = new TH1D(Form("_h_v2_pair_Lepphi_")            ,Form("Lep phi when dR is minimum")                ,24      ,-1*pi    ,1*pi); h_v2_pair_Lepphi->Sumw2();
   // distribution of idx
   h_v2_pair_Jet_AnJet_idx          = new TH2D(Form("_h_v2_pair_Jet_AnJet_idx_")          ,Form("idx scatter: Jet & AnJet after pairing")       ,10      ,0        ,10     ,10   ,0 ,10);    h_v2_pair_Jet_AnJet_idx->Sumw2();
   h_v2_pair_comp_Jet_pairbJet      = new TH2D(Form("_h_v2_pair_comp_Jet_pairbJet_")      ,Form("comparison for v1 & v2 Jet after pairing")     ,1000    ,0        ,1000   ,1000 ,0 ,1000);  h_v2_pair_comp_Jet_pairbJet->Sumw2();
   h_v2_pair_comp_AnJet_pairAnbJet  = new TH2D(Form("_h_v2_pair_comp_AnJet_pairAnbJet_")  ,Form("comparison for v1 & v2 AnJet after pairing")   ,1000    ,0        ,1000   ,1000 ,0 ,1000);  h_v2_pair_comp_AnJet_pairAnbJet->Sumw2();
   // invariant mass
   h_v2_pair_Jet_AnMu_Mass       = new TH1D(Form("_h_v2_pair_Jet_AnMu_Mass_")    ,Form("Jet & AnLep Invariant Mass when minimum dR")     ,1000 ,0.0  ,1000);                                          h_v2_pair_Jet_AnMu_Mass->Sumw2();
   h_v2_pair_AnJet_Mu_Mass       = new TH1D(Form("_h_v2_pair_AnJet_Mu_Mass_")    ,Form("AnJet & Lep Invariant Mass when minimum dR")       ,1000 ,0.0  ,1000);                                        h_v2_pair_AnJet_Mu_Mass->Sumw2();
   h_v2_pair_Jet_AnbJet_Mass     = new TH2D(Form("_h_v2_pair_Jet_AnbJet_Mass_")  ,Form("Jet+AnJet & pair_bJet+pair_AnbJet Invariant Mass when minimum dR")      , 1000, 0.0, 1000, 1000, 0.0, 1000);  h_v2_pair_Jet_AnbJet_Mass->Sumw2();
   // O3
   h_v2_CPO3_bfReco              = new TH1D(Form("_h_v2_CPO3_bfReco_")              , Form(""), 200   , -10    , 10);   h_v2_CPO3_bfReco->Sumw2();
   h_v2_CPO3_bfReco_ReRange      = new TH1D(Form("_h_v2_CPO3_bfReco_ReRange_"), Form(""), 40    , -2     , 2);    h_v2_CPO3_bfReco_ReRange->Sumw2();
   h_v2_CPO3_Pos                 = new TH1D(Form("_h_v2_CPO3_Pos_")           , Form(""), 200   , -10    , 10);   h_v2_CPO3_Pos->Sumw2();
   h_v2_CPO3_Neg                 = new TH1D(Form("_h_v2_CPO3_Neg_")           , Form(""), 200   , -10    , 10);   h_v2_CPO3_Neg->Sumw2();
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   /////////////////////////////
   /// To Check Evetn Weight ///
   /////////////////////////////

   h_LepbJetMass          = new TH1D(Form("h_LepbJetMass"), Form("Lep + bJet Invariant Mass"), 1000, 0, 1000); h_LepbJetMass->Sumw2();
   h_AnLepbJetMass        = new TH1D(Form("h_AnLepbJetMass"), Form("AnLep + bJet Invariant Mass"), 1000, 0, 1000); h_AnLepbJetMass->Sumw2();
   h2_TopMassVsLepBMass    = new TH2F(Form("TopMassVsLepBMass"), Form("Top Mass Vs (Lep+B) Mass"), 1000, 0, 1000, 1000,0,1000); h2_TopMassVsLepBMass->Sumw2();
   h2_AnTopMassVsLepBMass  = new TH2F(Form("AnTopMassVsLepBMass"), Form("Anti-Top Mass Vs (Lep+B) Mass"), 1000, 0, 1000, 1000,0,1000); h2_AnTopMassVsLepBMass->Sumw2();
   h2_AnLepBMassVsLepBMass     = new TH2F(Form("AnLepBMassVsLepBMass"), Form("(AntiLep+B) Mass Vs (Lep+B) Mass"), 1000, 0, 1000, 1000,0,1000); h2_AnLepBMassVsLepBMass->Sumw2();

   h_Diff_BBar_Pt     = new TH1D(Form("h_Diff_BBar_Pt"    ), Form("PT Difference between B and B bar "     ), 500, 0, 500); h_Diff_BBar_Pt->Sumw2();
   h_Diff_BBar_Energy = new TH1D(Form("h_Diff_BBar_Energy"), Form("Energy Difference between B and B bar"  ), 500, 0, 500); h_Diff_BBar_Energy->Sumw2();
   h_Diff_BBar_P      = new TH1D(Form("h_Diff_BBar_P"     ), Form("Momentum Difference between B and B bar"), 500, 0, 500); h_Diff_BBar_P->Sumw2();

   h_Diff_LepAnLep_Pt     = new TH1D(Form("h_Diff_LepAnLep_Pt"    ), Form("PT Difference between Lep. And AntiLep."      ), 500, 0, 500); h_Diff_LepAnLep_Pt->Sumw2();
   h_Diff_LepAnLep_Energy = new TH1D(Form("h_Diff_LepAnLep_Energy"), Form("Energy Difference between Lep. And AntiLep."  ), 500, 0, 500); h_Diff_LepAnLep_Energy->Sumw2();
   h_Diff_LepAnLep_P      = new TH1D(Form("h_Diff_LepAnLep_P"     ), Form("Momentum Difference between Lep. And AntiLep."), 500, 0, 500); h_Diff_LepAnLep_P->Sumw2();

   h_LepAnLepEngCheckO3 = new TH1D(Form("h_LepAnLepEngCheckO3"), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for O3"), 100, 0, 100); h_LepAnLepEngCheckO3->Sumw2();
   h_LepAnLepEngCheckOb = new TH1D(Form("h_LepAnLepEngCheckOb"), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for Ob"), 100, 0, 100); h_LepAnLepEngCheckOb->Sumw2();
   h_LepAnLepEngCheckO5 = new TH1D(Form("h_LepAnLepEngCheckO5"), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for O5"), 100, 0, 100); h_LepAnLepEngCheckO5->Sumw2();

   h_PileUpCheck  = new TH1D(Form("h_PileUpCheck"),  Form("NewPileUpMethod - OriginMethod"), 4000, -2, 2); h_PileUpCheck->Sumw2();
   //////////////////////////////////
   /// For CP-Violation Variables ///
   //////////////////////////////////
   for (int i =0; i < 13; ++i)
   {
      h_Reco_CPO_[i] = new TH1D(Form("h_Reco_CPO%d",i+1 ), Form("CPO%d",i+1   ), 200, -10, 10); h_Reco_CPO_[i]->Sumw2();
      h_Reco_CPO_ReRange_[i] = new TH1D(Form("h_Reco_CPO%d_ReRange",i+1 ), Form("CPO%d",i+1   ), 40, -2, 2); h_Reco_CPO_ReRange_[i]->Sumw2();
      h_LepAnLepEngCheck_[i] = new TH1D(Form("h_LepAnLepEngCheck_%d",i+1), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for CPO%d",i+1), 100, 0, 100); h_LepAnLepEngCheck_[i]->Sumw2();
   }
   ////////////////////////////////
   /// To Check O3 dilution ... ///
   ////////////////////////////////
   h_CPO3_Plus              = new TH1D(Form("h_CPO3_Plus"),              Form("CPO3 Positive"), 200, -10, 10); h_CPO3_Plus->Sumw2(); 
   h_CPO3_Minus             = new TH1D(Form("h_CPO3_Minus"),             Form("CPO3 Negative"), 200, -10, 10); h_CPO3_Minus->Sumw2(); 
   h_CPO3_JPRUp_Plus_Plus   = new TH1D(Form("h_CPO3_JPRUp_Plus_Plus"),   Form("CPO3_JPRUp Plus -> Plus"),   200, -10, 10); h_CPO3_JPRUp_Plus_Plus->Sumw2();                   
   h_CPO3_JPRUp_Plus_Minus  = new TH1D(Form("h_CPO3_JPRUp_Plus_Minus"),  Form("CPO3_JPRUp Plus -> Minus"),  200, -10, 10); h_CPO3_JPRUp_Plus_Minus->Sumw2();                   
   h_CPO3_JPRUp_Minus_Minus = new TH1D(Form("h_CPO3_JPRUp_Minus_Minus"), Form("CPO3_JPRUp Minus -> Minus"), 200, -10, 10); h_CPO3_JPRUp_Minus_Minus->Sumw2();                   
   h_CPO3_JPRUp_Minus_Plus  = new TH1D(Form("h_CPO3_JPRUp_Minus_Plus"),  Form("CPO3_JPRUp Minus -> Minus"), 200, -10, 10); h_CPO3_JPRUp_Minus_Plus->Sumw2();                   

   h_CPO3_JPRDown_Plus_Plus   = new TH1D(Form("h_CPO3_JPRDown_Plus_Plus"),   Form("CPO3_JPRDown Plus-> Plus"),   200, -10, 10); h_CPO3_JPRDown_Plus_Plus->Sumw2();
   h_CPO3_JPRDown_Plus_Minus  = new TH1D(Form("h_CPO3_JPRDown_Plus_Minus"),  Form("CPO3_JPRDown Plus-> Minus"),  200, -10, 10); h_CPO3_JPRDown_Plus_Minus->Sumw2();
   h_CPO3_JPRDown_Minus_Minus = new TH1D(Form("h_CPO3_JPRDown_Minus_Minus"), Form("CPO3_JPRDown Minus-> Minus"), 200, -10, 10); h_CPO3_JPRDown_Minus_Minus->Sumw2();
   h_CPO3_JPRDown_Minus_Plus  = new TH1D(Form("h_CPO3_JPRDown_Minus_Plus"),  Form("CPO3_JPRDown Minus-> Minus"), 200, -10, 10); h_CPO3_JPRDown_Minus_Plus->Sumw2();


   h_CPO3reco_evtweight             = new TH1D(Form("h_CPO3reco_evtweight"),             Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight->Sumw2();
   h_CPO3reco_evtweight_EtaVariUp   = new TH1D(Form("h_CPO3reco_evtweight_EtaVariUp"),   Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_EtaVariUp->Sumw2();
   h_CPO3reco_evtweight_EtaVariDown = new TH1D(Form("h_CPO3reco_evtweight_EtaVariDown"), Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_EtaVariDown->Sumw2();
   h_CPO3reco_evtweight_PhiVariUp   = new TH1D(Form("h_CPO3reco_evtweight_PhiVariUp"),   Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_PhiVariUp->Sumw2();
   h_CPO3reco_evtweight_PhiVariDown = new TH1D(Form("h_CPO3reco_evtweight_PhiVariDown"), Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_PhiVariDown->Sumw2();

   h_CPO3_EtaVariUp   = new TH1D(Form("h_CPO3_EtaVariUp"),   Form("CPO3_EtaVariUp"),   200, -10, 10); h_CPO3_EtaVariUp->Sumw2(); 
   h_CPO3_EtaVariUp_Plus_Plus   = new TH1D(Form("h_CPO3_EtaVariUp_Plus_Plus"),   Form("CPO3_EtaVariUp Plus-> Plus"),   200, -10, 10); h_CPO3_EtaVariUp_Plus_Plus->Sumw2();
   h_CPO3_EtaVariUp_Plus_Minus  = new TH1D(Form("h_CPO3_EtaVariUp_Plus_Minus"),  Form("CPO3_EtaVariUp Plus-> Minus"),  200, -10, 10); h_CPO3_EtaVariUp_Plus_Minus->Sumw2();
   h_CPO3_EtaVariUp_Minus_Minus = new TH1D(Form("h_CPO3_EtaVariUp_Minus_Minus"), Form("CPO3_EtaVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariUp_Minus_Minus->Sumw2();
   h_CPO3_EtaVariUp_Minus_Plus  = new TH1D(Form("h_CPO3_EtaVariUp_Minus_Plus"),  Form("CPO3_EtaVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariUp_Minus_Plus->Sumw2();

   h_CPO3_EtaVariDown   = new TH1D(Form("h_CPO3_EtaVariDown"),   Form("CPO3_EtaVariDown"),   200, -10, 10); h_CPO3_EtaVariDown->Sumw2(); 
   h_CPO3_EtaVariDown_Plus_Plus   = new TH1D(Form("h_CPO3_EtaVariDown_Plus_Plus"),   Form("CPO3_EtaVariDown Plus-> Plus"),   200, -10, 10); h_CPO3_EtaVariDown_Plus_Plus->Sumw2();
   h_CPO3_EtaVariDown_Plus_Minus  = new TH1D(Form("h_CPO3_EtaVariDown_Plus_Minus"),  Form("CPO3_EtaVariDown Plus-> Minus"),  200, -10, 10); h_CPO3_EtaVariDown_Plus_Minus->Sumw2();
   h_CPO3_EtaVariDown_Minus_Minus = new TH1D(Form("h_CPO3_EtaVariDown_Minus_Minus"), Form("CPO3_EtaVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariDown_Minus_Minus->Sumw2();
   h_CPO3_EtaVariDown_Minus_Plus  = new TH1D(Form("h_CPO3_EtaVariDown_Minus_Plus"),  Form("CPO3_EtaVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariDown_Minus_Plus->Sumw2();

   h_CPO3_PhiVariUp   = new TH1D(Form("h_CPO3_PhiVariUp"),   Form("CPO3_PhiVariUp"),   200, -10, 10); h_CPO3_PhiVariUp->Sumw2(); 
   h_CPO3_PhiVariUp_Plus_Plus   = new TH1D(Form("h_CPO3_PhiVariUp_Plus_Plus"),   Form("CPO3_PhiVariUp Plus-> Plus"),   200, -10, 10); h_CPO3_PhiVariUp_Plus_Plus->Sumw2();
   h_CPO3_PhiVariUp_Plus_Minus  = new TH1D(Form("h_CPO3_PhiVariUp_Plus_Minus"),  Form("CPO3_PhiVariUp Plus-> Minus"),  200, -10, 10); h_CPO3_PhiVariUp_Plus_Minus->Sumw2();
   h_CPO3_PhiVariUp_Minus_Minus = new TH1D(Form("h_CPO3_PhiVariUp_Minus_Minus"), Form("CPO3_PhiVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariUp_Minus_Minus->Sumw2();
   h_CPO3_PhiVariUp_Minus_Plus  = new TH1D(Form("h_CPO3_PhiVariUp_Minus_Plus"),  Form("CPO3_PhiVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariUp_Minus_Plus->Sumw2();

   h_CPO3_PhiVariDown   = new TH1D(Form("h_CPO3_PhiVariDown"),   Form("CPO3_PhiVariDown"),   200, -10, 10); h_CPO3_PhiVariDown->Sumw2();
   h_CPO3_PhiVariDown_Plus_Plus   = new TH1D(Form("h_CPO3_PhiVariDown_Plus_Plus"),   Form("CPO3_PhiVariDown Plus-> Plus"),   200, -10, 10); h_CPO3_PhiVariDown_Plus_Plus->Sumw2();
   h_CPO3_PhiVariDown_Plus_Minus  = new TH1D(Form("h_CPO3_PhiVariDown_Plus_Minus"),  Form("CPO3_PhiVariDown Plus-> Minus"),  200, -10, 10); h_CPO3_PhiVariDown_Plus_Minus->Sumw2();
   h_CPO3_PhiVariDown_Minus_Minus = new TH1D(Form("h_CPO3_PhiVariDown_Minus_Minus"), Form("CPO3_PhiVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariDown_Minus_Minus->Sumw2();
   h_CPO3_PhiVariDown_Minus_Plus  = new TH1D(Form("h_CPO3_PhiVariDown_Minus_Plus"),  Form("CPO3_PhiVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariDown_Minus_Plus->Sumw2();

   h_bTagWeight   = new TH1D(Form("h_bTagWeight"),   Form("bTag Weight "),   200, -2, 2); h_bTagWeight->Sumw2();
}

void ssb_analysis::DeclareHistosSyst(int index_)
{
   for (int j =0; j <10; j++) // Cutflow Loop //
   {
      // Cut Flow histogram //
      h_cf_sys_NLeptons[index_][j]       = new TH1D(Form("_h_cf_sys_NLeptons_%s_%d_"      ,v_SystFullName[index_].Data(), j),Form("Num. Lepton %s (%s)"                ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 20  , 0    , 20  ); h_cf_sys_NLeptons[index_][j]->Sumw2();         
      h_cf_sys_Lep1pt[index_][j]         = new TH1D(Form("_h_cf_sys_Lep1pt_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Leading Lepton pT %s (%s)"          ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_Lep1pt[index_][j]->Sumw2();           
      h_cf_sys_Lep1phi[index_][j]        = new TH1D(Form("_h_cf_sys_Lep1phi_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Leading Lepton phi %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_Lep1phi[index_][j]->Sumw2();          
      h_cf_sys_Lep1eta[index_][j]        = new TH1D(Form("_h_cf_sys_Lep1eta_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Leading Lepton eta %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50  , -2.5 , 2.5 ); h_cf_sys_Lep1eta[index_][j]->Sumw2();          
      h_cf_sys_Lep2pt[index_][j]         = new TH1D(Form("_h_cf_sys_Lep2pt_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Second Leading Lepton pT %s (%s)"   ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_Lep2pt[index_][j]->Sumw2();           
      h_cf_sys_Lep2phi[index_][j]        = new TH1D(Form("_h_cf_sys_Lep2phi_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Second Leading Lepton phi %s (%s)"  ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_Lep2phi[index_][j]->Sumw2();      
      h_cf_sys_Lep2eta[index_][j]        = new TH1D(Form("_h_cf_sys_Lep2eta_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Second Leading Lepton eta %s (%s)"  ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50  , -2.5 , 2.5 ); h_cf_sys_Lep2eta[index_][j]->Sumw2();       
      h_cf_sys_dilep_inv_mass[index_][j] = new TH1D(Form("_h_cf_sys_dilep_inv_mass_%s_%d_",v_SystFullName[index_].Data(), j),Form("Invariant Mass of Dilepton %s (%s)" ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_dilep_inv_mass[index_][j]->Sumw2();
      h_cf_sys_Jet1pt[index_][j]         = new TH1D(Form("_h_cf_sys_Jet1pt_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Leading Jet pT %s (%s)"             ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_Jet1pt[index_][j]->Sumw2();       
      h_cf_sys_Jet1phi[index_][j]        = new TH1D(Form("_h_cf_sys_Jet1phi_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Leading Jet phi %s (%s)"            ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_Jet1phi[index_][j]->Sumw2();       
      h_cf_sys_Jet1eta[index_][j]        = new TH1D(Form("_h_cf_sys_Jet1eta_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Leading Jet eta %s (%s)"            ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50  , -2.5 , 2.5 ); h_cf_sys_Jet1eta[index_][j]->Sumw2();       
      h_cf_sys_Jet2pt[index_][j]         = new TH1D(Form("_h_cf_sys_Jet2pt_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Second Leading Jet pT %s (%s)"      ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_Jet2pt[index_][j]->Sumw2();        
      h_cf_sys_Jet2phi[index_][j]        = new TH1D(Form("_h_cf_sys_Jet2phi_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Second Leading Jet phi %s (%s)"     ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_Jet2phi[index_][j]->Sumw2();       
      h_cf_sys_Jet2eta[index_][j]        = new TH1D(Form("_h_cf_sys_Jet2eta_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Second Leading Jet eta %s (%s)"     ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50  , -2.5 , 2.5 ); h_cf_sys_Jet2eta[index_][j]->Sumw2();      
      h_cf_sys_NJets[index_][j]          = new TH1D(Form("_h_cf_sys_NJets_%s_%d_"         ,v_SystFullName[index_].Data(), j),Form("Num. of Jets %s (%s)"               ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 20  , 0    , 20  ); h_cf_sys_NJets[index_][j]->Sumw2();        
      h_cf_sys_metpt[index_][j]          = new TH1D(Form("_h_cf_sys_metpt_%s_%d_"         ,v_SystFullName[index_].Data(), j),Form("MET %s (%s)"                        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_metpt[index_][j]->Sumw2();        
      h_cf_sys_metphi[index_][j]         = new TH1D(Form("_h_cf_sys_metphi_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("MET Phi %s (%s)"                    ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_metphi[index_][j]->Sumw2();       
      h_cf_sys_Nbjets[index_][j]         = new TH1D(Form("_h_cf_sys_Nbjets_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Num. of b-jets %s (%s)"             ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 20  , 0    , 20  ); h_cf_sys_Nbjets[index_][j]->Sumw2();        
      h_cf_sys_NPV[index_][j]            = new TH1D(Form("_h_cf_sys_NPV_%s_%d_"           ,v_SystFullName[index_].Data(), j),Form("Num. of Primary Vertex %s (%s)"     ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 100 , 0    , 100 ); h_cf_sys_NPV[index_][j]->Sumw2();

      h_sys_Lep1pt[index_][j]  = new TH1D(Form("h_sys_Lep1pt_%s_%d"   ,v_SystFullName[index_].Data(), j), Form("Leading Lepton pT %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Lep1pt[index_][j]->Sumw2(); 
      h_sys_Lep2pt[index_][j]  = new TH1D(Form("h_sys_Lep2pt_%s_%d_"  ,v_SystFullName[index_].Data(), j), Form("Second Leading Lepton pT %s (%s)" ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Lep2pt[index_][j]->Sumw2();
      h_sys_Lep1eta[index_][j] = new TH1D(Form("h_sys_Lep1eta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Leading Lepton Eta    %s (%s)"    ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Lep1eta[index_][j]->Sumw2();
      h_sys_Lep2eta[index_][j] = new TH1D(Form("h_sys_Lep2eta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Second Leading Lepton Eta %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Lep2eta[index_][j]->Sumw2();
      h_sys_Lep1phi[index_][j] = new TH1D(Form("h_sys_Lep1phi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Leading Lepton Phi %s (%s)"       ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Lep1phi[index_][j]->Sumw2();
      h_sys_Lep2phi[index_][j] = new TH1D(Form("h_sys_Lep2phi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Second Leading Lepton Phi %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Lep2phi[index_][j]->Sumw2();
   
      h_sys_Muonpt[index_][j]  = new TH1D(Form("h_sys_Muonpt_%s_%d_"  ,v_SystFullName[index_].Data(), j), Form("Muon pT %s (%s)"         ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Muonpt[index_][j]->Sumw2();
      h_sys_Elecpt[index_][j]  = new TH1D(Form("h_sys_Elecpt_%s_%d_"  ,v_SystFullName[index_].Data(), j), Form("Electron pT %s (%s)"     ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Elecpt[index_][j]->Sumw2();
      h_sys_Muoneta[index_][j] = new TH1D(Form("h_sys_Muoneta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Muon Eta %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Muoneta[index_][j]->Sumw2();
      h_sys_Eleceta[index_][j] = new TH1D(Form("h_sys_Eleceta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Electron Eta %s (%s)"    ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Eleceta[index_][j]->Sumw2();
      h_sys_Muonphi[index_][j] = new TH1D(Form("h_sys_Muonphi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Muon Phi %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Muonphi[index_][j]->Sumw2();
      h_sys_Elecphi[index_][j] = new TH1D(Form("h_sys_Elecphi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Electron Phi %s (%s)"    ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Elecphi[index_][j]->Sumw2();
   
      h_sys_Jet1pt[index_][j]  = new TH1D(Form("h_sys_Jet1pt_%s_%d_",  v_SystFullName[index_].Data(), j), Form("Leading Jet pT %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Jet1pt[index_][j]->Sumw2();
      h_sys_Jet2pt[index_][j]  = new TH1D(Form("h_sys_Jet2pt_%s_%d_",  v_SystFullName[index_].Data(), j), Form("Second Leading Jet pT %s (%s)" ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Jet2pt[index_][j]->Sumw2();
      h_sys_Jet1eta[index_][j] = new TH1D(Form("h_sys_Jet1eta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Leading Jet Eta %s (%s)"       ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Jet1eta[index_][j]->Sumw2();
      h_sys_Jet2eta[index_][j] = new TH1D(Form("h_sys_Jet2eta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Second Leading Jet Eta %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Jet2eta[index_][j]->Sumw2();
      h_sys_Jet1phi[index_][j] = new TH1D(Form("h_sys_Jet1phi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Leading Jet Phi %s (%s)"       ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Jet1phi[index_][j]->Sumw2();
      h_sys_Jet2phi[index_][j] = new TH1D(Form("h_sys_Jet2phi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Second Leading Jet Phi %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Jet2phi[index_][j]->Sumw2();
    
     
      h_sys_METpt[index_][j]  = new TH1D(Form("h_sys_METpt_%s_%d_"  ,v_SystFullName[index_].Data(), j), Form("MET pT %s (%s)" ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 200, 0.0, 200); h_sys_METpt[index_][j]->Sumw2();
      h_sys_METphi[index_][j] = new TH1D(Form("h_sys_METphi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("MET Phi %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_METphi[index_][j]->Sumw2();
   
      h_sys_DiLepMass[index_][j] = new TH1D(Form("h_sys_DiLepMass_%s_%d_"     ,v_SystFullName[index_].Data(), j),Form("Di-Lepton Invariant Mass %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 300, 0.0, 300); h_sys_DiLepMass[index_][j]->Sumw2();
      h_sys_Num_PV[index_][j]    = new TH1D(Form("h_sys_Num_PV_%s_%d_"        ,v_SystFullName[index_].Data(), j),     Form("Num of Primary Vertex after %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 100, 0.0, 100); h_sys_Num_PV[index_][j]->Sumw2();
      h_sys_Num_Jets[index_][j]  = new TH1D(Form("h_sys_Num_Jets_%s_%d_"      ,v_SystFullName[index_].Data(), j), Form("Num. of Jets after %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 20, 0.0, 20); h_sys_Num_Jets[index_][j]->Sumw2();
      h_sys_Num_bJets[index_][j] = new TH1D(Form("h_sys_Num_bJets_%s_%d_"     ,v_SystFullName[index_].Data(), j),Form("Num. of b Jets after %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 20, 0.0, 20); h_sys_Num_bJets[index_][j]->Sumw2();
   } // Cut Flow LOOP
   for (int j =0; j < 13; ++j)
   {
      h_sys_Reco_CPO_[index_][j]         = new TH1D(Form("h_sys_Reco_CPO%d_%s"        ,j+1 ,v_SystFullName[index_].Data() ), Form("CPO%d (%s)",j+1,v_SystFullName[index_].Data()   ), 200, -10, 10); h_sys_Reco_CPO_[index_][j]->Sumw2();
      h_sys_Reco_CPO_ReRange_[index_][j] = new TH1D(Form("h_sys_Reco_CPO%d_%s_ReRange",j+1 ,v_SystFullName[index_].Data() ), Form("CPO%d (%s)",j+1,v_SystFullName[index_].Data()   ), 40, -2, 2); h_sys_Reco_CPO_ReRange_[index_][j]->Sumw2();
   }

   h_sys_Top1Mass_[index_]     = new TH1D(Form("h_sys_Top1Mass_%s"     ,v_SystFullName[index_].Data() ), Form("Top1 Mass (%s)"     ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top1Mass_[index_]->Sumw2(); 
   h_sys_Top1pt_[index_]       = new TH1D(Form("h_sys_Top1pt_%s"       ,v_SystFullName[index_].Data() ), Form("Top1 pt (%s)"       ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top1pt_[index_]->Sumw2(); 
   h_sys_Top1Rapidity_[index_] = new TH1D(Form("h_sys_Top1Rapidity_%s" ,v_SystFullName[index_].Data() ), Form("Top1 Rapidity (%s)" ,v_SystFullName[index_].Data() ), 100, -5, 5); h_sys_Top1Rapidity_[index_]->Sumw2(); 
   h_sys_Top1phi_[index_]      = new TH1D(Form("h_sys_Top1phi_%s"      ,v_SystFullName[index_].Data() ), Form("Top1 phi (%s)"      ,v_SystFullName[index_].Data() ), 24, -1*pi, pi); h_sys_Top1phi_[index_]->Sumw2(); 
   h_sys_Top1Energy_[index_]   = new TH1D(Form("h_sys_Top1Energy_%s"   ,v_SystFullName[index_].Data() ), Form("Top1 Energy (%s)"   ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top1Energy_[index_]->Sumw2(); 
   h_sys_Top2Mass_[index_]     = new TH1D(Form("h_sys_Top2Mass_%s"     ,v_SystFullName[index_].Data() ), Form("Top2 Mass (%s)"     ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top2Mass_[index_]->Sumw2();
   h_sys_Top2pt_[index_]       = new TH1D(Form("h_sys_Top2pt_%s"       ,v_SystFullName[index_].Data() ), Form("Top2 pt (%s)"       ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top2pt_[index_]->Sumw2(); 
   h_sys_Top2Rapidity_[index_] = new TH1D(Form("h_sys_Top2Rapidity_%s" ,v_SystFullName[index_].Data() ), Form("Top2 Rapidity (%s)" ,v_SystFullName[index_].Data() ), 100, -5, 5); h_sys_Top2Rapidity_[index_]->Sumw2(); 
   h_sys_Top2phi_[index_]      = new TH1D(Form("h_sys_Top2phi_%s"      ,v_SystFullName[index_].Data() ), Form("Top2 phi (%s)"      ,v_SystFullName[index_].Data() ), 24, -1*pi, pi); h_sys_Top2phi_[index_]->Sumw2(); 
   h_sys_Top2Energy_[index_]   = new TH1D(Form("h_sys_Top2Energy_%s"   ,v_SystFullName[index_].Data() ), Form("Top2 Energy (%s)"   ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top2Energy_[index_]->Sumw2(); 

   h_sys_TopMass_[index_]       = new TH1D(Form("h_sys_TopMass_%s"       ,v_SystFullName[index_].Data() ), Form("Top Mass (%s)"       ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_TopMass_[index_]->Sumw2(); 
   h_sys_Toppt_[index_]         = new TH1D(Form("h_sys_Toppt_%s"         ,v_SystFullName[index_].Data() ), Form("Top pt (%s)"         ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Toppt_[index_]->Sumw2(); 
   h_sys_TopRapidity_[index_]   = new TH1D(Form("h_sys_TopRapidity_%s"   ,v_SystFullName[index_].Data() ), Form("Top Rapidity (%s)"   ,v_SystFullName[index_].Data() ), 100, -5, 5); h_sys_TopRapidity_[index_]->Sumw2(); 
   h_sys_Topphi_[index_]        = new TH1D(Form("h_sys_Topphi_%s"        ,v_SystFullName[index_].Data() ), Form("Top phi (%s)"        ,v_SystFullName[index_].Data() ), 24, -1*pi, pi); h_sys_Topphi_[index_]->Sumw2(); 
   h_sys_TopEnergy_[index_]     = new TH1D(Form("h_sys_TopEnergy_%s"     ,v_SystFullName[index_].Data() ), Form("Top Energy (%s)"     ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_TopEnergy_[index_]->Sumw2(); 
   h_sys_AnTopMass_[index_]     = new TH1D(Form("h_sys_AnTopMass_%s"     ,v_SystFullName[index_].Data() ), Form("AnTop Mass (%s)"     ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_AnTopMass_[index_]->Sumw2();
   h_sys_AnToppt_[index_]       = new TH1D(Form("h_sys_AnToppt_%s"       ,v_SystFullName[index_].Data() ), Form("AnTop pt (%s)"       ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_AnToppt_[index_]->Sumw2(); 
   h_sys_AnTopRapidity_[index_] = new TH1D(Form("h_sys_AnTopRapidity_%s" ,v_SystFullName[index_].Data() ), Form("AnTop Rapidity (%s)" ,v_SystFullName[index_].Data() ), 100, -5, 5); h_sys_AnTopRapidity_[index_]->Sumw2(); 
   h_sys_AnTopphi_[index_]      = new TH1D(Form("h_sys_AnTopphi_%s"      ,v_SystFullName[index_].Data() ), Form("AnTop phi (%s)"      ,v_SystFullName[index_].Data() ), 24, -1*pi, pi); h_sys_AnTopphi_[index_]->Sumw2(); 
   h_sys_AnTopEnergy_[index_]   = new TH1D(Form("h_sys_AnTopEnergy_%s"   ,v_SystFullName[index_].Data() ), Form("AnTop Energy (%s)"   ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_AnTopEnergy_[index_]->Sumw2(); 

}

void ssb_analysis::End()
{
   fout->Write();
   fout->Close();
   /////////////////////////////////
   /// For All-in-One Systematic ///
   /////////////////////////////////
   if (isAllSyst == true)
   {
      for (int i = 0; i < v_SystFullName.size(); ++i )
      {
         a_fout[i]->Write();
         a_fout[i]->Close();
      }
   }
}

void ssb_analysis::SetInputFileName( char *inname )
{
   FileName_ = inname;
}
void ssb_analysis::SetOutputFileName(char *outname)
{   
   outfile = outname;
   TString rm_rootout = TString(outfile).ReplaceAll(".root","");
   v_outName.clear();
   if (isAllSyst == true)
   {
      for (int i =0; i < v_SystFullName.size(); ++i)
      { TString on = Form("%s_%s.root",rm_rootout.Data(),v_SystFullName[i].Data());  
        v_outName.push_back(on);  
      }
   }
}

void ssb_analysis::GetTotalEvent()
{
   totalEvent = 0;
   cout << " FileName_ ? at Get TotalEvent () "<< FileName_ << endl;
   if ( TString( CenOfE ).Contains( "8TeV" ) )
   {
      if      ( TString(FileName_).Contains( "TTJets_FullLept" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 1 );}
      else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { totalEvent = SSBConfReader->GetNumber("TotalEvent", 2 );}
      else if ( TString(FileName_).Contains( "DYJetsToLL_M_50" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 3 );}
      else if ( TString(FileName_).Contains( "TTJets_HadronicMG" ) )   { totalEvent = SSBConfReader->GetNumber("TotalEvent", 4 );}
      else if ( TString(FileName_).Contains( "TTJets_SemiLeptMG" ) )   { totalEvent = SSBConfReader->GetNumber("TotalEvent", 5 );}
      else if ( TString(FileName_).Contains( "T_tW-channel" ) )        { totalEvent = SSBConfReader->GetNumber("TotalEvent", 6 );}
      else if ( TString(FileName_).Contains( "Tbar_tW-channel" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 7 );}
      else if ( TString(FileName_).Contains( "WW" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 8 );}
      else if ( TString(FileName_).Contains( "WZ" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 9 );}
      else if ( TString(FileName_).Contains( "ZZ" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 10);}
      else if ( TString(FileName_).Contains( "WJetsToLNu" ) )          { totalEvent = SSBConfReader->GetNumber("TotalEvent", 11);}
      else if ( TString(FileName_).Contains( "Data" ) )                { totalEvent = SSBConfReader->GetNumber("TotalEvent", 12);}
      else { cout << " File Name Error !! at GetTotalEvent()  " << endl; }
   }
   else if ( TString( CenOfE ).Contains( "13TeV" ) )
   {
      cout << "CENOFE works well !! " << endl;
      if      ( TString(FileName_).Contains( "TTJets" ) )              { 
          if      (TString(FileName_).Contains( "FSRUp" ) )      { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 1 ); }
          else if (TString(FileName_).Contains( "FSRDown" ) )    { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 2 ); }
          else if (TString(FileName_).Contains( "ISRUp" ) )      { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 3 ); }
          else if (TString(FileName_).Contains( "ISRDown" ) )    { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 4 ); }
          else if (TString(FileName_).Contains( "TuneUp" ) )     { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 5 ); }
          else if (TString(FileName_).Contains( "TuneDown" ) )   { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 6 ); }
          else if (TString(FileName_).Contains( "Herwig" ) )     { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 7 ); }
          else if (TString(FileName_).Contains( "EvtGen" ) )     { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 8 ); }
          else if (TString(FileName_).Contains( "hdampUp" ) )    { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 9 ); }
          else if (TString(FileName_).Contains( "hdampDown" ) )  { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 10 ); }
          else if (TString(FileName_).Contains( "GluoneMoveCRTune_erdON" ) ) { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 11 ); }
          else if (TString(FileName_).Contains( "GluoneMoveCRTune" ) )       { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 12 ); }
          else if (TString(FileName_).Contains( "QCDCRTune" ) )              { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 13 ); }
          else if (TString(FileName_).Contains( "erdON" ) )                  { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 14 ); }
          else if (TString(FileName_).Contains( "TopMass1735" ) )            { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 15 ); }
          else if (TString(FileName_).Contains( "TopMass1715" ) )            { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 16 ); }
          else if (TString(FileName_).Contains( "TopMass1755" ) )            { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 17 ); }
          else if (TString(FileName_).Contains( "TopMass1695" ) )            { totalEvent = SSBConfReader->GetNumber("TTJets_Sys", 18 ); }
          else { totalEvent = SSBConfReader->GetNumber("TotalEvent", 1 ); }// Defalut is Powheg 
      }
      else if ( TString(FileName_).Contains( "WJetsToLNu" ) ){
         if ( TString(FileName_).Contains( "WJetsToLNu_Mad" ) ) { totalEvent = SSBConfReader->GetNumber("TotalEvent", 14 ); }
         else { totalEvent = SSBConfReader->GetNumber("TotalEvent", 2 );} // WJetsToLNu
      }
      else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { totalEvent = SSBConfReader->GetNumber("TotalEvent", 3  );}
      else if ( TString(FileName_).Contains( "DYJetsToLL_M_50" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 4  );}
      else if ( TString(FileName_).Contains( "ST_tW_top" ) )           { totalEvent = SSBConfReader->GetNumber("TotalEvent", 5  );}
      else if ( TString(FileName_).Contains( "ST_tW_antitop" ) )       { totalEvent = SSBConfReader->GetNumber("TotalEvent", 6  );}
      else if ( TString(FileName_).Contains( "WW" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 7  );}
      else if ( TString(FileName_).Contains( "WZ" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 8  );}
      else if ( TString(FileName_).Contains( "ZZ" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 9  );}
      else if ( TString(FileName_).Contains( "TTbar_WJetToLNu" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 10 );}// TTW
      else if ( TString(FileName_).Contains( "TTbar_WQQ" ) )           { totalEvent = SSBConfReader->GetNumber("TotalEvent", 11 );}// TTW
      else if ( TString(FileName_).Contains( "TTbar_ZToLLNuNu" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 12 );}// TTZ
      else if ( TString(FileName_).Contains( "TTbar_ZQQ" ) )           { totalEvent = SSBConfReader->GetNumber("TotalEvent", 13 );}// TTZ
      else if ( TString(FileName_).Contains( "Z1Jet_LL" ) )            { totalEvent = SSBConfReader->GetNumber("TotalEvent", 15 );}
      else if ( TString(FileName_).Contains( "Z1Jet_EE" ) )            { totalEvent = SSBConfReader->GetNumber("TotalEvent", 16 );}
      else if ( TString(FileName_).Contains( "Z1Jet_TauTau" ) )        { totalEvent = SSBConfReader->GetNumber("TotalEvent", 17 );}
      else if ( TString(FileName_).Contains( "Data" ) )                { totalEvent = SSBConfReader->GetNumber("TotalEvent", 18 );}
      else { cout << " File Name Error !! at GetTotalEvent()  " << endl; }
   }
   else { cout << "CenOfE Error at GetTotalEvent() " << endl; }
   cout << "totalEvent ? " << totalEvent << endl;
}

// MC scale factor function
void ssb_analysis::MCSF()
{

//   double lumi = 19.6*1000;
//   double lumi = Lumi*1000;
//   double lumi = 40.24;
//   double lumi = 2.11*1000;
   double lumi = Lumi/1000000;
   double br   = 0.0159;
  // cout << "Lumi ? " << Lumi/1000000 << endl;
   cout << "lumi ? " << lumi<< endl;
   if ( TString(UsingTotEnv).Contains( "False" ) || TString(UsingTotEnv).Contains( "false" ) ) 
   {
      
      if ( TString( CenOfE ).Contains( "8TeV" ) ) // for 8 TeV
      {
         if      ( TString(FileName_).Contains( "TTJets_FullLept"     ) ) { mc_sf_ = 25.3*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { mc_sf_ = 860.5*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 3532.8*  lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "TTJets_HadronicMG"   ) ) { mc_sf_ = 106.9*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "TTJets_SemiLeptMG"   ) ) { mc_sf_ = 103*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "T_tW-channel"        ) ) { mc_sf_ = 11.2*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "Tbar_tW-channel"     ) ) { mc_sf_ = 11.2*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WW"                  ) ) { mc_sf_ = 5.8*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WZ"                  ) ) { mc_sf_ = 22.4*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "ZZ"                  ) ) { mc_sf_ = 9.0*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WJetsToLNu"          ) ) { mc_sf_ = 37509.0* lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "Data"                ) ) { mc_sf_ = 1;                                }
         else { cout << " File Name Error !! at MCSF() " << endl; }
      }
      else if ( TString( CenOfE ).Contains( "13TeV" ) ) // for 13 TeV
      {
         if      ( TString(FileName_).Contains( "TTJets"              ) ) { mc_sf_ = 831.76*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WJetsToLNu"          ) ) { mc_sf_ = 61526*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { mc_sf_ = 18810*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 5941.0*   lumi / NtupletotalEvent; }
         //else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 6025.2*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "ST_tW_top"           ) ) { mc_sf_ = 35.6*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "ST_tW_antitop"       ) ) { mc_sf_ = 35.6*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WW"                  ) ) { mc_sf_ = 118.7*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WZ"                  ) ) { mc_sf_ = 65.9*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "ZZ"                  ) ) { mc_sf_ = 31.8*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "TTbar_WJetToLNu"     ) ) { mc_sf_ = 0.2043*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "TTbar_WQQ"           ) ) { mc_sf_ = 0.4062*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "TTbar_ZToLLNuNu"     ) ) { mc_sf_ = 0.2529*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "TTbar_ZQQ"           ) ) { mc_sf_ = 0.5297*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "Z1Jet"               ) ) { mc_sf_ = 1921.8*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "Data"                ) ) { mc_sf_ = 1;                                 }
         else { cout << " File Name Error !! at MCSF()" << endl;}
      }
      else { cout << "Center Of Energy Error in MCSF()" << endl;  }
   }
   else if ( TString(UsingTotEnv).Contains( "True" ) || TString(UsingTotEnv).Contains( "true" ) )
   {
      if ( TString( CenOfE ).Contains( "8TeV" )  )
      {
         if      ( TString(FileName_).Contains( "TTJets_FullLept"     ) ) { mc_sf_ = 25.3*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { mc_sf_ = 860.5*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 3532.8*  lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "TTJets_HadronicMG"   ) ) { mc_sf_ = 106.9*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "TTJets_SemiLeptMG"   ) ) { mc_sf_ = 103*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "T_tW-channel"        ) ) { mc_sf_ = 11.2*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "Tbar_tW-channel"     ) ) { mc_sf_ = 11.2*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WW"                  ) ) { mc_sf_ = 5.8*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WZ"                  ) ) { mc_sf_ = 22.4*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "ZZ"                  ) ) { mc_sf_ = 9.0*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WJetsToLNu"          ) ) { mc_sf_ = 37509.0* lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "Data"                ) ) { mc_sf_ = 1;                          }
         else { cout << " File Name Error !! " << endl; }
      }
      else if ( TString(CenOfE).Contains("13TeV") )
      {
         if      ( TString(FileName_).Contains( "TTJets"              ) ) { mc_sf_ = 831.76*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WJetsToLNu"          ) ) { mc_sf_ = 61526*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { mc_sf_ = 18810*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 5941.0*   lumi / totalEvent; }
         //else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 6025.2*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "ST_tW_top"           ) ) { mc_sf_ = 35.6*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "ST_tW_antitop"       ) ) { mc_sf_ = 35.6*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WW"                  ) ) { mc_sf_ = 118.7*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WZ"                  ) ) { mc_sf_ = 65.9*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "ZZ"                  ) ) { mc_sf_ = 31.8*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "TTbar_WJetToLNu"     ) ) { mc_sf_ = 0.2043*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "TTbar_WQQ"           ) ) { mc_sf_ = 0.4062*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "TTbar_ZToLLNuNu"     ) ) { mc_sf_ = 0.2529*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "TTbar_ZQQ"           ) ) { mc_sf_ = 0.5297*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "Z1Jet"               ) ) { mc_sf_ = 1921.8*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "Data"                ) ) { mc_sf_ = 1;                                 }
         else { cout << " File Name Error !! at MCSF()" << endl;}

      }
      else { cout << "Center Of Energy Error in MCSF()" << endl; }
   } 
   else {cout << "MCSF error " << endl; }
   cout << "mc_sf_ : " << mc_sf_ << endl;
}

// Apply MC SF To Event //
void ssb_analysis::MCSFApply()
{
   evt_weight_beforemcsf_ =1; // Initailize evt_weight_beforemcsf_ //
   evt_weight_beforemcsf_ = evt_weight_; // keep event weight // 
   if ( !TString(FileName_).Contains( "Data") ){ evt_weight_ = evt_weight_*mc_sf_; } // apply MC scale factor // 
   else {evt_weight_ = 1;}
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      { 
         if ( !TString(FileName_).Contains( "Data") )
         {
            v_SystEvt[i]=evt_weight_beforemcsf_*mc_sf_;
         }
         else {
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   }
}
void ssb_analysis::GenWeightApply()
{
   double genweight = 1.0;
   evt_weight_beforegenweight_ =1; // Initailize evt_weight_beforegenweight_ //
   evt_weight_beforegenweight_ = evt_weight_; // keep event weight //
   if ( !TString(FileName_).Contains( "Data") ){
      if (Gen_EventWeight > 0.0){genweight =1;}
      else {genweight =-1;} 
      evt_weight_ = evt_weight_*genweight;
   }
   else {evt_weight_ = 1;}
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      { 
         if ( !TString(FileName_).Contains( "Data") )
         {
            v_SystEvt[i]=evt_weight_beforegenweight_*genweight;
         }
         else {
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   }
}
void ssb_analysis::PDFWeightApply()
{
   double pdfweight = 1.0;
   evt_weight_beforepdfweight_ = evt_weight_;
   if ( !TString(FileName_).Contains("Data"))
   { 
      if ( TString(FileName_).Contains("TTJets_") ){
         if ( PDFSys == 0 ) { pdfweight = 1.0; }// central //
         else if ( PDFSys > 0 )  
         {
             int index_pdf = PDFSys + 8;
             pdfweight = LHE_Weight->at(index_pdf)/LHE_Central;
             //cout << "LHE_Central : "<< LHE_Central <<" LHE_Id["<< index_pdf <<"] : " << LHE_Id->at(index_pdf) << endl;
         }
         else {cout << "Somethig Wrong in the PDF sys ... " << endl;}
      }
      else { pdfweight =1; }
      evt_weight_ = evt_weight_*pdfweight;
   }
   else { evt_weight_ = 1; }
   if ( isAllSyst == true )
   {
      double pdfweight_allsys = 1.0;
      for( int i =0; i < v_SystFullName.size(); ++i )
      { 
         pdfweight_allsys = 1.0;
         if ( !TString(FileName_).Contains( "Data") )
         {
            if (TString(FileName_).Contains("TTJets_")){
               /// PDF Index ///
               if (v_SystFullName.at(i).Contains("PDF_")){
                  TString SysPDF_ = v_SystFullName[i];
                  SysPDF_ = SysPDF_.ReplaceAll("PDF_","");
                  /*cout 
                  << "SysPDF_ : " << SysPDF_ << " SysPDF_ Int_t : " << SysPDF_.Atoi() 
                  << " LHE_Id : "<< LHE_Id->at(SysPDF_.Atoi() + 8)
                  << " EvtWeight : "<< LHE_Weight->at(SysPDF_.Atoi() + 8)/LHE_Central
                  << endl;*/
                  pdfweight_allsys = LHE_Weight->at(SysPDF_.Atoi() + 8)/LHE_Central;
               }
               else {pdfweight_allsys = 1;}
            } 

            v_SystEvt[i]=v_SystEvt[i]*pdfweight_allsys;
         }
         else {
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   }

}
void ssb_analysis::FactRenoApply()
{
   double factrenosys = 1.0;
   evt_weight_beforefactrenoweight_ =1; // Initailize evt_weight_beforefactrenoweight_ //
   evt_weight_beforefactrenoweight_ = evt_weight_; // keep event weight //
   if (LHE_Weight->size() == 0 && !TString(FileName_).Contains("Data") ) {cout << "There is no LHE_Weight!!!" <<endl;  return;}
   if (TString(FileName_).Contains("TTJets")) {
      if (FactRenoSys == 0) { factrenosys = 1; }
      else if (FactRenoSys == 1){ factrenosys = LHE_Weight->at(1)/LHE_Central; }
      else if (FactRenoSys == 2){ factrenosys = LHE_Weight->at(2)/LHE_Central; }
      else if (FactRenoSys == 3){ factrenosys = LHE_Weight->at(3)/LHE_Central; }
      else if (FactRenoSys == 4){ factrenosys = LHE_Weight->at(4)/LHE_Central; }
      else if (FactRenoSys == 5){ factrenosys = LHE_Weight->at(5)/LHE_Central; }
      else if (FactRenoSys == 6){ factrenosys = LHE_Weight->at(6)/LHE_Central; }
      else if (FactRenoSys == 7){ factrenosys = LHE_Weight->at(7)/LHE_Central; }
      else if (FactRenoSys == 8){ factrenosys = LHE_Weight->at(8)/LHE_Central; }
      else { cout << "Check out your FactRenoSys " << FactRenoSys << endl;  }
   }
   else { factrenosys =1; }
   evt_weight_ = evt_weight_*(factrenosys);
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      { 
         if ( !TString(FileName_).Contains( "Data") )
         {
            if ( !TString(FileName_).Contains( "TTJets") ){
               v_SystEvt[i]=v_SystEvt[i];
            }
            else {
               if (v_SystFullName[i].Contains("FactReno_1")){ v_SystEvt[i]=v_SystEvt[i]*LHE_Weight->at(1)/LHE_Central;}
               else if (v_SystFullName[i].Contains("FactReno_2")){v_SystEvt[i]=v_SystEvt[i]*LHE_Weight->at(2)/LHE_Central;}
               else if (v_SystFullName[i].Contains("FactReno_3")){v_SystEvt[i]=v_SystEvt[i]*LHE_Weight->at(3)/LHE_Central;}
               else if (v_SystFullName[i].Contains("FactReno_4")){v_SystEvt[i]=v_SystEvt[i]*LHE_Weight->at(4)/LHE_Central;}
               else if (v_SystFullName[i].Contains("FactReno_5")){v_SystEvt[i]=v_SystEvt[i]*LHE_Weight->at(5)/LHE_Central;}
               else if (v_SystFullName[i].Contains("FactReno_6")){v_SystEvt[i]=v_SystEvt[i]*LHE_Weight->at(6)/LHE_Central;}
               else if (v_SystFullName[i].Contains("FactReno_7")){v_SystEvt[i]=v_SystEvt[i]*LHE_Weight->at(7)/LHE_Central;}
               else if (v_SystFullName[i].Contains("FactReno_8")){v_SystEvt[i]=v_SystEvt[i]*LHE_Weight->at(8)/LHE_Central;}
               else {v_SystEvt[i]=v_SystEvt[i];}
            }
         }
         else { // Case of Data
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   }
}
void ssb_analysis::FragmentApply()
{
   double fragmentsys = 1.0;
   evt_weight_beforefragmentweight_ =1; // Initailize evt_weight_beforefragmentweight_ //
   evt_weight_beforefragmentweight_ = evt_weight_; // keep event weight //
   if (TString(FileName_).Contains("TTJets")) {
      if (FragmentSys == "nominal") { fragmentsys = 1.0; }
      else if (FragmentSys == "up"){ fragmentsys = Frag_Up_Weight; }
      else if (FragmentSys == "central"){ fragmentsys = Frag_Cen_Weight; }
      else if (FragmentSys == "down"){ fragmentsys = Frag_Down_Weight; }
      else if (FragmentSys == "peterson"){ fragmentsys = Frag_Peterson_Weight; }
      else { cout << "Check out your FragmentSys " << FragmentSys << endl;  }
   }
   else { fragmentsys =1; }
   evt_weight_ = evt_weight_*(fragmentsys);
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      { 
         if ( !TString(FileName_).Contains( "Data") )
         {
            if ( !TString(FileName_).Contains( "TTJets") ){
               v_SystEvt[i]=v_SystEvt[i];
            }
            else {
               if (v_SystFullName[i].Contains("FragmentCentral"))      { v_SystEvt[i]=v_SystEvt[i]*Frag_Cen_Weight;}
               else if (v_SystFullName[i].Contains("FragmentUp"))      { v_SystEvt[i]=v_SystEvt[i]*Frag_Up_Weight;}
               else if (v_SystFullName[i].Contains("FragmentDown"))    { v_SystEvt[i]=v_SystEvt[i]*Frag_Down_Weight;}
               else if (v_SystFullName[i].Contains("FragmentPeterson")){ v_SystEvt[i]=v_SystEvt[i]*Frag_Peterson_Weight;}
               else {v_SystEvt[i]=v_SystEvt[i];}
            }
         }
         else { // Case of Data
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   }
}
/// Systematic for Decay tables
void ssb_analysis::DecayTableApply()
{
   double dectabsys = 1.0;
   evt_weight_beforedectabweight_ =1; // Initailize evt_weight_beforedectabweight_ //
   evt_weight_beforedectabweight_ = evt_weight_; // keep event weight //
   if (TString(FileName_).Contains("TTJets")) {
      if (DecayTableSys == "central") { dectabsys = 1.0; }
      else if (DecayTableSys == "up"){ dectabsys = Semilep_BrUp_Weight; }
      else if (DecayTableSys == "down"){ dectabsys = Semilep_BrDown_Weight; }
      else { cout << "Check out your DecayTableSys " << DecayTableSys << endl;  }
   }
   else { dectabsys =1; }
   evt_weight_ = evt_weight_*(dectabsys);
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      { 
         if ( !TString(FileName_).Contains( "Data") )
         {
            if ( !TString(FileName_).Contains( "TTJets") ){
               v_SystEvt[i]=v_SystEvt[i];
            }
            else {
               if (v_SystFullName[i].Contains("DecayTableUp"))      { v_SystEvt[i]=v_SystEvt[i]*Semilep_BrUp_Weight;}
               else if (v_SystFullName[i].Contains("DecayTableDown"))    { v_SystEvt[i]=v_SystEvt[i]*Semilep_BrDown_Weight;}
               else {v_SystEvt[i]=v_SystEvt[i];}
            }
         }
         else { // Case of Data
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   }   
}

// Apply Trigger SF To Event //
void ssb_analysis::TriggerSFApply()
{
   double triggersf_ =1;
   evt_weight_beforeTrigger_ =1; // Initailize evt_weight_beforeTrigger_ //
   evt_weight_beforeTrigger_ = evt_weight_; // keep event weight // 
   if ( !TString(FileName_).Contains( "Data") )
    {  //evt_weight_ = evt_weight_*triggersf_;
/*      if ( TString(Decaymode).Contains( "dielec" ) ) { triggersf_ = 0.996; }
      else if (TString(Decaymode).Contains( "muel" )  )  { triggersf_ = 0.990;}
      else if (TString(Decaymode).Contains( "dimu" )  )  { triggersf_ = 0.996;}*/
      if ( TString(Decaymode).Contains( "dielec" ) ) { triggersf_ =SSBEffcal->TrigDiElec_Eff(Lep1,Lep2,TrigSFSys); }
      else if (TString(Decaymode).Contains( "muel" )  )  { //triggersf_ = SSBEffcal->TrigMuElec_Eff(TElectron,TMuon); 
         triggersf_ = SSBEffcal->TrigMuElec_Eff(Lep1,Lep2,TrigSFSys);
      }
      else if (TString(Decaymode).Contains( "dimu" )  )  { triggersf_ = SSBEffcal->TrigDiMuon_Eff(Lep1,Lep2,TrigSFSys);}
      else { cout << "CHECK OUT TRIGGER EFF"<< endl;}
      evt_weight_ = evt_weight_*triggersf_;
   }// apply Trigger scale factor //
   else { evt_weight_ = evt_weight_; }
   if ( isAllSyst == true )
   {
      double trigsfcen_ = 1;
      double trigsfup_ = 1;
      double trigsfdn_ = 1;
      if ( TString(Decaymode).Contains( "dimu" ) )  
      { 
         trigsfcen_ = SSBEffcal->TrigDiMuon_Eff(Lep1,Lep2,"central"); 
         trigsfup_  = SSBEffcal->TrigDiMuon_Eff(Lep1,Lep2,"up"); 
         trigsfdn_  = SSBEffcal->TrigDiMuon_Eff(Lep1,Lep2,"down"); 
      }
      else if ( TString(Decaymode).Contains( "dielec" ) ) 
      {
         trigsfcen_ = SSBEffcal->TrigDiElec_Eff(Lep1,Lep2,"central");
         trigsfup_  = SSBEffcal->TrigDiElec_Eff(Lep1,Lep2,"up"); 
         trigsfdn_  = SSBEffcal->TrigDiElec_Eff(Lep1,Lep2,"down"); 
      }
      else if ( TString(Decaymode).Contains( "muel" ) ) 
      {
         trigsfcen_ = SSBEffcal->TrigMuElec_Eff(Lep1,Lep2,"central");
         trigsfup_  = SSBEffcal->TrigMuElec_Eff(Lep1,Lep2,"up"); 
         trigsfdn_  = SSBEffcal->TrigMuElec_Eff(Lep1,Lep2,"down"); 
      }
      else  {cout << "CHECK OUT TRIGGER EFF ALL SYS"<< endl; }

      for( int i =0; i < v_SystFullName.size(); ++i )
      { 
         if ( !TString(FileName_).Contains( "Data") )
         {
            /// Get TRIGGER SF w SYS //
            if ( v_SystFullName[i].Contains("TrigSF") ){
               if (v_SystFullName[i].Contains("TrigSFUp") ) {v_SystEvt[i] = v_SystEvt[i]*trigsfup_;}
               else if (v_SystFullName[i].Contains("TrigSFDown") ) {v_SystEvt[i] = v_SystEvt[i]*trigsfdn_;}
               else {cout << "CHECK OUT TRIGGER SF SYSTEMATIC !!! " << v_SystFullName[i]<< endl;}
            } 
            else { v_SystEvt[i] = v_SystEvt[i]*trigsfcen_; }
         }
         else {
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   } 
}
void ssb_analysis::PileUpReWeightApply()
{
   evt_weight_beforePileup_ = 1;
   evt_weight_beforePileup_ = evt_weight_; // keep event weight // 
   double puweight_ = 1.;
   double pu_weight_central = puweight->weight( PileUp_Count_Intime );
   double pu_weight_up = puweightup->weight( PileUp_Count_Intime );
   double pu_weight_dn = puweightdn->weight( PileUp_Count_Intime );

   if ( !TString(FileName_).Contains( "Data") )
   {
//      puweight_ = puweight->weight( PileUp_Count_Intime );
      if (TString(PileUpSys).Contains("central") ) { puweight_   = pu_weight_central;}
      else if (TString(PileUpSys).Contains("up") ) { puweight_   = pu_weight_up; }
      else if (TString(PileUpSys).Contains("down") ) { puweight_ = pu_weight_dn; }
      else {  
      cout << "PileUp sys Error ... Defalut is Weight_PileUp ... : " << PileUpSys << endl;
   }
   evt_weight_ = evt_weight_*puweight_; } // apply PileUpReweight //
   else {evt_weight_ = 1;}
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      {
         if( !TString(FileName_).Contains("Data") )
         {
            if ( v_SystFullName[i].Contains("PileUpUp"))
            {
               v_SystEvt[i] = v_SystEvt[i]*pu_weight_up;
            }
            else if ( v_SystFullName[i].Contains("PileUpDown"))
            {
               v_SystEvt[i] = v_SystEvt[i]*pu_weight_dn;
            }
            else {
               v_SystEvt[i]=v_SystEvt[i]*pu_weight_central;
            }
         }
         else { // Not for Data // 
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      }
   }
}
void ssb_analysis::L1PreFireApply()
{
   evt_weight_beforeL1PreFire_ = 1;
   evt_weight_beforeL1PreFire_ = evt_weight_; // keep event weight // 
   double l1prefire_ = 1.;
   double l1prefire_central = L1_PreFire_Central;
   double l1prefire_up      = L1_PreFire_Up;
   double l1prefire_dn      = L1_PreFire_Down;

   if ( !TString(FileName_).Contains( "Data") )
   {
      if (TString(L1PreFireSys).Contains("central") ) { l1prefire_   = l1prefire_central;}
      else if (TString(L1PreFireSys).Contains("up") ) { l1prefire_   = l1prefire_up; }
      else if (TString(L1PreFireSys).Contains("down") ) { l1prefire_ = l1prefire_dn; }
      else if (TString(L1PreFireSys).Contains("none") ) { l1prefire_ = 1.0; }
      else {  
      cout << "PileUp sys Error ... Defalut is Weight_PileUp ... : " << L1PreFireSys << endl;
   }
   evt_weight_ = evt_weight_*l1prefire_; } // apply PileUpReweight //
   else {evt_weight_ = 1;}
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      {
         if( !TString(FileName_).Contains("Data") )
         {
            if ( v_SystFullName[i].Contains("L1PreFireUp"))
            {
               v_SystEvt[i] = v_SystEvt[i]*l1prefire_up;
            }
            else if ( v_SystFullName[i].Contains("L1PreFireDown"))
            {
               v_SystEvt[i] = v_SystEvt[i]*l1prefire_dn;
            }
            else {
               v_SystEvt[i]=v_SystEvt[i]*l1prefire_central;
            }
         }
         else { // Not for Data // 
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      }
   }
}
void ssb_analysis::Weight()
{
   double puweight_ = 1;
   double k_fac_;
   double genweight = 1;
   double pdfweight = 1;

/*   /////////////////////////
   /// To Get PDF Weight ///
   /////////////////////////
   if ( !TString(FileName_).Contains("Data") )
   {
      if (PDFSys == -1 ){ pdfweight =1; }
      else if ( PDFSys == 0 ) { pdfweight = PDFWeight_Cent->at(1); }
      else if ( PDFSys > 0 )  
      {
         if      ( PDFSys%2 == 1 ) { pdfweight = PDFWeight_Var1_Up->at(PDFSys/2); } //odd number is PDFUp
         else if ( PDFSys%2 == 0 ) { pdfweight = PDFWeight_Var1_Down->at(PDFSys/2 -1); } //even number is PDFDown
         else {cout << "Somethig Wrong in ther PDF sys ... " << endl;}
      }
   }*/
   evt_weight_ = 1;
   
   ///////////////////////
   // PileUp Systematic //
   ///////////////////////
   puweight_ = puweight->weight( PileUp_Count_Intime );
   if (TString(PileUpSys).Contains("central") ) { puweight_   = puweight->weight( PileUp_Count_Intime ); }
   else if (TString(PileUpSys).Contains("up") ) { puweight_   = puweightup->weight( PileUp_Count_Intime ); }
   else if (TString(PileUpSys).Contains("down") ) { puweight_ = puweightdn->weight( PileUp_Count_Intime ); }
   else {  
   cout << "PileUp sys Error ... Defalut is Weight_PileUp ... : " << PileUpSys << endl;
   }
   if ( TString(FileName_).Contains( "TTJets"   ) )
   {  
      k_fac_ = 1;
      evt_weight_ = k_fac_*mc_sf_*puweight_*lep_eff*pdfweight;
   }
//   { k_fac_ = 252.89/235.2; evt_weight_ = k_fac_*mc_sf_; }
   else if ( TString(FileName_).Contains( "DYJetsToLL" ) || 
             TString(FileName_).Contains( "WW"         ) ||
             TString(FileName_).Contains( "WZ"         ) ||
             TString(FileName_).Contains( "ZZ"         ) ||
             TString(FileName_).Contains( "WJetsToLNu" ) ||
             TString(FileName_).Contains( "tW" )   )
   {
      k_fac_ =1; 
      if ( TString(FileName_).Contains( "DYJetsToLL" ) || TString(FileName_).Contains( "WJetsToLNu" ) ) 
      {
         if (Gen_EventWeight > 0){genweight =1;}
         else {genweight =-1;}
      }
      evt_weight_ = k_fac_*mc_sf_*puweight_*lep_eff*genweight*pdfweight; 
   }
   else if ( TString(FileName_).Contains( "Data" ) )
   {
      k_fac_ =1; evt_weight_ = 1; puweight_=1; 
   }
   else {cout << "Event-Weight error !!" << endl; evt_weight_ =0; }
}
void ssb_analysis::NumPVCount()
{
   num_pv = 0;
   //cout << " start NumPVCount "<< endl;
   for (int i =0; i < Filter_PV->size(); ++i)
   {
      //cout <<  " Filter_PV->at(i) : " <<  Filter_PV->at(i) << endl;
      if ( Filter_PV->at(i) == true ) {num_pv++;}
   }
}

bool ssb_analysis::EveRun()
{
   bool everun = false;
   if (TString(FileName_).Contains( "Data" ) ) // Only Data
   {
      if ( TString(FileName_).Contains( "Jul" ) ) //Jul-Data
      {
         if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) || TString(Decaymode).Contains( "muel" )  )
         {
            if (Info_RunNumber >= 246908 && Info_RunNumber < 251604 ) { everun = true; }
         }
         else { cout << "YOU DON'T USE THIS FUNCTION , CHECH DECAY MODE" << endl;}
      }
      else 
      {
         if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) || TString(Decaymode).Contains( "muel" ) )
         {
            if (Info_RunNumber >= 251604 && Info_RunNumber <= 251883 ) { everun = true; }
         }
         else { cout << "YOU DON'T USE THIS FUNCTION , CHECH DECAY MODE" << endl;}
      }//prompt
   }
   else { everun = true;} //MC
   return everun;
}
bool ssb_analysis::ChannelIndex()
{
   bool chan_cut = true;
   int decay_idx = 0;
   if ( TString(Decaymode).Contains( "dielec" ) ) decay_idx = 22;
   if ( TString(Decaymode).Contains( "muel" ) )   decay_idx = 24;
   if ( TString(Decaymode).Contains( "dimuon" ) ) decay_idx = 26;

   if ( Channel_Idx != decay_idx && ( TString(FileName_).Contains("other") == false ) && TString(FileName_).Contains("TTJets") ) {chan_cut=false;} // only TTbar -> Di-Lepton Signal.
   if ( Channel_Idx == decay_idx && TString(FileName_).Contains("other") ) {chan_cut = false;} // only TTbar -> Di-Lepton Signal.
   return chan_cut;
}
bool ssb_analysis::METFilterAPP()
{
   bool metfilt_ = true;
   for ( int i = 0; i < METFilter_Name->size(); ++i )
   { 
      TString METFiltName = METFilter_Name->at(i);
      //cout << "METFiltName : " << METFilter_Name->at(i) << endl;
      if ( !(TString(FileName_).Contains("Data")) && 
            TString(METFiltName).Contains("Flag_eeBadScFilter")  ) {continue;}
      for( int j =0; j < v_METFilterName.size(); ++j )
      {
         if ( TString( METFiltName ).Contains(v_METFilterName.at(j)) )
         {
            if ( !(METFilter_isPass->at(i)) ){ metfilt_ = false; }
         } 
      }
   }
   for ( int i = 0; i < METFilterAdd_Name->size(); ++i )
   {
      TString METFiltName = METFilterAdd_Name->at(i);
      if ( !(TString(FileName_).Contains("Data")) && 
            TString(METFiltName).Contains("Flag_eeBadScFilter")  ) {continue;}
      for( int j =0; j < v_METFilterName.size(); ++j )
      {
         if ( TString( METFiltName ).Contains(v_METFilterName.at(j)) )
         {
            if ( !( METFilterAdd_isPass->at(i) ) ){ metfilt_ = false; }
         } 
      }
   }
   return metfilt_;
   //return true;
}

// Trigger Requirement Function
bool ssb_analysis::Trigger()
{
   ////////////////////////////////////////////////////////////////////////////////
   // single muon trigger                                                        //
   // HLT_IsoMu*, HLT_IsoTkMu* ; 뒤의 숫자는 pT threshold                           //
   // double muon trigger                                                        //
   // HLT_Mu(1)_TrkIsoVVL_Mu(2), HLt_Mu(1)_TrkIsoVVL_TkMu(2) ;                   //
   // (1): 1st muon pT, (2) 2nd muon pT, TrkIsoVVL: track isolation requirement  //
   ////////////////////////////////////////////////////////////////////////////////

   /// Variable for Trigger Function
   int ptrigindex;
   bool trigpass;

   ptrigindex =0;
   trigpass = false;
   TString trgName = ""; 
//   if (!TString(FileName_).Contains( "Data" ) ) { trigpass = true; return trigpass;}

   // if data name is DoubleMuon, exclude single muon trigger
   //bool isDoubleMuData = TString(FileName_).Contains("DoubleMuon");

   for(int i =0; i < Trigger_Name->size(); i++)
   {
      //cout << "Ntuple Triggers : " <<  Trigger_Name->at(i) << endl; 
      for(int j = 0; j < trigName.size(); j++)
      {
         trgName = trigName[j];
         
         //if( isDoubleMuData && (trgName.Contains("HLT_IsoMu") || trgName.Contains("HLT_IsoTkMu")) ) continue;
         
         if(TString(Decaymode).Contains( "dimuon" ) )
         {
            if( TString(FileName_).Contains( "Run2016H" ) ) 
            { 
               if(trgName.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") ||
                  trgName.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") ) {continue;}
            }
            else
            {
               if(trgName.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ||
                  trgName.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") ) {continue;}
            }
         }
         if(TString(Decaymode).Contains( "muel" ) )
         {
            if( TString(FileName_).Contains( "Run2016H" ) ) 
            { 
               if(trgName.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") ||
                  trgName.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) {continue;}
            }
            else
            {
               if(trgName.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ||
                  trgName.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) {continue;}
            }
         }

         if( TString( Trigger_Name->at(i) ).Contains( trigName.at(j) ) )//TString clone 
         {
            //cout << "trigName.at(j) ?" << trigName.at(j) << endl;
            if ( ( Trigger_isPass->at(i)  ) && 
                !( Trigger_isError->at(i) ) && 
                 ( Trigger_isRun->at(i) )      ) 
            { ptrigindex = ptrigindex+1; }
         }
      }
   }
   if ( ptrigindex > 0 ) { trigpass = true; }
   return trigpass;
}
// Function of Muon Rocheser Correction //
TLorentzVector* ssb_analysis::ApplyRocCor(TLorentzVector* tmp,int index_)
{
   TLorentzVector* applied = new TLorentzVector();
   double muoncorsf =  1.0;

   if( FileName_.Contains("Data") )
   {
      /// For DATA 
      int charge = Muon_Charge->at(index_); 
      if ( MuEnSys != "none")
      {
         /// Get Central Value ///
         muoncorsf = ssbmucor->kScaleDT(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(),0,0); 
         /// Calc Unc ///
         double totalunc = 0;
         if( MuEnSys != "central" ) {
            // Cal Stat. Unc. //
            TH1D* stat = new TH1D("stat", "", 1000, 0., 10.);
            for (int i = 0; i < 100; ++i) stat->Fill(ssbmucor->kScaleDT(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(),1,i) );
            double statunc = stat->GetRMS() / stat->GetMean();
            delete stat;
            double CorDm, FitDm, deviation, maxDeviation = 1.;
            for(auto i=0; i<5; ++i) {
               CorDm = ssbmucor->kScaleDT(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(),  4, i);
               FitDm = ssbmucor->kScaleDT(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(),  5, i); 
               deviation = (fabs(CorDm-1) > fabs(FitDm-1)) ? CorDm : FitDm;
               if( fabs(deviation-1) > fabs(maxDeviation-1) ) maxDeviation = deviation;
            } 
            //cout << "maxDeviation :" << maxDeviation << endl;
            totalunc = sqrt( pow(statunc,2) + pow(abs(deviation-muoncorsf),2) );
            //cout << "totalunc : " << totalunc<< endl;
            if ( MuEnSys == "up")       {totalunc =  1*totalunc;} 
            else if ( MuEnSys == "down"){totalunc = -1*totalunc;} 
            else {muoncorsf = 1; totalunc =0;cout << "Check Out YOU ApplyRocCor Option !!" << endl; } 
         }
         else {totalunc = 0;}
         muoncorsf = muoncorsf+totalunc;
      }
      else { muoncorsf = 1; }
      //cout << "muoncorsf : " << muoncorsf << endl;
   }
   else {
   /// For MC
      double u1 = Muon_rand1->at(index_);
      double u2 = Muon_rand2->at(index_);
      int nTkLayer = Muon_trackerLayers->at(index_); 
      int charge = Muon_Charge->at(index_); 
      TLorentzVector* genMu = (TLorentzVector*)GenMuon->At(index_); 
      if ( MuEnSys != "none"){
         /// Get Central Value ///
         if (genMu->Pt() > 0.0 && genMu->Eta() == 0){
            muoncorsf = ssbmucor->kScaleFromGenMC(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(), nTkLayer, genMu->Pt(), u1, 0, 0);
         }
         else {
            muoncorsf = ssbmucor->kScaleAndSmearMC(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(), nTkLayer, u1, u2, 0, 0);
         }
         /// Calc Unc ///
         double totalunc = 0;
         if( MuEnSys != "central" ) {
            // Cal Stat. Unc. //
            TH1D* stat = new TH1D("stat", "", 1000, 0., 10.);
            if (genMu->Pt() > 0.0 && genMu->Eta() == 0){
               for (int i = 0; i < 100; ++i) stat->Fill(ssbmucor->kScaleFromGenMC(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(), nTkLayer, genMu->Pt(), u1, 1, i) );
            }
            else {
               muoncorsf = ssbmucor->kScaleAndSmearMC(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(), nTkLayer, u1, u2, 0, 0);
               for (int i = 0; i < 100; ++i) stat->Fill(ssbmucor->kScaleAndSmearMC(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(), nTkLayer, u1, u2, 1, i));
            }
            double statunc = stat->GetRMS() / stat->GetMean();
            delete stat;
            double CorDm, FitDm, deviation, maxDeviation = 1.;
            for(auto i=0; i<5; ++i) {
               if (genMu->Pt() > 0.0 && genMu->Eta() == 0){
                  CorDm = ssbmucor->kScaleFromGenMC(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(), nTkLayer, genMu->Pt(), u1, 4, i);
                  FitDm = ssbmucor->kScaleFromGenMC(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(), nTkLayer, genMu->Pt(), u1, 5, i); 
               }
               else {
                  CorDm = ssbmucor->kScaleAndSmearMC(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(), nTkLayer, u1, u2, 4, i);
                  FitDm = ssbmucor->kScaleAndSmearMC(charge, tmp->Pt(), tmp->Eta(), tmp->Phi(), nTkLayer, u1, u2, 5, i);
               }
               deviation = (fabs(CorDm-1) > fabs(FitDm-1)) ? CorDm : FitDm;
               if( fabs(deviation-1) > fabs(maxDeviation-1) ) maxDeviation = deviation;
            } 
            //cout << "maxDeviation :" << maxDeviation << endl;
            totalunc = sqrt( pow(statunc,2) + pow(abs(deviation-muoncorsf),2) );
            //cout << "totalunc : " << totalunc<< endl;
            if ( MuEnSys == "up")       {totalunc =  1*totalunc;} 
            else if ( MuEnSys == "down"){totalunc = -1*totalunc;} 
            else {muoncorsf = 1; totalunc =0;cout << "Check Out YOU ApplyRocCor Option !!" << endl; } 
         }
         else {totalunc = 0;}
         muoncorsf = muoncorsf+totalunc;
      }
      else { muoncorsf = 1; }
      //cout << "muoncorsf : " << muoncorsf << endl;
   }
   if (muoncorsf != muoncorsf) {muoncorsf = 1;}
   applied->SetPtEtaPhiM(muoncorsf*tmp->Pt(),tmp->Eta(),tmp->Phi(),tmp->M());
   //cout << "tmp pt : " << tmp->Pt() << " applied pt : " << applied->Pt() << endl; 
   return applied;
}
void ssb_analysis::CorrectedMuonCollection()
{
   //cout << "Star CorrectedMuonCollec"<< endl;
   v_CorMuons.clear();
   TLorentzVector* tmp_;
   for (int i =0 ; i < Muon->GetEntriesFast(); ++i )
   {
      tmp_ = (TLorentzVector*)Muon->At(i);
      tmp_ = ApplyRocCor(tmp_,i);
      v_CorMuons.push_back(tmp_);
   }
}
void ssb_analysis::MakeElecCollection()
{
   //cout << "Star MakeElecCollection"<< endl;
   v_Elec.clear();
   TLorentzVector* tmp_;
   for (int i =0 ; i < Elec->GetEntriesFast(); ++i )
   {
      tmp_ = (TLorentzVector*)Elec->At(i);
      tmp_ = ApplyElecSCSM(tmp_,i);
      v_Elec.push_back(tmp_);
   }

}
TLorentzVector* ssb_analysis::ApplyElecSCSM( TLorentzVector* tmp,int index_)
{
   if (FileName_.Contains("Data")) {return tmp;}
   if (EleScSmSys == "Central")    {return tmp;}
   TLorentzVector* systl_ = new TLorentzVector();
   double scsmfac_ = 1.0;
   if (EleScSmSys == "UpUp"){scsmfac_ = Elec_ScSmUpUp->at(index_);}
   else if (EleScSmSys == "UpDown"){scsmfac_ = Elec_ScSmUpDown->at(index_);}
   else if (EleScSmSys == "DownUp"){scsmfac_ = Elec_ScSmDownUp->at(index_);}
   else if (EleScSmSys == "DownDown"){scsmfac_ = Elec_ScSmDownDown->at(index_);}
   else {cout << "Check Out Your EleScSmSys Option !!! " << EleScSmSys << endl;}
   //cout << "scsmfac_ : " << scsmfac_ << endl;
   systl_->SetPtEtaPhiM(scsmfac_*tmp->Pt(),tmp->Eta(),tmp->Phi(),tmp->M());
   return systl_;
}
// Function for Lepton Selection
void ssb_analysis::LeptonSelector()
{
   v_lepton_idx.clear();
   v_muon_idx.clear();
   v_electron_idx.clear();

   // dummy vector //
   v_lep_idx_temp.clear();
   v_muon_idx_temp.clear();
   v_electron_idx_temp.clear();

   /// Muon Correction ///
   CorrectedMuonCollection();
   /// Electron Class for Central & Syst.
   MakeElecCollection();
   // Case of Dimuon //
   if (TString(Decaymode).Contains("dimuon")) 
   { 
      Int_t nmu = Muon->GetEntriesFast();
      for (int i =0; i < nmu; ++i )
      {
         Mu_lep_sel = v_CorMuons.at(i);
         if ( ( v_lepton_iso->at(i) > lepisocut )   ||
              (*Mu_lep_sel).Pt()    < lep_pt        || 
              fabs( (*Mu_lep_sel).Eta() ) > lep_eta || 
              ( v_lepton_Id->at(i) == false )          ){continue;}         
         v_lep_idx_temp.push_back( i );

         if (v_lep_idx_temp.size()==1 && v_lepton_idx.size() == 0 ){  v_lepton_idx.push_back(i);}
         if (v_lep_idx_temp.size()>0)
         {
            if ( Muon_Charge->at( v_lep_idx_temp.at(0) ) != Muon_Charge->at( i ) ) { v_lepton_idx.push_back(i); }
         }
      }
   }
   // Case of Dielectron //
   else if (TString(Decaymode).Contains("dielec")) 
   { 
      Int_t nel = Elec->GetEntriesFast();

      for (int i =0; i < nel; ++i )
      {
         Ele_lep_sel = v_Elec.at(i);
         if ( (*Ele_lep_sel).Pt() != (*Ele_lep_sel).Pt() ){continue;}
         if ( (*Ele_lep_sel).Pt() < lep_pt )          {continue;}
         if ( fabs( (*Ele_lep_sel).Eta() ) > lep_eta ){continue;}   

         // Requirement of Electron Isolation with Electron Eta 
         //if (Elec_SCB_Medium->at(i) == false) {continue;}    
         //if (v_electron_Id->at(i) == false) {continue;}    
         if (v_lepton_Id->at(i) == false) {continue;}    
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.4442 && fabs( Elec_Supercluster_Eta->at(i) ) < 1.566 ) {continue;} 
         //if ( Elec_ChargeId_GsfCtfPx->at(i) == false ){continue;}//
         if ( Elec_ChargeId_GsfPx->at(i) == false ){continue;}//
         if ( fabs( Elec_Supercluster_Eta->at(i) ) <= 1.479 ) 
         { 
            if (v_lepton_iso->at(i) > 0.0695 ) {continue;}
            if ( fabs( Elec_Track_GsfdZ->at(i) ) > 0.10 || fabs( Elec_Track_GsfdXY->at(i) ) > 0.05 )  {continue;} 
         }
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.479 )  
         {  
            if (v_lepton_iso->at(i) > 0.0821 ) {continue;}
            if (fabs(Elec_Track_GsfdZ->at(i)) > 0.20 || fabs(Elec_Track_GsfdXY->at(i)) > 0.10 )  {continue;} 
         }

         v_lep_idx_temp.push_back( i );

         if ( v_lep_idx_temp.size()==1 && v_lepton_idx.size() == 0 ){  v_lepton_idx.push_back(i); }
         if ( v_lep_idx_temp.size() > 0 )
         {
            if ( Elec_Charge->at( v_lep_idx_temp.at(0) ) != Elec_Charge->at( i ) ) { v_lepton_idx.push_back(i); }
         }

      }
   }
   // Case of MuEl //
   else if (TString(Decaymode).Contains("muel"))
   {

      Int_t nmu = Muon->GetEntriesFast();
      for (int i =0; i < nmu; ++i )
      {
         Mu_lep_sel = v_CorMuons.at(i);
         if ( ( v_mlep_iso->at(i) > muonisocut ) || (*Mu_lep_sel).Pt() < muon_pt || fabs( (*Mu_lep_sel).Eta() ) > muon_eta || v_muon_Id->at(i) == false  ){continue;}         
         //if ( ( v_mlep_iso->at(i) > muonisocut ) || (*Mu_lep_sel).Pt() < 35 || fabs( (*Mu_lep_sel).Eta() ) > muon_eta || v_muon_Id->at(i) == false  ){continue;}         
         v_muon_idx.push_back( i );
      }

      Int_t nel = Elec->GetEntriesFast();
      for (int i =0; i < nel; ++i )
      {
         Ele_lep_sel = v_Elec.at(i);

         ////
         if ( (*Ele_lep_sel).Pt() < elec_pt ){continue;}
         //if ( (*Ele_lep_sel).Pt() < 35 ){continue;}

         if ( fabs( (*Ele_lep_sel).Eta() ) > elec_eta ){continue;}

         //if (Elec_SCB_Medium->at(i) == false) {continue;}    
         if (v_electron_Id->at(i) == false) {continue;}    
         if ( Elec_Conversion->at(i) == false ){continue;}
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.4442 && fabs( Elec_Supercluster_Eta->at(i) ) < 1.566 ) {continue;}
         if ( Elec_ChargeId_GsfPx->at(i) == false ){continue;}        
         if ( fabs( Elec_Supercluster_Eta->at(i) ) <= 1.479 ) 
         { 
            if ( fabs( Elec_Track_GsfdZ->at(i) ) > 0.10 || fabs( Elec_Track_GsfdXY->at(i) ) > 0.05 )  {continue;} 
         }
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.479 )  
         {  
            if (fabs(Elec_Track_GsfdZ->at(i)) > 0.20 || fabs(Elec_Track_GsfdXY->at(i)) > 0.10 )  {continue;} 
         }
         v_electron_idx_temp.push_back( i );

      }

      if (v_muon_idx.size() > 0 && v_electron_idx_temp.size() > 0 )
      {  
         for ( int i = 0; i < v_electron_idx_temp.size(); ++i )
         {
            if ( Muon_Charge->at( v_muon_idx.at(0) ) != Elec_Charge->at( v_electron_idx_temp.at(i) ) ) { v_electron_idx.push_back(v_electron_idx_temp.at(i) ); }
         }
      }
      
      if ( v_electron_idx.size() > 0 && v_muon_idx.size() > 0 ) 
      { 
         if ( Muon_Charge->at( v_muon_idx.at(0) ) == Elec_Charge->at( v_electron_idx.at(0) ) ) { cout << "MuEl seclection wrong !! " << endl; } 
      }

   }
   else { cout << "Lepton Selection error" << endl;}


   ///////////////////////////////////////////
   /// Electron Selection For Jet Cleaning /// 
   ///////////////////////////////////////////
   v_ele_idx_cand.clear();
   v_mu_idx_cand.clear();
   Int_t inele = Elec->GetEntriesFast();
   Int_t nmu   = Muon->GetEntriesFast();

   for (int ine = 0; ine <inele; ine++ )
   {  
      Ele_jet_clean = v_Elec.at(ine); 

      if (
           Elec_SCB_Veto->at(ine) == true &&
           Ele_jet_clean->Pt() > el_pt_jetcl &&
           fabs( Ele_jet_clean->Eta() ) < el_eta_jetcl &&
           Elec_Conversion->at(ine) == true
         )
      {
//         v_ele_idx_cand.push_back( ine );
         if ( fabs( Elec_Supercluster_Eta->at(ine) ) <= 1.479 )
         {
            if (
                 v_el_iso_jcl->at( ine ) < 0.175 &&
                 fabs( Elec_Track_GsfdXY->at(ine) ) < 0.05 &&
                 fabs( Elec_Track_GsfdZ->at(ine) ) < 0.10 
               )
            { v_ele_idx_cand.push_back( ine ); }

         } 
         else 
         {
            if (
                 v_el_iso_jcl->at( ine ) < 0.159 &&           
                 fabs( Elec_Track_GsfdXY->at(ine) ) < 0.10 &&
                 fabs( Elec_Track_GsfdXY->at(ine) ) < 0.20 
               )
            { v_ele_idx_cand.push_back( ine ); }
        
         }
      }

   }
     
   // Muon selectron for JetCleaning
   for ( int inmu =0; inmu < nmu; ++inmu )
   {
      Mu_jet_clean = v_CorMuons.at(inmu);
      if ( 
           (*Mu_jet_clean).Pt() > mu_pt_jetcl         && //Pt cut 
           fabs((*Mu_jet_clean).Eta()) < mu_eta_jetcl && //Eta cut
           v_mu_iso_jcl->at(inmu) < mu_iso_jetcl      && //Isolation cut...
           v_mu_Id_jcl->at(inmu) == true                 // Id cut ...
         )
      { v_mu_idx_cand.push_back( inmu ); }
   }
}
bool ssb_analysis::NumIsoLeptons()//YOU SHOULD REQUIRE THIS FUNCTION AFTER LEPTONSELETOR //
{
   bool numLeptons = true;
   if ( TString(Decaymode).Contains( "dimuon" ) ||
        TString(Decaymode).Contains( "dielec" )    ){  if ( v_lepton_idx.size() <= 1 ) {numLeptons=false;} } //Requiring  Two isolated Lepton .

   else if ( TString(Decaymode).Contains( "muel" ) ){ if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1 ) {numLeptons=false; } else {}}
   else if ( TString(Decaymode).Contains( "muonJet" ) ){if( v_lepton_idx.size() != 1){numLeptons=false;} }
   else { cout << "?? something wrong " << endl; }
   return numLeptons;
}
bool ssb_analysis::LeptonsPtAddtional()//YOU SHOULD REQUIRE THIS FUNCTION AFTER NumIsoLeptons //
{
   bool lepptadd = false;
   if ( TString(Decaymode).Contains( "dimuon" ) ||
        TString(Decaymode).Contains( "dielec" ) ||
        TString(Decaymode).Contains( "muel" )      ){  if ( Lep1->Pt() > 25 && Lep2->Pt() > 20 ) {lepptadd=true;} } //Pt of Leading Lepton should be over than 25 GeV and Second Leading Lepton Pt should be over thand 20 GeV.

   else if ( TString(Decaymode).Contains( "muonJet" ) ){if( v_lepton_idx.size() == 1){lepptadd=true;} }
   else { cout << "?? something wrong " << endl; }  
   return lepptadd;
}
bool ssb_analysis::MuVeto()
{
   bool muveto;
   muveto = false;
   int vetolep = 0;

   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      for ( int ilumuon = 0; ilumuon < Muon_Count; ++ilumuon )
      {
         if ( ilumuon != v_lepton_idx.at(0) && ilumuon != v_lepton_idx.at(1) )
         {
            if ( Muon_isLoose->at( ilumuon ) == true && Muon_PFIsodBeta04->at( ilumuon ) < 0.25 )
            {   
               vetoMu = v_CorMuons.at(ilumuon);
               if ( vetoMu->Pt() > 20 && fabs(vetoMu->Eta()) < 2.4 ) { vetolep++; }
            }
         }
      }

   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      // 3rd muon veto
      TLorentzVector *vetoMu = new TLorentzVector();
      for ( int ilumuon = 0; ilumuon < Muon_Count; ++ilumuon )
      {
         if ( Muon_isLoose->at( ilumuon ) == true && Muon_PFIsodBeta04->at( ilumuon ) < 0.25 )
         {   
            vetoMu = v_CorMuons.at(ilumuon);
            if ( vetoMu->Pt() > 20 && fabs(vetoMu->Eta()) < 2.4 ) { vetolep++; }
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      for ( int ilumuon = 0; ilumuon < Muon_Count; ++ilumuon )
      {        
         if ( ilumuon != v_muon_idx.at(0) )
         {             
            if ( Muon_isLoose->at( ilumuon ) == true && Muon_PFIsodBeta04->at( ilumuon ) < 0.25 )
            {   //vetolep++;
               vetoMu = v_CorMuons.at(ilumuon);
               if ( vetoMu->Pt() > 20 && fabs(vetoMu->Eta()) < 2.4 ) { vetolep++; }
            }  
         }        
      }
   }
   else { cout << "Something Wrong At MuVeto. Check Decaymode-" << endl;}

   if ( vetolep == 0 ){ muveto = true; }
   return muveto;
}
bool ssb_analysis::ElVeto()
{
   bool elveto;
   elveto = false;
   int vetolep = 0;
   if ( TString(Decaymode).Contains( "dimuon" ) )// case of dimuon ...
   {
      for ( int iele = 0; iele < Elec_Count; ++iele )
      {
         vetoEle = v_Elec.at(iele);
         if(  vetoEle->Pt() != vetoEle->Pt()) {continue;}
         if (
             Elec_SCB_Veto->at(iele) == true &&
             vetoEle->Pt() > 20                 &&
             fabs( vetoEle->Eta() ) < 2.5       &&
             Elec_Conversion->at(iele) == true     )
         {
          vetolep++; 
         }

      }
   } // dimuon requriement...
   else if( TString(Decaymode).Contains( "dielec" ) ) 
   {
      // 3rd Electron veto ..
      for ( int iele = 0; iele < Elec_Count; ++iele )
      {
         if ( iele != v_lepton_idx.at(0) && iele != v_lepton_idx.at(1) )
         {
            vetoEle = v_Elec.at(iele);
            if (
                 Elec_SCB_Veto->at(iele) == true &&
                 vetoEle->Pt()            > 20      &&
                 fabs( vetoEle->Eta() )   < 2.5     &&
                 Elec_Conversion->at(iele) == true     )
            { 
               if ( fabs( Elec_Supercluster_Eta->at(iele) ) <= 1.479 )
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.175 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.05 &&
                       fabs( Elec_Track_GsfdZ->at(iele) ) < 0.10 
                     ){ vetolep++; }
               
               } 
               else 
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.159 &&           
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.10 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.20 
                     ){ vetolep++; }
               }
            }
         }
      }

   } // case of dielectron...

   else if ( TString(Decaymode).Contains( "muel" ) ) 
   {
      for ( int iele = 0; iele < Elec_Count; ++iele )
      {
         if ( iele != v_electron_idx.at(0) )
         {
            vetoEle = v_Elec.at(iele);
            if(  vetoEle->Pt() != vetoEle->Pt()) {continue;}

            if (
                 Elec_SCB_Veto->at(iele) == true &&
                 vetoEle->Pt()            > 20   &&
                 fabs( vetoEle->Eta() )   < 2.5  &&
                 Elec_Conversion->at(iele) == true 
               )
            {
               if ( fabs( Elec_Supercluster_Eta->at(iele) ) <= 1.479 )
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.175 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.05 &&
                       fabs( Elec_Track_GsfdZ->at(iele) ) < 0.10 
                     )
                  { vetolep++; }
               
               } 
               else 
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.159 &&           
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.10 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.20 
                     )
                  { vetolep++; }
               }              
            }
         }
      }
   }
   else {cout << "EleVto check !!! " << endl;}
   if ( vetolep == 0 ){ elveto = true; }

   return elveto;
}
bool ssb_analysis::ThirdLeptonVeto()
{
   bool third_veto = true;
   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" )    )
   { 
      if ( v_lepton_idx.size() <= 1 ) {third_veto = false;;}
      //if ( v_lepton_idx.size() != 2 ) {third_veto = false;;}

      if( TString(Decaymode).Contains( "dimuon" ) ) 
      { 
         if ( Muon_Charge->at( v_lepton_idx.at(0) ) == Muon_Charge->at( v_lepton_idx.at(1) ) ) {third_veto = false;;}
         if ( MuVeto() == false ) {third_veto = false;;}
         if ( ElVeto() == false ) {third_veto = false;;}
      } 
      else if( TString(Decaymode).Contains( "dielec" ) ) 
      {
         if ( Elec_Charge->at( v_lepton_idx.at(0) ) == Elec_Charge->at( v_lepton_idx.at(1) ) ) {third_veto = false;;}
         if ( MuVeto() == false ) {third_veto = false;;}
         if ( ElVeto() == false ) {third_veto = false;;}
      }
   } //Requiring  Two isolated Lepton .

   else if ( TString(Decaymode).Contains( "muel" ) )
   { 
      if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1 ) { third_veto = false;; } 
      if ( Muon_Charge->at( v_muon_idx.at(0) ) == Elec_Charge->at( v_electron_idx.at(0) ) ) { third_veto = false;; }
      if ( MuVeto() == false ) {third_veto = false;;}
      if ( ElVeto() == false ) {third_veto = false;;}
   }

   else { cout << "something wrong at ThirdLeptonVeto veto Function !!!!" << endl; }
   return third_veto;
}
void ssb_analysis::LeptonOrder()
{
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      if (v_lepton_idx.size() > 1)
      {
         Lep1 = v_CorMuons.at( v_lepton_idx.at(0) );
         Lep2 = v_CorMuons.at( v_lepton_idx.at(1) );
         if (Muon_Charge->at( v_lepton_idx.at(0) ) < 0 )
         {
            Lep   = v_CorMuons.at( v_lepton_idx.at(0) ); 
            AnLep = v_CorMuons.at( v_lepton_idx.at(1) );
            // If muon(Lep1)'s charge in the first lepton idex is less then zero, this muon will be lepton that has minus charge
            // at the moment, another muon(Lep2) will be anti lepton that has plus charge
         }
         else
         {
            Lep   = v_CorMuons.at( v_lepton_idx.at(1) );
            AnLep = v_CorMuons.at( v_lepton_idx.at(0) );
         }
      }
   }
 
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      if (v_lepton_idx.size() > 1)
      {
         Lep1 = v_Elec.at( v_lepton_idx.at(0) );
         Lep2 = v_Elec.at( v_lepton_idx.at(1) );
         if (Elec_Charge->at( v_lepton_idx.at(0) ) < 0 )
         {
            Lep   = v_Elec.at( v_lepton_idx.at(0) );
            AnLep = v_Elec.at( v_lepton_idx.at(1) );
         }
         else
         {
            Lep   = v_Elec.at( v_lepton_idx.at(1) );
            AnLep = v_Elec.at( v_lepton_idx.at(0) );
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      if ( v_muon_idx.size() > 0 && v_electron_idx.size() > 0 )
      {
         TMuon = v_CorMuons.at( v_muon_idx.at(0) );
         TElectron = v_Elec.at( v_electron_idx.at(0) );
         if ((*TMuon).Pt() > (*TElectron).Pt())
         {
            Lep1 = v_CorMuons.at( v_muon_idx.at(0) );
            Lep2 = v_Elec.at( v_electron_idx.at(0) );
         }
         else
         {
            Lep1 = v_Elec.at( v_electron_idx.at(0) );
            Lep2 = v_CorMuons.at( v_muon_idx.at(0) );
         }
         if ( Muon_Charge->at( v_muon_idx.at(0) ) < 0 )
         {
            Lep   = v_CorMuons.at( v_muon_idx.at(0) );
            AnLep = v_Elec.at( v_electron_idx.at(0) );
         }
         else 
         {
            AnLep = v_CorMuons.at( v_muon_idx.at(0) );
            Lep   = v_Elec.at( v_electron_idx.at(0) );
         }
      }
   }
   else {cout << "Lepton TLorentzVector Error" << endl;}
}

void ssb_analysis::LeptonGetSF()
{

   lep_eff = 1.0;
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      lep_eff = SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,LepIdSFSys,LepIsoSFSys,LepTrackSFSys); 
   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      //lep_eff = SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,LepIdSFSys,LepIsoSFSys);
      lep_eff = SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,LepIdSFSys,LepRecoSFSys);// LepRecoSFSys is for Lepton //
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      //lep_eff = 1.0;
      lep_eff = SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,LepIdSFSys,LepIsoSFSys,LepTrackSFSys,LepRecoSFSys);//LepRecoSFSys is for Electron
   }
   else { lep_eff = 1.0; cout << "LeptonGetSF error !!!!" << endl; }
}
void ssb_analysis::LeptonSFApply()
{

   double lep_sf_cent = 0.0;

   evt_weight_beforeLepsf_ = 1.0;
   evt_weight_beforeLepsf_ = evt_weight_; // keep event weight 
   if ( !TString(FileName_).Contains( "Data") )
   {
      LeptonGetSF(); // Getting Lepton SF //
      evt_weight_ = evt_weight_*lep_eff;
   }
   else {evt_weight_ = 1;}

   /// For All-in-One Systematic ///
   if ( isAllSyst == true )
   {
      /// Get Central Lepton SF for each channel
      if ( !TString(FileName_).Contains( "Data") )
      {
         if ( TString(Decaymode).Contains( "dimuon" ) )
         {
            lep_sf_cent = SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","central","central"); // lep1 lep2 idsys isosys tracksys 
         }
         else if ( TString(Decaymode).Contains( "dielec" ) )
         {
            lep_sf_cent = SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"central","central");// LepRecoSFSys is for Lepton //
         }
         else if ( TString(Decaymode).Contains( "muel" ) )
         {
            lep_sf_cent = SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"central","central","central","central");//LepID, LepIso, Tracksys, ReconSys
         }
         else {lep_sf_cent = 0.0;cout << "Error In LeptonSFApply" << endl;}
      }
      for( int i =0; i < v_SystFullName.size(); ++i )
      {
         if ( !TString(FileName_).Contains( "Data") ){
            if ( v_SystFullName[i].Contains("LepIDUp"))
            {
               if ( TString(Decaymode).Contains( "dimuon" ) ) 
               {
                  //v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"up","central"); // lep1 lep2 idsys isosys 
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"up","central","central"); // lep1 lep2 idsys isosys tracksys 
               }
               else if ( TString(Decaymode).Contains( "dielec" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"up","central");// LepRecoSFSys is for Lepton //
               }
               else if ( TString(Decaymode).Contains( "muel" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"up","central","central","central");//LepID, LepIso, Tracksys, ReconSys
               }
               else {cout << "LeptonSFApply error for All-in-One Syst.!!!!" << endl;}
            }
            else if ( v_SystFullName[i].Contains("LepIDDown"))
            {
               if ( TString(Decaymode).Contains( "dimuon" ) ) 
               {
                  //v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"down","central"); // lep1 lep2 idsys isosys 
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"down","central","central"); // lep1 lep2 idsys isosys tracksys 
               }
               else if ( TString(Decaymode).Contains( "dielec" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"down","central");// LepRecoSFSys is for Lepton //
               }
               else if ( TString(Decaymode).Contains( "muel" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"down","central","central","central");//LepID, LepIso, Tracksys, ReconSys
               }
               else {cout << "LeptonSFApply error for All-in-One Syst.!!!!" << endl;}

            }
            else if ( v_SystFullName[i].Contains("LepIsoUp"))
            {
               if ( TString(Decaymode).Contains( "dimuon" ) ) 
               {
                  //v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","up"); // lep1 lep2 idsys isosys 
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","up","central"); // lep1 lep2 idsys isosys tracksys
               }
               else if ( TString(Decaymode).Contains( "dielec" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"central","central");// LepRecoSFSys is for Lepton //
               }
               else if ( TString(Decaymode).Contains( "muel" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"central","up","central","central");//LepID, LepIso, Tracksys, ReconSys
               }
               else {cout << "LeptonSFApply error for All-in-One Syst.!!!!" << endl;}
            }
            else if ( v_SystFullName[i].Contains("LepIsoDown"))
            {
               if ( TString(Decaymode).Contains( "dimuon" ) ) 
               {
                  //v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","down"); // lep1 lep2 idsys isosys 
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","down","central"); // lep1 lep2 idsys isosys tracksys
               }
               else if ( TString(Decaymode).Contains( "dielec" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"central","central");// LepRecoSFSys is for Lepton //
               }
               else if ( TString(Decaymode).Contains( "muel" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"central","down","central","central");//LepID, LepIso, Tracksys, ReconSys
               }
               else {cout << "LeptonSFApply error for All-in-One Syst.!!!!" << endl;}

            }

            else {
               /// At here, you will use leptons sf as central value ///
               v_SystEvt[i] = v_SystEvt[i]*lep_sf_cent;
            }
         }
         else { // Data 
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i]; 
//         cout << "v_SystFullName[" << i << "] : " << v_SystFullName[i] << " : " << v_SystEvt[i] << endl;
      }
   }
}
// Jet Selection and Jet cleaning 
void ssb_analysis::JetSelector()
{
   v_jet_idx.clear();
   v_jetup_idx.clear();
   v_jetdn_idx.clear();
   v_jetresup_idx.clear();
   v_jetresdn_idx.clear();

   v_jet_TL.clear();
   v_jetup_TL.clear();
   v_jetdn_TL.clear();
   v_jetresup_TL.clear();
   v_jetresdn_TL.clear();


   /// For JER 
   v_jetdpt_res.clear();    // dpt of JER for All-in-One syst.
   v_jetdpt_jesup.clear();  // dpt of JERUp for All-in-One syst.
   v_jetdpt_jesdn.clear();  // dpt of JERDn for All-in-One syst.
   v_jetdpt_resup.clear();  // dpt of JERUp for All-in-One syst.
   v_jetdpt_resdn.clear();  // dpt of JERDn for All-in-One syst.

   // Delta px JER Variation // 
   v_jetdpx_res.clear();    // dpx of JER for All-in-One syst.
   v_jetdpx_jesup.clear();  // dpx of JERUp for All-in-One syst.
   v_jetdpx_jesdn.clear();  // dpx of JERDn for All-in-One syst.
   v_jetdpx_resup.clear();  // dpx of JERUp for All-in-One syst.
   v_jetdpx_resdn.clear();  // dpx of JERDn for All-in-One syst.

   // Delta py JER Variation // 
   v_jetdpy_res.clear();    // dpy of JER for All-in-One syst.
   v_jetdpy_jesup.clear();  // dpy of JERUp for All-in-One syst.
   v_jetdpy_jesdn.clear();  // dpy of JERDn for All-in-One syst.
   v_jetdpy_resup.clear();  // dpy of JERUp for All-in-One syst.
   v_jetdpy_resdn.clear();  // dpy of JERDn for All-in-One syst.


   Int_t injet = Jet->GetEntriesFast();

   TJet = new TLorentzVector();

   TLorentzVector *JetUp =new TLorentzVector();
   TLorentzVector *JetDown = new TLorentzVector();

   TLorentzVector *JetCen = new TLorentzVector();
   TLorentzVector *JetUpPre = new TLorentzVector();
   TLorentzVector *JetDownPre = new TLorentzVector();
   TLorentzVector* TJetResUp;// = new TLorentzVector();
   TLorentzVector* TJetResDn;// = new TLorentzVector();

   TLorentzVector *Jetpre = new TLorentzVector();

   double jetpt_  = -999; 
   double jeteta_ = -999; 
   double jetphi_ = -999; 
   double jetenergy_ = -999; 
   for (int i =0; i < injet; i++)
   {

      TJet = new TLorentzVector();
      Jetpre = new TLorentzVector();
      Jetpre = (TLorentzVector*)Jet->At(i);
      //////////////////////////////
      /// JES Systematic Method1 ///
      //////////////////////////////
      if      ( TString(JetEnSys).Contains("JetEnNorm"     ) ) { Jetpre = (TLorentzVector*)Jet->At(i); }
      else if ( TString(JetEnSys).Contains("JetEnShiftedUp") )
      {
         jetpt_     = Jetpre->Pt();
         jeteta_    = Jetpre->Eta();
         jetphi_    = Jetpre->Phi();
         jetenergy_ = Jetpre->Energy();
         /// Apply TJet
         Jetpre->SetPtEtaPhiE(jetpt_*Jet_EnShiftedUp->at(i),jeteta_,jetphi_,jetenergy_*Jet_EnShiftedUp->at(i));
      }
      else if ( TString(JetEnSys).Contains("JetEnShiftedDown") )
      {
         jetpt_     = Jetpre->Pt();
         jeteta_    = Jetpre->Eta();
         jetphi_    = Jetpre->Phi();
         jetenergy_ = Jetpre->Energy();

         Jetpre->SetPtEtaPhiE(jetpt_*Jet_EnShiftedDown->at(i),jeteta_,jetphi_,jetenergy_*Jet_EnShiftedDown->at(i));
      }
      else { cout << "Check out Jet Systematic Config ... Default is Norminal Jet ..."<< endl; Jetpre = (TLorentzVector*)Jet->At(i); }

      //////////////////////////////
      /// JER Systematic Method1 ///
      //////////////////////////////
      if      ( TString(JetResSys).Contains("JetResNorm"     ) ) { TJet = JERSmearing(Jetpre,i,"Cen"); }
      else if ( TString(JetResSys).Contains("JetResShiftedUp") )
      {
         TJet = JERSmearing(Jetpre,i,"Up");
      }
      else if ( TString(JetResSys).Contains("JetResShiftedDown") )
      {
         TJet = JERSmearing(Jetpre,i,"Down");
      }
      else { cout << "Check out Jet Systematic Config ... Default is Norminal Jet ..."<< endl; TJet = JERSmearing(Jetpre,i,"Cen"); }
    
/*      /// Jet Pt Resolution Dilution Study ///
      if (JetPtResDil != "None" &&  JetPtResDil != "none" ) {
         /// You should check your systematic options, when you study Dilution... 
         //  These part won't be correct, if you use other systematic options... 
         //  So we have to check your systematic options and dilution study options //
         TJet = ApplyJetPtRes(Jetpre,i,JetPtResDil);
      }*/
      //cout << "After Pt Res Jet Pt : " << TJet->Pt() << endl;    
      if ( JetCleaning(TJet) == true ) {
        /* if (JetPhiResDil != "None" &&  JetPhiResDil != "none" ) 
         {
         //   cout << "JetPhiResDil : " << JetPhiResDil << endl;
            TJet = ApplyJetPhiRes(TJet,i,JetPhiResDil);
         }*/
         //cout << "After Phi Res Jet Pt : " << TJet->Pt() << endl;
         v_jetdpt_res.push_back(TJet->Pt()-Jetpre->Pt());  
         v_jetdpx_res.push_back(TJet->Px()-Jetpre->Px());  
         v_jetdpy_res.push_back(TJet->Py()-Jetpre->Py());  
         ////////////////////////////////
         /// Apply Requirement of Jet ///
         ////////////////////////////////
         if ( (*TJet).Pt() > jet_pt && fabs( (*TJet).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            //jet_pt: 30, jet_eta; 2.4, jet_id: PFLoose in analysis_config.config
            v_jet_idx.push_back( i ); 
            v_jet_TL.push_back( TJet ); 
            if ( (*TJet).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }
      ////////////////////////////////////////
      /// For JES Systematic with Method 2 ///
      ////////////////////////////////////////

      JetUp =new TLorentzVector();
      JetDown = new TLorentzVector();

      JetCen = new TLorentzVector();
      JetUpPre = new TLorentzVector();
      JetDownPre = new TLorentzVector();
      JetCen = (TLorentzVector*)Jet->At(i); 

      JetUpPre->SetPtEtaPhiE(JetCen->Pt()*Jet_EnShiftedUp->at(i),JetCen->Eta(),JetCen->Phi(),JetCen->Energy()*Jet_EnShiftedUp->at(i));
      JetDownPre->SetPtEtaPhiE(JetCen->Pt()*Jet_EnShiftedDown->at(i),JetCen->Eta(),JetCen->Phi(),JetCen->Energy()*Jet_EnShiftedDown->at(i));
            
      JetUp = JERSmearing(JetUpPre,i,"Cen");
      if (JetCleaning(JetUp) == true ) { 
         v_jetdpt_jesup.push_back(JetUp->Pt()-JetUpPre->Pt());  
         v_jetdpx_jesup.push_back(JetUp->Px()-JetUpPre->Px());  
         v_jetdpy_jesup.push_back(JetUp->Py()-JetUpPre->Py());  
         if ( (*JetUp).Pt() > jet_pt && fabs( (*JetUp).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            v_jetup_idx.push_back( i ); 
            v_jetup_TL.push_back( JetUp ); 
            if ( (*JetUp).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }     
      JetDown = JERSmearing(JetDownPre,i,"Cen");
      if (JetCleaning(JetDown) == true ) { 
         v_jetdpt_jesdn.push_back(JetDown->Pt()-JetDownPre->Pt());  
         v_jetdpx_jesdn.push_back(JetDown->Px()-JetDownPre->Px());  
         v_jetdpy_jesdn.push_back(JetDown->Py()-JetDownPre->Py());  
         if ( (*JetDown).Pt() > jet_pt && fabs( (*JetDown).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            v_jetdn_idx.push_back( i ); 
            v_jetdn_TL.push_back( JetDown ); 
            if ( (*JetDown).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }

      ////////////////////////////////////////
      /// For JER Systematic with Method 2 ///
      ////////////////////////////////////////

      TJetResUp = new TLorentzVector();
      TJetResDn = new TLorentzVector();

      TJetResUp = JERSmearing(Jetpre,i,"Up");
      if ( JetCleaning(TJetResUp) == true )
      {
         v_jetdpt_resup.push_back(TJetResUp->Pt()-Jetpre->Pt());  
         v_jetdpx_resup.push_back(TJetResUp->Px()-Jetpre->Px());  
         v_jetdpy_resup.push_back(TJetResUp->Py()-Jetpre->Py());  
         
         if ( (*TJetResUp).Pt() > jet_pt && fabs( (*TJetResUp).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            v_jetresup_idx.push_back( i ); 
            v_jetresup_TL.push_back( TJetResUp ); 
            if ( (*TJetResUp).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }

      TJetResDn = JERSmearing(Jetpre,i,"Down");
      if ( JetCleaning(TJetResDn) == true )
      {
         v_jetdpt_resdn.push_back(TJetResDn->Pt()-Jetpre->Pt());  
         v_jetdpx_resdn.push_back(TJetResDn->Px()-Jetpre->Px());  
         v_jetdpy_resdn.push_back(TJetResDn->Py()-Jetpre->Py());  
         
         if ( (*TJetResDn).Pt() > jet_pt && fabs( (*TJetResDn).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
         {
            v_jetresdn_idx.push_back( i ); 
            v_jetresdn_TL.push_back( TJetResDn ); 
            if ( (*TJetResDn).Pt() < 30. ){cout << "Something wrong " << endl;}
         }
      }
   }

//   JetCleaning(v_jet_idx_cand, v_jet_TL_cand, v_jet_idx, v_jet_TL);
   if ( TString(FileName_).Contains( "Data" ) )
   {
      v_jetdn_idx.clear();
      v_jetup_idx.clear();
      v_jetresdn_idx.clear();
      v_jetresup_idx.clear();

      v_jetdn_TL.clear();
      v_jetup_TL.clear();
      v_jetresdn_TL.clear();
      v_jetresup_TL.clear();

      v_jetdn_idx = v_jet_idx;
      v_jetup_idx = v_jet_idx;
      v_jetresdn_idx = v_jet_idx;
      v_jetresup_idx = v_jet_idx;

      v_jetdn_TL = v_jet_TL;
      v_jetup_TL = v_jet_TL;
      v_jetresdn_TL = v_jet_TL;
      v_jetresup_TL = v_jet_TL;
   }
}
void ssb_analysis::JetCleaning(std::vector<int> v_jidx_cand, std::vector<TLorentzVector*>v_jtl_cand, std::vector<int>& v_jidx, std::vector<TLorentzVector*>&v_jtl )
{
   bool overlapjetlep = false;
   v_jidx.clear();
   v_jtl.clear();
   
   TLorentzVector* TmpJet = new TLorentzVector();
   for ( int i =0; i < v_jidx_cand.size(); ++i ) // Jet Loop (before jet cleanning )
   {
      l_jet = 0;
      overlapjetlep = false;
      TmpJet = new TLorentzVector();
      TmpJet = v_jtl_cand.at(i);
      // electron
      if ( v_ele_idx_cand.size() > 0 )
      {
         for ( int j = 0; j < v_ele_idx_cand.size(); ++j )
         {
            selected_Ele_jet_clean = new TLorentzVector();
            selected_Ele_jet_clean = v_Elec.at( v_ele_idx_cand.at(j) );
            if ( (*selected_Ele_jet_clean).DeltaR( (*TmpJet) ) < 0.4){ overlapjetlep = true; }
         }
      }

      // muon
      if ( v_mu_idx_cand.size() > 0 )
      {  
         for ( int k = 0; k < v_mu_idx_cand.size(); ++k )
         {
            selected_Mu_jet_clean = new TLorentzVector();
            selected_Mu_jet_clean = v_CorMuons.at( v_mu_idx_cand.at(k) );
            if ( (*selected_Mu_jet_clean).DeltaR( (*TmpJet) ) < 0.4){ overlapjetlep = true; }
         }
      }
      if ( overlapjetlep != true  ) 
      {
         v_jidx.push_back( v_jidx_cand.at(i) );
         v_jtl.push_back( v_jtl_cand.at(i) );
      }
   }
}
bool ssb_analysis::JetCleaning(TLorentzVector* jet_)
{
   bool overlapjetlep = false;
   // electron
   for ( int j = 0; j < v_ele_idx_cand.size(); ++j )
   {
      selected_Ele_jet_clean = new TLorentzVector();
      selected_Ele_jet_clean = v_Elec.at( v_ele_idx_cand.at(j) );
      if ( (*selected_Ele_jet_clean).DeltaR( (*jet_) ) < 0.4){ overlapjetlep = true; }
   }
   // muon
   for ( int k = 0; k < v_mu_idx_cand.size(); ++k )
   {
      selected_Mu_jet_clean = new TLorentzVector();
      selected_Mu_jet_clean = v_CorMuons.at( v_mu_idx_cand.at(k) );
      if ( (*selected_Mu_jet_clean).DeltaR( (*jet_) ) < 0.4){ overlapjetlep = true; }
   }
   return !overlapjetlep; 
}
void ssb_analysis::JetDefiner()
{
   Jet1 = new TLorentzVector();
   Jet2 = new TLorentzVector();

   if (v_jet_TL.size() >=1)
   {
      Jet1 = v_jet_TL[0];
      if (v_jet_TL.size() > 1){ Jet2 = v_jet_TL[1]; }
   }
   if(isAllSyst == true){ 
      Jet1Up = new TLorentzVector();
      Jet2Up = new TLorentzVector();
      Jet1Dn = new TLorentzVector();
      Jet2Dn = new TLorentzVector();
      Jet1JERUp = new TLorentzVector();
      Jet2JERUp = new TLorentzVector();
      Jet1JERDn = new TLorentzVector();
      Jet2JERDn = new TLorentzVector();

      if (v_jetup_TL.size() >= 1) {
         Jet1Up = v_jetup_TL[0];
         if (v_jetup_TL.size() > 1) { Jet2Up = v_jetup_TL[1]; }
      }
      if ( v_jetdn_TL.size() >= 1 ) { 
         Jet1Dn = v_jetdn_TL[0]; 
         if ( v_jetdn_TL.size() > 1 ) Jet2Dn = v_jetdn_TL[1]; 
      }
      if (v_jetresup_TL.size() >= 1) {
         Jet1JERUp = v_jetresup_TL[0];
         if (v_jetresup_TL.size() > 1) { Jet2JERUp = v_jetresup_TL[1]; }
      }
      if ( v_jetresdn_TL.size() >= 1 ) { 
         Jet1JERDn = v_jetresdn_TL[0]; 
         if ( v_jetresdn_TL.size() > 1 ) Jet2JERDn = v_jetresdn_TL[1]; 
      }
   }
}
void ssb_analysis::METDefiner()
{
 
   Met = new TLorentzVector();
   MetJESUp = new TLorentzVector();
   MetJESDn = new TLorentzVector();
   MetJERUp = new TLorentzVector();
   MetJERDn = new TLorentzVector();
   /// For Data, We don't need to apply JES systematic ///   
   if (TString(FileName_).Contains( "Data" ))
   {
      Met =  (TLorentzVector*)METMUEGCleanCor->At(0); 
      MetJESDn = Met; MetJESUp = Met;
      MetJERDn = Met; MetJERUp = Met;
   }
   else {
      Met =  (TLorentzVector*)METMUCleanCor->At(0);
      if ( TString(JetEnSys).Contains("JetEnNorm"     ) )
      {
         Met = Met;
      }
      else if ( TString(JetEnSys).Contains("JetEnShiftedUp") )
      {
         Met->SetPtEtaPhiE(METMUCleanCor_JetEnShiftedUp_PT->at(0) ,0, METMUCleanCor_JetEnShiftedUp_Phi->at(0),0);
      }
      else if ( TString(JetEnSys).Contains("JetEnShiftedDown") )
      {
         Met->SetPtEtaPhiE(METMUCleanCor_JetEnShiftedDown_PT->at(0) ,0, METMUCleanCor_JetEnShiftedDown_Phi->at(0),0);
      }
      else if ( TString(MetSys).Contains("None") ){
         Met = Met;
      }
      else if ( TString(MetSys).Contains("Up") )
      { cout << "Case Up" << endl;
         Met->SetPtEtaPhiE(METMUCleanCor_UnclusteredEnShiftedUp_PT->at(0) ,0, METMUCleanCor_UnclusteredEnShiftedUp_Phi->at(0),0);
      }
      else if ( TString(MetSys).Contains("Down") )
      {  cout << "Case Down" << endl;
         Met->SetPtEtaPhiE(METMUCleanCor_UnclusteredEnShiftedDown_PT->at(0) ,0, METMUCleanCor_UnclusteredEnShiftedDown_Phi->at(0),0);
      }
      else 
      {
         Met = Met;
         cout << "Check out Jet Systematic Config ... Default is Norminal MET ..."<< endl;      
      }
      Met = METSmear( v_jetdpt_res, v_jetdpx_res, v_jetdpy_res, Met );
      // apply MET Smear //
      if (isAllSyst == true)
      {
         MetJESUp->SetPtEtaPhiE(METMUCleanCor_JetEnShiftedUp_PT->at(0) ,0, METMUCleanCor_JetEnShiftedUp_Phi->at(0),0);
         MetJESDn->SetPtEtaPhiE(METMUCleanCor_JetEnShiftedDown_PT->at(0) ,0, METMUCleanCor_JetEnShiftedDown_Phi->at(0),0);
         
         MetJERUp = (TLorentzVector*)METMUCleanCor->At(0);
         MetJERDn = (TLorentzVector*)METMUCleanCor->At(0);
         
         MetJESUp = METSmear( v_jetdpt_jesup, v_jetdpx_jesup, v_jetdpy_jesup, MetJESUp );
         MetJESDn = METSmear( v_jetdpt_jesdn, v_jetdpx_jesdn, v_jetdpy_jesdn, MetJESDn );
         /// Jer MET Smear ///
         MetJERUp = METSmear( v_jetdpt_resup, v_jetdpx_resup, v_jetdpy_resup, MetJERUp );
         MetJERDn = METSmear( v_jetdpt_resdn, v_jetdpx_resdn, v_jetdpy_resdn, MetJERDn );
      }
   }
}

TLorentzVector* ssb_analysis::METSmear(std::vector<double> v_dpt, std::vector<double> v_dpx, std::vector<double> v_dpy, TLorentzVector* met_ )
{
   if ( dojer&!TString(FileName_).Contains( "Data" ) ){ 
      TLorentzVector* smearmet_ = new TLorentzVector();
      double sum_dpt =0.;
      double sum_dpx =0.;
      double sum_dpy =0.;
      // check dpt
      for (int i = 0; i< v_dpt.size(); ++i )
      {   
         sum_dpt += v_dpt[i];
         sum_dpx += v_dpx[i];
         sum_dpy += v_dpy[i];
      }
  
      smearmet_->SetPxPyPzE(met_->Px()-sum_dpx,met_->Py()-sum_dpy,0,0);
      return smearmet_;
   }
   else { return met_;}
}
// DiMuon Mass Cut : step 1
bool ssb_analysis::DiLeptonMassCut()
{
   dimu_masscut = false;

   if ( ((*Lep1)+(*Lep2)).M() > 20 ){ dimu_masscut = true; }

   return dimu_masscut;
}

// ZVeto Cut : step 2
bool ssb_analysis::ZVetoCut()
{
   zvetocut = false;

   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) ) 
   {
      if ( ((*Lep1)+(*Lep2)).M() <= 76 || ((*Lep1)+(*Lep2)).M() >= 106 ){ zvetocut = true; }
   }
   else if ( TString(Decaymode).Contains( "muel" ) ){ zvetocut = true; }
   else {cout << "ZVeto Error !!" << endl;}

   return zvetocut;
}

// Num.Jet Cut : Step 3
bool ssb_analysis::NumJetCut(std::vector<int> v_jets)
{
   numjetcut = false;
   //if ( v_jets.size() == 2 ){ numjetcut = true; }
   if ( v_jets.size() >= 2 ){ numjetcut = true; }
   return numjetcut;
}

//MET Cut : step 4
bool ssb_analysis::METCut(TLorentzVector* met)
{
   metcut = false;
   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) )
   {
      if (met->Pt() > 40) { metcut =true; }
   }
   else if ( TString(Decaymode).Contains( "muel" ) ){ metcut =true; }
   else {cout << "METCut Error !!" << endl;}
   return metcut;
}
void ssb_analysis::BDsicApply()
{
   v_bjet_idx.clear();
   v_bjetup_idx.clear();
   v_bjetdn_idx.clear();
   v_bjetresup_idx.clear();
   v_bjetresdn_idx.clear();
   v_bjet_TL.clear();
   v_bjetup_TL.clear();
   v_bjetdn_TL.clear();
   v_bjetresup_TL.clear();
   v_bjetresdn_TL.clear();
   nbtagged =0;
   for (int ije = 0; ije < v_jet_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jet_idx.at(ije) ) > bdisccut ) { nbtagged++;v_bjet_idx.push_back( v_jet_idx.at(ije) ); v_bjet_TL.push_back( v_jet_TL.at(ije) ); }
   }
   /// For JES Up 
   for (int ije = 0; ije < v_jetup_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetup_idx.at(ije) ) > bdisccut ) { v_bjetup_idx.push_back( v_jetup_idx.at(ije) ); v_bjetup_TL.push_back( v_jetup_TL.at(ije) ); }
   }
   /// For JES Down 
   for (int ije = 0; ije < v_jetdn_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetdn_idx.at(ije) ) > bdisccut ) { v_bjetdn_idx.push_back( v_jetdn_idx.at(ije) ); v_bjetdn_TL.push_back( v_jetdn_TL.at(ije) ); }
   }
   /// For JER Up 
   for (int ije = 0; ije < v_jetresup_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetresup_idx.at(ije) ) > bdisccut ) { v_bjetresup_idx.push_back( v_jetresup_idx.at(ije) ); v_bjetresup_TL.push_back( v_jetresup_TL.at(ije) ); }
   }
   /// For JER Down 
   for (int ije = 0; ije < v_jetresdn_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetresdn_idx.at(ije) ) > bdisccut ) { v_bjetresdn_idx.push_back( v_jetresdn_idx.at(ije) ); v_bjetresdn_TL.push_back( v_jetresdn_TL.at(ije) ); }
   }

}
//B Jet Cut : step 5
bool ssb_analysis::BJetCut( std::vector<int>v_bjets )
{
   bjetcut = false;
   if (v_bjets.size() >= 1) { bjetcut = true; }
   return bjetcut;
}
void ssb_analysis::BTaggigSFApply()
{
   evt_weight_beforeBtag_ = 1;
   evt_weight_beforeBtag_ = evt_weight_;
   if( TString(FileName_).Contains( "Data" ) ) {return;}
   std::vector<double> v_Jet_pT;
   std::vector<double> v_Jet_eta;
   std::vector<double> v_Jet_bDisc;
   std::vector<int> v_Jet_flav;
   v_Jet_pT.clear();
   v_Jet_eta.clear();
   v_Jet_bDisc.clear();
   v_Jet_flav.clear();
   TLorentzVector* jet_ = new TLorentzVector();
   for ( int ijet = 0; ijet < v_jet_idx.size(); ++ijet )
   {
      int idx_ = v_jet_idx[ijet];
      jet_ = (TLorentzVector*)v_jet_TL[ijet];
      v_Jet_pT.push_back(jet_->Pt()); 
      v_Jet_eta.push_back(jet_->Eta()); 
      v_Jet_bDisc.push_back(Jet_bDisc->at(idx_)); 
      v_Jet_flav.push_back(abs(Jet_HadronFlavour->at(idx_))); 
   } 

   if( !TString(FileName_).Contains( "Data" ) ) {evt_weight_ = SSBEffcal->Btagging_EvenWeight( v_Jet_pT, v_Jet_eta, v_Jet_bDisc, bdisccut,v_Jet_flav )*evt_weight_;}
   else { evt_weight_ = 1; }
   if ( isAllSyst == true )
   {
      double b_sf_cent = 1.0;
      //b_sf_cent = SSBEffcal->Btagging_EvenWeightSys(v_Jet_pT,v_Jet_eta,v_Jet_bDisc,bdisccut,"Central");
      b_sf_cent = SSBEffcal->Btagging_EvenWeightSys( v_Jet_pT, v_Jet_eta, v_Jet_bDisc, bdisccut,v_Jet_flav,"Central");
      for( int i =0; i < v_SystFullName.size(); ++i )
      {
         if ( !TString(FileName_).Contains( "Data") ){
            if ( v_SystFullName[i].Contains("Jet")  ){continue;}
            else if ( v_SystFullName[i].Contains("BTagSF") || v_SystFullName[i].Contains("BTagEff") )
            {
               //v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->Btagging_EvenWeightSys(v_Jet_pT, v_Jet_eta, v_Jet_bDisc, bdisccut,v_SystFullName[i]);
               v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->Btagging_EvenWeightSys( v_Jet_pT, v_Jet_eta, v_Jet_bDisc, bdisccut,v_Jet_flav,v_SystFullName[i]);
            }
            else {
               v_SystEvt[i] = v_SystEvt[i]*b_sf_cent;
            }
         }
         else { // Data 
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i]; 
      }
   }
}
void ssb_analysis::BTaggigSFApplyJESR(std::vector<int>v_jidx, std::vector<TLorentzVector*>v_jtl, int sysidx_)
{
   if( TString(FileName_).Contains( "Data" ) ) {return;}
   std::vector<double> v_Jet_pT;
   std::vector<double> v_Jet_eta;
   std::vector<double> v_Jet_bDisc;
   std::vector<int> v_Jet_flav;
   v_Jet_pT.clear();
   v_Jet_eta.clear();
   v_Jet_bDisc.clear();
   v_Jet_flav.clear();
   for ( int ij = 0; ij < v_jidx.size(); ++ij )
   {
      int idx_ = v_jidx[ij];
      TLorentzVector* jidx_ = v_jtl[ij];
      v_Jet_pT.push_back(jidx_->Pt()); 
      v_Jet_eta.push_back(jidx_->Eta()); 
      v_Jet_bDisc.push_back(Jet_bDisc->at(idx_)); 
      v_Jet_flav.push_back(abs(Jet_HadronFlavour->at(idx_))); 
   }
   if ( isAllSyst == true )
   {
      v_SystEvt[sysidx_] = v_SystEvt[sysidx_]*SSBEffcal->Btagging_EvenWeightSys(v_Jet_pT,v_Jet_eta,v_Jet_bDisc,bdisccut,v_Jet_flav,"Central"); 
      m_Syst_EvtW[ v_SystFullName[sysidx_] ] = v_SystEvt[sysidx_]; 
   }
}
void ssb_analysis::BTaggigSFApplyDilu(){
   evt_weight_beforeBtag_ = 1;
   evt_weight_beforeBtag_ = evt_weight_;
   if( TString(FileName_).Contains( "Data" ) ) {return;}
   std::vector<double> v_Jet_pT;
   std::vector<double> v_Jet_eta;
   std::vector<double> v_Jet_bDisc;
   std::vector<int> v_Jet_flav;
   v_Jet_pT.clear();
   v_Jet_eta.clear();
   v_Jet_bDisc.clear();
   v_Jet_flav.clear();
   TLorentzVector* jet_ = new TLorentzVector();
   for ( int ijet = 0; ijet < v_jet_idx.size(); ++ijet )
   {
      int idx_ = v_jet_idx[ijet];
      jet_ = (TLorentzVector*)v_jet_TL[ijet];
      v_Jet_pT.push_back(jet_->Pt()); 
      v_Jet_eta.push_back(jet_->Eta()); 
      v_Jet_bDisc.push_back(Jet_bDisc->at(idx_)); 
      v_Jet_flav.push_back(abs(Jet_HadronFlavour->at(idx_))); 
   } 

   if( !TString(FileName_).Contains( "Data" ) ) {evt_weight_ = SSBEffcal->Btagging_EvenWeight( v_Jet_pT, v_Jet_eta, v_Jet_bDisc, bdisccut,v_Jet_flav )*evt_weight_;}

}
void ssb_analysis::BJetDefiner(std::vector<int> v_jets, std::vector<int> v_bjets)
{
   //////////////////
   // Define BJets //
   //////////////////
   bcandi = -9;
   if( v_bjet_TL.size() > 1)
   {
      bJet1 = v_bjet_TL.at(0); 
      bJet2 = v_bjet_TL.at(1);
   }

   else if (  v_bjet_idx.size() > 0 )
   {
      bJet1 = v_bjet_TL.at(0);
      for ( int k = 0; k < v_jet_idx.size(); ++k)
      {
         if ( v_jet_idx.at(k) != v_bjet_idx.at(0) ){ bJet2 =v_jet_TL[k]; break; }
      }
   }
 
   // Bjet TLorentzVector
   else if ( v_bjet_idx.size() == 0 )
   {
      if ( v_jet_idx.size() > 1 )
      {
         bJet1 = v_jet_TL[0];
         bJet2 = v_jet_TL[1];
      }
   }
   else { cout << "you should check out nbjet cut!! " << endl; }
}
void ssb_analysis::BJetDefiner()
{
   //////////////////
   // Define BJets //
   //////////////////
   /// For Central Jet ///
   if( v_bjet_TL.size() > 1)
   {
      bJet1 = v_bjet_TL.at(0); 
      bJet2 = v_bjet_TL.at(1);
   }
   else if ( v_bjet_TL.size() > 0 )
   {  
      bJet1 = v_bjet_TL.at(0);
      for ( int i = 0; i < v_jet_idx.size(); ++i )
      {
         if ( v_jet_idx[i] != v_bjet_idx[0] ){ bJet2 = v_jet_TL[i]; break; }
      }
   }
   else if ( v_bjet_TL.size() == 0 )
   {
      bJet1 = Jet1;
      bJet2 = Jet2;
   }
   else { cout << "you should check out nbjet cut!! Error about v_bjet_TL !! " << endl; }
   /// For JES ShiftedUp Jet ///
   if( v_bjetup_TL.size() > 1)
   {
      bJet1Up = v_bjetup_TL.at(0); 
      bJet2Up = v_bjetup_TL.at(1);
   }
   else if ( v_bjetup_TL.size() > 0 )
   {  
      bJet1Up = v_bjetup_TL.at(0);
      for ( int i = 0; i < v_jetup_idx.size(); ++i )
      {
         if ( v_jetup_idx[i] != v_bjetup_idx[0] ){ bJet2Up = v_jetup_TL[i]; break; }
      }
   }
   else if ( v_bjetup_TL.size() == 0 )
   {
      bJet1Up = Jet1Up;
      bJet2Up = Jet2Up;
   }
   else { cout << "you should check out nbjet cut!! Error about v_bjetup_TL !! " << endl; }
   if (isAllSyst == true) 
   { 
      if( v_bjetdn_TL.size() > 1)
      {
         bJet1Dn = v_bjetdn_TL.at(0); 
         bJet2Dn = v_bjetdn_TL.at(1);
      }
      else if ( v_bjetdn_TL.size() > 0 )
      {  
         bJet1Dn = v_bjetdn_TL.at(0);
         for ( int i = 0; i < v_jetdn_idx.size(); ++i )
         {
            if ( v_jetdn_idx[i] != v_bjetdn_idx[0] ){ bJet2Dn = v_jetdn_TL[i]; break; }
         }
      }
      else if ( v_bjetdn_TL.size() == 0 )
      {
         bJet1Dn = Jet1Dn;
         bJet2Dn = Jet2Dn;
      }
      
      else { cout << "you should check out nbjet cut!! Error about v_bjetdn_TL !! " << endl; }
      
      /// For JER ShiftedUp Jet ///
      if( v_bjetresup_TL.size() > 1)
      {
         bJet1JERUp = v_bjetresup_TL.at(0); 
         bJet2JERUp = v_bjetresup_TL.at(1);
      }
      else if ( v_bjetresup_TL.size() > 0 )
      {  
         bJet1JERUp = v_bjetresup_TL.at(0);
         for ( int i = 0; i < v_jetresup_idx.size(); ++i )
         {
            if ( v_jetresup_idx[i] != v_bjetresup_idx[0] ){ bJet2JERUp = v_jetresup_TL[i]; break; }
         }
      }
      else if ( v_bjetresup_TL.size() == 0 )
      {
         bJet1JERUp = Jet1JERUp;
         bJet2JERUp = Jet2JERUp;
      }
      else { cout << "you should check out nbjet cut!! Error about v_bjetup_TL !! " << endl; }
      /// For JER ShiftedDown Jet ///
      if( v_bjetresdn_TL.size() > 1)
      {
         bJet1JERDn = v_bjetresdn_TL.at(0); 
         bJet2JERDn = v_bjetresdn_TL.at(1);
      }
      else if ( v_bjetresdn_TL.size() > 0 )
      {  
         bJet1JERDn = v_bjetresdn_TL.at(0);
         for ( int i = 0; i < v_jetresdn_idx.size(); ++i )
         {
            if ( v_jetresdn_idx[i] != v_bjetresdn_idx[0] ){ bJet2JERDn = v_jetresdn_TL[i]; break; }
         }
      }
      else if ( v_bjetresdn_TL.size() == 0 )
      {
         bJet1JERDn = Jet1JERDn;
         bJet2JERDn = Jet2JERDn;
      }
      else { cout << "you should check out nbjet cut!! Error about v_bjetdn_TL !! " << endl; }
   }
}
/////////////////////
//Double batagging //
/////////////////////
bool ssb_analysis::DoubleBtag( std::vector<int>v_bjets )
{
   doublebtag = false;
   if ( v_bjets.size() >=2 ) {doublebtag = true;}
   return doublebtag;
}
//Kinematic Solver to reconstruct top
void ssb_analysis::KinSol()
{

   xconstraint = 0; yconstraint  = 0;
 
   ksolweight_ = -1.0;
   /// initialize array
   //cout << "At KinSol -- " << "Lep->Px() : " << Lep->Px() << " AnLep->Px() : " << AnLep->Px() << " bJet1->Px() : " << bJet1->Px() << " bJet2->Px() : " << bJet2->Px() << " Met->Px() : " << Met->Px()<< endl;
   xconstraint = Lep->Px() + AnLep->Px() + bJet1->Px() + bJet2->Px() + Met->Px();
   yconstraint = Lep->Py() + AnLep->Py() + bJet1->Py() + bJet2->Py() + Met->Py();

   ssbflsolver->SetConstraints( xconstraint , yconstraint );  

   TtFullLepKinSolver::NeutrinoSolution nuSol1  = ssbflsolver->getNuSolution( (*Lep),   (*AnLep), (*bJet1), (*bJet2) );
   TtFullLepKinSolver::NeutrinoSolution nuSol2  = ssbflsolver->getNuSolution( (*Lep),   (*AnLep), (*bJet2), (*bJet1) );


   if ( nuSol1.weight >nuSol2.weight )
   {
      ksolweight_ = nuSol1.weight;
      (*AnNu)     = nuSol1.neutrino;
      (*Nu)       = nuSol1.neutrinoBar;
      (*bJet)     = (*bJet2);
      (*AnbJet)   = (*bJet1);
   }
   else 
   {  
      ksolweight_ = nuSol2.weight;
      (*AnNu)     = nuSol2.neutrino; 
      (*Nu)       = nuSol2.neutrinoBar;
      (*bJet)     = (*bJet1);
      (*AnbJet)   = (*bJet2);
   }
   //(*Top)      = (*Lep)   + (*AnbJet) + (*AnNu);
   //(*AnTop)    = (*AnLep) + (*bJet)   + (*Nu);
   (*AnTop)     = (*Lep)   + (*AnbJet) + (*AnNu);
   (*Top)       = (*AnLep) + (*bJet)   + (*Nu);
   (*W1)        = (*Lep)   + (*AnNu);
   (*W2)        = (*AnLep) + (*Nu);


}
void ssb_analysis::KinSolSys( TLorentzVector* lep, TLorentzVector* alep, TLorentzVector* bj1, TLorentzVector* bj2, TLorentzVector* met,TString Sys_ )
{

   ksolweight_ = -1.0;
   /// initialize array

   xconstraint = lep->Px() + alep->Px() + bj1->Px() + bj2->Px() + met->Px();
   yconstraint = lep->Py() + alep->Py() + bj1->Py() + bj2->Py() + met->Py();

   ssbflsolver->SetConstraints( xconstraint , yconstraint );  

   TtFullLepKinSolver::NeutrinoSolution nuSol1  = ssbflsolver->getNuSolution( (*lep),   (*alep), (*bj1), (*bj2) );
   TtFullLepKinSolver::NeutrinoSolution nuSol2  = ssbflsolver->getNuSolution( (*lep),   (*alep), (*bj2), (*bj1) );

   if ( nuSol1.weight >nuSol2.weight )
   {
      ksolweight_ = nuSol1.weight;
      (*AnNu)     = nuSol1.neutrino;
      (*Nu)       = nuSol1.neutrinoBar;
      (*bJet)     = (*bj2);
      (*AnbJet)   = (*bj1);
   }
   else 
   {  
      ksolweight_ = nuSol2.weight;
      (*AnNu)     = nuSol2.neutrino; 
      (*Nu)       = nuSol2.neutrinoBar;
      (*bJet)     = (*bj1);
      (*AnbJet)   = (*bj2);
   }
   (*AnTop)     = (*lep)  + (*AnbJet) + (*AnNu);
   (*Top)       = (*alep) + (*bJet)   + (*Nu);
   (*W1)        = (*lep)  + (*AnNu);
   (*W2)        = (*alep) + (*Nu);
}
void ssb_analysis::Bjorken( TLorentzVector *lep, TLorentzVector *alep,TLorentzVector *bj,TLorentzVector *abj,TLorentzVector *nu,TLorentzVector *anu )
{
   x1_bj = -999;
   x2_bj = -999;
   x3_bj = -999;
   x1_bj = ( lep->Pz() + alep->Pz() + bj->Pz() + abj->Pz() + nu->Pz() + anu->Pz() 
           + lep->E()  + alep->E()  + bj->E()  + abj->E()  + nu->E()  + anu->E()  )/(2*4*1000);
   x2_bj = ( lep->E()  + alep->E()  + bj->E()  + abj->E()  + nu->E()  + anu->E() 
           - lep->Pz() - alep->Pz() - bj->Pz() - abj->Pz() - nu->Pz() - anu->Pz() )/(2*4*1000);
   x3_bj = ( x1_bj + x2_bj )*(4*1000);
   
   
}

void ssb_analysis::LepAnLepMisCharge()
{
   laecO3 = 100; laecOb = 100; laecO5 = 100;
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      if ( Lep->Pt() < 300 )
      {
         if ( AnLep->Pt() < 300 ){ laecO3 = 1.5; laecOb = 1.5; laecO5 = 1.5; }
         else                    { laecO3 = 3.5; laecOb = 3.5; laecO5 = 3.5; }
      }
      else 
      {
         if ( AnLep->Pt() < 300 ){ laecO3 = 5.5; laecOb = 5.5; laecO5 = 5.5; }
         else                    { laecO3 = 7.5; laecOb = 7.5; laecO5 = 7.5; }
      }
   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      if( Lep->Pt() >20 && Lep->Pt() < 50 )
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { laecO3 = 1.5; laecOb = 1.5; laecO5 = 1.5; }
               else                               { laecO3 = 3.5; laecOb = 3.5; laecO5 = 3.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 5.5; laecOb = 5.5; laecO5 = 5.5; }
               else                                { laecO3 = 7.5; laecOb = 7.5; laecO5 = 7.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 11.5; laecOb = 11.5; laecO5 = 11.5; }
               else                                { laecO3 = 13.5; laecOb = 13.5; laecO5 = 13.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 15.5; laecOb = 15.5; laecO5 = 15.5; }
               else                                { laecO3 = 17.5; laecOb = 17.5; laecO5 = 17.5; }

            } //
         }
      }
      else // Lep >= 50GeV 
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { laecO3 = 21.5; laecOb = 21.5; laecO5 = 21.5; }
               else                               { laecO3 = 23.5; laecOb = 23.5; laecO5 = 23.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 25.5; laecOb = 25.5; laecO5 = 25.5; }
               else                                { laecO3 = 27.5; laecOb = 27.5; laecO5 = 27.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 31.5; laecOb = 31.5; laecO5 = 31.5; }
               else                                { laecO3 = 33.5; laecOb = 33.5; laecO5 = 33.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 35.5; laecOb = 35.5; laecO5 = 35.5; }
               else                                { laecO3 = 37.5; laecOb = 37.5; laecO5 = 37.5; }

            } //
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      if (TMuon->Pt() < 300 )
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 1.5; laecOb = 1.5; laecO5 = 1.5; } // Barrel
            else                                    { laecO3 = 3.5; laecOb = 3.5; laecO5 = 3.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 5.5; laecOb = 5.5; laecO5 = 5.5; } // Barrel
            else                                    { laecO3 = 7.5; laecOb = 7.5; laecO5 = 7.5; } // Endcap
         }
      } 
      else 
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 11.5; laecOb = 11.5; laecO5 = 11.5; } // Barrel
            else                                    { laecO3 = 13.5; laecOb = 13.5; laecO5 = 13.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 15.5; laecOb = 15.5; laecO5 = 15.5; } // Barrel
            else                                    { laecO3 = 17.5; laecOb = 17.5; laecO5 = 17.5; } // Endcap
         }

      }

   }

}
void ssb_analysis::NewLepAnLepMisCharge()
{
   int disc_leptons = 0;
   int value_mis = 0;
   v_MisCharge.clear();
   if ( v_recocp_O.size() != 13) {return;} /// After getting CPVarialbes, You can execute this fucntion ...///
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      if ( Lep->Pt() < 300 )
      {
         if ( AnLep->Pt() < 300 ){ disc_leptons = 1.5;}
         else                    { disc_leptons = 3.5;}
      }
      else 
      {
         if ( AnLep->Pt() < 300 ){ disc_leptons = 5.5;}
         else                    { disc_leptons = 7.5;}
      }
   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      if( Lep->Pt() >20 && Lep->Pt() < 50 )
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { disc_leptons = 1.5; }
               else                               { disc_leptons = 3.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 5.5; }
               else                                { disc_leptons = 7.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 11.5; }
               else                                { disc_leptons = 13.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 15.5; }
               else                                { disc_leptons = 17.5; }

            } //
         }
      }
      else // Lep >= 50GeV 
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { disc_leptons = 21.5; }
               else                               { disc_leptons = 23.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 25.5; }
               else                                { disc_leptons = 27.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 31.5; }
               else                                { disc_leptons = 33.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { disc_leptons = 35.5; }
               else                                { disc_leptons = 37.5; }

            } //
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      if (TMuon->Pt() < 300 )
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { disc_leptons = 1.5; } // Barrel
            else                                    { disc_leptons = 3.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { disc_leptons = 5.5; } // Barrel
            else                                    { disc_leptons = 7.5; } // Endcap
         }
      } 
      else 
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { disc_leptons = 11.5; } // Barrel
            else                                    { disc_leptons = 13.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { disc_leptons = 15.5; } // Barrel
            else                                    { disc_leptons = 17.5; } // Endcap
         }

      }

   }
   /// Check Up Positive and Negative !!!
   for ( int i = 0; i < v_recocp_O.size(); ++i )
   {
      if (v_recocp_O[i] > 0.0000) {
         value_mis = disc_leptons;
      }
      else { value_mis = 50+disc_leptons; }
      v_MisCharge.push_back(value_mis);
   }


}
double ssb_analysis::JetPtResol(TLorentzVector *jet)
{
   double jpr_;
   double jpt_;
   double jeta_;
   jpt_ = jet->Pt();
   jeta_ = jet->Eta();
   jpr_ = ssbcpviol->getJPR(jeta_,jpt_);
   return jpr_;
}

double ssb_analysis::CalTopPtRewight(TString Opt_)
{
   double topptreweight =1.0;
   if (Opt_.Contains("None")){return topptreweight;}
   else if (Opt_.Contains("Apply"))
   {
      if ( TString(FileName_).Contains( "TTJets" ) )
      {
         TLorentzVector* genTop   = new TLorentzVector();
         TLorentzVector* genAnTop = new TLorentzVector();
         if (GenTop->GetEntriesFast() >0 && GenAnTop->GetEntriesFast() >0 ){
            genTop   = (TLorentzVector*)GenTop->At(0);
            genAnTop = (TLorentzVector*)GenAnTop->At(0);
         }
         else {
         // For Herwig Case //
            genTop = FindGenPar(6);
            genAnTop = FindGenPar(-6);
         }
         double a = 0.0615;
         double b = -0.0005;
         double toppt = genTop->Pt();
         double antoppt = genAnTop->Pt();
         double topsf = exp(a+b*toppt);
         double antopsf = exp(a+b*antoppt);
         topptreweight = sqrt(topsf*antopsf);
      }
      else {topptreweight = 1.0;}
   }
   else {return 1.0;}
   return topptreweight;
}
void ssb_analysis::TopPtReweightApply()
{
   double evt_weight_beforetopptweight_ = evt_weight_;
   double topptreweight =1.0;
   if ( TString(FileName_).Contains( "TTJets" ) )
   {
      topptreweight = CalTopPtRewight(TopPtSys);
   }
   else {topptreweight =1;}
   evt_weight_ = evt_weight_*topptreweight;

   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      {  
         double topptrew_allsys =1; 
         if ( !TString(FileName_).Contains("Data") )
         {
            if ( TString(v_SystFullName[i]).Contains("TopPt") && TString(FileName_).Contains("TTJets") ){ v_SystEvt[i] = v_SystEvt[i]*CalTopPtRewight("Apply");}
            else { v_SystEvt[i]=v_SystEvt[i];}
         }
         else {
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   }
}
TLorentzVector* ssb_analysis::JERSmearing(TLorentzVector *jet, int idx_, TString op_)
{
   double jerfrac_ = 1.;
   if ( op_.Contains( "Cen") == true ) { jerfrac_ = Jet_EnergyResolution_SF->at(idx_); }
   else if ( op_ == "Up") { jerfrac_ = Jet_EnergyResolution_SFUp->at(idx_); }
   else if ( op_ == "Down") { jerfrac_ = Jet_EnergyResolution_SFDown->at(idx_); }
   else {cout << "Check out your JERSmeraing option !! " << endl;jerfrac_ =1.0; }

   TLorentzVector* jetsmerd = new TLorentzVector();
   if ( dojer&!TString(FileName_).Contains( "Data" ) )
   {
      jetsmerd->SetPtEtaPhiE(jet->Pt()*jerfrac_,
                             jet->Eta(),jet->Phi(),
                            jet->Energy()*jerfrac_);
   }
   else {
      jetsmerd = jet;
   }
   //cout << "before smeared jet pt : " << jet->Pt()  << " Eta : " << jet->Eta() << " phi " << jet->Phi() << " Energy : " << jet->Energy() << endl;
   //cout << "after smeared jetsmerd pt : " << jetsmerd->Pt()  << " Eta : " << jetsmerd->Eta() << " phi " << jetsmerd->Phi() << " Energy : " << jetsmerd->Energy() << endl;
   return  jetsmerd;
}
TLorentzVector* ssb_analysis::ApplyJetPtRes(TLorentzVector *jet, int idx_, TString op_)
{
   double jerfrac_ = 1.;
   if ( op_.Contains( "None") == true ) { return jet; }
   else if ( op_ == "Up") { jerfrac_ = 1+Jet_EnergyResolution_MC->at(idx_); }
   else if ( op_ == "Down") { jerfrac_ = 1-Jet_EnergyResolution_MC->at(idx_); }
   else {cout << "Check out your JERSmeraing option !! " << endl;jerfrac_ =1.0; }
   //cout << "jerfrac_ : " << jerfrac_ << endl;
   TLorentzVector* variedjet = new TLorentzVector();
   if ( !TString(FileName_).Contains( "Data" ) )
   {
      variedjet->SetPtEtaPhiE(jet->Pt()*jerfrac_,
                             jet->Eta(),jet->Phi(),
                            jet->Energy()*jerfrac_);
   }
   else {
      variedjet = jet;
   }
   //cout << "before smeared jet pt : " << jet->Pt()  << " Eta : " << jet->Eta() << " phi " << jet->Phi() << " Energy : " << jet->Energy() << endl;
   //cout << "after smeared variedjet pt : " << variedjet->Pt()  << " Eta : " << variedjet->Eta() << " phi " << variedjet->Phi() << " Energy : " << variedjet->Energy() << endl;
   return  variedjet;
}
TLorentzVector* ssb_analysis::ApplyJetPhiRes(TLorentzVector *jet, int idx_, TString op_)
{
   double jerfrac_ = 1.;
   if ( op_.Contains( "None") == true ) { return jet; }
   else if ( op_ == "Up") { jerfrac_ = 1+Jet_PhiResolution_MC->at(idx_); }
   else if ( op_ == "Down") { jerfrac_ = 1-Jet_PhiResolution_MC->at(idx_); }
   else {cout << "Check out your JERSmeraing option !! " << endl;jerfrac_ =1.0; }
   TLorentzVector* variedjet = new TLorentzVector();
   if ( !TString(FileName_).Contains( "Data" ) )
   {
      variedjet->SetPtEtaPhiE(jet->Pt(),
                             jet->Eta(),jet->Phi()*jerfrac_,
                            jet->Energy());
   }
   else {
      variedjet = jet;
   }
   //cout << "before smeared jet pt : " << jet->Pt()  << " Eta : " << jet->Eta() << " phi " << jet->Phi() << " Energy : " << jet->Energy() << endl;
   //cout << "after smeared variedjet pt : " << variedjet->Pt()  << " Eta : " << variedjet->Eta() << " phi " << variedjet->Phi() << " Energy : " << variedjet->Energy() << endl;
   return  variedjet;
}
void ssb_analysis::ApplyJetPtPhiDilution()
{
   /// Find bJet1 & bJet2 in v_jet_TL !!!
   bool isPtDilu = false;
   bool isPhiDilu = false;
   //if ( JetPtPhiDil == false ){cout << "No excute ApplyJetPtPhiDilution JetPtPhiDil false " << endl;return;}
   if (bAndbBarDil == true) {return;}
   if (JetPtResDil != "None" && JetPtResDil != "none") {isPtDilu = true;}
   if (JetPhiResDil != "None" && JetPhiResDil != "none") {isPhiDilu = true;}
   if (isPtDilu && isPhiDilu) {cout << "You should check out your options of pt & phi dilution"<< endl;return;}
   if (isAllSyst&& (isPtDilu ==true ||isPhiDilu == true )) {cout << "You should turn off your Allsys options "<< endl;return;}
   if ( ( isPtDilu == false ) && (isPhiDilu==false) ){cout << "No excute ApplyJetPtPhiDilution" << endl;return;}
   bool findbjet1 = false;
   bool findbjet2 = false;
   int bjet1ind = -1;
   int bjet2ind = -1;
   int njetidx = 0;

   double bj1pt = bJet1->Pt();
   double bj1px = bJet1->Px();
   double bj1py = bJet1->Py();

   double bj2pt = bJet2->Pt();
   double bj2px = bJet2->Px();
   double bj2py = bJet2->Py();

   std::vector<double> v_diff_pt;
   std::vector<double> v_diff_px;
   std::vector<double> v_diff_py;

   TLorentzVector* empty = new TLorentzVector();
   TString b1op_ = "";
   TString b2op_ = "";

   if (FindSameObj(empty,bJet1) ==true ) {return;}
   if (FindSameObj(empty,bJet2) ==true ) {return;}
   for (std::vector<TLorentzVector*>::iterator it = v_jet_TL.begin(); it != v_jet_TL.end(); ++it)
   {
      //cout << " v_jet_TL !! " <<endl;
      if (FindSameObj(bJet1,(*it))) { findbjet1 = true; bjet1ind = njetidx; };
      if (FindSameObj(bJet2,(*it))) { findbjet2 = true; bjet2ind = njetidx; };
      njetidx++;
   }
   /// Apply JetPt-1sig.Vari !!
   //cout << "JetPtResDil : " << JetPtResDil << endl;
   if (isPtDilu){
      if ( JetPtResDil == "UpUp" || JetPtResDil == "upup") { b1op_ = "Up"; b2op_ = "Up"; }
      else if ( JetPtResDil == "UpDown" || JetPtResDil == "updown") { b1op_ = "Up"; b2op_ = "Down"; }
      else if ( JetPtResDil == "DownUp" || JetPtResDil == "downup") { b1op_ = "Down"; b2op_ = "Up"; }
      else if ( JetPtResDil == "DownDown" || JetPtResDil == "downdown") { b1op_ = "Down"; b2op_ = "Down"; }
      else { b1op_ = "None"; b2op_ = "None`"; }
      bJet1 = ApplyJetPtRes(bJet1,v_jet_idx[bjet1ind], b1op_);
      bJet2 = ApplyJetPtRes(bJet2,v_jet_idx[bjet2ind], b2op_); 
   }
   if (isPhiDilu){
      if ( JetPhiResDil == "UpUp" || JetPhiResDil == "upup") { b1op_ = "Up"; b2op_ = "Up"; }
      else if ( JetPhiResDil == "UpDown" || JetPhiResDil == "updown") { b1op_ = "Up"; b2op_ = "Down"; }
      else if ( JetPhiResDil == "DownUp" || JetPhiResDil == "downup") { b1op_ = "Down"; b2op_ = "Up"; }
      else if ( JetPhiResDil == "DownDown" || JetPhiResDil == "downdown") { b1op_ = "Down"; b2op_ = "Down"; }
      else { b1op_ = "None"; b2op_ = "None`"; }
      bJet1 = ApplyJetPhiRes(bJet1,v_jet_idx[bjet1ind], b1op_);
      bJet2 = ApplyJetPhiRes(bJet2,v_jet_idx[bjet2ind], b2op_); 
   }
   /// Changed the component of v_jet_TL !!
   v_jet_TL[bjet1ind] = bJet1;
   v_jet_TL[bjet2ind] = bJet2;

   v_diff_pt.push_back(bJet1->Pt() - bj1pt);
   v_diff_pt.push_back(bJet2->Pt() - bj2pt);
   v_diff_px.push_back(bJet1->Px() - bj1px);
   v_diff_px.push_back(bJet2->Px() - bj2px);
   v_diff_py.push_back(bJet1->Py() - bj1py);
   v_diff_py.push_back(bJet2->Py() - bj2py);

   /////////////////////
   /// Calculate MET /// 
   /////////////////////
   Met = METSmear(v_diff_pt,v_diff_px,v_diff_py,Met);
}

bool ssb_analysis::FindSameObj(TLorentzVector* ref, TLorentzVector* tar)
{
   bool sameobj = false;
   if (ref->DeltaR(*tar) == 0 && ref->Pt() == tar->Pt() && ref->Eta() == tar->Eta()){sameobj = true;}
   return sameobj;
}

/*void ssb_analysis::ApplyJetPtRes(TString op_)
{
   double bjerfrac1_ = 1.;
   double bjerfrac2_ = 1.;
   //cout << "Check out v_bjet size "<< v_bjet_idx.size() << endl;
   if (v_bjet_idx.size() <2){cout << "Not 2 bjets... " << endl; return;}
   int idxb1 = v_bjet_idx[0]; 
   int idxb2 = v_bjet_idx[1];

   double bj1pt = bJet1->Pt();
   double bj1eta = bJet1->Eta();
   double bj1phi = bJet1->Phi();
   double bj1energy = bJet1->Energy();
   double bj1px = bJet1->Px();
   double bj1py = bJet1->Py();

   double bj2pt = bJet2->Pt();
   double bj2eta = bJet2->Eta();
   double bj2phi = bJet2->Phi();
   double bj2energy = bJet2->Energy();
   double bj2px = bJet2->Px();
   double bj2py = bJet2->Py();
   std::vector<double> v_diff_pt;
   std::vector<double> v_diff_px;
   std::vector<double> v_diff_py;
   if ( op_.Contains( "None") == true ) { return; }
   else if ( op_ == "Up") 
   {
      bjerfrac1_ = 1+Jet_EnergyResolution_MC->at(idxb1); 
      bjerfrac2_ = 1+Jet_EnergyResolution_MC->at(idxb2); 
   }
   else if ( op_ == "Down") {
      bjerfrac1_ = 1-Jet_EnergyResolution_MC->at(idxb1); 
      bjerfrac2_ = 1-Jet_EnergyResolution_MC->at(idxb2); 
   }
   else {cout << "Check out your JERSmeraing option !! " << endl;bjerfrac1_ =1.0; bjerfrac2_ =1.0; }
   ////////////////////////////
   /// Change bJet1 & bJet2 ///
   ////////////////////////////
   bJet1->SetPtEtaPhiE(bj1pt*bjerfrac1_,bj1eta,bj1phi,bj1energy*bjerfrac1_);
   bJet2->SetPtEtaPhiE(bj2pt*bjerfrac2_,bj2eta,bj2phi,bj2energy*bjerfrac2_);

   cout 
   << "before apply pt res bJet1 pt : " << bj1pt << " eta : "  << bj1eta << " phi : " << bj1phi << " energy : " << bj1energy  
   << endl;
   cout 
   << "after apply pt res bJet1 pt : " << bJet1->Pt() << " eta : "  << bJet1->Eta() << " phi : " << bJet1->Phi() << " energy : " << bJet1->Energy()
   << endl;
   v_diff_pt.push_back(bJet1->Pt() - bj1pt);
   v_diff_px.push_back(bJet1->Px() - bj1px);
   v_diff_py.push_back(bJet1->Py() - bj1py);
   /////////////////////
   /// Calculate MET /// 
   /////////////////////
   Met = METSmear(v_diff_pt,v_diff_px,v_diff_py,Met);
}*/
void ssb_analysis::MakeVecforTL()
{
}
void ssb_analysis::ClearVectors()
{
   v_lepton_Id->clear();
   v_lepton_iso->clear();
   v_jet_Id->clear();
   v_mlep_iso->clear();
   v_elec_iso->clear();
   v_muon_Id->clear();  // for MuEl channel
   v_electron_Id->clear(); // for ElEl and MuEl channels

   /// User variables for JetCleanning cut
   v_mu_iso_jcl->clear();
   v_mu_Id_jcl->clear();
   v_el_iso_jcl->clear();
   v_el_Id_jcl->clear();
}
void ssb_analysis::ReadDupleList()
{
   if (!FileName_.Contains("Data")){ return; } 
   if (FileName_.Contains("Data_Double")){ return; } 
   if (FileName_.Contains("Data_MuonElec")){ return; } 
   /// Duplicate Event ///
   std::map<int,std::vector<int> > m_test;
   TString evtfileName = FileName_;
   char bar = '_';
/*   if (Decaymode == "dimu" ) {evtfileName.ReplaceAll("Single","Double");}
   else if ( Decaymode == "dielec" ){ evtfileName.ReplaceAll("Single","Double"); }
   else if ( Decaymode == "muel" ){ evtfileName.ReplaceAll("Single","Double"); }
   else { return;}*/
   if (Decaymode.Contains("dimu") ) {evtfileName="Data_DoubleMuon_Run2016";}
   else if ( Decaymode.Contains("dielec") ){ evtfileName="Data_DoubleEG_Run2016"; }
   else if ( Decaymode.Contains("muel") ){ evtfileName="Data_MuonEG_Run2016"; }
   else { return;}
   int NumRunFiles_ = 1;
   if (FileName_.Contains("Run2016B")){ NumRunFiles_ = 176; evtfileName.ReplaceAll("Run2016","Run2016B"); } 
   else if (FileName_.Contains("Run2016C"))   { NumRunFiles_ = 58;  evtfileName.ReplaceAll("Run2016","Run2016C"); } 
   else if (FileName_.Contains("Run2016D"))   { NumRunFiles_ = 98;  evtfileName.ReplaceAll("Run2016","Run2016D"); } 
   else if (FileName_.Contains("Run2016E"))   { NumRunFiles_ = 83;  evtfileName.ReplaceAll("Run2016","Run2016E"); } 
   else if (FileName_.Contains("Run2016F"))   { NumRunFiles_ = 61;  evtfileName.ReplaceAll("Run2016","Run2016F"); } 
   else if (FileName_.Contains("Run2016G"))   { NumRunFiles_ = 143; evtfileName.ReplaceAll("Run2016","Run2016G"); } 
   else if (FileName_.Contains("Run2016HV2")) { NumRunFiles_ = 155; evtfileName.ReplaceAll("Run2016","Run2016HV2"); } 
   else if (FileName_.Contains("Run2016HV3")) { NumRunFiles_ = 5;   evtfileName.ReplaceAll("Run2016","Run2016HV3"); }
   else {  cout << "Check Out Your InputFiles !!" << endl;  return;  } 

   FILE *XSectionInfo;
   char LumiNum[1000];
   char EvtNum[1000];
   char RunNum[1000];
   char col1[1000];
   int lumi_number;
   int evt_number;
   int run_number;
   dilepevt.clear();
//   int testlines =0;
   for (int i = 1; i < NumRunFiles_+1; ++i)
   {
      ostringstream ostr;
      ostr.str("");
      ostr << i;
      string ostr_ = "_" + ostr.str() + ".txt";
      string evt_file = evtfileName.Data()  + ostr_;
      string slash = "/";
      string tttt= evtfileName.Data();
//      string evt_filePath = "../EventList/" + tttt + slash + evt_file;
      //string evt_filePath = "/d3/scratch/sha/Analyses/SSB/MyAnalysis/25ns_Anlysis/Run2016/RemoveDuplicatedEvt/EventList/" + tttt + slash + evt_file;
      //string evt_filePath = "/d0/scratch/sha/Analyses/SSB/MyAnalysis/25ns_Anlysis/Run2016/RemoveDuplicatedEvt/EventList/" + tttt + slash + evt_file;
      string evt_filePath = "/u/user/gcho/TopPhysics/CPV/CMSSW_8_0_26_patch1/src/TOP-18-007/SSBAnalysis/AnalysisCode/25ns_Anlysis/Run2016/RemoveDuplicatedEvt/EventList/" + tttt + slash + evt_file;
      cout << "evt_filePath : " << evt_filePath << endl;
      XSectionInfo = fopen(evt_filePath.c_str(),"r");                                                                               

      if (XSectionInfo !=NULL){
         while (fscanf(XSectionInfo, "%s %s %d %s %d %s %d \n",RunNum,col1,&run_number, LumiNum,&lumi_number, EvtNum,& evt_number ) != EOF)
         {  //testlines++;
            dilepevt[run_number][lumi_number][evt_number]=1;
            //273730 Info_Luminosity: 1993 Info_EventNumber: -2023987825
            //if (run_number == 273730 && lumi_number == 1993 && evt_number == -2023987825  ) { cout << "CHECK ~!!!!!!:" << run_number<< " : "<< lumi_number << " : " <<evt_number  << endl;}
         }// End of ReadFile
         fclose(XSectionInfo);
      }
      else {
         delete XSectionInfo;
      }
   }// End of loop of fileList //
   // Check out Dilepton List 
//   cout << " testlines : " << testlines << endl;
/*   int test2lines =0.0; 
   /// For Debug ///
   for (map<int, map<int, map<int, int> > >::iterator iter = dilepevt.begin(); iter != dilepevt.end(); ++iter )
   {  //cout << "iter first : " << iter->first << endl;
      for ( map<int, map<int, int> >::iterator itsub = iter->second.begin(); itsub != iter->second.end(); ++itsub)
      {  
         for (map<int, int>::iterator ivec = itsub->second.begin(); ivec != itsub->second.end(); ++ivec)
         {
           test2lines++;
           //cout << "iter->first : " << iter->first  << "itsub->first : " << itsub->first << " ivec->first : " << ivec->first << endl; 
         }
      }
   }

   cout << " test2lines: " << test2lines << endl;*/
 
   cout << "End of Read Dilep Evt" << endl;
}
bool ssb_analysis::RMDuplEvt(int run_, int lumi_, int evt_num_)
{
   /// Duplicate Event ///
   bool vetoduple = false;
   if (dilepevt[run_][lumi_][evt_num_] == 0) {vetoduple =true;} // No Duplicated Event //
//   else {cout << "dilepevt[run_][lumi_][evt_num_] : " << dilepevt[run_][lumi_][evt_num_] << endl; }
//   cout << "Remove Duplicate Event " << endl;
   return vetoduple; //if vetoduple is false(1), this event is duplicated
}
void ssb_analysis::SetGenLepAnLep()
{
   GenLep   = new TLorentzVector(); 
   GenAnLep = new TLorentzVector();
   for (int i = 0; i < GenPar_pdgId->size(); ++i)
   {
      if ( GenPar_pdgId->at(i) == 2212 ) {continue;}
      TLorentzVector* par = (TLorentzVector*)GenPar->At(i);
      if ( TString(Decaymode).Contains( "dielec" ) ) {
         if ( GenPar_pdgId->at(i) == 11 ) { GenLep = par; }
         if ( GenPar_pdgId->at(i) == -11 ) { GenAnLep = par; }
      }
      else if ( TString(Decaymode).Contains( "dimuon" ) ) {
         if ( GenPar_pdgId->at(i) == 13 ) { GenLep = par; }
         if ( GenPar_pdgId->at(i) == -13 ) { GenAnLep = par; }
      }
      else if ( TString(Decaymode).Contains( "muel" ) ) {
         if ( GenPar_pdgId->at(i) == 11 ) { GenLep = par; }
         if ( GenPar_pdgId->at(i) == -11 ) { GenAnLep = par; }
         if ( GenPar_pdgId->at(i) == 13 ) { GenLep = par; }
         if ( GenPar_pdgId->at(i) == -13 ) { GenAnLep = par; }
      }
      else { cout << "Check Your Decay Mode !!!" << endl;return;}
   }
//   cout << "GenLep pt : " << GenLep->Pt() << endl;
//   cout << "GenAnLep pt : " << GenAnLep->Pt() << endl;

}
TLorentzVector* ssb_analysis::FindGenPar(int pdgid_)
{
    TLorentzVector* tarpar_ = new TLorentzVector();
    bool find_par = false;
    for (int i = 0; i < GenPar_pdgId->size(); ++i)
    {
       if ( GenPar_pdgId->at(i) == pdgid_ ) { tarpar_ = (TLorentzVector*)GenPar->At(i); find_par = true;} 
    }
    if (find_par == false)cout << "There is no target partice !! pdgId : " << pdgid_  << endl;
    return tarpar_;
}

void ssb_analysis::BandBbarJetDiff()
{
   /// This Function is for dilution factor study of difference of response between b and bbar jet ///
   /// Fraction of DeltaR is 0.001.
   if ( TString(FileName_).Contains( "Data" ) ){return;}
   if ( bAndbBarDil == false ){return;}
//   cout << "bAndbBarDil : " << bAndbBarDil  << endl;
   double frac = 0.0011;
   //double bjpx_,bjpy_,bjpz_,bje_; 
   //double anbjpx_,anbjpy_,anbjpz_,anbje_;
   double bjpt_,bjeta_,bjphi_,bje_; 
   double anbjpt_,anbjeta_,anbjphi_,anbje_;
 
   double fac_bj   = 1-0.5*0.0011;
   double fac_anbj = 1+0.5*0.0011;
   // Get BJet 4 - momentum
   /*bjpx_ = bJet->Px();
   bjpy_ = bJet->Py();
   bjpz_ = bJet->Pz();*/
   bjpt_ = bJet->Pt();
   bjeta_ = bJet->Eta();
   bjphi_ = bJet->Phi();
   bje_ = bJet->Energy();

   // Get AnBJet 4 - momentum
   /*anbjpx_ = AnbJet->Px();
   anbjpy_ = AnbJet->Py();
   anbjpz_ = AnbJet->Pz();*/
   anbjpt_ = AnbJet->Pt();
   anbjeta_ = AnbJet->Eta();
   anbjphi_ = AnbJet->Phi();
   anbje_ = AnbJet->Energy();

   /*cout << " ----------------- before bAndbbarJet Diff---------------------- " << endl;
   cout << "bJet px : " << bJet->Px() << " py : " << bJet->Py() << " pz : " <<  bJet->Pz() << " energy : " << bJet->Energy() 
   << "AnbJet px : " << AnbJet->Px() << " py : " << AnbJet->Py() << " pz : " <<  AnbJet->Pz() << " energy : " << AnbJet->Energy() << endl;*/

   // Apply Diff. of Response
   //bJet->SetPxPyPzE(bjpx_*fac_bj, bjpy_*fac_bj , bjpz_*fac_bj, bje_*fac_bj);
   bJet->SetPtEtaPhiE(bjpt_*fac_bj, bjeta_, bjphi_, bje_*fac_bj);
   //AnbJet->SetPxPyPzE(anbjpx_*fac_anbj,anbjpy_*fac_anbj,anbjpz_*fac_anbj,anbje_*fac_anbj);
   AnbJet->SetPtEtaPhiE(anbjpt_*fac_anbj,anbjeta_,anbjphi_,anbje_*fac_anbj);

   // For Debug 
   /*cout << " ----------------- after bAndbbarJet Diff---------------------- " << endl;
   cout << " fac_bj ? : " << fac_bj << " fac_anbj : " << fac_anbj << endl;
   cout << "bJet px : " << bJet->Px() << " py : " << bJet->Py() << " pz : " <<  bJet->Pz() << " energy : " << bJet->Energy() 
   << "AnbJet px : " << AnbJet->Px() << " py : " << AnbJet->Py() << " pz : " <<  AnbJet->Pz() << " energy : " << AnbJet->Energy() << endl;*/


   // Apply Top 
   (*Top)   = (*AnLep) + (*bJet)   + (*Nu);
   (*AnTop) = (*Lep)   + (*AnbJet) + (*AnNu);

}
void ssb_analysis::SetUpKINObs()
{
   isKinSol=false;
   v_leptons_VLV.clear();
   v_jets_VLV.clear();
   v_bjets_VLV.clear();
   v_lepidx_KIN.clear();
   v_anlepidx_KIN.clear();
   v_jetidx_KIN.clear();
   v_bjetidx_KIN.clear();
   v_btagging_KIN.clear();
   /// lepton ///
   v_leptons_VLV.push_back(common::TLVtoLV(*Lep));
   v_lepidx_KIN.push_back(0);
   v_leptons_VLV.push_back(common::TLVtoLV(*AnLep));
   v_anlepidx_KIN.push_back(1);
   /*cout << "jet info. for kinematic solver"  << endl;
   cout << v_jet_idx.size() << " " << v_jet_TL.size() << endl;*/
   const KinematicReconstruction* kinematicReconstruction(0); 
   kinematicReconstruction = new KinematicReconstruction(1, true);

   const LV met_LV = common::TLVtoLV(*Met);

   for (int ijet = 0; ijet < v_jet_idx.size(); ++ijet)
   {
      int idx_jet = v_jet_idx[ijet];
      if(Jet_bDisc->at(idx_jet)  > bdisccut ){
         v_bjets_VLV.push_back(common::TLVtoLV(*v_jet_TL[ijet]));
         //v_bjetidx_KIN.push_back(v_jet_idx[ijet]);
         v_bjetidx_KIN.push_back(ijet);
      }
      v_jets_VLV.push_back(common::TLVtoLV(*v_jet_TL[ijet]));
      v_jetidx_KIN.push_back(ijet);
      v_btagging_KIN.push_back(Jet_bDisc->at(idx_jet));
   }
   KinematicReconstructionSolutions kinematicReconstructionSolutions = kinematicReconstruction->solutions(v_lepidx_KIN,v_anlepidx_KIN, v_jetidx_KIN,  v_bjetidx_KIN,  v_leptons_VLV, v_jets_VLV, v_btagging_KIN, met_LV);
   //cout << "Num Sol : " << kinematicReconstructionSolutions.numberOfSolutions() << endl;
   //cout << "MET ? " << met_LV.pt() << endl;
   if (kinematicReconstructionSolutions.numberOfSolutions())
   {
      isKinSol= true;
      LV top1 = kinematicReconstructionSolutions.solution().top();
      LV top2 = kinematicReconstructionSolutions.solution().antiTop();
      LV bjet1 = kinematicReconstructionSolutions.solution().bjet();
      LV bjet2 = kinematicReconstructionSolutions.solution().antiBjet();
      LV neutrino1 = kinematicReconstructionSolutions.solution().neutrino();
      LV neutrino2 = kinematicReconstructionSolutions.solution().antiNeutrino();
      //kinematicReconstructionSolutions.solution().print();
      //Top = new TLorentzVector(common::LVtoTLV(top1));   
      Top       = new TLorentzVector(common::LVtoTLV(top1));
      AnTop     = new TLorentzVector(common::LVtoTLV(top2));
      bJet      = new TLorentzVector(common::LVtoTLV(bjet1));
      AnbJet    = new TLorentzVector(common::LVtoTLV(bjet2));
      Nu        = new TLorentzVector(common::LVtoTLV(neutrino1));
      AnNu      = new TLorentzVector(common::LVtoTLV(neutrino2));

      (*W1)        = (*Lep)  + (*AnNu);
      (*W2)        = (*AnLep)  + (*Nu);

   }
   //delete 
   delete kinematicReconstruction;
}
//void ssb_analysis::SetUpKINObsSyst()
void ssb_analysis::SetUpKINObsSyst(vector<int> v_jetsys_idx, vector<TLorentzVector*> v_jetsys_TL, TLorentzVector* metsys)
{
   isKinSol=false;
   v_leptons_VLV.clear();
   v_jets_VLV.clear();
   v_bjets_VLV.clear();
   v_lepidx_KIN.clear();
   v_anlepidx_KIN.clear();
   v_jetidx_KIN.clear();
   v_bjetidx_KIN.clear();
   v_btagging_KIN.clear();
   /// lepton ///
   v_leptons_VLV.push_back(common::TLVtoLV(*Lep));
   v_lepidx_KIN.push_back(0);
   v_leptons_VLV.push_back(common::TLVtoLV(*AnLep));
   v_anlepidx_KIN.push_back(1);
   //cout << "jet info. for kinematic solver"  << endl;
   //cout << v_jet_idx.size() << " " << v_jetsys_TL.size() << endl;
   const KinematicReconstruction* kinematicReconstruction(0); 
   kinematicReconstruction = new KinematicReconstruction(1, true);

   const LV met_LV = common::TLVtoLV(*metsys);

   for (int ijet = 0; ijet < v_jetsys_idx.size(); ++ijet)
   {
      int idx_jet = v_jetsys_idx[ijet];
      if(Jet_bDisc->at(idx_jet)  > bdisccut ){
         v_bjets_VLV.push_back(common::TLVtoLV(*v_jetsys_TL[ijet]));
         //v_bjetidx_KIN.push_back(v_jetsys_idx[ijet]);
         v_bjetidx_KIN.push_back(ijet);
      }
      v_jets_VLV.push_back(common::TLVtoLV(*v_jetsys_TL[ijet]));
      v_jetidx_KIN.push_back(ijet);
      v_btagging_KIN.push_back(Jet_bDisc->at(idx_jet));
   }
   KinematicReconstructionSolutions kinematicReconstructionSolutions = kinematicReconstruction->solutions(v_lepidx_KIN,v_anlepidx_KIN, v_jetidx_KIN,  v_bjetidx_KIN,  v_leptons_VLV, v_jets_VLV, v_btagging_KIN, met_LV);
   if (kinematicReconstructionSolutions.numberOfSolutions())
   {
      isKinSol= true;
      LV top1 = kinematicReconstructionSolutions.solution().top();
      LV top2 = kinematicReconstructionSolutions.solution().antiTop();
      LV bjet1 = kinematicReconstructionSolutions.solution().bjet();
      LV bjet2 = kinematicReconstructionSolutions.solution().antiBjet();
      LV neutrino1 = kinematicReconstructionSolutions.solution().neutrino();
      LV neutrino2 = kinematicReconstructionSolutions.solution().antiNeutrino();
      //Top = new TLorentzVector(common::LVtoTLV(top1));   
      Top       = new TLorentzVector(common::LVtoTLV(top1));
      AnTop     = new TLorentzVector(common::LVtoTLV(top2));
      bJet      = new TLorentzVector(common::LVtoTLV(bjet1));
      AnbJet    = new TLorentzVector(common::LVtoTLV(bjet2));
      Nu        = new TLorentzVector(common::LVtoTLV(neutrino1));
      AnNu      = new TLorentzVector(common::LVtoTLV(neutrino2));

      (*W1)        = (*Lep)  + (*AnNu);
      (*W2)        = (*AnLep)  + (*Nu);
      //cout << "Top pt " << Top->Pt() << endl;
   }
   delete kinematicReconstruction;
}