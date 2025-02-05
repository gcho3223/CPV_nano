#define ssb_analysis_cxx

#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "TMath.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

#include "./../interface/ssb_analysis.hpp"
#include "./../CommonTools.hpp"

using namespace std;

//void ssb_analysis::SetNanoAODHeader(Long64_t currentEntry) {
void ssb_analysis::SetNanoAODHeader() {


    std::string inputFileName(FileName_.Data());

    if (inputFileName.find("2016") != std::string::npos) {
        is2016 = true;
        is2017 = false;
        isRun3 = false;
        // 예제 트리거, 실제로 존재하는 경우만 설정합니다.
        if (fChain->GetBranch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ")) {
            
            map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"] = 0;
            fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"]);
        }

    } else if (inputFileName.find("2017") != std::string::npos) {
        is2016 = false;
        is2017 = true;
        isRun3 = false;
        
        #include "../analysis/2017/NanoAODTree.h" 
        if (fChain->GetBranch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")) {
            map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"] = 0;
            fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"]);
        }

    } else if (inputFileName.find("Run3") != std::string::npos) {
        is2016 = false;
        is2017 = false;
        isRun3 = true;
        #include "../analysis/Run3/NanoAODTree.h" 
        
        if (fChain->GetBranch("HLT_IsoMu24")) {
            map_TrigPathName_passed["HLT_IsoMu24"] = 0;
            fChain->SetBranchAddress("HLT_IsoMu24", &map_TrigPathName_passed["HLT_IsoMu24"]);
        }

    } else {
        std::cerr << "Error: Could not determine NanoAOD header version from file name." << std::endl;
        exit(1);
    }
}

void ssb_analysis::SetConfig()
{

}


void ssb_analysis::Loop( char *logfile )
{

   //GetTotalEvent();
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
   //cout << "Ntuple Total Event Check !! " << NtupletotalEvent << endl;

   //MCSF();
   cout << "MC_SF Check !! " << mc_sf_ << endl;

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
      //MakeVecforTL();

      // initailizing TLorentzVector
      //TLVInitial();

      evt_weight_ = 1;

      //GetVariables();
      //MCSFApply();
      //////////////////////

      ///////////////////////////////
      /// ** Dilepton Analysis ** ///
      ///////////////////////////////
      if ( TString(Decaymode).Contains( "dielec" ) ||
           TString(Decaymode).Contains( "dimuon" ) || 
           TString(Decaymode).Contains( "muel" )     ){

         /*if (fChain->GetBranch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ")) {
            
            map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"] = 0;
            fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"]);
         }*/
         //cout << "fCurrent : " << fCurrent << endl;
         //Trigger(nb);
         //Trigger();
         //if (map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"] != HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ) {cout << "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoV error !" << map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"] << " " << HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ << endl;}
         cout << "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ " <<HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ  << endl;       

      }
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


}//end Loop function

void ssb_analysis::GetNtupleTotalEvent( unsigned int totevent )// Not Used Function
{
}

void ssb_analysis::Start()
{
   //fout = new TFile(Form("output/%s",outfile),"RECREATE");
   if (strcmp(outdir, "None") != 0 ) {
      //fout = new TFile(Form("%s/%s",outdir,outfile),"RECREATE");
      fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr/%s/%s",outdir,outfile),"RECREATE");
   }
   else {
      fout = new TFile(Form("output/%s",outfile),"RECREATE");
   }
   cout << "fout - getname : " << fout->GetName() << endl;
//   else if      ( genLoopon == 1 ){ fout = new TFile(Form("output/%s",outfile),"UPDATE");}
   //if      ( genLoopon == 0 ){ fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/SSB_CPviolation/output/%s",outfile),"RECREATE");}
   //else if ( genLoopon == 1 ){ fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/SSB_CPviolation/output/%s",outfile),"UPDATE"  );}
//   else {cout << "genLoopon error" << endl;}
   fout->cd("");

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();
}

void ssb_analysis::DeclareHistos()
{
}

void ssb_analysis::End()
{
   fout->Write();
   fout->Close();
}

void ssb_analysis::SetInputFileName( char *inname )
{
   FileName_ = inname;
   char unsco_ = '_';
   Size_t unscoIndex = FileName_.Last(unsco_);
   FileName_.Remove(unscoIndex, FileName_.Length());
   if (FileName_.Contains("/")){
      unscoIndex = FileName_.Last('/');
      FileName_.Remove(0,unscoIndex+1);
   }
   //cout << "FileName_ : " << FileName_ << endl;
   
}
void ssb_analysis::SetOutputFileName(char *outname, char *sedir)
{   
   outfile = outname;
   outdir = sedir;
}

/*void ssb_analysis::MCSF()
{
   if (FileName_.Contains("Data")||FileName_.Contains("Single")||FileName_.Contains("EG")){ mc_sf_ = 1.; return; }
   /// Open Xsec Tables ///
   FILE *xsecs_;
   char sampleName[1000];
   double xsec_ = -1.;
   double br_ = -1.; 
   int totalevt_  = -1.; 
   int positive_  = -1.; 
   int negative_  = -1.; 
   int posi_nega_ = -1.; 
   string xsec_dir= "./xsecAndsample/";
   string xsec_filePath = xsec_dir+ XsecTable_.Data();
   //cout << "xsec_filePath : " << xsec_filePath << endl;
   /// SampleName | TotalEvt | Positive+Negative | Xsection | Branching Fraction |
   xsecs_ = fopen(xsec_filePath.c_str(),"r");
   map<string, int> m_sam_totalevt;
   map<string, double> m_sam_xsec;
   map<string, double> m_sam_br;
   map<string, int> m_sam_positive;
   map<string, int> m_sam_negative;
   map<string, int> m_sam_posi_nega;
   if (xsecs_!=NULL) 
   { 
      //cout << "Load Xsection Table!" << endl;
      while (fscanf(xsecs_, "%s %d %d %d %d %lf %lf\n", sampleName, &totalevt_, &positive_, &negative_, &posi_nega_, &xsec_, &br_ ) != EOF)
      {
         //cout 
         //<< "sampleName : " << sampleName << " totalevt_ : " << totalevt_
         //<< " positive_ " << positive_ << " negative_ : " << negative_ 
         //<< " posi_nega_ " << posi_nega_ << " xsec_ : " << xsec_ 
         //<< " br_ " << br_ 
         //<< endl;
         m_sam_totalevt[sampleName] = totalevt_;
         m_sam_positive[sampleName] = positive_;
         m_sam_negative[sampleName] = negative_;
         m_sam_posi_nega[sampleName] = posi_nega_;
         m_sam_xsec[sampleName] = xsec_;
         m_sam_br[sampleName] = br_;
      }
      fclose(xsecs_);
   }
   else {cout << "No xsec_filePath !!!" << xsec_filePath << endl;return;} 
   //cout << "Lumi : " << Lumi << endl;
   double lumi = Lumi/1000000;
   auto it = m_sam_xsec.find(FileName_.Data());
   if (it !=  m_sam_xsec.end()){
      cout << "SK Key " << FileName_.Data() << " found in the map."<< endl;
      mc_sf_ = (m_sam_xsec[FileName_.Data()]*m_sam_br[FileName_.Data()]*lumi)/m_sam_posi_nega[FileName_.Data()];
   }
   else {
      mc_sf_ =1.;
      cout << "Key " << FileName_.Data() << " not found in the map. mc sf is 1" << mc_sf_ << endl;
   }
   return;
}*/

bool ssb_analysis::Trigger() {
    bool trigpass = false;

    // Define TString variables to determine decay mode and run type

    // Determine if the sample is Data or MC
    bool isData = FileName_.Contains("Data_");

    // Determine the run period if it's 2018
    TString runPeriod = "";
    if (FileName_.Contains("2018")) {
        runPeriod = "2018";
    }
    else if (FileName_.Contains("2016PreVFP")) {
        runPeriod = "2016PreVFP";
    }
    else if (FileName_.Contains("2016PostVFP")) {
        runPeriod = "2016PostVFP";
    }
    else if (FileName_.Contains("2017")) {
        runPeriod = "2017";
    }
    else {
       std::cout << "check out runperiod!!! in Trigger function ! FileName_ : " << FileName_ << std::endl;
    }

    // MC case: if any trigger in TriggerList passes, set trigpass = true
    if (!isData) {
       if (PassTrigger(isData, runPeriod)) { // Assuming PassTrigger is a method that checks a specific trigger
           trigpass = true;
       }
    }
    // Data case: further specify based on FileName_
    else {

    }

    return trigpass;

}
/*bool ssb_analysis::Trigger(Long64_t nb) {
    bool trigpass = false;

    // Define TString variables to determine decay mode and run type

    // Determine if the sample is Data or MC
    bool isData = FileName_.Contains("Data_");

    // Determine the run period if it's 2018
    int runPeriod = 0;
    if (FileName_.Contains("2018")) {
        runPeriod = 2018;
    }

    // Define Trigger Lists
    std::vector<TString> SingleLeptonTriggerList = {"HLT_IsoMu24"};
    std::vector<TString> DiLeptonTriggerList = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"};
    std::vector<TString> TriggerList = {"HLT_IsoMu24","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"}; // Add all triggers used in the analysis

    // MC case: if any trigger in TriggerList passes, set trigpass = true
    if (!isData) {
        for (const auto& trigger : TriggerList) {
            if (PassTrigger(trigger, nb)) { // Assuming PassTrigger is a method that checks a specific trigger
                trigpass = true;
                break;
            }
        }
    }
    // Data case: further specify based on FileName_
    else {
        if (runPeriod == 2018) {
            if (FileName_.Contains("Single")) {
                // If a single lepton trigger passes and no dilepton triggers pass, set trigpass = true
                bool singleLeptonPass = false;
                for (const auto& trigger : SingleLeptonTriggerList) {
                    if (PassTrigger(trigger, nb)) {
                        singleLeptonPass = true;
                        break;
                    }
                }
                
                bool dileptonPass = false;
                for (const auto& trigger : DiLeptonTriggerList) {
                    if (PassTrigger(trigger, nb)) {
                        dileptonPass = true;
                        break;
                    }
                }
                
                // Pass condition for SingleLepton triggers in 2018
                if (singleLeptonPass && !dileptonPass) {
                    trigpass = true;
                }
            }
            // Double case: if any trigger in TriggerList passes, set trigpass = true
            else if (FileName_.Contains("Double")) {
                for (const auto& trigger : TriggerList) {
                    if (PassTrigger(trigger, nb)) {
                        trigpass = true;
                        break;
                    }
                }
            }
        }
    }

    return trigpass;

}*/
/*bool ssb_analysis::PassTrigger(TString triggerName, Long64_t currentEntry) {
   if (!fChain->GetListOfBranches()->FindObject(triggerName)) {
      std::cerr << "Trigger branch : " << triggerName << " Not in Tree" << std::endl;
      return false;
   }
   else {
      bool triggerStatus = false;
        fChain->SetBranchAddress(triggerName, &triggerStatus);
        // 현재 이벤트에 대한 트리거 상태를 가져옴
        fChain->GetEntry(currentEntry);
        //cout << "triggerName : "  << triggerName  << " triggerStatus : " << triggerStatus << endl;
        return triggerStatus;
   }
}*/
//bool ssb_analysis::PassTrigger() {
bool ssb_analysis::PassTrigger(bool isdata_, TString runPeriod_) {
   //if (map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"] ){std::cout << "trigger : " << map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"] << std::endl;}
   if (!isdata_){
      /// For MC single + dilepton trigger 
      if (runPeriod_.Contains("2016PreVFP")){
 
         int selectrigpass;
         cout << "map_TrigPathName_passed " << map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"] << endl;
         selectrigpass = map_TrigPathName_passed["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"];  
/*         int selectrigpass = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
         + HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ
         + HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL
         + HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
         + HLT_IsoMu24
         + HLT_IsoTkMu24;*/
         if (selectrigpass > 0){
            std::cout << "selectrigpass : " << selectrigpass << std::endl;
         }
      }
   }
   return true; 
}

