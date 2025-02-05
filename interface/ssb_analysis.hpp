#ifndef _ssb_analysis_

#define _ssb_analysis_
  
#include <set>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraphErrors.h>
#include "TLorentzVector.h"
#include "TEnv.h"
  
//#include "./../analysis/NanoAODTree.h"
#include "./../analysis/2016/NanoAODTree.h"
#include "./../analysis/2017/NanoAODTree.h"
#include "./../analysis/Run3/NanoAODTree.h"

// Textreader
#include "./../TextReader/TextReader.hpp"
/*
#include "./../interface/ssb_eff.hpp"
#include "./../interface/ssb_cpviol.hpp"
#include "./../kinsol/TtFullLepKinSolver.hpp"
#include "./../interface/BTagCalibrationStandalone.h"
#include "./../interface/LumiReWeighting.h"

//#include "./../RoccoR/RoccoR.h"
#include "./../RoccoR/RoccoR.h"

#include "./../KinSolv/analysisUtils.h"
#include "./../KinSolv/KinematicReconstruction.h"
#include "./../KinSolv/KinematicReconstructionSolution.h"
*/
using namespace std;

class ssb_analysis : public NanoAODTree 
{
   public:
      //declare functions
      //ssb_analysis(TTree *tree=0);
      ssb_analysis(TTree *tree=0, string ConfName="");
      virtual ~ssb_analysis();

      //basic frame
      void SetNanoAODHeader();
      virtual void Loop( char *logfile );
      void GetNtupleTotalEvent( unsigned int totevent );
      void Start();
      void End();

      //user define functions
      void SetInputFileName( char *inname );
      void SetOutputFileName(char *outname, char *sedir);
      void DeclareHistos();
      void SetConfig();

      bool Trigger();
      //bool Trigger();
      //bool Trigger(Long64_t nb);
      //bool PassTrigger();
      bool PassTrigger(bool isdata_, TString runPeriod_);
      //bool PassTrigger(TString triggerName, Long64_t currentEntry);

   private:
      //put variables that you want
      char *outfile;
      char *outdir;
      TFile *fout;
      std::vector<TString>v_outName;
      unsigned int NtupletotalEvent;
      unsigned int totalEvent;

      double evt_weight_;
      TString Decaymode;
      TString FileName_; 

      TextReader *SSBConfReader;
      double mc_sf_;
      double pi;

      std::map<std::string, Bool_t> map_TrigPathName_passed;
      bool is2016;
      bool is2017;
      bool isRun3;

   public:

      //declare histograms

      TH1D *h_PileUp;
      TH1D *h_PileUp_Up;
      TH1D *h_PileUp_Down;
      TH1D *h_NumEl;
      TH1D *h_NumMu;
      TH1D *h_Num_PV_BeforePreSel;
      TH1D *h_Num_PV_AfterMetFilter;
      TH1D *h_Num_PV_AfterTrigger;
};
#endif

#ifdef ssb_analysis_cxx

ssb_analysis::ssb_analysis(TTree *tree, string confName)
{
   if (tree == 0)
   {
      printf("ERROR: Can't find any input tree.\n");
   }
   Init(tree);
    
   // Text Reader from Jaehoon.
   SSBConfReader = new TextReader();
   string confDir = "./configs/";
   string confpath = "";
   confpath = confDir+confName;
   SSBConfReader->ReadFile(confpath);
   SSBConfReader->ReadVariables();
   SSBConfReader->PrintoutVariables();

   Decaymode    = SSBConfReader->GetText( "Channel" ); // Channel

   pi = TMath::Pi();

}

ssb_analysis::~ssb_analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   delete fout;
   delete SSBConfReader;
}

#endif
   

