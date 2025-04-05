#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector2.h>
#include <TH1D.h>
#include <TFile.h>

#include "./../CommonTools.hpp"
#include "./../TextReader/TextReader.hpp"


class Analysis {
public:
    /// Constructor and destructor ///
    Analysis(TChain *chain, std::string inputName, std::string seDirName, std::string outputName, const std::string &branchListFile, const std::string &configFile, int NumEvt);
    ~Analysis();
    void InitBranches(const std::string &branchListFile);
    void Loop();
    void Start();
    void DeclareHistos();
    // void SetVariables();
    // TLorentzVector createLorentzVector(float pt, float eta, float phi, float mass);

    static TextReader *SSBConfReader; //analysis, setvariables
    static TString FileName_; //trigger, MCSF, MCSFApply
    static TString RunPeriod; //setvariables, setobjectvariable, trigger
    static TString Decaymode; //setvariables, setobjectvariable, trigger, leptonselector, leptonorder, numisoleptons, leptonsptadditional, dileptonmasscut
    static TString XsecTable_; //setvariables, MCSF

    /// Maps for dynamic branch storage ///
    static std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> boolSingles; //initbranches, setvariables
    static std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UInt_t>>> uintSingles; //initbranches // Add uintSingles for UInt_t types
    static std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Float_t>>> floatSingles; //initbranches, setvariables
    static std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>> floatVectors; //initbranches, setobjectvariable, leptonselector
    static std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Bool_t>>> boolVectors; //initbranches, setvariables, setobjectvariable
    static std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Int_t>>> intVectors; //initbranches, setvariables, setobjectvariable, leptonselector, leptonorder
    static std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> noiseFilters; //SetVariables, METFilterAPP

    static double Lumi; //MCSF, setvariables

    static TLorentzVector Lep1;
    static TLorentzVector Lep2;
    static TLorentzVector Lep;
    static TLorentzVector AnLep;
    static TLorentzVector Met;
    static TLorentzVector Jet1;
    static TLorentzVector Jet2;
    static TLorentzVector bJet1;
    static TLorentzVector bJet2;

private:
    TString SetInputFileName( std::string inname );
    void NumPVCount();

    ///////////////
    // variables //
    ///////////////
    /// TTreeReader and TChain ///
    TChain *chain;
    TTreeReader fReader; 
    std::string outfile;
    std::string outdir;
    TFile *fout;
    
    TString cutflowName[10];

    /// Maps for dynamic branch storage ///
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Int_t>>> intSingles; //initbranches, loop
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UChar_t>>> ucharSingles; //initbranches // Add uintSingles for UInt_t types
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UInt_t>>> uintVectors; //initbranches
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UChar_t>>> ucharVectors; //initbranches

    int NumEvt; //loop
    int num_pv; //loop
    int maxEvt;
    // muons //
    TLorentzVector TMuon;
    TLorentzVector TElectron;

    // electrons //

    // jets //
    TTreeReaderArray<Float_t>* jets_btag;

    //LepIdType   v_muon_Id;
    //LepIdType   v_elec_Id;

    double pi;

    //Kinetic Variables                                                                                                               
    int    n_elep_Id;
    double el_id_1;
    double el_id_2;
    
    int    n_elec_Id;

    //Jetcleannig Variables
    double mu_pt_jetcl;
    double mu_eta_jetcl;
    double mu_iso_jetcl;

    double el_pt_jetcl;
    double el_eta_jetcl;
    double el_iso_jetcl;
    int    n_el_id_jetcl;
    double el_id_jetcl_1;
    double el_id_jetcl_2;

    // DeclareHistos //
    TH1D *h_Num_PV_beforePreSel;
    TH1D *h_Num_PV_afterMETFilter;
    TH1D *h_Num_PV_afterTrigger;
    TH1D *h_passedTriggers;
    TH2D *h_triggerEvent;
    // debug: selected muon information after object selection //
    TH1D *h_AllSelMu_multiplisity_wW;
    TH1D *h_AllSelMu_multiplisity_woW;
    TH1D *h_AllSelMu_pt;
    TH1D *h_AllSelMu_eta;
    TH1D *h_AllSelMu_phi;
    TH1D *h_AllSelMu_iso;
    // debug: all selected electrons
    TH1D *h_AllSelVetoEle_multiplisity_wW;
    TH1D *h_AllSelVetoEle_multiplisity_woW;
    TH1D *h_AllSelVetoEle_pt;
    TH1D *h_AllSelVetoEle_eta;
    TH1D *h_AllSelVetoEle_phi;
    TH1D *h_AllSelVetoEle_iso;

    TH1D *h_Num_PV[10];
    
    TH1D *h_Lep1pt[10];
    TH1D *h_Lep2pt[10];
    TH1D *h_Lep1eta[10];
    TH1D *h_Lep2eta[10];
    TH1D *h_Lep1phi[10];
    TH1D *h_Lep2phi[10];
    TH1D *h_DiLepMass[10];

    TH1D *h_Jet1pt[10];
    TH1D *h_Jet2pt[10];
    TH1D *h_Jet1eta[10];
    TH1D *h_Jet2eta[10];
    TH1D *h_Jet1phi[10];
    TH1D *h_Jet2phi[10];
    TH1D *h_bJet1pt[10];
    TH1D *h_bJet2pt[10];
    TH1D *h_bJet1eta[10];
    TH1D *h_bJet2eta[10];
    TH1D *h_bJet1phi[10];
    TH1D *h_bJet2phi[10];
    TH1D *h_Num_Jets[10];
    TH1D *h_Num_Jets_puid[10];
    TH1D *h_Num_bJets[10];

    TH1D *h_HT[10]; 

    TH1D *h_METpt[10]; 
    TH1D *h_METphi[10]; 

    /// Observables for CPV ///
    TH1D *h_Reco_CPO_[13];
    TH1D *h_Reco_CPO_ReRange_[13];
};

#endif // ANALYSIS_H

