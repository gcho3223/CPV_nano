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

// Common tools 
#include "./../CommonTools.hpp"
// Textreader
#include "./../TextReader/TextReader.hpp"


class Analysis {
public:
    // Constructor and destructor
Analysis(TChain *chain, std::string inputName, std::string seDirName, std::string outputName, const std::string &branchListFile, const std::string &configFile, int NumEvt);

    //~Analysis() = default;
    ~Analysis();

    // Set variable function
    void SetVariables();
    // Event loop function
    void Loop();
    // create TLorentzVector //
    TLorentzVector createLorentzVector(float pt, float eta, float phi, float mass);
    void Start();
    void DeclareHistos();
private:
    // TTreeReader and TChain
    TChain *chain;
    std::string outfile;
    std::string outdir;
    TTreeReader fReader;
    TFile *fout;
    //TextReader from Jaehoon.
    TextReader *SSBConfReader;
    int NumEvt; //

    double Lumi;
    TString FileName_;
    TString RunPeriod;
    TString Decaymode;
    TString XsecTable_;

    TString cutflowName[10];

    int num_pv;
    // Maps for dynamic branch storage
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> boolSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Int_t>>> intSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UInt_t>>> uintSingles; // Add uintSingles for UInt_t types
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UChar_t>>> ucharSingles; // Add uintSingles for UInt_t types
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Float_t>>> floatSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>> floatVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Bool_t>>> boolVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Int_t>>> intVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UInt_t>>> uintVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UChar_t>>> ucharVectors;

    // Initialize TTreeReader with branches from the list
    void InitBranches(const std::string &branchListFile);
    std::string removeSubstring(std::string &str, const std::string &keyword);
    //std::unique_ptr<TTreeReaderValue<bool>> DeepCopy(const std::unique_ptr<TTreeReaderValue<bool>>& src);
    template <typename T>std::unique_ptr<TTreeReaderValue<T>> DeepCopy(const std::unique_ptr<TTreeReaderValue<T>>& src);
    bool METFilterAPP();
    bool Trigger();
    bool SelTrigger(std::vector<std::string> v_sel);
    //bool ;
    //TString SetInputFileName( char *inname );
    TString SetInputFileName( std::string inname );
    void SetObjectVariable();
    void MCSF();
    void MCSFApply();
    void NumPVCount();
    void LeptonSelector(); //(muon & electron)
    void LeptonOrder(); // Lepton (muon & electron)
    bool ThirdLeptonVeto(); // Lepton (muon & electron)
    //void CorrectedMuonCollection(); // Muon  
    void MakeMuonCollection(); // Muon  
    void MakeElecCollection(); // Electron 
    void JetSelector(); // Jet
    void MakeJetCollection(); // Jet 
    bool JetCleaning(TLorentzVector* jet_);
    TLorentzVector JERSmearing(TLorentzVector* jet, int idx_, TString op_);
    void METDefiner();

    //
    bool NumIsoLeptons();
    bool LeptonsPtAddtional();
    bool DiLeptonMassCut();

    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> triggerList;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> noiseFilters;
    std::vector<int> v_lep_idx_temp;
    std::vector<int> v_muon_idx_temp;
    std::vector<int> v_electron_idx_temp;

    TLorentzVector Lep1; // 
    TLorentzVector Lep2; //
    TLorentzVector Lep; // 
    TLorentzVector AnLep; //
    TLorentzVector Met; // 

    // 
    TTreeReaderArray<float>* leptons_iso;
    TTreeReaderArray<bool>   *leptons_Id;

    // muons //
    TTreeReaderArray<Float_t>* muons_pt;
    TTreeReaderArray<Float_t>* muons_eta;
    TTreeReaderArray<Float_t>* muons_phi;
    TTreeReaderArray<Float_t>* muons_M;
    TTreeReaderArray<float>* muons_iso;
    TTreeReaderArray<float>* muonsveto_iso;
    TTreeReaderArray<bool>   *muons_Id;
    TTreeReaderArray<bool>   *muonsveto_Id;
    std::vector<TLorentzVector> muons;
    std::vector<TLorentzVector> muonsveto;
    TLorentzVector TMuon;
    TLorentzVector TElectron;

    // electrons //
    TTreeReaderArray<Float_t>* elecs_pt;
    TTreeReaderArray<Float_t>* elecs_eta;
    TTreeReaderArray<Float_t>* elecs_phi;
    TTreeReaderArray<Float_t>* elecs_M;
    TTreeReaderArray<float>*   elecs_iso;
    TTreeReaderArray<float>*   elecsveto_iso;
    TTreeReaderArray<Int_t>*   elecs_scbId;
    TTreeReaderArray<Bool_t>*  elecs_mvaId;
    TTreeReaderArray<Int_t>*   elecsveto_scbId;
    TTreeReaderArray<Bool_t>*  elecsveto_mvaId;
    std::vector<TLorentzVector> elecs;
    std::vector<TLorentzVector> elecsveto;

    // jets //
    TTreeReaderArray<Float_t>* jets_pt;
    TTreeReaderArray<Float_t>* jets_eta;
    TTreeReaderArray<Float_t>* jets_phi;
    TTreeReaderArray<Float_t>* jets_M;
    TTreeReaderArray<Bool_t>*  jets_Id;
    TTreeReaderArray<Float_t>* jets_btag;
    std::vector<TLorentzVector> jets;
    bool dojer;

    // MET //
    TTreeReaderValue<Float_t>* met_pt;
    TTreeReaderValue<Float_t>* met_phi;

    //LepIdType   v_muon_Id;
    //LepIdType   v_elec_Id;

    // index vector for object (muon, electron, jet)
    std::vector<int> v_lepton_idx; // Indecies of lepton
    std::vector<int> v_muon_idx; // Indecies of muon
    std::vector<int> v_electron_idx; // Indecies of electron
    std::vector<int> v_jet_idx;  // Indecies of jet

    //JetCleanning information from config.
    TString muiso_type;
    TString muId;
    TString eliso_type;
    TString elId;

    double pi;

    //to get trigger information from config.
    int num_dleptrig;
    int num_sleptrig;
    std::vector<std::string> DLtrigName; // trigger 
    std::vector<std::string> SLtrigName; // trigger 
    std::vector<std::string> trigName; // trigger 

    // lepton variables
    TString lepisotype;
    TString lepId;

    //Muon Type for MuEl - channel.
    TString MuonIsoType;
    TString MuonId;
    TString ElecId;
    TString ElecIsoType;
    TString JetId;
    TString JetbTag;
    TString JetEnSys;
    TString JetResSys;
    TString BTagSFSys;
    TString BTagEffSys;
    // MET Systematic ...
    TString MetSys;
    
    // Lepton Systematic
    TString LepIdSFSys;
    TString LepIsoSFSys;
    TString LepRecoSFSys;
    TString LepTrackSFSys;
    
    // PileUp Systematic ...
    std::string PileUpMCFile;
    std::string PileUpDATAFile;
    TString PileUpSys;
    TString L1PreFireSys;

    // PDF Systematic ... Getting the Num of PDF Set
    int PDFSys;
    // FactReno Systematic ... Getting the Num of FactReno Set
    int FactRenoSys;
    // Fragment Systematic ... Getting the Num of Fragment Set
    TString FragmentSys;
    // DecayTable Systematic ... Getting the Num of DecayTable Set
    TString DecayTableSys;
    // TopPtReweight Sys
    TString TopPtSys;


    int eleid_scbcut;
    int elevetoid_scbcut;

    /// Event weights ...
    double mc_sf_;
    double evt_weight_; 
    double evt_weight_beforemcsf_;
    double evt_weight_beforegenweight_;
    double evt_weight_beforefactrenoweight_;
    double evt_weight_beforedectabweight_;
    double evt_weight_beforefragmentweight_;
    double evt_weight_beforepdfweight_;
    double evt_weight_beforePileup_;
    double evt_weight_beforeL1PreFire_;
    double evt_weight_beforeTrigger_;
    double evt_weight_beforeLepsf_;
    double evt_weight_beforeBtag_;
    double lep_eff;

    //Kinetic Variables                                                                                                               
    double lep_pt;
    double lep_eta;
    double lepisocut;                                                                                                                 
    int    n_elep_Id;                                                                                                                 
    double el_id_1;                                                                                                                   
    double el_id_2;                                                                                                                   
    
    double muon_pt;
    double muon_eta;
    double muonisocut;                                                                                                                
    double elec_pt;
    double elec_eta;
    double elecisocut;
    int    n_elec_Id;


    double jet_pt;
    double jet_eta;
    int    jet_id;
    float  bdisccut;
    double met_cut;

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


    TH1D *h_Lep1pt[10];
    TH1D *h_Lep2pt[10];
    TH1D *h_Lep1eta[10];
    TH1D *h_Lep2eta[10];
    TH1D *h_Lep1phi[10];
    TH1D *h_Lep2phi[10];
    TH1D *h_Jet1pt[10];
    TH1D *h_Jet2pt[10];
    TH1D *h_Jet1eta[10];
    TH1D *h_Jet2eta[10];
    TH1D *h_Jet1phi[10];
    TH1D *h_Jet2phi[10];
    TH1D *h_HT[10]; 
    TH1D *h_METpt[10]; 
    TH1D *h_METphi[10]; 
    TH1D *h_DiLepMass[10]; 
    TH1D *h_Num_PV[10];
    TH1D *h_Num_Jets[10];
    TH1D *h_Num_bJets[10];

    TH1D *h_Reco_CPO_[13];
    TH1D *h_Reco_CPO_ReRange_[13];
};

#endif // ANALYSIS_H

