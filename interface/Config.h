#ifndef CONFIG_H
#define CONFIG_H

#include "Analysis.h"

class Config {
public:
    Config();
    ~Config();

    static void SetVariables();
    static void SetObjectVariable();
    static TLorentzVector createLorentzVector(float pt, float eta, float phi, float mass);
    
    // Lepton //
    static TTreeReaderArray<float>* leptons_iso;
    static TTreeReaderArray<bool>*  leptons_Id;
    // Muon //
    static TTreeReaderArray<bool>*  muonsveto_Id;
    static TTreeReaderArray<float>* muons_iso;
    static TTreeReaderArray<float>* muonsveto_iso;
    static TTreeReaderArray<bool>*  muons_Id;
    static double muon_pt;
    static double muon_eta;
    // Electron //
    static TTreeReaderArray<Bool_t>*  elecs_mvaId;
    static TTreeReaderArray<Bool_t>*  elecsveto_mvaId;
    static TTreeReaderArray<float>*   elecs_iso;
    static TTreeReaderArray<float>*   elecsveto_iso;
    static TTreeReaderArray<Int_t>*   elecs_scbId;
    static TTreeReaderArray<Int_t>*   elecsveto_scbId;
    static double elec_pt;
    static double elec_eta;
    //Muon Type for MuEl - channel.
    static TString MuonIsoType;
    static TString MuonId;
    static TString ElecIsoType;
    static TString ElecId;
    //Jet                                                                                         
    static TString JetId;
    static TString JetPUID;
    static TString JetbTag;
    static TString BTagSFSys;
    static TString BTagEffSys;
    static int    jet_id;
    static int    jet_puid;
    static double jet_pt;
    static double jet_eta;

private:
    static std::string removeSubstring(std::string &str, const std::string &keyword);
    template <typename T>static std::unique_ptr<TTreeReaderValue<T>> DeepCopy(const std::unique_ptr<TTreeReaderValue<T>>& src);
    
    //JetCleanning information from config.
    static TString muiso_type;
    static TString muId;
    static TString eliso_type;
    static TString elId;
    
    //to get trigger information from config.
    static int num_dleptrig;
    static int num_sleptrig;
    
    // lepton variables
    static TString lepisotype;
    static TString lepId;
    
    // MET Systematic ...
    static TString MetSys;

    // Lepton Systematic
    static TString LepIdSFSys;
    static TString LepIsoSFSys;
    static TString LepRecoSFSys;
    static TString LepTrackSFSys;

    // PileUp Systematic ...
    static std::string PileUpMCFile;
    static std::string PileUpDATAFile;
    static TString PileUpSys;
    static TString L1PreFireSys;

    // PDF Systematic ... Getting the Num of PDF Set
    static int PDFSys;
    // FactReno Systematic ... Getting the Num of FactReno Set
    static int FactRenoSys;
    // Fragment Systematic ... Getting the Num of Fragment Set
    static TString FragmentSys;
    // DecayTable Systematic ... Getting the Num of DecayTable Set
    static TString DecayTableSys;
    // TopPtReweight Sys
    static TString TopPtSys;

    static int elevetoid_scbcut;

};

#endif // CONFIG_H