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

    static TString lepId;
    // Lepton for dimuon & dielectron channel//
    static TTreeReaderArray<float>* leptons_iso;
    static TTreeReaderArray<bool>*  leptons_Id;
    
    // Muon for dimuon & emu channel //
    static TTreeReaderArray<bool>*  muons_Id;
    static TTreeReaderArray<float>* muons_iso;
    static double muon_pt;
    static double muon_eta;
    // Electron for dielectron and emu channel //
    static TTreeReaderArray<Bool_t>*  elecs_mvaId;
    static TTreeReaderArray<float>*   elecs_iso;
    static TTreeReaderArray<Int_t>*   elecs_scbId;
    static double elec_pt;
    static double elec_eta;
    // veto muon for jet cleaning and 3rd lepton veto //
    // veto muon //
    static TTreeReaderArray<bool>*  vetomuon_id;
    static TTreeReaderArray<float>* vetomuon_iso;
    static float vetomu_iso_cut;
    static float vetomu_pt;
    static float vetomu_eta;
    // electron veto //
    static TTreeReaderArray<Bool_t>*  vetoelec_mvaId;
    static TTreeReaderArray<float>*   vetoelec_iso;
    static TTreeReaderArray<Int_t>*   vetoelec_scbId;
    static float vetoelec_iso_cut;
    static float vetoelec_pt;
    static float vetoelec_eta;
    //Jet                                                                                         
    static int    jet_id;
    static int    jet_puid;
    static double jet_pt;
    static double jet_eta;
    // Electron veto?
    static int elevetoid_scbcut;
    
private:
    static std::string removeSubstring(std::string &str, const std::string &keyword);
    template <typename T>static std::unique_ptr<TTreeReaderValue<T>> DeepCopy(const std::unique_ptr<TTreeReaderValue<T>>& src);
    
    ///// Read informations from config. /////
    //to get trigger information from config.
    static int num_dleptrig;
    static int num_sleptrig;
    // lepton variables
    static TString lepisotype;
    // Muon Type for MuEl - channel.
    static TString MuonIsoType;
    static TString MuonId;
    static TString ElecIsoType;
    static TString ElecId;
    // Jet
    static TString JetId;
    static TString JetPUID;
    static TString JetbTag;
    static TString BTagSFSys;
    static TString BTagEffSys;
    //JetCleanning information from config.
    static TString vetomu_id;
    static TString vetomu_isoid;
    static TString vetoelec_id;
    static TString vetoelec_isoid;
    // MET Systematic ...
    static TString MetSys;
    // Lepton Systematic
    static TString LepIdSFSys;
    static TString LepIsoSFSys;
    static TString LepRecoSFSys;
    static TString LepTrackSFSys;
    // SF & weight & Systematics
    // pile up
    static TString PileUpSys;
    static std::string PileUpMCFile;
    static std::string PileUpDATAFile;
    // L1PreFire
    static TString L1PreFireSys;
    // PDF Systematic ... Getting the Num of PDF Set
    static int PDFSys;
    // FactReno Systematic ... Getting the Num of FactReno Set
    static int FactRenoSys;
    // Fragment
    static TString FragmentSys;
    // DecayTable
    static TString DecayTableSys;
    // TopPtReweight
    static TString TopPtSys;
};

#endif // CONFIG_H