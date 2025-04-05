#ifndef LEPTON_H
#define LEPTON_H

#include "Analysis.h"

class Lepton
{
public:
    Lepton();
    ~Lepton();

    static void LeptonSelector();
    static void LeptonOrder();
    static bool NumIsoLeptons();
    static bool LeptonsPtAddtional();
    static bool ThirdLeptonVeto(std::vector<int> vLep1Idx, std::vector<int> vLep2Idx = std::vector<int>());
    static bool DiLeptonMassCut();
    static bool ZMassVeto();
    // Lepton //    
    static double lep_pt;
    static double lep_eta;
    static double lepisocut;
    // index vector for object
    static std::vector<int> v_lepton_idx;
    static std::vector<int> v_muon_idx;
    static std::vector<int> v_electron_idx;
    // Muon //
    static TTreeReaderArray<Float_t>* muons_pt;
    static TTreeReaderArray<Float_t>* muons_eta;
    static TTreeReaderArray<Float_t>* muons_phi;
    static TTreeReaderArray<Float_t>* muons_M;
    static std::vector<TLorentzVector> muonsveto;
    static std::vector<TLorentzVector> muons;
    static double muonisocut;
    // Electron //
    static TTreeReaderArray<Float_t>* elecs_pt;
    static TTreeReaderArray<Float_t>* elecs_eta;
    static TTreeReaderArray<Float_t>* elecs_phi;
    static TTreeReaderArray<Float_t>* elecs_M;
    static std::vector<TLorentzVector> elecsveto;
    static std::vector<double> elecsveto_isoValue;
    static std::vector<TLorentzVector> elecs;
    static double elecisocut;
    // Z mass veto //    
    static bool zmassveto;
    static int eleid_scbcut;

private:
    static void MakeMuonCollection();
    static void MakeElecCollection();
    //static bool ThirdLeptonVeto(); // Lepton (muon & electron)
    //static void CorrectedMuonCollection(); // Muon  


    static std::vector<int> v_lep_idx_temp;
    static std::vector<int> v_muon_idx_temp;
    static std::vector<int> v_electron_idx_temp;

    
};

#endif // LEPTON_H
