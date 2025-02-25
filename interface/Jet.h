#ifndef JET_H
#define JET_H

#include "Analysis.h"
#include "Lepton.h"

class Jet {
public:
    Jet();
    ~Jet();

    static void JetSelector();
    static void JetDefiner();
    static bool nJet(std::vector<int> v_jet);

    //Muon Type for MuEl - channel.
    static TString JetEnSys;
    static TString JetResSys;

    // Jet //
    static TTreeReaderArray<Float_t>* jets_pt;
    static TTreeReaderArray<Float_t>* jets_eta;
    static TTreeReaderArray<Float_t>* jets_phi;
    static TTreeReaderArray<Float_t>* jets_M;
    static TTreeReaderArray<Float_t>* jets_btag;
    static TTreeReaderArray<Int_t>*  jets_id;
    static TTreeReaderArray<Int_t>*  jets_puid;
    static std::vector<double> v_jetdpt_res;
    static std::vector<double> v_jetdpx_res;
    static std::vector<double> v_jetdpy_res;
    
    static bool dojer;
    
    // for debugging Jet cleaning and selection
    static Int_t injet;
    static Int_t cleanedjet;
    static Int_t passedjet;

    static bool njet;
private:
    static void MakeJetCollection();
    static bool JetCleaning(TLorentzVector* jet_);
    static TLorentzVector JERSmearing(TLorentzVector* jet, int idx_, TString op_);
};

#endif // JET_H
