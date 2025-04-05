#include "../interface/Lepton.h"
#include "../interface/Config.h"
//-------------------------------Initialize variables-------------------------------//
/// Lepton ///
std::vector<int> Lepton::v_lep_idx_temp;
std::vector<int> Lepton::v_lepton_idx;
double Lepton::lep_pt;
double Lepton::lep_eta;
double Lepton::lepisocut;
/// Muon ////
std::vector<int> Lepton::v_muon_idx_temp;
std::vector<int> Lepton::v_muon_idx;
std::vector<TLorentzVector> Lepton::muons;
TTreeReaderArray<Float_t>* Lepton::muons_pt;
TTreeReaderArray<Float_t>* Lepton::muons_eta;
TTreeReaderArray<Float_t>* Lepton::muons_phi;
TTreeReaderArray<Float_t>* Lepton::muons_M;
double Lepton::muonisocut;

std::vector<TLorentzVector> Lepton::muonsveto;

/// Electron ///
std::vector<int> Lepton::v_electron_idx_temp;
std::vector<int> Lepton::v_electron_idx;
std::vector<TLorentzVector> Lepton::elecs;
TTreeReaderArray<Float_t>* Lepton::elecs_pt;
TTreeReaderArray<Float_t>* Lepton::elecs_eta;
TTreeReaderArray<Float_t>* Lepton::elecs_phi;
TTreeReaderArray<Float_t>* Lepton::elecs_M;
std::vector<TLorentzVector> Lepton::elecsveto;
int Lepton::eleid_scbcut;
std::vector<double> Lepton::elecsveto_isoValue;
double Lepton::elecisocut;
/// Z mass veto ///
bool Lepton::zmassveto;

namespace {
bool passIsolation(float iso, float isoCut) {return iso <= isoCut;}
bool passKinematicCuts(float _pt, float _eta, float ptCut, float etaCut) {return _pt >= ptCut && fabs(_eta) <= etaCut;}
bool passId(bool id) {return id;}
bool elecSCBId(int id, int idcut) {return id >= idcut;}
bool elecCharge(int id) {return id > 0;}
}
void Lepton::LeptonSelector()
{
    //std::cout << "||------------------------Lepton Selector------------------------||" << std::endl;

    // Clear output vectors
    v_lep_idx_temp.clear();
    v_lepton_idx.clear();
    v_muon_idx_temp.clear();
    v_muon_idx.clear();
    v_electron_idx.clear();
    v_electron_idx_temp.clear();
    elecsveto.clear();
    elecsveto_isoValue.clear();
    muonsveto.clear();

    // Prepare corrected collections
    MakeMuonCollection();
    MakeElecCollection();

    // Lambda functions for common checks
    //auto passIsolation = [](float iso, float isoCut)
    //{return iso <= isoCut;};

    //auto passKinematicCuts = [](float _pt, float _eta, float ptCut, float etaCut)
    //{return _pt >= ptCut && fabs(_eta) <= etaCut;};

    //auto passId = [](bool id) {return id;};

    //auto elecSCBId = [](int id, int idcut) {return id >= idcut;};
    //auto elecCharge = [](int id) {return id > 0;};

    /////////////////
    // Dimuon case //
    /////////////////
    if(TString(Analysis::Decaymode).Contains("dimuon"))
    {
        Int_t nmu = Config::muons_iso->GetSize();
        for(int i = 0; i < nmu; ++i)
        {
            //std::cout << "nmu[" << i << "]: " << "muons_iso: " << Config::muons_iso->At(i) << ", lepisocut: " << lepisocut << " | muons_pt: " << muons_pt->At(i) << ", lep_pt: " << lep_pt << " | muons_eta: " << muons_eta->At(i) << ", lep_eta: " << lep_eta << " | muons_Id: " << Config::muons_Id->At(i) << std::endl;
            if(!passIsolation(Config::muons_iso->At(i), lepisocut) || 
                !passKinematicCuts(muons_pt->At(i), muons_eta->At(i), lep_pt, lep_eta) || 
                !passId(Config::muons_Id->At(i)))
            {
                //std::cout << "muon[" << i << "] is not passed" << std::endl;
                continue;
            }
            v_lep_idx_temp.push_back(i);
            //v_lep_idx_temp.push_back(i); // I think this is copy mistake

            if(v_lep_idx_temp.size() == 1 && v_lepton_idx.empty())
            {v_lepton_idx.push_back(i);}
            else if(v_lep_idx_temp.size() > 0 && (*Analysis::intVectors["Muon_charge"])[v_lep_idx_temp[0]] != (*Analysis::intVectors["Muon_charge"])[i]) // push_back opposite charged muon
            {v_lepton_idx.push_back(i);}
        }
        // Assign final indices to v_muon_idx
        v_muon_idx = v_lepton_idx;
        //std::cout << "|| Dimuon case: size of muon index(v_muon_idx.size()) : " << v_muon_idx.size() << std::endl;
        //for(int i : v_muon_idx)
        //{std::cout << "v_muon_idx[" << i << "]: " << muons_pt->At(i) << ", " << muons_eta->At(i) << ", " << muons_phi->At(i) << ", " << muons_M->At(i) << std::endl;}
    }
    /////////////////////
    // Dielectron case //
    /////////////////////
    else if(TString(Analysis::Decaymode).Contains("dielec"))
    {
        //Int_t nel = Elec->GetEntriesFast();
        Int_t nel = elecs_pt->GetSize();
        for(int i = 0; i < nel; ++i)
        {
            //std::cout << "nel[" << i << "]: " << " | elecs_pt: " << elecs_pt->At(i) << ", lep_pt: " << lep_pt << " | elecs_eta: " << elecs_eta->At(i) << ", lep_eta: " << lep_eta << " | elecs_Id: " << Config::elecs_scbId->At(eleid_scbcut: " << eleid_scbcut << std::endl;
            //std::cout << "Electron_deltaEtaSC: " << (*Analysis::floatVectors["Electron_deltaEtaSC"])[i] << "fabs: " << fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) << "| elecCharge: " << elecCharge((*Analysis::intVectors["Electron_tightCharge"])[i]) << std::endl;
            //Ele_lep_sel = elecs.at(i);
            if(!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), lep_pt, lep_eta) ||
                !elecSCBId(Config::elecs_scbId->At(i), eleid_scbcut) ||
                (fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) > 1.4442 &&
                 fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) < 1.566) ||
                !(elecCharge((*Analysis::intVectors["Electron_tightCharge"])[i])))
            {continue;}
            v_lep_idx_temp.push_back(i);
            
            if(v_lep_idx_temp.size() == 1 && v_lepton_idx.empty())
            {
                v_lepton_idx.push_back(i);
                //v_electron_idx.push_back(i); // I think this is copy mistake
            }
            else if(v_lep_idx_temp.size() > 0 &&
                (*Analysis::intVectors["Electron_charge"])[v_lep_idx_temp[0]] != (*Analysis::intVectors["Electron_charge"])[i]) // push_back opposite charged electron
            {
                v_lepton_idx.push_back(i);
                //v_electron_idx.push_back(i); // I think this is copy mistake
            }
        }
        // Assign final indices to v_electron_idx
        v_electron_idx = v_lepton_idx;
        //std::cout << "|| Dielectron case: size of electron index(v_electron_idx.size()) : " << v_electron_idx.size() << std::endl;
    }
    ///////////////
    // MuEl case //
    ///////////////
    else if (TString(Analysis::Decaymode).Contains("muel"))
    {
        ///////////////////
        // Process muons //
        ///////////////////
        //Int_t nmu = Muon->GetEntriesFast();
        Int_t nmu = Config::muons_iso->GetSize();
        for (int i = 0; i < nmu; ++i)
        {
            if (!passIsolation(Config::muons_iso->At(i), muonisocut) || 
                //!passKinematicCuts(*Mu_lep_sel, muon_pt, muon_eta) || 
                //!passKinematicCuts(muons_pt->At(i), muons_eta->At(i), lep_pt, lep_eta) || 
                !passKinematicCuts(muons_pt->At(i), muons_eta->At(i), Config::muon_pt, Config::muon_eta) || 
                !passId(Config::muons_Id->At(i)))
            {continue;}
            v_muon_idx.push_back(i);
        }
        ///////////////////////
        // Process electrons //
        ///////////////////////
        Int_t nel = elecs_pt->GetSize();
        for (int i = 0; i < nel; ++i)
        {
          //if (!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), lep_pt, lep_eta) || 
          if (!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), Config::elec_pt, Config::elec_eta) || 
          !elecSCBId(Config::elecs_scbId->At(i), eleid_scbcut) || 
          (fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) > 1.4442 && 
          fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) < 1.566) || 
          !(elecCharge((*Analysis::intVectors["Electron_tightCharge"])[i])) || // Added missing logical operator
          !(*Analysis::boolVectors["Electron_convVeto"])[i])
          {continue;}
          v_electron_idx_temp.push_back(i);
        }
        /// muon and electron charge check ///
        if(!v_muon_idx.empty() && !v_electron_idx_temp.empty())
        {
            for(int i : v_electron_idx_temp)
            {
                // Added missing closing parenthesis
                if ((*Analysis::intVectors["Muon_charge"])[v_muon_idx[0]] != (*Analysis::intVectors["Electron_charge"])[i]) // I think this is correct
                //if ((*Analysis::intVectors["Muon_charge"])[v_lep_idx_temp[0]] != (*Analysis::intVectors["Electron_charge"])[i])
                {v_electron_idx.push_back(i);}
            }
        }
        //////////////////////////////
        // Validate charge matching //
        //////////////////////////////
        if(!v_muon_idx.empty() && !v_electron_idx.empty() &&
          (*Analysis::intVectors["Muon_charge"])[v_muon_idx[0]] == (*Analysis::intVectors["Electron_charge"])[v_electron_idx[0]])
        {
           if((*Analysis::intVectors["Muon_charge"])[v_muon_idx_temp[0]] != (*Analysis::intVectors["Electron_charge"])[v_electron_idx[0]]) // I think this is not needed
           {std::cerr << "MuEl selection wrong!" << std::endl;}
        }
    }
    else {std::cerr << "Lepton Selection error" << std::endl;}

    /////////////////////////////
    // Muon Selection for veto //
    /////////////////////////////
    //std::cout << "||--Muon/Electron for veto--                              ||" << std::endl;
    Int_t nmu = Config::vetomuon_iso->GetSize();

    for(int imu = 0; imu < nmu; ++imu)
    {
        //std::cout << "nmu2[" << imu << "]: " << "muons_iso: " << Config::muons_iso->At(imu) << ", lepisocut: " << lepisocut << " | muons_pt: " << muons_pt->At(imu) << ", lep_pt: " << lep_pt << " | muons_eta: " << muons_eta->At(imu) << ", lep_eta: " << lep_eta << " | muons_Id: " << Config::muons_Id->At(imu) << std::endl;
        //if(!passIsolation(Config::muons_iso->At(imu), lepisocut) || 
        if(!passIsolation(Config::muons_iso->At(imu), Config::vetomu_iso_cut) || 
           //!passKinematicCuts(muons_pt->At(imu), muons_eta->At(imu), lep_pt, lep_eta) ||
           !passKinematicCuts(muons_pt->At(imu), muons_eta->At(imu), Config::vetomu_pt, Config::vetomu_eta) ||
           !passId(Config::muons_Id->At(imu)))
        {/*std::cout << "muon[" << imu << "] is not passed" << std::endl;*/ continue;}
        muonsveto.push_back(muons.at(imu));
    }
    //std::cout << "|| " << muonsveto.size() << " muons are seclected" << std::endl;
    /////////////////////////////////
    // Electron Selection for veto //
    /////////////////////////////////
    if(Config::vetoelec_iso == nullptr)
    {
        //std::cout << "Error: Config::vetoelec_iso is null!" << std::endl;
        return;
    }
    //std::cout << "Debug! Electron seclection: elecs_iso: " << Config::vetoelec_iso->GetSize() << ", elecs.size: " << elecs.size() << std::endl;
    Int_t nel = Config::vetoelec_iso->GetSize();
    if(nel != elecs.size())
    {
        std::cerr << "Error: Size mismatch! nel: " << nel << ", elecs.size: " << elecs.size() << std::endl;
        return;
    }
    //if(!elecs_pt || !elecs_eta || !Analysis::floatVectors["Electron_deltaEtaSC"] || !Analysis::intVectors["Electron_tightCharge"])
    if(!elecs_pt || !elecs_eta || !Analysis::floatVectors["Electron_deltaEtaSC"] || !Analysis::intVectors["Electron_tightCharge"])
    {
        if(!elecs_pt) std::cout << "Error: elecs_pt is null!" << std::endl;
        if(!elecs_eta) std::cout << "Error: elecs_eta is null!" << std::endl;
        if(!Analysis::floatVectors["Electron_deltaEtaSC"]) std::cout << "Error: floatVectors['Electron_deltaEtaSC'] is null!" << std::endl;
        if(!Analysis::intVectors["Electron_tightCharge"]) std::cout << "Error: intVectors['Electron_tightCharge'] is null!" << std::endl;
        return;
    }
    for(int iel = 0; iel < nel; ++iel)
    {
        //std::cout << "|| Processing electron index: " << iel << std::endl;
        // 벡터 및 크기 확인
        //if(!passIsolation(Config::vetoelec_iso->At(iel), Config::vetoelec_iso_cut) || !elecs_pt || !elecs_eta || !Config::vetoelec_scbId || !Analysis::floatVectors["Electron_deltaEtaSC"] || !Analysis::intVectors["Electron_tightCharge"])
        if(!passIsolation(Config::vetoelec_iso->At(iel), Config::vetoelec_iso_cut) || !passKinematicCuts(elecs_pt->At(iel), elecs_eta->At(iel), Config::vetoelec_pt, Config::vetoelec_eta) || !Config::vetoelec_scbId || !Analysis::floatVectors["Electron_deltaEtaSC"] || !Analysis::intVectors["Electron_tightCharge"])
        {
            if(!passIsolation(Config::vetoelec_iso->At(iel), Config::vetoelec_iso_cut)) {/*std::cout << "Error: elecs_iso is not passed!: vetoelec_iso->At(iel): " << Config::vetoelec_iso->At(iel) << ", vetoelec_iso_cut: " << Config::vetoelec_iso_cut << std::endl;*/ return;}
            //if(!elecs_pt) {std::cerr << "Error: elecs_pt is nullptr!: elecs_pt: " << elecs_pt << std::endl; return;}
            //if(!elecs_eta) {std::cerr << "Error: elecs_eta is nullptr!: elecs_eta: " << elecs_eta << std::endl; return;}
            if(!passKinematicCuts(elecs_pt->At(iel), elecs_eta->At(iel), Config::vetoelec_pt, Config::vetoelec_eta)) {/*std::cout << "Error: elecs_pt&eta is not passed!: vetoelec_pt->At(iel): " << Config::vetoelec_pt->At(iel) << ", vetoelec_eta->At(iel): " << Config::vetoelec_eta->At(iel) << std::endl;*/ return;}
            if(!Config::vetoelec_scbId) {std::cerr << "Error: vetoelec_scbId is nullptr!" << std::endl; return;}
            if(!Analysis::floatVectors["Electron_deltaEtaSC"]) {std::cerr << "Error: floatVectors[\"Electron_deltaEtaSC\"] is nullptr!" << std::endl; return;}
            if(!Analysis::intVectors["Electron_tightCharge"]) {std::cerr << "Error: intVectors[\"Electron_tightCharge\"] is nullptr!" << std::endl; return;}
            return;
        }
        if(elecs_pt->GetSize() <= iel || elecs_eta->GetSize() <= iel ||
            Config::vetoelec_scbId->GetSize() <= iel ||
            (*Analysis::floatVectors["Electron_deltaEtaSC"]).GetSize() <= iel ||
            (*Analysis::intVectors["Electron_tightCharge"]).GetSize() <= iel)
        {
            std::cout << "Error: Out-of-bounds access detected at index " << iel << std::endl;
            continue;
        }
        /////////////////////////
        // Condition Debugging //
        /////////////////////////
        if(!passIsolation(Config::vetoelec_iso->At(iel), Config::vetoelec_iso_cut)) // isolation check
        {std::cout << "|| Failed isolation check for index " << iel << ", iso: " << Config::vetoelec_iso->At(iel) << ", iso_cut: " << Config::vetoelec_iso_cut << std::endl; continue;}
        if(!passKinematicCuts(elecs_pt->At(iel), elecs_eta->At(iel), Config::vetoelec_pt, Config::vetoelec_eta)) // kinematic cuts check
        {/*std::cout << "|| Failed kinematic cuts for index " << iel << ", elecs_pt->At(iel): " << elecs_pt->At(iel) << ", lep_pt: " << lep_pt << ", elecs_eta->At(iel): " << elecs_eta->At(iel) << ", lep_eta: " << lep_eta << std::endl;*/ continue;}
        if(!elecSCBId(Config::vetoelec_scbId->At(iel), eleid_scbcut)) // scbId check
        {/*std::cout << "|| Failed SCB ID for index " << iel << " Config::vetoelec_scbId->At(iel): " << Config::vetoelec_scbId->At(iel) << ", eleid_scbcut: " << eleid_scbcut << std::endl;*/ continue;}
        double combinedEta = (*Analysis::floatVectors["Electron_deltaEtaSC"])[iel] + elecs_eta->At(iel);
        if(fabs(combinedEta) > 1.4442 && fabs(combinedEta) < 1.566) // eta range check
        {/*std::cout << "|| Failed eta range for index " << iel << ", combinedEta: " << combinedEta << std::endl;*/ continue;}
        if(!elecCharge((*Analysis::intVectors["Electron_tightCharge"])[iel])) // charge check
        {/*std::cout << "|| Failed charge check for index " << iel << ", (*Analysis::intVectors[\"Electron_tightCharge\"])[iel]: " << (*Analysis::intVectors["Electron_tightCharge"])[iel] << std::endl;*/ continue;}
        /////////////////////////////
        // Add final selected elec //
        /////////////////////////////
        elecsveto.push_back(elecs.at(iel));
        elecsveto_isoValue.push_back(Config::vetoelec_iso->At(iel));
        //std::cout << "|| Final selected electron index: " << iel << std::endl;
    }
    //std::cout << "|| " << elecsveto.size() << " electrons are seclected" << std::endl;
    //std::cout << "|| elecsveto_idx: " << elecsveto_idx.size() << std::endl;
    //for(int i = 0; i < elecsveto_idx.size(); i++)
    //{std::cout << "|| elecsveto_idx[" << i << "]: " << elecsveto_idx[i] << std::endl;}
} // end of LeptonSelector //

void Lepton::LeptonOrder()
{
    // Helper function for assigning leptons and debugging
    //std::cout << "||------------------------Lepton Ordering------------------------||" << std::endl;
    auto assignLeptons = [&](const std::vector<TLorentzVector>& leptons,
                             const std::string& chargeKey,
                             const std::vector<int>& indices,
                             int idx1, int idx2)
    {
        if(indices.size() <= std::max(idx1, idx2))
        {
            std::cerr << "Error: indices size (" << indices.size() 
                      << ") is smaller than required index." << std::endl;
            return;
        }
        //std::cout << "|| Size of indices: " << indices.size()<< std::endl;
        //std::cout << "|| Size of leptons: " << leptons.size()<< std::endl;
        Analysis::Lep1 = leptons.at(indices[idx1]);
        Analysis::Lep2 = leptons.at(indices[idx2]);
        //std::cout << "|| End Lep1 & Lep2 assignment !  "<< std::endl;

        ///////////////////////////
        // Debugging information //
        ///////////////////////////
        //std::cerr << "Assigning Lep1 and Lep2:\n"
        //          << "  Lep1: Pt = " << Analysis::Lep1.Pt() << ", Eta = " << Analysis::Lep1.Eta() << "\n"
        //          << "  Lep2: Pt = " << Analysis::Lep2.Pt() << ", Eta = " << Analysis::Lep2.Eta() << std::endl;

        //std::cout << "Debug! point 1: assignLeptons in LeptonOrder"<< std::endl;
        if((*Analysis::intVectors[chargeKey])[indices[idx1]] < 0)
        {
            Analysis::Lep = leptons.at(indices[idx1]);
            Analysis::AnLep = leptons.at(indices[idx2]);
        }
        else
        {
            Analysis::Lep = leptons.at(indices[idx2]);
            Analysis::AnLep = leptons.at(indices[idx1]);
        }
        //std::cout << "Debug! point 2: assignLeptons in LeptonOrder"<< std::endl;
        ///////////////////////////////
        // Additional debugging info //
        ///////////////////////////////
        //std::cerr << "Assigning Lep and AnLep:\n"
        //          << "  Lep: Charge < 0, Pt = " << Analysis::Lep.Pt() << "\n"
        //          << "  AnLep: Charge > 0, Pt = " << Analysis::AnLep.Pt() << std::endl;
        //std::cout << "Debug! point 3: assignLeptons in LeptonOrder"<< std::endl;
    }; // end of assignLeptons //
    //std::cout << "|| Finished assignLeptons" << std::endl;
    ///////////////////////
    // Dimuon decay mode //
    ///////////////////////
    if(TString(Analysis::Decaymode).Contains("dimuon"))
    {
        //std::cout << "||--Decay mode: Dimuon--                                   ||" << std::endl;
        if(v_muon_idx.size() > 1)
        {
            //std::cout << "Debug! dimuon decay mode point 1" << std::endl;
            assignLeptons(muons, "Muon_charge", v_muon_idx, 0, 1);
            //std::cout << "Lep1: " << Analysis::Lep1.Pt() << ", " << Analysis::Lep1.Eta() << ", " << Analysis::Lep1.Phi() << ", " << Analysis::Lep1.M() << ", " << (*Analysis::intVectors["Muon_charge"])[v_muon_idx[0]] << std::endl;
            //std::cout << "Lep2: " << Analysis::Lep2.Pt() << ", " << Analysis::Lep2.Eta() << ", " << Analysis::Lep2.Phi() << ", " << Analysis::Lep2.M() << ", " << (*Analysis::intVectors["Muon_charge"])[v_muon_idx[1]] << std::endl;
        }
        else
        {
            //std::cerr << "Lepton TLorentzVector Error: v_muon_idx is empty or too small for Decaymode = dimuon" 
            //          << " (size: " << v_muon_idx.size() << ")" << std::endl;
            return;
        }
        //std::cout << "Debug! dimuon decay mode point 2" << std::endl;
    }
    ///////////////////////////
    // Dielectron decay mode //
    ///////////////////////////
    else if(TString(Analysis::Decaymode).Contains("dielec"))
    {
        //std::cout << "||--Decay mode: Dielectron--                               ||" << std::endl;
        if(v_electron_idx.size() > 1)
        {assignLeptons(elecs, "Electron_charge", v_electron_idx, 0, 1);}
        else
        {
            //std::cerr << "Lepton TLorentzVector Error: v_electron_idx is empty or too small for Decaymode = dielec" 
            //          << " (size: " << v_electron_idx.size() << ")" << std::endl;
            return;
        }
    }
    //////////////////////////////
    // Muon-Electron decay mode //
    //////////////////////////////
    else if(TString(Analysis::Decaymode).Contains("muel"))
    {
        //std::cout << "||--Decay mode: Muon-Electron--                            ||" << std::endl;
        if(v_muon_idx.size() > 0 && v_electron_idx.size() > 0)
        {
            if(muons.at(v_muon_idx.at(0)).Pt() > elecs.at(v_electron_idx.at(0)).Pt())
            {
                Analysis::Lep1 = muons.at(v_muon_idx.at(0));
                Analysis::Lep2 = elecs.at(v_electron_idx.at(0));
            }
            else
            {
                Analysis::Lep1 = elecs.at(v_electron_idx.at(0));
                Analysis::Lep2 = muons.at(v_muon_idx.at(0));
            }
            if((*Analysis::intVectors["Muon_charge"])[v_muon_idx[0]] < 0)
            {
                Analysis::Lep = muons.at(v_muon_idx.at(0));
                Analysis::AnLep = elecs.at(v_electron_idx.at(0));
            }
            else
            {
                Analysis::Lep = elecs.at(v_electron_idx.at(0));
                Analysis::AnLep = muons.at(v_muon_idx.at(0));
            }
        }
        else
        {
            std::cerr << "Lepton TLorentzVector Error: v_muon_idx size = " 
                      << v_muon_idx.size() << ", v_electron_idx size = " 
                      << v_electron_idx.size() << " for Decaymode = muel" << std::endl;
            return;
        }
    }
    ///////////////////////
    // Invalid Decaymode //
    ///////////////////////
    else
    {std::cerr << "Lepton TLorentzVector Error: Decaymode = " << Analysis::Decaymode << std::endl;}
} // end of LeptonOrder //

bool Lepton::LeptonsPtAddtional()
{
    bool lepptadd = false;
    if( TString(Analysis::Decaymode).Contains( "dimuon" ) ||
         TString(Analysis::Decaymode).Contains( "dielec" ) ||                                                                                      
         TString(Analysis::Decaymode).Contains( "muel" ) )
    { if( Analysis::Lep1.Pt() > 25 && Analysis::Lep2.Pt() > 20 ) {lepptadd=true;} } //Pt of Leading Lepton should be over than 25 GeV and Second Leading Lepton Pt should be over than 20 GeV.
    else if( TString(Analysis::Decaymode).Contains( "muonJet" ) )
    { if( v_lepton_idx.size() == 1){lepptadd=true;} }
    else { std::cerr << "?? something wrong " << std::endl; }                                                                                      
    return lepptadd;
} // end of LeptonsPtAddtional //

/// step 1 : Di-Lepton mass > 20 GeV ///
bool Lepton::DiLeptonMassCut()
{        
   bool dimu_masscut = false;
   if( ((Analysis::Lep1)+(Analysis::Lep2)).M() > 20 ){ dimu_masscut = true; }
   return dimu_masscut;
} // end of DiLeptonMassCut //

/// step 2 : Z-mass veto < 76 GeV or > 106 GeV ///
bool Lepton::ZMassVeto()
{
    zmassveto = false;
    if( TString(Analysis::Decaymode).Contains("dimuon") ||
        TString(Analysis::Decaymode).Contains("dielec") )
    {
        if( ((Analysis::Lep1)+(Analysis::Lep2)).M() < 76 || ((Analysis::Lep1)+(Analysis::Lep2)).M() > 106 )
        {zmassveto = true;}
    }
    return zmassveto;
} // end of ZMassVeto //

bool Lepton::NumIsoLeptons()
{
    bool numLeptons = true;
    if( TString(Analysis::Decaymode).Contains( "dimuon" ) || TString(Analysis::Decaymode).Contains( "dielec" ))
    { if( v_lepton_idx.size() <= 1 ) {numLeptons=false;} } //Requiring  Two isolated Lepton
    else if( TString(Analysis::Decaymode).Contains( "muel" ) )
    {
        if(v_muon_idx.size() < 1 || v_electron_idx.size() < 1 ) {numLeptons=false; }
        else {}
    }
    else if( TString(Analysis::Decaymode).Contains( "muonJet" ) )
    { if( v_lepton_idx.size() != 1){numLeptons=false;} }
    else { std::cerr << "?? something wrong " << std::endl; }
    return numLeptons;
} // end of NumIsoLeptons //

/////////////////////////////////////////////
// Make Muon and Electron Collection logic //
/////////////////////////////////////////////
void Lepton::MakeMuonCollection()
{
    // make muon collection tlorentzvectors//
    //std::cout << "||--Muon Collection--                                      ||" << std::endl;
    muons.clear();
    for(int imu = 0; imu <muons_pt->GetSize(); ++imu)
    {
        //std::cout << "Debug! MuonCollection: imu " << imu << " muons_pt->At(imu) : " << muons_pt->At(imu) << " muons_eta->At(imu) : " << muons_eta->At(imu) << " muons_phi->At(imu) : " << muons_phi->At(imu) << " muons_M->At(imu)  : "<< muons_M->At(imu)  << std::endl; 
        muons.push_back(Config::createLorentzVector(muons_pt->At(imu), muons_eta->At(imu), muons_phi->At(imu), muons_M->At(imu) ));
    }
    //std::cout << "|| size of muons.size : " << muons.size() << std::endl; 
    return;
} // end of MakeMuonCollection //

void Lepton::MakeElecCollection()
{
    // Electron collection logic
    //std::cout << "||--Electron Collection--                                  ||" << std::endl;
    elecs.clear();
    for(int iele = 0; iele < elecs_pt->GetSize(); ++iele)
    {
        //std::cout << "Debug! ElecCollection: iele " << iele << " elecs_pt->At(iele) : " << elecs_pt->At(iele) << " elecs_eta->At(iele) : " << elecs_eta->At(iele) << " elecs_phi->At(iele) : " << elecs_phi->At(iele) << " elecs_M->At(iele)  : "<< elecs_M->At(iele)  << std::endl; 
        elecs.push_back(Config::createLorentzVector(elecs_pt->At(iele), elecs_eta->At(iele), elecs_phi->At(iele), elecs_M->At(iele) )); 
    }
    //std::cout << "|| size of elecs.size : " << elecs.size() << std::endl; 
    return;
} // end of MakeElecCollection //

bool Lepton::ThirdLeptonVeto(std::vector<int> vLep1Idx, std::vector<int> vLep2Idx)
{
    /* Description of ThirdLeptonVeto
    1. candidate leptons are already selected, save their idices in v_lepton_idx by LeptonSelector()
    2. So, we just need to check it from the 3rd lepton in v_lepton_idx (loop can be started from 2)
    --- About channels ---
    (1) dimuon and dielectron:
        - check the charge of the first and second leptons, they should be opposite
        - check the 3rd lepton's pt, eta, iso cut, and id if the size of v_lepton_idx is greater than 3(>=3)
    (2) muon-electron:
        - check the size of v_muon_idx and v_electron_idx, they should be greater than 0
        - check the 3rd lepton's pt, eta, iso cut, and id if the size of each v_lepton_idx is greater than 2(>=2)
    --- procedure ---
    (1) if the size of v_lepton_idx is less than 1(<=1), that event is vetoed
    (2) the size of v_lepton_idx is greater than 3(>=3), do checking again the charge of the leptons (idx 0, 1)
    (3) check the pt, eta, iso cut, and id from the 3rd lepton(idx 2)
    (4) a lepton is passed, then it is regarded as a candidate lepton(3rd lepton) and the event will be vetoed
    v_lepton_idx에 기준을 통과한 렙톤들이 들어가 있음
    그럼 v_lepton_idx 사이즈만큼 루프를 돌면서
    처음 2개의 렙톤(idx번호가 0, 1)에 대해서 먼저 다시한번 확인하고 charge까지 확인
    그 다음 3번째 렙톤이 있으면 pt, eta, iso, id를 확인해서 candidate 렙톤인 경우 그 이벤트는 veto
    */
    bool thirdLepVeto = false;
    if(vLep1Idx.size() <= 1) {thirdLepVeto = false;}
    else if(vLep1Idx.size() >= 3)
    {
        if(TString(Analysis::Decaymode).Contains("dimuon"))
        {
            if((*Analysis::intVectors["Muon_charge"])[vLep1Idx[0]] == (*Analysis::intVectors["Muon_charge"])[vLep1Idx[1]]) {thirdLepVeto = true;}
            for(int i=2; i<vLep1Idx.size(); i++)
            {
                if(passKinematicCuts(muons_pt->At(vLep1Idx[i]), muons_eta->At(vLep1Idx[i]), lep_pt, lep_eta) &&
                   passIsolation(Config::muons_iso->At(vLep1Idx[i]), lepisocut) &&
                   passId(Config::muons_Id->At(vLep1Idx[i])))
                {
                    thirdLepVeto = true;
                    //std::cout << "|| 3rd lepton vetoed" << std::endl;
                    //std::cout << "|| pt: " << muons_pt->At(vLep1Idx[i]) << ", lep_pt: " << lep_pt << std::endl;
                    //std::cout << "|| eta: " << muons_eta->At(vLep1Idx[i]) << ", lep_eta: " << lep_eta << std::endl;
                    //std::cout << "|| iso: " << Config::muons_iso->At(vLep1Idx[i]) << ", lepisocut: " << lepisocut << std::endl;
                    //std::cout << "|| id: " << Config::muons_Id->At(vLep1Idx[i]) << ", Config id: " << Config::lepId << std::endl;
                }
            }
        }
        else if(TString(Analysis::Decaymode).Contains("dielec"))
        {
            if((*Analysis::intVectors["Electron_charge"])[vLep1Idx[0]] == (*Analysis::intVectors["Electron_charge"])[vLep1Idx[1]]) {thirdLepVeto = true;}
            for(int i=2; i<vLep1Idx.size(); i++)
            {
                if( passKinematicCuts(elecs_pt->At(vLep1Idx[i]), elecs_eta->At(vLep1Idx[i]), lep_pt, lep_eta) &&
                    elecSCBId(Config::elecs_scbId->At(vLep1Idx[i]), eleid_scbcut) &&
                    (fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[vLep1Idx[i]] + elecs_eta->At(vLep1Idx[i])) > 1.4442 &&
                    fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[vLep1Idx[i]] + elecs_eta->At(vLep1Idx[i])) < 1.566) &&
                    elecCharge((*Analysis::intVectors["Electron_tightCharge"])[vLep1Idx[i]]) )
                {thirdLepVeto = true;}
            }
        }
        else if(TString(Analysis::Decaymode).Contains("muel"))
        {
            if(vLep1Idx.size() < 1 || vLep2Idx.size() < 1 ) {thirdLepVeto = true;}
            if((*Analysis::intVectors["Muon_charge"])[vLep1Idx[0]] == (*Analysis::intVectors["Electron_charge"])[vLep2Idx[0]]) {thirdLepVeto = true;}
            for(int i=2; i<vLep1Idx.size(); i++)
            {
                if(passKinematicCuts(muons_pt->At(vLep1Idx[i]), muons_eta->At(vLep1Idx[i]), lep_pt, lep_eta) &&
                   passIsolation(Config::muons_iso->At(vLep1Idx[i]), lepisocut) &&
                   passId(Config::muons_Id->At(vLep1Idx[i])))
                {thirdLepVeto = true;}
            }
            for(int i=2; i<vLep2Idx.size(); i++)
            {
                if( passKinematicCuts(elecs_pt->At(vLep2Idx[i]), elecs_eta->At(vLep2Idx[i]), lep_pt, lep_eta) &&
                    elecSCBId(Config::elecs_scbId->At(vLep2Idx[i]), eleid_scbcut) &&
                    (fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[vLep2Idx[i]] + elecs_eta->At(vLep2Idx[i])) > 1.4442 &&
                    fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[vLep2Idx[i]] + elecs_eta->At(vLep2Idx[i])) < 1.566) &&
                    elecCharge((*Analysis::intVectors["Electron_tightCharge"])[vLep2Idx[i]]) )
                {thirdLepVeto = true;}
            }
        }
        else {std::cout << " 3rd lepton veto is not implemented for this decay mode: " << Analysis::Decaymode << std::endl;}
    }
    return thirdLepVeto;
} // end of ThirdLeptonVeto //