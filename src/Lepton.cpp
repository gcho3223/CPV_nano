#include "../interface/Lepton.h"
#include "../interface/Config.h"
//////////////////////////////////////////////////////////////////////////////////////
//                              Initialize variables                                //
//////////////////////////////////////////////////////////////////////////////////////
/// Lepton ///
double Lepton::lep_pt;
double Lepton::lep_eta;
double Lepton::lepisocut;
std::vector<int> Lepton::v_lep_idx_temp;
std::vector<int> Lepton::v_lepton_idx;
/// Muon ////
std::vector<TLorentzVector> Lepton::muons;
TTreeReaderArray<Float_t>* Lepton::muons_pt;
TTreeReaderArray<Float_t>* Lepton::muons_eta;
TTreeReaderArray<Float_t>* Lepton::muons_phi;
TTreeReaderArray<Float_t>* Lepton::muons_M;
std::vector<TLorentzVector> Lepton::muonsveto;
double Lepton::muonisocut;
std::vector<int> Lepton::v_muon_idx_temp;
std::vector<int> Lepton::v_muon_idx;
/// Electron ///
std::vector<TLorentzVector> Lepton::elecs;
TTreeReaderArray<Float_t>* Lepton::elecs_pt;
TTreeReaderArray<Float_t>* Lepton::elecs_eta;
TTreeReaderArray<Float_t>* Lepton::elecs_phi;
TTreeReaderArray<Float_t>* Lepton::elecs_M;
std::vector<TLorentzVector> Lepton::elecsveto;
double Lepton::elecisocut;
std::vector<int> Lepton::v_electron_idx_temp;
std::vector<int> Lepton::v_electron_idx;
/// Z mass veto ///
bool Lepton::zmassveto;

void Lepton::LeptonSelector()
{
    std::cout << "=============================================================" << std::endl;  
    std::cout << "||                     Lepton Selector                     ||" << std::endl;
    std::cout << "=============================================================" << std::endl;

    // Clear output vectors
    v_lepton_idx.clear();
    v_muon_idx.clear();
    v_electron_idx.clear();
    v_lep_idx_temp.clear();
    v_muon_idx_temp.clear();
    v_electron_idx_temp.clear();
    elecsveto.clear();
    muonsveto.clear();

    // Prepare corrected collections
    MakeMuonCollection();
    MakeElecCollection();

    // Lambda functions for common checks
    auto passKinematicCuts = [](float _pt, float _eta, float ptCut, float etaCut)
    {return _pt >= ptCut && fabs(_eta) <= etaCut;};

    auto passIsolation = [](float iso, float isoCut)
    {return iso <= isoCut;};

    auto passId = [](bool id) {return id;};

    auto elecSCBId = [](int id, int idcut) {return id >= idcut;};
    auto elecCharge = [](int id) {return id > 0;};

    /////////////////
    // Dimuon case //
    /////////////////
    if(TString(Analysis::Decaymode).Contains("dimuon"))
    {
        Int_t nmu = Config::muons_iso->GetSize();
        for(int i = 0; i < nmu; ++i)
        {
            if(!passIsolation(Config::muons_iso->At(i), lepisocut) || 
                !passKinematicCuts(muons_pt->At(i), muons_eta->At(i),lep_pt, lep_eta) || 
                !passId(Config::muons_Id->At(i)))
            {continue;}
            v_lep_idx_temp.push_back(i);
            //v_lep_idx_temp.push_back(i); // I think this is copy mistake

            if(v_lep_idx_temp.size() == 1 && v_lepton_idx.empty()) {v_lepton_idx.push_back(i);}
            else if(v_lep_idx_temp.size() > 0 && (*Analysis::intVectors["Muon_charge"])[v_lep_idx_temp[0]] != (*Analysis::intVectors["Muon_charge"])[i]) // push_back opposite charged muon
            {v_lepton_idx.push_back(i);}
        }
        // Assign final indices to v_muon_idx
        v_muon_idx = v_lepton_idx;
        std::cout << "|| Dimuon case: size of muon index(v_muon_idx.size()) : " << v_muon_idx.size() << std::endl;
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
            //Ele_lep_sel = elecs.at(i);
            if(!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), lep_pt, lep_eta) ||
                !elecSCBId(Config::elecs_scbId->At(i), Analysis::eleid_scbcut) ||
                (fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) > 1.4442 &&
                 fabs((*Analysis::floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) < 1.566) ||
                !(elecCharge((*Analysis::intVectors["Electron_tightCharge"])[i])))
            {continue;}
            v_lep_idx_temp.push_back(i);
            
            if(v_lep_idx_temp.size() == 1 && v_lepton_idx.empty())
            {
                v_lepton_idx.push_back(i);
                v_electron_idx.push_back(i); // I think this is copy mistake
            }
            else if(v_lep_idx_temp.size() > 0 &&
                (*Analysis::intVectors["Electron_charge"])[v_lep_idx_temp[0]] != (*Analysis::intVectors["Electron_charge"])[i]) // push_back opposite charged electron
            {
                v_lepton_idx.push_back(i);
                v_electron_idx.push_back(i); // I think this is copy mistake
            }
        }
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
                !passKinematicCuts(muons_pt->At(i), muons_eta->At(i),lep_pt, lep_eta) || 
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
          if (!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), lep_pt, lep_eta) || 
          !elecSCBId(Config::elecs_scbId->At(i), Analysis::eleid_scbcut) || 
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
                //if ((*Analysis::intVectors["Muon_charge"])[v_muon_idx[0]] != (*Analysis::intVectors["Electron_charge"])[i]) // I think this is correct
                if ((*Analysis::intVectors["Muon_charge"])[v_lep_idx_temp[0]] != (*Analysis::intVectors["Electron_charge"])[i])
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

    std::cout << "||--Muon/Electron Selection--                              ||" << std::endl;
    ////////////////////
    // Muon Selection //
    ////////////////////
    Int_t nmu = Config::muonsveto_iso->GetSize();

    for (int imu = 0; imu < nmu; ++imu)
    {
        if (!passIsolation(Config::muons_iso->At(imu), lepisocut) || 
            !passKinematicCuts(muons_pt->At(imu), muons_eta->At(imu), lep_pt, lep_eta) ||
            !passId(Config::muons_Id->At(imu)))
        {continue;}
        muonsveto.push_back(muons.at(imu));
    }
    std::cout << "|| " << muonsveto.size() << " muons are seclected" << std::endl;

    if (Config::elecsveto_iso == nullptr)
    {
        std::cerr << "Error: Config::elecsveto_iso is null!" << std::endl;
        return; // 또는 continue;
    }
    ////////////////////////
    // Electron Selection //
    ////////////////////////
    //std::cout << "Debug! Electron seclection: elecs_iso: " << Config::elecsveto_iso->GetSize() << ", elecs.size: " << elecs.size() << std::endl;
    Int_t nel = Config::elecsveto_iso->GetSize();
    if (nel != elecs.size())
    {
        std::cerr << "Error: Size mismatch! nel: " << nel << ", elecs.size: " << elecs.size() << std::endl;
        return; // 또는 continue;
    }
    if (!elecs_pt || !elecs_eta || !Analysis::floatVectors["Electron_deltaEtaSC"] || !Analysis::intVectors["Electron_tightCharge"])
    {
        if (!elecs_pt) std::cerr << "Error: elecs_pt is null!" << std::endl;
        if (!elecs_eta) std::cerr << "Error: elecs_eta is null!" << std::endl;
        if (!Analysis::floatVectors["Electron_deltaEtaSC"]) std::cerr << "Error: floatVectors['Electron_deltaEtaSC'] is null!" << std::endl;
        if (!Analysis::intVectors["Electron_tightCharge"]) std::cerr << "Error: intVectors['Electron_tightCharge'] is null!" << std::endl;
        return;
    }
    for (int iel = 0; iel < nel; ++iel)
    {
        //std::cout << "|| Processing electron index: " << iel << std::endl;
        // 벡터 및 크기 확인
        if (!elecs_pt || !elecs_eta || !Config::elecsveto_scbId || !Analysis::floatVectors["Electron_deltaEtaSC"] || !Analysis::intVectors["Electron_tightCharge"])
        {
            if (!elecs_pt) std::cerr << "Error: elecs_pt is nullptr!" << std::endl;
            if (!elecs_eta) std::cerr << "Error: elecs_eta is nullptr!" << std::endl;
            if (!Config::elecsveto_scbId) std::cerr << "Error: elecsveto_scbId is nullptr!" << std::endl;
            if (!Analysis::floatVectors["Electron_deltaEtaSC"]) std::cerr << "Error: floatVectors[\"Electron_deltaEtaSC\"] is nullptr!" << std::endl;
            if (!Analysis::intVectors["Electron_tightCharge"]) std::cerr << "Error: intVectors[\"Electron_tightCharge\"] is nullptr!" << std::endl;
            if (!elecs_pt || !elecs_eta || !Config::elecsveto_scbId || !Analysis::floatVectors["Electron_deltaEtaSC"] || !Analysis::intVectors["Electron_tightCharge"])
            {
                std::cerr << "Error: Null pointer detected for variables!" << std::endl;
                return;
            }
            std::cerr << "Error: Null pointer detected for variables!" << std::endl;
            return;
        }
        if (elecs_pt->GetSize() <= iel || elecs_eta->GetSize() <= iel ||
            Config::elecsveto_scbId->GetSize() <= iel ||
            (*Analysis::floatVectors["Electron_deltaEtaSC"]).GetSize() <= iel ||
            (*Analysis::intVectors["Electron_tightCharge"]).GetSize() <= iel)
        {
            std::cerr << "Error: Out-of-bounds access detected at index " << iel << std::endl;
            continue;
        }
        /////////////////////////
        // Condition Debugging //
        /////////////////////////
        if (!passKinematicCuts(elecs_pt->At(iel), elecs_eta->At(iel), lep_pt, lep_eta))
        {
            std::cout << "|| Failed kinematic cuts for index " << iel << std::endl;
            continue;
        }
        if (!elecSCBId(Config::elecsveto_scbId->At(iel), Analysis::eleid_scbcut))
        {
            std::cout << "|| Failed SCB ID for index " << iel << std::endl;
            continue;
        }
        double combinedEta = (*Analysis::floatVectors["Electron_deltaEtaSC"])[iel] + elecs_eta->At(iel);
        if (fabs(combinedEta) > 1.4442 && fabs(combinedEta) < 1.566)
        {
            std::cout << "|| Failed eta range for index " << iel << std::endl;
            continue;
        }
        if (!elecCharge((*Analysis::intVectors["Electron_tightCharge"])[iel]))
        {
            std::cout << "|| Failed charge check for index " << iel << std::endl;
            continue;
        }
        /////////////////////////////
        // Add final selected elec //
        /////////////////////////////
        elecsveto.push_back(elecs.at(iel));
        //std::cout << "|| Final selected electron index: " << iel << std::endl;
    }
    std::cout << "|| " << elecsveto.size() << " electrons are seclected" << std::endl;
} // end of LeptonSelector //

void Lepton::LeptonOrder()
{
    // Helper function for assigning leptons and debugging
    std::cout << "=============================================================" << std::endl;  
    std::cout << "||                     Lepton Ordering                     ||" << std::endl;
    std::cout << "=============================================================" << std::endl;
    auto assignLeptons = [&](const std::vector<TLorentzVector>& leptons,
                             const std::string& chargeKey,
                             const std::vector<int>& indices,
                             int idx1, int idx2)
    {
        if (indices.size() <= std::max(idx1, idx2))
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
        if ((*Analysis::intVectors[chargeKey])[indices[idx1]] < 0)
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
    if (TString(Analysis::Decaymode).Contains("dimuon"))
    {
        std::cout << "||--Decay mode: Dimuon--                                   ||" << std::endl;
        if (v_muon_idx.size() > 1)
        {
            //std::cout << "Debug! dimuon decay mode point 1" << std::endl;
            assignLeptons(muons, "Muon_charge", v_muon_idx, 0, 1);
        }
        else
        {
            //std::cerr << "Lepton TLorentzVector Error: v_muon_idx is empty or too small for Decaymode = dimuon" 
            //          << " (size: " << v_muon_idx.size() << ")" << std::endl;
            return;
        }
        //std::cout << "Debug! dimuon decay mode point 2" << std::endl;
    }
    else if (TString(Analysis::Decaymode).Contains("dielec"))
    {
        ///////////////////////////
        // Dielectron decay mode //
        ///////////////////////////
        std::cout << "||--Decay mode: Dielectron--                               ||" << std::endl;
        if (v_electron_idx.size() > 1)
        {assignLeptons(elecs, "Electron_charge", v_electron_idx, 0, 1);}
        else
        {
            //std::cerr << "Lepton TLorentzVector Error: v_electron_idx is empty or too small for Decaymode = dielec" 
            //          << " (size: " << v_electron_idx.size() << ")" << std::endl;
            return;
        }
    }
    else if (TString(Analysis::Decaymode).Contains("muel"))
    {
        //////////////////////////////
        // Muon-Electron decay mode //
        //////////////////////////////
        std::cout << "||--Decay mode: Muon-Electron--                            ||" << std::endl;
        if (v_muon_idx.size() > 0 && v_electron_idx.size() > 0)
        {
            if (muons.at(v_muon_idx.at(0)).Pt() > elecs.at(v_electron_idx.at(0)).Pt())
            {
                Analysis::Lep1 = muons.at(v_muon_idx.at(0));
                Analysis::Lep2 = elecs.at(v_electron_idx.at(0));
            }
            else
            {
                Analysis::Lep1 = elecs.at(v_electron_idx.at(0));
                Analysis::Lep2 = muons.at(v_muon_idx.at(0));
            }
            if ((*Analysis::intVectors["Muon_charge"])[v_muon_idx[0]] < 0)
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
    if ( TString(Analysis::Decaymode).Contains( "dimuon" ) ||
         TString(Analysis::Decaymode).Contains( "dielec" ) ||                                                                                      
         TString(Analysis::Decaymode).Contains( "muel" ) )
    { if( Analysis::Lep1.Pt() > 25 && Analysis::Lep2.Pt() > 20 ) {lepptadd=true;} } //Pt of Leading Lepton should be over than 25 GeV and Second Leading Lepton Pt should be over than 20 GeV.
    else if ( TString(Analysis::Decaymode).Contains( "muonJet" ) )
    { if( v_lepton_idx.size() == 1){lepptadd=true;} }
    else { std::cerr << "?? something wrong " << std::endl; }                                                                                      
    return lepptadd;
} // end of LeptonsPtAddtional //

/// step 1 : Di-Lepton mass > 20 GeV ///
bool Lepton::DiLeptonMassCut()
{        
   bool dimu_masscut = false;
   if ( ((Analysis::Lep1)+(Analysis::Lep2)).M() > 20 ){ dimu_masscut = true; }
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
        if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1 ) {numLeptons=false; }
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
    for (int imu = 0; imu <muons_pt->GetSize(); ++imu)
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
    for (int iele = 0; iele < elecs_pt->GetSize(); ++iele)
    {
        //std::cout << "Debug! ElecCollection: iele " << iele << " elecs_pt->At(iele) : " << elecs_pt->At(iele) << " elecs_eta->At(iele) : " << elecs_eta->At(iele) << " elecs_phi->At(iele) : " << elecs_phi->At(iele) << " elecs_M->At(iele)  : "<< elecs_M->At(iele)  << std::endl; 
        elecs.push_back(Config::createLorentzVector(elecs_pt->At(iele), elecs_eta->At(iele), elecs_phi->At(iele), elecs_M->At(iele) )); 
    }
    //std::cout << "|| size of elecs.size : " << elecs.size() << std::endl; 
    return;
} // end of MakeElecCollection //