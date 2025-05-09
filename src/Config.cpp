#include "../interface/Config.h"
#include "../interface/Lepton.h"
#include "../interface/Jet.h"
#include "../interface/Trigger.h"
#include "../interface/MET.h"
#include "../interface/MCSaleFactor.h"
//-------------------------------Initialize variables-------------------------------//
TTreeReaderArray<bool>* Config::vetomuon_id;
TTreeReaderArray<float>* Config::muons_iso;
TTreeReaderArray<float>* Config::vetomuon_iso;
TTreeReaderArray<bool>* Config::muons_Id;
TTreeReaderArray<float>* Config::elecs_iso;
TTreeReaderArray<float>* Config::vetoelec_iso;
TTreeReaderArray<Int_t>* Config::elecs_scbId;
TTreeReaderArray<Int_t>* Config::vetoelec_scbId;
TTreeReaderArray<float>* Config::leptons_iso;
TTreeReaderArray<bool>* Config::leptons_Id;
TTreeReaderArray<Bool_t>* Config::elecs_mvaId;
TTreeReaderArray<Bool_t>* Config::vetoelec_mvaId;
TString Config::vetomu_isoid;
TString Config::vetomu_id;
float Config::vetomu_pt;
float Config::vetomu_eta;
TString Config::vetoelec_isoid;
TString Config::vetoelec_id;
float Config::vetoelec_pt;
float Config::vetoelec_eta;
int Config::num_dleptrig;
int Config::num_sleptrig;
TString Config::lepisotype;
TString Config::lepId;
TString Config::MuonIsoType;
float Config::vetomu_iso_cut;
TString Config::MuonId;
TString Config::ElecIsoType;
float Config::vetoelec_iso_cut;
TString Config::ElecId;
TString Config::JetId;
TString Config::JetPUID;
TString Config::JetbTag;
TString Config::BTagSFSys;
TString Config::BTagEffSys;
TString Config::MetSys;
TString Config::LepIdSFSys;
TString Config::LepIsoSFSys;
TString Config::LepRecoSFSys;
TString Config::LepTrackSFSys;
std::string Config::PileUpMCFile;
std::string Config::PileUpDATAFile;
TString Config::PileUpSys;
TString Config::L1PreFireSys;
int Config::PDFSys;
int Config::FactRenoSys;
TString Config::FragmentSys;
TString Config::DecayTableSys;
TString Config::TopPtSys;
int Config::elevetoid_scbcut;
double Config::muon_pt;
double Config::muon_eta;
double Config::elec_pt;
double Config::elec_eta;
double Config::jet_pt;
double Config::jet_eta;
int Config::jet_id;
int Config::jet_puid;

void Analysis::InitBranches(const std::string &branchListFile)
{
    std::cout << "||------------------------InitBranches------------------------||" << std::endl;
    std::cout << "|| branchListFile : " << branchListFile << std::endl;
    std::ifstream infile(branchListFile);
    if (!infile) {throw std::runtime_error("Error: Could not open branch list file " + branchListFile);}

    std::string line;
    while(std::getline(infile, line))
    {
        /// Skip empty lines ///
        if(line.empty()) continue;

        std::istringstream iss(line);
        std::string branchName, objectType, dataType, varType;

        /// Parse the line for branch name, object type, data type, and variable type (vector/single) ///
        if(std::getline(iss, branchName, ',') &&
           std::getline(iss, objectType, ',') &&
           std::getline(iss, dataType, ',') &&
           std::getline(iss, varType, ','))
        {
            /// Trim whitespace from each part (both ends) ///
            auto trim = [](std::string &s)
            {
                s.erase(s.find_last_not_of(" \t\n\r") + 1);
                s.erase(0, s.find_first_not_of(" \t\n\r"));
            };

            trim(branchName);
            trim(objectType);
            trim(dataType);
            trim(varType);

            /// Initialize based on data type and variable type ///
            if(dataType == "Bool_t" && varType == "single")
            {Analysis::boolSingles[branchName] = std::make_unique<TTreeReaderValue<Bool_t>>(fReader, branchName.c_str());}
            else if(dataType == "Int_t" && varType == "single")
            {intSingles[branchName] = std::make_unique<TTreeReaderValue<Int_t>>(fReader, branchName.c_str());}
            else if(dataType == "UInt_t" && varType == "single")
            {Analysis::uintSingles[branchName] = std::make_unique<TTreeReaderValue<UInt_t>>(fReader, branchName.c_str());}
            else if(dataType == "Float_t" && varType == "single")
            {Analysis::floatSingles[branchName] = std::make_unique<TTreeReaderValue<Float_t>>(fReader, branchName.c_str());}
            else if(dataType == "UChar_t" && varType == "single") // New case for UChar_t single value
            { ucharSingles[branchName] = std::make_unique<TTreeReaderValue<UChar_t>>(fReader, branchName.c_str());}
            else if(dataType == "Bool_t" && varType == "vector")
            {Analysis::boolVectors[branchName] = std::make_unique<TTreeReaderArray<Bool_t>>(fReader, branchName.c_str());}
            else if(dataType == "Int_t" && varType == "vector")
            {Analysis::intVectors[branchName] = std::make_unique<TTreeReaderArray<Int_t>>(fReader, branchName.c_str());}
            else if(dataType == "UInt_t" && varType == "vector")
            {uintVectors[branchName] = std::make_unique<TTreeReaderArray<UInt_t>>(fReader, branchName.c_str());}
            else if(dataType == "Float_t" && varType == "vector")
            {Analysis::floatVectors[branchName] = std::make_unique<TTreeReaderArray<Float_t>>(fReader, branchName.c_str());}
            else if(dataType == "UChar_t" && varType == "vector") // New case for UChar_t vector
            {ucharVectors[branchName] = std::make_unique<TTreeReaderArray<UChar_t>>(fReader, branchName.c_str());}
            else
            {
                std::cerr << "Warning: Unsupported data type or format in branch list: "
                          << dataType << ", " << varType << " for branch " << branchName << std::endl;
            }

        }
        else {std::cerr << "Warning: Invalid branch format in list: " << line << std::endl;}
    }
} // end of InitBranches //

void Config::SetVariables()
{
    //std::cout << "||------------------------Set Variables------------------------||" << std::endl;

    Analysis::Lumi       = Analysis::SSBConfReader->GetNumber( "Luminosity" );
    Analysis::RunPeriod  = Analysis::SSBConfReader->GetText( "RunRange" );
    Analysis::Decaymode  = Analysis::SSBConfReader->GetText( "Channel" ); // Channel
    Analysis::XsecTable_ = Analysis::SSBConfReader->GetText( "XSecTablesName" );

    /// Set Trigger List ///
    num_dleptrig = Analysis::SSBConfReader->Size( "dileptrigger" );
    num_sleptrig = Analysis::SSBConfReader->Size( "singleleptrigger" );
    /// Set Lepton ID and Iso type ///
    lepId        = Analysis::SSBConfReader->GetText( "ID_type" );
    lepisotype   = Analysis::SSBConfReader->GetText( "Iso_type" );
    /// Set dilepton trigger ///
    //std::cout << "||------------------------Dilepton trigger list------------------------||" << std::endl;
    for(int i =0; i < num_dleptrig; ++i)
    {
        //std::cout << Analysis::SSBConfReader->GetText("dileptrigger",i+1) << std::endl;
        std::string tmptrg = Analysis::SSBConfReader->GetText("dileptrigger",i+1);
        Trigger::DLtrigName.push_back( removeSubstring( tmptrg, "_v") );
        Trigger::trigName.push_back( removeSubstring( tmptrg, "_v")  );
    }
    /// Set single lepton trigger ///
    //std::cout << "||------------------------Single lepton trigger list------------------------||" << std::endl;
    for(int i =0; i < num_sleptrig; ++i)
    {
        //std::cout << Analysis::SSBConfReader->GetText("singleleptrigger",i+1) << std::endl;
        std::string tmptrg = Analysis::SSBConfReader->GetText("singleleptrigger",i+1);
        Trigger::SLtrigName.push_back( removeSubstring(tmptrg,"_v") );
        Trigger::trigName.push_back( removeSubstring(tmptrg,"_v") );
    }
    /// Set trigger list ///
    for(int i = 0; i < Trigger::trigName.size(); ++i)
    {Trigger::triggerList[Trigger::trigName[i]] = DeepCopy<bool>(Analysis::boolSingles[Trigger::trigName[i]]);}
    /// Set Noise filter(MET, events filter) ///
    //std::cout << "||------------------------Noise filter list------------------------||" << std::endl;
    for(int i = 0; i < Analysis::SSBConfReader->Size("METFilters"); ++i)
    {
       //std::cout << "METFilters: " << SSBConfReader->GetText("METFilters",i+1) << std::endl;
       std::string tmpnoisefl = Analysis::SSBConfReader->GetText("METFilters",i+1);
       //std::cout << "METFilters: " <<  tmpnoisefl << std::endl;
       //std::cout << boolSingles[tmpnoisefl] << std::endl;
       if(Analysis::boolSingles[tmpnoisefl] ==NULL) {std::cout << "Error!!!" << tmpnoisefl << std::endl;}
       Analysis::noiseFilters[tmpnoisefl] = DeepCopy<bool>(Analysis::boolSingles[tmpnoisefl]);
       //std::cout << "test "<< std::endl;
    }

    /// Lep Info ( Di-Muon or Di-Electron ) ///
    lepisotype  = Analysis::SSBConfReader->GetText( "Iso_type" );
    lepId       = Analysis::SSBConfReader->GetText( "Lep_ID"   );
    /// Muon Info for Only MuEl channel ///
    MuonIsoType = Analysis::SSBConfReader->GetText( "MuonIso_type" );
    MuonId      = Analysis::SSBConfReader->GetText( "Muon_ID"      );
    /// Electron Info Only for MuEl channel ///
    ElecIsoType = Analysis::SSBConfReader->GetText( "ElecIso_type" );
    ElecId      = Analysis::SSBConfReader->GetText( "Elec_ID"      );  
    /// Jet Info ///
    JetId       = Analysis::SSBConfReader->GetText( "Jet_ID"   );
    JetPUID     = Analysis::SSBConfReader->GetText( "Jet_PUID" );
    JetbTag     = Analysis::SSBConfReader->GetText( "Jet_btag" );

    /// Kinematic variables for Object ///
    /// Di-Muon & Di-Electron channel///
    Lepton::lep_pt      = Analysis::SSBConfReader->GetNumber( "Lep_pt"  );
    Lepton::lep_eta     = Analysis::SSBConfReader->GetNumber( "Lep_eta" );
    Lepton::lepisocut   = Analysis::SSBConfReader->GetNumber( "Iso_cut" );
    /// Muon-Electron channel ///
    muon_pt             = Analysis::SSBConfReader->GetNumber( "Muon_pt"     );
    muon_eta            = Analysis::SSBConfReader->GetNumber( "Muon_eta"    );
    Lepton::muonisocut  = Analysis::SSBConfReader->GetNumber( "MuonIso_cut" );
    elec_pt             = Analysis::SSBConfReader->GetNumber( "Elec_pt"     );
    elec_eta            = Analysis::SSBConfReader->GetNumber( "Elec_eta"    );
    Lepton::elecisocut  = Analysis::SSBConfReader->GetNumber( "ElecIso_cut" );
    /// MET ///
    MET::met_cut        = Analysis::SSBConfReader->GetNumber( "MET_cut" );
    /// Jet ///
    jet_pt              = Analysis::SSBConfReader->GetNumber( "Jet_pt"  );
    jet_eta             = Analysis::SSBConfReader->GetNumber( "Jet_eta" );
    /// Jet Cleaning for All channels///
    vetomu_isoid        = Analysis::SSBConfReader->GetText( "Muiso_for_Jetcleanning"      );
    vetomu_id           = Analysis::SSBConfReader->GetText( "Muid_for_Jetcleanning"       );
    vetomu_iso_cut      = Analysis::SSBConfReader->GetNumber( "Muisocut_for_Jetcleanning" );
    vetomu_pt           = Analysis::SSBConfReader->GetNumber( "Mupt_for_Jetcleanning"     );
    vetomu_eta          = Analysis::SSBConfReader->GetNumber( "Mueta_for_Jetcleanning"    );
    vetoelec_isoid      = Analysis::SSBConfReader->GetText( "Eliso_for_Jetcleanning"      );
    vetoelec_id         = Analysis::SSBConfReader->GetText( "Elid_for_Jetcleanning"       );
    vetoelec_iso_cut    = Analysis::SSBConfReader->GetNumber( "Elisocut_for_Jetcleanning" );
    vetoelec_pt         = Analysis::SSBConfReader->GetNumber( "Elpt_for_Jetcleanning"     );
    vetoelec_eta        = Analysis::SSBConfReader->GetNumber( "Eleta_for_Jetcleanning"    );
    /// Jet Energy Resolution ///
    Jet::dojer          = Analysis::SSBConfReader->GetBool("DoJER");
    Jet::JetEnSys       = Analysis::SSBConfReader->GetText("JESsys"); //Getting Jet Energy Systematic Type ...
    Jet::JetResSys      = Analysis::SSBConfReader->GetText("JERsys"); //Getting Jet Smearing Systematic Type ... 
    /// BTagging ///
    BTagSFSys           = Analysis::SSBConfReader->GetText( "BtaggingSFSys"  );
    BTagEffSys          = Analysis::SSBConfReader->GetText( "BtaggingEffSys" );
    /// MET ///
    MetSys              = Analysis::SSBConfReader->GetText(   "METsys"       ); //Getting Jet Energy Systematic Type ...
    /// Systematic ///
    LepIdSFSys          = Analysis::SSBConfReader->GetText(   "LepIDSFSys"         );
    LepIsoSFSys         = Analysis::SSBConfReader->GetText(   "LepIsoSFSys"        );
    LepRecoSFSys        = Analysis::SSBConfReader->GetText(   "LepRecoSFSys"       );
    LepTrackSFSys       = Analysis::SSBConfReader->GetText(   "LepTrackSFSys"      );
    PDFSys              = Analysis::SSBConfReader->GetNumber( "PDFSys"             );
    FactRenoSys         = Analysis::SSBConfReader->GetNumber( "FactRenoSys"        );
    FragmentSys         = Analysis::SSBConfReader->GetText(   "FragmentSys"        );
    DecayTableSys       = Analysis::SSBConfReader->GetText(   "DecayTableSys"      );
    PileUpMCFile        = Analysis::SSBConfReader->GetText(   "PileUpMCFileName"   );
    PileUpDATAFile      = Analysis::SSBConfReader->GetText(   "PileUpDATAFileName" );
    PileUpSys           = Analysis::SSBConfReader->GetText(   "PileUpSys"          );
    L1PreFireSys        = Analysis::SSBConfReader->GetText(   "L1PreFireSys"       );
    TopPtSys            = Analysis::SSBConfReader->GetText(   "TopPtSys"           );
} // end of SetVariables //

void Config::SetObjectVariable()
{
    //std::cout << "||------------------------Set Object Variables------------------------||" << std::endl;
    /// Lepton Kinematic ///
    // muon //
    Lepton::muons_pt  = Analysis::floatVectors["Muon_pt"].get(); 
    Lepton::muons_eta = Analysis::floatVectors["Muon_eta"].get(); 
    Lepton::muons_phi = Analysis::floatVectors["Muon_phi"].get();
    Lepton::muons_M   = Analysis::floatVectors["Muon_mass"].get();
    // electron //
    Lepton::elecs_pt  = Analysis::floatVectors["Electron_pt"].get();
    Lepton::elecs_eta = Analysis::floatVectors["Electron_eta"].get(); 
    Lepton::elecs_phi = Analysis::floatVectors["Electron_phi"].get();
    Lepton::elecs_M   = Analysis::floatVectors["Electron_mass"].get();

    /////////////////////
    // Case of Di-Muon //
    /////////////////////
    if(TString(Analysis::Decaymode).Contains("dimuon"))
    {
        if(TString(lepisotype).Contains("PFIsodbeta03"))
        {
            if(Analysis::floatVectors["Muon_pfRelIso03_all"] == nullptr)
            {
                std::cerr << "Dimuon: Not found Muon_pfRelIso03_all in NanoAOD..." << std::endl;
                return;
            }
            leptons_iso = Analysis::floatVectors["Muon_pfRelIso03_all"].get();
            muons_iso   = Analysis::floatVectors["Muon_pfRelIso03_all"].get();
        }
        else if(TString(lepisotype).Contains("PFIsodbeta04"))
        {
            if(Analysis::floatVectors["Muon_pfRelIso04_all"] == nullptr)
            {
                std::cerr << "Dimuon: Not found Muon_pfRelIso04_all in NanoAOD..." << std::endl;
                return;
            }
            leptons_iso = Analysis::floatVectors["Muon_pfRelIso04_all"].get();
            muons_iso   = Analysis::floatVectors["Muon_pfRelIso04_all"].get();
        }
        else
        {
            std::cerr << "Muon Iso type Error" << std::endl;
            return;
        }
        if(leptons_iso != nullptr)
        {
            //std::cout << "|| At ObjSetup leptons_iso, muons_iso size: " << leptons_iso->GetSize() << ", " << muons_iso->GetSize() << std::endl;
        }
        else {std::cerr << "Error: leptons_iso is null after assignment!" << std::endl;}
        
        /// Muon ID ///
        if     (TString(lepId).Contains("Loose" )) { leptons_Id = Analysis::boolVectors["Muon_looseId"].get();  muons_Id = Analysis::boolVectors["Muon_looseId"].get(); }
        else if(TString(lepId).Contains("Medium")) { leptons_Id = Analysis::boolVectors["Muon_mediumId"].get(); muons_Id = Analysis::boolVectors["Muon_mediumId"].get();}
        else if(TString(lepId).Contains("Tight" )) { leptons_Id = Analysis::boolVectors["Muon_tightId"].get();  muons_Id = Analysis::boolVectors["Muon_tightId"].get(); }
        else { std::cout << "Muon ID Error" << std::endl; }

        if(leptons_Id != nullptr)
        {
            //std::cout << "|| leptons_Id/muons_Id size: " << leptons_Id->GetSize() << std::endl;
        }
        else {std::cerr << "Error: leptons_Id is null after assignment!" << std::endl;}

    }
    /////////////////////////
    // Case of Di-Electron //
    /////////////////////////
    else if(TString(Analysis::Decaymode).Contains("dielec"))
    {
        /// Electron iso type ///
        if(TString(lepisotype).Contains("PFIsoRho03"))
        {
            if(Analysis::floatVectors["Electron_pfRelIso03_all"] == nullptr)
            {
                std::cerr << "Dielectron: Not found Electron_pfRelIso03_all in NanoAOD..." << std::endl;
                return;
            }
            leptons_iso = Analysis::floatVectors["Electron_pfRelIso03_all"].get();
            elecs_iso   = Analysis::floatVectors["Electron_pfRelIso03_all"].get();
        }
        else if(TString(lepisotype).Contains("PFIsoRho04"))
        {
            if(Analysis::floatVectors["Electron_pfRelIso04_all"] == nullptr)
            {
                std::cerr << "Dielectron: Not found Electron_pfRelIso04_all in NanoAOD..." << std::endl;
                return;
            }
            leptons_iso = Analysis::floatVectors["Electron_pfRelIso04_all"].get();
            elecs_iso   = Analysis::floatVectors["Electron_pfRelIso04_all"].get();
        }
        else {std::cerr << "Electron Iso type Error" << std::endl;}
        if(leptons_iso != nullptr)
        {
            //std::cout << "|| At ObjSetup leptons_iso, elecs_iso size: " << leptons_iso->GetSize() << ", " << elecs_iso->GetSize() << std::endl;
        }
        else {std::cerr << "Error: leptons_iso is null after assignment!" << std::endl;}

        /// Electron ID ///
        elecs_scbId = Analysis::intVectors["Electron_cutBased"].get();
        if(TString(lepId).Contains("SCBLoose"))
        {
            //elecIdVariant = intVectors["Electron_cutBased"].get();
            Lepton::eleid_scbcut = 2;
        }
        else if(TString(lepId).Contains("SCBMedium")) {Lepton::eleid_scbcut = 3;}
        else if(TString(lepId).Contains("SCVTight"))  {Lepton::eleid_scbcut = 4;}
        else if(TString(lepId).Contains("SCBVeto"))   {Lepton::eleid_scbcut = 1;}
        else if(TString(lepId).Contains("MVALoose"))
        {
            elecs_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
            Lepton::eleid_scbcut = 2;
        }
        else if(TString(lepId).Contains("MVAMedium"))
        {
            elecs_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WP90"].get();
            Lepton::eleid_scbcut = 3;
        }
        else if(TString(lepId).Contains("MVATight"))
        {
            elecs_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WP80"].get();
            Lepton::eleid_scbcut = 4;
        }
        else if(TString(lepId).Contains("MVAVeto"))
        {
            elecs_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
            Lepton::eleid_scbcut = 1;
        }
        else {std::cerr << "Electron ID Error" << std::endl;}

        if(elecs_mvaId != nullptr)
        {
            //std::cout << "|| elecs_mvaId size: " << elecs_mvaId->GetSize() << std::endl;
        }
        else {std::cerr << "Error: elecs_mvaId is null after assignment!" << std::endl;}
    }
    //////////////////
    // Case of MuEl //
    //////////////////
    else if(TString(Analysis::Decaymode).Contains("muel"))
    {
        /// Muon iso type ///
        if     (TString(MuonIsoType).Contains("PFIsodbeta03")) {muons_iso = Analysis::floatVectors["Muon_pfRelIso03_all"].get();}
        else if(TString(MuonIsoType).Contains("PFIsodbeta04")) {muons_iso = Analysis::floatVectors["Muon_pfRelIso04_all"].get();}
        else {std::cerr << "Muon-electron: Not found Muon Iso type" << std::endl;}
        /// Muon ID ///
        if     (TString(MuonId).Contains("Loose")) {muons_Id = Analysis::boolVectors["Muon_looseId"].get();}
        else if(TString(MuonId).Contains("Tight")) {muons_Id = Analysis::boolVectors["Muon_tightId"].get();}
        else {std::cerr << "Muon-electron: Not found Muon ID type" << std::endl;}

        if(TString(ElecIsoType).Contains("PFIsoRho03"))
        {elecs_iso = Analysis::floatVectors["Electron_pfRelIso03_all"].get();}
        else if(TString(ElecIsoType).Contains("PFIsoRho04"))
        {elecs_iso = Analysis::floatVectors["Electron_pfRelIso04_all"].get();}
        else {std::cerr << "Muon-electron: Not found Electron Iso type" << std::endl;}

        /// Electron ID ///
        elecs_scbId = Analysis::intVectors["Electron_cutBased"].get();
        if(TString(ElecId).Contains("SCBLoose"))
        {
            //elecIdVariant = intVectors["Electron_cutBased"].get();
            Lepton::eleid_scbcut = 2;
        }
        else if(TString(ElecId).Contains("SCBMedium")) {Lepton::eleid_scbcut = 3;}
        else if(TString(ElecId).Contains("SCVTight"))  {Lepton::eleid_scbcut = 4;}
        else if(TString(ElecId).Contains("SCBVeto"))   {Lepton::eleid_scbcut = 1;}
        else if(TString(ElecId).Contains("MVALoose"))
        {
            elecs_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
            Lepton::eleid_scbcut = 2;
        }
        else if(TString(ElecId).Contains("MVAMedium"))
        {
            elecs_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WP90"].get();
            Lepton::eleid_scbcut = 3;
        }
        else if(TString(ElecId).Contains("MVATight"))
        {
            elecs_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WP80"].get();
            Lepton::eleid_scbcut = 4;
        }
        else if(TString(ElecId).Contains("MVAVeto"))
        {
            elecs_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
            Lepton::eleid_scbcut = 1;
        }
        else {std::cerr << "Muon-electron: Not found Electron ID type" << std::endl;}
    }
    /////////////////////////////
    /// Default Error Handling //
    /////////////////////////////
    else {std::cerr << "Error: Unknown Decaymode!" << std::endl;}

    ////////////////////////////////////////////////////////////////////
    // Set leptons for veto lepton (jet cleaning & third lepton veto) //
    ////////////////////////////////////////////////////////////////////
    // Muon information //
    //////////////////////
    if(TString(vetomu_isoid).Contains("PFIsodbeta03"))
    {
        if(Analysis::floatVectors["Muon_pfRelIso03_all"] == nullptr)
        {
            std::cerr << "Error: Muon_pfRelIso03_all branch not initialized!" << std::endl;
            return;
        }
        vetomuon_iso = Analysis::floatVectors["Muon_pfRelIso03_all"].get();
    }
    else if(TString(vetomu_isoid).Contains("PFIsodbeta04"))
    {
        if(Analysis::floatVectors["Muon_pfRelIso04_all"] == nullptr)
        {
            std::cerr << "Error: Muon_pfRelIso04_all branch not initialized!" << std::endl;
            return;
        }
        vetomuon_iso = Analysis::floatVectors["Muon_pfRelIso04_all"].get();
    }
    else
    {
        std::cerr << "Muon Iso type Error" << std::endl;
        return;
    }
    if(vetomuon_iso != nullptr)
    {
        //std::cout << "|| vetomuon_iso size: " << vetomuon_iso->GetSize() << std::endl;
    }
    else {std::cerr << "Error: vetomuon_iso is null after assignment!" << std::endl;}
    /// Muon ID ///
    if     (TString(vetomu_id).Contains( "Loose"  ) ) { vetomuon_id = Analysis::boolVectors["Muon_looseId"].get(); }
    else if(TString(vetomu_id).Contains( "Medium" ) ) { vetomuon_id = Analysis::boolVectors["Muon_mediumId"].get();}
    else if(TString(vetomu_id).Contains( "Tight"  ) ) { vetomuon_id = Analysis::boolVectors["Muon_tightId"].get(); }
    else{ std::cout << "Muon ID Error" << std::endl; }
    
    if(vetomuon_id != nullptr)
    {
        //std::cout << "|| vetomuon_id size: " << vetomuon_id->GetSize() << std::endl;
    }
    else {std::cerr << "Error: vetomuon_id is null after assignment!" << std::endl;}

    //////////////////////////
    // Electron information //
    //////////////////////////
    vetoelec_scbId = Analysis::intVectors["Electron_cutBased"].get();
    if(!Analysis::intVectors["Electron_cutBased"].get()) {std::cerr <<"Error:  Electron_cutBased !! is nullptr! ub Electron information" <<std::endl;}
    if(!vetoelec_scbId) {std::cerr << "Error: vetoelec_scbId is nullptr! ub Electron information " << std::endl;}

    if(TString(vetoelec_id).Contains("SCBLoose"))
    {
        //elecIdVariant = intVectors["Electron_cutBased"].get();
        elevetoid_scbcut = 2;
    }
    else if(TString(vetoelec_id).Contains("SCBMedium")) {elevetoid_scbcut = 3;}
    else if(TString(vetoelec_id).Contains("SCVTight"))  {elevetoid_scbcut = 4;}
    else if(TString(vetoelec_id).Contains("SCBVeto"))   {elevetoid_scbcut = 1;}
    else if(TString(vetoelec_id).Contains("MVALoose"))
    {
        vetoelec_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
        elevetoid_scbcut = 2;
    }
    else if(TString(vetoelec_id).Contains("MVAMedium"))
    {
        vetoelec_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WP90"].get();
        elevetoid_scbcut = 3;
    }
    else if(TString(vetoelec_id).Contains("MVATight"))
    {
        vetoelec_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WP80"].get();
        elevetoid_scbcut = 4;
    }
    else if(TString(vetoelec_id).Contains("MVAVeto"))
    {
        vetoelec_mvaId = Analysis::boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
        elevetoid_scbcut = 1;
    }
    else {std::cerr << "Electron ID in veto selection Error " << vetoelec_id << std::endl;}

    if(TString(vetoelec_isoid).Contains("PFIsoRho03"))
    {
        if(Analysis::floatVectors.find("Electron_pfRelIso03_all") != Analysis::floatVectors.end())
        {
            auto* ptr = Analysis::floatVectors["Electron_pfRelIso03_all"].get();
            if(ptr)
            {
                vetoelec_iso = ptr;
                //std::cout << "|| in sele.... vetoelec_iso .. " << std::endl;
            }
            else {std::cerr << "Error: 'Electron_pfRelIso03_all' is a null unique_ptr." << std::endl;}
        }
        else {std::cerr << "Error: 'Electron_pfRelIso03_all' not found in floatVectors." << std::endl;}
    }
    else if(TString(vetoelec_isoid).Contains("PFIsoRho04")) {std::cerr << "No PFIsoRho04 in NanoAOD..." << std::endl;}
    else {std::cerr << "Electron Iso type Error" << std::endl;}

    if(vetoelec_iso != nullptr)
    {
        //std::cout << "|| vetoelec_iso size: " << vetoelec_iso->GetSize() << std::endl;
    }
    else {std::cerr << "Error: vetoelec_iso is null after assignment!" << std::endl;}

    /////////////////////
    // Jet information //
    /////////////////////
    Jet::jets_pt    = Analysis::floatVectors["Jet_pt"].get(); 
    Jet::jets_eta   = Analysis::floatVectors["Jet_eta"].get(); 
    Jet::jets_phi   = Analysis::floatVectors["Jet_phi"].get();
    Jet::jets_M     = Analysis::floatVectors["Jet_mass"].get();
    Jet::jets_id    = Analysis::intVectors["Jet_jetId"].get();
    Jet::jets_puid  = Analysis::intVectors["Jet_puId"].get();
    Jet::jets_btag  = Analysis::floatVectors["Jet_btagDeepB"].get();

    if(Analysis::RunPeriod.Contains("2016"))
    {   //Ref: https://twiki.cern.ch/twiki/bin/view/CMS/JetID
        if     (JetId == "PFLoose")        { jet_id = 1; } // Loose ID
        else if(JetId == "PFTight")        { jet_id = 3; } // Loose & Tight ID
        else if(JetId == "PFLooseLepVeto") { jet_id = 5; } // Loose + Tight + TightLeptonVeto
        else if(JetId == "PFTightLepVeto") { jet_id = 7; } // Loose + Tight + TightLeptonVeto
        else { std::cout << "|| Jet condition error for 2016!" << std::endl; }
    } 
    else if(Analysis::RunPeriod.Contains("2017") || Analysis::RunPeriod.Contains("2018"))
    {
        if     (JetId == "PFTight")        { jet_id = 2; } // Tight ID
        else if(JetId == "PFTightLepVeto") { jet_id = 6; } // Tight + TightLeptonVeto
        else { std::cout << "|| Jet condition error for 2017/2018!" << std::endl; }
    } 
    else if(Analysis::RunPeriod.Contains("2022") || Analysis::RunPeriod.Contains("2023"))
    {
        if     (JetId == "PFTight")        { jet_id = 2; } // Tight ID
        else if(JetId == "PFTightLepVeto") { jet_id = 6; } // Tight + TightLeptonVeto
        else { std::cout << "|| Jet condition error for 2022/2023!" << std::endl; }
    } 
    else {std::cout << "RunPeriod not recognized!" << std::endl;}
    /// btagging WP ///
    if     ( TString(JetbTag).Contains( "CSVL"  ) )    { Jet::btagcut = 0.244; }
    else if( TString(JetbTag).Contains( "CSVM"  ) )    { Jet::btagcut = 0.679; }
    else if( TString(JetbTag).Contains( "CSVT"  ) )    { Jet::btagcut = 0.898; }
    else if( TString(JetbTag).Contains( "CISVL" ) )    { Jet::btagcut = 0.423; }
    else if( TString(JetbTag).Contains( "CISVM" ) )    { Jet::btagcut = 0.814; }
    else if( TString(JetbTag).Contains( "CISVT" ) )    { Jet::btagcut = 0.941; }
    else if( TString(JetbTag).Contains( "pfCSVV2L" ) ) { Jet::btagcut = 0.5426; }
    else if( TString(JetbTag).Contains( "pfCSVV2M" ) ) { Jet::btagcut = 0.8484; }
    else if( TString(JetbTag).Contains( "pfCSVV2T" ) ) { Jet::btagcut = 0.9535; }
    else if( TString(JetbTag).Contains( "deepCSVL" ) ) { Jet::btagcut = 0.2027; } // only for 2016PreVFP
    else if( TString(JetbTag).Contains( "deepCSVM" ) ) { Jet::btagcut = 0.6001; } // only for 2016PreVFP
    else if( TString(JetbTag).Contains( "deepCSVT" ) ) { Jet::btagcut = 0.8819; } // only for 2016PreVFP
    else if( TString(JetbTag).Contains( "deepJetL" ) ) { Jet::btagcut = 0.0508; } // only for 2016PreVFP
    else if( TString(JetbTag).Contains( "deepJetM" ) ) { Jet::btagcut = 0.2598; } // only for 2016PreVFP
    else if( TString(JetbTag).Contains( "deepJetT" ) ) { Jet::btagcut = 0.6502; } // only for 2016PreVFP
    else {std::cout << "|| bscriminator error !!" << std::endl;}
    /// jet PU ID under 50 GeV ///
    if     (JetPUID == "Fail")   {jet_puid = 0;} // fail all PU ID
    else if(JetPUID == "Loose")  {jet_puid = 1;} // pass loose ID
    else if(JetPUID == "Medium") {jet_puid = 3;} // pass loose & medium ID
    else if(JetPUID == "Tight")  {jet_puid = 7;} // pass all PU IDs
    else {std::cout << "|| Jet PU ID Error" << std::endl;}
    /// MET ///
    MET::met_pt    = Analysis::floatSingles["MET_pt"].get(); 
    MET::met_phi   = Analysis::floatSingles["MET_phi"].get(); 
    //std::cout << " met_pt : " << met_pt << std::endl; 
    //std::cout << " met_phi : " << met_phi << std::endl;

    /// Gen weight ///
    MCSaleFactor::Gen_weight = Analysis::floatSingles["genWeight"].get();
    //Config::Pileup_nPU = Config::floatSingles["Pileup_nPU"].get();
} // end of SetObjectVariable //

std::string Config::removeSubstring(std::string &str, const std::string &keyword)
{
    size_t pos = str.find(keyword);  // find specific keyword ! 
    if(pos != std::string::npos) {str.erase(pos);}
    //std::cout << "removeSubstring: " << str << std::endl;
    return str; 
} // end of removeSubstring //

TLorentzVector Config::createLorentzVector(float pt, float eta, float phi, float mass)
{
    //std::cout << "createLorentzVector " << std::endl;
    /// Check validity of inputs ///
    if(pt < 0 || std::isnan(pt) || std::isnan(eta) || std::isnan(phi) || std::isnan(mass))
    {throw std::runtime_error("Invalid inputs for TLorentzVector");}
    /// Normalize phi to the range [-π, π] ///
    phi = TVector2::Phi_mpi_pi(phi);
    /// Create the TLorentzVector ///
    TLorentzVector lv;// = new TLorentzVector();
    lv.SetPtEtaPhiM(static_cast<double>(pt),
                    static_cast<double>(eta),
                    static_cast<double>(phi),
                    static_cast<double>(mass));
    /// Ensure E^2 >= P^2 (physical validity check) ///
    /*if (lv.E() < lv.P()) {
        std::cout << "pt : " << pt << " eta: " << eta << " phi :" << phi << " mass " << mass << std::endl;
        throw std::runtime_error("Invalid LorentzVector: E < P");
    }*/
    return lv;
} // end of createLorentzVector //

template <typename T>
std::unique_ptr<TTreeReaderValue<T>> Config::DeepCopy(const std::unique_ptr<TTreeReaderValue<T>>& src)
{return std::make_unique<TTreeReaderValue<T>>(*src);} // end of DeepCopy //