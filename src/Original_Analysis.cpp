#include "../interface/Analysis.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Constructor: initialize TTreeReader with TChain and branch list file
Analysis::Analysis(TChain *inputChain, std::string inputName, std::string seDirName, std::string outputName, const std::string &branchListFile, const std::string &configFile, int NumEvt= -1)
    : chain(inputChain), fReader(inputChain), NumEvt(NumEvt), outdir(seDirName), outfile(outputName){
    if (!chain) {
        throw std::runtime_error("Error: Invalid TChain pointer!");
    }
    FileName_ = SetInputFileName(inputName);
    std::cout << "FileName_ : " << FileName_ << std::endl;
    // Load Configuration files //
    std::cout << "configFile : " << configFile << std::endl;  
    std::string confDir = "./configs/";
    std::string confpath = "";
    confpath = confDir+configFile;
    SSBConfReader = new TextReader();
    SSBConfReader->ReadFile(confpath);
    SSBConfReader->ReadVariables();
    SSBConfReader->PrintoutVariables();
    // Initialize branches based on branch list file
    InitBranches(branchListFile);
    cutflowName[0] = "Step_0";
    cutflowName[1] = "Step_1" ;
    cutflowName[2] = "Step_2";
    cutflowName[3] = "Step_3";
    cutflowName[4] = "Step_4";
    cutflowName[5] = "bTagged Jet >= 1";
    cutflowName[6] = "bTagged Jet >= 2";
    cutflowName[7] = "bTagged Jet == 2";
    cutflowName[8] = "Top-Recon.";
    cutflowName[9] = "Top-Pt-Rewight";
 
    pi = TMath::Pi();
    Start();
}

Analysis::~Analysis() {
    // Safely delete the TFile object

    if (fout) {
        fout->Write();
        std::cout << "fout write !!" <<  fout->GetName() << std::endl;
        fout->Close(); // Close the file before deleting
        delete fout;
        fout = nullptr;
        std::cout << "fout successfully deleted." << std::endl;
    }

    // Safely delete the TextReader object
    if (SSBConfReader) {
        delete SSBConfReader;
        SSBConfReader = nullptr;
        std::cout << "SSBConfReader successfully deleted." << std::endl;
    }

    // Add other cleanup as needed for dynamically allocated objects
    std::cout << "Analysis destructor completed." << std::endl;
}

void Analysis::InitBranches(const std::string &branchListFile) {
    std::cout << "branchListFile : " << branchListFile << std::endl;
    std::ifstream infile(branchListFile);
    if (!infile) {
        throw std::runtime_error("Error: Could not open branch list file " + branchListFile);
    }

    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string branchName, objectType, dataType, varType;

        // Parse the line for branch name, object type, data type, and variable type (vector/single)
        if (std::getline(iss, branchName, ',') &&
            std::getline(iss, objectType, ',') &&
            std::getline(iss, dataType, ',') &&
            std::getline(iss, varType, ',')) {

            // Trim whitespace from each part (both ends)
            auto trim = [](std::string &s) {
                s.erase(s.find_last_not_of(" \t\n\r") + 1);
                s.erase(0, s.find_first_not_of(" \t\n\r"));
            };

            trim(branchName);
            trim(objectType);
            trim(dataType);
            trim(varType);

            // Initialize based on data type and variable type
            if (dataType == "Bool_t" && varType == "single") {
                boolSingles[branchName] = std::make_unique<TTreeReaderValue<Bool_t>>(fReader, branchName.c_str());
            } else if (dataType == "Int_t" && varType == "single") {
                intSingles[branchName] = std::make_unique<TTreeReaderValue<Int_t>>(fReader, branchName.c_str());
            } else if (dataType == "UInt_t" && varType == "single") {
                uintSingles[branchName] = std::make_unique<TTreeReaderValue<UInt_t>>(fReader, branchName.c_str());
            } else if (dataType == "Float_t" && varType == "single") {
                floatSingles[branchName] = std::make_unique<TTreeReaderValue<Float_t>>(fReader, branchName.c_str());
            } else if (dataType == "UChar_t" && varType == "single") { // New case for UChar_t single value
                ucharSingles[branchName] = std::make_unique<TTreeReaderValue<UChar_t>>(fReader, branchName.c_str());
            } else if (dataType == "Bool_t" && varType == "vector") {
                boolVectors[branchName] = std::make_unique<TTreeReaderArray<Bool_t>>(fReader, branchName.c_str());
            } else if (dataType == "Int_t" && varType == "vector") {
                intVectors[branchName] = std::make_unique<TTreeReaderArray<Int_t>>(fReader, branchName.c_str());
            } else if (dataType == "UInt_t" && varType == "vector") {
                uintVectors[branchName] = std::make_unique<TTreeReaderArray<UInt_t>>(fReader, branchName.c_str());
            } else if (dataType == "Float_t" && varType == "vector") {
                floatVectors[branchName] = std::make_unique<TTreeReaderArray<Float_t>>(fReader, branchName.c_str());
            } else if (dataType == "UChar_t" && varType == "vector") { // New case for UChar_t vector
                ucharVectors[branchName] = std::make_unique<TTreeReaderArray<UChar_t>>(fReader, branchName.c_str());
            } else {
                std::cerr << "Warning: Unsupported data type or format in branch list: "
                          << dataType << ", " << varType << " for branch " << branchName << std::endl;
            }

        } else {
            std::cerr << "Warning: Invalid branch format in list: " << line << std::endl;
        }
    }
}



void Analysis::SetVariables() {
    std::cout << "Set varibles " << std::endl;

    Lumi = SSBConfReader->GetNumber( "Luminosity" );
    RunPeriod = SSBConfReader->GetText( "RunRange" );
    Decaymode = SSBConfReader->GetText( "Channel" ); // Channel
    XsecTable_ = SSBConfReader->GetText( "XSecTablesName" );
    /// Set Trigger List ///
    num_dleptrig = SSBConfReader->Size( "dileptrigger" );
    num_sleptrig = SSBConfReader->Size( "singleleptrigger" );

    ///
    lepId = SSBConfReader->GetText( "ID_type" );
    lepisotype = SSBConfReader->GetText( "Iso_type" );

    for(int i =0; i < num_dleptrig; ++i)
    {
        std::cout << SSBConfReader->GetText("dileptrigger",i+1) << std::endl;
        std::string tmptrg = SSBConfReader->GetText("dileptrigger",i+1);
        DLtrigName.push_back( removeSubstring( tmptrg, "_v") );
        trigName.push_back( removeSubstring( tmptrg, "_v")  );
    }

    for(int i =0; i < num_sleptrig; ++i)
    {
        std::cout << SSBConfReader->GetText("singleleptrigger",i+1) << std::endl;
        std::string tmptrg = SSBConfReader->GetText("singleleptrigger",i+1);
        SLtrigName.push_back( removeSubstring(tmptrg,"_v") );
        trigName.push_back( removeSubstring(tmptrg,"_v") );
    }

    for (int i = 0; i < trigName.size(); ++i){
        triggerList[trigName[i]] =DeepCopy<bool>( boolSingles[trigName[i]]);
    }

    /// Set Noise filter(MET, events filter) ///
    for (int i = 0; i < SSBConfReader->Size("METFilters"); ++i){
       //std::cout << "METFilters: " << SSBConfReader->GetText("METFilters",i+1) << std::endl;
       std::string tmpnoisefl = SSBConfReader->GetText("METFilters",i+1);
       //std::cout << "METFilters: " <<  tmpnoisefl << std::endl;
       //std::cout << boolSingles[tmpnoisefl] << std::endl;
       if(boolSingles[tmpnoisefl] ==NULL) {std::cout << "Error!!!" << tmpnoisefl << std::endl;}
       noiseFilters[tmpnoisefl] = DeepCopy<bool>( boolSingles[tmpnoisefl]);
       //std::cout << "test "<< std::endl;
    }

    // Lep Infor ( Di Muon or Di Electron )
    lepisotype = SSBConfReader->GetText( "Iso_type" );
    lepId      = SSBConfReader->GetText( "Lep_ID"   );
    
    // Muon Infor for Only MuEl channel
    MuonIsoType = SSBConfReader->GetText( "MuonIso_type" );
    MuonId      = SSBConfReader->GetText( "Muon_ID"      );
    
    // Electron Infor Only for MuEl channel
    ElecIsoType = SSBConfReader->GetText( "ElecIso_type" );
    ElecId      = SSBConfReader->GetText( "Elec_ID"      );  
    // Jet Infor
    JetId      = SSBConfReader->GetText("Jet_ID");
    JetbTag    = SSBConfReader->GetText("Jet_btag");
 
    
    // Kinematic variables for Object //(not MuEl)
    lep_pt     = SSBConfReader->GetNumber( "Lep_pt"  );
    lep_eta    = SSBConfReader->GetNumber( "Lep_eta" );
    lepisocut  = SSBConfReader->GetNumber( "Iso_cut" );
    jet_pt     = SSBConfReader->GetNumber( "Jet_pt"  );
    jet_eta    = SSBConfReader->GetNumber( "Jet_eta" );
    met_cut    = SSBConfReader->GetNumber( "MET_cut" );
 
    // Kinematic variables for Object (MuEl)
    muon_pt     = SSBConfReader->GetNumber( "Muon_pt"     );
    muon_eta    = SSBConfReader->GetNumber( "Muon_eta"    );
    muonisocut  = SSBConfReader->GetNumber( "MuonIso_cut" );
    elec_pt     = SSBConfReader->GetNumber( "Elec_pt"     );
    elec_eta    = SSBConfReader->GetNumber( "Elec_eta"    );
    elecisocut  = SSBConfReader->GetNumber( "ElecIso_cut" );

    muiso_type = SSBConfReader->GetText( "Muiso_for_Jetcleanning" );
    muId       = SSBConfReader->GetText( "Muid_for_Jetcleanning"  );
    eliso_type = SSBConfReader->GetText( "Eliso_for_Jetcleanning" );
    elId       = SSBConfReader->GetText( "Elid_for_Jetcleanning"  );

    dojer = SSBConfReader->GetBool("DoJER");
    JetEnSys    = SSBConfReader->GetText("JESsys"); // Getting Jet Energy Systematic Type ...
    JetResSys   = SSBConfReader->GetText("JERsys"); // Getting Jet Smearing Systematic Type ... 
       
    BTagSFSys    = SSBConfReader->GetText( "BtaggingSFSys"        );
    BTagEffSys   = SSBConfReader->GetText( "BtaggingEffSys"       );
       
    MetSys    = SSBConfReader->GetText("METsys"); // Getting Jet Energy Systematic Type ...
    LepIdSFSys   = SSBConfReader->GetText( "LepIDSFSys"           );
    LepIsoSFSys  = SSBConfReader->GetText( "LepIsoSFSys"          );
    LepRecoSFSys  = SSBConfReader->GetText( "LepRecoSFSys"        );
    LepTrackSFSys  = SSBConfReader->GetText( "LepTrackSFSys"      );
    PDFSys         = SSBConfReader->GetNumber( "PDFSys"           );
    FactRenoSys    = SSBConfReader->GetNumber( "FactRenoSys"      );
    FragmentSys  = SSBConfReader->GetText( "FragmentSys"          );
    DecayTableSys  = SSBConfReader->GetText( "DecayTableSys"      );
    PileUpMCFile = SSBConfReader->GetText( "PileUpMCFileName"     );
    PileUpDATAFile = SSBConfReader->GetText( "PileUpDATAFileName" );
    PileUpSys    = SSBConfReader->GetText( "PileUpSys"            );
    L1PreFireSys    = SSBConfReader->GetText( "L1PreFireSys"            );
    TopPtSys    = SSBConfReader->GetText(  "TopPtSys"             );
 
    ///
    //std::cout << "triggerList : " << triggerList.size()<< std::endl;
    //SetObjectVariable();
}

void Analysis::SetObjectVariable() {
    std::cout << "!!! SetObjectVariable start!!!" << std::endl;
    // Lepton Kinematic //
    // muon //
    muons_pt  = floatVectors["Muon_pt"].get(); 
    muons_eta = floatVectors["Muon_eta"].get(); 
    muons_phi = floatVectors["Muon_phi"].get();
    muons_M   = floatVectors["Muon_mass"].get();

    // electron //
    elecs_pt  = floatVectors["Electron_pt"].get();
    elecs_eta = floatVectors["Electron_eta"].get(); 
    elecs_phi = floatVectors["Electron_phi"].get();
    elecs_M   = floatVectors["Electron_mass"].get();

    ///////////////////////
    /// Case of Di-Muon ///
    ///////////////////////
    if (TString(Decaymode).Contains("dimuon")) {
        if (TString(lepisotype).Contains("PFIsodbeta03")) {
            if (floatVectors["Muon_pfRelIso03_all"] == nullptr) {
                std::cerr << "Error: Muon_pfRelIso03_all branch not initialized!" << std::endl;
                return;
            }
            leptons_iso = floatVectors["Muon_pfRelIso03_all"].get();
            muons_iso = floatVectors["Muon_pfRelIso03_all"].get();
        } else if (TString(lepisotype).Contains("PFIsodbeta04")) {
            if (floatVectors["Muon_pfRelIso04_all"] == nullptr) {
                std::cerr << "Error: Muon_pfRelIso04_all branch not initialized!" << std::endl;
                return;
            }
            leptons_iso = floatVectors["Muon_pfRelIso04_all"].get();
            muons_iso = floatVectors["Muon_pfRelIso04_all"].get();
        } else {
            std::cerr << "Muon Iso type Error" << std::endl;
            return;
        }

        if (leptons_iso != nullptr) {
            std::cout << "At ObjSetup leptons_iso size: " << leptons_iso->GetSize() << std::endl;
        } else {
            std::cerr << "Error: leptons_iso is null after assignment!" << std::endl;
        }
        
        if (muons_iso != nullptr) {
            std::cout << "muons_iso size: " << muons_iso->GetSize() << std::endl;
        } else {
            std::cerr << "Error: muons_iso is null after assignment!" << std::endl;
        }


        /// Muon ID
        if      (TString(lepId).Contains( "Loose"  ) )  { leptons_Id = boolVectors["Muon_looseId"].get();  muons_Id = boolVectors["Muon_looseId"].get(); }
        else if (TString(lepId).Contains( "Medium"  ) ) { leptons_Id = boolVectors["Muon_mediumId"].get(); muons_Id = boolVectors["Muon_mediumId"].get();}
        else if (TString(lepId).Contains( "Tight"  ) )  { leptons_Id = boolVectors["Muon_tightId"].get();  muons_Id = boolVectors["Muon_tightId"].get(); }
        else { std::cout << "Muon ID Error" << std::endl; }

        if (leptons_Id != nullptr) {
            std::cout << "leptons_Id size: " << leptons_Id->GetSize() << std::endl;
        } else {
            std::cerr << "Error: leptons_Id is null after assignment!" << std::endl;
        }
        
        if (muons_Id != nullptr) {
            std::cout << "muons_Id size: " << muons_Id->GetSize() << std::endl;
        } else {
            std::cerr << "Error: muons_Id is null after assignment!" << std::endl;
        }
    }

    ///////////////////////////
    /// Case of Di-Electron ///
    ///////////////////////////
    else if (TString(Decaymode).Contains("dielec")) {
        /// Electron iso type
        if (TString(lepisotype).Contains("PFIsoRho03")) {
            leptons_iso = floatVectors["Electron_pfRelIso03_all"].get();
            elecs_iso = floatVectors["Electron_pfRelIso03_all"].get();
        } else if (TString(lepisotype).Contains("PFIsoRho04")) {
            leptons_iso = floatVectors["Electron_pfRelIso03_all"].get();
            elecs_iso = floatVectors["Electron_pfRelIso03_all"].get();
            std::cerr << "No PFIsoRho04 in NanoAOD..." << std::endl;
        } else {
            std::cerr << "Electron Iso type Error" << std::endl;
        }

        /// Electron ID
        elecs_scbId = intVectors["Electron_cutBased"].get();
        if (TString(lepId).Contains("SCBLoose")) {
            //elecIdVariant = intVectors["Electron_cutBased"].get();
            eleid_scbcut = 2;
        } else if (TString(lepId).Contains("SCBMedium")) {
            eleid_scbcut = 3;
        } else if (TString(lepId).Contains("SCVTight")) {
            eleid_scbcut = 4;
        } else if (TString(lepId).Contains("SCBVeto")) {
            eleid_scbcut = 1;
        } else if (TString(lepId).Contains("MVALoose")) {
            elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
            eleid_scbcut = 2;
        } else if (TString(lepId).Contains("MVAMedium")) {
            elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP90"].get();
            eleid_scbcut = 3;
        } else if (TString(lepId).Contains("MVATight")) {
            elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP80"].get();
            eleid_scbcut = 4;
        } else if (TString(lepId).Contains("MVAVeto")) {
            elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
            eleid_scbcut = 1;
        } else {
            std::cerr << "Electron ID Error" << std::endl;
        }

    }

    ////////////////////
    /// Case of MuEl ///
    ////////////////////
    else if (TString(Decaymode).Contains("muel")) {
        /// Muon iso type
        if (TString(MuonIsoType).Contains("PFIsodbeta03")) {
            muons_iso = floatVectors["Muon_pfRelIso03_all"].get();
        } else if (TString(MuonIsoType).Contains("PFIsodbeta04")) {
            muons_iso = floatVectors["Muon_pfRelIso04_all"].get();
        } else {
            std::cerr << "Muon Iso type Error in MuEl case" << std::endl;
        }

        /// Muon ID
        if (TString(MuonId).Contains("Loose")) {
            muons_Id = boolVectors["Muon_tightId"].get();
        } else if (TString(MuonId).Contains("Tight")) {
            muons_Id = boolVectors["Muon_isTight"].get();
        } else {
            std::cerr << "Muon ID Error in MuEl case" << std::endl;
        }

        if (TString(ElecIsoType).Contains("PFIsoRho03")) {
            elecs_iso = floatVectors["Electron_pfRelIso03_all"].get();
        } else if (TString(ElecIsoType).Contains("PFIsoRho04")) {
            elecs_iso = floatVectors["Electron_pfRelIso03_all"].get();
            std::cerr << "No PFIsoRho04 in NanoAOD..." << std::endl;
        } else {
            std::cerr << "Electron Iso type Error" << std::endl;
        }

        /// Electron ID
        elecs_scbId = intVectors["Electron_cutBased"].get();
        if (TString(ElecId).Contains("SCBLoose")) {
            //elecIdVariant = intVectors["Electron_cutBased"].get();
            eleid_scbcut = 2;
        } else if (TString(ElecId).Contains("SCBMedium")) {
            eleid_scbcut = 3;
        } else if (TString(ElecId).Contains("SCVTight")) {
            eleid_scbcut = 4;
        } else if (TString(ElecId).Contains("SCBVeto")) {
            eleid_scbcut = 1;
        } else if (TString(ElecId).Contains("MVALoose")) {
            elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
            eleid_scbcut = 2;
        } else if (TString(ElecId).Contains("MVAMedium")) {
            elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP90"].get();
            eleid_scbcut = 3;
        } else if (TString(ElecId).Contains("MVATight")) {
            elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP80"].get();
            eleid_scbcut = 4;
        } else if (TString(ElecId).Contains("MVAVeto")) {
            elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
            eleid_scbcut = 1;
        } else {
            std::cerr << "Electron ID Error" << std::endl;
        }

    }

    /////////////////////////////
    /// Default Error Handling //
    /////////////////////////////
    else {
        std::cerr << "Error: Unknown Decaymode!" << std::endl;
    }

    //////////////////////////////////////////////////////////////////////
    /// Set leptons for veto lepton (jet cleaning & third lepton veto) ///
    //////////////////////////////////////////////////////////////////////
    ////////////////////////
    /// Muon information ///
    ////////////////////////
    if (TString(muiso_type).Contains("PFIsodbeta03")) {
        if (floatVectors["Muon_pfRelIso03_all"] == nullptr) {
            std::cerr << "Error: Muon_pfRelIso03_all branch not initialized!" << std::endl;
            return;
        }
        muonsveto_iso = floatVectors["Muon_pfRelIso03_all"].get();
    } else if (TString(muiso_type).Contains("PFIsodbeta04")) {
        if (floatVectors["Muon_pfRelIso04_all"] == nullptr) {
            std::cerr << "Error: Muon_pfRelIso04_all branch not initialized!" << std::endl;
            return;
        }
        muonsveto_iso = floatVectors["Muon_pfRelIso04_all"].get();
    } else {
        std::cerr << "Muon Iso type Error" << std::endl;
        return;
    }

    
    if (muonsveto_iso != nullptr) {
        std::cout << "muonsveto_iso size: " << muonsveto_iso->GetSize() << std::endl;
    } else {
        std::cerr << "Error: muonsveto_iso is null after assignment!" << std::endl;
    }


    /// Muon ID
    if      (TString(muId).Contains( "Loose"  ) )  { muonsveto_Id = boolVectors["Muon_looseId"].get(); }
    else if (TString(muId).Contains( "Medium"  ) ) { muonsveto_Id = boolVectors["Muon_mediumId"].get();}
    else if (TString(muId).Contains( "Tight"  ) )  { muonsveto_Id = boolVectors["Muon_tightId"].get(); }
    else { std::cout << "Muon ID Error" << std::endl; }
    
    if (muonsveto_Id != nullptr) {
        std::cout << "muonsveto_Id size: " << muonsveto_Id->GetSize() << std::endl;
    } else {
        std::cerr << "Error: muonsveto_Id is null after assignment!" << std::endl;
    }


    ////////////////////////////
    /// Electron information ///
    ////////////////////////////
    elecsveto_scbId = intVectors["Electron_cutBased"].get();
    if (!intVectors["Electron_cutBased"].get()) {std::cerr <<"Error:  Electron_cutBased !! is nullptr! ub Electron information" <<std::endl;}
    if (!elecsveto_scbId) std::cerr << "Error: elecsveto_scbId is nullptr! ub Electron information " << std::endl;

    if (TString(elId).Contains("SCBLoose")) {
        //elecIdVariant = intVectors["Electron_cutBased"].get();
        elevetoid_scbcut = 2;
    } else if (TString(elId).Contains("SCBMedium")) {
        elevetoid_scbcut = 3;
    } else if (TString(elId).Contains("SCVTight")) {
        elevetoid_scbcut = 4;
    } else if (TString(elId).Contains("SCBVeto")) {
        elevetoid_scbcut = 1;
    } else if (TString(elId).Contains("MVALoose")) {
        elecsveto_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
        elevetoid_scbcut = 2;
    } else if (TString(elId).Contains("MVAMedium")) {
        elecsveto_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP90"].get();
        elevetoid_scbcut = 3;
    } else if (TString(elId).Contains("MVATight")) {
        elecsveto_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP80"].get();
        elevetoid_scbcut = 4;
    } else if (TString(elId).Contains("MVAVeto")) {
        elecsveto_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
        elevetoid_scbcut = 1;
    } else {
        std::cerr << "Electron ID in veto selection Error " << elId << std::endl;
    }

    if (TString(eliso_type).Contains("PFIsoRho03")) {
        if (floatVectors.find("Electron_pfRelIso03_all") != floatVectors.end()) {
            auto* ptr = floatVectors["Electron_pfRelIso03_all"].get();
            if (ptr) {
                elecsveto_iso = ptr;
                std::cout << "in sele.... elecsveto_iso .. " << std::endl;
            } else {
                std::cerr << "Error: 'Electron_pfRelIso03_all' is a null unique_ptr." << std::endl;
            }
        } else {
            std::cerr << "Error: 'Electron_pfRelIso03_all' not found in floatVectors." << std::endl;
        }
    } else if (TString(eliso_type).Contains("PFIsoRho04")) {
        std::cerr << "No PFIsoRho04 in NanoAOD..." << std::endl;
    } else {
        std::cerr << "Electron Iso type Error" << std::endl;
    }

    if (elecsveto_iso != nullptr) {
        std::cout << "elecsveto_iso size: " << elecsveto_iso->GetSize() << std::endl;
    } else {
        std::cerr << "Error: elecsveto_iso is null after assignment!" << std::endl;
    }


    //////////////////////
    /// Set Jet object ///
    //////////////////////
     
    jets_pt  = floatVectors["Jet_pt"].get(); 
    jets_eta = floatVectors["Jet_eta"].get(); 
    jets_phi = floatVectors["Jet_phi"].get();
    jets_M   = floatVectors["Jet_mass"].get();
    jets_Id = boolVectors["Jet_jetId"].get();

    if (RunPeriod.Contains("2016")) {
        if      (JetId == "PFLoose") { jet_id = 1; } // Loose ID
        else if (JetId == "PFTight") { jet_id = 3; } // Loose + Tight ID
        else if (JetId == "PFLooseLepVeto") { jet_id = 5; } // Loose + Tight + TightLeptonVeto
        else if (JetId == "PFTightLepVeto") { jet_id = 7; } // Loose + Tight + TightLeptonVeto
        else { std::cout << "Jet condition error for 2016!" << std::endl; }
    } 
    else if (RunPeriod.Contains("2017") || RunPeriod.Contains("2018")) {
        if      (JetId == "PFTight") { jet_id = 2; } // Tight ID
        else if (JetId == "PFTightLepVeto") { jet_id = 6; } // Tight + TightLeptonVeto
        else { std::cout << "Jet condition error for 2017/2018!" << std::endl; }
    } 
    else if (RunPeriod.Contains("2022") || RunPeriod.Contains("2023")) {
        if      (JetId == "PFTight") { jet_id = 2; } // Tight ID
        else if (JetId == "PFTightLepVeto") { jet_id = 6; } // Tight + TightLeptonVeto
        else { std::cout << "Jet condition error for 2022/2023!" << std::endl; }
    } 
    else {
        std::cout << "RunPeriod not recognized!" << std::endl;
    }

    /// btagging WP ///
    if      ( TString(JetbTag).Contains( "CSVL"  ) )    { bdisccut = 0.244; }
    else if ( TString(JetbTag).Contains( "CSVM"  ) )    { bdisccut = 0.679; }
    else if ( TString(JetbTag).Contains( "CSVT"  ) )    { bdisccut = 0.898; }
    else if ( TString(JetbTag).Contains( "CISVL" ) )    { bdisccut = 0.423; }
    else if ( TString(JetbTag).Contains( "CISVM" ) )    { bdisccut = 0.814; }
    else if ( TString(JetbTag).Contains( "CISVT" ) )    { bdisccut = 0.941; }
    else if ( TString(JetbTag).Contains( "pfCSVV2L" ) ) { bdisccut = 0.5426; }
    else if ( TString(JetbTag).Contains( "pfCSVV2M" ) ) { bdisccut = 0.8484; }
    else if ( TString(JetbTag).Contains( "pfCSVV2T" ) ) { bdisccut = 0.9535; }
    else if ( TString(JetbTag).Contains( "deepCSVL" ) ) { bdisccut = 0.2027; }
    else if ( TString(JetbTag).Contains( "deepCSVM" ) ) { bdisccut = 0.6001; }
    else if ( TString(JetbTag).Contains( "deepCSVT" ) ) { bdisccut = 0.8819; }
    else if ( TString(JetbTag).Contains( "deepJetL" ) ) { bdisccut = 0.0508; }
    else if ( TString(JetbTag).Contains( "deepJetM" ) ) { bdisccut = 0.2598; }
    else if ( TString(JetbTag).Contains( "deepJetT" ) ) { bdisccut = 0.6502; }
    else { std::cout << "bscriminator error !!" << std::endl; }

    /// MET ///
    met_pt  = floatSingles["MET_pt"].get(); 
    met_phi  = floatSingles["MET_phi"].get(); 
//    std::cout << " met_pt : " << met_pt << std::endl; 
//    std::cout << " met_phi : " << met_phi << std::endl; 

    std::cout << "End of SetObjectVariable !" << std::endl;
}



// Event loop function
void Analysis::Loop() {
    if (!chain) {
        std::cerr << "Error: TChain is null!" << std::endl;
        return;
    }

    Long64_t nEntries = chain->GetEntries();
    fReader.Restart(); // Reset the reader to the beginning
    evt_weight_ = 1;
    MCSFApply();
    for (Long64_t ientry = 0; ientry < ((NumEvt == -1 || NumEvt > nEntries) ? nEntries : NumEvt); ++ientry) {
    // 루프 내용


        if (!fReader.Next()) { // Move to the next entry
            std::cerr << "Error: Failed to read entry " << ientry << std::endl;
            break;
        }
        SetObjectVariable(); //
        LeptonSelector(); 

        LeptonOrder();
        JetSelector(); 
        METDefiner();

        //if (i > 100) break; //%lld supports Long64_t
        if (ientry % 10000 == 0) {
            printf("Event %lld\n", ientry); //%lld supports Long64_t
        }


        // Noise (MET) Filter //
        if ( METFilterAPP() == false ) {continue;}
        // Trigger Requirement //
        if ( Trigger() == false ) {continue;}
  
        // Good Primary Vertex Selection //
        if (intSingles["PV_npvsGood"] && **intSingles["PV_npvsGood"] < 1) {
            continue;
        }
        std::cout << "**intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] <<std::endl;
        std::cout << "evt_weight_ : " << evt_weight_ << std::endl;
        FillHisto( h_Num_PV[0]    , **intSingles["PV_npvsGood"] , evt_weight_ );
        FillHisto( h_DiLepMass[0], ( (Lep1)+(Lep2) ).M(), evt_weight_ );
        FillHisto( h_Lep1pt[0] , (Lep1).Pt()  , evt_weight_ );
        FillHisto( h_Lep1eta[0], (Lep1).Eta() , evt_weight_ );
        FillHisto( h_Lep1phi[0], (Lep1).Phi() , evt_weight_ );
        FillHisto( h_Lep2pt[0] , (Lep2).Pt()  , evt_weight_ );
        FillHisto( h_Lep2eta[0], (Lep2).Eta() , evt_weight_ );
        FillHisto( h_Lep2phi[0], (Lep2).Phi() , evt_weight_ );
        FillHisto( h_METpt[0]   , Met.Pt()  , evt_weight_ );
        FillHisto( h_METphi[0]  , Met.Phi()  , evt_weight_ );
        FillHisto( h_Num_Jets[0]  , v_jet_idx.size(), evt_weight_ );
        //FillHisto( h_Num_bJets[1], nbtagged, evt_weight_ );

        if (LeptonsPtAddtional() == false ) {continue;}
        if ( DiLeptonMassCut() == false ) {continue;}
        /// Step 1 ///
        FillHisto( h_DiLepMass[1], ( (Lep1)+(Lep2) ).M(), evt_weight_ );
        FillHisto( h_Num_PV[1], num_pv, evt_weight_ );
        FillHisto( h_Lep1pt[1] , (Lep1).Pt()  , evt_weight_ );
        FillHisto( h_Lep1eta[1], (Lep1).Eta() , evt_weight_ );
        FillHisto( h_Lep1phi[1], (Lep1).Phi() , evt_weight_ );
        FillHisto( h_Lep2pt[1] , (Lep2).Pt()  , evt_weight_ );
        FillHisto( h_Lep2eta[1], (Lep2).Eta() , evt_weight_ );
        FillHisto( h_Lep2phi[1], (Lep2).Phi() , evt_weight_ );
        FillHisto( h_METpt[1]   , Met.Pt()  , evt_weight_ );
        FillHisto( h_METphi[1]  , Met.Phi()  , evt_weight_ );
        FillHisto( h_Num_Jets[1]  , v_jet_idx.size(), evt_weight_ );
 
        
        //FillHisto( h_Num_bJets[1], nbtagged, evt_weight_ );
       

    }// end of event iteration //
    std::cout << "End Loop !!" << std::endl;
}


template <typename T>
std::unique_ptr<TTreeReaderValue<T>> Analysis::DeepCopy(const std::unique_ptr<TTreeReaderValue<T>>& src) {                       
    return std::make_unique<TTreeReaderValue<T>>(*src);                                                                             
}

std::string Analysis::removeSubstring(std::string &str, const std::string &keyword) {
    size_t pos = str.find(keyword);  // find specific keyword ! 
    if (pos != std::string::npos) {
        str.erase(pos);  // 
    }
    std::cout << "str : " << str << std::endl;
    return str; 
}

bool Analysis::METFilterAPP()
{
   bool metfilt_ = true;
   for ( const auto &pair : noiseFilters )
   { 
      TString METFiltName = pair.first;
//      if ( !(TString(FileName_).Contains("Data")) && 
//            TString(METFiltName).Contains("Flag_eeBadScFilter")  ) {continue;}
            if ( !(**pair.second) ){ metfilt_ = false; }
   }
   return metfilt_;
   //return true;
}

// Trigger Requirement Function
bool Analysis::SelTrigger(std::vector<std::string> v_sel)
{
    std::string trgName = "";
    int ptrigindex = 0;
    bool passtrig_ = false;

    for (int j = 0; j < v_sel.size(); j++)
    {
        trgName = v_sel[j];

        // trgName cout in triggerList 
        auto it = triggerList.find(trgName);
        if (it != triggerList.end() && it->second) {
            if (**(it->second)) {ptrigindex++;}
        } else {
            std::cerr << "Error: Trigger " << trgName << " not found in triggerList." << std::endl;
        }
    }

    if (ptrigindex > 0) {
        passtrig_ = true;
    }
    return passtrig_;
}

bool Analysis::Trigger()
{
   /// Variable for Trigger Function
   bool singleTrig_ = false;
   bool doubleTrig_ = false;
   bool ispassselTrig_ = false;
   bool ispassvetoTrig_ = false;
   bool trigpass = false;
   std::vector<std::string> seltrigName;
   std::vector<std::string> vetotrigName;
 
   if (!TString(FileName_).Contains( "Data" ) ) { 
       //trigpass = true; return trigpass;
       ispassvetoTrig_ = false; 
       seltrigName = trigName;
       trigpass = SelTrigger(seltrigName);
   }
   else {
      // Set veto trigger & selected trigger //
      // Channel Index //
      if ( RunPeriod.Contains("2018") )/// Only 2018, SingleEG and Double EG combined 
      {
         if (TString(Decaymode).Contains("dimuon")){ // Dimuon // 
            if ( TString(FileName_).Contains( "Single") ) {
               seltrigName = SLtrigName;
               vetotrigName = DLtrigName;
               ispassselTrig_ = SelTrigger(seltrigName);
               ispassvetoTrig_ = SelTrigger(vetotrigName);
            }
            else if( TString(FileName_).Contains( "Double") ) {
               seltrigName = DLtrigName;
               vetotrigName = SLtrigName;
               ispassselTrig_ = SelTrigger(seltrigName);
               ispassvetoTrig_ = SelTrigger(vetotrigName);
            }
            else { std::cout << "Check out FileName_ in Trigger ()" << std::endl;}
         }
         else if(TString(Decaymode).Contains("muel")) {
            if ( TString(FileName_).Contains( "MuonEG") ) {
               seltrigName = DLtrigName;
               vetotrigName = SLtrigName;
               ispassselTrig_ = SelTrigger(seltrigName);
               ispassvetoTrig_ = SelTrigger(vetotrigName);
            }
            else if( TString(FileName_).Contains( "EGamma") ) {
               seltrigName = SLtrigName;
               vetotrigName = DLtrigName;
               ispassselTrig_ = SelTrigger(seltrigName);
               ispassvetoTrig_ = SelTrigger(vetotrigName);
            }
            else { std::cout << "Check out FileName_ in Trigger ()" << std::endl;}
         }
         else if(TString(Decaymode).Contains("dielec")) { 
            //seltrigName = DLtrigName;
            //vetotrigName = SLtrigName;
            if( TString(FileName_).Contains( "EGamma") ) {
               singleTrig_ = false;
               doubleTrig_ = false;
               singleTrig_ = SelTrigger(SLtrigName);
               doubleTrig_ = SelTrigger(DLtrigName);
               if (singleTrig_&&(!doubleTrig_)) {trigpass = true;}
               else if ((!singleTrig_)&&(doubleTrig_)) {trigpass = true;}
               else {trigpass = false;}
               ispassvetoTrig_ = false;
            }
         }
         else { std::cout << "Check out Decaymode in 2018-Trigger()" << std::endl; }
      }
      else {/// 2016 (AVP, NonAPV) && 2017 
         if ( TString(FileName_).Contains( "Single") ) {
            seltrigName = SLtrigName;
            vetotrigName = DLtrigName;
            //std::cout << "selected!! " << std::endl;
            ispassselTrig_ = SelTrigger(seltrigName);
            //std::cout << "ispassselTrig_ : " << ispassselTrig_ << std::endl;
            //std::cout << "veto !! " << std::endl;
            ispassvetoTrig_ = SelTrigger(vetotrigName);
            //std::cout << "ispassvetoTrig_ : " << ispassvetoTrig_ << std::endl;
            trigpass = ispassselTrig_;
         }
         else if( TString(FileName_).Contains( "Double") || TString(FileName_).Contains( "MuonEG")) {
            //std::cout << "SK : FileName_ "  << std::endl;
            seltrigName = DLtrigName;
            vetotrigName = SLtrigName;
            //std::cout << "SK : seltrigName:" << seltrigName  << std::endl;
            //std::cout << "SK : vetotrigName:" << vetotrigName  << std::endl;
            ispassselTrig_ = SelTrigger(trigName);
            //ispassvetoTrig_ = SelTrigger(vetotrigName);
            trigpass = ispassselTrig_;
         }
         else { std::cout << "Check out FileName_ in Trigger ()" << std::endl;}

      }

      if (ispassvetoTrig_ == true) {trigpass = false;}
   }
   //cout << "ispassselTrig_ : " << ispassselTrig_ << " ispassvetoTrig_ : " << ispassvetoTrig_ << "trigger : " << trigpass << endl; 
   return trigpass;
}
TString Analysis::SetInputFileName( std::string inname )
{  
   TString inputName = inname;
   char unsco_ = '_';
   Size_t unscoIndex = inputName.Last(unsco_);
   inputName.Remove(unscoIndex, inputName.Length());
   if (inputName.Contains("/")){
      unscoIndex = inputName.Last('/');
      inputName.Remove(0,unscoIndex+1);
   }
   //cout << "inputName : " << inputName << endl;
   return inputName; 
}

void Analysis::MCSF()
{
    if (FileName_.Contains("Data")||FileName_.Contains("Single")||FileName_.Contains("EG")){ mc_sf_ = 1.; return; }
    /// Open Xsec Tables ///
    FILE *xsecs_;
    char sampleName[1000];
    double xsec_ = -1.;
    double br_ = -1.; 
    int totalevt_  = -1.; 
    int positive_  = -1.; 
    int negative_  = -1.; 
    int posi_nega_ = -1.; 
    std::string xsec_dir= "./xsecAndsample/";
    std::string xsec_filePath = xsec_dir+ XsecTable_.Data();
    //cout << "xsec_filePath : " << xsec_filePath << std::endl;
    /// SampleName | TotalEvt | Positive+Negative | Xsection | Branching Fraction |
    xsecs_ = fopen(xsec_filePath.c_str(),"r");
    std::map<std::string, int> m_sam_totalevt;
    std::map<std::string, double> m_sam_xsec;
    std::map<std::string, double> m_sam_br;
    std::map<std::string, int> m_sam_positive;
    std::map<std::string, int> m_sam_negative;
    std::map<std::string, int> m_sam_posi_nega;
    if (xsecs_!=NULL) 
    { 
       //cout << "Load Xsection Table!" << std::endl;
       while (fscanf(xsecs_, "%s %d %d %d %d %lf %lf\n", sampleName, &totalevt_, &positive_, &negative_, &posi_nega_, &xsec_, &br_ ) != EOF)
       {
          std::cout 
          << "sampleName : " << sampleName << " totalevt_ : " << totalevt_
          << " positive_ " << positive_ << " negative_ : " << negative_ 
          << " posi_nega_ " << posi_nega_ << " xsec_ : " << xsec_ 
          << " br_ " << br_ 
          << std::endl;
          m_sam_totalevt[sampleName] = totalevt_;
          m_sam_positive[sampleName] = positive_;
          m_sam_negative[sampleName] = negative_;
          m_sam_posi_nega[sampleName] = posi_nega_;
          m_sam_xsec[sampleName] = xsec_;
          m_sam_br[sampleName] = br_;
       }
       fclose(xsecs_);
    }
    else {std::cout << "No xsec_filePath !!!" << xsec_filePath << std::endl;return;} 
    std::cout << "Lumi : " << Lumi << std::endl;
    double lumi = Lumi/1000000;
    auto it = m_sam_xsec.find(FileName_.Data());
    if (it !=  m_sam_xsec.end()){
        std::cout << "SK Key " << FileName_.Data() << " found in the std::map."<< std::endl;
        mc_sf_ = (m_sam_xsec[FileName_.Data()]*m_sam_br[FileName_.Data()]*lumi)/m_sam_posi_nega[FileName_.Data()];
        std::cout << "mc_sf_ " << mc_sf_ << std::endl;
    }
    else {
        mc_sf_ =1.;
        std::cout << "Key " << FileName_.Data() << " not found in the std::map. mc sf is 1" << mc_sf_ << std::endl;
    }
    return;
}


// Apply MC SF To Event //
void Analysis::MCSFApply()
{
    MCSF();
    evt_weight_beforemcsf_ =1; // Initailize evt_weight_beforemcsf_ //
    evt_weight_beforemcsf_ = evt_weight_; // keep event weight //
     
    if ( !TString(FileName_).Contains( "Data") ){ evt_weight_ = evt_weight_*mc_sf_; std::cout <<  "evt_weight_ : "<< evt_weight_ << " mc_sf_ : "  << mc_sf_  << std::endl; 
    } // apply MC scale factor // 
    else {evt_weight_ = 1;}
    std::cout << "SK MC SF : " << evt_weight_<< std::endl;
}

// Safely create a TLorentzVector
TLorentzVector Analysis::createLorentzVector(float pt, float eta, float phi, float mass) {
    //std::cout << "createLorentzVector " << std::endl;
    // Check validity of inputs
    if (pt < 0 || std::isnan(pt) || std::isnan(eta) || std::isnan(phi) || std::isnan(mass)) {
        throw std::runtime_error("Invalid inputs for TLorentzVector");
    }

    // Normalize phi to the range [-π, π]
    phi = TVector2::Phi_mpi_pi(phi);

    // Create the TLorentzVector
    TLorentzVector lv;// = new TLorentzVector();
    lv.SetPtEtaPhiM(static_cast<double>(pt),
                    static_cast<double>(eta),
                    static_cast<double>(phi),
                    static_cast<double>(mass));

    // Ensure E^2 >= P^2 (physical validity check)
    /*if (lv.E() < lv.P()) {
        std::cout << "pt : " << pt << " eta: " << eta << " phi :" << phi << " mass " << mass << std::endl;
        throw std::runtime_error("Invalid LorentzVector: E < P");
    }*/

    return lv;
}

void Analysis::LeptonSelector()
{
    // Clear output vectors
    v_lepton_idx.clear();
    v_muon_idx.clear();
    v_electron_idx.clear();
    v_lep_idx_temp.clear();
    v_muon_idx_temp.clear();
    v_electron_idx_temp.clear();

    // Prepare corrected collections
    MakeMuonCollection();
    MakeElecCollection();

    // Lambda functions for common checks
    auto passKinematicCuts = [](float _pt, float _eta, float ptCut, float etaCut) {
      return _pt >= ptCut && fabs(_eta) <= etaCut;
    };

    auto passIsolation = [](float iso, float isoCut) {
      return iso <= isoCut;
    };

    auto passId = [](bool id) {
      return id;
    };

    auto elecSCBId = [](int id, int idcut) {
      return id >= idcut; 
    };
    auto elecCharge = [](int id) {
      return id > 0; 
    };

    // Dimuon case
    if (TString(Decaymode).Contains("dimuon")) {
        Int_t nmu = muons_iso->GetSize();
        for (int i = 0; i < nmu; ++i) {

            if (!passIsolation(muons_iso->At(i), lepisocut) || 
                !passKinematicCuts(muons_pt->At(i), muons_eta->At(i),lep_pt, lep_eta) || 
                !passId(muons_Id->At(i))) {
                continue;
            }

            v_lep_idx_temp.push_back(i);
            v_lep_idx_temp.push_back(i);
            if (v_lep_idx_temp.size() == 1 && v_lepton_idx.empty()) {
                v_lepton_idx.push_back(i);
            } else if (v_lep_idx_temp.size() > 0 &&
                        (*intVectors["Muon_charge"])[v_lep_idx_temp[0]] != (*intVectors["Muon_charge"])[i]){

                v_lepton_idx.push_back(i);
            }
        }
        // Assign final indices to v_muon_idx
        v_muon_idx = v_lepton_idx;

    }
    // Dielectron case
    else if (TString(Decaymode).Contains("dielec")) {
        //Int_t nel = Elec->GetEntriesFast();
        Int_t nel = elecs_pt->GetSize();
        for (int i = 0; i < nel; ++i) {
            //Ele_lep_sel = elecs.at(i);


            if (!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), lep_pt, lep_eta) ||
                !elecSCBId(elecs_scbId->At(i), eleid_scbcut) ||
                (fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) > 1.4442 &&
                 fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) < 1.566) ||
                !(elecCharge((*intVectors["Electron_tightCharge"])[i]))) {
                continue;
            }

            v_lep_idx_temp.push_back(i);
            if (v_lep_idx_temp.size() == 1 && v_lepton_idx.empty()) {
                v_lepton_idx.push_back(i);
                v_electron_idx.push_back(i);
            } else if (v_lep_idx_temp.size() > 0 &&
              ((*intVectors["Electron_charge"])[v_lep_idx_temp[0]] != (*intVectors["Electron_charge"])[i])){
                v_lepton_idx.push_back(i);
                v_electron_idx.push_back(i);
            }
        }
    }
    // MuEl case
    else if (TString(Decaymode).Contains("muel")) {
        // Process muons
        //Int_t nmu = Muon->GetEntriesFast();
        Int_t nmu = muons_iso->GetSize();
        for (int i = 0; i < nmu; ++i) {

            if (!passIsolation(muons_iso->At(i), muonisocut) || 
                //!passKinematicCuts(*Mu_lep_sel, muon_pt, muon_eta) || 
                !passKinematicCuts(muons_pt->At(i), muons_eta->At(i),lep_pt, lep_eta) || 
                !passId(muons_Id->At(i))) {
                continue;
            }

            v_muon_idx.push_back(i);
        }

        // Process electrons

        Int_t nel = elecs_pt->GetSize();
        for (int i = 0; i < nel; ++i) {
        
          if (!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), lep_pt, lep_eta) || 
          !elecSCBId(elecs_scbId->At(i), eleid_scbcut) || 
          (fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) > 1.4442 && 
          fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) < 1.566) || 
          !(elecCharge((*intVectors["Electron_tightCharge"])[i])) || // Added missing logical operator
          !(*boolVectors["Electron_convVeto"])[i]) { 
              continue;
          }   
            v_electron_idx_temp.push_back(i);
        }


        if (!v_muon_idx.empty() && !v_electron_idx_temp.empty()) {
            for (int i : v_electron_idx_temp) {
                if ((*intVectors["Muon_charge"])[v_lep_idx_temp[0]] != (*intVectors["Electron_charge"])[i]) { // Added missing closing parenthesis
                    v_electron_idx.push_back(i);
                }
            }
        }

        // Validate charge matching
        if (!v_muon_idx.empty() && !v_electron_idx.empty() &&
           (*intVectors["Muon_charge"])[v_muon_idx[0]] == (*intVectors["Electron_charge"])[v_electron_idx[0]]) {
           if ((*intVectors["Muon_charge"])[v_muon_idx_temp[0]] != (*intVectors["Electron_charge"])[v_electron_idx[0]]) {
              std::cerr << "MuEl selection wrong!" << std::endl;
           }
        }
    }
    else {
        std::cerr << "Lepton Selection error" << std::endl;
    }
    std::cout << " ---sk 333333--- " << std::endl;
    // Muon Selection
    Int_t nmu = muonsveto_iso->GetSize();

    for (int imu = 0; imu < nmu; ++imu) {
        if (!passIsolation(muons_iso->At(imu), lepisocut) || 
            !passKinematicCuts(muons_pt->At(imu), muons_eta->At(imu), lep_pt, lep_eta) ||
            !passId(muons_Id->At(imu))) {
            continue;
        }

        muonsveto.push_back(muons.at(imu));

    }
    std::cout << " ---sk 444444--- " << std::endl;

    if (elecsveto_iso == nullptr) {
        std::cerr << "Error: elecsveto_iso is null!" << std::endl;
        return; // 또는 continue;
    }

    // Electron Selection
    std::cout << "elecs_iso --- " << elecsveto_iso->GetSize() << ", elecs.size: " << elecs.size() << std::endl;

    Int_t nel = elecsveto_iso->GetSize();
    if (nel != elecs.size()) {
        std::cerr << "Error: Size mismatch! nel: " << nel << ", elecs.size: " << elecs.size() << std::endl;
        return; // 또는 continue;
    }

    std::cout << "---sk 55555---" << std::endl;

    if (!elecs_pt || !elecs_eta || !floatVectors["Electron_deltaEtaSC"] || !intVectors["Electron_tightCharge"]) {
        if (!elecs_pt) std::cerr << "Error: elecs_pt is null!" << std::endl;
        if (!elecs_eta) std::cerr << "Error: elecs_eta is null!" << std::endl;
        if (!floatVectors["Electron_deltaEtaSC"]) std::cerr << "Error: floatVectors['Electron_deltaEtaSC'] is null!" << std::endl;
        if (!intVectors["Electron_tightCharge"]) std::cerr << "Error: intVectors['Electron_tightCharge'] is null!" << std::endl;
        return;
    }


    for (int iel = 0; iel < nel; ++iel) {
        std::cout << "Processing electron index: " << iel << std::endl;

        // 벡터 및 크기 확인
        if (!elecs_pt || !elecs_eta || !elecsveto_scbId ||
            !floatVectors["Electron_deltaEtaSC"] || !intVectors["Electron_tightCharge"]) {
            if (!elecs_pt) std::cerr << "Error: elecs_pt is nullptr!" << std::endl;
            if (!elecs_eta) std::cerr << "Error: elecs_eta is nullptr!" << std::endl;
            if (!elecsveto_scbId) std::cerr << "Error: elecsveto_scbId is nullptr!" << std::endl;
            if (!floatVectors["Electron_deltaEtaSC"]) std::cerr << "Error: floatVectors[\"Electron_deltaEtaSC\"] is nullptr!" << std::endl;
            if (!intVectors["Electron_tightCharge"]) std::cerr << "Error: intVectors[\"Electron_tightCharge\"] is nullptr!" << std::endl;

            if (!elecs_pt || !elecs_eta || !elecsveto_scbId ||
                !floatVectors["Electron_deltaEtaSC"] || !intVectors["Electron_tightCharge"]) {
                std::cerr << "Error: Null pointer detected for variables!" << std::endl;
                return;
            }



            std::cerr << "Error: Null pointer detected for variables!" << std::endl;
            return;
        }

        if (elecs_pt->GetSize() <= iel || elecs_eta->GetSize() <= iel ||
            elecsveto_scbId->GetSize() <= iel ||
            (*floatVectors["Electron_deltaEtaSC"]).GetSize() <= iel ||
            (*intVectors["Electron_tightCharge"]).GetSize() <= iel) {
            std::cerr << "Error: Out-of-bounds access detected at index " << iel << std::endl;
            continue;
        }

        // 조건 디버깅
        if (!passKinematicCuts(elecs_pt->At(iel), elecs_eta->At(iel), lep_pt, lep_eta)) {
            std::cout << "Failed kinematic cuts for index " << iel << std::endl;
            continue;
        }

        if (!elecSCBId(elecsveto_scbId->At(iel), eleid_scbcut)) {
            std::cout << "Failed SCB ID for index " << iel << std::endl;
            continue;
         }

         double combinedEta = (*floatVectors["Electron_deltaEtaSC"])[iel] + elecs_eta->At(iel);
         if (fabs(combinedEta) > 1.4442 && fabs(combinedEta) < 1.566) {
             std::cout << "Failed eta range for index " << iel << std::endl;
             continue;
         }
         
         if (!elecCharge((*intVectors["Electron_tightCharge"])[iel])) {
             std::cout << "Failed charge check for index " << iel << std::endl;
             continue;
         }
         
         // 최종 선택된 elec 추가
         elecsveto.push_back(elecs.at(iel));
         std::cout << "Selected electron index: " << iel << std::endl;
    }


std::cout << "---sk 66666---" << std::endl;


}


void Analysis::LeptonOrder() {
    // Helper function for assigning leptons and debugging
    std::cout << "test !! LeptonOrder! " << std::endl;
    auto assignLeptons = [&](const std::vector<TLorentzVector>& leptons,
                             const std::string& chargeKey,
                             const std::vector<int>& indices,
                             int idx1, int idx2) {
        if (indices.size() <= std::max(idx1, idx2)) {
            std::cerr << "Error: indices size (" << indices.size() 
                      << ") is smaller than required index." << std::endl;
            return;
        }
        std::cout << "assignLeptons !!!!"<< std::endl;
        std::cout << "indices ? " << indices.size()<< std::endl;
        std::cout << "leptons ? " << leptons.size()<< std::endl;
        Lep1 = leptons.at(indices[idx1]);
        Lep2 = leptons.at(indices[idx2]);
        std::cout << "end Lep1 Lep2 assignment !  "<< std::endl;

        // Debugging information
        std::cerr << "Assigning Lep1 and Lep2:\n"
                  << "  Lep1: Pt = " << Lep1.Pt() << ", Eta = " << Lep1.Eta() << "\n"
                  << "  Lep2: Pt = " << Lep2.Pt() << ", Eta = " << Lep2.Eta() << std::endl;

        std::cout << "sk4-1"<< std::endl;
        if ((*intVectors[chargeKey])[indices[idx1]] < 0) {
            Lep = leptons.at(indices[idx1]);
            AnLep = leptons.at(indices[idx2]);
        } else {
            Lep = leptons.at(indices[idx2]);
            AnLep = leptons.at(indices[idx1]);
        }
        std::cout << "sk4-2"<< std::endl;
        
        // Additional debugging information
        std::cerr << "Assigning Lep and AnLep:\n"
                  << "  Lep: Charge < 0, Pt = " << Lep.Pt() << "\n"
                  << "  AnLep: Charge > 0, Pt = " << AnLep.Pt() << std::endl;
        std::cout << "sk4-3"<< std::endl;
    }; // end of assignLeptons //

    std::cout << "sk1 " << std::endl;

    // Handle dimuon decay mode
    if (TString(Decaymode).Contains("dimuon")) {
    std::cout << "sk2 " << std::endl;

        if (v_muon_idx.size() > 1) {
    std::cout << "sk3 " << std::endl;
            assignLeptons(muons, "Muon_charge", v_muon_idx, 0, 1);
        } /*else {
            std::cerr << "Lepton TLorentzVector Error: v_muon_idx is empty or too small for Decaymode = dimuon" 
                      << " (size: " << v_muon_idx.size() << ")" << std::endl;
            return;
        }*/
    std::cout << "sk4 " << std::endl;
    }
    else if (TString(Decaymode).Contains("dielec")) {
    // Handle dielectron decay mode
        if (v_electron_idx.size() > 1) {
            assignLeptons(elecs, "Electron_charge", v_electron_idx, 0, 1);
        } else {
            std::cerr << "Lepton TLorentzVector Error: v_electron_idx is empty or too small for Decaymode = dielec" 
                      << " (size: " << v_electron_idx.size() << ")" << std::endl;
            return;
        }
    }
    else if (TString(Decaymode).Contains("muel")) {
    // Handle muon-electron decay mode
        if (v_muon_idx.size() > 0 && v_electron_idx.size() > 0) {
            if (muons.at(v_muon_idx.at(0)).Pt() > elecs.at(v_electron_idx.at(0)).Pt()) {
                Lep1 = muons.at(v_muon_idx.at(0));
                Lep2 = elecs.at(v_electron_idx.at(0));
            } else {
                Lep1 = elecs.at(v_electron_idx.at(0));
                Lep2 = muons.at(v_muon_idx.at(0));
            }

            if ((*intVectors["Muon_charge"])[v_muon_idx[0]] < 0) {
                Lep = muons.at(v_muon_idx.at(0));
                AnLep = elecs.at(v_electron_idx.at(0));
            } else {
                Lep = elecs.at(v_electron_idx.at(0));
                AnLep = muons.at(v_muon_idx.at(0));
            }
        } else {
            std::cerr << "Lepton TLorentzVector Error: v_muon_idx size = " 
                      << v_muon_idx.size() << ", v_electron_idx size = " 
                      << v_electron_idx.size() << " for Decaymode = muel" << std::endl;
            return;
        }
    }
    // Handle invalid Decaymode
    else {
        std::cerr << "Lepton TLorentzVector Error: Decaymode = " << Decaymode << std::endl;
    }
}


void Analysis::MakeMuonCollection() {
    // make muon collection tlorentzvectors//
    std::cout << "Muon Correction !!" << std::endl;
    muons.clear();
    for (int imu = 0; imu < muons_pt->GetSize(); ++imu){
        //std::cout << "imu " << imu << " muons_pt->At(imu) : " << muons_pt->At(imu) << " muons_eta->At(imu) : " << muons_eta->At(imu) << " muons_phi->At(imu) : " << muons_phi->At(imu) << " muons_M->At(imu)  : "<< muons_M->At(imu)  << std::endl; 
        muons.push_back(createLorentzVector(muons_pt->At(imu), muons_eta->At(imu), muons_phi->At(imu), muons_M->At(imu) )); 
    }
    std::cout << "end of MakeMuonCollection !" << std::endl;
    std::cout << "size of muons.size : " << muons.size() << std::endl; 
    return;
}


void Analysis::MakeElecCollection() {
    // Electron collection logic
    elecs.clear();
    for (int iele = 0; iele < elecs_pt->GetSize(); ++iele){
        elecs.push_back(createLorentzVector(elecs_pt->At(iele), elecs_eta->At(iele), elecs_phi->At(iele), elecs_M->At(iele) )); 
    }
 
    std::cout << "end of MakeElecCollection !" << std::endl;
    std::cout << "size of elecs.size : " << elecs.size() << std::endl; 
    return;
}

void Analysis::MakeJetCollection() {
    // Jet collection logic
    jets.clear();
    for (int ijet = 0; ijet < jets_pt->GetSize(); ++ijet){
        jets.push_back(createLorentzVector(jets_pt->At(ijet), jets_eta->At(ijet), jets_phi->At(ijet), jets_M->At(ijet) )); 
    } 
    return;
}

void Analysis::JetSelector()
{
    MakeJetCollection();
    v_jet_idx.clear();

    Int_t injet = jets_pt->GetSize();

    TLorentzVector JetCalib;
    TLorentzVector JetRaw;

    double jetpt_ = -999;
    double jeteta_ = -999;
    double jetphi_ = -999;
    double jetenergy_ = -999;

    for (int i = 0; i < injet; i++)
    {
        //JetRaw = *(static_cast<TLorentzVector*>(jets.at(i)));
        JetRaw = jets.at(i);


        //////////////////////////////
        /// JES Systematic Method1 ///
        //////////////////////////////
        if (TString(JetEnSys).Contains("JetEnNorm"))
        {
            // 그대로 유지
        }
        else if (TString(JetEnSys).Contains("JetEnShiftedUp"))
        {
            jetpt_ = JetRaw.Pt();
            jeteta_ = JetRaw.Eta();
            jetphi_ = JetRaw.Phi();
            jetenergy_ = JetRaw.Energy();
            //JetRaw.SetPtEtaPhiE(jetpt_ * Jet_EnShiftedUp->at(i), jeteta_, jetphi_, jetenergy_ * Jet_EnShiftedUp->at(i));
        }
        else if (TString(JetEnSys).Contains("JetEnShiftedDown"))
        {
            jetpt_ = JetRaw.Pt();
            jeteta_ = JetRaw.Eta();
            jetphi_ = JetRaw.Phi();
            jetenergy_ = JetRaw.Energy();
            //JetRaw.SetPtEtaPhiE(jetpt_ * Jet_EnShiftedDown->at(i), jeteta_, jetphi_, jetenergy_ * Jet_EnShiftedDown->at(i));
        }
        else
        {
            std::cout << "Check out Jet Systematic Config ... Default is Norminal Jet ..." << std::endl;
        }

        //////////////////////////////
        /// JER Systematic Method1 ///
        //////////////////////////////
        JetCalib = (JERSmearing(&JetRaw, i, JetResSys));

        if (JetCleaning(&JetCalib))
        {
            //v_jetdpt_res.push_back(JetCalib.Pt() - JetRaw.Pt());
            //v_jetdpx_res.push_back(JetCalib.Px() - JetRaw.Px());
            //v_jetdpy_res.push_back(JetCalib.Py() - JetRaw.Py());

            ////////////////////////////////
            /// Apply Requirement of Jet ///
            ////////////////////////////////
            if (JetCalib.Pt() > jet_pt && fabs(JetCalib.Eta()) < jet_eta && JetId[i] >= jet_id)
            {
                v_jet_idx.push_back(i);
                //v_jet_TL.push_back(std::make_shared<TLorentzVector>(JetCalib));
                if (JetCalib.Pt() < 30.0)
                {
                    std::cout << "Something wrong " << std::endl;
                }
            }
        }
    }
}

bool Analysis::JetCleaning(TLorentzVector* jet_)
{
    for (const auto& ele : elecsveto)
    {
        if (ele.DeltaR(*jet_) < 0.4)
            return false;
    }
    for (const auto& mu : elecsveto)
    {
        if (mu.DeltaR(*jet_) < 0.4)
            return false;
    }
    return true;
}

TLorentzVector Analysis::JERSmearing(TLorentzVector* jet, int idx_, TString op_)
{
    double jerfrac_ = 1.0;

    if (op_.Contains("Norm")){
        //jerfrac_ = Jet_EnergyResolution_SF->at(idx_);
    }
    else if (op_ == "Up"){std::cout << "no op " << std::endl;
        //jerfrac_ = Jet_EnergyResolution_SFUp->at(idx_);
    }
    else if (op_ == "Down"){std::cout << "no op" << std::endl;
        //jerfrac_ = Jet_EnergyResolution_SFDown->at(idx_);
    }
    else{
        std::cout << "Check out your JERSmearing option !!" << "op_ : "<< op_<< std::endl;
    }
    if (dojer && !TString(FileName_).Contains("Data"))
    {
        return TLorentzVector(
            jet->Px() * jerfrac_,
            jet->Py() * jerfrac_,
            jet->Pz() * jerfrac_,
            jet->Energy() * jerfrac_
        );
    }

    return *jet;
}
bool Analysis::NumIsoLeptons()//YOU SHOULD REQUIRE THIS FUNCTION AFTER LEPTONSELETOR //
{        
    bool numLeptons = true;
    if ( TString(Decaymode).Contains( "dimuon" ) ||
         TString(Decaymode).Contains( "dielec" )    ){  if ( v_lepton_idx.size() <= 1 ) {numLeptons=false;} } //Requiring  Two isolated Lepton .
 
    else if ( TString(Decaymode).Contains( "muel" ) ){ if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1 ) {numLeptons=false; } else {}}
    else if ( TString(Decaymode).Contains( "muonJet" ) ){if( v_lepton_idx.size() != 1){numLeptons=false;} }
    else { std::cerr << "?? something wrong " << std::endl; }
    return numLeptons;
}

bool Analysis::LeptonsPtAddtional()//YOU SHOULD REQUIRE THIS FUNCTION AFTER NumIsoLeptons //                                        
{     
    bool lepptadd = false;                                                                                                               
    if ( TString(Decaymode).Contains( "dimuon" ) ||
         TString(Decaymode).Contains( "dielec" ) ||                                                                                      
         TString(Decaymode).Contains( "muel" )      ){  if ( Lep1.Pt() > 25 && Lep2.Pt() > 20 ) {lepptadd=true;} } //Pt of Leading Lepton should be over than 25 GeV and Second Leading Lepton Pt should be over thand 20 GeV.
       
    else if ( TString(Decaymode).Contains( "muonJet" ) ){if( v_lepton_idx.size() == 1){lepptadd=true;} }
    else { std::cerr << "?? something wrong " << std::endl; }                                                                                      
    return lepptadd;
}
void Analysis::Start()
{
   //fout = new TFile(Form("output/%s",outfile),"RECREATE");
   if (strcmp(outdir.c_str(), "None") != 0 ) {
      fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr/%s/%s", outdir.c_str(), outfile.c_str()), "RECREATE");

   }
   else {
      //fout = new TFile(Form("output/%s",outfile),"RECREATE");
      fout = new TFile(Form("output/%s", outfile.c_str()), "RECREATE");
   }
   std::cout << "fout - getname : " << fout->GetName() << std::endl;
   fout->cd("");

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();
}

void Analysis::DeclareHistos()
{
 
   for (int i =0 ; i < 10 ; i++)
   {
      h_Lep1pt[i]  = new TH1D(Form("h_Lep1pt_%d" ,i), Form("Leading Lepton pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep1pt[i]->Sumw2(); 
      h_Lep2pt[i]  = new TH1D(Form("h_Lep2pt_%d" ,i), Form("Second Leading Lepton pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep2pt[i]->Sumw2();
      h_Lep1eta[i] = new TH1D(Form("h_Lep1eta_%d",i), Form("Leading Lepton Eta    %s"    ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep1eta[i]->Sumw2();
      h_Lep2eta[i] = new TH1D(Form("h_Lep2eta_%d",i), Form("Second Leading Lepton Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep2eta[i]->Sumw2();
      h_Lep1phi[i] = new TH1D(Form("h_Lep1phi_%d",i), Form("Leading Lepton Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Lep1phi[i]->Sumw2();
      h_Lep2phi[i] = new TH1D(Form("h_Lep2phi_%d",i), Form("Second Leading Lepton Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Lep2phi[i]->Sumw2();
   
      h_Jet1pt[i]  = new TH1D(Form("h_Jet1pt_%d", i), Form("Leading Jet pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet1pt[i]->Sumw2();
      h_Jet2pt[i]  = new TH1D(Form("h_Jet2pt_%d", i), Form("Second Leading Jet pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet2pt[i]->Sumw2();
      h_Jet1eta[i] = new TH1D(Form("h_Jet1eta_%d",i), Form("Leading Jet Eta %s"       ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet1eta[i]->Sumw2();
      h_Jet2eta[i] = new TH1D(Form("h_Jet2eta_%d",i), Form("Second Leading Jet Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet2eta[i]->Sumw2();
      h_Jet1phi[i] = new TH1D(Form("h_Jet1phi_%d",i), Form("Leading Jet Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Jet1phi[i]->Sumw2();
      h_Jet2phi[i] = new TH1D(Form("h_Jet2phi_%d",i), Form("Second Leading Jet Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Jet2phi[i]->Sumw2();
    
     
      h_METpt[i]  = new TH1D(Form("h_METpt_%d",i), Form("MET pT %s" ,cutflowName[i].Data()), 200, 0.0, 200); h_METpt[i]->Sumw2();
      h_METphi[i] = new TH1D(Form("h_METphi_%d",i), Form("MET Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_METphi[i]->Sumw2();

      h_DiLepMass[i] = new TH1D(Form("h_DiLepMass_%d",i),Form("Di-Lepton Invariant Mass %s",cutflowName[i].Data()), 300, 0.0, 300); h_DiLepMass[i]->Sumw2();
      h_Num_PV[i]    = new TH1D(Form("h_Num_PV_%d",i),     Form("Num of Primary Vertex after %s",cutflowName[i].Data()), 100, 0.0, 100); h_Num_PV[i]->Sumw2();
      h_Num_Jets[i]  = new TH1D(Form("h_Num_Jets_%d",i), Form("Num. of Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_Jets[i]->Sumw2();
      h_Num_bJets[i] = new TH1D(Form("h_Num_bJets_%d",i),Form("Num. of b Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_bJets[i]->Sumw2();
    }
}


void Analysis::METDefiner()
{
//   Met.SetPtEtaPhiM(,0,,0);//MET_phi MET_pt // 
    //Met.SetPtEtaPhiM(met_pt, 0, met_phi, 0); // MET_phi MET_pt
    //Met.SetPtEtaPhiM(*met_pt, 0, *met_phi, 0); // MET_phi MET_pt
    //Met.SetPtEtaPhiM(static_cast<double>(*met_pt), 0, static_cast<double>(*met_phi), 0); // MET_phi MET_pt
    Met.SetPtEtaPhiM(static_cast<double>(**met_pt), 0, static_cast<double>(**met_phi), 0); // MET_phi MET_pt
}

bool Analysis::DiLeptonMassCut()
{        
   bool dimu_masscut = false;

   if ( ((Lep1)+(Lep2)).M() > 20 ){ dimu_masscut = true; }
         
   return dimu_masscut;
}        
      
