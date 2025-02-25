#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "../interface/Analysis.h"
#include "../interface/Config.h"
#include "../interface/MET.h"
#include "../interface/Jet.h"
#include "../interface/Lepton.h"
#include "../interface/Trigger.h"
#include "../interface/MCSaleFactor.h"

//////////////////////////////////////////////////////////////////////////////////////
//                              Initialize variables                                //
//////////////////////////////////////////////////////////////////////////////////////
TextReader *Analysis::SSBConfReader;
TString Analysis::FileName_;
TString Analysis::RunPeriod;
TString Analysis::Decaymode;
TString Analysis::XsecTable_;
std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> Analysis::boolSingles;
std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UInt_t>>> Analysis::uintSingles;
std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Float_t>>> Analysis::floatSingles;
std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>> Analysis::floatVectors;
std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Bool_t>>> Analysis::boolVectors;
std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Int_t>>> Analysis::intVectors;
std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> Analysis::noiseFilters;
double Analysis::Lumi;
TLorentzVector Analysis::Lep1 = TLorentzVector();
TLorentzVector Analysis::Lep2 = TLorentzVector();
TLorentzVector Analysis::Lep = TLorentzVector();
TLorentzVector Analysis::AnLep = TLorentzVector();
TLorentzVector Analysis::Met = TLorentzVector();
TLorentzVector Analysis::Jet1 = TLorentzVector();
TLorentzVector Analysis::Jet2 = TLorentzVector();
std::vector<TLorentzVector> Analysis::jets;
std::vector<int> Analysis::v_jet_idx;
std::vector<int> Analysis::v_jet_idx_puid;
std::vector<std::shared_ptr<TLorentzVector>> Analysis::v_jet_TL;
//std::vector<int> Analysis::v_jet_TL;
std::vector<int> Analysis::v_jet_Id;
TTreeReaderValue<Float_t>* Analysis::met_pt;
TTreeReaderValue<Float_t>* Analysis::met_phi;
int Analysis::eleid_scbcut;
float Analysis::bdisccut;
double Analysis::met_cut;
double Analysis::mc_sf_;
double Analysis::evt_weight_;
double Analysis::evt_weight_beforemcsf_;

// Constructor: initialize TTreeReader with TChain and branch list file
Analysis::Analysis(TChain *inputChain, std::string inputName, std::string seDirName, std::string outputName, const std::string &branchListFile, const std::string &configFile, int NumEvt= -1)
    : chain(inputChain), fReader(inputChain), NumEvt(NumEvt), outdir(seDirName), outfile(outputName)
{
    if (!chain) {
        throw std::runtime_error("Error: Invalid TChain pointer!");
    }
    FileName_ = SetInputFileName(inputName);
    std::cout << "=============================================================" << std::endl;  
    std::cout << "||                      --Analysis--                       ||" << std::endl;
    std::cout << "|| FileName_ : " << FileName_ << std::endl;
    std::cout << "|| configFile : " << configFile << std::endl;  
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
    cutflowName[8] = "Top-Reco";
    cutflowName[9] = "Top-Pt-Rewight";
 
    pi = TMath::Pi();
    Start();
}

Analysis::~Analysis()
{
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
//////////////////////////////////////////////////////////////////////////////////////
//                              Event loop function                                 //
//////////////////////////////////////////////////////////////////////////////////////
void Analysis::Loop()
{
    if(!chain)
    {
        std::cerr << "Error: TChain is null!" << std::endl;
        return;
    }
    Long64_t nEntries = chain->GetEntries();
    fReader.Restart(); // Reset the reader to the beginning
    evt_weight_ = 1;
    MCSaleFactor::MCSFApply();

    ////////////////
    // loop start //
    ////////////////
    for(Long64_t ientry = 0; ientry < ((NumEvt == -1 || NumEvt > nEntries) ? nEntries : NumEvt); ++ientry)
    {
        // start of event loop
        if(!fReader.Next()) { // Move to the next entry
            std::cerr << "Error: Failed to read entry " << ientry << std::endl;
            break;
        }
        //if (i > 100) break; //%lld supports Long64_t
        if(ientry % 10000 == 0) {printf("Event %lld\n", ientry);} //%lld supports Long64_t
        //std::cout << "----------------------- Event: " << ientry << " ----------------------- " << std::endl;
        //////////////////////////////////////////////////////////////////
        // Step 0-1 : Set Object Variables                              //
        Config::SetObjectVariable();
        Lepton::LeptonSelector(); 
        Lepton::LeptonOrder();
        Jet::JetSelector(); 
        MET::METDefiner();

        // check the vertex distribution before pre-selection //
        //std::cout << "Debug @step 0-1: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
        FillHisto(h_Num_PV_beforePreSel, **intSingles["PV_npvsGood"], evt_weight_);
        if(MET::METFilterAPP() == true) {FillHisto(h_Num_PV_afterMETFilter, **intSingles["PV_npvsGood"], evt_weight_);}
        if(Trigger::Trig() == true)     {FillHisto(h_Num_PV_afterTrigger, **intSingles["PV_npvsGood"], evt_weight_);}

        // Step 0-2 : Apply Noise (MET) Filter                          //
        if( MET::METFilterAPP() == false ) {continue;}
        // Step 0-3 : Apply Trigger Requirement                         //
        if( Trigger::Trig() == false ) {continue;}
        // Step 0-4 : Apply Good Primary Vertex Selection               //
        if(intSingles["PV_npvsGood"] && **intSingles["PV_npvsGood"] < 1) {continue;}
        //////////////////////////////////////////////////////////////////
        //std::cout << "Debug @step 0-4: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
        //std::cout << "evt_weight_ : " << evt_weight_ << std::endl;
        //std::cout << "Debug @Step 0: v_jet_idx size : " << Analysis::v_jet_idx.size() << std::endl;
        FillHisto(h_Num_PV[0],          **intSingles["PV_npvsGood"],  evt_weight_);
        FillHisto(h_DiLepMass[0],       ( Lep1 + Lep2 ).M(),          evt_weight_);
        FillHisto(h_Lep1pt[0],          Lep1.Pt(),                    evt_weight_);
        FillHisto(h_Lep1eta[0],         Lep1.Eta(),                   evt_weight_);
        FillHisto(h_Lep1phi[0],         Lep1.Phi(),                   evt_weight_);
        FillHisto(h_Lep2pt[0],          Lep2.Pt(),                    evt_weight_);
        FillHisto(h_Lep2eta[0],         Lep2.Eta(),                   evt_weight_);
        FillHisto(h_Lep2phi[0],         Lep2.Phi(),                   evt_weight_);
        FillHisto(h_METpt[0],           Met.Pt(),                     evt_weight_);
        FillHisto(h_METphi[0],          Met.Phi(),                    evt_weight_);
        FillHisto(h_Num_Jets[0],        v_jet_idx.size(),             evt_weight_);
        FillHisto(h_Num_Jets_puid[0],   v_jet_idx_puid.size(),        evt_weight_);
        //////////////////////////////////////////////////////////////////
        // Third lepton veto                                            //
        //////////////////////////////////////////////////////////////////
        //if(Lepton::ThirdLeptonVeto() == false ) {continue;}
        if(Lepton::LeptonsPtAddtional() == false ) {continue;}
        //////////////////////////////////////////////////////////////////
        // Step 1 : Lepton requirement & Di-Lepton mass (> 20 GeV)cut   //
        //////////////////////////////////////////////////////////////////
        //std::cout << "Debug @step 1: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
        if(Lepton::DiLeptonMassCut() == false ) {continue;}
 
        FillHisto(h_DiLepMass[1],       ( Lep1 + Lep2 ).M(),          evt_weight_); // You can check Z mass region before vetoing Z mass
        FillHisto(h_Num_PV[1],          **intSingles["PV_npvsGood"],  evt_weight_);
        FillHisto(h_Lep1pt[1],          Lep1.Pt(),                    evt_weight_);
        FillHisto(h_Lep1eta[1],         Lep1.Eta(),                   evt_weight_);
        FillHisto(h_Lep1phi[1],         Lep1.Phi(),                   evt_weight_);
        FillHisto(h_Lep2pt[1],          Lep2.Pt(),                    evt_weight_);
        FillHisto(h_Lep2eta[1],         Lep2.Eta(),                   evt_weight_);
        FillHisto(h_Lep2phi[1],         Lep2.Phi(),                   evt_weight_);
        FillHisto(h_METpt[1],           Met.Pt(),                     evt_weight_);
        FillHisto(h_METphi[1],          Met.Phi(),                    evt_weight_);
        FillHisto(h_Num_Jets[1],        v_jet_idx.size(),             evt_weight_);
        FillHisto(h_Num_Jets_puid[1],   v_jet_idx_puid.size(),        evt_weight_);

        //////////////////////////////////////////////////////////////////
        // Step 2 : Z-mass veto < 76 GeV or > 106 GeV                   //
        //////////////////////////////////////////////////////////////////
        //std::cout << "Debug @step 2: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
        if(Lepton::ZMassVeto() == false) {continue;}

        FillHisto(h_DiLepMass[2],       ( Lep1 + Lep2 ).M(),            evt_weight_);
        FillHisto(h_Num_PV[2],          **intSingles["PV_npvsGood"],    evt_weight_);
        FillHisto(h_Lep1pt[2],          Lep1.Pt(),                      evt_weight_);
        FillHisto(h_Lep1eta[2],         Lep1.Eta(),                     evt_weight_);
        FillHisto(h_Lep1phi[2],         Lep1.Phi(),                     evt_weight_);
        FillHisto(h_Lep2pt[2],          Lep2.Pt(),                      evt_weight_);
        FillHisto(h_Lep2eta[2],         Lep2.Eta(),                     evt_weight_);
        FillHisto(h_Lep2phi[2],         Lep2.Phi(),                     evt_weight_);
        FillHisto(h_METpt[2],           Met.Pt(),                       evt_weight_);
        FillHisto(h_METphi[2],          Met.Phi(),                      evt_weight_);
        FillHisto(h_Num_Jets[2],        v_jet_idx.size(),               evt_weight_);
        FillHisto(h_Num_Jets_puid[2],   v_jet_idx_puid.size(),          evt_weight_);

        //////////////////////////////////////////////////////////////////
        // Step 3 : nJet >= 2                                           //
        //////////////////////////////////////////////////////////////////
        //std::cout << "Debug @step 3: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
        Jet::JetDefiner();
        if(Jet::nJet(v_jet_idx) == true)
        {
            FillHisto(h_DiLepMass[3],       ( Lep1 + Lep2 ).M(),            evt_weight_);
            FillHisto(h_Num_PV[3],          **intSingles["PV_npvsGood"],    evt_weight_);
            FillHisto(h_Lep1pt[3],          Lep1.Pt(),                      evt_weight_);
            FillHisto(h_Lep1phi[3],         Lep1.Phi(),                     evt_weight_);
            FillHisto(h_Lep1eta[3],         Lep1.Eta(),                     evt_weight_);
            FillHisto(h_Lep2pt[3],          Lep2.Pt(),                      evt_weight_);
            FillHisto(h_Lep2eta[3],         Lep2.Eta(),                     evt_weight_);
            FillHisto(h_Lep2phi[3],         Lep2.Phi(),                     evt_weight_);
            FillHisto(h_METpt[3],           Met.Pt(),                       evt_weight_);
            FillHisto(h_METphi[3],          Met.Phi(),                      evt_weight_);
            FillHisto(h_Jet1pt[3],          Jet1.Pt(),                      evt_weight_);
            FillHisto(h_Jet1eta[3],         Jet1.Eta(),                     evt_weight_);
            FillHisto(h_Jet1phi[3],         Jet1.Phi(),                     evt_weight_);
            FillHisto(h_Jet2pt[3],          Jet2.Pt(),                      evt_weight_);
            FillHisto(h_Jet2eta[3],         Jet2.Eta(),                     evt_weight_);
            FillHisto(h_Jet2phi[3],         Jet2.Phi(),                     evt_weight_);
            FillHisto(h_Num_Jets[3],        v_jet_idx.size(),               evt_weight_);
            FillHisto(h_Num_Jets_puid[3],   v_jet_idx_puid.size(),          evt_weight_);

            //////////////////////////////////////////////////////////////////
            // Step 4 : MET cut                                             //
            // dimuon or dielec channel : MET > 40 GeV                      //
            // muel channel : pass this step                                //
            //////////////////////////////////////////////////////////////////
            if(MET::METcut(&Met) == true)
            {
                //std::cout << "Debug @step 4: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
                FillHisto(h_DiLepMass[4],       ( Lep1 + Lep2 ).M(),            evt_weight_);
                FillHisto(h_Num_PV[4],          **intSingles["PV_npvsGood"],    evt_weight_);
                FillHisto(h_Lep1pt[4],          Lep1.Pt(),                      evt_weight_);
                FillHisto(h_Lep1phi[4],         Lep1.Phi(),                     evt_weight_);
                FillHisto(h_Lep1eta[4],         Lep1.Eta(),                     evt_weight_);
                FillHisto(h_Lep2pt[4],          Lep2.Pt(),                      evt_weight_);
                FillHisto(h_Lep2eta[4],         Lep2.Eta(),                     evt_weight_);
                FillHisto(h_Lep2phi[4],         Lep2.Phi(),                     evt_weight_);
                FillHisto(h_METpt[4],           Met.Pt(),                       evt_weight_);
                FillHisto(h_METphi[4],          Met.Phi(),                      evt_weight_);
                FillHisto(h_Jet1pt[4],          Jet1.Pt(),                      evt_weight_);
                FillHisto(h_Jet1eta[4],         Jet1.Eta(),                     evt_weight_);
                FillHisto(h_Jet1phi[4],         Jet1.Phi(),                     evt_weight_);
                FillHisto(h_Jet2pt[4],          Jet2.Pt(),                      evt_weight_);
                FillHisto(h_Jet2eta[4],         Jet2.Eta(),                     evt_weight_);
                FillHisto(h_Jet2phi[4],         Jet2.Phi(),                     evt_weight_);
                FillHisto(h_Num_Jets[4],        v_jet_idx.size(),               evt_weight_);
                FillHisto(h_Num_Jets_puid[4],   v_jet_idx_puid.size(),          evt_weight_);
            }// end of step 4: MET cut //
        }// end of step 3:nJet >= 2
    }// end of event iteration //
    std::cout << "End Loop !!" << std::endl;
} // endl of loop //

TString Analysis::SetInputFileName(std::string inname)
{  
   TString inputName = inname;
    char unsco_ = '_';
    Size_t unscoIndex = inputName.Last(unsco_);
    inputName.Remove(unscoIndex, inputName.Length());
    if(inputName.Contains("/"))
    {
      unscoIndex = inputName.Last('/');
      inputName.Remove(0,unscoIndex+1);
    }
    //cout << "inputName : " << inputName << endl;
    return inputName;
} // end of SetInputFileName //

void Analysis::Start()
{
    //fout = new TFile(Form("output/%s",outfile),"RECREATE");
    if(strcmp(outdir.c_str(), "None") != 0 )
    {fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr/%s/%s", outdir.c_str(), outfile.c_str()), "RECREATE");}
    else
    {
       //fout = new TFile(Form("output/%s",outfile),"RECREATE");
       fout = new TFile(Form("output/%s", outfile.c_str()), "RECREATE");
    }
    std::cout << "||                       --Start--                         ||" << std::endl;
    std::cout << "|| fout - getname : " << fout->GetName() << std::endl;
    fout->cd("");

    TDirectory *dir = gDirectory;
    dir->cd();

    DeclareHistos();
} // end of start

void Analysis::DeclareHistos()
{
    ////////////////////////////////////////////
    /// kinematics of objects for each steps ///
    ////////////////////////////////////////////
    h_Num_PV_beforePreSel   = new TH1D(Form("h_Num_PV_beforePreSel"),   Form("Num. of Primary Vertex before Pre-Selection"),    100, 0.0, 100); h_Num_PV_beforePreSel->Sumw2();
    h_Num_PV_afterMETFilter = new TH1D(Form("h_Num_PV_afterMETFilter"), Form("Num. of Primary Vertex after MET Filter"),        100, 0.0, 100); h_Num_PV_afterMETFilter->Sumw2();
    h_Num_PV_afterTrigger   = new TH1D(Form("h_Num_PV_afterTrigger"),   Form("Num. of Primary Vertex after Trigger"),           100, 0.0, 100); h_Num_PV_afterTrigger->Sumw2();
    for(int i =0 ; i < 10 ; i++)
    {
        h_Num_PV[i]    = new TH1D(Form("h_Num_PV_%d",i),     Form("Num of Primary Vertex after %s",cutflowName[i].Data()), 100, 0.0, 100); h_Num_PV[i]->Sumw2();

        h_Lep1pt[i]         = new TH1D(Form("h_Lep1pt_%d" ,i),       Form("Leading Lepton pT %s"        ,   cutflowName[i].Data()), 250, 0.0, 250); h_Lep1pt[i]->Sumw2(); 
        h_Lep2pt[i]         = new TH1D(Form("h_Lep2pt_%d" ,i),       Form("Second Leading Lepton pT %s" ,   cutflowName[i].Data()), 250, 0.0, 250); h_Lep2pt[i]->Sumw2();
        h_Lep1eta[i]        = new TH1D(Form("h_Lep1eta_%d",i),       Form("Leading Lepton Eta    %s"    ,   cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep1eta[i]->Sumw2();
        h_Lep2eta[i]        = new TH1D(Form("h_Lep2eta_%d",i),       Form("Second Leading Lepton Eta %s",   cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep2eta[i]->Sumw2();
        h_Lep1phi[i]        = new TH1D(Form("h_Lep1phi_%d",i),       Form("Leading Lepton Phi %s"       ,   cutflowName[i].Data()), 24, -1*pi, pi); h_Lep1phi[i]->Sumw2();
        h_Lep2phi[i]        = new TH1D(Form("h_Lep2phi_%d",i),       Form("Second Leading Lepton Phi %s",   cutflowName[i].Data()), 24, -1*pi, pi); h_Lep2phi[i]->Sumw2();
        h_DiLepMass[i]      = new TH1D(Form("h_DiLepMass_%d",i),     Form("Di-Lepton Invariant Mass %s",    cutflowName[i].Data()), 300, 0.0, 300); h_DiLepMass[i]->Sumw2();
   
        h_Jet1pt[i]         = new TH1D(Form("h_Jet1pt_%d", i),       Form("Leading Jet pT %s"        ,      cutflowName[i].Data()), 250, 0.0, 250); h_Jet1pt[i]->Sumw2();
        h_Jet2pt[i]         = new TH1D(Form("h_Jet2pt_%d", i),       Form("Second Leading Jet pT %s" ,      cutflowName[i].Data()), 250, 0.0, 250); h_Jet2pt[i]->Sumw2();
        h_Jet1eta[i]        = new TH1D(Form("h_Jet1eta_%d",i),       Form("Leading Jet Eta %s"       ,      cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet1eta[i]->Sumw2();
        h_Jet2eta[i]        = new TH1D(Form("h_Jet2eta_%d",i),       Form("Second Leading Jet Eta %s",      cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet2eta[i]->Sumw2();
        h_Jet1phi[i]        = new TH1D(Form("h_Jet1phi_%d",i),       Form("Leading Jet Phi %s"       ,      cutflowName[i].Data()), 24, -1*pi, pi); h_Jet1phi[i]->Sumw2();
        h_Jet2phi[i]        = new TH1D(Form("h_Jet2phi_%d",i),       Form("Second Leading Jet Phi %s",      cutflowName[i].Data()), 24, -1*pi, pi); h_Jet2phi[i]->Sumw2();
        h_Num_Jets[i]       = new TH1D(Form("h_Num_Jets_%d",i),      Form("Num. of Jets after %s",          cutflowName[i].Data()), 20, 0.0, 20);   h_Num_Jets[i]->Sumw2();
        h_Num_bJets[i]      = new TH1D(Form("h_Num_bJets_%d",i),     Form("Num. of b Jets after %s",        cutflowName[i].Data()), 20, 0.0, 20);   h_Num_bJets[i]->Sumw2();
        h_Num_Jets_puid[i]  = new TH1D(Form("h_Num_Jets_puid_%d",i), Form("Num. of Jets after PU ID %s",    cutflowName[i].Data()), 20, 0.0, 20);   h_Num_Jets_puid[i]->Sumw2();
    
        h_METpt[i]  = new TH1D(Form("h_METpt_%d",i), Form("MET pT %s" ,cutflowName[i].Data()), 200, 0.0, 200); h_METpt[i]->Sumw2();
        h_METphi[i] = new TH1D(Form("h_METphi_%d",i), Form("MET Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_METphi[i]->Sumw2();
    }
} // end of DeclareHistos //