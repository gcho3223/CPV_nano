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

//-------------------------------Initialize variables-------------------------------//
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
TLorentzVector Analysis::bJet1 = TLorentzVector();
TLorentzVector Analysis::bJet2 = TLorentzVector();

// Constructor: initialize TTreeReader with TChain and branch list file
Analysis::Analysis(TChain *inputChain, std::string inputName, std::string seDirName, std::string outputName, const std::string &branchListFile, const std::string &configFile, int NumEvt= -1)
    : chain(inputChain), fReader(inputChain), NumEvt(NumEvt), outdir(seDirName), outfile(outputName)
{
    if (!chain) {
        throw std::runtime_error("Error: Invalid TChain pointer!");
    }
    FileName_ = SetInputFileName(inputName);
    std::cout << "||------------------------Analysis------------------------||" << std::endl;
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
    // Apply MC SF //
    MCSaleFactor::evt_weight_ = 1;
    MCSaleFactor::evt_weight = 1;
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
        MCSaleFactor::evt_weight = MCSaleFactor::evt_weight_; //initialze evt_weight applied MC SF
        Config::SetObjectVariable();

        /// SF & weight
        MCSaleFactor::GenWeight();
        //MCSaleFactor::PileupWeight();

        /// object ///
        Lepton::LeptonSelector();
        Lepton::LeptonOrder();
        Jet::JetSelector(); 
        Jet::JetDefiner();
        MET::METDefiner();
        Jet::bTaggingJet();
        Jet::bTaggingJetDefiner();

        // check the vertex distribution before pre-selection //
        //std::cout << "Debug @step 0-1: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
        FillHisto(h_Num_PV_beforePreSel, **intSingles["PV_npvsGood"], MCSaleFactor::evt_weight);
        // check the selected muon //
        FillHisto(h_AllSelMu_multiplisity_wW,  Lepton::v_muon_idx.size(), MCSaleFactor::evt_weight);
        FillHisto(h_AllSelMu_multiplisity_woW, Lepton::v_muon_idx.size());
        for(int i=0; i<Lepton::v_muon_idx.size(); i++)
        {
            FillHisto(h_AllSelMu_pt,    Lepton::muons[Lepton::v_muon_idx[i]].Pt(),    MCSaleFactor::evt_weight);
            FillHisto(h_AllSelMu_eta,   Lepton::muons[Lepton::v_muon_idx[i]].Eta(),   MCSaleFactor::evt_weight);
            FillHisto(h_AllSelMu_phi,   Lepton::muons[Lepton::v_muon_idx[i]].Phi(),   MCSaleFactor::evt_weight);
            FillHisto(h_AllSelMu_iso,   Config::muons_iso->At(Lepton::v_muon_idx[i]), MCSaleFactor::evt_weight);
        }
        FillHisto(h_AllSelVetoEle_multiplisity_wW,  Lepton::elecsveto.size(),   MCSaleFactor::evt_weight);
        FillHisto(h_AllSelVetoEle_multiplisity_woW, Lepton::elecsveto.size());
        for(int i=0; i<Lepton::elecsveto.size(); i++)
        {
            FillHisto(h_AllSelVetoEle_pt,     Lepton::elecsveto[i].Pt(),        MCSaleFactor::evt_weight);
            FillHisto(h_AllSelVetoEle_eta,    Lepton::elecsveto[i].Eta(),       MCSaleFactor::evt_weight);
            FillHisto(h_AllSelVetoEle_phi,    Lepton::elecsveto[i].Phi(),       MCSaleFactor::evt_weight);
            FillHisto(h_AllSelVetoEle_iso,    Lepton::elecsveto_isoValue[i],    MCSaleFactor::evt_weight);
        }
        // Step 0-2 : Apply Noise (MET) Filter                          //
        if(MET::METFilterAPP() == false) {continue;}
        FillHisto(h_Num_PV_afterMETFilter, **intSingles["PV_npvsGood"], MCSaleFactor::evt_weight);
        // Step 0-3 : Apply Trigger Requirement                         //
        if(Trigger::Trig() == false) {continue;}
        FillHisto(h_Num_PV_afterTrigger, **intSingles["PV_npvsGood"], MCSaleFactor::evt_weight);
        //for(const auto &trigNamePassed : Trigger::passedTriggersThisEvent)
        //{
        //    auto it = std::find(Trigger::trigName.begin(), Trigger::trigName.end(), trigNamePassed);
        //    if(it != Trigger::trigName.end())
        //    {
        //        int trigIndex = std::distance(Trigger::trigName.begin(), it);
        //        h_triggerEvent->Fill(ientry, trigIndex);
        //    }
        //}
        //Trigger::passedTriggersThisEvent.clear();
        // Step 0-4 : Apply Good Primary Vertex Selection               //
        if(intSingles["PV_npvsGood"] && **intSingles["PV_npvsGood"] < 1) {continue;}
        //////////////////////////////////////////////////////////////////
        //std::cout << "Debug @step 0-4: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
        //std::cout << "Debug @Step 0: v_jet_idx size : " << Jet::v_jet_idx.size() << std::endl;
        FillHisto(h_Num_PV[0],          **intSingles["PV_npvsGood"],  MCSaleFactor::evt_weight);
        FillHisto(h_DiLepMass[0],       ( Lep1 + Lep2 ).M(),          MCSaleFactor::evt_weight);
        FillHisto(h_Lep1pt[0],          Lep1.Pt(),                    MCSaleFactor::evt_weight);
        FillHisto(h_Lep1eta[0],         Lep1.Eta(),                   MCSaleFactor::evt_weight);
        FillHisto(h_Lep1phi[0],         Lep1.Phi(),                   MCSaleFactor::evt_weight);
        FillHisto(h_Lep2pt[0],          Lep2.Pt(),                    MCSaleFactor::evt_weight);
        FillHisto(h_Lep2eta[0],         Lep2.Eta(),                   MCSaleFactor::evt_weight);
        FillHisto(h_Lep2phi[0],         Lep2.Phi(),                   MCSaleFactor::evt_weight);
        FillHisto(h_METpt[0],           Met.Pt(),                     MCSaleFactor::evt_weight);
        FillHisto(h_METphi[0],          Met.Phi(),                    MCSaleFactor::evt_weight);
        FillHisto(h_Num_Jets[0],        Jet::v_jet_idx.size(),        MCSaleFactor::evt_weight);
        
        ////////////////////////////////////////////////////////////////////////////////
        // Step 1                                                                     //
        // (1) Lepton requirement & Di-Lepton mass (> 20 GeV)cut                      //
        // (2) 3rd Lepton veto                                                        //
        // (3) LeptonsPtAddtional(leading lepton > 25GeV, sub-leading lepton > 20GeV) //
        ////////////////////////////////////////////////////////////////////////////////
        if(Lepton::NumIsoLeptons() == false) {continue;}
        if( ((TString(Decaymode).Contains("dimuon"))    && Lepton::ThirdLeptonVeto(Lepton::v_muon_idx) == true ) ||
            (TString(Decaymode).Contains("dielectron")) && Lepton::ThirdLeptonVeto(Lepton::v_electron_idx) == true ||
            (TString(Decaymode).Contains("muel"))       && Lepton::ThirdLeptonVeto(Lepton::v_muon_idx, Lepton::v_electron_idx) == true
        ) {continue;}
        if(Lepton::LeptonsPtAddtional() == false ) {continue;}
        if(Lepton::DiLeptonMassCut() == false ) {continue;}
        //std::cout << "Debug @step 1: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
 
        FillHisto(h_DiLepMass[1],       ( Lep1 + Lep2 ).M(),          MCSaleFactor::evt_weight); // You can check Z mass region before vetoing Z mass
        FillHisto(h_Num_PV[1],          **intSingles["PV_npvsGood"],  MCSaleFactor::evt_weight);
        FillHisto(h_Lep1pt[1],          Lep1.Pt(),                    MCSaleFactor::evt_weight);
        FillHisto(h_Lep1eta[1],         Lep1.Eta(),                   MCSaleFactor::evt_weight);
        FillHisto(h_Lep1phi[1],         Lep1.Phi(),                   MCSaleFactor::evt_weight);
        FillHisto(h_Lep2pt[1],          Lep2.Pt(),                    MCSaleFactor::evt_weight);
        FillHisto(h_Lep2eta[1],         Lep2.Eta(),                   MCSaleFactor::evt_weight);
        FillHisto(h_Lep2phi[1],         Lep2.Phi(),                   MCSaleFactor::evt_weight);
        FillHisto(h_METpt[1],           Met.Pt(),                     MCSaleFactor::evt_weight);
        FillHisto(h_METphi[1],          Met.Phi(),                    MCSaleFactor::evt_weight);
        FillHisto(h_Num_Jets[1],        Jet::v_jet_idx.size(),        MCSaleFactor::evt_weight);

        //////////////////////////////////////////////////////////////////
        // Step 2 : Z-mass veto < 76 GeV or > 106 GeV                   //
        //////////////////////////////////////////////////////////////////
        //std::cout << "Debug @step 2: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
        if(Lepton::ZMassVeto() == false) {continue;}

        FillHisto(h_DiLepMass[2],       ( Lep1 + Lep2 ).M(),            MCSaleFactor::evt_weight);
        FillHisto(h_Num_PV[2],          **intSingles["PV_npvsGood"],    MCSaleFactor::evt_weight);
        FillHisto(h_Lep1pt[2],          Lep1.Pt(),                      MCSaleFactor::evt_weight);
        FillHisto(h_Lep1eta[2],         Lep1.Eta(),                     MCSaleFactor::evt_weight);
        FillHisto(h_Lep1phi[2],         Lep1.Phi(),                     MCSaleFactor::evt_weight);
        FillHisto(h_Lep2pt[2],          Lep2.Pt(),                      MCSaleFactor::evt_weight);
        FillHisto(h_Lep2eta[2],         Lep2.Eta(),                     MCSaleFactor::evt_weight);
        FillHisto(h_Lep2phi[2],         Lep2.Phi(),                     MCSaleFactor::evt_weight);
        FillHisto(h_METpt[2],           Met.Pt(),                       MCSaleFactor::evt_weight);
        FillHisto(h_METphi[2],          Met.Phi(),                      MCSaleFactor::evt_weight);
        FillHisto(h_Num_Jets[2],        Jet::v_jet_idx.size(),          MCSaleFactor::evt_weight);

        //////////////////////////////////////////////////////////////////
        // Step 3 : nJet >= 2                                           //
        //////////////////////////////////////////////////////////////////
        //std::cout << "Debug @step 3: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
        if(Jet::nJet(Jet::v_jet_idx) == true)
        {
            FillHisto(h_DiLepMass[3],       ( Lep1 + Lep2 ).M(),            MCSaleFactor::evt_weight);
            FillHisto(h_Num_PV[3],          **intSingles["PV_npvsGood"],    MCSaleFactor::evt_weight);
            FillHisto(h_Lep1pt[3],          Lep1.Pt(),                      MCSaleFactor::evt_weight);
            FillHisto(h_Lep1phi[3],         Lep1.Phi(),                     MCSaleFactor::evt_weight);
            FillHisto(h_Lep1eta[3],         Lep1.Eta(),                     MCSaleFactor::evt_weight);
            FillHisto(h_Lep2pt[3],          Lep2.Pt(),                      MCSaleFactor::evt_weight);
            FillHisto(h_Lep2eta[3],         Lep2.Eta(),                     MCSaleFactor::evt_weight);
            FillHisto(h_Lep2phi[3],         Lep2.Phi(),                     MCSaleFactor::evt_weight);
            FillHisto(h_METpt[3],           Met.Pt(),                       MCSaleFactor::evt_weight);
            FillHisto(h_METphi[3],          Met.Phi(),                      MCSaleFactor::evt_weight);
            FillHisto(h_Jet1pt[3],          Jet1.Pt(),                      MCSaleFactor::evt_weight);
            FillHisto(h_Jet1eta[3],         Jet1.Eta(),                     MCSaleFactor::evt_weight);
            FillHisto(h_Jet1phi[3],         Jet1.Phi(),                     MCSaleFactor::evt_weight);
            FillHisto(h_Jet2pt[3],          Jet2.Pt(),                      MCSaleFactor::evt_weight);
            FillHisto(h_Jet2eta[3],         Jet2.Eta(),                     MCSaleFactor::evt_weight);
            FillHisto(h_Jet2phi[3],         Jet2.Phi(),                     MCSaleFactor::evt_weight);
            FillHisto(h_Num_Jets[3],        Jet::v_jet_idx.size(),          MCSaleFactor::evt_weight);

            //////////////////////////////////////////////////////////////////
            // Step 4 : MET cut                                             //
            // dimuon or dielec channel : MET > 40 GeV                      //
            // muel channel : pass this step                                //
            //////////////////////////////////////////////////////////////////
            if(MET::METcut(&Met) == true)
            {
                //std::cout << "Debug @step 4: **intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] << std::endl;
                FillHisto(h_DiLepMass[4],       ( Lep1 + Lep2 ).M(),            MCSaleFactor::evt_weight);
                FillHisto(h_Num_PV[4],          **intSingles["PV_npvsGood"],    MCSaleFactor::evt_weight);
                FillHisto(h_Lep1pt[4],          Lep1.Pt(),                      MCSaleFactor::evt_weight);
                FillHisto(h_Lep1phi[4],         Lep1.Phi(),                     MCSaleFactor::evt_weight);
                FillHisto(h_Lep1eta[4],         Lep1.Eta(),                     MCSaleFactor::evt_weight);
                FillHisto(h_Lep2pt[4],          Lep2.Pt(),                      MCSaleFactor::evt_weight);
                FillHisto(h_Lep2eta[4],         Lep2.Eta(),                     MCSaleFactor::evt_weight);
                FillHisto(h_Lep2phi[4],         Lep2.Phi(),                     MCSaleFactor::evt_weight);
                FillHisto(h_METpt[4],           Met.Pt(),                       MCSaleFactor::evt_weight);
                FillHisto(h_METphi[4],          Met.Phi(),                      MCSaleFactor::evt_weight);
                FillHisto(h_Jet1pt[4],          Jet1.Pt(),                      MCSaleFactor::evt_weight);
                FillHisto(h_Jet1eta[4],         Jet1.Eta(),                     MCSaleFactor::evt_weight);
                FillHisto(h_Jet1phi[4],         Jet1.Phi(),                     MCSaleFactor::evt_weight);
                FillHisto(h_Jet2pt[4],          Jet2.Pt(),                      MCSaleFactor::evt_weight);
                FillHisto(h_Jet2eta[4],         Jet2.Eta(),                     MCSaleFactor::evt_weight);
                FillHisto(h_Jet2phi[4],         Jet2.Phi(),                     MCSaleFactor::evt_weight);
                FillHisto(h_Num_Jets[4],        Jet::v_jet_idx.size(),          MCSaleFactor::evt_weight);

                //////////////////////////////////////////////////////////////////
                // Step 5 : b-tagging jet cut                                   //
                // 1 or more b-tagged jet                                       //
                //////////////////////////////////////////////////////////////////
                if(Jet::bTaggingJetCut(Jet::v_bjet_idx) == true)
                {
                    FillHisto(h_DiLepMass[5],       ( Lep1 + Lep2 ).M(),            MCSaleFactor::evt_weight);
                    FillHisto(h_Num_PV[5],          **intSingles["PV_npvsGood"],    MCSaleFactor::evt_weight);
                    FillHisto(h_Lep1pt[5],          Lep1.Pt(),                      MCSaleFactor::evt_weight);
                    FillHisto(h_Lep1phi[5],         Lep1.Phi(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_Lep1eta[5],         Lep1.Eta(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_Lep2pt[5],          Lep2.Pt(),                      MCSaleFactor::evt_weight);
                    FillHisto(h_Lep2eta[5],         Lep2.Eta(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_Lep2phi[5],         Lep2.Phi(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_METpt[5],           Met.Pt(),                       MCSaleFactor::evt_weight);
                    FillHisto(h_METphi[5],          Met.Phi(),                      MCSaleFactor::evt_weight);
                    FillHisto(h_Jet1pt[5],          Jet1.Pt(),                      MCSaleFactor::evt_weight);
                    FillHisto(h_Jet1eta[5],         Jet1.Eta(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_Jet1phi[5],         Jet1.Phi(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_Jet2pt[5],          Jet2.Pt(),                      MCSaleFactor::evt_weight);
                    FillHisto(h_Jet2eta[5],         Jet2.Eta(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_Jet2phi[5],         Jet2.Phi(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_bJet1pt[5],         bJet1.Pt(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_bJet1eta[5],        bJet1.Eta(),                    MCSaleFactor::evt_weight);
                    FillHisto(h_bJet1phi[5],        bJet1.Phi(),                    MCSaleFactor::evt_weight);
                    FillHisto(h_bJet2pt[5],         bJet2.Pt(),                     MCSaleFactor::evt_weight);
                    FillHisto(h_bJet2eta[5],        bJet2.Eta(),                    MCSaleFactor::evt_weight);
                    FillHisto(h_bJet2phi[5],        bJet2.Phi(),                    MCSaleFactor::evt_weight);
                    FillHisto(h_Num_Jets[5],        Jet::v_jet_idx.size(),          MCSaleFactor::evt_weight);
                    FillHisto(h_Num_bJets[5],       Jet::v_bjet_idx.size(),         MCSaleFactor::evt_weight);
                }
            }// end of step 4: MET cut //
        }// end of step 3:nJet >= 2
        int nStoredTriggers = Trigger::triggerPassCount.size();
        for(auto &pair : Trigger::triggerPassCount)
        {
            int xbin = h_passedTriggers->GetXaxis()->FindBin(pair.first.c_str());
            h_passedTriggers->SetBinContent(xbin, pair.second);
        }
        for(int i=0; i<Trigger::trigName.size(); i++)
        {
            h_passedTriggers->GetXaxis()->SetBinLabel(i+1, Trigger::trigName[i].c_str());
            //h_triggerEvent->GetYaxis()->SetBinLabel(i+1, Trigger::trigName[i].c_str());
        }
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

    // DeclareHistos() 호출 전에 maxEvt 설정
    maxEvt = (NumEvt == -1 || NumEvt > chain->GetEntries()) ? chain->GetEntries() : NumEvt;
    
    DeclareHistos();
} // end of start

void Analysis::DeclareHistos()
{
    ////////////////////////////////////////////
    /// kinematics of objects for each steps ///
    ////////////////////////////////////////////
    h_Num_PV_beforePreSel   = new TH1D(Form("h_Num_PV_beforePreSel"),   Form("Num. of Primary Vertex before Pre-Selection;# of PV;Entries"),    100, 0.0, 100); h_Num_PV_beforePreSel->Sumw2();
    h_Num_PV_afterMETFilter = new TH1D(Form("h_Num_PV_afterMETFilter"), Form("Num. of Primary Vertex after MET Filter;# of PV;Entries"),        100, 0.0, 100); h_Num_PV_afterMETFilter->Sumw2();
    h_Num_PV_afterTrigger   = new TH1D(Form("h_Num_PV_afterTrigger"),   Form("Num. of Primary Vertex after Trigger;# of PV;Entries"),           100, 0.0, 100); h_Num_PV_afterTrigger->Sumw2();
    h_passedTriggers        = new TH1D(Form("h_passedTriggers"),        Form("Num. of Triggers passed;;Entries"),                               6, 0, 6);       h_passedTriggers->Sumw2();
    
    /// Debug: pt, eta, phi, iso value for selected muon after object selection
    //after object selection
    h_AllSelMu_multiplisity_wW       = new TH1D(Form("h_AllSelMu_multiplisity_wW"),  Form("Selected muon multiplicity @ initial event w/ W;# of muons;Entries"), 5, 0.0, 5);       h_AllSelMu_multiplisity_wW->Sumw2();
    h_AllSelMu_multiplisity_woW      = new TH1D(Form("h_AllSelMu_multiplisity_woW"), Form("Selected muon multiplicity @ initial event wo W;# of muons;Entries"), 5, 0.0, 5);       h_AllSelMu_multiplisity_woW->Sumw2();
    h_AllSelMu_pt                    = new TH1D(Form("h_AllSelMu_pt"),               Form("All selected muon pT @ initial event;pT [GeV];Entries"),              1000, 0.0, 1000); h_AllSelMu_pt->Sumw2();
    h_AllSelMu_eta                   = new TH1D(Form("h_AllSelMu_eta"),              Form("All selected muon eta @ initial event;#eta;Entries"),                 50, -2.5, 2.5);   h_AllSelMu_eta->Sumw2();
    h_AllSelMu_phi                   = new TH1D(Form("h_AllSelMu_phi"),              Form("All selected muon phi @ initial event;#phi;Entries"),                 24, -1*pi, pi);   h_AllSelMu_phi->Sumw2();
    h_AllSelMu_iso                   = new TH1D(Form("h_AllSelMu_iso"),              Form("All selected muon iso @ initial event;Isolation;Entries"),            100, 0.0, 1);     h_AllSelMu_iso->Sumw2();
    // debug: all selected electrons
    h_AllSelVetoEle_multiplisity_wW  = new TH1D(Form("h_AllSelVetoEle_multiplisity_wW"),  Form("Selected electron multiplicity @ initial event w/ W;# of electrons;Entries"),   5, 0.0, 5);       h_AllSelVetoEle_multiplisity_wW->Sumw2();
    h_AllSelVetoEle_multiplisity_woW = new TH1D(Form("h_AllSelVetoEle_multiplisity_woW"), Form("Selected electron multiplicity @ initial event wo W;# of electrons;Entries"),   5, 0.0, 5);       h_AllSelVetoEle_multiplisity_woW->Sumw2();
    h_AllSelVetoEle_pt               = new TH1D(Form("h_AllSelVetoEle_pt"),               Form("All selected electron pT @ initial event;pT [GeV];Entries"),                    1000, 0.0, 1000); h_AllSelVetoEle_pt->Sumw2();
    h_AllSelVetoEle_eta              = new TH1D(Form("h_AllSelVetoEle_eta"),              Form("All selected electron eta @ initial event;#eta;Entries"),                       50, -2.5, 2.5);   h_AllSelVetoEle_eta->Sumw2();
    h_AllSelVetoEle_phi              = new TH1D(Form("h_AllSelVetoEle_phi"),              Form("All selected electron phi @ initial event;#phi;Entries"),                       24, -1*pi, pi);   h_AllSelVetoEle_phi->Sumw2();
    h_AllSelVetoEle_iso              = new TH1D(Form("h_AllSelVetoEle_iso"),              Form("All selected electron iso @ initial event;Isolation;Entries"),                  100, 0.0, 1);     h_AllSelVetoEle_iso->Sumw2();

    for(int i =0 ; i < 10 ; i++)
    {
        h_Num_PV[i]        = new TH1D(Form("h_Num_PV_%d",i),       Form("Num. of Primary Vertex after %s;# of PV;Entries",              cutflowName[i].Data()), 100, 0.0, 100);   h_Num_PV[i]->Sumw2();

        h_Lep1pt[i]        = new TH1D(Form("h_Lep1pt_%d" ,i),      Form("Leading Lepton pT %s;pT [GeV];Entries",                        cutflowName[i].Data()), 1000, 0.0, 1000); h_Lep1pt[i]->Sumw2(); 
        h_Lep1eta[i]       = new TH1D(Form("h_Lep1eta_%d",i),      Form("Leading Lepton Eta %s;#eta;Entries",                           cutflowName[i].Data()), 50, -2.5, 2.5);   h_Lep1eta[i]->Sumw2();
        h_Lep1phi[i]       = new TH1D(Form("h_Lep1phi_%d",i),      Form("Leading Lepton Phi %s;#phi;Entries",                           cutflowName[i].Data()), 24, -1*pi, pi);   h_Lep1phi[i]->Sumw2();
        h_Lep2pt[i]        = new TH1D(Form("h_Lep2pt_%d" ,i),      Form("Sub Leading Lepton pT %s;pT [GeV];Entries",                    cutflowName[i].Data()), 1000, 0.0, 1000); h_Lep2pt[i]->Sumw2();
        h_Lep2eta[i]       = new TH1D(Form("h_Lep2eta_%d",i),      Form("Sub Leading Lepton Eta %s;#eta;Entries",                       cutflowName[i].Data()), 50, -2.5, 2.5);   h_Lep2eta[i]->Sumw2();
        h_Lep2phi[i]       = new TH1D(Form("h_Lep2phi_%d",i),      Form("Sub Leading Lepton Phi %s;#phi;Entries",                       cutflowName[i].Data()), 24, -1*pi, pi);   h_Lep2phi[i]->Sumw2();
        h_DiLepMass[i]     = new TH1D(Form("h_DiLepMass_%d",i),    Form("Dilepton Invariant Mass %s;M [GeV];Entries",                   cutflowName[i].Data()), 300, 0.0, 300);   h_DiLepMass[i]->Sumw2();
   
        h_Jet1pt[i]        = new TH1D(Form("h_Jet1pt_%d", i),      Form("Leading Jet pT %s;pT [GeV];Entries",                           cutflowName[i].Data()), 1000, 0.0, 1000); h_Jet1pt[i]->Sumw2();
        h_Jet1eta[i]       = new TH1D(Form("h_Jet1eta_%d",i),      Form("Leading Jet Eta %s;#eta;Entries",                              cutflowName[i].Data()), 50, -2.5, 2.5);   h_Jet1eta[i]->Sumw2();
        h_Jet1phi[i]       = new TH1D(Form("h_Jet1phi_%d",i),      Form("Leading Jet Phi %s;#phi;Entries",                              cutflowName[i].Data()), 24, -1*pi, pi);   h_Jet1phi[i]->Sumw2();
        h_Jet2pt[i]        = new TH1D(Form("h_Jet2pt_%d", i),      Form("Sub Leading Jet pT %s;pT [GeV];Entries",                       cutflowName[i].Data()), 1000, 0.0, 1000); h_Jet2pt[i]->Sumw2();
        h_Jet2eta[i]       = new TH1D(Form("h_Jet2eta_%d",i),      Form("Sub Leading Jet Eta %s;#eta;Entries",                          cutflowName[i].Data()), 50, -2.5, 2.5);   h_Jet2eta[i]->Sumw2();
        h_Jet2phi[i]       = new TH1D(Form("h_Jet2phi_%d",i),      Form("Sub Leading Jet Phi %s;#phi;Entries",                          cutflowName[i].Data()), 24, -1*pi, pi);   h_Jet2phi[i]->Sumw2();
        h_Num_Jets[i]      = new TH1D(Form("h_Num_Jets_%d",i),     Form("Num. of Jets after %s;# of Jets;Entries",                      cutflowName[i].Data()), 20, 0.0, 20);     h_Num_Jets[i]->Sumw2();

        h_bJet1pt[i]       = new TH1D(Form("h_bJet1pt_%d", i),     Form("Leading b-tagging Jet pT %s;pT [GeV];Entries",                 cutflowName[i].Data()), 1000, 0.0, 1000); h_bJet1pt[i]->Sumw2();
        h_bJet1eta[i]      = new TH1D(Form("h_bJet1eta_%d",i),     Form("Leading b-tagging Jet Eta %s;#eta;Entries",                    cutflowName[i].Data()), 50, -2.5, 2.5);   h_bJet1eta[i]->Sumw2();
        h_bJet1phi[i]      = new TH1D(Form("h_bJet1phi_%d",i),     Form("Leading b-tagging Jet Phi %s;#phi;Entries",                    cutflowName[i].Data()), 24, -1*pi, pi);   h_bJet1phi[i]->Sumw2();
        h_bJet2pt[i]       = new TH1D(Form("h_bJet2pt_%d", i),     Form("Sub Leading b-tagging Jet pT %s;pT [GeV];Entries",             cutflowName[i].Data()), 1000, 0.0, 1000); h_bJet2pt[i]->Sumw2();
        h_bJet2eta[i]      = new TH1D(Form("h_bJet2eta_%d",i),     Form("Sub Leading b-tagging Jet Eta %s;#eta;Entries",                cutflowName[i].Data()), 50, -2.5, 2.5);   h_bJet2eta[i]->Sumw2();
        h_bJet2phi[i]      = new TH1D(Form("h_bJet2phi_%d",i),     Form("Sub Leading b-tagging Jet Phi %s;#phi;Entries",                cutflowName[i].Data()), 24, -1*pi, pi);   h_bJet2phi[i]->Sumw2();
        h_Num_bJets[i]     = new TH1D(Form("h_Num_bJets_%d",i),    Form("Num. of b-tagging Jets after %s;# of b-tagging Jets;Entries",  cutflowName[i].Data()), 20, 0.0, 20);     h_Num_bJets[i]->Sumw2();
    
        h_METpt[i]         = new TH1D(Form("h_METpt_%d",i),        Form("MET pT %s;pT [GeV];Entries",                                   cutflowName[i].Data()), 1000, 0.0, 1000); h_METpt[i]->Sumw2();
        h_METphi[i]        = new TH1D(Form("h_METphi_%d",i),       Form("MET Phi %s;#phi;Entries",                                      cutflowName[i].Data()), 24, -1*pi, pi);   h_METphi[i]->Sumw2();
    }
} // end of DeclareHistos //