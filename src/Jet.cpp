#include "../interface/Jet.h"
#include "../interface/Config.h"
//-------------------------------Initialize variables-------------------------------//
std::vector<TLorentzVector> Jet::jets;
TString Jet::JetEnSys;
TString Jet::JetResSys;
TTreeReaderArray<Float_t>* Jet::jets_pt;
TTreeReaderArray<Float_t>* Jet::jets_eta;
TTreeReaderArray<Float_t>* Jet::jets_phi;
TTreeReaderArray<Float_t>* Jet::jets_M;
TTreeReaderArray<Float_t>* Jet::jets_btag;
//TTreeReaderArray<Bool_t>* Jet::jets_id;
TTreeReaderArray<Int_t>* Jet::jets_id;
TTreeReaderArray<Int_t>* Jet::jets_puid;
std::vector<int> Jet::v_jet_Id;
std::vector<int> Jet::v_jet_idx;
std::vector<int> Jet::v_bjet_idx;
std::vector<std::shared_ptr<TLorentzVector>> Jet::v_jet_TL;
std::vector<std::shared_ptr<TLorentzVector>> Jet::v_bjet_TL;
std::vector<double> Jet::v_jetdpt_res;
std::vector<double> Jet::v_jetdpx_res;
std::vector<double> Jet::v_jetdpy_res;
bool Jet::dojer;

Int_t Jet::injet;
Int_t Jet::cleanedjet;
Int_t Jet::passedjet;

bool Jet::njet;

double Jet::maxpt;
double Jet::subpt;

bool Jet::btaggingjet;
float Jet::btagcut;

int Jet::nbtagjet;

double Jet::maxbpt;
double Jet::subbpt;

void Jet::JetSelector()
{
    //std::cout << "||------------------------Jet Selector------------------------||" << std::endl;
    MakeJetCollection();
    v_jet_TL.clear();
    v_jet_idx.clear();

    v_jetdpt_res.clear();
    v_jetdpx_res.clear();
    v_jetdpy_res.clear();

    injet = jets_pt->GetSize();
    cleanedjet = 0;
    passedjet = 0;
    //std::cout << "|| Debug! JetSelector: injet : " << injet << std::endl;

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
        ////////////////////////////
        // JES Systematic Method1 //
        ////////////////////////////
        if(TString(JetEnSys).Contains("JetEnNorm"))
        {/* 그대로 유지*/}
        else if(TString(JetEnSys).Contains("JetEnShiftedUp"))
        {
            jetpt_ = JetRaw.Pt();
            jeteta_ = JetRaw.Eta();
            jetphi_ = JetRaw.Phi();
            jetenergy_ = JetRaw.Energy();
            //JetRaw.SetPtEtaPhiE(jetpt_ * Jet_EnShiftedUp->at(i), jeteta_, jetphi_, jetenergy_ * Jet_EnShiftedUp->at(i));
        }
        else if(TString(JetEnSys).Contains("JetEnShiftedDown"))
        {
            jetpt_ = JetRaw.Pt();
            jeteta_ = JetRaw.Eta();
            jetphi_ = JetRaw.Phi();
            jetenergy_ = JetRaw.Energy();
            //JetRaw.SetPtEtaPhiE(jetpt_ * Jet_EnShiftedDown->at(i), jeteta_, jetphi_, jetenergy_ * Jet_EnShiftedDown->at(i));
        }
        else
        {std::cout << "Check out Jet Systematic Config ... Default is Norminal Jet ..." << std::endl;}

        ////////////////////////////
        // JER Systematic Method1 //
        ////////////////////////////
        //std::cout << " || " << i << "th jet : " << JetRaw.Pt() << std::endl;
        JetCalib = (JERSmearing(&JetRaw, i, JetResSys));
        if(JetCleaning(&JetCalib))
        {
            cleanedjet++;
            v_jetdpt_res.push_back(JetCalib.Pt() - JetRaw.Pt());
            v_jetdpx_res.push_back(JetCalib.Px() - JetRaw.Px());
            v_jetdpy_res.push_back(JetCalib.Py() - JetRaw.Py());

            //////////////////////////////
            // Apply Requirement of Jet //
            //////////////////////////////
            //std::cout << "|| JetCalib.Eta() : " << fabs(JetCalib.Eta()) << ", jet_eta: " << Config::jet_eta << std::endl;
            //std::cout << "|| Jets_Id[" << i << "] : " << jets_id->At(i) << ", jet_id: " << Config::jet_id << std::endl;
            if(JetCalib.Pt() > Config::jet_pt && fabs(JetCalib.Eta()) < Config::jet_eta && jets_id->At(i) >= Config::jet_id)
            {
                //std::cout << " || " << i << "th jet (passed): " << JetCalib.Pt() << std::endl;
                if(JetCalib.Pt() < 50.0) // apply jet PU ID under 50 GeV
                {
                    if(jets_puid->At(i) >= Config::jet_puid)
                    {
                        v_jet_idx.push_back(i);
                        v_jet_TL.push_back(std::make_shared<TLorentzVector>(JetCalib));
                        //std::cout << "||---------------" << i << "th selected jet PU id ---------------||" << std::endl;
                        //std::cout << "|| JetCalib.Pt() : " << JetCalib.Pt() << std::endl;
                        //passedjet++; // for debug
                    }
                }
                else
                {
                    v_jet_idx.push_back(i);
                    v_jet_TL.push_back(std::make_shared<TLorentzVector>(JetCalib));
                    //std::cout << "||---------------" << i << "th selected jet Debug Info---------------||" << std::endl;
                    //std::cout << "|| JetCalib.Pt() : " << JetCalib.Pt() << ", jet_pt: " << Config::jet_pt << std::endl;
                    //passedjet++; // for debug
                }
                if(JetCalib.Pt() < 30.0) {std::cout << "pT of JetCalib is too small " << std::endl;}
            }
        }
    }
    //for(int i=0; i<v_jet_TL.size(); i++)
    //{std::cout << "Debug! JetSelector: v_jet_TL[" << i << "] : " << v_jet_TL[i]->Pt() << std::endl;}
    //std::cout << "|| Debug! JetSelector: v_jet_idx size : " << v_jet_idx.size() << ", v_jet_TL size : " << v_jet_TL.size() << std::endl;
    //std::cout << "|| Debug! JetSelector: injet : " << injet << ", cleanedjet : " << cleanedjet << ", passedjet : " << passedjet << std::endl;
} // end of JetSelector //

void Jet::MakeJetCollection()
{
    // Jet collection logic
    jets.clear();
    for (int ijet = 0; ijet < jets_pt->GetSize(); ++ijet)
    {jets.push_back(Config::createLorentzVector(jets_pt->At(ijet), jets_eta->At(ijet), jets_phi->At(ijet), jets_M->At(ijet) ));} 
    return;
} // end of MakeJetCollection //

bool Jet::JetCleaning(TLorentzVector* jet_)
{
    for (const auto& ele : Lepton::elecsveto) // remove jet if dR is samller than 0.4 because it is too close to electron
    {
        if (ele.DeltaR(*jet_) < 0.4) {return false;}
    }
    for (const auto& mu : Lepton::muonsveto)
    {
        if (mu.DeltaR(*jet_) < 0.4) {return false;}
    }
    return true;
} // end of JetCleaning //

TLorentzVector Jet::JERSmearing(TLorentzVector* jet, int idx_, TString op_)
{
    double jerfrac_ = 1.0;
    if (op_.Contains("Norm")) {/*jerfrac_ = Jet_EnergyResolution_SF->At(idx_)*/;}
    else if (op_ == "Up")
    {
        std::cout << "Up option" << std::endl;
        //jerfrac_ = Jet_EnergyResolution_SFUp->At(idx_);
    }
    else if (op_ == "Down")
    {
        std::cout << "Down option" << std::endl;
        //jerfrac_ = Jet_EnergyResolution_SFDown->At(idx_);
    }
    else
    {std::cout << "Check out your JERSmearing option !!" << "op_ : "<< op_<< std::endl;}
    if (dojer && !TString(Analysis::FileName_).Contains("Data"))
    {
        return TLorentzVector(
            jet->Px() * jerfrac_,
            jet->Py() * jerfrac_,
            jet->Pz() * jerfrac_,
            jet->Energy() * jerfrac_
        );
    }
    return *jet;
} // end of JERSmearing //
//////////////////////////
/// step 3 : nJet >= 2 ///
//////////////////////////
bool Jet::nJet(std::vector<int> v_jet)
{
    njet = false;
    if (v_jet.size() >= 2) {njet = true;}
    return njet;
}
/// find leading and subleading jet cut ///
void Jet::JetDefiner()
{
    //std::cout << " ************************************************************************ " << std::endl;
    //std::cout << "Debug! JetDefiner: v_jet_TL.size() : " << Analysis::v_jet_TL.size() << std::endl;
//
    maxpt = 0;
    subpt = 0;
    Analysis::Jet1 = TLorentzVector(); // Leading Jet
    Analysis::Jet2 = TLorentzVector(); // Subleading Jet

    if(v_jet_TL.size() >= 2)
    {
        for(int ijet=0; ijet<v_jet_TL.size(); ijet++)
        {
            //std::cout << "Debug! JetDefiner: ijet : " << ijet << ", pt: " << v_jet_TL[ijet]->Pt() << std::endl;
            if(v_jet_TL[ijet]->Pt() > maxpt)
            {/// get leading jet
                maxpt = v_jet_TL[ijet]->Pt();
                Analysis::Jet1 = *(v_jet_TL[ijet]);
            }
            else if(v_jet_TL[ijet]->Pt() <= Analysis::Jet1.Pt() && v_jet_TL[ijet]->Pt() > subpt)
            {/// get sub leading jet
                subpt = v_jet_TL[ijet]->Pt();
                Analysis::Jet2 = *(v_jet_TL[ijet]);
            }
        }
    }
    else {/*std::cout << "Debug! JetDefiner: the number of jets is less than 1 in v_jet_TL" << std::endl;*/ return;}
    //for(int ijet=0; ijet<Analysis::v_jet_TL.size(); ijet++)
    //{std::cout << "Debug! JetDefiner: v_jet_TL[" << ijet << "] : " << Analysis::v_jet_TL[ijet]->Pt() << std::endl;}
    //std::cout << "Debug! JetDefiner: Jet1 : " << Analysis::Jet1.Pt() << ", Jet2 : " << Analysis::Jet2.Pt() << std::endl;
    //std::cout << " ************************************************************************ " << std::endl;
}
//////////////////////////////////
/// step 5 : b-tagging jet cut ///
//////////////////////////////////
bool Jet::bTaggingJetCut(std::vector<int> v_bjet_idx)
{
    btaggingjet = false;
    if(v_bjet_idx.size() >= 1) {btaggingjet = true;}
    return btaggingjet;
}
/// find b-tagging jet ///
void Jet::bTaggingJet() 
{
    v_bjet_idx.clear();
    v_bjet_TL.clear();

    nbtagjet = 0;

    for(int ijet=0; ijet<v_jet_idx.size(); ijet++)
    {
        if(jets_btag->At(v_jet_idx[ijet]) > btagcut)
        {
            nbtagjet++;
            v_bjet_idx.push_back(v_jet_idx[ijet]);
            v_bjet_TL.push_back(v_jet_TL[ijet]);
        }
    }
}
/// find leading and subleading b-tagging jet ///
void Jet::bTaggingJetDefiner() 
{
    maxbpt = 0;
    subbpt = 0;
    Analysis::bJet1 = TLorentzVector();
    Analysis::bJet2 = TLorentzVector();

    if(v_jet_TL.size() >= 2)
    {
        /// case 1: v_bjet_TL.size() == 0 ///
        if(v_bjet_TL.size() == 0)
        {
            //std::cout << "Debug! case1 : v_bjet_TL.size() == 0" << std::endl;
            Analysis::bJet1 = Analysis::Jet1;
            Analysis::bJet2 = Analysis::Jet2;
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::Jet1 : " << Analysis::Jet1.Pt() << std::endl;
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::Jet2 : " << Analysis::Jet2.Pt() << std::endl;
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::bJet1 : " << Analysis::bJet1.Pt() << std::endl;
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::bJet2 : " << Analysis::bJet2.Pt() << std::endl;
        }
        /// case 2: v_bjet_TL.size() == 1 ///
        else if(v_bjet_TL.size() == 1)
        {
            //std::cout << "Debug! case2 : v_bjet_TL.size() == 1" << std::endl;
            Analysis::bJet1 = *(v_bjet_TL.at(0));
            if(Analysis::bJet1.Pt() == Analysis::Jet1.Pt()) {Analysis::bJet2 = Analysis::Jet2;}
            else {Analysis::bJet2 = Analysis::Jet1;}
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::Jet1 : " << Analysis::Jet1.Pt() << std::endl;
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::Jet2 : " << Analysis::Jet2.Pt() << std::endl;
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::bJet1 : " << Analysis::bJet1.Pt() << std::endl;
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::bJet2 : " << Analysis::bJet2.Pt() << std::endl;
        }
        /// case 3: v_bjet_TL.size() >= 2 ///
        else if(v_bjet_TL.size() >= 2)
        {
            //std::cout << "Debug! case3 : v_bjet_TL.size() > 1" << std::endl;
            for(int ibjet=0; ibjet<v_bjet_TL.size(); ibjet++)
            {
                if(v_bjet_TL[ibjet]->Pt() > maxbpt)
                {
                    maxbpt = v_bjet_TL[ibjet]->Pt();
                    Analysis::bJet1 = *(v_bjet_TL[ibjet]);
                }
                else if(v_bjet_TL[ibjet]->Pt() < Analysis::bJet1.Pt() && v_bjet_TL[ibjet]->Pt() > subbpt)
                {
                    subbpt = v_bjet_TL[ibjet]->Pt();
                    Analysis::bJet2 = *(v_bjet_TL[ibjet]);
                }
            }
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::bJet1 : " << Analysis::bJet1.Pt() << std::endl;
            //std::cout << "Debug! bTaggingJetDefiner: Analysis::bJet2 : " << Analysis::bJet2.Pt() << std::endl;
        }
    }
    else {return;}
    //std::cout << "Debug! bTaggingJetDefiner: Analysis::bJet1 : " << Analysis::bJet1.Pt() << ", Analysis::bJet2 : " << Analysis::bJet2.Pt() << std::endl;
}
