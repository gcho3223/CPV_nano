#ifndef HISTOLIST_H
#define HISTOLIST_H
///////////////////////////////////////////////////
//              Histogram list                   //
///////////////////////////////////////////////////
#include <vector>
#include <TString.h>

vector<TString> hist1D =
{
    "h_Num_PV_beforePreSel", "h_Num_PV_afterMETFilter", "h_Num_PV_afterTrigger", "h_passedTriggers",
    // All selected muons
    "h_AllSelMu_multiplisity_wW", "h_AllSelMu_multiplisity_woW", 
    "h_AllSelMu_pt", "h_AllSelMu_eta", "h_AllSelMu_phi", "h_AllSelMu_iso",
    // All selected veto electrons
    "h_AllSelVetoEle_multiplisity_wW", "h_AllSelVetoEle_multiplisity_woW",
    "h_AllSelVetoEle_pt", "h_AllSelVetoEle_eta", "h_AllSelVetoEle_phi", "h_AllSelVetoEle_iso",

    "h_Num_PV_0",   "h_Num_PV_1",   "h_Num_PV_2",   "h_Num_PV_3",   "h_Num_PV_4",   "h_Num_PV_5",   //"h_Num_PV_6",   "h_Num_PV_7",   "h_Num_PV_8",   "h_Num_PV_9",
    
    "h_Lep1pt_0",   "h_Lep1pt_1",   "h_Lep1pt_2",   "h_Lep1pt_3",   "h_Lep1pt_4",   "h_Lep1pt_5",   //"h_Lep1pt_6",   "h_Lep1pt_7",   "h_Lep1pt_8",   "h_Lep1pt_9",
    "h_Lep1eta_0",  "h_Lep1eta_1",   "h_Lep1eta_2", "h_Lep1eta_3",  "h_Lep1eta_4",  "h_Lep1eta_5",  //"h_Lep1eta_6",  "h_Lep1eta_7",  "h_Lep1eta_8",  "h_Lep1eta_9",
    "h_Lep1phi_0",  "h_Lep1phi_1",  "h_Lep1phi_2",  "h_Lep1phi_3",  "h_Lep1phi_4",  "h_Lep1phi_5",  //"h_Lep1phi_6",  "h_Lep1phi_7",  "h_Lep1phi_8",  "h_Lep1phi_9",
    "h_Lep2pt_0",   "h_Lep2pt_1",   "h_Lep2pt_2",   "h_Lep2pt_3",   "h_Lep2pt_4",   "h_Lep2pt_5",   //"h_Lep2pt_6",   "h_Lep2pt_7",   "h_Lep2pt_8",   "h_Lep2pt_9",
    "h_Lep2eta_0",  "h_Lep2eta_1",  "h_Lep2eta_2",  "h_Lep2eta_3",  "h_Lep2eta_4",  "h_Lep2eta_5",  //"h_Lep2eta_6",  "h_Lep2eta_7",  "h_Lep2eta_8",  "h_Lep2eta_9",
    "h_Lep2phi_0",  "h_Lep2phi_1",  "h_Lep2phi_2",  "h_Lep2phi_3",  "h_Lep2phi_4",  "h_Lep2phi_5",  //"h_Lep2phi_6",  "h_Lep2phi_7",  "h_Lep2phi_8",  "h_Lep2phi_9",
    "h_DiLepMass_0","h_DiLepMass_1","h_DiLepMass_2","h_DiLepMass_3","h_DiLepMass_4","h_DiLepMass_5",//"h_DiLepMass_6","h_DiLepMass_7","h_DiLepMass_8","h_DiLepMass_9"//,

    "h_Jet1pt_0",   "h_Jet1pt_1",   "h_Jet1pt_2",   "h_Jet1pt_3",   "h_Jet1pt_4",   "h_Jet1pt_5",   //"h_Jet1pt_6",   "h_Jet1pt_7",   "h_Jet1pt_8",   "h_Jet1pt_9",
    "h_Jet1eta_0",  "h_Jet1eta_1",  "h_Jet1eta_2",  "h_Jet1eta_3",  "h_Jet1eta_4",  "h_Jet1eta_5",  //"h_Jet1eta_6",  "h_Jet1eta_7",  "h_Jet1eta_8",  "h_Jet1eta_9",
    "h_Jet1phi_0",  "h_Jet1phi_1",  "h_Jet1phi_2",  "h_Jet1phi_3",  "h_Jet1phi_4",  "h_Jet1phi_5",  //"h_Jet1phi_6",  "h_Jet1phi_7", "h_Jet1phi_8",  "h_Jet1phi_9",
    "h_Jet2pt_0",   "h_Jet2pt_1",   "h_Jet2pt_2",   "h_Jet2pt_3",   "h_Jet2pt_4",   "h_Jet2pt_5",   //"h_Jet2pt_6",   "h_Jet2pt_7",   "h_Jet2pt_8",   "h_Jet2pt_9",
    "h_Jet2eta_0",  "h_Jet2eta_1",  "h_Jet2eta_2",  "h_Jet2eta_3",  "h_Jet2eta_4",  "h_Jet2eta_5",  //"h_Jet2eta_6",  "h_Jet2eta_7",  "h_Jet2eta_8",  "h_Jet2eta_9",
    "h_Jet2phi_0",  "h_Jet2phi_1",  "h_Jet2phi_2",  "h_Jet2phi_3",  "h_Jet2phi_4",  "h_Jet2phi_5",  //"h_Jet2phi_6",  "h_Jet2phi_7",  "h_Jet2phi_8",  "h_Jet2phi_9",
    "h_Num_Jets_0", "h_Num_Jets_1", "h_Num_Jets_2", "h_Num_Jets_3", "h_Num_Jets_4", "h_Num_Jets_5", //"h_Num_Jets_6", "h_Num_Jets_7", "h_Num_Jets_8", "h_Num_Jets_9",

    "h_bJet1pt_0",  "h_bJet1pt_1",  "h_bJet1pt_2",  "h_bJet1pt_3",  "h_bJet1pt_4",  "h_bJet1pt_5",  //"h_bJet1pt_6",  "h_bJet1pt_7",  "h_bJet1pt_8",  "h_bJet1pt_9",
    "h_bJet1eta_0", "h_bJet1eta_1", "h_bJet1eta_2", "h_bJet1eta_3", "h_bJet1eta_4", "h_bJet1eta_5", //"h_bJet1eta_6", "h_bJet1eta_7", "h_bJet1eta_8", "h_bJet1eta_9",
    "h_bJet1phi_0", "h_bJet1phi_1", "h_bJet1phi_2", "h_bJet1phi_3", "h_bJet1phi_4", "h_bJet1phi_5", //"h_bJet1phi_6", "h_bJet1phi_7", "h_bJet1phi_8", "h_bJet1phi_9",
    "h_bJet2pt_0",  "h_bJet2pt_1",  "h_bJet2pt_2",  "h_bJet2pt_3",  "h_bJet2pt_4",  "h_bJet2pt_5",  //"h_bJet2pt_6",  "h_bJet2pt_7",  "h_bJet2pt_8",  "h_bJet2pt_9",
    "h_bJet2eta_0", "h_bJet2eta_1", "h_bJet2eta_2", "h_bJet2eta_3", "h_bJet2eta_4", "h_bJet2eta_5", //"h_bJet2eta_6", "h_bJet2eta_7", "h_bJet2eta_8", "h_bJet2eta_9",
    "h_bJet2phi_0", "h_bJet2phi_1", "h_bJet2phi_2", "h_bJet2phi_3", "h_bJet2phi_4", "h_bJet2phi_5", //"h_bJet2phi_6", "h_bJet2phi_7", "h_bJet2phi_8", "h_bJet2phi_9",
    "h_Num_bJets_0","h_Num_bJets_1","h_Num_bJets_2","h_Num_bJets_3","h_Num_bJets_4","h_Num_bJets_5",//"h_Num_bJets_6","h_Num_bJets_7","h_Num_bJets_8","h_Num_bJets_9",

    "h_METpt_0",    "h_METpt_1",    "h_METpt_2",    "h_METpt_3",    "h_METpt_4",    "h_METpt_5",    //"h_METpt_6",    "h_METpt_7",    "h_METpt_8",    "h_METpt_9",
    "h_METphi_0",   "h_METphi_1",   "h_METphi_2",   "h_METphi_3",   "h_METphi_4",   "h_METphi_5",   //"h_METphi_6",  "h_METphi_7",   "h_METphi_8",    "h_METphi_9",

    //"h_HT_0",       "h_HT_1",       "h_HT_2",       "h_HT_3",       "h_HT_4",       "h_HT_5",       "h_HT_6",       "h_HT_7",       "h_HT_8",       "h_HT_9",

    ///////////////////////
    // O3 after top reco //
    ///////////////////////
    //"h_Reco_CPO3",
};

vector<TString> hist2D =
{
	
};

#endif // HISTOLIST_H