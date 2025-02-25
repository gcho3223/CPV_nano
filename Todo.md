# *To do*

## step 1
- Lepton requirement & Di-Lepton mass (> 20 GeV)cut

## step 2
- Z-mass veto < 76 GeV or > 106 GeV

## step 3
- nJet >= 2

## Sample
- UL2016PreVFP/TTbar_Signal_1.list

# Check tigger list for dilepton channel & SF, ID & WP, calcuate b-tagging efficiency
- ref
    - pdf: https://indico.cern.ch/event/995560/contributions/4189577/attachments/2174069/3671077/Dilepton_TriggerSF_TOPPAG.pdf
    - pdf: file:///C:/Users/suici/OneDrive/Desktop/TrigEffAN_temp.pdf
    - twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Trigger_Efficiency_Scale_Factors
    - b-tagging: https://btv-wiki.docs.cern.ch/PerformanceCalibration/fixedWPSFRecommendations/#b-tagging-efficiencies-in-simulation

## 1. Trigger list for dilepton channel & SF
### 2016
- emu
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL       | Data runs(B~G) & MC |
| HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL      | Data runs(B~G) & MC |
| HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ    | Data runs H         |
| HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ   | Data runs H         |
| HLT_Ele27_WPTight_Gsf                                | Data & MC           |
| HLT_IsoMu24                                          | Data & MC           |
| HLT_IsoTkMu24                                        | Data & MC           |
- ee
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_Ele23_Ele12_caloIdL_TrackIdL_IsoVL_DZ            | Data & MC           |
| HLT_DoubleEle33_CaloIdL_MW                           | Data & MC           |
| HLT_Ele27_WPTight_Gsf                                | Data & MC           |
- mumu
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL                     | Data runs(B~G) & MC | confirmed with Config file
| HLT_Mu17_TrkIsoVVL_TkMu8_TrackIdVVL                  | Data runs(B~G) & MC | confirmed with Config file
| HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ                  | Data runs H         | confirmed with Config file
| HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ                | Data runs H         | confirmed with Config file
| HLT_IsoMu24                                          | Data & MC           | confirmed with Config file
| HLT_IsoTkMu24                                        | Data & MC           | confirmed with Config file

### MET Filter
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_PFMET300                                         | Data & MC           |
| HLT_MET200                                           | Data & MC           |
| HLT_PFHT300_PFMET110                                 | Data & MC           |
| HLT_PFMET170_HBHECleaned                             | Data & MC           |
| HLT_PFMET(NoMu)120_PFMHT(NoMu)120_IDTight            | Data & MC           |

### 2017
- emu
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ    | Data & MC |
| HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL      | Data & MC |
| HLT_Ele35_WPTight_Gsf                                | Data & MC |
| HLT_IsoMu27                                          | Data & MC |
- ee
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_Ele23_Ele12_caloIdL_TrackIdL_IsoVL               | Data & MC           |
| HLT_DoubleEle33_CaloIdL_MW                           | Data & MC           |
| HLT_Ele35_WPTight_Gsf                                | Data & MC           |
- mumu
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass80           | Data & MC           |
| HLT_IsoMu17                                          | Data & MC           |


### 2018
- emu
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ    | Data & MC           |
| HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL      | Data & MC           |
| HLT_Ele32_WPTight_Gsf                                | Data & MC           |
| HLT_IsoMu24                                          | Data & MC           |
- ee
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_Ele23_Ele12_caloIdL_TrackIdL_IsoVL               | Data & MC           |
| HLT_DoubleEle25_CaloIdL_MW                           | Data & MC           |
| HLT_Ele32_WPTight_Gsf                                | Data & MC           |
- mumu
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8          | Data & MC           |
| HLT_IsoMu24                                          | Data & MC           |

### 2017 & 2018
|                   Name                               |        Dataset      |
|------------------------------------------------------|---------------------|
| HLT_PFMET200_HBHECleaned                             | Data & MC           |
| HLT_PFMET200_HBHE_BeamHaloCleaned                    | Data & MC           |
| HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned             | Data & MC           |
| HLT_PFMET(NoMu)120_PFMHT(NoMu)120_IDTight            | Data & MC           |
| HLT_PFMET(NoMu)120_PFMHT(NoMu)120_IDTight_PFHT60     | Data & MC           |
| HLT_PFHT500_PFMET100_PFMHT100_IDTight                | Data & MC           |
| HLT_PFHT700_PFMET85_PFMHT85_IDTight                  | Data & MC           |
| HLT_PFHT800_PFMET75_PFMHT75_IDTight                  | Data & MC           |

## from twiki(TOPTrigger)
### 2016
|                   Name                               |        Dataset      |    Channel     |
|------------------------------------------------------|---------------------|----------------|
| HLT_Iso(Tk)Mu24_v4                                   | Data & MC           |   SingleMuon   |
| HLT_Iso(Tk)Mu22_eta2p1_v4                            | Data & MC           |   SingleMuon   |
| HLT_Ele32_eta2p1_WPTight_Gsf_v8                      | Data & MC           | SingleElectron |
| HLT_Ele27_WPTight_Gsf_v7                             | Data & MC           | SingleElectron |
| HLT_Ele25_eta2p1_WPTight_Gsf_v7                      | Data & MC           | SingleElectron |
| HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9 or| Data & MC           | dilepton: emu  |
| HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4| Data & MC           | dilepton: emu  |
| HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9 or | Data & MC           | dilepton: e,mu |
| HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3 or| Data & MC           | dilepton: e,mu |
| HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4| Data & MC           | dilepton: e,mu |
| HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 or            | Data & MC           | dilepton: mu   |
| HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6             | Data & MC           | dilepton: mu   |
| HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9         | Data & MC           | dilepton: e    |
| HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v8             | Data & MC           | dilepton: e    |
| HLT_PFHT450_SixJet40_BTagCSV_p056_v5 or              | Data & MC           | fully hadronic |
| HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_v5                                                  |
### 2017
|                   Name                                              |        Dataset         |    Channel     |
|---------------------------------------------------------------------|------------------------|----------------|
| HLT_Ele35/38/40_WPTight_Gsf_v, HLT_Ele32_WPTight_Gsf_L1DoubleEG_v   | Data & MC              |       e        |
| HLT_IsoMu27_v OR HLT_IsoMu24_eta2p1_v                               | Data run(A~D) & MC     |      mu        |
| HLT_IsoMu27_v                                                       | Data run(E~F) & MC     |      mu        |
| HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v                                | Data & MC              |    e + HT      |
| HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v            | Data & MC              |    e + jet     |
| HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ                                 | Data run(A~B) & MC     |  double muon   |
| HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8/Mass3p8_v                 | Data run(from C) & MC  |  double muon   |
| HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL(_DZ)_v                       | Data & MC              | double electron|
| HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL(_DZ)_v or           | Data & MC              | muon + electron|
| HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v or             | Data & MC              | double electron|
| HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v                 | Data & MC              | double electron|
| HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v                           | Data run(A~B) & MC     | fully hadronic |
| HLT_PFHT430_SixJet40_BTagCSV_p080_v                                 -------------------------------------------
| HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v                        
| HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v                              | Data run(C~F) & MC     | fully hadronic |
| HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v                    | Data run(D) & MC       | fully hadronic |
### 2018
|                   Name                                        |   Dataset     |    Channel     |
|---------------------------------------------------------------|---------------|----------------|
| HLT_Ele32_WPTight_Gsf_v                                       | Data & MC     |       e        |
| HLT_IsoMu24_v                                                 | Data & MC     |      mu        |
| HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v                          | Data & MC     |    e + HT      |
| HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v      | Data & MC     |    e + jet     |
| HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v                 | Data & MC     |  double muon   |
| HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL(_DZ)_v                 | Data & MC     | double electron|
| HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL(_DZ)_v or     | Data & MC     | muon + electron|
| HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v or       | Data & MC     | muon + electron|
| HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v           | Data & MC     | muon + electron|
| HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5_v or                 | Data & MC     | fully hadronic |
| HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v              | Data & MC     | fully hadronic |

## 2. Object ID & WP
### electrons
    - pT > 15 GeV
    - |etaSuperCluster| < 2.4
    - Barrel-endcap gap removed 1.4442 < |etaSuperCluster| < 1.566
    - GSF electrons passing the POG recommended **VID-v2 tight** cut based working point
### muons
    - pT > 15 GeV (pT > 20 GeV)
    - |eta| < 2.4
    - Global muon passing the POG recommended **TightMuon** ID criteria
    - PF-based combined relative isolation (**lepPfIso < 0.15**)
### jets
    - Anti-kT R=0.4
    - pT > 30 GeV
    - |eta| < 2.4
    - Tight Lepton Veto PF Jet ID
    - Jet objects within dR(jet, lepton) < 0.4 overlap of selected leptons are not counted

## 3. Jet_puId apply to jet pT under 50 GeV
- ref1: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
- ref2: https://indico.cern.ch/event/1034671/contributions/4345150/attachments/2246420/3828130/PUID_UL16_UpdateWithSummer20Samples.pdf
- Jet_puId flag
    - puId == 0 means 000: fail all PU ID
    - puId == 1 means 001: pass loose ID, fail medium & tight
    - puId == 3 means 011: pass loose & medium ID, fail tight
    - puId == 7 means 111: pass all PU IDs (loose, medium, tight)
- (USE THIS) Loose  : 99(95)% efficiency for |η| < 2.5 (> 2.5) prompt jets
- Medium : 90% efficiency for prompt jets
- Tight  : 80% efficiency for prompt jets

## 4. b-tagging efficiency
