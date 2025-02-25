ANCode/
├─ configs/
├─ interface/
│  ├─ Analysis.h
│  ├─ Config.h
│  ├─ Jet.h
│  ├─ Lepton.h
│  ├─ MCScaleFactor.h
│  ├─ MET.h
│  ├─ Trigger.h
├─ src/
│  ├─ Analysis.cc
│  ├─ Config.cc
│  ├─ Jet.cc
│  ├─ Lepton.cc
│  ├─ MCScaleFactor.cc
│  ├─ MET.cc
│  ├─ Trigger.cc
└─ main_ssb.cpp

<ConfigHandler>
InitBranches
SetVariables
SetObjectVariable

<MET>
METFilterAPP
METDefiner

<Trigger>
SelTrigger
Trigger

<MCSF>
MCSF
MCSFApply

<Lepton>
LeptonSelector
LeptonOrder
MakeMuonCollection
MakeElecCollection
NumIsoLeptons
LeptonsPtAddtional

<Jet>
MakeJetCollection
JetSelector
JetCleaning
JERSmearing

<Start>
Start

