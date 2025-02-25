# SSB NanoAOD Analysis Code

This branch contains a package for analyzing CMS NanoAOD data.
The package is modulized for their function.

## Package Structure and Code Summary
### Package Structure
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
### **Main Components**
- **`branchlist/`**: Contains lists of branches that are used in the analysis, ensuring compatibility with different versions of NanoAOD.
- **`configs/`**: Contains various configuration files for different datasets and analysis setups.
- **`input/`**: Stores input dataset lists, specifying ROOT files to be analyzed.
- **`interface/`**: Contains the header files for the analysis classes.
  - **`Analysis.h`**: Defines the analysis class, including member functions and variables used for event processing.
  - **`Config.h`**: Contains the configuration for the analysis.
  - **`Jet.h`**: Contains the jet class.
  - **`Lepton.h`**: Contains the lepton class.
  - **`MCScaleFactor.h`**: Contains the MC scale factor class.
  - **`MET.h`**: Contains the MET class.
  - **`Trigger.h`**: Contains the trigger class.
- **`src/`**: Analysis classes
  - **`Analysis.cpp`**: Analysis loop
  - **`Config.cpp`**: Configuration for reading branches, set variables & objects, create Lorentzvector 
  - **`Jet.cpp`**: contains jet selection, jet smearing and cleaning, jet difiner, discriminate the number of jets in the event
  - **`Lepton.cpp`**: contains lepton selection, odering leptons w.r.t pT, lepton pT additional, di-lepton mass cut, Z-mass veto, number of iso leptons, make lepton collection for muon and electron
  - **`MCScaleFactor.cpp`**: read MC scale factor table and apply to the event
  - **`MET.cpp`**: MET definer and filter application
  - **`Trigger.cpp`**: Trigger selection
- **`main_ssb.cpp`**: The main execution file that initializes the analysis and handles input data.
- **`CommonTools.cpp/hpp`**: Provides utility functions used across multiple parts of the analysis.
- **`TextReader/TextReader.cpp/hpp`**: Handles reading text-based inputs for configurations or dataset lists.
- **`output/`**: Contains processed ROOT files after running the analysis.
- **`xsecAndsample/`**: Holds cross-section and sample information for different datasets.

## Running the Analysis
### Execution Format
Based on `run_ssb_check.sh`, the program should be executed with the following arguments:
```sh
./ssb_analysis <runPeriod> <StudyName> <Channels> <input_list>
```
Where:
- **`runPeriod`**: Specifies the data-taking period (e.g., `UL2016PreVFP`, `UL2017`, `UL2018`).
- **`StudyName`**: Defines the study version (e.g., `Testv1`).
- **`Channels`**: The physics channel being analyzed:
  - `MuMu` → Dimuon channel
  - `ElEl` → Dielectron channel
  - `MuEl` → Muon-electron channel
- **`input_list`**: The dataset file (e.g., `input/UL2016PreVFP/TTbar_Signal_1.list`).

## Branchlist Directory and Its Role
The `branchlist/` directory contains text files that list the specific branches read from NanoAOD files. These files ensure compatibility with different NanoAOD versions by specifying which branches to extract. Each line in `branchlist/branch_list.txt` follows this format:
```
<Branch Name>, <Category>, <Data Type>, <Structure>
```
Where:
- **Branch Name**: The name of the branch in the ROOT file.
- **Category**: Indicates the type of physics object (e.g., `trigger`, `muon`, `electron`).
- **Data Type**: The type of data stored in the branch (e.g., `Bool_t`, `UInt_t`, `Float_t`).
- **Structure**: Specifies whether it is a `single` value or a `vector`.