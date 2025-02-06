# SSB NanoAOD Analysis Code

This repository contains code for analyzing CMS NanoAOD data using the SSB framework. The analysis processes event data and extracts relevant physics observables.

## Installation and Setup

### Cloning the Repository
To get started, clone the repository from GitHub and check out the required branch:
```sh
git clone --branch Run2_ULSummer20_v1 https://github.com/physicist87/SSBNanoAODANCode.git
cd SSBNanoAODANCode
```

### Compilation Instructions
To compile the analysis program, use the provided `Makefile`:
```sh
make -f Makefile
```
This will generate the executable **`ssb_analysis`**.

## Repository Structure and Code Summary

### **Main Components**
- **`main_ssb.cpp`**: The main execution file that initializes the analysis and handles input data.
- **`src/Analysis.cpp`**: Implements core analysis functions, including event selection and variable calculations.
- **`interface/Analysis.h`**: Defines the analysis class, including member functions and variables used for event processing.
- **`CommonTools.cpp/hpp`**: Provides utility functions used across multiple parts of the analysis.
- **`TextReader/TextReader.cpp/hpp`**: Handles reading text-based inputs for configurations or dataset lists.
- **`configs/`**: Contains various configuration files for different datasets and analysis setups.
- **`input/`**: Stores input dataset lists, specifying ROOT files to be analyzed.
- **`output/`**: Contains processed ROOT files after running the analysis.
- **`xsecAndsample/`**: Holds cross-section and sample information for different datasets.
- **`branchlist/`**: Contains lists of branches that are used in the analysis, ensuring compatibility with different versions of NanoAOD.

## Running the Analysis
### Input Files
Input file lists are stored in the `input/` directory. Example:
- `input/UL2016PreVFP/TTbar_Signal_1.list`

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

### Example Execution Commands
For a dimuon analysis in UL2016 PreVFP:
```sh
./ssb_analysis UL2016PreVFP Testv1 MuMu input/UL2016PreVFP/TTbar_Signal_1.list
```
For a dielectron analysis in UL2017:
```sh
./ssb_analysis UL2017 StudyX ElEl input/UL2016PreVFP/TTbar_Signal_1.list
```
For a muon-electron analysis in UL2018:
```sh
./ssb_analysis UL2018 AnalysisA MuEl input/UL2016PreVFP/TTbar_Signal_1.list
```

### Running a Check
To validate the setup and check if the program is correctly configured, use:
```sh
bash run_ssb_check.sh
```
Alternatively, execute it directly:
```sh
./run_ssb_check.sh
```

### Checking the Output
The processed ROOT files will be stored in the `output/` directory. Example:
```sh
ls output/Testv1/UL2016PreVFP/MuMu/
```
Expected output:
```
UL2016PreVFP_TTbar_Signal_1.root
TTbar_Signal_1.root
```

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

### **Example Entries**
```
HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, trigger, Bool_t, single
HLT_IsoMu24, trigger, Bool_t, single
nMuon, muon, UInt_t, single
Muon_pfRelIso03_all, muon, Float_t, vector
Muon_eta, muon, Float_t, vector
Muon_mass, muon, Float_t, vector
```
This structure ensures that only the necessary branches are accessed, optimizing performance and reducing memory usage.

## Notes
- Ensure that ROOT is properly installed and configured before running the analysis.
- Modify `Makefile` if necessary to match your system's compiler settings.
- Use the configuration files in `configs/` to customize your analysis workflow.

## Contact
For questions or contributions, please open an issue or contact the maintainers.

---
This README provides essential instructions for setting up and running the SSB NanoAOD analysis. If you need additional details, feel free to modify and expand it!

