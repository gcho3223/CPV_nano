#include "../interface/MET.h"
//-------------------------------Initialize variables-------------------------------//
TTreeReaderValue<Float_t>* MET::met_pt;
TTreeReaderValue<Float_t>* MET::met_phi;
bool MET::metcut;
double MET::met_cut;

void MET::METDefiner()
{
//   Met.SetPtEtaPhiM(,0,,0);//MET_phi MET_pt // 
   //Met.SetPtEtaPhiM(Analysis::met_pt, 0, Analysis::met_phi, 0); // MET_phi MET_pt
   //Met.SetPtEtaPhiM(*met_pt, 0, *met_phi, 0); // MET_phi MET_pt
   //Met.SetPtEtaPhiM(static_cast<double>(*met_pt), 0, static_cast<double>(*met_phi), 0); // MET_phi MET_pt
   Analysis::Met.SetPtEtaPhiM(static_cast<double>(**(MET::met_pt)), 0, static_cast<double>(**(MET::met_phi)), 0);
} // end of METDefiner

bool MET::METFilterAPP()
{
   bool metfilt_ = true;
   for ( const auto &pair : Analysis::noiseFilters )
   { 
      TString METFiltName = pair.first;
//      if ( !(TString(FileName_).Contains("Data")) && 
//            TString(METFiltName).Contains("Flag_eeBadScFilter")  ) {continue;}
      if ( !(**pair.second) ){ metfilt_ = false; }
   }
   return metfilt_;
   //return true;
} // end of METFilterAPP

/// step 4 : MET > 40 GeV ///
bool MET::METcut(TLorentzVector* met)
{
   MET::metcut = false;
   if(TString(Analysis::Decaymode).Contains("dimuon") || TString(Analysis::Decaymode).Contains("dielec"))
   { if(met->Pt() > 40) {MET::metcut = true;} }
   else if(TString(Analysis::Decaymode).Contains("muel")) {MET::metcut = true;}
   else {std::cerr << "step 4(MET cut) Error Please check" << std::endl;}
   return MET::metcut;
}
