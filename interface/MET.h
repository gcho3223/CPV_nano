#ifndef MET_H
#define MET_H

#include "Analysis.h"

class MET {
public:
    MET();
    ~MET();

    static void METDefiner();
    static bool METFilterAPP();
    static bool METcut(TLorentzVector* met);

    static TTreeReaderValue<Float_t>* met_pt;
    static TTreeReaderValue<Float_t>* met_phi;
    static double met_cut;
    static bool metcut;
};

#endif // MET_H
