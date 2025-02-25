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

    static bool metcut;
};

#endif // MET_H
