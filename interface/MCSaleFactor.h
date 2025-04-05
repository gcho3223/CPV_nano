#ifndef MCSaleFactor_H
#define MCSaleFactor_H

#include "Analysis.h"

class MCSaleFactor {
public:
    MCSaleFactor();
    ~MCSaleFactor();

    static void MCSF();
    static void MCSFApply();
    static void GenWeight();
    //static void PileupWeight();

    static double evt_weight_;
    static double evt_weight;
    static double evt_weight_gW_;
    // MC SF
    static double evt_weight_beforemcsf_;
    static double mc_sf_;
    // Gen weight
    static TTreeReaderValue<Float_t>* Gen_weight;
    static double evt_weight_beforegenweight_;
    static double genweight;
    // FactRen weight
    static double evt_weight_beforefactrenoweight_;
    // Decay table weight
    static double evt_weight_beforedectabweight_;
    // Fragment weight
    static double evt_weight_beforefragmentweight_;
    // PDF weight
    static double evt_weight_beforepdfweight_;
    // Pileup weight
    static double evt_weight_beforePileup_;
    static double puweight;
    // L1PreFire weight 
    static double evt_weight_beforeL1PreFire_;
    // Trigger weight
    static double evt_weight_beforeTrigger_;
    // Lepsf weight
    static double evt_weight_beforeLepsf_;
    // Btag weight
    static double evt_weight_beforeBtag_;
    // Lep_eff
    static double lep_eff;
};

#endif // MCSaleFactor_H
