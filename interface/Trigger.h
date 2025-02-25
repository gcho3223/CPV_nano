#ifndef TRIGGER_H
#define TRIGGER_H

#include "Analysis.h"

class Trigger {
public:
    Trigger();
    ~Trigger();

    static bool Trig();
    static bool SelTrig(std::vector<std::string> v_sel);

    /// Maps for dynamic branch storage ///
    static std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> triggerList; //seltrigger

    static std::vector<std::string> DLtrigName; 
    static std::vector<std::string> SLtrigName; 
    static std::vector<std::string> trigName;
};

#endif // TRIGGER_H
