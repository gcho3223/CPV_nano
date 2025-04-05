#include "../interface/Trigger.h"
//////////////////////////////////////////////////////////////////////////////////////
//                              Initialize variables                                //
//////////////////////////////////////////////////////////////////////////////////////
std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> Trigger::triggerList;
std::vector<std::string> Trigger::DLtrigName;
std::vector<std::string> Trigger::SLtrigName;
std::vector<std::string> Trigger::trigName;

std::unordered_map<std::string, int> Trigger::triggerPassCount;
std::vector<std::string> Trigger::passedTriggersThisEvent;

bool Trigger::Trig()
{
    //std::cout << "||------------------------Trigger------------------------||" << std::endl;
    /// Variable for Trigger Function
    bool singleTrig_ = false;
    bool doubleTrig_ = false;
    bool ispassselTrig_ = false;
    bool ispassvetoTrig_ = false;
    bool trigpass = false;
    std::vector<std::string> seltrigName;
    std::vector<std::string> vetotrigName;
    std::vector<std::string> filteredTrig;
 
    if(!TString(Analysis::FileName_).Contains("Data")) // MC //
    { 
        // Debugging line
        //std::cout << "FileName_ : " << Analysis::FileName_ << std::endl;
        //for(int i = 0; i < seltrigName.size(); i++) {std::cout << "seltrigName["<< i << "] : " << seltrigName[i] << std::endl;}
        ////trigpass = true; return trigpass; //original
        
        //ispassvetoTrig_ = false; //original
        //seltrigName = trigName; //original
        //trigpass = SelTrig(seltrigName); //original
        
        //---------------------------------------------------//
        //if(trigpass == true) {std::cout << "trigpass is true, PASSED!!" << std::endl;}
        //else {std::cout << "trigpass is false, FAILED!!" << std::endl;}
        
        ispassvetoTrig_ = false; 
        seltrigName = trigName;
        // print trigger list ///
        //std::cout << "seltrigName : ";
        //for(const auto &name:seltrigName) {std::cout << name << " ";}
        //std::cout << std::endl;
        for(const auto &trigger : seltrigName)
        {
            bool skipTrigger = false;
            if(TString(Analysis::Decaymode).Contains("dimuon"))
            {
                if(!TString(Analysis::RunPeriod).Contains("2016H")) // these trigger are only applied to Data 2016H dimuon channel!
                {
                    if(TString(trigger).Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ") ||
                       TString(trigger).Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"))
                    {
                        skipTrigger = true;
                        //std::cout << "skipTrigger: " << trigger << std::endl;
                    }
                }
                else {skipTrigger = false;}
            }
            if(!skipTrigger) {filteredTrig.push_back(trigger);}
        }
        trigpass = SelTrig(filteredTrig);
        //std::cout << "trigpass : " << trigpass << std::endl;
    }
    else // Data //
    {
        // Set veto trigger & selected trigger //
        // Channel Index //
        if(Analysis::RunPeriod.Contains("2018"))/// Only 2018, SingleEG and Double EG combined 
        {
            if(TString(Analysis::Decaymode).Contains("dimuon")) // Dimuon // 
            { 
                if(TString(Analysis::FileName_).Contains("Single"))
                {
                    seltrigName = SLtrigName;
                    vetotrigName = DLtrigName;
                    ispassselTrig_ = SelTrig(seltrigName);
                    ispassvetoTrig_ = SelTrig(vetotrigName);
                }
                else if(TString(Analysis::FileName_).Contains("Double"))
                {
                    seltrigName = DLtrigName;
                    vetotrigName = SLtrigName;
                    ispassselTrig_ = SelTrig(seltrigName);
                    ispassvetoTrig_ = SelTrig(vetotrigName);
                }
                else {std::cout << "Check out FileName_ in Trigger ()" << std::endl;}
            }
            else if(TString(Analysis::Decaymode).Contains("dielec"))
            { 
                //seltrigName = DLtrigName;
                //vetotrigName = SLtrigName;
                if(TString(Analysis::FileName_).Contains("EGamma"))
                {
                    singleTrig_ = false;
                    doubleTrig_ = false;
                    singleTrig_ = SelTrig(SLtrigName);
                    doubleTrig_ = SelTrig(DLtrigName);
                    if (singleTrig_&&(!doubleTrig_)) {trigpass = true;}
                    else if ((!singleTrig_)&&(doubleTrig_)) {trigpass = true;}
                    else {trigpass = false;}
                    ispassvetoTrig_ = false;
                }
            }
            else if(TString(Analysis::Decaymode).Contains("muel"))
            {
                if(TString(Analysis::FileName_).Contains("MuonEG"))
                {
                    seltrigName = DLtrigName;
                    vetotrigName = SLtrigName;
                    ispassselTrig_ = SelTrig(seltrigName);
                    ispassvetoTrig_ = SelTrig(vetotrigName);
                }
                else if(TString(Analysis::FileName_).Contains("EGamma"))
                {
                    seltrigName = SLtrigName;
                    vetotrigName = DLtrigName;
                    ispassselTrig_ = SelTrig(seltrigName);
                    ispassvetoTrig_ = SelTrig(vetotrigName);
                }
                else {std::cout << "Check out FileName_ in Trigger ()" << std::endl;}
            }
            else {std::cout << "Check out Decaymode in 2018-Trigger()" << std::endl;}
        } // end of 2018 //
        else /// 2016 (AVP, NonAPV) && 2017 ///
        {
            if(TString(Analysis::FileName_).Contains("Single"))
            {
                seltrigName = SLtrigName;
                vetotrigName = DLtrigName;
                ispassselTrig_ = SelTrig(seltrigName);
                ispassvetoTrig_ = SelTrig(vetotrigName);
                trigpass = ispassselTrig_;
                
                // Debugging line
                //std::cout << "Single Trigger!" << std::endl;
                //std::cout << "seltrigName: ";
                //for(const auto &name:seltrigName) {std::cout << name << " ";}
                //std::cout << std::endl;
                //std::cout << "vetotrigName: ";
                //for(const auto &name:vetotrigName) {std::cout << name << " ";}
                //std::cout << std::endl;
                //std::cout << "ispassselTrig_ : " << ispassselTrig_ << std::endl;
                //std::cout << "ispassvetoTrig_ : " << ispassvetoTrig_ << std::endl;
            }
            else if(TString(Analysis::FileName_).Contains("Double") || TString(Analysis::FileName_).Contains( "MuonEG"))
            {
                seltrigName = DLtrigName;
                vetotrigName = SLtrigName;
                ispassselTrig_ = SelTrig(trigName);
                //ispassvetoTrig_ = SelTrig(vetotrigName);
                trigpass = ispassselTrig_;
                
                // Debugging line
                //std::cout << "Double/MuonEG Trigger!" << std::endl;
                //std::cout << "seltrigName: ";
                //for(const auto &name:seltrigName) {std::cout << name << " ";}
                //std::cout << std::endl;
                //std::cout << "vetotrigName: ";
                //for(const auto &name:vetotrigName) {std::cout << name << " ";}
                //std::cout << std::endl;
                //std::cout << "ispassselTrig_ : " << ispassselTrig_ << std::endl;
                //std::cout << "ispassvetoTrig_ : " << ispassvetoTrig_ << std::endl;
            }
            else {std::cout << "Check out FileName_ in Trigger ()" << std::endl;}
        } // end of 2016 (AVP, NonAPV) && 2017 //
        if (ispassvetoTrig_ == true) {trigpass = false;}
    }
   //cout << "ispassselTrig_ : " << ispassselTrig_ << " ispassvetoTrig_ : " << ispassvetoTrig_ << "trigger : " << trigpass << endl; 
   return trigpass;
} // end of Trigger //

bool Trigger::SelTrig(std::vector<std::string> v_sel)
{
    std::string trgName = "";
    int ptrigindex = 0;
    bool passtrig_ = false;

    for(int j = 0; j < v_sel.size(); j++)
    {
        trgName = v_sel[j];
        /// trgName cout in triggerList ///
        auto it = triggerList.find(trgName);
        if(it != triggerList.end() && it->second)
        {
            if(**(it->second))
            {
                ptrigindex++;
                triggerPassCount[trgName]++;
                passedTriggersThisEvent.push_back(trgName);
                //std::cout << "trgName : " << trgName << " is passed" << std::endl;
            }
            else {//std::cout << "trgName : " << trgName << " is failed" << std::endl;
            }
        }
        else {std::cerr << "Error: Trigger " << trgName << " not found in triggerList." << std::endl;}
    }
    if (ptrigindex > 0) {passtrig_ = true;}
    return passtrig_;
} // end of SelTrig //