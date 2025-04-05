#include "../interface/MCSaleFactor.h"
#include "../interface/Config.h"
//-------------------------------Initialize variables-------------------------------//
double MCSaleFactor::evt_weight_;
double MCSaleFactor::evt_weight;
double MCSaleFactor::evt_weight_gW_;
double MCSaleFactor::evt_weight_beforemcsf_;
double MCSaleFactor::mc_sf_;
double MCSaleFactor::evt_weight_beforegenweight_;
double MCSaleFactor::genweight;
TTreeReaderValue<Float_t>* MCSaleFactor::Gen_weight;

void MCSaleFactor::MCSFApply()
{
    //std::cout << "||------------------------MC Scale Factor Apply------------------------||" << std::endl;
    MCSF();
    evt_weight_beforemcsf_ =1; // Initailize evt_weight_beforemcsf_ //
    evt_weight_beforemcsf_ = evt_weight_; // keep event weight //
    
    if ( !TString(Analysis::FileName_).Contains( "Data") )
    {
        evt_weight_ = evt_weight_ * mc_sf_;
        //std::cout <<  "|| MCSFApply evt_weight_ : "<< evt_weight_ << ", mc_sf_ : "  << mc_sf_  << std::endl;
    } // apply MC scale factor // 
    else {evt_weight_ = 1;}
    //std::cout << "|| MC SF(event weight) for Data: " << evt_weight_<< std::endl;
} // end of MCSFApply //

void MCSaleFactor::MCSF()
{
    if (Analysis::FileName_.Contains("Data")||Analysis::FileName_.Contains("Single")||Analysis::FileName_.Contains("EG")){ mc_sf_ = 1.; return; }
    /// Open Xsec Tables ///
    FILE *xsecs_;
    char sampleName[1000];
    double xsec_ = -1.;
    double br_ = -1.; 
    int totalevt_  = -1.; 
    int positive_  = -1.; 
    int negative_  = -1.; 
    int posi_nega_ = -1.; 
    std::string xsec_dir= "./xsecAndsample/";
    std::string xsec_filePath = xsec_dir+ Analysis::XsecTable_.Data();
    //cout << "xsec_filePath : " << xsec_filePath << std::endl;
    /// SampleName | TotalEvt | Positive+Negative | Xsection | Branching Fraction |
    xsecs_ = fopen(xsec_filePath.c_str(),"r");
    std::map<std::string, int> m_sam_totalevt;
    std::map<std::string, double> m_sam_xsec;
    std::map<std::string, double> m_sam_br;
    std::map<std::string, int> m_sam_positive;
    std::map<std::string, int> m_sam_negative;
    std::map<std::string, int> m_sam_posi_nega;

    if (xsecs_!=NULL) 
    {
        std::cout << "||------------------------Load Xsection Table------------------------||" << std::endl;
        //cout << "Load Xsection Table!" << std::endl;
        while (fscanf(xsecs_, "%s %d %d %d %d %lf %lf\n", sampleName, &totalevt_, &positive_, &negative_, &posi_nega_, &xsec_, &br_ ) != EOF)
        {
            std::cout 
            << "|| sampleName : " << sampleName << " totalevt_ : " << totalevt_
            << " positive_ " << positive_ << " negative_ : " << negative_ 
            << " posi_nega_ " << posi_nega_ << " xsec_ : " << xsec_ 
            << " br_ " << br_ 
            << std::endl;
            m_sam_totalevt[sampleName] = totalevt_;
            m_sam_positive[sampleName] = positive_;
            m_sam_negative[sampleName] = negative_;
            m_sam_posi_nega[sampleName] = posi_nega_;
            m_sam_xsec[sampleName] = xsec_;
            m_sam_br[sampleName] = br_; //braching ratio
        }
        fclose(xsecs_);
    }
    else {std::cout << "No xsec_filePath !!!" << xsec_filePath << std::endl;return;}

    std::cout << "|| Lumi : " << Analysis::Lumi << std::endl;
    double lumi = Analysis::Lumi/1000000;
    
    auto it = m_sam_xsec.find(Analysis::FileName_.Data());
    if (it !=  m_sam_xsec.end())
    {
        std::cout << "|| Key " << Analysis::FileName_.Data() << " found in the std::map."<< std::endl;
        mc_sf_ = (m_sam_xsec[Analysis::FileName_.Data()]*m_sam_br[Analysis::FileName_.Data()]*lumi)/m_sam_posi_nega[Analysis::FileName_.Data()];
        std::cout << "|| mc_sf_ " << mc_sf_ << std::endl;
    }
    else
    {
        mc_sf_ =1.;
        std::cout << "|| Key " << Analysis::FileName_.Data() << " not found in the std::map. mc sf is 1" << mc_sf_ << std::endl;
    }
    return;
} // end of MCSF //

void MCSaleFactor::GenWeight()
{
    genweight = 1.0;
    evt_weight_beforemcsf_ = 1;
    evt_weight_beforemcsf_ = evt_weight_;
    //evt_weight_gW_ = evt_weight_;
    if(!TString(Analysis::FileName_).Contains("Data")) // MC sample //
    {
        if(*Gen_weight->Get() > 0.0) {genweight = 1;}
        else {genweight = -1;}
        evt_weight = evt_weight * genweight;
        //std::cout << "|| p2 Gen_weight : " << *Gen_weight->Get() << ", genweight: " << genweight << ", evt_weight_: " << evt_weight_ << ", evt_weight: " << evt_weight << std::endl;
    }
    else {evt_weight_ = 1;} // Data //
}

//void MCSaleFactor::PileupWeight()
//{
//    Config::puweight = 1.0;
//    Config::evt_weight_beforePileup_ = 1;
//    Config::evt_weight_beforePileup_ = Config::evt_weight_;
//    if(!TString(Config::FileName_).Contains("Data")) // MC sample //
//    {if(TString(Config::PileUpSys).Contains("central")) {Config::puweight = Pileup_nPU;}}
//    Config::evt_weight_ = Config::evt_weight_ * Config::puweight;
//    std::cout << "|| Pileup_nPU : " << Pileup_nPU << ", puweight : " << Config::puweight << ", evt_weight_ : " << Config::evt_weight_ << std::endl;
//}