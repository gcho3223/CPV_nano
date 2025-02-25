#include "../interface/MCSaleFactor.h"

void MCSaleFactor::MCSFApply()
{
    std::cout << "=============================================================" << std::endl;
    std::cout << "||                MC Scale Factor Apply                    ||" << std::endl;
    std::cout << "=============================================================" << std::endl;
    MCSF();
    Analysis::evt_weight_beforemcsf_ =1; // Initailize evt_weight_beforemcsf_ //
    Analysis::evt_weight_beforemcsf_ = Analysis::evt_weight_; // keep event weight //
     
    if ( !TString(Analysis::FileName_).Contains( "Data") )
    {
        Analysis::evt_weight_ = Analysis::evt_weight_ * Analysis::mc_sf_;
        std::cout <<  "|| evt_weight_ : "<< Analysis::evt_weight_ << ", mc_sf_ : "  << Analysis::mc_sf_  << std::endl;
    } // apply MC scale factor // 
    else {Analysis::evt_weight_ = 1;}
    std::cout << "|| MC SF(event weight) for Data: " << Analysis::evt_weight_<< std::endl;
} // end of MCSFApply //

void MCSaleFactor::MCSF()
{
    if (Analysis::FileName_.Contains("Data")||Analysis::FileName_.Contains("Single")||Analysis::FileName_.Contains("EG")){ Analysis::mc_sf_ = 1.; return; }
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
        std::cout << "||-------------------Load Xsection Table-------------------||" << std::endl;
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
        Analysis::mc_sf_ = (m_sam_xsec[Analysis::FileName_.Data()]*m_sam_br[Analysis::FileName_.Data()]*lumi)/m_sam_posi_nega[Analysis::FileName_.Data()];
        std::cout << "|| mc_sf_ " << Analysis::mc_sf_ << std::endl;
    }
    else
    {
        Analysis::mc_sf_ =1.;
        std::cout << "|| Key " << Analysis::FileName_.Data() << " not found in the std::map. mc sf is 1" << Analysis::mc_sf_ << std::endl;
    }
    return;
} // end of MCSF //