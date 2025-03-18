#include "DrawKinematics.h"

CPVari O3Asym(TH1D *hist)
{
    CPVari tmp;

    double posCount = 0.0;
    double posErr = 0.0;
    double fpos = 0.0;
    double pos_w2 = 0.0;
    double negCount = 0.0;
    double negErr = 0.0;
    double fneg = 0.0;
    double neg_w2 = 0.0;
    double ref_val = 20.496633;
    double norfac = 0.0000;
    double total = 0.0000;
    double Asym = 0.0000;
    double AsymErr = 0.0000;
    double xbin = 0;

    xbin = hist->GetNbinsX();
    //hist->Rebin(10);
    for (int i = 1; i <= hist->GetNbinsX(); i++)
	{
   		double binCenter = hist->GetBinCenter(i);
   		double binContent = hist->GetBinContent(i);
   		if (binCenter > 0)
        {
            posCount += binContent;
            posErr += hist->GetBinError(i);
            pos_w2 += hist->GetSumw2()->At(i); //Sumw2: sum of squares of weights
        }
   		else if (binCenter < 0)
        {
            negCount += binContent;
            negErr += hist->GetBinError(i);
            neg_w2 += hist->GetSumw2()->At(i);
        }
	}
    //norfac = ref_val/xsec;
    total = posCount + negCount;
    fpos = posCount / total;
    fneg = negCount / total;
    Asym = (posCount - negCount) / total;
    AsymErr = sqrt(
                    pow(2*negCount/pow(posCount+negCount,2),2)*pos_w2
                    + pow((-2)*posCount/pow(posCount+negCount,2),2)*neg_w2
                    );

    tmp.asym_ = Asym;
    tmp.asym_err_ = AsymErr;
    tmp.num_p_ = posCount;
    tmp.num_p_err_ = posErr;
    tmp.f_p_ = fpos;
    tmp.num_m_ = negCount;
    tmp.num_m_err_ = negErr;
    tmp.f_m_ = fneg;
    tmp.num_total_ = total;

    return tmp;
}