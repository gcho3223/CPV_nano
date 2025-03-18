#include "DrawKinematics.h"

void DrawOverflowBin(TH1D *his, double min, double max)
{
    double contentUnder = 0;
    double contentOver  = 0;
    
    for (int i=0; i<his->GetNbinsX()-1; i++)
	{
        if (his->GetBinCenter(i+1) < min)
		{contentUnder += his->GetBinContent(i+1);}
        if (his->GetBinCenter(i+1) > max)
		{contentOver += his->GetBinContent(i+1);}
    }
    for (int i=0; i<his->GetNbinsX()-1; i++)
	{
        if (his->GetBinCenter(i+1) < min && his->GetBinCenter(i+2) > min) his->SetBinContent(i+2, contentUnder);
        if (his->GetBinCenter(i+1) < max && his->GetBinCenter(i+2) > max) his->SetBinContent(i+1, contentOver);
    }
    his->GetXaxis()->SetRangeUser(min, max);
}