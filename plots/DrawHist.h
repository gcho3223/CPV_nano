#ifndef DRAWHIST_H
#define DRAWHIST_H

#include <string>
#include <cmath>
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLatex.h"

using std::string;

void DrawHist1D(TH1D *hist, string sample, string version, string savepath);
void DrawHist2D(TH2D *hist, string sample, string version, string savepath);

void DrawOverflowBin(TH1D* his, double min, double max);

struct CPVari
{
    double asym_ = 0.;
    double asym_err_ = -999.;
    double num_p_ = 0.;
    double num_p_err_ = -999.;
    double f_p_ = 0.;
    double num_m_ = 0.;
    double num_m_err_ = -999.;
    double f_m_ = 0.;
    double num_total_ = 0.;
};

#endif // DRAWHIST_H