//////////////////////////////////////////////////////////////////////////////////
//	SinglePlotter.C Description													//
//	root -l -b -q SinglePlotter.C'("version","run_period","channel")'			//
//	version: version name														//
//	run_period: run period														//
//	channel: channel name														//
//////////////////////////////////////////////////////////////////////////////////
#include "CPVari.cc"
#include "DrawKinematics.cc"
#include "DrawOverflowBin.cc"

void Plotter_Kinematics(string version, string period, string channel)
{
	vector<string> sample = {
        "TTbar_Signal_1"
    };

	string savepath = Form("./Job_Version/%s/%s/%s",version.c_str(),period.c_str(),channel.c_str()); // add sample name after channel

	for(int isn=0; isn<sample.size(); isn++)
	{
		TString fin;
		/// file read ///
		fin = Form("../output/%s/%s/%s/%s_test.root",version.c_str(),period.c_str(),channel.c_str(),sample[isn].c_str());
		TFile *f = new TFile(fin,"READ");
		if(f==NULL){cout << "File is not found........." << endl;}
		else {cout << fin.Data() << endl;}
		/// draw histogram!! ///
		for(int ih1=0; ih1<hist1D.size(); ih1++)
		{
			TH1D *hist1d = (TH1D*)f->Get(Form("%s",hist1D[ih1].Data()));
			if(hist1d == NULL){cout << "hist1d: " << hist1D[ih1].Data() << " is not found" << endl;}
			DrawKinematics(hist1d,sample[isn],version,savepath); // draw kinematics for pt, eta, phi, mass, O3, Num of jets, Num of b-jets,...
		}
		//for(int ih2=0; ih2<hist2D.size(); ih2++)
		//{
		//	TH2D *hist2d = (TH2D*)f->Get(Form("%s",hist2D[ih2].Data()));
		//	if(hist2d == NULL){cout << "hist2d: " << hist2D[ih2].Data() << " is not found" << endl;}
		//	Draw2DHist(hist2d,version,sample[isn],fout3d,iscpv);
		//}
		delete f;
	}
}