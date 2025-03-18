//////////////////////////////////////////////////////////////////////////////////
//	SinglePlotter.C Description													//
//	root -l -b -q SinglePlotter.C'("version","run_period","channel")'			//
//	version: version name														//
//	run_period: run period														//
//	channel: channel name														//
//////////////////////////////////////////////////////////////////////////////////
#include "DrawOverflowBin.cc"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLatex.h"

void Plotter_Single(string version, string period, string channel)
{
	void DrawHist1D(TH1D *hist1d, string sample, string version, string savepath);
	void DrawHist2D(TH2D *hist2d, string sample, string version, string savepath);

	vector<string> sample = {
        "TTbar_Signal_1"
    };
	vector<string> hist1D = {"h_passedTriggers"};
	vector<string> hist2D = {"h_triggerEvent"};

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
			TH1D *hist1d = (TH1D*)f->Get(Form("%s",hist1D[ih1].c_str()));
			if(hist1d == NULL){cout << "hist1d: " << hist1D[ih1].c_str() << " is not found" << endl;}
			DrawHist1D(hist1d,sample[isn],version,savepath); // draw kinematics for pt, eta, phi, mass, O3, Num of jets, Num of b-jets,...
		}
		for(int ih2=0; ih2<hist2D.size(); ih2++)
		{
			TH2D *hist2d = (TH2D*)f->Get(Form("%s",hist2D[ih2].c_str()));
			if(hist2d == NULL){cout << "hist2d: " << hist2D[ih2].c_str() << " is not found" << endl;}
			DrawHist2D(hist2d,sample[isn],version,savepath);
		}
		delete f;
	}
}

void DrawHist1D(TH1D *hist1d, string sample, string version, string savepath)
{
	TLegend *leg;

	TString hname = hist1d->GetName();

	bool isStat = false;
	//double xmin = 0;
	//double xmax = 100;
	//double ymin = 0;
	//double ymax = 100;
	double TitleOffset = 1.5;
	double xTitleOffset = 1.5;
	double yTitleOffset = 1.5;
	double xLabelSize = 0.023;
	double yLabelSize = 0.023;

	TCanvas *c1 = new TCanvas("c1",hname,800,800);
	c1->cd();
	c1->SetMargin(0.12,0.12,0.12,0.12); // left, right, top, bottom
	hist1d->Draw("hist");
	hist1d->SetTitle("Passed triggers");
	hist1d->SetTitleOffset(TitleOffset);
	hist1d->GetYaxis()->SetTitle("Entries");
	hist1d->GetYaxis()->SetTitleOffset(yTitleOffset);
	hist1d->GetXaxis()->SetLabelSize(xLabelSize);
	if(!isStat) {hist1d->SetStats(0);}
	
	c1->SaveAs(Form("%s/%s.pdf",savepath.c_str(),hist1d->GetName()));
}

void DrawHist2D(TH2D *hist2d, string sample, string version, string savepath)
{
	TLegend *leg;

	TString hname = hist2d->GetName();

	bool isRerange = true;
	bool isStat = false;
	double xmin = 0;
	double xmax = 100;
	double ymin = 0;
	double ymax = 100;
	double TitleOffset = 1.5;
	double xTitleOffset = 1.5;
	double yTitleOffset = 1.5;
	double xTitleSize = 0.025;
	double yTitleSize = 0.025;
	double xLabelSize = 0.023;
	double yLabelSize = 0.02;

	TCanvas *c2 = new TCanvas("c2",hname,800,800);
	c2->cd();
	c2->SetMargin(0.28,0.05,0.12,0.12); // left, right, top, bottom
	hist2d->Draw("col0");
	hist2d->SetTitle("passed triggers per event");
	hist2d->SetTitleOffset(TitleOffset);
	//hist2d->GetYaxis()->SetTitle("Entries");
	hist2d->GetYaxis()->SetTitleOffset(yTitleOffset);
	hist2d->GetYaxis()->SetLabelSize(yLabelSize);
	hist2d->GetYaxis()->LabelsOption(">");
	hist2d->GetXaxis()->SetTitle("Entries");
	hist2d->GetXaxis()->SetTitleOffset(xTitleOffset);
	hist2d->GetXaxis()->SetTitleSize(xTitleSize);
	hist2d->GetXaxis()->SetLabelSize(xLabelSize);

	if(!isStat) {hist2d->SetStats(0);}

	c2->SaveAs(Form("%s/%s.pdf",savepath.c_str(),hist2d->GetName()));

	if(isRerange)
	{
		hist2d->GetXaxis()->SetRangeUser(xmin,xmax);
		c2->SaveAs(Form("%s/%s_re.pdf",savepath.c_str(),hist2d->GetName()));
	}
}