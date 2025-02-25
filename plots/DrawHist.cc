#include "Histolist.h"
#include "DrawHist.h"

void DrawHist1D(TH1D *hist, string sample, string version, string savepath)
{
	void DrawOverflowBin(TH1D *his, double min, double max);
	CPVari O3Asym(TH1D *hist);
    CPVari tmp;

	TLegend *leg;

	TString hname = hist->GetName();

	string step;
	for(int i=0; i<10; i++)
	{
		if(hname.Contains(Form("_%d",i)))
		{step = Form("step%d",i); break;}
	}
	savepath = Form("%s/%s",savepath.c_str(),step.c_str());
	
	/// canvas ///
	TCanvas *c = new TCanvas("c",hname,800,800);
	c->cd();
	c->SetMargin(0.12,0.12,0.12,0.12);
	hist->Draw("hist");
	hist->SetTitle("");
	hist->GetYaxis()->SetTitle("Entries");
	hist->GetYaxis()->SetTitleOffset(1.5);
    /// Primary vertex ///
	if(hname.Contains("PV"))\
	{
		hist->GetXaxis()->SetRangeUser(0,100);
		//savepath = Form("%s/%s/PV",savepath.c_str());
	}
	/// pT ///
	if(hname.Contains("pt"))
	{
		//hist->Rebin(5);
		hist->GetXaxis()->SetTitle("GeV");
		DrawOverflowBin(hist, 0, 200);
		savepath = Form("%s/pT",savepath.c_str());
	}
	/// eta ///
	if(hname.Contains("eta"))
	{
		//hist->Rebin(2);
		hist->GetXaxis()->SetTitle("#eta");
		//DrawOverflowBin(hist, -2.5, 2.5);
		savepath = Form("%s/eta",savepath.c_str());
	}
	/// phi ///
	if(hname.Contains("phi"))
	{
		//hist->Rebin(4);
		hist->GetXaxis()->SetTitle("#phi");
		//DrawOverflowBin(hist, -3.14, 3.14);
		savepath = Form("%s/phi",savepath.c_str());
	}
	/// invariant mass ///
	if(hname.Contains("Mass"))
	{
		//hist->Rebin(5);
		hist->GetXaxis()->SetTitle("Invariant Mass (GeV)");
		//DrawOverflowBin(hist, 0, 300);
	}
    /// number of jets ///
	if(hname.Contains("Num_Jets") || hname.Contains("Num_bJets"))
	{
		hist->GetXaxis()->SetTitle("# of Jets");
        hist->GetXaxis()->SetRangeUser(0, 20);
		savepath = Form("%s/Num_Jets",savepath.c_str());
	}
	/// O3 ///
	if(hname.Contains("O3"))
	{
		c->SetLogy();
		//hist->Rebin(10);

		tmp = O3Asym(hist);

		hist->GetXaxis()->SetTitle("O_{3}/m_{t}^{4}");
		hist->SetStats(0);
		leg = new TLegend(.55,.7,.93,.87);
		leg->SetTextSize(0.023);
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
		if(sample == "Data_DoubleMuon" || sample == "Data_SingleMuon" || sample == "Data")
		{leg->AddEntry(hist, Form("Event: %0.f",tmp.num_total_),"");}
		else
		{leg->AddEntry(hist, Form("Event: %.3f",tmp.num_total_),"");}
		leg->AddEntry(hist, Form("Mean: %.4f",hist->GetMean()),"");
		leg->AddEntry(hist, Form("Std Dev: %.4f",hist->GetRMS()),"");
		if(sample == "DYJetsToLL_M_10To50" || sample == "DYJetsToLL_M_50" || sample == "ST_tW_antitop" || sample == "ST_tW_top" ||
		sample == "TTJets_Signal" || sample == "TTJets_others" ||
		sample == "TTbar_WJetToLNu" || sample == "TTbar_WQQ" || sample == "TTbar_ZQQ" || sample == "TTbar_ZToLLNuNu" ||
		sample == "WW" || sample == "WZ" || sample == "ZZ" || sample == "WJetsToLNu" ||
		sample == "DY" || sample == "ST" || sample == "TTV" || sample == "Diboson")
		{
			leg->AddEntry(hist, Form("N (O_{3} < 0): %.4f",tmp.num_m_),"");
			leg->AddEntry(hist, Form("N (O_{3} > 0): %.4f",tmp.num_p_),"");
		}
		else
		{
			leg->AddEntry(hist, Form("N (O_{3} < 0): %0.f",tmp.num_m_),"");
			leg->AddEntry(hist, Form("N (O_{3} > 0): %0.f",tmp.num_p_),"");
		}
		leg->AddEntry(hist, Form("A: %.4f",tmp.asym_),"");
		leg->Draw();

		//DrawOverflowBin(hist, -2, 2);
		savepath = Form("%s/O3",savepath.c_str());
	}

	/// save canvas ///
	gSystem->mkdir(Form("%s/",savepath.c_str()),kTRUE);
	c->SaveAs(Form("%s/%s.pdf",savepath.c_str(),hist->GetName()));

	delete c;
}
void DrawHist2D(TH2D *hist, string sample, string version, string savepath)
{
	TString hname = hist->GetName();

	TLegend *leg2d;

	TCanvas *c = new TCanvas("c",hname,800,800);

	c->cd();
	c->SetMargin(0.12,0.12,0.12,0.12);
	hist->Draw("colz");
	hist->SetTitle("");
	/// invariant mass ///
	if(hname.Contains("comp") || hname.Contains("Mass"))
	{
		hist->RebinX(10);
		hist->RebinY(10);
		hist->GetXaxis()->SetRangeUser(0,1000);
		hist->GetYaxis()->SetRangeUser(0,1000);
	}
	if(hname.Contains("dR") || hname.Contains("idx"))
	{
		hist->GetXaxis()->SetRangeUser(0,10);
		hist->GetYaxis()->SetRangeUser(0,10);
	}
	
	hist->GetYaxis()->SetTitleOffset(1.5);
	//hist->GetXaxis()->SetTitle("#DeltaR_{b jet}");
	//hist->GetYaxis()->SetTitle("#DeltaR_{#bar{b} jet}");

	// check the number of deltaR>3, deltaR<3
	double deltaR_under3 = hist->Integral(5,30,5,30);
	double deltaR_over3 = hist->Integral()- deltaR_under3;
	double frac_dR_under3 = deltaR_under3/hist->Integral();
	double frac_dR_over3 = deltaR_over3/hist->Integral();
	//fout3d << "Sample name: " << sample.c_str() << endl;
	//fout3d << "entries: " << hist->GetEntries() << ", Event: " << hist->Integral() << endl;
	//fout3d << "DeltaR < 3: " << deltaR_under3 << ", DeltaR > 3: " << deltaR_over3 << endl;
	//fout3d << "fraction: DeltaR < 3: " << deltaR_under3/hist->Integral() << " , DeltaR > 3: " << deltaR_over3/hist->Integral() << endl;
	//fout3d << endl;

	leg2d = new TLegend(.15,.7,.43,.8);
	leg2d->SetTextSize(0.023);
    leg2d->SetBorderSize(0);
    leg2d->SetFillColor(0);
    leg2d->SetFillStyle(0);
	leg2d->AddEntry(hist, Form("DeltaR < 3: %.3f",frac_dR_under3),"");
	leg2d->AddEntry(hist, Form("DeltaR > 3: %.3f",frac_dR_over3),"");
	//leg2d->Draw();

	// save canvas
	c->SaveAs(Form("%s/%s.pdf",savepath.c_str(),hist->GetName()));

	delete c;
}