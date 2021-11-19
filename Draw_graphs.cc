#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"

#include <iostream>
#include <vector>
using namespace std;
using std::cout;
using std::endl;

void SetCanvas(TCanvas* c);
void SetHist(TH1* h,TString titx, TString tity);
void SavePNGandEPS(TCanvas** arrayCan, Int_t nobr);

Int_t Draw(){
	
	map<Int_t, Int_t> trigRanges = {{6, 7}, {12, 20}, {20, 30}};
	
	map<Int_t, vector<TH1D*>> histomap = {{6, {}}, {12, {}}, {20, {}}};
	
	//MAP OF HISTOGRAMS
	//KEY: TT = 
	//		6			12			20
	//VECTOR CONSTITUENTS:
	//0 = XSECT SPECTRUM WITH TT
	//1 = PP HISTOGRAM IN COUNTS
	//2 = OO HISTOGRAM IN COUNTS
	//3 = PP HISTOGRAM NORMALIZED PER TRIGGER TRACKS INTEGRAL
	//4 = OO HISTOGRAM NORMALIZED PER TRIGGER TRACKS INTEGRAL

	vector<TH1D*> DRhisto;
	
	//VECTOR OF DELTA_RECOIL HISTOGRAMS
	//0 = DELTA RECOIL FOR PP AT 6-7 && 12-20
	//1 = DELTA RECOIL FOR OO AT 6-7 && 12-20
	//2 = RATIO OF DELTA RECOIL FOR PP AND OO AT 6-7 && 12-20
	
	TString name = "";
	
	for (const auto& item : trigRanges){
		Int_t ttLow = item.first;
		name = Form("../TT-%d_%d/Res_tt_%d.root", ttLow, trigRanges.at(ttLow), ttLow);
		TFile *f11 = TFile::Open(name,"READ");
		if(!f11) return 20; 

		//READ Xsect HISTOGRAM 
		name = Form("hTT_MB_PartLevel_tt_%d", ttLow);
		TH1D*  RJSp_Xs = (TH1D*) f11->Get(name);
		histomap.at(ttLow).push_back(RJSp_Xs->Clone());
		if(!RJSp_Xs) return 21; 
		f11->Close();

		name = Form("../TT-%d_%d/Res_pois_tt_%d.root", ttLow);
		TFile *f12 = TFile::Open(name,"READ");
		if(!f12) return 20; 

		//READ counts pp HISTOGRAM 
		name = Form("fhRecoilJetPt_pp_MB_PartLevel_tt_%d", ttLow);   
		TH1D*  RJSp_Cpp = (TH1D*) f12->Get(name);
		histomap.at(ttLow).push_back(RJSp_Cpp->Clone());
		if(!RJSp_Cpp) return 21; 

		//READ counts OO HISTOGRAM 
		name = Form("fhRecoilJetPt_OO_MB_PartLevel_tt_%d", ttLow);   
		TH1D*  RJSp_COO = (TH1D*) f12->Get(name);
		histomap.at(ttLow).push_back(RJSp_COO->Clone());
		if(!RJSp_COO) return 21; 

		//READ normalized pp HISTOGRAM 
		name = Form("fhRecoilJetPt_pp_MB_PartLevel_tt_%d", ttLow);   
		TH1D*  RJSp_Npp = (TH1D*) f12->Get(name);
		histomap.at(ttLow).push_back(RJSp_Npp->Clone());
		if(!RJSp_Npp) return 21; 

		//READ normalized OO HISTOGRAM 
		name = Form("fhRecoilJetPt_OO_MB_PartLevel_tt_%d", ttLow);   
		TH1D*  RJSp_NOO = (TH1D*) f12->Get(name);
		histomap.at(ttLow).push_back(RJSp_NOO->Clone());
		if(!RJSp_NOO) return 21; 
		f12->Close();
	}

	name = Form("../DeltaRecoil_tt_6-7_tt_12-20.root");
	TFile *f13 = TFile::Open(name,"READ");
	if(!f13) return 20; 
	
	//READ Delta recoil 6-7 && 12-20 HISTOGRAM for pp 
	name = Form("Delta_RecoilJetPt_pp_MB_PartLevel_tt_6-7_tt_12-20");
	TH1D*  DRJSp_pp_612 = (TH1D*) f13->Get(name);
	DRhisto.push_back(DRJSp_pp_612->Clone());
	if(!DRJSp_pp_612) return 21; 
	
	//READ Delta recoil 6-7 && 12-20 HISTOGRAM for OO
	name = Form("Delta_RecoilJetPt_OO_MB_PartLevel_tt_6-7_tt_12-20");
	TH1D*  DRJSp_OO_612 = (TH1D*) f13->Get(name);
	DRhisto.push_back(DRJSp_OO_612->Clone());
	if(!DRJSp_OO_612) return 21; 
	
	f13->Close();
	
	
	

		//----------------------------------------------------------
		/*gStyle->SetOptTitle(0);
		gStyle->SetOptStat(0);
		TCanvas *c[10];
		TLegend *leg;
		Int_t io = 0;

		//DRAW CAVAS WHERE ONE COMAPARES BOTH HISTOGRAMS
		c[io] = new TCanvas(Form("c%d",io), "JetCrossSections",0,0,600,500 ); 
		SetCanvas((TCanvas*) c[io]);
		c[io]->cd()->SetLogy(); 

		leg = new TLegend(0.5,0.6,0.95,0.95," ","brNDC");  //this will be legend
		leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.04);
		leg->AddEntry((TObject*) 0, "pp #sqrt{#it{s}} = 13 TeV","");
		leg->AddEntry((TObject*) 0, "anti-#it{k}_{T} #it{R} = 0.4","");

		SetHist((TH1*) hard, "#it{p}_{T,ch jet} (GeV/#it{c})", "cross section (mb)");
		hard->SetLineColor(2);
		soft->SetLineColor(4);

		leg->AddEntry((TObject*) hard, "hard bins","l");
		leg->AddEntry((TObject*) soft, "minimum bias","l");

		hard->Draw(); //here you draw the first histo
		soft->Draw("same"); //here you 

		leg->Draw();
		io++;


		//Print out commands for figure printout in terminal
		SavePNGandEPS((TCanvas**) c, io);*/
	
	TCanvas* c = new TCanvas("canv1", "canvas", 0., 0., 600., 800.);
	c->Divide(5,5);
	c->cd(1);
	int i = 1;
	for (const auto& item: histomap) {
		for (const auto& graph: item.second){
			c->cd(i);
			graph->Draw();
			i++;
		}
	}
	
	
		return 11;
}

//_____________________________________________________________________


void SetHist(TH1* h,TString titx, TString tity){

   h->GetXaxis()->SetTitle(titx.Data());
   h->GetYaxis()->SetTitle(tity.Data());
   h->GetXaxis()->SetTitleSize(0.05);
   h->GetYaxis()->SetTitleSize(0.05);
   h->GetXaxis()->SetLabelSize(0.04);
   h->GetYaxis()->SetLabelSize(0.04);

   h->GetYaxis()->SetTitleOffset(1.3);
   h->GetXaxis()->SetTitleOffset(1.2);

}

//_____________________________________________________________________

void SetCanvas(TCanvas* c){
   c->SetLeftMargin(0.17);
   c->SetBottomMargin(0.14);
   c->SetRightMargin(0.05);
   c->SetTopMargin(0.05);
   c->SetTickx();
   c->SetTicky();
}


//_____________________________________________________________________


void SavePNGandEPS(TCanvas** arrayCan, Int_t nobr){

   cout<<endl<<endl<<endl;
   //PNG IMAGES
   gROOT->ProcessLine("TImage *img = TImage::Create();");
   for(int j=0; j<nobr;j++){
      TString nameBase =arrayCan[j]->GetTitle();
      if(nameBase.BeginsWith("x")) continue;
      TString namePNG = nameBase + ".png";
      TString OutImg = "img->FromPad(" + (TString) (arrayCan[j]->GetName()) + "); img->WriteImage(\"./Results/" + namePNG.Data() + "\");";
	  gROOT->ProcessLine(OutImg);
   }
/*
   //EPS IMAGES
   for(int j=0; j<nobr;j++){
      TString nameBase =arrayCan[j]->GetTitle();
      if(nameBase.BeginsWith("x")) continue;
      TString nameEPS = nameBase + ".eps";
      cout<<arrayCan[j]->GetName()<<"->SaveAs(\""<<nameEPS.Data()<<"\");"<<endl;
   }
*/
}


