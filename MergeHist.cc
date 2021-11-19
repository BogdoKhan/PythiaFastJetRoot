#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TF1.h"

#include <vector>
#include <map>
#include <iostream>
using std::cout;
using std::cin;
using std::endl;

using namespace std; 

int GetSummarySpectrum(const Int_t& seed, const Int_t& TTLow);
void Scaling(Int_t trigRangeLow);
void DoPoisSmearing(Int_t trigRangeLow);
void MakeDeltaRecoilSp(const Int_t& trigRangeLow_1, const Int_t& trigRangeLow_2);

//int GetSummarySpectrum(const Int_t& seed, const Int_t& TTLow)
//void Scaling(Int_t trigRangeLow)
//void DoPoisSmearing(Int_t trigRangeLow)
//void MakeDeltaRecoilSp(const Int_t& trigRangeLow_1, const Int_t& trigRangeLow_2)

TString name = "";

enum {kMB, kFIT, kTRIG};  //kMB = pythia min bias,   kFIT =  event triggered with FIT coincidence
TString trg[]={"MB","FIT"}; 
enum {kOO, kpp, ktr};
TString RType[] = {"OO", "pp"};

TH1D* fhTTHPartLevel[kTRIG];
TH1D* fhRecoilJetPtTTH_PartLevel[kTRIG];
TH2D* fhRecoilJetPhiTTH_PartLevel[kTRIG]; 

void Sum(TFile *f1, const Int_t& seed, const Int_t& TTLow, const Int_t& ptval){
	Int_t cislo = seed;
	Int_t ttl = TTLow;
	Int_t pthmin = ptval;
	
	for(Int_t ig = kMB; ig< kTRIG; ig++){
		name = Form("hTT_%s_PartLevel_seed%d_tt_%d_pt%d", trg[ig].Data(), cislo, ttl, pthmin);
		TH1D* htt = (TH1D*) f1->Get(name);
		if(!htt){ 
			cout << "Wrong pT for this file!" << endl;
			return ;
				 }

		name = Form("fhRecoilJetPt_%s_PartLevel_seed%d_tt_%d_pt%d", trg[ig].Data(), cislo, ttl, pthmin);
		TH1D* hpt = (TH1D*) f1->Get(name);
		if(!hpt){ 
			cout << "Wrong pT for this file!" << endl;
			return ;
				 }

		name = Form("fhRecoilJetPhi_%s_PartLevel_seed%d_tt_%d_pt%d", trg[ig].Data(),cislo, ttl, pthmin);
		TH2D* hphi = (TH2D*) f1->Get(name);
		if(!hphi){ 
			cout << "Wrong pT for this file!" << endl;
			return ;
				 }

		name = Form("fHistXsections_%s_%d_tt_%d_pt%d", trg[ig].Data(), cislo, ttl, pthmin);
		TProfile* hxs = (TProfile*) f1->Get(name);
		if(!hxs){ 
			cout << "Wrong pT for this file!" << endl;
			return ;
				 }
		
		name = Form("fHistTrials_%s_%d_tt_%d_pt%d", trg[ig].Data(), cislo, ttl, pthmin);
		TH1F* hntr = (TH1F*) f1->Get(name);
		if(!hntr){ 
			cout << "Wrong pT for this file!" << endl;
			return ;
				 }

		Double_t xsection = hxs->GetMean(2);
		Double_t ntrials =  hntr->Integral();
		Double_t weight = 1;
		if (ntrials != 0) weight = xsection/ntrials;

		fhTTHPartLevel[ig]->Add(htt, weight);
		fhRecoilJetPtTTH_PartLevel[ig]->Add(hpt, weight);
		fhRecoilJetPhiTTH_PartLevel[ig]->Add(hphi, weight);
		delete htt;
		delete hpt;
		delete hphi;
		delete hxs;
		delete hntr;
	}
	cout << "File was added" << endl;
}

int GetSummarySpectrum(const Int_t& seed, const Int_t& TTLow){
	
	vector<Int_t> pthmin = {4, 11, 21, 36, 56, 84, 117, 156, 200, 249};  //hard bin

	Double_t pTbins3[]   = {-20,-15,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,100,120,140,160,180,200};
	const Int_t npTbins3 = sizeof(pTbins3)/sizeof(Double_t)-1;

	const Int_t narrPhi=128;
	Double_t arrPhi[narrPhi+1];
	Double_t p = TMath::TwoPi()/narrPhi;

	for(Int_t i=0; i<=narrPhi; i++){ 
		arrPhi[i] = i*p;
	}
	
	Int_t cislo = seed;
	Int_t trigRangeLow = TTLow;
	
	//cout << "Enter seed, TTLow: " << endl;
	//cin >> cislo >> trigRangeLow;

	for(Int_t ig = kMB; ig<kTRIG; ig++){
		name = Form("hTT_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		fhTTHPartLevel[ig] = new TH1D(name.Data(),name.Data(), 100, 0, 100);

		name = Form("fhRecoilJetPt_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		fhRecoilJetPtTTH_PartLevel[ig] = new TH1D(name.Data(), name.Data(), 200, -20, 180);

		name = Form("fhRecoilJetPhi_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		fhRecoilJetPhiTTH_PartLevel[ig] = new TH2D(name.Data(), name.Data(), npTbins3, pTbins3, narrPhi, arrPhi);
	}      


	TString filename = "";
	while (filename != "stop") {
		cout << "Enter filename: " << endl;
		cin >> filename;
		if (filename == "stop") {
			break;
		}
		
		std::string fname = static_cast<std::string>(filename);
		Int_t pos = fname.find("HB") + 2;
		Int_t endpos = fname.find("_", pos);
		std::string pt = fname.substr(pos, endpos - pos);
		Int_t ptval = std::stoi(pt);
		
		TFile *file = TFile::Open(filename,"READ");
	   	if(!file){
			cout << "Not a file!" << endl;
			break;
		}
		Sum (file, cislo, trigRangeLow, ptval);
	}
	
	name = Form("Res_tt_%d.root", trigRangeLow);
	TFile *outFile = new TFile(name, "RECREATE");
	outFile->cd();
		
	for(Int_t ig = kMB; ig<kTRIG; ig++){
		fhTTHPartLevel[ig]->Write();
		fhRecoilJetPtTTH_PartLevel[ig]->Write();
		fhRecoilJetPhiTTH_PartLevel[ig]->Write();
	}
	outFile->Close();
	
	for(Int_t igm = kMB; igm<kTRIG; igm++){
		delete fhTTHPartLevel[igm];
		delete fhRecoilJetPtTTH_PartLevel[igm];
		delete fhRecoilJetPhiTTH_PartLevel[igm];
	}
	delete outFile;
	return 0;
}

template<typename T>
T Scaler(T histo, TString& name, const Double_t& factor){
	histo->SetName(name);
	histo->Scale(factor);
	return histo;
}

void Scaling(Int_t trigRangeLow){
	name = Form("Res_tt_%d.root", trigRangeLow);
	TFile* f1 = TFile::Open(name, "UPDATE");
	name = Form("Res_scaled_tt_%d.root", trigRangeLow);
	TFile *outFile = new TFile(name, "RECREATE");
	outFile->cd();
	//Make counts for known luminosities: let us have 4 nb^-1 -> 4 nb^-1 = 4e6 mb^-1;
	//pp-luminosities: 3 pb^-1 = 3e9 mb^-1
	//OO-luminosities: 0.5 nb^-1; 1 nb^-1; 2nb^-1; 4nb^-1;
	//Xsection for OO-collisions: A(O) = 16; Xs(OO) = Xs(pp)*A^2;
	Double_t luminosity = 4e6; //4 nb^-1 -> 4e6 mb^-1;
	Double_t luminosity_pp = 3e9;
	Int_t A = 16; //oxygen mass in amu
	Double_t factorOO = A * A * luminosity;
	
	for(Int_t ig = kMB; ig< kTRIG; ig++){
		name = Form("hTT_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		TH1D* htt = (TH1D*) f1->Get(name);
		TH1D* htt_pp = (TH1D*) htt->Clone();
		name = Form("hTT_OO_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		htt = Scaler(htt, name, factorOO);
		htt->Write();

		name = Form("hTT_pp_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		htt_pp = Scaler(htt_pp, name, luminosity_pp);
		htt_pp->Write();


		name = Form("fhRecoilJetPt_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		TH1D* hpt = (TH1D*) f1->Get(name);
		TH1D* hpt_pp = (TH1D*) hpt->Clone();
		name = Form("fhRecoilJetPt_OO_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		hpt = Scaler(hpt, name, factorOO);
		hpt->Write();

		name = Form("fhRecoilJetPt_pp_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		hpt_pp = Scaler(hpt_pp, name, luminosity_pp);
		hpt_pp->Write();


		name = Form("fhRecoilJetPhi_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		TH2D* hphi = (TH2D*) f1->Get(name);
		TH2D* hphi_pp = (TH2D*) hphi->Clone();
		name = Form("fhRecoilJetPhi_OO_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		hphi = Scaler(hphi, name, factorOO);
		hphi->Write();

		name = Form("fhRecoilJetPhi_pp_%s_PartLevel_tt_%d", trg[ig].Data(), trigRangeLow);
		hphi_pp = Scaler(hphi_pp, name, luminosity_pp);
		hphi_pp->Write();

	}
	outFile->Close();
}

TH1D* MakePois(TFile* f1, const TString& name) {
	TRandom3 gener;
	TH1D* histo = (TH1D*) f1->Get(name);
	TH1D* Pois_hist = (TH1D*) histo->Clone();
	Int_t nbins = Pois_hist->GetNbinsX();
	for (Int_t in = 1; in <= nbins; in++){
		Int_t binc = Pois_hist->GetBinContent(in);
		binc = gener.TRandom::Poisson(binc);
		Pois_hist->SetBinContent(in, binc);
		Pois_hist->SetBinError(in, TMath::Sqrt(binc));
	}
	cout << "Poisson smearing for " << name << " is done!\n";
	return Pois_hist;
}

void DoPoisSmearing(Int_t trigRangeLow){
	name = Form("Res_scaled_tt_%d.root", trigRangeLow);
	TFile* f1 = TFile::Open(name, "UPDATE");
	name = Form("Res_pois_tt_%d.root", trigRangeLow);
	TFile *outFile = new TFile(name, "RECREATE");
	outFile->cd();
	
	TH1D* Pois_hist_TT[kTRIG]; 
	TH1D* Pois_hist_RJ[kTRIG]; 
	
	Int_t nbins = 0;
	TRandom3 gener;
	
	for(Int_t ig = kMB; ig< kTRIG; ig++){
		for(Int_t ij = kOO; ij< ktr; ij++){
			name = Form("hTT_%s_%s_PartLevel_tt_%d", RType[ij].Data(), trg[ig].Data(), trigRangeLow);
			Pois_hist_TT[ig] = MakePois(f1, name);
			Pois_hist_TT[ig]->Write();

			name = Form("fhRecoilJetPt_%s_%s_PartLevel_tt_%d", RType[ij].Data(), trg[ig].Data(), trigRangeLow);
			Pois_hist_RJ[ig] = MakePois(f1, name);
			Pois_hist_RJ[ig]->Write();
		
			name = Form("fhRecoilJetPhi_%s_%s_PartLevel_tt_%d", RType[ij].Data(), trg[ig].Data(), trigRangeLow);
			TH2D* hphi = (TH2D*) f1->Get(name);
			hphi->Write();
		}
	}
	outFile->Close();
}

vector<TH1D*> MakeVectorHist(vector<TH1D*>& hist_vector, TString name, const Int_t& trigRangeLow_1){
	TFile* f1 = TFile::Open(name, "UPDATE");
	
	vector<Double_t> integr_1;
	integr_1.resize(2);

	//get scaled recoil jet spectra 
	//get integral from trigger track spectra
	for(Int_t ig = kMB; ig< kTRIG; ig++){
		for(Int_t ij = kOO; ij< ktr; ij++){
			name = Form("hTT_%s_%s_PartLevel_tt_%d", RType[ij].Data(), trg[ig].Data(), trigRangeLow_1);
			TH1D* htt = (TH1D*) f1->Get(name);
			if(!htt){ 
				cout << "Wrong TT_Range for this file!" << trg[ig].Data() << 1 << endl;
				return {};
					 }
			integr_1[ig] = htt->Integral();

			//scale recoil jet spectrum to 1/integral
			name = Form("fhRecoilJetPt_%s_%s_PartLevel_tt_%d", RType[ij].Data(), trg[ig].Data(), trigRangeLow_1);
			TH1D* hpt = (TH1D*) f1->Get(name);
			if(!hpt){ 
				cout << "Wrong TT_Range for this file!" << trg[ig].Data() << 2 << endl;
				return {};
					 }
			TH1D* Scaled_hpt = (TH1D*) hpt->Clone();
			Scaled_hpt->Scale(1/integr_1[ig]);
			name = Form("Scaled_RecoilJetPt_%s_%s_PartLevel_tt_%d", RType[ij].Data(), trg[ig].Data(), trigRangeLow_1);
			Scaled_hpt->SetName(name);
			Scaled_hpt->Write();
			hist_vector.push_back(Scaled_hpt);

			delete htt;
			delete hpt;
			delete Scaled_hpt;
		}
	}
	f1->Close();
	return hist_vector;
}

Double_t appr_func(Double_t* x, Double_t* par) {
	Double_t slp = 1;
	if (par[1] != 0) { slp = par[1]; }
	Double_t fitval = par[0]*TMath::Exp(-x[0]/slp);
	return fitval;
}

void MakeDeltaRecoilSp(const Int_t& trigRangeLow_1, const Int_t& trigRangeLow_2){
	map<Int_t, Int_t> trigRanges = {{6, 7}, {12, 20}, {20, 30}};
	
	vector<TH1D*> hist_vector = {};
	
	name = Form("../TT-%d_%d/Res_pois_tt_%d.root", trigRangeLow_1, trigRanges.at(trigRangeLow_1), trigRangeLow_1);
	MakeVectorHist(hist_vector, name, trigRangeLow_1);
		
	name = Form("../TT-%d_%d/Res_pois_tt_%d.root", trigRangeLow_2, trigRanges.at(trigRangeLow_2), trigRangeLow_2);
	MakeVectorHist(hist_vector, name, trigRangeLow_2);
	
	
	name = Form("../DeltaRecoil_tt_%d-%d_tt_%d-%d.root", trigRangeLow_1, trigRanges.at(trigRangeLow_1), trigRangeLow_2, trigRanges.at(trigRangeLow_2));	
	TFile *outFile = new TFile(name, "RECREATE");
	outFile->Close();
		
	name = Form("../TT-%d_%d/Res_pois_tt_%d.root", trigRangeLow_1, trigRanges.at(trigRangeLow_1), trigRangeLow_1);
	TFile* f1 = TFile::Open(name, "READ");
	name = Form("../TT-%d_%d/Res_pois_tt_%d.root", trigRangeLow_2, trigRanges.at(trigRangeLow_2), trigRangeLow_2);
	TFile* f2 = TFile::Open(name, "READ");
	
	for(Int_t ig = kMB; ig< kTRIG; ig++){
		for(Int_t ij = kOO; ij< ktr; ij++){
			name = Form("Scaled_RecoilJetPt_%s_%s_PartLevel_tt_%d", RType[ij].Data(), trg[ig].Data(), trigRangeLow_1);
			TH1D* ScRJ_1 = (TH1D*) f1->Get(name);
			name = Form("Scaled_RecoilJetPt_%s_%s_PartLevel_tt_%d", RType[ij].Data(), trg[ig].Data(), trigRangeLow_2);
			TH1D* ScRJ_2 = (TH1D*) f2->Get(name);

			name = Form("Delta_RecoilJetPt_%s_%s_PartLevel_tt_%d-%d_tt_%d-%d", RType[ij].Data(),
						trg[ig].Data(),
						trigRangeLow_1,
						trigRanges.at(trigRangeLow_1),
						trigRangeLow_2,
						trigRanges.at(trigRangeLow_2));

			TH1D* DeltaRecoil = (TH1D*) ScRJ_1->Clone();
			DeltaRecoil->Add(ScRJ_1, ScRJ_2, -1, 1);
			DeltaRecoil->SetName(name);

			name = Form("../DeltaRecoil_tt_%d-%d_tt_%d-%d.root", trigRangeLow_1, trigRanges.at(trigRangeLow_1), trigRangeLow_2, trigRanges.at(trigRangeLow_2));
			outFile = TFile::Open(name, "UPDATE");
			DeltaRecoil->Write();
			outFile->Close();

			delete ScRJ_1;
			delete ScRJ_2;
			delete DeltaRecoil;
		}
	}
	f1->Close();
	f2->Close();
	outFile->Close();
	
	TF1* func = new TF1("expoconst", appr_func, -20,200, 2);
	func->SetParameters(1, 1);
	//TCanvas* c1 = new TCanvas("c1", "canvas", 800, 600);
	
	
	name = Form("../DeltaRecoil_tt_%d-%d_tt_%d-%d.root", trigRangeLow_1, trigRanges.at(trigRangeLow_1), trigRangeLow_2, trigRanges.at(trigRangeLow_2));
	outFile = TFile::Open(name, "UPDATE");
	name = Form("Delta_RecoilJetPt_%s_%s_PartLevel_tt_%d-%d_tt_%d-%d", RType[0].Data(),
						trg[0].Data(),
						trigRangeLow_1,
						trigRanges.at(trigRangeLow_1),
						trigRangeLow_2,
						trigRanges.at(trigRangeLow_2));
	TH1D* DRJ_OO = (TH1D*) outFile->Get(name);
	DRJ_OO->SetTitle(name);
	DRJ_OO->SetMarkerStyle(kDot);
	DRJ_OO->SetMarkerColor(4);
	
	TH1D* gs = (TH1D*) DRJ_OO->Clone();
	gs->Draw();
	DRJ_OO->Fit("expoconst");
	
	name = Form("Delta_RecoilJetPt_%s_%s_PartLevel_tt_%d-%d_tt_%d-%d", RType[1].Data(),
						trg[0].Data(),
						trigRangeLow_1,
						trigRanges.at(trigRangeLow_1),
						trigRangeLow_2,
						trigRanges.at(trigRangeLow_2));
	TH1D* DRJ_pp = (TH1D*) outFile->Get(name);
	DRJ_pp->SetTitle(name);
	DRJ_pp->SetMarkerStyle(kDot);
	DRJ_pp->SetMarkerColor(2);
	
	THStack* DRJ_stack = new THStack("Histogram stack","");
	DRJ_stack->Add(DRJ_OO);
	DRJ_stack->Add(DRJ_pp);
	DRJ_stack->Write();
	outFile->Close();
	
	delete func;
	delete DRJ_stack;
}



