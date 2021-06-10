#include "TH1D.h"
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

#include <vector>
#include <iostream>
using std::cout;
using std::cin;
using std::endl;

enum {kMB, kFIT, kTRIG};  //kMB = pythia min bias,   kFIT =  event triggered with FIT coincidence
TString trg[]={"MB","FIT"};  

TString name = "";

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

int pois(const Int_t& seed, const Int_t& TTLow){
	
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
		name = Form("hTT_%s_PartLevel_seed%d_tt_%d", trg[ig].Data(), cislo, trigRangeLow);
		fhTTHPartLevel[ig] = new TH1D(name.Data(),name.Data(), 100, 0, 100);

		name = Form("fhRecoilJetPt_%s_PartLevel_seed%d_tt_%d", trg[ig].Data(), cislo, trigRangeLow);
		fhRecoilJetPtTTH_PartLevel[ig] = new TH1D(name.Data(), name.Data(), 200, -20, 180);

		name = Form("fhRecoilJetPhi_%s_PartLevel_seed%d_tt_%d", trg[ig].Data(),cislo, trigRangeLow);
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
	
	TFile *outFile = new TFile("Res.root", "RECREATE");
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


void Scaling(TString fname, Int_t cislo, Int_t trigRangeLow){
	TFile* f1 = TFile::Open(fname, "UPDATE");
	TFile *outFile = new TFile("Res_scaled.root", "RECREATE");
	outFile->cd();
	//Make counts for known luminosities: let us have 4 nb^-1 -> 4 nb^-1 = 4e6 mb^-1;
	//pp-luminosities: 3 pb^-1 = 3e9 mb^-1
	//OO-luminosities: 0.5 nb^-1; 1 nb^-1; 2nb^-1; 4nb^-1;
	//Xsection for OO-collisions: A(O) = 16; Xs(OO) = Xs(pp)*A^2;
	Int_t luminosity = 4e6; //4 nb^-1 -> 4e6 mb^-1;
	Int_t A = 16; //oxygen mass in amu
	Int_t factorOO = A * A * luminosity;

	
	TH1D* Pois_hist[kTRIG]; 
	
	Int_t nbins = 0;
	TRandom3 gener;
	
	for(Int_t ig = kMB; ig< kTRIG; ig++){
		
		name = Form("hTT_%s_PartLevel_seed%d_tt_%d", trg[ig].Data(), cislo, trigRangeLow);
		TH1D* htt = (TH1D*) f1->Get(name);
		htt->Scale(factorOO);
		htt->Write();

		name = Form("fhRecoilJetPt_%s_PartLevel_seed%d_tt_%d", trg[ig].Data(), cislo, trigRangeLow);
		TH1D* hpt = (TH1D*) f1->Get(name);
		hpt->Scale(factorOO);
		hpt->Write();
		
		Pois_hist[ig] = (TH1D*) hpt->Clone();
		nbins = Pois_hist[ig]->GetNbinsX();
		for (Int_t in = 0; in < nbins; in++){
			Int_t binc = Pois_hist[ig]->GetBinContent(in);
			binc = gener.TRandom::Poisson(binc);
			Pois_hist[ig]->SetBinContent(in, binc);
			Pois_hist[ig]->SetBinError(in, TMath::Sqrt(binc));
		}
		name = Form("Pois_fhRecoilJetPt_%s_PartLevel_seed%d_tt_%d", trg[ig].Data(), cislo, trigRangeLow);
		Pois_hist[ig]->SetName(name.Data());
		Pois_hist[ig]->Write();

		name = Form("fhRecoilJetPhi_%s_PartLevel_seed%d_tt_%d", trg[ig].Data(),cislo, trigRangeLow);
		TH2D* hphi = (TH2D*) f1->Get(name);
		hphi->Scale(factorOO);
		hphi->Write();
	}
	
	
}




