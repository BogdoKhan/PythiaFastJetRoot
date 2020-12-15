#include "TROOT.h"

void dummy(){
	TFile* MyFile = new TFile("PP7_pythia_ANTIKT04_tune14_c3866185.root","READ");
	gROOT->ProcessLine("gDirectory->ls()");
	int i;
	cin >> i;
	if (i == 1){
		gROOT->ProcessLine("hTrackPt->Draw()");
	}
	else if (i == 2) {
		gROOT->ProcessLine("hJetPt->Draw()");
	}
}