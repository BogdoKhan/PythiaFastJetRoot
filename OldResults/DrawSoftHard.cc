#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include <iostream>
#include <string>

using namespace std;

int main(){
	TCanvas* c1 = new TCanvas;
	gPad->SetLogy();
	
	string str = "";
	cout << "Enter whether histogram you would like to see: \n"
		<< "1 - Hard bin, pT;\n" << "2 - Soft bin, pT;\n" 
		<< "3 - both of them in the same canvas\n"
		<< ".q - quit the script\n";
	cin >> str;
	
	while (str != ".q"){
		if (str == "1"){
			TFile* f = new TFile("./Results/Hard_PYAkT04_tune14_c5092.root");
			f->cd();
			TH1D* hHard = (TH1D*) f->Get("hTrackPtsum");
			hHard->Draw();
			c1->Update();
		}
		else if (str == "2"){	
			TFile* a = new TFile("./Results/Soft_PYAkT04_tune14_c5092.root");
			a->cd();
			TH1D* hSoft = (TH1D*) a->Get("hTrackPt");
			hSoft->Draw();
			c1->Update();
		}
		else if (str == "3"){
			TFile* f = new TFile("./Results/Hard_PYAkT04_tune14_c5092.root");
			f->cd();
			TH1D* hHard = (TH1D*) f->Get("hTrackPtsum");
			hHard->Draw();
			c1->Update();

			TFile* a = new TFile("./Results/Soft_PYAkT04_tune14_c5092.root");
			a->cd();
			TH1D* hSoft = (TH1D*) a->Get("hTrackPt");
			hSoft->Draw("same");
			c1->Update();
		}
		else {
			cout << "Wrong input! Type .q to exit" << endl;
		}
		cin >> str;
	}
	return 0;
}