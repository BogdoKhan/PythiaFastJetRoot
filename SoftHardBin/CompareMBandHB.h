#include "TROOT.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TF1.h"
#include "TLine.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include <stdio.h>
#include <stdlib.h>
#include "TColor.h"
#include <iostream>
#include <string>
using namespace std;

using std::cout;
using std::endl;

void SetCanvas(TCanvas* c);
void SetHist(TH1* h,TString titx, TString tity);

void SavePNGandEPS(TCanvas** arrayCan, Int_t nobr);

Int_t CompareMBandHB(vector<string> argv);