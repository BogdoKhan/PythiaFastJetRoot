/*#include "TH1D.h"
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

using std::cout;
using std::endl;

void SetCanvas(TCanvas* c);
void SetHist(TH1* h,TString titx, TString tity);

void SavePNGandEPS(TCanvas** arrayCan, Int_t nobr);*/
#include "VisualizeResults.h"


Int_t VisRes(vector<string> argv){
	TString seed = argv[1];
   TFile *fSoft = new TFile("./Results/HardSoftBins/Soft_PYAkT04_tune14_c" + seed + ".root", "READ");
   if(!fSoft) return 1;
   TH1D* hTrackPtSoft = (TH1D*) fSoft->Get("hTrackPt");
   if(!hTrackPtSoft) return 2;
   hTrackPtSoft->SetLineColor(2);
   hTrackPtSoft->Rebin(2);

   TFile *fHard = new TFile("./Results/HardSoftBins/Hard_PYAkT04_tune14_c" + seed + ".root", "READ");
   if(!fHard) return 1;
   TH1D* hTrackPtHard = (TH1D*) fHard->Get("hTrackPtsum");
   if(!hTrackPtHard) return 2;
   hTrackPtHard->SetLineColor(4);
   hTrackPtHard->Rebin(2);

   int io=0;
   TCanvas *c[500];
   TLegend *leg=NULL;
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);

   c[io] = new TCanvas(Form("c%d",io),"PtTrackSpectra",600,500);
   SetCanvas((TCanvas*) c[io]);
   c[io]->SetLogy();
   SetHist((TH1D*) hTrackPtSoft, "p_{T}^{track} (GeV/c)", "cross section (mb)");

   hTrackPtSoft->SetMinimum(1e-6);
   hTrackPtSoft->Draw();
   hTrackPtHard->Draw("same");

   //legend
   leg = new TLegend(0.5,0.75,0.9,0.95," ","brNDC");
   leg->SetTextSize(0.045);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->AddEntry((TObject*) hTrackPtSoft, "minimum bias","l");
   leg->AddEntry((TObject*) hTrackPtHard, "hard bins","l");
   leg->Draw();

   io++;
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //Ratio  Hard bins / MinBias 

   TH1D *hRatio = (TH1D*) hTrackPtHard->Clone("hRatio");
   hRatio->Divide((TH1D*)hTrackPtSoft);
   c[io] = new TCanvas(Form("c%d",io),"RatioPtTrackSpectra",600,500);
   SetCanvas((TCanvas*) c[io]);
   SetHist((TH1D*) hRatio, "p_{T}^{track} (GeV/c)", "hard bins / min bias");
   hRatio->Draw();

   TLine *l =new TLine(0,1,50,1);
   l->SetLineStyle(2);
   l->Draw();

   io++;




   //Print on screen a lines which copied and pasted to teminal makes png figure
   SavePNGandEPS((TCanvas**) c, io);

   return 11;
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

void SetHist(TH1* h,TString titx, TString tity){

   h->GetXaxis()->SetTitle(titx.Data());
   h->GetYaxis()->SetTitle(tity.Data());
   h->GetXaxis()->SetTitleSize(0.05);
   h->GetYaxis()->SetTitleSize(0.05);
   h->GetXaxis()->SetLabelSize(0.04);
   h->GetYaxis()->SetLabelSize(0.04);

   h->GetYaxis()->SetTitleOffset(1.5);
   h->GetXaxis()->SetTitleOffset(1.2);

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
	TString OutImg = "img->FromPad(" + (TString) arrayCan[j]->GetName() + "); img->WriteImage(\"" + "./Results/HardSoftBins/" + namePNG.Data() + "\");";
	gROOT->ProcessLine(OutImg);
   }
 
 }

