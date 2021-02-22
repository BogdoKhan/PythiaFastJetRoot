#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"

#include <iostream>
using std::cout;
using std::endl;

void SetCanvas(TCanvas* c);
void SetHist(TH1* h,TString titx, TString tity);
void SavePNGandEPS(TCanvas** arrayCan, Int_t nobr);

Int_t Draw(){

   //OPEN FIRST ROOT FILE	
   TFile *f1 = new TFile("./Results/PP7_pth_akT_kT04_tune14_c5024.root","READ");
   if(!f1) return 20; 
   //READ HISTOGRAM 
   TH1D*  hard = (TH1D*) f1->Get("hard_hJetXsect");
   //TH1D*  hard = (TH1D*) f1->Get("hard_hJetPt");
   if(!hard) return 21; 

   //READ HISTOGRAM 
   TH1D*  soft = (TH1D*) f1->Get("soft_hJetPt");
   //TH1D*  soft = (TH1D*) f2->Get("soft_hJetPt");
   if(!soft) return 21; 

   //----------------------------------------------------------
   gStyle->SetOptTitle(0);
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
   SavePNGandEPS((TCanvas**) c, io);
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


