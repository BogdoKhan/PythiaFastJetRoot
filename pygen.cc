// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

//THREAD-SAFE RANDOM NUMBER GENERATOR FOR PYTHIA SEED

#if defined (_MSC_VER)  // Visual studio
    #define thread_local __declspec( thread )
#elif defined (__GCC__) // GCC
    #define thread_local __thread
#endif

#include <random>
#include <time.h>
#include <thread>
//OTHER LIBRARIES

#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"  
#include <ctime>
#include <future>
#include <vector>
#include <sys/stat.h>
#include <ctime>
#include <cstring>
#include <string>
#include <mutex>
#include <functional>
#include <fstream>
#include "TROOT.h"


using namespace std;
using namespace Pythia8;
//___________________________________________________________________
//               FUNCTIONS


bool EtaCut(fastjet::PseudoJet fjJet, double etaMin, double etaMax);
void MakeDir(const string& dirpath, const string& dirname);
//___________________________________________________________________

//_________________________________________________________________________ 

void MakeDir(const string& dirpath, const string& dirname){
	string path = dirpath + dirname;
	const char* folder = path.c_str();
	//folder = "./Results";
	struct stat sb;

	if (stat(folder, &sb) != 0 && !S_ISDIR(sb.st_mode)) {
		mkdir(folder, 0755);
	} 
}

//RANDOM INT NUMBER GENERATOR WITH THREAD SAFETY
//WARNING - POTENTIAL MEMORY LEAK: NEED TO STORE/DELETE GENERATORS

int intRand(const int & min, const int & max, vector<mt19937*>& gens) {
	mutex addGen_excl;
    static thread_local mt19937* generator = nullptr;
    if (!generator) generator = new mt19937(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    uniform_int_distribution<int> distribution(min, max);
	lock_guard<mutex> lock_Gen(addGen_excl);
	gens.push_back(generator);
    return distribution(*generator);
}

enum {kMB, kFIT, kTRIG};  //kMB = pythia min bias,   kFIT =  event triggered with FIT coincidence
TString trg[]={"MB","FIT"};  

class TH1Coll {
public:
	//vector<TH1D> h1t(2) = [fhTTH_mb],[fhTTH_FIT]
	//vector<TH1D> h1h(2) = [fhRec-pT_mb], [fhRec-pt_FIT]
	//vector<TH2D> h2(2) = [fhRec-phi_mb],[fhRec-phi_FIT]
	//vector<TProfile> hXsect(2) = [Xsect_mb],[Xsect_FIT]
	//vector<TH1F> hNtrials(2) = [Ntrials_mb], [Ntrials_FIT]
	
	vector<TH1D> h1t;
	vector<TH1D> h1h;
	vector<TH2D> h2;
	vector<TProfile> hXsect;
	vector<TH1F> hNtrials;
	
	TH1D* fhTTH_PartLevel[kTRIG];
	TH1D* fhRecoilJetPtTTH_PartLevel[kTRIG];
	TH2D* fhRecoilJetPhiTTH_PartLevel[kTRIG]; 
	TProfile* hXsect_part[kTRIG];
	TH1F* hNtrials_part[kTRIG];
	
	TH1Coll():
	_ttl(-1), _tth(-1), _ptmin(-1), _ptmax(-1), _nseed(0)
	{
		MakeCol();
	}
	
	TH1Coll(const int& ttl,
			const int& tth,
			const int& ptmin,
			const int& ptmax,
		   	const int& nseed) :
	_ttl(ttl), _tth(tth),
	_ptmin(ptmin), _ptmax(ptmax),
	_nseed(nseed){
		MakeCol();
	}

	~TH1Coll() {
		/*for(Int_t ig = kMB; ig<kTRIG; ig++){
			delete[] fhTTH_PartLevel[ig];
			delete[] fhRecoilJetPtTTH_PartLevel[ig];
			delete[] fhRecoilJetPhiTTH_PartLevel[ig];
		}*/
	}
	
	int sttl() const {
		return _ttl;
	}
	
	int stth() const {
		return _tth;
	}
	
	int sptmin() const {
		return _ptmin;
	}
	
	int sptmax() const {
		return _ptmax;
	}
	
	int snseed() const {
		return _nseed;
	}
	
   
      
private:
	int _ttl;
	int _tth;
	int _ptmin;
	int _ptmax;
	int _nseed;
	
	void MakeCol() {
		TString name;
		Double_t pTbins3[]   = {-20,-15,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,100,120,140,160,180,200};
		const Int_t npTbins3 = sizeof(pTbins3)/sizeof(Double_t)-1;

		const Int_t narrPhi=128;
		Double_t arrPhi[narrPhi+1];
		Double_t p = TMath::TwoPi()/narrPhi;
		for(Int_t i=0; i<=narrPhi; i++) arrPhi[i] = i*p;


		for(Int_t ig = kMB; ig<kTRIG; ig++){
			name = Form("hTT_%s_PartLevel_seed_%d_t_%d_p%d", trg[ig].Data(),_nseed, _ttl, _ptmin);
			fhTTH_PartLevel[ig] = new TH1D(name.Data(),name.Data(), 100, 0, 100);
			fhTTH_PartLevel[ig]->Sumw2();
			h1t.push_back(*fhTTH_PartLevel[ig]);
			
			//HJET SPECTRA
			name = Form("fhRecoilJetPt_%s_PartLevel_seed_%d_t_%d_p%d", trg[ig].Data(), _nseed, _ttl, _ptmin);
			fhRecoilJetPtTTH_PartLevel[ig] = new TH1D(name.Data(), name.Data(), 200, -20, 180);
			fhRecoilJetPtTTH_PartLevel[ig]->Sumw2();
			h1h.push_back(*fhRecoilJetPtTTH_PartLevel[ig]);
			
			//HJET_Phi
			name = Form("fhRecoilJetPhi_%s_PartLevel_seed_%d_t_%d_p%d", trg[ig].Data(), _nseed, _ttl, _ptmin);
			fhRecoilJetPhiTTH_PartLevel[ig] = new TH2D(name.Data(), name.Data(), npTbins3, pTbins3, narrPhi, arrPhi);
			fhRecoilJetPhiTTH_PartLevel[ig]->Sumw2();
			h2.push_back(*fhRecoilJetPhiTTH_PartLevel[ig]);
			
			//Profile for Xsection
			name = Form("hXsect_%s_PartLevel_seed_%d_t_%d_p%d", trg[ig].Data(), _nseed, _ttl, _ptmin);
			hXsect_part[ig] = new TProfile(name.Data(), name.Data(), 1, 0, 1);
   			hXsect_part[ig]->GetYaxis()->SetTitle("xsection");
			hXsect.push_back(*hXsect_part[ig]);
			
			//Add Ntrials
			name = Form("hNtrials_%s_PartLevel_seed_%d_t_%d_p%d", trg[ig].Data(), _nseed, _ttl, _ptmin);
  			hNtrials_part[ig]= new TH1F(name.Data(), name.Data(), 1, 0, 1);
   			hNtrials_part[ig]->GetYaxis()->SetTitle("trials");
			hNtrials.push_back(*hNtrials_part[ig]);
		}
	}
};

TH1Coll operator+(const TH1Coll& lhs, const TH1Coll& rhs){
	/*TH1Coll* thcnew = new TH1Coll(lhs.sttl(), lhs.stth(),
								  lhs.sptmin(), lhs.sptmax(),
								  lhs.snseed());*/
	TH1Coll thcnew = lhs;
	for (size_t i = 0; i <= 1; i++) {
		thcnew.h1t[i] = lhs.h1t[i] + rhs.h1t[i];
		thcnew.h1h[i] = lhs.h1h[i] + rhs.h1h[i];
		thcnew.h2[i].Add(&(lhs.h2[i]), &(rhs.h2[i]),1,1);
		thcnew.hXsect[i].Add(&(lhs.hXsect[i]), &(rhs.hXsect[i]),1,1);
		thcnew.hNtrials[i].Add(&(lhs.hNtrials[i]), &(rhs.hNtrials[i]),1,1);
	}
	return thcnew;
}


TH1Coll Run(int argc, vector<string> argv) {


   Int_t cislo = -1;                 //unique number for each file
   Int_t tune  = -1;                 //pythia tune
   Int_t trigRangeLow  = 6;
   Int_t trigRangeHigh = 7;
   Int_t ptHatMin  = -1;
   Int_t ptHatMax  = -1;
   Float_t jr=-1;
	vector<mt19937*> GenStorage;

   if(argc!=8){  
      cout<<"Usage:"<<endl<<"./pygen <PythiaTune> <Number> <TTlow> <TT high>  <pTHat min> <pTHat max> <jetR>"<<endl;
	   TH1Coll* dummy = new TH1Coll;
      return *dummy;
   }
   tune  = stoi(argv[1]);
   cislo = stoi(argv[2]);
   trigRangeLow  = stoi(argv[3]);
   trigRangeHigh = stoi(argv[4]);
   ptHatMin      = stoi(argv[5]);
   ptHatMax      = stoi(argv[6]);
   jr            = stof(argv[7]);
	
	TH1Coll* EvCollection = new TH1Coll(trigRangeLow, trigRangeHigh, ptHatMin, ptHatMax, cislo);
	


   Int_t nEvent=(Int_t) 1e1;
   //Int_t nEvent=(Int_t) 2e5;
  
   TString name; 
   //______________________________________

   //random for single TT selection
   TRandom3*  random = new TRandom3(cislo);
   random->SetSeed(0);
	int randseed = intRand(0,100000000, GenStorage);

   //__________________________________________________________________________
   //                        ANALYSIS SETTINGS

   double jetParameterR   = jr; //jet R
   double trackLowPtCut   = 0.15; //GeV
   double trackEtaCut     = 0.9;
   Float_t fkDeltaPhiCut  = TMath::Pi()-0.6;// rad
   //__________________________________________________________________________
   //                        PYTHIA SETTINGS

   // Generator. Process selection. LHC initialization. Histogram.
   Pythia pythia;
   pythia.readString("Beams:idA = 2212"); //beam 1 proton
   pythia.readString("Beams:idB = 2212"); //beam 2 proton
   pythia.readString("Beams:eCM = 6350.");
   //p+pb
   pythia.readString(Form("Tune:pp = %d",tune));  //tune 1-13    5=defaulr TUNE4C,  6=Tune 4Cx, 7=ATLAS MB Tune A2-CTEQ6L1

   pythia.readString("Random:setSeed = on");
   //pythia.readString(Form("Random:seed = %d",cislo));
   pythia.readString(Form("Random:seed = %d",randseed));
	
	std::ofstream ofs;
	ofs.open ("test.txt", std::ofstream::out | std::ofstream::app);
	ofs << randseed << endl;
	ofs.close();

   //pythia.readString("SoftQCD:inelastic = on");
   pythia.readString("HardQCD:all = on");
   //pythia.readString("SoftQCD:all = on");

   if(ptHatMin<1 || ptHatMax <1){     
      pythia.readString("PhaseSpace:pTHatMin = 0."); // <<<<<<<<<<<<<<<<<<<<<<<
   }else{
      name = Form("PhaseSpace:pTHatMin = %f", (Float_t) ptHatMin);
      pythia.readString(name.Data()); 
      name = Form("PhaseSpace:pTHatMax = %f", (Float_t) ptHatMax);
      pythia.readString(name.Data()); 
   }
//   pythia.settings.addWord("ANALYSIS:outputFile"   ,"output.root" );
  // string outputFile = pythia.settings.word("ANALYSIS:outputFile");

   pythia.readString("310:mayDecay  = off"); //K0s
   pythia.readString("3122:mayDecay = off"); //labda0
   pythia.readString("3112:mayDecay = off"); //sigma-
   pythia.readString("3212:mayDecay = off"); //sigma0
   pythia.readString("3222:mayDecay = off"); //sigma+
   pythia.readString("3312:mayDecay = off"); //xi-
   pythia.readString("3322:mayDecay = off"); //xi+
   pythia.readString("3334:mayDecay = off"); //omega-



   pythia.init();
 
   //___________________________________________________ 
   //                      FAST JET  SETTINGS

   double jetEtaMin = - trackEtaCut + jetParameterR; //signal jet eta range
   double jetEtaMax = - jetEtaMin;


   fastjet::Strategy strategy = fastjet::Best;
   fastjet::RecombinationScheme recombScheme = fastjet::BIpt_scheme;

   fastjet::JetDefinition *jetDefAKT = NULL;
   fastjet::JetDefinition *jetDefKT = NULL;

   jetDefAKT = new fastjet::JetDefinition(fastjet::antikt_algorithm, 
                                              jetParameterR, 
                                              recombScheme, 
                                              strategy);

   jetDefKT = new fastjet::JetDefinition(fastjet::kt_algorithm, 
                                              jetParameterR, 
                                              recombScheme, 
                                              strategy);
 

 
   fastjet::GhostedAreaSpec ghostareaspec(trackEtaCut, 1, 0.01); //ghost max rap, repeat, ghostarea default 0.01
   fastjet::AreaType areaType = fastjet::active_area_explicit_ghosts;
   fastjet::AreaDefinition *areaDef = new fastjet::AreaDefinition(areaType, ghostareaspec);

   // Fastjet input
   std::vector<fastjet::PseudoJet> fjInputsAKT;
   std::vector<fastjet::PseudoJet> fjInputsKT;

   std::vector<Int_t> fTrigTracksGen; //index of trigger tracks candidates
   //___________________________________________________ 
   //FIT DETECTOR
   Float_t FIT_A_low = 2.2;
   Float_t FIT_A_up  = 5.0;
 
   Float_t FIT_C_low = -3.4;
   Float_t FIT_C_up  = -2.3;
   //___________________________________________________ 
   //                HISTOGRAMS
   //enum {kMB, kFIT, kTRIG};  //kMB = pythia min bias,   kFIT =  event triggered with FIT coincidence
   //TString trg[]={"MB","FIT"};  

   //TT SPECTRA
	TH1D* fhTTH_PartLevel[kTRIG];
	TH1D* fhRecoilJetPtTTH_PartLevel[kTRIG];
	TH2D* fhRecoilJetPhiTTH_PartLevel[kTRIG]; 
	TProfile* hXsect_part[kTRIG];
	TH1F* hNtrials_part[kTRIG];

   Double_t pTbins3[]   = {-20,-15,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,100,120,140,160,180,200};
   const Int_t npTbins3 = sizeof(pTbins3)/sizeof(Double_t)-1;
   
   const Int_t narrPhi=128;
   Double_t arrPhi[narrPhi+1];
   Double_t p = TMath::TwoPi()/narrPhi;
   for(Int_t i=0; i<=narrPhi; i++) arrPhi[i] = i*p;
  

   for(Int_t ig = kMB; ig<kTRIG; ig++){
		name = Form("hTT_%s_PartLevel_seed%d_tt_%d_pt%d", trg[ig].Data(), cislo, trigRangeLow, ptHatMin);
		fhTTH_PartLevel[ig] = new TH1D(name.Data(),name.Data(), 100, 0, 100);
		fhTTH_PartLevel[ig]->Sumw2();

		//HJET SPECTRA
		name = Form("fhRecoilJetPt_%s_PartLevel_seed%d_tt_%d_pt%d", trg[ig].Data(), cislo, trigRangeLow, ptHatMin);
		fhRecoilJetPtTTH_PartLevel[ig] = new TH1D(name.Data(), name.Data(), 200, -20, 180);
		fhRecoilJetPtTTH_PartLevel[ig]->Sumw2();

		name = Form("fhRecoilJetPhi_%s_PartLevel_seed%d_tt_%d_pt%d", trg[ig].Data(),cislo, trigRangeLow, ptHatMin);
		fhRecoilJetPhiTTH_PartLevel[ig] = new TH2D(name.Data(), name.Data(), npTbins3, pTbins3, narrPhi, arrPhi);
		fhRecoilJetPhiTTH_PartLevel[ig]->Sumw2();

		name = Form("fHistXsections_%s_%d_tt_%d_pt%d", trg[ig].Data(), cislo, trigRangeLow, ptHatMin);
		hXsect_part[ig] = new TProfile(name.Data(), name.Data(), 1, 0, 1);
		hXsect_part[ig]->GetYaxis()->SetTitle("xsection");

		name = Form("fHistTrials_%s_%d_tt_%d_pt%d", trg[ig].Data(), cislo, trigRangeLow, ptHatMin);
		hNtrials_part[ig] = new TH1F(name.Data(), name.Data(), 1, 0, 1);
		hNtrials_part[ig]->GetYaxis()->SetTitle("trials");
   }      


   //PT HARD
   //___________________________________________________ 
   Bool_t bAtLeastOneTrigger = kFALSE; //was there a high pT particle?
   Bool_t bFIT = kFALSE; //was the even triggered by FIT coincidence?
   
   Double_t ptTriggGen=-1.0;

   Double_t dphi;
   Double_t ptjetcorr; 


   Int_t  idxtrig; 

   Double_t rho  = 0.;
   Int_t nJetAcc = 0;
   
  
   Double_t ptLJ = -1.;
   Double_t ptSJ = -1.;
   unsigned int iLJ = -1;
   unsigned int iSJ = -1;
   Double_t frhovec[999];

   Int_t fFIT_A = 0; 
   Int_t fFIT_C = 0;
   Bool_t bfill[2];

   //___________________________________________________ 
   // Begin event loop. Generate event. Skip if error. List first one.
   for(int iEvent = 0; iEvent < nEvent; iEvent++){
      if(!pythia.next()) continue;

      //          INITIALIZATION
      fjInputsAKT.clear();
      fjInputsKT.clear();
      fTrigTracksGen.clear();
      bAtLeastOneTrigger = kFALSE;
      bFIT               = kFALSE;

      fFIT_A = 0; 
      fFIT_C = 0;
      //______________________________
      // Fill charged jets and search for trigger tracks
  

      for(Int_t i = 0; i < pythia.event.size(); ++i){
         if(pythia.event[i].isFinal() && pythia.event[i].isCharged()){

            if(FIT_C_low < pythia.event[i].eta() && pythia.event[i].eta() < FIT_C_up) fFIT_C++;
            if(FIT_A_low < pythia.event[i].eta() && pythia.event[i].eta() < FIT_A_up) fFIT_A++;

            if(TMath::Abs(pythia.event[i].eta()) > trackEtaCut) continue;      //eta cut
            

            if(trigRangeLow <= ((Float_t) pythia.event[i].pT()) && 
                               ((Float_t) pythia.event[i].pT()) < trigRangeHigh){
                   fTrigTracksGen.push_back(i);
                   bAtLeastOneTrigger = kTRUE;
            }
         }
      }

      //--------------------------------------------------------
      if(!bAtLeastOneTrigger) continue; //skip events without high pT trigger


      //CONINCIDENCE CONDITION ?       
      if(fFIT_A > 0 && fFIT_C > 0) bFIT = kTRUE; 
      bfill[kMB]  = 1;
      bfill[kFIT] = bFIT;


      //FILL JET CONSTITUENTS
      for(Int_t i = 0; i < pythia.event.size(); ++i){
         if(pythia.event[i].isFinal() && pythia.event[i].isCharged()){

            if(TMath::Abs(pythia.event[i].eta()) > trackEtaCut) continue;      //eta cut

            if(pythia.event[i].pT() < trackLowPtCut) continue;

            fjInputsAKT.push_back( fastjet::PseudoJet(pythia.event[i].px(),
                                                   pythia.event[i].py(),
                                                   pythia.event[i].pz(),
                                                   pythia.event[i].pAbs()));

            fjInputsKT.push_back( fastjet::PseudoJet(pythia.event[i].px(),
                                                   pythia.event[i].py(),
                                                   pythia.event[i].pz(),
                                                   pythia.event[i].pAbs()));

         }
      }
       //fjInputs.size(); 
      //--------------------------------------------------------
      vector<fastjet::PseudoJet> inclusiveJetsAKT;
      vector<fastjet::PseudoJet> inclusiveJetsKT;
      fastjet::ClusterSequenceArea clustSeqAKT(fjInputsAKT, *jetDefAKT, *areaDef);
      fastjet::ClusterSequenceArea clustSeqKT(fjInputsKT, *jetDefKT, *areaDef);

      inclusiveJetsAKT = clustSeqAKT.inclusive_jets(0.15); //list of jets pT > 150 MeV
      inclusiveJetsKT  = clustSeqKT.inclusive_jets(0.);   //list of jets pT > 0 

      //______________________________________________________________________________
      //  CALCULATE RHO

      rho = 0.;
      nJetAcc = 0;
     
      
      ptLJ=-1;
      ptSJ=-1;
      iLJ = -1;
      iSJ = -1;

      //Exclude 2 leading jets
      for(unsigned int ijet = 0; ijet < inclusiveJetsKT.size(); ijet++){
         fastjet::PseudoJet fjJet = inclusiveJetsKT.at(ijet);
         //vector < fastjet::PseudoJet > constit = clustSeqAKT.constituents(fjJet);
         if(!EtaCut(fjJet, jetEtaMin, jetEtaMax)) continue;
         if(fjJet.pt()   < 0.150) continue;

         if(fjJet.pt() > ptLJ){
            ptSJ  = ptLJ;
            iSJ   = iLJ;
            
            ptLJ  = fjJet.pt();
            iLJ   = ijet;
         }else if(fjJet.pt() > ptSJ){
            ptSJ  = fjJet.pt();
            iSJ   = ijet;
         }
      }

      for(unsigned int ijet = 0; ijet < inclusiveJetsKT.size(); ijet++){
         fastjet::PseudoJet fjJet = inclusiveJetsKT.at(ijet);
         //vector < fastjet::PseudoJet > constit = clustSeqAKT.constituents(fjJet);
         if(!EtaCut(fjJet, jetEtaMin, jetEtaMax)) continue;

         if(ijet==iLJ) continue; //skip two leading kT jets
         if(ijet==iSJ) continue;

         frhovec[nJetAcc]  = fjJet.pt()/fjJet.area();
         nJetAcc++;
      }

      if(nJetAcc>0){
         rho = TMath::Median(nJetAcc, frhovec);
      }


      //______________________________________________________________________________
      //               ASSIGN TRIGGER HADRON
      Int_t nn   =(Int_t)(fTrigTracksGen.size());
      Int_t rnd  = random->Integer(nn); //0 to ntriggers-1
      idxtrig  = fTrigTracksGen[rnd];


      if( idxtrig >=0 ){ 

         if(pythia.event[idxtrig].isFinal() &&
            pythia.event[idxtrig].isCharged() &&
            TMath::Abs(pythia.event[idxtrig].eta()) < trackEtaCut){

            ptTriggGen = pythia.event[idxtrig].pT(); //TT  pT

			for(Int_t ig=0; ig< kTRIG; ig++){
				   if(!bfill[ig]) continue;
				   fhTTH_PartLevel[ig]->Fill(ptTriggGen);
			}

            //Count jets and trigger-jet pairs at MC  generator level
            for(unsigned int ijet = 0; ijet < inclusiveJetsAKT.size(); ijet++){
               fastjet::PseudoJet fjJet = inclusiveJetsAKT.at(ijet);
               //vector < fastjet::PseudoJet > constit = clustSeqAKT.constituents(fjJet);
               if(!EtaCut(fjJet, jetEtaMin, jetEtaMax)) continue;
               if(fjJet.pt()   < 0.150) continue;
   
               dphi = TVector2::Phi_0_2pi(pythia.event[idxtrig].phi() - fjJet.phi());
               ptjetcorr = fjJet.pt()  - rho * fjJet.area();

			   for(Int_t ig=0; ig< kTRIG; ig++){
					if(!bfill[ig]) continue; 
					fhRecoilJetPhiTTH_PartLevel[ig]->Fill(ptjetcorr, dphi);
			   }

			   if(TMath::Abs((Float_t) TVector2::Phi_mpi_pi(dphi)) < fkDeltaPhiCut) continue;
			   for(Int_t ig=0; ig< kTRIG; ig++){
				  	if(!bfill[ig]) continue;
				  	fhRecoilJetPtTTH_PartLevel[ig]->Fill(ptjetcorr);
			   }

            }//jet loop
         }//etacut
      }//TT
   }// End of event loop.
   for(Int_t ig=0; ig< kTRIG; ig++){
		if(!bfill[ig]) continue;
		hXsect_part[ig]->Fill(0.5,pythia.info.sigmaGen());
		hNtrials_part[ig]->Fill(0.5, pythia.info.nAccepted());
   }

   //____________________________________________________
   //          SAVE OUTPUT
 
   for(Int_t ig=0; ig< kTRIG; ig++){
      //TT SPECTRA
	   (*EvCollection).h1t[ig] = *fhTTH_PartLevel[ig];
      //HJET SPECTRA
	   (*EvCollection).h1h[ig] = *fhRecoilJetPtTTH_PartLevel[ig];
	   (*EvCollection).h2[ig] = *fhRecoilJetPhiTTH_PartLevel[ig];
	   //Xsect & Ntrials spectra
	   (*EvCollection).hXsect[ig] = *hXsect_part[ig];
	   (*EvCollection).hNtrials[ig] = *hNtrials_part[ig];
   }
   for(Int_t ig=0; ig< kTRIG; ig++){
	   delete fhTTH_PartLevel[ig];
	   delete fhRecoilJetPtTTH_PartLevel[ig];
	   delete fhRecoilJetPhiTTH_PartLevel[ig];
	   delete hXsect_part[ig];
	   delete hNtrials_part[ig];
  }

   pythia.stat();
	
	for (auto& item: GenStorage) {
		delete item;
	}
	
   return *EvCollection;
}

//_________________________________________________________________________ 
//_________________________________________________________________________ 
//_________________________________________________________________________ 
//_________________________________________________________________________ 
bool EtaCut(fastjet::PseudoJet fjJet, double etaMin, double etaMax) {
   if(fjJet.eta() > etaMax || fjJet.eta() < etaMin){
      return false;
   }else{
      return true;
   }
}



int main (int argc, char* argv[]){
	
	std::ofstream ofs;
	ofs.open ("test.txt", std::ofstream::out | std::ofstream::trunc);
	ofs.close();
	
	ROOT::EnableThreadSafety();
	
	Int_t tune = -1;
	Int_t rndseed = -1;
	Float_t jr = -1;
	Int_t nthr = -1;
	TString name = "";
	//vector<mt19937> GeneratorsStorage;
	
   if(argc!=5){  
      cout<<"Usage:"<<endl<<"./pygen <PythiaTune> <Number> <jetR> <NumberOfThreads>"<<endl;
      return 0;
   }
	
	tune = atoi(argv[1]);
	rndseed = atoi(argv[2]);
	Int_t cislo = atoi(argv[2]);
	jr = atof(argv[3]);
	nthr = atoi(argv[4]);
	
	vector<Int_t> ttl = {6, 12, 20};     //TT bin
	vector<Int_t> tth = {7, 20, 30};    //TT bin 
	vector<Int_t> pthmin = {4, 11, 21, 36, 56, 84, 117, 156, 200, 249};  //hard bin
	vector<Int_t> pthmax = {11, 21, 36, 56, 84, 117, 156, 200, 249, 1000}; //hard bin
	
	//jr = 0.4;     //jet R
		
	vector<future<TH1Coll>> futs;
	
	string seedpath = "./Results/" + to_string(rndseed) + "/";
	char* nseed = argv[2];
	
	MakeDir("./", "Results");
	MakeDir("./Results/", to_string(rndseed));
	
	vector<TH1Coll> Col_fspace;
	for (size_t itt = 0; itt < ttl.size(); itt++){
		vector<TH1Coll> Collection;
		
		string tt_path = "TT-" + to_string(ttl.at(itt)) + "_" + to_string(tth.at(itt));
		MakeDir(seedpath, tt_path);
		
		for (size_t ipt = 0; ipt < pthmin.size(); ipt++){
			TH1Coll* hcol = new TH1Coll(ttl.at(itt), tth.at(itt), pthmin.at(ipt), pthmax.at(ipt), rndseed*51);
			
			ofs.open ("test.txt", std::ofstream::out | std::ofstream::app);
			ofs << tt_path << "; pT-" << pthmin.at(ipt) << "_" << pthmax.at(ipt) << endl;
			ofs.close();

			for (Int_t i = 0; i < nthr; i++){
				rndseed += i;
				futs.push_back(
					async([itt, ipt, i, tune, rndseed, ttl, tth, pthmin, pthmax, jr]{
						vector<string> vec = {"",
											   to_string(tune), 
										   to_string(rndseed), 
										   to_string(ttl.at(itt)),
										   to_string(tth.at(itt)),
										   to_string(pthmin.at(ipt)),
										   to_string(pthmax.at(ipt)),
										   to_string(jr)};
								return Run(8, vec);
							}
						)
					);
			}
			
			mutex M_sum_of_items; //creates mutual exclusion for threads during the summation of ready data
			for (auto& item : futs) {
					if (item.valid()){
						lock_guard<mutex> lock_sum(M_sum_of_items); //thread is locked to use only one set of data and avoid "data race"
						TH1Coll collitem = item.get();
						*hcol = *hcol + collitem;
					}
			}
			Collection.push_back(*hcol);
			
//OUTPUT

			TString tag = Form("PP6p35_ANTIKT%02d_TT%d_%d",
			   TMath::Nint(jr*10), 
			   ttl.at(itt), tth.at(itt));
			
			if(pthmin.at(ipt)<0 || pthmax.at(ipt) <0){
				name = Form("./Results/%s/%s_tune%d_MB_c%s.root", nseed, tag.Data(), tune, nseed);
			}else{
				name = Form("./Results/%s/TT-%s_%s/%s_HB%d_%d_tune%d_c%s.root", nseed, 
							to_string(ttl.at(itt)).c_str(), to_string(tth.at(itt)).c_str(), 
							tag.Data(), pthmin.at(ipt), pthmax.at(ipt), tune, nseed);
			}

			TFile* outFile = new TFile(name.Data(), "RECREATE");
			outFile->cd();
			
			for(Int_t ig=0; ig< kTRIG; ig++){
				name = Form("hTT_%s_PartLevel_seed%d_tt_%d_pt%d", trg[ig].Data(), cislo, ttl.at(itt), pthmin.at(ipt));
				(*hcol).h1t[ig].SetName(name.Data());
				(*hcol).h1t[ig].Write();
				
				name = Form("fhRecoilJetPt_%s_PartLevel_seed%d_tt_%d_pt%d", trg[ig].Data(), cislo, ttl.at(itt), pthmin.at(ipt));
				(*hcol).h1h[ig].SetName(name.Data());
				(*hcol).h1h[ig].Write();

				name = Form("fhRecoilJetPhi_%s_PartLevel_seed%d_tt_%d_pt%d", trg[ig].Data(),cislo, ttl.at(itt), pthmin.at(ipt));
				(*hcol).h2[ig].SetName(name.Data());
				(*hcol).h2[ig].Write();

				name = Form("fHistXsections_%s_%d_tt_%d_pt%d", trg[ig].Data(), cislo, ttl.at(itt), pthmin.at(ipt));
				(*hcol).hXsect[ig].SetName(name.Data());
				(*hcol).hXsect[ig].Write();
				
				name = Form("fHistTrials_%s_%d_tt_%d_pt%d", trg[ig].Data(), cislo, ttl.at(itt), pthmin.at(ipt));
				(*hcol).hNtrials[ig].SetName(name.Data());
				(*hcol).hNtrials[ig].Write();
				
			}
			 outFile->Close();
		}
		
	}
	return 0;
}
