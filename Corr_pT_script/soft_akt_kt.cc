// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

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
#include "TFile.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TProfile.h"
#include "TList.h"
#include "TVector3.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"  
#include <ctime>

#include <sys/stat.h>
#include "aux_func.h"

using namespace Pythia8;

//___________________________________________________________________

int soft(int argc, vector<string> argv) {
   //at the moment the macro expects 2 arguments   : random seed   and jet radius


   Int_t seed  = -1;     //dummy initialization of the unique initial random seed number for each file.    
   Int_t tune  = 14;     //pythia Monash tune. pythia tunes define procecess that are considered in collisions  
   Int_t nEvent= 1e5;    //the number of events which will be processed. if you need more change this number and recompile
   if(argc!=3){  
      cout<<"Usage:"<<endl<<"./pygen <Seed> <jetR>"<<endl;    //just testing whether the number of arguments is correct
      return 0;
   }
   seed   = stoi(argv[1]);    // initialization of the random seed from the first argument
   //__________________________________________________________________________
   //                        ANALYSIS SETTINGS

   Double_t jetParameterR   = (Double_t) stof(argv[2]); //initialization of jet cone radius from the second argument
   Double_t trackEtaCut     = 0.9; //psedorapidity  cut on accepted tracks
                                   //try look yourself what pseudorapidity means
				   //roughly speaking it is rapidity for massless particles
				   //at high  energies all particles can be considered massles (80% of all produced particles are pions)
				   //pseudorapidity is related to polar angle, ie. the angle w.r.t. beam axis
				   //pseudorapidity 0  ..... direction perpendicular to beam axis
                                   //pseudorapidity -> infinity   if you approach beam axis
				   //the |pseudorapidity| < 0.9  selects region of polar angles covered by the ALICE central barrel detectors
				   //you can try to calculate wich coverage this is in degrees 

   TString name;  //auxiliary string
   //__________________________________________________________________________
   //        CONFIGURATION OF  PYTHIA

   // Generator. Process selection. LHC initialization. Histogram.
   Pythia pythia;
   //COLLISION SYSTEM  what are you going to collide  = two protons with center of mass energy 13 TeV 
   pythia.readString("Beams:idA = 2212"); //beam 1 2212  is  proton PDG code  (you can try to google out what are PDG codes of other particles, e.g.  photon)
   pythia.readString("Beams:idB = 2212"); //beam 2 p

   //CENTER OF MASS ENERGY of collision 
   pythia.readString("Beams:eCM = 13000."); // in GeV

   //PYTHIA TUNE
	//pythia.readString(Form("Tune:ee = %d",7));
   pythia.readString(Form("Tune:pp = %d",tune));  //tune 1-14    5=defaulr TUNE4C,  6=Tune 4Cx, 7=ATLAS MB Tune A2-CTEQ6L1  14=Monash
	

   //SEED OF PSEUDORANDOM GENERATOR   ...  try to run pythia with the same initial seed, you will see that you will get identical results
   //              when running you should therefore make sure that you are changing the random seed
   //              if you will set the seed to 0  pythia will generate initial seed randomly by itself
   //
   pythia.readString("Random:setSeed = on");
   pythia.readString(Form("Random:seed = %d",seed));

   //QCD PROCESSES THAT WILL BE INCLUDED
   //pythia.readString("HardQCD:all = on");             //for hard bin configuration -ie you generate events with certain cutoff on energy that was exchanged 
   pythia.readString("SoftQCD:inelastic = on");     //for minimum bias configuration -ie you generate any event

   //SKIP PROCESSES WITH Q2 < 5 GeV/c  
   //pythia.readString("PhaseSpace:pTHatMin = 10."); // <<<<<<<<<<<<<<<<<<<<<<< this is the minimum cutoff.    the mimimum cutoff allowed by pythia is 3 
   //pythia.readString("PhaseSpace:pTHatMax = 50."); // <<<<<<<<<<<<<<<<<<<<<<< this is the maxium cutoff  
	
   //SWITCH OFF DECAYS TO SECONDARY PARTICLES which would decay via weak interaction
   pythia.readString("310:mayDecay  = off"); //K0s
   pythia.readString("3122:mayDecay = off"); //labda0
   pythia.readString("3112:mayDecay = off"); //sigma-
   //pythia.readString("3212:mayDecay = off"); //sigma0
   pythia.readString("3222:mayDecay = off"); //sigma+
   pythia.readString("3312:mayDecay = off"); //xi-
   pythia.readString("3322:mayDecay = off"); //xi0
   pythia.readString("3334:mayDecay = off"); //omega-


	
   pythia.init();
 
   //___________________________________________________ 
   //          FAST JET CONFIGURATION 


   fastjet::Strategy strategy = fastjet::Best;
   //PT RECOMBINATION SCHEME  ... defines a way how the fourmomenta of particles will be combined to jets
   fastjet::RecombinationScheme recombScheme = fastjet::BIpt_scheme;

   //Define clustering algorithms
	fastjet::JetDefinition *jetDefAKTCh = NULL; 
	fastjet::JetDefinition *jetDefKTCh = NULL;

	jetDefAKTCh = new fastjet::JetDefinition(fastjet::antikt_algorithm,
                                                 jetParameterR,  
                                                 recombScheme, 
                                                 strategy);
	jetDefKTCh = new fastjet::JetDefinition(fastjet::kt_algorithm,
                                                 jetParameterR,  
                                                 recombScheme, 
                                                 strategy);

   //PARAMETERS THAT GO INTO THE MEASUREMENT OF ACTIVE JET AREAS. 
   //fast jet is going to inject to event  many infinitesimally soft particles (energy 10^-99)
   //these particles do not influence your jets (infrared and colinear safety of antikt)
   //by counting how many ghost particles were found in the jet
   //he will measure the jet area
   fastjet::GhostedAreaSpec ghostareaspec(trackEtaCut, 1, 0.005);       
   fastjet::AreaType areaType = fastjet::active_area_explicit_ghosts;
   fastjet::AreaDefinition *areaDef = new fastjet::AreaDefinition(areaType, ghostareaspec);

   // FASTJET INPUT  : the algorithmus needs as an input list of particle fourmomenta verctors
   // in steps he merges these fourvecotors together and gets jets
	std::vector<fastjet::PseudoJet> fjInputsAKT;
	std::vector<fastjet::PseudoJet> fjInputsKT;

  
   //___________________________________________________ 
   //                HISTOGRAMS

	//Histograms for the calculation of Xsection
	
	TString Xsstring = "fHistXsection";
	TProfile* fHistXsection = new TProfile(Xsstring, "fHistXsection", 1, 0, 1);
	fHistXsection->GetYaxis()->SetTitle("xsection");
	
	TString HTrials = "fHistTrials";
	TH1F* fHistTrials = new TH1F(HTrials, "fHistTrials", 1, 0, 1);
	fHistTrials->GetYaxis()->SetTitle("trials");
	
	//Histograms for JETS
   
	//pT distribution of jets, kT clustering
	TH1D* fJetPt_KT;
	name = "soft_hJetPt_kT"; 
	fJetPt_KT = new TH1D(name.Data(), "Jet pT distribution for kT clustering sequence", 100, 0, 50.0);
	fJetPt_KT->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
	fJetPt_KT->GetYaxis()->SetTitle("counts");
	fJetPt_KT->Sumw2();
	
	//pT distribution of jets, anti-kT clustering
	TH1D* fJetPt_AKT;
	name = "soft_hJetPt_anti-kT"; 
	fJetPt_AKT = new TH1D(name.Data(), "Jet pT distribution for anti-kT clustering sequence", 100, 0, 50.0);
	fJetPt_AKT->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
	fJetPt_AKT->GetYaxis()->SetTitle("counts");
	fJetPt_AKT->Sumw2();
	
	//pT distribution of jets, corrected values
	TH1D* fJetPt_SUB;
	name = "soft_hJetPt"; 
	fJetPt_SUB = new TH1D(name.Data(), "Jet pT distribution with the subtraction of underlying events", 140, -20.0, 50.0);
	fJetPt_SUB->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
	fJetPt_SUB->GetYaxis()->SetTitle("counts");
	fJetPt_SUB->Sumw2();
	



   //___________________________________________________ 
   // BEGIN EVENT LOOP. Generate event. Skip if error.
   for(int iEvent = 0; iEvent < nEvent; iEvent++){
      if(!pythia.next()) continue;

	   //reset the list of particles
		fjInputsAKT.resize(0);   
		fjInputsKT.resize(0);

      //______________________________
      // LOOP OVER CHARGED TRACKS 

      for(Int_t i = 0; i < pythia.event.size(); i++){
         if(!pythia.event[i].isFinal()) continue;  //select FINAL state PARTICLES ONLY  (after electromagnetic and strong decays)
	                                          // if you would comment this line out you would see that
						  //  pythia would give you also  quarks, color strings, resonances
						  //  pythia event contains all history of the scattering process
						  //  starting from the initial scatterd quarks till the final state particles
						  //  in detector you measure final state particles only 
         if(!pythia.event[i].isCharged()) continue; //select CHARGED PARTICLES   - neutral particles do not make tracks in detector
          //TASK: here should be also some cut on pseudorapidity of ALICE central barrel, try to implement is
		 //      plot pseudorapity and azimuthal distribution of tracks with pT > 150 MeV
		 if(pythia.event[i].pT() > 0.15 && TMath::Abs(pythia.event[i].eta()) < trackEtaCut) {
			 
			 fjInputsAKT.push_back( fastjet::PseudoJet(pythia.event[i].px(),         //make a list of charged particle 4-vectors from which he makes a jet
                                                   pythia.event[i].py(),
                                                   pythia.event[i].pz(),
                                                   pythia.event[i].pAbs()));
			 fjInputsKT.push_back( fastjet::PseudoJet(pythia.event[i].px(),         //make a list of charged particle 4-vectors from which he makes a jet
									   pythia.event[i].py(),
									   pythia.event[i].pz(),
									   pythia.event[i].pAbs()));
		 }				
      }
      
/*      //__________________________________________________
      //BUILD a CHARGED JET
      for(Int_t i = 0; i < pythia.event.size(); i++){
         if(pythia.event[i].isFinal() && pythia.event[i].isCharged()&& pythia.event[i].pT() > 10.0 && TMath::Abs(pythia.event[i].eta())< 0.9){
            fjInputs.push_back( fastjet::PseudoJet(pythia.event[i].px(),         //make a list of charged particle 4-vectors from which he makes a jet
                                                   pythia.event[i].py(),
                                                   pythia.event[i].pz(),
                                                   pythia.event[i].pAbs()));

            //(fjInputs.at(fjInputs.size()-1)).set_user_index(i); 

            cout<<"TRK PT="<<pythia.event[i].pT()<<
                  "ETA ="<< pythia.event[i].eta()<<
                  "PHI ="<< pythia.event[i].phi()<<"  "<<i<<endl;

         }
      }
*/
      //--------------------------------------------------------
		vector<fastjet::PseudoJet> inclusiveAKTJetsCh;
		fastjet::ClusterSequenceArea clustSeq_Sig_AKT(fjInputsAKT, *jetDefAKTCh, *areaDef);
		vector<fastjet::PseudoJet> inclusiveKTJetsCh;
		fastjet::ClusterSequenceArea clustSeq_Sig_KT(fjInputsKT, *jetDefKTCh, *areaDef);

		inclusiveAKTJetsCh = clustSeq_Sig_AKT.inclusive_jets(0.15); 
		inclusiveKTJetsCh = clustSeq_Sig_KT.inclusive_jets(0.); 

		  //_______________________________________________________________________________________

		  //_______________________________________________________________________________________
		  //LOOP OVER jets and print their CONSTITUENTS
		double rhosum = 0.;
		for(unsigned int ijet = 0; ijet < inclusiveKTJetsCh.size(); ijet++){ //loop over all full jets, kT
			fastjet::PseudoJet fjJet_KT = inclusiveKTJetsCh.at(ijet);		  
			if(!(TMath::Abs(fjJet_KT.eta()) < (trackEtaCut))) continue; 
			rhosum += (fjJet_KT.pt() / fjJet_KT.area());
			fJetPt_KT->Fill((Double_t) fjJet_KT.pt()); 
		}
		double rho = rhosum / inclusiveKTJetsCh.size();

		for(unsigned int aktjet = 0; aktjet < inclusiveAKTJetsCh.size(); aktjet++){ //loop over all full jets, anti-kT + subtraction
			fastjet::PseudoJet fjJet_AKT = inclusiveAKTJetsCh.at(aktjet);		  
			if(!(TMath::Abs(fjJet_AKT.eta()) < (trackEtaCut))) continue; 
			fJetPt_AKT->Fill((Double_t) fjJet_AKT.pt()); 
			//subtract underlying events
			Double_t corr_pt = (Double_t) fjJet_AKT.pt() - rho * (Double_t) fjJet_AKT.area();
			fJetPt_SUB->Fill(corr_pt); 
		}

   }// End of event loop.

	//Fill histograms to obtain cross sections
	fHistXsection->Fill(0.5,pythia.info.sigmaGen());
	fHistTrials->Fill(0.5, pythia.info.nAccepted());
	Double_t xsection = fHistXsection->GetMean(2);
	Double_t ntrials =  fHistTrials->Integral();
	Double_t weight = xsection/ntrials;
	
	fJetPt_SUB->Scale(weight);
	
	
   //____________________________________________________
   //          SAVE OUTPUT
	MakeDir("./", "Results");
	
	TString tag = Form("hard_PP7_pth_akT_kT%02d", TMath::Nint(jetParameterR*10) );

	TFile* outFile = new TFile(Form("./Results/%s_tune%d_c%d.root",tag.Data(), tune, seed), "UPDATE");
	outFile->cd();

	fJetPt_SUB->Write();

   	outFile->Close();


	pythia.stat();
   return 0;
}



