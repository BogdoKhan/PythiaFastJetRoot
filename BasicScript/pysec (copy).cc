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



using namespace Pythia8;


//___________________________________________________________________

int main(int argc, char* argv[]) {
   //at the moment the macro expects 2 arguments   : random seed   and jet radius

   Int_t seed  = -1;     //dummy initialization of the unique initial random seed number for each file.    
   Int_t tune  = 14;     //pythia Monash tune. pythia tunes define procecess that are considered in collisions  
   Int_t nEvent= 1e1;    //the number of events which will be processed. if you need more change this number and recompile
   if(argc!=3){  
      cout<<"Usage:"<<endl<<"./pygen <Seed> <jetR>"<<endl;    //just testing whether the number of arguments is correct
      return 0;
   }
   seed   = atoi(argv[1]);    // initialization of the random seed from the first argument
   //__________________________________________________________________________
   //                        ANALYSIS SETTINGS

   Double_t jetParameterR   = (Double_t) atof(argv[2]); //initialization of jet cone radius from the second argument
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
   pythia.readString("Beams:idA = 2212"); //beam 1 proton    2212  is  proton PDG code  (you can try to google out what are PDG codes of other particles, e.g.  photon)
   pythia.readString("Beams:idB = 2212"); //beam 2 proton

   //CENTER OF MASS ENERGY of collision 
   pythia.readString("Beams:eCM = 13000."); // in GeV

   //PYTHIA TUNE
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
   //pythia.readString("PhaseSpace:pTHatMin = 20."); // <<<<<<<<<<<<<<<<<<<<<<< this is the minimum cutoff.    the mimimum cutoff allowed by pythia is 3 
   // pythia.readString("PhaseSpace:pTHatMax = 5."); // <<<<<<<<<<<<<<<<<<<<<<< this is the maxium cutoff  at the moment commented out

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

   //DEFINE ANTIKT ALGORITHM  
   fastjet::JetDefinition *jetDefAKTCh = NULL; 

   jetDefAKTCh = new fastjet::JetDefinition(fastjet::antikt_algorithm,
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
   std::vector<fastjet::PseudoJet> fjInputs;

  
   //___________________________________________________ 
   //                HISTOGRAMS

   
   //PT DISTRIBUTION OF TRACKS
   TH1D* fTrackPt;
   name = "hTrackPt";
   fTrackPt = new TH1D(name.Data(),"Track pT distribution", 50, 0.0, 50.0);
   fTrackPt->GetXaxis()->SetTitle("track p_{T} [GeV/c]");
   fTrackPt->GetYaxis()->SetTitle("counts");
   fTrackPt->Sumw2(); //tells to root to properly handle statistical errors
   
   //PT DISTRIBUTIONS OF JETS
   TH1D* fJetPt;
   name = "hJetPt"; 
   fJetPt = new TH1D(name.Data(), "Jet pT distribution", 50, 0, 50.0);
   fJetPt->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
   fJetPt->GetYaxis()->SetTitle("counts");
   fJetPt->Sumw2();




   //___________________________________________________ 
   // BEGIN EVENT LOOP. Generate event. Skip if error.
   for(int iEvent = 0; iEvent < nEvent; iEvent++){
      if(!pythia.next()) continue;

      //here you start with event 
      fjInputs.resize(0);  //reset the list of particles 

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
         fTrackPt->Fill(pythia.event[i].pT()); //FILL HISTOGRAM WITH TRACK PT
      }
      
      //__________________________________________________
      //BUILD a CHARGED JET
      for(Int_t i = 0; i < pythia.event.size(); i++){
         if(pythia.event[i].isFinal() && pythia.event[i].isCharged()&& TMath::Abs(pythia.event[i].eta())<0.9){
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

      //--------------------------------------------------------
      vector<fastjet::PseudoJet> inclusiveJetsCh;
      fastjet::ClusterSequenceArea clustSeq_Sig(fjInputs, *jetDefAKTCh, *areaDef);

      inclusiveJetsCh = clustSeq_Sig.inclusive_jets(0.150); //the lowerst accepted  jet should have pT > 150 MeV
                                                            // the lowest pT that can be measured by ALICE
      //_______________________________________________________________________________________
      //LOOP OVER JETS
      for(unsigned int ijet = 0; ijet < inclusiveJetsCh.size(); ijet++){
         fastjet::PseudoJet fjJet = inclusiveJetsCh.at(ijet);
    
         fJetPt->Fill((Double_t) fjJet.pt()); //FILL HISTOGRAM WITH JET PT

      }//jet loop

      //_______________________________________________________________________________________
      //LOOP OVER jets and print their CONSTITUENTS
      for(unsigned int ijet = 0; ijet < inclusiveJetsCh.size(); ijet++){ //loop over all full jets
          cout<<"JET ................ "<<ijet<<endl;
          fastjet::PseudoJet fjJet = inclusiveJetsCh.at(ijet);
          if(fjJet.pt()<0.150) continue; 
            vector<fastjet::PseudoJet> constituents =  clustSeq_Sig.constituents(inclusiveJetsCh[ijet]); //for jet get list of constituents
            for(unsigned int ict = 0; ict < constituents.size(); ict++){
               if(constituents.at(ict).pt() < 0.15) continue;
               cout<<"JETC pt="<<constituents.at(ict).pt() 
                   <<" eta ="<<constituents.at(ict).eta() 
                   <<" phi ="<<constituents.at(ict).phi()<<endl;
                   //<<" lab ="<<constituents.at(ict).user_index()<<endl;
            }
     }
cout<<"=============================================="<<endl;

   }// End of event loop.

   //____________________________________________________
   //          SAVE OUTPUT

   TString tag = Form("PP7_pythia_ANTIKT%02d", TMath::Nint(jetParameterR*10) );

   TFile* outFile = new TFile(Form("%s_tune%d_c%d.root",tag.Data(), tune, seed), "RECREATE");
   outFile->cd();

   fTrackPt->Write();
   fJetPt->Write();

   outFile->Close();


   pythia.stat();
   return 0;
}



