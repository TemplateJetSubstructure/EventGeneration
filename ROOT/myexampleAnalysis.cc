#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>


#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "myexampleAnalysis.h"
#include "myTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "Pythia8/Pythia.h"

using namespace std;

// Constructor 
myexampleAnalysis::myexampleAnalysis(){
    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    fOutName = "test.root";
    tool = new myTools();

    // jet def 
    m_jet_def               = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.8);
    m_jet_def_ca  = new fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.8);

    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis End " << endl;
}

// Destructor 
myexampleAnalysis::~myexampleAnalysis(){
    delete tool;
    delete m_jet_def;
    delete m_jet_def_ca;
}

// Begin method
void myexampleAnalysis::Begin(){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("EventTree", "Event Tree for myexample");
    
   DeclareBranches();
   ResetBranches();
   

   return;
}

// End
void myexampleAnalysis::End(){
    
    tT->Write();
    tF->Close();
    return;
}

// Analyze
void myexampleAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8){
    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if (!pythia8->next()) return;
    //pythia8->event.list();

    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches 
    ResetBranches();
    
    // new event-----------------------
    fTEventNumber = ievt;
    std::vector <fastjet::PseudoJet>           particlesForJets;
    std::vector <fastjet::PseudoJet>           bhadrons;
    std::vector <fastjet::PseudoJet>           chadrons;
    std::vector <fastjet::PseudoJet>           HS;
    
    // Particle loop -----------------------------------------------------------
    for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){

      fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e());
      p.set_user_info(new MyUserInfo(pythia8->event[ip].id(),ip,pythia8->event[ip].charge()));
      
      if (tool->IsBHadron(pythia8->event[ip].id())){
        bhadrons.push_back(p);
      }
      if (fabs(pythia8->event[ip].id()) < 40){
        HS.push_back(p);
      }

      // particles for jets --------------
      if (!pythia8->event[ip].isFinal() )      continue;
      if (fabs(pythia8->event[ip].id())  ==12) continue;
      if (fabs(pythia8->event[ip].id())  ==13) continue;
      if (fabs(pythia8->event[ip].id())  ==14) continue;
      if (fabs(pythia8->event[ip].id())  ==16) continue;

      particlesForJets.push_back(p);
      
    } // end particle loop -----------------------------------------------

    fastjet::ClusterSequence cs(particlesForJets, *m_jet_def);
    vector<fastjet::PseudoJet> myJets = fastjet::sorted_by_pt(cs.inclusive_jets(25.0));
    fastjet::contrib::SoftDrop sd(0.0, 0.1); //beta = 0, zcut = 0.1
    fastjet::contrib::Nsubjettiness nSub1KT(1, fastjet::contrib::Njettiness::wta_kt_axes, 1., 0.8, 0.8);
    fastjet::contrib::Nsubjettiness nSub2KT(2, fastjet::contrib::Njettiness::wta_kt_axes, 1., 0.8, 0.8);

    for (unsigned int ij = 0; ij < myJets.size(); ij++){   
      if(fTNJetsSmallRFilled == MaxNJetSmallR) {cout << "Warning: More than " << MaxNJetSmallR << " small R jets" << endl; continue;}      
      if (myJets[ij].pt() > 500){
	fastjet::ClusterSequence cs2(myJets[ij].constituents(), *m_jet_def_ca);
	vector<fastjet::PseudoJet> myJets2 = sorted_by_pt(cs2.inclusive_jets(0.));
	fastjet::PseudoJet sdJet = sd(myJets2[0]);
	double mg = sdJet.m();
	double tau1 = nSub1KT(myJets[ij]);
	double tau2 = nSub2KT(myJets[ij]);
	double tau21 = tau1 > 0 ? tau2/tau1 : 0.;
	fTJsmallPt       [fTNJetsSmallRFilled] = myJets[ij].pt();
        fTJsmallEta      [fTNJetsSmallRFilled] = myJets[ij].eta();
        fTJsmallPhi      [fTNJetsSmallRFilled] = myJets[ij].phi();
        fTJsmallM        [fTNJetsSmallRFilled] = myJets[ij].m();     
	fTJsmallMg       [fTNJetsSmallRFilled] = mg;
	fTJsmalltau21    [fTNJetsSmallRFilled] = tau21;
	double Emax = -1;
        int flav = 0;
        for (unsigned int i=0; i < HS.size(); i++){
          if (HS[i].delta_R(myJets[ij]) > 0.4) continue;
          if (HS[i].e() > Emax){
            flav = HS[i].user_info<MyUserInfo>().pdg_id();
            Emax = HS[i].e();
          }
	}
	fTJsmalltype     [fTNJetsSmallRFilled] = flav;
	fTNJetsSmallRFilled++;
      }
    }
    tT->Fill();
      
    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent End " << endl;
    return;
}



// declate branches
void myexampleAnalysis::DeclareBranches(){
   
   // Event Properties 
   tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
   
   // smallR jets
   tT->Branch("NJetsFilledSmallR",         &fTNJetsSmallRFilled,       "NJetsFilledSmallR/I");
   tT->Branch("JsmallPt",                  &fTJsmallPt,                "JsmallPt[NJetsFilledSmallR]/F");
   tT->Branch("JsmallEta",                 &fTJsmallEta,               "JsmallEta[NJetsFilledSmallR]/F");
   tT->Branch("JsmallPhi",                 &fTJsmallPhi,               "JsmallPhi[NJetsFilledSmallR]/F");
   tT->Branch("JsmallM",                   &fTJsmallM,                 "JsmallM[NJetsFilledSmallR]/F");
   tT->Branch("JsmallMg",                   &fTJsmallMg,                 "JsmallMg[NJetsFilledSmallR]/F");
   tT->Branch("Jsmalltau21",                   &fTJsmalltau21,                 "Jsmalltau21[NJetsFilledSmallR]/F");
   tT->Branch("Jsmalltype",                   &fTJsmalltype,                 "Jsmalltype[NJetsFilledSmallR]/I");
   tT->GetListOfBranches()->ls();
    
   return;
}


// resets vars
void myexampleAnalysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;

      fTNJetsSmallRFilled=0;
      for (int iP=0; iP < MaxNJetSmallR; ++iP){
          fTJsmallPt      [iP]= -999;
          fTJsmallPhi     [iP]= -999;
          fTJsmallEta     [iP]= -999;
          fTJsmallM       [iP]= -999;
	  fTJsmallMg      [iP]= -999;
	  fTJsmalltau21   [iP]= -999;
	  fTJsmalltype    [iP]= -999;
      }
}
