#ifndef myexampleTOOLS_H
#define myexampleTOOLS_H 

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  

#include "Pythia8/Pythia.h"

#include "myFastJetBase.h"

using namespace std;

class myTools {
    private:
        int m_test;



    public:
        myTools();
        
        // methods
	double width(fastjet::PseudoJet jet);
	int ntrack(fastjet::PseudoJet jet,double pT);
        double JetCharge(fastjet::PseudoJet jet,double kappa);
	double JetTrackMass(fastjet::PseudoJet jet,int which);
	bool IsBHadron(int pdgId);
	bool IsCHadron(int pdgId);
	bool Btag(fastjet::PseudoJet jet,vector<fastjet::PseudoJet> bhadrons,vector<fastjet::PseudoJet> chadrons,double jetrad,double b, double c, double uds);
	bool BosonMatch(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> Bosons, double jetrad, int BosonID);
	bool IsIsolated(Pythia8::Particle* particle, Pythia8::Pythia* pythia8, float rel_iso, float conesize);
	fastjet::PseudoJet Add(fastjet::PseudoJet jet);
	fastjet::PseudoJet Drop(fastjet::PseudoJet jet);
	fastjet::PseudoJet Angles(fastjet::PseudoJet jet);
	fastjet::PseudoJet Scale(fastjet::PseudoJet jet);

};

#endif

