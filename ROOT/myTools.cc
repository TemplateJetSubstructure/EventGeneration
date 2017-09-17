#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  

#include "myTools.h"
#include "myFastJetBase.h"

#include "TRandom3.h"
#include "TError.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Constructor 
myTools::myTools(){
    m_test = 0;
}

fastjet::PseudoJet myTools::Add(fastjet::PseudoJet jet){
  fastjet::PseudoJet out(0.,0.,0.,0.);
  TRandom3 *rand = new TRandom3(0);
  double myang = rand->Uniform(0,2*3.141);
  double myrad = rand->Uniform(0,1);
  myrad = sqrt(myrad);
  double eta = jet.eta()+0.4*myrad*cos(myang);
  double phi = jet.phi()+0.4*myrad*sin(myang);
  fastjet::PseudoJet out2(0.,0.,0.,0.);
  out2.reset_momentum_PtYPhiM(0.5,eta,phi,0.);
  out = jet;
  out+=out2;
  return out;
}

fastjet::PseudoJet myTools::Angles(fastjet::PseudoJet jet){
  fastjet::PseudoJet out(0.,0.,0.,0.);
  TRandom3 *rand = new TRandom3(0);
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    double phi_smear = rand->Gaus(jet.constituents()[i].phi(),0.005);
    double eta_smear = rand->Gaus(jet.constituents()[i].eta(),0.005);
    fastjet::PseudoJet out2;
    //out2.reset_momentum_PtYPhiM(jet.constituents()[i].e()/cosh(eta_smear),eta_smear,phi_smear,0.);
    out2.reset_momentum_PtYPhiM(jet.constituents()[i].pt(),eta_smear,phi_smear,0.);
    out+=out2;
  }
  return out;
}

fastjet::PseudoJet myTools::Drop(fastjet::PseudoJet jet){
  fastjet::PseudoJet out(0.,0.,0.,0.);
  TRandom3 *rand = new TRandom3(0);
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    double r = 0.25*exp(-2*jet.constituents()[i].e());
    double flip = rand->Uniform(0,1);
    if (( r < flip) && (jet.constituents()[i].e() <2.5))
      continue;
    out+=jet.constituents()[i];
  }
  return out;
}

fastjet::PseudoJet myTools::Scale(fastjet::PseudoJet jet){
  fastjet::PseudoJet out(0.,0.,0.,0.);
  TRandom3 *rand = new TRandom3(0);
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    double fact = 0.05*(1+1.5/jet.pt());
    double fact_smear = rand->Gaus(1,fact);
    out+=jet.constituents()[i]*fact_smear;
  }
  return out;
}

double myTools::width(fastjet::PseudoJet jet){
  double out=0.;
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    out+= (jet.constituents()[i].pt()/jet.pt())*pow(jet.constituents()[i].delta_R(jet),2);
  }
  return sqrt(out);
}

int myTools::ntrack(fastjet::PseudoJet jet,double pT){
  int out=0.;
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    if (fabs(jet.constituents()[i].user_info<MyUserInfo>().charge()) > 0.1 && jet.constituents()[i].pt() > pT){
      out++;
    }
  }
  return out;
}

double myTools::JetTrackMass(fastjet::PseudoJet jet, int which){
  TLorentzVector trackjet = TLorentzVector(0.,0.,0.,0.);
  int ntracks=0;
  TRandom3 *rand = new TRandom3(0);
  TLorentzVector myjet = TLorentzVector(jet.px(),jet.py(),jet.pz(),jet.e());
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    if (fabs(jet.constituents()[i].user_info<MyUserInfo>().charge()) > 0.1){
      //if (jet.constituents()[i].pt() < 0.5) continue;
      TLorentzVector hold = TLorentzVector(jet.constituents()[i].px(),jet.constituents()[i].py(),jet.constituents()[i].pz(),jet.constituents()[i].e());
      
      double phi = hold.Phi();
      double eta = hold.Eta();
      double pt = hold.Pt();

      //if (fabs(eta) > 2.5){
      //continue;
      //}
      //if (pt < 0.5) continue;

      if (which < 0){
	
	//Efficiency
	
	double drop = 0;
	if (pt < 1){
	  if (fabs(eta) <0.9){
	    drop = 0.9;
	  }
	  else if(fabs(eta) < 1.4){
	    drop = 0.9;
	  }
	  else if(fabs(eta) < 2.5){
	    drop = 0.85;
	  }
	}
	else if (pt < 30){
	  if (fabs(eta) <0.9){
	    drop = 0.98;
	  }
	  else if(fabs(eta) < 1.4){
	    drop = 0.93;
	  }
	  else if(fabs(eta) < 2.5){
	    drop = 0.9;
	  }
	}
	else {
	  if (fabs(eta) <0.9){
	    drop = 0.95;
	  }
	  else if(fabs(eta) < 1.4){
	    drop = 0.85;
	  }
	  else if(fabs(eta) < 2.5){
	    drop = 0.8;
	  }
	}
	
	double pi = 4*atan(1.);
	drop = 2.*drop*atan(myjet.DeltaR(hold)/0.002)/pi;
	
	if (rand->Uniform(0,1) > drop){
	  continue;
	}
	
	//Resolutions 
	
	double phismear = 0.003;
	if (pt < 1){
	  if (fabs(eta) <0.9){
	    phismear = 0.003;
	  }
	  else if(fabs(eta) < 1.4){
	    phismear = 0.004;
	  }
	  else if(fabs(eta) < 2.5){
	    phismear = 0.007;
	  }
	}
	else if (pt < 5){
	  if (fabs(eta) < 0.9){
	    phismear = 0.002;
	  }
	  else if (fabs(eta) < 1.4){
	    phismear = 0.003;
	  }
	  else if (fabs(eta) < 2.5){
	    phismear = 0.004;
	  }	  
	}
	else if (pt < 20){
	  if (fabs(eta) < 0.9){
	    phismear = 0.0005;
	  }
	  else if (fabs(eta) < 1.4){
	    phismear = 0.0006;
	  }
	  else if (fabs(eta) < 2.5){
	    phismear = 0.0008;
	  }
	}
	else {
	  if (fabs(eta) < 0.9){
	    phismear = 0.0002;
	  }
	  else if (fabs(eta) < 1.4){
	    phismear = 0.0002;
	  }
	  else if (fabs(eta) < 2.5){
	    phismear = 0.0004;
	  }
	}
	
	double pTsmear = 0.02;
	
	if (pt < 1){
	  if (fabs(eta) < 0.9){
	    pTsmear = 0.02;
	  }
	  else if (fabs(eta) < 1.4){
	    pTsmear = 0.02;
	  }
	  else if (fabs(eta) < 2.5){
	    pTsmear = 0.03;
	  }
	} 
	else if (pt < 10){
	  if (fabs(eta) < 0.9){
	    pTsmear = 0.01;
	  }
	  else if (fabs(eta) < 1.4){
	    pTsmear = 0.02;
	  }
	  else if (fabs(eta) < 2.5){
	    pTsmear = 0.03;
	  }
	} 
	else if (pt < 30){
	  if (fabs(eta) < 0.9){
	    pTsmear = 0.015;
	  }
	  else if (fabs(eta) < 1.4){
	    pTsmear = 0.025;
	  }
	  else if (fabs(eta) < 2.5){
	    pTsmear = 0.04;
	  }
	} 
	else if (pt < 50){
	  if (fabs(eta) < 0.9){
	    pTsmear = 0.02;
	  }
	  else if (fabs(eta) < 1.4){
	    pTsmear = 0.04;
	  }
	  else if (fabs(eta) < 2.5){
	    pTsmear = 0.08;
	  }
	} 
	else if (pt < 100){
	  if (fabs(eta) < 0.9){
	    pTsmear = 0.04;
	  }
	  else if (fabs(eta) < 1.4){
	    pTsmear = 0.07;
	  }
	  else if (fabs(eta) < 2.5){
	    pTsmear = 0.2;
	  }
	} 
	else {
	  if (fabs(eta) < 0.9){
	    pTsmear = 0.08;
	  }
	  else if (fabs(eta) < 1.4){
	    pTsmear = 0.14;
	  }
	  else if (fabs(eta) < 2.5){
	    pTsmear = 0.4;
	  }
	} 
	
	phi = rand->Gaus(phi,phismear);
	eta = rand->Gaus(eta,phismear);
	pt = rand->Gaus(pt,pTsmear*pt);
	hold.SetPtEtaPhiM(pt,eta,phi,hold.M());
      }	
      if (fabs(eta) > 2.5){
        continue;
      }
      if (pt < 0.5) continue;
      trackjet+=hold;
      ntracks+=1;
    }
  }
  delete rand;
  if (abs(which)==1) return trackjet.M();
  if (abs(which)==2) return trackjet.Pt();
  if (abs(which)==3) return ntracks;
}

double myTools::JetCharge(fastjet::PseudoJet jet,double kappa){
  //Returns the jet charge with weighting factor kappa
  double charge=0.;
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    charge+=jet.constituents()[i].user_info<MyUserInfo>().charge()*pow(jet.constituents()[i].pt(),kappa);
  }
  return charge/pow(jet.pt(),kappa);
}

bool myTools::IsBHadron(int pdgId){
  int abs_pdgId = abs(pdgId);
  int abs_pdgId_mod10k = (abs(pdgId)%10000);
  if( (abs_pdgId_mod10k>=500 && abs_pdgId_mod10k<600) /*mesons*/  ||
      (abs_pdgId>=5000      && abs_pdgId<6000)      /*baryons*/   )
    return true;

  return false;
}

bool myTools::IsCHadron(int pdgId){
  int abs_pdgId = abs(pdgId);
  int abs_pdgId_mod10k = (abs(pdgId)%10000);
  if( (abs_pdgId_mod10k>=400 && abs_pdgId_mod10k<500) /*mesons*/  ||
      (abs_pdgId>=4000      && abs_pdgId<5000)      /*baryons*/   )
    return true;

  return false;
}

bool myTools::Btag(fastjet::PseudoJet jet,vector<fastjet::PseudoJet> bhadrons,vector<fastjet::PseudoJet> chadrons,double jetrad,double b,double c,double uds){

  //TRandom3 *rand = new TRandom3(0);

  int foundb=0;
  int foundc=0;
  
  for (unsigned int i=0; i<bhadrons.size(); i++){
    if (bhadrons[i].pt() < 5) continue;
    if (bhadrons[i].delta_R(jet)<jetrad){
      foundb=1;
      return true; 
    }
  }

  /*
  for (unsigned int i=0; i<chadrons.size(); i++){
    if (chadrons[i].delta_R(jet)<jetrad){
      if (foundb!=1) foundc=1;
    }
  }

  if (foundb==1){
    double flip = rand->Uniform(0.,1.);
    if (flip < b){
      delete rand;
      return true;
    }
  }
  if (foundc==1){
    double flip= rand->Uniform(0.,1.);
    if (flip < 1./c){
      delete rand;
      return true;
    }
  }
  double flip= rand->Uniform(0.,1.);
  if (flip < 1./uds){
    delete rand;
    return true;
  }
  
  delete rand;
  */
  return false;
}

bool myTools::BosonMatch(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> Bosons, double jetrad, int BosonID ){

  for (unsigned int i=0; i<Bosons.size(); i++){
      if (Bosons[i].user_info<MyUserInfo>().pdg_id() != BosonID) continue;
      if (Bosons[i].delta_R(jet)<jetrad){
        return true;
      }
  }
  return false;
}

bool myTools::IsIsolated(Pythia8::Particle* particle, Pythia8::Pythia* pythia8, float rel_iso, float ConeSize){
    float sumpT=0;
    fastjet::PseudoJet part(particle->px(), particle->py(), particle->pz(),particle->e() );
    for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){
        if (!pythia8->event[ip].isFinal() )      continue;
        if (fabs(pythia8->event[ip].id())  ==12) continue;
        if (fabs(pythia8->event[ip].id())  ==14) continue;
        if (fabs(pythia8->event[ip].id())  ==16) continue;
        if (pythia8->event[ip].pT()       < 0.5) continue;
        if (&pythia8->event[ip] == particle    ) continue; //same particle
        fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e() );
        if(p.delta_R(part)>ConeSize)             continue;
        sumpT+=p.pt();
    }
    if(sumpT/part.pt()>rel_iso) return false;
    else return true;
}
