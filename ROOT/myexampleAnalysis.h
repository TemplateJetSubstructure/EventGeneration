#ifndef  myexampleAnalysis_H
#define  myexampleAnalysis_H

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"

#include "myTools.h"
#include "myFastJetBase.h"
#include "Pythia8/Pythia.h"

using namespace std;

class myexampleAnalysis{
    private:
        int  ftest;
        int  fDebug;
        string fOutName;

        TFile *tF;
        TTree *tT;
        myTools *tool;
	
	float fTMETx;
	float fTMETy;
	float fTLepx;
	float fTLepy;
	float fTLepz;

	float fTlep1x;
	float fTlep1y;
	float fTlep1z;
	int fTlep1id;
	float fTEventMass;

	float fTlep2x;
	float fTlep2y;
	float fTlep2z;
	int fTlep12d;

	float fTjetm;
	float fTjetpt;

        // Tree Vars ---------------------------------------
        int              fTEventNumber;
        int              fTPassBoost;
        int              fTPassResolve;
	int fTwhichfat;
        // Els
        static const int MaxNEles    = 5;
        int              fTNElesFilled;
        float            fTElesPt     [MaxNEles];
        float            fTElesEta    [MaxNEles];
        float            fTElesPhi    [MaxNEles];
        float            fTElesE      [MaxNEles];
        int              fTElesCharge [MaxNEles];
        // Muons
        static const int MaxNMuons    = 5;
        int              fTNMuonsFilled;
        float            fTMuonsPt     [MaxNMuons];
        float            fTMuonsEta    [MaxNMuons];
        float            fTMuonsPhi    [MaxNMuons];
        float            fTMuonsE      [MaxNMuons];
        int              fTMuonsCharge [MaxNMuons];
        // Bosons 
        static const int MaxNBosons    = 5;
        int              fTNBosonsFilled;
        float            fTBosonPt [MaxNBosons];
        float            fTBosonEta[MaxNBosons];
        float            fTBosonPhi[MaxNBosons];
        float            fTBosonM  [MaxNBosons];
        int              fTBosonID [MaxNBosons];

        static const int MaxNJetSmallR = 20;
        int              fTNJetsSmallRFilled;
        float            fTJsmallPt        [MaxNJetSmallR];
        float            fTJsmallEta       [MaxNJetSmallR];
        float            fTJsmallPhi       [MaxNJetSmallR];
        float            fTJsmallM         [MaxNJetSmallR];
        float            fTJsmallCharge    [MaxNJetSmallR];
	float            fTJsmallTrackpT    [MaxNJetSmallR];
	float            fTJsmallTrackMass    [MaxNJetSmallR];
        float            fTJsmallTrackpTR    [MaxNJetSmallR];
        float            fTJsmallTrackMassR    [MaxNJetSmallR];
	int              fTJsmallBtag      [MaxNJetSmallR];
	int              fTJsmallCtag      [MaxNJetSmallR];
	int              fTJsmallFtag      [MaxNJetSmallR];
        int              fTJsmallWPicked   [MaxNJetSmallR];
	int              fTJsmallBPicked   [MaxNJetSmallR];
	float            fTJsmallMg         [MaxNJetSmallR];
	float            fTJsmalltau21         [MaxNJetSmallR];
	int            fTJsmalltype         [MaxNJetSmallR];
	int fTsmallntrack[MaxNJetSmallR];
	int fTsmallntrackR[MaxNJetSmallR];
	
	float            fTJsmallScale        [MaxNJetSmallR];
	float            fTJsmallAngle        [MaxNJetSmallR];
	float            fTJsmallDrop        [MaxNJetSmallR];
	float            fTJsmallAdd        [MaxNJetSmallR];

	float            fTJsmallScaleP        [MaxNJetSmallR];
        float            fTJsmallAngleP        [MaxNJetSmallR];
        float            fTJsmallDropP        [MaxNJetSmallR];
        float            fTJsmallAddP        [MaxNJetSmallR];

	float fTtopmass; 
	float fTtopmass2;

        static const int MaxNJetLargeR = 20;
        int              fTNJetsLargeRFilled;
        float            fTJlargeRPt        [MaxNJetLargeR];
        float            fTJlargeREta       [MaxNJetLargeR];
        float            fTJlargeRPhi       [MaxNJetLargeR];
        float            fTJlargeRM         [MaxNJetLargeR];
        float            fTJlargeRMungroomed[MaxNJetLargeR];
        float            fTJlargeRCharge    [MaxNJetLargeR];
        int              fTJlargeRBtag      [MaxNJetLargeR];
        int              fTJlargeWplusMatch [MaxNJetLargeR];
        int              fTJlargeWminusMatch[MaxNJetLargeR];
        int              fTJlargeZpMatch    [MaxNJetLargeR];
        int              fTJlargeHpMatch    [MaxNJetLargeR];
	
	int fTmode;

        fastjet::JetDefinition     *m_jet_def;
        fastjet::JetDefinition     *m_jet_def_ca;



    public:
        myexampleAnalysis ();
        ~myexampleAnalysis ();
        
        void Begin();
        void AnalyzeEvent(int iEvt, Pythia8::Pythia *pythia8);
        void End();
        void DeclareBranches();
        void ResetBranches();
        void Debug(int debug){
            fDebug = debug;
        }
        void SetOutName(string outname){
            fOutName = outname;
        }
        
};

#endif

