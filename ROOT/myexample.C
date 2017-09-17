#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "TString.h"
#include "TSystem.h"
#include "TError.h"
//#include "TPythia8.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

#include "fastjet/PseudoJet.hh"  
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

#include "myTools.h"
#include "myexampleAnalysis.h"

#include "boost/program_options.hpp"

using std::cout;
using std::endl;
using std::string;
using std::map;
using namespace std;
namespace po = boost::program_options;

int getSeed(int seed){                                                                                                                                               
  if (seed > -1) return seed;                                                                                                                                      \
  int timeSeed = time(NULL);                                                                                                                                       \
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729);                                                                                                         \
}

int main(int argc, char* argv[]){
    // argument parsing  ------------------------
    cout << "Called as: ";
    for(int ii=0; ii<argc; ++ii){
        cout << argv[ii] << " ";
    }
    cout << endl;

    // agruments 
    int nEvents = 0;
    int fDebug  = 0;
    string outName = "test.root";

    po::options_description desc("Allowed options");
    desc.add_options() //100000
      ("help", "produce help message")
      ("NEvents", po::value<int>(&nEvents)->default_value(10000) ,    "Number of Events ")
      ("Debug",   po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
      ("OutFile", po::value<string>(&outName)->default_value("test.root"), "output file name")
      ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    }
    //------

   Pythia8::Pythia* pythia8b = new Pythia8::Pythia();

   int    seed      =-1;
   seed = getSeed(seed);

   // Configure and initialize pythia                                                                                                                                
   pythia8b->readString("Random:setSeed = on");
   std::stringstream ss; ss << "Random:seed = " << seed;
   cout << ss.str() << endl;
   pythia8b->readString(ss.str());

   //Sometimes, we may want to use a particular seed in order to reproduce results.
   //pythia8b->readString("Random:setSeed = on");
   //pythia8b->readString("Random:seed = 798631097");
   
   myexampleAnalysis * analysis1 = new myexampleAnalysis();

   pythia8b->readString("HardQCD:all  = on");
   pythia8b->readString("PhaseSpace:pTHatMin  = 400");
   pythia8b->readString("PhaseSpace:pTHatMax  = 800");
   
   pythia8b->readString("Beams:idA = 2212");
   pythia8b->readString("Beams:idB = 2212");
   pythia8b->readString("Beams:eCM = 13000");
   pythia8b->init();//2212 /* p */, 2212 /* p */, 13000. /* TeV */);
   
   analysis1->SetOutName(outName);
   analysis1->Begin();
   cout << "running on " << nEvents << " events " << endl;
   for (Int_t iev = 0; iev < nEvents; iev++) {
     if (iev%100==0) cout << iev << " " << nEvents << endl;
     analysis1->AnalyzeEvent(iev, pythia8b);
   }
   analysis1->End();
   pythia8b->stat(); 
   delete analysis1;
   
   // that was it
   delete pythia8b;
   return 0;
}
