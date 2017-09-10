# EventGeneration

UFO instructions (from Felix Yu)

1.  Untar the files and put them in the "models/" subdirectory of a MadGraph5 installation.
2.  At the command line, run "./bin/mg5" from the MadGraph5 installation directory.
3.  Import the Zprime model.  The command is:
mg5> import model Zprime_UFO

4.  Generate your signal Feynman graphs.  For Zprimes, this is the step that distinguishes resonances below and above the top pair threshold.  For Zprimes below the top pair threshold, use
mg5> generate p p > zpl, zpl > j j
mg5> add process p p > zpl, zpl > b b~

while for Zprimes above 350 GeV, use
mg5> generate p p > zph, zph > j j
mg5> add process p p > zph, zph > b b~

5.  Output the graphs to a directory:
mg5> output Dijet_signal

6.  You can either "launch" the mg5 event generation interface, or work from the command line to generate events.  In either case, you have to edit the param_card.dat to adjust the coupling, denoted by "gz" in the "Block frblock" as well as the resonance mass "MZpl", or "MZph", in the "Block mass".  Then, edit the run_card.dat for the appropriate center of mass energy and number of events, and finally, you can execute the event generation.  From the command line, I use
$ cd Dijet_signal/
$ ./bin/generate_events 2 4 test_14TeV
where "2" denotes a multicore system, "4" denotes the number of cores to use, and "test_14TeV" is the tag name of the event run.

You should then see your events in the "Dijet_signal/Events/test_14TeV/" directory, which you can parse and analyze.
