import sys
import argparse
import time

from ROOT import *
gRandom.SetSeed(int(time.time()))

def myText(x, y, text, color = 1):
    l = TLatex()
    l.SetTextSize(0.025)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x,y,text)
    pass

'''
This is the information inside the tree:

   TTreeReaderValue<Int_t> EventNumber = {fReader, "EventNumber"};
   TTreeReaderValue<Int_t> NJetsFilledSmallR = {fReader, "NJetsFilledSmallR"};
   TTreeReaderArray<Float_t> JsmallPt = {fReader, "JsmallPt"};
   TTreeReaderArray<Float_t> JsmallEta = {fReader, "JsmallEta"};
   TTreeReaderArray<Float_t> JsmallPhi = {fReader, "JsmallPhi"};
   TTreeReaderArray<Float_t> JsmallM = {fReader, "JsmallM"};
   TTreeReaderArray<Float_t> JsmallMg = {fReader, "JsmallMg"};
   TTreeReaderArray<Float_t> Jsmalltau21 = {fReader, "Jsmalltau21"};
   TTreeReaderArray<Int_t> Jsmalltype = {fReader, "Jsmalltype"};

run with "--template <source>" to change the template source where
<source> = ('random', 'subl', 'lead' (default))
'''
parser = argparse.ArgumentParser(description="closure tests for jet substructure templates")
parser.add_argument('--template', type=str, choices=('random', 'subl', 'lead'), default='lead',
                    help="template source ('lead' = leading jet; 'subl' = subleading jet;"
                    "'random' = random jet)")
parser.add_argument('--limit', type=int, default=0, help="stop after this many events")
args = parser.parse_args()

myfile = TFile("in.root")
mytree = myfile.Get("EventTree")

pThisto = TH1F("","",20,500,1000)
pThisto2 = TH1F("","",20,500,1000)
pThistoq = TH1F("","",20,500,1000)
pThistog = TH1F("","",20,500,1000)

etahisto = TH1F("","",20,-4,4)
etahisto2 = TH1F("","",20,-4,4)

eta_pT_histo = TH2F("", "", 20, 500, 1000, 20, -4, 4)
eta_pT_histo2 = TH2F("", "",  20, 500, 1000, 20, -4, 4)

templ_jetstr = '({})'.format(args.template)
masshisto = TH1D("","Template jet mass {};m (MeV);".format(templ_jetstr),20,0,300)
masshistoq = TH1D("","",20,0,300)
masshistog = TH1D("","",20,0,300)

test_jetstr = '(random)' if args.template == 'random' else '(lead)' if args.template == 'subl' else '(subl)'
masshisto2 = TH1D("","Test jet mass {};m (MeV)".format(test_jetstr),20,0,300)
masshistoq2 = TH1D("","",20,0,300)
masshistog2 = TH1D("","",20,0,300)

masshisto_built = TH1D("","",20,0,300)
masshistoq_built = TH1D("","",20,0,300)
masshistog_built = TH1D("","",20,0,300)

masshisto_built2 = TH1D("","",20,0,300)
masshistoq_built2 = TH1D("","",20,0,300)
masshistog_built2 = TH1D("","",20,0,300)

template_q = TH2D("","",20,500,1000,20,0,300)
template_g = TH2D("","",20,500,1000,20,0,300)
templates_q = [TH1D("","",20,0,300) for _ in range(22)]
templates_g = [TH1D("","",20,0,300) for _ in range(22)]

# random selection of 1's and zeros. We need this because we have to use the opposite jet for testing
# ie, we need the same random number the second time we loop through
randomjets = [gRandom.Integer(2) for _ in range(mytree.GetEntries())]
random_template_jet = lambda evt_num: randomjets[evt_num]
random_test_jet = lambda evt_num: 0 if randomjets[evt_num] == 1 else 1

# utility functions
index_for_template_jet = lambda evt_num: 1 if args.template == 'subl' else random_template_jet(evt_num) if args.template == 'random' else 0
index_for_test_jet = lambda evt_num: 1 if index_for_template_jet(evt_num) == 0 else 0
enough_jets = lambda evt: evt.NJetsFilledSmallR > 1 or (evt.NJetsFilledSmallR == 1 and args.template not in ('random', 'subl'))

def fill_templates(i, evt):
    jet_i = index_for_template_jet(i)
    template, templates = (template_g, templates_g) if evt.Jsmalltype[jet_i] == 21 else (template_q, templates_q)
    template.Fill(evt.JsmallPt[jet_i], evt.JsmallM[jet_i])
    templates[pThisto.GetXaxis().FindBin(mytree.JsmallPt[jet_i])].Fill(evt.JsmallM[jet_i])

def fill_histos(i, evt):
    template_jet_i = index_for_template_jet(i)
    test_jet_i = index_for_test_jet(i)

    if evt.NJetsFilledSmallR > test_jet_i:
        masshisto2_qg = masshistog2 if evt.Jsmalltype[test_jet_i] == 21 else masshistoq2
        masshisto2.Fill(mytree.JsmallM[test_jet_i])
        masshisto2_qg.Fill(mytree.JsmallM[test_jet_i])
        pThisto2.Fill(mytree.JsmallPt[test_jet_i])
        etahisto2.Fill(mytree.JsmallEta[test_jet_i])
        eta_pT_histo.Fill(mytree.JsmallPt[test_jet_i], mytree.JsmallEta[test_jet_i])
    if evt.NJetsFilledSmallR > template_jet_i:
        masshisto_qg = masshistog if evt.Jsmalltype[template_jet_i] == 21 else masshistoq
        masshisto.Fill(mytree.JsmallM[template_jet_i])
        masshisto_qg.Fill(mytree.JsmallM[template_jet_i])
        pThisto.Fill(mytree.JsmallPt[template_jet_i])
        etahisto.Fill(mytree.JsmallEta[template_jet_i])
        eta_pT_histo2.Fill(mytree.JsmallPt[template_jet_i], mytree.JsmallEta[template_jet_i])

for i in range(min(args.limit or mytree.GetEntries(), mytree.GetEntries())):
    mytree.GetEntry(i)
    if enough_jets(mytree):
        fill_templates(i, mytree)
        fill_histos(i, mytree)

nreps = 5
def fillrandom(templ, pT, *histos):
    """
    Fill the histos with the same @nreps random values drawn from the template in @templ corresponding to @pT
    """
    for _ in range(nreps):
        random_mass = templ[pThisto.GetXaxis().FindBin(pT)].GetRandom()
        for h in histos:
            h.Fill(random_mass)

for i in range(min(args.limit or mytree.GetEntries(), mytree.GetEntries())):
    mytree.GetEntry(i)
    if enough_jets(mytree):
        template_jet_i = 1 if args.template == 'subl' else random_template_jet(i) if args.template == 'random' else 0
        test_jet_i = 0 if template_jet_i == 1 else 1

        # the template jet always exists, or we would've skipped the event
        if mytree.Jsmalltype[template_jet_i] == 21:
            fillrandom(templates_g, mytree.JsmallPt[template_jet_i], masshisto_built, masshistog_built)
        else:
            fillrandom(templates_q, mytree.JsmallPt[template_jet_i], masshisto_built, masshistoq_built)

        # The test jet may not exist
        if mytree.NJetsFilledSmallR > test_jet_i:
            if mytree.Jsmalltype[test_jet_i] == 21:
                fillrandom(templates_g, mytree.JsmallPt[test_jet_i], masshisto_built2, masshistog_built2)
            else:
                fillrandom(templates_q, mytree.JsmallPt[test_jet_i], masshisto_built2, masshistoq_built2)


c = TCanvas("a","a",500,500)
def output(masshisto, masshistoq, masshistog, masshisto_built, masshistoq_built, masshistog_built, pThisto, etahisto, suffix):
    gStyle.SetOptStat(0)
    gPad.SetLeftMargin(0.15)
    gPad.SetTopMargin(0.05)
    masshisto.Draw()
    masshisto.GetYaxis().SetTitle("arbitrary units")
    masshisto.GetXaxis().SetTitle("m_{MMDT} [GeV]")
    masshisto.GetYaxis().SetTitleOffset(1.8)
    masshisto.GetXaxis().SetTitleOffset(1.3)
    masshisto.SetLineColor(1)
    masshisto.GetYaxis().SetRangeUser(0,masshisto.GetMaximum()*1.2)

    masshistoq.SetLineColor(2)
    masshistoq.SetLineStyle(3)
    masshistoq.Draw("same")

    masshistog.SetLineColor(4)
    masshistog.SetLineStyle(3)
    masshistog.Draw("same")

    masshisto_built.SetLineColor(1)
    masshisto_built.SetMarkerStyle(20)
    masshisto_built.SetMarkerSize(1)
    masshisto_built.Scale(masshisto.Integral()/masshisto_built.Integral())
    masshisto_built.Draw("samep")

    masshistoq_built.SetLineColor(2)
    masshistoq_built.SetLineStyle(3)
    masshistoq_built.SetMarkerColor(2)
    masshistoq_built.SetMarkerStyle(20)
    masshistoq_built.SetMarkerSize(1)
    masshistoq_built.Scale(masshistoq.Integral()/masshistoq_built.Integral())
    masshistoq_built.Draw("samep")

    masshistog_built.SetLineColor(2)
    masshistog_built.SetLineStyle(3)
    masshistog_built.SetMarkerColor(4)
    masshistog_built.SetMarkerStyle(20)
    masshistog_built.SetMarkerSize(1)
    masshistog_built.Scale(masshistog.Integral()/masshistog_built.Integral())
    masshistog_built.Draw("samep")


    leg = TLegend(.45, .4, .85, .7)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry(masshisto,"Total","L")
    leg.AddEntry(masshistoq,"quarks","L")
    leg.AddEntry(masshistog,"gluons","L")
    leg.AddEntry(masshisto_built,"from templates, total","ep")
    leg.AddEntry(masshistoq_built,"from templates, quarks","ep")
    leg.AddEntry(masshistog_built,"from templates, gluons","ep")
    leg.Draw()

    myText(0.2,0.9,"#scale[1.5]{#bf{Pythia 8.226 QCD dijets}}")

    c.Print("masshisto_%s.pdf" % suffix)


pThisto.SetLineColor(2)
pThisto.Scale(1.0/pThisto.Integral())
pThisto.Draw()
pThisto2.SetLineColor(4)
pThisto2.Scale(1.0/pThisto2.Integral())
pThisto2.Draw("same")
c.Print("pTspectrum (templ = {}).pdf".format(args.template))

etahisto.SetLineColor(2)
etahisto.Scale(1.0/etahisto.Integral())
etahisto.Draw()
etahisto2.SetLineColor(4)
etahisto2.Scale(1.0/etahisto2.Integral())
etahisto2.Draw("same")
c.Print("eta (templ = {}).pdf".format(args.template))

eta_pT_histo.Draw("lego2")
c.Print("eta_pt template {}.pdf".format(templ_jetstr))

eta_pT_histo2.Draw("lego2")
c.Print("eta_pt test {}.pdf".format(test_jetstr))

eta_pT_histo.Scale(eta_pT_histo2.Integral() / eta_pT_histo.Integral())
diff_histo = eta_pT_histo2 - eta_pT_histo
diff_histo.Draw("lego2")
c.Print("eta_pT_diff (templ = {}).pdf".format(args.template))

output(masshisto, masshistoq, masshistog, masshisto_built, masshistoq_built, masshistog_built, pThisto, etahisto, "template {}".format(templ_jetstr))
output(masshisto2, masshistoq2, masshistog2, masshisto_built2, masshistoq_built2, masshistog_built2, pThisto2, etahisto2, "test {}".format(test_jetstr))



