from ROOT import *

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
'''

myfile = TFile("test.root")
mytree = myfile.Get("EventTree")

pThisto = TH1F("","",20,500,1000)
pThisto2 = TH1F("","",20,500,1000)
pThistoq = TH1F("","",20,500,1000)
pThistog = TH1F("","",20,500,1000)

etahisto = TH1F("","",20,-5,5)
etahisto2 = TH1F("","",20,-5,5)

masshisto = TH1D("","",20,0,300)
masshistoq = TH1D("","",20,0,300)
masshistog = TH1D("","",20,0,300)

masshisto2 = TH1D("","",20,0,300)
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
templates_q = {}
templates_g = {}
for i in range(0,22):
    templates_q[i] = TH1D("","",20,0,300)
    templates_g[i] = TH1D("","",20,0,300)

for i in range(mytree.GetEntries()):
    mytree.GetEntry(i)
    if (mytree.NJetsFilledSmallR > 0):
        pThisto.Fill(mytree.JsmallPt[0])
        etahisto.Fill(mytree.JsmallEta[0])
        masshisto.Fill(mytree.JsmallM[0])
        if (mytree.Jsmalltype[0]==21):
            masshistog.Fill(mytree.JsmallM[0])
            template_g.Fill(mytree.JsmallPt[0],mytree.JsmallM[0])
            templates_g[pThisto.GetXaxis().FindBin(mytree.JsmallPt[0])].Fill(mytree.JsmallM[0])
        else:
            masshistoq.Fill(mytree.JsmallM[0])
            template_q.Fill(mytree.JsmallPt[0],mytree.JsmallM[0])
            templates_q[pThisto.GetXaxis().FindBin(mytree.JsmallPt[0])].Fill(mytree.JsmallM[0])
            pass
        pass
    if (mytree.NJetsFilledSmallR > 1):
        pThisto2.Fill(mytree.JsmallPt[1])
        etahisto2.Fill(mytree.JsmallEta[1])
        masshisto2.Fill(mytree.JsmallM[1])
        if (mytree.Jsmalltype[1]==21):
            masshistog2.Fill(mytree.JsmallM[1])
        else:
            masshistoq2.Fill(mytree.JsmallM[1])
    pass

nreps = 5
def fillrandom(templ, pT, *histos):
    for _ in range(nreps):
        random_mass = templ[pThisto.GetXaxis().FindBin(pT)].GetRandom()
        for h in histos:
            h.Fill(random_mass)
for i in range(mytree.GetEntries()):
    mytree.GetEntry(i)
    if (mytree.NJetsFilledSmallR > 0):
        if (mytree.Jsmalltype[0]==21):
            fillrandom(templates_g, mytree.JsmallPt[0], masshisto_built, masshistog_built)
        else:
            fillrandom(templates_q, mytree.JsmallPt[0], masshisto_built, masshistoq_built)
    if (mytree.NJetsFilledSmallR > 1):
        if (mytree.Jsmalltype[1]==21):
            fillrandom(templates_g, mytree.JsmallPt[1], masshisto_built2, masshistog_built2)
        else:
            fillrandom(templates_q, mytree.JsmallPt[1], masshisto_built2, masshistoq_built2)
    pass

def output(masshisto, masshistoq, masshistog, masshisto_built, masshistoq_built, masshistog_built, pThisto, etahisto, suffix):
    c = TCanvas("a","a",500,500)
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

    pThisto.Draw()
    c.Print("pTspectrum_%s.pdf" % suffix)

    etahisto.Draw()
    c.Print("eta_%s.pdf" % suffix)

output(masshisto, masshistoq, masshistog, masshisto_built, masshistoq_built, masshistog_built, pThisto, etahisto, "leading")
output(masshisto2, masshistoq2, masshistog2, masshisto_built2, masshistoq_built2, masshistog_built2, pThisto2, etahisto2, "subleading")
