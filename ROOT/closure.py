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
pThistoq = TH1F("","",20,500,1000)
pThistog = TH1F("","",20,500,1000)

masshisto = TH1D("","",20,0,300)
masshistoq = TH1D("","",20,0,300)
masshistog = TH1D("","",20,0,300)

masshisto_built = TH1D("","",20,0,300)

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
    pass

nreps = 5
for i in range(mytree.GetEntries()):
    mytree.GetEntry(i)
    if (mytree.NJetsFilledSmallR > 0):
        if (mytree.Jsmalltype[0]==21):
            for j in range(nreps):
                masshisto_built.Fill(templates_g[pThisto.GetXaxis().FindBin(mytree.JsmallPt[0])].GetRandom())
                pass
        else:
            for j in range(nreps):
                masshisto_built.Fill(templates_q[pThisto.GetXaxis().FindBin(mytree.JsmallPt[0])].GetRandom())
                pass
            pass
    pass
        
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

leg = TLegend(.5,.4,.85,.7)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(masshisto,"Total","L")
leg.AddEntry(masshistoq,"quarks","L")
leg.AddEntry(masshistog,"gluons","L")
leg.AddEntry(masshisto_built,"from templates","ep")
leg.Draw()

myText(0.2,0.9,"#scale[1.5]{#bf{Pythia 8.226 QCD dijets}}")

c.Print("masshisto.pdf")

pThisto.Draw()
c.Print("pTspectrum.pdf")
