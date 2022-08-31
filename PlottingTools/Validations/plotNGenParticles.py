import sys
import os
import ROOT
import numpy as np
from tdrstyle import *
setTDRStyle()
import CMS_lumi as CMS_lumi
from ROOT import gStyle, TH1F, TCanvas, TLegend, TLatex, TEfficiency, SetOwnership, TH2F
from ROOT import kRed, kBlue, kGreen, kMagenta, kOrange, kCyan, kBlack

kcolors = [kRed, kBlue, kGreen+2, kMagenta+2, kOrange, kCyan+1, kBlack]
markers = [24, 25, 26, 27, 32, 28, 21]


iPeriod = 0
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12

treename = "syncTree"

#_______________________________________________________________________________
def drawText(text, x=0.17, y=0.35, font_size=0.):
    tex = TLatex(x, y,text)
    if font_size > 0.:
      tex.SetTextSize(font_size)
      tex.SetTextSize(0.05)
      tex.SetNDC()
      tex.Draw()
      return tex

def plotNGenParticles(rfile, variable, xbins, xtitle, cut, text, plotsuffix, plotfolder):
    nbin = xbins[0]; xmin = xbins[1]; xmax = xbins[2]

    c1 = ROOT.TCanvas("c1","c1",800,600)
    #c = newCanvas()
    tch = ROOT.TChain(treename)
    tch.Add(rfile)
    hist= ROOT.TH1F("hist","hist",xbins[0], xbins[1], xbins[2])
    tch.Draw(variable+">>hist", cut)
    hist.SetMaximum(hist.GetBinContent(hist.GetMaximumBin()) * 1.4)
    hist.GetXaxis().SetLabelSize(0.05)
    hist.GetYaxis().SetLabelSize(0.05)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle("Events")
    hist.SetMarkerSize(1.5)
    hist.Draw("histtext")
    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)

    txt = drawText(text, 0.2,0.9,0.04)

    c1.SaveAs(plotfolder+"genParticles_HHbbWW_"+plotsuffix+".pdf")
    del hist,c1,tch


def plotNGenLepPt(rfile, variables, xbins, xtitle, cuts, legnames, text, plotsuffix, plotfolder):
    nbin = xbins[0]; xmin = xbins[1]; xmax = xbins[2]

    c1 = ROOT.TCanvas("c1","c1",800,600)
    #c = newCanvas()
    tch = ROOT.TChain(treename)
    tch.Add(rfile)
    hist = []
    ymax = 0.0
    for i,var in enumerate(variables):
        cut  = cuts[i]
        hist.append( ROOT.TH1F("hist%d"%i,"hist%d"%i,xbins[0], xbins[1], xbins[2]))
        tch.Draw(var+">>hist%d"%i, cut)
        entries = hist[i].GetEntries()
        bincontent = hist[i].GetBinContent(nbin)+hist[i].GetBinContent(nbin+1)
        hist[i].SetBinContent(nbin, bincontent)
        hist[i].Scale(1.0/entries)
        hist[i].SetLineColor(kcolors[i])
        if ymax < hist[i].GetBinContent(hist[i].GetMaximumBin()):
            ymax = hist[i].GetBinContent(hist[i].GetMaximumBin())
    base  = TH1F("base","base",nbin,xmin,xmax)
    base.SetMinimum(0)
    base.SetMaximum(ymax*1.2)
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
    base.GetXaxis().SetTitleSize(0.05)
    base.GetYaxis().SetTitleSize(0.05)
    base.GetXaxis().SetTitle(xtitle)
    base.GetYaxis().SetTitle("Normalized to unity")
    base.Draw("")
    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)

    leg = TLegend(0.45,0.66,.9,0.7+len(legnames)*0.05, "", "brNDC");
    leg.SetHeader("HH#rightarrow bbWW")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)

    for i,h in enumerate(hist):
        h.Draw("histsame")
        leg.AddEntry(h, legnames[i],"l")
    leg.Draw("same");

    txt = drawText(text, 0.2,0.9,0.04)

    c1.SaveAs(plotfolder+"genParticles_HHbbWW_"+plotsuffix+".pdf")
    del hist,c1,tch,txt,leg




SLdict = {
        "nGenX":"num of X,X#rightarrow HH",
        "nGenHiggs":"num of Higgs,X#rightarrow HH",
        "nGenW1FromHiggs":"num of W1,H#rightarrow W1W2",
        "nGenW2FromHiggs":"num of W2,H#rightarrow W1W2",
        "nGenBQuarkFromHiggs": "num of bquarks, H#rightarrow bb",
        "nGenLepFromW1FromHiggs":"num of ele/mu,W1#rightarrow l#nu",
        "nGenLepFromTauFromW1FromHiggs":"num of ele/mu,W1#rightarrow #tau#nu#rightarrow l+3#nu",
        "nGenQuarkFromW2FromHiggs" : "num of quark(abs(pdgid)=1,2,3,4,5), W2#rightarrow qq",
        }
DLdict = {
        "nGenX":"num of X,X#rightarrow HH",
        "nGenHiggs":"num of Higgs,X#rightarrow HH",
        "nGenW1FromHiggs":"num of W1,H#rightarrow W1W2",
        "nGenW2FromHiggs":"num of W2,H#rightarrow W1W2",
        "nGenBQuarkFromHiggs": "num of bquarks, H#rightarrow bb",
        "nGenLepFromW1FromHiggs":"num of ele/mu,W1#rightarrow l#nu",
        "nGenLepFromW2FromHiggs":"num of ele/mu,W2#rightarrow l#nu",
        "nGenLepFromTauFromW1FromHiggs":"num of ele/mu,W1#rightarrow #tau#nu#rightarrow l+3#nu",
        "nGenLepFromTauFromW2FromHiggs":"num of ele/mu,W2#rightarrow #tau#nu#rightarrow l+3#nu",
        }
def plotAllNGen(plotter):
    xbins = [3,0,3]
    cut = "1"
    Plotfolder = os.path.join(plotter.baseDir,"genParticlesPlots/")
    if not os.path.isdir(Plotfolder): os.mkdir(Plotfolder)
    rfile = plotter.inputFile
    channel = "DL_"+plotter.shortName
    thisdict = DLdict
    if plotter.isSL: thisdict = SLdict
    if plotter.isSL: channel = "SL_"+plotter.shortName
    year = plotter.runyear

    leppts = ["GenLepFromW1FromHiggs_pt","GenLepFromTauFromW1FromHiggs_pt"]
    legnames = ["ele/mu from W1#rightarrow ele/mu+#nu","ele/mu from W1#rightarrow#tau#nu#rightarrow ele/mu+3#nu"]
    lepcuts = ["nGenLepFromW1FromHiggs>0","nGenLepFromTauFromW1FromHiggs>0"]
    plotsuffix =  "%s_%d_"%(channel,year)+"_".join(leppts)
    ptbins = [50,0,100.0]
    text = "Run%d,%s"%(year, channel)
    plotNGenLepPt(rfile, leppts, ptbins, "genlep pT [GeV]", lepcuts, legnames, text, plotsuffix, Plotfolder)
    for variable in thisdict.keys():
        xtitle = thisdict[variable]
        if "Nonres" in channel: xtitle = xtitle.replace("X","GG")
        plotsuffix = "%s_%d_"%(channel,year)+variable
        #plotNGenParticles(rfile, variable, xbins, xtitle, cut, text, plotsuffix, Plotfolder)

