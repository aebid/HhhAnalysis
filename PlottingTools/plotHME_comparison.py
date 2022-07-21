import ROOT
from ROOT import kRed, kBlue, kGreen, kMagenta, kOrange, kCyan, kBlack
import numpy as np
import sys
import os
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.argv.append( '-b' )
sys.argv.append( '-q' )


kcolors = [kRed, kBlue, kGreen+2, kMagenta+2, kOrange, kCyan+1, kBlack]
markers = [24, 25, 26, 27, 32, 28, 21]
rfiledir = "/eos/user/t/tahuang/ReRunHME_Florian_2016_20220720/"
rfiledir = "/eos/user/t/tahuang/ReRunHME_Florian_2016_20220720_it50000_gp0p18/"
Plotfolder = rfiledir+"compareHMEPlots/"
if not os.path.exists(Plotfolder):
    os.system("mkdir "+Plotfolder)
treename = "Events"

### plot HME from one file
def plotHME(rfile, legnames, variables, xbins, cut, text, plotsuffix):
    nbin = xbins[0]; xmin = xbins[1]; xmax = xbins[2]

    histlist = []

    xleg = 0.6
    if "800" in text or "850" in text or "900" in text or "1000" in text:
      xleg = 0.15
    legend = ROOT.TLegend(xleg,0.65,xleg+0.3,0.7+len(legnames)*.045); 
    c1 = ROOT.TCanvas("c1","c1",800,600)
    ROOT.gPad.SetGrid(1,1)
    total_env = None
    hs = ROOT.THStack("Shapecomparison", "Reco HME distribution")
    chain = ROOT.TChain(treename)
    chain.Add(rfile)
    for i,legname in enumerate(legnames):
        variable = variables[i]
        #weight = "total_weight"
        #finalcut = cut+"&& %s >= 250.0"%variable
	finalcut = cut
        histlist.append(ROOT.TH1F("hist%d"%i,"hist%d"%i,xbins[0], xbins[1], xbins[2]))

        chain.Draw(variable + ">> hist%d"%i, finalcut)
        integral = histlist[i].Integral()
        mean = histlist[i].GetMean()
        entries = histlist[i].GetEntries()
        histlist[i].SetBinContent(1, histlist[i].GetBinContent(0)+histlist[i].GetBinContent(1))
        histlist[i].SetBinContent(nbin, histlist[i].GetBinContent(nbin+1)+histlist[i].GetBinContent(nbin))
        histlist[i].SetLineColor(kcolors[i])
        histlist[i].SetLineWidth(2)
        histlist[i].SetStats(0)
        legend.AddEntry(histlist[i], legname+": %d events"%integral,"l")
        hs.Add(histlist[i])

    ##GetEntries() return underflow+Integral(1, Nbins)+overflow
    legend.SetHeader(text+" total:%d events"%histlist[0].GetEntries())
    #c1.SetLogy()
    hs.Draw("nostackhist")
    legend.Draw("same")
    hs.GetHistogram().GetXaxis().SetTitle("Reco HME [GeV]")
    hs.GetHistogram().GetYaxis().SetTitle("Events")
    c1.SaveAs(Plotfolder+"HMEdistribution_"+plotsuffix+".pdf")
    c1.SaveAs(Plotfolder+"HMEdistribution_"+plotsuffix+".C")

"""
def plotDY2D(DY, xvar, yvar, xbins, ybins, cut, plotsuffix):
    nbins = xbins[0]; xmin = xbins[1]; xmax = xbins[2]
    xbins = []
    binwidth = (xmax-xmin)/nbins
    for i in range(0, nbins+1):
        xbins.append(xmin + i*binwidth)
    xbins = np.asarray(xbins)
    ybins = np.asarray(ybins)
    histname = DY+"_"+xvar+"_"+yvar
    hist = ROOT.TH2F(histname, histname, len(xbins)-1, xbins, len(ybins) -1, ybins)
    hist.SetStats(0)

    c1 = ROOT.TCanvas("c1","c1",800,600)
    rfile = bfile_dict[DY]["path"]
    weight = bfile_dict[DY]["weight"]
    #weight = "total_weight"
    xsec,event_weight_sum = get_xsection_eventweightsum_file(rfile, treename) 
    #event_weight_sum = bfile_dict[dy]["event_weight_sum"]
    relative_weight = TotalLumi*1000*xsec/event_weight_sum
    #weight = "total_weight"
    tch = ROOT.TChain(treename)
    tch.Add(rfile)
    finalcut = cut+"*%f"%relative_weight +"* %s "%weight
    finalcut = cut
    tch.Draw("%s:%s"%(yvar, xvar) + ">>"+histname, finalcut,"colz")
    hist.GetXaxis().SetTitle(xvar)
    hist.GetYaxis().SetTitle(yvar)
    #hist.SetMinimum(-15)
    #hist.SetMaximum(15)

    c1.SaveAs(Plotfolder+histname+"_"+plotsuffix+"_noMaxMin.pdf")



def plotDY_datadriven(DY, variable, xbins, cut, plotsuffix):
    nbin = xbins[0]; xmin = xbins[1]; xmax = xbins[2]
    #h1 = ROOT.TH1F("hist1","hist1",xbins[0], xbins[1], xbins[2])

    histlist = []
    tchlist = []


    legend = ROOT.TLegend(0.65,0.7,0.88,0.7+4*.05); 
    c1 = ROOT.TCanvas("c1","c1",800,600)

    samplename = None; channel = "ee"
    if "isMuMu" in cut:
        samplename = "DoubleMuon"
        channel = "#mu#mu"
    elif "isElEl" in cut:
        samplename = "DoubleEG"

    hs = ROOT.THStack("Shapecomparison", "Drell-Yan distribution in %s channel"%channel)
    for i,dy in enumerate(DY):
        #print("keys in full_local_samplelist ",full_local_samplelist[dy].keys())
        if "Data" not in dy:
            #samplename = full_local_samplelist[dy].keys()[0]
            samplename = list(full_local_samplelist[dy].keys())[0] ## for python3
        rfile = full_local_samplelist[dy][samplename]["path"]
        print("samplename ", samplename, " rfile ",rfile)
        #weight = bfile_dict[dy]["weight"]
        #weight = "total_weight"
        weight = "dy_Mbtag_weight"
        relative_weight = 1.0
        if "Data" not in dy:
            xsec,event_weight_sum = get_xsection_eventweightsum_file(rfile, treename) 
            relative_weight = TotalLumi*1000*xsec/event_weight_sum
            weight = "dy_Mbtag_weight*sample_weight*event_reco_weight"
            print("dy ",dy, " xsec ",xsec, " weight sum ", event_weight_sum," relative weight ",relative_weight)
        #weight = "total_weight"
        histlist.append(ROOT.TH1F("hist%d"%i,"hist%d"%i,xbins[0], xbins[1], xbins[2]))
        tchlist.append(ROOT.TChain(treename))
        tchlist[i].Add(rfile)
        finalcut = cut+"*%f"%relative_weight +"* %s "%weight
        finalcut = cut
        tchlist[i].Draw(variable + ">> hist%d"%i, finalcut)

        #integral = histlist[i].Integral()
        mean = histlist[i].GetMean()
        entries = histlist[i].GetEntries()
        print("mean ",mean)
        histlist[i].SetLineColor(kcolors[i])
        histlist[i].SetLineWidth(2)
        histlist[i].SetStats(0)
        #legend.AddEntry(histlist[i], dy+": %d events"%abs(histlist[i].Integral()),"l")
        legend.AddEntry(histlist[i], dy+": %d events"%histlist[i].Integral(),"l")
        hs.Add(histlist[i])
        #if i == 0:
        #    histlist[i].SetTitle("Drell-Yan distribution, no btagging")
        #    histlist[i].GetXaxis().SetTitle(variable)
        #    histlist[i].Draw()
        #else:
        #    histlist[i].Draw("same")


    c1.SetLogy()
    hs.Draw("nostackhist")
    legend.Draw("same")
    hs.GetHistogram().GetXaxis().SetTitle(variable)
    hs.GetHistogram().GetYaxis().SetTitle("Events")
    c1.SaveAs(Plotfolder+"DY_datadriven_"+variable+"_"+plotsuffix+".pdf")

def plotDY2D_datadriven(DY, xvar, yvar, xbins, ybins, cut, plotsuffix):
    nbins = xbins[0]; xmin = xbins[1]; xmax = xbins[2]
    xbins = []
    binwidth = (xmax-xmin)/nbins
    for i in range(0, nbins+1):
        xbins.append(xmin + i*binwidth)
    xbins = np.asarray(xbins)
    ybins = np.asarray(ybins)
    histname = DY+"datadriven_"+xvar+"_"+yvar
    hist = ROOT.TH2F(histname, histname, len(xbins)-1, xbins, len(ybins) -1, ybins)
    hist.SetStats(0)

    c1 = ROOT.TCanvas("c1","c1",800,600)
    samplename = None; channel = "ee"
    if "isMuMu" in cut:
        samplename = "DoubleMuon"
        channel = "#mu#mu"
    elif "isElEl" in cut:
        samplename = "DoubleEG"

    if "Data" not in DY:
        samplename = list(full_local_samplelist[DY].keys())[0]

    rfile = full_local_samplelist[DY][samplename]["path"]
    weight = "dy_Mbtag_weight"
    relative_weight = 1.0
    if "Data" not in DY:
        xsec,event_weight_sum = get_xsection_eventweightsum_file(rfile, treename) 
        relative_weight = TotalLumi*1000*xsec/event_weight_sum
        weight = "dy_Mbtag_weight*sample_weight*event_reco_weight"
        print("dy ",DY, " xsec ",xsec, " weight sum ", event_weight_sum," relative weight ",relative_weight)
    #weight = "total_weight"
    tch = ROOT.TChain(treename)
    tch.Add(rfile)
    finalcut = cut+"*%f"%relative_weight +"* %s "%weight
    finalcut = cut
    tch.Draw("%s:%s"%(yvar, xvar) + ">>"+histname, finalcut,"colztext")
    hist.GetXaxis().SetTitle(xvar)
    hist.GetYaxis().SetTitle(yvar)
    hist.SetTitle(DY+", "+channel+" channel")
    hist.SetMinimum(0)
    #hist.SetMaximum(15)

    c1.SaveAs(Plotfolder+histname+"_"+plotsuffix+"_noMaxMin.pdf")


#plotsuffix = "nobtag_to_Mbtag_20201209"
#plotsuffix = "nobtag_20201209_2btaggingcut_ElEl"
#plotsuffix = "nobtag_20201209_ElEl"
plotsuffix = "2MbtagfromWeight_20210112_MuMu"
DYprocess = ["DYToLL2J","DYToLL1J","DYToLL0J","DYM10to50"]
DYprocess_untagged = ["Data_untagged", "TT_untagged"]

#DYprocess = ["DYToLL1J","DYToLL0J","DYM10to50"]
thecut = "((isMuEl || isElMu) && ll_M<91-15 && ll_M>12 && hme_h2mass_reco>250)"
thecut = "(ll_M<91-15 && ll_M>12 && hme_h2mass_reco>250 && nnout_MTandMT2_MJJ_M400>0.2)"
#thecut = "(ll_M<91-15 && ll_M>12 && hme_h2mass_reco>250)"
suffix = "nnoutcut0p2_20210112_noweight"
#suffix = "20210112_noweight"
#cut = "(isMuMu && ll_M<91-15 && ll_M>12 && hme_h2mass_reco>250 && jet1_cMVAv2 > 0.4432 && jet2_cMVAv2 > 0.4432)"
#cut = "(isElEl && ll_M<91-15 && ll_M>12 && hme_h2mass_reco>250 && jet1_cMVAv2 > 0.4432 && jet2_cMVAv2 > 0.4432)"
for ch in ["MuMu", "ElEl"]:
    #for ch in ["MuMu", "ElEl","ElMu"]:
    if ch == "MuMu":
        cut = "("+thecut +" && isMuMu)"
    elif ch == "ElEl":
        cut = "("+thecut +" && isElEl)"
    elif ch == "ElMu":
        cut = "("+thecut +" && (isElMu || isMuEl))"
    plotsuffix=suffix+ch
    #plotDY(DYprocess, "nnout_MTandMT2_MJJ_M400", [20, 0.0, 1.0], cut, plotsuffix)
    #plotDY(DYprocess, "hme_h2mass_reco", [20, 250, 1250.0],  cut, plotsuffix)
    #plotDY_datadriven(DYprocess_untagged, "nnout_MTandMT2_MJJ_M400", [20, 0.0, 1.0], cut, plotsuffix)
    #plotDY_datadriven(DYprocess_untagged, "hme_h2mass_reco", [20, 250, 1250.0],  cut, plotsuffix)
    #plotDY_datadriven(DYprocess_untagged, "dy_Mbtag_weight", [40, 0.0, 0.04],  cut, plotsuffix)


"""
#variables = ["hme_mass_peak","hme_mass_peak_divSol","hme_mass_peak_boosted","hme_mass_peak_boosted_divSol"]
#legnames =  ["2 subjets, weight case", "2 subjets, average case","fatjet, weight case","fatjet, average case"]
variables = ["HME", "hme_mass_peak","hme_mass_peak_gravity"]
legnames = ["Florian implementaion", "Tao implementaion", "Tao implementation, weighted"]

#masslist = [1000, 1250, 1500, 1750, 2000, 2500, 3000, 700, 800, 900]
#masslist = [260, 300, 400, 500, 600, 750, 800, 900]

masslist = [250, 260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650, 700,750, 800, 850, 900, 1000] 
for mass in masslist:
    rfile = rfiledir+ "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_HME_ijob*.root"%(mass)
    text = "Radion M=%d GeV"%mass
    cut  = "1"
    #xbins = [90,500.0, 3500.0] ## for boosted
    xbins = [60,220.0, 1420.0]
    plotsuffix = "Radion_M%d_weight_average"%mass
    plotHME(rfile, legnames, variables, xbins, cut, text, plotsuffix)
