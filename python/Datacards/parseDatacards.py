import os
import ROOT
import re
import numpy as np
import datetime
import Datacards
import sys


processlist = 	["TTbar","SingleTop","Drell_Yan","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV","Signal"]
statslist =     ["CMS_eff_b_heavy","CMS_eff_b_light","CMS_pu", "CMS_pdf", "CMS_eff_trigger","CMS_eff_e","CMS_eff_mu","CMS_iso_mu","QCDscale"]


def createShapes_1D(masslist, filesuffix, inputdir, prefix_out, outdir):
    for mass in masslist:
	processnames = ["TT","sT","DY","data_untagged","TT_untagged","sT_untagged","ttV","VV", "RadionM{}".format(mass)]
        outdir_mass = os.path.join(outdir, "M%d"%mass)
	thisfilesuffix = filesuffix
	#if not parametricDNN:
	#    thisfilesuffix = filesuffix.replace("dedicatedDNN", "dedicatedDNN%d"%mass)
        infile = os.path.join(inputdir, "%s_SignalM%d.root"%(thisfilesuffix, mass))

	print("filesuffix ", filesuffix," thisfilesuffix ", thisfilesuffix," infile ", infile)
	infile = ROOT.TFile.Open(infile)
        channels = ["MuMu","MuEl","ElEl"]
	if not os.path.exists(outdir_mass):
	    os.system("mkdir "+outdir_mass)
        for channel in channels:
            outfile_name = os.path.join(outdir_mass, prefix_out+"_M%d_%s_autoMC.root"%(mass, channel))
	    outfile  = ROOT.TFile.Open(outfile_name, "recreate")

	    histname_data_obs = "data_obs_{ch}_M{m}".format(ch = channel, m = mass)
	    hist_data_obs = infile.Get(histname_data_obs)
	    hist_clone_data_obs = hist_data_obs.Clone("data_obs")
	    outfile.Write()
	    del hist_clone_data_obs
	    for i, process in enumerate(processlist):
	      if (channel == "MuMu" or channel == "ElEl") and process == "Drell_Yan":
		continue
	      if channel == "MuEl" and ("untagged" in process):
		continue
	      histname_nominal = processnames[i]+"_"+channel
	      hist_nominal = infile.Get(histname_nominal)
	      if "Signal" in process:
		hist_nominal.Scale(1e-3/5.0)
	      hist_clone_nominal = hist_nominal.Clone(processlist[i])
	      outfile.Write()
	      del hist_clone_nominal
	      if "data" in process:
		continue
	      for stat in statslist:
		histname_up = processnames[i]+"_"+channel+"_"+stat+"_up"
		histname_down = processnames[i]+"_"+channel+"_"+stat+"_down"
		histname_sys = process+"_"+stat
		hist_up = infile.Get(histname_up)
		hist_down = infile.Get(histname_down)
		if "Signal" in process:
		  hist_up.Scale(1e-3/5.0)
		  hist_down.Scale(1e-3/5.0)
		if "QCDscale" in histname_sys and "untagged" in histname_sys:
		  histname_sys = process+"_"+stat+process.replace("_untagged", "")
		elif "QCDscale" in histname_sys and "untagged" not in histname_sys:
		  histname_sys = process+"_"+stat+process
		if "CMS_eff_trigger" in histname_sys:
		  histname_sys = process+"_"+stat+"_"+channel
		hist_clone_up = hist_up.Clone(histname_sys+"Up")
		hist_clone_down = hist_down.Clone(histname_sys+"Down")
		outfile.Write()
		del hist_clone_up
		del hist_clone_down
	    outfile.Close()

def writeintoworkspace_1D(masslist,filesuffix, variable, xmin, xmax, inputdir, prefix_out, outdir):
    parametricDNN = True
    if "dedicatedDNN" in filesuffix:
       parametricDNN = False
    #filesuffix = "Hhh_FinalBGYield_xsec1pb_NN_nnout_MTonly"
    #variable = "DNN"
    #xmin = 0.12; xmax = 2.76
    for mass in masslist:
        outdir_mass = os.path.join(outdir, "M%d"%mass)
	thisfilesuffix = filesuffix
	#if not parametricDNN:
	#    thisfilesuffix = filesuffix.replace("dedicatedDNN", "dedicatedDNN%d"%mass)
	print("filesuffix ", filesuffix," thisfilesuffix ", thisfilesuffix)
        infile = os.path.join(inputdir, "%s_SignalM%d.root"%(thisfilesuffix, mass))
        os.system("mkdir -p "+outdir_mass)
        channels = ["MuMu","MuEl","ElEl"]
        for channel in channels:
            outfile = os.path.join(outdir_mass, prefix_out+"_M%d_%s_shapes.root"%(mass, channel))
            command = 'root -b -q "writeworkspace1D.C(%d, \\"%s\\", \\"%s\\", %f, %f,  \\"%s\\", \\"%s\\")" >> workspace1D.log'%(mass, channel, variable, xmin, xmax, infile, outfile)
            print "writeintoworkspace command ",command
            os.system(command)
        ### example
        # root -b -q 'writeworkspace1D.C (400, "MuMu","NN", 0.0, 2.4, "HHbbWW_20180625_NNoutput_MjjCR_test/Hhh_FinalBGYield_xsec1pb_NN_nnout_MTonly_SignalM400.root","test.root")'
        ###

def writeintoworkspace_2D(masslist, filesuffix, variable, xmin, xmax, yvariable, ymin, ymax, inputdir, prefix_out, outdir):
    parametricDNN = True
    if "dedicatedDNN" in filesuffix:
       parametricDNN = False
    #filesuffix = "Hhh_FinalBGYield_xsec1pb_NNvsHME_nnout_MTonly"
    #variable = "DNN"
    #xmin = 0.12; xmax = 2.76
    #yvariable = "HME"
    #ymin = 250; ymax = 1200
    rootws_fname = "root2ws_%s_2D.sh"%filesuffix
    rootws = open(rootws_fname,"write")
    for mass in masslist:
        outdir_mass = os.path.join(outdir, "M%d"%mass)
	thisfilesuffix = filesuffix
	if not parametricDNN:
	    thisfilesuffix = filesuffix+"%d"%mass
        infile = os.path.join(inputdir, "%s_SignalM%d.root"%(thisfilesuffix, mass))
        #infile = os.path.join(inputdir, "%s_SignalM%d.root"%(filesuffix, mass))
        os.system("mkdir -p "+outdir_mass)
        channels = ["MuMu","MuEl","ElEl"]
        for channel in channels:
            outfile = os.path.join(outdir_mass, prefix_out+"_M%d_%s_shapes.root"%(mass, channel))
            command = 'root -b -q "writeworkspace2D.C(%d, \\"%s\\", \\"%s\\", %f, %f, \\"%s\\", %f, %f,  \\"%s\\", \\"%s\\")" >> workspace2D.log'%(mass, channel, variable, xmin, xmax, yvariable, ymin, ymax,  infile, outfile)
            #print "writeintoworkspace command ",command
            #os.system(command)
	    rootws.write(command+"\n")
        ### example
    os.system("chmod 775 "+rootws_fname)
    #os.system("source "+rootws_fname)



def generatedatacard(rootdir, masspoints, prefix, addobservation, autoMC):
    print("generate data card for ", rootdir)
    #channels = ["MuMu","MuEl","ElEl"]
    channels = ["MuMu","ElEl","MuEl"]
    allchannels = ["MuMu","ElEl","MuEl","MuMu_ElEl_MuEl"]
    processes = ["TTbar","SingleTop","Drell_Yan","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV","Signal"]
    datacard = Datacards.Datacards()
    for mass in masspoints:
        channels_rootfile = {}
        dir_mass = os.path.join(rootdir, "M%d"%mass)
        if not os.path.isdir(dir_mass):
            print "wrong directory for root files including shapes ", dir_mass
            sys.exit()
	for channel in allchannels:
            outfile = os.path.join(dir_mass, prefix+"_M{mass}_{channel}.dat".format(mass = mass, channel = channel))
            outfile2 = os.path.join(dir_mass, prefix+"_M{mass}_{channel}_autoMCthresh10.dat".format(mass = mass, channel = channel))
            wsfile = os.path.join(dir_mass, prefix+"_M{mass}_{channel}_combine_workspace.root".format(mass = mass, channel = channel))
            wsfile2 = os.path.join(dir_mass, prefix+"_M{mass}_{channel}_combine_workspace_autoMCthresh10.root".format(mass = mass, channel = channel))
	    #GGToX0ToHHTo2B2L2Nu_M260_ElEl_signalxsec1fb_tminus1
            limitfile  = os.path.join(dir_mass, prefix+"_M{mass}_{channel}_signalxsec1fb_tminus1.log".format(mass = mass, channel = channel))
            limitfile2 = os.path.join(dir_mass, prefix+"_M{mass}_{channel}_signalxsec1fb_tminus1_autoMCthresh10.log".format(mass = mass, channel = channel))
	    #os.system("mv "+outfile +" "+outfile2)
	    #os.system("mv "+ wsfile +" "+wsfile2)
	    #os.system("mv "+ limitfile +" "+limitfile2)

        for channel in channels:
            thisrootfile = os.path.join(dir_mass, prefix+"_M{mass}_{channel}_autoMC.root".format(mass = mass, channel = channel))
            channels_rootfile[channel] = thisrootfile
            if not os.path.isfile(thisrootfile):
                print "rootfile does not exist ", thisrootfile
                sys.exit()
            outfile = os.path.join(dir_mass, prefix+"_M{mass}_{channel}.dat".format(mass = mass, channel = channel))

            #print "channel ",channel," outfile ",outfile
            datacard.generateDatacard( outfile, thisrootfile, channel, mass, addobservation, autoMC)
        outfile_all = os.path.join(dir_mass, prefix+"_M{mass}_MuMu_ElEl_MuEl.dat".format(mass = mass))
        datacard.generateDatacard_multichannels(outfile_all, channels_rootfile, channels, mass, addobservation, autoMC)
        #outfile_MuMu_MuEl = os.path.join(dir_mass, prefix+"_M{mass}_MuMu_MuEl.dat".format(mass = mass))
        #datacard.generateDatacard_multichannels(outfile_MuMu_MuEl, channels_rootfile, ["MuMu","MuEl"], mass, addobservation)

        



def alldatacards(rootdir, model, masspoints, addobservation, outdir, analysistype):
    #os.system("cp %s %s"%(rootfile, outdir))
    channels = ["MuMu","MuEl","ElEl"]
    processes = ["signal","TT","DY","sT"]
    for mass in masspoints:
        outfile = os.path.join(outdir, "%s_M%d.txt"%(analysistype, mass))
        suffix = "M%d"%mass
        i = 0
        #if mass<=400:
        #    i = 11
        #suffix = "RadionM%d_WP%i"%(mass, i) 
        rootfile = os.path.join(rootdir, "HME_%s_2D_bin10_M%d_workspace.root"%(model, mass))
        generatedatacard(outfile, rootfile, channels, processes, mass, suffix, addobservation, analysistype)



def alldatacards_HME(rootfile, masspoints, addobservation, outdir, analysistype):
    os.system("cp %s %s"%(rootfile, outdir))
    channels = ["MuMu","MuEl","ElEl"]
    processes = ["signal","TT","DY","sT"]
    #for mass in masspoints:
    for mass in masspoints:
      #for i in range(0, 12):
        i = 11
        if mass>400:
            i = 0
        suffix = "RadionM%d_WP%d"%(mass,i)
        outfile = os.path.join(outdir, "%s_%s.txt"%(analysistype, suffix))
        generatedatacard(outfile, rootfile, channels, processes, mass, suffix, addobservation, analysistype)


def extranumber(s):
    num = []
    for t in s.split():
        try:
            num.append(float(t))
        except ValueError:
            pass

    return num


def submitdatacards(cardname, method, mass):
    logfile = "expectedlimits_tmp.log"
    os.system("combine {cardname}  -M {method} -m {mass} -t -1 -s -1 > {logfile}".format(cardname = cardname, method = method, mass = mass, logfile = logfile))
    #os.system("combine {cardname}  -M {method} -t 100 -s -1 > {logfile}".format(cardname = cardname, method = method, logfile = logfile))
    logopen = open(logfile, "read")
    sigma_xsec = 1.0
    percents = [2.5, 16.0, 50.0, 84.0, 97.5]
    limits = {}
    for line in logopen:
        if not line.startswith("Expected"):
            continue
        #print "line ",line
        for percent in percents:
            pattern = "Expected%5.1f%%: r < "%(percent)
            #print " percent ",percent, " pattern ",pattern
            if re.match(pattern, line):
                limits[percent] = float(line.split(" ")[-1].rstrip()) * sigma_xsec
    #print "limits ",limits
    return limits



def submitdatacards_v2(cardname, method):
    logfile = "expectedlimits_tmp.log"
    os.system("combine {cardname}  -M {method} -t 100  -m {mass} -s -1 > {logfile}".format(cardname = cardname, method = method, mass = mass, logfile = logfile))
    logopen = open(logfile, "read")
    sigma_xsec = 1.0
    #percents = [2.5, 16.0, 50.0, 84.0, 97.5]
    limits = {}
    for line in logopen:
        if line.startswith ("median expected limit: "):
            nums = extranumber(line)
            limits[50.0]  = nums[0] * sigma_xsec
        elif line.startswith("   68% expected band :"):
            nums = extranumber(line)
            limits[16.0] = nums[0] * sigma_xsec
            limits[84.0] = nums[1] * sigma_xsec
        elif line.startswith("   95% expected band :"):
            nums = extranumber(line)
            limits[2.50] = nums[0] * sigma_xsec
            limits[97.50] = nums[1] * sigma_xsec
        else :
            pass
              
    print "limits ",limits
    return limits

def submitdatacards_v3(cardname, mass, method):
    logfile = "expectedlimits_tmp.log"
    outdir = "shapeM%d_%s"%(mass, method)
    os.system("mkdir -p "+outdir)
    os.system("combine {cardname}  -M {method} -m {mass} -t -1 -s -1 --plots --out {outdir} > {logfile}".format(cardname = cardname, method = method, mass = mass, outdir = outdir, logfile = logfile))

def runCombineTools(cardnameprefix, masspoints, method, logfile):
    for mass in masspoints:
        cardname = cardnameprefix + "_M%d.txt"%mass
        os.system("combine {cardname}  -M {method} -t 10000 -s -1 >> {logfile}".format(cardname = cardname, method = method, logfile = logfile))

def runImpacts(cardnameprefix, xpoints, logfile):
    for x in xpoints:
        cardname = cardnameprefix + "_M%d.txt"%x
        wsname = cardnameprefix + "_M%d.root"%x
        samplename = cardnameprefix + "_M%d"%x
        impacfile = samplename+"_impacts"
        os.system("text2workspace.py {cardname} -o {wsname} >> {logfile}".format(cardname = cardname, wsname = wsname, logfile = logfile))
        os.system("combineTool.py -M Impacts -d {wsname} -m 125 --doInitialFit >> {logfile}".format(wsname = wsname, logfile = logfile))
        os.system("combineTool.py -M Impacts -d {wsname} -m 125 --doFits --parallel 4 >> {logfile}".format(wsname = wsname, logfile = logfile))
        os.system("combineTool.py -M Impacts -d {wsname} -m 125 -o {impactfile}.json >> {logfile}".format(wsname = wsname, impactfile = impacfile,logfile = logfile))
        os.system("plotImpacts.py -i {impactfile}.json -o {impactfile} >> {logfile}".format(impactfile = impacfile, logfile = logfile))

def getlimits(logfile, mass):
    return True
    
def CombineLimitplots(filelist, histlist, masspoints, xtitle, legends, text, plotname):

    #colors = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kMagenta+2, ROOT.kOrange, ROOT.kBlack]
    #markers = [20, 20,  21,22,23,24]
    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kMagenta+2, ROOT.kCyan+2, ROOT.kBlack]
    markers = [20,  21,22,23,24]
    tfilelist = []
    for f in filelist:
        tfilelist.append( ROOT.TFile(f, "READ"))


    c1 = ROOT.TCanvas("c1","c1",700, 800)
    #c1 = ROOT.TCanvas("c1","c1",600, 800)

    c1.SetLogy()
    c1.SetGridx()  
    c1.SetGridy()  
    c1.SetTickx()  
    c1.SetTicky()  
    minx = min(masspoints)*0.8
    maxx = max(masspoints)*1.3
    b1 = ROOT.TH1F("b2","b2", len(masspoints)*2, minx, maxx)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
    b1.SetTitle(" ")
    b1.GetXaxis().SetTitle(xtitle)
    b1.GetYaxis().SetTitle("95% CL limit on #sigma(gg#rightarrow X^{spin 0} #rightarrow HH) #times #it{Br}(HH #rightarrow b#bar{b}l#nul#nu) [fb]")
    b1.GetYaxis().SetTitleOffset(1.4)
    #yhigh = max(twosigma_up)*1.2
    #ylow = min(twosigma_low)*.9
    yhigh = 2000.0 #2000.0
    ylow = 1.0#1.0
    b1.GetYaxis().SetRangeUser(ylow, yhigh)
    b1.SetStats(0)
    b1.Draw()


    leg0 = ROOT.TLegend(0.65,0.8,0.9,0.90)
    leg0.SetFillColor(ROOT.kWhite)
    leg0.SetTextFont(42)
    #le0g.AddEntry(grdata,"observed","pl")

    #leg = ROOT.TLegend(0.52,0.48,0.9,0.52+0.05*len(tfilelist))
    leg = ROOT.TLegend(0.56,0.48,0.9,0.5+0.05*len(tfilelist))
    #leg.SetHeader("DNN inputs: Kinematic+")
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextFont(42)
    leg.SetTextSize(.028)
    #leg.AddEntry(grdata,"observed","pl")
    for i, tf in enumerate(tfilelist):
        tf.cd()
        gcentralname = histlist[i]+"_central"
        if "data" in histlist[i]:
            gcentralname = histlist[i]
        g_central = tf.Get( gcentralname )
        g_central.SetLineColor(colors[i])
        g_central.SetMarkerColor(colors[i])
        g_central.SetMarkerStyle(markers[i])
        if i == 0:
            g_onesigma = tf.Get(histlist[i]+"_onesigma")
            g_twosigma = tf.Get(histlist[i]+"_twosigma")
            g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
            g_twosigma.SetLineStyle(2)
            g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
            g_onesigma.SetLineStyle(2)
            leg0.AddEntry(g_central,"Expected 95% upper limit","l")
            leg0.AddEntry(g_onesigma,"1 std. deviation, Expected","f")
            leg0.AddEntry(g_twosigma,"2 std. deviation, Expected","f")
            g_twosigma.Draw("fe3same")
            g_onesigma.Draw("fe3same")

        #if i ==0:
        #    g_central.Draw("lsame")
        #    thisleg = leg.AddEntry(g_central, legends[i],"l")
        #else:
        g_central.Draw("lpsame")
        thisleg = leg.AddEntry(g_central, legends[i],"pl")
        thisleg.SetTextColor(colors[i])
    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+" "+"35.87 fb^{-1}(13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.2,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    #tex1.Draw("same")
    
    leg0.Draw("same")
    leg.Draw("same")
    c1.SaveAs(plotname+"_95Upperlmit_Nodata_combined.pdf")
    c1.SaveAs(plotname+"_95Upperlmit_Nodata_combined.C")
    img = ROOT.TImage.Create();
    img.FromPad(c1)
    img.WriteImage(plotname+"_95Upperlmit_Nodata_combined.png")
    
def makeBrazilPlot(masspoints, central, onesigma_up, onesigma_low, twosigma_up, twosigma_low, xtitle, text, plotname):
    fakeerrors = [0.0]*len(masspoints)
    #g_onesigma = TGraphAsymmErrors(len(masspoints),  np.array(masspoints),  np.array(central), np.array(fakeerrors), np.array(fakeerrors), np.array(onesigma_low), np.array(onesigma_up))
    #g_twosigma = TGraphAsymmErrors(len(masspoints),  np.array(masspoints),  np.array(central), np.array(fakeerrors), np.array(fakeerrors), np.array(twosigma_low), np.array(twosigma_up))
    c1 = ROOT.TCanvas("c1","c1",600, 800)
    outfilename = plotname+".root"
    tfile = ROOT.TFile(outfilename,"UPDATE")
    
    c1.SetLogy()
    c1.SetGridx()  
    c1.SetGridy()  
    c1.SetTickx()  
    c1.SetTicky()  
    onesigma_low.reverse()
    twosigma_low.reverse()
    onesigma_all = onesigma_up + onesigma_low
    twosigma_all = twosigma_up + twosigma_low
    masspoints_all = masspoints + list(reversed(masspoints))
    masspoints_f =  np.array(masspoints)+0.0
    print "allXpoints ",masspoints_all," onesigma ",onesigma_all," float masspoints ",masspoints_f
    g_central = ROOT.TGraph(len(masspoints),  np.array(masspoints)+0.0,  np.array(central))
    g_onesigma = ROOT.TGraph(len(masspoints)*2,  np.array(masspoints_all)+0.0,  np.array(onesigma_all))
    g_twosigma = ROOT.TGraph(len(masspoints)*2,  np.array(masspoints_all)+0.0,  np.array(twosigma_all))


    #g_twosigma.SetFillColor(ROOT.kYellow)
    g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
    g_twosigma.SetLineStyle(2)
    g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
    g_onesigma.SetLineStyle(2)
    g_central.SetLineWidth(2)
    g_central.SetLineStyle(7)
    g_central.SetLineColor(9)
    g_central.SetMarkerStyle(20)
    g_central.SetMarkerSize(1)
    g_central.SetMarkerColor(9)

    #b1 = ROOT.TH1F("b2","b2",14, 250.0, 950.0)
    minx = min(masspoints)*0.8
    maxx = max(masspoints)*1.3
    b1 = ROOT.TH1F("b2","b2", len(masspoints)*2, minx, maxx)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
    b1.SetTitle(" ")
    b1.GetYaxis().SetTitle("95% C.L. limit on #sigma(gg#rightarrow X #rightarrow HH) X #bf{#it{Br}}(HH #rightarrow b#bar{b}l#nul#nu) (fb)")
    b1.GetXaxis().SetTitle(xtitle)
    #yhigh = max(twosigma_up)*1.2
    #ylow = min(twosigma_low)*.9
    yhigh = 2000.0 #10000.0
    ylow = 1.0#10.0
    b1.GetYaxis().SetRangeUser(ylow, yhigh)
    b1.SetStats(0)

    b1.Draw()
    g_twosigma.Draw("fe3same")
    g_onesigma.Draw("fe3same")
    g_central.Draw("lpsame")
    samplename =  plotname.split("/")[-1]
    g_central.SetName("%s_central"%samplename)
    g_onesigma.SetName("%s_onesigma"%samplename)
    g_twosigma.SetName("%s_twosigma"%samplename)
    g_central.Write()
    g_onesigma.Write()
    g_twosigma.Write()

    
    leg = ROOT.TLegend(0.65,0.8,0.9,0.90)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextFont(42)
    #leg.AddEntry(grdata,"observed","pl")
    leg.AddEntry(g_central,"Expected 95% upper limit","l")
    leg.AddEntry(g_onesigma,"1 std. deviation, Expected","f")
    leg.AddEntry(g_twosigma,"2 std. deviation, Expected","f")
    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+" "+"35.87 fb^{-1}(13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.2,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    tex1.Draw("same")
    
    leg.Draw("same")
    c1.SaveAs(plotname+"_95Upperlmit_Nodata_tmp.pdf")
    c1.SaveAs(plotname+"_95Upperlmit_Nodata_tmp.C")
    img = ROOT.TImage.Create();
    img.FromPad(c1)
    img.WriteImage(plotname+"_95Upperlmit_Nodata_tmp.png")
    #c1.SaveAs(plotname+"_95Upperlmit_Nodata_tmp.png")

def Limitplots(cardnameprefix, masspoints, method, text, plotname):
    onesigma_up = []
    twosigma_up = []
    onesigma_low = []
    twosigma_low = []
    central = []
    method = "Asymptotic"
    """
    for mass in masspoints:
        cardname = cardnameprefix + "_M%d.txt"%mass
        submitdatacards_v3(cardname, mass, 'MaxLikelihoodFit')
    """ 
    for mass in masspoints:
        cardname = cardnameprefix + "_M%d.txt"%mass
        limits = submitdatacards(cardname, method, mass)
        print "cardname ",cardname," limits ",limits,"\n"
        central.append(limits[50.0])
        twosigma_low.append(limits[2.5])
        onesigma_low.append(limits[16.0])
        onesigma_up.append(limits[84.0])
        twosigma_up.append(limits[97.5])

    xtitle = "Mass [GeV]"
    makeBrazilPlot(masspoints, central, onesigma_up, onesigma_low, twosigma_up, twosigma_low, xtitle, text, plotname)

def GetNLimits(cardnameprefix, mass, Ntimes, method, text, plotname):
    method = "Asymptotic"
    cardname = cardnameprefix + "_M%d.txt"%mass
    limits = submitdatacards(cardname, method, mass)
    hist = ROOT.TH1F("expectedlimits_M%d"%mass, "expectedlimits_M%d"%mass,100, limits[2.5], limits[97.5])
    for i in range(Ntimes):
        limits = submitdatacards_v2(cardname, method, mass)
        print "cardname ",cardname," limits ",limits,"\n"
        hist.Fill(limits[50.0])
        os.system("sleep 1")
    c1 = ROOT.TCanvas("c1","c1",800, 600)
    hist.Draw("hist")
    hist.SetLineWidth(2)
    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+" "+"35.87 fb^{-1}(13 TeV),2016")
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.2,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    tex1.Draw("same")
    
    c1.SaveAs(plotname+"_95Upperlmit_centralexpected_Nodata_tmp.pdf")
    c1.SaveAs(plotname+"_95Upperlmit_centralexpected_Nodata_tmp.C")
    


def Limitplots_HME(cardnameprefix, xpoints, mass,  method, text, plotname):
    onesigma_up = []
    twosigma_up = []
    onesigma_low = []
    twosigma_low = []
    central = []
    method = "Asymptotic"
    """
    for mass in masspoints:
        cardname = cardnameprefix + "_M%d.txt"%mass
        submitdatacards_v3(cardname, mass, 'MaxLikelihoodFit')
    """ 
    for x in xpoints:
        suffix = "RadionM%d_WP%d"%(mass, x)
        cardname = cardnameprefix + "_%s.txt"%(suffix)
        limits = submitdatacards(cardname, method, mass)
        print "cardname ",cardname," limits ",limits,"\n"
        central.append(limits[50.0])
        twosigma_low.append(limits[2.5])
        onesigma_low.append(limits[16.0])
        onesigma_up.append(limits[84.0])
        twosigma_up.append(limits[97.5])

    xtitle = "WorkingPoint (NN cut)"
    makeBrazilPlot(xpoints, central, onesigma_up, onesigma_low, twosigma_up, twosigma_low, xtitle, text, plotname)


def parseDatacards_NNcut1D(inputdir, nnout, outdir_prefix):
    parametricDNN = True
    autoMC = True
    if "dedicatedDNN" in nnout:
       parametricDNN = False
    masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
    #masslist = [400]
    NNout = "DNN";
    NNcutlist = [0.0, 0.04,  0.12,  0.20,  0.28, 0.36, 0.40,  0.48,  0.56, 0.60,  0.72]
    NNcutlist = [0.0]
    NoMjjbins =False
    if NoMjjbins:
       outdir_prefix = outdir_prefix+"_MjjcutSonly"
       #outdir_prefix = outdir_prefix+"_MjjMerged"
    #nnout = "nnout_MTonly"
    for nncut in NNcutlist:
        
	#nncutsuffix = "nnstep0p04_nncut0p%s"%(str(nncut)[2:])
	nncutsuffix = "nnstep0p05_nncut0p%s"%(str(nncut)[2:])
	#nncutsuffix = "nnstep0p033_nncut0p%s"%(str(nncut)[2:])
	#nncutsuffix = "nnstep0p025_nncut0p%s"%(str(nncut)[2:])
	#nncutsuffix = "nnstep0p1_nncut0p%s"%(str(nncut)[2:])


	#nncutsuffix = "nnstep0p1_nncut0p%s_HME1p1"%(str(nncut)[2:])
	nncutsuffix = "nnstep0p1_nncut0p%s_HME1p15"%(str(nncut)[2:])
	#nncutsuffix = "nnstep0p1_nncut0p%s_HME1p2"%(str(nncut)[2:])
	#nncutsuffix = "nnstep0p033_nncut0p%s_HME1p15"%(str(nncut)[2:])
	#nncutsuffix = "nnstep0p2_nncut0p%s_HME1p15"%(str(nncut)[2:])
	outdir = outdir_prefix+"_"+nnout+"_"+nncutsuffix+"_limits/"
	os.system("mkdir -p "+outdir)
	NNoutMin = nncut;  NNoutMax = 3.0-nncut*2
	if NoMjjbins:
	    NNoutMax = 1.0-nncut*2 ## no Mjj macro bin
	## 1D
	#filesuffix = "Hhh_FinalBGYield_xsec1pb_NN_%s_%s"%( nnout, nncutsuffix)
	## 2D
	filesuffix = "Hhh_FinalBGYield_xsec1pb_NNvsHME_%s_%s"%( nnout, nncutsuffix)
	
	#if "dedicated" in nnout:
	#    mass = int(nnout[-3:])
	#    masslist = [mass]
	
	#writeintoworkspace_1D(masslist, filesuffix, NNout, NNoutMin, NNoutMax, inputdir, prefix_out, outdir)
        createShapes_1D(masslist, filesuffix, inputdir, prefix_out, outdir)

    #writeintoworkspace_2D(masslist,inputdir, prefix_out, outdir)
	generatedatacard(outdir, masslist, prefix_out, True, autoMC)

def parseDatacards_HME1D(inputdir, nnout, outdir_prefix):
    autoMC = True
    filesuffix = "Hhh_FinalBGYield_xsec1pb_HME_%s"%(nnout)
    outdir = outdir_prefix+"_"+nnout+"_limits/"
    os.system("mkdir -p "+outdir)
    masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
    createShapes_1D(masslist, filesuffix, inputdir, prefix_out, outdir)
    generatedatacard(outdir, masslist, prefix_out, True, autoMC)



def parseDatacards_NNcut2D(inputdir, nnout, outdir_prefix):
    masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
    NNout = "DNN";
    NNcutlist = [0.0, 0.04,  0.12,  0.20,  0.28, 0.36, 0.40,  0.48,  0.56, 0.60,  0.72]
    NNcutlist = [0.0, 0.04,  0.12,  0.20]
    NNcutlist = [0.0]

    HMEout = "HME"; HMEMin = 250.0; HMEMax=1200.0
    #nnout = "nnout_MTonly"
    for nncut in NNcutlist:
	nncutsuffix = "nnstep0p04_nncut0p%s"%(str(nncut)[2:])
	outdir = outdir_prefix+"_NNvsHME_"+nnout+"_"+nncutsuffix+"_HMEbinsv4_limits2D/"
	#nncutsuffix = "nnstep0p05_nncut0p%s"%(str(nncut)[2:])
	#nncutsuffix = "nnstep0p033_nncut0p%s"%(str(nncut)[2:])
	#nncutsuffix = "nnstep0p025_nncut0p%s"%(str(nncut)[2:])
	nncutsuffix = "nnstep0p1_nncut0p%s"%(str(nncut)[2:])
	nncutsuffix = "nnstep0p02_nncut0p%s"%(str(nncut)[2:])
	nncutsuffix = "nnstep0p2_nncut0p%s"%(str(nncut)[2:])
	nncutsuffix = "nnstep0p01_nncut0p%s"%(str(nncut)[2:])
	outdir = outdir_prefix+"_NNvsHME_"+nnout+"_"+nncutsuffix+"_HMEbinsv4_1p15_limits2D/"
	#outdir = outdir_prefix+"_NNvsHME_"+nnout+"_"+nncutsuffix+"_HMEbinsv4_1p1_limits2D/"
	#outdir = outdir_prefix+"_NNvsHME_"+nnout+"_"+nncutsuffix+"_HMEbinsv4_1p05_limits2D/"

	os.system("mkdir -p "+outdir)
	NNoutMin = nncut;  NNoutMax = 3.0-nncut*2
	##HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D/Hhh_FinalBGYield_xsec1pb_NNvsHME_nnout_MTonly_nnstep0p04_nncut0p72_
	filesuffix = "Hhh_FinalBGYield_xsec1pb_NNvsHME_%s_%s"%( nnout, nncutsuffix)
	#writeintoworkspace_1D(masslist, filesuffix, NNout, NNoutMin, NNoutMax, inputdir, prefix_out, outdir)

	#writeintoworkspace_2D(masslist,filesuffix, NNout, NNoutMin, NNoutMax, HMEout, HMEMin, HMEMax, inputdir, prefix_out, outdir)
	#generatedatacard(outdir, masslist, prefix_out, True)

prefix_out = "GGToX0ToHHTo2B2L2Nu"
suffix = 'Radion_1D_nnoutMTandMJJ_xsec1pb'
output_suffix = '{:%Y-%m-%d}_{}'.format(datetime.date.today(), suffix)
#datacarddir = "Datacards/"

#os.system("mkdir -p "+datacarddir)
inputdir = "HHbbWW_20180625_NNoutput_MjjCR_test/"
inputdir = "HHbbWW_20200331_NNoutput_MjjCR/"
inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy1D/"
#inputdir = "HHbbWW_20200608_NNoutput_MjjCR_NNcutstudy1D_dedicatedDNN/"
#inputdir = "HHbbWW_20200608_NNoutput_NOMjjCR_NNcutstudy1D_dedicatedDNN/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy1D_binsize005/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy1D_binsize0p033/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy1D_binsize0p025/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy1D_binsize0p10/"
inputdir = "HHbbWW_20200807_NNoutput_MjjCR_NNcutstudy1D/"
#inputdir = "HHbbWW_20200623_NNoutput_NoMjjCR_NNcutstudy1D/"
#inputdir = "HHbbWW_20200623_NNoutput_NoMjjCRMjjcut_NNcutstudy1D/"
inputdir = "HHbbWW_20200807_test_MjjSignalonly/"
inputdir = "HHbbWW_20200812_NNoutput_MergedMjj_NNcutstudy1D_MCstat/"
#inputdir = "HHbbWW_20200812_NNoutput_NoMjjCR_NNcutstudy1D_MCstat/"
inputdir = "HHbbWW_20200812_NNoutput_MjjCR_NNcutstudy1D_MCstat/"

#inputdir = "HHbbWW_20200814_NNvsHME_linearized1D_MjjCR_NNoutp1_HMEbinsv4_1p15/"
#inputdir = "HHbbWW_20200814_NNvsHME_linearized1D_MjjcutSonly_NNoutp05_HMEbinsv4_1p15/"
#inputdir = "HHbbWW_20200814_NNvsHME_linearized1D_MjjMerged_NNoutp05_HMEbinsv4_1p15/"
inputdir = "HHbbWW_20200814_NNvsHME_linearized1D_MjjCR_HMEbinsv4_1p15/"
#inputdir = "HHbbWW_20200814_NNvsHME_linearized1D_MjjcutSonly_HMEbinsv4_1p15/"
#inputdir = "HHbbWW_20200814_NNvsHME_linearized1D_MjjMerged_NNoutp05_HMEbinsv4_1p15/"
#inputdir = "HHbbWW_20200814_NNvsHME_linearized1D_MjjMerged_HMEbinsv4_1p15/"

inputdir = "HHbbWW_20200814_NNvsHME_linearized1D_MjjCR_HMEbinsv4/"
inputdir = "HHbbWW_20201006_NNvsHME_linearized1D_MjjCR_HMEbinsv4/"


outdir = prefix_out+inputdir
method = "Asymptotic"
masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
#masslist = [260, 270, 300, 350, 400, 450, 500]
#masslist = [400]
nnout = "nnout_MTonly"
prefix_out = prefix_out+"_NNvsHME_MCstat"
#parseDatacards_NNcut1D(inputdir, nnout, prefix_out)
#parseDatacards_NNcut1D(inputdir, "nnout_MTandMT2", prefix_out)
parseDatacards_NNcut1D(inputdir, "nnout_MTandMT2_MJJ", prefix_out)
### HME, anti DNN cut
#### HME, anti DNN cut
inputdir = "HHbbWW_20200821_HME_HMEbinsv4_1p05/"
prefix_out = prefix_out + "_HME"
#parseDatacards_HME1D(inputdir, "nnout_MTandMT2cut0p8", prefix_out)
#parseDatacards_HME1D(inputdir, "nnout_MTandMT2_MJJcut0p8", prefix_out)
#parseDatacards_HME1D(inputdir, "nnout_MTandMT2_MJJcut0p95", prefix_out)
#parseDatacards_NNcut1D(inputdir, "nnout_MTandMT2_HMEMJJ_dedicatedDNN", prefix_out)
#parseDatacards_NNcut1D(inputdir, "nnout_MTandMT2_HME_dedicatedDNN", prefix_out)
#parseDatacards_NNcut1D(inputdir, "nnout_MTandMT2_HMEMJJ_dedicatedDNN750", prefix_out)
#parseDatacards_NNcut1D(inputdir, "nnout_MTandMT2_HME_dedicatedDNN750", prefix_out)
inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D/"
inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv2/"
inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv3/"
inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p1/"
inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p15/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p2/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p15_NNout0p05/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p15_NNout0p033/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p1_NNout0p025/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p15_NNout0p025/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p05/"
inputdir  = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p15_NNout0p1/"
#inputdir  = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p1_NNout0p02/"
inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p15_NNout0p02/"
inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p15_NNout0p2/"
#inputdir = "HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy2D_HMEbinsv4_1p15_NNout0p01/"
#parseDatacards_NNcut2D(inputdir, nnout, prefix_out)
#parseDatacards_NNcut2D(inputdir, "nnout_MTandMT2", prefix_out)
#parseDatacards_NNcut2D(inputdir, "nnout_MTandMT2_MJJ", prefix_out)
#masslist = [260, 270, 400, 750, 900]
#masslist = [400]
filesuffix = "Hhh_FinalBGYield_xsec1pb_NN_nnout_MTonly"
#alldatacards(rootfile, masslist, False, outdir)
cardnameprefix = os.path.join(outdir, "shape")
plotname = os.path.join(outdir, "Radion")
#masslist = [260, 500, 650,  900]
#Limitplots(cardnameprefix, masslist, method, plotname)
#modellist = ['MTonly','MT2only','MTandMT2','MTandMT2_MJJ','MTandMT2_HME','MTandMT2_HMEMJJ', 'MTandMJJ']
#modellist = ["MTonly","MTandMT2", "MTandMT2_MJJ"]
modellist = ["MTonly","MTandMT2"]
#modellist = ["MTonly"]
