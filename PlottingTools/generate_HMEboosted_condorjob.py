import os
import sys 
sys.argv.append( '-b' )


#masslist = [1000, 1250, 1500, 1750, 2000, 2500, 3000, 700, 800, 900, 750, 850, 650, 600]
masslist = [250, 260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650, 700,750, 800, 850, 900, 1000] 
#masslist = [2000]

iterations = 10000 ## default=10k, opt: 10k, 50k, 100k
nbins_rebin = 1 ## default=10, opt:1,10,20,50
metRes = 25.2##default, opt: 50
gravitypercent= 0.18 #0.18default, opt: 0.18, 0.3, 0.4
njobs = 80
bashscript = "runHME_boosted.sh"

def getEntries(filename, treename="Events"):
    import ROOT
    tfile = ROOT.TFile(filename)
    tree = tfile.Get(treename)
    return tree.GetEntries()

def generate_run_HME(masslist, workdir, inputdir, outdir, njobs):
    if not os.path.exists(workdir):
	os.system("mkdir "+workdir)
	os.system("mkdir "+workdir+"/output")
	os.system("mkdir "+workdir+"/error")
	os.system("mkdir "+workdir+"/log")
    if not os.path.exists(outdir):
	os.system("mkdir "+outdir)
    os.system("cp "+bashscript + " " +workdir)
    fname_all = workdir+"condor_RunHME_allMass_boosted.sh"
    script_all = open(fname_all, "write")
    script_all.write("#!/bin/bash\n")
    #script_all.write("eval `scramv1 runtime -sh`\n")
    for mass in masslist:
       inputfile = inputdir + "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d.root"%mass
       nTotal = getEntries(inputfile)
       nEv_per_job = nTotal/njobs + 1
       for ijob in range(njobs):
	   nStart = nEv_per_job*ijob
	   nEnd = nEv_per_job*(ijob+1)
	   if nEnd > nTotal: nEnd = nTotal
	   outputfile = outdir + "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_HME_ijob%d.root"%(mass, ijob)
	   batchfname = workdir + "Batch_GluGluToRadionToHHTo2B2VTo2L2Nu_HME_M%d_ijob%d.cmd"%(mass, ijob)
	   print("Mass ", mass, " totalevt ", nTotal," nstart ", nStart, " nEnd ", nEnd, " out ", outputfile)
	   script_all.write("condor_submit "+batchfname + "\n")
	   batchscript = open(batchfname, "write")
	   batchscript.write("""universe              = vanilla 
executable            = {script}
arguments             = {arg1} {arg2} {arg3} {arg4} {arg5} {arg6}
output                = output/{suffix}.$(ClusterId).$(ProcId).out
error                 = error/{suffix}.$(ClusterId).$(ProcId).err
log                   = log/{suffix}.$(ClusterId).log
request_memory        = 4000M                                                                                                                        
+JobFlavour           = "tomorrow"
Notification          = Complete
notify_user           = taohuang@email.tamu.edu
queue
       """.format(script = bashscript, suffix = "GluGluToRadionToHHTo2B2VTo2L2NuM%dijob%d"%(mass, ijob), arg1=inputfile, arg2=outputfile, arg3=nStart, arg4=nEnd, arg5=iterations, arg6=gravitypercent))
    os.system("chmod 775 "+fname_all)

inputdir = "/eos/user/t/tahuang/Florian_2016_ForReHME/"
newfolder = "ReRunHME_Florian_2016_20220720_it10000_test/"
outdir = "/eos/user/t/tahuang/"+newfolder
workdir = "/afs/cern.ch/work/t/tahuang/HHAnalysis/CMSSW_10_2_0/src/HhhAnalysis/HeavyMassEstimator/test/"+newfolder


generate_run_HME(masslist, workdir, inputdir, outdir, njobs)
