#!/usr/bin/python
import glob
import os
import sys
sys.path.append('/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD/')
import ROOT
from localSamplelist import *

def getTotalEvents(filename, treename= "Friends"):
    h1 = ROOT.TH1F("h1","h1",10,0,10)
    ch = ROOT.TChain(treename)
    ch.Add(filename)
    ch.Draw("1>>h1","ll_M<76")
    return h1.GetEntries()
#benchmarks =  ["ttV","Wjet","sT","VV","DY","TT","Data"]
benchmarks =  ["ttV","Wjet","sT","VV","DY","TT","Data"]
#benchmarks = ["ttV","Wjet","sT","VV","DY","TT","Data","radion","graviton"]
#benchmarks = ["radion_M400","graviton_M400"]
signals = []
masses = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650,750, 800, 900]
for mass in masses:
	signals.append("RadionM%d"%mass)
	#benchmarks.append("graviton_M%d"%mass)
#alliters = [100, 1000, 10000, 100000, 1000000]
benchmarks = signals
#benchmarks = signals + benchmark
#benchmarks.append("radion_M750")
#benchmarks = ["TT"]
#benchmarks.append("sT_top")
#benchmarks.append("sT_antitop")
#benchmarks.append("DYM10to50")#69k
#benchmarks.append("DYToLL0J")#838k
#benchmarks.append("DYToLL1J")#5095k
#benchmarks.append("DYToLL2J")#15265k
#
DYestimation =False

alliters = [10000]
jobdir = "20220717_HHbbWW_newHME_10k_v2"
if DYestimation:
   benchmarks = ["sT","TT","Data"]
   jobdir = jobdir+"_DYestimation"
#jobdir = "20180205_10k_ALLSamples_Louvain"### write  all samples to same dir
submitscript_name = "Singalonly"
outputdir = "/eos/user/t/tahuang/%s/"%(jobdir)
totaljobs = 40
os.system("mkdir -p %s" % jobdir)
os.system("mkdir -p %s" % outputdir)
CMSSW_BASE = os.environ['CMSSW_BASE']
# kepperror > /dev/null
#drop error &> /dev/null
def generateJobs():
    submitscript = open("submitall%s_%s.sh"%(jobdir, submitscript_name),"w")
    for job in benchmarks:
      if DYestimation:
         job = job + untagged_suffix
      for dataname in full_local_samplelist[job].keys():
          if DYestimation:
             datname = dataname + untagged_suffix
	  for iters in alliters:
	    totaljobs = 10
	    filename = full_local_samplelist[job][dataname]['path']
	    #events_per_job = 1000
	    #if not(os.path.isfile(filename)):
	    #    print "dataname ",dataname," filename ",filename," SKIP!! "
	    #    continue
	    #statinfo = os.stat(filename)
	    ##totaljobs = int(statinfo.st_size*1.0/1e6)+1
	    #totalevents = getTotalEvents(filename, "Friends")
	    #totaljobs = int(totalevents*1.0/events_per_job)+1
	    #if "Radion" in job or job == "Data":
	    #    totaljobs = totaljobs *10
	    #print("job ",job," dataname ",dataname," totaljobs ",totaljobs," file size (M) %.2f"%(statinfo.st_size*1.0/1e6)," totevents ",totalevents," file ",filename)
	    #for ijob in range(0, totaljobs):
	    for ijob in range(0, totaljobs):
		jobscript = open("{0}/Send_PlotterProducer_{1}_{2}_{3}.cmd".format(jobdir, dataname, ijob, iters), "w")
		jobscript.write("""universe              = vanilla 
executable            = runhme.sh
arguments             = {njobs} {ijob} {jobtype} {outfolder}
output                = output/newhme.$(ClusterId).$(ProcId).out
error                 = error/newhme.$(ClusterId).$(ProcId).err
log                   = log/newhme.$(ClusterId).log
request_cpus          = 2
request_memory        = 4000M
+JobFlavour           = "workday"
Notification          = Complete
queue
""".format(jobtype=job, dataname=dataname, njobs=totaljobs, ijob=ijob, jobdir = jobdir, outfolder=outputdir))
		jobscript.close()

		#os.system("cat %s/Send_PlotterProducer_%s_%d_%d.cmd"%(jobdir, job, ijob, iters))
		submitscript.write("condor_submit {0}/Send_PlotterProducer_{1}_{2}_{3}.cmd\n".format(jobdir, dataname, ijob, iters))
    submitscript.close()
    os.system("chmod +x submitall%s_%s.sh"%(jobdir, submitscript_name))

def checkoutput():
    resubmitscript = open("submitall%s_%s_resubmit.sh"%(jobdir, submitscript_name),"w")
    resubmitscript.write("""#!/bin/bash
cd $CMSSW_BASE/src/HhhAnalysis/PlottingTools/
	""")
    for job in benchmarks:
      if job == "TT":
          continue
      for dataname in full_local_samplelist[job].keys():
	  totaljobs = 40
	  filename = full_local_samplelist[job][dataname]['path']
	  events_per_job = 1000
	  if DYestimation:
	      events_per_job = 2000
	      if dataname == "MuonEG":
	          continue
	      filename = untagged_samplelist[job][dataname]["path"]
	  if not(os.path.isfile(filename)):
	      print "dataname ",dataname," filename ",filename," SKIP!! "
	      continue
	  totalevents = getTotalEvents(filename, "Friends")
	  totaljobs = int(totalevents*1.0/events_per_job)+1
	  if "Radion" in job or job == "Data":
	      totaljobs = totaljobs *2
	  allout = os.path.join(outputdir, dataname+"_ijob*.root")
	  nevent_HME = getTotalEvents(allout, "Friends")
	  ratio = float(nevent_HME)/float(totalevents)
	  outfiles = glob.glob(allout)
          iters  =  alliters[0]
	  if len(outfiles)  != totaljobs or abs(ratio -1.0)>0.01:
	      print "error!!! ",dataname, " some output files are missing, expected outfiles ", totaljobs," nevents  ",totalevents," final out ", len(outfiles), " nevent_HME ",nevent_HME
	      for ijob in range(0, totaljobs):
		  outfile =  os.path.join(outputdir, dataname+"_ijob%d.root"%ijob)
		  if not os.path.isfile(outfile):
		      print "sbatch {0}/Send_PlotterProducer_{1}_{2}_{3}.cmd".format(jobdir, dataname, ijob, iters)
		      resubmitscript.write("""
sbatch {0}/Send_PlotterProducer_{1}_{2}_{3}.cmd""".format(jobdir, dataname, ijob, iters))
    resubmitscript.close()
    os.system("chmod +x submitall%s_%s_resubmit.sh"%(jobdir, submitscript_name))


def mergeoutput():
    for job in benchmarks:
      for dataname in full_local_samplelist[job].keys():
	  finalout = os.path.join(outputdir, dataname+"_HME_Friends.root")
	  allout = os.path.join(outputdir, dataname+"_ijob*.root")
	  os.system("hadd -f "+finalout+" "+allout)
	  filename = full_local_samplelist[job][dataname]['path']
	  nevent_HME = getTotalEvents(finalout)
	  nevent_expected = getTotalEvents(filename)
	  ratio = float(nevent_HME)/float(nevent_expected)
	  if abs(ratio-1.0) > 0.01:
	      print "error!!! nevent expected ", nevent_expected," after HME ",nevent_HME," diff ",abs(ratio-1.0)



#python PlotterProducer.py -b HaddNo {jobtype}
#os.system("./submitallHME.sh")
generateJobs()
#checkoutput()
#mergeoutput()
