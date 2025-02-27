#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

#from HhhAnalysis.NanoAOD.HHbbWWProducer import *
from HHbbWWProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import *
from  PhysicsTools.NanoAODTools.postprocessing.modules.jme.mht import *
#from  PhysicsTools.NanoAODTools.postprocessing.examples.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *

files = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAOD/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/00000/82738046-5711-E811-ADAF-0CC47A4D75F0.root"]
filesTTbar= [
'root://cms-xrd-global.cern.ch//store/group/cmst3/group/nanoAOD/NanoTestProd006/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer17MiniAOD-92X-NanoCrabProd006/171006_155430/0000/nanolzma_1.root',
]
fileMuMu = ["root://cms-xrd-global.cern.ch//store/data/Run2016B/MuonEG/NANOAOD/05Feb2018_ver2-v1/00000/48512742-610C-E811-B1A6-0CC47A1DF800.root"]
filesTTbar = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIFall17_NanoAOD.root"]
filesSignal = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/GravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV_NanoAOD.root"]

selection='''(Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
 Sum$(Muon_pt > 20) >= 2 ||
 Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
 Sum$(Muon_pt > 20 && Muon_tightId) >= 1 ||
 (Sum$(Muon_pt > 20) == 0 && Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) == 0 && MET_pt > 150 ) ) 
 &&  Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) >= 2 
'''

selectionALL='''Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
 Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
 Sum$(Jet_pt > 40 && Jet_jetId) >= 4   || 
Sum$(Jet_pt *(abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) > 160  || 
MET_pt > 100  || Sum$(Muon_pt > 20 && Muon_tightId) >= 1
'''
mht_hh = lambda : mhtProducer( lambda j : j.pt > 20,
                            lambda mu : mu.pt > 5 and mu.pfRelIso04_all < 0.4,
                            lambda el : el.pt > 5 and el.pfRelIso03_all < 0.4 )

#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[jecUncertAll_cppOut(),jetmetUncertaintiesAll(),btagSFProducer("cmva"),vhbb()],provenance=True)
#p.run()

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

p=PostProcessor(".",fileMuMu,"1","keep_and_drop.txt",[mht_hh(),  HHbbWWProducer(True)],provenance=True, jsonInput=runsAndLumis())
#p=PostProcessor(".",inputFiles(),selection.replace('\n',' '),"keep_and_drop.txt",modules=[puWeight(),jetmetUncertaintiesAll(),mhtVHbb(),btagSFProducer("cmva"),hhbbWW()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
p.run()

print "DONE"
os.system("ls -lR")
