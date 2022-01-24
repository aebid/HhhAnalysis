#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import warnings

from RunConfiguration import *

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

#sys.path.append('/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD')
#sys.path.append('/afs/cern.ch/work/t/tahuang/HHAnalysis/CMSSW_10_2_0/src/HhhAnalysis/python/NanoAOD')
sys.path.append(os. getcwd())
from countHistogramProducer import *
from genParticleProducer import *
#from HHbbWWProducer import *
from Devin_sync_producer import *


#from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
#btagSF2016_cMVA = lambda : btagSFProducer("2016",  algo = 'cmva', sfFileName='cMVAv2_Moriond17_B_H.csv')## file is under NanoAODTools
#btagSF2017_deepCSV = lambda : btagSFProducer("2017", algo = 'deepcsv', sfFileName='DeepCSV_94XSF_V3_B_F.csv')

#def btagSFyear(year):
#    return{
#	2016: btagSF2016_cMVA(),
#	2017: btagSF2017_deepCSV()
#    }[year]

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
#def jetmetUncertaintiesyear(year):
#    """jet uncertainties and recalibration for MC"""
#    return{
#	2016: jetmetUncertainties2016(),
#	2017: jetmetUncertainties2017()
#    }[year]

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetRecalib import *

from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import *

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.mht import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
#from  PhysicsTools.NanoAODTools.postprocessing.examples.puWeightProducer import *
def puWeightyear(year):
    return {
      2016: puAutoWeight_2016(),
      2017: puAutoWeight_2017(),
      2018: puAutoWeight_2018()
    }[year]

###mhtProducer(jetSelection, muonSelection, electronSelection), adding the selected obj and filling MHT_pt and MHT_phi
mht_hh = lambda : mhtProducer( lambda j : j.pt > 20 and abs(j.eta) < 2.4,
                            lambda mu : mu.pt > 10 and abs(mu.eta) < 2.4,
                            lambda el : el.pt > 10 and abs(el.eta) < 2.5 )

#outputdir = "/fdata/hepx/store/user/taohuang/HH_NanoAOD/"
#filesTTbar = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIFall17_NanoAOD.root"]
#filesSignal = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/GravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV_NanoAOD_2.root"]
#filesTTbar2017 = ["/fdata/hepx/store/mc/RunIIFall17NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/2E9913D6-BAA9-E811-8ABE-0CC47A4DEF3E.root"]
#fileSignal2017 = ["/fdata/hepx/store/mc/RunIIFall17NanoAOD/GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph_correctedcfg/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/70B64607-EEB2-E811-B1D7-A0369FD0B268.root"]
#filesDY1J2017 = ["/fdata/hepx/store/mc/RunIIFall17NanoAODv4/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/90000/62B5DB6B-F42C-D648-A6EF-5F26495985BF.root"]
#file2017Test = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD2017/5E621211-8B42-E811-9903-001E67F8FA2E.root"]
outputdir   = "/eos/user/t/tahuang/2021NtupleProduceTest/"
#file2016Test = ["/eos/user/t/tahuang/NanoAODSamples_2019/myNanoProdMc2016_NANO_brazos_20190807.root"]
#file2017Test = ["/eos/user/t/tahuang/NanoAODSamples_2019/myNanoProdMc2017_NANO_brazos.root"]
#file2018Test = ["/eos/user/t/tahuang/NanoAODSamples_2019/myNanoProdMc2018_NANO_10_20190805Test.root"]

file2016Test = ["/eos/user/t/tahuang/NanoAODSamples_2021/test_m300_nanohadd.root"]
#file2016Test = ["/eos/user/d/daebi/hhBBww_datasets/2016/signal/m270/NanoAODproduction_2016_cfg_NANO_2.root"]

#modules = [ puWeightyear(Runyear), countHistogramAll_2016(), jetmetUncertaintiesyear(Runyear), btagSFyear(Runyear),  mht_hh(), HHbbWWProducer(True, verbose = 1) ]
modules = [ puWeightyear(Runyear),  HHbbWWProducer(True, verbose = 1) ]
#modules = [puWeightyear(Runyear)]
#p=PostProcessor(outputdir, fileSignal2017,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True)
p=PostProcessor(outputdir, file2016Test,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True)
#p=PostProcessor(outputdir, filesTTbar2017,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True)
#p=PostProcessor(outputdir, filesDY1J2017,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True)
print "run Postprocessor here "
p.run()


