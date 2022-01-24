#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import warnings

#from RunConfiguration import *
import __builtin__
#__builtin__.Runyear = 2016
#Runyear = 2016
__builtin__.Runyear = 2017
#Runyear = 2017

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

<<<<<<< Updated upstream
#sys.path.append('/afs/cern.ch/work/t/tahuang/HHAnalysis/CMSSW_10_2_0/src/HhhAnalysis/python/NanoAOD')
sys.path.append('/afs/cern.ch/work/d/daebi/diHiggs/CMSSW_10_2_26/src/HhhAnalysis/python/NanoAOD')
=======
sys.path.append('/afs/cern.ch/work/t/tahuang/HHAnalysis/CMSSW_10_2_0/src/HhhAnalysis/python/NanoAOD')
#sys.path.append('/afs/cern.ch/work/d/daebi/diHiggs/CMSSW_10_2_0/src/HhhAnalysis/python/NanoAOD')
>>>>>>> Stashed changes
#from countHistogramProducer import *
from genParticleProducer import *
#from HHbbWWProducer_sync import *
#from Devin_sync_producer import *
from HHbbWWProducer_v2 import *

from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
btagSF2016_cMVA = lambda : btagSFProducer("2016",  algo = 'cmva', sfFileName='cMVAv2_Moriond17_B_H.csv')## file is under NanoAODTools
btagSF2017_deepCSV = lambda : btagSFProducer("2017", algo = 'deepcsv', sfFileName='DeepCSV_94XSF_V3_B_F.csv')
btagSF2018_deepCSV = lambda : btagSFProducer("2018", algo = 'deepcsv', sfFileName='DeepCSV_102XSF_V1.csv')
def btagSFyear(year):
    return{
	2016: btagSF2016_cMVA(),
	2017: btagSF2017_deepCSV(),
	2018: btagSF2018_deepCSV()
    }[year]

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
def jetmetUncertaintiesyear(year):
    """jet uncertainties and recalibration for MC"""
    return{
	2016: jetmetUncertainties2016(),
	2017: jetmetUncertainties2017(),
	2018: jetmetUncertainties2018()
    }[year]

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


#outputdir   = "tao_oldsample_mytest/"
#outputdir   = "add_cov_test/"
#outputdir   = "mysample_test_m{mass}/".format(mass = mass)
#outputdir = "sync/"
#file2016Test = ["add_cov/NANOAODSIM_M-260_1.root"]
#file2016Test = ["tao_oldsample/myNanoProdMc2016_NANO_brazos_20190807.root"]
outputdir   = "/eos/user/t/tahuang/2021NtupleProduceTest/"
#file2016Test = ["/eos/user/t/tahuang/NanoAODSamples_2019/myNanoProdMc2016_NANO_brazos_20190807.root"]
#file2017Test = ["/eos/user/t/tahuang/NanoAODSamples_2019/myNanoProdMc2017_NANO_brazos.root"]
#file2018Test = ["/eos/user/t/tahuang/NanoAODSamples_2019/myNanoProdMc2018_NANO_10_20190805Test.root"]

file2016Test = ["/eos/user/t/tahuang/Sync_2021/sync_2016_m750.root"]
file2017Test = ["/eos/user/t/tahuang/Sync_2021/sync_2017_m750.root"]
file2018Test = ["tao_oldsample/myNanoProdMc2018_NANO_10_20190805Test.root"]

#file2016Test = ["test/test_m300_nanohadd.root"]
#file2016Test = ["sync/sync_2016_m750.root"]
#file2017Test = ["tao_oldsample/myNanoProdMc2017_NANO_brazos.root"]
#file2018Test = ["tao_oldsample/myNanoProdMc2018_NANO_10_20190805Test.root"]

#modules = [ puWeightyear(Runyear), countHistogramAll_2016(), jetmetUncertaintiesyear(Runyear), btagSFyear(Runyear),  mht_hh(), HHbbWWProducer(True, verbose = 1) ]
modules = [ puWeightyear(Runyear),   HHbbWWProducer(True, verbose = 1) ]
#p=PostProcessor(outputdir, filesTTbar2017,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True)
#p=PostProcessor(outputdir, filesDY1J2017,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True)
##note: make minor change to postprocessor.py to change the default output TTree name:
## https://github.com/tahuang1991/nanoAOD-tools/blob/master/python/postprocessing/framework/postprocessor.py
p=PostProcessor(outputdir, file2017Test,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True, outtreeName="syncTree")
print "run Postprocessor here "
p.run()

