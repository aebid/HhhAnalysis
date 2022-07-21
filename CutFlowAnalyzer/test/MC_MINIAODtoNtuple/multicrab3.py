from CRABClient.UserUtilities import config
username = "tahuang"
storesite = "/eos/user/t/tahuang/"
outfolder = "MC_Hhh_analysis_Run2016_20220720"
config = config()

config.section_("General")
config.General.workArea     = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['out_ana.root']
config.JobType.psetName     = "runDiHiggsWWAnalyzer.py"
config.JobType.pyCfgParams  = ["sampleType=100"]
config.JobType.inputFiles   = ["Files/EfficienciesAndSF_BCDEFGH_Tracking.root","Files/EfficienciesAndSF_BCDEF_Tracking.root","Files/EfficienciesAndSF_GH_ISO.root","Files/EfficienciesAndSF_BCDEF_ID.root","Files/EfficienciesAndSF_BCDEF_trigger.root","Files/EfficienciesAndSF_GH_Tracking.root","Files/EfficienciesAndSF_BCDEF_ISO.root","Files/EfficienciesAndSF_GH_ID.root","Files/EfficienciesAndSF_GH_trigger.root"]
#config.JobType.inputFiles   = ["Files/"]

config.section_("Data")
config.Data.inputDBS        = 'global'
config.Data.splitting       = 'FileBased'
config.Data.unitsPerJob     = 1
config.Data.outLFNDirBase   = '/store/user/%s/'%(username) + outfolder
config.Data.publication     = False

config.section_("Site")
config.Site.storageSite     =  'T3_CH_CERNBOX' ## or T3_US_FNALLPC, 'T3_US_TAMU'
# Submit or write the Plotter files
OnlySubmitCRAB = True

import os
import glob
import operator
def findNewestDir(directory):
  os.chdir(directory)
  dirs = {}
  for dir in glob.glob('*'):
    if os.path.isdir(dir):
      dirs[dir] = os.path.getctime(dir)
  lister = sorted(dirs.iteritems(), key=operator.itemgetter(1))
  return lister[-1][0]

datasets  = []; NumSample = []; sampleN_short = []
doTT=False; doDY=False; doVV=False; doSingleT=False; doWjets=False; dottV=False; doRadion=True; doGravition=False;
#doTT=True; doDY=True; doVV=True; doSingleT=True; doWjets=True; dottV=True; doRadion=True; doGravition=True;
# SIGNAL
if doRadion:
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('100'); sampleN_short.append('Rad_260')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM') 
  NumSample.append('101'); sampleN_short.append('Rad_270')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('102'); sampleN_short.append('Rad_300')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('103'); sampleN_short.append('Rad_350')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('104'); sampleN_short.append('Rad_400')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('105'); sampleN_short.append('Rad_450')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('106'); sampleN_short.append('Rad_500')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('107'); sampleN_short.append('Rad_550')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('108'); sampleN_short.append('Rad_600')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('109'); sampleN_short.append('Rad_650')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('110'); sampleN_short.append('Rad_750')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('111'); sampleN_short.append('Rad_800')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('112'); sampleN_short.append('Rad_900')
  datasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('113'); sampleN_short.append('Rad_1000')

if doGravition:
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM') 
  NumSample.append('114'); sampleN_short.append('Grav_260')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('115'); sampleN_short.append('Grav_270')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('116'); sampleN_short.append('Grav_300')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('117'); sampleN_short.append('Grav_350')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('118'); sampleN_short.append('Grav_400')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('119'); sampleN_short.append('Grav_450')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('120'); sampleN_short.append('Grav_500')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('121'); sampleN_short.append('Grav_550')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('122'); sampleN_short.append('Grav_600')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('123'); sampleN_short.append('Grav_650')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-700_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('124'); sampleN_short.append('Grav_700')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('125'); sampleN_short.append('Grav_800')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('126'); sampleN_short.append('Grav_900')
  datasets.append('/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('127'); sampleN_short.append('Grav_1000')
# TT
if doTT:
  datasets.append('/TTTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM')
  NumSample.append('13'); sampleN_short.append('TT')
# DY
if doDY:
  datasets.append('/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('14'); sampleN_short.append('DY')
  datasets.append('/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('15'); sampleN_short.append('DY')
  datasets.append('/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('16'); sampleN_short.append('DY')
  datasets.append('/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('17'); sampleN_short.append('DY')
# VV
if doVV:
  datasets.append('/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('18'); sampleN_short.append('VV')
  datasets.append('/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('19'); sampleN_short.append('VV')
  datasets.append('/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('20'); sampleN_short.append('VV')
  datasets.append('/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('21'); sampleN_short.append('VV')
  datasets.append('/WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('22'); sampleN_short.append('VV')
  datasets.append('/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('23'); sampleN_short.append('VV')
  datasets.append('/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('24'); sampleN_short.append('VV')
  datasets.append('/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM')
  NumSample.append('25'); sampleN_short.append('VV')
  datasets.append('/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('26'); sampleN_short.append('VV')
# Single-t
if doSingleT:
  datasets.append('/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('27'); sampleN_short.append('sT')
  datasets.append('/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('28'); sampleN_short.append('sT')
  datasets.append('/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('29'); sampleN_short.append('sT')
  datasets.append('/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('30'); sampleN_short.append('sT')
  datasets.append('/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('31'); sampleN_short.append('sT')
# W + Jets
if doWjets:
  datasets.append('/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('32'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
  NumSample.append('33'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
  NumSample.append('34'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('35'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('36'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('37'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('38'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('39'); sampleN_short.append('Wjet')
# tt + V
if dottV:
  datasets.append('/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('40'); sampleN_short.append('ttV')
  datasets.append('/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
  NumSample.append('41'); sampleN_short.append('ttV')
  datasets.append('/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('42'); sampleN_short.append('ttV')
  datasets.append('/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/MINIAODSIM')
  NumSample.append('43'); sampleN_short.append('ttV')

plotter_f = open("for_plotter.py",'w')
check_f   = open("check_crab.sh",'w')
check_f.write("#!/bin/bash\n")
if __name__ == '__main__':
  from CRABAPI.RawCommand import crabCommand
  from multiprocessing import Process
  from CRABClient.ClientExceptions import ClientException
  from httplib import HTTPException
  def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))
  i=0
  args = []; args.append("-99")
  configlist = []
  lastSampleShort="NOSAMPLESHORT"
  for dataset in datasets:
    # To check the status
    check_f.write("crab status -d crab_projects/crab_"+dataset.split('/')[1]+"\n")
    # To plot easily the datasets
    if not OnlySubmitCRAB:
      sampleN   = dataset.split('/')[1].split('_')[0]+dataset.split('/')[1].split('_')[1]
      path      = storesite + outfolder + dataset.split('/')[1] + "/crab_" + dataset.split('/')[1] + "/"
      NewestDir = findNewestDir(path)
      path      = path + NewestDir
      newSample = True
      if(sampleN_short[i]==lastSampleShort): newSample = False
      if((newSample and lastSampleShort!="NOSAMPLESHORT") or i==int(len(NumSample)-1) ):
        plotter_f.write('  this_cat      = "cat HADD/' + sampleN_short[i-1] + '_*' ' > HADD/' + sampleN_short[i-1] + '.txt"\n')
        plotter_f.write('  this_hadd     = "hadd -T -f -k ' + storesite + outfolder + oldDataset + '/crab_' + oldDataset + '.root @HADD/' + sampleN_short[i-1] + '.txt"\n')
        plotter_f.write('  this_NtotPath = " ' + storesite + outfolder + oldDataset + '/crab_' + oldDataset + '.root"\n')
      if(newSample):
        plotter_f.write('if( whichSample == "' + sampleN_short[i] + '" ):\n')
        plotter_f.write('  Find_str.append("find ' + path + ' | grep root | grep -v failed > HADD/' + sampleN_short[i] + "_" + sampleN + '.txt")\n')
      lastSampleShort = sampleN_short[i]
    if OnlySubmitCRAB:
      #os.system("cp runDiHiggsWWAnalyzer.py runDiHiggsWWAnalyzer_%s.py"%sampleN_short[i])
      #config.JobType.psetName  = "runDiHiggsWWAnalyzer_%s.py"%sampleN_short[i]
      config.Data.inputDataset = dataset
      config.General.requestName = dataset.split('/')[1]
      args[0] = "sampleType=%s"%NumSample[i]
      config.JobType.pyCfgParams = args
      #crabCommand('submit', config = config)
      #p = Process(target=crabCommand('submit', config=config),  args=())
      p = Process(target=submit, args=(config,))
      p.start()
      p.join()
    i=i+1
    oldDataset = dataset.split('/')[1]
print "bash check_crab.sh > RESULTS.txt"
