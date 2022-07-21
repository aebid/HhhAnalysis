import FWCore.ParameterSet.Config as cms
import os
import sys

process = cms.Process("DiHiggsAnalyzerSL")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
### DL sample
#/store/mc/RunIISummer16MiniAODv3/GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/70000/3A6D44E6-913D-E911-AD45-AC1F6B1AF224.root
### SL sample
#/store/mc/RunIISummer16MiniAODv3/GluGluToRadionToHHTo2B2WToLNu2J_M-270_narrow_13TeV-madgraph/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/270000/60F34423-2E3F-E911-B900-0CC47A6C17FC.root
#


process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
        #'file:/eos/uscms/store/user/tahuang/DiHiggs/out_sim.root'
	#'file:/eos/uscms/store/user/tahuang/DiHiggs/out_miniaod.root'
	#'file:/fdata/hepx/store/user/tahuang/TEST_LOCALLY/Run2016D-23Sep2016_MINIAOD_12B2DEA9-B68C-E611-99A4-0CC47A1DF810.root'
	#'file:/fdata/hepx/store/user/tahuang/TEST_LOCALLY/DYJETS_7A385961-C6D9-E611-85B2-0025905B85BC.root'
	#'file:/fdata/hepx/store/user/tahuang/TEST_LOCALLY/GravitonM400-02DF7FEF-74D9-E611-956A-02163E013746.root'
	#'file:/fdata/hepx/store/user/tahuang/TEST_LOCALLY/RadionM400-3217F073-EF25-E611-862F-A0369F7FC210.root'
	#'/store/data/Run2017C/DoubleMuon/MINIAOD/PromptReco-v3/000/300/742/00000/00CE998C-717E-E711-AD9B-02163E019CB5.root'
	#'/store/data/Run2017C/DoubleMuon/MINIAOD/PromptReco-v2/000/300/636/00000/FC587F5C-E57D-E711-830B-02163E01A3C2.root'
	'/store/mc/RunIISummer16MiniAODv3/GluGluToRadionToHHTo2B2WToLNu2J_M-270_narrow_13TeV-madgraph/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/270000/60F34423-2E3F-E911-B900-0CC47A6C17FC.root'
    )
)

from HhhAnalysis.MCProduction.InputFileHelpers import *
#inputdir = ['/fdata/hepx/store/user/tahuang/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_1M/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_1M/170329_023747/0000/']
#inputdir = ['/fdata/hepx/store/user/tahuang/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_10k/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_10k/170330_023219/0000/']
#process = useInputDir(process, inputdir)
nEvents = 10000
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(nEvents) 
)

process.MessageLogger = cms.Service("MessageLogger", 
    destinations = cms.untracked.vstring("cout"), 
    cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))
)

process.eventCounterFilter = cms.EDFilter("EventCounterFilter")

import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
#print hlt," hlt.triggerResultsFilter ",hlt.triggerResultsFilter
# accept if any path succeeds (explicit OR)
"""
process.hltfilter = hlt.triggerResultsFilter.clone(
	HLTPaths = ('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*','HLT_Mu17*'),
	l1tResults = '',#not use L1t results
	throw = cms.bool(False) 
"""
triggerPaths = cms.vstring( #'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*',
			    # 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*',
			    # 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*',
			    # 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*',
			     'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_*',##added 20170928
	)
process.hltfilter = cms.EDFilter( "TriggerResultsFilter",
        triggerConditions = triggerPaths,
	hltResults = cms.InputTag( "TriggerResults","","HLT"),
	#l1tResults = cms.InputTag( "hltGtDigis" ),
	l1tResults = cms.InputTag( "" ),
	l1tIgnoreMask = cms.bool( False ),
	l1techIgnorePrescales = cms.bool( False ),
	daqPartitions = cms.uint32( 1 ),
	throw = cms.bool(True)    
)


muonPOGSFdir = os.getenv( "CMSSW_BASE" ) +"/src/HhhAnalysis/CutFlowAnalyzer/test/MuonEffAndSF_2016Data/"
process.DiHiggsWWBBAna = cms.EDAnalyzer('DiHiggsWWBBSLAnalyzer',
  verbose = cms.untracked.int32(0),
  #enum {Data = 0,TTbar, DYJets, DY0Jets, DY1Jets, DY2Jets, ZZTo2L2Q, ZZTo2L2Nu, ZZTo4L, WWToLNuQQ, WWTo2L2Nu, WZTo2L2Q, WZTo1L3Nu, WZTo1L1Nu2Q, WZTo3LNu, ST_tchannel_top, ST_tchannel_antitop, ST_schannel, ST_tW_antitop, ST_tW_top, WJetsToLNu, WJetsToLNu_HT100To200, WJetsToLNu_HT200To400, WJetsToLNu_HT400To600, WJetsToLNu_HT600To800, WJetsToLNu_HT800To1200, WJetsToLNu_HT1200To2500, WJetsToLNu_HT2500ToInf, TTWJetsToQQ, TTWJetsToLNu, TTZToQQ, TTZToLLNuNu };//add other background
  #enum {Rad_260=100, Rad_270, Rad_300, Rad_350, Rad_400, Rad_450, Rad_500, Rad_550, Rad_600, Rad_650, Rad_750, Rad_800, Rad_900,  Rad_1000, Grav_260, Grav_270, Grav_300, Grav_350, Grav_400, Grav_450, Grav_500, Grav_550, Grav_600, Grav_650, Grav_750, Grav_800, Grav_900,  Grav_1000};
  SampleType = cms.untracked.int32(100),
  #############gen level
  #genParticles = cms.InputTag("genParticles"),
  genParticles = cms.InputTag("prunedGenParticles"),#minAOD
  genjets = cms.InputTag("slimmedGenJets"),
  #genjets = cms.InputTag("ak4GenJetsNoNu"),
  #genjets = cms.InputTag("ak4GenJets"),

  ##trigger matching
  #doTriggerMatching = cms.bool(True),
  #TriggerResults = cms.InputTag("TriggerResults","","HLT"),
  ##TriggerObjects = cms.InputTag("selectedPatTrigger"),
  #TriggerObjects = cms.InputTag("slimmedPatTrigger"),
 
  #hltPaths = triggerPaths,
  #deltaPtRel_trigger = cms.untracked.double(.5),
  #deltaR_trigger  = cms.untracked.double(.1),

  ##reco 
  ##muons = cms.InputTag("cleanPatPFMuonsTriggerMatch"),
  #muons = cms.InputTag("slimmedMuons"),
  ##2016data: Run BCDEF use 2016Medium, GH use Medium
  #mu_id = cms.untracked.string("2016Medium"),
  #mu_PFIso = cms.untracked.double(0.15),#tight iso
  #triggerSFFile = cms.string(muonPOGSFdir+"EfficienciesAndSF_BCDEF_trigger.root"),
  #isoSFFile = cms.string(muonPOGSFdir+"EfficienciesAndSF_BCDEF_ISO.root"),
  #idSFFile = cms.string(muonPOGSFdir+"EfficienciesAndSF_BCDEF_ID.root"),
  #trackingSFFile = cms.string(muonPOGSFdir+"EfficienciesAndSF_BCDEFGH_Tracking.root"),
  #triggerSFhist = cms.string("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"),
  #isoSFhist = cms.string("TightISO_MediumID_pt_eta/abseta_pt_ratio"),
  #idSFhist = cms.string("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio"),
  #trackingSFhist = cms.string("ratio_eff_eta3_dr030e030_corr"),

  #electrons = cms.InputTag("slimmedElectrons"),
  #jets = cms.InputTag("slimmedJets"),
  mets = cms.InputTag("slimmedMETs"),
  #primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  ##beamSpot = cms.InputTag("offlineBeamSpot"),
  ##tracks = cms.InputTag("generalTracks"),
  ##TrackRefitter = cms.InputTag("TrackRefitter"),
  #Traj = cms.InputTag("TrackRefitter"),

  debug = cms.untracked.bool(False),
  onlyGenLevel = cms.bool(False),

  #simulation = cms.bool(True),
  runMMC = cms.bool(False)
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("out_ana_sl.root")
    #fileName = cms.string("file:/fdata/hepx/store/user/taohuang/DiHiggsAnalysisSample/out_ann_radion_M400_20160411.root")
)

process.phlt = cms.Path(process.hltfilter)
process.pDiHiggsWWBBAna = cms.Path(
  #process.eventCounterFilter*
  #process.hltfilter*
  process.DiHiggsWWBBAna
)

process.pdump = cms.Path(process.dump)

process.schedule = cms.Schedule(process.pDiHiggsWWBBAna)
