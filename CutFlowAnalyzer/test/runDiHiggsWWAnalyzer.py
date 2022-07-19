import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os
import sys

process = cms.Process("DiHiggsAnalyzer")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
options = VarParsing('analysis')
options.register("sampleType", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int)
options.parseArguments()
#print("sys.argv ", sys.argv, " sampleType ", options.sampleType)

process.source = cms.Source("PoolSource",
  secondaryFileNames = cms.untracked.vstring(),
  fileNames = cms.untracked.vstring(
  #'file:/fdata/hepx/store/user/lpernie/TEST_LOCALLY/DYJETS_7A385961-C6D9-E611-85B2-0025905B85BC.root'
  #'/store/mc/RunIISpring16MiniAODv2/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/60000/0005201C-D41B-E611-8E37-002481E0D398.root'
  '/store/mc/RunIISummer16MiniAODv2/GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/F23BE62D-6EC1-E611-8C86-0CC47A57CB62.root'
  )
)


process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(100) 
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
triggerPaths = cms.vstring( 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*',
			     'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*',
			     'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*',
			     'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*',
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

scalefactor_file = "Files/"
#scalefactor_file = ""

process.DiHiggsWWBBAna = cms.EDAnalyzer('DiHiggsWWBBAnalyzer',
  verbose = cms.untracked.int32(0),
  #enum {Data = 0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, TTbar, DYJets, DY0Jets, DY1Jets, DY2Jets, ZZTo2L2Q, ZZTo2L2Nu, ZZTo4L, WWToLNuQQ, WWTo2L2Nu, WZTo2L2Q, WZTo1L3Nu, WZTo1L1Nu2Q, WZTo3LNu, ST_tchannel_top, ST_tchannel_antitop, ST_schannel, ST_tW_antitop, ST_tW_top, WJetsToLNu, WJetsToLNu_HT100To200, WJetsToLNu_HT200To400, WJetsToLNu_HT400To600, WJetsToLNu_HT600To800, WJetsToLNu_HT800To1200, WJetsToLNu_HT1200To2500, WJetsToLNu_HT2500ToInf, TTWJetsToQQ, TTWJetsToLNu, TTZToQQ, TTZToLLNuNu };//add other background
  SampleType = cms.untracked.int32(options.sampleType),
  #genParticles = cms.InputTag("genParticles"),
  genParticles = cms.InputTag("prunedGenParticles"),#minAOD
  genjets = cms.InputTag("slimmedGenJets"),
  #genjets = cms.InputTag("ak4GenJetsNoNu"),
  #genjets = cms.InputTag("ak4GenJets"),

  #trigger matching
  doTriggerMatching = cms.bool(True),
  hltPaths = triggerPaths,
  deltaPtRel_trigger = cms.untracked.double(.5),
  deltaR_trigger  = cms.untracked.double(.1),

  #reco 
  #muons = cms.InputTag("cleanPatPFMuonsTriggerMatch"),
  muons           = cms.InputTag("slimmedMuons"),
  #2016data: Run BCDEF use 2016Medium, GH use Medium
  mu_id           = cms.untracked.string("2016Medium"),
  mu_PFIso        = cms.untracked.double(0.15),#tight iso
  triggerSFFile   = cms.string(scalefactor_file + "EfficienciesAndSF_BCDEF_trigger.root"),
  isoSFFile       = cms.string(scalefactor_file + "EfficienciesAndSF_BCDEF_ISO.root"),
  idSFFile        = cms.string(scalefactor_file + "EfficienciesAndSF_BCDEF_ID.root"),
  trackingSFFile  = cms.string(scalefactor_file + "EfficienciesAndSF_BCDEFGH_Tracking.root"),
  triggerSFhist   = cms.string("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"),
  isoSFhist       = cms.string("TightISO_MediumID_pt_eta/abseta_pt_ratio"),
  idSFhist        = cms.string("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio"),
  trackingSFhist  = cms.string("ratio_eff_eta3_dr030e030_corr"),

  electrons       = cms.InputTag("slimmedElectrons"),
  jets            = cms.InputTag("slimmedJets"),
  mets            = cms.InputTag("slimmedMETs"),
  beamSpot        = cms.InputTag("offlineBeamSpot"),
  triggerEvent    = cms.InputTag("patTriggerEvent"),
  tracks          = cms.InputTag("generalTracks"),
  TriggerResults  = cms.InputTag("TriggerResults","","HLT"),
  TriggerObjects  = cms.InputTag("selectedPatTrigger"),
  TrackRefitter   = cms.InputTag("TrackRefitter"),
  primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  Traj            = cms.InputTag("TrackRefitter"),

  debug = cms.untracked.bool(False),
  onlyGenLevel = cms.bool(False),
  runMMC = cms.bool(False)
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("out_ana.root")
)

process.phlt = cms.Path(process.hltfilter)
process.pDiHiggsWWBBAna = cms.Path(
  process.eventCounterFilter*
  process.hltfilter*
  process.DiHiggsWWBBAna
)

process.pdump = cms.Path(process.dump)

process.schedule = cms.Schedule(process.pDiHiggsWWBBAna)
