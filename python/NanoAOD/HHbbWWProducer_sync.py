import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
print "ROOT version ",ROOT.gROOT.GetVersion()
from math import sqrt,cos
import random
#random.randint(0,100)

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..
import sys
sys.path.append('/afs/cern.ch/work/t/tahuang/HHAnalysis/CMSSW_10_2_0/src/HhhAnalysis/python/NanoAOD')

import __builtin__

if hasattr(__builtin__, "Runyear"):
    Runyear = __builtin__.Runyear
else:
    import RunConfiguration as RunConfig
    Runyear = RunConfig.Runyear

#from RunConfiguration import *
if Runyear == 2016:
    import POGRecipesRun2016 as POGRecipesRun2
elif Runyear == 2017:
    import POGRecipesRun2017 as POGRecipesRun2
elif Runyear == 2018:
    import POGRecipesRun2018 as POGRecipesRun2
else:
    sys.exit("wrong run year value: %d"%Runyear)


import Run2DiLeptonTrigger as Run2Trigger

print "HHbbWWProducer, import finished here, Runyear ", Runyear

errorcolor1 = '\x1b[1;31m'
errorcolor2 = '\x1b[0m'

class HHbbWWProducer(Module):
    ## data or MC, which L1 trigger, HLT?
    ###kwargs: triggertype, verbose, run_lumi
    def __init__(self,isMC, **kwargs):
	print "init HHbbWWProduer"
	self.writeHistFile=True
        self.isMC = isMC ## two mode: data or MC
        self.looseEvtSelection = True
        self.noBtaggingCut = True
	print "kwargs ",kwargs
	self.triggertype  = ''##"DoubleMuon, DoubleEG, MuonEG"
	self.DYestimation = False
	self.CheckBtaggingEff = False
	self.deltaR_trigger_reco = 0.1; self.deltaPtRel_trigger_reco = 0.5
	self.verbose = 3
	self.muonEta = 2.4; self.EGEta = 2.5
	self.leadingMuonPt = {"DoubleMuon": 20.0, "MuonEG":25.0}
	self.subleadingMuonPt = {"DoubleMuon": 10.0, "MuonEG":10.0}
	self.leadingEGPt = {"DoubleEG": 25.0, "MuonEG":25.0}
	self.subleadingEGPt = {"DoubleEG": 15.0, "MuonEG":15.0}
	self.jetPt = 20; self.jetEta = 2.4
	self.deltaR_j_l = 0.3 #jet,lepton seperation
	self.maxnjets = 20

	self.run_lumi = None
	self.__dict__.update((key, kwargs[key]) for key in kwargs.keys())
	self.CheckBtaggingEff = (self.CheckBtaggingEff and self.DYestimation and self.isMC)
        print "self.run_lumi  ",self.run_lumi," trigger type ",self.triggertype," verbose ",self.verbose, " DYestimation ",self.DYestimation

	self.addGenToTree = False
	self.addSystematic = False

	self.ll_sys_branchname = ["Isosf","IDsf","trgsf","trackingsf","HLTsafeIDsf"]
	self.jjbtag_sys_branchname = ["jjbtag_light","jjbtag_heavy"]
	self.lep1_branchname_sysvalue = {}
	self.lep2_branchname_sysvalue = {}
	self.jjbtag_branchname_sysvalue = {}


	### for debug
	self.nEv_DoubleEG = 0
	self.nEv_DoubleMuon = 0
	self.nEv_MuonEG = 0
	#self.lepSFmanager = POGRecipesRun2.LeptonSFManager(True)

    def beginJob(self, histFile=None,histDirName=None):
	print "BeginJob "
	#Module.beginJob(self,histFile,histDirName)
	#self.addObject(self.h_cutflow)
	#self.addObject(self.h_cutflow_weight)
	#self.addObject(self.h_cutflowlist["DoubleMuon"])
	#self.addObject(self.h_cutflowlist["DoubleEG"])
	#self.addObject(self.h_cutflowlist["MuonEG"])

    def endJob(self):
        pass


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	print "BeginFiles "
	self.h_eventcounter = ROOT.TH1F("h_raweventcounter","h_raweventcounter", 20, 0, 20.0)
	self.h_cutflow = ROOT.TH1F("h_cutflow","h_cutflow", 20, 0, 20.0)
	self.h_cutflow_weight = ROOT.TH1F("h_cutflow_weight","h_cutflow_weight", 200, 0, 200.0)
	self.h_cutflowlist = {"DoubleMuon":ROOT.TH1F("h_cutflow_DoubleMuon","h_cutflow_DoubleMuon", 20, 0, 20.0), 
			      "DoubleEG":ROOT.TH1F("h_cutflow_DoubleEG","h_cutflow_DoubleEG", 20, 0, 20.0), 
			      "MuonEG":ROOT.TH1F("h_cutflow_MuonEG","h_cutflow_MuonEG", 20, 0, 20.0)}

	self.h_cutflowlist["DoubleMuon"].SetLineColor(ROOT.kRed)
	self.h_cutflowlist["DoubleEG"].SetLineColor(ROOT.kBlue)
	self.h_cutflowlist["MuonEG"].SetLineColor(ROOT.kBlack)
	
        self.out = wrappedOutputTree
        self.out.branch("ak4Jet1_idx",  "I");
        self.out.branch("ak4Jet1_pt",  "F");
        self.out.branch("ak4Jet1_E",  "F");
        self.out.branch("ak4Jet1_eta",  "F");
        self.out.branch("ak4Jet1_phi",  "F");
        self.out.branch("ak4Jet1_cMVAv2",  "F");
        self.out.branch("ak4Jet1_CSV",  "F");
        self.out.branch("ak4Jet2_idx",  "I");
        self.out.branch("ak4Jet2_pt",  "F");
        self.out.branch("ak4Jet2_E",  "F");
        self.out.branch("ak4Jet2_eta",  "F");
        self.out.branch("ak4Jet2_phi",  "F");
        self.out.branch("ak4Jet2_cMVAv2",  "F");
        self.out.branch("ak4Jet2_CSV",  "F");
        self.out.branch("ak4Jet3_idx",  "I");
        self.out.branch("ak4Jet3_pt",  "F");
        self.out.branch("ak4Jet3_E",  "F");
        self.out.branch("ak4Jet3_eta",  "F");
        self.out.branch("ak4Jet3_phi",  "F");
        self.out.branch("ak4Jet3_cMVAv2",  "F");
        self.out.branch("ak4Jet3_CSV",  "F");
        self.out.branch("ak4Jet4_idx",  "I");
        self.out.branch("ak4Jet4_pt",  "F");
        self.out.branch("ak4Jet4_E",  "F");
        self.out.branch("ak4Jet4_eta",  "F");
        self.out.branch("ak4Jet4_phi",  "F");
        self.out.branch("ak4Jet4_cMVAv2",  "F");
        self.out.branch("ak4Jet4_CSV",  "F");

        self.out.branch("ak8Jet1_pt",            "F");
        self.out.branch("ak8Jet1_E",             "F");
        self.out.branch("ak8Jet1_eta",           "F");
        self.out.branch("ak8Jet1_phi",           "F");
        self.out.branch("ak8Jet1_msoftdrop",     "F");
        self.out.branch("ak8Jet1_tau1",          "F");
        self.out.branch("ak8Jet1_tau2",          "F");
        self.out.branch("ak8Jet1_btagHbb",       "F");
	self.out.branch("ak8Jet1_subjet0_pt",    "F");
	self.out.branch("ak8Jet1_subjet0_eta",   "F");
	self.out.branch("ak8Jet1_subjet0_phi",   "F");
	self.out.branch("ak8Jet1_subjet0_CSV",   "F");   
	self.out.branch("ak8Jet1_subjet1_pt",    "F");
	self.out.branch("ak8Jet1_subjet1_eta",   "F");
	self.out.branch("ak8Jet1_subjet1_phi",   "F");
	self.out.branch("ak8Jet1_subjet1_CSV",   "F");   
        self.out.branch("ak8Jet2_pt",            "F");
        self.out.branch("ak8Jet2_E",             "F");
        self.out.branch("ak8Jet2_eta",           "F");
        self.out.branch("ak8Jet2_phi",           "F");
        self.out.branch("ak8Jet2_msoftdrop",     "F");
        self.out.branch("ak8Jet2_btagHbb",       "F");
        self.out.branch("ak8Jet2_tau1",          "F");
        self.out.branch("ak8Jet2_tau2",          "F");
	self.out.branch("ak8Jet2_subjet0_pt",    "F");
	self.out.branch("ak8Jet2_subjet0_eta",   "F");
	self.out.branch("ak8Jet2_subjet0_phi",   "F");
	self.out.branch("ak8Jet2_subjet0_CSV",   "F");   
	self.out.branch("ak8Jet2_subjet1_pt",    "F");
	self.out.branch("ak8Jet2_subjet1_eta",   "F");
	self.out.branch("ak8Jet2_subjet1_phi",   "F");
	self.out.branch("ak8Jet2_subjet1_CSV",   "F");   

        self.out.branch("ak8lsJet1_pt",            "F");
        self.out.branch("ak8lsJet1_E",             "F");
        self.out.branch("ak8lsJet1_eta",           "F");
        self.out.branch("ak8lsJet1_phi",           "F");
        self.out.branch("ak8lsJet1_msoftdrop",     "F");
        self.out.branch("ak8lsJet1_tau1",          "F");
        self.out.branch("ak8lsJet1_tau2",          "F");
        self.out.branch("ak8lsJet1_btagHbb",       "F");
	self.out.branch("ak8lsJet1_subjet0_pt",    "F");
	self.out.branch("ak8lsJet1_subjet0_eta",   "F");
	self.out.branch("ak8lsJet1_subjet0_phi",   "F");
	self.out.branch("ak8lsJet1_subjet0_CSV",   "F");   
	self.out.branch("ak8lsJet1_subjet1_pt",    "F");
	self.out.branch("ak8lsJet1_subjet1_eta",   "F");
	self.out.branch("ak8lsJet1_subjet1_phi",   "F");
	self.out.branch("ak8lsJet1_subjet1_CSV",   "F");   
        self.out.branch("ak8lsJet2_pt",            "F");
        self.out.branch("ak8lsJet2_E",             "F");
        self.out.branch("ak8lsJet2_eta",           "F");
        self.out.branch("ak8lsJet2_phi",           "F");
        self.out.branch("ak8lsJet2_msoftdrop",     "F");
        self.out.branch("ak8lsJet2_btagHbb",       "F");
        self.out.branch("ak8lsJet2_tau1",          "F");
        self.out.branch("ak8lsJet2_tau2",          "F");
	self.out.branch("ak8lsJet2_subjet0_pt",    "F");
	self.out.branch("ak8lsJet2_subjet0_eta",   "F");
	self.out.branch("ak8lsJet2_subjet0_phi",   "F");
	self.out.branch("ak8lsJet2_subjet0_CSV",   "F");   
	self.out.branch("ak8lsJet2_subjet1_pt",    "F");
	self.out.branch("ak8lsJet2_subjet1_eta",   "F");
	self.out.branch("ak8lsJet2_subjet1_phi",   "F");
	self.out.branch("ak8lsJet2_subjet1_CSV",   "F");   


	#self.out.branch("lepstype",  "F")
        self.out.branch("mu1_pt",  "F");
        self.out.branch("mu1_TrgObjfilterbits",  "I");
        self.out.branch("mu1_E",  "F");
        self.out.branch("mu1_eta",  "F");
        self.out.branch("mu1_phi",  "F");
        self.out.branch("mu1_pdgid",  "I");
        self.out.branch("mu1_charge",  "F");
        self.out.branch("mu1_sip3D",  "F");
        self.out.branch("mu1_miniRelIso",  "F");
        self.out.branch("mu1_dxy",  "F");
        self.out.branch("mu1_dxyAbs",  "F");
        self.out.branch("mu1_dz",  "F");
        #self.out.branch("mu1_miniRelIso_chg",  "F");
        self.out.branch("mu1_segmentCompatibility",  "F");
        self.out.branch("mu1_leptonMVA",  "F");
	self.out.branch("mu1_miniRelIsoCharged",        "F");
	#self.out.branch("mu1_miniRelIsoNeutral",        "F");
        self.out.branch("mu2_pt",  "F");
        self.out.branch("mu2_E",  "F");
        self.out.branch("mu2_eta",  "F");
        self.out.branch("mu2_phi",  "F");
        self.out.branch("mu2_pdgid",  "I");
        self.out.branch("mu2_TrgObjfilterbits",  "I");
        self.out.branch("mu2_charge",  "F");
        self.out.branch("mu2_sip3D",  "F");
        self.out.branch("mu2_miniRelIso",  "F");
        self.out.branch("mu2_dxy",  "F");
        self.out.branch("mu2_dxyAbs",  "F");
        self.out.branch("mu2_dz",  "F");
        #self.out.branch("mu2_miniRelIso_chg",  "F");
        self.out.branch("mu2_segmentCompatibility",  "F");
        self.out.branch("mu2_leptonMVA",  "F");
	self.out.branch("mu2_miniRelIsoCharged",        "F");
	#self.out.branch("mu2_miniRelIsoNeutral",        "F");

        self.out.branch("ele1_pt",  "F");
        self.out.branch("ele1_TrgObjfilterbits",  "I");
        self.out.branch("ele1_E",  "F");
        self.out.branch("ele1_eta",  "F");
        self.out.branch("ele1_phi",  "F");
        self.out.branch("ele1_pdgid",  "I");
        self.out.branch("ele1_charge",  "F");
        self.out.branch("ele1_sip3D",  "F");
        self.out.branch("ele1_miniRelIso",  "F");
        self.out.branch("ele1_dxy",  "F");
        self.out.branch("ele1_dxyAbs",  "F");
        self.out.branch("ele1_dz",  "F");
        #self.out.branch("ele1_miniRelIso_chg",  "F");
        self.out.branch("ele1_ntMVAeleID",  "F");
        self.out.branch("ele1_leptonMVA",  "F");
	self.out.branch("ele1_miniRelIsoCharged",        "F");
	#self.out.branch("ele1_miniRelIsoNeutral",        "F");
	self.out.branch("ele1_mvaIdFall17noIso",         "F");
        self.out.branch("ele2_pt",  "F");
        self.out.branch("ele2_E",  "F");
        self.out.branch("ele2_eta",  "F");
        self.out.branch("ele2_phi",  "F");
        self.out.branch("ele2_pdgid",  "I");
        self.out.branch("ele2_TrgObjfilterbits",  "I");
        self.out.branch("ele2_charge",  "F");
        self.out.branch("ele2_sip3D",  "F");
        self.out.branch("ele2_miniRelIso",  "F");
        self.out.branch("ele2_dxy",  "F");
        self.out.branch("ele2_dxyAbs",  "F");
        self.out.branch("ele2_dz",  "F");
        #self.out.branch("ele2_miniRelIso_chg",  "F");
        self.out.branch("ele2_leptonMVA",  "F");
        self.out.branch("ele2_ntMVAeleID",  "F");
	self.out.branch("ele2_miniRelIsoCharged",        "F");
	self.out.branch("ele2_miniRelIsoNeutral",        "F");
	self.out.branch("ele2_mvaIdFall17noIso",         "F");

	self.out.branch("PFMET",  "F")
	self.out.branch("PFMETphi",  "F")
	self.out.branch("n_presel_mu", "I")
	self.out.branch("n_presel_ele", "I")
	self.out.branch("n_presel_ak4Jet", "I")
	self.out.branch("n_presel_ak8Jet", "I")
	self.out.branch("n_presel_ak8lsJet", "I")
	"""
	self.out.branch("nJetsL",  "F")
	self.out.branch("jj_M",  "F")
	self.out.branch("el_hltsafeid", "F")
	self.out.branch("llreco_weight",  "F")
	self.out.branch("llid_weight",  "F")
	self.out.branch("lliso_weight",  "F")
	self.out.branch("lltrigger_weight",  "F")
	self.out.branch("ht_jets",  "F")
	self.out.branch("ht",  "F")
	self.out.branch("ll_M",  "F")
	self.out.branch("llmet_M",  "F")
	self.out.branch("llmetjj_MT2",  "F")
	self.out.branch("llmetjj_M",  "F")
	self.out.branch("lljj_M",  "F")
	self.out.branch("cosThetaStar",  "F")
	self.out.branch("ll_DR_l_l",  "F")
	self.out.branch("jj_DR_j_j",  "F")
	self.out.branch("llmetjj_DPhi_ll_jj",  "F")
	self.out.branch("llmetjj_DPhi_ll_met",  "F")
	self.out.branch("llmetjj_DPhi_llmet_jj",  "F")
	self.out.branch("ll_pt",  "F")
	self.out.branch("ll_eta",  "F")
	self.out.branch("jj_pt",  "F")
	self.out.branch("jj_eta",  "F")
	self.out.branch("llmetjj_minDR_l_j",  "F")
	self.out.branch("llmetjj_MTformula",  "F")
	self.out.branch("ll_DPhi_l_l",  "F")
	self.out.branch("ll_DEta_l_l",  "F")
	self.out.branch("event_pu_weight",  "F")
	self.out.branch("event_lep_weight",  "F")
	self.out.branch("event_btag_weight",  "F")
	self.out.branch("jjbtag_heavy",  "F")
	self.out.branch("jjbtag_light",  "F")
	"""
        self.out.branch("PU_weight",  "F");
        self.out.branch("MC_weight",  "F");
	self.out.branch("pu",  "F")
	self.out.branch("run",  "I")
	self.out.branch("ls",  "I")##lumiblocks
	self.out.branch("event",  "I")
	#self.out.branch("event_weight",  "F")
	#self.out.branch("DY_BDT_flat",  "F")
	#self.out.branch("dy_nobtag_to_btagM_weight",  "F")
	#self.out.branch("mt2",  "F")
	#self.out.branch("mt2_bb",  "F")
	#self.out.branch("mt2_ll",  "F")
	#self.out.branch("event_number",  "I")

	"""
	if self.isMC and self.addSystematic:
	    for name in self.ll_sys_branchname:
		self.out.branch("lep1%s"%name,"F")
		self.out.branch("lep1%s_up"%name,"F")
		self.out.branch("lep1%s_down"%name,"F")
		self.out.branch("lep2%s"%name,"F")
		self.out.branch("lep2%s_up"%name,"F")
		self.out.branch("lep2%s_down"%name,"F")
	    self.out.branch("ll_Tallintrgeff","F")
	    self.out.branch("ll_Tallintrgeff_up","F")
	    self.out.branch("ll_Tallintrgeff_down","F")
	    for name in self.jjbtag_sys_branchname:
		self.out.branch("%s"%name,"F")
		self.out.branch("%s_up"%name,"F")
		self.out.branch("%s_down"%name,"F")

	    self.out.branch("event_pu_weight_up",  "F")
	    self.out.branch("event_pu_weight_down",  "F")
	    self.out.branch("event_pdf_weight_up",  "F")
	    self.out.branch("event_pdf_weight_down",  "F")
	    #self.out.branch("LHEScaleWeight", "F", lenVar = "nLHEScaleWeight")
	    for i in range(0, 9):
		self.out.branch("LHEScaleWeight_%d"%i, "F")
	if self.DYestimation or self.addGenToTree:	
	    self.out.branch("genjet1_partonFlavour",  "I");
	    self.out.branch("genjet2_partonFlavour",  "I");
	    self.out.branch("genmet_pt",  "F");
	    self.out.branch("genmet_phi",  "F");
	##how to add gen information???
	if self.isMC and self.addGenToTree:
	    self.out.branch("genb1_pt",  "F");##quark level
	    self.out.branch("genjet1_pt",  "F");
	    self.out.branch("genjet1_E",  "F");
	    self.out.branch("genjet1_eta",  "F");
	    self.out.branch("genjet1_phi",  "F");
	    self.out.branch("genb2_pt",  "F");
	    self.out.branch("genjet2_pt",  "F");
	    self.out.branch("genjet2_E",  "F");
	    self.out.branch("genjet2_eta",  "F");
	    self.out.branch("genjet2_phi",  "F");
	    #self.out.branch("genW1_pt",  "F");
	    #self.out.branch("genW1_mass",  "F");
	    #self.out.branch("genW2_pt",  "F");
	    #self.out.branch("genW2_mass",  "F");
	    #self.out.branch("gennu1_pt",  "F");##nu1 from W1
	    #self.out.branch("gennu1_eta",  "F");
	    #self.out.branch("gennu1_phi",  "F");
	    #self.out.branch("gennu2_pt",  "F");
	    #self.out.branch("gennu2_eta",  "F");
	    #self.out.branch("gennu2_phi",  "F");
	    self.out.branch("genl1_pt",  "F");##genl1 from lep1, leading 
	    self.out.branch("genl1_E",  "F");
	    self.out.branch("genl1_eta",  "F");
	    self.out.branch("genl1_phi",  "F");
	    self.out.branch("genl1_id",  "I");
	    self.out.branch("genl2_pt",  "F");
	    self.out.branch("genl2_E",  "F");
	    self.out.branch("genl2_eta",  "F");
	    self.out.branch("genl2_phi",  "F");
	    self.out.branch("genl2_id",  "I");
	    self.out.branch("dR_genb1_genjet1",  "F");
	    self.out.branch("dR_genb1_genjet2",  "F");
	    self.out.branch("dR_genb2_genjet1",  "F");
	    self.out.branch("dR_genb2_genjet2",  "F");
	    self.out.branch("dR_genl1_lepFromW1",  "F");
	    self.out.branch("dR_genl1_lepFromW2",  "F");
	    self.out.branch("dR_genl2_lepFromW1",  "F");
	    self.out.branch("dR_genl2_lepFromW2",  "F");
	    self.out.branch("met_diNuetrino_pt",  "F");
	    self.out.branch("met_diNuetrino_phi",  "F");
	    if self.CheckBtaggingEff:
		self.out.branch("alljets_pt", "F", n=self.maxnjets)
		self.out.branch("alljets_eta", "F", n=self.maxnjets)
		self.out.branch("alljets_cMVAv2", "F", n=self.maxnjets)
		self.out.branch("alljets_DeepCSV", "F", n=self.maxnjets)
		#self.out.branch("alljets_partonFlavour", "I", n=self.maxnjets)
		#self.out.branch("alljets_hadronFlavour", "I", n=self.maxnjets)
		self.out.branch("alljets_genpartonFlavour", "I", n=self.maxnjets)
	    """
	
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	self.h_eventcounter.Write()
	self.h_cutflow.Write()
	self.h_cutflowlist["DoubleMuon"].Write()
	self.h_cutflowlist["DoubleEG"].Write()
	self.h_cutflowlist["MuonEG"].Write()
	self.h_cutflow_weight.Write()
	#inputTree.SetBranchStatus("Jet_btagSF_*", 0)
	#c1 = ROOT.TCanvas("h_cutflow","h_cutflow")
	#self.h_cutflowlist["MuonEG"].Draw("hist")
	#self.h_cutflowlist["DoubleEG"].Draw("histsame")
	#self.h_cutflowlist["DoubleMuon"].Draw("histsame")
	#c1.Write()
	#c2 = ROOT.TCanvas("h_cutflow_weight","h_cutflow_weight")
	#c2.cd()
	#self.h_cutflow_weight.Draw("hist")
	#c2.Write()

    #def goldenJason(self, run, lumi):
    #    run_str = '%d'%run
    #    if run_str in self.run_lumi.keys():
    #        alllumis = self.run_lumi[run_str]
    #        for lumirange in alllumis:
    #            if lumi >= lumirange[0] and lumi <= lumirange[1]:
    #                return True
    #    else:
    #        return False


    def findingLeptonPairs(self, leptons_mu, leptons_el):

	leptonpairs = []
	alltriggertypes = []
	if self.isMC:
	    alltriggertypes = ["DoubleMuon", "DoubleEG", "MuonEG"]
	else:
	    alltriggertypes = [self.triggertype]


	for triggertype in alltriggertypes:
	    if triggertype == "DoubleMuon" and len(leptons_mu) >= 2:
		nmuons = len(leptons_mu)
		for imu1 in range(0, nmuons):
		    for imu2 in range(imu1+1, nmuons):
			if leptons_mu[imu1].charge * leptons_mu[imu2].charge < 0 and leptons_mu[imu1].pt > self.leadingMuonPt["DoubleMuon"] and leptons_mu[imu2].pt > self.subleadingMuonPt["DoubleMuon"]:
			    leptonpairs.append(( leptons_mu[imu1],  leptons_mu[imu2] ))

	    elif triggertype == "DoubleEG" and len(leptons_el) >= 2:
		nelectrons = len(leptons_el)
		for iel1 in range(0, nelectrons):
		    for iel2 in range(iel1+1, nelectrons):
			if leptons_el[iel1].charge * leptons_el[iel2].charge < 0 and leptons_el[iel1].pt > self.leadingEGPt["DoubleEG"] and leptons_el[iel2].pt > self.subleadingEGPt["DoubleEG"]:
			    leptonpairs.append(( leptons_el[iel1],  leptons_el[iel2] ))

	    elif triggertype == "MuonEG" and len(leptons_mu) >= 1 and len(leptons_el) >= 1:
		nmuons = len(leptons_mu)
		nelectrons = len(leptons_el)
		for imu1 in range(0, nmuons):
		    for iel1 in range(0, nelectrons):
			if leptons_mu[imu1].charge * leptons_el[iel1].charge > 0:
			    continue
			if (leptons_mu[imu1].pt >= self.leadingMuonPt["MuonEG"] and leptons_el[iel1].pt >= self.subleadingEGPt["MuonEG"]) or (leptons_mu[imu1].pt >= self.subleadingMuonPt["MuonEG"] and leptons_el[iel1].pt >= self.leadingEGPt["MuonEG"]):
			    if leptons_mu[imu1].pt > leptons_el[iel1].pt:
				leptonpairs.append(( leptons_mu[imu1],  leptons_el[iel1] ))
			    else:
				leptonpairs.append(( leptons_el[iel1],  leptons_mu[imu1] ))

	if self.verbose > 2:
	    for leps in leptonpairs:
		print "Lepton pair leading ",leps[0].pdgId, " pt ",leps[0].pt," subleading ",leps[1].pdgId," pt ",leps[1].pt
        return leptonpairs


    def pt(self, jet, isMC):
        ## the MC has JER smearing applied which has output branch Jet_pt_smeared which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if isMC:
            return jet.pt_smeared
        else:
            return jet.pt		    
    
    def met(self, met, isMC):
        ## the MC has JER smearing applied which has output branch met_[pt/phi]_smeared which should be compared 
        ## with data branch MET_[pt/phi]. This essentially aliases the two branches to one common variable.
        if isMC:
            return (met.pt_smeared,met.phi_smeared)
        else:
            return (met.pt,met.phi)

    def printLepton(self, lep, string):
	print string," pdgid ",lep.pdgId," pt ",lep.pt, " eta ",lep.eta, " phi ",lep.phi
    def printTrigObjs(self, objs, string):
	print string
	for obj in objs:
	    print "trigger obj id ",obj.id, " pt ",obj.pt," l1pt ",obj.l1pt, " eta ",obj.eta, " phi ", obj.phi," filterbits ",obj.filterBits
	   
    def matchTriggerObjwithHLT(self, path, trigobjs):
	##triggerobj filterBits: 
	##old one in 2018:
	##1 = CaloIdL_TrackIdL_IsoVL, 2 = WPLoose, 4 = WPTight, 8 = OverlapFilter PFTau for Electron (PixelMatched e/gamma); 
	##1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau for Muon
	##new one in 2019:
	##electron: 1 = CaloIdL_TrackIdL_IsoVL, 2 = 1e (WPTight), 4 = 1e (WPLoose), 8 = OverlapFilter PFTau, 16 = 2e, 32 = 1e-1mu, 64 = 1e-1tau, 128 = 3e, 256 = 2e-1mu, 512 = 1e-2mu
	##Muon: 1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = 1mu, 16 = 2mu, 32 = 1mu-1e, 64 = 1mu-1tau, 128 = 3mu, 256 = 2mu-1e, 512 =1mu-2e 

	alllegs = path.split('_')[1:-1]
	leg_trigobj_map  = {};
        leg_trigobj_map['leg1'] = []
        leg_trigobj_map['leg2'] = []
	foundleg1 = False; foundleg2 = False
	testobjs = []
	legpts = []
	for x in alllegs:
	    legobjmap = None
	    if "Mu" in x:
	    	pt = 8.0 ## HLT obj pt
	    	if "TkMu" not in x:
		    pt = int(x[2:])
                else: 
		    pt = int(x[4:])
                legpts.append(pt)
		#matchedobjs = list(filter(lambda obj: obj.id == 13 and obj.pt >= pt and obj.l1pt > pt/2.0 and (obj.filterBits & 3) > 0, trigobjs))
		matchedobjs = list(filter(lambda obj: obj.id == 13 and obj.pt >= pt and (obj.filterBits & (1 | 2 | 16 | 32)) > 0, trigobjs))
	        if self.verbose > 3: 
		    testobjs = list(filter(lambda obj: obj.id == 13, trigobjs))
		    self.printTrigObjs(testobjs, "testobjs Mu legpt%d "%pt)


	        if not(foundleg1):
		    leg_trigobj_map['leg1'] = matchedobjs
		    foundleg1 = True
		elif not (foundleg2):
		    leg_trigobj_map['leg2'] = matchedobjs
		    foundleg2 = True
		else:
		    print errorcolor1+"error in matching HLT to trigger objects ", path, errorcolor2
	    elif "Ele" in x:
		pt = int(x[3:])
                legpts.append(pt)
		#matchedobjs = list(filter(lambda obj: obj.id == 11 and obj.pt >= pt and (obj.filterBits & 1) > 0, trigobjs))
		matchedobjs = list(filter(lambda obj: obj.id == 11 and obj.pt >= pt and (obj.filterBits & (1 | 16 | 32)) > 0 , trigobjs))
	        if not(foundleg1):
		    leg_trigobj_map['leg1'] = matchedobjs
		    foundleg1 = True
		elif not (foundleg2):
		    leg_trigobj_map['leg2'] = matchedobjs
		    foundleg2 = True
		else:
		    print errorcolor1+"error in matching HLT to trigger objects ", path,errorcolor2

	if self.verbose > 3: 
	    print "firedpath ", path," legpts ",legpts, " matched l1objs, leg1 ",len(leg_trigobj_map['leg1']), " leg2 ",len(leg_trigobj_map['leg2'])


	return leg_trigobj_map 
         
	#selectedObjs = []
	#for x in alllegs:
	#    if "Mu" in x:
	#    	pt = 8.0
	#    	if "TkMu" not in x:
	#	    pt = int(x[2:])
	#	for obj in trigobjs:
	#	    #if obj.id == 13 and (obj.filterBits & 1)>0 and obj.pt >= pt and (obj.l1pt >= 4.0 or obj.l1pt_2 >= 4.0):
	#	    #if obj.id == 13 and obj.pt >= pt:
	#	    #if obj.id == 13 and obj.pt >= pt and obj.l1pt >= 4.0 and (obj.filterBits & 3) > 0:
	#	    if obj.id == 13 and obj.pt >= pt and (obj.filterBits & 3) > 0:
	#	    	selectedObjs.append(obj)
	#	    elif  obj.id == 13 and  self.verbose >3 :
	#	        print "path ",path, " trig pt ",pt, " trig obj pt ",obj.pt, " l1pt ",obj.l1pt," filterbits ",obj.filterBits
	#    elif "Ele" in x:
	#	pt = int(x[3:])
	#	for obj in trigobjs:
	#	    #if obj.id == 11 and (obj.filterBits & 1)>0 and obj.pt >= pt and (obj.l1pt >= 10.0 or obj.l1pt_2 >= 10.0):##mini pt for L1 seed? 10GeV?
	#	    #if obj.id == 11 and obj.pt >= pt:##mini pt for L1 seed? 10GeV?
	#	    if obj.id == 11 and obj.pt >= pt and  (obj.filterBits & 1) > 0:##mini pt for L1 seed? 10GeV?
	#	    	selectedObjs.append(obj)
	#	    elif  obj.id == 11 and  self.verbose >3 :
	#	        print "path ",path, " trig pt ",pt, " trig obj pt ",obj.pt, " l1pt ",obj.l1pt," filterbits ",obj.filterBits

	#return set(selectedObjs)

    def findTriggerType(self, leptons):
	if abs(leptons[0].pdgId) == 13 and abs(leptons[1].pdgId) == 13:
	    return "DoubleMuon"
	elif abs(leptons[0].pdgId) == 11 and  abs(leptons[1].pdgId) == 11:
	    return "DoubleEG"
	elif (abs(leptons[0].pdgId) == 13 and  abs(leptons[1].pdgId) == 11) or (abs(leptons[0].pdgId) == 11 and  abs(leptons[1].pdgId) == 13):
	    return "MuonEG"
	else:
	    print "in findTriggerType, lepton id is 11 or 13, error!!!"
	    return ""

    def findfiredPaths(self, event, allHLTPaths):
	return filter(lambda path : hasattr(event, path) and getattr(event,  path, False), allHLTPaths)

    def matchHLTPath(self, event, run, trigobjs, leptons, matchedTrigobjs = []):
	""" check which HLT is fired """
	##whether cut on l1pt for triggerobjects?
	## DoubleMu, L1 seed: DoubleMu_11_4 or 12_5
	## MuonEG: L1 seed: Mu20_EG10, L1_Mu5_IsoEG18..
	##DoubleEG: L1_SingleIsoEG22, L1_DoubleEG_15_10...

	trigobjs.sort(key=lambda x: x.pt,reverse=True)##sort trigobjs by pt increasing order
	l1t = self.findTriggerType(leptons)	
        if not  self.isMC and l1t != self.triggertype:
	     print  errorcolor1+"error!!! for data, the lepton's type is not the same as trigger type "+errorcolor2
	     return False

        HLTpaths, L1Pts = Run2Trigger.findHLTPathsAndL1Pts(Runyear, self.triggertype, run, self.isMC)

	def l1ptcut(path, leg1Objs, leg2Objs):
	    """check whether trigger objects matched to two legs could pass l1pt cuts or not"""
	    if len(leg1Objs)  == 0 and len(leg2Objs) == 0:
	        return False,list(),list()
	    #l1ptcut = l1ptcutForPath(path)
	    l1ptcut_list = L1Pts[path] #

	    for l1ptcut in l1ptcut_list:
		leg1Objs_l1ptcut = list(filter(lambda obj : obj.l1pt >= l1ptcut[0], leg1Objs))
		leg2Objs_l1ptcut = list(filter(lambda obj : obj.l1pt >= l1ptcut[1], leg2Objs))
		passl1pt_double = (len(leg1Objs_l1ptcut)  >= 1 and len(leg2Objs_l1ptcut) >= 1 and not(leg1Objs_l1ptcut == leg2Objs_l1ptcut and len(leg1Objs_l1ptcut) == 1))
                passl1pt_single = (l1ptcut[1] < 0.1 and len(leg1Objs_l1ptcut)  >= 1) or (l1ptcut[0] < 0.1 and len(leg2Objs_l1ptcut)  >= 1) ## single seed
                if passl1pt_double or passl1pt_single:
		    if passl1pt_single and self.verbose > 3:
		        print "pass L1 pt cuts with single seed !! l1ptcut_list ",l1ptcut_list
		    return True,leg1Objs_l1ptcut, leg2Objs_l1ptcut
	    return False,list(),list()


        
        allfiredHLTs = self.findfiredPaths(event, HLTpaths)
        if self.verbose > 3:
	    print "all fired paths ",allfiredHLTs
	    self.printLepton(leptons[0], "leading lep"); self.printLepton(leptons[1], "subleading lep")

        matchTolep1 = lambda tobj : tobj.id == abs(leptons[0].pdgId) and deltaR(tobj.eta, tobj.phi, leptons[0].eta, leptons[0].phi) < self.deltaR_trigger_reco and abs(tobj.pt-leptons[0].pt)/leptons[0].pt < self.deltaPtRel_trigger_reco 
        matchTolep2 = lambda tobj : tobj.id == abs(leptons[1].pdgId) and deltaR(tobj.eta, tobj.phi, leptons[1].eta, leptons[1].phi) < self.deltaR_trigger_reco and abs(tobj.pt-leptons[1].pt)/leptons[1].pt < self.deltaPtRel_trigger_reco 
        for  path in allfiredHLTs:
	    leg_trigobj_map = self.matchTriggerObjwithHLT(path, trigobjs)
            goodLeg1objs_lep1 = list(filter(matchTolep1, leg_trigobj_map['leg1']))
            goodLeg2objs_lep1 = list(filter(matchTolep1, leg_trigobj_map['leg2']))
            goodLeg1objs_lep2 = list(filter(matchTolep2, leg_trigobj_map['leg1']))
            goodLeg2objs_lep2 = list(filter(matchTolep2, leg_trigobj_map['leg2']))
	    if self.verbose > 3: 
		self.printTrigObjs( leg_trigobj_map['leg1'], "triggerobj matched to leg1")
		self.printTrigObjs( leg_trigobj_map['leg2'], "triggerobj matched to leg2")
		print "goodLeg1objs_lep1 ",len(goodLeg1objs_lep1), " goodLeg2objs_lep1 ",len(goodLeg2objs_lep1)," goodLeg1objs_lep2 ",len(goodLeg1objs_lep2)," goodLeg2objs_lep2 ",len(goodLeg2objs_lep2)
	         
	    
	    ##add l1pt cut:
	    fired_option1, leg1list_opt1, leg2list_opt1 = l1ptcut(path, goodLeg1objs_lep1, goodLeg2objs_lep2)
	    fired_option2, leg1list_opt2, leg2list_opt2 = l1ptcut(path, goodLeg1objs_lep2, goodLeg2objs_lep1)
            if fired_option1  or fired_option2:
	        if fired_option1:  matchedTrigobjs = [leg1list_opt1, leg2list_opt1]
		if fired_option2:  matchedTrigobjs = [leg2list_opt2, leg1list_opt2]
	        if self.verbose > 3:
		    print "passing L1pt cuts, find fired HLT path finally! matched trig objs ",matchedTrigobjs
    		return True
	        
	return False


    def trigger_reco_matching(self, lepton, trigobjs):
	""" match reco object to triggger obj """
	matchedObj = None; matched_dR = self.deltaR_trigger_reco; matched_dPtRel = self.deltaPtRel_trigger_reco
	for tobj in trigobjs:
	    if abs(lepton.pdgId) != tobj.id:
		continue
	    dR = deltaR(tobj.eta, tobj.phi, lepton.eta, lepton.phi)
	    dPtRel = abs(tobj.pt-lepton.pt)/lepton.pt 
	    ##check id, dR, dpt
	    if dR < self.deltaR_trigger_reco and dPtRel < self.deltaPtRel_trigger_reco:
		if self.verbose > 3 and matchedObj != None:
		    print "trigger obj l1 pt ", tobj.l1pt," l2pt ", tobj.l2pt," pt ",tobj.pt, " eta ",tobj.eta, " bits ", tobj.filterBits, " offlep pt ",lepton.pt," eta ",lepton.eta," id ",lepton.pdgId," got a second match trig obj: with trig obj pt ",matchedObj.pt," eta ",matchedObj.eta," phi ",matchedObj.phi
		if matchedObj == None or (matchedObj != None and dPtRel < matched_dPtRel): 
		    matchedObj = tobj; matched_dR = dR; matched_dPtRel = dPtRel
	return matchedObj

		    


    def fillCutFlow(self, cutflow_bin, weight):
	#cutflow_bin += 1.0
	if self.verbose > 4:
	    print "final filling cutflow ",cutflow_bin," triggertype ",self.triggertype," weight ",weight
	self.h_cutflow_weight.Fill(weight)##used to mornitor the weight
	while cutflow_bin > 0:
	    self.h_cutflow.Fill( cutflow_bin, weight)
	    self.h_cutflowlist[self.triggertype].Fill( cutflow_bin, weight)
	    cutflow_bin = cutflow_bin - 1

    def fillCutFlow_leptons(self, cutflow_bin, weight, leptonpairs):
	allL1triggertype = []
	if self.isMC:
	    allL1triggertype = [self.findTriggerType(x) for x in leptonpairs]
	else:
	    allL1triggertype = [self.triggertype]
	
	#only fill for first case
	#print "allL1triggertype ",allL1triggertype
	#l1t = allL1triggertype[0]
	while cutflow_bin > 0:
	    self.h_cutflow.Fill( cutflow_bin, weight)
            for l1t in allL1triggertype:
		self.h_cutflowlist[l1t].Fill( cutflow_bin, weight)
	    cutflow_bin = cutflow_bin - 1
	

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
	#hlts = Collection(event, "HLT")

        self.h_eventcounter.Fill(1)	
        ### PV	
	PV = Object(event, "PV")
	event_pu_weight = 1.0; event_pu_weight_up = 1.0; event_pu_weight_down = 1.0
	LHEPdfweight_stddev = 0.0; LHEScaleWeight = None
	pu = PV.npvsGood
        run = getattr(event,"run", False)	
	luminosityBlock = getattr(event, "luminosityBlock", False)
	ievent =  getattr(event,"event", False)
	RhoFastjetCentralCalo = getattr(event, "fixedGridRhoFastjetCentralCalo", 0.0)
	event_reco_weight = 1.0 ## for pu_weight*btag_SF*lepSF
	sample_weight = 1.0
	genweight = 1.0
	if self.isMC:
	    event_pu_weight = event.puWeight
	    event_pu_weight_up = event.puWeightUp/event_pu_weight
	    event_pu_weight_down = event.puWeightDown/event_pu_weight
	    nvtx = int(getattr(event, "Pileup_nTrueInt"))
	    if self.verbose > 3:
		print "pugood ", pu, " nvtx ",nvtx, " event.puWeight ",event.puWeight," up ",event.puWeightUp," down ",event.puWeightDown
	    genweight = getattr(event, "genWeight", 1)
	    sample_weight = genweight * event.puWeight
	    
	    if self.addSystematic and hasattr(event, "LHEPdfWeight"):
		LHEPdfWeight   = event.LHEPdfWeight
	        LHEPdfWeight_sum = sum(LHEPdfWeight)
		nLHEPdfWeight = len(LHEPdfWeight)
		LHEPdfWeight_ave = LHEPdfWeight_sum/nLHEPdfWeight
		LHEPdfweight_tmplist = [(x - LHEPdfWeight_ave)*(x-LHEPdfWeight_ave) for x in LHEPdfWeight]
		LHEPdfweight_stddev = sqrt(sum(LHEPdfweight_tmplist)/(len(LHEPdfWeight)-1.0))
    		#print "sum ",sum(LHEPdfWeight)," ave ",LHEPdfWeight_ave," LHEPdfweight_stddev ",LHEPdfweight_stddev
	    if self.addSystematic and hasattr(event, "LHEScaleWeight"):
		LHEScaleWeight = list(event.LHEScaleWeight)
    		#print "LHEScaleWeight ",LHEScaleWeight, " type ", type(LHEScaleWeight)

	if self.verbose > 1:
	    print "run ",run," luminosityBlock ",luminosityBlock," ievent ",ievent," sample_weight ",sample_weight
	cutflow_bin = 0
	##to get even weight sum
	##first time to fill cutflow histogram at the begining 
	self.h_cutflowlist["DoubleMuon"].Fill( cutflow_bin, sample_weight )
	self.h_cutflowlist["DoubleEG"].Fill( cutflow_bin, sample_weight )
	self.h_cutflowlist["MuonEG"].Fill( cutflow_bin, sample_weight )
	self.h_cutflow.Fill( cutflow_bin, sample_weight)

	### MET
        met = Object(event, "MET")
	metPt = met.pt; metPhi = met.phi
	if hasattr(event, "MET_pt_nom"):
	    metPt = getattr(event, "MET_pt_nom")
	    metPhi = getattr(event, "MET_phi_nom")
	    #print "use calibrated MET pt ",metPt," phi ",metPhi," old pt ",met.pt," phi ",met.phi

        electrons = list(Collection(event, "Electron"))
        muons = list(Collection(event, "Muon"))
	muon_pt_corrected = None 
	if hasattr(event, "Muon_pt_corrected"):
	    muon_pt_corrected = getattr(event, "Muon_pt_corrected")
	    for muon in muons:
	        imu = muons.index(muon)
    		#print "imu ",imu, " pt ",muon.pt," corrected pt ",muon_pt_corrected[imu]
		muon.pt = muon_pt_corrected[imu]
	
        jets = list(Collection(event, "Jet"))
	mht = Object(event, "MHT")
	#jet_btagSF = None; jet_btagSF_up = None; jet_btagSF_down = None
	#if self.isMC:
	#    jet_btagSF      = event.Jet_btagSF
	#    jet_btagSF_up   = event.Jet_btagSF_up
	#    jet_btagSF_down = event.Jet_btagSF_down
	#

        #####################################
        ## di-leps selection
        #####################################



	###cutstep: dilepton pt, and eta
	for mu in muons:
	    if  self.verbose > 3:
	        print "Muon id ",mu.pdgId, " pt ",mu.pt," eta ",mu.eta
	lepton_muons = list(filter(lambda x : POGRecipesRun2.muonPreselection(x), muons))
   
	for el in electrons:
	    if  self.verbose > 3:
	        print "Electron id ", el.pdgId, " pt ", el.pt," eta ", el.eta
	## check eta and whether it passes the conversion veto(not from photon conversion), convVeto = True: good electron
	lepton_electrons_pre = list(filter(lambda x :  POGRecipesRun2.electronPreselection(x), electrons))

	def ele_muon_dR(ele):
	    for mu in lepton_muons:
	        if deltaR(ele.eta, ele.phi, mu.eta, mu.phi) < 0.3:
	            return False
	    return True

        lepton_muons.sort(key=lambda x:x.pt,reverse=True)	
	lepton_electrons = [x for x in lepton_electrons_pre if ele_muon_dR(x)]
        lepton_electrons.sort(key=lambda x:x.pt,reverse=True)	
        #for el in lepton_electrons:
	#    print "Electron id ", el.pdgId, " pt ", el.pt," eta ", el.eta
        
        
        #####################################
        ## di-jets selection
        #####################################
        # no recalibration-2019/06/13
	#if hasattr(event, "Jet_pt_nom"):
	#    jet_pt_corrected = getattr(event, "Jet_pt_nom")
	#    for jet in jets:
        #        ijet = jets.index(jet)
        #        #print "Jet pt from reCalibration : ",jet_pt_corrected[ijet]," old pt ",jet.pt
        #        jets[ijet].pt = jet_pt_corrected[ijet]
	#
	###cut: JetPt, JetEta
        
        jets_pre = [x for x in jets if  POGRecipesRun2.ak4jetPreselection(x)]
	#print "#jets ",len(jets)," #jets_pre ",len(jets_pre)
	def jet_lep_dR(jet):
	    for mu in lepton_muons:
	        if deltaR(jet.eta, jet.phi, mu.eta, mu.phi) < 0.4:
	            return False
	    for ele in lepton_electrons:
	        if deltaR(jet.eta, jet.phi, ele.eta, ele.phi) < 0.4:
	            return False
	    return True

	jets_clean = [x for x in jets_pre if jet_lep_dR(x)]
	        

	#### select final two jets: jets with maximum pT , for DY estimation 
        ## Run2016: sort alljets with a cMVAv2 descreasing order 
	## Run2017: use DeepCSV
	#if Runyear == 2016:
	#    jets_clean.sort(key=lambda x:x.btagCMVA, reverse=True)	
        #elif Runyear == 2017 or Runyear== 2018:
	#    jets_clean.sort(key=lambda x:x.btagDeepB, reverse=True)	
        #else:
	#    sys.exit(errorcolor1+"wrong run year value: %d"%Runyear+errorcolor2)

        ak8jets = list(Collection(event, "FatJet"))
	ak8subjets = list(Collection(event, "SubJet"))
	def fatjet_subjetcut(jet):
	    subjet1_idx = jet.subJetIdx1
	    subjet2_idx = jet.subJetIdx2
	    #print "subjet1_idx ",subjet1_idx," subjet2_idx ",subjet2_idx," njet ",len(ak8subjets)
	    if subjet1_idx < 0 or subjet2_idx < 0 or subjet1_idx >= len(ak8subjets)  or subjet2_idx >= len(ak8subjets):
		return False
	    subjet1 = ak8subjets[subjet1_idx]
	    subjet2 = ak8subjets[subjet2_idx]
	    subjet_pt_btagging_pass = False
	    #subjet_pt = subjet1.pt if subjet1.btagDeepB > subjet2.btagDeepB else subjet2.pt
	    #subjet_pt_btagging_pass = (subjet1.btagDeepB > 0.4941 or subjet2.btagDeepB > 0.4941) and subjet_pt > 30
	    #if subjet1.btagDeepFlavB < subjet2.btagDeepFlavB:
	    if POGRecipesRun2.Ak8subjetMediumBtagging(subjet1) and subjet1.pt > 30:
	        subjet_pt_btagging_pass = True
	    if POGRecipesRun2.Ak8subjetMediumBtagging(subjet2) and subjet2.pt > 30:
	        subjet_pt_btagging_pass = True
	    return subjet1.pt > 20 and abs(subjet1.eta)<2.4 and subjet2.pt>20 and abs(subjet2.eta)<2.4 and subjet_pt_btagging_pass

	def fajet_subjetCSVsum(jet):  
	    subjet1_idx = jet.subJetIdx1
	    subjet2_idx = jet.subJetIdx2
	    subjet1 = ak8subjets[subjet1_idx]
	    subjet2 = ak8subjets[subjet2_idx]
	    return subjet1.btagDeepB+subjet2.btagDeepB
	

	def fatjet_lep_dR(jet):
	    for mu in lepton_muons:
	        if deltaR(jet.eta, jet.phi, mu.eta, mu.phi) < 0.8:
	            return False
	    for ele in lepton_electrons:
	        if deltaR(jet.eta, jet.phi, ele.eta, ele.phi) < 0.8:
	            return False
	    return True


	ak8jets_pre = [x for x in ak8jets if POGRecipesRun2.ak8jetPreselection(x)]
	ak8jets_clean = [x for x in ak8jets_pre if fatjet_lep_dR(x) and fatjet_subjetcut(x)] 
	
        ak8jets_CSVsum_clean = ak8jets_clean
	ak8jets_CSVsum_clean.sort(key=lambda x: fajet_subjetCSVsum(x), reverse=True)
	
        mismatchEvents = []## for debugging
        ak8lsjets = list(Collection(event, "FatJetAK8LSLoose"))
	ak8lssubjets = list(Collection(event, "SubJetAK8LSLoose"))
	def fatjetak8ls_subjetcut(jet):
	    subjet1_idx = jet.subJetIdx1
	    subjet2_idx = jet.subJetIdx2
	    if ievent in mismatchEvents:
		print luminosityBlock," ",ievent,"\t AK8jet pt ",jet.pt," eta ",jet.eta ," phi ",jet.phi," subjet1_idx ",subjet1_idx," subjet2_idx ",subjet2_idx," njet ",len(ak8lssubjets)
	    if subjet1_idx < 0 or subjet2_idx < 0 or subjet1_idx >= len(ak8lssubjets)  or subjet2_idx >= len(ak8lssubjets):
		return False
	    subjet1 = ak8lssubjets[subjet1_idx]
	    subjet2 = ak8lssubjets[subjet2_idx]
	    if ievent in mismatchEvents:
	        print ievent," subjet1 pt ",subjet1.pt," eta ",subjet1.eta," subjet2 pt ",subjet2.pt," eta ",subjet2.eta
	  
	    return subjet1.pt > 20 and abs(subjet1.eta)<2.4 and subjet2.pt>20 and abs(subjet2.eta)<2.4

	def fatjetak8ls_lep_dR(jet):
	    ##dR_lep_fajet > 0.1
	    #for mu in lepton_muons:
	    #    if deltaR(jet.eta, jet.phi, mu.eta, mu.phi) < 0.1:
	    #        return False
	    #for ele in lepton_electrons:
	    #    if deltaR(jet.eta, jet.phi, ele.eta, ele.phi) < 0.1:
	    #        return False

	    ##min(dR_lep_fajet) < 1.2
	    #for mu in lepton_muons:
	    #    if deltaR(jet.eta, jet.phi, mu.eta, mu.phi) < 1.2:
	    #        return True
	    #for ele in lepton_electrons:
	    #    if deltaR(jet.eta, jet.phi, ele.eta, ele.phi) < 1.2:
	    #        return True
	    #return False
	    dR_lep_jet = []
	    #dR_lep_subjet = []
	    mindR_lep_jet = 999.0
	    mindR_lep_subjet = 999.0
	    subjet1_idx = jet.subJetIdx1
	    subjet2_idx = jet.subJetIdx2
	    #print "subjet1_idx ",subjet1_idx," subjet2_idx ",subjet2_idx," njet ",len(jets)
	    if subjet1_idx < 0 or subjet2_idx < 0 or subjet1_idx >= len(ak8lssubjets)  or subjet2_idx >= len(ak8lssubjets):
		return False
	    subjet1 = ak8lssubjets[subjet1_idx]
	    subjet2 = ak8lssubjets[subjet2_idx]
	    for mu in lepton_muons:
	        thisdR_lep_jet = deltaR(jet.eta, jet.phi, mu.eta, mu.phi)
		if ievent in mismatchEvents:
		    print "\t\t muon pt ",mu.pt, " eta ",mu.eta," phi ",mu.phi, " dR ",thisdR_lep_jet
	        dR_lep_jet.append(thisdR_lep_jet)
	        if thisdR_lep_jet < mindR_lep_jet:
	            mindR_lep_jet = thisdR_lep_jet
	            if thisdR_lep_jet < 1.2: 
	                mindR_lep_subjet = min(deltaR(subjet1.eta, subjet1.phi, mu.eta, mu.phi), deltaR(subjet2.eta, subjet2.phi, mu.eta, mu.phi))
	    for ele in lepton_electrons:
	        thisdR_lep_jet = deltaR(jet.eta, jet.phi, ele.eta, ele.phi)
		if ievent in mismatchEvents:
		    print "\t\t electron pt ", ele.pt, " eta ", ele.eta," phi ", ele.phi," dR ",thisdR_lep_jet
	        dR_lep_jet.append(thisdR_lep_jet)
	        if thisdR_lep_jet < mindR_lep_jet :
	            mindR_lep_jet = thisdR_lep_jet
	            if thisdR_lep_jet < 1.2: 
	                mindR_lep_subjet = min (deltaR(subjet1.eta, subjet1.phi, ele.eta, ele.phi), deltaR(subjet2.eta, subjet2.phi, ele.eta, ele.phi))

	    if ievent in mismatchEvents:
		print ievent," dR_lep_jet ",dR_lep_jet," mindR_lep_jet ",mindR_lep_jet," mindR_lep_subjet ",mindR_lep_subjet
	    return mindR_lep_jet < 1.2 and  mindR_lep_subjet > 0.1
            

	def fatjetak8ls_bjetclean(jet):
	    #Hbbfatjet = False
	    #for fatjet in ak8jets_clean:
	    #    ##require fatjet passing the loose Hbb btagging cut 
	    #    if fatjet.btagHbb>0.5803:
	    #        Hbbfatjet = True
	    #    if deltaR(jet.eta, jet.phi, fatjet.eta, fatjet.phi) < 1.6 and fatjet.btagHbb>0.5803:
	    #        return False
            #if Hbbfatjet:
	    #    return True 
	    if ievent in mismatchEvents:
	        print "#ak8jets ",len(ak8jets_CSVsum_clean)," #ak4jets ",len(jets_clean)," Ak8lsjet pt ",jet.pt," eta ",jet.eta," phi ",jet.phi
	    if len(ak8jets_CSVsum_clean) >= 1:
	        fatjet = ak8jets_CSVsum_clean[0]
	        if deltaR(jet.eta, jet.phi, fatjet.eta, fatjet.phi) < 1.6:
		    return False
		else:
		    return True
	    jets_btagging = jets_clean
	    ## which b-tagging method we should use
	    jets_btagging.sort(key=lambda x:x.btagDeepFlavB,reverse=True)
            
	    if ievent in mismatchEvents:
		for j in jets_btagging:
		    dR =  deltaR(jet.eta, jet.phi, j.eta, j.phi)
		    print ievent," ak4jet index ",jets_btagging.index(j), " pt ",j.pt," eta ",j.eta," phi ",j.phi," CVS ",j.btagDeepB," dR(ak4j,ak8lsj) ",dR
	    if len(jets_btagging) >= 2:
	        if deltaR(jet.eta, jet.phi, jets_btagging[0].eta, jets_btagging[0].phi) < 1.2 or deltaR(jet.eta, jet.phi, jets_btagging[1].eta, jets_btagging[1].phi) < 1.2:
	            return False
		else: 
		    return True
	        
	    return True



	if ievent in mismatchEvents:
	    print "event ",ievent," ls ",luminosityBlock," all AK8jet "
	    for jet in ak8lsjets:
		print "\t AK8jet pt ",jet.pt," eta ",jet.eta ," phi ",jet.phi," jetid ",jet.jetId, " tau1 ",jet.tau1," tau2 ",jet.tau2," tau2/tau1 ",jet.tau2/jet.tau1 
	ak8lsjets_pre = [x for x in ak8lsjets if  POGRecipesRun2.ak8lsjetPreselection(x)]
	ak8lsjets_prev2 = [x for x in ak8lsjets_pre if fatjetak8ls_subjetcut(x) and fatjetak8ls_lep_dR(x)] 
	ak8lsjets_clean = [x for x in ak8lsjets_prev2 if fatjetak8ls_bjetclean(x)]
	if ievent in mismatchEvents:
	    print ievent," #ak8lsjets ",len(ak8lsjets)," # ak8lsjets_pre ",len(ak8lsjets_pre)," # ak8lsjets_prev2 ",len(ak8lsjets_prev2)," # ak8lsjets_clean ",len(ak8lsjets_clean)



        jets_clean.sort(key=lambda x:x.pt,reverse=True)	
        ak8jets_clean.sort(key=lambda x:x.pt,reverse=True)	
        ak8lsjets_clean.sort(key=lambda x:x.pt,reverse=True)	




	"""
        #jet selection, old
        #jetsForHiggs = [x for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and x.Pt>20 and abs(x.eta)<2.5]
        #if (len(jetsForHiggs) < 2): return False
        self.out.fillBranch("jet1_idx",  hJidx[0]);
        self.out.fillBranch("jet1_pt",  hj1_p4.Pt());
        self.out.fillBranch("jet1_E",  hj1_p4.E());
        self.out.fillBranch("jet1_eta",  hj1_p4.Eta());
        self.out.fillBranch("jet1_phi",  hj1_p4.Phi());
        self.out.fillBranch("jet1_cMVAv2",  jet1.btagCMVA);
        self.out.fillBranch("jet1_CSV",  jet1.btagDeepB);
        self.out.fillBranch("jet2_idx",  hJidx[1]);
        self.out.fillBranch("jet2_pt", hj2_p4.Pt());
        self.out.fillBranch("jet2_E",  hj2_p4.E());
        self.out.fillBranch("jet2_eta", hj2_p4.Eta());
        self.out.fillBranch("jet2_phi",  hj2_p4.Phi());
        self.out.fillBranch("jet2_cMVAv2",  jet2.btagCMVA);
        self.out.fillBranch("jet2_CSV",  jet2.btagDeepB);
	#self.out.fillBranch("isElEl", isElEl) # 0 or 1
	#self.out.fillBranch("isElMu",  isElMu) # 0 or 1, mu_pt<el_pt
	#self.out.fillBranch("isMuEl",  isMuEl) # 0 or 1
	#self.out.fillBranch("isMuMu",  isMuMu) # 0 or 1
	#self.out.fillBranch("isSF",  isSF)
	#self.out.fillBranch("lepstype", Lepstype)
        self.out.fillBranch("lep1_pt",  lep1_p4.Pt());
        self.out.fillBranch("lep1_E",  lep1_p4.E());
        self.out.fillBranch("lep1_eta",  lep1_p4.Eta());
        self.out.fillBranch("lep1_phi",  lep1_p4.Phi());
        self.out.fillBranch("lep1_pdgid", leptons[0].pdgId );
        self.out.fillBranch("lep1_id",  lep1_id);
        self.out.fillBranch("lep1_iso",  lep1_iso);
        self.out.fillBranch("lep1_TrgObjfilterbits",  lep1_TrgObjfilterbits);
        self.out.fillBranch("lep1_E",  lep1_p4.E());
        self.out.fillBranch("lep2_pt",  lep2_p4.Pt());
        self.out.fillBranch("lep2_E",   lep2_p4.E());
        self.out.fillBranch("lep2_eta",  lep2_p4.Eta());
        self.out.fillBranch("lep2_phi",  lep2_p4.Phi());
        self.out.fillBranch("lep2_pdgid",  leptons[1].pdgId);
        self.out.fillBranch("lep2_id",  lep2_id);
        self.out.fillBranch("lep2_iso",  lep2_iso);
        self.out.fillBranch("lep2_TrgObjfilterbits",  lep2_TrgObjfilterbits);
	self.out.fillBranch("ht",  ht)
	self.out.fillBranch("ht_jets",  ht_jets)
	self.out.fillBranch("nJetsL",  nJetsL)
	self.out.fillBranch("jj_M",  jj_M)
	self.out.fillBranch("el_hltsafeid", llHLTsafeIDsf)
	self.out.fillBranch("llreco_weight",  lltrackingsf)
	self.out.fillBranch("llid_weight",  llIDsf)
	self.out.fillBranch("lliso_weight",  llIsosf)
	self.out.fillBranch("lltrigger_weight",  lltrgsf)
	self.out.fillBranch("ll_M",  ll_M)
	self.out.fillBranch("llmet_M",  llmet_M)
	#self.out.fillBranch("llmetjj_MT2",  "F")
	self.out.fillBranch("llmetjj_M", llmetjj_M)
	self.out.fillBranch("lljj_M",  lljj_M)
	self.out.fillBranch("cosThetaStar",  cosThetaStar)
	self.out.fillBranch("ll_DR_l_l",  ll_DR_l_l)
	self.out.fillBranch("jj_DR_j_j",  jj_DR_j_j)
	self.out.fillBranch("llmetjj_DPhi_ll_jj",  llmetjj_DPhi_ll_jj)
	self.out.fillBranch("llmetjj_DPhi_ll_met",  llmetjj_DPhi_ll_met)
	self.out.fillBranch("llmetjj_DPhi_llmet_jj",  llmetjj_DPhi_llmet_jj)
	self.out.fillBranch("ll_pt",  ll_pt)
	self.out.fillBranch("ll_eta",  ll_eta)
	self.out.fillBranch("jj_pt",  jj_pt)
	self.out.fillBranch("jj_eta",  jj_eta)
	self.out.fillBranch("llmetjj_minDR_l_j",  llmetjj_minDR_l_j)
	self.out.fillBranch("llmetjj_MTformula",  llmetjj_MTformula)
	self.out.fillBranch("ll_DPhi_l_l",  ll_DPhi_l_l)
	self.out.fillBranch("pu",  pu)
	self.out.fillBranch("event_pu_weight",  event_pu_weight)
	self.out.fillBranch("event_lep_weight",  event_lep_weight)
	self.out.fillBranch("event_btag_weight",  hJets_BtagSF[0] * hJets_BtagSF[1])
        self.out.fillBranch("sample_weight",  sample_weight)
        self.out.fillBranch("genweight",  genweight)
	self.out.fillBranch("event_reco_weight",  event_reco_weight)
	self.out.fillBranch("total_weight",  event_reco_weight * sample_weight)
	self.out.fillBranch("event_run",  run)
	self.out.fillBranch("event_lumiblock",  luminosityBlock)
	#self.out.fillBranch("event_weight",  "F")
	#self.out.fillBranch("DY_BDT_flat",  "F")
	#self.out.fillBranch("dy_nobtag_to_btagM_weight",  "F")
	#self.out.fillBranch("mt2",  "F")
	#self.out.fillBranch("mt2_bb",  "F")
	#self.out.fillBranch("mt2_ll",  "F")
	#self.out.fillBranch("event_number",  ievent)
	if self.isMC and self.addSystematic:
	    for name in self.ll_sys_branchname:
		self.out.fillBranch("lep1%s"%name,        self.lep1_branchname_sysvalue[name][0])
		self.out.fillBranch("lep1%s_up"%name,     self.lep1_branchname_sysvalue[name][1])
		self.out.fillBranch("lep1%s_down"%name,   self.lep1_branchname_sysvalue[name][2])
		self.out.fillBranch("lep2%s"%name,        self.lep2_branchname_sysvalue[name][0])
		self.out.fillBranch("lep2%s_up"%name,     self.lep2_branchname_sysvalue[name][1])
		self.out.fillBranch("lep2%s_down"%name,   self.lep2_branchname_sysvalue[name][2])
	    self.out.fillBranch("event_pu_weight_up",    event_pu_weight_up)
	    self.out.fillBranch("event_pu_weight_down",  event_pu_weight_down)
	    self.out.fillBranch("event_pdf_weight_up",   1.0 + LHEPdfweight_stddev)
	    self.out.fillBranch("event_pdf_weight_down", 1.0 - LHEPdfweight_stddev)
	    #self.out.fillBranch("LHEScaleWeight", LHEScaleWeight)
	    #print "len(LHEScaleWeight) ",len(LHEScaleWeight)
	    for i in range(0, len(LHEScaleWeight)):
		self.out.fillBranch("LHEScaleWeight_%d"%i, LHEScaleWeight[i])

	if self.isMC:
	    genmet = Object(event, "GenMET")
	    genParticles = Collection(event, "GenPart")
	    genjets = Collection(event, "GenJet")
	    nGenPart = getattr(event, "nGenPart", False)
	    nGenJet = getattr(event, "nGenJet", False)

	    nGenW1FromHiggs=0; nGenW2FromHiggs=0; nGenLepFromW1FromHiggs=0; nGenLepFromW2FromHiggs=0; nGenNuFromW1FromHiggs =0; nGenNuFromW2FromHiggs=0
	    nGenBQuarkFromHiggs=0
	    if hasattr(event, "nGenX"):
	        #nGenX   =  getattr(event, "nGenX", 0) 
		nGenW1FromHiggs = getattr(event, "nGenW1FromHiggs", 0) 
		nGenW2FromHiggs = getattr(event, "nGenW2FromHiggs", 0) 
		nGenLepFromW1FromHiggs = getattr(event, "nGenLepFromW1FromHiggs", 0)
		nGenNuFromW1FromHiggs = getattr(event, "nGenNuFromW1FromHiggs", 0)
		nGenLepFromW2FromHiggs = getattr(event, "nGenLepFromW2FromHiggs", 0)
		nGenNuFromW2FromHiggs = getattr(event, "nGenNuFromW2FromHiggs", 0)
		nGenBQuarkFromHiggs = getattr(event, "nGenBQuarkFromHiggs", 0)
	    if nGenW1FromHiggs  == 1 and nGenW2FromHiggs  == 1:
	        GenW1FromHiggs_mass = getattr(event, "GenW1FromHiggs_mass")
	        GenW1FromHiggs_pt = getattr(event, "GenW1FromHiggs_pt")
	        GenW2FromHiggs_mass = getattr(event, "GenW2FromHiggs_mass")
	        GenW2FromHiggs_pt = getattr(event, "GenW2FromHiggs_pt")
	        #self.out.fillBranch("genW1_pt",    GenW1FromHiggs_pt[0])
	        #self.out.fillBranch("genW1_mass",  GenW1FromHiggs_mass[0])
	        #self.out.fillBranch("genW2_pt",    GenW2FromHiggs_pt[0])
	        #self.out.fillBranch("genW2_mass",  GenW2FromHiggs_mass[0])
		if self.verbose > 3:
		  print "Gen W1 from Higgs, mass ", GenW1FromHiggs_mass[0], " pt ",GenW1FromHiggs_pt[0],\
		    "Gen W2 from Higgs, mass ",GenW2FromHiggs_mass[0], " pt ",GenW2FromHiggs_pt[0]
	    genNuFromW1_p4 = ROOT.TLorentzVector()
	    genNuFromW2_p4 = ROOT.TLorentzVector()
	    genLepFromW1_p4 = ROOT.TLorentzVector()
	    genLepFromW2_p4 = ROOT.TLorentzVector()
            dR_genl1_lepFromW1 = 99.0            
            dR_genl1_lepFromW2 = 99.0            
            dR_genl2_lepFromW1 = 99.0            
            dR_genl2_lepFromW2 = 99.0            
            foundLepsfromWs = False
	    if nGenLepFromW1FromHiggs == 1 and nGenNuFromW1FromHiggs == 1 and nGenLepFromW2FromHiggs == 1 and nGenNuFromW2FromHiggs == 1:
	        foundLepsfromWs = True
	        GenLepFromW1FromHiggs_pt   = getattr(event, "GenLepFromW1FromHiggs_pt")
	        GenLepFromW1FromHiggs_eta  = getattr(event, "GenLepFromW1FromHiggs_eta")
	        GenLepFromW1FromHiggs_phi  = getattr(event, "GenLepFromW1FromHiggs_phi")
	        GenNuFromW1FromHiggs_pt   = getattr(event, "GenNuFromW1FromHiggs_pt")
	        GenNuFromW1FromHiggs_eta  = getattr(event, "GenNuFromW1FromHiggs_eta")
	        GenNuFromW1FromHiggs_phi  = getattr(event, "GenNuFromW1FromHiggs_phi")
	        GenLepFromW2FromHiggs_pt   = getattr(event, "GenLepFromW2FromHiggs_pt")
	        GenLepFromW2FromHiggs_eta  = getattr(event, "GenLepFromW2FromHiggs_eta")
	        GenLepFromW2FromHiggs_phi  = getattr(event, "GenLepFromW2FromHiggs_phi")
	        GenNuFromW2FromHiggs_pt   = getattr(event, "GenNuFromW2FromHiggs_pt")
	        GenNuFromW2FromHiggs_eta  = getattr(event, "GenNuFromW2FromHiggs_eta")
	        GenNuFromW2FromHiggs_phi  = getattr(event, "GenNuFromW2FromHiggs_phi")
		genLepFromW1_p4.SetPtEtaPhiM(GenLepFromW1FromHiggs_pt[0], GenLepFromW1FromHiggs_eta[0], GenLepFromW1FromHiggs_phi[0], 0.0)
		genNuFromW1_p4.SetPtEtaPhiM(GenNuFromW1FromHiggs_pt[0], GenNuFromW1FromHiggs_eta[0], GenNuFromW1FromHiggs_phi[0], 0.0)
		genLepFromW2_p4.SetPtEtaPhiM(GenLepFromW2FromHiggs_pt[0], GenLepFromW2FromHiggs_eta[0], GenLepFromW2FromHiggs_phi[0], 0.0)
		genNuFromW2_p4.SetPtEtaPhiM(GenNuFromW2FromHiggs_pt[0], GenNuFromW2FromHiggs_eta[0], GenNuFromW2FromHiggs_phi[0], 0.0)
                diNu_p4 = genNuFromW1_p4+genNuFromW2_p4
		met_diNu_pt = diNu_p4.Pt()
		met_diNu_phi = diNu_p4.Phi()
		if  self.addGenToTree:
		    self.out.fillBranch("met_diNuetrino_pt",  met_diNu_pt);
		    self.out.fillBranch("met_diNuetrino_phi",  met_diNu_phi);
		    #self.out.fillBranch("gennu1_pt",  GenNuFromW1FromHiggs_pt[0]) #nu1 from W1
		    #self.out.fillBranch("gennu1_eta", GenNuFromW1FromHiggs_eta[0]) 
		    #self.out.fillBranch("gennu1_phi", GenNuFromW1FromHiggs_phi[0]) 
		    #self.out.fillBranch("gennu2_pt",  GenNuFromW2FromHiggs_pt[0]) 
		    #self.out.fillBranch("gennu2_eta", GenNuFromW2FromHiggs_eta[0]) 
		    #self.out.fillBranch("gennu2_phi", GenNuFromW2FromHiggs_phi[0]) 
		if self.verbose > 3:
                  print "GenLepFromW1 pt ",GenLepFromW1FromHiggs_pt[0], " eta ",GenLepFromW1FromHiggs_eta[0], " phi ",GenLepFromW1FromHiggs_phi[0],\
		      " GenLepFromW2 pt ",GenLepFromW2FromHiggs_pt[0], " eta ",GenLepFromW2FromHiggs_eta[0], " phi ",GenLepFromW2FromHiggs_phi[0]
		  print "MET from two nuetrinos pt ",met_diNu_pt, " phi ",met_diNu_phi
		  print "genmet_pt ",genmet.pt," phi ",genmet.phi

	    genb1FromHiggs_p4 = ROOT.TLorentzVector()
	    genb2FromHiggs_p4 = ROOT.TLorentzVector()
            dR_genb1_genjet1 = 99.0
            dR_genb1_genjet2 = 99.0
            dR_genb2_genjet1 = 99.0
            dR_genb2_genjet2 = 99.0
            foundBQuarkFromHiggs = False
	    if nGenBQuarkFromHiggs == 2 :
	        foundBQuarkFromHiggs = True
	        GenBQuarkFromHiggs_pt = getattr(event, "GenBQuarkFromHiggs_pt")
	        GenBQuarkFromHiggs_mass = getattr(event, "GenBQuarkFromHiggs_mass")
	        GenBQuarkFromHiggs_eta = getattr(event, "GenBQuarkFromHiggs_eta")
	        GenBQuarkFromHiggs_phi = getattr(event, "GenBQuarkFromHiggs_phi")
		genb1FromHiggs_p4.SetPtEtaPhiM(GenBQuarkFromHiggs_pt[0], GenBQuarkFromHiggs_eta[0], GenBQuarkFromHiggs_phi[0], GenBQuarkFromHiggs_mass[0])
		genb2FromHiggs_p4.SetPtEtaPhiM(GenBQuarkFromHiggs_pt[1], GenBQuarkFromHiggs_eta[1], GenBQuarkFromHiggs_phi[1], GenBQuarkFromHiggs_mass[1])
		if  self.addGenToTree:
		    self.out.fillBranch("genb1_pt",  GenBQuarkFromHiggs_pt[0]) #b1
		    self.out.fillBranch("genb2_pt",  GenBQuarkFromHiggs_pt[1]) #b2
		    #self.out.fillBranch("genb1_eta",  GenBQuarkFromHiggs_eta[0]) #b1
		    #self.out.fillBranch("genb1_phi",  GenBQuarkFromHiggs_phi[0]) #b1
		    #self.out.fillBranch("genb2_eta",  GenBQuarkFromHiggs_eta[1]) #b2
		    #self.out.fillBranch("genb2_phi",  GenBQuarkFromHiggs_phi[1]) #b2
		if self.verbose > 3:
		    print "GEn B quark from Higgs pt ",GenBQuarkFromHiggs_pt, " eta ",GenBQuarkFromHiggs_eta, " phi ",GenBQuarkFromHiggs_phi
	    
	    lep1_genindex = leptons[0].genPartIdx
	    lep2_genindex = leptons[1].genPartIdx
	    if lep1_genindex >= 0 and lep1_genindex < nGenPart and self.addGenToTree:
		genl1 = genParticles[lep1_genindex]
		if self.verbose > 3:
		    print "lep1 pt ",genl1.pt," eta ",genl1.eta, " phi ",genl1.phi
		genl1_p4 = ROOT.TLorentzVector(); genl1_p4.SetPtEtaPhiM(genl1.pt, genl1.eta, genl1.phi, genl1.mass)
		self.out.fillBranch("genl1_pt",  genl1_p4.Pt());
		self.out.fillBranch("genl1_E",  genl1_p4.E());
		self.out.fillBranch("genl1_eta",  genl1_p4.Eta());
		self.out.fillBranch("genl1_phi",  genl1_p4.Phi());
		self.out.fillBranch("genl1_id",  abs(genl1.pdgId));
		if foundLepsfromWs:
		    dR_genl1_lepFromW1 = genl1_p4.DeltaR(genLepFromW1_p4)
		    dR_genl1_lepFromW2 = genl1_p4.DeltaR(genLepFromW2_p4)
	    if lep2_genindex >= 0 and lep2_genindex < nGenPart and self.addGenToTree:
		genl2 = genParticles[lep2_genindex]
		if self.verbose > 3:
		    print "lep2 pt ",genl2.pt," eta ",genl2.eta, " phi ",genl2.phi
		genl2_p4 = ROOT.TLorentzVector(); genl2_p4.SetPtEtaPhiM(genl2.pt, genl2.eta, genl2.phi, genl2.mass)
		self.out.fillBranch("genl2_pt",  genl2_p4.Pt());
		self.out.fillBranch("genl2_E",  genl2_p4.E());
		self.out.fillBranch("genl2_eta",  genl2_p4.Eta());
		self.out.fillBranch("genl2_phi",  genl2_p4.Phi());
		self.out.fillBranch("genl2_id",  abs(genl2.pdgId));
		if foundLepsfromWs:
		    dR_genl2_lepFromW1 = genl2_p4.DeltaR(genLepFromW1_p4)
		    dR_genl2_lepFromW2 = genl2_p4.DeltaR(genLepFromW2_p4)
	    ## jet matching  
            genjet1_partonflavour = 0; genjet2_partonflavour = 0
	    jjbtag_heavy = 1.0; jjbtag_light = 1.0
    	    if jet1.genJetIdx >= 0 and jet1.genJetIdx < nGenJet:
		genjet1 = genjets[jet1.genJetIdx]
		genjet1_p4 = ROOT.TLorentzVector(); 
		if self.verbose > 3:
		    print "genjet1 pt ",genjet1.pt," eta ",genjet1.eta, " phi ",genjet1.phi
		genjet1_p4.SetPtEtaPhiM(genjet1.pt, genjet1.eta, genjet1.phi, genjet1.mass)
		genjet1_partonflavour = flavour(genjet1.partonFlavour)
    		if self.addGenToTree:
		    self.out.fillBranch("genjet1_pt", genjet1_p4.Pt());
		    self.out.fillBranch("genjet1_E",  genjet1_p4.E());
		    self.out.fillBranch("genjet1_eta",  genjet1_p4.Eta());
		    self.out.fillBranch("genjet1_phi",  genjet1_p4.Phi());
		if self.DYestimation or self.addGenToTree:	
		    self.out.fillBranch("genjet1_partonFlavour",  genjet1_partonflavour);
		if foundBQuarkFromHiggs:
		    dR_genb1_genjet1 = genjet1_p4.DeltaR(genb1FromHiggs_p4)
		    dR_genb2_genjet1 = genjet1_p4.DeltaR(genb2FromHiggs_p4)
		    
    	    if jet2.genJetIdx >= 0 and jet2.genJetIdx < nGenJet :
		genjet2 = genjets[jet2.genJetIdx]
		genjet2_p4 = ROOT.TLorentzVector(); 
		if self.verbose > 3:
		    print "genjet2 pt ",genjet2.pt," eta ",genjet2.eta, " phi ",genjet2.phi
		genjet2_partonflavour = flavour(genjet2.partonFlavour)
		genjet2_p4.SetPtEtaPhiM(genjet2.pt, genjet2.eta, genjet2.phi, genjet2.mass)
    		if self.addGenToTree:
		    self.out.fillBranch("genjet2_pt", genjet2_p4.Pt());
		    self.out.fillBranch("genjet2_E",  genjet2_p4.E());
		    self.out.fillBranch("genjet2_eta",  genjet2_p4.Eta());
		    self.out.fillBranch("genjet2_phi",  genjet2_p4.Phi());
		if self.DYestimation or self.addGenToTree:	
		    self.out.fillBranch("genjet2_partonFlavour",  genjet2_partonflavour);

		if foundBQuarkFromHiggs:
		    dR_genb1_genjet2 = genjet2_p4.DeltaR(genb1FromHiggs_p4)
		    dR_genb2_genjet2 = genjet2_p4.DeltaR(genb2FromHiggs_p4) 
	    #self.jjbtag_sys_branchname = ["jjbtag_light","jjbtag_heavy"]
            self.jjbtag_branchname_sysvalue["jjbtag_heavy"] = [1.0, 1.0, 1.0]
            self.jjbtag_branchname_sysvalue["jjbtag_light"] = [1.0, 1.0, 1.0]
	    if genjet1_partonflavour > 0:
	        jjbtag_heavy = jjbtag_heavy * hJets_BtagSF[0] 
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][0] = jjbtag_heavy
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][1] = self.jjbtag_branchname_sysvalue["jjbtag_heavy"][1]* hJets_BtagSF_up[0]/hJets_BtagSF[0]
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][2] = self.jjbtag_branchname_sysvalue["jjbtag_heavy"][2]* hJets_BtagSF_down[0]/hJets_BtagSF[0]
	    else:
	        jjbtag_light = jjbtag_light * hJets_BtagSF[0] 
		self.jjbtag_branchname_sysvalue["jjbtag_light"][0] = jjbtag_light
		self.jjbtag_branchname_sysvalue["jjbtag_light"][1] = self.jjbtag_branchname_sysvalue["jjbtag_light"][1]* hJets_BtagSF_up[0]/hJets_BtagSF[0]
		self.jjbtag_branchname_sysvalue["jjbtag_light"][2] = self.jjbtag_branchname_sysvalue["jjbtag_light"][2]* hJets_BtagSF_down[0]/hJets_BtagSF[0]
	    if genjet2_partonflavour > 0:
	        jjbtag_heavy = jjbtag_heavy * hJets_BtagSF[1] 
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][0] = jjbtag_heavy
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][1] = self.jjbtag_branchname_sysvalue["jjbtag_heavy"][1]* hJets_BtagSF_up[1]/hJets_BtagSF[1]
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][2] = self.jjbtag_branchname_sysvalue["jjbtag_heavy"][2]* hJets_BtagSF_down[1]/hJets_BtagSF[1]
	    else:
	        jjbtag_light = jjbtag_light * hJets_BtagSF[1] 
		self.jjbtag_branchname_sysvalue["jjbtag_light"][0] = jjbtag_light
		self.jjbtag_branchname_sysvalue["jjbtag_light"][1] = self.jjbtag_branchname_sysvalue["jjbtag_light"][1]* hJets_BtagSF_up[1]/hJets_BtagSF[1]
		self.jjbtag_branchname_sysvalue["jjbtag_light"][2] = self.jjbtag_branchname_sysvalue["jjbtag_light"][2]* hJets_BtagSF_down[1]/hJets_BtagSF[1]
	    self.out.fillBranch("jjbtag_heavy", jjbtag_heavy);
	    self.out.fillBranch("jjbtag_light", jjbtag_light);
	    if self.addSystematic:
		for name in self.jjbtag_branchname_sysvalue.keys():
		    self.out.fillBranch(name,   self.jjbtag_branchname_sysvalue[name][0])
		    self.out.fillBranch(name+"_up",   self.jjbtag_branchname_sysvalue[name][1])
		    self.out.fillBranch(name+"_down", self.jjbtag_branchname_sysvalue[name][2])


	    #genflavour = lambda jet : (jet.genJetIdx < nGenJet and jet.genJetIdx >=0)*flavour(genjets[abs(jet.genJetIdx)].partonFlavour)
            if self.CheckBtaggingEff:
		alljets_genpartonFlavour = []
		for jet in alljets:
		    if jet.genJetIdx < nGenJet and jet.genJetIdx >=0:
			alljets_genpartonFlavour.append(flavour(genjets[abs(jet.genJetIdx)].partonFlavour))
		    else:
			alljets_genpartonFlavour.append(0)
			
		alljets_genpartonFlavour = resize_alljets(alljets_genpartonFlavour)
		alljets_partonFlavour = [jet.partonFlavour for jet in alljets]
		alljets_hadronFlavour = [jet.hadronFlavour for jet in alljets]
		alljets_partonFlavour = resize_alljets(alljets_partonFlavour)
		alljets_hadronFlavour = resize_alljets(alljets_hadronFlavour)
		self.out.fillBranch("alljets_pt", alljets_pt)
		self.out.fillBranch("alljets_eta", alljets_eta)
		self.out.fillBranch("alljets_cMVAv2", alljets_cMVAv2)
		self.out.fillBranch("alljets_DeepCSV", alljets_DeepCSV)
		self.out.fillBranch("alljets_genpartonFlavour", alljets_genpartonFlavour)
		#self.out.fillBranch("alljets_partonFlavour", alljets_partonFlavour)
		#self.out.fillBranch("alljets_hadronFlavour", alljets_hadronFlavour)

	    if self.DYestimation or self.addGenToTree:	
		self.out.fillBranch("jet1_partonFlavour",  flavour(jet1.partonFlavour));
		self.out.fillBranch("jet1_hadronFlavour",  flavour(jet1.hadronFlavour));
		self.out.fillBranch("jet2_partonFlavour",  flavour(jet2.partonFlavour));
		self.out.fillBranch("jet2_hadronFlavour",  flavour(jet2.hadronFlavour));
		self.out.fillBranch("genmet_pt",  genmet.pt);
		self.out.fillBranch("genmet_phi", genmet.phi);

    	    if  self.addGenToTree:
		self.out.fillBranch("dR_genb1_genjet1",  dR_genb1_genjet1);
		self.out.fillBranch("dR_genb1_genjet2",  dR_genb1_genjet2);
		self.out.fillBranch("dR_genb2_genjet1",  dR_genb2_genjet1);
		self.out.fillBranch("dR_genb2_genjet2",  dR_genb2_genjet2);
		self.out.fillBranch("dR_genl1_lepFromW1",  dR_genl1_lepFromW1);
		self.out.fillBranch("dR_genl1_lepFromW2",  dR_genl1_lepFromW2);
		self.out.fillBranch("dR_genl2_lepFromW1",  dR_genl2_lepFromW1);
		self.out.fillBranch("dR_genl2_lepFromW2",  dR_genl2_lepFromW2);


	"""
	def fillBranches_muon(leptons):
	    if (len(leptons) >= 1):
		lep1 = leptons[0]
		lep1_p4 = ROOT.TLorentzVector()
		lep1_p4.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, leptons[0].mass)
		self.out.fillBranch("mu1_pt",                     lep1.pt);
		self.out.fillBranch("mu1_E",                      lep1_p4.E());
		self.out.fillBranch("mu1_eta",                    lep1.eta);
		self.out.fillBranch("mu1_phi",                    lep1.phi);
		self.out.fillBranch("mu1_pdgid",                  lep1.pdgId);
		self.out.fillBranch("mu1_charge",                 lep1.charge);
		self.out.fillBranch("mu1_sip3D",                  lep1.sip3d);
		self.out.fillBranch("mu1_miniRelIso",             lep1.miniPFRelIso_all);
		self.out.fillBranch("mu1_miniRelIsoCharged",      lep1.miniPFRelIso_chg);
		self.out.fillBranch("mu1_dxy",                    lep1.dxy);
		self.out.fillBranch("mu1_dz",                     lep1.dz);
		self.out.fillBranch("mu1_segmentCompatibility",   lep1.segmentComp);
		self.out.fillBranch("mu1_leptonMVA",              lep1.mvaTTH);
		self.out.fillBranch("mu1_dxyAbs",                 abs(lep1.dxy));
	    else:
		self.out.fillBranch("mu1_pt",                     -10000.0);
		self.out.fillBranch("mu1_E",                      -10000.0);
		self.out.fillBranch("mu1_eta",                    -10000.0);
		self.out.fillBranch("mu1_phi",                    -10000.0);
		self.out.fillBranch("mu1_pdgid",                  -100000);
		self.out.fillBranch("mu1_charge",                 -10000.0);
		self.out.fillBranch("mu1_sip3D",                  -10000.0);
		self.out.fillBranch("mu1_miniRelIso",             -10000.0);
		self.out.fillBranch("mu1_miniRelIsoCharged",      -10000.0);
		self.out.fillBranch("mu1_dxy",                    -10000.0);
		self.out.fillBranch("mu1_dxyAbs",                 -10000.0);
		self.out.fillBranch("mu1_dz",                     -10000.0);
		self.out.fillBranch("mu1_segmentCompatibility",   -10000.0);
		self.out.fillBranch("mu1_leptonMVA",              -10000.0);
	    
	    if (len(leptons) >= 2):
		lep2 = leptons[1]
		lep2_p4 = ROOT.TLorentzVector()
		lep2_p4.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, leptons[1].mass)
		self.out.fillBranch("mu2_pt",                     lep2.pt);
		self.out.fillBranch("mu2_E",                      lep2_p4.E());
		self.out.fillBranch("mu2_eta",                    lep2.eta);
		self.out.fillBranch("mu2_phi",                    lep2.phi);
		self.out.fillBranch("mu2_pdgid",                  lep2.pdgId);
		self.out.fillBranch("mu2_charge",                 lep2.charge);
		self.out.fillBranch("mu2_sip3D",                  lep2.sip3d);
		self.out.fillBranch("mu2_miniRelIso",             lep2.miniPFRelIso_all);
		self.out.fillBranch("mu2_miniRelIsoCharged",      lep2.miniPFRelIso_chg);
		self.out.fillBranch("mu2_dxy",                    lep2.dxy);
		self.out.fillBranch("mu2_dz",                     lep2.dz);
		self.out.fillBranch("mu2_segmentCompatibility",   lep2.segmentComp);
		self.out.fillBranch("mu2_leptonMVA",              lep2.mvaTTH);
		self.out.fillBranch("mu2_dxyAbs",                 abs(lep2.dxy));
	    else:
		self.out.fillBranch("mu2_pt",                     -10000.0);
		self.out.fillBranch("mu2_E",                      -10000.0);
		self.out.fillBranch("mu2_eta",                    -10000.0);
		self.out.fillBranch("mu2_phi",                    -10000.0);
		self.out.fillBranch("mu2_pdgid",                  -100000);
		self.out.fillBranch("mu2_charge",                 -10000.0);
		self.out.fillBranch("mu2_sip3D",                  -10000.0);
		self.out.fillBranch("mu2_miniRelIso",             -10000.0);
		self.out.fillBranch("mu2_miniRelIsoCharged",      -10000.0);
		self.out.fillBranch("mu2_dxy",                    -10000.0);
		self.out.fillBranch("mu2_dxyAbs",                 -10000.0);
		self.out.fillBranch("mu2_dz",                     -10000.0);
		self.out.fillBranch("mu2_segmentCompatibility",   -10000.0);
		self.out.fillBranch("mu2_leptonMVA",              -10000.0);

	def fillBranches_electron(leptons):
	    if (len(leptons) >= 1):
		lep1 = leptons[0]
		lep1_p4 = ROOT.TLorentzVector()
		lep1_p4.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, leptons[0].mass)
		self.out.fillBranch("ele1_pt",                    lep1.pt);
		self.out.fillBranch("ele1_E",                     lep1_p4.E());
		self.out.fillBranch("ele1_eta",                   lep1.eta);
		self.out.fillBranch("ele1_phi",                   lep1.phi);
		self.out.fillBranch("ele1_pdgid",                 lep1.pdgId);
		self.out.fillBranch("ele1_charge",                lep1.charge);
		self.out.fillBranch("ele1_sip3D",                 lep1.sip3d);
		self.out.fillBranch("ele1_miniRelIso",            lep1.miniPFRelIso_all);
		self.out.fillBranch("ele1_miniRelIsoCharged",     lep1.miniPFRelIso_chg);
		self.out.fillBranch("ele1_dxy",                   lep1.dxy);
		self.out.fillBranch("ele1_dz",                    lep1.dz);
		self.out.fillBranch("ele1_leptonMVA",             lep1.mvaTTH);
		self.out.fillBranch("ele1_dxyAbs",                abs(lep1.dxy));
		self.out.fillBranch("ele1_ntMVAeleID",            0.0);##not available ?
	    else:
		self.out.fillBranch("ele1_pt",                    -10000.0);
		self.out.fillBranch("ele1_E",                     -10000.0);
		self.out.fillBranch("ele1_eta",                   -10000.0);
		self.out.fillBranch("ele1_phi",                   -10000.0);
		self.out.fillBranch("ele1_pdgid",                 -100000);
		self.out.fillBranch("ele1_charge",                -10000.0);
		self.out.fillBranch("ele1_sip3D",                 -10000.0);
		self.out.fillBranch("ele1_miniRelIso",            -10000.0);
		self.out.fillBranch("ele1_miniRelIsoCharged",     -10000.0);
		self.out.fillBranch("ele1_dxy",                   -10000.0);
		self.out.fillBranch("ele1_dxyAbs",                -10000.0);
		self.out.fillBranch("ele1_dz",                    -10000.0);
		self.out.fillBranch("ele1_ntMVAeleID",            -10000.0);
		self.out.fillBranch("ele1_leptonMVA",             -10000.0);
	    
	    if (len(leptons) >= 2):
		lep2 = leptons[1]
		lep2_p4 = ROOT.TLorentzVector()
		lep2_p4.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, leptons[1].mass)
		self.out.fillBranch("ele2_pt",                    lep2.pt);
		self.out.fillBranch("ele2_E",                     lep2_p4.E());
		self.out.fillBranch("ele2_eta",                   lep2.eta);
		self.out.fillBranch("ele2_phi",                   lep2.phi);
		self.out.fillBranch("ele2_pdgid",                 lep2.pdgId);
		self.out.fillBranch("ele2_charge",                lep2.charge);
		self.out.fillBranch("ele2_sip3D",                 lep2.sip3d);
		self.out.fillBranch("ele2_miniRelIso",            lep2.miniPFRelIso_all);
		self.out.fillBranch("ele2_miniRelIsoCharged",     lep2.miniPFRelIso_chg);
		self.out.fillBranch("ele2_dxy",                   lep2.dxy);
		self.out.fillBranch("ele2_dz",                    lep2.dz);
		self.out.fillBranch("ele2_leptonMVA",             lep2.mvaTTH);
		self.out.fillBranch("ele2_dxyAbs",                abs(lep2.dxy));
		self.out.fillBranch("ele2_ntMVAeleID",            0.0);##not available ?
	    else:
		self.out.fillBranch("ele2_pt",                    -10000.0);
		self.out.fillBranch("ele2_E",                     -10000.0);
		self.out.fillBranch("ele2_eta",                   -10000.0);
		self.out.fillBranch("ele2_phi",                   -10000.0);
		self.out.fillBranch("ele2_pdgid",                 -100000);
		self.out.fillBranch("ele2_charge",                -10000.0);
		self.out.fillBranch("ele2_sip3D",                 -10000.0);
		self.out.fillBranch("ele2_miniRelIso",            -10000.0);
		self.out.fillBranch("ele2_miniRelIsoCharged",     -10000.0);
		self.out.fillBranch("ele2_dxy",                   -10000.0);
		self.out.fillBranch("ele2_dxyAbs",                -10000.0);
		self.out.fillBranch("ele2_dz",                    -10000.0);
		self.out.fillBranch("ele2_ntMVAeleID",            -10000.0);
		self.out.fillBranch("ele2_leptonMVA",             -10000.0);

	def fillBranches_jet(jets):
            if (len(jets) >= 1):
		jet1 = jets[0]
		jet1_p4 = ROOT.TLorentzVector()
		jet1_p4.SetPtEtaPhiM(jets[0].pt, jets[0].eta, jets[0].phi, jets[0].mass)
		self.out.fillBranch("ak4Jet1_pt",                 jet1.pt);
		self.out.fillBranch("ak4Jet1_E",                  jet1_p4.E());
		self.out.fillBranch("ak4Jet1_eta",                jet1.eta);
		self.out.fillBranch("ak4Jet1_phi",                jet1.phi);
		self.out.fillBranch("ak4Jet1_cMVAv2",             jet1.btagCMVA);
		#self.out.fillBranch("ak4Jet1_CSV",                jet1.btagDeepB);
		self.out.fillBranch("ak4Jet1_CSV",                jet1.btagDeepFlavB);
	    else:
		self.out.fillBranch("ak4Jet1_pt",                 -10000.0);
		self.out.fillBranch("ak4Jet1_E",                  -10000.0);
		self.out.fillBranch("ak4Jet1_eta",                -10000.0);
		self.out.fillBranch("ak4Jet1_phi",                -10000.0);
		self.out.fillBranch("ak4Jet1_cMVAv2",             -10000.0);
		self.out.fillBranch("ak4Jet1_CSV",                -10000.0);

            if (len(jets) >= 2):
		jet2 = jets[1]
		jet2_p4 = ROOT.TLorentzVector()
		jet2_p4.SetPtEtaPhiM(jets[1].pt, jets[1].eta, jets[1].phi, jets[1].mass)
		self.out.fillBranch("ak4Jet2_pt",                 jet2.pt);
		self.out.fillBranch("ak4Jet2_E",                  jet2_p4.E());
		self.out.fillBranch("ak4Jet2_eta",                jet2.eta);
		self.out.fillBranch("ak4Jet2_phi",                jet2.phi);
		self.out.fillBranch("ak4Jet2_cMVAv2",             jet2.btagCMVA);
	       #self.out.fillBranch("ak4Jet2_CSV",                jet2.btagDeepB);
		self.out.fillBranch("ak4Jet2_CSV",                jet2.btagDeepFlavB);
	    else:
		self.out.fillBranch("ak4Jet2_pt",                 -10000.0);
		self.out.fillBranch("ak4Jet2_E",                  -10000.0);
		self.out.fillBranch("ak4Jet2_eta",                -10000.0);
		self.out.fillBranch("ak4Jet2_phi",                -10000.0);
		self.out.fillBranch("ak4Jet2_cMVAv2",             -10000.0);
		self.out.fillBranch("ak4Jet2_CSV",                -10000.0);

            if (len(jets) >= 3):
		jet3 = jets[2]
		jet3_p4 = ROOT.TLorentzVector()
		jet3_p4.SetPtEtaPhiM(jets[2].pt, jets[2].eta, jets[2].phi, jets[2].mass)
		self.out.fillBranch("ak4Jet3_pt",                 jet3.pt);
		self.out.fillBranch("ak4Jet3_E",                  jet3_p4.E());
		self.out.fillBranch("ak4Jet3_eta",                jet3.eta);
		self.out.fillBranch("ak4Jet3_phi",                jet3.phi);
		self.out.fillBranch("ak4Jet3_cMVAv2",             jet3.btagCMVA);
	       #self.out.fillBranch("ak4Jet3_CSV",                jet3.btagDeepB);
		self.out.fillBranch("ak4Jet3_CSV",                jet3.btagDeepFlavB);
	    else:
		self.out.fillBranch("ak4Jet3_pt",                 -10000.0);
		self.out.fillBranch("ak4Jet3_E",                  -10000.0);
		self.out.fillBranch("ak4Jet3_eta",                -10000.0);
		self.out.fillBranch("ak4Jet3_phi",                -10000.0);
		self.out.fillBranch("ak4Jet3_cMVAv2",             -10000.0);
		self.out.fillBranch("ak4Jet3_CSV",                -10000.0);
            if (len(jets) >= 4):
		jet4 = jets[3]
		jet4_p4 = ROOT.TLorentzVector()
		jet4_p4.SetPtEtaPhiM(jets[3].pt, jets[3].eta, jets[3].phi, jets[3].mass)
		self.out.fillBranch("ak4Jet4_pt",                 jet4.pt);
		self.out.fillBranch("ak4Jet4_E",                  jet4_p4.E());
		self.out.fillBranch("ak4Jet4_eta",                jet4.eta);
		self.out.fillBranch("ak4Jet4_phi",                jet4.phi);
		self.out.fillBranch("ak4Jet4_cMVAv2",             jet4.btagCMVA);
	       #self.out.fillBranch("ak4Jet4_CSV",                jet4.btagDeepB);
		self.out.fillBranch("ak4Jet4_CSV",                jet4.btagDeepFlavB);
	    else:
		self.out.fillBranch("ak4Jet4_pt",                 -10000.0);
		self.out.fillBranch("ak4Jet4_E",                  -10000.0);
		self.out.fillBranch("ak4Jet4_eta",                -10000.0);
		self.out.fillBranch("ak4Jet4_phi",                -10000.0);
		self.out.fillBranch("ak4Jet4_cMVAv2",             -10000.0);
		self.out.fillBranch("ak4Jet4_CSV",                -10000.0);



	def fillBranches_ak8jet(ak8jets, ak8subjets):
            if (len(ak8jets) >= 1):
		ak8jet1 = ak8jets[0]
		ak8jet1_p4 = ROOT.TLorentzVector()
		ak8jet1_p4.SetPtEtaPhiM(ak8jets[0].pt, ak8jets[0].eta, ak8jets[0].phi, ak8jets[0].mass)
		subjet1 = ak8subjets[ak8jet1.subJetIdx1]
		subjet2 = ak8subjets[ak8jet1.subJetIdx2]
		self.out.fillBranch("ak8Jet1_pt",                 ak8jet1.pt);
		self.out.fillBranch("ak8Jet1_E",                  ak8jet1_p4.E());
		self.out.fillBranch("ak8Jet1_eta",                ak8jet1.eta);
		self.out.fillBranch("ak8Jet1_phi",                ak8jet1.phi);
		self.out.fillBranch("ak8Jet1_msoftdrop",          ak8jet1.msoftdrop);
		self.out.fillBranch("ak8Jet1_btagHbb",            ak8jet1.btagHbb);
		self.out.fillBranch("ak8Jet1_tau1",               ak8jet1.tau1);
		self.out.fillBranch("ak8Jet1_tau2",               ak8jet1.tau2);
		self.out.fillBranch("ak8Jet1_subjet0_pt",         subjet1.pt);
		self.out.fillBranch("ak8Jet1_subjet0_eta",        subjet1.eta);
		self.out.fillBranch("ak8Jet1_subjet0_phi",        subjet1.phi);
		self.out.fillBranch("ak8Jet1_subjet0_CSV",        subjet1.btagDeepB);   
		#self.out.fillBranch("ak8Jet1_subjet0_CSV",        subjet1.btagDeepFlavB);   
		self.out.fillBranch("ak8Jet1_subjet1_pt",         subjet2.pt);
		self.out.fillBranch("ak8Jet1_subjet1_eta",        subjet2.eta);
		self.out.fillBranch("ak8Jet1_subjet1_phi",        subjet2.phi);
		self.out.fillBranch("ak8Jet1_subjet1_CSV",        subjet2.btagDeepB);
		#self.out.fillBranch("ak8Jet1_subjet1_CSV",        subjet2.btagDeepFlavB);
	    else:
		self.out.fillBranch("ak8Jet1_pt",                 -10000.0);
		self.out.fillBranch("ak8Jet1_E",                  -10000.0);
		self.out.fillBranch("ak8Jet1_eta",                -10000.0);
		self.out.fillBranch("ak8Jet1_phi",                -10000.0);
		self.out.fillBranch("ak8Jet1_msoftdrop",          -10000.0);
		self.out.fillBranch("ak8Jet1_btagHbb",            -10000.0);
		self.out.fillBranch("ak8Jet1_tau1",               -10000.0);
		self.out.fillBranch("ak8Jet1_tau2",               -10000.0);
		self.out.fillBranch("ak8Jet1_subjet0_pt",         -10000.0);
		self.out.fillBranch("ak8Jet1_subjet0_eta",        -10000.0);
		self.out.fillBranch("ak8Jet1_subjet0_phi",        -10000.0);
		self.out.fillBranch("ak8Jet1_subjet0_CSV",        -10000.0);   
		self.out.fillBranch("ak8Jet1_subjet1_pt",         -10000.0);
		self.out.fillBranch("ak8Jet1_subjet1_eta",        -10000.0);
		self.out.fillBranch("ak8Jet1_subjet1_phi",        -10000.0);
		self.out.fillBranch("ak8Jet1_subjet1_CSV",        -10000.0);   
	    
            if (len(ak8jets) >= 2):
		ak8jet2 = ak8jets[1]
		ak8jet2_p4 = ROOT.TLorentzVector()
		ak8jet2_p4.SetPtEtaPhiM(ak8jets[1].pt, ak8jets[1].eta, ak8jets[1].phi, ak8jets[1].mass)
		subjet1 = ak8subjets[ak8jet2.subJetIdx1]
		subjet2 = ak8subjets[ak8jet2.subJetIdx2]
		self.out.fillBranch("ak8Jet2_pt",                 ak8jet2.pt);
		self.out.fillBranch("ak8Jet2_E",                  ak8jet2_p4.E());
		self.out.fillBranch("ak8Jet2_eta",                ak8jet2.eta);
		self.out.fillBranch("ak8Jet2_phi",                ak8jet2.phi);
		self.out.fillBranch("ak8Jet2_msoftdrop",          ak8jet2.msoftdrop);
		self.out.fillBranch("ak8Jet2_btagHbb",            ak8jet2.btagHbb);
		self.out.fillBranch("ak8Jet2_tau1",               ak8jet2.tau1);
		self.out.fillBranch("ak8Jet2_tau2",               ak8jet2.tau2);
		self.out.fillBranch("ak8Jet2_subjet0_pt",         subjet1.pt);
		self.out.fillBranch("ak8Jet2_subjet0_eta",        subjet1.eta);
		self.out.fillBranch("ak8Jet2_subjet0_phi",        subjet1.phi);
		self.out.fillBranch("ak8Jet2_subjet0_CSV",        subjet1.btagDeepB);   
		#self.out.fillBranch("ak8Jet2_subjet0_CSV",        subjet1.btagDeepFlavB);   
		self.out.fillBranch("ak8Jet2_subjet1_pt",         subjet2.pt);
		self.out.fillBranch("ak8Jet2_subjet1_eta",        subjet2.eta);
		self.out.fillBranch("ak8Jet2_subjet1_phi",        subjet2.phi);
		self.out.fillBranch("ak8Jet2_subjet1_CSV",        subjet2.btagDeepB);
		#self.out.fillBranch("ak8Jet2_subjet1_CSV",        subjet2.btagDeepFlavB);
	    else:
		self.out.fillBranch("ak8Jet2_pt",                 -10000.0);
		self.out.fillBranch("ak8Jet2_E",                  -10000.0);
		self.out.fillBranch("ak8Jet2_eta",                -10000.0);
		self.out.fillBranch("ak8Jet2_phi",                -10000.0);
		self.out.fillBranch("ak8Jet2_msoftdrop",          -10000.0);
		self.out.fillBranch("ak8Jet2_tau1",               -10000.0);
		self.out.fillBranch("ak8Jet2_tau2",               -10000.0);
		self.out.fillBranch("ak8Jet2_btagHbb",            -10000.0);
		self.out.fillBranch("ak8Jet2_subjet0_pt",         -10000.0);
		self.out.fillBranch("ak8Jet2_subjet0_eta",        -10000.0);
		self.out.fillBranch("ak8Jet2_subjet0_phi",        -10000.0);
		self.out.fillBranch("ak8Jet2_subjet0_CSV",        -10000.0);   
		self.out.fillBranch("ak8Jet2_subjet1_pt",         -10000.0);
		self.out.fillBranch("ak8Jet2_subjet1_eta",        -10000.0);
		self.out.fillBranch("ak8Jet2_subjet1_phi",        -10000.0);
		self.out.fillBranch("ak8Jet2_subjet1_CSV",        -10000.0);   
	    
	def fillBranches_ak8lsjet(ak8lsjets, ak8lssubjets):
            if (len(ak8lsjets) >= 1):
		ak8lsjet1 = ak8lsjets[0]
		ak8lsjet1_p4 = ROOT.TLorentzVector()
		ak8lsjet1_p4.SetPtEtaPhiM(ak8lsjets[0].pt, ak8lsjets[0].eta, ak8lsjets[0].phi, ak8lsjets[0].mass)
		subjet1 = ak8lssubjets[ak8lsjet1.subJetIdx1]
		subjet2 = ak8lssubjets[ak8lsjet1.subJetIdx2]
		self.out.fillBranch("ak8lsJet1_pt",                 ak8lsjet1.pt);
		self.out.fillBranch("ak8lsJet1_E",                  ak8lsjet1_p4.E());
		self.out.fillBranch("ak8lsJet1_eta",                ak8lsjet1.eta);
		self.out.fillBranch("ak8lsJet1_phi",                ak8lsjet1.phi);
		self.out.fillBranch("ak8lsJet1_msoftdrop",          ak8lsjet1.msoftdrop);
		self.out.fillBranch("ak8lsJet1_btagHbb",            ak8lsjet1.btagHbb);
		self.out.fillBranch("ak8lsJet1_tau1",               ak8lsjet1.tau1);
		self.out.fillBranch("ak8lsJet1_tau2",               ak8lsjet1.tau2);
		self.out.fillBranch("ak8lsJet1_subjet0_pt",         subjet1.pt);
		self.out.fillBranch("ak8lsJet1_subjet0_eta",        subjet1.eta);
		self.out.fillBranch("ak8lsJet1_subjet0_phi",        subjet1.phi);
 	        self.out.fillBranch("ak8lsJet1_subjet0_CSV",        subjet1.btagDeepB);   
		#self.out.fillBranch("ak8lsJet1_subjet0_CSV",        subjet1.btagDeepFlavB);   
		self.out.fillBranch("ak8lsJet1_subjet1_pt",         subjet2.pt);
		self.out.fillBranch("ak8lsJet1_subjet1_eta",        subjet2.eta);
		self.out.fillBranch("ak8lsJet1_subjet1_phi",        subjet2.phi);
	        self.out.fillBranch("ak8lsJet1_subjet1_CSV",        subjet2.btagDeepB);
		#self.out.fillBranch("ak8lsJet1_subjet1_CSV",        subjet2.btagDeepFlavB);
	    else:
		self.out.fillBranch("ak8lsJet1_pt",                 -10000.0);
		self.out.fillBranch("ak8lsJet1_E",                  -10000.0);
		self.out.fillBranch("ak8lsJet1_eta",                -10000.0);
		self.out.fillBranch("ak8lsJet1_phi",                -10000.0);
		self.out.fillBranch("ak8lsJet1_msoftdrop",          -10000.0);
		self.out.fillBranch("ak8lsJet1_btagHbb",            -10000.0);
		self.out.fillBranch("ak8lsJet1_tau1",               -10000.0);
		self.out.fillBranch("ak8lsJet1_tau2",               -10000.0);
		self.out.fillBranch("ak8lsJet1_subjet0_pt",         -10000.0);
		self.out.fillBranch("ak8lsJet1_subjet0_eta",        -10000.0);
		self.out.fillBranch("ak8lsJet1_subjet0_phi",        -10000.0);
		self.out.fillBranch("ak8lsJet1_subjet0_CSV",        -10000.0);   
		self.out.fillBranch("ak8lsJet1_subjet1_pt",         -10000.0);
		self.out.fillBranch("ak8lsJet1_subjet1_eta",        -10000.0);
		self.out.fillBranch("ak8lsJet1_subjet1_phi",        -10000.0);
		self.out.fillBranch("ak8lsJet1_subjet1_CSV",        -10000.0);   
	    
            if (len(ak8lsjets) >= 2):
		ak8lsjet2 = ak8lsjets[1]
		ak8lsjet2_p4 = ROOT.TLorentzVector()
		ak8lsjet2_p4.SetPtEtaPhiM(ak8lsjets[1].pt, ak8lsjets[1].eta, ak8lsjets[1].phi, ak8lsjets[1].mass)
		subjet1 = ak8lssubjets[ak8lsjet2.subJetIdx1]
		subjet2 = ak8lssubjets[ak8lsjet2.subJetIdx2]
		self.out.fillBranch("ak8lsJet2_pt",                 ak8lsjet2.pt);
		self.out.fillBranch("ak8lsJet2_E",                  ak8lsjet2_p4.E());
		self.out.fillBranch("ak8lsJet2_eta",                ak8lsjet2.eta);
		self.out.fillBranch("ak8lsJet2_phi",                ak8lsjet2.phi);
		self.out.fillBranch("ak8lsJet2_msoftdrop",          ak8lsjet2.msoftdrop);
		self.out.fillBranch("ak8lsJet2_btagHbb",            ak8lsjet2.btagHbb);
		self.out.fillBranch("ak8lsJet2_tau1",               ak8lsjet2.tau1);
		self.out.fillBranch("ak8lsJet2_tau2",               ak8lsjet2.tau2);
		self.out.fillBranch("ak8lsJet2_subjet0_pt",         subjet1.pt);
		self.out.fillBranch("ak8lsJet2_subjet0_eta",        subjet1.eta);
		self.out.fillBranch("ak8lsJet2_subjet0_phi",        subjet1.phi);
	        self.out.fillBranch("ak8lsJet2_subjet0_CSV",        subjet1.btagDeepB);   
		#self.out.fillBranch("ak8lsJet2_subjet0_CSV",        subjet1.btagDeepFlavB);   
		self.out.fillBranch("ak8lsJet2_subjet1_pt",         subjet2.pt);
		self.out.fillBranch("ak8lsJet2_subjet1_eta",        subjet2.eta);
		self.out.fillBranch("ak8lsJet2_subjet1_phi",        subjet2.phi);
	        self.out.fillBranch("ak8lsJet2_subjet1_CSV",        subjet2.btagDeepB);
		#self.out.fillBranch("ak8lsJet2_subjet1_CSV",        subjet2.btagDeepFlavB);
	    else:
		self.out.fillBranch("ak8lsJet2_pt",                 -10000.0);
		self.out.fillBranch("ak8lsJet2_E",                  -10000.0);
		self.out.fillBranch("ak8lsJet2_eta",                -10000.0);
		self.out.fillBranch("ak8lsJet2_phi",                -10000.0);
		self.out.fillBranch("ak8lsJet2_msoftdrop",          -10000.0);
		self.out.fillBranch("ak8lsJet2_tau1",               -10000.0);
		self.out.fillBranch("ak8lsJet2_tau2",               -10000.0);
		self.out.fillBranch("ak8lsJet2_btagHbb",            -10000.0);
		self.out.fillBranch("ak8lsJet2_subjet0_pt",         -10000.0);
		self.out.fillBranch("ak8lsJet2_subjet0_eta",        -10000.0);
		self.out.fillBranch("ak8lsJet2_subjet0_phi",        -10000.0);
		self.out.fillBranch("ak8lsJet2_subjet0_CSV",        -10000.0);   
		self.out.fillBranch("ak8lsJet2_subjet1_pt",         -10000.0);
		self.out.fillBranch("ak8lsJet2_subjet1_eta",        -10000.0);
		self.out.fillBranch("ak8lsJet2_subjet1_phi",        -10000.0);
		self.out.fillBranch("ak8lsJet2_subjet1_CSV",        -10000.0);   
	    
	        
	        

        met_p4 = ROOT.TLorentzVector()## mass=0, eta=0
        met_p4.SetPtEtaPhiM(metPt,0.,metPhi,0.) # only use met vector to derive transverse quantities
	self.out.fillBranch("n_presel_mu", len(lepton_muons))
	self.out.fillBranch("n_presel_ele", len(lepton_electrons))
	self.out.fillBranch("n_presel_ak4Jet", len(jets_clean))
	self.out.fillBranch("n_presel_ak8Jet", len(ak8jets_clean))
	self.out.fillBranch("n_presel_ak8lsJet", len(ak8lsjets_clean))
	self.out.fillBranch("PFMET",  met_p4.Pt())
	self.out.fillBranch("PFMETphi", met_p4.Phi())
	self.out.fillBranch("pu",  pu)
	self.out.fillBranch("PU_weight",  event_pu_weight)
        self.out.fillBranch("MC_weight",  genweight)
	self.out.fillBranch("run",  run)
	self.out.fillBranch("event",  ievent)
	self.out.fillBranch("ls",  luminosityBlock)

	if ievent == 191977:
	    for jet in jets_clean:
	        print "event id ",ievent, " jet pt ",jet.pt," eta ",jet.eta," phi ",jet.phi 

	fillBranches_muon(lepton_muons)
	fillBranches_electron(lepton_electrons)
	fillBranches_jet(jets_clean)
	fillBranches_ak8jet(ak8jets_clean, ak8subjets)
	fillBranches_ak8lsjet(ak8lsjets_clean, ak8lssubjets)
	

        return True
                

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

#hhbbWW = lambda : HHbbWWProducer(True, "") 
#hhbbWW_data = lambda x : HHbbWWProducer(False, x) 
