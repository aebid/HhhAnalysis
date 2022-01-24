import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
print "ROOT version ", ROOT.gROOT.GetVersion()
from math import sqrt, cos
import copy
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.output import FriendOutput

from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc...
import sys
#sys.path.append('...')

import __builtin__

if hasattr(__builtin__, "Runyear"):
    Runyear = __builtin__.Runyear
else:
    import RunConfiguration as RunConfig
    Runyear = RunConfig.Runyear

import POGRecipesRun2_allyears as POGRecipesRun2


if not (Runyear == 2016 or Runyear == 2017 or Runyear == 2018):
    sys.exit("Wrong run year value: {runyear}".format(runyear = Runyear))

import Run2DiLeptonTrigger as Run2Trigger

print("HHbbWWProducer, import finished here, Runyear ", Runyear)

errorcolor1 = '\x1b[1;31m'
errorcolor2 = '\x1b[0m'

class HHbbWWProducer(Module):
    ## data or MC, which L! trigger, HLT?
    ###kwargs: triggertype, verbose, run_lumi
    def __init__(self, isMC, **kwargs):
        print "init HHbbWWProducer"
        self.writeHistFile = True
        self.isMC = isMC
        #self.runyear = Runyear

        #Self variables
        self.ievent = -1
        self.luminosityBlock = -1
        self.run = -1
        self.muons_pre = []
        self.muons_fakeable = []
        self.muons_tight = []
        self.electrons_pre = []
        self.electrons_cleaned = []
        self.electrons_fakeable = []
        self.electrons_tight = []
        self.jets = []
        self.jets_pre = []
        self.jets_clean = []
        self.jets_btagged_medium = []
        self.jets_btagged_loose = []
        self.ak8jets_pre = []
        self.ak8jets_clean = []
        self.ak8jets_btagged = []
        self.ak8subjets = []
        self.met = -1
        self.PU_weight = -1
        self.MC_weight = -1
        self.HLT = -1
        self.flag = -1
        self.taus = []

    def beginJob(self, histFile=None, histDirName=None):
        print "BeginJob "

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        print "BeginFiles "
        self.out = wrappedOutputTree
        self.Single_Signal = copy.deepcopy(wrappedOutputTree)
        self.Single_Signal._tree.SetName("syncTree_hhbb1l_SR")
        self.Single_Signal._file = wrappedOutputTree._file
        self.Single_Signal._intree = wrappedOutputTree._intree
        self.Single_Signal._tree.SetDirectory(wrappedOutputTree._file)

        self.Single_Fake = copy.deepcopy(wrappedOutputTree)
        self.Single_Fake._tree.SetName("syncTree_hhbb1l_Fake")
        self.Single_Fake._file = wrappedOutputTree._file
        self.Single_Fake._intree = wrappedOutputTree._intree
        self.Single_Fake._tree.SetDirectory(wrappedOutputTree._file)

        self.Double_Signal = copy.deepcopy(wrappedOutputTree)
        self.Double_Signal._tree.SetName("syncTree_hhbb2l_SR")
        self.Double_Signal._file = wrappedOutputTree._file
        self.Double_Signal._intree = wrappedOutputTree._intree
        self.Double_Signal._tree.SetDirectory(wrappedOutputTree._file)

        self.Double_Fake = copy.deepcopy(wrappedOutputTree)
        self.Double_Fake._tree.SetName("syncTree_hhbb2l_Fake")
        self.Double_Fake._file = wrappedOutputTree._file
        self.Double_Fake._intree = wrappedOutputTree._intree
        self.Double_Fake._tree.SetDirectory(wrappedOutputTree._file)

        self.addbranches(self.out)
        self.addbranches(self.Single_Signal)
        self.addbranches(self.Single_Fake)
        self.addbranches(self.Double_Signal)
        self.addbranches(self.Double_Fake)

    def addbranches(self, out):
        out.branch("event", "I");		#Event Number
        out.branch("ls", "I");			#Lumi Section Number
        out.branch("run", "I");			#Run Number
        out.branch("n_presel_mu", "I");		#Number of muons passing the preselection
        out.branch("n_fakeablesel_mu", "I");	#Number of muons passing the fakeable selection
        out.branch("n_mvasel_mu", "I");		#Number of muons passing the tight selection
        out.branch("n_presel_ele", "I");	#Number of electrons passing the preselection
        out.branch("n_fakeablesel_ele", "I");	#Number of electrons passing the fakeable selection
        out.branch("n_mvasel_ele", "I");	#Number of electrons passing the tight selection
        out.branch("n_presel_ak4Jet", "I");	#Number of AK4 Jets passing the preselection and cleaning
        out.branch("n_presel_ak4JetVBF", "I");	#Number of AK4 Jets passing the VBF selection (only preselection if inclusive) and cleaning
        out.branch("n_presel_ak8Jet", "I");	#Number of AK8 Jets passing the preselection and cleaning (+btagging)
        out.branch("n_presel_ak8lsJet", "I");	#Number of AK8LS Jets passing the preselection
        out.branch("n_loose_ak4BJet", "I");	#Number of AK4 b-jets passing the preselection and loose b-tagging WP
        out.branch("n_medium_ak4BJet", "I");	#Number of AK4 b-jets passing the preselection and medium b-tagging WP
        out.branch("is_ee", "I");		#True if two selected leptons are electrons (only in the DL category)
        out.branch("is_mm", "I");		#True if two selected leptons are muons (only in the DL category)
        out.branch("is_em", "I");		#True if two selected leptons have opposite flavor (only in the DL category)
        out.branch("is_boosted", "I");		#True if the event falls into the boosted category (SL, DL)
        out.branch("is_semiboosted", "I");	#True if the event falls into the boosted category (SL)
        out.branch("is_resolved", "I");		#True if the event falls into the resolved category (SL, DL)

        ######## Muons ########			#2 leading preselected muons (in conept!!!)
        for i in [1, 2]:
          out.branch("mu{i}_pt".format(i = i), "F");
          out.branch("mu{i}_conept".format(i = i), "F");		#Corrected pT: pT of the muon if lepMVA > 0.90 and passes Medium ID, else 0.90*pT (associated jet)
          out.branch("mu{i}_eta".format(i = i), "F");
          out.branch("mu{i}_phi".format(i = i), "F");
          out.branch("mu{i}_E".format(i = i), "F");
          out.branch("mu{i}_charge".format(i = i), "I");
          out.branch("mu{i}_miniRelIso".format(i = i), "F");	#Relative isolation variable computed with the miniIso approach
          out.branch("mu{i}_PFRelIso04".format(i = i), "F");	#PF relative isolation dR = 0.4, total (with rho*EA PU corrections)
          out.branch("mu{i}_jetNDauChargedMVASel".format(i = i), "I");	#Number of charged constituents in the associated jet
          out.branch("mu{i}_jetPtRel".format(i = i), "F");	#Relative pT of the muon wrt the associated jet using the LepAware approach
          out.branch("mu{i}_jetRelIso".format(i = i), "F");	#Relative isolation of the muon wrt the associated jet using the LepAware approach
          out.branch("mu{i}_jetDeepJet".format(i = i), "F");	#DeepJet value of the associated jet
          out.branch("mu{i}_sip3D".format(i = i), "F");
          out.branch("mu{i}_dxy".format(i = i), "F");		#computed with innerTrack
          out.branch("mu{i}_dxyAbs".format(i = i), "F");
          out.branch("mu{i}_dz".format(i = i), "F");		#computed with innerTrack
          out.branch("mu{i}_segmentCompatibility".format(i = i), "F");
          out.branch("mu{i}_leptonMVA".format(i = i), "F");
          out.branch("mu{i}_mediumID".format(i = i), "I");	#Flag if muon passes HIP-safe medium PF muon ID
          out.branch("mu{i}_dpt_div_pt".format(i = i), "F");	#Relative error on the track pt (computed with best muon track)
          out.branch("mu{i}_isfakeablesel".format(i = i), "I");	#Flag if muon passes fakeable selections
          out.branch("mu{i}_ismvasel".format(i = i), "I");	#Flag if muon passes mva-based selections
          out.branch("mu{i}_isGenMatched".format(i = i), "F");

        ######## Electrons ########		#2 leading preselected electrons (in conept!!!)
        for i in [1, 2]:
          out.branch("ele{i}_pt".format(i = i), "F");
          out.branch("ele{i}_conept".format(i = i), "F");	#Corrected pT: pT of the electron if lepMVA > 0.90, else 0.90*pT (associated jet)
          out.branch("ele{i}_eta".format(i = i), "F");
          out.branch("ele{i}_phi".format(i = i), "F");
          out.branch("ele{i}_E".format(i = i), "F");
          out.branch("ele{i}_charge".format(i = i), "I");
          out.branch("ele{i}_miniRelIso".format(i = i), "F");	#Relative isolation variable computed with the miniIso approach
          out.branch("ele{i}_PFRelIso04".format(i = i), "F");	#PF relative isolation dR = 0.4, total (with rho*EA PU corrections)
          out.branch("ele{i}_jetNDauChargedMVASel".format(i = i), "I");	#Number of charged constituents in the associated jet
          out.branch("ele{i}_jetPtRel".format(i = i), "F");	#Relative pT of the electron wrt the associated jet using the LepAware approach
          out.branch("ele{i}_jetRelIso".format(i = i), "F");	#pT ratio of the electron wrt the associated jet using the LepAware approach
          out.branch("ele{i}_jetDeepJet".format(i = i), "F");	#DeepJet value of the associated jet
          out.branch("ele{i}_sip3D".format(i = i), "F");
          out.branch("ele{i}_dxy".format(i = i), "F");		#computed with innerTrack
          out.branch("ele{i}_dxyAbs".format(i = i), "F");
          out.branch("ele{i}_dz".format(i = i), "F");		#computed with innerTrack
          out.branch("ele{i}_ntMVAeleID".format(i = i), "F");	#non-triggering POG MVA electron ID discriminator
          out.branch("ele{i}_leptonMVA".format(i = i), "F");
          out.branch("ele{i}_passesConversionVeto".format(i = i), "I");
          out.branch("ele{i}_nMissingHits".format(i = i), "I");
          out.branch("ele{i}_sigmaEtaEta".format(i = i), "F");	#\sigma_{i\eta i\eta}
          out.branch("ele{i}_HoE".format(i = i), "F");		#H/E
          out.branch("ele{i}_OoEminusOoP".format(i = i), "F");	#1/E - 1/P
          out.branch("ele{i}_isfakeablesel".format(i = i), "I");	#flag if electorn passes fakeable selection
          out.branch("ele{i}_ismvasel".format(i = i), "I");	#flag if electron passes tight selection
          out.branch("ele{i}_isGenMatched".format(i = i), "I");

        ######## AK4 Jets ########             #4 leading preselected AK4 jets
        for i in [1, 2, 3, 4]:
          out.branch("ak4Jet{i}_pt".format(i = i), "F");	#pT of leading (highest pT) AK4 jet passing selection criteria for AK4 jets
          out.branch("ak4Jet{i}_eta".format(i = i), "F");
          out.branch("ak4Jet{i}_phi".format(i = i), "F");
          out.branch("ak4Jet{i}_E".format(i = i), "F");
          out.branch("ak4Jet{i}_CSV".format(i = i), "F");	#DeepFlavB score
          out.branch("ak4Jet{i}_btagSF".format(i = i), "F");	#b-tagging SF

        ######## AK4 VBF Jets ########             #2 leading VBF jets in inclusive analysis OR 2 selected VBF jets in the event category
        for i in [1, 2]:
          out.branch("ak4JetVBF{i}_pt".format(i = i), "F");        #pT of leading (highest pT) AK4 jet passing selection criteria for AK4 jets
          out.branch("ak4JetVBF{i}_eta".format(i = i), "F");
          out.branch("ak4JetVBF{i}_phi".format(i = i), "F");
          out.branch("ak4JetVBF{i}_E".format(i = i), "F");
          out.branch("ak4JetVBF{i}_CSV".format(i = i), "F");       #DeepFlavB score
          out.branch("ak4JetVBF{i}_btagSF".format(i = i), "F");    #b-tagging SF

        ######## AK8 Jets ########             #2 leading preselected + cleaned + btagged AK8 jets
        for i in [1, 2]:
          out.branch("ak8Jet{i}_pt".format(i = i), "F");	#pT of leading (highest pT) AK8 jet passing selection criteria for AK8 jets
          out.branch("ak8Jet{i}_eta".format(i = i), "F");
          out.branch("ak8Jet{i}_phi".format(i = i), "F");
          out.branch("ak8Jet{i}_E".format(i = i), "F");
          out.branch("ak8Jet{i}_msoftdrop".format(i = i), "F");	#softdrop mass
          out.branch("ak8Jet{i}_tau1".format(i = i), "F");	#N-subjettiness (1 axis)
          out.branch("ak8Jet{i}_tau2".format(i = i), "F");	#N-subjettiness (2 axis)
          for j in [0, 1]:
            out.branch("ak8Jet{i}_subjet{j}_pt".format(i = i, j = j), "F");
            out.branch("ak8Jet{i}_subjet{j}_eta".format(i = i, j = j), "F");
            out.branch("ak8Jet{i}_subjet{j}_phi".format(i = i, j = j), "F");
            out.branch("ak8Jet{i}_subjet{j}_CSV".format(i = i, j = j), "F");	#DeepB score

        ######## MET ########             #Missing transverse momentum
        out.branch("PFMET", "F");	#MET pT
        out.branch("PFMETphi", "F");       #MET phi

        ######## HME ########             #HH mass reconstructed by HME algorithm (arXiv:1701.04442)
        out.branch("HME", "F");		#HH mass computed for leading and subleading fakeable lepton and the two AK4 jets with the highest DeepB score

        ######## Event weights ########
        out.branch("PU_weight", "F");	#Pileup weight
        out.branch("PU_jetID_SF", "F");	#SF due to PU jet ID cut on the selected AK4 jets
        out.branch("MC_weight", "F");	#MC generator weight
        out.branch("topPt_wgt", "F");	#top pT reweighting SF
        out.branch("btag_SF", "F");	#b-tagging SF (w/o the corrective ratio)
        out.branch("trigger_SF", "F");	#trigger SF
        out.branch("lepton_IDSF", "F");	#lepton ID SF derived from selected tight leptons
        out.branch("lepton_IDSF_recoToLoose", "F");	#reco-to-loose lepton ID SF derived from selected tight leptons
        out.branch("lepton_IDSF_looseToTight", "F");	#loose-to-tight lepton ID SF derived from selected tight leptons
        out.branch("L1prefire", "F");	#L1 prefiring weight
        out.branch("fakeRate", "F");	#event-level jet -> lepton fake rate (applied only in fake CR)
        out.branch("vbf_m_jj", "F");	#mass of the VBF jet pair (filled only in the event categories)
        out.branch("vbf_dEta_jj", "F");	#difference in eta of the VBF jet pair (filled only in the event categories)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.Single_Signal.write()
        self.Single_Fake.write()
        self.Double_Signal.write()
        self.Double_Fake.write()


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        self.ievent = -1
        self.luminosityBlock = -1
        self.run = -1
        self.muons_pre = []
        self.muons_fakeable = []
        self.muons_tight = []
        self.electrons_pre = []
        self.electrons_cleaned = []
        self.electrons_fakeable = []
        self.electrons_tight = []
        self.jets = []
        self.jets_pre = []
        self.jets_clean = []
        self.jets_btagged_medium = []
        self.jets_btagged_loose = []
        self.ak8jets_pre = []
        self.ak8jets_clean = []
        self.ak8jets_btagged = []
        self.ak8subjets = []
        self.met = -1
        self.PU_weight = -1
        self.MC_weight = -1
        self.HLT = -1
        self.flag = -1
        self.taus = []




        self.ievent =  getattr(event,"event", False)
        self.luminosityBlock = getattr(event, "luminosityBlock", False)
        self.run = getattr(event, "run", False)


        self.met = Object(event, "MET")
        self.PU_weight = event.puWeight
        self.MC_weight = getattr(event, "genWeight", 1)

        electrons = list(Collection(event, "Electron"))
        muons = list(Collection(event, "Muon"))
        muon_pt_corrected = None
        if hasattr(event, "Muon_pt_corrected"):
          muon_pt_corrected = getattr(event, "Muon_pt_corrected")
          print("Corrected Muon pT")
          for muon in muons:
            imu = muons.index(muon)
            muon.pt = muon_pt_corrected[imu]
        self.jets = list(Collection(event, "Jet"))
        mht = Object(event, "MHT")
        ak8jets = list(Collection(event, "FatJet"))
        self.ak8subjets = list(Collection(event, "SubJet"))

        self.muons_pre = [x for x in muons if self.muonPreselection(x)]

        self.electrons_pre = [x for x in electrons if self.electronPreselection(x)]
        self.electrons_cleaned = [x for x in self.electrons_pre if self.leptonCleaning(x)]

        self.muons_pre.sort(key=lambda x:self.conept(x), reverse=True)
        self.electrons_cleaned.sort(key=lambda x:self.conept(x), reverse=True)

        self.muons_fakeable = [x for x in self.muons_pre if self.muonFakeable(x)]
        self.muons_tight = [x for x in self.muons_fakeable if self.muonTight(x)]

        self.electrons_fakeable = [x for x in self.electrons_cleaned if self.electronFakeable(x)]
        self.electrons_tight = [x for x in self.electrons_fakeable if self.electronTight(x)]

        self.jets_pre = [x for x in self.jets if self.ak4jetPreselection(x)]
        self.jets_clean = [x for x in self.jets_pre if self.jetCleaning(x, 0.4)]
        self.jets_btagged_medium = [x for x in self.jets_clean if self.ak4jetBtagging(x, "medium")]
        self.jets_btagged_loose = [x for x in self.jets_clean if self.ak4jetBtagging(x, "loose")]

        self.ak8jets_pre = [x for x in ak8jets if self.ak8jetPreselection(x)]
        self.ak8jets_clean = [x for x in self.ak8jets_pre if self.jetCleaning(x, 0.8)]
        self.ak8jets_btagged = [x for x in self.ak8jets_clean if self.ak8jetBtagging(x)]


        self.HLT = Object(event, "HLT")
        self.flag = Object(event, "Flag")
        self.taus = list(Collection(event, "Tau"))
        #category = "None"
        single_category = self.single_lepton()
        double_category = self.double_lepton()



        self.fillBranches(self.out)
        if (("Single" in single_category) and ("Signal" in single_category)):
          self.fillBranches(self.Single_Signal)
          self.Single_Signal.fill()
        if (("Single" in single_category) and ("Fake" in single_category)):
          self.fillBranches(self.Single_Fake)
          self.Single_Fake.fill()
        if (("Double" in double_category) and ("Signal" in double_category)):
          self.fillBranches(self.Double_Signal)
          self.Double_Signal.fill()
        if (("Double" in double_category) and ("Fake" in double_category)):
          self.fillBranches(self.Double_Fake)
          self.Double_Fake.fill()

        return True

    def conept(self, lep):
      """
      #https://github.com/CERN-PH-CMG/cmgtools-lite/blob/f8a34c64a4489d94ff9ac4c0d8b0b06dad46e521/TTHAnalysis/python/tools/conept.py#L74
      if (abs(lep.pdgId) != 11 and abs(lep.pdgId) != 13): return lep.pt
      if (abs(lep.pdgId) != 13 or lep.mediumMuonId > 0) and lep.mvaTTH > 0.90: return lep.pt
      return lep.pt #Currently do not have jetPtRatiov2, looking for a fix
      #else: return 0.90* lep.pt / lep.jetPtRatiov2
      """
      #Fix taken from Florian
      #https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/BaseHHtobbWW.py#L951-L962Le
      if (abs(lep.pdgId) != 11 and abs(lep.pdgId) != 13): return lep.pt
      elif (abs(lep.pdgId) == 11 and lep.mvaTTH > 0.30): return lep.pt
      elif (abs(lep.pdgId) == 13 and lep.mediumId and lep.mvaTTH > 0.50): return lep.pt
      else: return 0.9 * lep.pt * (1.0 + lep.jetRelIso)



    #Object Selection https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/-/blob/master/Legacy/objects.md
    def get_jet_from_lepton(self, lep):
      jetId = lep.jetIdx
      if jetId < 0 or jetId > len(self.jets):
          return 0
      return self.jets[jetId]

    def jetDeepJetUpper(self, muon):
      jet_pt = 0.9*muon.pt*(1+muon.jetRelIso)
      min_pt = 20; max_pt = 45
      wploose = [0.0613, 0.0521, 0.0494]
      wpmedium = [0.3093, 0.3033, 0.2770]
      x = min(max(0, jet_pt - min_pt)/(max_pt - min_pt), 1)
      return x*wploose[Runyear-2016] + (1-x)*wpmedium[Runyear-2016]

    def muonPreselection(self, muon):
      #pT >= 5, abs(eta) <= 2.4, abs(dxy) <= 500 um, abs(dz) <= 1mm, miniIso <= 0.4, sip3D <= 8, looseID
      return abs(muon.eta)<=2.4 and muon.pt>=5 and abs(muon.dxy)<=0.05 and abs(muon.dz)<=0.1 and muon.miniPFRelIso_all<=0.4 and muon.sip3d<=8 and muon.looseId

    def muonFakeable(self, muon):
      #Preselection, lepton cone-pT >= 10, jetDeepJet <= medium WP (per year), if leptonMVA <= 0.5: JetRelIso <= 0.5 and jetDeepJet upper cut obtained with interpolation
      jetDeepJet_MedWP = [0.3093, 0.3033, 0.2770]
      jet = self.get_jet_from_lepton(muon)
      #JetRelIso < 0.8 *** Changed https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/BaseHHtobbWW.py#L1033
      #return self.conept(muon) >= 10 and jet.btagDeepFlavB <= jetDeepJet_MedWP[Runyear - 2016] and (muon.mvaTTH > 0.5 or (muon.jetRelIso < 0.5 and jet.btagDeepFlavB <= self.jetDeepJetUpper(muon)))
      #If no associated jet, skip jet cuts *** 
      if jet == 0: return self.conept(muon) >= 10 and (muon.mvaTTH > 0.5 or muon.jetRelIso < 0.8)
      else: return self.conept(muon) >= 10 and jet.btagDeepFlavB <= jetDeepJet_MedWP[Runyear - 2016] and (muon.mvaTTH > 0.5 or (muon.jetRelIso < 0.8 and jet.btagDeepFlavB <= self.jetDeepJetUpper(muon)))

    def muonTight(self, muon):
      #Fakeable object selection, leptonMVA >= 0.5, Medium ID
      return (muon.mvaTTH >= 0.50 and muon.mediumId)

    def leptonCleaning(self, ele):
      #Preselected electrons are removed when they overlap (within cone size R = 0.3) with a preselected muon
      for mu in self.muons_pre:
        if deltaR(ele.eta, ele.phi, mu.eta, mu.phi) < 0.3:
          return False
      return True

    def electronPreselection(self, ele):
      #pt >= 7, abs(eta) <= 2.5, abs(dxy) <= 500um, abs(dz) <= 1mm, miniIso <= 0.4, sip3d <= 8, Loose working point, Fall17 v2 noIso MVA-based electron ID, number of missing inner hits <= 1
      return abs(ele.eta)<=2.5 and ele.pt>=7.0 and abs(ele.dxy)<=0.05 and abs(ele.dz)<=0.1 and ele.miniPFRelIso_all<=0.4 and ele.sip3d<=8 and ele.mvaFall17V2noIso_WPL and ele.lostHits<=1

    def electronFakeable(self, ele):
      #preselection, lepton cone-pT >= 10, (abs(eta) <= 1.479 and sieie <= 0.011), (1.479 < abs(eta) < 2.5 and sieie <= 0.030), hoe <= 0.10, -0.04 <= eInvMinusPInv, jetDeepJet <= medium WP, (leptonMVA <= 0.80 and JetRelIso <= 0.7 and Fall17V2noIso_WP80), number of missing inner hits = 0, convVeto
      eta = ele.eta + ele.deltaEtaSC
      #if not ((abs(eta) > 1.479 and abs(eta) < 2.5 and ele.sieie <= 0.030) or (abs(eta) <= 1.479 and ele.sieie <= 0.011)): *** you lose 1 event by including the eta < 2.5 cut, not sure if this is a big deal or not
      if not ((abs(eta) > 1.479 and ele.sieie <= 0.030) or (abs(eta) <= 1.479 and ele.sieie <= 0.011)):
        return False
      wpmedium = [0.3093, 0.3033, 0.2770]
      wptight = [0.7221, 0.7489, 0.7264]
      jet = self.get_jet_from_lepton(ele)
      if jet != 0:
      #jet.btagDeepFlavB <= wpmedium[Runyear-2016]) *** Changed https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/BaseHHtobbWW.py#L1069-L1071
        if ele.mvaTTH < 0.3:
          if not (jet.btagDeepFlavB <= wptight[Runyear-2016]): return False
        else:
          if not (jet.btagDeepFlavB <= wpmedium[Runyear-2016]): return False
      #ele.mvaTTH <= 0.30 *** Changed https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/BaseHHtobbWW.py#L1067
      #ele.mvaFall17V2noIso_WP80 *** Changed https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/BaseHHtobbWW.py#L1067
      if (ele.mvaTTH < 0.30 and not (ele.jetRelIso < 0.7 and ele.mvaFall17V2noIso_WP90)):
        return False
      return self.conept(ele) >= 10 and ele.hoe <= 0.10 and ele.eInvMinusPInv >= -0.04 and ele.lostHits == 0 and ele.convVeto

    def electronTight(self, ele):
      #Fakeable object selection, leptonMVA >= 0.80
      return ele.mvaTTH >= 0.30

    def jetCleaning(self, jet, dR): 
      #Not all fakeables, only leading in SL or DL!!!
      #AK4 jets are removed if they overlap with fakeable muons or electrons within dR < 0.4: AK8 dR < 0.8
      muons_fakeable = self.muons_fakeable; electrons_fakeable = self.electrons_fakeable
      leptons_fakeable = muons_fakeable + electrons_fakeable; leptons_fakeable.sort(key=lambda x:self.conept(x), reverse=True)

      for i in range(min(len(leptons_fakeable), 2)):
        #Single vs Double channels. Florian's code shows only need 1/2 for SL/DL channel (min(length, 2)) and we can check channels by len(leptons_fakeable) later
        lep = leptons_fakeable[i]
        if deltaR(jet.eta, jet.phi, lep.eta, lep.phi) < dR:
          return False
      return True

      #for mu in muons_fakeable:
      #  if deltaR(jet.eta, jet.phi, mu.eta, mu.phi) < 0.4:
      #    return False
      #for ele in electrons_fakeable:
      #  if deltaR(jet.eta, jet.phi, ele.eta, ele.phi) < 0.4:
      #    return False
      #return True

    def ak4jetPreselection(self, jet):
      #PF jet ID: 2016 - loose, 2017 - tight, 2018 - tight, pt >= 25, abs(eta) < 2.4, Jet PU ID (loose WP for pt < 50)
      badevent = False
      #if self.ievent in [8558, 9513, 11223, 11640, 11999, 185867, 186482, 195813, 64706, 67492, 74559, 150026, 150962, 79914, 96359, 145368, 146583, 165085, 189430, 190426]:
      #  print "Bad event found, ", self.ievent
      #  badevent = True
      if badevent: print "eta = ", jet.eta, " pt = ", jet.pt, " ID = ", jet.jetId
      if not (jet.pt > 50 or jet.puId & 4): #Bit operators!!!
        return False
      if badevent: print "Pass puID"
      PFJetID = [1, 2, 2] # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
      return (abs(jet.eta) <= 2.4 and jet.pt >= 25 and jet.jetId & PFJetID[Runyear-2016])

    def ak4jetBtagging(self, jet, wp):
      #The pfDeepFlavour (DeepJet) algorithm is used.
      #Gitlab gives tight/med/loose for each runyear, but does not explain which one to use
      wploose = [0.0614, 0.0521, 0.0494]
      wpmedium = [0.3093, 0.3033, 0.2770]
      wptight = [0.7221, 0.7489, 0.7264]
      if wp == "loose": return jet.btagDeepFlavB > wploose[Runyear-2016]
      if wp == "medium": return jet.btagDeepFlavB > wpmedium[Runyear-2016]
      if wp == "tight": return jet.btagDeepFlavB > wptight[Runyear-2016]

    def ak8jetPreselection(self, jet):
      #PF jet ID: 2016 - loose, 2017 - tight, 2018 - tight, pt >= 200, abs(eta) <= 2.4, two subjets each pt >= 20 and abs(eta) <= 2.4, 30 < msoftdrop < 210 GeV, tau2/tau1 <= 0.75
      ak8subjets = self.ak8subjets
      subjet1_idx = jet.subJetIdx1
      subjet2_idx = jet.subJetIdx2
      if subjet1_idx < 0 or subjet2_idx < 0 or subjet1_idx >= len(ak8subjets) or subjet2_idx >= len(ak8subjets):
        return False
      subjet1 = ak8subjets[subjet1_idx]
      subjet2 = ak8subjets[subjet2_idx]
      if not ((subjet1.pt >= 20 and abs(subjet1.eta) <= 2.4) and (subjet2.pt >= 20 and abs(subjet2.eta) <= 2.4)):
        return False
      PFJetID = [1, 2, 2] # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
      return (jet.jetId & PFJetID[Runyear-2016] and jet.pt >= 200 and abs(jet.eta) <= 2.4 and jet.msoftdrop >= 30 and jet.msoftdrop <= 210 and jet.tau2/jet.tau1 <= 0.75)

    def ak8jetBtagging(self, jet):
      #The DeepCSV b-tagging algorithm is used. The b-tagging algorithm is applied on the subjets.
      ak8subjets = self.ak8subjets
      subjet1_idx = jet.subJetIdx1
      subjet2_idx = jet.subJetIdx2
      if subjet1_idx < 0 or subjet2_idx < 0 or subjet1_idx >= len(self.ak8subjets) or subjet2_idx >= len(ak8subjets):
        return False
      subjet1 = ak8subjets[subjet1_idx]
      subjet2 = ak8subjets[subjet2_idx]
      wploose = [0.2217, 0.1522, 0.1241]
      wpmedium = [0.6321, 0.4941, 0.4184]
      wptight = [0.8953, 0.8001, 0.7527]
      return ((subjet1.btagDeepB > wpmedium[Runyear-2016] and subjet1.pt >= 30) or (subjet2.btagDeepB > wpmedium[Runyear-2016] and subjet2.pt >= 30))

    #Event Selection
    def single_lepton(self):
      #Pass MET filters
      #At least 1 fakeable lepton
      #If the leading cone-pT lepton is e (mu), pass single e (mu) trigger
      #cone-pt > 32 (25) for e (mu)
      #Invariant mass of each pair of preselected leptons (electrons NOT cleaned) must be greater than 12 GeV
      #Not more than 1 tight lepton - tight should be same as highest cone pT fakeable
      #Tau veto: no tau passing pt > 20, abs(eta) < 2.3, abs(dxy) <= 1000, abs(dz) <= 0.2, decay modes = {0, 1, 2, 10, 11}, and byMediumDeepTau2017v2VSjet, byVLooseDeepTau2017v2VSmu, byVVVLooseDeepTau2017v2VSe. Taus overlapping with fakeable electrons or fakeable muons within dR < 0.3 are not considered for the tau veto
      #No pair of same-flavor, opposite sign preselected leptons within 10 GeV of the Z mass
      #At least 1 medium btag (that can be on a AK8 jet): (#selJetsAK8_b >= 1 || #b-medium >= 1)
      #Minimal number of jets to construct an Hbb and admit an hadronic W with a missing jet: (#selJetsAK8_b == 0 && #selJetsAK4 >= 3) || (#selJetsAK8_b >= 1 && nJet_that_not_bb >= 1)
      muons_fakeable = self.muons_fakeable; electrons_fakeable = self.electrons_fakeable;
      muons_tight = self.muons_tight; electrons_tight = self.electrons_tight;
      ak8jets_btagged = self.ak8jets_btagged; jets_btagged = self.jets_btagged_medium;
      jets_pre = self.jets_pre;

      fake_leptons = muons_fakeable + electrons_fakeable
      fake_leptons.sort(key=lambda x:self.conept(x), reverse=True)
      tight_leptons = muons_tight + electrons_tight
      tight_leptons.sort(key=lambda x:self.conept(x), reverse=True)


      if not (len(fake_leptons) >= 1 and self.met_filters(self.flag)): return "None" 
      leading_lepton = fake_leptons[0]
      lep_type = ""
      if (leading_lepton in muons_fakeable): lep_type = "muon"
      if (leading_lepton in electrons_fakeable): lep_type = "ele"
      HLT = self.HLT
      if not ((lep_type == "muon" and self.single_muon_trigger(HLT)) or (lep_type == "ele" and self.single_electron_trigger(HLT))): return "None"
      if not ((lep_type == "muon" and self.conept(leading_lepton) > 25) or (lep_type == "ele" and self.conept(leading_lepton) > 32)): return "None"
      #if not (self.invar_mass_check()): return "None"
      if not ((len(tight_leptons) == 0) or (len(tight_leptons) == 1 and tight_leptons[0] == leading_lepton)): return "None"
      if not (self.tau_veto() and self.Zmass_and_invar_mass_cut()): return "None"
      if not (len(ak8jets_btagged) >= 1 or len(jets_btagged) >= 1): return "None"
      Jet_that_not_bb = []
      if len(ak8jets_btagged) > 0:
        #for i in range(len(ak8jets_btagged)): #Maybe 'leading ak8' isn't by pt, but by btag?
        #ak8jets_btagged.sort(key=lambda x:x.btagDeepB, reverse=True) #Nope, btag didn't change anything
        #  if (len(ak8jets_btagged) >= 2): print i, ak8jets_btagged[i].pt, ak8jets_btagged[i].btagDeepB
        Jet_that_not_bb = [x for x in jets_pre if deltaR(x.eta, x.phi, ak8jets_btagged[0].eta, ak8jets_btagged[0].phi) > 1.2]
      if not ((len(ak8jets_btagged) == 0 and len(jets_pre) >= 3) or (len(ak8jets_btagged) >= 1 and len(Jet_that_not_bb) >= 1)): return "None"
      category_string = "Single"
      if len(ak8jets_btagged) >= 1:
        category_string += "_HbbFat_WjjRes"
        if len(Jet_that_not_bb) >= 2: category_string += "_allReco"
        else: category_string += "_MissJet"
      else:
        category_string += "_Res"
        if len(jets_pre) >= 4:
          category_string += "_allReco"
          if len(jets_btagged) >= 1: category_string += "_2b"
          else: category_string += "_1b"
        else:
          category_string += "_MissWJet"
          if len(jets_btagged) >= 1: category_string += "_2b"
          else: category_string += "_1b"
      if len(tight_leptons) == 1: category_string += "_Signal"
      if len(tight_leptons) == 0: category_string += "_Fake"
      return category_string


    def double_lepton(self):
      #Pass MET filters
      #At least 2 fakeable leptons (choose the leading 2 in cone pT in the following)
      #if both fakeable leptons are electrons, the event needs to pass either the single electron or the double electron trigger; if both fakeable leptons are muons, the event needs to pass either the single muon or the double muon trigger; if one fakeable lepton is an electron and the other fakeable lepton is a muon, the event needs to pass either the single electron or the single muon or the muon+electron trigger
      #cone pT > 25 for the leading fakeable lepton and cone pT > 15 GeV for the subleading fakeable lepton
      #The 2 fakeable leptons must have opposite charge
      #No pair of same-flavor, opposite-sign preselected leptons within 10 GeV of the Z mass
      #The event needs to pass the selection in either the boosted or the resolved category:
      #In order for the event to pass the selection in the boosted category, it needs to contain at least one b-tagged Ak8 jet (see object definition)
      #In order for the event to pass the selection in the resolved category, it needs to contain at least 2 AK4 jets, of which at least one passes the medium working-point of the DeepJet b-tagging algorithm
      #The two categories are mutually exclusive. The higher priority is given to the boosted category, i.e. events satisfying the criteria for the boosted category are not considered for the resolved category.
      #The leading and subleading fakeable lepton both pass the tight lepton selection criteria
      #In MC, require MC matching of the leading and subleading fakelable lepton
      #Either the leading fakeable lepton, the subleading fakeable lepton or both fail the tight lepton selection criteria
      #In MC, require MC matching of the leading and subleading fakelable lepton

      badevent = False
      #if self.ievent in [8129, 8903, 184406, 186589, 186735, 187684, 187823, 192682, 193720, 194123, 194144, 194341, 194581, 194834, 195610, 195849, 64695, 64935, 66324, 66755, 67338, 67675, 67869, 72159, 72360, 72506, 73311, 73470, 73791, 74030, 75256, 75592, 148296, 148454, 149800, 149915, 150546, 150552, 150582, 76062, 76413, 76528, 76772, 77641, 78214, 78484, 78999, 79290, 79410, 84294, 84377, 84547, 84955, 85579, 85810, 85982, 86440, 87766, 97426, 97466, 97585, 97950, 98826, 99091, 144722, 89935, 89987, 91151, 91325, 91446, 91921, 146040, 146545, 147249, 147320, 147476, 147895, 151222, 151541, 151612, 151912, 151993, 165342, 165343, 166106, 166467, 167215, 16773, 167952, 188771, 188829, 189366, 189473, 189762, 189911, 189947, 190427, 190454, 190589, 191350, 191998]: #Him not in me
      #if self.ievent in [8009, 8096, 8271, 8358, 8423, 8474, 8485, 8489, 8500, 8641, 8805, 8904, 9166, 9474, 9544, 9783, 9792, 9817, 9907, 9893, 9896, 10009, 10089, 10114, 10202, 10255, 10362, 10441, 10628, 10730, 10713, 10767, 10785, 10962, 11265, 11321, 11474, 11497, 11538, 11612, 11695, 11723, 11912, 11990, 184018, 184049, 184078, 184110, 184336, 184482, 184493, 184605, 184883, 184888, 184928, 184941, 184910, 184966, 185114, 185356, 185532, 185581, 185630, 185671, 185736, 185822, 185819, 185832, 186113, 186225, 186288, 186318, 186376, 186448, 186524, 186586, 186705, 186731, 186894, 186981, 187133, 187261, 187350, 187398, 187534, 187680, 187739, 187918, 187959, 187965, 164098, 164143, 164296, 164330, 164355, 164358, 164429, 164454, 164463, 164545, 164695, 192061, 192117, 192172, 192214, 192419, 192428, 192568, 192639, 192687, 192685, 192711, 192822, 192832, 192976, 193313, 193369, 193677, 193694, 193736, 193739, 193852, 193878, 193879, 194149, 194261, 194239, 194247, 194275, 194300, 194371, 194440, 194457, 194465, 194470, 194844, 194858, 194923, 194965, 194986, 195243, 195263, 195450, 195470, 195483, 195969, 64019, 64030, 64116, 64287, 64491, 64727, 64889, 64897, 64940, 65020, 65057, 65096, 65226, 65276, 65306, 65341, 65349, 65378, 65406, 65436, 65573, 65695, 65699, 66078, 66077, 66276, 66349, 66546, 66598, 66652, 66720, 66766, 66799, 66817, 66854, 66973, 66976, 67086, 67539, 67683, 67733, 67787, 67813, 67933, 72023, 72101, 72143, 72145, 72200, 72221, 72267, 72317, 72402, 72614, 72767, 72838, 72911, 72981, 73018, 73128, 73182, 73271, 73569, 73605, 73619, 73730, 73809, 73972, 74064, 74117, 74103, 74161, 74293, 74390, 74394, 74454, 74577, 74717, 74724, 74758, 74809, 74821, 74840, 74913, 74980, 75098, 75121, 75147, 75375, 75388, 75438, 75465, 75511, 75571, 75623, 75823, 75813, 75835, 75875, 75989, 148143, 148329, 148355, 148446, 148527, 148535, 148628, 148669, 148696, 148716, 148918, 148993, 149292, 149331, 149433, 149554, 149597, 149600, 149859, 149908, 149910, 149931, 150199, 150216, 150328, 150323, 150543, 150650, 150673, 150785, 150922, 76010, 76099, 76126, 76153, 76195, 76328, 76441, 76474, 76541, 76610, 76721, 77031, 77190, 77280, 77409, 77483, 77682, 77812, 77817, 77836, 77870, 77922, 78030, 78201, 78257, 78414, 78538, 78754, 78792, 78846, 78903, 78966, 79083, 79170, 79176, 79209, 79231, 79262, 79295, 79318, 79361, 79414, 79472, 79474, 79560, 79652, 79667, 79684, 79682, 79714, 79721, 79748, 79865, 79987, 84099, 84128, 84152, 84151, 84194, 84414, 84413, 84531, 84542, 84597, 84605, 84675, 84782, 84820, 84841, 84865, 85139, 85146, 85141, 85173, 85251, 85313, 85393, 85462, 85472, 85589, 85697, 85699, 85710, 85968, 86091, 86090, 86186, 86232, 86316, 86408, 86434, 86585, 86636, 86974, 86982, 87062, 87313, 87345, 87382, 87410, 87437, 87497, 87517, 87602, 87741, 87753, 87752, 87844, 96022, 96093, 96088, 96168, 96257, 96455, 96472, 96537, 96572, 96721, 96722, 96738, 96821, 97184, 97296, 97298, 97531, 97551, 97617, 97724, 97743, 97852, 97843, 97918, 98005, 98076, 98217, 98220, 98267, 98406, 98464, 98544, 98596, 98649, 98676, 98740, 98756, 98817, 98861, 98856, 99014, 99048, 99046, 99123, 99131, 99161, 99308, 99343, 99403, 99550, 99690, 99708, 99960, 144039, 144310, 144526, 144595, 144608, 144604, 144664, 144675, 144920, 89051, 89074, 89173, 89179, 89301, 89516, 89647, 89776, 89847, 90021, 90179, 90354, 90655, 90966, 90995, 91211, 91376, 91440, 91448, 91550, 91668, 91881, 91883, 91898, 91933, 145128, 145306, 145340, 145360, 145431, 145527, 145662, 145741, 145755, 145786, 145903, 145920, 146070, 146103, 146332, 146358, 146379, 146389, 146399, 146562, 146585, 146588, 146615, 146698, 146763, 146772, 146804, 146818, 146941, 147044, 147061, 147181, 147417, 147410, 147459, 147472, 147544, 147553, 147719, 147764, 147840, 147851, 147890, 147946, 88227, 88234, 88317, 88336, 88411, 88446, 88593, 88701, 88736, 88772, 88815, 151021, 151060, 151143, 151150, 151173, 151229, 151386, 151553, 151718, 151944, 165041, 165059, 165097, 165069, 165188, 165244, 165345, 165360, 165505, 165539, 165608, 166007, 166257, 166337, 166502, 166516, 166784, 166797, 166918, 166940, 167149, 167146, 167155, 167309, 167328, 167400, 167478, 167504, 167548, 167556, 167569, 167925, 167957, 188047, 188125, 188189, 188179, 188303, 188368, 188404, 188666, 188672, 188688, 188700, 188880, 189099, 189076, 189090, 189198, 189201, 189347, 189684, 189760, 189848, 189896, 189940, 189950, 190051, 190074, 190062, 190144, 190271, 190385, 190361, 190478, 190658, 190694, 190679, 190714, 190795, 190827, 190877, 190904, 190897, 190963, 190965, 191202, 191360, 191460, 191547, 191557, 191800, 191799, 191960]:

      #if self.ievent in [8263, 11473, 185206, 192910, 64356, 67396, 150373, 84084, 86649, 87316, 89796, 90402, 91337, 91649, 145810, 146329, 167232, 167784, 189906]: #Him not me, fake
      #if self.ievent in [8044, 8064, 8645, 9019, 11351, 11530, 185022, 185234, 187378, 164481, 164569, 192058, 193346, 193860, 194005, 194259, 194568, 194668, 194662, 195726, 195895, 64664, 66435, 66668, 66793, 67244, 67455, 72218, 73521, 73549, 74489, 74496, 75080, 75352, 75603, 148392, 149679, 149906, 149979, 150635, 150934, 76002, 76847, 77449, 79557, 79976, 84719, 86914, 86926, 87319, 87376, 97109, 97414, 97593, 97839, 98493, 98816, 99223, 90477, 90827, 145074, 145366, 146223, 147002, 88024, 88080, 151007, 151655, 165643, 166672, 167133, 167523, 167724, 167706, 167935, 188238, 189180]: #Me not him, fake

      #  print "This is the bad event ", self.ievent
      #  badevent = True

      muons_fakeable = self.muons_fakeable; electrons_fakeable = self.electrons_fakeable;
      muons_tight = self.muons_tight; electrons_tight = self.electrons_tight;
      ak8jets_btagged = self.ak8jets_btagged; jets_btagged = self.jets_btagged_medium;
      jets_pre = self.jets_pre;

      fake_leptons = muons_fakeable + electrons_fakeable
      fake_leptons.sort(key=lambda x:self.conept(x), reverse=True)
      tight_leptons = muons_tight + electrons_tight
      tight_leptons.sort(key=lambda x:self.conept(x), reverse=True)
      if badevent:
        print "This is the bad event, starting"
      if not (len(fake_leptons) >= 2 and self.met_filters(self.flag)): return "None"
      leading_lepton = fake_leptons[0]
      subleading_lepton = fake_leptons[1]
      lep_type = ""
      if (leading_lepton in muons_fakeable):
        if (subleading_lepton in muons_fakeable): lep_type = "MuMu"
        if (subleading_lepton in electrons_fakeable): lep_type = "ElMu"
      if (leading_lepton in electrons_fakeable):
        if (subleading_lepton in muons_fakeable): lep_type = "ElMu"
        if (subleading_lepton in electrons_fakeable): lep_type = "ElEl"
      if badevent:
        print "This is the bad event, lep_type ", lep_type
      HLT = self.HLT
      if badevent:
        print "This is the bad event, got HLT"
      if not ((lep_type == "MuMu" and (self.double_muon_trigger(HLT) or self.single_muon_trigger(HLT))) or (lep_type == "ElEl" and (self.double_electron_trigger(HLT) or self.single_electron_trigger(HLT))) or ((lep_type == "ElMu" and (self.muon_electron_trigger(HLT) or self.single_muon_trigger(HLT) or self.single_electron_trigger(HLT))))): return "None"
      if badevent:
        print "This is the bad event, pass triggers"
      if not (self.conept(leading_lepton) > 25 and self.conept(subleading_lepton) > 15 and leading_lepton.charge != subleading_lepton.charge): return "None"
      #if not (self.invar_mass_check()): return "None"

      if badevent:
        print "This is the bad event, pass conept"

      if not ((len(tight_leptons) <= 2)): return "None"
      if not (self.Zmass_and_invar_mass_cut()): return "None"

      if badevent:
        print "This is the bad event, pass mass cut"

      category_string = "Double"
      if len(ak8jets_btagged) >= 1: category_string += "_Boosted"
      elif len(jets_pre) >= 2 and len(jets_btagged) >= 1: category_string += "_Resolved"
      if ((leading_lepton in tight_leptons) and (subleading_lepton in tight_leptons)): category_string += "_Signal"
      else: category_string += "_Fake"
      if (("Boosted" not in category_string) and ("Resolved" not in category_string)): return "None"

      if badevent:
        print "This is the bad event, string is ", category_string

      return category_string


    def tau_veto(self):
      #Tau veto: no tau passing pt>20, abs(eta) < 2.3, abs(dxy) <= 1000, abs(dz) <= 0.2, "decayModeFindingNewDMs", decay modes = {0, 1, 2, 10, 11}, and "byMediumDeepTau2017v2VSjet", "byVLooseDeepTau2017v2VSmu", "byVVVLooseDeepTau2017v2VSe". Taus overlapping with fakeable electrons or fakeable muons within dR < 0.3 are not considered for the tau veto
      #False -> Gets Removed : True -> Passes veto
      for i in self.taus:
        if (i.pt > 20 and abs(i.eta) < 2.3 and abs(i.dxy) <= 1000 and abs(i.dz) <= 0.2 and i.idDecayModeNewDMs and i.decayMode in [0,1,2,10,11] and i.idDeepTau2017v2p1VSjet >= 16 and i.idDeepTau2017v2p1VSmu >= 1 and i.idDeepTau2017v2p1VSe >= 1):
          for fake in (self.muons_fakeable + self.electrons_fakeable):
            if deltaR(fake.eta, fake.phi, i.eta, i.phi) > 0.3:
              return False
      return True

    def Zmass_and_invar_mass_cut(self):
      #No pair of same-flavor, opposite-sign preselected leptons within  10GeV of the Z mass
      #The invariant mass of each pair of preselected leptons (electrons not cleaned wrt muons) must be greater than 12 GeV
      muons_pre = self.muons_pre; electrons_pre = self.electrons_pre;
      pre_leptons = muons_pre + electrons_pre
      pre1 = ROOT.TLorentzVector(); pre2 = ROOT.TLorentzVector(); pre_pair = ROOT.TLorentzVector();
      Z_mass = 91.1876 #Z_mass = 80   #Florian's code doesn't have it at 80? Am I dumb? He also cuts out invar_mass < 12?
      for i in pre_leptons:
        pre1.SetPtEtaPhiM(i.pt, i.eta, i.phi, i.mass)
        for j in pre_leptons:
          if i == j: continue
          if not ((i in muons_pre and j in muons_pre) or (i in electrons_pre and j in electrons_pre)): continue
          if (i.charge == j.charge): continue
          pre2.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)
          pre_pair = pre1 + pre2
          if (abs(pre_pair.M() - Z_mass) < 10 or pre_pair.M() < 12): #Include 12 here to remove separate invar mass check
            return False
      return True

    def single_electron_trigger(self, HLT):
      if Runyear == 2016:
        return (HLT.Ele27_WPTight_Gsf or HLT.Ele25_eta2p1_WPTight_Gsf or HLT.Ele27_eta2p1_WPLoose_Gsf)
      elif Runyear == 2017:
        return (HLT.Ele35_WPTight_Gsf or HLT.Ele32_WPTight_Gsf)
      elif Runyear == 2018:
        return (HLT.Ele32_WPTight_Gsf)
      else:
        return False

    def single_muon_trigger(self, HLT):
      if Runyear == 2016:
        return (HLT.IsoMu22 or HLT.IsoTkMu22 or HLT.IsoMu22_eta2p1 or HLT.IsoTkMu22_eta2p1 or HLT.IsoMu24 or HLT.IsoTkMu24)
      elif Runyear == 2017:
        return (HLT.IsoMu24 or HLT.IsoMu27)
      elif Runyear == 2018:
        return (HLT.IsoMu24 or HLT.IsoMu27)
      else:
        return False

    def double_electron_trigger(self, HLT):
      if Runyear == 2016:
        return (HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)
      elif Runyear == 2017:
        return (HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL)
      elif Runyear == 2018:
        return (HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL)
      else:
        return False

    def double_muon_trigger(self, HLT):
      if Runyear == 2016:
        return (HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL or HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or HLT.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL or HLT.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ)
      elif Runyear == 2017:
        return (HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)
      elif Runyear == 2018:
        return (HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)
      else:
        return False

    def muon_electron_trigger(self, HLT):
      if Runyear == 2016:
        return (HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL or HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or HLT.Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL or HLT.Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ)
      elif Runyear == 2017:
        return (HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)
      elif Runyear == 2018:
        return (HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL)
      else:
        return False

    def met_filters(self, flag):
      #if not (((Runyear == 2017 or Runyear == 2018) and flag.ecalBadCalibReducedMINIAODFilter) or (Runyear == 2016)): #ecalBadCalibRediucedMINIAODFilter not in nanoAOD, using Flag_ecalBadCalibFilterV2 for now? Even florian doesn't have this figured out https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/METScripts.py#L25
      #  return False
      if not ((flag.eeBadScFilter and not self.isMC) or (self.isMC)):
        return False
      return (flag.goodVertices and flag.globalSuperTightHalo2016Filter and flag.HBHENoiseFilter and flag.HBHENoiseIsoFilter and flag.EcalDeadCellTriggerPrimitiveFilter and flag.BadPFMuonFilter)


    def fillBranches(self, out):
        out.fillBranch("event", self.ievent);
        out.fillBranch("ls", self.luminosityBlock);
        out.fillBranch("run", self.run);
        out.fillBranch("n_presel_mu", len(self.muons_pre));
        out.fillBranch("n_fakeablesel_mu", len(self.muons_fakeable));
        out.fillBranch("n_mvasel_mu", len(self.muons_tight));
        out.fillBranch("n_presel_ele", len(self.electrons_cleaned));
        out.fillBranch("n_fakeablesel_ele", len(self.electrons_fakeable));
        out.fillBranch("n_mvasel_ele", len(self.electrons_tight));
        out.fillBranch("n_presel_ak4Jet", len(self.jets_clean));
        out.fillBranch("n_presel_ak4JetVBF", -999);
        out.fillBranch("n_presel_ak8Jet", len(self.ak8jets_btagged));
        #out.fillBranch("n_presel_ak8lsJet", );
        out.fillBranch("n_loose_ak4BJet", len(self.jets_btagged_loose));
        out.fillBranch("n_medium_ak4BJet", len(self.jets_btagged_medium));
        out.fillBranch("is_ee", -999);
        out.fillBranch("is_mm", -999);
        out.fillBranch("is_em", -999);
        out.fillBranch("is_boosted", -999);
        out.fillBranch("is_semiboosted", -999);
        out.fillBranch("is_resolved", -999);

        ######## Muons ########
        for i in [1, 2]:
          if (i <= len(self.muons_pre)):
            mu = self.muons_pre[i-1]
            pt = mu.pt; conept = self.conept(mu); eta = mu.eta; phi = mu.phi; mass = mu.mass;
            mu_p4 = ROOT.TLorentzVector(); mu_p4.SetPtEtaPhiM(pt, eta, phi, mass)
            E = mu_p4.E(); charge = mu.charge;
            miniRelIso = mu.miniPFRelIso_all; PFRelIso04 = -999; jetNDauChargedMVASel = -999; jetPtRel = -999;
            jetRelIso = -999; jetDeepJet = -999; sip3D = mu.sip3d; dxy = mu.dxy; dxyAbs = abs(mu.dxy); dz = mu.dz;
            segmentCompatibility = mu.segmentComp; leptonMVA = mu.mvaTTH; mediumID = -999; dpt_div_pt = -999;
            isfakeablesel = (mu in self.muons_fakeable); ismvasel = (mu in self.muons_tight); isGenMatched = -999;
          else:
            value = -999
            pt = value; conept = value; eta = value; phi = value; mass = value;
            E = value; charge = value;
            miniRelIso = value; PFRelIso04 = value; jetNDauChargedMVASel = value; jetPtRel = value;
            jetRelIso = value; jetDeepJet = value; sip3D = value; dxy = value; dxyAbs = value; dz = value;
            segmentCompatibility = value; leptonMVA = value; mediumID = value; dpt_div_pt = value;
            isfakeablesel = value; ismvasel = value; isGenMatched = value;


          out.fillBranch("mu{i}_pt".format(i = i), pt);
          out.fillBranch("mu{i}_conept".format(i = i), conept);
          out.fillBranch("mu{i}_eta".format(i = i), eta);
          out.fillBranch("mu{i}_phi".format(i = i), phi);
          out.fillBranch("mu{i}_E".format(i = i), E);
          out.fillBranch("mu{i}_charge".format(i = i), charge);
          out.fillBranch("mu{i}_miniRelIso".format(i = i), miniRelIso);
          out.fillBranch("mu{i}_PFRelIso04".format(i = i), PFRelIso04);
          out.fillBranch("mu{i}_jetNDauChargedMVASel".format(i = i), jetNDauChargedMVASel);
          out.fillBranch("mu{i}_jetPtRel".format(i = i), jetPtRel);
          out.fillBranch("mu{i}_jetRelIso".format(i = i), jetRelIso);
          out.fillBranch("mu{i}_jetDeepJet".format(i = i), jetDeepJet);
          out.fillBranch("mu{i}_sip3D".format(i = i), sip3D);
          out.fillBranch("mu{i}_dxy".format(i = i), dxy);
          out.fillBranch("mu{i}_dxyAbs".format(i = i), dxyAbs);
          out.fillBranch("mu{i}_dz".format(i = i), dz);
          out.fillBranch("mu{i}_segmentCompatibility".format(i = i), segmentCompatibility);
          out.fillBranch("mu{i}_leptonMVA".format(i = i), leptonMVA);
          out.fillBranch("mu{i}_mediumID".format(i = i), mediumID);
          out.fillBranch("mu{i}_dpt_div_pt".format(i = i), dpt_div_pt);
          out.fillBranch("mu{i}_isfakeablesel".format(i = i), isfakeablesel);
          out.fillBranch("mu{i}_ismvasel".format(i = i), ismvasel);
          out.fillBranch("mu{i}_isGenMatched".format(i = i), isGenMatched);

        ######## Electrons ########
        for i in [1, 2]:
          if (i <= len(self.electrons_cleaned)):
            ele = self.electrons_cleaned[i-1]
	    pt = ele.pt; conept = self.conept(ele); eta = ele.eta; phi = ele.phi; mass = ele.mass;
            ele_p4 = ROOT.TLorentzVector(); ele_p4.SetPtEtaPhiM(pt, eta, phi, mass)
            E = ele_p4.E(); charge = ele.charge;
            miniRelIso = ele.miniPFRelIso_all; PFRelIso04 = -999; jetNDauChargedMVASel = -999; jetPtRel = -999;
            jetRelIso = -999; jetDeepJet = -999; sip3D = ele.sip3d; dxy = ele.dxy; dxyAbs = abs(ele.dxy); dz = ele.dz;
            ntMVAeleID = -999; leptonMVA = ele.mvaTTH; passesConversionVeto = ele.convVeto; nMissingHits = ele.lostHits;
            sigmaEtaEta = ele.sieie; HoE = ele.hoe; OoEminusOoP = ele.eInvMinusPInv;
            isfakeablesel = (ele in self.electrons_fakeable); ismvasel = (ele in self.electrons_tight); isGenMatched = -999;
          else:
            value = -9999
            pt = value; conept = value; eta = value; phi = value; mass = value;
            E = value; charge = value;
            miniRelIso = value; PFRelIso04 = value; jetNDauChargedMVASel = value; jetPtRel = value;
            jetRelIso = value; jetDeepJet = value; sip3D = value; dxy = value; dxyAbs = value; dz = value;
            ntMVAeleID = value; leptonMVA = value; passesConversionVeto = value; nMissingHits = value;
            sigmaEtaEta = value; HoE = value; OoEminusOoP = value;
            isfakeablesel = value; ismvasel = value; isGenMatched = value;

          out.fillBranch("ele{i}_pt".format(i = i), pt);
          out.fillBranch("ele{i}_conept".format(i = i), conept);
          out.fillBranch("ele{i}_eta".format(i = i), eta);
          out.fillBranch("ele{i}_phi".format(i = i), phi);
          out.fillBranch("ele{i}_E".format(i = i), E);
          out.fillBranch("ele{i}_charge".format(i = i), charge);
          out.fillBranch("ele{i}_miniRelIso".format(i = i), miniRelIso);
          out.fillBranch("ele{i}_PFRelIso04".format(i = i), PFRelIso04);
          out.fillBranch("ele{i}_jetNDauChargedMVASel".format(i = i), jetNDauChargedMVASel);
          out.fillBranch("ele{i}_jetPtRel".format(i = i), jetPtRel);
          out.fillBranch("ele{i}_jetRelIso".format(i = i), jetRelIso);
          out.fillBranch("ele{i}_jetDeepJet".format(i = i), jetDeepJet);
          out.fillBranch("ele{i}_sip3D".format(i = i), sip3D);
          out.fillBranch("ele{i}_dxy".format(i = i), dxy);
          out.fillBranch("ele{i}_dxyAbs".format(i = i), dxyAbs);
          out.fillBranch("ele{i}_dz".format(i = i), dz);
          out.fillBranch("ele{i}_ntMVAeleID".format(i = i), ntMVAeleID);
          out.fillBranch("ele{i}_leptonMVA".format(i = i), leptonMVA);
          out.fillBranch("ele{i}_passesConversionVeto".format(i = i), passesConversionVeto);
          out.fillBranch("ele{i}_nMissingHits".format(i = i), nMissingHits);
          out.fillBranch("ele{i}_sigmaEtaEta".format(i = i), sigmaEtaEta);
          out.fillBranch("ele{i}_HoE".format(i = i), HoE);
          out.fillBranch("ele{i}_OoEminusOoP".format(i = i), OoEminusOoP);
          out.fillBranch("ele{i}_isfakeablesel".format(i = i), isfakeablesel);
          out.fillBranch("ele{i}_ismvasel".format(i = i), ismvasel);
          out.fillBranch("ele{i}_isGenMatched".format(i = i), isGenMatched);

        ######## AK4 Jets ########
        for i in [1, 2, 3, 4]:
          if (i <= len(self.jets_clean)):
            ak4jet = self.jets_clean[i-1]
            pt = ak4jet.pt; eta = ak4jet.eta; phi = ak4jet.phi; mass = ak4jet.mass;
            ak4jet_p4 = ROOT.TLorentzVector(); ak4jet_p4.SetPtEtaPhiM(pt, eta, phi, mass);
            E = ak4jet_p4.E(); CSV = ak4jet.btagDeepFlavB; btagSF = -999;
          else:
            value = -99999
            pt = value; eta = value; phi = value; mass = value;
            E = value; CSV = value; btagSF = value;

          out.fillBranch("ak4Jet{i}_pt".format(i = i), pt);
          out.fillBranch("ak4Jet{i}_eta".format(i = i), eta);
          out.fillBranch("ak4Jet{i}_phi".format(i = i), phi);
          out.fillBranch("ak4Jet{i}_E".format(i = i), E);
          out.fillBranch("ak4Jet{i}_CSV".format(i = i), CSV);
          out.fillBranch("ak4Jet{i}_btagSF".format(i = i), btagSF);

        ######## AK4 VBF Jets ########
        for i in [1, 2]:
          out.fillBranch("ak4JetVBF{i}_pt".format(i = i), -999);
          out.fillBranch("ak4JetVBF{i}_eta".format(i = i), -999);
          out.fillBranch("ak4JetVBF{i}_phi".format(i = i), -999);
          out.fillBranch("ak4JetVBF{i}_E".format(i = i), -999);
          out.fillBranch("ak4JetVBF{i}_CSV".format(i = i), -999);
          out.fillBranch("ak4JetVBF{i}_btagSF".format(i = i), -999);

        ######## AK8 Jets ########
        for i in [1, 2]:
          if (i <= len(self.ak8jets_clean)):
            ak8jet = self.ak8jets_clean[i-1]
            pt = ak8jet.pt; eta = ak8jet.eta; phi = ak8jet.phi; mass = ak8jet.mass;
            ak8jet_p4 = ROOT.TLorentzVector(); ak8jet_p4.SetPtEtaPhiM(pt, eta, phi, mass);
            E = ak8jet_p4.E(); msoftdrop = ak8jet.msoftdrop; tau1 = ak8jet.tau1; tau2 = ak8jet.tau2;
            ak8subjet1 = self.ak8subjets[ak8jet.subJetIdx1]; ak8subjet2 = self.ak8subjets[ak8jet.subJetIdx2];
            sub0_pt = ak8subjet1.pt; sub0_eta = ak8subjet1.eta; sub0_phi = ak8subjet1.phi; sub0_btagDeepB = ak8subjet1.btagDeepB;
            sub1_pt = ak8subjet2.pt; sub1_eta = ak8subjet2.eta; sub1_phi = ak8subjet2.phi; sub1_btagDeepB = ak8subjet2.btagDeepB;

          else:
            value = -9999
            pt = value; eta = value; phi = value; mass = value;
            E = value; msoftdrop = value; tau1 = value; tau2 = value;
            sub0_pt = value; sub0_eta = value; sub0_phi = value; sub0_btagDeepB = value;
            sub1_pt = value; sub1_eta = value; sub1_phi = value; sub1_btagDeepB = value;

          out.fillBranch("ak8Jet{i}_pt".format(i = i), pt);
          out.fillBranch("ak8Jet{i}_eta".format(i = i), eta);
          out.fillBranch("ak8Jet{i}_phi".format(i = i), phi);
          out.fillBranch("ak8Jet{i}_E".format(i = i), E);
          out.fillBranch("ak8Jet{i}_msoftdrop".format(i = i), msoftdrop);
          out.fillBranch("ak8Jet{i}_tau1".format(i = i), tau1);
          out.fillBranch("ak8Jet{i}_tau2".format(i = i), tau2);

          out.fillBranch("ak8Jet{i}_subjet0_pt".format(i = i), sub0_pt);
          out.fillBranch("ak8Jet{i}_subjet0_eta".format(i = i), sub0_eta);
          out.fillBranch("ak8Jet{i}_subjet0_phi".format(i = i), sub0_phi);
          out.fillBranch("ak8Jet{i}_subjet0_CSV".format(i = i), sub0_btagDeepB);
          out.fillBranch("ak8Jet{i}_subjet1_pt".format(i = i), sub1_pt);
          out.fillBranch("ak8Jet{i}_subjet1_eta".format(i = i), sub1_eta);
          out.fillBranch("ak8Jet{i}_subjet1_phi".format(i = i), sub1_phi);
          out.fillBranch("ak8Jet{i}_subjet1_CSV".format(i = i), sub1_btagDeepB);


        
        ######## MET ########
        met_p4 = ROOT.TLorentzVector()
        met_p4.SetPtEtaPhiM(self.met.pt, 0., self.met.phi, 0.)
        out.fillBranch("PFMET", met_p4.Pt());
        out.fillBranch("PFMETphi", met_p4.Phi());

        ######## HME ########
        out.fillBranch("HME", -999);

        ######## Event Weights ########
        out.fillBranch("PU_weight", self.PU_weight);
        out.fillBranch("PU_jetID_SF", -999);
        out.fillBranch("MC_weight", self.MC_weight);
        out.fillBranch("topPt_wgt", -999);
        out.fillBranch("btag_SF", -999);
        out.fillBranch("trigger_SF", -999);
        out.fillBranch("lepton_IDSF", -999);
        out.fillBranch("lepton_IDSF_recoToLoose", -999);
        out.fillBranch("lepton_IDSF_looseToTight", -999);
        out.fillBranch("L1prefire", -999);
        out.fillBranch("fakeRate", -999);
        out.fillBranch("vbf_m_jj", -999);
        out.fillBranch("vbf_dEta_jj", -999);


