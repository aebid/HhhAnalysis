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
        self.jets_btagged = []
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
        self.jets_btagged = []
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
        self.jets_clean = [x for x in self.jets_pre if self.ak4jetCleaning(x)]
        self.jets_btagged = [x for x in self.jets_clean if self.ak4jetBtagging(x)]

        self.ak8jets_pre = [x for x in ak8jets if self.ak8jetPreselection(x)]
        self.ak8jets_clean = [x for x in self.ak8jets_pre if self.ak8jetCleaning(x)]
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
      if self.get_jet_from_lepton(muon) == 0:
        return False
      jetDeepJet_MedWP = [0.3093, 0.3033, 0.2770]
      jet = self.get_jet_from_lepton(muon)
      #JetRelIso < 0.8 *** Changed https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/BaseHHtobbWW.py#L1033
      return self.conept(muon) >= 10 and jet.btagDeepFlavB <= jetDeepJet_MedWP[Runyear - 2016] and (muon.mvaTTH > 0.5 or (muon.jetRelIso < 0.8 and jet.btagDeepFlavB <= self.jetDeepJetUpper(muon)))

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
      if not ((abs(ele.eta) > 1.479 and abs(ele.eta) < 2.5 and ele.sieie <= 0.030) or (abs(ele.eta) <= 1.479 and ele.sieie <= 0.011)):
        return False
      if self.get_jet_from_lepton(ele) == 0:
        return False
      wpmedium = [0.3093, 0.3033, 0.2770]
      if not (self.get_jet_from_lepton(ele).btagDeepFlavB <= wpmedium[Runyear-2016]):
        return False
      #ele.mvaTTH <= 0.30 *** Changed https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/BaseHHtobbWW.py#L1067
      #ele.mvaFall17V2noIso_WP80 *** Changed https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/BaseHHtobbWW.py#L1067
      if (ele.mvaTTH < 0.30 and not (ele.jetRelIso < 0.7 and ele.mvaFall17V2noIso_WP90)):
        return False
      return self.conept(ele) >= 10 and ele.hoe <= 0.10 and ele.eInvMinusPInv >= -0.04 and ele.lostHits == 0 and ele.convVeto

    def electronTight(self, ele):
      #Fakeable object selection, leptonMVA >= 0.80
      return ele.mvaTTH >= 0.30

    def ak4jetPreselection(self, jet):
      #PF jet ID: 2016 - loose, 2017 - tight, 2018 - tight, pt >= 25, abs(eta) < 2.4, Jet PU ID (loose WP for pt < 50)
      if (jet.pt < 50 and not jet.puId >= 4):
        return False
      PFJetID = [0, 2, 2]
      return (abs(jet.eta) <= 2.4 and jet.pt >= 25 and jet.jetId >= PFJetID[Runyear-2016])

    def ak4jetCleaning(self, jet):
      #AK4 jets are removed if they overlap with fakeable muons or electrons within dR< 0.4
      muons_fakeable = self.muons_fakeable; electrons_fakeable = self.electrons_fakeable
      for mu in muons_fakeable:
        if deltaR(jet.eta, jet.phi, mu.eta, mu.phi) < 0.4:
          return False
      for ele in electrons_fakeable:
        if deltaR(jet.eta, jet.phi, ele.eta, ele.phi) < 0.4:
          return False
      return True

    def ak4jetBtagging(self, jet):
      #The pfDeepFlavour (DeepJet) algorithm is used.
      #Gitlab gives tight/med/loose for each runyear, but does not explain which one to use
      wploose = [0.0614, 0.0521, 0.0494]
      wpmedium = [0.3093, 0.3033, 0.2770]
      wptight = [0.7221, 0.7489, 0.7264]
      return jet.btagDeepFlavB >= wpmedium[Runyear-2016]

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
      PFJetID = [0, 2, 2]
      return (jet.jetId >= PFJetID[Runyear-2016] and jet.pt >= 200 and abs(jet.eta) <= 2.4 and jet.msoftdrop > 30 and jet.msoftdrop < 210 and jet.tau2/jet.tau1 <= 0.75)

    def ak8jetCleaning(self, jet):
      #AK8 jets are removed if they overlap with fakeable muons or electrons within dR < 0.8
      muons_fakeable = self.muons_fakeable; electrons_fakeable = self.electrons_fakeable;
      for mu in muons_fakeable:
        if deltaR(jet.eta, jet.phi, mu.eta, mu.phi) < 0.8:
          return False
      for ele in electrons_fakeable:
        if deltaR(jet.eta, jet.phi, ele.eta, ele.phi) < 0.8:
          return False
      return True

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
      return ((subjet1.btagDeepB > wpmedium[Runyear-2016] and subjet1.pt > 30) or (subjet2.btagDeepB > wpmedium[Runyear-2016] and subjet2.pt > 30))

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
      ak8jets_btagged = self.ak8jets_btagged; jets_btagged = self.jets_btagged;
      jets_pre = self.jets_pre;

      fake_leptons = muons_fakeable + electrons_fakeable
      fake_leptons.sort(key=lambda x:self.conept(x), reverse=True)
      tight_leptons = muons_tight + electrons_tight
      tight_leptons.sort(key=lambda x:self.conept(x), reverse=True)
      if not (len(fake_leptons) > 0 and self.met_filters(self.flag)): return "None"
      leading_lepton = fake_leptons[0]
      if not ((leading_lepton in muons_fakeable and self.single_muon_trigger(self.HLT)) or (leading_lepton in electrons_fakeable and self.single_electron_trigger(self.HLT))): return "None"
      if not ((leading_lepton in muons_fakeable and self.conept(leading_lepton) > 25) or (leading_lepton in electrons_fakeable and self.conept(leading_lepton) > 32)): return "None"
      if not (self.invar_mass_check()): return "None"
      if not ((len(tight_leptons) == 0) or (len(tight_leptons) == 1 and tight_leptons[0] == leading_lepton)): return "None"
      if not (self.tau_veto() and self.Zmass_cut()): return "None"
      if not (len(ak8jets_btagged) >= 1 or len(jets_btagged) >= 1): return "None"
      Jet_that_not_bb = []
      if len(ak8jets_btagged) > 0:
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

      muons_fakeable = self.muons_fakeable; electrons_fakeable = self.electrons_fakeable;
      muons_tight = self.muons_tight; electrons_tight = self.electrons_tight;
      ak8jets_btagged = self.ak8jets_btagged; jets_btagged = self.jets_btagged;
      jets_pre = self.jets_pre;

      fake_leptons = muons_fakeable + electrons_fakeable
      fake_leptons.sort(key=lambda x:self.conept(x), reverse=True)
      tight_leptons = muons_tight + electrons_tight
      tight_leptons.sort(key=lambda x:self.conept(x), reverse=True)
      if not (len(fake_leptons) >= 2 and self.met_filters(self.flag)): return "None"
      leading_lepton = fake_leptons[0]
      subleading_lepton = fake_leptons[1]
      HLT = self.HLT
      if not ((leading_lepton in muons_fakeable and subleading_lepton in muons_fakeable and (self.double_muon_trigger(HLT) or self.single_muon_trigger(HLT))) or (leading_lepton in electrons_fakeable and subleading_lepton in electrons_fakeable and (self.double_electron_trigger(HLT) or self.single_electron_trigger(HLT))) or ((leading_lepton in muons_fakeable and subleading_lepton in electrons_fakeable) or ((leading_lepton in electrons_fakeable and subleading_lepton in muons_fakeable)) and self.muon_electron_trigger(HLT))): return "None"
      if not (self.conept(leading_lepton) > 25 and self.conept(subleading_lepton) > 15 and leading_lepton.charge != subleading_lepton.charge): return "None"
      if not (self.invar_mass_check()): return "None"

      if not ((len(tight_leptons) <= 2)): return "None"
      if not (self.Zmass_cut()): return "None"
      category_string = "Double"
      if len(ak8jets_btagged) >= 1: category_string += "_Boosted"
      elif len(jets_pre) >= 2 and len(jets_btagged) >= 1: category_string += "_Resolved"
      if ((leading_lepton in tight_leptons) and (subleading_lepton in tight_leptons)): category_string += "_Signal"
      else: category_string += "_Fake"
      if (("Boosted" not in category_string) and ("Resolved" not in category_string)): return "None"

      return category_string


    def invar_mass_check(self):
      #The invariant mass of each pair of preselected leptons (electrons not cleaned wrt muons) must be greater than 12 GeV
      lep1_p4 = ROOT.TLorentzVector(); lep2_p4 = ROOT.TLorentzVector(); pair_p4 = ROOT.TLorentzVector();
      pre_leptons = self.electrons_pre + self.muons_pre
      for i in pre_leptons:
        lep1_p4.SetPtEtaPhiM(i.pt, i.eta, i.phi, i.mass)
        for j in pre_leptons:
          if i == j: continue
          lep2_p4.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)
          pair_p4 = lep1_p4 + lep2_p4
          pair_mass = pair_p4.M()
          if (pair_mass < 12):
            return False
      return True

    def tau_veto(self):
      #Tau veto: no tau passing pt>20, abs(eta) < 2.3, abs(dxy) <= 1000, abs(dz) <= 0.2, "decayModeFindingNewDMs", decay modes = {0, 1, 2, 10, 11}, and "byMediumDeepTau2017v2VSjet", "byVLooseDeepTau2017v2VSmu", "byVVVLooseDeepTau2017v2VSe". Taus overlapping with fakeable electrons or fakeable muons within dR < 0.3 are not considered for the tau veto
      #False -> Gets Removed : True -> Passes veto
      for i in self.taus:
        if (i.pt > 20 and abs(i.eta) < 2.3 and abs(i.dxy) <= 1000 and abs(i.dz) <= 0.2 and i.idDecayModeNewDMs and i.decayMode in [0,1,2,10,11] and i.idDeepTau2017v2p1VSjet >= 16 and i.idDeepTau2017v2p1VSmu >= 1 and i.idDeepTau2017v2p1VSe >= 1):
          for fake in (self.muons_fakeable + self.electrons_fakeable):
            if deltaR(fake.eta, fake.phi, i.eta, i.phi) > 0.3:
              return False
      return True

    def Zmass_cut(self):
      #No pair of same-flavor, opposite-sign preselected leptons within  10GeV of the Z mass
      muons_pre = self.muons_pre; electrons_pre = self.electrons_pre;
      pre_leptons = muons_pre + electrons_pre
      pre1 = ROOT.TLorentzVector(); pre2 = ROOT.TLorentzVector(); pre_pair = ROOT.TLorentzVector();
      Z_mass = 80
      for i in pre_leptons:
        for j in pre_leptons:
          if i == j: continue
          if not ((i in muons_pre and j in muons_pre) or (i in electrons_pre and j in electrons_pre)): continue
          if (i.charge == j.charge): continue
          pre1.SetPtEtaPhiM(i.pt, i.eta, i.phi, i.mass)
          pre2.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)
          pre_pair = pre1 + pre2
          if (abs(pre_pair.M() - Z_mass) < 10):
            return False
      return True

    def single_electron_trigger(self, HLT):
      if Runyear == 2016:
        return (HLT.Ele27_WPTight_Gsf and HLT.Ele25_eta2p1_WPTight_Gsf and HLT.Ele27_eta2p1_WPLoose_Gsf)
      elif Runyear == 2017:
        return (HLT.Ele35_WPTight_Gsf and HLT.Ele32_WPTight_Gsf)
      elif Runyear == 2018:
        return (HLT.Ele32_WPTight_Gsf)
      else:
        return False

    def single_muon_trigger(self, HLT):
      if Runyear == 2016:
        return (HLT.IsoMu22 and HLT.IsoTkMu22 and HLT.IsoMu22_eta2p1 and HLT.IsoTkMu22_eta2p1 and HLT.IsoMu24 and HLT.IsoTkMu24)
      elif Runyear == 2017:
        return (HLT.IsoMu24 and HLT.IsoMu27)
      elif Runyear == 2018:
        return (HLT.IsoMu24 and HLT.IsoMu27)
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
        return (HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL and HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ and HLT.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL and HLT.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ)
      elif Runyear == 2017:
        return (HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 and HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)
      elif Runyear == 2018:
        return (HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)
      else:
        return False

    def muon_electron_trigger(self, HLT):
      if Runyear == 2016:
        return (HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL and HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ and HLT.Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL and HLT.Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ)
      elif Runyear == 2017:
        return (HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ and HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ and HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL and HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)
      elif Runyear == 2018:
        return (HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ and HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ and HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL)
      else:
        return False

    def met_filters(self, flag):
      if not (((Runyear == 2017 or Runyear == 2018) and flag.ecalBadCalibReducedMINIAODFilter) or (Runyear == 2016)):
        return False
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
        out.fillBranch("n_loose_ak4BJet", -999);
        out.fillBranch("n_medium_ak4BJet", -999);
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


        """
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
        """

