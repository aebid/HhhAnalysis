import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
import ROOT
from coffea.nanoevents.methods import vector


import object_selection
import event_selection

class EventProcess():
    def __init__(self, inputFile, sampleName, isMC, Runyear, isSL, debug=0):
        self.fname = inputFile
        self.sampleName = sampleName
        self.isMC  = isMC
        self.Runyear = Runyear
        self.isSL = isSL
        self.debug  = debug
        print("Starting NanoAOD processing")


        #fname = "sync_2016_m750.root"
        events = NanoEventsFactory.from_root(self.fname, schemaclass=NanoAODSchema.v7).events()
        self.nLeps = 1 #Single lepton or Di Lepton channels
        if not self.isSL: self.nLeps = 2
        debug = 0

        self.events = events
        self.muons = ak.pad_none(events.Muon, 1)
        self.electrons = ak.pad_none(events.Electron, 1)
        self.taus = ak.pad_none(events.Tau, 1)
        self.ak4_jets = ak.pad_none(events.Jet, 1)
        self.ak8_jets = ak.pad_none(events.FatJet, 1)
        self.ak8_subjets = ak.pad_none(events.SubJet, 1)
        self.HLT = events.HLT
        self.flag = events.Flag

        self.jetDeepJet_WP_loose  = [0.0613, 0.0521, 0.0494]
        self.jetDeepJet_WP_medium = [0.3093, 0.3033, 0.2770]
        self.jetDeepJet_WP_tight  = [0.7221, 0.7489, 0.7264]
        self.PFJetID = [1, 2, 2]

        if debug > 0:
            print("Muons: ",       self.muons)
            print("Electrons: ",   self.electrons)
            print("Taus: ",        self.taus)
            print("AK4 Jets: ",    self.ak4_jets)
            print("AK8 Jets: ",    self.ak8_jets)
            print("AK8 SubJets: ", self.ak8_subjets)
            print("HLT: ",         self.HLT)
        #from object_selection import object_selection
        #from event_selection import SL_selection,DL_selection

    def add_conept(self):
        return object_selection.add_conept(self)
    def link_jets(self):
        return object_selection.link_jets(self)
    def muon_selection(self):
        return object_selection.muon_selection(self)
    def electron_selection(self):
        return object_selection.electron_selection(self)
    def ak4_jet_selection(self):
        return object_selection.ak4_jet_selection(self)
    def ak8_jet_selection(self):
        return object_selection.ak8_jet_selection(self)

    def single_lepton_category(self):
        return event_selection.single_lepton_category(self)


"""
### in  nanoAOD_processing.py

from bbWWProcessor_coffea import EventProcess


fname = "sync_2016_m750.root"
nLeps = 1 #Single lepton or Di Lepton channels
isSL = True
debug = 0

Runyear = 2016
isMC = True


eventProcess = EventProcess(fname, "RadionDLM750", isMC, Runyear, isSL, debug)

#eventProcess.cone_pt()
#print(eventProcess.muons)
#print(eventProcess.muons.conept)

##object selection
#eventProcess.object_selection()
##SL
#eventProcess.SL_selection()
"""
