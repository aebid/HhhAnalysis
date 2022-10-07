import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
import ROOT
from coffea.nanoevents.methods import vector




class EventProcess():
    def __init__(self, inputFile, sampleName, isMC, runYear, isSL, debug=0):
        self.fname = inputFile
        self.sampleName = samplleName
        self.isMC  = isMC
        self.runYear = runYear
        self.isSL = isSL
        self.debug  = dehug
        print("Starting NanoAOD processing")


        #fname = "sync_2016_m750.root"
        events = NanoEventsFactory.from_root(self.fname, schemaclass=NanoAODSchema.v7).events()
        self.nLeps = 1 #Single lepton or Di Lepton channels
        if not self.isSL: self.nLeps = 2
        debug = 0


        self.muons = ak.pad_none(events.Muon, 1)
        self.electrons = ak.pad_none(events.Electron, 1)
        self.taus = ak.pad_none(events.Tau, 1)
        self.ak4_jets = ak.pad_none(events.Jet, 1)
        self.ak8_jets = ak.pad_none(events.FatJet, 1)
        self.ak8_subjets = ak.pad_none(events.SubJet, 1)
        self.HLT = events.HLT
        self.flag = events.Flag


        if debug > 0:
            print("Muons: ",       self.muons)
            print("Electrons: ",   self.electrons)
            print("Taus: ",        self.taus)
            print("AK4 Jets: ",    self.ak4_jets)
            print("AK8 Jets: ",    self.ak8_jets)
            print("AK8 SubJets: ", self.ak8_subjets)
            print("HLT: ",         self.HLT)
        from object_selection import object_selection
        from event_selection import SL_selection,DL_selection



### in  nanoAOD_processing.py

from bbWWProcessor_coffea import EventProcess


fname = "sync_2016_m750.root"
nLeps = 1 #Single lepton or Di Lepton channels
isSL = True
debug = 0

Runyear = 2016
isMC = True


eventProcess = EventProcess(fname, "RadionDLM750", isMC, Runyear, isSL, debug)
##objec selection
eventProcess.object_selection()
##SL
eventProcess.SL_selection()
