### in  nanoAOD_processing.py

from bbWWProcessor_coffea import EventProcess
import awkward as ak


import time
startTime = time.time()


fname = "sync_2016_m750.root"
nLeps = 1 #Single lepton or Di Lepton channels
isSL = True
debug = 1

Runyear = 2016
isMC = True


eventProcess = EventProcess(fname, "RadionDLM750", isMC, Runyear, isSL, debug)

eventProcess.add_conept()
eventProcess.link_jets()
eventProcess.muon_selection()
eventProcess.electron_selection()
eventProcess.ak4_jet_selection()
eventProcess.ak8_jet_selection()



eventProcess.single_lepton_category()
eventProcess.double_lepton_category()

if debug:
    eventProcess.print_object_selection()
    eventProcess.print_event_selection()

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))


eventProcess.create_df()

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
