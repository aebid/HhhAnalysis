### in  nanoAOD_processing.py

from bbWWProcessor_coffea import EventProcess
import awkward as ak


import time
startTime = time.time()


fname = "sync_2016_m750.root"
nLeps = 1 #Single lepton or Di Lepton channels
isSL = True
debug = 0

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

##object selection
#eventProcess.object_selection()
##SL
#eventProcess.SL_selection()





print("Events with 1 object comparison (my new coffea value) // (my old nanoAOD value [Tallinn value if different])")
print("Muons preselected: ", ak.sum(ak.any(eventProcess.muons.preselected, axis=1)), " // 93605")
print("Muons fakeable: ", ak.sum(ak.any(eventProcess.muons.fakeable, axis=1)), " // 81978")
print("Muons tight: ", ak.sum(ak.any(eventProcess.muons.tight, axis=1)), " // 78340")


print("Electrons preselected: ", ak.sum(ak.any(eventProcess.electrons.preselected, axis=1)), " // NA")
print("Electrons cleaned: ", ak.sum(ak.any(eventProcess.electrons.cleaned, axis=1)), " // 75430")
print("Electrons fakeable: ", ak.sum(ak.any(eventProcess.electrons.fakeable, axis=1)), " // 58833")
print("Electrons tight: ", ak.sum(ak.any(eventProcess.electrons.tight, axis=1)), " // 56395")


print("AK4 Jets preselected: ", ak.sum(ak.any(eventProcess.ak4_jets.preselected, axis=1)), " // NA")
print("AK4 Jets cleaned: ", ak.sum(ak.any(eventProcess.ak4_jets.cleaned, axis=1)), " // 144403(144446)")
print("AK4 Jets loose Btag: ", ak.sum(ak.any(eventProcess.ak4_jets.loose_btag, axis=1)), " // NA")
print("AK4 Jets medium Btag: ", ak.sum(ak.any(eventProcess.ak4_jets.medium_btag, axis=1)), " // NA")

print("AK8 Jets preselected: ", ak.sum(ak.any(eventProcess.ak8_jets.preselected, axis=1)), " // 77501")
print("AK8 Jets cleaned: ", ak.sum(ak.any(eventProcess.ak8_jets.cleaned, axis=1)), " // 69384")
print("AK8 Jets Btag: ", ak.sum(ak.any(eventProcess.ak8_jets.btag, axis=1)), " // 54065(53678)")


print("N events: ", len(eventProcess.events))
print("N single events: ", ak.sum(eventProcess.events.single_lepton))



executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
