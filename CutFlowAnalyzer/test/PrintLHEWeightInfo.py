from DataFormats.FWLite import Events, Handle
events = Events("/eos/cms/store/mc/RunIIFall17DRPremix/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/00000/F2C3E077-6650-E811-879E-FA163EEB9AB3.root")
for event in events: 
  break
handle = Handle("LHEEventProduct")
event.getByLabel("externalLHEProducer", handle)
lhe = handle.product()
for i in lhe.weights(): 
  print i.id, i.wgt
