import ROOT


runyear = 2017
myfilelist = ["sync/sync_2016_m750_Friend.root", "sync/sync_2017_m750_Friend.root"]
tallinfilelist = ["/afs/cern.ch/user/k/kaehatah/public/sync_legacy_bbww/era_2016/sync_bbww_Tallinn_2016_v21.root", "/afs/cern.ch/user/k/kaehatah/public/sync_legacy_bbww/era_2017/sync_bbww_Tallinn_2017_v13.root"]



for runyear in [2016, 2017]:
  f = ROOT.TFile(myfilelist[runyear-2016])
  tree = f.Get("syncTree")
  tree_1l_SR = f.Get("syncTree_hhbb1l_SR")
  tree_1l_Fake = f.Get("syncTree_hhbb1l_Fake")
  tree_2l_SR = f.Get("syncTree_hhbb2l_SR")
  tree_2l_Fake = f.Get("syncTree_hhbb2l_Fake")

  n_mu = tree.GetEntries("n_presel_mu")
  n_mu_fake = tree.GetEntries("n_fakeablesel_mu")
  n_mu_tight = tree.GetEntries("n_mvasel_mu")
  n_ele = tree.GetEntries("n_presel_ele")
  n_ele_fake = tree.GetEntries("n_fakeablesel_ele")
  n_ele_tight = tree.GetEntries("n_mvasel_ele")
  n_ak4 = tree.GetEntries("n_presel_ak4Jet")
  n_ak8 = tree.GetEntries("n_presel_ak8Jet")

  n_2l_SR = tree_2l_SR.GetEntries()
  n_2l_Fake = tree_2l_Fake.GetEntries()
  n_1l_SR = tree_1l_SR.GetEntries()
  n_1l_Fake = tree_1l_Fake.GetEntries()



  f2 = ROOT.TFile(tallinfilelist[runyear-2016])
  tree2 = f2.Get("syncTree")
  tree_1l_SR_Tall = f2.Get("syncTree_hhbb1l_SR")
  tree_1l_Fake_Tall = f2.Get("syncTree_hhbb1l_Fake")
  tree_2l_SR_Tall = f2.Get("syncTree_hhbb2l_SR")
  tree_2l_Fake_Tall = f2.Get("syncTree_hhbb2l_Fake")

  n_mu_Tall = tree2.GetEntries("n_presel_mu")
  n_mu_fake_Tall = tree2.GetEntries("n_fakeablesel_mu")
  n_mu_tight_Tall = tree2.GetEntries("n_mvasel_mu")
  n_ele_Tall = tree2.GetEntries("n_presel_ele")
  n_ele_fake_Tall = tree2.GetEntries("n_fakeablesel_ele")
  n_ele_tight_Tall = tree2.GetEntries("n_mvasel_ele")
  n_ak4_Tall = tree2.GetEntries("n_presel_ak4Jet")
  n_ak8_Tall = tree2.GetEntries("n_presel_ak8Jet")

  n_2l_SR_Tall = tree_2l_SR_Tall.GetEntries()
  n_2l_Fake_Tall = tree_2l_Fake_Tall.GetEntries()
  n_1l_SR_Tall = tree_1l_SR_Tall.GetEntries()
  n_1l_Fake_Tall = tree_1l_Fake_Tall.GetEntries()





  print "Runyear = ", runyear
  print "Category 	: My Result 	: Tallinn"
  print "1 pre mu 	: ", n_mu, " 	: ", n_mu_Tall
  print "1 fake mu 	: ", n_mu_fake, " 	: ", n_mu_fake_Tall
  print "1 tight mu 	: ", n_mu_tight, " 	: ", n_mu_tight_Tall
  print "1 pre ele 	: ", n_ele, " 	: ", n_ele_Tall
  print "1 fake ele 	: ", n_ele_fake, " 	: ", n_ele_fake_Tall
  print "1 tight ele 	: ", n_ele_tight, " 	: ", n_ele_tight_Tall
  print "1 pre ak4 	: ", n_ak4, " 	: ", n_ak4_Tall
  print "1 pre ak8 	: ", n_ak8, " 	: ", n_ak8_Tall

  print "n 2l SR 	: ", n_2l_SR, " 	: ", n_2l_SR_Tall
  print "n 2l Fake	: ", n_2l_Fake, " 	: ", n_2l_Fake_Tall
  print "n 1l SR 	: ", n_1l_SR, " 	: ", n_1l_SR_Tall
  print "n 1l Fake 	: ", n_1l_Fake, " 	: ", n_1l_Fake_Tall


  my_events = []; tal_events = []
  #mytree = tree_2l_Fake; taltree = tree_2l_Fake_Tall
  #mytree = tree_2l_SR; taltree = tree_2l_SR_Tall
  #mytree = tree_1l_Fake; taltree = tree_1l_Fake_Tall
  #mytree = tree_1l_SR; taltree = tree_1l_SR_Tall
  mytree = tree; taltree = tree2

  if runyear == 2017:
    """
    for i in range(taltree.GetEntries()):
      taltree.GetEntry(i); tal_events.append(taltree.event)
    for i in range(mytree.GetEntries()):
      mytree.GetEntry(i); my_events.append(mytree.event)
    """
    for i in range(taltree.GetEntries()):
      taltree.GetEntry(i); 
      if taltree.n_presel_ak4Jet > 0: tal_events.append(taltree.event)
    for i in range(mytree.GetEntries()):
      mytree.GetEntry(i); 
      if mytree.n_presel_ak4Jet > 0: my_events.append(mytree.event)
  print "AK4 jet count"

  hnm = [x for x in tal_events if x not in my_events]
  nmh = [x for x in my_events if x not in tal_events]
  #print "Single Lepton Fake Region"
  print "Him not in me ", len(hnm), hnm
  print "Me not in him ", len(nmh), nmh
