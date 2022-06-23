# Installing the Package:

Follow these instructions to checkout CMSSW_10_2_0 branch   
```
mkdir myFolder   
cmsrel CMSSW_10_2_0   
cd CMSSW_10_2_0/src   
cmsenv
git clone https://github.com/tahuang1991/nanoAOD-tools.git PhysicsTools/NanoAODTools   
git clone https://github.com/tahuang1991/HhhAnalysis.git   
cd HhhAnalysis
scram b -j 12   
```
tahuang1991/nanoAOD-tools is built based on cms-nanoAOD/nanoAOD-tools with minor change

# Produce NanoAOD sample v7
https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv7

# Running on data or simulation using NanoAOD sample to produce Ntuple with Dilepton and Single lepton selections:   

MC:
```
cd HhhAnalysis/python/NanoAOD
cmsRun postproc_local.py
```
# NanoAOD data sheet: 
https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html

# Notes on DiHiggsWWBBAnalyzer.cc with MiniAOD sample
This analyzer is targeted to run on DiHiggs, TTbar, DY samples.     
Sampletype is used to mark which sample it is running on.   
If you are using data, set Sampletype=0 and analyzer will skip gen-level selection   
If you are using signal sample, e.g. B6, you can set Sampletype=0, or Sampletype=6 which will check gen-level information   
If you are using TTbar, you can set Sampletype=0, or Sampletype=13 which will check gen-level information   
If you are using DY sample, you can set Sampletype=0, or Sampletype=14 which will check gen-level information( not implemented yet)   
The event will be saved only if it pass preselection at reco levle but we also sometimes keep track of gen-level information   
 

# Samples Location
https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/-/blob/master/Legacy/datasets.md
