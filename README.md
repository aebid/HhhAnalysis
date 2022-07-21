# Installing the Package:

Follow these instructions to checkout CMSSW_10_2_0 branch and add HeavyMassEstimator as submodule
```
mkdir myFolder   
cmsrel CMSSW_10_2_0   
cd CMSSW_10_2_0/src   
cmsenv
git clone https://github.com/tahuang1991/nanoAOD-tools.git PhysicsTools/NanoAODTools   
git clone https://github.com/tahuang1991/HhhAnalysis.git   
cd HhhAnalysis
git submoudle https://github.com/tahuang1991/HeavyMassEstimator
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
This analyzer is targeted to run on Data and CM samples like DiHiggs, TTbar, Drell-yan etc. 

Sampletype is used to mark which sample it is running on      
 - If you are using data, set Sampletype=0 and analyzer will skip gen-level selection   
 - If you are using signal sample, e.g. Radion260, you can set Sampletype=0, or Sampletype=100 which will check gen-level information   
 - If you are using TTbar, you can set Sampletype=0, or Sampletype=1 which will check gen-level information   
 - If you are using DY sample, you can set Sampletype=0, or Sampletype=2 which will check gen-level information( not implemented yet) 

The event will be saved only if it passes preselection at reco level or all generator levle particles are found  
 
# Run python version HME under PlottingTools


  - PlottingTools/runHME_HHbbWW_boosted.py is used to read in events from TTree and compute the HME for each event
  ```
  python runHME_HHbbWW_boosted.py -i infile.root -o outfile.root -nStart 0 -nEnd 100 -it 10000  -gp 0.18
  ```
  and the bash script version is following
  ```
  ./runHME_boosted.sh infile.root outfile.root 0 100 10000 0.18
  ```
  - PlottingTools/generateHMEJobs_HHbbWW_condor.py is to generate splitted condor jobs for one sample to run HME 


# Samples Location
https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/-/blob/master/Legacy/datasets.md
