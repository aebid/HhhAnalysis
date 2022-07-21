# The folder is used to make MC vs Data plots and Compuate HME and MT2

## Make the Data vs MC plots:

1) Adjust "PlotterProducer.py" with the correct input files:
  a) in HhhAnalysis/CutFlowAnalyzer/test/DATA_MINIAODtoNtuple you have a file called multicrab3.py
  b) in 'multicrab3.py' set OnlySubmitCRAB to False
  c) "python multicrab3.py" will create a file called 'for_plotter.py'
  d) copy the content of for_plotter.py in PlotterProducer.py (under # Data)
  e) Redo these steps for the file "multicrab3.py" in HhhAnalysis/CutFlowAnalyzer/test/MC_MINIAODtoNtuple 
2) source Send_PlotterProducer.sh
3) Once the jobs are done you can run: 'python FinalPlotter.py -b'

## HME
To Checkout python version Heavy mass estimator(HME) by
```
git submodule add https://github.com/tahuang1991/HeavyMassEstimator
```
Python scripts to compute HME include runHME_HHbbWW.py, runHME_HHbbWW_boosted.py, runHME_HHbbWW_general.py etc and the scripts to generate condor/slrm batch jobs include generateHMEJobs.py, generateHMEJobs_HHbbWW.py, generateHMEJobs_HHZZBB.py, generateHMEJobs_HHbbWW_condor.py etc. bash scripts like runHME_boosted.sh is used in condor job to run runHME_HHbbWW_boosted.py

Here is the example to run runHME_HHbbWW_boosted.py
  ```
  python runHME_HHbbWW_boosted.py -i infile.root -o outfile.root -nStart 0 -nEnd 100 -it 10000  -gp 0.18
  ```
Above exmaple is to run HME from eventid=0 to eventid=99 from infile.root. The number iteration is 10000 and the method to get average reco mass around peak uses gp=0.18. Finally it would write TTree with reco HME as new branch into outfile.root



## MT2
mt2_* files are used to compute the MT2 variable. The Swig is used to convert C++ version MT2 code into python version 

