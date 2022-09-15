import sys
import os

from ROOT import TFile, TDirectory, TTree, gROOT
from plotNGenParticles import *

## run quiet mode
sys.argv.append( '-b' )
gROOT.SetBatch(1)

class HHbbWWPlotter():
  def __init__(self, inputFile, year, baseDir, sample, sampleShortName, isSL=True, isData=False, treename="syncTree"):
    self.inputFile = inputFile
    self.baseDir = baseDir
    if not os.path.isdir(self.baseDir): 
        os.mkdir(self.baseDir)
    self.ext = ".pdf"
    self.isSL = isSL
    self.isData = isData
    self.shortName = sampleShortName
    self.legend = sampleShortName
    self.sample = sample
    self.runyear = year
    self.file = TFile.Open(self.inputFile)
    self.tree = (self.file).Get(treename)
  def setLegend(self, leg):
    self.legend = leg
  def getTree(self):
    return self.tree




##non-res, LO 
#rfile2017DL_nonres="GluGluToHHTo2B2VTo2L2Nu_node_10_run2017_2ACB1EB15515_Friend.root"
#rfile2016DL_nonres="GluGluToHHTo2B2VTo2L2Nu_node_11_run2016_52CA4A8443D3_Friend.root"
#rfile2017SL_nonres="GluGluToHHTo2B2WToLNu2J_node_11_run2017_AC7A17C0E21B_Friend.root"
#rfile2016SL_nonres="GluGluToHHTo2B2WToLNu2J_node_6_run2016_cfg_NANO_5_Friend.root"
##non-res, NLO 
datapath = "/Users/taohuang/Documents/DiHiggs/20220819_HHbbWW_genPart/"
rfile2016SL_nonres = datapath+"GluGluToHHTo2B2VLNu2J_node_cHHH1_run2016_EC732DAFA0D8_Friend.root"
rfile2017SL_nonres = datapath+"GluGluToHHTo2B2VLNu2J_node_cHHH0_run2017_AB1DFE7A91C9_Friend.root"
rfile2016DL_nonres = datapath+"GluGluToHHTo2B2VTo2L2Nu_node_cHHH1_run2016_7DF388FAD247_Friend.root"
rfile2017DL_nonres = datapath+"GluGluToHHTo2B2VTo2L2Nu_node_cHHH1_run2017_CF0FDD252EF0_Friend.root"
rfile2017DLtest    = datapath+"GluGluToHHTo2B2VTo2L2Nu_node_10_run2017_C1F513F6A155_testfile_Friend.root"
rfile2016DL_res    = datapath+"sync_2016_m750_Friend.root"
rfile2017DL_res    = datapath+"sync_2017_m750_Friend.root"
rfile2016SL_res    = datapath+"NanoAOD_SingleLepton_M500_run2016_Friend.root"
rfile2017SL_res    = datapath+"NanoAODproduction_2017_Singlelepton_RadionM500_Friend.root"
##2016SL non-resnonance

rfile2016SL_nonres_plt = HHbbWWPlotter(rfile2016SL_nonres, 2016, "./","GluGluToHHTo2B2VLNu2J_node_cHHH1","Nonres_cHHH1", True)
rfile2017SL_nonres_plt = HHbbWWPlotter(rfile2017SL_nonres, 2017, "./","GluGluToHHTo2B2VLNu2J_node_cHHH0","Nonres_cHHH0", True)
rfile2016DL_nonres_plt = HHbbWWPlotter(rfile2016DL_nonres, 2016, "./","GluGluToHHTo2B2VTo2L2J_node_cHHH1","Nonres_cHHH1",False)
rfile2017DL_nonres_plt = HHbbWWPlotter(rfile2017DL_nonres, 2017, "./","GluGluToHHTo2B2VTo2L2J_node_cHHH1","Nonres_cHHH1",False)
rfile2016SL_res_plt = HHbbWWPlotter(rfile2016SL_res, 2016, "./","GluGluToRadionToHHTo2B2VLNu2J_M500","RadionM500", True)
rfile2017SL_res_plt = HHbbWWPlotter(rfile2017SL_res, 2017, "./","GluGluToRadionToHHTo2B2VLNu2J_M500","RadionM500", True)
rfile2016DL_res_plt = HHbbWWPlotter(rfile2016DL_res, 2016, "./","GluGluToRadionToHHTo2B2VTo2L2J_M750","RadionM750",False)
rfile2017DL_res_plt = HHbbWWPlotter(rfile2017DL_res, 2017, "./","GluGluToRadionToHHTo2B2VTo2L2J_M750","RaidonM750",False)

plotAllNGen(rfile2016SL_nonres_plt)
plotAllNGen(rfile2017SL_nonres_plt)
plotAllNGen(rfile2016DL_nonres_plt)
plotAllNGen(rfile2017DL_nonres_plt)
plotAllNGen(rfile2016SL_res_plt)
plotAllNGen(rfile2017SL_res_plt)
plotAllNGen(rfile2016DL_res_plt)
plotAllNGen(rfile2017DL_res_plt)
