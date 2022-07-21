import os, sys, random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
from math import *
import argparse
import numpy as np

sys.path.append(os.path.join(os.path.abspath('..'), 'HeavyMassEstimator/python'))
from HeavyMassEstimatorHHZZBB import *

import Samplelist 
execfile("start.py")
execfile("functions.py")
#Creating folders and parameters
def AddFiles(chain, inputDir):
    theInputFiles = []
    for d in range(len(inputDir)):
        my_dir = inputDir[d]
        if os.path.isdir(my_dir):
	    ls = os.listdir(my_dir)
	    print "mydir ", my_dir, " allfile ",ls
	    for x in ls:
	       if x.endswith('.root'):
		    print my_dir + x	       
		    chain.Add(my_dir + x)
	    #chain.Add(my_dir + x for x in ls if x.endswith('.root'))
        elif os.path.isfile(my_dir):
	    chain.Add(my_dir)
	    continue
        elif d==len(inputDir)-1:
	    print "ERROR: No input files were selected"
	    exit()
doTest = False
doHME = True
#refPDF = ROOT.TFile("REFPDFPU40_ZZBB.root","READ")
#onshellZmasspdf = refPDF.Get("onshellZmasspdf")
#offshellZmasspdf = refPDF.Get("offshellZmasspdf")
##onshellZmasspdf.Print(); offshellZmasspdf.Print()
#recobjetrescalec1pdfPU40 = refPDF.Get("recobjetrescalec1pdfPU40")

tree_name="DiHiggsWWBBAna/evtree"
benchmarks = ["ttV", "Wjet", "singTop", "VV", "DY", "TTbar", "Data"]
#whichSample = "Graviton_0419"

parser = argparse.ArgumentParser(description='runHME')
parser.add_argument("-n", "--njobs", dest="njobs", type=int, default=100, help="total splitted jobs. [Default: 100]")
parser.add_argument("-i", "--ijob", dest="ijob", type=int, default=0, help="ith job in splitting. [Default: 0]")
parser.add_argument("-o", "--outputdir", dest="outputdir", default="./", help="output directory [Default: ./]")
parser.add_argument("-jt", "--jobtype", dest="jobtype", default="radion", help="sample type to run. [Default: radion]")
parser.add_argument("-it", "--iterations", dest="iterations",type=int, default=10000, help="iteration number used HME [Default: 100000]")
args = parser.parse_args()
whichSample = args.jobtype
makeHadd = "HaddNo"
"""
#print "Executing: python", sys.argv[0] , "-b", sys.argv[2], sys.argv[3], "(Arg1=makeHadd: HaddYes or HaddNo /|\ Arg2=whichSample: TT, DY, VV, sT, Wjet, ttV, Data)"
#makeHadd = sys.argv[2]
#whichSample = sys.argv[3]
if( makeHadd!="HaddYes" and makeHadd!="HaddNo" ): print "WARNING! 1st argument has to be HaddYes or HaddNo"; sys.exit()
#if( whichSample!="TT" and whichSample!="DY" and whichSample!="VV" and whichSample!="sT" and whichSample!="Wjet" and whichSample!="ttV" and whichSample!="Data" ):  print "WARNING! 2nd argument have to be TT, DY, VV, singTop, Wjet, ttV, or Data"; sys.exit()
######################################
## define shell commands
######################################
Find_str = []; this_cat = ""; this_hadd = ""; this_NtotPath = ""
for i in range(len(Samplelist.outAnalist[whichSample])):
  Find_str.append("find %s | grep root | grep -v failed > HADD/%s_%d.txt"%(Samplelist.outAnalist[whichSample][i], Samplelist.sampleFullName[whichSample], i))
    #Find_str.append("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg/170405_172215 | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt")
this_cat = "cat HADD/%s_* > HADD/%s.txt"%(Samplelist.sampleFullName[whichSample], whichSample)
if whichSample != "DATA":
    this_NtotPath = "/fdata/hepx/store/user/%s/%s/crab_%s.root"%(user,Samplelist.sampleFullName[whichSample], Samplelist.sampleFullName[whichSample])
    this_hadd     = "hadd -T -f -k %s @HADD/%s.txt"%(this_NtotPath, whichSample)
    outpath = "/fdata/hepx/store/user/%s/%s"%(user,Samplelist.sampleFullName[whichSample])
    os.system("mkdir -p %s" % outpath)

######################################
# Running commands
######################################
for this_find in Find_str:
  print "this find ",this_find
  os.system(this_find)

os.system(this_cat)
print "this_hadd ", this_hadd
if (makeHadd=="HaddYes" and whichSample!="Data"): os.system(this_hadd)
"""
######################################
#this_NtotPath  = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/src/HhhAnalysis/CutFlowAnalyzer/test/out_ana_ZZBBJJ.root"
#this_NtotPath = "/fdata/hepx/store/user/tahuang/MC_Hhh_analysis/GluGluToRadionToHHTo2B2ZTo2L2J_M-400_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2ZTo2L2J_M-400_narrow_13TeV-madgraph-v2/170626_081846/0000/combined/out_ana.root"
#this_NtotPath = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/src/HhhAnalysis/CutFlowAnalyzer/test/out_ana_10k_graviton.root"
#this_NtotPath = "/fdata/hepx/store/user/taohuang/DiHiggsAnalysisSample/out_ann_%s_20160411.root"%whichSample
this_NtotPath = Samplelist.outAnalist[Samplelist.sampleFullName[whichSample]]
print "this_NtotPath ",this_NtotPath
#if( whichSample != "Data" ):
#  Ntot_path = this_NtotPath
#  MyFile =  ROOT.TFile.Open(Ntot_path,"read");
#  h_prehlt  =  ROOT.TH1F(MyFile.Get("TriggerResults/hevent_filter")); nTOT_prehlt = h_prehlt.GetBinContent(2)
#  h_posthlt =  ROOT.TH1F(MyFile.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt = h_posthlt.GetBinContent(2);
#else:
#  nTOT_prehlt = 0.; nTOT_posthlt = 0.
nTOT_prehlt = 1.; nTOT_posthlt = 1.0
TCha = ROOT.TChain(tree_name)
#with open(this_cat.split(">")[1].split(" ")[1],"r") as f:
#  for line in f:
#    if not line.isspace():
#      TCha.Add(str(line[:-1]))
#TCha.Add("/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/src/HhhAnalysis/CutFlowAnalyzer/test/out_ana_10k_graviton.root")
#TCha.Add(this_NtotPath)
AddFiles(TCha, this_NtotPath)
print whichSample, "TChain has", TCha.GetEntries(), "entries."
TotalEv = TCha.GetEntries()
nStart = TotalEv/args.njobs*args.ijob
nEnd = TotalEv/args.njobs*(args.ijob+1)
if args.ijob == args.njobs-1:
    nEnd = -1
#f = ROOT.TFile('/fdata/hepx/store/user/%s/Hhh_For_Plotting/'%user + whichSample + '_ijob%d.root'%args.ijob,'recreate'); f.cd()
f = ROOT.TFile(args.outputdir + whichSample + '_ijob%d.root'%args.ijob,'recreate'); f.cd()
print "output ",args.outputdir + whichSample + '_ijob%d.root'%args.ijob

ptbinSF= [20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0, 200.0, 300.]
etabinSF = [0, 0.9, 1.2, 2.1, 2.4]
myptbinSF = np.asarray(ptbinSF)
myetabinSF = np.asarray(etabinSF)

# PRESELECTION
# Weights
h_pre_Nev_preHLT         = ROOT.TH1F("h_pre_Nev_preHLT","",1,-0.5,0.5);          h_pre_Nev_preHLT.GetXaxis().SetTitle("#Events (pre HLT)"); h_pre_Nev_preHLT.SetBinContent(1,nTOT_prehlt)
h_pre_Nev_posHLT         = ROOT.TH1F("h_pre_Nev_posHLT","",1,-0.5,0.5);          h_pre_Nev_posHLT.GetXaxis().SetTitle("#Events (post HLT)"); h_pre_Nev_posHLT.SetBinContent(1,nTOT_posthlt)
h_pre_XsecBr             = ROOT.TH1F("h_pre_XsecBr","",1,-0.5,0.5);              h_pre_XsecBr.GetXaxis().SetTitle("XsecBr");
h_pre_muon1_triggerSF    = ROOT.TH2F("h_pre_muon1_triggerSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon1_triggerSF.GetXaxis().SetTitle("#mu 1 Trigger SF");
h_pre_muon1_isoSF        = ROOT.TH2F("h_pre_muon1_isoSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF); h_pre_muon1_isoSF.GetXaxis().SetTitle("#mu 1 Iso SF");
h_pre_muon1_idSF         = ROOT.TH2F("h_pre_muon1_idSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon1_idSF.GetXaxis().SetTitle("#mu 1 ID SF");
h_pre_muon1_trackingSF   = ROOT.TH1F("h_pre_muon1_trackingSF","",50,-2.4,2.4); h_pre_muon1_trackingSF.GetXaxis().SetTitle("#mu 1 tracking SF");
h_pre_muon1_SF_bg2       = ROOT.TH2F("h_pre_muon1_SF_bg2","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon1_SF_bg2.GetXaxis().SetTitle("#mu SF bg2");
h_pre_muon1_SF_bg1       = ROOT.TH1F("h_pre_moun1_SF_bg1","",50,-2.4,2.4); h_pre_muon1_SF_bg1.GetXaxis().SetTitle("#mu SF bg1");
h_pre_muon2_triggerSF    = ROOT.TH2F("h_pre_muon2_triggerSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon2_triggerSF.GetXaxis().SetTitle("#mu 1 Trigger SF");
h_pre_muon2_isoSF        = ROOT.TH2F("h_pre_muon2_isoSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);      h_pre_muon2_isoSF.GetXaxis().SetTitle("#mu 1 Iso SF");
h_pre_muon2_idSF         = ROOT.TH2F("h_pre_muon2_idSF","",len(etabinSF)-1,myetabinSF, len(ptbinSF)-1, myptbinSF);       h_pre_muon2_idSF.GetXaxis().SetTitle("#mu 1 ID SF");
h_pre_muon2_trackingSF   = ROOT.TH1F("h_pre_muon2_trackingSF","",50,-2.4,2.4); h_pre_muon2_trackingSF.GetXaxis().SetTitle("#mu 1 tracking SF");
h_pre_muon2_SF_bg2       = ROOT.TH2F("h_pre_muon2_SF_bg2","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon2_SF_bg2.GetXaxis().SetTitle("#mu SF bg2");
h_pre_muon2_SF_bg1       = ROOT.TH1F("h_pre_moun2_SF_bg1","",50,-2.4,2.4); h_pre_muon1_SF_bg1.GetXaxis().SetTitle("#mu SF bg1");
# Regression variables
h_pre_numOfVertices      = ROOT.TH1F("h_pre_numOfVertices","",50,0.,50.);      h_pre_numOfVertices.GetXaxis().SetTitle("");
h_pre_b1jet_mt           = ROOT.TH1F("h_pre_b1jet_mt","",50,0.,100.);          h_pre_b1jet_mt.GetXaxis().SetTitle("");
h_pre_b1jet_leadTrackPt  = ROOT.TH1F("h_pre_b1jet_leadTrackPt","",50,0.,100.); h_pre_b1jet_leadTrackPt.GetXaxis().SetTitle("");
h_pre_b1jet_leptonPtRel  = ROOT.TH1F("h_pre_b1jet_leptonPtRel","",50,0.,10.);  h_pre_b1jet_leptonPtRel.GetXaxis().SetTitle("");
h_pre_b1jet_leptonPt     = ROOT.TH1F("h_pre_b1jet_leptonPt","",50,0.,100.);    h_pre_b1jet_leptonPt.GetXaxis().SetTitle("");
h_pre_b1jet_leptonDeltaR = ROOT.TH1F("h_pre_b1jet_leptonDeltaR","",50,0.,5.);  h_pre_b1jet_leptonDeltaR.GetXaxis().SetTitle("");
h_pre_b1jet_neHEF        = ROOT.TH1F("h_pre_b1jet_neHEF","",50,0.,2.);         h_pre_b1jet_neHEF.GetXaxis().SetTitle("");
h_pre_b1jet_neEmEF       = ROOT.TH1F("h_pre_b1jet_neEmEF","",50,0.,2.);        h_pre_b1jet_neEmEF.GetXaxis().SetTitle("");
h_pre_b1jet_vtxNtracks   = ROOT.TH1F("h_pre_b1jet_vtxNtracks","",50,0.,50.);   h_pre_b1jet_vtxNtracks.GetXaxis().SetTitle("");
h_pre_b1jet_vtxPt        = ROOT.TH1F("h_pre_b1jet_vtxPt","",50,0.,100.);       h_pre_b1jet_vtxPt.GetXaxis().SetTitle("");
h_pre_b1jet_vtxMass      = ROOT.TH1F("h_pre_b1jet_vtxMass","",50,0.,50.);      h_pre_b1jet_vtxMass.GetXaxis().SetTitle("");
h_pre_b1jet_vtx3DSig     = ROOT.TH1F("h_pre_b1jet_vtx3DSig","",50,0.,5.);      h_pre_b1jet_vtx3DSig.GetXaxis().SetTitle("");
h_pre_b1jet_vtx3DVal     = ROOT.TH1F("h_pre_b1jet_vtx3DVal","",50,0.,5.);      h_pre_b1jet_vtx3DVal.GetXaxis().SetTitle("");

# Muons Gen
h_gen_pre_MU1_pt         = ROOT.TH1F("h_gen_pre_MU1_pt","",50,10.,350.);    h_gen_pre_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]");    
h_gen_pre_MU2_pt         = ROOT.TH1F("h_gen_pre_MU2_pt","",50,10.,350.);    h_gen_pre_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); 
h_gen_pre_MU1_eta        = ROOT.TH1F("h_gen_pre_MU1_eta","",50,-3.,3.);     h_gen_pre_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta");          
h_gen_pre_MU2_eta        = ROOT.TH1F("h_gen_pre_MU2_eta","",50,-3.,3.);     h_gen_pre_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta");       
h_gen_pre_mass_l1l2      = ROOT.TH1F("h_gen_pre_mass_l1l2","",40,20.,150.); h_gen_pre_mass_l1l2.GetXaxis().SetTitle("m(l,l)");               
h_gen_pre_dR_l1l2        = ROOT.TH1F("h_gen_pre_dR_l1l2","",50,0.,5.);      h_gen_pre_dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)");           
# Muons Reco
h_pre_MU1_pt             = ROOT.TH1F("h_pre_MU1_pt","",50,10.,350.);    h_pre_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]");    
h_pre_MU2_pt             = ROOT.TH1F("h_pre_MU2_pt","",50,10.,350.);    h_pre_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); 
h_pre_MU1_eta            = ROOT.TH1F("h_pre_MU1_eta","",50,-3.,3.);     h_pre_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta");          
h_pre_MU2_eta            = ROOT.TH1F("h_pre_MU2_eta","",50,-3.,3.);     h_pre_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta");       
h_pre_mass_l1l2          = ROOT.TH1F("h_pre_mass_l1l2","",40,20.,150.); h_pre_mass_l1l2.GetXaxis().SetTitle("m(l,l)");               
h_pre_dR_l1l2            = ROOT.TH1F("h_pre_dR_l1l2","",50,0.,5.);      h_pre_dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)");           
# B-jets Gen
h_gen_pre_J1_pt          = ROOT.TH1F("h_gen_pre_J1_pt","",50,20.,350.);     h_gen_pre_J1_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]");     
h_gen_pre_J2_pt          = ROOT.TH1F("h_gen_pre_J2_pt","",50,20.,350.);     h_gen_pre_J2_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]");  
h_gen_pre_J1_eta         = ROOT.TH1F("h_gen_pre_J1_eta","",50,-3.,3.);      h_gen_pre_J1_eta.GetXaxis().SetTitle("Lead. jet #eta");           
h_gen_pre_J2_eta         = ROOT.TH1F("h_gen_pre_J2_eta","",50,-3.,3.);      h_gen_pre_J2_eta.GetXaxis().SetTitle("Sublead. jet #eta");        
h_gen_pre_mass_b1b2      = ROOT.TH1F("h_gen_pre_mass_b1b2","",40,40.,400.); h_gen_pre_mass_b1b2.GetXaxis().SetTitle("m(j,j)");      
h_gen_pre_dR_b1b2        = ROOT.TH1F("h_gen_pre_dR_b1b2","",50,0.,6.);      h_gen_pre_dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)");           
# B-jets Reco
h_pre_J1_pt              = ROOT.TH1F("h_pre_J1_pt","",50,20.,350.);     h_pre_J1_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]");     
h_pre_J2_pt              = ROOT.TH1F("h_pre_J2_pt","",50,20.,350.);     h_pre_J2_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]");  
h_pre_J1_eta             = ROOT.TH1F("h_pre_J1_eta","",50,-3.,3.);      h_pre_J1_eta.GetXaxis().SetTitle("Lead. jet #eta");           
h_pre_J2_eta             = ROOT.TH1F("h_pre_J2_eta","",50,-3.,3.);      h_pre_J2_eta.GetXaxis().SetTitle("Sublead. jet #eta");        
h_pre_mass_b1b2          = ROOT.TH1F("h_pre_mass_b1b2","",40,40.,400.); h_pre_mass_b1b2.GetXaxis().SetTitle("m(j,j)");      
h_pre_dR_b1b2            = ROOT.TH1F("h_pre_dR_b1b2","",50,0.,6.);      h_pre_dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)");           
# Mix Gen level
h_gen_pre_met_pt         = ROOT.TH1F("h_gen_pre_met_pt","",50,10.,400.);    h_gen_pre_met_pt.GetXaxis().SetTitle("MET [GeV]");                
h_gen_pre_mass_trans     = ROOT.TH1F("h_gen_pre_mass_trans","",50,0.,250.); h_gen_pre_mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]");     
h_gen_pre_dR_l1l2b1b2    = ROOT.TH1F("h_gen_pre_dR_l1l2b1b2","",50,0.,6.);  h_gen_pre_dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)");     
h_gen_pre_dR_minbl       = ROOT.TH1F("h_gen_pre_dR_minbl","",50,0.,6.);     h_gen_pre_dR_minbl.GetXaxis().SetTitle("min #Delta R(lj)");     
h_gen_pre_dphi_llmet     = ROOT.TH1F("h_gen_pre_dphi_llmet","",50,-4.,4.);  h_gen_pre_dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)");
# Mix Reco
h_pre_met_pt             = ROOT.TH1F("h_pre_met_pt","",50,10.,400.);    h_pre_met_pt.GetXaxis().SetTitle("MET [GeV]");                
h_pre_mass_trans         = ROOT.TH1F("h_pre_mass_trans","",50,0.,250.); h_pre_mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]");     
h_pre_dR_l1l2b1b2        = ROOT.TH1F("h_pre_dR_l1l2b1b2","",50,0.,6.);  h_pre_dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)");     
h_pre_dR_minbl           = ROOT.TH1F("h_pre_dR_minbl","",50,0.,6.);     h_pre_dR_minbl.GetXaxis().SetTitle("min #Delta R(lj)");     
h_pre_dphi_llmet         = ROOT.TH1F("h_pre_dphi_llmet","",50,-4.,4.);  h_pre_dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)");
# CLEANING CUT
# Muons Gen
h_gen_cc_MU1_pt          = ROOT.TH1F("h_gen_cc_MU1_pt","",50,10.,350.);    h_gen_cc_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]");    
h_gen_cc_MU2_pt          = ROOT.TH1F("h_gen_cc_MU2_pt","",50,10.,350.);    h_gen_cc_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); 
h_gen_cc_MU1_eta         = ROOT.TH1F("h_gen_cc_MU1_eta","",50,-3.,3.);     h_gen_cc_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta");          
h_gen_cc_MU2_eta         = ROOT.TH1F("h_gen_cc_MU2_eta","",50,-3.,3.);     h_gen_cc_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta");       
h_gen_cc_mass_l1l2       = ROOT.TH1F("h_gen_cc_mass_l1l2","",40,20.,150.); h_gen_cc_mass_l1l2.GetXaxis().SetTitle("m(l,l)");             
h_gen_cc_dR_l1l2         = ROOT.TH1F("h_gen_cc_dR_l1l2","",50,0.,5.);      h_gen_cc_dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)");           
# gen_B-jets Gen
h_gen_cc_J1_pt           = ROOT.TH1F("h_gen_cc_J1_pt","",50,20.,350.);     h_gen_cc_J1_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]");     
h_gen_cc_J2_pt           = ROOT.TH1F("h_gen_cc_J2_pt","",50,20.,350.);     h_gen_cc_J2_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]");  
h_gen_cc_J1_eta          = ROOT.TH1F("h_gen_cc_J1_eta","",50,-3.,3.);      h_gen_cc_J1_eta.GetXaxis().SetTitle("Lead. jet #eta");           
h_gen_cc_J2_eta          = ROOT.TH1F("h_gen_cc_J2_eta","",50,-3.,3.);      h_gen_cc_J2_eta.GetXaxis().SetTitle("Sublead. jet #eta");        
h_gen_cc_mass_b1b2       = ROOT.TH1F("h_gen_cc_mass_b1b2","",40,40.,400.); h_gen_cc_mass_b1b2.GetXaxis().SetTitle("m(j,j)");              
h_gen_cc_dR_b1b2         = ROOT.TH1F("h_gen_cc_dR_b1b2","",50,0.,6.);      h_gen_cc_dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)");           
# gen_Mix Gen
h_gen_cc_met_pt          = ROOT.TH1F("h_gen_cc_met_pt","",50,10.,400.);    h_gen_cc_met_pt.GetXaxis().SetTitle("MET [GeV]");                
h_gen_cc_mass_trans      = ROOT.TH1F("h_gen_cc_mass_trans","",50,0.,250.); h_gen_cc_mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]");   
h_gen_cc_dR_l1l2b1b2     = ROOT.TH1F("h_gen_cc_dR_l1l2b1b2","",50,0.,6.);  h_gen_cc_dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)");
h_gen_cc_dR_minbl        = ROOT.TH1F("h_gen_cc_dR_minbl","",50,0.,6.);     h_gen_cc_dR_l1l2b1b2.GetXaxis().SetTitle("min #Delta R(lj)");
h_gen_cc_dphi_llmet      = ROOT.TH1F("h_gen_cc_dphi_llmet","",50,-4.,4.);  h_gen_cc_dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)");

# Muons Reco
h_cc_MU1_pt              = ROOT.TH1F("h_cc_MU1_pt","",50,10.,350.);    h_cc_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]");    
h_cc_MU2_pt              = ROOT.TH1F("h_cc_MU2_pt","",50,10.,350.);    h_cc_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); 
h_cc_MU1_eta             = ROOT.TH1F("h_cc_MU1_eta","",50,-3.,3.);     h_cc_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta");          
h_cc_MU2_eta             = ROOT.TH1F("h_cc_MU2_eta","",50,-3.,3.);     h_cc_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta");       
h_cc_mass_l1l2           = ROOT.TH1F("h_cc_mass_l1l2","",40,20.,150.); h_cc_mass_l1l2.GetXaxis().SetTitle("m(l,l)");             
h_cc_dR_l1l2             = ROOT.TH1F("h_cc_dR_l1l2","",50,0.,5.);      h_cc_dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)");           
# B-jets Reco
h_cc_J1_pt               = ROOT.TH1F("h_cc_J1_pt","",50,20.,350.);     h_cc_J1_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]");     
h_cc_J2_pt               = ROOT.TH1F("h_cc_J2_pt","",50,20.,350.);     h_cc_J2_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]");  
h_cc_J1_eta              = ROOT.TH1F("h_cc_J1_eta","",50,-3.,3.);      h_cc_J1_eta.GetXaxis().SetTitle("Lead. jet #eta");           
h_cc_J2_eta              = ROOT.TH1F("h_cc_J2_eta","",50,-3.,3.);      h_cc_J2_eta.GetXaxis().SetTitle("Sublead. jet #eta");        
h_cc_mass_b1b2           = ROOT.TH1F("h_cc_mass_b1b2","",40,40.,400.); h_cc_mass_b1b2.GetXaxis().SetTitle("m(j,j)");              
h_cc_dR_b1b2             = ROOT.TH1F("h_cc_dR_b1b2","",50,0.,6.);      h_cc_dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)");           
# Mix Reco
h_cc_met_pt              = ROOT.TH1F("h_cc_met_pt","",50,10.,400.);    h_cc_met_pt.GetXaxis().SetTitle("MET [GeV]");                
h_cc_mass_trans          = ROOT.TH1F("h_cc_mass_trans","",50,0.,250.); h_cc_mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]");   
h_cc_dR_l1l2b1b2         = ROOT.TH1F("h_cc_dR_l1l2b1b2","",50,0.,6.);  h_cc_dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)");
h_cc_dR_minbl            = ROOT.TH1F("h_cc_dR_minbl","",50,0.,6.);  h_cc_dR_l1l2b1b2.GetXaxis().SetTitle("min #Delta R(lj)");
h_cc_dphi_llmet          = ROOT.TH1F("h_cc_dphi_llmet","",50,-4.,4.);  h_cc_dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)");

h_h2mass_weight_gen      = ROOT.TH1F("h_h2mass_weight_gen","",100,200.,1200.);  h_h2mass_weight_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight1_gen     = ROOT.TH1F("h_h2mass_weight1_gen","",100,200.,1200.);  h_h2mass_weight1_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_reco     = ROOT.TH1F("h_h2mass_weight_reco","",100,200.,1200.);  h_h2mass_weight_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight1_reco    = ROOT.TH1F("h_h2mass_weight1_reco","",100,200.,1200.);  h_h2mass_weight1_reco.GetXaxis().SetTitle("HME reco mass [GeV]");

nEv = 0
for ev in TCha:
  if nEv < nStart or (nEnd > 0 and nEv >= nEnd):
      #print "nEv ",nEv," nstart ",nStart," nEnd ",nEnd
      nEv = nEv +1
      continue
  if (doTest and nEv%10 == 0 ):
      print "nEv ",nEv
  elif (nEv%10000 == 0):
      print "nEv ",nEv
  
  if (doTest and nEv>=10000):
      break
  """
  if ev.htoZZ:
      print "ZZBB event "
      if not ev.findAllGenParticles:
        print "not Find all Gen particles"
  elif ev.htoWW:
      print "WWBB event"
      continue
  else:
      print "NOTZZBBWW event"

  """
  # CUTS
  met_pt 	      = ev.genmet_neutrinos_pt 
  met_px 	      = ev.genmet_neutrinos_px 
  met_py 	      = ev.genmet_neutrinos_py
  MET_cut             = True#(ev.met_pt>20)
  MuMu_cut            = (ev.muon1_pt>20 and fabs(ev.muon1_eta)<2.4 and ev.muon2_pt>10 and fabs(ev.muon2_eta)<2.4 and ev.mass_l1l2>12)
  B1B2_cut            = (ev.b1jet_pt>20 and fabs(ev.b1jet_eta)<2.4 and ev.b2jet_pt>20 and fabs(ev.b2jet_eta)<2.4)
  preselection        = (MET_cut and MuMu_cut and B1B2_cut)
  mt_cleancut         = (ev.mass_trans>10)
  NoZ_cleancut        = (ev.mass_l1l2<(91-15))
  mJJ_cleancut        = (ev.mass_b1b2>30)
  dR_lljj_cleancut    = (ev.dR_l1l2b1b2>0.2 and ev.dR_l1l2b1b2<4)
  dR_ll_cleancut      = (ev.dR_l1l2<3.8)
  dR_jj_cleancut      = (ev.dR_b1b2<3.8)
  cleaning_cuts       = (mt_cleancut and mJJ_cleancut and dR_lljj_cleancut and dR_ll_cleancut and dR_jj_cleancut)
  if (preselection):
      print "eventid ",nEv," pass pre cut"
  if (cleaning_cuts ):
      print "eventid ",nEv," pass CC cut"
  
  #if (preselection and ev.findAllGenParticles and doHME):
  if (ev.findAllGenParticles and doHME):
      mu1_p4 	      = ROOT.TLorentzVector(ev.mu1_px, ev.mu1_py, ev.mu1_pz, ev.mu1_energy)
      mu2_p4 	      = ROOT.TLorentzVector(ev.mu2_px, ev.mu2_py, ev.mu2_pz, ev.mu2_energy)
      b1jet_p4 	      = ROOT.TLorentzVector(ev.b1_px, ev.b1_py, ev.b1_pz, ev.b1_energy)
      b2jet_p4 	      = ROOT.TLorentzVector(ev.b2_px, ev.b2_py, ev.b2_pz, ev.b2_energy)
      #print "ev_nu1_px ",ev.nu1_px," py ",ev.nu1_py," nu2_px ",ev.nu2_px," py ",ev.nu2_py
      met_vec2        = ROOT.TVector2(ev.nu1_px + ev.nu2_px, ev.nu1_py + ev.nu2_py)
      hme_gen = HeavyMassEstimator()
      hme_gen.setKinematic(mu1_p4, mu2_p4, b1jet_p4, b2jet_p4, met_vec2)
      #hme_gen.setonshellZmasspdf(onshellZmasspdf)
      #hme_gen.setoffshellZmasspdf(offshellZmasspdf)
      #hme_gen.setrecobjetrescalec1pdf(recobjetrescalec1pdfPU40)
      #hme_gen.showKinematic()
      hme_gen.setIterations(args.iterations)
      hme_gen.runHME()
      if hme_gen.hme_h2Mass.GetEntries()>0:
	  print "Gen Level most probable reco mass ",hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin())," entries ",hme_gen.hme_h2Mass.GetEntries()," true mass ",ev.h2tohh_mass
	  hme_gen.hme_h2Mass.SetName("hme_h2Mass_ev%d_genlevel"%nEv)
          hme_gen.hme_h2MassWeight1.SetName("hme_h2MassWeight1_ev%d_genlevel"%nEv)
	  if nEv%100 == 0:
	      hme_gen.hme_h2Mass.Write()
	      hme_gen.hme_h2MassWeight1.Write()
          h_h2mass_weight_gen.Fill(hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin()))
          h_h2mass_weight1_gen.Fill(hme_gen.hme_h2MassWeight1.GetXaxis().GetBinCenter(hme_gen.hme_h2MassWeight1.GetMaximumBin()))
  
  if (preselection and doHME and ev.findAllGenParticles):
      muon1_p4 	  = ROOT.TLorentzVector(ev.muon1_px, ev.muon1_py, ev.muon1_pz, ev.muon1_energy)
      muon2_p4 	  = ROOT.TLorentzVector(ev.muon2_px, ev.muon2_py, ev.muon2_pz, ev.muon2_energy)
      b1jet_p4 	  = ROOT.TLorentzVector(ev.b1jet_px, ev.b1jet_py, ev.b1jet_pz, ev.b1jet_energy)
      b2jet_p4 	  = ROOT.TLorentzVector(ev.b2jet_px, ev.b2jet_py, ev.b2jet_pz, ev.b2jet_energy)
      met_vec2    = ROOT.TVector2(met_px, met_py)
      hme = HeavyMassEstimator()
      hme.setKinematic(muon1_p4, muon2_p4, b1jet_p4, b2jet_p4, met_vec2)
      #hme.setonshellZmasspdf(onshellZmasspdf)
      #hme.setoffshellZmasspdf(offshellZmasspdf)
      #hme.setrecobjetrescalec1pdf(recobjetrescalec1pdfPU40)
      #hme.showKinematic()
      hme.setIterations(args.iterations)
      hme.runHME()
      if hme.hme_h2Mass.GetEntries()>0:
	  print "Reco Level most probable reco mass ",hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())," entries ",hme.hme_h2Mass.GetEntries()," true mass ",ev.h2tohh_mass
	  hme.hme_h2Mass.SetName("hme_h2Mass_ev%d_recolevel"%nEv)
          hme.hme_h2MassWeight1.SetName("hme_h2MassWeight1_ev%d_recolevel"%nEv)
	  if nEv%100 == 0:
	      hme.hme_h2Mass.Write()
	      hme.hme_h2MassWeight1.Write()
          h_h2mass_weight_reco.Fill(hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin()))
          h_h2mass_weight1_reco.Fill(hme.hme_h2MassWeight1.GetXaxis().GetBinCenter(hme.hme_h2MassWeight1.GetMaximumBin()))
  # WEIGHTS
  weight = 1.
  if( whichSample!="Data" ):
    # Cross section (We only weights to 1fb-1, the real normalization for xsec will be placed into FinalPlotter.py)
    xSec = ev.XsecBr
    if(nEv==0): h_pre_XsecBr.SetBinContent(1, xSec)
    Lumi    = 1. * 1000 # Convert from 1fb-1 to 1pb-1
    weight = weight * (Lumi*xSec)/nTOT_prehlt
    # Scale factors
    weight = weight * ev.muon1_pogSF * ev.muon2_pogSF
  # Minimal Selection
  if( preselection ):
    # SF
    h_pre_muon1_triggerSF.Fill(abs(ev.muon1_eta), ev.muon1_pt, ev.muon1_triggerSF)
    h_pre_muon1_isoSF.Fill(abs(ev.muon1_eta), ev.muon1_pt, ev.muon1_isoSF)
    h_pre_muon1_idSF.Fill(abs(ev.muon1_eta), ev.muon1_pt, ev.muon1_idSF)
    h_pre_muon1_trackingSF.Fill(ev.muon1_eta, ev.muon1_trackingSF)
    h_pre_muon1_SF_bg1.Fill(ev.muon1_eta)
    h_pre_muon1_SF_bg2.Fill(abs(ev.muon1_eta), ev.muon1_pt)
    h_pre_muon2_triggerSF.Fill(abs(ev.muon2_eta), ev.muon2_pt, ev.muon2_triggerSF)
    h_pre_muon2_isoSF.Fill(abs(ev.muon2_eta), ev.muon2_pt, ev.muon2_isoSF)
    h_pre_muon2_idSF.Fill(abs(ev.muon2_eta), ev.muon2_pt, ev.muon2_idSF)
    h_pre_muon2_trackingSF.Fill(ev.muon2_eta, ev.muon2_trackingSF)
    h_pre_muon2_SF_bg1.Fill(ev.muon2_eta)
    h_pre_muon2_SF_bg2.Fill(abs(ev.muon2_eta), ev.muon2_pt)
    # Regression Variables
    h_pre_numOfVertices.Fill( ev.numOfVertices, weight )   
    h_pre_b1jet_mt.Fill( ev.b1jet_mt, weight )
    h_pre_b1jet_leadTrackPt.Fill( ev.b1jet_leadTrackPt, weight )
    h_pre_b1jet_leptonPtRel.Fill( ev.b1jet_leptonPtRel, weight )
    h_pre_b1jet_leptonPt.Fill( ev.b1jet_leptonPt, weight )
    h_pre_b1jet_leptonDeltaR.Fill( ev.b1jet_leptonDeltaR, weight )
    h_pre_b1jet_neHEF.Fill( ev.b1jet_neHEF, weight )
    h_pre_b1jet_neEmEF.Fill( ev.b1jet_neEmEF, weight )
    h_pre_b1jet_vtxNtracks.Fill( ev.b1jet_vtxNtracks, weight )
    h_pre_b1jet_vtxPt.Fill( ev.b1jet_vtxPt, weight )
    h_pre_b1jet_vtxMass.Fill( ev.b1jet_vtxMass, weight )
    h_pre_b1jet_vtx3DSig.Fill( ev.b1jet_vtx3DSig, weight )
    h_pre_b1jet_vtx3DVal.Fill( ev.b1jet_vtx3DVal, weight )
    # Kinematic Variables
    ####Gen Level
    if (ev.findAllGenParticles):
	h_gen_pre_MU1_pt.Fill( ev.mu1_pt, weight )
	h_gen_pre_MU2_pt.Fill( ev.mu2_pt, weight )
	h_gen_pre_MU2_eta.Fill( ev.mu2_eta, weight )
	h_gen_pre_MU1_eta.Fill( ev.mu1_eta, weight )
	h_gen_pre_mass_l1l2.Fill( ev.mass_genl1l2, weight )
	h_gen_pre_dR_l1l2.Fill( ev.dR_genl1l2, weight )
	h_gen_pre_J1_pt.Fill( ev.b1_pt, weight )
	h_gen_pre_J2_pt.Fill( ev.b2_pt, weight )
	h_gen_pre_J1_eta.Fill( ev.b1_eta, weight )
	h_gen_pre_J2_eta.Fill( ev.b2_eta, weight )
	h_gen_pre_mass_b1b2.Fill( ev.mass_genb1b2, weight )
	h_gen_pre_dR_b1b2.Fill( ev.dR_genb1b2, weight )
	h_gen_pre_met_pt.Fill( ev.genmet_pt, weight )
	h_gen_pre_mass_trans.Fill( ev.mass_gentrans, weight )
	h_gen_pre_dR_l1l2b1b2.Fill( ev.dR_genl1l2b1b2, weight )
	h_gen_pre_dR_minbl.Fill( ev.dR_genminbl, weight )
	h_gen_pre_dphi_llmet.Fill( ev.dphi_genllmet, weight )
    #### RECO level
    h_pre_MU1_pt.Fill( ev.muon1_pt, weight )
    h_pre_MU2_pt.Fill( ev.muon2_pt, weight )
    h_pre_MU2_eta.Fill( ev.muon2_eta, weight )
    h_pre_MU1_eta.Fill( ev.muon1_eta, weight )
    h_pre_mass_l1l2.Fill( ev.mass_l1l2, weight )
    h_pre_dR_l1l2.Fill( ev.dR_l1l2, weight )
    h_pre_J1_pt.Fill( ev.b1jet_pt, weight )
    h_pre_J2_pt.Fill( ev.b2jet_pt, weight )
    h_pre_J1_eta.Fill( ev.b1jet_eta, weight )
    h_pre_J2_eta.Fill( ev.b2jet_eta, weight )
    h_pre_mass_b1b2.Fill( ev.mass_b1b2, weight )
    h_pre_dR_b1b2.Fill( ev.dR_b1b2, weight )
    h_pre_met_pt.Fill( ev.met_pt, weight )
    h_pre_mass_trans.Fill( ev.mass_trans, weight )
    h_pre_dR_l1l2b1b2.Fill( ev.dR_l1l2b1b2, weight )
    h_pre_dR_minbl.Fill( ev.dR_minbl, weight )
    h_pre_dphi_llmet.Fill( ev.dphi_llmet, weight )
    if( cleaning_cuts ):
      # Kinematic Variables
      if ev.findAllGenParticles:
	  h_gen_cc_MU1_pt.Fill( ev.mu1_pt, weight )
	  h_gen_cc_MU2_pt.Fill( ev.mu2_pt, weight )
	  h_gen_cc_MU2_eta.Fill( ev.mu2_eta, weight )
	  h_gen_cc_MU1_eta.Fill( ev.mu1_eta, weight )
	  h_gen_cc_mass_l1l2.Fill( ev.mass_genl1l2, weight )
	  h_gen_cc_dR_l1l2.Fill( ev.dR_genl1l2, weight )
	  h_gen_cc_J1_pt.Fill( ev.b1_pt, weight )
	  h_gen_cc_J2_pt.Fill( ev.b2_pt, weight )
	  h_gen_cc_J1_eta.Fill( ev.b1_eta, weight )
	  h_gen_cc_J2_eta.Fill( ev.b2_eta, weight )
	  h_gen_cc_mass_b1b2.Fill( ev.mass_genb1b2, weight )
	  h_gen_cc_dR_b1b2.Fill( ev.dR_genb1b2, weight )
	  h_gen_cc_met_pt.Fill( ev.genmet_pt, weight )
	  h_gen_cc_mass_trans.Fill( ev.mass_gentrans, weight )
	  h_gen_cc_dR_l1l2b1b2.Fill( ev.dR_genl1l2b1b2, weight )
	  h_gen_cc_dR_minbl.Fill( ev.dR_genminbl, weight )
	  h_gen_cc_dphi_llmet.Fill( ev.dphi_genllmet, weight )
      ###RECO 
      h_cc_MU1_pt.Fill( ev.muon1_pt, weight )
      h_cc_MU2_pt.Fill( ev.muon2_pt, weight )
      h_cc_MU2_eta.Fill( ev.muon2_eta, weight )
      h_cc_MU1_eta.Fill( ev.muon1_eta, weight )
      h_cc_mass_l1l2.Fill( ev.mass_l1l2, weight )
      h_cc_dR_l1l2.Fill( ev.dR_l1l2, weight )
      h_cc_J1_pt.Fill( ev.b1jet_pt, weight )
      h_cc_J2_pt.Fill( ev.b2jet_pt, weight )
      h_cc_J1_eta.Fill( ev.b1jet_eta, weight )
      h_cc_J2_eta.Fill( ev.b2jet_eta, weight )
      h_cc_mass_b1b2.Fill( ev.mass_b1b2, weight )
      h_cc_dR_b1b2.Fill( ev.dR_b1b2, weight )
      h_cc_met_pt.Fill( ev.met_pt, weight )
      h_cc_mass_trans.Fill( ev.mass_trans, weight )
      h_cc_dR_l1l2b1b2.Fill( ev.dR_l1l2b1b2, weight )
      h_cc_dR_minbl.Fill( ev.dR_minbl, weight )
      h_cc_dphi_llmet.Fill( ev.dphi_llmet, weight )
  nEv = nEv +1

h_pre_muon1_triggerSF.Divide(h_pre_muon1_SF_bg2)
h_pre_muon1_isoSF.Divide(h_pre_muon1_SF_bg2)
h_pre_muon1_idSF.Divide(h_pre_muon1_SF_bg2)
h_pre_muon1_trackingSF.Divide(h_pre_muon1_SF_bg1)
h_pre_muon2_triggerSF.Divide(h_pre_muon2_SF_bg2)
h_pre_muon2_isoSF.Divide(h_pre_muon2_SF_bg2)
h_pre_muon2_idSF.Divide(h_pre_muon2_SF_bg2)
h_pre_muon2_trackingSF.Divide(h_pre_muon2_SF_bg1)
#normalize1D(h_pre_muon1_trackingSF)
# Writing histograms
f.cd()
h_pre_XsecBr.Write()
h_pre_muon1_triggerSF.Write()
h_pre_muon1_isoSF.Write()
h_pre_muon1_idSF.Write()
h_pre_muon1_trackingSF.Write()
h_pre_muon2_triggerSF.Write()
h_pre_muon2_isoSF.Write()
h_pre_muon2_idSF.Write()
h_pre_muon2_trackingSF.Write()
h_pre_Nev_preHLT.Write()
h_pre_Nev_posHLT.Write()
h_pre_numOfVertices.Write()  
h_pre_b1jet_mt.Write()
h_pre_b1jet_leadTrackPt.Write()
h_pre_b1jet_leptonPtRel.Write()
h_pre_b1jet_leptonPt.Write()
h_pre_b1jet_leptonDeltaR.Write()
h_pre_b1jet_neHEF.Write()
h_pre_b1jet_neEmEF.Write()
h_pre_b1jet_vtxNtracks.Write()
h_pre_b1jet_vtxPt.Write()
h_pre_b1jet_vtxMass.Write()
h_pre_b1jet_vtx3DSig.Write()
h_pre_b1jet_vtx3DVal.Write()
h_pre_MU1_pt.Write()
h_pre_MU2_pt.Write()
h_pre_MU1_eta.Write()
h_pre_MU2_eta.Write()
h_pre_mass_l1l2.Write()
h_pre_dR_l1l2.Write()
h_pre_J1_pt.Write()
h_pre_J2_pt.Write()
h_pre_J1_eta.Write()
h_pre_J2_eta.Write()
h_pre_mass_b1b2.Write()
h_pre_dR_b1b2.Write()
h_pre_met_pt.Write()
h_pre_mass_trans.Write()
h_pre_dR_l1l2b1b2.Write()
h_pre_dphi_llmet.Write()
# Cleaning Cuts
h_cc_MU1_pt.Write()
h_cc_MU2_pt.Write()
h_cc_MU2_eta.Write()
h_cc_MU1_eta.Write()
h_cc_mass_l1l2.Write()
h_cc_dR_l1l2.Write()
h_cc_J1_pt.Write()
h_cc_J2_pt.Write()
h_cc_J1_eta.Write()
h_cc_J2_eta.Write()
h_cc_mass_b1b2.Write()
h_cc_dR_b1b2.Write()
h_cc_met_pt.Write()
h_cc_mass_trans.Write()
h_cc_dR_l1l2b1b2.Write()
h_cc_dphi_llmet.Write()
h_h2mass_weight_gen.Write()
h_h2mass_weight1_gen.Write()
h_h2mass_weight_reco.Write()
h_h2mass_weight1_reco.Write()
f.Close()
print "Done."
