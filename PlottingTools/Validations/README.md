# Ntuples Validations

The list of plots to validate the ntuples before next step

## Run genParticleProducer.py to find gen particles from decay chain
Instruction to just run modules in genParticleProducer.py to find gen particles from decay chain
```
cd CMSSW_X_X/src
cd HhhAnalysis/python/NanoAOD
python postproc_genpart.py
```
Make sure that the input file is replaced by the one you want to process
### Descriptions of gen information in Ntuples

The follow is the Gen particle information filled in Ntuples by running genHH module in genParticleProducer.py


| branch name |  data type | description |
| ------------- |:-------------:|:-----|
| nGenX    | int  | number of gen X  |
| GenX_status | int | statu code of gen X. if more than one X is found, it is a vector. same for similar branches |
| GenX_phi    | float | global phi of gen X |
| GenX_eta    | float | global eta of gen X |
| GenX_mass    | float | mass of gen X |
| GenX_idx    | int | index of gen X in genParticles |
| GenX_pt    | float | transverse momentum of gen X |
| GenX_statusFlags   | int | status flag of gen X |
| GenX_pdgId    | int | pdgId of gen X |
| GenX_charge    | int | charge of gen X |
| nGenHiggs   | int | num of gen Higgs, from X->HH |
| GenHiggs_xx | xx | xxx |
| nGenW1FromHiggs   | int | num of gen W1, from X->HH->bbW1W2  |
| GenW1FromHiggs_xx | xx | xxx |
| nGenW2FromHiggs   | int | num of gen W2, from X->HH->bbW1W2. in DL channel W2->lv and in SL, W2->qq  |
| GenW2FromHiggs_xx | xx | xxx |
| nGenLepFromW1FromHiggs   | int | num of gen lepton (ele/muon), from W1->lv with W1 from H->WW  |
| GenLepFromW1FromHiggs_xx | xx | xxx |
| nGenNuFromW1FromHiggs   | int | num of gen neutrino (ele/muon neutrino), from  W1->lv with W1 from H->WW  |
| GenNuFromW1FromHiggs_xx | xx | xxx |
| nGenTauNuFromW1FromHiggs   | int | num of gen tau neutrino (tau neutrino), from W1->Tau+v with W1 from H->WW  |
| GenTauNuFromW1FromHiggs_xx | xx | xxx |
| nGenLepFromTauFromW1FromHiggs   | int | num of gen lepton (ele/muon) from tau decay, from  W1->lv->Tau+v->l+v+v+v with with W1 from H->WW   |
| GenLepFromTauFromW1FromHiggs_xx | xx | xxx |
| nGenNuFromTauFromW1FromHiggs   | int | num of gen (ele/muon/tau neutrino), from W1->Tau+v->l+v+v+v with W1 from H->WW  |
| GenNuFromTauFromW1FromHiggs_xx | xx | xxx |
| nGenLepFromW2FromHiggs   | int | num of gen lepton (ele/muon), W2->lv with W2 from H->WW |
| GenLepFromW2FromHiggs_xx | xx | xxx |
| nGenNuFromW2FromHiggs   | int | num of gen neutrino (ele/muon neutrino), from W2->lv with W2 from H->WW |
| GenNuFromW2FromHiggs_xx | xx | xxx |
| nGenTauNuFromW2FromHiggs   | int | num of gen tau neutrino (tau neutrino), from W2->Tau+v  with W2 from H->WW |
| GenTauNuFromW2FromHiggs_xx | xx | xxx |
| nGenLepFromTauFromW2FromHiggs   | int | num of gen lepton (ele/muon) from tau decay, from W2->lv->Tau+v->l+v+v+v with W2 from H->WW  |
| GenLepFromTauFromW2FromHiggs_xx | xx | xxx |
| nGenNuFromTauFromW2FromHiggs   | int | num of gen (ele/muon/tau neutrino), from W2->Tau+v->l+v+v+v with W2 from H->WW |
| GenNuFromTauFromW2FromHiggs_xx | xx | xxx |
| nGenQuarkFromW2FromHiggs   | int | num of gen quark with abs(pdgid)=1,2,3,4,5, from W2->qq with W2 from H->WW |
| GenQuarkFromTauFromW2FromHiggs_xx | xx | xxx |
| nGenBQuarkFromHiggs   | int | num of gen b and bbar quarks, from X->HH->bbW1W2  |
| GenBQuarkFromHiggs_xx | xx | xxx |


The follow is the Gen particle information filled in Ntuples by running genTTbar module in genParticleProducer.py


| branch name |  data type | description |
| ------------- |:-------------:|:-----|
| nGenTop    | int  | number of gen top/antitop  |
| GenTop_xx | xx | xxx |
| nGenBQuarkFromTop    | int  | number of gen b and bbar from top decay  |
| GenBQuarkFromTop_xx | xx | xxx |
| nGenWFromTop    | int  | number of gen W from top decay  |
| GenWFromTop_xx | xx | xxx |
| nGenLepFromWFromTop    | int  | number of gen lepton(ele/muon) from W decay with W from t->Wb |
| GenLepFromWFromTop_xx | xx | xxx |
| nGenNuFromWFromTop    | int  | number of gen leptonic neutrino from W decay with W from t->Wb |
| GenNuFromWFromTop_xx | xx | xxx |
| nGenLepFromTauFromTop    | int  | number of gen lepton(ele/muon) from Tau decay with Tau from t->Wb->Tau+nu+b |
| GenLepFromTauFromTop_xx | xx | xxx |
| nGenQaurkFromWFromTop    | int  | number of gen quark with abs(pdgId)=1,2,3,4,5 from W->qq with W from t->Wb |
| GenQuarkFromWFromTop_xx | xx | xxx |

Attentions to the MC samples:
 - for SL and DL signal samples, W can decay into ele/muon and tau, with one third W going to tau and tau decaying into ele/muon+neutrinos
 - for DL signal smaples, it contains that H->ZZ with one Z->ll and other Z->nu+nu
 
Here is the Higgs decay branch ratio that matters in this analysis, with Higgs mass= 125 GeV. The table is taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR

| Higgs decay channel | branch ratio |
| ------------- |:-------------:|
| H->bb  | 5.792E-01 |
| H->ww  | 2.170E-01 |
| H->zz  | 2.667E-02 |

## Generator level, for Signal and Major backgrounds

### single lepton (SL) channel, resonant signal 
- distribution of gen particle number in one event: X,Higgs,W/b quark from Higgs decay, lepton from W, quark from W decay
  - nGenX, nGenHiggs,nGenWFromHiggs, nGenBQuarkFromHiggs, nGenLepFromW1FromHiggs, nGenQuarkFromW2FromHiggs etc
- eta, pt, mass distributions of gen particles in X->HH->bbWW->bblvqq
- phi and pt, px, py distributions of genMET
- 2D distributions of X mass and X mass reconstructed from HH, Higgs mass and Higgs mass recontructed from bb and WW, W mass and W mass reconstructed from lv and qq,  Higgs mass and Higgs mass reconstructed from lvqq, X mass and X mass reconstructed from bblvqq
- 2D distribution of W mass with one heavier W as x-axis and one lighter W as y-axis

### double lepton (DL) channel, resonant signal
- distribution of gen particle number in one event: X,Higgs,W/b quark from Higgs decay, lepton from W1/W2
  - nGenX, nGenHiggs,nGenWFromHiggs, nGenBQuarkFromHiggs, nGenLepFromW1FromHiggs, nGenLepFromW2FromHiggs etc
- eta, pt mass distributions of gen particles in X->HH->bbWW->bblvlv
- phi and pt, px, py distributions of genMET
- 2D distributions of X mass and X mass reconstructed from HH, Higgs mass and Higgs mass recontructed from bb and WW, W mass and W mass reconstructed from lv,  Higgs mass and Higgs mass reconstructed from lvlv, X mass and X mass reconstructed from bblvlv
- 2D distribution of W mass with one heavier W as x-axis and one lighter W as y-axis

### single lepton (SL) channel, top pair production 
- distribution of gen particle number in one event: top,W/b quark from top decay, lepton from W, quark from W decay
  - nGenTop, nGenWFromTop, nGenBQuarkFromTop, nGenLepFromWFromTop etc
- eta, pt, mass distributions of gen particles in tt->bWbW->blvbqq
- invariant mass of top pair
- phi and pt, px, py distributions of genMET
- 2D distributions of top mass and top mass reconstructed from bW, W mass and W mass reconstructed from lv and qq,  top mass and top mass constructed from blv or bqq
- 2D distribution of W mass with one heavier W as x-axis and one lighter W as y-axis

### double lepton (DL) channel, top pair production 
- distribution of gen particle number in one event: top,W/b quark from top decay, lepton from W, quark from W decay
  - nGenTop, nGenWFromTop, nGenBQuarkFromTop, nGenLepFromWFromTop etc
- eta, pt, mass distributions of gen particles in tt->bWbW->blvblv
- invariant mass of top pair
- phi and pt, px, py distributions of genMET
- 2D distributions of top mass and top mass reconstructed from bW, W mass and W mass reconstructed from lv,  top mass and top mass constructed from blv
- 2D distribution of W mass with one heavier W as x-axis and one lighter W as y-axis

Besides, here is the list of high-level kinematics to calculate:
- invariant mass, pT of dilepton for DL channel by adding up the Lorentz vectors of two leptons 
- invariant mass, pT of dilightquarks for SL channel by adding up the Lorentz vectors of two quarks from W 
- pT of dibquarks by adding up the  Lorentz vectors of two b quarks from Higgs 
- delta R between two leptons for DL, between two quarks from W for SL and between b quarks from Higgs
- minimum delta R between lepton and quark
- delta phi between two leptons (lepton) system and two bquarks system and delta phi between two leptons (lepton) system and genMET for DL(SL) channel 
- transverse mass (MT) for DL channel


see examples to calculate above variables here: https://github.com/tahuang1991/HhhAnalysis/blob/CMSSW_10_2_0/CutFlowAnalyzer/plugins/DiHiggsWWBBSLAnalyzer.cc#L1580

## Genjets
 Match Genjets to quarks in above particle decay chains by deltaR method
 - eta and pt distributions of all genjets that matched to each quarks, including light quarks and b quarks
 - delta R of gen quark and gen jet distributions
 - W Mass reconstructed gen jets matched light quarks from W and Higgs mass reconstructed from gen jets matched to b quarks from Higgs decay. Make it 1D when plotting it alone and 2D when comparing with gen level true particle mass
 - bquark_pt/bgenjet_pt for leading genjet and sub leading genjet that are matched to b quarks, with deltaR(bquark, bgenjet) < 0.3
 
 The calculate the high-level variables with quarks replaced by matched genjets
 
## Reco level, both MC(signals + backgrounds) and data
 here is the list of high-level kinematics with reco level objects to calculate (If AK8jet is selected, then do it for subjet):
- eta, pt distributions of selected reco objects
- MET significance and MET resolution on px, py by doing (recomet_px-genmet_px) and (recomet_py-genmet_py)
- lepton MVA score
- b-tagging discriminator distribution of bjets
- bquark_pt/bjet_pt for leading bjet and sub leading bjet that are matched to b quarks, with deltaR(bquark, bjet) < 0.3
- invariant mass, pT of dilepton for DL channel by adding up the Lorentz vectors of two leptons 
- invariant mass, pT of dijets for SL channel by adding up the Lorentz vectors of two jets from W
- invariant mass, pT of dijets by adding up the  Lorentz vectors of two b jets from Higgs 
- delta R between two leptons for DL, between two quarks from W for SL and between b jets from Higgs
- minimum delta R between lepton and jet
- delta phi between two leptons (lepton) system and two bjets system and delta phi between two leptons (lepton) system and MET for DL(SL) channel 
- transverse mass (MT) for DL channel
- Stransverse Mass (MT2) for DL channel: https://www.hep.phy.cam.ac.uk/~lester/mt2/
- HME for both channels 



## Plotting Style 

Plots should be made  with title, label, stat box and text with explanation etc
 - Text: sample type (Signal with M=? or TTbar ), run condition (Run2016/2017/2018) 
 - X-Title: variable name 
 - Stat box: at least with total number of events, underflow and overflow
 - the plots to validate ntuples at reco level should try to use standard CMS style
   - CMS Publish guide: https://twiki.cern.ch/twiki/bin/view/CMS/Internal/PubGuidelines
   - CMS figure guide: https://twiki.cern.ch/twiki/bin/view/CMS/Internal/FigGuidelines
   
   To use the ROOT .C macros in pyROOT, the simplest way is:
```
from ROOT import gROOT

# Set TDR styles
gROOT.LoadMacro("tdrstyle.C")
gROOT.ProcessLine("setTDRStyle();")

# Add CMS text
gROOT.LoadMacro("CMS_lumi.C")
gROOT.ProcessLine("CMS_lumi(c1);")
```
 
