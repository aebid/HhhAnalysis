# Ntuples Validations

The list of plots to validate the ntuples before next step

## Generator level, for Signal and Major backgrounds

### single lepton (SL) channel, resonance signal 
- eta, pt, mass distributions of gen particles in X->HH->bbWW->bblvqq
- phi and pt, px, py distributions of genMET
- 2D distributions of X mass and X mass reconstructed from HH, Higgs mass and Higgs mass recontructed from bb and WW, W mass and W mass reconstructed from lv and qq,  Higgs mass and Higgs mass reconstructed from lvqq, X mass and X mass reconstructed from bblvqq
- 2D distribution of W mass with one heavier W as x-axis and one lighter W as y-axis

### double lepton (DL) channel, resonance signal
- eta, pt mass distributions of gen particles in X->HH->bbWW->bblvlv
- phi and pt, px, py distributions of genMET
- 2D distributions of X mass and X mass reconstructed from HH, Higgs mass and Higgs mass recontructed from bb and WW, W mass and W mass reconstructed from lv,  Higgs mass and Higgs mass reconstructed from lvlv, X mass and X mass reconstructed from bblvlv
- 2D distribution of W mass with one heavier W as x-axis and one lighter W as y-axis

### single lepton (SL) channel, top pair production 
- eta, pt, mass distributions of gen particles in tt->bWbw->blvbqq
- invariant mass of top pair
- phi and pt, px, py distributions of genMET
- 2D distributions of top mass and top mass reconstructed from bW, W mass and W mass reconstructed from lv and qq,  top mass and top mass constructed from blv or bqq
- 2D distribution of W mass with one heavier W as x-axis and one lighter W as y-axis

### double lepton (DL) channel, top pair production 
- eta, pt, mass distributions of gen particles in tt->bWbw->blvblv
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
- transverse mass for DL channel


see examples to calculate above variables here: https://github.com/tahuang1991/HhhAnalysis/blob/CMSSW_10_2_0/CutFlowAnalyzer/plugins/DiHiggsWWBBSLAnalyzer.cc#L1580

## Genjets
 Match Genjets to quarks in above particle decay chains by deltaR method
 - eta and pt distributions of all genjets that matched to each quarks, including light quarks and b quarks
 - deltR of gen quark and gen jet distributions
 - W Mass reconstructed gen jets matched light quarks from W and Higgs mass reconstructed from gen jets matched to b quarks from Higgs decay. Make it 1D when plotting it alone and 2D when comparing with gen level true particle mass
 
 The calculate the high-level variables with quarks replaced by matched genjets
 
## Reco level
 here is the list of high-level kinematics with reco level objects to calculate:
- eta, pt distributions of selected reco objects
- MET significance 
- lepton MVA score
- b-tagging discriminator distribution of bjets
- invariant mass, pT of dilepton for DL channel by adding up the Lorentz vectors of two leptons 
- invariant mass, pT of dijets for SL channel by adding up the Lorentz vectors of two jets from W
- invariant amss, pT of dijets by adding up the  Lorentz vectors of two b jets from Higgs 
- delta R between two leptons for DL, between two quarks from W for SL and between b jets from Higgs
- minimum delta R between lepton and jet
- delta phi between two leptons (lepton) system and two bjets system and delta phi between two leptons (lepton) system and MET for DL(SL) channel 
- transverse mass for DL channel
- Stransverse Mass for DL channel: https://www.hep.phy.cam.ac.uk/~lester/mt2/
- HME for both channels 
