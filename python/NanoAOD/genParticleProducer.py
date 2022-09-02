import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import math
import logging


sign = lambda x: int(math.copysign(1, x) if x != 0 else 0)


statusFlagsMap = {
  # comments taken from:
  # DataFormats/HepMCCandidate/interface/GenParticle.h
  # PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h
  #
  # nomenclature taken from:
  # PhysicsTools/NanoAOD/python/genparticles_cff.py
  #
  #TODO: use this map in other gen-lvl particle selectors as well
  # GenLepFromTauFromTop -> isDirectPromptTauDecayProduct &&
  #                         isDirectHardProcessTauDecayProduct &&
  #                         isLastCopy &&
  #                         ! isDirectHadronDecayProduct
  # GenLepFromTau -> isDirectTauDecayProduct (or isDirectPromptTauDecayProduct?) &&
  #                  isLastCopy &&
  #                  ! isDirectHadronDecayProduct
  #                  (&& maybe isHardProcessTauDecayProduct?)
  #
  # GenLepFromTop -> isPrompt &&
  #                  isHardProcess &&
  #                  (isLastCopy || isLastCopyBeforeFSR) &&
  #                  ! isDirectHadronDecayProduct
  #
  # Not sure if to choose (isLastCopy or isLastCopyBeforeFSR) or just isFirstCopy:
  # GenWZQuark, GenHiggsDaughters, GenVbosons
  #
  # Have no clue what exactly to require from GenTau
  #
  #
  'isPrompt'                           : 0,  # any decay product NOT coming from hadron, muon or tau decay
  'isDecayedLeptonHadron'              : 1,  # a particle coming from hadron, muon, or tau decay
                                             # (does not include resonance decays like W,Z,Higgs,top,etc)
                                             # equivalent to status 2 in the current HepMC standard
  'isTauDecayProduct'                  : 2,  # a direct or indirect tau decay product
  'isPromptTauDecayProduct'            : 3,  # a direct or indirect decay product of a prompt tau
  'isDirectTauDecayProduct'            : 4,  # a direct tau decay product
  'isDirectPromptTauDecayProduct'      : 5,  # a direct decay product from a prompt tau
  'isDirectHadronDecayProduct'         : 6,  # a direct decay product from a hadron
  'isHardProcess'                      : 7,  # part of the hard process
  'fromHardProcess'                    : 8,  # the direct descendant of a hard process particle of the same pdg id
  'isHardProcessTauDecayProduct'       : 9,  # a direct or indirect decay product of a tau from the hard process
  'isDirectHardProcessTauDecayProduct' : 10, # a direct decay product of a tau from the hard process
  'fromHardProcessBeforeFSR'           : 11, # the direct descendant of a hard process particle of the same pdg id
                                             # for outgoing particles the kinematics are those before QCD or QED FSR
  'isFirstCopy'                        : 12, # the first copy of the particle in the chain with the same pdg id
  'isLastCopy'                         : 13, # the last copy of the particle in the chain with the same pdg id
                                             # (and therefore is more likely, but not guaranteed,
                                             # to carry the final physical momentum)
  'isLastCopyBeforeFSR'                : 14, # the last copy of the particle in the chain with the same pdg id
                                             # before QED or QCD FSR (and therefore is more likely,
                                             # but not guaranteed, to carry the momentum after ISR;
                                             # only really makes sense for outgoing particles
}

def printGenParticles(genparts, genParticles, string=''):
    print("debugging ",string," num of genParts to debug ", len(genparts))
    genParticles = list(genParticles)
    for genp in genparts:
      idx = None
      if hasattr(genp, 'idx'):
        idx = genp.idx
      print("---- id ",genp.pdgId," idx ",idx," pt ",genp.pt," eta ",genp.eta," mass ",genp.mass, " status ",genp.status, " its motheridx ",genp.genPartIdxMother)
      if genp.genPartIdxMother >0:
        genp_mother = genParticles[genp.genPartIdxMother]
        print("\t its  mother id ",genp_mother.pdgId, " mass ",genp_mother.mass, " status ",genp_mother.status)
      if idx > 0 :
        genp_daughters = filter(lambda genPart: genPart.genPartIdxMother == idx, genParticles)
        for genpd in genp_daughters: 
          print("\t its  daughter id ",genpd.pdgId, " mass ",genpd.mass, " status ",genpd.status)

class MassTable:
  def __init__(self):
    self.pdgTable = ROOT.TDatabasePDG()

  def getMass(self, mass, pdgId):
    if mass > 10. or (pdgId == 22 and mass > 1.) or abs(pdgId) == 24 or pdgId == 23:
      return mass
    else:
      genParticleInstance = self.pdgTable.GetParticle(pdgId)
      if not genParticleInstance:
        # Since most of the common low-mass particles are defined in ROOT's PDG table,
        # and that it's more than likely we don't need such generator-level information,
        # we can safely set the masses of such particles to 0 GeV
        logging.debug("Setting the mass to 0 GeV for a particle with PDG id of %d" % pdgId)
        return 0.
      return genParticleInstance.Mass()

  def getCharge(self, pdgId):
    genParticleInstance = self.pdgTable.GetParticle(pdgId)
    if not genParticleInstance:
      # It's more than likely that we don't need to know the charges of generator-level particles
      # that are not defined in ROOT's PDG id table. Therefore, we assign neutral charges to
      # these particles.
      logging.debug("Setting the charge to neutral for a particle with PDG id of %d" % pdgId)
      return 0
    return sign(genParticleInstance.Charge())


class GenPartAux:
  def __init__(self, genPart, idx, massTable):
    self.pt               = genPart.pt
    self.eta              = genPart.eta
    self.phi              = genPart.phi
    self.mass             = massTable.getMass(genPart.mass, genPart.pdgId)
    self.pdgId            = genPart.pdgId
    self.charge           = massTable.getCharge(genPart.pdgId)
    self.status           = genPart.status
    self.statusFlags      = genPart.statusFlags
    self.genPartIdxMother = genPart.genPartIdxMother
    self.idx              = idx

  def __str__(self):
    return("pt = %.3f eta = %.3f phi = %.3f mass = %.3f pdgId = %i charge = %i status = %i " \
           "statusFlags = %i mom = %i idx = %i" % \
      (self.pt, self.eta, self.phi, self.mass, self.pdgId, self.charge, self.status, \
       self.statusFlags, self.genPartIdxMother, self.idx))

  def __repr__(self):
    return self.__str__()

  def chechIf(self, condition):
    assert(condition in statusFlagsMap)
    return (self.statusFlags & (1 << statusFlagsMap[condition]) != 0)


class SelectionOptions:
  SAVE_TAU                      = 0
  SAVE_LEPTONIC_TAU             = 1
  SAVE_HADRONIC_TAU             = 2
  SAVE_LEPTON_FROM_TAU          = 3
  SAVE_LEPTONIC_NU_FROM_TAU     = 4
  SAVE_TAU_NU_FROM_LEPTONIC_TAU = 5
  SAVE_TAU_NU_FROM_HADRONIC_TAU = 6

  SAVE_TOP                               = 9
  SAVE_W_FROM_TOP                        = 10
  SAVE_BQUARK_FROM_TOP                   = 11
  SAVE_LEPTON_FROM_W_FROM_TOP            = 12
  SAVE_LEPTONIC_NU_FROM_W_FROM_TOP       = 13
  SAVE_TAU_FROM_TOP                      = 14
  SAVE_TAU_NU_FROM_TOP                   = 15
  SAVE_LEPTON_FROM_TAU_FROM_TOP          = 16
  SAVE_LEPTON_NU_FROM_TAU_FROM_TOP       = 17
  SAVE_TAU_NU_FROM_LEPTONIC_TAU_FROM_TOP = 18
  SAVE_TAU_NU_FROM_HADRONIC_TAU_FROM_TOP = 19
  SAVE_NU_FROM_TAU_FROM_TOP              = 20
  SAVE_QUARK_FROM_W_FROM_TOP             = 21



class SelectionXToHHTobbWWOptions:
  SAVE_HIG_From_X                            = 0
  SAVE_BQUARK_FROM_HIG                       = 1
  SAVE_W1_FROM_HIG                           = 2
  SAVE_W2_FROM_HIG                           = 3
  SAVE_LEPTON_FROM_W1_FROM_HIG               = 4
  SAVE_LEPTONIC_NU_FROM_W1_FROM_HIG          = 5
  SAVE_LEPTON_FROM_W2_FROM_HIG               = 6
  SAVE_LEPTONIC_NU_FROM_W2_FROM_HIG          = 7
  SAVE_QUARK_FROM_W2_FROM_HIG                = 8
  SAVE_TAU_NU_FROM_W1_FROM_HIG               = 9
  SAVE_LEPTON_FROM_TAU_FROM_W1_FROM_HIG      = 10
  SAVE_NU_FROM_TAU_FROM_W1_FROM_HIG          = 11
  SAVE_TAU_NU_FROM_W2_FROM_HIG               = 12
  SAVE_LEPTON_FROM_TAU_FROM_W2_FROM_HIG      = 13
  SAVE_NU_FROM_TAU_FROM_W2_FROM_HIG          = 14
  SAVE_TAU_FROM_WZ_FROM_HIG                  = 15



def genLeptonSelection(genParticles):
  return filter(lambda genPart: abs(genPart.pdgId) in [11, 13] and genPart.status == 1, genParticles)

def genPromptLeptonSelection(genParticles):
  return filter(
    lambda genLepton:
      genLepton.chechIf('isLastCopy') and
      not genLepton.chechIf('isDirectHadronDecayProduct') and
      (
        genLepton.chechIf('isPrompt') or
        genLepton.chechIf('isDirectPromptTauDecayProduct')
      ),
    genLeptonSelection(genParticles)
  )




def genHiggsSelection(genParticles):
  return filter(
    lambda genPart:
      genPart.pdgId == 25 and \
      (genParticles[genPart.genPartIdxMother].pdgId != 25 if genPart.genPartIdxMother >= 0 else True),
    genParticles
  )

def genRadionGravitonSelection(genParticles): ## X->HH
  XCandidates = []
  pdgIds =[35, 39]## add potential id for X
     
  XCandidates = filter (lambda genPart: abs(genPart.pdgId) in pdgIds, genParticles)
  #printGenParticles(XCandidates, genParticles, "genRadionGravitonSelection before final filter")

  return filter(
    lambda genPart:
      (genParticles[genPart.genPartIdxMother].pdgId != genPart.pdgId if genPart.genPartIdxMother >= 0 else True),
    XCandidates
  )

def genHiggsDaughtersSelection(genParticles):
  return filter(
    lambda genPart:
      genPart.pdgId != 25 and \
      (genParticles[genPart.genPartIdxMother].pdgId == 25 if genPart.genPartIdxMother >= 0 else False),
    genParticles
  )




def genWZquarkSelection(genParticles):
  return filter(
    lambda genPart:
      abs(genPart.pdgId) in [1, 2, 3, 4, 5, 6] and genPart.genPartIdxMother >= 0 and \
      abs(genParticles[genPart.genPartIdxMother].pdgId) in [23, 24],
    genParticles
  )

def genVbosonSelection(genParticles):
  return filter(
    lambda genPart:
      abs(genPart.pdgId) in [23, 24] and genPart.genPartIdxMother >= 0 and \
      genParticles[genPart.genPartIdxMother].pdgId != genPart.pdgId,
    genParticles
  )

def genNuSelection(genParticles):
  return filter(lambda genPart: abs(genPart.pdgId) in [12, 14, 16], genParticles)

def genTopSelection(genParticles, choice, enable_consistency_checks = True):
  genTopMotherMap = {}
  genWMotherMap = {}
  for genPart in genParticles:
    if abs(genPart.pdgId) == 6:
      genTopMotherMap[genPart.idx] = genPart.genPartIdxMother # daughter -> mother
    elif abs(genPart.pdgId) == 24:
      genWMotherMap[genPart.genPartIdxMother] = genPart.idx # mother -> daughter

  genTopCandidates = {}
  for genTopIdx in genTopMotherMap:
    if genTopIdx not in genTopMotherMap.values():
      genTopCandidates[genTopIdx] = []

  if choice == SelectionOptions.SAVE_TOP:
    return map(lambda genTopIdx: genParticles[genTopIdx], genTopCandidates)

  for genPart in genParticles:
    if genPart.genPartIdxMother in genTopCandidates:
      genTopCandidates[genPart.genPartIdxMother].append(genPart.idx)

  # top always decays into W + q, where q = d, s, b with b the most common one
  genBquarkFromTop               = []
  genWFromTop                    = []
  genLepsFromWfromTop            = []
  genNusFromWfromTop             = []
  genTausFromTop                 = []
  genNusTauFromTop               = []
  genLepsFromTauFromTop          = []
  genNuLepFromTauFromTop         = []
  genNuTauFromLeptonicTauFromTop = []
  genNuTauFromHadronicTauFromTop = []
  genQuarkFromWfromTop           = []

  for genTopCandidateIdx, genTopCandidateDaughterIdxs in genTopCandidates.items():
    genTopCandidate = genParticles[genTopCandidateIdx]
    if len(genTopCandidateDaughterIdxs) != 2:
      ##no intermediate state in NanoAOd samples
      raise ValueError("Invalid number of top (%s) decay products (%s): %i" % \
        (genTopCandidate, ', '.join(map(lambda idx: genParticles[idx], genTopCandidateDaughterIdxs)), len(genTopCandidateDaughterIdxs))
      )

    genWsfromTop = [ genParticles[idx] for idx in genTopCandidateDaughterIdxs if genParticles[idx].pdgId == 24 * sign(genTopCandidate.pdgId) ]
    genQsfromTop = [ genParticles[idx] for idx in genTopCandidateDaughterIdxs if sign(genTopCandidate.pdgId) * genParticles[idx].pdgId in [1, 3, 5] ]

    if len(genWsfromTop) != 1:
      raise ValueError("Not exactly 1 W boson found from top (%s) decay: %s" % \
        (genTopCandidate, ', '.join(map(lambda idx: str(genParticles[idx]), genTopCandidateDaughterIdxs)))
      )
    if len(genQsfromTop) != 1:
      raise ValueError("Not exactly 1 quark found from top (%s) decay: %s" % \
        (genTopCandidate, ', '.join(map(lambda idx: str(genParticles[idx]), genTopCandidateDaughterIdxs)))
      )

    if abs(genQsfromTop[0].pdgId) == 5:
      genBquarkFromTop.extend(genQsfromTop)
    if choice == SelectionOptions.SAVE_BQUARK_FROM_TOP:
      continue # skip the next part

    genWfromTop = genWsfromTop[0]
    # now search for leptons and/or taus which are descendant of the W boson
    # however, the W might ,,travel'' until it decays, i.e. its immediate daughter can be a single W boson
    # that's why we have to loop over the generator level particle collection and find out which W is the last one
    # in the top decay chain
    while genWfromTop.idx in genWMotherMap:
      genWfromTop = genParticles[genWMotherMap[genWfromTop.idx]]

    # let's look at W's decay products; possibilities include:
    # 1) leptonic: W -> l vl
    # 2) tauonic:  W -> tau vtau
    # 3) hadronic: W -> q q'
    if choice == SelectionOptions.SAVE_W_FROM_TOP:
      genWFromTop.append(genWfromTop)
      continue # skip the next part


    genWfromTopDaughters = filter(lambda genPart: genPart.genPartIdxMother == genWfromTop.idx, genParticles)
    if len(genWfromTopDaughters) != 2:
      raise ValueError("Invalid number (%i) of W (%s) daughters from top (%s) decay: %s" % \
        (len(genWfromTopDaughters), genWfromTop, genTopCandidate, ', '.join(map(str, genWfromTopDaughters)))
      )

    if any(map(lambda genPart: abs(genPart.pdgId) in [11, 13], genWfromTopDaughters)):
      lepsFromWfromTop = filter(lambda genPart: -sign(genWfromTop.pdgId) * genPart.pdgId in [11, 13], genWfromTopDaughters)
      if len(lepsFromWfromTop) != 1:
        raise ValueError("Inconsistent W (%s) decay products from top (%s) decay: %s" % \
          (genWfromTop, genTopCandidate, ', '.join(map(str, genWfromTopDaughters)))
        )
      genLepFromWfromTop = lepsFromWfromTop[0]
      nusLepFromWfromTop = filter(lambda genPart: genPart.pdgId != sign(genWfromTop.pdgId) * (abs(genLepFromWfromTop.pdgId) + 1), genWfromTopDaughters)
      if len(nusLepFromWfromTop) != 1:
        raise ValueError("Inconsistent W (%s) decay products from top (%s) decay: %s" % \
          (genWfromTop, genTopCandidate, ', '.join(map(str, genWfromTopDaughters)))
        )
      genLepsFromWfromTop.extend(lepsFromWfromTop)
      genNusFromWfromTop.extend(nusLepFromWfromTop)
    elif any(map(lambda genPart: abs(genPart.pdgId) == 15, genWfromTopDaughters)):
      genTausFromWfromTop   = filter(lambda genPart: genPart.pdgId == -sign(genWfromTop.pdgId) * 15, genWfromTopDaughters)
      genNusTauFromWfromTop = filter(lambda genPart: genPart.pdgId ==  sign(genWfromTop.pdgId) * 16, genWfromTopDaughters)

      if len(genTausFromWfromTop) != 1 or len(genNusTauFromWfromTop) != 1:
        raise ValueError("Inconsistent W (%s) tauonic decay products from top (%s) decay: %s" % \
          (genWfromTop, genTopCandidate, ', '.join(map(str, genWfromTopDaughters)))
        )

      genTausFromTop.extend(genTausFromWfromTop)
      genNusTauFromTop.extend(genNusTauFromWfromTop)

      if choice not in [
          SelectionOptions.SAVE_LEPTON_FROM_TAU_FROM_TOP,
          SelectionOptions.SAVE_LEPTON_NU_FROM_TAU_FROM_TOP,
          SelectionOptions.SAVE_TAU_NU_FROM_LEPTONIC_TAU_FROM_TOP,
          SelectionOptions.SAVE_TAU_NU_FROM_HADRONIC_TAU_FROM_TOP,
        ]:
        continue

      genTauFromWfromTop = genTausFromWfromTop[0]
      # let's check tau decay products from t -> Wb, W -> tau vtau
      # map all mother taus to their daughters and pick the one we chose
      genTauMotherMap = {}
      for genPart in genParticles:
        if abs(genPart.pdgId) == 15:
          genTauMotherMap[genPart.genPartIdxMother] = genPart.idx

      genTauFromWfromTop_toDecay = genTauFromWfromTop
      while genTauFromWfromTop_toDecay.idx in genTauMotherMap:
        genTauFromWfromTop_toDecay = genParticles[genTauMotherMap[genTauFromWfromTop_toDecay.idx]]

      # now check tau decay products
      tauFromWfromTopDaughters = filter(lambda genPart: genPart.genPartIdxMother == genTauFromWfromTop_toDecay.idx, genParticles)
      if any(map(lambda genPart: abs(genPart.pdgId) in [11, 13], tauFromWfromTopDaughters)):
        # leptonic tau decay, record lepton, tau neutrino and lepton neutrino
        nusTauFromTauFromWfromTop = filter(lambda genPart: genPart.pdgId == 16 * sign(genTauFromWfromTop_toDecay.pdgId), tauFromWfromTopDaughters)
        if len(nusTauFromTauFromWfromTop) != 1:
          raise ValueError("Not enough tau neutrinos in leptonic tau (%s) decay from W (%s) from top (%s): %s"\
            (genTauFromWfromTop_toDecay, genWfromTop, genTopCandidate, ', '.join(map(str, tauFromWfromTopDaughters)))
          )
        lepsFromTauFromTauFromWfromTop = filter(lambda genPart: abs(genPart.pdgId) in [11, 13], tauFromWfromTopDaughters)
        if len(lepsFromTauFromTauFromWfromTop) != 1:
          raise ValueError("Too many leptons in leptonic tau (%s) decay from W (%s) from top (%s): %s"\
            (genTauFromWfromTop_toDecay, genWfromTop, genTopCandidate, ', '.join(map(str, tauFromWfromTopDaughters)))
          )
        nusLepFromTauFromWfromTop = filter(
          lambda genPart: genPart.pdgId == -sign(lepsFromTauFromTauFromWfromTop[0].pdgId) * (abs(lepsFromTauFromTauFromWfromTop[0].pdgId) + 1),
          tauFromWfromTopDaughters
        )
        if len(nusLepFromTauFromWfromTop) != 1:
          raise ValueError("Not enough lepton neutrinos in leptonic tau (%s) decay from W (%s) from top (%s): %s"\
            (genTauFromWfromTop_toDecay, genWfromTop, genTopCandidate, ', '.join(map(str, tauFromWfromTopDaughters)))
          )
        genLepsFromTauFromTop.extend(lepsFromTauFromTauFromWfromTop)
        genNuLepFromTauFromTop.extend(nusLepFromTauFromWfromTop)
        genNuTauFromLeptonicTauFromTop.extend(nusTauFromTauFromWfromTop)
      else:
        # hadronic tau decay, record the tau neutrinos only
        nusTauFromTauFromWfromTop = filter(lambda genPart: genPart.pdgId == 16 * sign(genTauFromWfromTop_toDecay.pdgId), tauFromWfromTopDaughters)
        if len(nusTauFromTauFromWfromTop) != 1:
          raise ValueError("Invalid hadronic tau (%s) decay products from W (%s) decay in top (%s) decay: %s" % \
            (genTauFromWfromTop_toDecay, genWfromTop, genTopCandidate, ', '.join(map(str, tauFromWfromTopDaughters)))
          )
        genNuTauFromHadronicTauFromTop.extend(nusTauFromTauFromWfromTop)
    elif all(map(lambda genPart: abs(genPart.pdgId) in [1, 2, 3, 4, 5], genWfromTopDaughters)):
      # hadronic case
      if enable_consistency_checks:
        genWfromTopDaughters_pdgIdSorted = list(sorted(genWfromTopDaughters, key = lambda genPart: abs(genPart.pdgId), reverse = True))
        if not ((genWfromTopDaughters_pdgIdSorted[0].pdgId == -5 * sign(genWfromTop.pdgId) and genWfromTopDaughters_pdgIdSorted[1].pdgId ==  4 * sign(genWfromTop.pdgId)) or \
                (genWfromTopDaughters_pdgIdSorted[0].pdgId == -5 * sign(genWfromTop.pdgId) and genWfromTopDaughters_pdgIdSorted[1].pdgId ==  2 * sign(genWfromTop.pdgId)) or \
                (genWfromTopDaughters_pdgIdSorted[0].pdgId ==  4 * sign(genWfromTop.pdgId) and genWfromTopDaughters_pdgIdSorted[1].pdgId == -3 * sign(genWfromTop.pdgId)) or \
                (genWfromTopDaughters_pdgIdSorted[0].pdgId ==  4 * sign(genWfromTop.pdgId) and genWfromTopDaughters_pdgIdSorted[1].pdgId == -1 * sign(genWfromTop.pdgId)) or \
                (genWfromTopDaughters_pdgIdSorted[0].pdgId == -3 * sign(genWfromTop.pdgId) and genWfromTopDaughters_pdgIdSorted[1].pdgId ==  2 * sign(genWfromTop.pdgId)) or \
                (genWfromTopDaughters_pdgIdSorted[0].pdgId ==  2 * sign(genWfromTop.pdgId) and genWfromTopDaughters_pdgIdSorted[1].pdgId == -1 * sign(genWfromTop.pdgId))):
          raise ValueError("Invalid hadronic W (%s) decay products from top (%s): %s" % \
            (genWfromTop, genTopCandidate, ', '.join(map(str, genWfromTopDaughters)))
          )
      genQuarkFromWfromTop.extend(genWfromTopDaughters)
    else:
      raise ValueError("Invalid W (%s) daughters from top (%s) decay: %s" % \
        (genWfromTop, genTopCandidate, ', '.join(map(str, genWfromTopDaughters)))
      )

  if choice == SelectionOptions.SAVE_BQUARK_FROM_TOP:
    return genBquarkFromTop
  if choice == SelectionOptions.SAVE_W_FROM_TOP:
    return genWFromTop
  if choice == SelectionOptions.SAVE_LEPTON_FROM_W_FROM_TOP:
    return genLepsFromWfromTop
  if choice == SelectionOptions.SAVE_LEPTONIC_NU_FROM_W_FROM_TOP:
    return genNusFromWfromTop
  if choice == SelectionOptions.SAVE_TAU_FROM_TOP:
    return genTausFromTop
  if choice == SelectionOptions.SAVE_TAU_NU_FROM_TOP:
    return genNusTauFromTop
  if choice == SelectionOptions.SAVE_TAU_NU_FROM_HADRONIC_TAU_FROM_TOP:
    return genNuTauFromHadronicTauFromTop
  if choice == SelectionOptions.SAVE_LEPTON_FROM_TAU_FROM_TOP:
    return genLepsFromTauFromTop
  if choice == SelectionOptions.SAVE_TAU_NU_FROM_LEPTONIC_TAU_FROM_TOP:
    return genNuTauFromLeptonicTauFromTop
  if choice == SelectionOptions.SAVE_LEPTON_NU_FROM_TAU_FROM_TOP:
    return genNuLepFromTauFromTop
  if choice == SelectionOptions.SAVE_QUARK_FROM_W_FROM_TOP:
    return genQuarkFromWfromTop

  raise ValueError("Invalid selection option: %i" % choice)

def genTauSelection(genParticles, choice, enable_consistency_checks = False):
  genTauMotherMap = {}
  for genPart in genParticles:
    if abs(genPart.pdgId) == 15:
      genTauMotherMap[genPart.idx] = genPart.genPartIdxMother

  genTauCandidates = {}
  for genTauIdx in genTauMotherMap:
    if genTauIdx not in genTauMotherMap.values():
      genTauCandidates[genTauIdx] = { 'daughters' : [], 'isLeptonic' : False }

  for genPart in genParticles:
    if genPart.genPartIdxMother in genTauCandidates:
      genTauCandidates[genPart.genPartIdxMother]['daughters'].append(genPart)
      genTauCandidates[genPart.genPartIdxMother]['isLeptonic'] = genTauCandidates[genPart.genPartIdxMother]['isLeptonic'] or \
                                                                 abs(genPart.pdgId) in [11, 13]

  # assert that the decay products of the leptonic taus are consistent
  if enable_consistency_checks:
    for genTauIdx in genTauCandidates:
      if genTauCandidates[genTauIdx]['isLeptonic']:
        genTauCurrent   = genParticles[genTauIdx]
        genTauDaughters = genTauCandidates[genTauIdx]['daughters']

        if len(genTauDaughters) != 3:
          raise ValueError("Not enough daughters in leptonic tau (%s) decay: %s" % \
            (genTauCurrent, ', '.join(map(str, genTauDaughters)))
          )
        genTauDaughterLep   = None
        genTauDaughterNuLep = None
        genTauDaughterNuTau = None
        for daughter in genTauDaughters:
          if abs(daughter.pdgId) in [11, 13]:
            genTauDaughterLep = daughter
          elif abs(daughter.pdgId) in [12, 14]:
            genTauDaughterNuLep = daughter
          elif abs(daughter.pdgId) == 16:
            genTauDaughterNuTau = daughter

        if genTauDaughterNuLep.pdgId is None:
          raise ValueError("Could not find lepton nu from leptonic tau (%s) decay (daughters: %s)" % \
            (genTauCurrent, ', '.join(map(str, genTauDaughters)))
          )
        if genTauDaughterNuTau.pdgId is None:
          raise ValueError("Could not find tau nu from leptonic tau (%s) decay (daughters: %s)" % \
            (genTauCurrent, ', '.join(map(str, genTauDaughters)))
          )

        if sign(genTauCurrent.pdgId) != sign(genTauDaughterLep.pdgId):
          raise ValueError("Wrong signs in tau (%s) and lepton charges (%s)" % (genTauCurrent, genTauDaughterLep))
        if sign(genTauDaughterLep.pdgId) * (abs(genTauDaughterLep.pdgId) + 1) != -genTauDaughterNuLep.pdgId:
          raise ValueError("Inconsistent pdgIds b/w lepton (%s) and lepton nu (%s) in leptonic tau decay" % \
            (genTauDaughterLep, genTauDaughterNuLep)
          )

  if choice == SelectionOptions.SAVE_TAU:
    return map(lambda genTauIdx: genParticles[genTauIdx], genTauCandidates)
  elif choice in [
      SelectionOptions.SAVE_HADRONIC_TAU,
      SelectionOptions.SAVE_TAU_NU_FROM_HADRONIC_TAU,
    ]:
    genHadronicTauIdxs = filter(lambda genTauIdx: not genTauCandidates[genTauIdx]['isLeptonic'], genTauCandidates)
    if choice == SelectionOptions.SAVE_HADRONIC_TAU:
      return map(lambda genLeptonicTauIdx: genParticles[genLeptonicTauIdx], genHadronicTauIdxs)

    genHadronicTauDaughterArrays = map(
      lambda genLeptonicTauIdx: genTauCandidates[genLeptonicTauIdx]['daughters'], genHadronicTauIdxs
    )
    genNusTauFromHadTaus = []
    for genHadronicTauDaughters in genHadronicTauDaughterArrays:
      for genHadronicTauDaughter in genHadronicTauDaughters:
        if abs(genHadronicTauDaughter.pdgId) == 16:
          genNusTauFromHadTaus.append(genHadronicTauDaughter)
    if choice == SelectionOptions.SAVE_TAU_NU_FROM_HADRONIC_TAU:
      return genNusTauFromHadTaus
  elif choice in [
      SelectionOptions.SAVE_LEPTONIC_TAU,
      SelectionOptions.SAVE_LEPTON_FROM_TAU,
      SelectionOptions.SAVE_LEPTONIC_NU_FROM_TAU,
      SelectionOptions.SAVE_TAU_NU_FROM_LEPTONIC_TAU,
    ]:
    genLeptonicTauIdxs = filter(lambda genTauIdx: genTauCandidates[genTauIdx]['isLeptonic'], genTauCandidates)
    if choice == SelectionOptions.SAVE_LEPTONIC_TAU:
      return map(lambda genLeptonicTauIdx: genParticles[genLeptonicTauIdx], genLeptonicTauIdxs)

    genLeptonicTauDaughterArrays = map(
      lambda genLeptonicTauIdx: genTauCandidates[genLeptonicTauIdx]['daughters'], genLeptonicTauIdxs
    )
    genLeptonsFromTaus, genNusLepFromTaus, genNusTauFromLepTaus = [], [], []
    for genLeptonicTauDaughters in genLeptonicTauDaughterArrays:
      for genLeptonicTauDaughter in genLeptonicTauDaughters:
        if abs(genLeptonicTauDaughter.pdgId) in [11, 13]:
          genLeptonsFromTaus.append(genLeptonicTauDaughter)
        elif abs(genLeptonicTauDaughter.pdgId) in [12, 14]:
          genNusLepFromTaus.append(genLeptonicTauDaughter)
        elif abs(genLeptonicTauDaughter.pdgId) == 16:
          genNusTauFromLepTaus.append(genLeptonicTauDaughter)

    if choice == SelectionOptions.SAVE_LEPTON_FROM_TAU:
      return genLeptonsFromTaus
    if choice == SelectionOptions.SAVE_LEPTONIC_NU_FROM_TAU:
      return genNusLepFromTaus
    if choice == SelectionOptions.SAVE_TAU_NU_FROM_LEPTONIC_TAU:
      return genNusTauFromLepTaus
  else:
    raise ValueError("Choice %i not implemented" % choice)

def findLastDecendantWithSameId(cand, genParticles, allowEmission = False):
  nextgeneration = filter(lambda genPart: genPart.pdgId == cand.pdgId and genPart.genPartIdxMother == cand.idx, genParticles)
  allnextgenerations = filter(lambda genPart: genPart.genPartIdxMother == cand.idx, genParticles)
  if len(nextgeneration) == 0: 
    return cand 
  if len(nextgeneration) == 1:
    if allowEmission or len(allnextgenerations)  == 1:
          return findLastDecendantWithSameId(nextgeneration[0], genParticles, allowEmission)
    elif not(allowEmission) and len(allnextgenerations) > 1:
      return cand  
  else:
      ## example: tau->3tau
      printGenParticles([cand], genParticles, "Cand in findLastDecenDant")
      printGenParticles(nextgeneration, genParticles,"findLastDecendantError!!")
      return cand
      raise ValueError("Invalid number of nextgeneration(%s) of cand (%s) in findLastDecendantWithSameId %i:" % \
	(', '.join(map(lambda genPart: str(genPart), nextgeneration)), str(cand), len(nextgeneration))
      )


def genXToHHTobbWWSelection(genParticles, choice, enable_consistency_checks = True):
  ## resonance case: X->HH,X pdgid = 35 for radion and 39 for graviton
  ##nonresonance case: gg->HH, maybe ignore gg
  genBQuarkFromHIG = []
  genW1FromHIG = None
  genW2FromHIG = None
  genLepFromW1FromHIG        = None
  genNuFromW1FromHIG         = None
  genLepFromW2FromHIG        = None
  genNuFromW2FromHIG         = None
  debugXToHHTobbWW           = False
  #XCandidates  = genRadionGravitonSelection(genParticles)
  XCandidates = []
  pdgIds =[35, 39]## add potential id for X
     
  XCandidates = filter (lambda genPart: abs(genPart.pdgId) in pdgIds, genParticles)
  XCandidatesidx = [X.idx  for X in XCandidates]
  finalXCandidates = filter(lambda X : X.genPartIdxMother not in XCandidatesidx, XCandidates)
  #if nonresonance and len(finalXCandidates) == 2://Dihiggs from GG
  #  Xid1 =  finalXCandidates[0].genPartIdxMother
  #  Xid2 =  finalXCandidates[1].genPartIdxMother
  #  if Xid1 == Xid2: 
  #    finalXCandidates = [genParticles[Xid1]]
  #    XCandidatesidx = [Xid1] 
  #  else:
  #    finalXCandidates = []
  #    XCandidatesidx = []
     
  if len(finalXCandidates) != 1: 
    if debugXToHHTobbWW: 
      print("Warning!! Xcandidates is not one ",len(finalXCandidates))
      printGenParticles(finalXCandidates, genParticles, "genXToHHTobbWWSelection")
    ##avoid printout for non-signal MC
    #return []
  
  ##step1 find X->HH
  genHiggs = filter(lambda genPart: genPart.pdgId == 25, genParticles)
  genHiggsidx = [Higgs.idx for Higgs in genHiggs]
  genHiggsFromX = []
  if len(finalXCandidates) != 1: ## use non-res logic 
    genHiggsFromX = filter(lambda genPart: genPart.genPartIdxMother not in genHiggsidx, genHiggs)
  else:
    genHiggsFromX = filter(lambda genPart: genPart.genPartIdxMother in XCandidatesidx, genHiggs)


  if len(genHiggsFromX) != 2:
    raise ValueError("Invalid number of X->HH (%s) decay products (%s): %i" % \
      (XCandidates[0], ', '.join(map(lambda genPart: str(genPart), genHiggsFromX)), len(genHiggsFromX))
    )

  if choice == SelectionXToHHTobbWWOptions.SAVE_HIG_From_X:
    return genHiggsFromX
  ##step2 find HH->bbWW
  ##possible combination: HH->bbbb; HH->WWWW;  HH->bbbb?
  ##possible HH->bbZZ
  ## no intermediate state in NanoAOD sample ?
  genBQuarkFromHiggs = filter(lambda genPart: abs(genPart.pdgId) == 5 and genPart.genPartIdxMother in genHiggsidx, genParticles)
  genWFromHiggs = filter(lambda genPart: abs(genPart.pdgId) == 24 and genPart.genPartIdxMother in genHiggsidx, genParticles)
  genZFromHiggs = filter(lambda genPart: abs(genPart.pdgId) == 23 and genPart.genPartIdxMother in genHiggsidx, genParticles)
  isHHbbWW = len(genBQuarkFromHiggs) == 2 and len(genWFromHiggs) == 2 

  if choice == SelectionXToHHTobbWWOptions.SAVE_BQUARK_FROM_HIG:
    if len(genBQuarkFromHiggs)  != 2 or genBQuarkFromHiggs[0].genPartIdxMother != genBQuarkFromHiggs[1].genPartIdxMother:
      if debugXToHHTobbWW:
	      print("Invalid number of H->bb decay products (%s): %i" % \
          ( ', '.join(map(lambda genPart: str(genPart), genBQuarkFromHiggs)), len(genBQuarkFromHiggs)))
      return []
    return genBQuarkFromHiggs

  if len(genWFromHiggs)  != 2 or genWFromHiggs[0].genPartIdxMother != genWFromHiggs[1].genPartIdxMother:
    if len(genZFromHiggs) != 2 and debugXToHHTobbWW:
      print("Invalid number of H->bb (%s):%d or H->WW (%s):%d in HH (%s) decay" % \
        (', '.join(map(lambda genPart: str(genPart), genBQuarkFromHiggs)), len(genBQuarkFromHiggs), \
         ', '.join(map(lambda genPart: str(genPart), genWFromHiggs)), len(genWFromHiggs), \
         ', '.join(map(lambda genPart: str(genPart), genHiggsFromX))))
    #HH sample also include HH->ZZbb
    #printGenParticles(genHiggsFromX, genParticles, "debugging H->bbWW") 
    #printGenParticles(genZFromHiggs, genParticles, "debugging H->ZZ") 
    return []
      
  ##step3 find W->lv, (only to lepton)
  def findLepandNufromW(genW):
    genW = findLastDecendantWithSameId(genW, genParticles, True)
    genWDaughters = filter(lambda genPart: genPart.genPartIdxMother == genW.idx, genParticles)
    daughters = {}
    daughters['Lep'] =  []
    daughters['Nu']  =  []
    if len(genWDaughters) != 2:
      raise ValueError("Invalid number (%i) of W (%s) daughters decay: %s" % \
        (len(genWDaughters), genW, ', '.join(map(str, genWDaughters)))
      )

    if any(map(lambda genPart: abs(genPart.pdgId) in [11, 13], genWDaughters)):
      lepsFromW = filter(lambda genPart: -sign(genW.pdgId) * genPart.pdgId in [11, 13], genWDaughters)
      if len(lepsFromW) != 1:
        raise ValueError("Inconsistent W (%s) decay products decay: %s" % \
          (genW, ', '.join(map(str, genWDaughters)))
        )
      genLepFromW = lepsFromW[0]
      nusLepFromW = filter(lambda genPart: genPart.pdgId != sign(genW.pdgId) * (abs(genLepFromW.pdgId) + 1), genWDaughters)
      if len(nusLepFromW) != 1:
        raise ValueError("Inconsistent W (%s) decay products decay: %s" % \
          (genW, ', '.join(map(str, genWDaughters)))
        )
      daughters['Lep'].append(genLepFromW)
      daughters['Nu'].append(nusLepFromW[0])

    return daughters

  ##step3 find W->q1q2, (only to quarks, hadronic case)
  def findQuarksfromW(genW):
    genW = findLastDecendantWithSameId(genW, genParticles, True)
    genWDaughters = filter(lambda genPart: genPart.genPartIdxMother == genW.idx, genParticles)
    if len(genWDaughters) != 2:
      raise ValueError("Invalid number (%i) of W (%s) daughters decay: %s" % \
        (len(genWDaughters), genW, ', '.join(map(str, genWDaughters)))
      )
    daughters = filter(lambda genPart: abs(genPart.pdgId) in [1,2,3,4,5], genWDaughters) 
    if enable_consistency_checks and len(daughters) == 2:
      genWDaughters_pdgIdSorted = list(sorted(genWDaughters, key = lambda genPart: abs(genPart.pdgId), reverse = True))
      if not ((genWDaughters_pdgIdSorted[0].pdgId == -5 * sign(genW.pdgId) and genWDaughters_pdgIdSorted[1].pdgId ==  4 * sign(genW.pdgId)) or \
              (genWDaughters_pdgIdSorted[0].pdgId == -5 * sign(genW.pdgId) and genWDaughters_pdgIdSorted[1].pdgId ==  2 * sign(genW.pdgId)) or \
              (genWDaughters_pdgIdSorted[0].pdgId ==  4 * sign(genW.pdgId) and genWDaughters_pdgIdSorted[1].pdgId == -3 * sign(genW.pdgId)) or \
              (genWDaughters_pdgIdSorted[0].pdgId ==  4 * sign(genW.pdgId) and genWDaughters_pdgIdSorted[1].pdgId == -1 * sign(genW.pdgId)) or \
              (genWDaughters_pdgIdSorted[0].pdgId == -3 * sign(genW.pdgId) and genWDaughters_pdgIdSorted[1].pdgId ==  2 * sign(genW.pdgId)) or \
              (genWDaughters_pdgIdSorted[0].pdgId ==  2 * sign(genW.pdgId) and genWDaughters_pdgIdSorted[1].pdgId == -1 * sign(genW.pdgId))):
        raise ValueError("Invalid hadronic W (%s) decay products from Higgs: %s" % \
          (genW, ', '.join(map(str, genWDaughters)))
        )

    return daughters
          
  if choice in [ SelectionXToHHTobbWWOptions.SAVE_TAU_FROM_WZ_FROM_HIG ]:
    genTau = filter(lambda genPart: abs(genPart.pdgId) == 15 and genPart.genPartIdxMother >0 \
     and abs(genParticles[genPart.genPartIdxMother].pdgId) in [23,24], genParticles)
    return genTau
  
  genWFromHiggs.sort(key=lambda x:x.mass,reverse=True)	
  genW1FromHIG = None; genW2FromHIG = None
  ##W1->lv; W2->lv (for DL) or qq (for SL)
  if len(findQuarksfromW(genWFromHiggs[0])) == 2: ## W->qq
    genW1FromHIG  = genWFromHiggs[1]
    genW2FromHIG  = genWFromHiggs[0] 
  else:
    genW1FromHIG  = genWFromHiggs[0]
    genW2FromHIG  = genWFromHiggs[1]
  isHHbbWWlvqq = len(findQuarksfromW(genW2FromHIG)) == 2 and findLepandNufromW(genW1FromHIG)['Lep'] != [] 
  isHHbbWWlvlv = findLepandNufromW(genW2FromHIG)['Lep'] != [] and findLepandNufromW(genW1FromHIG)['Lep'] != [] 
  ## HHbbWW DL/SL sample usually include W->tau+nu
  if isHHbbWW and not (isHHbbWWlvqq or isHHbbWWlvlv) and debugXToHHTobbWW:
    print("WARNING!!: this event did not fall into single lepton category nor double lepton category")
    printGenParticles(genWFromHiggs, genParticles, "debugging Ws(from Higgs->ww) decays")
    printGenParticles([findLastDecendantWithSameId(genW1FromHIG, genParticles, True)], genParticles,"debugging w1")
    printGenParticles([findLastDecendantWithSameId(genW2FromHIG, genParticles, True)], genParticles,"debugging w2")
    
       
  if choice == SelectionXToHHTobbWWOptions.SAVE_W1_FROM_HIG:
      return [genW1FromHIG]
  if choice == SelectionXToHHTobbWWOptions.SAVE_W2_FROM_HIG:
      return [genW2FromHIG]

  if choice in [ SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_W1_FROM_HIG ,\
    SelectionXToHHTobbWWOptions.SAVE_LEPTONIC_NU_FROM_W1_FROM_HIG ]:
    W1Daughters = findLepandNufromW(genW1FromHIG)
    if choice == SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_W1_FROM_HIG:
        return W1Daughters['Lep']
    else:
        return W1Daughters['Nu']
  elif choice in [ SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_W2_FROM_HIG ,\
    SelectionXToHHTobbWWOptions.SAVE_LEPTONIC_NU_FROM_W2_FROM_HIG ]:
    W2Daughters = findLepandNufromW(genW2FromHIG)
    if choice == SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_W2_FROM_HIG:
        return W2Daughters['Lep']
    else:
        return W2Daughters['Nu']
  elif choice in [ SelectionXToHHTobbWWOptions.SAVE_QUARK_FROM_W2_FROM_HIG]:
    return findQuarksfromW(genW2FromHIG) 
  elif choice in [ SelectionXToHHTobbWWOptions.SAVE_TAU_NU_FROM_W1_FROM_HIG ]:
    # tau_neutrino from W1->tau+nu
    genWDaughters = filter(lambda genPart: abs(genPart.pdgId) == 16 and genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genW1FromHIG, genParticles, True).idx, genParticles)
    return genWDaughters 
  elif choice in [ SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_TAU_FROM_W1_FROM_HIG ]:
    # ele/mu from W1->tau+nu->ele/mu+Nu+Nu+Nu
    genWDaughters = filter(lambda genPart: abs(genPart.pdgId) == 15 and genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genW1FromHIG, genParticles, True).idx, genParticles)
    if len(genWDaughters) == 0 : return []
    genTauDaughters = filter(lambda genPart: abs(genPart.pdgId) in [11,13] and  genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genWDaughters[0], genParticles, True).idx, genParticles)
    return genTauDaughters 
  elif choice in [ SelectionXToHHTobbWWOptions.SAVE_NU_FROM_TAU_FROM_W1_FROM_HIG ]:
    # two Nu from tau decay in W1->tau+nu->ele/mu+nu+Nu+Nu
    genWDaughters = filter(lambda genPart: abs(genPart.pdgId) == 15 and genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genW1FromHIG, genParticles, True).idx, genParticles)
    if len(genWDaughters) == 0 : return []
    genTauDaughters = filter(lambda genPart: abs(genPart.pdgId) in [12,14,16] and  genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genWDaughters[0], genParticles, True).idx, genParticles)
    return genTauDaughters 
  elif choice in [ SelectionXToHHTobbWWOptions.SAVE_TAU_NU_FROM_W2_FROM_HIG ]:
    #Tau_nu from W2->tau+nu
    genWDaughters = filter(lambda genPart: abs(genPart.pdgId) == 16 and genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genW2FromHIG, genParticles, True).idx, genParticles)
    return genWDaughters
  elif choice in [ SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_TAU_FROM_W2_FROM_HIG ]:
    # ele/mu from W2->tau+nu->ele/mu+nu+Nu+Nu
    genWDaughters = filter(lambda genPart: abs(genPart.pdgId) == 15 and genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genW2FromHIG, genParticles, True).idx, genParticles)
    if len(genWDaughters) == 0 : return []
    genTauDaughters = filter(lambda genPart: abs(genPart.pdgId) in [11,13] and  genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genWDaughters[0], genParticles, True).idx, genParticles)
    return genTauDaughters 
  elif choice in [ SelectionXToHHTobbWWOptions.SAVE_NU_FROM_TAU_FROM_W2_FROM_HIG ]:
    # two neutrinos from Tau decay in W2->tau+nu->ele/mu+nu+Nu+Nu
    genWDaughters = filter(lambda genPart: abs(genPart.pdgId) == 15 and genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genW2FromHIG, genParticles, True).idx, genParticles)
    if len(genWDaughters) == 0 : return []
    genTauDaughters = filter(lambda genPart: abs(genPart.pdgId) in [12,14,16] and  genPart.genPartIdxMother \
      == findLastDecendantWithSameId(genWDaughters[0], genParticles, True).idx, genParticles)
    return genTauDaughters 
  else:
    raise ValueError("Choice %i not implemented in  XToHHTobbWW" % choice)
    return []



class genParticleProducer(Module):

  def __init__(self, genEntry, verbose = False):
    self.massTable = MassTable()
    self.branchLenNames  = {}
    self.selections      = {}
    self.branchBaseNames = []

    self.genBranches = {
        "pt"          : "F",
        "eta"         : "F",
        "phi"         : "F",
        "mass"        : "F",
        "pdgId"       : "I",
        "charge"      : "I",
        "status"      : "I",
        "statusFlags" : "I",
	"idx"         : "I",
      }

    for branchBaseName, selection in genEntry.items():
      self.branchBaseNames.append(branchBaseName)
      self.selections[branchBaseName]     = selection
      self.branchLenNames[branchBaseName] = "n%s" % branchBaseName

    if verbose:
      logging.getLogger().setLevel(logging.DEBUG)

  def beginJob(self):
    pass

  def endJob(self):
    pass

  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    for branchBaseName in self.branchBaseNames:
      for branchName, branchType in self.genBranches.items():
        self.out.branch(
          "%s_%s" % (branchBaseName, branchName),
          branchType,
          lenVar = self.branchLenNames[branchBaseName]
        )

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    genParticles  = map(lambda genPartIdx: GenPartAux(genPartIdx[1], genPartIdx[0], self.massTable), enumerate(Collection(event, "GenPart")))
    allgenParticles = Collection(event, "GenPart")
    #print(" type(genParticles) ",type(genParticles)  ," type(allgenParticles) ",type(allgenParticles))

    for branchBaseName in self.branchBaseNames:
      gen_arr = self.selections[branchBaseName](genParticles)
      gen_arr = list(sorted(gen_arr, key = lambda genPart: genPart.pt, reverse = True)) # sort by pT
      #printGenParticles(gen_arr, allgenParticles,branchBaseName)

      for branchName, branchType in self.genBranches.items():
        self.out.fillBranch(
          "%s_%s" % (branchBaseName, branchName),
          map(lambda genPart: getattr(genPart, branchName), gen_arr)
        )
    return True


genLeptonEntry                      = ("GenLep",                         genPromptLeptonSelection)
genLeptonAllEntry                   = ("GenLepAll",                      genLeptonSelection)
genXEntry                           = ("GenX",                       genRadionGravitonSelection)
genHiggsEntry                       = ("GenHiggs",                       genHiggsSelection)
genHiggsDaughtersEntry              = ("GenHiggsDaughters",              genHiggsDaughtersSelection)
genNuEntry                          = ("GenNu",                          genNuSelection)
genWZquarkEntry                     = ("GenWZQuark",                     genWZquarkSelection)
genVbosonEntry                      = ("GenVbosons",                     genVbosonSelection)
genBQuarkFromHIGEntry               = ("GenBQuarkFromHiggs",             (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_BQUARK_FROM_HIG)))
genW1FromHIGEntry                   = ("GenW1FromHiggs",                 (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_W1_FROM_HIG)))
genW2FromHIGEntry                   = ("GenW2FromHiggs",                 (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_W2_FROM_HIG)))
genLepFromW1FromHIGEntry            = ("GenLepFromW1FromHiggs",          (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_W1_FROM_HIG)))
genNuFromW1FromHIGEntry             = ("GenNuFromW1FromHiggs",           (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_LEPTONIC_NU_FROM_W1_FROM_HIG)))
genLepFromW2FromHIGEntry            = ("GenLepFromW2FromHiggs",          (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_W2_FROM_HIG)))
genNuFromW2FromHIGEntry             = ("GenNuFromW2FromHiggs",           (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_LEPTONIC_NU_FROM_W2_FROM_HIG)))
genTauNuFromW1FromHIGEntry          = ("GenTauNuFromW1FromHiggs",        (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_TAU_NU_FROM_W1_FROM_HIG)))
genLepFromTauFromW1FromHIGEntry     = ("GenLepFromTauFromW1FromHiggs",   (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_TAU_FROM_W1_FROM_HIG)))
genNuFromTauFromW1FromHIGEntry      = ("GenNuFromTauFromW1FromHiggs",    (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_NU_FROM_TAU_FROM_W1_FROM_HIG)))
genTauNuFromW2FromHIGEntry          = ("GenTauNuFromW2FromHiggs",        (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_TAU_NU_FROM_W2_FROM_HIG)))
genLepFromTauFromW2FromHIGEntry     = ("GenLepFromTauFromW2FromHiggs",   (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_LEPTON_FROM_TAU_FROM_W2_FROM_HIG)))
genNuFromTauFromW2FromHIGEntry      = ("GenNuFromTauFromW2FromHiggs",    (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_NU_FROM_TAU_FROM_W2_FROM_HIG)))
genQuarkFromW2FromHIGEntry          = ("GenQuarkFromW2FromHiggs",        (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_QUARK_FROM_W2_FROM_HIG)))
genTauFromWZFromHIGEntry            = ("GenTauFromWZFromHiggs",           (lambda genParticles : genXToHHTobbWWSelection(genParticles, SelectionXToHHTobbWWOptions.SAVE_TAU_FROM_WZ_FROM_HIG)))
genTauEntry                         = ("GenTau",                         (lambda genParticles : genTauSelection(genParticles, SelectionOptions.SAVE_TAU)))
genLeptonicTauEntry                 = ("GenLeptonicTau",                 (lambda genParticles : genTauSelection(genParticles, SelectionOptions.SAVE_LEPTONIC_TAU)))
genHadronicTauEntry                 = ("GenHadronicTau",                 (lambda genParticles : genTauSelection(genParticles, SelectionOptions.SAVE_HADRONIC_TAU)))
genLepFromTauEntry                  = ("GenLepFromTau",                  (lambda genParticles : genTauSelection(genParticles, SelectionOptions.SAVE_LEPTON_FROM_TAU)))
genNuLepFromTauEntry                = ("GenNuLepFromTau",                (lambda genParticles : genTauSelection(genParticles, SelectionOptions.SAVE_LEPTONIC_NU_FROM_TAU)))
genNuTauFromLeptonicTauEntry        = ("GenNuTauFromLeptonicTau",        (lambda genParticles : genTauSelection(genParticles, SelectionOptions.SAVE_TAU_NU_FROM_LEPTONIC_TAU)))
genNuTauFromHadronicTauEntry        = ("GenNuTauFromHadronicTau",        (lambda genParticles : genTauSelection(genParticles, SelectionOptions.SAVE_TAU_NU_FROM_HADRONIC_TAU)))
genTopEntry                         = ("GenTop",                         (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_TOP)))
genLepFromTopEntry                  = ("GenLepFromWFromTop",             (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_LEPTON_FROM_W_FROM_TOP)))
genNuFromTopEntry                   = ("GenNuFromWFromTop",              (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_LEPTONIC_NU_FROM_W_FROM_TOP)))
genTauFromTopEntry                  = ("GenTauFromTop",                  (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_TAU_FROM_TOP)))
genNuTauFromTopEntry                = ("GenNuTauFromTop",                (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_TAU_NU_FROM_TOP)))
genNuFromHadronicTauFromTopEntry    = ("GenNuFromHadronicTauFromTop",    (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_TAU_NU_FROM_HADRONIC_TAU_FROM_TOP)))
genNuFromLeptonicTauFromTopEntry    = ("GenNuFromLeptonicTauFromTop",    (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_TAU_NU_FROM_LEPTONIC_TAU_FROM_TOP)))
genNuTauFromLeptonicTauFromTopEntry = ("GenNuTauFromLeptonicTauFromTop", (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_LEPTON_NU_FROM_TAU_FROM_TOP)))
genLepFromTauFromTopEntry           = ("GenLepFromTauFromTop",           (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_LEPTON_FROM_TAU_FROM_TOP)))
genQuarkFromWFromTopEntry           = ("GenQuarkFromWFromTop",           (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_QUARK_FROM_W_FROM_TOP)))
genBQuarkFromTopEntry               = ("GenBQuarkFromTop",               (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_BQUARK_FROM_TOP)))
genWFromTopEntry                    = ("GenWFromTop",               (lambda genParticles : genTopSelection(genParticles, SelectionOptions.SAVE_W_FROM_TOP)))

# provide these variables as the 2nd arguments to the import option for the nano_postproc.py script
genLepton                      = lambda : genParticleProducer(dict([genLeptonEntry]))                      # all prompt stable leptons
genLeptonAll                   = lambda : genParticleProducer(dict([genLeptonAllEntry]))                   # all stable leptons
genHiggs                       = lambda : genParticleProducer(dict([genHiggsEntry]))                       # all Higgs (first in the decay chain)
genHiggsDaughters              = lambda : genParticleProducer(dict([genHiggsDaughtersEntry]))              # all Higgs daughters
genBQuarkFromHiggs             = lambda : genParticleProducer(dict([genBQuarkFromHIGEntry]))
genW1FromHiggs                 = lambda : genParticleProducer(dict([genW1FromHIGEntry]))
genW2FromHiggs                 = lambda : genParticleProducer(dict([genW2FromHIGEntry]))
genLepFromW1FromHiggs          = lambda : genParticleProducer(dict([genLepFromW1FromHIGEntry]))
genNuFromW1FromHiggs           = lambda : genParticleProducer(dict([genNuFromW1FromHIGEntry]))
genLepFromW2FromHiggs          = lambda : genParticleProducer(dict([genLepFromW2FromHIGEntry]))
genNuFromW2FromHiggs           = lambda : genParticleProducer(dict([genNuFromW2FromHIGEntry]))
genLepFromTauFromW1FromHiggs   = lambda : genParticleProducer(dict([genLepFromTauFromW1FromHIGEntry]))
genNuFromTauFromW1FromHiggs    = lambda : genParticleProducer(dict([genNuFromTauFromW1FromHIGEntry]))
genLepFromTauFromW2FromHiggs   = lambda : genParticleProducer(dict([genLepFromTauFromW2FromHIGEntry]))
genNuFromTauFromW2FromHiggs    = lambda : genParticleProducer(dict([genNuFromTauFromW2FromHIGEntry]))
genQuarksFromW2FromHiggs       = lambda : genParticleProducer(dict([genQuarkFromW2FromHIGEntry]))
genTauFromWZFromHiggs          = lambda : genParticleProducer(dict([genTauFromWZFromHIGEntry]))
genTau                         = lambda : genParticleProducer(dict([genTauEntry]))                         # all taus
genNu                          = lambda : genParticleProducer(dict([genNuEntry]))                          # all neutrinos
genWZquark                     = lambda : genParticleProducer(dict([genWZquarkEntry]))                     # all quarks coming from W or Z decay
genVboson                      = lambda : genParticleProducer(dict([genVbosonEntry]))                      # all W and Z bosons (first in the decay chain)
genLeptonicTau                 = lambda : genParticleProducer(dict([genLeptonicTauEntry]))                 # only taus (tau(l)) decaying leptonically: tau(l) -> l vl vtau
genHadronicTau                 = lambda : genParticleProducer(dict([genHadronicTauEntry]))                 # only taus (tau(h)) decaying hadronically: tau(h) -> tauh vtau (NB! NOT RECONSTRUCTED GEN HAD TAU!)
genLepFromTau                  = lambda : genParticleProducer(dict([genLepFromTauEntry]))                  # only leptons (l) coming from tau decay: tau -> l vl vtau
genNuLepFromTau                = lambda : genParticleProducer(dict([genNuLepFromTauEntry]))                # only lepton neutrinos (vl) coming from tau decay: tau -> l vl vtau
genNuTauFromLeptonicTau        = lambda : genParticleProducer(dict([genNuTauFromLeptonicTauEntry]))        # only tau neutrinos (vtau) coming from leptonic tau decay: tau -> l vl vtau
genNuTauFromHadronicTau        = lambda : genParticleProducer(dict([genNuTauFromHadronicTauEntry]))        # only tau neutrinos (vtau) coming from hadronic tau decay: tau -> tauh vtau
genTop                         = lambda : genParticleProducer(dict([genTopEntry]))                         # all tops
genLepFromTop                  = lambda : genParticleProducer(dict([genLepFromTopEntry]))                  # only leptons (l) from t -> W b, W -> l vl decay
genNuFromTop                   = lambda : genParticleProducer(dict([genNuFromTopEntry]))                   # only lepton neutrinos (vl) from t -> W b, W -> l vl decay
genTauFromTop                  = lambda : genParticleProducer(dict([genTauFromTopEntry]))                  # only taus (tau) from t -> W b, W -> tau vtau decay
genNuTauFromTop                = lambda : genParticleProducer(dict([genNuTauFromTopEntry]))                # only tau neutrinos (vtau) from t -> W b, W -> tau vtau decay
genNuFromHadronicTauFromTop    = lambda : genParticleProducer(dict([genNuFromHadronicTauFromTopEntry]))    # only tau neutrinos (vtau2) from t -> W b, W -> tau vtau, tau -> tauh vtau2 decay
genNuFromLeptonicTauFromTop    = lambda : genParticleProducer(dict([genNuFromLeptonicTauFromTopEntry]))    # only tau neutrinos (vtau2) from t -> W b, W -> tau vtau, tau -> l vl vtau2 decay
genNuTauFromLeptonicTauFromTop = lambda : genParticleProducer(dict([genNuTauFromLeptonicTauFromTopEntry])) # only lepton neutrinos (vl) from t -> W b, W -> tau vtau, tau -> l vl vtau2 decay
genLepFromTauFromTop           = lambda : genParticleProducer(dict([genLepFromTauFromTopEntry]))           # only leptons (l) from t -> W b, W -> tau vtau, tau -> l vl vtau2 decay
genQuarkFromWFromTop           = lambda : genParticleProducer(dict([genQuarkFromWFromTopEntry]))           # only quarks (q, q') from t -> W b, W -> q q' decay
genBQuarkFromTop               = lambda : genParticleProducer(dict([genBQuarkFromTopEntry]))                # only b quarks from t -> W b
genWFromTop                    = lambda : genParticleProducer(dict([genWFromTopEntry]))               # only W from t -> W b

genAll = lambda : genParticleProducer(dict([
    genLeptonEntry,
    genLeptonAllEntry,
    genXEntry,
    genHiggsEntry,
    #genHiggsDaughtersEntry,
    genBQuarkFromHIGEntry,
    genW1FromHIGEntry,
    genW2FromHIGEntry,
    genLepFromW1FromHIGEntry,
    genNuFromW1FromHIGEntry,
    genLepFromW2FromHIGEntry,
    genNuFromW2FromHIGEntry,
    genNuEntry,
    genWZquarkEntry,
    genVbosonEntry,
    genTauEntry,
    genLeptonicTauEntry,
    genHadronicTauEntry,
    genLepFromTauEntry,
    genNuLepFromTauEntry,
    genNuTauFromLeptonicTauEntry,
    genNuTauFromHadronicTauEntry,
    genTopEntry,
    genLepFromTopEntry,
    genNuFromTopEntry,
    genTauFromTopEntry,
    genNuTauFromTopEntry,
    genNuFromHadronicTauFromTopEntry,
    genNuFromLeptonicTauFromTopEntry,
    genNuTauFromLeptonicTauFromTopEntry,
    genLepFromTauFromTopEntry,
    genQuarkFromWFromTopEntry,
    genBQuarkFromTopEntry,
  ]))

genHHAndTTbar = lambda : genParticleProducer(dict([
    genXEntry,
    genHiggsEntry,
    genBQuarkFromHIGEntry,
    genW1FromHIGEntry,
    genW2FromHIGEntry,
    genTauFromWZFromHIGEntry,
    genLepFromW1FromHIGEntry,
    genNuFromW1FromHIGEntry,
    genLepFromW2FromHIGEntry,
    genNuFromW2FromHIGEntry,
    genQuarkFromW2FromHIGEntry,
    genTauNuFromW1FromHIGEntry,
    genTauNuFromW2FromHIGEntry,
    genLepFromTauFromW1FromHIGEntry,
    genNuFromTauFromW1FromHIGEntry,
    genLepFromTauFromW2FromHIGEntry,
    genNuFromTauFromW2FromHIGEntry,
    #genHiggsDaughtersEntry,
    #below is ttbar, for SL and DL
    genTopEntry,
    genBQuarkFromTopEntry,
    genWFromTopEntry,
    genLepFromTopEntry,
    genNuFromTopEntry,
    genLepFromTauFromTopEntry,
    genQuarkFromWFromTopEntry,
  ]))

##X->HH->bbWW->bblvlv or bbtauvtauv with tau ->lvv
genHH = lambda : genParticleProducer(dict([
    genXEntry,
    genHiggsEntry,
    genBQuarkFromHIGEntry,
    genW1FromHIGEntry,
    genW2FromHIGEntry,
    genTauFromWZFromHIGEntry,
    genLepFromW1FromHIGEntry,
    genNuFromW1FromHIGEntry,
    genLepFromW2FromHIGEntry,
    genNuFromW2FromHIGEntry,
    genTauNuFromW1FromHIGEntry,
    genTauNuFromW2FromHIGEntry,
    genLepFromTauFromW1FromHIGEntry,
    genNuFromTauFromW1FromHIGEntry,
    genLepFromTauFromW2FromHIGEntry,
    genNuFromTauFromW2FromHIGEntry,
    genQuarkFromW2FromHIGEntry,
  ]))

genTest = lambda : genParticleProducer(dict([
    genHiggsEntry,
]))

genTTbar = lambda : genParticleProducer(dict([
    genTopEntry,
    genBQuarkFromTopEntry,
    genWFromTopEntry,
    genLepFromTopEntry,
    genNuFromTopEntry,
    genLepFromTauFromTopEntry,
    genQuarkFromWFromTopEntry,
  ]))

