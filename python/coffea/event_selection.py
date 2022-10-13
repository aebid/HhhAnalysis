import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
import ROOT
from coffea.nanoevents.methods import vector


def single_lepton_category(EventProcess):
    #Pass MET filters
    #At least 1 fakeable lepton
    #If the leading cone-pT lepton is e (mu), pass single e (mu) trigger
    #cone-pt > 32 (25) for e (mu)
    #Invariant mass of each pair of preselected leptons (electrons NOT cleaned) must be greater than 12 GeV
    #Not more than 1 tight lepton - tight should be same as highest cone pT fakeable
    #Tau veto: no tau passing pt > 20, abs(eta) < 2.3, abs(dxy) <= 1000, abs(dz) <= 0.2, decay modes = {0, 1, 2, 10, 11}, and byMediumDeepTau2017v2VSjet, byVLooseDeepTau2017v2VSmu, byVVVLooseDeepTau2017v2VSe. Taus overlapping with fakeable electrons or fakeable muons within dR < 0.3 are not considered for the tau veto
    #No pair of same-flavor, opposite sign preselected leptons within 10 GeV of the Z mass
    #At least 1 medium btag (that can be on a AK8 jet): (#selJetsAK8_b >= 1 || #b-medium >= 1)
    #Minimal number of jets to construct an Hbb and admit an hadronic W with a missing jet: (#selJetsAK8_b == 0 && #selJetsAK4 >= 3) || (#selJetsAK8_b >= 1 && nJet_that_not_bb >= 1)

    cutflow_val = 0

    muons = EventProcess.muons
    electrons = EventProcess.electrons
    taus = EventProcess.taus
    ak4_jets = EventProcess.ak4_jets
    ak8_jets = EventProcess.ak8_jets
    flag = EventProcess.flag
    HLT = EventProcess.HLT
    events = EventProcess.events
    Runyear = EventProcess.Runyear
    isMC = EventProcess.isMC


    print("N events: ", len(events))
    leptons_preselected = ak.concatenate([electrons.mask[electrons.preselected], muons.mask[muons.preselected]], axis=1)
    leptons_fakeable = ak.concatenate([electrons.mask[electrons.fakeable], muons.mask[muons.fakeable]], axis=1)
    leptons_tight = ak.concatenate([electrons.mask[electrons.tight], muons.mask[muons.tight]], axis=1)

    leptons_preselected_sorted = leptons_preselected[ak.argsort(leptons_preselected.conept, axis=1, ascending=False)]
    leptons_fakeable_sorted = leptons_fakeable[ak.argsort(leptons_fakeable.conept, axis=1, ascending=False)]
    leptons_tight_sorted = leptons_tight[ak.argsort(leptons_tight.conept, axis=1, ascending=False)]

    leading_leptons = ak.firsts(leptons_fakeable_sorted, axis=1)



    #Require at least 1 fakeable (or tight) lepton
    one_fakeable_lepton = ak.sum(leptons_fakeable_sorted.fakeable, axis=1) >= 1

    #Require MET filters
    MET_filters = (
        (
            (flag.eeBadScFilter) & (isMC == 0) | (isMC == 1)
        ) &
        (flag.goodVertices) & (flag.globalSuperTightHalo2016Filter) & (flag.HBHENoiseFilter) &
        (flag.HBHENoiseIsoFilter) & (flag.EcalDeadCellTriggerPrimitiveFilter) & (flag.BadPFMuonFilter)
    )

    single_step1_mask = (one_fakeable_lepton & MET_filters)


    events["single_lepton"] = ak.where(
        single_step1_mask,
        True,
        False
    )

    print("N single events step1: ", ak.sum(events.single_lepton))
    cutflow_val += 1

    #Require no more than 1 tight lepton
    one_tight_lepton = ak.sum(leptons_tight_sorted.tight, axis=1) <= 1
    #one_tight_lepton = ( (ak.sum(muons.tight, axis=1) + ak.sum(electrons.tight, axis=1)) <= 1)

    single_step2_mask = one_tight_lepton

    events["single_lepton"] = ak.where(
        events["single_lepton"] & single_step2_mask,
        True,
        False
    )

    print("N single events step2: ", ak.sum(events.single_lepton))
    cutflow_val += 1



    #Z mass and invariant mass cuts
    #No pair of same-flavor, opposite-sign preselected leptons within 10GeV of the Z mass (91.1876)
    #Invariant mass of each pair of preselected leptons (electrons not cleaned) must be greater than 12 GeV
    #leptons_preselected = ak.concatenate([electrons.mask[electrons.preselected], muons.mask[muons.preselected]], axis=1)
    lep_pairs_for_Zmass_and_Invarmass = ak.combinations(leptons_preselected, 2)
    first_leps, second_leps = ak.unzip(lep_pairs_for_Zmass_and_Invarmass)

    lep1_lorentz_vec = ak.zip(
        {
            "pt": ak.fill_none(first_leps.pt, 0.0),
            "eta": ak.fill_none(first_leps.eta, 0.0),
            "phi": ak.fill_none(first_leps.phi, 0.0),
            "mass": ak.fill_none(first_leps.mass, 0.0),
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior,
    )

    lep2_lorentz_vec = ak.zip(
        {
            "pt": ak.fill_none(second_leps.pt, 0.0),
            "eta": ak.fill_none(second_leps.eta, 0.0),
            "phi": ak.fill_none(second_leps.phi, 0.0),
            "mass": ak.fill_none(second_leps.mass, 0.0),
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior,
    )

    Zmass_and_invar_mass_cut = ak.any(
        (
            (ak.is_none(first_leps, axis=1) | ak.is_none(second_leps, axis=1))
            | #OR
            (
                ( (lep1_lorentz_vec + lep2_lorentz_vec).mass > 12.0) &
                (
                    (ak.fill_none(first_leps.charge, 0) == ak.fill_none(second_leps.charge, 0))
                    | #OR
                    (abs(ak.fill_none(first_leps.pdgId, 0)) != abs(ak.fill_none(second_leps.pdgId, 0)))
                    | #OR
                    (abs((lep1_lorentz_vec + lep2_lorentz_vec).mass - 91.1876) > 10.0)
                )
            )
        ) == 0, axis = 1
    )

    single_step3_mask = (Zmass_and_invar_mass_cut == 0)

    events["single_lepton"] = ak.where(
        events["single_lepton"] & single_step3_mask,
        True,
        False
    )

    print("N single events step3: ", ak.sum(events.single_lepton))
    cutflow_val += 1


    #If the leading cone-pT lepton is e (mu), pass single e (mu) trigger


    if Runyear == 2016:
        electron_trigger_cuts = ((HLT.Ele27_WPTight_Gsf) | (HLT.Ele25_eta2p1_WPTight_Gsf) | (HLT.Ele27_eta2p1_WPLoose_Gsf))
        muon_trigger_cuts = ((HLT.IsoMu22) | (HLT.IsoTkMu22) | (HLT.IsoMu22_eta2p1) | (HLT.IsoTkMu22_eta2p1) | (HLT.IsoMu24) | (HLT.IsoTkMu24))
    elif Runyear == 2017:
        electron_trigger_cuts = ((HLT.Ele35_WPTight_Gsf) | (HLT.Ele32_WPTight_Gsf))
        muon_trigger_cuts = ((HLT.IsoMu24) | (HLT.IsoMu27))
    elif Runyear == 2018:
        electron_trigger_cuts = ((HLT.Ele32_WPTight_Gsf) & (Runyear == 2018))
        muon_trigger_cuts = ((HLT.IsoMu24) | (HLT.IsoMu27))
    else:
        electron_trigger_cuts = ak.zeros_like(HLT)
        muon_trigger_cuts = ak.zeros_like(HLT)

    single_step4_mask = ak.where(
        (abs(leading_leptons.pdgId) == 11),
            electron_trigger_cuts,
            muon_trigger_cuts
    )

    events["single_lepton"] = ak.where(
        events["single_lepton"] & single_step4_mask,
            True,
            False
    )

    print("N single events step4: ", ak.sum(events.single_lepton))
    cutflow_val += 1


    leading_conept_cut = ak.where(
        (abs(leading_leptons.pdgId) == 11),
            leading_leptons.conept >= 32.0,
            leading_leptons.conept >= 25.0
    )

    single_step5_mask = leading_conept_cut

    events["single_lepton"] = ak.where(
        events["single_lepton"] & single_step5_mask,
            True,
            False
    )

    print("N single events step5: ", ak.sum(events.single_lepton))
    cutflow_val += 1


    leading_lepton_MC_match = (leading_leptons.genPartFlav == 1) | (leading_leptons.genPartFlav == 15)

    single_step6_mask = leading_lepton_MC_match

    if isMC:

        events["single_lepton"] = ak.where(
            events["single_lepton"] & single_step6_mask,
                True,
                False
        )

    print("Is MC? Doing step 6: ", isMC)
    print("N single events step6: ", ak.sum(events.single_lepton))
    cutflow_val += 1



    #Tau veto: no tau passing pt>20, abs(eta) < 2.3, abs(dxy) <= 1000, abs(dz) <= 0.2, "decayModeFindingNewDMs", decay modes = {0, 1, 2, 10, 11}, and "byMediumDeepTau2017v2VSjet", "byVLooseDeepTau2017v2VSmu", "byVVVLooseDeepTau2017v2VSe". Taus overlapping with fakeable electrons or fakeable muons within dR < 0.3 are not considered for the tau veto
    #False -> Gets Removed : True -> Passes veto
    tau_veto_pairs = ak.cartesian([taus, leptons_fakeable_sorted], nested=True)
    taus_for_veto, leps_for_veto = ak.unzip(tau_veto_pairs)

    tau_veto_cleaning = ak.min(abs(taus_for_veto.delta_r(leps_for_veto)), axis=2) >= 0.3

    tau_veto_selection = (
    (taus.pt > 20) & (abs(taus.eta) < 2.3) & (abs(taus.dxy) <= 1000.0) & (abs(taus.dz) <= 0.2) & (taus.idDecayModeNewDMs) &
    (
        (taus.decayMode == 0) | (taus.decayMode == 1) | (taus.decayMode == 2) | (taus.decayMode == 10) | (taus.decayMode == 11)
    ) &
    (taus.idDeepTau2017v2p1VSjet >= 16) & (taus.idDeepTau2017v2p1VSmu >= 1) & (taus.idDeepTau2017v2p1VSe >= 1)
    )

    tau_veto = ak.any(tau_veto_cleaning & tau_veto_selection, axis=1) == 0

    single_step7_mask = tau_veto

    events["single_lepton"] = events.single_lepton & single_step7_mask

    print("N single events step7: ", ak.sum(events.single_lepton))
    cutflow_val += 1


    #Jet cuts
    #1 or more btagged ak8_jets or 1 or more btagged ak4_jets
    one_btagged_jet = (ak.sum(ak4_jets.medium_btag_single, axis=1) >= 1) | (ak.sum(ak8_jets.btag_single, axis=1) >= 1)

    single_step8_mask = one_btagged_jet

    events["single_lepton"] = events.single_lepton & single_step8_mask



    print("N single events step8: ", ak.sum(events.single_lepton))
    cutflow_val += 1



    #Count ak4 jets that are dR 1.2 away from btagged ak8 jets ***Is this only from leading ak8 jet or from all?***
    #Require either:
        # 0 ak8 btagged jets and 3 or more cleaned ak4 jets
        # or
        # 1 or more ak8 btagged jets and 1 or more cleaned ak4 jets 1.2dR away from an ak8 bjet

    ak8_jets_btag_single = ak8_jets.mask[ak8_jets.btag_single]
    ak8_jets_btag_single_sorted = ak8_jets_btag_single[ak.argsort(ak8_jets_btag_single.pt, axis=1, ascending=False)]

    clean_ak4_jets_btagged_ak8_jets = ak.cartesian([ak4_jets.mask[ak4_jets.cleaned_single], ak8_jets_btag_single_sorted], nested=True)
    clean_ak4_for_veto, btag_ak8_for_veto = ak.unzip(clean_ak4_jets_btagged_ak8_jets)

    ak4_jets.jets_that_not_bb = abs(clean_ak4_for_veto.delta_r(btag_ak8_for_veto)) > 1.2
    n_jets_that_not_bb = ak.firsts(ak.sum(ak4_jets.jets_that_not_bb, axis=1))

    jet_btag_veto = (
        ( (n_jets_that_not_bb >= 1) & (ak.sum(ak8_jets.btag_single, axis=1) >= 1) )
        |
        ( (ak.sum(ak8_jets.btag_single, axis=1) == 0) & (ak.sum(ak4_jets.cleaned_single, axis=1) >= 3) )
    )

    jet_btag_veto_1 = (n_jets_that_not_bb >= 1) & (ak.sum(ak8_jets.btag_single, axis=1) >= 1)
    jet_btag_veto_2 = (ak.sum(ak8_jets.btag_single, axis=1) == 0) & (ak.sum(ak4_jets.cleaned_single, axis=1) >= 3)

    single_step9_mask = jet_btag_veto

    events["single_lepton"] = events.single_lepton & single_step9_mask

    print("N single events step9: ", ak.sum(events.single_lepton))
    cutflow_val += 1


    """
      category_string = "Single"
      if len(ak8jets_btagged) >= 1:##boosted
        category_string += "_HbbFat_WjjRes"
        if len(Jet_that_not_bb) >= 2: category_string += "_allReco"
        else: category_string += "_MissJet"
      else:## resolved
        category_string += "_Res"
        if len(jets_clean) >= 4:
          category_string += "_allReco"
          if len(jets_btagged) >= 1: category_string += "_2b"
          else: category_string += "_1b"
        else:
          category_string += "_MissWJet"
          if len(jets_btagged) >= 1: category_string += "_2b"
          else: category_string += "_1b"
      #if len(tight_leptons) == 1: category_string += "_Signal"
      #if len(tight_leptons) == 0: category_string += "_Fake"
      if len(tight_leptons) == 1 and  tight_leptons[0] == leading_lepton: category_string += "_Signal"
      else: category_string += "_Fake"
      SF = self.get_lepton_SF(leading_lepton)
      self.lepton_IDSF_recoToLoose = SF[0]
      self.lepton_IDSF_looseToTight = SF[1]
      self.lepton_IDSF = self.lepton_IDSF_recoToLoose*self.lepton_IDSF_looseToTight

      self.get_trigger_eff_SF(lep_type, fake_leptons)

      if "Fake" in category_string:
        print("Is fake single, testing jet->lepton fakerate")
        self.get_jet_to_lepton_fake_rate_SF(leading_lepton)
      return category_string

    """

    EventProcess.single_cutflow.append(cutflow_val)
