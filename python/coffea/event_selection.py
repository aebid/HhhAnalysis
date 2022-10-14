import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
import ROOT
from coffea.nanoevents.methods import vector

def increment_cutflow(events, cut, cutflow_name):
    events[cutflow_name] = ak.where(
        cut,
            (events[cutflow_name] + 1),
            events[cutflow_name]
    )

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

    debug = EventProcess.debug

    events["single_cutflow"] = np.zeros_like(events.run)

    if debug: print("N events: ", len(events))
    leptons_preselected = ak.concatenate([electrons.mask[electrons.preselected], muons.mask[muons.preselected]], axis=1)
    leptons_fakeable = ak.concatenate([electrons.mask[electrons.fakeable], muons.mask[muons.fakeable]], axis=1)
    leptons_tight = ak.concatenate([electrons.mask[electrons.tight], muons.mask[muons.tight]], axis=1)

    leptons_preselected_sorted = leptons_preselected[ak.argsort(leptons_preselected.conept, axis=1, ascending=False)]
    leptons_fakeable_sorted = leptons_fakeable[ak.argsort(leptons_fakeable.conept, axis=1, ascending=False)]
    leptons_tight_sorted = leptons_tight[ak.argsort(leptons_tight.conept, axis=1, ascending=False)]

    leading_leptons = leptons_fakeable_sorted[:,0]


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

    if debug: print("N single events step1: ", ak.sum(events.single_lepton))
    increment_cutflow(events, events.single_lepton, "single_cutflow")


    #Require no more than 1 tight lepton
    one_tight_lepton = ak.sum(leptons_tight_sorted.tight, axis=1) <= 1
    #one_tight_lepton = ( (ak.sum(muons.tight, axis=1) + ak.sum(electrons.tight, axis=1)) <= 1)

    single_step2_mask = one_tight_lepton

    events["single_lepton"] = events.single_lepton & single_step2_mask

    if debug: print("N single events step2: ", ak.sum(events.single_lepton))
    increment_cutflow(events, events.single_lepton, "single_cutflow")


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


    events["single_lepton"] = events.single_lepton & single_step3_mask


    if debug: print("N single events step3: ", ak.sum(events.single_lepton))
    increment_cutflow(events, events.single_lepton, "single_cutflow")


    #If the leading cone-pT lepton is e (mu), pass single e (mu) trigger

    single_step4_mask = ak.where(
        ak.is_none(leading_leptons) == 0,
            ak.where(
                (abs(leading_leptons.pdgId) == 11),
                    EventProcess.electron_trigger_cuts,
                    EventProcess.muon_trigger_cuts
            ),
            False
    )

    events["single_lepton"] = events.single_lepton & single_step4_mask

    if debug: print("N single events step4: ", ak.sum(events.single_lepton))
    increment_cutflow(events, events.single_lepton, "single_cutflow")


    leading_conept_cut = ak.where(
        ak.is_none(leading_leptons) == 0,
            ak.where(
                (abs(leading_leptons.pdgId) == 11),
                    leading_leptons.conept >= 32.0,
                    leading_leptons.conept >= 25.0
            ),
            False
    )

    single_step5_mask = leading_conept_cut


    events["single_lepton"] = events.single_lepton & single_step5_mask


    if debug: print("N single events step5: ", ak.sum(events.single_lepton))
    increment_cutflow(events, events.single_lepton, "single_cutflow")


    leading_lepton_MC_match = ak.where(
        ak.is_none(leading_leptons) == 0,
            (leading_leptons.genPartFlav == 1) | (leading_leptons.genPartFlav == 15),
            False
    )

    single_step6_mask = leading_lepton_MC_match

    if isMC:
        events["single_lepton"] = events.single_lepton & single_step6_mask


    if debug: print("Is MC? Doing step 6: ", isMC)
    if debug: print("N single events step6: ", ak.sum(events.single_lepton))
    increment_cutflow(events, events.single_lepton, "single_cutflow")


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

    if debug: print("N single events step7: ", ak.sum(events.single_lepton))
    increment_cutflow(events, events.single_lepton, "single_cutflow")


    #Jet cuts
    #1 or more btagged ak8_jets or 1 or more btagged ak4_jets
    one_btagged_jet = (ak.sum(ak4_jets.medium_btag_single, axis=1) >= 1) | (ak.sum(ak8_jets.btag_single, axis=1) >= 1)

    single_step8_mask = one_btagged_jet

    events["single_lepton"] = events.single_lepton & single_step8_mask

    if debug: print("N single events step8: ", ak.sum(events.single_lepton))
    increment_cutflow(events, events.single_lepton, "single_cutflow")



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
    n_jets_that_not_bb = (ak.sum(ak4_jets.jets_that_not_bb, axis=1))[:,0]

    jet_btag_veto = (
        ( (n_jets_that_not_bb >= 1) & (ak.sum(ak8_jets.btag_single, axis=1) >= 1) )
        |
        ( (ak.sum(ak8_jets.btag_single, axis=1) == 0) & (ak.sum(ak4_jets.cleaned_single, axis=1) >= 3) )
    )

    single_step9_mask = jet_btag_veto

    events["single_lepton"] = events.single_lepton & single_step9_mask

    if debug: print("N single events step9: ", ak.sum(events.single_lepton))
    increment_cutflow(events, events.single_lepton, "single_cutflow")


    events["single_lepton_category"] = ak.where(
        events.single_lepton,
            ak.where(
                ak.sum(ak8_jets.btag_single, axis=1) >= 1, #Boosted Category
                    ak.where(
                        n_jets_that_not_bb >= 2,
                            "Single_HbbFat_WjjRes_AllReco",
                            "Single_HbbFat_WjjRes_MissJet"
                    ), #End Boosted Category
                    ak.where( #Resolved Category
                        ak.sum(ak4_jets.cleaned_single, axis=1) >= 4,
                            ak.where(
                                ak.sum(ak4_jets.medium_btag_single, axis=1) >= 1,
                                    "Single_Res_allReco_2b",
                                    "Single_Res_allReco_1b"
                            ),
                            ak.where(
                                ak.sum(ak4_jets.medium_btag_single, axis=1) >= 1,
                                    "Single_Res_MissWJet_2b",
                                    "Single_Res_MissWJet_1b"
                            ) #End Resolved Category
                    )
            ),
            False
    )


    events["single_lepton_signal_or_fake"] = ak.where(
        events.single_lepton,
            ak.where(
                (ak.sum(leptons_tight_sorted.tight, axis=1) == 1) & leading_leptons.tight,
                    "Signal",
                    "Fake"
            ),
            False
    )

    if debug: print("Single Lep Categories: ", events.single_lepton_category)
    if debug: print("Single Lep Signal or Fake: ", events.single_lepton_signal_or_fake)



def double_lepton_category(EventProcess):
    #Pass MET filters
    #At least 2 fakeable leptons (choose the leading 2 in cone pT in the following)
    #if both fakeable leptons are electrons, the event needs to pass either the single electron or the double electron trigger; if both fakeable leptons are muons, the event needs to pass either the single muon or the double muon trigger; if one fakeable lepton is an electron and the other fakeable lepton is a muon, the event needs to pass either the single electron or the single muon or the muon+electron trigger
    #cone pT > 25 for the leading fakeable lepton and cone pT > 15 GeV for the subleading fakeable lepton
    #The 2 fakeable leptons must have opposite charge
    #No pair of same-flavor, opposite-sign preselected leptons within 10 GeV of the Z mass
    #The event needs to pass the selection in either the boosted or the resolved category:
    #In order for the event to pass the selection in the boosted category, it needs to contain at least one b-tagged Ak8 jet (see object definition)
    #In order for the event to pass the selection in the resolved category, it needs to contain at least 2 AK4 jets, of which at least one passes the medium working-point of the DeepJet b-tagging algorithm
    #The two categories are mutually exclusive. The higher priority is given to the boosted category, i.e. events satisfying the criteria for the boosted category are not considered for the resolved category.
    #The leading and subleading fakeable lepton both pass the tight lepton selection criteria
    #In MC, require MC matching of the leading and subleading fakeable lepton
    #Either the leading fakeable lepton, the subleading fakeable lepton or both fail the tight lepton selection criteria
    #In MC, require MC matching of the leading and subleading fakeable lepton

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

    debug = EventProcess.debug

    events["double_cutflow"] = np.zeros_like(events.run)

    if debug: print("N events: ", len(events))
    leptons_preselected = ak.concatenate([electrons.mask[electrons.preselected], muons.mask[muons.preselected]], axis=1)
    leptons_fakeable = ak.concatenate([electrons.mask[electrons.fakeable], muons.mask[muons.fakeable]], axis=1)
    leptons_tight = ak.concatenate([electrons.mask[electrons.tight], muons.mask[muons.tight]], axis=1)

    leptons_preselected_sorted = leptons_preselected[ak.argsort(leptons_preselected.conept, axis=1, ascending=False)]
    leptons_fakeable_sorted = leptons_fakeable[ak.argsort(leptons_fakeable.conept, axis=1, ascending=False)]
    leptons_tight_sorted = leptons_tight[ak.argsort(leptons_tight.conept, axis=1, ascending=False)]

    leading_leptons = leptons_fakeable_sorted[:,0]
    subleading_leptons = leptons_fakeable_sorted[:,1]


    #Require at least 2 fakeable (or tight) lepton
    one_fakeable_lepton = ak.sum(leptons_fakeable_sorted.fakeable, axis=1) >= 2

    #Require MET filters
    MET_filters = (
        (
            (flag.eeBadScFilter) & (isMC == 0) | (isMC == 1)
        ) &
        (flag.goodVertices) & (flag.globalSuperTightHalo2016Filter) & (flag.HBHENoiseFilter) &
        (flag.HBHENoiseIsoFilter) & (flag.EcalDeadCellTriggerPrimitiveFilter) & (flag.BadPFMuonFilter)
    )

    double_step1_mask = (one_fakeable_lepton & MET_filters)

    events["double_lepton"] = double_step1_mask

    if debug: print("N double events step1: ", ak.sum(events.double_lepton))
    increment_cutflow(events, events.double_lepton, "double_cutflow")


    #Leading lepton cone pT > 25.0, subleading lepton cone pT > 15.0, leading lepton charge != subleading lepton charge
    cone_pt_cuts = (leading_leptons.conept > 25.0) & (subleading_leptons.conept > 15.0)
    charge_cuts = (leading_leptons.charge != subleading_leptons.charge)
    double_step2_mask = ak.where(
        (ak.is_none(leading_leptons) == 0) & (ak.is_none(subleading_leptons) == 0),
            (cone_pt_cuts & charge_cuts),
            False
    )

    events["double_lepton"] = events.double_lepton & double_step2_mask

    if debug: print("N double events step2: ", ak.sum(events.double_lepton))
    increment_cutflow(events, events.double_lepton, "double_cutflow")


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

    double_step3_mask = (Zmass_and_invar_mass_cut == 0)


    events["double_lepton"] = events.double_lepton & double_step3_mask

    if debug: print("N double events step3: ", ak.sum(events.double_lepton))
    increment_cutflow(events, events.double_lepton, "double_cutflow")


    #HLT cuts
    #If Leps are MuMu then pass MuMu or Mu trigger
    #If Leps are ElEl then pass ElEl or El trigger
    #If Leps are MuEl then pass MuEl or El or Mu trigger

    HLT_cuts = ak.where(
        (ak.is_none(leading_leptons) == 0) & (ak.is_none(subleading_leptons) == 0),
            ak.where(
                (abs(leading_leptons.pdgId) == 11) & (abs(subleading_leptons.pdgId) == 11), #ElEl
                    EventProcess.electron_trigger_cuts | EventProcess.double_electron_trigger_cuts,
                    ak.where(
                        (abs(leading_leptons.pdgId) == 13) & (abs(subleading_leptons.pdgId) == 13), #MuMu
                            EventProcess.muon_trigger_cuts | EventProcess.double_muon_trigger_cuts,
                            ak.where(
                                ((abs(leading_leptons.pdgId) == 11) & (abs(subleading_leptons.pdgId) == 13)) | ((abs(leading_leptons.pdgId) == 13) & (abs(subleading_leptons.pdgId) == 11)), #MuEl or ElMu
                                    EventProcess.muon_electron_trigger_cuts | EventProcess.muon_trigger_cuts | EventProcess.electron_trigger_cuts,
                                    False
                            )
                    )
            ),
            False
    )

    double_step4_mask = HLT_cuts

    events["double_lepton"] = events.double_lepton & double_step4_mask

    if debug: print("N double events step4: ", ak.sum(events.double_lepton))
    increment_cutflow(events, events.double_lepton, "double_cutflow")


    #MC match for leading and subleading lepton on MC
    leading_lepton_MC_match = ak.where(
        ak.is_none(leading_leptons) == 0,
            (leading_leptons.genPartFlav == 1) | (leading_leptons.genPartFlav == 15),
            False
    )

    subleading_lepton_MC_match = ak.where(
        ak.is_none(subleading_leptons) == 0,
            (subleading_leptons.genPartFlav == 1) | (subleading_leptons.genPartFlav == 15),
            False
    )

    double_step5_mask = leading_lepton_MC_match & subleading_lepton_MC_match


    if isMC:
        events["double_lepton"] = events.double_lepton & double_step5_mask


    if debug: print("Is MC? Doing step 5: ", isMC)
    if debug: print("N double events step5: ", ak.sum(events.double_lepton))
    increment_cutflow(events, events.double_lepton, "double_cutflow")


    #No more than 2 tight leptons
    n_tight_leptons = ak.sum(leptons_tight.tight, axis=1)

    double_step6_mask = n_tight_leptons <= 2

    events["double_lepton"] = events.double_lepton & double_step6_mask

    if debug: print("N double events step6: ", ak.sum(events.double_lepton))
    increment_cutflow(events, events.double_lepton, "double_cutflow")


    leading_ak8_jet_cleaned = ak8_jets[ak.argsort(ak8_jets.pt, axis=1, ascending=False)][:,0]
    events["double_lepton_category"] = ak.where(
        events.double_lepton,
            ak.where(
                (ak.sum(ak8_jets.cleaned_double, axis=1) >= 1) & leading_ak8_jet_cleaned.btag_double,
                    "Double_HbbFat",
                    ak.where(
                        (ak.sum(ak4_jets.cleaned_double, axis=1) >= 2) & (ak.sum(ak4_jets.medium_btag_double, axis=1) == 1),
                            "Double_Res_1b",
                            ak.where(
                                (ak.sum(ak4_jets.cleaned_double, axis=1) >= 2) & (ak.sum(ak4_jets.medium_btag_double, axis=1) >= 2),
                                    "Double_Res_2b",
                                    "False"
                            )
                    )
            ),
            False
    )


    events["double_lepton_signal_or_fake"] = ak.where(
        events.double_lepton,
            ak.where(
                (leading_leptons.tight) & (subleading_leptons.tight),
                    "Signal",
                    "Fake"
            ),
            False
    )

    if debug: print("double Lep Categories: ", events.double_lepton_category)
    if debug: print("double Lep Signal or Fake: ", events.double_lepton_signal_or_fake)
