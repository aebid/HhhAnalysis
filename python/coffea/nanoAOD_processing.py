import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
import ROOT
from coffea.nanoevents.methods import vector

import time
startTime = time.time()

print("Starting NanoAOD processing")


fname = "sync_2016_m750.root"
events = NanoEventsFactory.from_root(fname, schemaclass=NanoAODSchema.v7).events()
nLeps = 1 #Single lepton or Di Lepton channels
debug = 0

Runyear = 2016
isMC = True

muons = ak.pad_none(events.Muon, 1)
electrons = ak.pad_none(events.Electron, 1)
taus = ak.pad_none(events.Tau, 1)
ak4_jets = ak.pad_none(events.Jet, 1)
ak8_jets = ak.pad_none(events.FatJet, 1)
ak8_subjets = ak.pad_none(events.SubJet, 1)
HLT = events.HLT
flag = events.Flag

if debug > 0:
    print("Muons: ", muons)
    print("Electrons: ", electrons)
    print("Taus: ", taus)
    print("AK4 Jets: ", ak4_jets)
    print("AK8 Jets: ", ak8_jets)
    print("AK8 SubJets: ", ak8_subjets)
    print("HLT: ", HLT)
    print("Flag: ", flag)



muons["conept"] = ak.where(
    (abs(muons.pdgId) == 13) & (muons.mediumId) & (muons.mvaTTH > 0.50),
        muons.pt,
        0.9 * muons.pt * (1.0 + muons.jetRelIso)
)


electrons["conept"] = ak.where(
    (abs(electrons.pdgId) == 11) & (electrons.mvaTTH > 0.30),
        electrons.pt,
        0.9 * electrons.pt * (1.0 + electrons.jetRelIso)
)

muons.ak4_jets = ak4_jets[muons.jetIdx].mask[(muons.jetIdx >= 0) & (muons.jetIdx < ak.num(ak4_jets))]
electrons.ak4_jets = ak4_jets[electrons.jetIdx].mask[(electrons.jetIdx >= 0) & (electrons.jetIdx < ak.num(ak4_jets))]
ak8_jets.subjet1 = ak8_subjets[ak8_jets.subJetIdx1].mask[(ak8_jets.subJetIdx1 >= 0) & (ak8_jets.subJetIdx1 < ak.num(ak8_subjets))]
ak8_jets.subjet2 = ak8_subjets[ak8_jets.subJetIdx2].mask[(ak8_jets.subJetIdx2 >= 0) & (ak8_jets.subJetIdx2 < ak.num(ak8_subjets))]



if debug > 0:
    print("AK8 jets to subjet1: ", ak8_jets.subjet1)
    print("AK8 jets to subjet2: ", ak8_jets.subjet2)

jetDeepJet_WP_loose  = [0.0613, 0.0521, 0.0494]; jetDeepJet_WP_medium = [0.3093, 0.3033, 0.2770]; jetDeepJet_WP_tight  = [0.7221, 0.7489, 0.7264]
jetDeepJet_min_pt = 20.0; jetDeepJet_max_pt = 45.0

muons.jetDeepJet_Upper_x1 = ak.where(
    (muons.conept - jetDeepJet_min_pt) >= 0.0,
        muons.conept - jetDeepJet_min_pt,
        0.0
)

muons.jetDeepJet_Upper_x2 = ak.where(
    (muons.jetDeepJet_Upper_x1/(jetDeepJet_max_pt - jetDeepJet_min_pt)) <= 1.0,
        (muons.jetDeepJet_Upper_x1/(jetDeepJet_max_pt - jetDeepJet_min_pt)),
        1.0
)

muons.jetDeepJet_Upper = muons.jetDeepJet_Upper_x2 * jetDeepJet_WP_loose[Runyear - 2016] + (1 - muons.jetDeepJet_Upper_x2) * jetDeepJet_WP_medium[Runyear - 2016]




def do_object_selection():

    #Muons
    muon_preselection_mask = (
        (abs(muons.eta) <= 2.4) & (muons.pt >= 5.0) & (abs(muons.dxy) <= 0.05) &
        (abs(muons.dz) <= 0.1) & (muons.miniPFRelIso_all <= 0.4) &
        (muons.sip3d <= 8) & (muons.looseId)
    )

    muon_fakeable_mask = (
        (muons.conept >= 10.0) &
        ak.where(
            ak.is_none(muons.ak4_jets, axis=1),
                (
                    (muons.mvaTTH > 0.50)
                    | #OR
                    (muons.jetRelIso < 0.80)
                ),
                (muons.ak4_jets.btagDeepFlavB <= jetDeepJet_WP_medium[Runyear - 2016]) &
                (
                    (muons.mvaTTH > 0.50)
                    | #OR
                    (muons.jetRelIso < 0.80) & (muons.ak4_jets.btagDeepFlavB <= muons.jetDeepJet_Upper)
                )
            )
    )

    muon_tight_mask = (
        (muons.mvaTTH >= 0.50) & (muons.mediumId)
    )

    muons["preselected"] = ak.where(
        muon_preselection_mask,
            True,
            False
    )

    muons["fakeable"] = ak.where(
        muon_fakeable_mask & muons.preselected,
            True,
            False
    )

    muons["tight"] = ak.where(
        muon_tight_mask & muons.fakeable,
            True,
            False
    )



    #Electrons
    electron_preselection_mask = (
        (abs(electrons.eta) <= 2.5) & (electrons.pt >= 7.0)& (abs(electrons.dxy) <= 0.05) &
        (abs(electrons.dz) <= 0.1) & (electrons.miniPFRelIso_all <= 0.4) &
        (electrons.sip3d <= 8) & (electrons.mvaFall17V2noIso_WPL) & (electrons.lostHits <= 1)
    )

    electron_fakeable_mask = (
        (electrons.conept >= 10.0) & (electrons.hoe <= 0.10) & (electrons.eInvMinusPInv >= -0.04) &
        (electrons.lostHits == 0) & (electrons.convVeto) &
        (
            (abs(electrons.eta + electrons.deltaEtaSC) > 1.479) & (electrons.sieie <= 0.030)
            | #OR
            (abs(electrons.eta + electrons.deltaEtaSC) <= 1.479) & (electrons.sieie <= 0.011)
        ) &
        ak.where(
            ak.is_none(electrons.ak4_jets, axis=1),
                True,
                ak.where(
                    electrons.mvaTTH < 0.3,
                        electrons.ak4_jets.btagDeepFlavB <= jetDeepJet_WP_tight[Runyear - 2016],
                        electrons.ak4_jets.btagDeepFlavB <= jetDeepJet_WP_medium[Runyear - 2016]
                )
        ) &
        ak.where(
            electrons.mvaTTH < 0.30,
                (electrons.jetRelIso < 0.7) & (electrons.mvaFall17V2noIso_WP90),
                True
        )
    )

    electron_tight_mask = electrons.mvaTTH >= 0.30


    electrons["preselected"] = ak.where(
        electron_preselection_mask,
            True,
            False
    )

    ele_mu_pair_for_cleaning = ak.cartesian(
        [electrons.mask[electrons.preselected], muons.mask[muons.preselected]], nested=True
    )

    ele_for_cleaning, mu_for_cleaning = ak.unzip( ele_mu_pair_for_cleaning )

    electron_cleaning_dr = ak.where(
        (ak.is_none(mu_for_cleaning, axis=2) == 0) & (ak.is_none(ele_for_cleaning, axis=2) == 0),
            abs(ele_for_cleaning.delta_r(mu_for_cleaning)),
            electrons.preselected
    )

    electron_cleaning_mask = ak.min(electron_cleaning_dr, axis=2) > 0.30

    electrons["cleaned"] = ak.where(
        electron_cleaning_mask & electrons.preselected,
            True,
            False
    )

    electrons["fakeable"] = ak.where(
        electron_fakeable_mask & electrons.cleaned,
            True,
            False
    )


    electrons["tight"] = ak.where(
        electron_tight_mask & electrons.fakeable,
            True,
            False
    )



    #AK4 Jets
    PFJetID = [1, 2, 2]

    ak4_jet_preselection_mask = (
        (
            (ak4_jets.pt > 50.0) | ((ak4_jets.puId & 4) > 0)
        ) &
        (abs(ak4_jets.eta) <= 2.4) & (ak4_jets.pt >= 25.0) &
        ((ak4_jets.jetId & PFJetID[Runyear - 2016]) > 0)
    )


    ak4_jets_loose_btag_mask = ak4_jets.btagDeepFlavB > jetDeepJet_WP_loose[Runyear - 2016]

    ak4_jets_medium_btag_mask = ak4_jets.btagDeepFlavB > jetDeepJet_WP_medium[Runyear - 2016]


    ak4_jets["preselected"] = ak.where(
        ak4_jet_preselection_mask,
            True,
            False
    )

    leptons_fakeable = ak.concatenate([electrons.mask[electrons.fakeable], muons.mask[muons.fakeable]], axis=1)
    leptons_fakeable = leptons_fakeable[ak.argsort(leptons_fakeable.conept, ascending=False)]

    ak4_jet_lep_pair_for_cleaning = ak.cartesian([ak4_jets.mask[ak4_jets.preselected], leptons_fakeable], nested=True)
    ak4_jet_for_cleaning, lep_for_cleaning = ak.unzip( ak4_jet_lep_pair_for_cleaning )

    ak4_jet_cleaning_dr = ak.where(
        (ak.is_none(lep_for_cleaning, axis=2) == 0) & (ak.is_none(ak4_jet_for_cleaning, axis=2) == 0),
            abs(ak4_jet_for_cleaning.delta_r(lep_for_cleaning)),
            ak4_jets.preselected
    )

    ak4_jet_cleaning_mask = ak.min(ak4_jet_cleaning_dr, axis=2) > 0.40

    ak4_jets["cleaned"] = ak.where(
        ak4_jet_cleaning_mask & ak4_jets.preselected,
            True,
            False
    )

    ak4_jets["loose_btag"] = ak.where(
        ak4_jets_loose_btag_mask & ak4_jets.cleaned,
            True,
            False
    )

    ak4_jets["medium_btag"] = ak.where(
        ak4_jets_medium_btag_mask & ak4_jets.cleaned,
            True,
            False
    )



    #AK8 Jets
    ak8_jets.tau2overtau1 = ak.where(
        (ak.is_none(ak8_jets, axis=1) == 0) & (ak.is_none(ak8_jets.tau2, axis=1) == 0) & (ak.is_none(ak8_jets.tau1, axis=1) == 0),
            ak8_jets.tau2 / ak8_jets.tau1,
            10.0
    )

    ak8_jet_preselection_mask = ak.where(
        ak.is_none(ak8_jets, axis=1) == 0,
            (ak.is_none(ak8_jets.subjet1) == 0) & (ak.is_none(ak8_jets.subjet2) == 0) &
            (ak8_jets.subjet1.pt >= 20) & (abs(ak8_jets.subjet1.eta) <= 2.4) &
            (ak8_jets.subjet2.pt >= 20) & (abs(ak8_jets.subjet2.eta) <= 2.4) &
            ( (ak8_jets.jetId & PFJetID[Runyear - 2016]) > 0) & (ak8_jets.pt >= 200) &
            (abs(ak8_jets.eta) <= 2.4) & (ak8_jets.msoftdrop >= 30) & (ak8_jets.msoftdrop <= 210) &
            (ak8_jets.tau2 / ak8_jets.tau1 <= 0.75),
            False
    )

    ak8_btagDeepB_WP_loose  = [0.2217, 0.1522, 0.1241]; ak8_btagDeepB_WP_medium = [0.6321, 0.4941, 0.4184]; ak8_btagDeepB_WP_tight  = [0.8953, 0.8001, 0.7527]

    ak8_jet_btag_mask = (
        (ak8_jets.subjet1.btagDeepB > ak8_btagDeepB_WP_medium[Runyear - 2016]) & (ak8_jets.subjet1.pt >= 30)
        | #OR
        (ak8_jets.subjet2.btagDeepB > ak8_btagDeepB_WP_medium[Runyear - 2016]) & (ak8_jets.subjet2.pt >= 30)
    )

    ak8_jets["preselected"] = ak.where(
        ak8_jet_preselection_mask,
            True,
            False
    )

    ak8_jet_lep_pair_for_cleaning = ak.cartesian([ak8_jets.mask[ak8_jets.preselected], leptons_fakeable], nested=True)
    ak8_jet_for_cleaning, lep_for_cleaning = ak.unzip( ak8_jet_lep_pair_for_cleaning )

    ak8_jet_cleaning_dr = ak.where(
        (ak.is_none(lep_for_cleaning, axis=2) == 0) & (ak.is_none(ak8_jet_for_cleaning, axis=2) == 0),
            abs(ak8_jet_for_cleaning.delta_r(lep_for_cleaning)),
            ak8_jets.preselected
    )

    ak8_jet_cleaning_mask = ak.min(ak8_jet_cleaning_dr, axis=2) > 0.80

    ak8_jets["cleaned"] = ak.where(
        ak8_jet_cleaning_mask & ak8_jets.preselected,
            True,
            False
    )

    ak8_jets["btag"] = ak.where(
        ak8_jet_btag_mask & ak8_jets.cleaned,
            True,
            False
    )

    return



def do_single_lepton_category():
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



    #Tau veto: no tau passing pt>20, abs(eta) < 2.3, abs(dxy) <= 1000, abs(dz) <= 0.2, "decayModeFindingNewDMs", decay modes = {0, 1, 2, 10, 11}, and "byMediumDeepTau2017v2VSjet", "byVLooseDeepTau2017v2VSmu", "byVVVLooseDeepTau2017v2VSe". Taus overlapping with fakeable electrons or fakeable muons within dR < 0.3 are not considered for the tau veto
    #False -> Gets Removed : True -> Passes veto
    tau_veto_pairs = ak.cartesian([taus, leptons_fakeable_sorted])
    taus_for_veto, leps_for_veto = ak.unzip(tau_veto_pairs)

    """
    tau_veto_cleaning = abs(taus_for_veto.delta_r(leps_for_veto)) >= 0.3
    tau_veto_selection = (
    (taus_for_veto.pt > 20) & (abs(taus_for_veto.eta) < 2.3) & (abs(taus_for_veto.dxy) <= 1000.0) & (abs(taus_for_veto.dz) <= 0.2) & (taus_for_veto.idDecayModeNewDMs) &
    (
        (taus_for_veto.decayMode == 0) | (taus_for_veto.decayMode == 1) | (taus_for_veto.decayMode == 2) | (taus_for_veto.decayMode == 10) | (taus_for_veto.decayMode == 11)
    ) &
    (taus_for_veto.idDeepTau2017v2p1VSjet >= 16) & (taus_for_veto.idDeepTau2017v2p1VSmu >= 1) & (taus_for_veto.idDeepTau2017v2p1VSe >= 1)
    )


    tau_veto = ak.any(tau_veto_cleaning & tau_veto_selection, axis=1) == 0
    """


    """
    electron_cleaning_dr = ak.where(
        (ak.is_none(mu_for_cleaning, axis=2) == 0) & (ak.is_none(ele_for_cleaning, axis=2) == 0),
            abs(ele_for_cleaning.delta_r(mu_for_cleaning)),
            electrons.preselected
    )
    electron_cleaning_mask = ak.min(electron_cleaning_dr, axis=2) > 0.30
    """

    print(abs(taus_for_veto.delta_r(leps_for_veto)))
    print("taus eta", taus.eta)
    print("Taus phi", taus.phi)
    print("Leps eta", leptons_fakeable_sorted.eta)
    print("Leps phi", leptons_fakeable_sorted.phi)



    print(abs(taus_for_veto.delta_r(leps_for_veto)))

    tau_veto_cleaning = abs(taus_for_veto.delta_r(leps_for_veto)) >= 0.3


    tau_veto_selection = (
    (taus_for_veto.pt > 20) & (abs(taus_for_veto.eta) < 2.3) & (abs(taus_for_veto.dxy) <= 1000.0) & (abs(taus_for_veto.dz) <= 0.2) & (taus_for_veto.idDecayModeNewDMs) &
    (
        (taus_for_veto.decayMode == 0) | (taus_for_veto.decayMode == 1) | (taus_for_veto.decayMode == 2) | (taus_for_veto.decayMode == 10) | (taus_for_veto.decayMode == 11)
    ) &
    (taus_for_veto.idDeepTau2017v2p1VSjet >= 16) & (taus_for_veto.idDeepTau2017v2p1VSmu >= 1) & (taus_for_veto.idDeepTau2017v2p1VSe >= 1)
    )

    print(tau_veto_selection)

    tau_veto = ak.any(tau_veto_cleaning & tau_veto_selection, axis=1) == 0
    print(tau_veto)




    print("Start the check")
    print(tau_veto_cleaning[14])
    print(ak.any(tau_veto_cleaning, axis=1)[14])
    print(ak.any(tau_veto_selection, axis=1)[14])
    print(tau_veto[14])

    print(events.single_lepton)
    print(tau_veto)

    single_step7_mask = tau_veto

    events["single_lepton"] = ak.where(
        events["single_lepton"] & single_step7_mask,
            True,
            False
    )

    print("N single events step7: ", ak.sum(events.single_lepton))



    """ OLD FOR LOOP VERSION
      for i in self.taus:
        tau_lepton_overlap = False
        if (i.pt > 20.0 and abs(i.eta) < 2.3 and abs(i.dxy) <= 1000.0 and abs(i.dz) <= 0.2 and i.idDecayModeNewDMs and i.decayMode in [0,1,2,10,11] and i.idDeepTau2017v2p1VSjet >= 16 and i.idDeepTau2017v2p1VSmu >= 1 and i.idDeepTau2017v2p1VSe >= 1):
          if self.debug >2 :
            self.printObject(i, "Tau for tau_veto")
          for fake in (fakeable_leptons):
            if self.debug >2 :
              self.printObject(fake, "lepton for tau_veto")
              print("dR ",deltaR(fake.eta, fake.phi, i.eta, i.phi))
            if deltaR(fake.eta, fake.phi, i.eta, i.phi) < 0.3:
              tau_lepton_overlap = True
              break
          if tau_lepton_overlap:
            continue
          else:
            if self.debug >2 :
              print("Veto this event by tau!!")
            return False
      return True
    """


    return




#Object Selection

do_object_selection()
do_single_lepton_category()

print("Events with 1 object comparison (my new coffea value) // (my old nanoAOD value [Tallinn value if different])")
print("Muons preselected: ", ak.sum(ak.any(muons.preselected, axis=1)), " // 93605")
print("Muons fakeable: ", ak.sum(ak.any(muons.fakeable, axis=1)), " // 81978")
print("Muons tight: ", ak.sum(ak.any(muons.tight, axis=1)), " // 78340")


print("Electrons preselected: ", ak.sum(ak.any(electrons.preselected, axis=1)), " // NA")
print("Electrons cleaned: ", ak.sum(ak.any(electrons.cleaned, axis=1)), " // 75430")
print("Electrons fakeable: ", ak.sum(ak.any(electrons.fakeable, axis=1)), " // 58833")
print("Electrons tight: ", ak.sum(ak.any(electrons.tight, axis=1)), " // 56395")

print("AK4 Jets preselected: ", ak.sum(ak.any(ak4_jets.preselected, axis=1)), " // NA")
print("AK4 Jets cleaned: ", ak.sum(ak.any(ak4_jets.cleaned, axis=1)), " // 144403(144446)")
print("AK4 Jets loose Btag: ", ak.sum(ak.any(ak4_jets.loose_btag, axis=1)), " // NA")
print("AK4 Jets medium Btag: ", ak.sum(ak.any(ak4_jets.medium_btag, axis=1)), " // NA")

print("AK8 Jets preselected: ", ak.sum(ak.any(ak8_jets.preselected, axis=1)), " // 77501")
print("AK8 Jets cleaned: ", ak.sum(ak.any(ak8_jets.cleaned, axis=1)), " // 69384")
print("AK8 Jets Btag: ", ak.sum(ak.any(ak8_jets.btag, axis=1)), " // 54065(53678)")

print("N events: ", len(events))
print("N single events: ", ak.sum(events.single_lepton))




#It looks like each branch adds about 10 seconds to execution time
"""
df = ak._v2.to_rdataframe(
    {
        "muon_pt": muons.pt,
        "muon_conept": muons.conept,
        "muon_eta": muons.eta,
        "muon_phi": muons.phi,
    }
)
"""
"""
        "muon_E": muons.pt,
        "muon_charge": muons.charge,
        "muon_miniRelIso": muons.miniPFRelIso_all,
        "muon_PFRelIso04": muons.pfRelIso04_all,
        "muon_jetNDauChargedMVASel": muons.pt,
        "muon_jetPtRel": muons.jetPtRelv2,
        "muon_jetRelIso": muons.jetRelIso,
        "muon_jetDeepJet": muons.ak4_jets.btagDeepFlavB,
        "muon_sip3D": muons.sip3d,
        "muon_dxy": muons.dxy,
        "muon_dxyAbs": abs(muons.dxy),
        "muon_dz": muons.dz,
        "muon_segmentCompatibility": muons.segmentComp,
        "muon_leptonMVA": muons.mvaTTH,
        "muon_mediumID": muons.mediumId,
        "muon_dpt_div_pt": muons.ptErr / muons.pt,
        "muon_preselected": muons.preselected,
        "muon_fakeable": muons.fakeable,
        "muon_tight": muons.tight,
        "muon_isGenMatched": muons.pt,
    }
)
"""

#df.Describe().Print()

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
