import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
import ROOT
from coffea.nanoevents.methods import vector


def add_conept(EventProcess):
    muons = EventProcess.muons
    electrons = EventProcess.electrons
    Runyear = EventProcess.Runyear
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

    return


def link_jets(EventProcess):
    muons = EventProcess.muons
    electrons = EventProcess.electrons
    ak4_jets = EventProcess.ak4_jets
    ak8_jets = EventProcess.ak8_jets
    ak8_subjets = EventProcess.ak8_subjets
    Runyear = EventProcess.Runyear
    jetDeepJet_WP_loose = EventProcess.jetDeepJet_WP_loose
    jetDeepJet_WP_medium = EventProcess.jetDeepJet_WP_medium
    jetDeepJet_WP_tight  = EventProcess.jetDeepJet_WP_tight

    #Linking as a muons["ak4_jets"] causes the muon object to return None for some reason? Must do a soft-link with muons.ak4_jets
    muons.ak4_jets = ak4_jets[muons.jetIdx].mask[(muons.jetIdx >= 0)]
    electrons.ak4_jets = ak4_jets[electrons.jetIdx].mask[(electrons.jetIdx >= 0)]
    ak8_jets.subjet1 = ak8_subjets[ak8_jets.subJetIdx1].mask[(ak8_jets.subJetIdx1 >= 0)]
    ak8_jets.subjet2 = ak8_subjets[ak8_jets.subJetIdx2].mask[(ak8_jets.subJetIdx2 >= 0)]

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


def muon_selection(EventProcess):
    #Muons
    muons = EventProcess.muons
    Runyear = EventProcess.Runyear
    jetDeepJet_WP_medium = EventProcess.jetDeepJet_WP_medium
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

    muon_tight_mask = ((muons.mvaTTH >= 0.50) & (muons.mediumId))

    muons["preselected"] = muon_preselection_mask

    muons["fakeable"] = muon_fakeable_mask & muons.preselected

    muons["tight"] = muon_tight_mask & muons.fakeable

    muons["MC_Match"] = (muons.genPartFlav == 1) | (muons.genPartFlav == 15)


def electron_selection(EventProcess):
    #Electrons
    electrons = EventProcess.electrons
    muons = EventProcess.muons
    Runyear = EventProcess.Runyear
    jetDeepJet_WP_medium = EventProcess.jetDeepJet_WP_medium
    jetDeepJet_WP_tight = EventProcess.jetDeepJet_WP_tight

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

    electrons["MC_Match"] = (electrons.genPartFlav == 1) | (electrons.genPartFlav == 15)


def ak4_jet_selection(EventProcess):
    #AK4 Jets
    ak4_jets = EventProcess.ak4_jets
    muons = EventProcess.muons
    electrons = EventProcess.electrons
    Runyear = EventProcess.Runyear
    jetDeepJet_WP_loose = EventProcess.jetDeepJet_WP_loose
    jetDeepJet_WP_medium = EventProcess.jetDeepJet_WP_medium
    PFJetID = EventProcess.PFJetID

    ak4_jet_preselection_mask = (
        (
            (ak4_jets.pt > 50.0) | ((ak4_jets.puId & 4) > 0)
        ) &
        (abs(ak4_jets.eta) <= 2.4) & (ak4_jets.pt >= 25.0) &
        ((ak4_jets.jetId & PFJetID[Runyear - 2016]) > 0)
    )


    ak4_jets_loose_btag_mask = ak4_jets.btagDeepFlavB > jetDeepJet_WP_loose[Runyear - 2016]

    ak4_jets_medium_btag_mask = ak4_jets.btagDeepFlavB > jetDeepJet_WP_medium[Runyear - 2016]

    ak4_jets["preselected"] = ak4_jet_preselection_mask

    leptons_fakeable = ak.concatenate([electrons.mask[electrons.fakeable], muons.mask[muons.fakeable]], axis=1)
    leptons_fakeable = leptons_fakeable[ak.argsort(leptons_fakeable.conept, ascending=False)]

    ak4_jet_lep_pair_for_cleaning = ak.cartesian([ak4_jets.mask[ak4_jets.preselected], leptons_fakeable], nested=True)
    ak4_jet_for_cleaning, lep_for_cleaning = ak.unzip( ak4_jet_lep_pair_for_cleaning )

    ak4_jet_cleaning_dr_all = ak.fill_none(abs(ak4_jet_for_cleaning.delta_r(lep_for_cleaning)), True)
    ak4_jet_cleaning_mask_all = ak.min(ak4_jet_cleaning_dr_all, axis=2) > 0.40

    ak4_jet_cleaning_dr_single = ak.fill_none(abs(ak4_jet_for_cleaning.delta_r(lep_for_cleaning)), True)[:,:,0:1]
    ak4_jet_cleaning_mask_single = ak.min(ak4_jet_cleaning_dr_single, axis=2) > 0.40

    ak4_jet_cleaning_dr_double = ak.fill_none(abs(ak4_jet_for_cleaning.delta_r(lep_for_cleaning)), True)[:,:,0:2]
    ak4_jet_cleaning_mask_double = ak.min(ak4_jet_cleaning_dr_double, axis=2) > 0.40



    ak4_jets["cleaned_all"] = ak4_jet_cleaning_mask_all & ak4_jets.preselected
    ak4_jets["cleaned_single"] = ak4_jet_cleaning_mask_single & ak4_jets.preselected
    ak4_jets["cleaned_double"] = ak4_jet_cleaning_mask_double & ak4_jets.preselected

    ak4_jets["loose_btag_all"] = ak4_jets_loose_btag_mask & ak4_jets.cleaned_all
    ak4_jets["loose_btag_single"] = ak4_jets_loose_btag_mask & ak4_jets.cleaned_single
    ak4_jets["loose_btag_double"] = ak4_jets_loose_btag_mask & ak4_jets.cleaned_double

    ak4_jets["medium_btag_all"] = ak4_jets_medium_btag_mask & ak4_jets.cleaned_all
    ak4_jets["medium_btag_single"] = ak4_jets_medium_btag_mask & ak4_jets.cleaned_single
    ak4_jets["medium_btag_double"] = ak4_jets_medium_btag_mask & ak4_jets.cleaned_double


def ak8_jet_selection(EventProcess):
    #AK8 Jets
    ak8_jets = EventProcess.ak8_jets
    muons = EventProcess.muons
    electrons = EventProcess.electrons
    Runyear = EventProcess.Runyear
    PFJetID = EventProcess.PFJetID

    ak8_jets.tau2overtau1 = ak.where(
        (ak.is_none(ak8_jets, axis=1) == 0) & (ak.is_none(ak8_jets.tau2, axis=1) == 0) & (ak.is_none(ak8_jets.tau1, axis=1) == 0),
            ak8_jets.tau2 / ak8_jets.tau1,
            10.0
    )

    ak8_jet_preselection_mask = (
        (ak.is_none(ak8_jets.subjet1) == 0) & (ak.is_none(ak8_jets.subjet2) == 0) &
        (ak8_jets.subjet1.pt >= 20) & (abs(ak8_jets.subjet1.eta) <= 2.4) &
        (ak8_jets.subjet2.pt >= 20) & (abs(ak8_jets.subjet2.eta) <= 2.4) &
        ( (ak8_jets.jetId & PFJetID[Runyear - 2016]) > 0) & (ak8_jets.pt >= 200) &
        (abs(ak8_jets.eta) <= 2.4) & (ak8_jets.msoftdrop >= 30) & (ak8_jets.msoftdrop <= 210) &
        (ak8_jets.tau2 / ak8_jets.tau1 <= 0.75)
    )

    ak8_btagDeepB_WP_loose  = [0.2217, 0.1522, 0.1241]; ak8_btagDeepB_WP_medium = [0.6321, 0.4941, 0.4184]; ak8_btagDeepB_WP_tight  = [0.8953, 0.8001, 0.7527]

    ak8_jet_btag_mask = (
        (ak8_jets.subjet1.btagDeepB > ak8_btagDeepB_WP_medium[Runyear - 2016]) & (ak8_jets.subjet1.pt >= 30)
        | #OR
        (ak8_jets.subjet2.btagDeepB > ak8_btagDeepB_WP_medium[Runyear - 2016]) & (ak8_jets.subjet2.pt >= 30)
    )

    ak8_jets["preselected"] = ak8_jet_preselection_mask


    leptons_fakeable = ak.concatenate([electrons.mask[electrons.fakeable], muons.mask[muons.fakeable]], axis=1)
    leptons_fakeable = leptons_fakeable[ak.argsort(leptons_fakeable.conept, ascending=False)]

    ak8_jet_lep_pair_for_cleaning = ak.cartesian([ak8_jets.mask[ak8_jets.preselected], leptons_fakeable], nested=True)
    ak8_jet_for_cleaning, lep_for_cleaning = ak.unzip( ak8_jet_lep_pair_for_cleaning )

    ak8_jet_cleaning_dr_all = ak.fill_none(abs(ak8_jet_for_cleaning.delta_r(lep_for_cleaning)), True)
    ak8_jet_cleaning_mask_all = ak.min(ak8_jet_cleaning_dr_all, axis=2) > 0.80

    ak8_jet_cleaning_dr_single = ak8_jet_cleaning_dr_all[:,:,0:1]
    ak8_jet_cleaning_mask_single = ak.min(ak8_jet_cleaning_dr_single, axis=2) > 0.80

    ak8_jet_cleaning_dr_double = ak8_jet_cleaning_dr_all[:,:,0:2]
    ak8_jet_cleaning_mask_double = ak.min(ak8_jet_cleaning_dr_double, axis=2) > 0.80

    ak8_jets["cleaned_all"] = ak8_jet_cleaning_mask_all & ak8_jets.preselected
    ak8_jets["cleaned_single"] = ak8_jet_cleaning_mask_single & ak8_jets.preselected
    ak8_jets["cleaned_double"] = ak8_jet_cleaning_mask_double & ak8_jets.preselected

    ak8_jets["btag_all"] = ak8_jet_btag_mask & ak8_jets.cleaned_all
    ak8_jets["btag_single"] = ak8_jet_btag_mask & ak8_jets.cleaned_single
    ak8_jets["btag_double"] = ak8_jet_btag_mask & ak8_jets.cleaned_double




def test_selections(muons):
    muons["test"] = True
    return
