import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
import ROOT
from coffea.nanoevents.methods import vector


def add_conept(muons, electrons):
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

def link_jets(muons, electrons, ak4_jets, ak8_jets, ak8_subjets):
    muons.ak4_jets = ak4_jets[muons.jetIdx].mask[(muons.jetIdx >= 0) & (muons.jetIdx < ak.num(ak4_jets))]
    electrons.ak4_jets = ak4_jets[electrons.jetIdx].mask[(electrons.jetIdx >= 0) & (electrons.jetIdx < ak.num(ak4_jets))]
    ak8_jets.subjet1 = ak8_subjets[ak8_jets.subJetIdx1].mask[(ak8_jets.subJetIdx1 >= 0) & (ak8_jets.subJetIdx1 < ak.num(ak8_subjets))]
    ak8_jets.subjet2 = ak8_subjets[ak8_jets.subJetIdx2].mask[(ak8_jets.subJetIdx2 >= 0) & (ak8_jets.subJetIdx2 < ak.num(ak8_subjets))]

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


def muon_selection(muons):
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


def electron_selection(electrons, muons):
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



def ak4_jet_selection(ak4_jets, muons, electrons):
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


def ak8_jet_selection(ak8_jets):
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


def test_selections(muons):
    muons["test"] = True
    return
