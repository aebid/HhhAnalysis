import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
import ROOT
from coffea.nanoevents.methods import vector


def create_df(EventProcess):

    print("Starting create DF")

    muons = EventProcess.muons[ak.argsort(EventProcess.muons.conept, ascending=False)]
    electrons = EventProcess.electrons[ak.argsort(EventProcess.electrons.conept, ascending=False)]
    ak4_jets = EventProcess.ak4_jets[ak.argsort(EventProcess.ak4_jets.pt, ascending=False)]
    ak8_jets = EventProcess.ak8_jets[ak.argsort(EventProcess.ak8_jets.pt, ascending=False)]
    ak8_subjets = EventProcess.ak8_subjets
    events = EventProcess.events
    ls = EventProcess.ls
    run = EventProcess.run

    muon1 = muons[:,0]
    muon2 = muons[:,1]
    electron1 = electrons[:,0]
    electron2 = electrons[:,1]
    ak4_jet1 = ak4_jets[:,0]
    ak4_jet2 = ak4_jets[:,1]
    ak4_jet3 = ak4_jets[:,2]
    ak4_jet4 = ak4_jets[:,3]
    ak8_jets.subjet1 = ak8_subjets[ak8_jets.subJetIdx1].mask[(ak8_jets.subJetIdx1 >= 0)]
    ak8_jets.subjet2 = ak8_subjets[ak8_jets.subJetIdx2].mask[(ak8_jets.subJetIdx2 >= 0)]
    ak8_jet1 = ak8_jets[:,0]
    ak8_jet2 = ak8_jets[:,1]
    ak8_jet1_subjet1 = ak8_jets.subjet1[:,0]
    ak8_jet1_subjet2 = ak8_jets.subjet2[:,0]
    ak8_jet2_subjet1 = ak8_jets.subjet1[:,1]
    ak8_jet2_subjet2 = ak8_jets.subjet2[:,1]

    event_dict = {
        #Event level information
        'event': np.array(events.event),
        'ls': np.array(ls),
        'run': np.array(run),
        'n_presel_muons': np.array(ak.sum(muons.preselected, axis=1), dtype=np.int32),
        'n_fakeable_muons': np.array(ak.sum(muons.fakeable, axis=1), dtype=np.int32),
        'n_tight_muons': np.array(ak.sum(muons.tight, axis=1), dtype=np.int32),

        'n_presel_electrons': np.array(ak.sum(electrons.preselected, axis=1), dtype=np.int32),
        'n_fakeable_electrons': np.array(ak.sum(electrons.fakeable, axis=1), dtype=np.int32),
        'n_tight_electrons': np.array(ak.sum(electrons.tight, axis=1), dtype=np.int32),

        'n_presel_ak4_jets': np.array(ak.sum(ak4_jets.preselected, axis=1), dtype=np.int32),
        'n_cleaned_ak4_jets': np.array(ak.sum(ak4_jets.cleaned_all, axis=1), dtype=np.int32),
        'n_loose_btag_ak4_jets': np.array(ak.sum(ak4_jets.loose_btag_all, axis=1), dtype=np.int32),
        'n_medium_btag_ak4_jets': np.array(ak.sum(ak4_jets.medium_btag_all, axis=1), dtype=np.int32),

        'n_presel_ak8_jets': np.array(ak.sum(ak8_jets.preselected, axis=1), dtype=np.int32),
        'n_cleaned_ak8_jets': np.array(ak.sum(ak8_jets.cleaned_all, axis=1), dtype=np.int32),
        'n_btag_ak8_jets': np.array(ak.sum(ak8_jets.btag_all, axis=1), dtype=np.int32),

        'Single_HbbFat_WjjRes_AllReco': np.array(events.Single_HbbFat_WjjRes_AllReco, dtype=np.int32),
        'Single_HbbFat_WjjRes_MissJet': np.array(events.Single_HbbFat_WjjRes_MissJet, dtype=np.int32),
        'Single_Res_allReco_2b': np.array(events.Single_Res_allReco_2b, dtype=np.int32),
        'Single_Res_allReco_1b': np.array(events.Single_Res_allReco_1b, dtype=np.int32),
        'Single_Res_MissWJet_2b': np.array(events.Single_Res_MissWJet_2b, dtype=np.int32),
        'Single_Res_MissWJet_1b': np.array(events.Single_Res_MissWJet_1b, dtype=np.int32),
        'Single_Signal': np.array(events.Single_Signal, dtype=np.int32),
        'Single_Fake': np.array(events.Single_Fake, dtype=np.int32),
        'single_category_cutflow': np.array(events.single_cutflow, dtype=np.int32),
        'single_is_e': np.array(ak.fill_none(events.is_e, False), dtype=np.int32),
        'single_is_m': np.array(ak.fill_none(events.is_m, False), dtype=np.int32),

        'Double_HbbFat': np.array(events.Double_HbbFat, dtype=np.int32),
        'Double_Res_1b': np.array(events.Double_Res_1b, dtype=np.int32),
        'Double_Res_2b': np.array(events.Double_Res_2b, dtype=np.int32),
        'Double_Signal': np.array(events.Double_Signal, dtype=np.int32),
        'Double_Fake': np.array(events.Double_Fake, dtype=np.int32),
        'double_category_cutflow': np.array(events.double_cutflow, dtype=np.int32),
        'double_is_ee': np.array(ak.fill_none(events.is_ee, False), dtype=np.int32),
        'double_is_mm': np.array(ak.fill_none(events.is_mm, False), dtype=np.int32),
        'double_is_em': np.array(ak.fill_none(events.is_em, False), dtype=np.int32),
    }


    #df.Define('muon1_pt', np.array(ak.fill_none(muon1.pt, -999.0), dtype=np.float32))
    muon1_dict = {
        #Muon information
        'muon1_pt': np.array(ak.fill_none(muon1.pt, -999.0), dtype=np.float32),
        'muon1_conept': np.array(ak.fill_none(muon1.conept, -999.0), dtype=np.float32),
        'muon1_eta': np.array(ak.fill_none(muon1.eta, -999.0), dtype=np.float32),
        'muon1_phi': np.array(ak.fill_none(muon1.phi, -999.0), dtype=np.float32),
        'muon1_E': np.array(ak.fill_none(muon1.E, -999.0), dtype=np.float32),
        'muon1_charge': np.array(ak.fill_none(muon1.charge, -999.0), dtype=np.float32),
        'muon1_miniRelIso': np.array(ak.fill_none(muon1.miniPFRelIso_all, -999.0), dtype=np.float32),
        'muon1_PFRelIso04': np.array(ak.fill_none(muon1.pfRelIso04_all, -999.0), dtype=np.float32),
        #'muon1_jetNDauChargedMVASel':
        'muon1_jetPtRel': np.array(ak.fill_none(muon1.jetPtRelv2, -999.0), dtype=np.float32),
        'muon1_jetRelIso': np.array(ak.fill_none(muon1.jetRelIso, -999.0), dtype=np.float32),
        #'muon1_jetDeepJet': np.array(ak.fill_none(muon1.ak4_jets.btagDeepFlavB * muon1.has_ak4_jet, -999.0), dtype=np.float32),
        'muon1_sip3D': np.array(ak.fill_none(muon1.sip3d, -999.0), dtype=np.float32),
        'muon1_dxy': np.array(ak.fill_none(muon1.dxy, -999.0), dtype=np.float32),
        'muon1_dxyAbs': np.array(ak.fill_none(abs(muon1.dxy), -999.0), dtype=np.float32),
        'muon1_dz': np.array(ak.fill_none(muon1.dz, -999.0), dtype=np.float32),
        'muon1_segmentCompatibility': np.array(ak.fill_none(muon1.segmentComp, -999.0), dtype=np.float32),
        'muon1_leptonMVA': np.array(ak.fill_none(muon1.mvaTTH, -999.0), dtype=np.float32),
        'muon1_mediumID': np.array(ak.fill_none(muon1.mediumId, -999.0), dtype=np.float32),
        'muon1_dpt_dv_pt': np.array(ak.fill_none(muon1.ptErr/muon1.pt, -999.0), dtype=np.float32),
        'muon1_isfakeablesel': np.array(ak.fill_none(muon1.fakeable, -999.0), dtype=np.float32),
        'muon1_ismvasel': np.array(ak.fill_none(muon1.tight, -999.0), dtype=np.float32),
        'muon1_isGenMatched': np.array(ak.fill_none(muon1.MC_Match, -999.0), dtype=np.float32),
    }

    muon2_dict = {
        #Muon information
        'muon2_pt': np.array(ak.fill_none(muon2.pt, -999.0), dtype=np.float32),
        'muon2_conept': np.array(ak.fill_none(muon2.conept, -999.0), dtype=np.float32),
        'muon2_eta': np.array(ak.fill_none(muon2.eta, -999.0), dtype=np.float32),
        'muon2_phi': np.array(ak.fill_none(muon2.phi, -999.0), dtype=np.float32),
        'muon2_E': np.array(ak.fill_none(muon2.E, -999.0), dtype=np.float32),
        'muon2_charge': np.array(ak.fill_none(muon2.charge, -999.0), dtype=np.float32),
        'muon2_miniRelIso': np.array(ak.fill_none(muon2.miniPFRelIso_all, -999.0), dtype=np.float32),
        'muon2_PFRelIso04': np.array(ak.fill_none(muon2.pfRelIso04_all, -999.0), dtype=np.float32),
        #'muon1_jetNDauChargedMVASel':
        'muon2_jetPtRel': np.array(ak.fill_none(muon2.jetPtRelv2, -999.0), dtype=np.float32),
        'muon2_jetRelIso': np.array(ak.fill_none(muon2.jetRelIso, -999.0), dtype=np.float32),
        #'muon1_jetDeepJet': np.array(ak.fill_none(muon1.ak4_jets.btagDeepFlavB * muon1.has_ak4_jet, -999.0), dtype=np.float32),
        'muon2_sip3D': np.array(ak.fill_none(muon2.sip3d, -999.0), dtype=np.float32),
        'muon2_dxy': np.array(ak.fill_none(muon2.dxy, -999.0), dtype=np.float32),
        'muon2_dxyAbs': np.array(ak.fill_none(abs(muon2.dxy), -999.0), dtype=np.float32),
        'muon2_dz': np.array(ak.fill_none(muon2.dz, -999.0), dtype=np.float32),
        'muon2_segmentCompatibility': np.array(ak.fill_none(muon2.segmentComp, -999.0), dtype=np.float32),
        'muon2_leptonMVA': np.array(ak.fill_none(muon2.mvaTTH, -999.0), dtype=np.float32),
        'muon2_mediumID': np.array(ak.fill_none(muon2.mediumId, -999.0), dtype=np.float32),
        'muon2_dpt_dv_pt': np.array(ak.fill_none(muon2.ptErr/muon2.pt, -999.0), dtype=np.float32),
        'muon2_isfakeablesel': np.array(ak.fill_none(muon2.fakeable, -999.0), dtype=np.float32),
        'muon2_ismvasel': np.array(ak.fill_none(muon2.tight, -999.0), dtype=np.float32),
        'muon2_isGenMatched': np.array(ak.fill_none(muon2.MC_Match, -999.0), dtype=np.float32),
    }
    muon_dicts = muon1_dict | muon2_dict


    electron1_dict = {
        #Electron information
        'electron1_pt': np.array(ak.fill_none(electron1.pt, -999.0), dtype=np.float32),
        'electron1_conept': np.array(ak.fill_none(electron1.conept, -999.0), dtype=np.float32),
        'electron1_eta': np.array(ak.fill_none(electron1.eta, -999.0), dtype=np.float32),
        'electron1_phi': np.array(ak.fill_none(electron1.phi, -999.0), dtype=np.float32),
        'electron1_E': np.array(ak.fill_none(electron1.E, -999.0), dtype=np.float32),
        'electron1_charge': np.array(ak.fill_none(electron1.charge, -999.0), dtype=np.float32),
        'electron1_miniRelIso': np.array(ak.fill_none(electron1.miniPFRelIso_all, -999.0), dtype=np.float32),
        'electron1_PFRelIso04': np.array(ak.fill_none(electron1.pfRelIso03_all, -999.0), dtype=np.float32),
        #'electron1_jetNDauChargedMVASel':
        'electron1_jetPtRel': np.array(ak.fill_none(electron1.jetPtRelv2, -999.0), dtype=np.float32),
        'electron1_jetRelIso': np.array(ak.fill_none(electron1.jetRelIso, -999.0), dtype=np.float32),
        #'electron1_jetDeepJet': np.array(ak.fill_none(muon1.ak4_jets.btagDeepFlavB * muon1.has_ak4_jet, -999.0), dtype=np.float32),
        'electron1_sip3D': np.array(ak.fill_none(electron1.sip3d, -999.0), dtype=np.float32),
        'electron1_dxy': np.array(ak.fill_none(electron1.dxy, -999.0), dtype=np.float32),
        'electron1_dxyAbs': np.array(ak.fill_none(abs(electron1.dxy), -999.0), dtype=np.float32),
        'electron1_dz': np.array(ak.fill_none(electron1.dz, -999.0), dtype=np.float32),
        #'electron1_ntMVAeleID':
        'electron1_leptonMVA': np.array(ak.fill_none(electron1.mvaTTH, -999.0), dtype=np.float32),
        'electron1_passesConversionVeto': np.array(ak.fill_none(electron1.convVeto, -999.0), dtype=np.float32),
        'electron1_nMissingHits': np.array(ak.fill_none(electron1.lostHits, -999.0), dtype=np.int32),
        'electron1_sigmaEtaEta': np.array(ak.fill_none(electron1.sieie, -999.0), dtype=np.float32),
        'electron1_HoE': np.array(ak.fill_none(electron1.hoe, -999.0), dtype=np.float32),
        'electron1_OoEminusOoP': np.array(ak.fill_none(electron1.eInvMinusPInv, -999.0), dtype=np.float32),
        'electron1_isfakeablesel': np.array(ak.fill_none(electron1.fakeable, -999.0), dtype=np.float32),
        'electron1_ismvasel': np.array(ak.fill_none(electron1.tight, -999.0), dtype=np.float32),
        'electron1_isGenMatched': np.array(ak.fill_none(electron1.MC_Match, -999.0), dtype=np.float32),
    }


    electron2_dict = {
        #Electron information
        'electron2_pt': np.array(ak.fill_none(electron2.pt, -999.0), dtype=np.float32),
        'electron2_conept': np.array(ak.fill_none(electron2.conept, -999.0), dtype=np.float32),
        'electron2_eta': np.array(ak.fill_none(electron2.eta, -999.0), dtype=np.float32),
        'electron2_phi': np.array(ak.fill_none(electron2.phi, -999.0), dtype=np.float32),
        'electron2_E': np.array(ak.fill_none(electron2.E, -999.0), dtype=np.float32),
        'electron2_charge': np.array(ak.fill_none(electron2.charge, -999.0), dtype=np.float32),
        'electron2_miniRelIso': np.array(ak.fill_none(electron2.miniPFRelIso_all, -999.0), dtype=np.float32),
        'electron2_PFRelIso04': np.array(ak.fill_none(electron2.pfRelIso03_all, -999.0), dtype=np.float32),
        #'electron2_jetNDauChargedMVASel':
        'electron2_jetPtRel': np.array(ak.fill_none(electron2.jetPtRelv2, -999.0), dtype=np.float32),
        'electron2_jetRelIso': np.array(ak.fill_none(electron2.jetRelIso, -999.0), dtype=np.float32),
        #'electron2_jetDeepJet': np.array(ak.fill_none(muon1.ak4_jets.btagDeepFlavB * muon1.has_ak4_jet, -999.0), dtype=np.float32),
        'electron2_sip3D': np.array(ak.fill_none(electron2.sip3d, -999.0), dtype=np.float32),
        'electron2_dxy': np.array(ak.fill_none(electron2.dxy, -999.0), dtype=np.float32),
        'electron2_dxyAbs': np.array(ak.fill_none(abs(electron2.dxy), -999.0), dtype=np.float32),
        'electron2_dz': np.array(ak.fill_none(electron2.dz, -999.0), dtype=np.float32),
        #'electron2_ntMVAeleID':
        'electron2_leptonMVA': np.array(ak.fill_none(electron2.mvaTTH, -999.0), dtype=np.float32),
        'electron2_passesConversionVeto': np.array(ak.fill_none(electron2.convVeto, -999.0), dtype=np.float32),
        'electron2_nMissingHits': np.array(ak.fill_none(electron2.lostHits, -999.0), dtype=np.int32),
        'electron2_sigmaEtaEta': np.array(ak.fill_none(electron2.sieie, -999.0), dtype=np.float32),
        'electron2_HoE': np.array(ak.fill_none(electron2.hoe, -999.0), dtype=np.float32),
        'electron2_OoEminusOoP': np.array(ak.fill_none(electron2.eInvMinusPInv, -999.0), dtype=np.float32),
        'electron2_isfakeablesel': np.array(ak.fill_none(electron2.fakeable, -999.0), dtype=np.float32),
        'electron2_ismvasel': np.array(ak.fill_none(electron2.tight, -999.0), dtype=np.float32),
        'electron2_isGenMatched': np.array(ak.fill_none(electron2.MC_Match, -999.0), dtype=np.float32),
    }
    electron_dicts = electron1_dict | electron2_dict


    ak4_jet1_dict = {
        #AK4 Jet information
        'ak4_jet1_pt': np.array(ak.fill_none(ak4_jet1.pt, -999.0), dtype=np.float32),
        'ak4_jet1_eta': np.array(ak.fill_none(ak4_jet1.eta, -999.0), dtype=np.float32),
        'ak4_jet1_phi': np.array(ak.fill_none(ak4_jet1.phi, -999.0), dtype=np.float32),
        'ak4_jet1_E': np.array(ak.fill_none(ak4_jet1.E, -999.0), dtype=np.float32),
        'ak4_jet1_CSV': np.array(ak.fill_none(ak4_jet1.btagDeepFlavB, -999.0), dtype=np.float32),
        #'ak4_jet1_btagSF': np.array(ak.fill_none(ak4_jet1.btagSF, -999.0), dtype=np.float32),
    }

    ak4_jet2_dict = {
        #AK4 Jet information
        'ak4_jet2_pt': np.array(ak.fill_none(ak4_jet2.pt, -999.0), dtype=np.float32),
        'ak4_jet2_eta': np.array(ak.fill_none(ak4_jet2.eta, -999.0), dtype=np.float32),
        'ak4_jet2_phi': np.array(ak.fill_none(ak4_jet2.phi, -999.0), dtype=np.float32),
        'ak4_jet2_E': np.array(ak.fill_none(ak4_jet2.E, -999.0), dtype=np.float32),
        'ak4_jet2_CSV': np.array(ak.fill_none(ak4_jet2.btagDeepFlavB, -999.0), dtype=np.float32),
        #'ak4_jet2_btagSF': np.array(ak.fill_none(ak4_jet2.btagSF, -999.0), dtype=np.float32),
    }

    ak4_jet3_dict = {
        #AK4 Jet information
        'ak4_jet3_pt': np.array(ak.fill_none(ak4_jet3.pt, -999.0), dtype=np.float32),
        'ak4_jet3_eta': np.array(ak.fill_none(ak4_jet3.eta, -999.0), dtype=np.float32),
        'ak4_jet3_phi': np.array(ak.fill_none(ak4_jet3.phi, -999.0), dtype=np.float32),
        'ak4_jet3_E': np.array(ak.fill_none(ak4_jet3.E, -999.0), dtype=np.float32),
        'ak4_jet3_CSV': np.array(ak.fill_none(ak4_jet3.btagDeepFlavB, -999.0), dtype=np.float32),
        #'ak4_jet3_btagSF': np.array(ak.fill_none(ak4_jet3.btagSF, -999.0), dtype=np.float32),
    }

    ak4_jet4_dict = {
        #AK4 Jet information
        'ak4_jet4_pt': np.array(ak.fill_none(ak4_jet4.pt, -999.0), dtype=np.float32),
        'ak4_jet4_eta': np.array(ak.fill_none(ak4_jet4.eta, -999.0), dtype=np.float32),
        'ak4_jet4_phi': np.array(ak.fill_none(ak4_jet4.phi, -999.0), dtype=np.float32),
        'ak4_jet4_E': np.array(ak.fill_none(ak4_jet4.E, -999.0), dtype=np.float32),
        'ak4_jet4_CSV': np.array(ak.fill_none(ak4_jet4.btagDeepFlavB, -999.0), dtype=np.float32),
        #'ak4_jet4_btagSF': np.array(ak.fill_none(ak4_jet4.btagSF, -999.0), dtype=np.float32),
    }
    ak4_jet_dicts = ak4_jet1_dict | ak4_jet2_dict


    ak8_jet1_dict = {
        #AK8 Jet information
        'ak8_jet1_pt': np.array(ak.fill_none(ak8_jet1.pt, -999.0), dtype=np.float32),
        'ak8_jet1_eta': np.array(ak.fill_none(ak8_jet1.eta, -999.0), dtype=np.float32),
        'ak8_jet1_phi': np.array(ak.fill_none(ak8_jet1.phi, -999.0), dtype=np.float32),
        'ak8_jet1_E': np.array(ak.fill_none(ak8_jet1.E, -999.0), dtype=np.float32),
        'ak8_jet1_msoftdrop': np.array(ak.fill_none(ak8_jet1.msoftdrop, -999.0), dtype=np.float32),
        'ak8_jet1_tau1': np.array(ak.fill_none(ak8_jet1.tau1, -999.0), dtype=np.float32),
        'ak8_jet1_tau2': np.array(ak.fill_none(ak8_jet1.tau2, -999.0), dtype=np.float32),

        'ak8_jet1_subjet1_pt': np.array(ak.fill_none(ak8_jet1_subjet1.pt, -999.0), dtype=np.float32),
        'ak8_jet1_subjet1_eta': np.array(ak.fill_none(ak8_jet1_subjet1.eta, -999.0), dtype=np.float32),
        'ak8_jet1_subjet1_phi': np.array(ak.fill_none(ak8_jet1_subjet1.phi, -999.0), dtype=np.float32),
        'ak8_jet1_subjet1_CSV': np.array(ak.fill_none(ak8_jet1_subjet1.btagDeepB, -999.0), dtype=np.float32),
        'ak8_jet1_subjet2_pt': np.array(ak.fill_none(ak8_jet1_subjet2.pt, -999.0), dtype=np.float32),
        'ak8_jet1_subjet2_eta': np.array(ak.fill_none(ak8_jet1_subjet2.eta, -999.0), dtype=np.float32),
        'ak8_jet1_subjet2_phi': np.array(ak.fill_none(ak8_jet1_subjet2.phi, -999.0), dtype=np.float32),
        'ak8_jet1_subjet2_CSV': np.array(ak.fill_none(ak8_jet1_subjet2.btagDeepB, -999.0), dtype=np.float32),
    }

    ak8_jet2_dict = {
        #AK8 Jet information
        'ak8_jet2_pt': np.array(ak.fill_none(ak8_jet2.pt, -999.0), dtype=np.float32),
        'ak8_jet2_eta': np.array(ak.fill_none(ak8_jet2.eta, -999.0), dtype=np.float32),
        'ak8_jet2_phi': np.array(ak.fill_none(ak8_jet2.phi, -999.0), dtype=np.float32),
        'ak8_jet2_E': np.array(ak.fill_none(ak8_jet2.E, -999.0), dtype=np.float32),
        'ak8_jet2_msoftdrop': np.array(ak.fill_none(ak8_jet2.msoftdrop, -999.0), dtype=np.float32),
        'ak8_jet2_tau1': np.array(ak.fill_none(ak8_jet2.tau1, -999.0), dtype=np.float32),
        'ak8_jet2_tau2': np.array(ak.fill_none(ak8_jet2.tau2, -999.0), dtype=np.float32),

        'ak8_jet2_subjet1_pt': np.array(ak.fill_none(ak8_jet2_subjet1.pt, -999.0), dtype=np.float32),
        'ak8_jet2_subjet1_eta': np.array(ak.fill_none(ak8_jet2_subjet1.eta, -999.0), dtype=np.float32),
        'ak8_jet2_subjet1_phi': np.array(ak.fill_none(ak8_jet2_subjet1.phi, -999.0), dtype=np.float32),
        'ak8_jet2_subjet1_CSV': np.array(ak.fill_none(ak8_jet2_subjet1.btagDeepB, -999.0), dtype=np.float32),
        'ak8_jet2_subjet2_pt': np.array(ak.fill_none(ak8_jet2_subjet2.pt, -999.0), dtype=np.float32),
        'ak8_jet2_subjet2_eta': np.array(ak.fill_none(ak8_jet2_subjet2.eta, -999.0), dtype=np.float32),
        'ak8_jet2_subjet2_phi': np.array(ak.fill_none(ak8_jet2_subjet2.phi, -999.0), dtype=np.float32),
        'ak8_jet2_subjet2_CSV': np.array(ak.fill_none(ak8_jet2_subjet2.btagDeepB, -999.0), dtype=np.float32),
    }
    ak8_jet_dicts = ak8_jet1_dict | ak8_jet2_dict


    df = ROOT.RDF.MakeNumpyDataFrame(event_dict | muon_dicts | electron_dicts | ak4_jet_dicts | ak8_jet_dicts)

    print("Saving DF")

    df.Snapshot('coffea_output', 'my_tree_rdf.root')

    print("Finished DF")
