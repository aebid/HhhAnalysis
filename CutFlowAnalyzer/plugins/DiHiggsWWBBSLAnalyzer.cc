// Package:    DiHiggsWWBBSLAnalyzer
// Class:      DiHiggsWWBBSLAnalyzer
/**\class DiHiggsWWBBSLAnalyzer DiHiggsWWBBSLAnalyzer.cc DiHiggsWW/DiHiggsWWBBSLAnalyzer/plugins/DiHiggsWWBBSLAnalyzer.cc
 */
// Original Author:  tao huang
//         Created:  Wed, 26 Nov 2014 17:58:07 GMT
//
// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//POG
#include "HhhAnalysis/CutFlowAnalyzer/interface/POGRecipesRun2.h"
//headers from root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"

//typedef std::pair<float, float> EtaPhi;
using namespace reco;

//float dxy(const reco::Candidate *cand, const reco::Vertex *point){
//  return (-(cand->vx()-point->x())*cand->py()+(cand->vy()+point->y())*cand->px())/cand->pt();
//};
//float dz(const reco::Candidate *cand, const reco::Vertex *point){


//class MMC;
//#define WMass 80.385   // W mass
//#define SMHMass 125.03 // SM module higgs mass
class DiHiggsWWBBSLAnalyzer : public edm::EDAnalyzer {
  public:
    explicit DiHiggsWWBBSLAnalyzer(const edm::ParameterSet&);
    ~DiHiggsWWBBSLAnalyzer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    edm::ParameterSet cfg_;
    // Labels to access
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;        // reconstructed muons
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;        // reconstructed electron
    edm::EDGetTokenT<pat::JetCollection> jetToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genjetToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjectsToken_;
    //edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    //edm::EDGetTokenT<reco::TrackCollection> trackRefToken_;
    //edm::EDGetTokenT< std::vector<Trajectory> > trajToken_;
    //edm::EDGetTokenT<pat::TriggerEvent> triggerEventToken_;
    //edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
    edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
    bool debug_;
    enum {Data = 0, TTbar, DYJets, DY0Jets, DY1Jets, DY2Jets, ZZTo2L2Q, ZZTo2L2Nu, ZZTo4L, WWToLNuQQ, WWTo2L2Nu, WZTo2L2Q, WZTo1L3Nu, WZTo1L1Nu2Q, WZTo3LNu, ST_tchannel_top, ST_tchannel_antitop, ST_schannel, ST_tW_antitop, ST_tW_top, WJetsToLNu, WJetsToLNu_HT100To200, WJetsToLNu_HT200To400, WJetsToLNu_HT400To600, WJetsToLNu_HT600To800, WJetsToLNu_HT800To1200, WJetsToLNu_HT1200To2500, WJetsToLNu_HT2500ToInf, TTWJetsToQQ, TTWJetsToLNu, TTZToQQ, TTZToLLNuNu};//add other background
    enum {Rad_260=100, Rad_270, Rad_300, Rad_350, Rad_400, Rad_450, Rad_500, Rad_550, Rad_600, Rad_650, Rad_750, Rad_800, Rad_900,  Rad_1000, Grav_260, Grav_270, Grav_300, Grav_350, Grav_400, Grav_450, Grav_500, Grav_550, Grav_600, Grav_650, Grav_750, Grav_800, Grav_900,  Grav_1000};
    //enum {Graviton_260=200, Rad_270, Rad_300, Rad_350, Rad_400, Rad_450, Rad_500, Rad_550, Rad_600, Rad_650, Rad_750, Rad_800, Rad_900,  Rad_1000, Grav_260, Grav_270, Grav_300, Grav_350, Grav_400, Grav_450, Grav_500, Grav_550, Grav_600, Grav_650, Grav_750, Grav_800, Grav_900,  Grav_1000};
    enum {Rad_260_ZZbb=300 };

    int sampleType_;
    //bool runMCMatch;//select physics objects by matching them to gen information
    int jetSelectionAlgo_;
    int muonSelectionAlgo_;
    int jetCorrectionAlgo_;
    int metCorrectionAlgo_;
    float mu_eta_;
    float mu_PFIso_;
    std::string mu_id_;
    float el_eta_;
    float el_Iso_;
    float leadingpt_mumu_;//mu-mu events, leading
    float trailingpt_mumu_;//
    float leadingpt_muel_;//
    float trailingpt_muel_;//
    float leadingpt_elmu_;//
    float trailingpt_elmu_;//
    float leadingpt_elel_;
    float trailingpt_elel_;
    float jet_eta_;
    float jet_leadingpt_;
    float jet_trailingpt_;
    std::string bjetDiscrName_;
    float bjetDiscrCut_loose_;
    float bjetDiscrCut_medium_;
    float bjetDiscrCut_tight_;
    float met_;
    float jetleptonDeltaR_;
    float iterations_;
    //gen matching 
    float leptonsDeltaR_;//dR(gen, reco)
    float jetsDeltaR_;//gen matching 
    bool onlyGenLevel_;
    std::string triggerSFFile_;
    std::string isoSFFile_;
    std::string idSFFile_;
    std::string trackingSFFile_;
    std::string triggerSFhist_;
    std::string isoSFhist_;
    std::string idSFhist_;
    std::string trackingSFhist_;
    std::vector<std::string> hltPaths_;
    bool doTriggerMatching_;
    float deltaPtRel_trigger_;
    float deltaR_trigger_;
    // debuglevel constrol 
    int verbose_; 
    void print();
    void printHtoWWChain();
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  private:
    TH1F *hevent;
    TTree *evtree;
    TFile *output;
    void initBranches();
    void fillbranches(); 
    edm::Service< TFileService > fs;

    //----------branches of tree ---------------------------
    int ievent;
    //******************************
    //GEN Level
    //******************************
  private:
    const reco::Candidate* findLeptonDaughter(const reco::Candidate* );
    const reco::Candidate* findNuDaughter(const reco::Candidate* );
    const reco::Candidate* findLeptonDescendants(const reco::Candidate*, int& );
    const reco::Candidate* findNuDescendants(const reco::Candidate*, int& );
    //decendants and ancestor
    const reco::Candidate* stableDecendant(const reco::Candidate* cand, int id);
    //const reco::Candidate* stablehtoWWdecendant(Particle p, PdgId id);
    const reco::Candidate* findDecendant(const reco::Candidate* cand, int id, bool first=false);
    const reco::Candidate* findAncestor(const reco::Candidate* cand, int id, bool first=false);
    bool hasMother(const reco::Candidate* cand, int id);
    bool hasDaughter(const reco::Candidate* cand, int id);

    void printCandidate(const reco::Candidate* );
    void printallDecendants(const reco::Candidate* );
    void printallAncestors(const reco::Candidate* );
    void checkGenParticlesHbbWW2L2Nu(edm::Handle<reco::GenParticleCollection> genParticleColl);
    void checkGenParticlesHbbWW2J2Nu(edm::Handle<reco::GenParticleCollection> genParticleColl);
    void checkGenParticlesZZbb(edm::Handle<reco::GenParticleCollection> genParticleColl);
    void checkGenParticlesZZbb_2L2J(edm::Handle<reco::GenParticleCollection> genParticleColl);
    void checkGenParticlesTTbar(edm::Handle<reco::GenParticleCollection> genParticleColl);
    void checkGenParticlesDY(edm::Handle<reco::GenParticleCollection> genParticleColl);
    void matchGenJet2Parton(edm::Handle<std::vector<reco::GenJet>> genjetColl);
    void matchmuon2Gen();//match pat:::Muon to gen muon 
    void matchBjets2Gen();//match genjet to gen b and then match pat::Jet to genjet

  private:
    // ---------- Candidates in signal channel ---------------------------
    const reco::Candidate* lep1_W1_cand;
    const reco::Candidate* nu1_W1_cand;
    const reco::Candidate* lep2_W2_cand;
    const reco::Candidate* nu2_W2_cand;
    const reco::Candidate* q1_W2_cand;
    const reco::Candidate* q2_W2_cand;
    const reco::Candidate* lep1_htoWW_cand;
    const reco::Candidate* lep2_htoWW_cand;
    const reco::Candidate* b1_htobb_cand;
    const reco::Candidate* b2_htobb_cand;
    const reco::Candidate* Xtohh_cand;

    const reco::Candidate* lep1cand;
    const reco::Candidate* nu1cand;
    const reco::Candidate* lep2cand;
    const reco::Candidate* nu2cand;
    const reco::Candidate* q1cand;//light quark from W2
    const reco::Candidate* q2cand;//light quark from W2
    const reco::Candidate* b1cand;
    const reco::Candidate* b2cand;
    const reco::Candidate* w1cand;
    const reco::Candidate* w2cand;
    const reco::Candidate* htoWWcand;
    const reco::Candidate* htoBBcand;
    const reco::Candidate* Xtohhcand;
    const reco::Candidate* t1cand;
    const reco::Candidate* t2cand;
    const reco::GenJet* b1genjet;
    const reco::GenJet* b2genjet;
    const reco::GenJet* q1genjet;
    const reco::GenJet* q2genjet;

    //ZZbb
    const reco::Candidate* z1cand;
    const reco::Candidate* z2cand;
    const reco::Candidate* htoZZcand;

    POGRecipesRun2::MuonPOGSFManager* muonPOGSFs;
  private:
    TLorentzVector lep1_lorentz;
    TLorentzVector lep2_lorentz;
    TLorentzVector bbar_lorentz;
    TLorentzVector nu1_lorentz;
    TLorentzVector nu2_lorentz;
    TLorentzVector met_lorentz;
    TLorentzVector bjet_lorentz;
    TLorentzVector bbarjet_lorentz;
    TLorentzVector stableDecendantsLorentz(const reco::Candidate* cand); 
    TLorentzVector calculateMET(); 
    float XsecBr;
    //gen particles
    bool findAllGenParticles;
    float lep1_energy;
    float lep1_pt;
    float lep1_px;
    float lep1_py;
    float lep1_pz;
    float lep1_eta;
    float lep1_phi;
    int lep1_pdgid;
    float w1_energy;
    float w1_pt;
    float w1_eta;
    float w1_phi;
    float w1_px;
    float w1_py;
    float w1_pz;
    float w1_mass;
    float nu1_energy;
    float nu1_pt;
    float nu1_px;
    float nu1_py;
    float nu1_pz;
    float nu1_eta;
    float nu1_phi;
    int nu1_pdgid;
    bool Wtolep1nu1;

    float w2_energy;
    float w2_pt;
    float w2_eta;
    float w2_phi;
    float w2_px;
    float w2_py;
    float w2_pz;
    float w2_mass;
    //w2->q1 q2
    float q1_energy;
    float q1_pt;
    float q1_px;
    float q1_py;
    float q1_pz;
    float q1_eta;
    float q1_phi;
    int   q1_pdgid;
    float q2_energy;
    float q2_pt;
    float q2_px;
    float q2_py;
    float q2_pz;
    float q2_eta;
    float q2_phi;
    int   q2_pdgid;
    float lep2_energy;
    float lep2_pt;
    float lep2_px;
    float lep2_py;
    float lep2_pz;
    float lep2_eta;
    float lep2_phi;
    int   lep2_pdgid;
    float nu2_energy;
    float nu2_pt;
    float nu2_px;
    float nu2_py;
    float nu2_pz;
    float nu2_eta;
    float nu2_phi;
    int   nu2_pdgid;
    bool Wtolep2nu2;
    

    float htoWW_energy;
    float htoWW_pt;
    float htoWW_eta;
    float htoWW_phi;
    float htoWW_px;
    float htoWW_py;
    float htoWW_pz; 
    float htoWW_mass;
    bool htoWW;

    //ZZbb
    float z1_energy;
    float z1_pt;
    float z1_eta;
    float z1_phi;
    float z1_px;
    float z1_py;
    float z1_pz;
    float z1_mass;
    float z2_energy;
    float z2_pt;
    float z2_eta;
    float z2_phi;
    float z2_px;
    float z2_py;
    float z2_pz;
    float z2_mass;
    float htoZZ_energy;
    float htoZZ_pt;
    float htoZZ_eta;
    float htoZZ_phi;
    float htoZZ_px;
    float htoZZ_py;
    float htoZZ_pz; 
    float htoZZ_mass;
    bool htoZZ;

    float b1_energy;
    float b1_pt;
    float b1_eta;
    float b1_phi;
    float b1_px;
    float b1_py;
    float b1_pz;
    float b2_energy;
    float b2_pt;
    float b2_eta;
    float b2_phi;
    float b2_px;
    float b2_py;
    float b2_pz;

    float htobb_energy;
    float htobb_pt;
    float htobb_eta;
    float htobb_phi;
    float htobb_px;
    float htobb_py;
    float htobb_pz;
    float htobb_mass;
    //genjet
    float b1genjet_energy;
    float b1genjet_pt;
    float b1genjet_eta;
    float b1genjet_phi;
    float b1genjet_px;
    float b1genjet_py;
    float b1genjet_pz;
    float b1genjet_mass;
    float dR_b1genjet;
    float b2genjet_energy;
    float b2genjet_pt;
    float b2genjet_eta;
    float b2genjet_phi;
    float b2genjet_px;
    float b2genjet_py;
    float b2genjet_pz;
    float b2genjet_mass;
    float dR_b2genjet;
    bool hastwogenjets;

    float q1genjet_energy;
    float q1genjet_pt;
    float q1genjet_eta;
    float q1genjet_phi;
    float q1genjet_px;
    float q1genjet_py;
    float q1genjet_pz;
    float q1genjet_mass;
    float dR_q1genjet;
    float q2genjet_energy;
    float q2genjet_pt;
    float q2genjet_eta;
    float q2genjet_phi;
    float q2genjet_px;
    float q2genjet_py;
    float q2genjet_pz;
    float q2genjet_mass;
    float dR_q2genjet;
    bool hastwogenjets_wqq;

    float genmet_neutrinos_pt;
    float genmet_neutrinos_phi;
    float genmet_neutrinos_px;
    float genmet_neutrinos_py;
    float genmet_pt;
    float genmet_phi;
    float genmet_px;
    float genmet_py;

    float Xtohh_energy;
    float Xtohh_px;
    float Xtohh_py;
    float Xtohh_pz;
    float Xtohh_mass;
    //for ttbar
    bool ttbar;
    float t1_energy;
    float t1_px;
    float t1_py;
    float t1_pz;
    float t1_mass;
    float t2_energy;
    float t2_px;
    float t2_py;
    float t2_pz;
    float t2_mass;

    float dR_genbl;
    float dR_genb1l1;
    float dR_genb1l2;
    float dR_genb2l1;
    float dR_genb2l2;
    float dR_genl1l2;
    float dR_genl1l2b1b2;
    float dphi_genl1l2b1b2;
    float dR_genb1b2;
    float dR_genminbl;
    float mass_genl1l2;
    float mass_gennu1nu2;
    float energy_genl1l2;
    float pt_genl1l2;
    float phi_genl1l2;
    float eta_genl1l2;
    float mass_genb1b2;
    float energy_genb1b2;
    float pt_genb1b2;
    float phi_genb1b2;
    float eta_genb1b2;
    float dphi_genllbb;
    float dphi_genllmet;
    float mass_gentrans;
    //Wtoqq
    float mass_genq1q2;
    float energy_genq1q2;
    float pt_genq1q2;
    float phi_genq1q2;
    float eta_genq1q2;
    //reco leve
    unsigned int numOfVertices;
    int numberOfmuon1;
    int numberOfmuon2;
    float muon1_energy;
    float muon1_pt;
    float muon1_eta;
    float muon1_phi;
    float muon1_px;
    float muon1_py;
    float muon1_pz;
    float muon1_isoVar;
    float muon1_dxy;
    float muon1_dz;
    float muon1_genpt;
    float muon1_geneta;
    float muon1_genphi;
    //muonid loose, medium, tight
    int muon1_id;
    float dR_lep1;
    float muon1_triggerSF;
    float muon1_isoSF;
    float muon1_idSF;
    float muon1_trackingSF;
    float muon1_pogSF;
    float muon1_trigger_dR;
    float muon1_trigger_dPtRel;
    float muon2_energy;
    float muon2_pt;
    float muon2_eta;
    float muon2_phi;
    float muon2_px;
    float muon2_py;
    float muon2_pz;
    float muon2_isoVar;
    float muon2_dxy;
    float muon2_dz;
    float muon2_genpt;
    float muon2_geneta;
    float muon2_genphi;
    int muon2_id;
    float dR_lep2;
    float muon2_triggerSF;
    float muon2_isoSF;
    float muon2_idSF;
    float muon2_trackingSF;
    float muon2_pogSF;
    float muon2_trigger_dR;
    float muon2_trigger_dPtRel;
    bool hastwomuons;

    float dR_b1jet;
    float dR_b2jet;
    float b1jet_px;
    float b1jet_py;
    float b1jet_pz;
    float b1jet_eta;
    float b1jet_phi;
    float b1jet_pt;
    float b1jet_energy;
    float b1jet_mass;
    unsigned int b1jet_btag;
    float b1jet_bDiscVar;
    float b1jet_mt;
    float b1jet_leadTrackPt;
    float b1jet_leptonPdgId;
    float b1jet_leptonPhi;
    float b1jet_leptonEta;
    float b1jet_leptonPtRel;
    float b1jet_leptonPt;
    float b1jet_leptonDeltaR;
    float b1jet_neHEF;//neutralHardonEfraction
    float b1jet_neEmEF;//neutralEmEnergyFraction
    float b1jet_vtxMass;
    float b1jet_vtxNtracks;
    float b1jet_vtxPx;
    float b1jet_vtxPy;
    float b1jet_vtxPt;
    float b1jet_vtx3DSig;
    float b1jet_vtx3DVal;
    float b1jet_vtxPosX;
    float b1jet_vtxPosY;
    float b1jet_vtxPosZ;
    float b2jet_px;
    float b2jet_py;
    float b2jet_pz;
    float b2jet_eta;
    float b2jet_phi;
    float b2jet_pt;
    float b2jet_energy;
    float b2jet_mass;
    unsigned int b2jet_btag;
    float b2jet_bDiscVar;
    float b2jet_mt;
    float b2jet_leadTrackPt;
    float b2jet_leptonPdgId;
    float b2jet_leptonPhi;
    float b2jet_leptonEta;
    float b2jet_leptonPtRel;
    float b2jet_leptonPt;
    float b2jet_leptonDeltaR;
    float b2jet_neHEF;//neutralHardonEfraction
    float b2jet_neEmEF;//neutralEmEnergyFraction
    float b2jet_vtxMass;
    float b2jet_vtxNtracks;
    float b2jet_vtxPx;
    float b2jet_vtxPy;
    float b2jet_vtxPt;
    float b2jet_vtx3DSig;
    float b2jet_vtx3DVal;
    float b2jet_vtxPosX;
    float b2jet_vtxPosY;
    float b2jet_vtxPosZ;
    bool hastwojets;
    bool hasb1jet;
    bool hasb2jet;

    float met_pt;
    float met_phi;
    float met_px;
    float met_py;

    float HT;

    float minMass;
    float MINdR_bl;
    float dR_bl;
    float dR_b1l1;
    float dR_b1l2;
    float dR_b2l1;
    float dR_b2l2;
    float dR_l1l2;
    float dR_l1l2b1b2;
    float dphi_l1l2b1b2;
    float dR_b1b2;
    float dR_minbl;
    float mass_l1l2;
    float energy_l1l2;
    float pt_l1l2;
    float phi_l1l2;
    float eta_l1l2;
    float mass_b1b2;
    float energy_b1b2;
    float pt_b1b2;
    float phi_b1b2;
    float eta_b1b2;
    float dphi_llbb;
    float dphi_llmet;
    float mass_trans;
  private:
    //bool runMMC_;
    //MMC* thismmc;
    // MMC tree branches
};

DiHiggsWWBBSLAnalyzer::DiHiggsWWBBSLAnalyzer(const edm::ParameterSet& iConfig){
  //****************************************************************************
  //                 SET GEN LEVEL VARIABLES AND MATCHING                      
  //****************************************************************************
  genParticlesToken_    = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  genjetToken_          = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets"));
  jetsDeltaR_ = iConfig.getUntrackedParameter<double>("jetsDeltaR",0.4);
  //leptonsDeltaR_ = iConfig.getUntrackedParameter<double>("leptonsDeltaR",0.1);
  metToken_             = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));

  //****************************************************************************
  //                 SET RECO LEVEL VARIABLES AND COUNTERS                       
  //****************************************************************************
  /*
  muonToken_            = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  electronToken_        = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_             = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));

  triggerBitsToken_     = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"));
  triggerObjectsToken_  = consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerObjects"));
  //beamSpotToken_        = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  //triggerEventToken_    = consumes<pat::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerEvent"));
  //tracksToken_          = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  //trackRefToken_        = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("TrackRefitter"));
  //trajToken_            = consumes< std::vector<Trajectory> >(iConfig.getParameter<edm::InputTag>("Traj"));
  primaryVerticesToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"));

  mu_eta_               = iConfig.getUntrackedParameter<double>("mu_eta",2.4);
  mu_PFIso_             = iConfig.getUntrackedParameter<double>("mu_PFIso", 0.25);
  mu_id_                = iConfig.getUntrackedParameter<std::string>("mu_id","2016Medium");
  el_eta_               = iConfig.getUntrackedParameter<double>("el_eta",2.5);
  el_Iso_               = iConfig.getUntrackedParameter<double>("el_Iso",0.04);
  leadingpt_mumu_       = iConfig.getUntrackedParameter<double>("leadingpt_mumu",20);
  trailingpt_mumu_      = iConfig.getUntrackedParameter<double>("trailingpt_mumu",10);
  leadingpt_muel_       = iConfig.getUntrackedParameter<double>("leadingpt_muel",25);
  trailingpt_muel_      = iConfig.getUntrackedParameter<double>("trailingpt_muel",15);
  leadingpt_elmu_       = iConfig.getUntrackedParameter<double>("leadingpt_elmu",25);
  trailingpt_elmu_      = iConfig.getUntrackedParameter<double>("trailingpt_elmu",10);
  leadingpt_elel_       = iConfig.getUntrackedParameter<double>("leadingpt_elel",25);
  trailingpt_elel_      = iConfig.getUntrackedParameter<double>("trailingpt_elel",15);
  jet_eta_              = iConfig.getUntrackedParameter<double>("jet_eta",2.5);
  jet_leadingpt_        = iConfig.getUntrackedParameter<double>("jet_leadingpt",20);
  jet_trailingpt_       = iConfig.getUntrackedParameter<double>("jet_trailingpt",20);
  bjetDiscrName_        = iConfig.getUntrackedParameter<std::string>("bjetDiscrName","pfCombinedMVAV2BJetTags");
  bjetDiscrCut_loose_   = iConfig.getUntrackedParameter<double>("bjetDiscrCut_loose",0.5);
  bjetDiscrCut_medium_  = iConfig.getUntrackedParameter<double>("bjetDiscrCut_medium",0.7);
  bjetDiscrCut_tight_   = iConfig.getUntrackedParameter<double>("bjetDiscrCut_tight",0.9);
  jetleptonDeltaR_      = iConfig.getUntrackedParameter<double>("jetleptonDeltaR",0.3);
  hltPaths_             = iConfig.getParameter<std::vector<std::string>>("hltPaths");
  doTriggerMatching_         = iConfig.getParameter<bool>("doTriggerMatching");
  deltaPtRel_trigger_   = iConfig.getUntrackedParameter<double>("deltaPtRel_trigger",0.5);
  deltaR_trigger_       = iConfig.getUntrackedParameter<double>("deltaR_trigger",0.1);

  triggerSFFile_        = iConfig.getParameter<std::string>("triggerSFFile");
  isoSFFile_            = iConfig.getParameter<std::string>("isoSFFile");
  idSFFile_             = iConfig.getParameter<std::string>("idSFFile");
  trackingSFFile_       = iConfig.getParameter<std::string>("trackingSFFile");
  triggerSFhist_        = iConfig.getParameter<std::string>("triggerSFhist");
  isoSFhist_            = iConfig.getParameter<std::string>("isoSFhist");
  idSFhist_             = iConfig.getParameter<std::string>("idSFhist");
  trackingSFhist_       = iConfig.getParameter<std::string>("trackingSFhist");
  */
  sampleType_           = iConfig.getUntrackedParameter<int>("SampleType",0);
  onlyGenLevel_         = iConfig.getParameter<bool>("onlyGenLevel");
  debug_                = iConfig.getUntrackedParameter<bool>("debug",false);
  verbose_              = iConfig.getUntrackedParameter<int>("verbose",0);
  //runMMC_               = iConfig.getParameter<bool>("runMMC");
  /*
     iterations_ = iConfig.getUntrackedParameter<int>("iterations",100000);
     seed_ = iConfig.getParameter<int>("seed");
     RefPDFfile_ = iConfig.getParameter<std::string>("RefPDFfile");
     std::cout <<" RefPDFfile_ " << RefPDFfile_ << std::endl;
   */
  // initilize candidates pointer
  //now do what ever initialization is needed
  std::cout  <<"DiHiggs sampletype_ "<< sampleType_ << std::endl;
  ievent = 0;
  lep1_W1_cand = NULL;
  nu1_W1_cand = NULL;
  lep2_W2_cand = NULL;
  nu2_W2_cand = NULL;
  lep1_htoWW_cand = NULL;
  lep2_htoWW_cand = NULL;
  b1_htobb_cand = NULL;
  b2_htobb_cand = NULL;
  Xtohh_cand = NULL;
  /*if (not onlyGenLevel_ and sampleType_ > Data){
      std::vector<std::string> SFTypes;
      std::vector<std::string> SFHists;
      std::vector<std::string> SFFiles;
      SFTypes.push_back("Trigger"); SFTypes.push_back("ISO"); SFTypes.push_back("ID"); SFTypes.push_back("Tracking");
      SFFiles.push_back(triggerSFFile_); SFFiles.push_back(isoSFFile_); SFFiles.push_back(idSFFile_); SFFiles.push_back(trackingSFFile_);
      SFHists.push_back(triggerSFhist_); SFHists.push_back(isoSFhist_); SFHists.push_back(idSFhist_); SFHists.push_back(trackingSFhist_);
      muonPOGSFs = new POGRecipesRun2::MuonPOGSFManager(SFTypes, SFFiles, SFHists);
  }*/
}


void DiHiggsWWBBSLAnalyzer::initBranches(){
  XsecBr=-1;
  findAllGenParticles = false;
  lep1cand = NULL;
  nu1cand = NULL;
  lep2cand = NULL;
  nu2cand = NULL;
  q1cand = NULL;
  q2cand = NULL;
  w1cand =NULL;
  w2cand = NULL;
  z1cand =NULL;
  z2cand = NULL;
  b1cand = NULL;
  b2cand = NULL;
  htoWWcand = NULL;
  htoBBcand = NULL;
  Xtohhcand = NULL;
  t1cand = NULL;
  t2cand = NULL;

  lep1_energy = -999999;
  lep1_pt = -1;
  lep1_eta = -9;
  lep1_phi = -9;
  lep1_px = -999999;
  lep1_py = -999999;
  lep1_pz = -999999;
  lep1_pdgid = -99;
  nu1_energy = -99999;
  nu1_pt = -999999;
  nu1_eta = -999999;
  nu1_phi = -999999;
  nu1_px = -999999;
  nu1_py = -999999;
  nu1_pz = -999999;
  nu1_pdgid = -99;
  w1_energy = -999999;
  w1_pt = -1;
  w1_eta = -9;
  w1_phi = -9;
  w1_px = -999999;
  w1_py = -999999;
  w1_pz = -999999;
  w1_mass = -999999;

  q1_energy = -999999;
  q1_eta =-9;
  q1_phi =-9;
  q1_pt =-1;
  q1_px = -999999;
  q1_py = -999999;
  q1_pz = -999999;
  q1_pdgid = -99;
  q2_energy = -999999;
  q2_eta =-9;
  q2_phi =-9;
  q2_pt =-1;
  q2_px = -999999;
  q2_py = -999999;
  q2_pz = -999999;
  q2_pdgid = -99;
  lep2_energy = -999999;
  lep2_eta =-9;
  lep2_phi =-9;
  lep2_pt =-1;
  lep2_px = -999999;
  lep2_py = -999999;
  lep2_pz = -999999;
  lep2_pdgid = -99;
  nu2_energy = -1;
  nu2_eta =-9;
  nu2_phi =-9;
  nu2_pt =-1;
  nu2_px = -999999;
  nu2_py = -999999;
  nu2_pz = -999999;
  nu2_pdgid = -99;
  w2_energy = -99999;
  w2_eta =-9;
  w2_phi =-9;
  w2_pt =-1;
  w2_px = -999999;
  w2_py = -999999;
  w2_pz = -999999;
  w2_mass = -999999;

  htoWW_energy = -999999;
  htoWW_pt = -1;
  htoWW_px = -999999;
  htoWW_py = -999999;
  htoWW_pz = -999999;
  htoWW_mass = -999999;
  htoWW =false;

  //ZZbb
  z1_energy = -999999;
  z1_pt = -1;
  z1_eta = -9;
  z1_phi = -9;
  z1_px = -999999;
  z1_py = -999999;
  z1_pz = -999999;
  z1_mass = -999999;
  z2_energy = -1;
  z2_pt = -1;
  z2_eta = -9;
  z2_phi = -9;
  z2_px = -999999;
  z2_py = -999999;
  z2_pz = -999999;
  z2_mass = -999999;
  htoZZ_energy = -1;
  htoZZ_pt = -1;
  htoZZ_px = -999999;
  htoZZ_py = -999999;
  htoZZ_pz = -999999;
  htoZZ_mass = -999999;
  htoZZ = false;

  b1_energy = -999999;
  b1_pt = -1;
  b1_eta = -9;
  b1_phi = -9;
  b1_px = -999999.0;
  b1_py = -999999.0;
  b1_pz = -999999.0;
  b2_energy = -999999;
  b2_pt = -1;
  b2_eta = -9;
  b2_phi = -9;
  b2_px = -999999.0;
  b2_py = -999999.0;
  b2_pz = -999999.0;

  htobb_energy = -999999;
  htobb_px = -999999.0;
  htobb_py = -999999.0;
  htobb_pz = -999999.0;
  htobb_mass = -999999.0;

  Xtohh_energy = -999999;
  Xtohh_px = -999999.0;
  Xtohh_py = -999999.0;
  Xtohh_pz = -999999.0;
  Xtohh_mass = -999999.0;

  b1genjet_px = -999999.0;
  b1genjet_py = -999999.0;
  b1genjet_pz = -999999.0;
  b1genjet_eta = -9;
  b1genjet_phi = -9;
  b1genjet_pt  = -1;
  b1genjet_energy=-1;
  b2genjet_px = -999999.0;
  b2genjet_py = -999999.0;
  b2genjet_pz = -999999.0;
  b2genjet_eta=-9;
  b2genjet_phi=-9;
  b2genjet_pt=-1;
  b2genjet_energy=-1;
  dR_b1genjet=jetsDeltaR_;
  dR_b2genjet=jetsDeltaR_;
  hastwogenjets = false;

  q1genjet_px = -999999.0;
  q1genjet_py = -999999.0;
  q1genjet_pz = -999999.0;
  q1genjet_eta = -9;
  q1genjet_phi = -9;
  q1genjet_pt  = -1;
  q1genjet_energy=-1;
  q2genjet_px = -999999.0;
  q2genjet_py = -999999.0;
  q2genjet_pz = -999999.0;
  q2genjet_eta=-9;
  q2genjet_phi=-9;
  q2genjet_pt=-1;
  q2genjet_energy=-1;
  dR_q1genjet=jetsDeltaR_;
  dR_q2genjet=jetsDeltaR_;
  hastwogenjets_wqq = false;

  genmet_neutrinos_pt = -999999.;
  genmet_neutrinos_phi = -999999.;
  genmet_neutrinos_pt = -999999.;
  genmet_neutrinos_phi = -999999.;
  genmet_neutrinos_px = -999999.;
  genmet_neutrinos_py = -999999.;
  genmet_pt = -999999.;
  genmet_phi = -999999.;
  genmet_px = -999999.;
  genmet_py = -999999.;

  //ttar
  ttbar = false;
  t1_px = -999999.0;
  t1_py = -999999.0;
  t1_pz = -999999.0;
  t1_energy = -99999;
  t1_mass = -1;
  t2_px = -999999.0;
  t2_py = -999999.0;
  t2_pz = -999999.0;
  t2_energy = -99999;
  t2_mass = -1;

  dR_genbl=-1.0;
  dR_genb1l1=-1.0;
  dR_genb1l2=-1.0;
  dR_genb2l1=-1.0;
  dR_genb2l2=-1.0;
  dR_genl1l2=-1.0;
  dR_genl1l2b1b2=-1.0;
  dphi_genl1l2b1b2=-1.0;
  dR_genb1b2=-1.0;
  dR_genminbl = -1.0;
  mass_genl1l2 = -1.0;
  mass_gennu1nu2 = -1.0;
  energy_genl1l2 = -1;
  pt_genl1l2 = -1;
  phi_genl1l2 = -1;
  eta_genl1l2 = -1;
  mass_genb1b2 = -1.0;
  energy_genb1b2 = -1;
  pt_genb1b2 = -1;
  phi_genb1b2 = -9;
  eta_genb1b2 = -9;
  dphi_genllbb = -10;
  dphi_genllmet = -10;
  mass_gentrans = -1;

  //wtoqq
  mass_genq1q2 = -1.0;
  energy_genq1q2 = -1;
  pt_genq1q2 = -1;
  phi_genq1q2 = -9;
  eta_genq1q2 = -9;
  //reco level
  numOfVertices = 0;
  numberOfmuon1 = -1;
  numberOfmuon2 = -1;
  muon1_px = -999999.0;
  muon1_py = -999999.0;
  muon1_pz = -999999.0;
  muon1_eta = -9;
  muon1_phi = -9;
  muon1_pt = -1;
  muon1_energy = -1;
  muon1_genpt = -1;
  muon1_geneta = -9.0;
  muon1_genphi = -9.0;
  muon1_isoVar = 10.0;
  muon1_dxy = -999999.0;
  muon1_dz  = -999999.0;
  muon1_triggerSF = 1.0;
  muon1_isoSF = 1.0;
  muon1_idSF = 1.0;
  muon1_trackingSF = 1.0;
  muon1_pogSF = 1.0;
  muon1_trigger_dR = 999.0;
  muon1_trigger_dPtRel = 999.0;
  muon2_px = -999999.0;
  muon2_py = -999999.0;
  muon2_pz = -999999.0;
  muon2_eta = -9;
  muon2_phi = -9;
  muon2_pt = -1;
  muon2_energy = -1;
  muon2_genpt = -1;
  muon2_geneta = -9.0;
  muon2_genphi = -9.0;
  muon2_isoVar = 10.0;
  muon2_dxy = -999999.0;
  muon2_dz  = -999999.0;
  muon2_triggerSF = 1.0;
  muon2_isoSF = 1.0;
  muon2_idSF = 1.0;
  muon2_trackingSF = 1.0;
  muon2_pogSF = 1.0;
  muon2_trigger_dR = 999.0;
  muon2_trigger_dPtRel = 999.0;
  dR_lep1 = 2.0;
  dR_lep2 = 2.0;
  hastwomuons = false;

  dR_b1jet = jetsDeltaR_;
  dR_b2jet = jetsDeltaR_;
  b1jet_px = -999999.0;
  b1jet_py = -999999.0;
  b1jet_pz = -999999.0;
  b1jet_eta=-1;
  b1jet_phi=-9;
  b1jet_pt = -1;
  b1jet_energy=-1;
  b1jet_btag = 0;
  b1jet_bDiscVar = -1;
  b1jet_leadTrackPt = -999999.0;
  b1jet_mt  = -999999.0;
  b1jet_leptonPdgId = 0;
  b1jet_leptonPhi = -9.0;
  b1jet_leptonEta = -9.0;
  b1jet_leptonPtRel = -1.0;
  b1jet_leptonPt = -1.0;
  b1jet_leptonDeltaR = 0;
  b1jet_neHEF = -1.0;//neutralHardonEfraction
  b1jet_neEmEF = -1.0;//neutralEmEnergyFraction
  b1jet_vtxMass = -1.0;
  b1jet_vtxNtracks = -1;
  b1jet_vtxPx = -999999.0;
  b1jet_vtxPy = -999999.0;
  b1jet_vtxPt = -1.0;
  b1jet_vtx3DSig = -999999.0;
  b1jet_vtx3DVal = -999999.0;
  b1jet_vtxPosX  = -999999.0;
  b1jet_vtxPosY = -999999.0;
  b1jet_vtxPosZ = -999999.0;

  b2jet_px = -999999.0;
  b2jet_py = -999999.0;
  b2jet_pz = -999999.0;
  b2jet_eta=-9;
  b2jet_phi=-9;
  b2jet_pt=-1;
  b2jet_energy=-1;
  b2jet_btag = 0;
  b2jet_bDiscVar = -1.;
  b2jet_leadTrackPt = -1.;
  b2jet_mt  = -999999.0;
  b2jet_leptonPdgId = 0;
  b2jet_leptonPhi = -9.0;
  b2jet_leptonEta = -9.0;
  b2jet_leptonPtRel = -1.0;
  b2jet_leptonPt  = -1.0;
  b2jet_leptonDeltaR = -1.0;
  b2jet_neHEF = -1.0;//neutralHardonEfraction
  b2jet_neEmEF = -1.0;//neutralEmEnergyFraction
  b2jet_vtxMass = -1.;
  b2jet_vtxNtracks = 0;
  b2jet_vtxPx = -999999.0;
  b2jet_vtxPy = -999999.0;
  b2jet_vtxPt = -1.0;
  b2jet_vtx3DSig = -1.0;
  b2jet_vtx3DVal = -1.0;
  b2jet_vtxPosX = -999999.0;
  b2jet_vtxPosY = -999999.0;
  b2jet_vtxPosZ = -999999.0;
  hastwojets = false;
  hasb1jet=false;
  hasb2jet=false;

  met_pt = -999.;
  met_phi = -9.;
  met_px = -999999.0;
  met_py = -999999.0;

  HT = 0;

  minMass=999.;
  dR_bl=-1.0;
  dR_b1l1=-1.0;
  dR_b1l2=-1.0;
  dR_b2l1=-1.0;
  dR_b2l2=-1.0;
  MINdR_bl=-1.0;
  dR_l1l2=-1.0;
  dR_l1l2b1b2=-1.0;
  dphi_l1l2b1b2=-1.0;
  dR_b1b2=-1.0;
  dR_minbl = -1.0;
  dR_genl1l2 = -1;
  mass_l1l2 = -1.0;
  energy_l1l2 = -1;
  pt_l1l2 = -1;
  phi_l1l2 = -1;
  eta_l1l2 = -1;
  mass_b1b2 = -1.0;
  energy_b1b2 = -1;
  pt_b1b2 = -1;
  phi_b1b2 = -1;
  eta_b1b2 = -1;
  dphi_llbb = -10;
  dphi_llmet = -10;
  mass_trans = -1;
}

DiHiggsWWBBSLAnalyzer::~DiHiggsWWBBSLAnalyzer(){
}

//step 1 find gen particles if it is MC sample
//step 2 select pat objects: matched by gen or not, and apply clear-up cuts at the same time 
//step 3 run MMC on selected objects: gen level or reco level
void DiHiggsWWBBSLAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  initBranches(); 
  ievent++;
  if (debug_) 
      std::cout << "event  " << iEvent.id().event()<<" ievent "<< ievent << std::endl;
  //Compute weight
  float BR_h_bb   = 0.577;
  float BR_h_WW   = 0.215;
  float BR_W_lnu  = 0.3257;
  float BR_final = pow(BR_h_bb,2) + pow(BR_h_WW,2)*pow(BR_W_lnu,4) + 2*BR_h_bb*BR_h_WW*pow(BR_W_lnu,2);
  XsecBr=-1;
  if(sampleType_==TTbar)                    XsecBr = 87.31; 
  else if(sampleType_==DYJets)                   XsecBr = 18610;
  else if(sampleType_==DY0Jets)                  XsecBr = 4758.9;
  else if(sampleType_==DY1Jets)                  XsecBr = 929.1;
  else if(sampleType_==DY2Jets)                  XsecBr = 337.1;
  else if(sampleType_==ZZTo2L2Q)                 XsecBr = 3.22;
  else if(sampleType_==ZZTo2L2Nu)                XsecBr = 0.564;
  else if(sampleType_==ZZTo4L)                   XsecBr = 1.256;
  else if(sampleType_==WWToLNuQQ)                XsecBr = 49.997;
  else if(sampleType_==WWTo2L2Nu)                XsecBr = 12.178;
  else if(sampleType_==WZTo2L2Q)                 XsecBr = 5.595;
  else if(sampleType_==WZTo1L3Nu)                XsecBr = 3.033;
  else if(sampleType_==WZTo1L1Nu2Q)              XsecBr = 10.71;
  else if(sampleType_==WZTo3LNu)                 XsecBr = 4.42965;
  else if(sampleType_==ST_tchannel_top)          XsecBr = 136.02;
  else if(sampleType_==ST_tchannel_antitop)      XsecBr = 80.95;
  else if(sampleType_==ST_schannel)              XsecBr = 3.36;
  else if(sampleType_==ST_tW_antitop)            XsecBr = 19.5545;
  else if(sampleType_==ST_tW_top)                XsecBr = 19.5545;
  else if(sampleType_==WJetsToLNu)               XsecBr = 61526.7;
  else if(sampleType_==WJetsToLNu_HT100To200)    XsecBr = 1627.45;
  else if(sampleType_==WJetsToLNu_HT200To400)    XsecBr = 435.237;
  else if(sampleType_==WJetsToLNu_HT400To600)    XsecBr = 59.181;
  else if(sampleType_==WJetsToLNu_HT600To800)    XsecBr = 14.580;
  else if(sampleType_==WJetsToLNu_HT800To1200)   XsecBr = 6.656;
  else if(sampleType_==WJetsToLNu_HT1200To2500)  XsecBr = 1.608;
  else if(sampleType_==WJetsToLNu_HT2500ToInf)   XsecBr = 0.0389;
  else if(sampleType_==TTWJetsToQQ)              XsecBr = 0.4062;
  else if(sampleType_==TTWJetsToLNu)             XsecBr = 0.2043;
  else if(sampleType_==TTZToQQ)                  XsecBr = 0.5297;
  else if(sampleType_==TTZToLLNuNu)              XsecBr = 0.2529;
  else if(sampleType_>= Rad_260)                  XsecBr = 1.;// To be added

  //****************************************************************************
  //                GENERATOR LEVEL                       
  //****************************************************************************
  edm::Handle<reco::GenParticleCollection> genParticleColl;
  edm::Handle<std::vector<reco::GenJet>> genjetColl;
  try{
      iEvent.getByToken(genParticlesToken_, genParticleColl);
      iEvent.getByToken(genjetToken_, genjetColl);
  } catch (...){
  	std::cout <<"no Gen information "<< std::endl;
  }

  if (sampleType_>=Rad_260) {
    //checkGenParticlesHbbWW2L2Nu(genParticleColl);
    checkGenParticlesHbbWW2J2Nu(genParticleColl);
    //if (not findAllGenParticles)
	  //checkGenParticlesZZbb(genParticleColl);
    //if (findAllGenParticles) std::cout <<"find all gen particles "<< std::endl;
    //else std::cout <<"Not all gen particles "<< std::endl;
  }
  //else if (sampleType_>=Rad_260_ZZbb)                                          checkGenParticlesZZbb_2L2J(genParticleColl);
  ///else if (sampleType_>=Rad_260_ZZbb)                                          checkGenParticlesZZbb(genParticleColl);
  else if (sampleType_==TTbar)                                          checkGenParticlesTTbar(genParticleColl);
  else if (sampleType_>=DYJets and sampleType_<=DY2Jets)                checkGenParticlesDY(genParticleColl);
  if (sampleType_>Data and findAllGenParticles)  matchGenJet2Parton( genjetColl );

  if (findAllGenParticles) {
      fillbranches(); //fill Gen info into tree
      std::cout <<"find all gen particles "<< std::endl;
  }
 

  if (onlyGenLevel_ and sampleType_>Data){
        if (findAllGenParticles)
	    evtree->Fill();
	return;
  }
  
  //std::cout <<"reco level "<< std::endl;
  //****************************************************************************
  //                RECO LEVEL
  //****************************************************************************
  /*
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByToken(primaryVerticesToken_, primaryVertices);
  if (primaryVertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = primaryVertices->front();
  //std::cout <<"number of vertices "<< primaryVertices->size() << std::endl;
  numOfVertices = primaryVertices->size() ;
   */
  //****************************************************************************
  //                Triggering matching 
  //                1. find all trigger objects matched to hlt paths
  //                2. match reco objects to trigger objects 
  //****************************************************************************
  /*
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  std::vector<pat::TriggerObjectStandAlone> matchedTriggerObjects;
  //edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  try{
      iEvent.getByToken(triggerBitsToken_, triggerBits);
      iEvent.getByToken(triggerObjectsToken_, triggerObjects);
      //iEvent.getByToken(triggerPrescales_, triggerPrescales);
  } catch (...){
      std::cout <<"no Trigger results information "<< std::endl;
  }

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //std::cout << "\nTRIGGER NAMEs, size " << names.size() << " TRIGGER OBJECTS, size "<< (*triggerObjects).size() << std::endl;
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);
      //std::cout << "\tTrigger object:  id "<< obj.pdgId() << " pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      // Print trigger object collection and type
      //std::cout << "\t   Collection: " << obj.collection() << std::endl;
      for (auto path : hltPaths_)
      	if (obj.hasPathName(path)) {
	    matchedTriggerObjects.push_back(obj);
	    break;
	}
  }
  
  std::cout <<"matched trigger objects, size  "<< matchedTriggerObjects.size() << std::endl; 
  for (auto obj : matchedTriggerObjects){
      std::cout << "\tTrigger object:  id "<< obj.pdgId() << " pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      // Print trigger object collection and type
      std::cout << "\t   Collection: " << obj.collection() << std::endl;
  }
  if (doTriggerMatching_ and matchedTriggerObjects.size() <= 1){
      std::cout <<"doTriggerMatching_ but size of matched tigger objects is " << matchedTriggerObjects.size() << std::endl;
      std::cout <<"Skip this Event !!! "<< std::endl;
      return ;
  }*/

  //****************************************************************************
  //                MET
  //****************************************************************************

  //std::cout <<"MET "<< std::endl;
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
  if (sampleType_ > Data){
      const reco::GenMET *genmet = met.genMET();
      genmet_px = genmet->px(); genmet_py = genmet->py(); genmet_phi = genmet->phi(); genmet_pt = genmet->pt();
  }
  printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). MET with JES up/down: %.1f/%.1f\n", 
	  met.pt(), met.phi(), met.sumEt(),met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));
  met_px = met.px(); met_py = met.py(); met_phi = met.phi(); met_pt = met.pt();
  //****************************************************************************
  //                Di-Leptons selection
  //                Medium ID + loose PF-based combined relative isolation with detaB correction
  //****************************************************************************
  //
  /*
  edm::Handle<pat::MuonCollection> muons;
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(muonToken_, muons);
  iEvent.getByToken(electronToken_, electrons);
  std::vector<const reco::Candidate *> pleptons;
  std::vector<const reco::Candidate *> nleptons;
  //std::cout <<"diMuon "<< std::endl;
  for (const pat::Muon &mu : *muons) {
    if (doTriggerMatching_){
	bool triggermatching = false;
	for (auto obj : matchedTriggerObjects){
	    float dr = deltaR(obj.eta(), obj.phi(), mu.eta(), mu.phi());
	    float dPtRel = std::fabs(obj.pt()-mu.pt())/mu.pt();
	    if (obj.id(trigger::TriggerMuon) and dr < deltaR_trigger_ and dPtRel < deltaPtRel_trigger_){
		    triggermatching = true;
		    break;
	    }
	}
	if (not triggermatching) continue;
    }

    bool muonid = false;
    if (mu_id_ == "2016Medium") muonid = POGRecipesRun2::is2016MediumMuon(mu);
    else if (mu_id_ == "Medium") muonid = POGRecipesRun2::isMediumMuon(mu);
    else std::cout <<"mu ID is not correct: "<< mu_id_ << std::endl;
    float isoVar = POGRecipesRun2::MuonIsoPFbased(mu);
    if (fabs(mu.eta())<mu_eta_ and mu.pt()>10 and muonid and isoVar <= mu_PFIso_){ 
	//and fabs(mu.muonBestTrack()->dz(PV.position()))<0.1 and 
	//((mu.pt()>20 and fabs(mu.muonBestTrack()->dxy(PV.position()))<0.02) or 
	//(mu.pt()<20 and fabs(mu.muonBestTrack()->dxy(PV.position()))<0.01))){
	if (mu.charge()>0) pleptons.push_back(&mu);
	else if (mu.charge()<0) nleptons.push_back(&mu);
	if (debug_)
	    std::cout <<"get one muon passed selection eta "<< mu.eta() <<" pt "<< mu.pt()<<" isovar "<<  isoVar <<" charge "<< mu.charge()<< std::endl;
	const reco::GenParticle * genp = mu.genParticle();
	if (genp and debug_)
	  std::cout <<"matched genParticle: id "<< genp->pdgId()<<" px "<< genp->px() <<" py "<< genp->py()<<" pz "<< genp->pz() << std::endl;
    }
    if(debug_) printf("muon with pt %4.1f, charge %d, IsoVar %5.3f, dz(PV) %+5.3f, dxy(PV)%+5.3f, POG loose id %d, tight id %d\n",
	  mu.pt(), mu.charge(), isoVar, mu.muonBestTrack()->dz(PV.position()), fabs(mu.muonBestTrack()->dxy(PV.position())),mu.isLooseMuon(), mu.isTightMuon(PV));
  }

  //Ignore electron now
  for (const pat::Electron &el : *electrons) {
    if (el.pt()>1000000) continue;
  }


  const reco::Candidate * selectedPlep = NULL;
  const reco::Candidate * selectedNlep = NULL;
  const reco::Candidate * selectedleadinglep = NULL;
  const reco::Candidate * selectedsubleadinglep = NULL;
  TLorentzVector dilep_p4;
  float sumPt=0.0;
  for (const reco::Candidate *plep : pleptons) {
    for (const reco::Candidate *nlep : nleptons){
	//select lepton pairs with larger sumPt
	dilep_p4.SetPxPyPzE(plep->px()+nlep->px(), plep->py()+nlep->py(), plep->pz()+nlep->pz(), plep->energy()+nlep->energy());
	if(((plep->pt()>leadingpt_mumu_ and nlep->pt()>trailingpt_mumu_) or (nlep->pt()>leadingpt_mumu_ and plep->pt()>trailingpt_mumu_))and dilep_p4.M()>12 and (plep->pt()+nlep->pt())>sumPt){
	  selectedPlep = plep;
	  selectedNlep = nlep;
	  sumPt = plep->pt()+nlep->pt();
	  if (plep->pt() > nlep->pt()){
	      selectedleadinglep = plep;
	      selectedsubleadinglep = nlep;
	  }else {
	      selectedleadinglep = nlep;
	      selectedsubleadinglep = plep;
	  }
	}
    }
  }

  //bool hastwomuons = false;
  //two leptons and adding SFs
  if (sumPt>=30){
    muon1_px = selectedleadinglep->px(); muon1_py = selectedleadinglep->py(); muon1_pz = selectedleadinglep->pz(); 
    muon1_energy = selectedleadinglep->energy();
    muon1_pt = selectedleadinglep->pt(); muon1_eta = selectedleadinglep->eta(); muon1_phi = selectedleadinglep->phi();
    muon1_dxy = dxy(selectedleadinglep, &PV); muon1_dz = dz(selectedleadinglep, &PV);
    muon2_px = selectedsubleadinglep->px(); muon2_py = selectedsubleadinglep->py(); muon2_pz = selectedsubleadinglep->pz(); 
    muon2_energy = selectedsubleadinglep->energy();
    muon2_pt = selectedsubleadinglep->pt(); muon2_eta = selectedsubleadinglep->eta(); muon2_phi = selectedsubleadinglep->phi();
    muon2_dxy = dxy(selectedsubleadinglep, &PV); muon2_dz = dz(selectedsubleadinglep, &PV);
    for (auto obj : matchedTriggerObjects){
	float dr1 = deltaR(obj.eta(), obj.phi(), muon1_eta, muon1_phi);
	float dPtRel1 = std::fabs(obj.pt()-muon1_pt)/muon1_pt;
	if (obj.id(trigger::TriggerMuon) and dr1 < muon1_trigger_dR and dPtRel1 < muon1_trigger_dPtRel){
	    muon1_trigger_dR = dr1;
	    muon1_trigger_dPtRel = dPtRel1;
	}
	float dr2 = deltaR(obj.eta(), obj.phi(), muon2_eta, muon2_phi);
	float dPtRel2 = std::fabs(obj.pt()-muon2_pt)/muon2_pt;
	if (obj.id(trigger::TriggerMuon) and dr2 < muon2_trigger_dR and dPtRel2 < muon2_trigger_dPtRel){
	    muon2_trigger_dR = dr2;
	    muon2_trigger_dPtRel = dPtRel2;
	}
    }

    hastwomuons = true;
    if (sampleType_>Data){
	float dR_gen_reco11 = deltaR(lep1_eta, lep1_phi, muon1_eta, muon1_phi);
	float dR_gen_reco12 = deltaR(lep1_eta, lep1_phi, muon2_eta, muon2_phi);
	float dR_gen_reco21 = deltaR(lep2_eta, lep2_phi, muon1_eta, muon1_phi);
	float dR_gen_reco22 = deltaR(lep2_eta, lep2_phi, muon2_eta, muon2_phi);
	dR_lep1 = (dR_gen_reco11 < dR_gen_reco21) ? dR_gen_reco11:dR_gen_reco21;
	dR_lep2 = (dR_gen_reco12 < dR_gen_reco22) ? dR_gen_reco12:dR_gen_reco22;
	//apply SF for MC samplpes 
	//float triggerSF1 =  POGRecipesRun2::getMuonTriggerSF(std::abs(muon1_eta), muon1_pt, triggerSFFile_, triggerSFhist_);
	//float isoSF1 =  POGRecipesRun2::getMuonISOSF(std::abs(muon1_eta), muon1_pt, isoSFFile_, isoSFhist_);
	//float idSF1 =  POGRecipesRun2::getMuonIDSF(std::abs(muon1_eta), muon1_pt, idSFFile_, idSFhist_);
	//float trackingSF1 = POGRecipesRun2::getMuonTrackingSF(muon1_eta, trackingSFFile_, trackingSFhist_);
	float triggerSF1 = muonPOGSFs->getMuonTriggerMCSF(std::abs(muon1_eta), muon1_pt);
	float isoSF1 = muonPOGSFs->getMuonISOMCSF(std::abs(muon1_eta), muon1_pt);
	float idSF1 = muonPOGSFs->getMuonIDMCSF(std::abs(muon1_eta), muon1_pt);
	float trackingSF1 = muonPOGSFs->getMuonTrackingMCSF(muon1_eta);
	muon1_pogSF = triggerSF1*isoSF1*idSF1*trackingSF1;
	muon1_triggerSF = triggerSF1; muon1_isoSF = isoSF1; muon1_idSF = idSF1; muon1_trackingSF = trackingSF1;
	//float triggerSF2 =  POGRecipesRun2::getMuonTriggerSF(std::abs(muon2_eta), muon2_pt, triggerSFFile_, triggerSFhist_);
	//float isoSF2 =  POGRecipesRun2::getMuonISOSF(std::abs(muon2_eta), muon2_pt, isoSFFile_, isoSFhist_);
	//float idSF2 =  POGRecipesRun2::getMuonIDSF(std::abs(muon2_eta), muon2_pt, idSFFile_, idSFhist_);
	//float trackingSF2 = POGRecipesRun2::getMuonTrackingSF(muon2_eta, trackingSFFile_, trackingSFhist_);
	float triggerSF2 = muonPOGSFs->getMuonTriggerMCSF(std::abs(muon2_eta), muon2_pt);
	float isoSF2 = muonPOGSFs->getMuonISOMCSF(std::abs(muon2_eta), muon2_pt);
	float idSF2 = muonPOGSFs->getMuonIDMCSF(std::abs(muon2_eta), muon2_pt);
	float trackingSF2 = muonPOGSFs->getMuonTrackingMCSF(muon2_eta);
	muon2_pogSF = triggerSF2*isoSF2*idSF2*trackingSF2;
	muon2_triggerSF = triggerSF2; muon2_isoSF = isoSF2; muon2_idSF = idSF2; muon2_trackingSF = trackingSF2;

	printf("leadinglepton: pt %5.1f, eta %+4.2f, triggerSF %+2.4f, isoSF %+2.4f, idSF %+2.4f, trackingSF %+2.4f\n", muon1_pt, muon1_eta, triggerSF1, isoSF1, idSF1, trackingSF1);
	printf("subleadlepton: pt %5.1f, eta %+4.2f, triggerSF %+2.4f, isoSF %+2.4f, idSF %+2.4f, trackingSF %+2.4f\n", muon2_pt, muon2_eta, triggerSF2, isoSF2, idSF2, trackingSF2);
    }else {
	printf("leadinglepton: pt %5.1f, eta %+4.2f, pogSF %+2.4f\n", muon1_pt, muon1_eta, muon1_pogSF);
	printf("subleadlepton: pt %5.1f, eta %+4.2f, pogSF %+2.4f\n", muon2_pt, muon2_eta, muon2_pogSF);
    }
  }
  */

  //****************************************************************************
  //                Di-Jets selection
  //****************************************************************************
  /*
  //std::cout <<"diJet "<< std::endl;
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  std::vector<pat::Jet> allbjets;
  for (const pat::Jet &j : *jets) {
    HT +=  j.pt();//scalar sum of all jet pt 
    if (j.pt() < jet_trailingpt_ or fabs(j.eta()) > jet_eta_) continue;
    if (hastwomuons){
	TLorentzVector jet_p4(j.px(), j.py(), j.pz(), j.energy());    
	float dR1 = jet_p4.DeltaR(TLorentzVector(selectedPlep->px(), selectedPlep->py(), selectedPlep->pz(), selectedPlep->energy()));
	float dR2 = jet_p4.DeltaR(TLorentzVector(selectedNlep->px(), selectedNlep->py(), selectedNlep->pz(), selectedNlep->energy()));
	if (dR1 < jetleptonDeltaR_ or dR2 < jetleptonDeltaR_) continue;
    }	
    float bDiscVar = j.bDiscriminator(bjetDiscrName_);
    if (bDiscVar < bjetDiscrCut_loose_)  continue;
    allbjets.push_back(j);

    if (debug_){
	printf("Jet with pt %6.1f, eta %+4.2f, pileup mva disc %+.2f, btag CSV %.3f, CISV %.3f\n",
	  j.pt(),j.eta(), j.userFloat("pileupJetId:fullDiscriminant"), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")));
	//printf("Jet with virtex: vtxMass %+4.2f, vtxNtracks %.1f, vtxPt %+4.2f, vtx3DSig %+4.2f, vtx3DVal %+4.2f, vtxPosX %+4.2f, vtxPosY %+4.2f, vtxPosZ %+4.2f", j.userFloat("vtxMass"), j.userFloat("vtxNtracks"), sqrt(j.userFloat("vtxPx")*j.userFloat("vtxPx") + j.userFloat("vtxPy")*j.userFloat("vtxPy")), j.userFloat("vtx3DSig"), j.userFloat("vtx3DVal"), j.userFloat("vtxPosX"), j.userFloat("vtxPosY"), j.userFloat("vtxPosZ"));
    }
    const reco::GenParticle * genp = j.genParticle();
    if (genp and debug_)
	std::cout <<"matched genParticle: id "<< genp->pdgId()<<" px "<< genp->px() <<" py "<< genp->py()<<" pz "<< genp->pz() << std::endl;
  }

  // sort jets by pt
  std::sort(allbjets.begin(), allbjets.end(), [](pat::Jet& jet1, pat::Jet& jet2) { return jet1.pt() > jet2.pt(); });
  if (debug_) std::cout <<"allbjets size "<< allbjets.size() << std::endl;
  unsigned int jet1=0, jet2=0;
  int numOfMediumbtags  = 1;
  float diff_higgsmass = 9999;
  for (unsigned int i=0; i<allbjets.size(); i++){
    for (unsigned int j=i+1; j<allbjets.size(); j++){
	int mbtags = 0;
	float bDiscVar1 = allbjets[i].bDiscriminator(bjetDiscrName_);
	float bDiscVar2 = allbjets[j].bDiscriminator(bjetDiscrName_);
	 if (bDiscVar1 > bjetDiscrCut_medium_) mbtags++;
	if (bDiscVar2 > bjetDiscrCut_medium_) mbtags++;
	if (mbtags >=  1)
	  hastwojets = true;
	else continue;

	if (mbtags > numOfMediumbtags){//first priority: 2 medium btags
	  jet1 = i;
	  jet2 = j;
	  numOfMediumbtags = mbtags;
	}else if (mbtags == numOfMediumbtags){//second priority: invariant mass close to M_H
	  TLorentzVector dijet_p4(allbjets[i].px()+allbjets[j].px(), allbjets[i].py()+allbjets[j].py(), 
		allbjets[i].pz()+allbjets[j].pz(),allbjets[i].energy()+allbjets[j].energy());
	  if (fabs(dijet_p4.M()-125)<diff_higgsmass){
	    jet1 = i;
	    jet2 = j;
	    diff_higgsmass = fabs(dijet_p4.M()-125); 
	  }
	}
    }
  }


  if (allbjets.size() >= 2 and hastwojets){
    b1jet_px = allbjets[jet1].px(); b1jet_py = allbjets[jet1].py(); b1jet_pz = allbjets[jet1].pz(); b1jet_energy = allbjets[jet1].energy();
    b1jet_pt = allbjets[jet1].pt(); b1jet_eta = allbjets[jet1].eta(); b1jet_phi = allbjets[jet1].phi();
    b1jet_bDiscVar = allbjets[jet1].bDiscriminator(bjetDiscrName_);
    b1jet_mt  = allbjets[jet1].mt();//ptcorrection?
    auto daus1(allbjets[jet1].daughterPtrVector());
    std::sort(daus1.begin(), daus1.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // C++11
    bool leadinglepton1 = false;
    float leadinglepton1_px = 0; float leadinglepton1_py = 0; float leadinglepton1_pz = 0; float leadinglepton1_p = 0;
    for (unsigned int i = 0; i < daus1.size(); ++i) {
	const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus1[i]);
	if (i==0)
	  b1jet_leadTrackPt = cand.pt();
	if ((abs(cand.pdgId()) == 13 or abs(cand.pdgId()) == 11) and not(leadinglepton1)){
	  b1jet_leptonPdgId = cand.pdgId(); b1jet_leptonPt = cand.pt(); 
	  b1jet_leptonEta = cand.eta(); b1jet_leptonPhi = cand.phi();
	  leadinglepton1_px = cand.px(); 
	  leadinglepton1_py = cand.py(); 
	  leadinglepton1_pz = cand.pz(); 
	  leadinglepton1_p = cand.p(); 
	  leadinglepton1 = true;
	}
	if (leadinglepton1 and b1jet_leadTrackPt>0)
	  break;
    }
    float lepXj1 = (leadinglepton1_px*b1jet_px+leadinglepton1_py*b1jet_py+leadinglepton1_pz*b1jet_pz)/allbjets[jet1].p(); 
    float pTrel2_1 = leadinglepton1_p*leadinglepton1_p - lepXj1*lepXj1;
    if (leadinglepton1)
	b1jet_leptonPtRel = std::sqrt(pTrel2_1);
    else 
	b1jet_leptonPtRel = 0;
    b1jet_leptonDeltaR = deltaR(b1jet_leptonEta, b1jet_leptonPhi, b1jet_eta, b1jet_phi);
    //b1jet_neHEF = allbjets[jet1].neutralHadronEnergy()/(allbjets[jet1].p4()*allbjets[jet1].rawFactor()).energy();//neutralHardonEfraction
    //b1jet_neEmEF = allbjets[jet1].neutralEmEnergy()/(allbjets[jet1].p4()*allbjets[jet1].rawFactor()).energy();//neutralEmEnergyFraction
    b1jet_neHEF = allbjets[jet1].neutralHadronEnergyFraction();
    b1jet_neEmEF = allbjets[jet1].neutralEmEnergyFraction();
    //b1jet_vtxMass = allbjets[jet1].userFloat("vtxMass"); b1jet_vtxNtracks = allbjets[jet1].userFloat("vtxNtracks"); 
    //b1jet_vtxPx = allbjets[jet1].userFloat("vtxPx"); b1jet_vtxPy = allbjets[jet1].userFloat("vtxPy");
    //b1jet_vtxPt = sqrt(b1jet_vtxPx*b1jet_vtxPx + b1jet_vtxPy*b1jet_vtxPy);
    //b1jet_vtx3DSig = allbjets[jet1].userFloat("vtx3DSig"); b1jet_vtx3DVal = allbjets[jet1].userFloat("vtx3DVal");
    //b1jet_vtxPosX = allbjets[jet1].userFloat("vtxPosX"); b1jet_vtxPosY = allbjets[jet1].userFloat("vtxPosY"); b1jet_vtxPosZ = allbjets[jet1].userFloat("vtxPosZ");

    b2jet_px = allbjets[jet2].px(); b2jet_py = allbjets[jet2].py(); b2jet_pz = allbjets[jet2].pz(); b2jet_energy = allbjets[jet2].energy();
    b2jet_pt = allbjets[jet2].pt(); b2jet_eta = allbjets[jet2].eta(); b2jet_phi = allbjets[jet2].phi();
    b2jet_bDiscVar = allbjets[jet2].bDiscriminator(bjetDiscrName_);
    b2jet_mt  = allbjets[jet2].mt();//ptcorrection?
    auto daus2(allbjets[jet2].daughterPtrVector());
    std::sort(daus2.begin(), daus2.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // C++11
    bool leadinglepton2 = false;
    float leadinglepton2_px = 0; float leadinglepton2_py = 0; float leadinglepton2_pz = 0; float leadinglepton2_p = 0;
    for (unsigned int i = 0; i < daus2.size(); ++i) {
	const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus2[i]);
	if (i==0)
	  b2jet_leadTrackPt = cand.pt();
	if ((abs(cand.pdgId()) == 13 or abs(cand.pdgId()) == 11)  and not(leadinglepton2)){
	  b2jet_leptonPdgId = cand.pdgId(); b2jet_leptonPt = cand.pt(); 
	  b2jet_leptonEta = cand.eta(); b2jet_leptonPhi = cand.phi();
	  leadinglepton2_px = cand.px(); 
	  leadinglepton2_py = cand.py(); 
	  leadinglepton2_pz = cand.pz(); 
	  leadinglepton2_p = cand.p(); 
	  leadinglepton2 = true;
	}
	if (leadinglepton2 and b2jet_leadTrackPt>0)
	  break;
    }
    float lepXj2 = (leadinglepton2_px*b2jet_px+leadinglepton2_py*b2jet_py+leadinglepton2_pz*b2jet_pz)/allbjets[jet2].p(); 
    float pTrel2_2 = leadinglepton2_p*leadinglepton2_p - lepXj2*lepXj2;
    if (leadinglepton2)
	b2jet_leptonPtRel = std::sqrt(pTrel2_2);
    else 
	b2jet_leptonPtRel = 0;
    //std::cout <<"b1jet_leptonPtRel "<<b1jet_leptonPtRel <<" b2jet_leptonPtRel "<< b2jet_leptonPtRel << std::endl;
    b2jet_leptonDeltaR = deltaR(b2jet_leptonEta, b2jet_leptonPhi, b2jet_eta, b2jet_phi);
    b2jet_neHEF = allbjets[jet2].neutralHadronEnergyFraction();
    b2jet_neEmEF = allbjets[jet2].neutralEmEnergyFraction();
    //b2jet_neHEF = allbjets[jet2].neutralHadronEnergy()/(allbjets[jet2].p4()*allbjets[jet2].rawFactor()).energy();//neutralHardonEfraction
    //b2jet_neEmEF = allbjets[jet2].neutralEmEnergy()/(allbjets[jet2].p4()*allbjets[jet2].rawFactor()).energy();//neutralEmEnergyFraction
    //b2jet_vtxMass = allbjets[jet2].userFloat("vtxMass"); b2jet_vtxNtracks = allbjets[jet2].userFloat("vtxNtracks"); 
    //b2jet_vtxPx = allbjets[jet2].userFloat("vtxPx"); b2jet_vtxPy = allbjets[jet2].userFloat("vtxPy");
    //b2jet_vtxPt = sqrt(b2jet_vtxPx*b2jet_vtxPx + b2jet_vtxPy*b2jet_vtxPy);
    //b2jet_vtx3DSig = allbjets[jet2].userFloat("vtx3DSig"); b2jet_vtx3DVal = allbjets[jet2].userFloat("vtx3DVal");
    //b2jet_vtxPosX = allbjets[jet2].userFloat("vtxPosX"); b2jet_vtxPosY = allbjets[jet2].userFloat("vtxPosY"); b2jet_vtxPosZ = allbjets[jet2].userFloat("vtxPosZ");

    if (sampleType_>Data){
	float dR_gen_reco11 = deltaR(b1_eta, b1_phi, b1jet_eta, b1jet_phi);
	float dR_gen_reco12 = deltaR(b1_eta, b1_phi, b2jet_eta, b2jet_phi);
	float dR_gen_reco21 = deltaR(b2_eta, b2_phi, b1jet_eta, b1jet_phi);
	float dR_gen_reco22 = deltaR(b2_eta, b2_phi, b2jet_eta, b2jet_phi);
	dR_b1jet = (dR_gen_reco11 < dR_gen_reco21) ? dR_gen_reco11:dR_gen_reco21;
	dR_b2jet = (dR_gen_reco12 < dR_gen_reco22) ? dR_gen_reco12:dR_gen_reco22;
    }
    
    printf("Jet1: pt %5.1f, eta %+4.2f, mt %5.1f, btag_var %+1.2f, \n", b1jet_pt, b1jet_eta, b1jet_mt, b1jet_bDiscVar);
    printf("Jet2: pt %5.1f, eta %+4.2f, mt %5.1f, btag_var %+1.2f, \n", b2jet_pt, b2jet_eta, b2jet_mt, b2jet_bDiscVar);
  }

  if (hastwomuons and hastwojets) std::cout <<"Event has two muons and two bjets " << std::endl;
  else if (hastwomuons and !hastwojets) std::cout <<"Event has two muons BUT not two bjets " << std::endl;
  else if (!hastwomuons and hastwojets) std::cout <<"Event has two jets BUT not two muons " << std::endl;
  else if (!hastwomuons and !hastwojets) std::cout <<"Event does not have two jets nor two muons " << std::endl;

  if (hastwomuons and hastwojets){
    TLorentzVector Muon1_p4(muon1_px, muon1_py, muon1_pz, muon1_energy); 
    TLorentzVector Muon2_p4(muon2_px, muon2_py, muon2_pz, muon2_energy); 
    TLorentzVector b1jet_p4(b1jet_px, b1jet_py, b1jet_pz, b1jet_energy); 
    TLorentzVector b2jet_p4(b2jet_px, b2jet_py, b2jet_pz, b2jet_energy); 
    dR_bl   = (b1jet_p4.Pt()>b2jet_p4.Pt()) ? (b1jet_p4.DeltaR( (Muon1_p4.Pt()>Muon2_p4.Pt()) ? Muon1_p4 : Muon2_p4 )) : (b2jet_p4.DeltaR( (Muon1_p4.Pt()>Muon2_p4.Pt()) ? Muon1_p4 : Muon2_p4 ));
    dR_b1l1 = b1jet_p4.DeltaR(Muon1_p4);
    dR_b1l2 = b1jet_p4.DeltaR(Muon2_p4);
    dR_b2l1 = b2jet_p4.DeltaR(Muon1_p4);
    dR_b2l2 = b2jet_p4.DeltaR(Muon2_p4);
    dR_b1b2 = b1jet_p4.DeltaR(b2jet_p4);
    dR_l1l2 = Muon1_p4.DeltaR(Muon2_p4);
    dR_l1l2b1b2 = (Muon1_p4+Muon2_p4).DeltaR(b1jet_p4+b2jet_p4);
    dphi_l1l2b1b2 = fabs(TVector2::Phi_mpi_pi( (Muon1_p4+Muon2_p4).Phi()-(b1jet_p4+b2jet_p4).Phi() ));
    TLorentzVector ll_p4 = Muon1_p4+Muon2_p4;
    TLorentzVector bjets_p4 = b1jet_p4+b2jet_p4;
    dR_minbl = std::min(std::min(dR_b1l1,dR_b1l2),std::min(dR_b2l1,dR_b2l2));
    dphi_llbb = TVector2::Phi_mpi_pi(ll_p4.Phi()-bjets_p4.Phi());
    dphi_llmet = TVector2::Phi_mpi_pi(ll_p4.Phi()-met_phi);
    mass_l1l2 = ll_p4.M(); energy_l1l2 = ll_p4.Energy(); pt_l1l2 = ll_p4.Pt(); eta_l1l2 = ll_p4.Eta(); phi_l1l2 = ll_p4.Phi();
    mass_b1b2 = bjets_p4.M(); energy_b1b2 = bjets_p4.Energy(); pt_b1b2 = bjets_p4.Pt(); eta_b1b2 = bjets_p4.Eta(); phi_b1b2 = bjets_p4.Phi();
    mass_trans = sqrt(2*ll_p4.Pt()*met_pt*(1-cos(dphi_llmet)));
    //if (dR_b1l1 > jetleptonDeltaR_ and dR_b1l2 > jetleptonDeltaR_ and dR_b2l1 > jetleptonDeltaR_ and dR_b2l2 > jetleptonDeltaR_) hasdRljet =true;
    //MINdR_bl = dR_b1l1*(dR_b1l1<dR_b1l2 && dR_b1l1<dR_b2l1 && dR_b1l1<dR_b2l2) + dR_b2l1*(dR_b2l1<dR_b2l2 && dR_b2l1<dR_b1l1 && dR_b2l1<dR_b1l2) + dR_b1l2*(dR_b1l2<dR_b1l1 && dR_b1l2<dR_b2l1 && dR_b1l2<dR_b2l2) + dR_b2l2*(dR_b2l2<dR_b1l1 && dR_b2l2<dR_b1l2 && dR_b2l2<dR_b2l1);

    //MT2: In order to construct MT2 for either the t tbar -> bW bW system or our signal H -> h h -> bb WW,
    //we group each pair b_jet-lepton into an object: we then get "Particle A" and "Particle B" (each one given by a b_jet-lepton object),
    //whose Kinematics is to be fed into the MT2 variable.
    //There are two possible Lepton-Bquark pairings. We compute MT2 for both and pick the smallest value.

  }*/
  //if((hastwomuons and hastwojets) or findAllGenParticles) evtree->Fill();
  if(findAllGenParticles) evtree->Fill();
}

void DiHiggsWWBBSLAnalyzer::beginJob(){
  hevent = fs->make<TH1F>("hevent", "event counter",10,0,10);
  evtree = fs->make<TTree>("evtree", "evtree");
  //output = new TFile("output.root","recreate");
  // output->cd();
  evtree->Branch("ievent",&ievent);
  //evtree = new TTree("evtree","event tree");
  evtree->Branch("findAllGenParticles",&findAllGenParticles);
  evtree->Branch("XsecBr",&XsecBr, "XsecBr/F");
  evtree->Branch("lep1_px",&lep1_px, "lep1_px/F");
  evtree->Branch("lep1_py",&lep1_py, "lep1_py/F");
  evtree->Branch("lep1_pz",&lep1_pz, "lep1_pz/F");
  evtree->Branch("lep1_eta",&lep1_eta, "lep1_eta/F");
  evtree->Branch("lep1_phi",&lep1_phi, "lep1_phi/F");
  evtree->Branch("lep1_pt",&lep1_pt, "lep1_pt/F");
  evtree->Branch("lep1_energy",&lep1_energy, "lep1_energy/F");
  evtree->Branch("lep1_pdgid",&lep1_pdgid, "lep1_pdgid/I");
  evtree->Branch("nu1_px",&nu1_px, "nu1_px/F");
  evtree->Branch("nu1_py",&nu1_py, "nu1_py/F");
  evtree->Branch("nu1_pz",&nu1_pz, "nu1_pz/F");
  evtree->Branch("nu1_eta",&nu1_eta, "nu1_eta/F");
  evtree->Branch("nu1_phi",&nu1_phi, "nu1_phi/F");
  evtree->Branch("nu1_pt",&nu1_pt, "nu1_pt/F");
  evtree->Branch("nu1_energy",&nu1_energy, "nu1_energy/F");
  evtree->Branch("nu1_pdgid",&nu1_pdgid, "nu1_pdgid/I");
  evtree->Branch("lep2_px",&lep2_px, "lep2_px/F");
  evtree->Branch("lep2_py",&lep2_py, "lep2_py/F");
  evtree->Branch("lep2_pz",&lep2_pz, "lep2_pz/F");
  evtree->Branch("lep2_eta",&lep2_eta, "lep2_eta/F");
  evtree->Branch("lep2_phi",&lep2_phi, "lep2_phi/F");
  evtree->Branch("lep2_pt",&lep2_pt, "lep2_pt/F");
  evtree->Branch("lep2_energy",&lep2_energy, "lep2_energy/F");
  evtree->Branch("lep2_pdgid",&lep2_pdgid, "lep2_pdgid/I");
  evtree->Branch("nu2_px",&nu2_px, "nu2_px/F");
  evtree->Branch("nu2_py",&nu2_py, "nu2_py/F");
  evtree->Branch("nu2_pz",&nu2_pz, "nu2_pz/F");
  evtree->Branch("nu2_eta",&nu2_eta, "nu2_eta/F");
  evtree->Branch("nu2_phi",&nu2_phi, "nu2_phi/F");
  evtree->Branch("nu2_pt",&nu2_pt, "nu2_pt/F");
  evtree->Branch("nu2_energy",&nu2_energy, "nu2_energy/F");
  evtree->Branch("nu2_pdgid",&nu2_pdgid, "nu2_pdgid/I");
  //evtree->Branch("nu1and2_pt",&nu1and2_pt, "nu1and2_pt/F");
  //evtree->Branch("nu1and2_px",&nu1and2_px, "nu1and2_px/F");
  //evtree->Branch("nu1and2_py",&nu1and2_py, "nu1and2_py/F");
  //evtree->Branch("nu1and2_diBaxis_p",&nu1and2_diBaxis_p, "nu1and2_diBaxis_p/F");
  //evtree->Branch("nu1and2_diBaxis_t",&nu1and2_diBaxis_t, "nu1and2_diBaxis_t/F");
  evtree->Branch("q1_px",    &q1_px,     "q1_px/F");
  evtree->Branch("q1_py",    &q1_py,     "q1_py/F");
  evtree->Branch("q1_pz",    &q1_pz,     "q1_pz/F");
  evtree->Branch("q1_eta",   &q1_eta,    "q1_eta/F");
  evtree->Branch("q1_phi",   &q1_phi,    "q1_phi/F");
  evtree->Branch("q1_pt",    &q1_pt,     "q1_pt/F");
  evtree->Branch("q1_energy",&q1_energy, "q1_energy/F");
  evtree->Branch("q1_pdgid", &q1_pdgid,  "q1_pdgid/I");
  evtree->Branch("q2_px",    &q2_px,     "q2_px/F");
  evtree->Branch("q2_py",    &q2_py,     "q2_py/F");
  evtree->Branch("q2_pz",    &q2_pz,     "q2_pz/F");
  evtree->Branch("q2_eta",   &q2_eta,    "q2_eta/F");
  evtree->Branch("q2_phi",   &q2_phi,    "q2_phi/F");
  evtree->Branch("q2_pt",    &q2_pt,     "q2_pt/F");
  evtree->Branch("q2_energy",&q2_energy, "q2_energy/F");
  evtree->Branch("q2_pdgid", &q2_pdgid,  "q2_pdgid/I");

  evtree->Branch("w1_mass",&w1_mass, "w1_mass/F");
  evtree->Branch("w1_px",&w1_px, "w1_px/F");
  evtree->Branch("w1_py",&w1_py, "w1_py/F");
  evtree->Branch("w1_pz",&w1_pz, "w1_pz/F");
  evtree->Branch("w1_energy",&w1_energy, "w1_energy/F");
  evtree->Branch("w1_pt",&w1_pt, "w1_pt/F");
  evtree->Branch("w1_eta",&w1_eta, "w1_eta/F");
  evtree->Branch("w1_phi",&w1_phi, "w1_phi/F");
  evtree->Branch("w2_mass",&w2_mass, "w2_mass/F");
  evtree->Branch("w2_px",&w2_px, "w2_px/F");
  evtree->Branch("w2_py",&w2_py, "w2_py/F");
  evtree->Branch("w2_pz",&w2_pz, "w2_pz/F");
  evtree->Branch("w2_energy",&w2_energy, "w2_energy/F");
  evtree->Branch("w2_pt",&w2_pt, "w2_pt/F");
  evtree->Branch("w2_eta",&w2_eta, "w2_eta/F");
  evtree->Branch("w2_phi",&w2_phi, "w2_phi/F");

  evtree->Branch("htoWW_energy",&htoWW_energy,"htoWW_energy/F");
  evtree->Branch("htoWW_px",&htoWW_px,"htoWW_px/F");
  evtree->Branch("htoWW_py",&htoWW_py,"htoWW_px/F");
  evtree->Branch("htoWW_pz",&htoWW_pz,"htoWW_pz/F");
  evtree->Branch("htoWW_mass",&htoWW_mass,"htoWW_mass/F");
  //evtree->Branch("Wtolep1nu1",&Wtolep1nu1,"Wtolep1nu1/B");
  //evtree->Branch("Wtolep2nu2",&Wtolep2nu2,"Wtolep2nu2/B");
  evtree->Branch("htoWW",&htoWW,"htoWW/B");

  //ZZbb
  evtree->Branch("z1_mass",&z1_mass, "z1_mass/F");
  evtree->Branch("z1_px",&z1_px, "z1_px/F");
  evtree->Branch("z1_py",&z1_py, "z1_py/F");
  evtree->Branch("z1_pz",&z1_pz, "z1_pz/F");
  evtree->Branch("z1_energy",&z1_energy, "z1_energy/F");
  evtree->Branch("z1_pt",&z1_pt, "z1_pt/F");
  evtree->Branch("z1_eta",&z1_eta, "z1_eta/F");
  evtree->Branch("z1_phi",&z1_phi, "z1_phi/F");
  evtree->Branch("z2_mass",&z2_mass, "z2_mass/F");
  evtree->Branch("z2_px",&z2_px, "z2_px/F");
  evtree->Branch("z2_py",&z2_py, "z2_py/F");
  evtree->Branch("z2_pz",&z2_pz, "z2_pz/F");
  evtree->Branch("z2_energy",&z2_energy, "z2_energy/F");
  evtree->Branch("z2_pt",&z2_pt, "z2_pt/F");
  evtree->Branch("z2_eta",&z2_eta, "z2_eta/F");
  evtree->Branch("z2_phi",&z2_phi, "z2_phi/F");
  evtree->Branch("htoZZ_energy",&htoZZ_energy,"htoZZ_energy/F");
  evtree->Branch("htoZZ_px",&htoZZ_px,"htoZZ_px/F");
  evtree->Branch("htoZZ_py",&htoZZ_py,"htoZZ_px/F");
  evtree->Branch("htoZZ_pz",&htoZZ_pz,"htoZZ_pz/F");
  evtree->Branch("htoZZ_mass",&htoZZ_mass,"htoZZ_mass/F");
  evtree->Branch("htoZZ",&htoZZ,"htoZZ/B");

  evtree->Branch("b1_px",&b1_px, "b1_px/F");
  evtree->Branch("b1_py",&b1_py, "b1_py/F");
  evtree->Branch("b1_pz",&b1_pz, "b1_pz/F");
  evtree->Branch("b1_eta",&b1_eta, "b1_eta/F");
  evtree->Branch("b1_phi",&b1_phi, "b1_phi/F");
  evtree->Branch("b1_pt",&b1_pt, "b1_pt/F");
  evtree->Branch("b1_energy",&b1_energy, "b1_energy/F");
  evtree->Branch("b2_px",&b2_px, "b2_px/F");
  evtree->Branch("b2_py",&b2_py, "b2_py/F");
  evtree->Branch("b2_pz",&b2_pz, "b2_pz/F");
  evtree->Branch("b2_eta",&b2_eta, "b2_eta/F");
  evtree->Branch("b2_phi",&b2_phi, "b2_phi/F");
  evtree->Branch("b2_pt",&b2_pt, "b2_pt/F");
  evtree->Branch("b2_energy",&b2_energy, "b2_energy/F");
  evtree->Branch("htobb_px",&htobb_px, "htobb_px/F");
  evtree->Branch("htobb_py",&htobb_py, "htobb_py/F");
  evtree->Branch("htobb_pz",&htobb_pz, "htobb_pz/F");
  evtree->Branch("htobb_energy",&htobb_energy, "htobb_energy/F");
  evtree->Branch("htobb_mass",&htobb_mass, "htobb_mass/F");
  evtree->Branch("Xtohh_px",&Xtohh_px,"Xtohh_px/F");
  evtree->Branch("Xtohh_py",&Xtohh_py,"Xtohh_py/F");
  evtree->Branch("Xtohh_pz",&Xtohh_pz,"Xtohh_pz/F");
  evtree->Branch("Xtohh_energy",&Xtohh_energy,"Xtohh_energy/F");
  evtree->Branch("Xtohh_mass",&Xtohh_mass,"Xtohh_mass/F");

  evtree->Branch("b1genjet_px",&b1genjet_px, "b1genjet_px/F");
  evtree->Branch("b1genjet_py",&b1genjet_py, "b1genjet_py/F");
  evtree->Branch("b1genjet_pz",&b1genjet_pz, "b1genjet_pz/F");
  evtree->Branch("b1genjet_eta",&b1genjet_eta, "b1genjet_eta/F");
  evtree->Branch("b1genjet_phi",&b1genjet_phi, "b1genjet_phi/F");
  evtree->Branch("b1genjet_pt",&b1genjet_pt, "b1genjet_pt/F");
  evtree->Branch("b1genjet_energy",&b1genjet_energy, "b1genjet_energy/F");
  evtree->Branch("b1genjet_mass",&b1genjet_mass, "b1genjet_mass/F");
  evtree->Branch("b2genjet_px",&b2genjet_px, "b2genjet_px/F");
  evtree->Branch("b2genjet_py",&b2genjet_py, "b2genjet_py/F");
  evtree->Branch("b2genjet_pz",&b2genjet_pz, "b2genjet_pz/F");
  evtree->Branch("b2genjet_eta",&b2genjet_eta, "b2genjet_eta/F");
  evtree->Branch("b2genjet_phi",&b2genjet_phi, "b2genjet_phi/F");
  evtree->Branch("b2genjet_pt",&b2genjet_pt, "b2genjet_pt/F");
  evtree->Branch("b2genjet_energy",&b2genjet_energy, "b2genjet_energy/F");
  evtree->Branch("b2genjet_mass",&b2genjet_mass, "b2genjet_mass/F");
  evtree->Branch("dR_b1genjet", &dR_b1genjet,"dR_b1genjet/F");  
  evtree->Branch("dR_b2genjet", &dR_b2genjet,"dR_b2genjet/F");  
  evtree->Branch("hastwogenjets", &hastwogenjets,"hastwogenjets/B");  

  evtree->Branch("q1genjet_px"    ,&q1genjet_px,     "q1genjet_px/F");
  evtree->Branch("q1genjet_py"    ,&q1genjet_py,     "q1genjet_py/F");
  evtree->Branch("q1genjet_pz"    ,&q1genjet_pz,     "q1genjet_pz/F");
  evtree->Branch("q1genjet_eta"   ,&q1genjet_eta,    "q1genjet_eta/F");
  evtree->Branch("q1genjet_phi"   ,&q1genjet_phi,    "q1genjet_phi/F");
  evtree->Branch("q1genjet_pt"    ,&q1genjet_pt,     "q1genjet_pt/F");
  evtree->Branch("q1genjet_energy",&q1genjet_energy, "q1genjet_energy/F");
  evtree->Branch("q1genjet_mass"  ,&q1genjet_mass,   "q1genjet_mass/F");
  evtree->Branch("q2genjet_px"    ,&q2genjet_px,     "q2genjet_px/F");
  evtree->Branch("q2genjet_py"    ,&q2genjet_py,     "q2genjet_py/F");
  evtree->Branch("q2genjet_pz"    ,&q2genjet_pz,     "q2genjet_pz/F");
  evtree->Branch("q2genjet_eta"   ,&q2genjet_eta,    "q2genjet_eta/F");
  evtree->Branch("q2genjet_phi"   ,&q2genjet_phi,    "q2genjet_phi/F");
  evtree->Branch("q2genjet_pt"    ,&q2genjet_pt,     "q2genjet_pt/F");
  evtree->Branch("q2genjet_energy",&q2genjet_energy, "q2genjet_energy/F");
  evtree->Branch("q2genjet_mass"  ,&q2genjet_mass,   "q2genjet_mass/F");
  evtree->Branch("dR_q1genjet", &dR_q1genjet,"dR_q1genjet/F");  
  evtree->Branch("dR_q2genjet", &dR_q2genjet,"dR_q2genjet/F");  
  evtree->Branch("hastwogenjets_wqq", &hastwogenjets_wqq,"hastwogenjets_wqq/B");  


  evtree->Branch("genmet_neutrinos_pt",&genmet_neutrinos_pt,"genmet_neutrinos_pt/F");
  evtree->Branch("genmet_neutrinos_phi",&genmet_neutrinos_phi,"genmet_neutrinos_phi/F");
  evtree->Branch("genmet_neutrinos_px",&genmet_neutrinos_px,"genmet_neutrinos_px/F");
  evtree->Branch("genmet_neutrinos_py",&genmet_neutrinos_py,"genmet_neutrinos_py/F");
  evtree->Branch("genmet_pt",&genmet_pt,"genmet_pt/F");
  evtree->Branch("genmet_phi",&genmet_phi,"genmet_phi/F");
  evtree->Branch("genmet_px",&genmet_px,"genmet_px/F");
  evtree->Branch("genmet_py",&genmet_py,"genmet_py/F");

  evtree->Branch("t1_px",&t1_px,"t1_px/F");
  evtree->Branch("t1_py",&t1_py,"t1_py/F");
  evtree->Branch("t1_pz",&t1_pz,"t1_pz/F");
  evtree->Branch("t1_energy",&t1_energy,"t1_energy/F");
  evtree->Branch("t1_mass",&t1_mass,"t1_mass/F");
  evtree->Branch("t2_px",&t2_px,"t2_px/F");
  evtree->Branch("t2_py",&t2_py,"t2_py/F");
  evtree->Branch("t2_pz",&t2_pz,"t2_pz/F");
  evtree->Branch("t2_energy",&t2_energy,"t2_energy/F");
  evtree->Branch("t2_mass",&t2_mass,"t2_mass/F");

  evtree->Branch("dR_genbl",&dR_genbl, "dR_genbl/F");
  evtree->Branch("dR_genb1l1",&dR_genb1l1, "dR_genb1l1/F");
  evtree->Branch("dR_genb1l2",&dR_genb1l2, "dR_genb1l2/F");
  evtree->Branch("dR_genb2l1",&dR_genb2l1, "dR_genb2l1/F");
  evtree->Branch("dR_genb2l2",&dR_genb2l2, "dR_genb2l2/F");
  evtree->Branch("dR_genb1b2",&dR_genb1b2, "dR_genb1b2/F");
  evtree->Branch("dR_genl1l2",&dR_genl1l2, "dR_genl1l2/F");
  evtree->Branch("dR_genl1l2b1b2",&dR_genl1l2b1b2, "dR_genl1l2b1b2/F");
  evtree->Branch("dphi_genl1l2b1b2",&dphi_genl1l2b1b2, "dphi_genl1l2b1b2/F");
  evtree->Branch("dR_genminbl",&dR_genminbl, "dR_genminbl/F");
  evtree->Branch("mass_genb1b2",&mass_genb1b2, "mass_genb1b2/F");
  evtree->Branch("energy_genb1b2",&energy_genb1b2, "energy_genb1b2/F");
  evtree->Branch("pt_genb1b2",&pt_genb1b2, "pt_genb1b2/F");
  evtree->Branch("phi_genb1b2",&phi_genb1b2, "phi_genb1b2/F");
  evtree->Branch("eta_genb1b2",&eta_genb1b2, "eta_genb1b2/F");
  evtree->Branch("mass_genl1l2",&mass_genl1l2, "mass_genl1l2/F");
  evtree->Branch("mass_gennu1nu2",&mass_gennu1nu2, "mass_gennu1nu2/F");
  evtree->Branch("energy_genl1l2",&energy_genl1l2, "energy_genl1l2/F");
  evtree->Branch("pt_genl1l2",&pt_genl1l2, "pt_genl1l2/F");
  evtree->Branch("phi_genl1l2",&phi_genl1l2, "phi_genl1l2/F");
  evtree->Branch("eta_genl1l2",&eta_genl1l2, "eta_genl1l2/F");
  evtree->Branch("dphi_genllbb",&dphi_genllbb, "dphi_genllbb/F");
  evtree->Branch("dphi_genllmet",&dphi_genllmet, "dphi_genllmet/F");
  evtree->Branch("mass_gentrans",&mass_gentrans, "mass_gentrans/F");
  //Wtoqq
  evtree->Branch("mass_genq1q2",    &mass_genq1q2,   "mass_genq1q2/F");
  evtree->Branch("energy_genq1q2",&energy_genq1q2, "energy_genq1q2/F");
  evtree->Branch("pt_genq1q2",        &pt_genq1q2,     "pt_genq1q2/F");
  evtree->Branch("phi_genq1q2",      &phi_genq1q2,    "phi_genq1q2/F");
  evtree->Branch("eta_genq1q2",      &eta_genq1q2,    "eta_genq1q2/F");
  //reco level
  evtree->Branch("numOfVertices",&numOfVertices, "numOfVertices/I");
  evtree->Branch("muon1_px",&muon1_px, "muon1_px/F");
  evtree->Branch("muon1_py",&muon1_py, "muon1_py/F");
  evtree->Branch("muon1_pz",&muon1_pz, "muon1_pz/F");
  evtree->Branch("muon1_eta",&muon1_eta, "muon1_eta/F");
  evtree->Branch("muon1_phi",&muon1_phi, "muon1_phi/F");
  evtree->Branch("muon1_pt",&muon1_pt, "muon1_pt/F");
  evtree->Branch("muon1_energy",&muon1_energy, "muon1_energy/F");
  evtree->Branch("muon1_isoVar",&muon1_isoVar, "muon1_isoVar/F");
  evtree->Branch("muon1_triggerSF",&muon1_triggerSF, "muon1_triggerSF/F");
  evtree->Branch("muon1_isoSF",&muon1_isoSF, "muon1_isoSF/F");
  evtree->Branch("muon1_idSF",&muon1_idSF, "muon1_idSF/F");
  evtree->Branch("muon1_trackingSF",&muon1_trackingSF, "muon1_trackingSF/F");
  evtree->Branch("muon1_pogSF",&muon1_pogSF, "muon1_pogSF/F");
  evtree->Branch("muon1_trigger_dR",&muon1_trigger_dR, "muon1_trigger_dR/F");
  evtree->Branch("muon1_trigger_dPtRel",&muon1_trigger_dPtRel, "muon1_trigger_dPtRel/F");
  evtree->Branch("muon1_dxy",&muon1_dxy, "muon1_dxy/F");
  evtree->Branch("muon1_dz",&muon1_dz, "muon1_dz/F");
  evtree->Branch("muon2_px",&muon2_px, "muon2_px/F");
  evtree->Branch("muon2_py",&muon2_py, "muon2_py/F");
  evtree->Branch("muon2_pz",&muon2_pz, "muon2_pz/F");
  evtree->Branch("muon2_eta",&muon2_eta, "muon2_eta/F");
  evtree->Branch("muon2_phi",&muon2_phi, "muon2_phi/F");
  evtree->Branch("muon2_pt",&muon2_pt, "muon2_pt/F");
  evtree->Branch("muon2_energy",&muon2_energy, "muon2_energy/F");
  evtree->Branch("muon2_dxy",&muon1_dxy, "muon2_dxy/F");
  evtree->Branch("muon2_dz",&muon1_dz, "muon2_dz/F");
  evtree->Branch("muon2_isoVar",&muon2_isoVar, "muon2_isoVar/F");
  evtree->Branch("muon2_triggerSF",&muon2_triggerSF, "muon2_triggerSF/F");
  evtree->Branch("muon2_isoSF",&muon2_isoSF, "muon2_isoSF/F");
  evtree->Branch("muon2_idSF",&muon2_idSF, "muon2_idSF/F");
  evtree->Branch("muon2_trackingSF",&muon2_trackingSF, "muon2_trackingSF/F");
  evtree->Branch("muon2_pogSF",&muon2_pogSF, "muon2_pogSF/F");
  evtree->Branch("muon2_trigger_dR",&muon2_trigger_dR, "muon2_trigger_dR/F");
  evtree->Branch("muon2_trigger_dPtRel",&muon2_trigger_dPtRel, "muon2_trigger_dPtRel/F");
  evtree->Branch("dR_lep1",&dR_lep1, "dR_lep1/F");
  evtree->Branch("dR_lep2",&dR_lep2, "dR_lep2/F");
  evtree->Branch("hastwomuons",&hastwomuons, "hastwomuons/B");
  //evtree->Branch("hasmuon1",&hasmuon1, "hasmuon1/B");
  //evtree->Branch("hasmuon2",&hasmuon2, "hasmuon2/B");
  //evtree->Branch("hasRecomuon1",&hasRecomuon1, "hasRecomuon1/B");
  //evtree->Branch("hasRecomuon2",&hasRecomuon2, "hasRecomuon2/B");
  //evtree->Branch("muon1_hasgenMu",&muon1_hasgenMu, "muon1_hasgenMu/B");
  //evtree->Branch("muon2_hasgenMu",&muon2_hasgenMu, "muon2_hasgenMu/B");

  evtree->Branch("b1jet_px",&b1jet_px, "b1jet_px/F");
  evtree->Branch("b1jet_py",&b1jet_py, "b1jet_py/F");
  evtree->Branch("b1jet_pz",&b1jet_pz, "b1jet_pz/F");
  evtree->Branch("b1jet_eta",&b1jet_eta, "b1jet_eta/F");
  evtree->Branch("b1jet_phi",&b1jet_phi, "b1jet_phi/F");
  evtree->Branch("b1jet_pt",&b1jet_pt, "b1jet_pt/F");
  evtree->Branch("b1jet_energy",&b1jet_energy, "b1jet_energy/F");
  evtree->Branch("b1jet_mass",&b1jet_mass, "b1jet_mass/F");
  evtree->Branch("b1jet_btag",&b1jet_btag, "b1jet_btag/i");//unsigned int
  evtree->Branch("b1jet_bDiscVar",&b1jet_bDiscVar, "b1jet_bDiscVar/F");
  evtree->Branch("b1jet_mt",&b1jet_mt, "b1jet_mt/F");
  evtree->Branch("b1jet_leadTrackPt",&b1jet_leadTrackPt, "b1jet_leadTrackPt/F");
  evtree->Branch("b1jet_leptonPdgId",&b1jet_leptonPdgId, "b1jet_leptonPdgId/F");
  evtree->Branch("b1jet_leptonPhi",&b1jet_leptonPhi, "b1jet_leptonPhi/F");
  evtree->Branch("b1jet_leptonEta",&b1jet_leptonEta, "b1jet_leptonEta/F");
  evtree->Branch("b1jet_leptonPt",&b1jet_leptonPt, "b1jet_leptonPt/F");
  evtree->Branch("b1jet_leptonPtRel",&b1jet_leptonPtRel, "b1jet_leptonPtRel/F");
  evtree->Branch("b1jet_leptonDeltaR",&b1jet_leptonDeltaR, "b1jet_leptonDeltaR/F");
  evtree->Branch("b1jet_neHEF",&b1jet_neHEF, "b1jet_neHEF/F");
  evtree->Branch("b1jet_neEmEF",&b1jet_neEmEF, "b1jet_neEmEF/F");
  evtree->Branch("b1jet_vtxMass",&b1jet_vtxMass, "b1jet_vtxMass/F");
  evtree->Branch("b1jet_vtxNtracks",&b1jet_vtxNtracks, "b1jet_vtxNtracks/F");
  evtree->Branch("b1jet_vtxPx",&b1jet_vtxPx, "b1jet_vtxPx/F");
  evtree->Branch("b1jet_vtxPy",&b1jet_vtxPy, "b1jet_vtxPy/F");
  evtree->Branch("b1jet_vtxPt",&b1jet_vtxPt, "b1jet_vtxPt/F");
  evtree->Branch("b1jet_vtx3DSig",&b1jet_vtx3DSig, "b1jet_vtx3DSig/F");
  evtree->Branch("b1jet_vtx3DVal",&b1jet_vtx3DVal, "b1jet_vtx3DVal/F");
  evtree->Branch("b1jet_vtxPosX",&b1jet_vtxPosX, "b1jet_vtxPosX/F");
  evtree->Branch("b1jet_vtxPosY",&b1jet_vtxPosY, "b1jet_vtxPosY/F");
  evtree->Branch("b1jet_vtxPosZ",&b1jet_vtxPosZ, "b1jet_vtxPosZ/F");

  evtree->Branch("b2jet_px",&b2jet_px, "b2jet_px/F");
  evtree->Branch("b2jet_py",&b2jet_py, "b2jet_py/F");
  evtree->Branch("b2jet_pz",&b2jet_pz, "b2jet_pz/F");
  evtree->Branch("b2jet_eta",&b2jet_eta, "b2jet_eta/F");
  evtree->Branch("b2jet_phi",&b2jet_phi, "b2jet_phi/F");
  evtree->Branch("b2jet_pt",&b2jet_pt, "b2jet_pt/F");
  evtree->Branch("b2jet_energy",&b2jet_energy, "b2jet_energy/F");
  evtree->Branch("b2jet_mass",&b2jet_mass, "b2jet_mass/F");
  evtree->Branch("b2jet_btag",&b2jet_btag, "b2jet_btag/i");
  evtree->Branch("b2jet_bDiscVar",&b2jet_bDiscVar, "b2jet_bDiscVar/F");
  evtree->Branch("b2jet_mt",&b2jet_mt, "b2jet_mt/F");
  evtree->Branch("b2jet_leadTrackPt",&b2jet_leadTrackPt, "b2jet_leadTrackPt/F");
  evtree->Branch("b2jet_leptonPdgId",&b2jet_leptonPdgId, "b2jet_leptonPdgId/F");
  evtree->Branch("b2jet_leptonPhi",&b2jet_leptonPhi, "b2jet_leptonPhi/F");
  evtree->Branch("b2jet_leptonEta",&b2jet_leptonEta, "b2jet_leptonEta/F");
  evtree->Branch("b2jet_leptonPt",&b2jet_leptonPt, "b2jet_leptonPt/F");
  evtree->Branch("b2jet_leptonPtRel",&b2jet_leptonPtRel, "b2jet_leptonPtRel/F");
  evtree->Branch("b2jet_leptonDeltaR",&b2jet_leptonDeltaR, "b2jet_leptonDeltaR/F");
  evtree->Branch("b2jet_neHEF",&b2jet_neHEF, "b2jet_neHEF/F");
  evtree->Branch("b2jet_neEmEF",&b2jet_neEmEF, "b2jet_neEmEF/F");
  evtree->Branch("b2jet_vtxMass",&b2jet_vtxMass, "b2jet_vtxMass/F");
  evtree->Branch("b2jet_vtxNtracks",&b2jet_vtxNtracks, "b2jet_vtxNtracks/F");
  evtree->Branch("b2jet_vtxPx",&b2jet_vtxPx, "b2jet_vtxPx/F");
  evtree->Branch("b2jet_vtxPy",&b2jet_vtxPy, "b2jet_vtxPy/F");
  evtree->Branch("b2jet_vtxPt",&b2jet_vtxPt, "b2jet_vtxPt/F");
  evtree->Branch("b2jet_vtx3DSig",&b2jet_vtx3DSig, "b2jet_vtx3DSig/F");
  evtree->Branch("b2jet_vtx3DVal",&b2jet_vtx3DVal, "b2jet_vtx3DVal/F");
  evtree->Branch("b2jet_vtxPosX",&b2jet_vtxPosX, "b2jet_vtxPosX/F");
  evtree->Branch("b2jet_vtxPosY",&b2jet_vtxPosY, "b2jet_vtxPosY/F");
  evtree->Branch("b2jet_vtxPosZ",&b2jet_vtxPosZ, "b2jet_vtxPosZ/F");

  evtree->Branch("dR_b1jet", &dR_b1jet,"dR_b1jet/F");  
  evtree->Branch("dR_b2jet", &dR_b2jet,"dR_b2jet/F");  
  evtree->Branch("hastwojets",&hastwojets, "hastwojets/B");
  evtree->Branch("hasb1jet",&hasb1jet, "hasb1jet/B");
  evtree->Branch("hasb2jet",&hasb2jet, "hasb2jet/B");

  evtree->Branch("met_pt",&met_pt,"met_pt/F");
  evtree->Branch("met_phi",&met_phi,"met_phi/F");
  evtree->Branch("met_px",&met_px,"met_px/F");
  evtree->Branch("met_py",&met_py,"met_py/F");

  evtree->Branch("HT",&HT,"HT/F");

  evtree->Branch("dR_bl",&dR_bl, "dR_bl/F");
  evtree->Branch("dR_b1l1",&dR_b1l1, "dR_b1l1/F");
  evtree->Branch("dR_b1l2",&dR_b1l2, "dR_b1l2/F");
  evtree->Branch("dR_b2l1",&dR_b2l1, "dR_b2l1/F");
  evtree->Branch("dR_b2l2",&dR_b2l2, "dR_b2l2/F");
  evtree->Branch("dR_b1b2",&dR_b1b2, "dR_b1b2/F");
  evtree->Branch("MINdR_bl",&MINdR_bl, "MINdR_bl/F");
  evtree->Branch("dR_l1l2",&dR_l1l2, "dR_l1l2/F");
  evtree->Branch("dR_l1l2b1b2",&dR_l1l2b1b2, "dR_l1l2b1b2/F");
  evtree->Branch("dphi_l1l2b1b2",&dphi_l1l2b1b2, "dphi_l1l2b1b2/F");
  evtree->Branch("dR_minbl",&dR_minbl, "dR_minbl/F");
  evtree->Branch("dR_genl1l2",&dR_genl1l2, "dR_genl1l2/F");
  evtree->Branch("mass_b1b2",&mass_b1b2, "mass_b1b2/F");
  evtree->Branch("energy_b1b2",&energy_b1b2, "energy_b1b2/F");
  evtree->Branch("pt_b1b2",&pt_b1b2, "pt_b1b2/F");
  evtree->Branch("phi_b1b2",&phi_b1b2, "phi_b1b2/F");
  evtree->Branch("eta_b1b2",&eta_b1b2, "eta_b1b2/F");
  evtree->Branch("mass_l1l2",&mass_l1l2, "mass_l1l2/F");
  evtree->Branch("energy_l1l2",&energy_l1l2, "energy_l1l2/F");
  evtree->Branch("pt_l1l2",&pt_l1l2, "pt_l1l2/F");
  evtree->Branch("phi_l1l2",&phi_l1l2, "phi_l1l2/F");
  evtree->Branch("eta_l1l2",&eta_l1l2, "eta_l1l2/F");
  evtree->Branch("dphi_llbb",&dphi_llbb, "dphi_llbb/F");
  evtree->Branch("dphi_llmet",&dphi_llmet, "dphi_llmet/F");
  evtree->Branch("mass_trans",&mass_trans, "mass_trans/F");   
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiHiggsWWBBSLAnalyzer::endJob() {
  hevent->Fill(1, ievent);
  //   std::cout << "endJob, ievent  " << ievent << std::endl;
  //output->Write();
  // output->Close();
  /*
  // release space 
  delete mu_onshellW_lorentz;
  delete mu_offshellW_lorentz;
  delete nu_onshellW_lorentz;
  delete nu_offshellW_lorentz;
  delete onshellW_lorentz;
  delete offshellW_lorentz;
  delete htoWW_lorentz;
  delete htoBB_lorentz;
  delete Xtohh_lorentz;*/
  //delete jets_lorentz;
}

// ------------ method called to check singal genParticles   ------------
void DiHiggsWWBBSLAnalyzer::checkGenParticlesHbbWW2L2Nu(edm::Handle<reco::GenParticleCollection> genParticleColl){

  std::vector<reco::GenParticle*> b1Coll; 
  std::vector<reco::GenParticle*> b2Coll;
  std::vector<reco::GenParticle*> W1Coll;
  std::vector<reco::GenParticle*> W2Coll;
  std::vector<const reco::Candidate*> htoWWColl;
  std::vector<const reco::Candidate*> htoBBColl;

  std::cout <<"*********** start to check GenParticles for HHtobbWW to bb2L2Nu sample ***********"<< std::endl;
  bool Xtohh = false;
  for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
    //particle id, (electron12)(muon13),(b5),(W+24),(SM higgs25)
    // particle id  it->pdgId()
    //std::cout << "Gen paticles: id " << it->pdgId() << std::endl; 
    //if (it->pdgId() == 24 && hasDaughter(it->clone(), -13) )
    if (it->pdgId() == 24){
      const reco::Candidate* tmpw1 = it->mother();
      while (tmpw1->pdgId() == 24 && tmpw1->numberOfMothers() == 1) tmpw1 = tmpw1->mother();
      if (tmpw1->pdgId() == 25)  W1Coll.push_back(it->clone());
    }
    else if (it->pdgId() == -24){
      const reco::Candidate* tmpw2 = it->mother();
      while (tmpw2->pdgId() == -24 && tmpw2->numberOfMothers() == 1) tmpw2 = tmpw2->mother();
      if (tmpw2->pdgId() == 25)  W2Coll.push_back(it->clone());
    }
    else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 )
    {
      if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
      //   std::cout << "find bquark candidate" << std::endl;
      b1Coll.push_back(it->clone());
    }
    else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 )
    {
      if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
      // std::cout << "find bbarquark candidate" << std::endl;
      b2Coll.push_back(it->clone());
    }

  }// all Gen particles

  std::cout <<"size W1Coll "<< W1Coll.size() <<" W2Coll "<< W2Coll.size()<<" b1Coll "<< b1Coll.size()<<" b2Coll "<< b2Coll.size()<< std::endl;
  //htoWW
  if (W1Coll.size() && W2Coll.size()){
    for (auto W1_cand : W1Coll)
      for (auto W2_cand : W2Coll){
        const reco::Candidate* W1_mother = W1_cand->mother();
        const reco::Candidate* W2_mother = W2_cand->mother();
        while (W1_mother->pdgId() == -24) W1_mother = W1_mother->mother();
        while (W2_mother->pdgId() ==  24) W2_mother = W2_mother->mother();
        if (W1_mother == W2_mother && W1_mother->pdgId() == 25) {
          htoWWColl.push_back(W1_mother);
          htoWW = true;
          break;
        }
      }
  }
  //htoBB
  if (b1Coll.size() && b2Coll.size()){
    for(auto b1_cand : b1Coll)
      for (auto b2_cand : b2Coll) {
        const reco::Candidate* b1_mother = b1_cand->mother();
        const reco::Candidate* b2_mother = b2_cand->mother();
        if (b1_mother == b2_mother && b1_mother->pdgId() == 25) {
          htoBBColl.push_back(b1_mother);
          break;
        }
      }
  }
  //Xtohh
  if (htoWWColl.size() && htoBBColl.size()){
    for (auto htoWW_cand : htoWWColl){
      for (auto htoBB_cand : htoBBColl){
        const reco::Candidate* htoWW_mother = htoWW_cand->mother();
        const reco::Candidate* htoBB_mother = htoBB_cand->mother();
        while (htoWW_mother->pdgId() == 25)  htoWW_mother = htoWW_mother->mother();
        while (htoBB_mother->pdgId() == 25)  htoBB_mother = htoBB_mother->mother();
        //if (htoWW_mother == htoBB_mother && (htoBB_mother->pdgId()==99927 || htoBB_mother->pdgId()==99926)){ 
        //release the pdgId requirement to use other signal samples: graviton, radion
        if (htoWW_mother == htoBB_mother){ 
          Xtohhcand = htoWW_mother;
          htoWWcand = htoWW_cand;
          htoBBcand = htoBB_cand;
          Xtohh = true;
          break;
        }
      }
      if (Xtohh) break;
    }
    if (Xtohh)
      std::cout << "find X candidate, id "<< Xtohhcand->pdgId() <<" mass "<< Xtohhcand->mass() << std::endl;
    else std::cout <<"failed to find Xtohh "<< std::endl;
  }

  if (Xtohh){
    b1cand = findDecendant(htoBBcand, 5, false);
    b2cand = findDecendant(htoBBcand, -5, false);
    w1cand = findDecendant(htoWWcand, 24, false);
    w2cand = findDecendant(htoWWcand, -24, false);   
    if (debug_){
      if (hasDaughter(w1cand, -13) and hasDaughter(w2cand, 13)) std::cout <<" find two muons "<< std::endl;
      if (hasDaughter(w1cand, -11) and hasDaughter(w2cand, 11)) std::cout <<" find two electrons "<< std::endl;
      if ((hasDaughter(w1cand, -13) and hasDaughter(w2cand, 11)) or (hasDaughter(w1cand, -11) and hasDaughter(w2cand, 13))) std::cout <<" find two electron+muon "<< std::endl;
    }
    if ((hasDaughter(w1cand, -13) or hasDaughter(w1cand, 13) or hasDaughter(w1cand, -11) or hasDaughter(w1cand, 11))
      and (hasDaughter(w2cand, -13) or hasDaughter(w2cand, 13) or hasDaughter(w2cand, -11) or hasDaughter(w2cand, 11))){
      lep1cand = findLeptonDaughter(w1cand);
      nu1cand = findNuDaughter(w1cand);
      lep2cand = findLeptonDaughter(w2cand);
      nu2cand = findNuDaughter(w2cand);
      //make sure all candiates are in same frame
      while (lep1cand->numberOfDaughters()==1 and  lep1cand->daughter(0)->pdgId()==lep1cand->pdgId())
        lep1cand = lep1cand->daughter(0);
      while (nu1cand->numberOfDaughters()==1 and  nu1cand->daughter(0)->pdgId()==nu1cand->pdgId())
        nu1cand = nu1cand->daughter(0);
      while (lep2cand->numberOfDaughters()==1 and  lep2cand->daughter(0)->pdgId()==lep2cand->pdgId())
        lep2cand = lep2cand->daughter(0);
      while (nu2cand->numberOfDaughters()==1 and  nu2cand->daughter(0)->pdgId()==nu2cand->pdgId())
        nu2cand = nu2cand->daughter(0);
      findAllGenParticles = true;
    }else{
      findAllGenParticles = false;
      std::cout <<"failed to two leptons from  W decays "<< std::endl;
    }
    if (debug_){
      std::cout <<" w1 " ; printCandidate(w1cand);
      std::cout <<" w2 " ; printCandidate(w2cand);
      std::cout <<" b1 " ; printCandidate(b1cand);
      std::cout <<" b2 " ; printCandidate(b2cand);
      std::cout <<" lep1 " ; printCandidate(lep1cand);
      std::cout <<" nu1 " ; printCandidate(nu1cand);
      std::cout <<" lep2 " ; printCandidate(lep2cand);
      std::cout <<" nu2 " ; printCandidate(nu2cand);
    }
  }
  std::cout <<"*********** end in checking GenParticles for Signal sample ***********"<< std::endl;
}

// ------------ method called to check singal genParticles   ------------
void DiHiggsWWBBSLAnalyzer::checkGenParticlesZZbb(edm::Handle<reco::GenParticleCollection> genParticleColl){

  std::vector<reco::GenParticle*> b1Coll; 
  std::vector<reco::GenParticle*> b2Coll;
  std::vector<reco::GenParticle*> ZColl;
  std::vector<const reco::Candidate*> htoZZColl;
  std::vector<const reco::Candidate*> htoBBColl;

  std::cout <<"*********** start to check GenParticles for ZZBB sample ***********"<< std::endl;
  bool Xtohh = false;
  for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
    //particle id, (electron12)(muon13),(b5),(W+24),(SM higgs25), (Z 23)
    // particle id  it->pdgId()
    //std::cout << "Gen paticles: id " << it->pdgId() << std::endl; 
    //if (it->pdgId() == 24 && hasDaughter(it->clone(), -13) )
    if (it->pdgId() == 23){
	const reco::Candidate* tmpz = it->mother();
	//while (tmpz->pdgId() == 23 && tmpz->numberOfMothers() == 1) tmpz = tmpw1->mother();
	if (tmpz->pdgId() == 25)  ZColl.push_back(it->clone());
    }
    else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 )
    {
	if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
	//   std::cout << "find bquark candidate" << std::endl;
	b1Coll.push_back(it->clone());
    }
    else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 )
    {
	if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
	// std::cout << "find bbarquark candidate" << std::endl;
	b2Coll.push_back(it->clone());
    }

  }// all Gen particles

  std::cout <<"size ZColl "<< ZColl.size() <<" b1Coll "<< b1Coll.size()<<" b2Coll "<< b2Coll.size()<< std::endl;
  //htoWW
  if (ZColl.size()>=2){
    //for (auto Z1_cand : W1Coll)
    //	for (auto Z2_cand : W2Coll){
    for (unsigned i=0; i < ZColl.size()-1; i++){
	auto Z1_cand = ZColl.at(i);
	for (unsigned j=i+1; j < ZColl.size(); j++){
	  auto Z2_cand = ZColl.at(j);
	  const reco::Candidate* Z1_mother = Z1_cand->mother();
	  const reco::Candidate* Z2_mother = Z2_cand->mother();
	  while (Z1_mother->pdgId() == 23) Z1_mother = Z1_mother->mother();
	  while (Z2_mother->pdgId() == 23) Z2_mother = Z2_mother->mother();
	  if (Z1_mother == Z2_mother && Z1_mother->pdgId() == 25) {
	    htoZZColl.push_back(Z1_mother);
	    htoZZ = true;
	    break;
	  }
	}
    }
  }

  //htoBB
  if (b1Coll.size() && b2Coll.size()){
    for(auto b1_cand : b1Coll)
	for (auto b2_cand : b2Coll) {
	  const reco::Candidate* b1_mother = b1_cand->mother();
	  const reco::Candidate* b2_mother = b2_cand->mother();
	  if (b1_mother == b2_mother && b1_mother->pdgId() == 25) {
	    htoBBColl.push_back(b1_mother);
	    break;
	  }
	}
  }
  //Xtohh
  if (htoZZColl.size() && htoBBColl.size()){
    for (auto htoZZ_cand : htoZZColl){
	for (auto htoBB_cand : htoBBColl){
	  const reco::Candidate* htoZZ_mother = htoZZ_cand->mother();
	  const reco::Candidate* htoBB_mother = htoBB_cand->mother();
	  while (htoZZ_mother->pdgId() == 25)  htoZZ_mother = htoZZ_mother->mother();
	  while (htoBB_mother->pdgId() == 25)  htoBB_mother = htoBB_mother->mother();
	  //if (htoWW_mother == htoBB_mother && (htoBB_mother->pdgId()==99927 || htoBB_mother->pdgId()==99926)){ 
	  //release the pdgId requirement to use other signal samples: graviton, radion
	  if (htoZZ_mother == htoBB_mother){ 
	    Xtohhcand = htoZZ_mother;
	    htoZZcand = htoZZ_cand;
	    htoBBcand = htoBB_cand;
	    Xtohh = true;
	    break;
	  }
	}
	if (Xtohh) break;

    }
    if (Xtohh)
	std::cout << "find X candidate, id "<< Xtohhcand->pdgId() <<" mass "<< Xtohhcand->mass() << std::endl;
    else std::cout <<"failed to find Xtohh "<< std::endl;
  }

  int neutrino_id = 0;
  if (Xtohh){
    b1cand = findDecendant(htoBBcand, 5, false);
    b2cand = findDecendant(htoBBcand, -5, false);
    //not working here 
    const reco::Candidate* htoZZtmp = htoZZcand;
    while (htoZZtmp->numberOfDaughters() == 1 and htoZZtmp->daughter(0)->pdgId() == 25) htoZZtmp = htoZZtmp->daughter(0);
    for (unsigned i=0; i < htoZZtmp->numberOfDaughters(); i++){
    	const reco::Candidate* tmp = htoZZtmp->daughter(i);
	if (tmp->pdgId() == 23 and hasDaughter(tmp, -13) and hasDaughter(tmp, 13))
	    z1cand = tmp;
	//else if (tmp->pdgId() == 23 and hasDaughter(tmp, -11) and hasDaughter(tmp, 11))
	else if (tmp->pdgId() == 23 and hasDaughter(tmp, -12) and hasDaughter(tmp, 12)){
	    neutrino_id = 12 ;
	    z2cand = tmp;
	}else if (tmp->pdgId() == 23 and hasDaughter(tmp, -14) and hasDaughter(tmp, 14)){
	    neutrino_id = 14 ;
	    z2cand = tmp;
	}else if (tmp->pdgId() == 23 and hasDaughter(tmp, -16) and hasDaughter(tmp, 16)){
	    neutrino_id = 16 ;
	    z2cand = tmp;
	}

    }

    if (not (z1cand and z2cand)){
	std::cout <<"Failed to two Zs decaying mu mu and nu nu"<< std::endl;
	return ;
    }


    if (debug_){
	if (hasDaughter(z1cand, -13) and hasDaughter(z1cand, 13)) std::cout <<" find two muons "<< std::endl;
	if (hasDaughter(z1cand, -11) and hasDaughter(z1cand, 11)) std::cout <<" find two electrons "<< std::endl;
	//if ((hasDaughter(zcand, -13) and hasDaughter(w2cand, 11)) or (hasDaughter(w1cand, -11) and hasDaughter(w2cand, 13))) std::cout <<" find two electron+muon "<< std::endl;
    }
    if (hasDaughter(z1cand, -13) and hasDaughter(z1cand, 13) and z2cand){
	lep1cand = stableDecendant(z1cand, -13);
	lep2cand = stableDecendant(z1cand, 13);
	nu1cand = stableDecendant(z2cand, neutrino_id);
	nu2cand = stableDecendant(z2cand, neutrino_id*(-1));
	//make sure all candiates are in same frame
	while (lep1cand->numberOfDaughters()==1 and  lep1cand->daughter(0)->pdgId()==lep1cand->pdgId())
	  lep1cand = lep1cand->daughter(0);
	while (nu1cand->numberOfDaughters()==1 and  nu1cand->daughter(0)->pdgId()==nu1cand->pdgId())
	  nu1cand = nu1cand->daughter(0);
	while (lep2cand->numberOfDaughters()==1 and  lep2cand->daughter(0)->pdgId()==lep2cand->pdgId())
	  lep2cand = lep2cand->daughter(0);
	while (nu2cand->numberOfDaughters()==1 and  nu2cand->daughter(0)->pdgId()==nu2cand->pdgId())
	  nu2cand = nu2cand->daughter(0);
	findAllGenParticles = true;
	std::cout <<" lep1 " ; printCandidate(lep1cand);
	std::cout <<" nu1 " ; printCandidate(nu1cand);
	std::cout <<" lep2 " ; printCandidate(lep2cand);
	std::cout <<" nu2 " ; printCandidate(nu2cand);
    }else{
	findAllGenParticles = false;
	std::cout <<"failed to two muons from  Z decays "<< std::endl;
    }
    if (debug_){
	std::cout <<" z1 " ; printCandidate(z1cand);
	std::cout <<" z2 " ; printCandidate(z2cand);
	std::cout <<" b1 " ; printCandidate(b1cand);
	std::cout <<" b2 " ; printCandidate(b2cand);
    }
  }
  std::cout <<"*********** end in checking GenParticles for ZZbb sample ***********"<< std::endl;
}


// ------------ method called to check singal genParticles   ------------
void DiHiggsWWBBSLAnalyzer::checkGenParticlesHbbWW2J2Nu(edm::Handle<reco::GenParticleCollection> genParticleColl){

  std::vector<reco::GenParticle*> b1Coll; 
  std::vector<reco::GenParticle*> b2Coll;
  std::vector<reco::GenParticle*> W1Coll;
  std::vector<reco::GenParticle*> W2Coll;
  std::vector<const reco::Candidate*> htoWWColl;
  std::vector<const reco::Candidate*> htoBBColl;

  std::cout <<"*********** start to check GenParticles for HbbWW SL sample ***********"<< std::endl;
  bool Xtohh = false;
  for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
    //particle id, (electron12)(muon13),(b5),(W+24),(SM higgs25)
    // particle id  it->pdgId()
    //std::cout << "Gen paticles: id " << it->pdgId() << std::endl; 
    //if (it->pdgId() == 24 && hasDaughter(it->clone(), -13) )
    if (it->pdgId() == 24){
      const reco::Candidate* tmpw1 = it->mother();
      while (tmpw1->pdgId() == 24 && tmpw1->numberOfMothers() == 1) tmpw1 = tmpw1->mother();
      if (tmpw1->pdgId() == 25)  W1Coll.push_back(it->clone());
    }
    else if (it->pdgId() == -24){
      const reco::Candidate* tmpw2 = it->mother();
      while (tmpw2->pdgId() == -24 && tmpw2->numberOfMothers() == 1) tmpw2 = tmpw2->mother();
      if (tmpw2->pdgId() == 25)  W2Coll.push_back(it->clone());
    }
    else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 )
    {
      if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
      //   std::cout << "find bquark candidate" << std::endl;
      b1Coll.push_back(it->clone());
    }
    else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 )
    {
      if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
      // std::cout << "find bbarquark candidate" << std::endl;
      b2Coll.push_back(it->clone());
    }

  }// all Gen particles

  std::cout <<"size W1Coll "<< W1Coll.size() <<" W2Coll "<< W2Coll.size()<<" b1Coll "<< b1Coll.size()<<" b2Coll "<< b2Coll.size()<< std::endl;
  //htoWW, find the last copy
  if (W1Coll.size() && W2Coll.size()){
    for (auto W1_cand : W1Coll)
      for (auto W2_cand : W2Coll){
        const reco::Candidate* W1_mother = W1_cand->mother();
        const reco::Candidate* W2_mother = W2_cand->mother();
        while (W1_mother->pdgId() == -24) W1_mother = W1_mother->mother();
        while (W2_mother->pdgId() ==  24) W2_mother = W2_mother->mother();
        if (W1_mother == W2_mother && W1_mother->pdgId() == 25) {
          htoWWColl.push_back(W1_mother);
          htoWW = true;
          break;
        }
      }
  }
  //htoBB
  if (b1Coll.size() && b2Coll.size()){
    for(auto b1_cand : b1Coll)
      for (auto b2_cand : b2Coll) {
        const reco::Candidate* b1_mother = b1_cand->mother();
        const reco::Candidate* b2_mother = b2_cand->mother();
        if (b1_mother == b2_mother && b1_mother->pdgId() == 25) {
          htoBBColl.push_back(b1_mother);
          break;
        }
      }
  }
  //Xtohh
  if (htoWWColl.size() && htoBBColl.size()){
    for (auto htoWW_cand : htoWWColl){
      for (auto htoBB_cand : htoBBColl){
        const reco::Candidate* htoWW_mother = htoWW_cand->mother();
        const reco::Candidate* htoBB_mother = htoBB_cand->mother();
        while (htoWW_mother->pdgId() == 25)  htoWW_mother = htoWW_mother->mother();
        while (htoBB_mother->pdgId() == 25)  htoBB_mother = htoBB_mother->mother();
        //if (htoWW_mother == htoBB_mother && (htoBB_mother->pdgId()==99927 || htoBB_mother->pdgId()==99926)){ 
        //release the pdgId requirement to use other signal samples: graviton, radion
        if (htoWW_mother == htoBB_mother){ 
          Xtohhcand = htoWW_mother;
          htoWWcand = htoWW_cand;
          htoBBcand = htoBB_cand;
          Xtohh = true;
          break;
        }
      }
      if (Xtohh) break;

    }
    if (Xtohh)
      std::cout << "find X candidate, id "<< Xtohhcand->pdgId() <<" mass "<< Xtohhcand->mass() << std::endl;
    else std::cout <<"failed to find Xtohh "<< std::endl;
  }

  if (Xtohh){
    b1cand = findDecendant(htoBBcand, 5, false);
    b2cand = findDecendant(htoBBcand, -5, false);
    w1cand = findDecendant(htoWWcand, 24, false);
    w2cand = findDecendant(htoWWcand, -24, false);   
    if (debug_){
      if (hasDaughter(w1cand, -13) and hasDaughter(w2cand, 13)) std::cout <<" find two muons from H->ww"<< std::endl;
      else if(hasDaughter(w1cand, -13) or hasDaughter(w2cand, 13)) std::cout <<" find one muon from H->ww"<< std::endl; 
      if (hasDaughter(w1cand, -11) and hasDaughter(w2cand, 11)) std::cout <<" find two electrons from H->ww"<< std::endl;
      else if(hasDaughter(w1cand, -11) and hasDaughter(w2cand, 11)) std::cout <<" find one electron from H->ww"<< std::endl; 
      if ((hasDaughter(w1cand, -13) and hasDaughter(w2cand, 11)) or (hasDaughter(w1cand, -11) and hasDaughter(w2cand, 13))) std::cout <<" find two electron+muon from H->ww"<< std::endl;
    }
    if (hasDaughter(w1cand, -11) or hasDaughter(w1cand, 11) or hasDaughter(w1cand, -13) or hasDaughter(w1cand, 13)){
      lep1cand = findLeptonDaughter(w1cand);
      nu1cand = findNuDaughter(w1cand);
      //make sure all candiates are in same frame
      while (lep1cand->numberOfDaughters()==1 and  lep1cand->daughter(0)->pdgId()==lep1cand->pdgId())
        lep1cand = lep1cand->daughter(0);
      while (nu1cand->numberOfDaughters()==1 and  nu1cand->daughter(0)->pdgId()==nu1cand->pdgId())
        nu1cand = nu1cand->daughter(0);
      //while (lep2cand->numberOfDaughters()==1 and  lep2cand->daughter(0)->pdgId()==lep2cand->pdgId())
      //  lep2cand = lep2cand->daughter(0);
      //while (nu2cand->numberOfDaughters()==1 and  nu2cand->daughter(0)->pdgId()==nu2cand->pdgId())
      //  nu2cand = nu2cand->daughter(0);
      //W->q1q2
      if (w2cand->numberOfDaughters()==2){
        q1cand = w2cand->daughter(0);
        q2cand = w2cand->daughter(1);
        if (abs(q1cand->pdgId()) <= 4 and abs(q2cand->pdgId()) <= 4){
          std::cout <<"one W leptonically decay and the other W decays into light quarks" << std::endl;
        }else if (abs(q1cand->pdgId()) >=5 or abs(q2cand->pdgId()) >=5){
          std::cout <<"one W leptonically decay and the other W decays into b/t " << std::endl;
        }else if ((abs(q1cand->pdgId()) >= 11 and abs(q1cand->pdgId()) <= 16)){
          std::cout <<"one W leptonically decay and the other W decays into leptons " << std::endl;
        }else{
          std::cout <<"one W leptonically decay and the other W decays into unknowns" << std::endl;
          std::cout <<"w2 "; printCandidate(w2cand); 
          std::cout <<"q1 "; printCandidate(q1cand); 
          std::cout <<"q2 "; printCandidate(q2cand); 
        }
        while(q1cand->numberOfDaughters()==1 and  q1cand->daughter(0)->pdgId()==q1cand->pdgId())
          q1cand = q1cand->daughter(0);
        while(q2cand->numberOfDaughters()==1 and  q2cand->daughter(0)->pdgId()==q2cand->pdgId())
          q2cand = q2cand->daughter(0);
      }else{
        std::cout <<"one W leptonically decay and the num of daugters from the other W "<<w2cand->numberOfDaughters() << std::endl;
        return ;
      }
      findAllGenParticles = true;
    }else if (hasDaughter(w2cand, -11) or hasDaughter(w2cand, 11) or hasDaughter(w2cand, -13) or hasDaughter(w2cand, 13)){
      lep1cand = findLeptonDaughter(w2cand);
      nu1cand = findNuDaughter(w2cand);
      //make sure all candiates are in same frame
      while (lep1cand->numberOfDaughters()==1 and  lep1cand->daughter(0)->pdgId()==lep1cand->pdgId())
        lep1cand = lep1cand->daughter(0);
      while (nu1cand->numberOfDaughters()==1 and  nu1cand->daughter(0)->pdgId()==nu1cand->pdgId())
        nu1cand = nu1cand->daughter(0);
      //while (lep2cand->numberOfDaughters()==1 and  lep2cand->daughter(0)->pdgId()==lep2cand->pdgId())
      //  lep2cand = lep2cand->daughter(0);
      //while (nu2cand->numberOfDaughters()==1 and  nu2cand->daughter(0)->pdgId()==nu2cand->pdgId())
      //  nu2cand = nu2cand->daughter(0);
      //W->q1q2
      if (w1cand->numberOfDaughters()==2){
        q1cand = w1cand->daughter(0);
        q2cand = w1cand->daughter(1);
        if (abs(q1cand->pdgId()) <= 4 and abs(q2cand->pdgId()) <= 4){
          std::cout <<"one W leptonically decay and the other W decays into light quarks" << std::endl;
        }else if (abs(q1cand->pdgId()) >=5 or abs(q2cand->pdgId()) >=5){
          std::cout <<"one W leptonically decay and the other W decays into b/t " << std::endl;
        }else if ((abs(q1cand->pdgId()) >= 11 and abs(q1cand->pdgId()) <= 16)){
          std::cout <<"one W leptonically decay and the other W decays into leptons " << std::endl;
        }else{
          std::cout <<"one W leptonically decay and the other W decays into unknowns" << std::endl;
          std::cout <<"w2 "; printCandidate(w2cand); 
          std::cout <<"q1 "; printCandidate(q1cand); 
          std::cout <<"q2 "; printCandidate(q2cand); 
        }
        while(q1cand->numberOfDaughters()==1 and  q1cand->daughter(0)->pdgId()==q1cand->pdgId())
          q1cand = q1cand->daughter(0);
        while(q2cand->numberOfDaughters()==1 and  q2cand->daughter(0)->pdgId()==q2cand->pdgId())
          q2cand = q2cand->daughter(0);
      }else{
        std::cout <<"one W leptonically decay and the num of daugters from the other W "<<w1cand->numberOfDaughters() << std::endl;
        return ;
      }
      findAllGenParticles = true;
    }else{
      findAllGenParticles = false;
      std::cout <<"failed to two leptons from  W decay and two quarks from the other W "<< std::endl;
    }
    if (debug_){
      std::cout <<" w1 " ; printCandidate(w1cand);
      std::cout <<" w2 " ; printCandidate(w2cand);
      std::cout <<" b1 " ; printCandidate(b1cand);
      std::cout <<" b2 " ; printCandidate(b2cand);
      std::cout <<" lep1 " ; printCandidate(lep1cand);
      std::cout <<" nu1 " ; printCandidate(nu1cand);
      std::cout <<" quark1 " ; printCandidate(q1cand);
      std::cout <<" quark2 " ; printCandidate(q2cand);
    }
  }//end of Xtohh
  std::cout <<"*********** end in checking GenParticles for Signal sample ***********"<< std::endl;
}


// ------------ method called to check genParticles ZZbb->2L2J2B   ------------
void DiHiggsWWBBSLAnalyzer::checkGenParticlesZZbb_2L2J(edm::Handle<reco::GenParticleCollection> genParticleColl){

  std::vector<reco::GenParticle*> b1Coll; 
  std::vector<reco::GenParticle*> b2Coll;
  std::vector<reco::GenParticle*> ZColl;
  std::vector<const reco::Candidate*> htoZZColl;
  std::vector<const reco::Candidate*> htoBBColl;

  std::cout <<"*********** start to check GenParticles for ZZBB, ZZ->2L2J sample ***********"<< std::endl;
  bool Xtohh = false;
  for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
    //particle id, (electron12)(muon13),(b5),(W+24),(SM higgs25), (Z 23)
    // particle id  it->pdgId()
    //std::cout << "Gen paticles: id " << it->pdgId() << std::endl; 
    //if (it->pdgId() == 24 && hasDaughter(it->clone(), -13) )
    if (it->pdgId() == 23){
	const reco::Candidate* tmpz = it->mother();
	//while (tmpz->pdgId() == 23 && tmpz->numberOfMothers() == 1) tmpz = tmpw1->mother();
	if (tmpz->pdgId() == 25)  ZColl.push_back(it->clone());
    }
    else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 )
    {
	if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
	//   std::cout << "find bquark candidate" << std::endl;
	b1Coll.push_back(it->clone());
    }
    else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 )
    {
	if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
	// std::cout << "find bbarquark candidate" << std::endl;
	b2Coll.push_back(it->clone());
    }

  }// all Gen particles

  std::cout <<"size ZColl "<< ZColl.size() <<" b1Coll "<< b1Coll.size()<<" b2Coll "<< b2Coll.size()<< std::endl;
  //htoWW
  if (ZColl.size()>=2){
    //for (auto Z1_cand : W1Coll)
    //	for (auto Z2_cand : W2Coll){
    for (unsigned i=0; i < ZColl.size()-1; i++){
	auto Z1_cand = ZColl.at(i);
	for (unsigned j=i+1; j < ZColl.size(); j++){
	  auto Z2_cand = ZColl.at(j);
	  const reco::Candidate* Z1_mother = Z1_cand->mother();
	  const reco::Candidate* Z2_mother = Z2_cand->mother();
	  while (Z1_mother->pdgId() == 23) Z1_mother = Z1_mother->mother();
	  while (Z2_mother->pdgId() == 23) Z2_mother = Z2_mother->mother();
	  if (Z1_mother == Z2_mother && Z1_mother->pdgId() == 25) {
	    htoZZColl.push_back(Z1_mother);
	    htoZZ = true;
	    break;
	  }
	}
    }
  }

  //htoBB
  if (b1Coll.size() && b2Coll.size()){
    for(auto b1_cand : b1Coll)
	for (auto b2_cand : b2Coll) {
	  const reco::Candidate* b1_mother = b1_cand->mother();
	  const reco::Candidate* b2_mother = b2_cand->mother();
	  if (b1_mother == b2_mother && b1_mother->pdgId() == 25) {
	    htoBBColl.push_back(b1_mother);
	    break;
	  }
	}
  }
  //Xtohh
  if (htoZZColl.size() && htoBBColl.size()){
    for (auto htoZZ_cand : htoZZColl){
	for (auto htoBB_cand : htoBBColl){
	  const reco::Candidate* htoZZ_mother = htoZZ_cand->mother();
	  const reco::Candidate* htoBB_mother = htoBB_cand->mother();
	  while (htoZZ_mother->pdgId() == 25)  htoZZ_mother = htoZZ_mother->mother();
	  while (htoBB_mother->pdgId() == 25)  htoBB_mother = htoBB_mother->mother();
	  //if (htoWW_mother == htoBB_mother && (htoBB_mother->pdgId()==99927 || htoBB_mother->pdgId()==99926)){ 
	  //release the pdgId requirement to use other signal samples: graviton, radion
	  if (htoZZ_mother == htoBB_mother){ 
	    Xtohhcand = htoZZ_mother;
	    htoZZcand = htoZZ_cand;
	    htoBBcand = htoBB_cand;
	    Xtohh = true;
	    break;
	  }
	}
	if (Xtohh) break;

    }
    if (Xtohh)
	std::cout << "find X candidate, id "<< Xtohhcand->pdgId() <<" mass "<< Xtohhcand->mass() << std::endl;
    else std::cout <<"failed to find Xtohh "<< std::endl;
  }

  int quarkid = 0;
  if (Xtohh){
    b1cand = findDecendant(htoBBcand, 5, false);
    b2cand = findDecendant(htoBBcand, -5, false);
   
    const reco::Candidate* htoZZtmp = htoZZcand;
    while (htoZZtmp->numberOfDaughters() == 1 and htoZZtmp->daughter(0)->pdgId() == 25) htoZZtmp = htoZZtmp->daughter(0);
    for (unsigned i=0; i < htoZZtmp->numberOfDaughters(); i++){
    	const reco::Candidate* tmp = htoZZtmp->daughter(i);
	if (tmp->pdgId() == 23 and hasDaughter(tmp, -13) and hasDaughter(tmp, 13))
	    z1cand = tmp;
	//light quarks(not including b quark)
	else if (tmp->pdgId() == 23 and hasDaughter(tmp, -1) and hasDaughter(tmp, 1)){
	    quarkid = 1;
	    z2cand = tmp;
	}else if (tmp->pdgId() == 23 and hasDaughter(tmp, -2) and hasDaughter(tmp, 2)){
	    quarkid = 2;
	    z2cand = tmp;
	}else if (tmp->pdgId() == 23 and hasDaughter(tmp, -3) and hasDaughter(tmp, 3)){
	    quarkid = 3;
	    z2cand = tmp;
	}else if (tmp->pdgId() == 23 and hasDaughter(tmp, -4) and hasDaughter(tmp, 4)){
	    quarkid = 4;
	    z2cand = tmp;
	}else {
	    std::cout <<"Z->JJ, light jet is not from quarkid = 1,2,3,4 " <<std::endl;
	    z2cand = NULL;
	    return ;
	}

    }

    if (not (z1cand and z2cand)){
	std::cout <<"Failed to two Zs decaying mu mu and nu nu"<< std::endl;
	return ;
    }


    if (debug_){
	if (hasDaughter(z1cand, -13) and hasDaughter(z1cand, 13)) std::cout <<" find two muons "<< std::endl;
	if (hasDaughter(z1cand, -11) and hasDaughter(z1cand, 11)) std::cout <<" find two electrons "<< std::endl;
	std::cout <<"Z->JJ quark id "<< quarkid << std::endl;
	//if ((hasDaughter(zcand, -13) and hasDaughter(w2cand, 11)) or (hasDaughter(w1cand, -11) and hasDaughter(w2cand, 13))) std::cout <<" find two electron+muon "<< std::endl;
    }
    if (hasDaughter(z1cand, -13) and hasDaughter(z1cand, 13) and quarkid>0){
	lep1cand = stableDecendant(z1cand, -13);
	lep2cand = stableDecendant(z1cand, 13);
        nu1cand = findDecendant(z2cand, quarkid, true);
        nu2cand = findDecendant(z2cand, quarkid*(-1), true);
	//make sure all candiates are in same frame
	while (lep1cand->numberOfDaughters()==1 and  lep1cand->daughter(0)->pdgId()==lep1cand->pdgId())
	  lep1cand = lep1cand->daughter(0);
	//while (nu1cand->numberOfDaughters()==1 and  nu1cand->daughter(0)->pdgId()==nu1cand->pdgId())
	//  nu1cand = nu1cand->daughter(0);
	while (lep2cand->numberOfDaughters()==1 and  lep2cand->daughter(0)->pdgId()==lep2cand->pdgId())
	  lep2cand = lep2cand->daughter(0);
	//while (nu2cand->numberOfDaughters()==1 and  nu2cand->daughter(0)->pdgId()==nu2cand->pdgId())
	//  nu2cand = nu2cand->daughter(0);
	findAllGenParticles = true;
	std::cout <<" lep1 " ; printCandidate(lep1cand);
	std::cout <<" nu1 " ; printCandidate(nu1cand);
	std::cout <<" lep2 " ; printCandidate(lep2cand);
	std::cout <<" nu2 " ; printCandidate(nu2cand);
    }else{
	findAllGenParticles = false;
	std::cout <<"failed to two muons from  Z decays "<< std::endl;
    }
    if (debug_){
	std::cout <<" z1 " ; printCandidate(z1cand);
	std::cout <<" z2 " ; printCandidate(z2cand);
	std::cout <<" b1 " ; printCandidate(b1cand);
	std::cout <<" b2 " ; printCandidate(b2cand);
    }
  }
  std::cout <<"*********** end in checking GenParticles for ZZbb, ZZ->2L2J sample ***********"<< std::endl;
}




// ------------ method called to check ttbar genParticles   ------------
void DiHiggsWWBBSLAnalyzer::checkGenParticlesTTbar(edm::Handle<reco::GenParticleCollection> genParticleColl){

  std::vector<reco::GenParticle*> b1Coll; 
  std::vector<reco::GenParticle*> b2Coll;
  std::vector<reco::GenParticle*> W1Coll;
  std::vector<reco::GenParticle*> W2Coll;
  std::vector<const reco::Candidate*> tColl;
  std::vector<const reco::Candidate*> tbarColl;
  std::cout <<"*********** start to check GenParticles for TTbar sample ***********"<< std::endl;
  for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
    //particle id, (electron12)(muon13),(b5),(W+24),(SM higgs25)
    //particle id  it->pdgId()
    //std::cout << "Gen paticles: id " << it->pdgId() << std::endl; 
    //if (it->pdgId() == 24 && hasDaughter(it->clone(), -13) )
    if (it->pdgId() == 24 ){
	const reco::Candidate* tmpw1 = it->mother();
	//while (tmpw1->pdgId() == 24 && tmpw1->numberOfMothers() == 1) tmpw1 = tmpw1->mother();
	if (tmpw1->pdgId() == 6)  W1Coll.push_back(it->clone());
    }
    else if (it->pdgId() == -24 ){
	const reco::Candidate* tmpw2 = it->mother();
	//while (tmpw2->pdgId() == -24 && tmpw2->numberOfMothers() == 1) tmpw2 = tmpw2->mother();
	if (tmpw2->pdgId() == -6)  W2Coll.push_back(it->clone());
    }
    else if (it->pdgId() == 5 && it->mother()->pdgId() == 6 ){
	if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
	b1Coll.push_back(it->clone());
    }
    else if (it->pdgId() == -5 && it->mother()->pdgId() == -6 ){
	if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
	b2Coll.push_back(it->clone());
    }
  }// all Gen particles

  if(debug_) std::cout <<"size W1Coll "<< W1Coll.size() <<" W2Coll "<< W2Coll.size()<<" b1Coll "<< b1Coll.size()<<" b2Coll "<< b2Coll.size()<< std::endl;
  //t->bW+ 6->5,24
  if (W1Coll.size() && b1Coll.size()){
    for (auto W1_cand : W1Coll)
	for (auto b1_cand : b1Coll){
	  const reco::Candidate* b1_mother = b1_cand->mother();
	  const reco::Candidate* W1_mother = W1_cand->mother();
	  //        while (W1_mother->pdgId() == 24) W1_mother = W1_mother->mother();
	  if (W1_mother == b1_mother && W1_mother->pdgId() == 6) {
	    tColl.push_back(W1_mother);
	    //std::cout <<" find t->bW+" << std::endl;
	    break;
	  }
	}
  }
  //tbar->bbarW- -6->-5,-24
  if (W2Coll.size() && b2Coll.size()){
    for(auto W2_cand : W2Coll)
	for (auto b2_cand : b2Coll) {
	  const reco::Candidate* W2_mother = W2_cand->mother();
	  const reco::Candidate* b2_mother = b2_cand->mother();
	  //   while (W2_mother->pdgId() == -24) W2_mother = W2_mother->mother();
	  if (W2_mother == b2_mother && W2_mother->pdgId() == -6) {
	    tbarColl.push_back(W2_mother);
	    //std::cout <<" find tbar->bbarW- " << std::endl;
	    break;
	  }
	}
  }
  std::cout <<"tColl size " << tColl.size() <<" tbarColl " << tbarColl.size() << std::endl;
  //bool ttbar =false;
  if (tColl.size()==1 and tbarColl.size() ==1){
    //ttbar = true;
    t1cand = tColl.at(0);
    t2cand = tbarColl.at(0);
    //while (t1cand->numberOfMothers()==1 and (t1cand->mother())->pdgId()==6) t1cand = t1cand->mother();
    //while (t2cand->numberOfMothers()==1 and (t2cand->mother())->pdgId()==-6) t2cand = t2cand->mother();
    std::cout << "find t and tbar candidate " << std::endl;
    b1cand = findDecendant(t1cand, 5, false);
    w1cand = findDecendant(t1cand, 24, false);
    b2cand = findDecendant(t2cand, -5, false);
    w2cand = findDecendant(t2cand, -24, false);   
    //for (unsigned int i=0; i < b1cand->numberOfDaughters(); i++)
    //std::cout <<"candidate id "<< b1cand->pdgId()<<" daughter i "<< i <<" id "<<(b1cand->daughter(i))->pdgId()<< std::endl;
    if (hasDaughter(w1cand, -13) and hasDaughter(w2cand, 13)){
	lep1cand = findLeptonDaughter(w1cand);
	nu1cand = findNuDaughter(w1cand);
	lep2cand = findLeptonDaughter(w2cand);
	nu2cand = findNuDaughter(w2cand);
	//make sure all candiates are in same frame
	while (lep1cand->numberOfDaughters()==1 and  lep1cand->daughter(0)->pdgId()==lep1cand->pdgId())
	  lep1cand = lep1cand->daughter(0);
	while (nu1cand->numberOfDaughters()==1 and  nu1cand->daughter(0)->pdgId()==nu1cand->pdgId())
	  nu1cand = nu1cand->daughter(0);
	while (lep2cand->numberOfDaughters()==1 and  lep2cand->daughter(0)->pdgId()==lep2cand->pdgId())
	  lep2cand = lep2cand->daughter(0);
	while (nu2cand->numberOfDaughters()==1 and  nu2cand->daughter(0)->pdgId()==nu2cand->pdgId())
	  nu2cand = nu2cand->daughter(0);
	findAllGenParticles = true;
	ttbar = true;
	std::cout <<" lep1 " ; printCandidate(lep1cand);
	std::cout <<" nu1 " ; printCandidate(nu1cand);
	std::cout <<" lep2 " ; printCandidate(lep2cand);
	std::cout <<" nu2 " ; printCandidate(nu2cand);
    }else{
	findAllGenParticles = false;
	std::cout <<"failed to two leptons from  W decays "<< std::endl;
    }
    if (debug_){
	std::cout <<" w1 " ; printCandidate(w1cand);
	std::cout <<" w2 " ; printCandidate(w2cand);
	std::cout <<" b1 " ; printCandidate(b1cand);
	std::cout <<" b2 " ; printCandidate(b2cand);
    }
  }
  std::cout <<"*********** end in checking GenParticles for TTbar sample ***********"<< std::endl;
}

// ------------ method called to check DY genParticles   ------------
void DiHiggsWWBBSLAnalyzer::checkGenParticlesDY(edm::Handle<reco::GenParticleCollection> genParticleColl){

  std::vector<reco::GenParticle*> lept1Coll;
  std::vector<reco::GenParticle*> lept2Coll;
  std::vector<reco::GenParticle*> bJet1Coll;
  std::vector<reco::GenParticle*> bJet2Coll;
  std::vector<const reco::Candidate*> ZColl;

  std::cout <<"*********** start to check GenParticles for DY sample ***********"<< std::endl;
  for(reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
    if( debug_ and abs(it->pdgId()) == 13){ 
	std::cout<<"ID: "<<it->pdgId()<<" ->numOfmothers: "<<it->numberOfMothers()<<std::endl;
	//if( it->numberOfMothers()==0 ) std::cout<<"  DY: NO MOTHER!"<<std::endl;
	//if( it->numberOfMothers()==1 ) std::cout<<"  DY: Mother is: "<<it->mother()->pdgId()<<std::endl;
	//if( it->numberOfMothers()==2 ) std::cout<<"  DY: MORE THAN 2 MOTHERS. "<<(it->mother(0))->pdgId()<<" "<<(it->mother(1))->pdgId()<<std::endl;
	if( it->numberOfMothers()==1  and (it->mother()->pdgId()==22 or it->mother()->pdgId()==23))
	  std::cout <<"found muons from gamma or Z, id "<< it->mother()->pdgId()<<" muon status "<< it->status()<<" mother's pt "<< it->mother()->pt() <<" mass "<< it->mother()->mass() << std::endl;
    }
    if( it->pdgId() == -13 and (it->mother()->pdgId()==22 or it->mother()->pdgId()==23))      lept1Coll.push_back(it->clone());
    if( it->pdgId() == 13 and (it->mother()->pdgId()==22 or it->mother()->pdgId()==23))       lept2Coll.push_back(it->clone());
    if( fabs(it->pdgId()) == 5 )  bJet1Coll.push_back(it->clone());
    if( fabs(it->pdgId()) == -5 ) bJet2Coll.push_back(it->clone());
  }// all Gen particles
  //Fill Final quantites, only Z/gamma->ll are accepted 
  if(lept1Coll.size()>0 and lept2Coll.size()>0){
    for (auto lept1Cand : lept1Coll){
	for (auto lept2Cand : lept2Coll){
	  const reco::Candidate *tmpcand1 = lept1Cand->mother();
	  const reco::Candidate *tmpcand2 = lept2Cand->mother();
	  if (tmpcand1 == tmpcand2){
	    ZColl.push_back(tmpcand1);
	    TLorentzVector part(tmpcand1->px(), tmpcand1->py(), tmpcand1->pz(), tmpcand1->energy());
	    std::cout <<"find Z/gamma to ll , pdgId "<< tmpcand1->pdgId() <<" mass "<< part.M()<<" pt "<< part.Pt() << std::endl;
	  }
	}
    }
  }

  if (ZColl.size()>0){
    lep1cand = stableDecendant(ZColl[0], -13);
    lep2cand = stableDecendant(ZColl[0], 13);
  }

  //not care about where b comes from 
  if(bJet1Coll.size()>0){
    float minPT=5;//min pt cut 
    for (auto bJet1Cand : bJet1Coll){
	const reco::Candidate *tmpcand = bJet1Cand;
	TLorentzVector part(tmpcand->px(), tmpcand->py(), tmpcand->pz(), tmpcand->energy());
	if(part.Pt()>minPT){
	  b1cand = tmpcand;
	  minPT = part.Pt();
	}
    }
  }
  if(bJet2Coll.size()>0){
    float minPT=5;
    for (auto bJet2Cand : bJet2Coll){
	const reco::Candidate *tmpcand = bJet2Cand;
	TLorentzVector part(tmpcand->px(), tmpcand->py(), tmpcand->pz(), tmpcand->energy());
	if(part.Pt()>minPT){
	  b2cand = tmpcand;
	  minPT = part.Pt();
	}
    }
  }
  if(ZColl.size()>0 and bJet1Coll.size()>0 and bJet2Coll.size()>0) findAllGenParticles=true;
  //std::cout <<"size lept1Coll "<< lept1Coll.size()<<" size lept2Coll "<< lept2Coll.size()<<" bJet1Coll "<< bJet1Coll.size()<<" bJet2Coll "<< bJet2Coll.size()<<" -> "<<findAllGenParticles<< std::endl;
  std::cout <<"*********** end in checking GenParticles for DY sample ***********"<< std::endl;
}

void DiHiggsWWBBSLAnalyzer::matchGenJet2Parton(edm::Handle<std::vector<reco::GenJet>> genjetColl){
  TLorentzVector b1_p4(b1cand->px(), b1cand->py(), b1cand->pz(), b1cand->energy());
  TLorentzVector b2_p4(b2cand->px(), b2cand->py(), b2cand->pz(), b2cand->energy());
  TLorentzVector q1_p4(q1cand->px(), q1cand->py(), q1cand->pz(), q1cand->energy());
  TLorentzVector q2_p4(q2cand->px(), q2cand->py(), q2cand->pz(), q2cand->energy());
  // float dR_b1genjet = 9;
  // float dR_b2genjet = 9;

  bool hasb1genjet = false, hasb2genjet = false;
  bool hasq1genjet = false, hasq2genjet = false;
  for (const reco::GenJet& genjet: *genjetColl){
    TLorentzVector genjet_p4(genjet.px(), genjet.py(), genjet.pz(), genjet.energy());
    float dR_b1 = genjet_p4.DeltaR(b1_p4);
    float dR_b2 = genjet_p4.DeltaR(b2_p4);
    float dR_q1 = genjet_p4.DeltaR(q1_p4);
    float dR_q2 = genjet_p4.DeltaR(q2_p4);
    //std::cout <<"dR_b1 "<<dR_b1 <<" dR_b2 "<< dR_b2 <<" genjet px "<< genjet.px() << " py "<< genjet.py()<<" pz "<< genjet.pz()<<" E "<<genjet.energy() <<std::endl;
    if (dR_b1 < dR_b2 and dR_b1 < dR_b1genjet){
      b1genjet = &genjet;
      dR_b1genjet = dR_b1;
      hasb1genjet = true;
    }else if (dR_b1 > dR_b2 and dR_b2 < dR_b2genjet){
      b2genjet = &genjet;
      dR_b2genjet = dR_b2;
      hasb2genjet = true;
    }

    if (dR_q1 < dR_q2 and dR_q1 < dR_q1genjet){
      q1genjet = &genjet;
      dR_q1genjet = dR_q1;
      hasq1genjet = true;
    }else if (dR_q1 > dR_q2 and dR_q2 < dR_q2genjet){
      q2genjet = &genjet;
      dR_q2genjet = dR_q2;
      hasq2genjet = true;
    }
  }
  if (hasb1genjet and hasb2genjet){
    hastwogenjets = true;
    std::cout <<"dR_b1genjet "<< dR_b1genjet<< " dR_b2genjet "<< dR_b2genjet <<" genjet b1 pt  "<< b1genjet->pt() << " b2 pt "<< b2genjet->pt()<< std::endl;
    //std::cout <<" b1genjet: genparticle components num: "<< b1genjet->numberOfDaughters() << std::endl;
    //for (auto genp : b1genjet->getGenConstituents())
    //  std::cout <<"genp id "<< genp->pdgId()<<" px "<< genp->px()<<" py "<< genp->py()<<" pt "<< genp->pt() << std::endl;
  }
  if (hasq1genjet and hasq2genjet){
    hastwogenjets_wqq = true;
    std::cout <<"dR_q1genjet "<< dR_q1genjet<< " dR_q2genjet "<< dR_q2genjet <<" genjet q1 pt  "<< q1genjet->pt() << " q2 pt "<< q2genjet->pt()<< std::endl;
  }

}

//---------- method called to fill branches --------------------------------------------
void 
DiHiggsWWBBSLAnalyzer::fillbranches(){
  using namespace std;//include min, max
  if (findAllGenParticles){
    lep1_energy = lep1cand->energy();
    lep1_pt = lep1cand->pt();
    lep1_eta = lep1cand->eta();
    lep1_phi = lep1cand->phi();
    lep1_px = lep1cand->px();
    lep1_py = lep1cand->py();
    lep1_pz = lep1cand->pz();
    lep1_pdgid = lep1cand->pdgId();
    nu1_energy = nu1cand->energy();
    nu1_pt = nu1cand->pt();
    nu1_eta = nu1cand->eta();
    nu1_phi = nu1cand->phi();
    nu1_px = nu1cand->px();
    nu1_py = nu1cand->py();
    nu1_pz = nu1cand->pz();
    nu1_pdgid = nu1cand->pdgId();

    q1_energy = q1cand->energy();
    q1_pt     = q1cand->pt();
    q1_eta    = q1cand->eta();
    q1_phi    = q1cand->phi();
    q1_px     = q1cand->px();
    q1_py     = q1cand->py();
    q1_pz     = q1cand->pz();
    q1_pdgid  = q1cand->pdgId();
    q2_energy = q2cand->energy();
    q2_pt     = q2cand->pt();
    q2_eta    = q2cand->eta();
    q2_phi    = q2cand->phi();
    q2_px     = q2cand->px();
    q2_py     = q2cand->py();
    q2_pz     = q2cand->pz();
    q2_pdgid  = q2cand->pdgId();
    //lep2_energy = lep2cand->energy();
    //lep2_pt = lep2cand->pt();
    //lep2_eta = lep2cand->eta();
    //lep2_phi = lep2cand->phi();
    //lep2_px = lep2cand->px();
    //lep2_py = lep2cand->py();
    //lep2_pz = lep2cand->pz();
    //lep2_pdgid = lep2cand->pdgId();
    //nu2_energy = nu2cand->energy();
    //nu2_pt = nu2cand->pt();
    //nu2_px = nu2cand->px();
    //nu2_py = nu2cand->py();
    //nu2_pz = nu2cand->pz();
    //nu2_eta = nu2cand->eta();
    //nu2_phi = nu2cand->phi();
    //nu2_pdgid = nu2cand->pdgId();

    b1_energy = b1cand->energy();
    b1_pt = b1cand->pt();
    b1_eta = b1cand->eta();
    b1_phi = b1cand->phi();
    b1_px = b1cand->px();
    b1_py = b1cand->py();
    b1_pz = b1cand->pz();
    b2_energy = b2cand->energy();
    b2_pt = b2cand->pt();
    b2_eta = b2cand->eta();
    b2_phi = b2cand->phi();
    b2_px = b2cand->px();
    b2_py = b2cand->py();
    b2_pz = b2cand->pz();
    if (sampleType_ > Data and (htoWW or ttbar)){
      w1_energy = w1cand->energy();
      w1_pt = w1cand->pt();
      w1_eta = w1cand->eta();
      w1_phi = w1cand->phi();
      w1_px = w1cand->px();
      w1_py = w1cand->py();
      w1_pz = w1cand->pz();
      w1_mass = w1cand->mass();
      w2_energy = w2cand->energy();
      w2_pt = w2cand->pt();
      w2_eta = w2cand->eta();
      w2_phi = w2cand->phi();
      w2_px = w2cand->px();
      w2_py = w2cand->py();
      w2_pz = w2cand->pz();
      w2_mass = w2cand->mass();
    }
    if (sampleType_>=Rad_260 and htoWW){
      htoWW_energy = htoWWcand->energy();
      htoWW_px = htoWWcand->px();
      htoWW_py = htoWWcand->py();
      htoWW_pz = htoWWcand->pz();
      htoWW_mass = htoWWcand->mass();

      htobb_energy = htoBBcand->energy();
      htobb_px = htoBBcand->px();
      htobb_py = htoBBcand->py();
      htobb_pz = htoBBcand->pz();
      htobb_mass = htoBBcand->mass();

      Xtohh_energy = Xtohhcand->energy();
      Xtohh_px = Xtohhcand->px();
      Xtohh_py = Xtohhcand->py();
      Xtohh_pz = Xtohhcand->pz();
      Xtohh_mass = Xtohhcand->mass();
    }else if (sampleType_==TTbar){
      t1_energy = t1cand->energy();
      t1_px = t1cand->px();
      t1_py = t1cand->py();
      t1_pz = t1cand->pz();
      t2_energy = t2cand->energy();
      t2_px = t2cand->px();
      t2_py = t2cand->py();
      t2_pz = t2cand->pz();
    }else if (sampleType_ >= Rad_260_ZZbb and htoZZ){
      z1_energy = z1cand->energy();
      z1_pt = z1cand->pt();
      z1_eta = z1cand->eta();
      z1_phi = z1cand->phi();
      z1_px = z1cand->px();
      z1_py = z1cand->py();
      z1_pz = z1cand->pz();
      z1_mass = z1cand->mass();
      z2_energy = z2cand->energy();
      z2_pt = z2cand->pt();
      z2_eta = z2cand->eta();
      z2_phi = z2cand->phi();
      z2_px = z2cand->px();
      z2_py = z2cand->py();
      z2_pz = z2cand->pz();
      z2_mass = z2cand->mass();
        
      htobb_energy = htoBBcand->energy();
      htobb_px = htoBBcand->px();
      htobb_py = htoBBcand->py();
      htobb_pz = htoBBcand->pz();
      htobb_mass = htoBBcand->mass();

      Xtohh_energy = Xtohhcand->energy();
      Xtohh_px = Xtohhcand->px();
      Xtohh_py = Xtohhcand->py();
      Xtohh_pz = Xtohhcand->pz();
      Xtohh_mass = Xtohhcand->mass();
    }

  if (hastwogenjets){
    b1genjet_energy = b1genjet->energy();
    b1genjet_pt = b1genjet->pt();
    b1genjet_eta = b1genjet->eta();
    b1genjet_phi = b1genjet->phi();
    b1genjet_px = b1genjet->px();
    b1genjet_py = b1genjet->py();
    b1genjet_pz = b1genjet->pz();
    b2genjet_energy = b2genjet->energy();
    b2genjet_pt = b2genjet->pt();
    b2genjet_eta = b2genjet->eta();
    b2genjet_phi = b2genjet->phi();
    b2genjet_px = b2genjet->px();
    b2genjet_py = b2genjet->py();
    b2genjet_pz = b2genjet->pz();
  }
  if (hastwogenjets_wqq){
    q1genjet_energy = q1genjet->energy();
    q1genjet_pt     = q1genjet->pt();
    q1genjet_eta    = q1genjet->eta();
    q1genjet_phi    = q1genjet->phi();
    q1genjet_px     = q1genjet->px();
    q1genjet_py     = q1genjet->py();
    q1genjet_pz     = q1genjet->pz();
    q2genjet_energy = q2genjet->energy();
    q2genjet_pt     = q2genjet->pt();
    q2genjet_eta    = q2genjet->eta();
    q2genjet_phi    = q2genjet->phi();
    q2genjet_px     = q2genjet->px();
    q2genjet_py     = q2genjet->py();
    q2genjet_pz     = q2genjet->pz();
  }

    TLorentzVector lep1_p4(lep1cand->px(), lep1cand->py(), lep1cand->pz(), lep1cand->energy());
    //TLorentzVector lep2_p4(lep2cand->px(), lep2cand->py(), lep2cand->pz(), lep2cand->energy());
    //TLorentzVector diNeutrinos_p4(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py(), nu1cand->pz()+nu2cand->pz(), nu1cand->energy()+nu2cand->energy());
    //genmet_neutrinos_px = diNeutrinos_p4.Px();
    //genmet_neutrinos_py = diNeutrinos_p4.Py();
    //genmet_neutrinos_pt = diNeutrinos_p4.Pt();
    //genmet_neutrinos_phi = diNeutrinos_p4.Phi();
    //mass_gennu1nu2 = diNeutrinos_p4.M();

    TLorentzVector b1_p4(b1cand->px(), b1cand->py(), b1cand->pz(), b1cand->energy());
    TLorentzVector b2_p4(b2cand->px(), b2cand->py(), b2cand->pz(), b2cand->energy());
    TLorentzVector q1_p4(q1cand->px(), q1cand->py(), q1cand->pz(), q1cand->energy());
    TLorentzVector q2_p4(q2cand->px(), q2cand->py(), q2cand->pz(), q2cand->energy());
    //dR_genbl   = (b1_p4.Pt()>b2_p4.Pt()) ? (b1_p4.DeltaR( (lep1_p4.Pt()>lep2_p4.Pt()) ? lep1_p4 : lep2_p4 )) : (b2_p4.DeltaR( (lep1_p4.Pt()>lep2_p4.Pt()) ? lep1_p4 : lep2_p4 ));
    dR_genb1l1 = b1_p4.DeltaR(lep1_p4);
    //dR_genb1l2 = b1_p4.DeltaR(lep2_p4);
    dR_genb2l1 = b2_p4.DeltaR(lep1_p4);
    //dR_genb2l2 = b2_p4.DeltaR(lep2_p4);
    dR_genb1b2 = b1_p4.DeltaR(b2_p4);
    //dR_genl1l2 = lep1_p4.DeltaR(lep2_p4);
    //dR_genl1l2b1b2 = (lep1_p4+lep2_p4).DeltaR(b1_p4+b2_p4);
    //dphi_genl1l2b1b2 = TVector2::Phi_mpi_pi( (lep1_p4+lep2_p4).Phi()-(b1_p4+b2_p4).Phi() );
    //TLorentzVector genll_p4 = lep1_p4+lep2_p4;
    TLorentzVector genbb_p4 = b1_p4+b2_p4;
    TLorentzVector genqq_p4 = q1_p4+q2_p4;
    //dR_genminbl = min(min(dR_genb1l1,dR_genb1l2),min(dR_genb2l1,dR_genb2l2));
    //dphi_genllbb = TVector2::Phi_mpi_pi(genll_p4.Phi()-genbb_p4.Phi());
    //dphi_genllmet = TVector2::Phi_mpi_pi(genll_p4.Phi()-genmet_phi);
    //mass_genl1l2 = genll_p4.M(); energy_genl1l2 = genll_p4.Energy(); pt_genl1l2 = genll_p4.Pt(); eta_genl1l2 = genll_p4.Eta(); phi_genl1l2 = genll_p4.Phi();
    mass_genb1b2 = genbb_p4.M(); energy_genb1b2 = genbb_p4.Energy(); pt_genb1b2 = genbb_p4.Pt(); eta_genb1b2 = genbb_p4.Eta(); phi_genb1b2 = genbb_p4.Phi();
    mass_genq1q2 = genqq_p4.M(); energy_genq1q2 = genqq_p4.Energy(); pt_genq1q2 = genqq_p4.Pt(); eta_genq1q2 = genqq_p4.Eta(); phi_genq1q2 = genqq_p4.Phi();
    //mass_gentrans = sqrt(2*genll_p4.Pt()*genmet_pt*(1-cos(dphi_genllmet)));
  }
}


//-------------- method called to find a muon daughter for given candidate -------------------------
const reco::Candidate* 
DiHiggsWWBBSLAnalyzer::findLeptonDaughter(const reco::Candidate *Wcand){

  const reco::Candidate *tmpcand = NULL;
  int count = 0;
  for (unsigned int i = 0; i<Wcand->numberOfDaughters(); i++)
    if (abs(Wcand->daughter(i)->pdgId()) == 13 or abs(Wcand->daughter(i)->pdgId()) == 11) {
    count++;
    tmpcand = Wcand->daughter(i);
    }
  if (count != 1) std::cout << " this candidate has more one lep daughter " << std::endl;
  return tmpcand;
}


//------------ method called to find muon descendants for given candidate -------------------------
// return one muon descendant and the number of muon descendants by "count"
const reco::Candidate*
DiHiggsWWBBSLAnalyzer::findLeptonDescendants(const reco::Candidate *cand, int& count){

  const reco::Candidate * tmpcand = NULL;
  for (unsigned int i = 0; i<cand->numberOfDaughters(); i++){
    tmpcand = findLeptonDescendants(cand->daughter(i), count);
    if (abs(cand->daughter(i)->pdgId()) == 13 or abs(cand->daughter(i)->pdgId()) == 11) count++;
    if (tmpcand == NULL && cand->daughter(i)->pdgId() == cand->pdgId())  tmpcand = cand->daughter(i);                            
  }

  return tmpcand;
}


//-------------- method called to find a neutrino daughter for given candidate -------------------------
const reco::Candidate* 
DiHiggsWWBBSLAnalyzer::findNuDaughter(const reco::Candidate *Wcand){

  const reco::Candidate *tmpcand = NULL;
  int count = 0;
  for (unsigned int i = 0; i<Wcand->numberOfDaughters(); i++)
    if (abs(Wcand->daughter(i)->pdgId()) == 12 or abs(Wcand->daughter(i)->pdgId()) == 14) {
      count++;
      tmpcand = Wcand->daughter(i);
    }
  if (count != 1) std::cout << " this candidate has more one nu daughter " << std::endl;
  return tmpcand;
}

//------------ method called to find neutrino descendants for given candidate -------------------------
// return one neutrino descendant and the number of neutrno descendants by "count"
const reco::Candidate*
DiHiggsWWBBSLAnalyzer::findNuDescendants(const reco::Candidate *cand, int& count){

  const reco::Candidate * tmpcand = NULL;
  for (unsigned int i = 0; i<cand->numberOfDaughters(); i++){
    tmpcand = findLeptonDescendants(cand->daughter(i), count);
    if (abs(cand->daughter(i)->pdgId()) == 12 or abs(cand->daughter(i)->pdgId()) == 14) count++;
    if (tmpcand == NULL && cand->daughter(i)->pdgId() == cand->pdgId())  tmpcand = cand->daughter(i);                            
  }
  return tmpcand;
}


//------------- method called to find stable decendant with pdgid = id
const reco::Candidate* 
DiHiggsWWBBSLAnalyzer::stableDecendant(const reco::Candidate* cand, int id){
  const reco::Candidate* tmp = NULL;
  for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
    if ((cand->daughter(i))->pdgId() == id && (cand->daughter(i))->status() == 1)
	return 	tmp=cand->daughter(i);
    else if (stableDecendant(cand->daughter(i),id)) 
	return tmp=stableDecendant(cand->daughter(i),id);
  }   
  return tmp;
}


//------------ method called to find decendant with pdgid = id, 
//if first it true, then return the candidate closest to seed
//if first it false, then return the candidate farthest to seed
const reco::Candidate* 
DiHiggsWWBBSLAnalyzer::findDecendant(const reco::Candidate* cand, int id, bool first){
  /*const reco::Candidate* tmp1 = cand;
  while (tmp1->numberOfDaughters() == 1 and tmp1->daughter(0)->pdgId() == cand->pdgId()) tmp1 = tmp1->daughter(0);
  std::cout <<"In Findecentdant daughters of cand, id "<< cand->pdgId() << std::endl;
  for (unsigned int i=0; i < tmp1->numberOfDaughters(); i++){
      printCandidate(tmp1->daughter(i));
  }*/
 
  const reco::Candidate* tmp = NULL;
  for (unsigned int i=0; i < cand->numberOfDaughters(); i++){

    if ((cand->daughter(i))->pdgId() == id && first && cand->pdgId() != id)
	return 	tmp=cand->daughter(i);
    else if ((cand->daughter(i))->pdgId() == id && !first && !hasDaughter(cand->daughter(i), id)) 
	return  tmp=cand->daughter(i); // tmp does not has daughter with pdgid = id
    else if ((cand->daughter(i))->pdgId() == id && !first && (cand->daughter(i))->numberOfDaughters()>1) 
	return  tmp=cand->daughter(i);// tmp has more one daughters therefore it is final-states
    else if (findDecendant(cand->daughter(i),id, first)) 
	return tmp=findDecendant(cand->daughter(i),id);
  }
  return tmp;
}

//---------- method called to find a ancestor with pdgid = id, 
//if first is true, then return the candidate closest to seed
//if first is false, then return the candidate furthest to seed
const reco::Candidate*
DiHiggsWWBBSLAnalyzer::findAncestor(const reco::Candidate* cand, int id, bool first){

  const reco::Candidate* tmp = NULL;
  for (unsigned int i=0; i < cand->numberOfMothers(); i++){

    if ((cand->mother(i))->pdgId() == id && first && cand->pdgId() != id)
	return 	tmp=cand->mother(i);
    else if ((cand->mother(i))->pdgId() == id && !first && !hasMother(cand->mother(i), id)) 
	return  tmp=cand->mother(i);
    else if (findAncestor(cand->mother(i),id, first)) 
	return tmp=findAncestor(cand->mother(i),id, first);
  }
  return tmp;

}

//---------- method called to check whether cand has mother with pdgid = id -----------------------------------
bool
DiHiggsWWBBSLAnalyzer::hasMother(const reco::Candidate* cand, int id){

  for (unsigned int i=0; i < cand->numberOfMothers(); i++)
    if ((cand->mother(i))->pdgId() == id) return true;
  return false;

}

//-------- method called to check whether cand has daughter with pdgid = id ------------------------------------
bool
DiHiggsWWBBSLAnalyzer::hasDaughter(const reco::Candidate* cand, int id){

  //for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
  //	   std::cout <<"candidate id "<< cand->pdgId()<<" daughter i "<< i <<" id "<<(cand->daughter(i))->pdgId()<< std::endl;
  if (not(cand)) return false;
  //std::cout <<"cand id "<< cand->pdgId() <<" targetid "<< id << std::endl;
  for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
    if ((cand->daughter(i))->pdgId() == id) return true;
  return false;

}


//----------- method called to calculate total lorentz vector from b bar jets ---------------
//
TLorentzVector
DiHiggsWWBBSLAnalyzer::stableDecendantsLorentz(const reco::Candidate *cand){

  TLorentzVector bjets = TLorentzVector();
  for (unsigned i = 0; i < cand->numberOfDaughters(); i++){
    if ((cand->daughter(i))->status() == 1){
	TLorentzVector tmp((cand->daughter(i))->px(), (cand->daughter(i))->py(), 
	    (cand->daughter(i))->pz(),(cand->daughter(i))->energy());
	//printCandidate(cand->daughter(i));
	bjets = bjets+tmp; 
    }else bjets = bjets+stableDecendantsLorentz(cand->daughter(i));
  }

  return bjets;
}

//-------------  method called to calculate MET in simuation ------------------ 
TLorentzVector 
DiHiggsWWBBSLAnalyzer::calculateMET(){

  TLorentzVector METlorentz = TLorentzVector();
  TVector2 met_pxpy(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py());
  METlorentz.SetPxPyPzE(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py(),0,met_pxpy.Mod());

  return METlorentz;
}




// ------------ method called for printing additional information, useful for debugging  ------------
void
DiHiggsWWBBSLAnalyzer::print() {
  //todo: combine all necessary print out for debug 
  std::cout << "print() " << std::endl;
}

//------------- method called for printing h->WW->mumununu chain -------------------------
/*
   void 
   DiHiggsWWBBSLAnalyzer::printHtoWWChain(){

   float h_energy = lep1_energy+lep2_energy+nu1_energy+nu2_energy;
   float h_px = lep1_px+lep2_px+nu1_px+nu2_px;
   float h_py = lep1_py+lep2_py+nu1_py+nu2_py;
   float h_pz = lep1_pz+lep2_pz+nu1_pz+nu2_pz;

   TLorentzVector final_p4(h_px, h_py, h_pz, h_energy);

   float h_energy_bb = b1_energy+b2_energy; 
   float h_px_bb = b1_px+b2_px; 
   float h_py_bb = b1_py+b2_py;
   float h_pz_bb = b1_pz+b2_pz;
   TLorentzVector h_bb_p4(h_px_bb,h_py_bb,h_pz_bb,h_energy_bb);

   TLorentzVector X_final_tot(final_p4+h_bb_p4); 
   if (htoWW){

   TLorentzVector h_p4(htoWW_px, htoWW_py, htoWW_pz, htoWW_energy);
   TLorentzVector WW_p4(lep1_W1_cand->px()+lep2_W2_cand->px(),lep1_W1_cand->py()+lep2_W2_cand->py(), 
   lep1_W1_cand->pz()+lep2_W2_cand->pz(),lep1_W1_cand->energy()+lep2_W2_cand->energy());
   std::cout << "invariant mass from h_p4: " << h_p4.M() 
   << "	, from WW_p4 " << WW_p4.M() 
   << "	, from final_p4 " << final_p4.M() << std::endl;
   if (abs(WW_p4.M()-h_p4.M())>1)  std::cout << "h->WW, invariant mass reconstruction discrepancy ? " << std::endl; 
   std::cout <<  " H -> WW " << std::endl;
   const reco::Candidate* tmp_cand1 = NULL;
   const reco::Candidate* tmp_cand2 = NULL;

   for (unsigned int n = 0; n<lep1_htoWW_cand->numberOfDaughters(); n++){                                  
   if ((lep1_htoWW_cand->daughter(n))->pdgId()==-24) tmp_cand1 =  lep1_htoWW_cand->daughter(n);
   if ((lep1_htoWW_cand->daughter(n))->pdgId()== 24) tmp_cand2 =  lep1_htoWW_cand->daughter(n);
   if (n >= 2) std::cout << "h has more 2 daughters, id " << (lep1_htoWW_cand->daughter(n))->pdgId() << std::endl;
   }
   while (tmp_cand1 != lep1_W1_cand || tmp_cand2 != lep2_W2_cand){
   TLorentzVector W1_p4(tmp_cand1->px(),tmp_cand1->py(),tmp_cand1->pz(), tmp_cand1->energy());
   TLorentzVector W2_p4(tmp_cand2->px(),tmp_cand2->py(),tmp_cand2->pz(), tmp_cand2->energy());
   TLorentzVector tmp_WW_p4(W1_p4+W2_p4);
   std::cout <<"W- num of daughters "<< tmp_cand1->numberOfDaughters() << " W-mass " << W1_p4.M() << " W1 four momentum "; W1_p4.Print();
   std::cout <<"W+ num of daughters "<< tmp_cand2->numberOfDaughters() << " W+mass " << W2_p4.M() << " W2 four momentum "; W2_p4.Print();
   std::cout << "Total invariant mass " << tmp_WW_p4.M() <<" tmp_WW four momentum "; tmp_WW_p4.Print(); 
//if (tmp_cand1 != lep1_W1_cand) {
for (unsigned int i = 0; i<tmp_cand1->numberOfDaughters(); i++){
std::cout << " daughter of W- , id " << tmp_cand1->daughter(i)->pdgId() << "  status " << tmp_cand1->daughter(i)->status() <<std::endl; 
TLorentzVector dau_W1_p4(tmp_cand1->daughter(i)->px(),tmp_cand1->daughter(i)->py(),tmp_cand1->daughter(i)->pz(),tmp_cand1->daughter(i)->energy());
std::cout << " four momentum "; dau_W1_p4.Print();
if (tmp_cand1->daughter(i)->pdgId() == -24) tmp_cand1 = tmp_cand1->daughter(i);  
}
// }
//if (tmp_cand2 != lep2_W2_cand) {
for (unsigned int j = 0; j<tmp_cand2->numberOfDaughters(); j++){
std::cout << " daughter of W+ , id " << tmp_cand2->daughter(j)->pdgId() << "  status " << tmp_cand2->daughter(j)->status() <<std::endl; 
TLorentzVector dau_W2_p4(tmp_cand2->daughter(j)->px(),tmp_cand2->daughter(j)->py(),tmp_cand2->daughter(j)->pz(),tmp_cand2->daughter(j)->energy());
std::cout << " four momentum "; dau_W2_p4.Print();
if (tmp_cand2->daughter(j)->pdgId() == 24) tmp_cand2 = tmp_cand2->daughter(j);  
}
// }
}
std::cout << "WW -> mumununu" << std::endl;
while (tmp_cand1->status() != 1 || tmp_cand2->status() != 1){

std::cout << "-------------  begin of this loop ----------------------" << std::endl;

int size1 = tmp_cand1->numberOfDaughters();
int size2 = tmp_cand2->numberOfDaughters();
int muon1 = -1;
int muon2 = -1;
float px=0;
float py=0;
float pz=0;
float energy=0;
if (tmp_cand1->pdgId() == 13)  {
  //std::cout << "cand1 reaches final states, status of particle " << tmp_cand1->status() << std::endl;
  px += nu1_px;
  py += nu1_py;
  pz += nu1_pz;
  energy += nu1_energy;
}
std::cout << "cand1, id"<< tmp_cand1->pdgId() << " status " << tmp_cand1->status() <<" size of daughters " << size1 << std::endl; 
std::cout << " daughters of " << ((tmp_cand1->pdgId()==-24)?"W-  ":"muon- ") << std::endl;
for (int i = 0; i < size1; i++){
  const Candidate * d1 = tmp_cand1->daughter(i); 
  std::cout << "daughter id " << d1->pdgId() << "  status " << d1->status() << std::endl;
  if (d1->pdgId() == 13 ) muon1 = i;
  printCandidate(d1);
  px += d1->px();
  py += d1->py();
  pz += d1->pz();
  energy += d1->energy(); 
}
TLorentzVector cand1_lorentz(px, py, pz, energy); 
std::cout << " W- mass from W- Candidate " << w1_mass << " from mu-,nu " << cand1_lorentz.M() << std::endl;
if (muon1 != -1 && tmp_cand1->status() != 1) tmp_cand1 = tmp_cand1->daughter(muon1);
float px2 = 0.0;
float py2 = 0.0;
float pz2 = 0.0;
float energy2 = 0.0;
if (tmp_cand2->pdgId() == -13)  {
  //  std::cout << "cand2 reaches final states, status of particle "<< tmp_cand2->status() << std::endl;
  px2 += nu2_px;
  py2 += nu2_py;
  pz2 += nu2_pz;
  energy2 += nu2_energy;
}
std::cout << "cand2, id" << tmp_cand2->pdgId() <<" status " << tmp_cand2->status() <<" size of daughters " << size2 << std::endl; 
std::cout << " daughters of " << ((tmp_cand2->pdgId()==24)?"W+  ":"muon+ ") << std::endl;
for (int j = 0; j < size2; j++){
  const Candidate * d2 = tmp_cand2->daughter(j); 
  std::cout << "daughter id " << d2->pdgId() << "  status " << d2->status() << std::endl;
  if (d2->pdgId() == -13 ) muon2 = j;
  printCandidate(d2);
  px2 += d2->px();
  py2 += d2->py();
  pz2 += d2->pz();
  energy2 += d2->energy(); 
}
TLorentzVector cand2_lorentz(px2, py2, pz2, energy2); 
std::cout << " W+ mass from W+ Candidate " << w2_mass << " from mu+,nu " << cand2_lorentz.M() << std::endl;
if (muon2 != -1 && tmp_cand2->status() != 1) tmp_cand2 = tmp_cand2->daughter(muon2);
TLorentzVector tmp = cand1_lorentz+cand2_lorentz;
std::cout <<"Total px " << tmp.Px() << " py " << tmp.Py() << " pz " << tmp.Pz() << " E " << tmp.Energy() << std::endl;
std::cout << " invariant mass from daughters of WW " << tmp.M() << std::endl;  
std::cout << "For Next loop status of cand1 " << tmp_cand1->status()  << "	cand2 " << tmp_cand2->status() << std::endl; 
std::cout << "-------------  end of this loop ----------------------" << std::endl;
} 
if (findAllGenParticles){	
  std::cout <<"htoWW invariant mass " << final_p4.M() <<" four momentum " << std::endl; 
  final_p4.Print();
  std::cout <<"htobb invariant mass " << h_bb_p4.M() << " four momentum " << std::endl;
  h_bb_p4.Print();
  std::cout <<"Xtohh invariant mass " << X_final_tot.M() <<" total momentum " << std::endl;
  X_final_tot.Print();
} 

}//end if (htoWW)

}*/



//---------- method called to print candidates for debug ---------------------
void
DiHiggsWWBBSLAnalyzer::printCandidate(const reco::Candidate* cand){

  std::cout <<" Candidate id: "<< cand->pdgId() << " mass: " << cand->mass() <<" (P,E)= ("<< cand->px() <<", "<< cand->py()<<", "<< cand->pz()<<", "<< cand->energy() <<")" <<"(Pt,E) = ("<< cand->pt() <<", "<< cand->eta() <<", "<< cand->phi()<<", "<<cand->energy()<< ")" 
  <<" status: " << cand->status() << std::endl; 
  //<< (cand->isLastCopy() ? "lastCopy" : "NOTLastCopy") << (cand->isLastCopyBeforeFSR() ? " LastCopyBeforeFSR" : " NOTLastCopybeforeFSR")<< std::endl;

}

//--------- method called to print all decendants for cand -------------------
void 
DiHiggsWWBBSLAnalyzer::printallDecendants(const reco::Candidate* cand){

  if (cand->status() != 0 && cand->numberOfDaughters() > 0){
    std::cout << "******************  children of id "<< cand->pdgId() <<"      *********************" << std::endl;
    for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
	printCandidate(cand->daughter(i));
    std::cout << "***********************************************************" << std::endl;
    for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
	printallDecendants(cand->daughter(i));

  }
}


//--------- method called to print all Ancestors for cand -------------------
void 
DiHiggsWWBBSLAnalyzer::printallAncestors(const reco::Candidate* cand){

  if (cand->status() != 0 && cand->numberOfMothers() > 0){
    std::cout << "******************  mothers of id "<< cand->pdgId() <<"      *********************" << std::endl;
    for (unsigned int i=0; i < cand->numberOfMothers(); i++)
	printCandidate(cand->mother(i));
    std::cout << "***********************************************************" << std::endl;
    for (unsigned int i=0; i < cand->numberOfMothers(); i++)
	printallAncestors(cand->mother(i));

  }
}




// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiHiggsWWBBSLAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiHiggsWWBBSLAnalyzer);


