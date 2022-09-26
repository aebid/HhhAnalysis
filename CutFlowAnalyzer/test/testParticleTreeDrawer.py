import FWCore.ParameterSet.Config as cms

process = cms.Process("testParticle")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
ev1 = 86617
ev2 = 86617
process.source = cms.Source("PoolSource",
    eventsToProcess = cms.untracked.VEventRange("1:%d-1:%d"%(ev1, ev2)),
    #fileNames = cms.untracked.vstring('/store/relval/2008/5/20/RelVal-RelValTTbar-1211209682-FakeConditions-2nd/0000/08765709-5826-DD11-9CE8-000423D94700.root')
    fileNames = cms.untracked.vstring('file:/eos/user/t/tahuang/GluGluToRadionToHHTo2B2WToLNu2J_M-250_Run2016/6C2B2F7D-DF3C-E911-B863-B083FED00117.root')
)

process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
    #src = cms.InputTag("genParticles"),
    src = cms.InputTag("prunedGenParticles"),
    maxEventsToPrint  = cms.untracked.int32(1)
)

process.printTree2 = cms.EDAnalyzer("ParticleTreeDrawer",
    #src = cms.InputTag("genParticles"),
    src = cms.InputTag("prunedGenParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(False),
    printIndex  = cms.untracked.bool(False)
)

process.printEventNumber = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.printTree1*process.printTree2)
process.outpath = cms.EndPath(process.printEventNumber)


