# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: NanoAODproduction_2016_cfg -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --no_exec --conditions 102X_mcRun2_asymptotic_v8 --era Run2_2016,run2_nanoAOD_94X2016 --nThreads=8
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('NANO',eras.Run2_2016,eras.run2_nanoAOD_94X2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

high_precision = True
if high_precision:
  from  PhysicsTools.NanoAOD.common_cff import *
  process.genParticleTable.variables = cms.PSet(
      pt  = Var("pt", float, precision=23),
      phi = Var("phi", float, precision=23),
      eta  = Var("eta", float, precision=23),
      mass = Var("?!((abs(pdgId)>=1 && abs(pdgId)<=5) || (abs(pdgId)>=11 && abs(pdgId)<=16) || pdgId==21 || pdgId==111 || abs(pdgId)==211 || abs(pdgId)==421 || abs(pdgId)==411 || (pdgId==22 && mass<1))?mass:0", float,precision="?((abs(pdgId)==6 || abs(pdgId)>1000000) && statusFlags().isLastCopy())?20:8",doc="Mass stored for all particles with the exception of quarks (except top), leptons/neutrinos, photons with mass < 1 GeV, gluons, pi0(111), pi+(211), D0(421), and D+(411). For these particles, you can lookup the value from PDG."),
      pdgId  = Var("pdgId", int, doc="PDG id"),
      status  = Var("status", int, doc="Particle status. 1=stable"),
      genPartIdxMother = Var("?numberOfMothers>0?motherRef(0).key():-1", int, doc="index of the mother particle"),
      statusFlags = (Var(
          "statusFlags().isLastCopyBeforeFSR()                  * 16384 +"
          "statusFlags().isLastCopy()                           * 8192  +"
          "statusFlags().isFirstCopy()                          * 4096  +"
          "statusFlags().fromHardProcessBeforeFSR()             * 2048  +"
          "statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +"
          "statusFlags().isHardProcessTauDecayProduct()         * 512   +"
          "statusFlags().fromHardProcess()                      * 256   +"
          "statusFlags().isHardProcess()                        * 128   +"
          "statusFlags().isDirectHadronDecayProduct()           * 64    +"
          "statusFlags().isDirectPromptTauDecayProduct()        * 32    +"
          "statusFlags().isDirectTauDecayProduct()              * 16    +"
          "statusFlags().isPromptTauDecayProduct()              * 8     +"
          "statusFlags().isTauDecayProduct()                    * 4     +"
          "statusFlags().isDecayedLeptonHadron()                * 2     +"
          "statusFlags().isPrompt()                             * 1      ",
          int,
          doc=(
              "gen status flags stored bitwise, bits are: "
              "0 : isPrompt, "
              "1 : isDecayedLeptonHadron, "
              "2 : isTauDecayProduct, "
              "3 : isPromptTauDecayProduct, "
              "4 : isDirectTauDecayProduct, "
              "5 : isDirectPromptTauDecayProduct, "
              "6 : isDirectHadronDecayProduct, "
              "7 : isHardProcess, "
              "8 : fromHardProcess, "
              "9 : isHardProcessTauDecayProduct, "
              "10 : isDirectHardProcessTauDecayProduct, "
              "11 : fromHardProcessBeforeFSR, "
              "12 : isFirstCopy, "
              "13 : isLastCopy, "
              "14 : isLastCopyBeforeFSR, "
            )
          )
      ),
  )


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:NanoAODproduction_2016_cfg_PAT.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('NanoAODproduction_2016_cfg nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('NanoAODproduction_2016_cfg_NANO.root'),
    fakeNameForCrab = cms.untracked.bool(True),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mcRun2_asymptotic_v8', '')

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
