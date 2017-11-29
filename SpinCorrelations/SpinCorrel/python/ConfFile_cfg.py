import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
        'file:/eos/user/n/nfilipov/zg_sherpa/sherpa_EEG_MASTER_cff_py_GEN.root'
#        'file:/eos/user/n/nfilipov/zg_sherpa/sherpa_MMG_MASTER_cff_py_GEN.root'
#       'file:/eos/user/n/nfilipov/zg_sherpa/sherpa_TTG_MASTER_cff_py_GEN.root'
#        'file:/eos/user/n/nfilipov/zg_sherpa/sherpa_TTGSPIN_MASTER_cff_py_GEN.root'
    )
)

process.myProducerLabel = cms.EDProducer('SpinCorrel',
                                         genParticles = cms.InputTag("genParticles")
                                         )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('/eos/user/n/nfilipov/zg_sherpa/EEGout_genParticles.root'),
                               outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_genParticles_*_*'
        ),
)

  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
