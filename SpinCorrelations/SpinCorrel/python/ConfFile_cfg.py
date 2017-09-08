import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
        'file:/afs/cern.ch/work/n/nfilipov/public/forGabi/zg_sherpa/zg_sherpa/ttg/sherpa_TTG_MASTER_cff_py_GEN.root'
    )
)

process.myProducerLabel = cms.EDProducer('SpinCorrel',
                                         genParticles = cms.InputTag("genParticles")
                                         )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('myOutputFile.root'),
                               outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_genParticles_*_*'
        ),
)

  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
