from FastSimulation.Configuration.mixNoPU_cfi import *
from FastSimulation.Configuration.mixWithPU_cfi import *
# Gen Particles from Mixing Module
genParticlesFromMixingModule = cms.EDProducer("GenParticleProducer",
                                              saveBarCodes = cms.untracked.bool(True),
                                              useCrossingFrame = cms.untracked.bool(True),
                                              mix = cms.string("mixGenPU"),                          
                                              abortOnUnknownPDGCode = cms.untracked.bool(False)
                                              )
# Famos SimHits producer
from FastSimulation.EventProducer.FamosSimHits_cff import *
famosSimHits.GenParticleLabel = "genParticlesFromMixingModule" 
# PileupSummaryInfo
from SimGeneral.PileupInformation.AddPileupSummary_cfi import *
addPileupInfo.PileupMixingLabel = 'mixGenPU'
addPileupInfo.simHitLabel = 'famosSimHits'

famosMixing = cms.Sequence(
    mixGenPU+
    genParticlesFromMixingModule+
    addPileupInfo
    )
    
