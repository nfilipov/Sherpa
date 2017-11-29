import FWCore.ParameterSet.Config as cms

import FastSimulation.HighLevelTrigger.DummyModule_cfi
from FastSimulation.Tracking.GlobalPixelTracking_cff import *

#############################################
# Reconstruct tracks with pixel seeds
#############################################

# Take all pixel tracks for b tagging track reco (pTMin>1GeV, nHits>=8)
hltFastTrackMergerForB = cms.EDProducer("FastTrackMerger",
    SaveTracksOnly = cms.untracked.bool(True),
    TrackProducers = cms.VInputTag(cms.InputTag("globalPixelWithMaterialTracks"),
                                   cms.InputTag("globalPixelTrackCandidates")),
    ptMin = cms.untracked.double(1.0),
    minHits = cms.untracked.uint32(8)
)




###############################

hltBLifetimeRegionalCkfTrackCandidatesHbb = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
hltBLifetimeRegionalPixelSeedGeneratorHbbVBF = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
hltBLifetimeBTagIP3D1stTrkRegionalCkfTrackCandidatesJet20 = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
hltBLifetimeBTagIP3D1stTrkRegionalCkfTrackCandidatesJet20 = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
hltBLifetimeRegionalPixelSeedGeneratorbbPhiL1FastJet = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
hltBLifetimeBTagIP3D1stTrkRegionalPixelSeedGeneratorJet20HbbL1FastJet = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
hltBLifetimeDiBTagIP3D1stTrkRegionalPixelSeedGeneratorJet20HbbL1FastJet = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
hltBLifetimeFastRegionalPixelSeedGeneratorHbbVBF = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
hltBLifetimeRegionalPixelSeedGeneratorbbPhiL1FastJetFastPV = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
hltFastPixelBLifetimeRegionalPixelSeedGeneratorHbb = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()

hltBLifetimeRegionalCkfTrackCandidatesHbb = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesHbbVBF = cms.Sequence(globalPixelTracking)
hltBLifetimeBTagIP3D1stTrkRegionalCkfTrackCandidatesJet20Hbb = cms.Sequence(globalPixelTracking)
hltBLifetimeDiBTagIP3D1stTrkRegionalCkfTrackCandidatesJet20Hbb = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesbbPhi = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesbbPhiL1FastJet = cms.Sequence(globalPixelTracking)
hltBLifetimeBTagIP3D1stTrkRegionalCkfTrackCandidatesJet20HbbL1FastJet = cms.Sequence(globalPixelTracking)
hltBLifetimeDiBTagIP3D1stTrkRegionalCkfTrackCandidatesJet20HbbL1FastJet = cms.Sequence(globalPixelTracking)
hltBLifetimeFastRegionalCkfTrackCandidatesHbbVBF = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesbbPhiL1FastJetFastPV = cms.Sequence(globalPixelTracking)
hltFastPixelBLifetimeRegionalCkfTrackCandidatesHbb = cms.Sequence(globalPixelTracking)

hltBLifetimeRegionalCtfWithMaterialTracksHbb = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksHbbVBF = hltFastTrackMergerForB.clone()
hltBLifetimeBTagIP3D1stTrkRegionalCtfWithMaterialTracksJet20Hbb = hltFastTrackMergerForB.clone()
hltBLifetimeDiBTagIP3D1stTrkRegionalCtfWithMaterialTracksJet20Hbb = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksbbPhi = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJet = hltFastTrackMergerForB.clone()
hltBLifetimeBTagIP3D1stTrkRegionalCtfWithMaterialTracksJet20HbbL1FastJet = hltFastTrackMergerForB.clone()
hltBLifetimeDiBTagIP3D1stTrkRegionalCtfWithMaterialTracksJet20HbbL1FastJet = hltFastTrackMergerForB.clone()
hltBLifetimeFastRegionalCtfWithMaterialTracksHbbVBF = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV = hltFastTrackMergerForB.clone()
hltFastPixelBLifetimeRegionalCtfWithMaterialTracksHbb = hltFastTrackMergerForB.clone()

hltBLifetimeRegionalCkfTrackCandidates = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesSingleTop = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesEleJetSingleTop = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesIsoEleJetSingleTop = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesRA2b = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesRAzr = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesHbb = cms.Sequence(globalPixelTracking)
hltBLifetimeRegional3DCkfTrackCandidatesJet30Hbb = cms.Sequence(globalPixelTracking)
hltBLifetimeRegional3D1stTrkCkfTrackCandidatesJet20Hbb = cms.Sequence(globalPixelTracking)
hltBLifetimeRegional3DCkfTrackCandidatesJet30Hbb = cms.Sequence(globalPixelTracking)
hltBLifetimeBTagIP3D1stTrkRegionalCkfTrackCandidatesJet20Hbb = cms.Sequence(globalPixelTracking)
hltBLifetimeDiBTagIP3D1stTrkRegionalCkfTrackCandidatesJet20Hbb = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesbbPhi = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesGammaB = cms.Sequence(globalPixelTracking)

hltBLifetimeRegionalCtfWithMaterialTracks = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksSingleTop = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksEleJetSingleTop = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksIsoEleJetSingleTop = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksRA2b = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksRAzr = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksHbb = hltFastTrackMergerForB.clone()
hltBLifetimeRegional3DCtfWithMaterialTracksJet30Hbb = hltFastTrackMergerForB.clone()
hltBLifetimeRegional3D1stTrkCtfWithMaterialTracksJet20Hbb = hltFastTrackMergerForB.clone()
hltBLifetimeRegional3DCtfWithMaterialTracksJet30Hbb = hltFastTrackMergerForB.clone()
hltBLifetimeBTagIP3D1stTrkRegionalCtfWithMaterialTracksJet20Hbb = hltFastTrackMergerForB.clone()
hltBLifetimeDiBTagIP3D1stTrkRegionalCtfWithMaterialTracksJet20Hbb = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksbbPhi = hltFastTrackMergerForB.clone()
hltBLifetimeRegionalCtfWithMaterialTracksGammaB = hltFastTrackMergerForB.clone()


#############################################
# Reconstruct muons for MumuK
#############################################
import FWCore.ParameterSet.Config as cms

# Take all pixel-seeded tracks for b tagging track reco (pTMin>1GeV, nHits>=8) 
hltCtfWithMaterialTracksMumuk = cms.EDProducer("FastTrackMerger",
    SaveTracksOnly = cms.untracked.bool(True),
    TrackProducers = cms.VInputTag(cms.InputTag("globalPixelWithMaterialTracks"),
                                   cms.InputTag("globalPixelTrackCandidates")),
    ptMin = cms.untracked.double(3.0),
    minHits = cms.untracked.uint32(5)
)

# produce ChargedCandidates from tracks
hltMumukAllConeTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hltCtfWithMaterialTracksMumuk"),
    particleType = cms.string('mu-')
)

hltCkfTrackCandidatesMumuk = cms.Sequence(cms.SequencePlaceholder("HLTL3muonrecoSequence"))


#############################################
# Reconstruct muons for JPsiToMumu
#############################################

# Take all pixel-seeded tracks for b tagging track reco (pTMin>1GeV, nHits>=8) 
hltCtfWithMaterialTracksMumu = cms.EDProducer("FastTrackMerger",
    SaveTracksOnly = cms.untracked.bool(True),
    TrackProducers = cms.VInputTag(cms.InputTag("hltL3Muons")),
    ptMin = cms.untracked.double(3.0),
    minHits = cms.untracked.uint32(5)
)

# produce ChargedCandidates from tracks
hltMuTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hltCtfWithMaterialTracksMumu"),
    particleType = cms.string('mu-')
)

hltCkfTrackCandidatesMumu = cms.Sequence(cms.SequencePlaceholder("HLTL3muonrecoNocandSequence"))


