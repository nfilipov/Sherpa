import FWCore.ParameterSet.Config as cms
import os

source = cms.Source("EmptySource")

generator = cms.EDFilter("SherpaGeneratorFilter",
  maxEventsToPrint = cms.int32(0),
  filterEfficiency = cms.untracked.double(1.0),
  crossSection = cms.untracked.double(-1),
  SherpaProcess = cms.string('TTG'),
  SherpackLocation = cms.string('./'),
  SherpackChecksum = cms.string('20df6cce9ce0d95d43132ecac4086286'),
  FetchSherpack = cms.bool(False),
  SherpaPath = cms.string('./'),
  SherpaPathPiece = cms.string('./'),
  SherpaResultDir = cms.string('Result'),
  SherpaDefaultWeight = cms.double(1.0),
  SherpaParameters = cms.PSet(parameterSets = cms.vstring(
                             "MPI_Cross_Sections",
                             "Run"),
                              MPI_Cross_Sections = cms.vstring(
				" MPIs in Sherpa, Model = Amisic:",
				" semihard xsec = 43.3293 mb,",
				" non-diffractive xsec = 17.0318 mb with nd factor = 0.3142."
                                                  ),
                              Run = cms.vstring(
				" (run){",
				" EVENTS = 100;",
				" EVENT_MODE = HepMC;",
				" HEPMC2_GENEVENT_OUTPUT = hepmc;",
				" # avoid comix re-init after runcard modification",
				" WRITE_MAPPING_FILE 3;",
				" ME_SIGNAL_GENERATOR Comix Amegic LOOPGEN;",
				" EVENT_GENERATION_MODE Weighted;",
				" LOOPGEN:=BlackHat;",
				"}(run)",
				" (beam){",
				" BEAM_1 = 2212; BEAM_ENERGY_1 = 6500.;",
				" BEAM_2 = 2212; BEAM_ENERGY_2 = 6500.;",
				"}(beam)",
				" (processes){",
				" Process 93 93 -> 15 -15 22;",
				" Order (*,3);",
				" CKKW sqr(20./E_CMS);",
				" Print_Graphs MyGraphs;",
				" #Integration_Error 0.01;",
				" ME_Generator Amegic;",
				" End process;",
				"}(processes)",
				" (selector){",
				" Mass 15 -15 50 E_CMS;",
				" PT 22 10. E_CMS;",
				" PseudoRapidity 22 -3.0 3.0;",
				" PseudoRapidity 15 -3.0 3.0;",
				" PseudoRapidity -15 -3.0 3.0;",
				" DeltaR 22 15 0.6 100000;",
				" DeltaR 22 -15 0.6 100000;",
				"}(selector)",
				" (shower){",
				" CSS_EW_MODE = 1",
				"}(shower)",
				" (integration){",
				" FINISH_OPTIMIZATION = Off",
				"}(integration)",
				" (isr){",
				" PDF_LIBRARY     = LHAPDFSherpa",
				" PDF_SET         = CT10",
				" PDF_SET_VERSION = 0",
				"}(isr)",
				" (mi){",
				" MI_HANDLER = Amisic  # None or Amisic",
				"}(mi)"
                                                  ),
                             )
)

ProductionFilterSequence = cms.Sequence(generator)

