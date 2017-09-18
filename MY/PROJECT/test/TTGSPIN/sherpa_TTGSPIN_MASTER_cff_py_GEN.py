# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: MY/PROJECT/python/sherpa_TTGSPIN_MASTER_cff.py -s GEN -n 100 --no_exec --conditions auto:mc --eventcontent RAWSIM
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('MY/PROJECT/python/sherpa_TTGSPIN_MASTER_cff.py nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('/eos/user/n/nfilipov/zg_sherpa/sherpa_TTGSPIN_MASTER_cff_py_GEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDFilter("SherpaGeneratorFilter",
    SherpaProcess = cms.string('TTGSPIN'),
    SherpaParameters = cms.PSet(
        parameterSets = cms.vstring('MPI_Cross_Sections', 
            'Run'),
        Run = cms.vstring(' (run){', 
            ' EVENTS = 100;', 
            ' EVENT_MODE = HepMC;', 
            ' HEPMC2_GENEVENT_OUTPUT = hepmc;', 
            ' # avoid comix re-init after runcard modification', 
            ' WRITE_MAPPING_FILE 3;', 
            ' #HARD_DECAYS = 1', 
            ' HARD_SPIN_CORRELATIONS = 1', 
            ' SOFT_SPIN_CORRELATIONS = 1', 
            '}(run)', 
            ' (beam){', 
            ' BEAM_1 = 2212; BEAM_ENERGY_1 = 6500.;', 
            ' BEAM_2 = 2212; BEAM_ENERGY_2 = 6500.;', 
            '}(beam)', 
            ' (processes){', 
            ' # Process 93 93 -> 15 -15 22 93{1};', 
            ' Process 93 93 -> 15 -15 22;', 
            ' Order (*,3);', 
            ' CKKW sqr(20./E_CMS);', 
            ' Print_Graphs MyGraphs;', 
            ' # Integration_Error 0.01;', 
            ' End process;', 
            '}(processes)', 
            ' (selector){', 
            ' Mass 15 -15 30 E_CMS;', 
            ' PT 22 10. E_CMS;', 
            ' # PT 90 50. E_CMS;', 
            ' # PT -90 50. E_CMS;', 
            ' PseudoRapidity 22 -2.6 2.6;', 
            ' #GP PseudoRapidity 90 -3.0 3.0;', 
            ' #GP PseudoRapidity -90 -3.0 3.0;', 
            ' DeltaR 22 15 0.05 100000', 
            ' DeltaR 22 -15 0.05 100000', 
            '}(selector)', 
            ' (shower){', 
            ' CSS_EW_MODE = 1', 
            '}(shower)', 
            ' (integration){', 
            ' FINISH_OPTIMIZATION = Off', 
            '}(integration)', 
            ' (isr){', 
            ' PDF_LIBRARY     = LHAPDFSherpa', 
            ' PDF_SET         = CT10', 
            ' PDF_SET_VERSION = 0', 
            '}(isr)', 
            ' (me){', 
            ' # ME_SIGNAL_GENERATOR = Internal Comix', 
            ' ME_SIGNAL_GENERATOR = Amegic', 
            ' EVENT_GENERATION_MODE = Unweighted;', 
            ' # ME_SIGNAL_GENERATOR Comix Amegic LOOPGEN;', 
            ' # LOOPGEN:=BlackHat;', 
            '}(me)', 
            ' (mi){', 
            ' MI_HANDLER = Amisic  # None or Amisic', 
            '}(mi)'),
        MPI_Cross_Sections = cms.vstring(' MPIs in Sherpa, Model = Amisic:', 
            ' semihard xsec = 43.3293 mb,', 
            ' non-diffractive xsec = 17.0318 mb with nd factor = 0.3142.')
    ),
    filterEfficiency = cms.untracked.double(1.0),
    FetchSherpack = cms.bool(False),
    SherpackChecksum = cms.string('e7104ab53d6f45ecaf017503a5403dde'),
    SherpaResultDir = cms.string('Result'),
    SherpaPath = cms.string('./'),
    crossSection = cms.untracked.double(-1),
    maxEventsToPrint = cms.int32(0),
    SherpaPathPiece = cms.string('./'),
    SherpackLocation = cms.string('./'),
    SherpaDefaultWeight = cms.double(1.0)
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

