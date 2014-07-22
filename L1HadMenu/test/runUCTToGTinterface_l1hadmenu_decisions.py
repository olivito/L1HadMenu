import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("ReRunningL1")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,skipEvents = cms.untracked.uint32(0))

readFiles.extend( [
    "/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/02B79593-F47F-E311-8FF6-003048FFD796.root",
] )

# Make the framework shut up.
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS162_V2::All'

# Load emulation and RECO sequences
isMC = True
if not isMC:
    process.load("L1Trigger.UCT2015.emulation_cfi")
    print "Running on data!"     
else:
    process.load("L1Trigger.UCT2015.emulationMC_cfi")

# Load sequences
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("L1Trigger.UCT2015.uctl1extraparticles_cfi")

## jet seed: 10 by default
###process.UCT2015Producer.jetSeed = cms.uint32(10) 

## configuration for sums: 7 by default
process.UCT2015Producer.regionETCutForHT = cms.uint32(7)
## eta range:
# 0-21 corresponds to all eta
# 4-17 corresponds to |eta| < 3.0 (default)
# 5-16 corresponds to |eta| < 2.2
process.UCT2015Producer.minGctEtaForSums = cms.uint32(4)
process.UCT2015Producer.maxGctEtaForSums = cms.uint32(17)

process.l1hadmenu = cms.EDProducer( 'L1HadMenuDecisionProducer' ,
                              L1JetsCentInputTag = cms.InputTag("l1extraParticlesUCT","Central","L1CustomNtupleProc"),
                              L1JetsFwdInputTag = cms.InputTag("l1extraParticlesUCT","Forward","L1CustomNtupleProc"),
                              L1EtMissInputTag = cms.InputTag("l1extraParticlesUCT","MET","L1CustomNtupleProc"),
                              L1MHTInputTag = cms.InputTag("l1extraParticlesUCT","MHT","L1CustomNtupleProc"),
                              )

process.p1 = cms.Path(
    process.emulationSequence *
    process.uct2015L1Extra *
    process.l1hadmenu
       #  *process.YourFavoritePlottingRoutine  --> This ends at l1hadmenu decisions production, anything after is up to the analyst 
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("output.root"),
    outputCommands = cms.untracked.vstring('drop *',
          'keep *_l1hadmenu_*_ReRunningL1',
    ) 
)

process.outPath = cms.EndPath( process.output )


