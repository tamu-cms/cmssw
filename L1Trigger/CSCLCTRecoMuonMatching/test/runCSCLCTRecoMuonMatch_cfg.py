# Test of CSCLCTRecoMuonMatcher on data that carry both RAW and reco, e.g. the ZMu skim
# (/Muon0/Run2026D-ZMu-PromptReco-v1/RAW-RECO), or AOD with its RAW parent given as secondary files.
#
#   cmsRun runCSCLCTRecoMuonMatch_cfg.py inputFiles=<file> globalTag=<GT of the runs> [era=Run3_2026] \
#       [maxEvents=100] [secondaryInputFiles=<RAW parent of an AOD file>] [outputFile=...] \
#       [segmentSource=auto|standalone|matches] [maxDeltaHalfStrip=N] [maxDeltaWireGroup=N]
#
# Data LCTs are unpacked from rawDataCollector here (no data tier stores them). Only the CSC geometry is needed.
# Output: tree cscLCTRecoMuonMatchAnalyzer/lctMuon in outputFile, one entry per LCT with its matched muons (vector
# branches). Default windows: +-6 half-strips, +-2 wire groups (deliberately loose, to be tightened).
#
# Emulator LCTs: append to a copy of L1Trigger/CSCTriggerPrimitives/test/runCSCTriggerPrimitiveProducer_cfg.py
# (it sets process.schedule, has no TFileService, and hard-codes the data GT 140X_dataRun3_v3):
#   from Configuration.AlCa.GlobalTag import GlobalTag
#   process.GlobalTag = GlobalTag(process.GlobalTag, "<GT of the runs>", "")
#   process.load("L1Trigger.CSCLCTRecoMuonMatching.cscLCTRecoMuonMatchAnalyzer_cfi")
#   process.cscLCTRecoMuonMatchAnalyzer.lcts = "cscTriggerPrimitiveDigis"  # TMB level (L1TdeCSCTPG uses :MPCSORTED)
#   process.cscLCTRecoMuonMatchAnalyzer.matcher.lctBXWindow = <N>         # emulator LCTs span BX 5-11 around 8
#   process.TFileService = cms.Service("TFileService", fileName=cms.string("lctRecoMuonMatch_emu.root"))
#   process.lctMatchPath = cms.Path(process.cscLCTRecoMuonMatchAnalyzer)
#   process.schedule.append(process.lctMatchPath)
# and run it with unpack=True l1=True run3=True runCCLUTOTMB=True runME11ILT=True useGEMAlignment=True
# saveEdmOutput=False [use6BitGEMCSCBendingAngle=True] inputFiles=<file>: without run3=True the cfg emulates Run-2
# TMBs (era Run2_2018), and with the ME1/1 ILT on, empty alignment LUTs throw. The emulator configuration (one
# bend LUT width for all ME1/1 chambers while the firmware is 4- or 6-bit per chamber, the all-zero placeholder
# alignment, clctNplanesHitPattern) is not matched to the runs' firmware configuration: carry that as a caveat.
import importlib

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from L1Trigger.CSCLCTRecoMuonMatching.cscLCTRecoMuonMatchAnalyzer_cfi import (
    cscLCTRecoMuonMatchAnalyzer as _defaults,
)

options = VarParsing("analysis")
options.register("globalTag", "", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "global tag of the runs analysed (required)")
options.register("era", "Run3_2026", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "era module name in Configuration.Eras, e.g. Run3_2025 or Run3_2026")
options.register("segmentSource", _defaults.matcher.segmentSource.value(), VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "standalone | matches | auto")
options.register("maxDeltaHalfStrip", _defaults.matcher.maxDeltaHalfStrip.value(),
                 VarParsing.multiplicity.singleton, VarParsing.varType.int, "match window in key half-strips")
options.register("maxDeltaWireGroup", _defaults.matcher.maxDeltaWireGroup.value(),
                 VarParsing.multiplicity.singleton, VarParsing.varType.int, "match window in key wire groups")
options.setType("outputFile", VarParsing.varType.string)  # no automatic _numEvent<N> suffix
options.setDefault("outputFile", "lctRecoMuonMatch.root")
options.setDefault("maxEvents", 100)
options.parseArguments()

if not options.globalTag:
    raise RuntimeError("globalTag= is required: use the global tag of the runs being analysed")

era = getattr(importlib.import_module("Configuration.Eras.Era_%s_cff" % options.era), options.era)
process = cms.Process("LCTMUON", era)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = True
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, "")

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
    secondaryFileNames=cms.untracked.vstring(options.secondaryInputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))

process.load("EventFilter.CSCRawToDigi.cscUnpacker_cfi")  # muonCSCDigis, incl. MuonCSCCorrelatedLCTDigi
process.load("L1Trigger.CSCLCTRecoMuonMatching.cscLCTRecoMuonMatchAnalyzer_cfi")
process.cscLCTRecoMuonMatchAnalyzer.matcher.segmentSource = options.segmentSource
process.cscLCTRecoMuonMatchAnalyzer.matcher.maxDeltaHalfStrip = options.maxDeltaHalfStrip
process.cscLCTRecoMuonMatchAnalyzer.matcher.maxDeltaWireGroup = options.maxDeltaWireGroup

process.p = cms.Path(process.muonCSCDigis + process.cscLCTRecoMuonMatchAnalyzer)
