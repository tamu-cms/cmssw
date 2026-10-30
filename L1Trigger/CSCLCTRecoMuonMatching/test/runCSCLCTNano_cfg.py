# NanoAOD-format ntuple of CSC LCTs, reco muons, trigger information and (MC) simulation truth (cscLCTNano_cff), from
# inputs that carry both RAW and reco: the ZMu skim (/Muon0/Run2026D-ZMu-PromptReco-v1/RAW-RECO), AOD with its RAW
# parent given as secondary files, or MC GEN-SIM-RECO with its GEN-SIM-DIGI-RAW parent as secondary files.
#
#   cmsRun runCSCLCTNano_cfg.py inputFiles=<file> globalTag=<GT of the input> [era=Run3_2026] [geometry=<name>] \
#       [mc=True] [lcts=<module[:instance]>] [maxEvents=100] [secondaryInputFiles=<RAW parent>] [outputFile=...] \
#       [segmentSource=auto|standalone|matches] [maxDeltaHalfStrip=N] [maxDeltaWireGroup=N] [lctBXWindow=N] \
#       [hltProcess=HLT]
#
# Every event of the input is written (tag and probe is done on the ntuple). LCTs: data default = unpacked from
# rawDataCollector here (no data tier stores them); MC default = the emulator LCTs simCscTriggerPrimitiveDigis of the
# DIGI step (TMB level). geometry: empty = from the global tag (Run 3); for Phase 2 the name of a
# Configuration/Geometry Geometry<name>Reco_cff, e.g. ExtendedRun4D121. The match windows apply to the LCT-muon and
# (MC) the LCT-SimTrack links. lctBXWindow (default -1 = all LCTs): emulator LCTs span the readout window around BX 8,
# data LCTs carry a 1-bit BX and must keep -1.
#
# Emulator LCTs from data: append to a copy of L1Trigger/CSCTriggerPrimitives/test/runCSCTriggerPrimitiveProducer_cfg.py
# (it sets process.schedule, hard-codes the data GT 140X_dataRun3_v3 and has no NanoAOD output):
#   from Configuration.AlCa.GlobalTag import GlobalTag
#   process.GlobalTag = GlobalTag(process.GlobalTag, "<GT of the runs>", "")
#   process.load("L1Trigger.CSCLCTRecoMuonMatching.cscLCTNano_cff")
#   process.cscLCTRecoMuonTable.lcts = "cscTriggerPrimitiveDigis"  # TMB level (L1TdeCSCTPG uses :MPCSORTED)
#   process.cscLCTRecoMuonTable.matcher.lctBXWindow = <N>         # emulator LCTs span BX 5-11 around 8
#   and the NanoAODOutputModule EndPath of this file, appended to process.schedule;
# run it with unpack=True l1=True run3=True runCCLUTOTMB=True runME11ILT=True useGEMAlignment=True
# saveEdmOutput=False [use6BitGEMCSCBendingAngle=True] inputFiles=<file>: without run3=True the cfg emulates Run-2
# TMBs (era Run2_2018), and with the ME1/1 ILT on, empty alignment LUTs throw. The emulator configuration (one
# bend LUT width for all ME1/1 chambers while the firmware is 4- or 6-bit per chamber, the all-zero placeholder
# alignment, clctNplanesHitPattern) is not matched to the runs' firmware configuration: carry that as a caveat.
import importlib

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from L1Trigger.CSCLCTRecoMuonMatching.cscLCTRecoMuonTableProducer_cfi import cscLCTRecoMuonTableProducer as _defaults

options = VarParsing("analysis")
options.register("globalTag", "", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "global tag of the input (required)")
options.register("era", "Run3_2026", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "era module name in Configuration.Eras, e.g. Run3_2026 or Phase2C17I13M9")
options.register("geometry", "", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "empty: geometry from the global tag; else Configuration.Geometry.Geometry<geometry>Reco_cff")
options.register("mc", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "simulation: add GenPart, SimTrack and LCTSimTrackMatch")
options.register("lcts", "", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "LCT collection; empty: data muonCSCDigis:MuonCSCCorrelatedLCTDigi (unpacked here), "
                 "MC simCscTriggerPrimitiveDigis")
options.register("hltProcess", "HLT", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "process name of the HLT products (TriggerResults, hltTriggerSummaryAOD)")
options.register("segmentSource", _defaults.matcher.segmentSource.value(), VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "standalone | matches | auto")
options.register("maxDeltaHalfStrip", _defaults.matcher.maxDeltaHalfStrip.value(),
                 VarParsing.multiplicity.singleton, VarParsing.varType.int, "LCT-muon match window in key half-strips")
options.register("maxDeltaWireGroup", _defaults.matcher.maxDeltaWireGroup.value(),
                 VarParsing.multiplicity.singleton, VarParsing.varType.int, "LCT-muon match window in key wire groups")
options.register("lctBXWindow", _defaults.matcher.lctBXWindow.value(), VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "keep LCTs with |BX - 8| <= lctBXWindow; < 0: all (data: keep -1)")
options.setType("outputFile", VarParsing.varType.string)  # no automatic _numEvent<N> suffix
options.setDefault("outputFile", "cscLCTNano.root")
options.setDefault("maxEvents", 100)
options.parseArguments()

if not options.globalTag:
    raise RuntimeError("globalTag= is required: use the global tag of the input")
lcts = options.lcts or ("simCscTriggerPrimitiveDigis" if options.mc else "muonCSCDigis:MuonCSCCorrelatedLCTDigi")

era = getattr(importlib.import_module("Configuration.Eras.Era_%s_cff" % options.era), options.era)
process = cms.Process("CSCLCTNANO", era)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = True
if options.geometry:
    process.load("Configuration.Geometry.Geometry%sReco_cff" % options.geometry)
else:
    process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, "")

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
    secondaryFileNames=cms.untracked.vstring(options.secondaryInputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.load("L1Trigger.CSCLCTRecoMuonMatching.cscLCTNano_cff")
process.cscLCTRecoMuonTable.lcts = lcts
process.cscLCTRecoMuonTable.matcher.segmentSource = options.segmentSource
process.cscLCTRecoMuonTable.matcher.maxDeltaHalfStrip = options.maxDeltaHalfStrip
process.cscLCTRecoMuonTable.matcher.maxDeltaWireGroup = options.maxDeltaWireGroup
process.cscLCTRecoMuonTable.matcher.lctBXWindow = options.lctBXWindow
process.cscLCTSimTrackTable.maxDeltaHalfStrip = options.maxDeltaHalfStrip
process.cscLCTSimTrackTable.maxDeltaWireGroup = options.maxDeltaWireGroup
process.patTrigger.processName = options.hltProcess
for _module in (process.slimmedPatTrigger, process.unpackedPatTrigger):
    _module.triggerResults = cms.InputTag("TriggerResults", "", options.hltProcess)

tasks = [process.cscLCTNanoTask]
if lcts.split(":")[0] == "muonCSCDigis":
    process.load("EventFilter.CSCRawToDigi.cscUnpacker_cfi")
    tasks.append(cms.Task(process.muonCSCDigis))
if options.mc:
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")  # GenParticlePruner
    tasks.append(process.cscLCTNanoMCTask)

process.out = cms.OutputModule("NanoAODOutputModule",
    fileName=cms.untracked.string(options.outputFile),
    compressionAlgorithm=cms.untracked.string("LZMA"),
    compressionLevel=cms.untracked.int32(9),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep nanoaodFlatTable_*Table_*_*",
        "keep edmTriggerResults_TriggerResults__%s" % options.hltProcess,  # HLT_* bits
        "keep edmTriggerResults_l1bits__*",  # L1_* bits
    ))
process.e = cms.EndPath(process.out, cms.Task(*tasks))
