# NanoAOD-format ntuple of CSC LCTs, reco muons, trigger information and (MC) simulation truth (GemCscNano_cff), from
# Run-3 inputs that carry both RAW and reco: the ZMu skim (/Muon0/Run2026D-ZMu-PromptReco-v1/RAW-RECO), AOD with its RAW
# parent given as secondary files (not yet run), or MC GEN-SIM-RECO with its GEN-SIM-(DIGI-)RAW parent as secondary
# files.
#
#   cmsRun runGemCscNano_cfg.py inputFiles=<file> globalTag=<GT of the input> era=<Run-3 era, e.g. Run3_2026> [mc=True] \
#       [maxEvents=-1] [secondaryInputFiles=<RAW parent>] [outputFile=...] [muonCut=<reco::Muon cut; empty = all>]
#
# Events written: a stored muon (muonCut) with pT > 10 GeV in 0.9 < |eta| < 2.4, or in MC a final-state generator muon
# with pT > 10 GeV in that range (gemCscNanoMuonFilter, gemCscNanoGenFilter of the cff). LCTs: unpacked from
# rawDataCollector here (no data tier stores them; MC RAW is packed from the emulator's MPC-sorted LCTs), every digi.
# Only Run-3 eras are accepted.
#
# Emulator LCTs from data: append to a copy of L1Trigger/CSCTriggerPrimitives/test/runCSCTriggerPrimitiveProducer_cfg.py
# (it sets process.schedule, hard-codes the data GT 140X_dataRun3_v3 and has no NanoAOD output):
#   from Configuration.AlCa.GlobalTag import GlobalTag
#   process.GlobalTag = GlobalTag(process.GlobalTag, "<GT of the runs>", "")
#   process.load("L1Trigger.GemCscNano.GemCscNano_cff")
#   process.cscLCTRecoMuonTable.lcts = "cscTriggerPrimitiveDigis"  # TMB level (L1TdeCSCTPG uses :MPCSORTED)
#   from L1Trigger.GemCscNano.GemCscNano_cff import setGemCscNanoSelection
#   setGemCscNanoSelection(process)  # after any change of muonCut or gemCscNanoSelection
#   process.gemCscNanoMuonPath = cms.Path(process.gemCscNanoMuonFilter)  # the event preselection
#   process.out = <the NanoAODOutputModule of this file; its l1bits keep statement then names the process L1CSCTPG>
#   process.e = cms.EndPath(process.out, process.gemCscNanoTask)
#   process.schedule.extend([process.gemCscNanoMuonPath, process.e])
# (not run since the event preselection was added, 2026-11-02)
# (the base cfg's era is Run3; the HLT objects and bits are read from the process "HLT", the cff default);
# run it with unpack=True l1=True run3=True runCCLUTOTMB=True runME11ILT=True useGEMAlignment=True
# saveEdmOutput=False [use6BitGEMCSCBendingAngle=True] inputFiles=<file>: without run3=True the cfg emulates Run-2
# TMBs (era Run2_2018), and with the ME1/1 ILT on, empty alignment LUTs throw. The emulator configuration (one
# bend LUT width for all ME1/1 chambers while the firmware is 4- or 6-bit per chamber, the all-zero placeholder
# alignment, clctNplanesHitPattern) is not matched to the runs' firmware configuration: carry that as a caveat.
import importlib

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from L1Trigger.GemCscNano.cscLCTRecoMuonTableProducer_cfi import cscLCTRecoMuonTableProducer as _defaults

options = VarParsing("analysis")
options.register("globalTag", "", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "global tag of the input (required)")
options.register("era", "", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Run-3 era of the input, a module name in Configuration.Eras, e.g. Run3_2026 (required)")
options.register("mc", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "simulation: add GenPart, SimTrack and SimTrackCSCLayer, and the generator-muon event selection")
options.register("muonCut", _defaults.muonCut.value(), VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "reco::Muon selection: only these muons are stored (Muon and the per-muon tables)")
options.setType("outputFile", VarParsing.varType.string)  # no automatic _numEvent<N> suffix
options.setDefault("outputFile", "GemCscNano.root")
options.parseArguments()

if not options.globalTag:
    raise RuntimeError("globalTag= is required: use the global tag of the input")
if not options.era.startswith("Run3") or "FastSim" in options.era:
    raise RuntimeError("era=%s: GemCscNano is defined for Run-3 inputs (full simulation or data) only; era= is required"
                       % options.era)

era = getattr(importlib.import_module("Configuration.Eras.Era_%s_cff" % options.era), options.era)
process = cms.Process("GEMCSCNANO", era)

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

process.load("L1Trigger.GemCscNano.GemCscNano_cff")
process.load("EventFilter.CSCRawToDigi.cscUnpacker_cfi")
process.cscLCTRecoMuonTable.muonCut = options.muonCut
from L1Trigger.GemCscNano.GemCscNano_cff import setGemCscNanoSelection
setGemCscNanoSelection(process)

tasks = [process.gemCscNanoTask, cms.Task(process.muonCSCDigis)]
selectEvents = ["gemCscNanoMuonPath"]
process.gemCscNanoMuonPath = cms.Path(process.gemCscNanoMuonFilter)
if options.mc:
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")  # GenParticlePruner
    tasks.append(process.gemCscNanoMCTask)
    process.gemCscNanoGenPath = cms.Path(process.gemCscNanoGenFilter)
    selectEvents.append("gemCscNanoGenPath")

process.out = cms.OutputModule("NanoAODOutputModule",
    fileName=cms.untracked.string(options.outputFile),
    compressionAlgorithm=cms.untracked.string("LZMA"),
    compressionLevel=cms.untracked.int32(9),
    SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring(*selectEvents)),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep nanoaodFlatTable_*Table_*_*",
        "keep edmTriggerResults_TriggerResults__HLT",  # HLT_* bits
        "keep edmTriggerResults_l1bits__%s" % process.name_(),  # L1_* bits (not those an input RECO process kept)
    ))
process.e = cms.EndPath(process.out, cms.Task(*tasks))
