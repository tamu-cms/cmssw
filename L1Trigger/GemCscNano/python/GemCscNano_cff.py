# NanoAOD-format ntuple of CSC LCTs, ALCTs and CLCTs, reco muons, trigger information and simulation truth, one
# entry per event (Run 3). Events: those with a stored muon (cscLCTRecoMuonTable.muonCut) with pT > 10 GeV and
# 0.9 < |eta| < 2.4, or (MC) a final-state generator muon with pT > 10 GeV and 0.9 < |eta| < 2.4 (gemCscNanoMuonFilter,
# gemCscNanoGenFilter; the user's definition of 2026-11-02). Per event: every LCT, ALCT and CLCT; the reco muons
# passing muonCut (default pT > 10 GeV and the medium or tight cut-based ID). No LCT matching and no tag-and-probe
# selection is applied.
#   LCT, ALCT, CLCT, Muon,              cscLCTRecoMuonTable (CSCLCTRecoMuonTableProducer)
#   MuonCSCChamber,
#   MuonCSCSegment, MuonCSCLayer,
#   MuonCSCRecHitADC, MuonGEMChamber,
#   MuonGEMHitMatch
#   TrigObj                             the HLT muon objects (central muon selection, PAT trigger chain of miniAOD run
#                                       here on hltTriggerSummaryAOD) with one filter bit per single-muon path
#   L1Mu                                the central table, every BX of the uGMT readout
#   PV                                  primaryVertexTable: the central PV count columns and main PV, from RECO
#                                       offlinePrimaryVertices
#   TMBHeader, TMBHeaderRaw             tmbHeaderTable: the (O)TMB DAQ header of each chamber readout, words 0-26
#                                       (CSCTMBStatusDigi keeps no more); needs muonCSCDigis.UnpackStatusDigis = True
#   EMTFMuon                            emtfMuonTable: the EMTF output tracks as the uGMT receives them
#                                       (gmtStage2Digis:EMTF), every BX
#   EMTFHit, EMTFTrack, EMTFLCT         emtfDaqTable on the EMTF unpacker emtfStage2Digis (gemCscNanoDataTask: data
#                                       only, MC RAW has no EMTF readout): the EMTF DAQ hits, tracks and LCTs
#   HLT_* and L1_* bits                 the muon HLT paths and L1 seeds (gemCscNanoHLTBits, gemCscNanoL1Bits),
#                                       written by NanoAODOutputModule
#   GenPart, SimTrack, SimTrackCSCLayer gemCscNanoMCTask (MC only; GenPart and Pileup as in the central NanoAOD)
#   Pileup
# Module labels end in "Table" where NanoAODOutputModule must keep the product (keep nanoaodFlatTable_*Table_*_*).
# The PAT trigger, gen-particle and NanoAOD modules keep their central labels (patTrigger, prunedGenParticles,
# finalGenParticles, triggerObjectTable, l1bits, ...): this cff cannot be loaded into a job that also runs the PAT or
# NANO steps.
import FWCore.ParameterSet.Config as cms

from L1Trigger.GemCscNano.cscLCTRecoMuonTableProducer_cfi import cscLCTRecoMuonTableProducer
from EventFilter.L1TRawToDigi.emtfStage2Digis_cfi import emtfStage2Digis
from L1Trigger.GemCscNano.cscLCTSimTrackTableProducer_cfi import cscLCTSimTrackTableProducer
from L1Trigger.GemCscNano.emtfDaqTableProducer_cfi import emtfDaqTableProducer
from L1Trigger.GemCscNano.emtfMuonTableProducer_cfi import emtfMuonTableProducer
from L1Trigger.GemCscNano.tmbHeaderTableProducer_cfi import tmbHeaderTableProducer
from L1Trigger.GemCscNano.primaryVertexTableProducer_cfi import primaryVertexTableProducer
from L1Trigger.GemCscNano.triggerResultsSelector_cfi import triggerResultsSelector
from PhysicsTools.NanoAOD.genparticles_cff import finalGenParticles, genIso, genParticleTable
from PhysicsTools.NanoAOD.globals_cff import puTable as _puTable
from PhysicsTools.NanoAOD.l1trig_cff import l1MuTable as _l1MuTable
from PhysicsTools.NanoAOD.triggerObjects_cff import l1bits as _l1bits
from PhysicsTools.NanoAOD.triggerObjects_cff import mksel, unpackedPatTrigger
from PhysicsTools.NanoAOD.triggerObjects_cff import triggerObjectTable as _triggerObjectTable
from PhysicsTools.NanoAOD.vertices_cff import vertexTable as _vertexTable
from PhysicsTools.PatAlgos.slimming.genParticles_cff import (
    packedGenParticles,
    prunedGenParticles,
    prunedGenParticlesWithStatusOne,
)
from PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi import selectedPatTrigger
from PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi import slimmedPatTrigger
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger as _patTrigger

cscLCTRecoMuonTable = cscLCTRecoMuonTableProducer.clone()

# CSC acceptance of the event preselection: the user's definition (2026-11-02), the nominal CSC coverage. A PSet, so a
# cfg reads it from the process after process.load (importing this cff before cms.Process(name, era) exists would
# skip the era modifiers of every module it imports).
gemCscNanoSelection = cms.PSet(cscAcceptance=cms.string("abs(eta) > 0.9 && abs(eta) < 2.4"))


def gemCscNanoSelectionCuts(muonCut, cscAcceptance):
    """The cuts of gemCscNanoMuonFilter and gemCscNanoGenFilter for a muonCut (empty: all muons) and an acceptance."""
    return ("%spt > 10 && %s" % ("(%s) && " % muonCut if muonCut else "", cscAcceptance),
            "abs(pdgId) == 13 && status == 1 && pt > 10 && %s" % cscAcceptance)


def setGemCscNanoSelection(process):
    """Set both filter cuts from process.cscLCTRecoMuonTable.muonCut and process.gemCscNanoSelection: call it after
    changing either, or the event preselection keeps the default cuts."""
    process.gemCscNanoMuonFilter.cut, process.gemCscNanoGenFilter.cut = gemCscNanoSelectionCuts(
        process.cscLCTRecoMuonTable.muonCut.value(), process.gemCscNanoSelection.cscAcceptance.value())


_muonCut, _genCut = gemCscNanoSelectionCuts(cscLCTRecoMuonTable.muonCut.value(),
                                            gemCscNanoSelection.cscAcceptance.value())
gemCscNanoMuonFilter = cms.EDFilter("MuonSelector", src=cms.InputTag("muons"), cut=cms.string(_muonCut),
                                    filter=cms.bool(True))
gemCscNanoGenFilter = cms.EDFilter("GenParticleSelector", src=cms.InputTag("genParticles"), cut=cms.string(_genCut),
                                   filter=cms.bool(True))

# As miniAOD makes slimmedPatTrigger (switchOnTriggerStandAlone, packed path names), then the central TrigObj table.
patTrigger = _patTrigger.clone(onlyStandAlone=True, packTriggerPathNames=cms.bool(True))
# TrigObj: the HLT muon objects only (the central muon selection on the objects of saveTags filters), every one of
# them (not only those passing a bit), with one filter bit per single-muon path: the path accepted the event and the
# object is in its last filter (PATTriggerProducer: pathFired && lastFilter). Tag and probe on HLT_IsoMu24 uses bit 1.
GEMCSCNANO_MUON_PATHS = ["HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_CascadeMu100",
                         "HLT_HighPtTkMu100"]
triggerObjectTable = _triggerObjectTable.clone(selections=cms.PSet(Muon=_triggerObjectTable.selections.Muon.clone(
    skipObjectsNotPassingQualityBits=False,
    qualityBits=cms.VPSet(*[mksel("path('%s_v*', 1, 0)" % p, p + " accepted, object in its last filter")
                            for p in GEMCSCNANO_MUON_PATHS]))))

# HLT and L1 bits: the muon paths and seeds only. A name (after HLT_ or L1_, without _v<N>) is kept if it names a
# muon object (Mu or mu not followed by a lower-case letter, Muon, muon as in Dimuon, Csc, CSC, Cosmic, cosmic) and no
# other physics object: electron or e/gamma (Ele, EG not preceded by an upper-case letter: not the NEG of
# SingleMu22_BMTF_NEG), photon, jet, HT, MET, ETM, PNet, b tag, hadronic tau (Tau, Prong; Tau in a name with Tau3Mu or
# Tau3mu is the three-muon candidate and is kept).
_muonToken = "[Mm]u(?![a-z])|[Mm]uon|Csc|CSC|[Cc]osmic"
_otherObject = "Ele|(^|[^A-Z])EG|Photon|Pho|Jet|HT|MET|ETM|PNet|BTag|Btag|Prong|^(?!.*Tau3[Mm]u).*Tau"
gemCscNanoHLTBits = triggerResultsSelector.clone(src=cms.InputTag("TriggerResults", "", "HLT"), prefix="HLT_",
                                                 select=_muonToken, veto=_otherObject)
# the central l1bits (final decision), without the unprefireable-event bits (not muon seeds)
l1bits = _l1bits.clone(storeUnprefireableBits=False)
gemCscNanoL1Bits = triggerResultsSelector.clone(src="l1bits", prefix="L1_", select=_muonToken, veto=_otherObject)
# PV: the central NanoAOD PV columns computable on RECO, on offlinePrimaryVertices (miniAOD's offlineSlimmedPrimaryVertices
# is its copy, same order: PATVertexSlimmer.cc), with the central good-vertex cut.
primaryVertexTable = primaryVertexTableProducer.clone(src="offlinePrimaryVertices", goodPvCut=_vertexTable.goodPvCut,
                                                      name="PV")
# TMB header: the unpacker's CSCTMBStatusDigis (it makes them only with UnpackStatusDigis = True).
tmbHeaderTable = tmbHeaderTableProducer.clone(src=cms.InputTag("muonCSCDigis", "MuonCSCTMBStatusDigi"))
# EMTF: the output tracks as the uGMT reads them out (data and MC); the EMTF DAQ readout (data only).
emtfMuonTable = emtfMuonTableProducer.clone(src=cms.InputTag("gmtStage2Digis", "EMTF"))
emtfDaqTable = emtfDaqTableProducer.clone(src="emtfStage2Digis")
# Every BX of the uGMT readout: the producer clamps [minBX, maxBX] to the BXVector's range (SimpleFlatTableProducer.h).
l1MuTable = _l1MuTable.clone(minBX=-1000, maxBX=1000)

gemCscNanoTask = cms.Task(
    cscLCTRecoMuonTable,
    patTrigger,
    selectedPatTrigger,
    slimmedPatTrigger,
    unpackedPatTrigger,
    triggerObjectTable,
    gemCscNanoHLTBits,
    l1bits,
    gemCscNanoL1Bits,
    l1MuTable,
    primaryVertexTable,
    emtfMuonTable,
    tmbHeaderTable,
)

# Data only: the EMTF readout (FEDs 1384, 1385), not packed in MC RAW.
gemCscNanoDataTask = cms.Task(emtfStage2Digis, emtfDaqTable)

# MC: the central GenPart (the miniAOD pruning genParticles -> prunedGenParticlesWithStatusOne -> prunedGenParticles,
# then the NanoAOD finalGenParticles, with genIso on packedGenParticles), and the SimTracks with their links to GenPart
# and their CSC layers.
cscLCTSimTrackTable = cscLCTSimTrackTableProducer.clone()
# Pileup: the central table on RECO inputs: addPileupInfo (the input of miniAOD's slimmedAddPileupInfo, which keeps every
# BX's counts and the BX-0 details this table reads: PileupSummaryInfoSlimmer.cc; Run3Winter26 RECO has no slimmed copy)
# and offlinePrimaryVertices. NPUTablesProducer writes no columns when addPileupInfo is missing.
puTable = _puTable.clone(src="addPileupInfo", pvsrc="offlinePrimaryVertices")

gemCscNanoMCTask = cms.Task(
    prunedGenParticlesWithStatusOne,
    prunedGenParticles,
    packedGenParticles,
    finalGenParticles,
    genIso,
    genParticleTable,
    cscLCTSimTrackTable,
    puTable,
)
