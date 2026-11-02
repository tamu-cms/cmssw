# NanoAOD-format ntuple of CSC LCTs, reco muons, trigger information and simulation truth, one entry per event (Run 3).
# Events: those with a stored muon (cscLCTRecoMuonTable.muonCut) with pT > 10 GeV and 0.9 < |eta| < 2.4, or (MC) a
# final-state generator muon with pT > 10 GeV and 0.9 < |eta| < 2.4 (gemCscNanoMuonFilter, gemCscNanoGenFilter; the
# user's definition of 2026-11-02). Per event: every LCT; the reco muons passing muonCut (default pT > 10 GeV and the
# medium or tight cut-based ID). LCT matching and tag-and-probe selections are applied when the ntuple is analysed.
#   LCT, Muon, MuonCSCChamber,          cscLCTRecoMuonTable (CSCLCTRecoMuonTableProducer)
#   MuonCSCSegment, MuonCSCLayer,
#   MuonCSCRecHitADC, MuonGEMChamber,
#   MuonGEMHitMatch
#   TrigObj, L1Mu                       the central NanoAOD tables (TrigObj from the PAT trigger chain of miniAOD, run
#                                       here on hltTriggerSummaryAOD; L1Mu with every BX of the uGMT readout)
#   HLT_* and L1_* bits                 TriggerResults of the HLT process and the central l1bits, written by
#                                       NanoAODOutputModule
#   GenPart, SimTrack, SimTrackCSCLayer gemCscNanoMCTask (MC only; GenPart as in the central NanoAOD)
# Module labels end in "Table" where NanoAODOutputModule must keep the product (keep nanoaodFlatTable_*Table_*_*).
# The PAT trigger, gen-particle and NanoAOD modules keep their central labels (patTrigger, prunedGenParticles,
# finalGenParticles, triggerObjectTable, l1bits, ...): this cff cannot be loaded into a job that also runs the PAT or
# NANO steps.
import FWCore.ParameterSet.Config as cms

from EventFilter.OnlineMetaDataRawToDigi.tcdsRawToDigi_cfi import tcdsRawToDigi as _tcdsRawToDigi
from L1Trigger.GemCscNano.cscLCTRecoMuonTableProducer_cfi import cscLCTRecoMuonTableProducer
from L1Trigger.GemCscNano.cscLCTSimTrackTableProducer_cfi import cscLCTSimTrackTableProducer
from L1Trigger.L1TGlobal.simGtExtFakeProd_cfi import simGtExtUnprefireable
from PhysicsTools.NanoAOD.genparticles_cff import finalGenParticles, genIso, genParticleTable
from PhysicsTools.NanoAOD.l1trig_cff import l1MuTable as _l1MuTable
from PhysicsTools.NanoAOD.triggerObjects_cff import l1bits, triggerObjectTable, unpackedPatTrigger
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
cscAcceptance = gemCscNanoSelection.cscAcceptance.value()
gemCscNanoMuonFilter = cms.EDFilter(
    "MuonSelector",
    src=cms.InputTag("muons"),
    cut=cms.string("(%s) && pt > 10 && %s" % (cscLCTRecoMuonTable.muonCut.value(), cscAcceptance)),
    filter=cms.bool(True),
)
gemCscNanoGenFilter = cms.EDFilter(
    "GenParticleSelector",
    src=cms.InputTag("genParticles"),
    cut=cms.string("abs(pdgId) == 13 && status == 1 && pt > 10 && %s" % cscAcceptance),
    filter=cms.bool(True),
)

# As miniAOD makes slimmedPatTrigger (switchOnTriggerStandAlone, packed path names), then the central TrigObj table.
patTrigger = _patTrigger.clone(onlyStandAlone=True, packTriggerPathNames=cms.bool(True))
# The central l1bits with its unprefireable-event bits: simGtExtUnprefireable reads the TCDS record, unpacked here from
# RAW as in the central RawToDigi (Configuration/StandardSequences/python/RawToDigi_cff.py).
tcdsDigis = _tcdsRawToDigi.clone()
# Every BX of the uGMT readout: the producer clamps [minBX, maxBX] to the BXVector's range (SimpleFlatTableProducer.h).
l1MuTable = _l1MuTable.clone(minBX=-1000, maxBX=1000)

gemCscNanoTask = cms.Task(
    cscLCTRecoMuonTable,
    patTrigger,
    selectedPatTrigger,
    slimmedPatTrigger,
    unpackedPatTrigger,
    triggerObjectTable,
    tcdsDigis,
    simGtExtUnprefireable,
    l1bits,
    l1MuTable,
)

# MC: the central GenPart (the miniAOD pruning genParticles -> prunedGenParticlesWithStatusOne -> prunedGenParticles,
# then the NanoAOD finalGenParticles, with genIso on packedGenParticles), and the SimTracks with their links to GenPart
# and their CSC layers.
cscLCTSimTrackTable = cscLCTSimTrackTableProducer.clone()

gemCscNanoMCTask = cms.Task(
    prunedGenParticlesWithStatusOne,
    prunedGenParticles,
    packedGenParticles,
    finalGenParticles,
    genIso,
    genParticleTable,
    cscLCTSimTrackTable,
)
