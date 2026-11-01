# NanoAOD-format ntuple of CSC LCTs, reco muons, trigger information and simulation truth, one entry per event. Every
# event and every LCT is stored, and the reco muons passing cscLCTRecoMuonTable.muonCut (default pT > 10 GeV and the
# medium or tight cut-based ID); tag-and-probe selections are applied when the ntuple is analysed.
#   LCT, Muon, LCTMuonMatch             cscLCTRecoMuonTable (CSCLCTRecoMuonTableProducer)
#   TrigObj (muon objects only), L1Mu   the central NanoAOD tables; TrigObj needs the PAT trigger chain of miniAOD, run
#                                       here on hltTriggerSummaryAOD
#   HLT_* and L1_* bits                 TriggerResults of the HLT process and l1bits, written by NanoAODOutputModule
#   GenPart, SimTrack, LCTSimTrackMatch gemCscNanoMCTask (MC only)
# Module labels end in "Table" where NanoAODOutputModule must keep the product (keep nanoaodFlatTable_*Table_*_*).
# The PAT trigger and NanoAOD modules keep their central labels (patTrigger, triggerObjectTable, l1bits, ...): this cff
# cannot be loaded into a job that also runs the PAT or NANO steps.
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Modifier_phase2_trigger_cff import phase2_trigger
from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2GTtables_cff import l1tP2GTTrigConvert as _l1tP2GTTrigConvert
from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nanotables_cff import gmtTkMuTable, staMuTable
from L1Trigger.GemCscNano.cscLCTRecoMuonTableProducer_cfi import cscLCTRecoMuonTableProducer
from L1Trigger.GemCscNano.cscLCTSimTrackTableProducer_cfi import cscLCTSimTrackTableProducer
from PhysicsTools.NanoAOD.genparticles_cff import genParticleTable as _genParticleTable
from PhysicsTools.NanoAOD.l1trig_cff import l1MuTable
from PhysicsTools.NanoAOD.triggerObjects_cff import l1bits as _l1bits
from PhysicsTools.NanoAOD.triggerObjects_cff import mksel, unpackedPatTrigger
from PhysicsTools.NanoAOD.triggerObjects_cff import triggerObjectTable as _triggerObjectTable
from PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi import selectedPatTrigger
from PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi import slimmedPatTrigger
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger as _patTrigger

cscLCTRecoMuonTable = cscLCTRecoMuonTableProducer.clone()

# As miniAOD makes slimmedPatTrigger (switchOnTriggerStandAlone, packed path names), then the central TrigObj table
# restricted to its muon selection, plus bit 14: the object passed the last filter of HLT_IsoMu24 (the central bits
# 1 = Iso and 3 = 1mu also match the final filters of HLT_IsoMu20 and HLT_IsoMu27).
patTrigger = _patTrigger.clone(onlyStandAlone=True, packTriggerPathNames=cms.bool(True))
triggerObjectTable = _triggerObjectTable.clone(selections=cms.PSet(Muon=_triggerObjectTable.selections.Muon.clone()))
triggerObjectTable.selections.Muon.qualityBits.append(
    mksel("path('HLT_IsoMu24_v*', 1, 0)", "last filter of HLT_IsoMu24")
)
# The unprefireable-event bit needs simGtExtUnprefireable, which only miniAOD makes; without it the bit is always false.
l1bits = _l1bits.clone(storeUnprefireableBits=False)

gemCscNanoTask = cms.Task(
    cscLCTRecoMuonTable,
    patTrigger,
    selectedPatTrigger,
    slimmedPatTrigger,
    unpackedPatTrigger,
    triggerObjectTable,
    l1bits,
    l1MuTable,
)

# Phase 2: the central muon selection uses Run-3 HLT collections and filters. The Phase-2 menu (HLT_75e33) keeps its L3
# muons in hltPhase2L3MuonCandidates; the bits keep the central positions (1 = Iso, 3 = 1mu, 10 = 1mu (Mu50),
# 14 = last filter of the IsoMu24 path) with the final filters of HLT_IsoMu24_FromL1TkMuon and HLT_Mu50_FromL1TkMuon;
# only objects with one of these bits are stored. L1 bits from the Phase-2 GT and the Phase-2 GMT muons (L1gmtMuon,
# L1gmtTkMuon) are added; L1Mu and the TrigObj l1 columns stay those of the legacy GMT emulation.
_none = mksel("0", "None")
phase2_trigger.toModify(
    triggerObjectTable.selections.Muon,
    sel="type(83) && pt > 5 && coll('hltPhase2L3MuonCandidates')",
    qualityBits=cms.VPSet(
        _none,
        mksel("filter('hltL3crIso*IsoRegionalNewFiltered0p07*')", "Iso"),
        _none,
        mksel("filter('hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk')", "1mu"),
        _none,
        _none,
        _none,
        _none,
        _none,
        _none,
        mksel("filter('hltL3fL1TkSingleMu22L3Filtered50Q')", "1mu (Mu50)"),
        _none,
        _none,
        _none,
        mksel("path('HLT_IsoMu24_FromL1TkMuon', 1, 0)", "last filter of HLT_IsoMu24_FromL1TkMuon"),
    ),
)
phase2_trigger.toReplaceWith(l1bits, _l1tP2GTTrigConvert.clone())
phase2_trigger.toReplaceWith(gemCscNanoTask, gemCscNanoTask.copyAndAdd(staMuTable, gmtTkMuTable))

# MC: the generator muons with all their ancestors, written as the central GenPart table (without the miniAOD-only
# isolation and at full float precision), and the SimTracks with their links to GenPart and to the LCTs.
cscLCTGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src=cms.InputTag("genParticles"),
    select=cms.vstring("drop *", "++keep abs(pdgId) == 13"),
)
genParticleTable = _genParticleTable.clone(src="cscLCTGenParticles", externalVariables=cms.PSet())
for _var in ("pt", "eta", "phi", "mass"):
    getattr(genParticleTable.variables, _var).precision = -1
cscLCTSimTrackTable = cscLCTSimTrackTableProducer.clone(
    lctTable=cms.InputTag("cscLCTRecoMuonTable", "lct"),
    genPartAssociation="cscLCTGenParticles",
)

gemCscNanoMCTask = cms.Task(cscLCTGenParticles, genParticleTable, cscLCTSimTrackTable)
