# L1Trigger/GemCscNano

Association of CSC LCTs with reco muons through the muons' own CSC segments (no propagation, no ΔR), and a NanoAOD-format
ntuple of LCTs, reco muons, trigger information and simulation truth with the associations as index columns.

## Library

```cpp
#include "L1Trigger/GemCscNano/interface/CSCLCTRecoMuonMatcher.h"

// in the module constructor: matcher_(ps.getParameter<edm::ParameterSet>("matcher")) and, in fillDescriptions,
// CSCLCTRecoMuonMatcher::fillPSetDescription(matcherDesc); desc.add<edm::ParameterSetDescription>("matcher", matcherDesc);
const auto matches = matcher_.match(lcts, muons, cscGeometry);  // CSCCorrelatedLCTDigiCollection, reco::MuonCollection
for (const auto& lm : matches) {            // one entry per valid LCT
  // lm.chamber, *lm.lct
  for (const auto& mm : lm.muons) {         // matched muons, closest first (possibly none, possibly several)
    // muons[mm.muonIndex], mm.segment, mm.deltaHalfStrip, mm.deltaEighthStrip, mm.deltaWireGroup, mm.key ...
  }
}
```

BuildFile: `<use name="L1Trigger/GemCscNano"/>`.

- Each CSC segment of each muon is placed on the key layer (layer 3) in LCT units: key half-strip (0-based; ME1/1:
  ME1b 0-127, ME1a 128-223), eighth-strip (CCLUT convention, es = 4·hs at the half-strip centre), key wire group
  (0-based). The position comes from the segment's layer-3 rechit when it has a valid one, otherwise from the
  segment projected onto layer 3, always with the geometry of the segment's own ring.
- An LCT is matched to every muon with a segment in its chamber within |ΔHS| ≤ `maxDeltaHalfStrip` and
  |ΔWG| ≤ `maxDeltaWireGroup` (Δ = LCT − segment). Matching is many-to-many: an LCT can have several muons (e.g.
  muons sharing a segment), a muon can match several LCTs, also in one chamber. Per LCT, a muon keeps its closest
  segment. ME1/1 ring-1 (ME1b) and ring-4 (ME1a) segments both map to the ring-1 LCT chamber.
- `segmentSource`: `standalone` (segments of the standalone-muon fit), `matches` (the segment arbitrated to the
  muon in each crossed chamber), `auto` (standalone when available, otherwise matches).
- `lctBXWindow`: keep LCTs with |BX − 8| ≤ window; < 0 = no selection. Data LCTs carry a 1-bit BX, so keep it < 0
  for them; emulator LCTs span the readout window around BX 8.
- Muon selection is left to the caller of the library (all muons of the collection given to `match()` are matched);
  the ntuple producer passes only the muons of its `muonCut`.
- `Config` has no default values. The defaults of `fillPSetDescription` are ±6 half-strips and ±2 wire groups
  (deliberately loose, to be tightened from the residual distributions), `segmentSource` `auto` (a proposal of the
  matching study) and `lctBXWindow` −1 (no BX selection, as data LCTs need). Negative windows are rejected.
- The LCT pointers in the result refer to the LCT collection, the segment pointers to the event products reached
  through the muons (standalone-muon rechits or `cscSegments`); they are valid as long as those products are.
- ME1/1: the half-strip window spans the ME1b/ME1a numbering boundary (127/128), so an ME1b segment near half-strip
  127 can match an ME1a LCT near 128 and vice versa; compare `key.ring` with the LCT half-strip (< 128: ME1b) to keep
  same-ring pairs only.

Inputs need reco and raw in the same event (data LCTs are unpacked from RAW): e.g. the ZMu RAW-RECO skim, or AOD with
its RAW parent. MiniAOD cannot be used (no CSC segments; the segment refs in `slimmedMuons` are not usable).

## NanoAOD ntuple (`python/GemCscNano_cff.py`, `test/runGemCscNano_cfg.py`)

One `Events` entry per input event, written by `NanoAODOutputModule`. No event selection, and no LCT selection beyond
`lctBXWindow`; the reco muons are preselected by `muonCut` (cfg option of the same name; an empty cut keeps all muons). Tag and probe is applied when
the ntuple is analysed. Row indices are 0-based, −1 = none.

| Collection | Content | Producer |
|---|---|---|
| `LCT` | every LCT the matcher considers (valid, inside `lctBXWindow`), in collection order: chamber, hs, es, wg, quality, patterns, slope (`slope6` = `getSlopeEx(6)`), bend, bx, type, gemLayer, `nMuonMatch`, (MC) `nSimTrackMatch` | `CSCLCTRecoMuonTableProducer` |
| `Muon` | the reco muons passing `muonCut` (default `pt > 10 && (passed("CutBasedIdMedium") || passed("CutBasedIdTight"))`), in collection order; only they are matched: central NanoAOD names where the quantity is the same (pt, eta, phi, charge, isGlobal, isTracker, isStandalone, isPFcand, looseId, mediumId, tightId, pfRelIso04_all, dxy, dz, nStations), standalone-track kinematics (0 if none: `standalonePt == 0`), `nCSCSegments`, `nLCTMatch` | `CSCLCTRecoMuonTableProducer` |
| `LCTMuonMatch` | one row per matched LCT–muon pair, grouped by LCT, closest muon first: `lctIdx`, `muonIdx`, the segment (source, ring, hits, χ²/ndof, time), its key-layer position and the residuals `dHS`, `dES`, `dWG` | `CSCLCTRecoMuonTableProducer` |
| `MuonCSCChamber` | one row per CSC chamber match (`reco::Muon::matches()`) of each stored muon, in muon order: `muonIdx`, chamber (ME1/1: normally a ring-1 = ME1b and a ring-4 = ME1a row for each crossing, with nearly the same propagated state and edges; the ring crossed is the row with `keyInStripPlane`, neither beyond the outer ends of ME1/1); the track propagated by reco (inner track if any, else standalone; `TrackDetectorAssociator`) on the chamber surface: `x`, `y`, `dXdZ`, `dYdZ` and their errors, `edgeX`, `edgeY` (< 0: inside), `nDigisInRange`; that position carried to layer 3 in LCT units: `keyHS`, `keyES`, `keyWG`, `keyValid`, `keyInStripPlane`, `eta`, `phi`; the segments there: `nSegmentMatches`, `hasMatcherSegment` (the configured `segmentSource`), `hasStandaloneSegment`, `hasArbitratedSegment` | `CSCLCTRecoMuonTableProducer` |
| `MuonGEMChamber` | one row per GEM chamber match (`reco::Muon::matches()`, GEM subdetector) of each stored muon, in muon order: `muonIdx`, region, station (0 = GE0), ring, layer (0 = super-chamber), chamber; the track propagated by reco on the chamber surface: `x`, `y`, `dXdZ`, `dYdZ` and their errors, `edgeX`, `edgeY` (< 0: inside); the eta partition whose bounds hold that point (`roll`, 0 if none: outside the active area or a super-chamber row) with its `strip` and `pad` coordinates (−1 if none), `eta`, `phi`; `nHitMatches` (rows of `MuonGEMHitMatch`), `nSegmentMatches` (GEM segments, reco's segment path) | `CSCLCTRecoMuonTableProducer` |
| `MuonGEMHitMatch` | one row per GEM rechit reco matched to a `MuonGEMChamber` row (`reco::MuonGEMHitMatch`), grouped by chamber row: `muonGEMChamberIdx`, `x` (rechit local x in its eta partition), `dx` (x − the propagated x: reco's matching variable), `bx` | `CSCLCTRecoMuonTableProducer` |
| `TrigObj` | HLT muon objects with `filterBits`, as in the central NanoAOD (Run 3: its muon selection; Phase 2: `hltPhase2L3MuonCandidates`, bits 1 = Iso, 3 = 1mu, 10 = Mu50 set by the final filters of `HLT_IsoMu24_FromL1TkMuon` / `HLT_Mu50_FromL1TkMuon`, only objects with a bit are kept), plus bit 14 = the object passed the last filter of `HLT_IsoMu24_v*` (Phase 2: `HLT_IsoMu24_FromL1TkMuon`); the central bits 1 and 3 also match the final filters of IsoMu20 and IsoMu27 | PAT trigger chain + `triggerObjectTable` |
| `L1Mu` | L1 muons (legacy GMT, BX −2..2) | central `l1MuTable` |
| `L1gmtMuon`, `L1gmtTkMuon` | Phase 2 only: Phase-2 GMT standalone and tracker muons | `DPGAnalysis/Phase2L1TNanoAOD` |
| `HLT_*`, `L1_*` | HLT path bits of the HLT process; L1 bits (Run 3 from `gtStage2Digis`, Phase 2 from the Phase-2 GT) | `NanoAODOutputModule`, `l1bits` |
| `GenPart` | MC: generator muons and all their ancestors, central GenPart columns at full precision | `GenParticlePruner` + central `genParticleTable` |
| `SimTrack` | MC: muon SimTracks and every SimTrack with a CSC simhit (signal only): kinematics, vertex, `genPartIdx` (row in GenPart), `parentIdx`, `nCSCSimHits`, `nLCTMatch` | `CSCLCTSimTrackTableProducer` |
| `LCTSimTrackMatch` | MC: one row per LCT–SimTrack pair within the same windows as the muons (default ±6 HS / ±2 WG), grouped by LCT, closest first: `lctIdx`, `simTrackIdx`, the simhit key-layer position (the track's own layer-3 simhit, else the line through its simhits in the chamber) and the residuals | `CSCLCTSimTrackTableProducer` |

```bash
# data: ZMu RAW-RECO (LCTs unpacked here)
cmsRun L1Trigger/GemCscNano/test/runGemCscNano_cfg.py inputFiles=file:<file> \
  globalTag=<GT of the runs> era=Run3_2026 maxEvents=1000
# Run-3 MC: GEN-SIM-RECO with its GEN-SIM-DIGI-RAW parent (emulator LCTs simCscTriggerPrimitiveDigis)
cmsRun ... mc=True inputFiles=file:<GEN-SIM-RECO> secondaryInputFiles=file:<GEN-SIM-DIGI-RAW parent> \
  globalTag=<GT of the sample> era=Run3_2026
# central Run-3 MC whose GEN-SIM-RAW has no emulator LCTs (e.g. Run3Winter25 DY): unpack them from its RAW
cmsRun ... mc=True lcts=muonCSCDigis:MuonCSCCorrelatedLCTDigi inputFiles=file:<GEN-SIM-RECO> \
  secondaryInputFiles=file:<GEN-SIM-RAW parent> globalTag=<GT of the sample> era=Run3_2025
# Phase-2 MC: GEN-SIM-RECO of the Phase-2 RelVals carries simhits and LCTs
cmsRun ... mc=True inputFiles=file:<GEN-SIM-RECO> globalTag=<GT> era=Phase2C22I13M9 geometry=ExtendedRun4D121
```

The header of the cfg shows how to write the ntuple from emulator LCTs re-emulated from data. The cff keeps the
central labels of the PAT trigger and NanoAOD modules (`patTrigger`, `triggerObjectTable`, `l1bits`, ...), so it cannot be
loaded into a job that also runs the PAT or NANO steps. Table names have no
underscore and the link columns follow the NanoAOD `<target>Idx` convention; coffea's `NanoAODSchema` does not know the
link tables, so they are read as plain collections.

## Known limitations

- Smoke test on one ZMu RAW-RECO file of run 403567 (±2 HS / ±2 WG windows, provisional): `standalone` and `matches`
  pick the same segment in 95.8% of the chambers both have; in ME+1/1a and ME−1/1b mean ΔES is ≈ 0.5 eighth-strip
  lower than in ME+1/1b and ME−1/1a; about 0.35% of segments have no valid key position (no layer-3 rechit and a
  projection just outside the wire plane) and cannot match.
- GEM-only LCT types (CLCT-2GEM, ALCT-2GEM) take a coordinate from GEM: expect larger residuals for them.
- A crossed chamber without a segment gives no candidate. The matcher does not use the propagated position; the
  ntuple stores it per chamber in `MuonCSCChamber`, the input for an LCT-efficiency denominator. Two descriptions of
  "crossed" are available: the edges at the chamber surface (`edgeX < 0 && edgeY < 0`, the odd-layer wire plane) and
  the layer-3 point (`keyValid`, `keyInStripPlane`), and they do not always agree. In ME1/1 the edges and
  the propagated state are nearly the same in the ring-1 and ring-4 rows, so count one crossing per LCT chamber, take the
  ring from `keyInStripPlane`. Which definition the denominator uses is not decided.
- Phase 2: `L1Mu` and the `l1*` columns of `TrigObj` come from the legacy GMT emulation; the Phase-2 L1 muons are in
  `L1gmtMuon` / `L1gmtTkMuon`.
