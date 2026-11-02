# L1Trigger/GemCscNano

A NanoAOD-format ntuple of CSC LCTs, reco muons, the muons' CSC and GEM information, trigger information and
simulation truth, for the R&D of the LCT and GEM-CSC trigger algorithms. The ntuple stores inputs only: the LCT-muon
and LCT-SimTrack matching, with its windows and segment choices, runs when the ntuple is analysed (GEM-CSC-dev
`studies/2026-10-20_lct-recomuon-matching/nanoana/lct_matching.py`).

Inputs need reco and raw in the same event (data LCTs are unpacked from RAW): e.g. the ZMu RAW-RECO skim, or AOD with
its RAW parent. MiniAOD cannot be used (no CSC segments; the segment refs in `slimmedMuons` are not usable).

Library: `interface/CSCLayerPoint.h`, a point on a CSC layer in strip and wire-group coordinates (`fs`, `wg`) and a
straight line evaluated at a layer, used by the producers.

## NanoAOD ntuple (`python/GemCscNano_cff.py`, `test/runGemCscNano_cfg.py`)

One `Events` entry per input event, written by `NanoAODOutputModule`. No event selection and no LCT selection: every
LCT digi is stored; the reco muons are preselected by `muonCut` (cfg option of the same name; an empty cut keeps all
muons). Matching and tag and probe are applied when the ntuple is analysed. Row indices are 0-based, −1 = none.

| Collection | Content | Producer |
|---|---|---|
| `LCT` | every digi of the LCT collection, in collection order: chamber, `valid` (`isValid()`), hs, es, wg, quality, patterns, slope (`slope6` = `getSlopeEx(6)`), bend, bx, type, gemLayer | `CSCLCTRecoMuonTableProducer` |
| `Muon` | the reco muons passing `muonCut` (default `pt > 10 && (passed("CutBasedIdMedium") || passed("CutBasedIdTight"))`), in collection order (the doc string names the cut used): central NanoAOD names where the quantity is the same (pt, eta, phi, charge, isGlobal, isTracker, isStandalone, isPFcand, looseId, mediumId, tightId, pfRelIso04_all, dxy, dz, nStations), standalone-track kinematics (0 if none: `standalonePt == 0`) | `CSCLCTRecoMuonTableProducer` |
| `MuonCSCChamber` | one row per CSC chamber match (`reco::Muon::matches()`) of each stored muon, in muon order: `muonIdx`, chamber (ME1/1: normally a ring-1 = ME1b and a ring-4 = ME1a row for each crossing, with nearly the same propagated state and edges; `MuonCSCLayer_propInStripPlane` tells which ring's strip plane the track crosses); the track propagated by reco (inner track if any, else standalone; `TrackDetectorAssociator`) on the chamber surface: `x`, `y`, `dXdZ`, `dYdZ` and their errors, `edgeX`, `edgeY` (< 0: inside), `nDigisInRange`, `eta`, `phi` of that point; `nSegmentMatches` | `CSCLCTRecoMuonTableProducer` |
| `MuonCSCSegment` | one row per CSC segment of a stored muon: source 0 = a valid 4D CSC hit of its standalone track (`outerTrack()->recHits()` order), source 1 = a segment match of one of its CSC chamber matches (`matches()`, `segmentMatches` order) with reco's arbitration `mask`, link `muonCSCChamberIdx`; chamber ids (ring), nHits, chi2, ndof, time, position and direction (chamber frame) with errors | `CSCLCTRecoMuonTableProducer` |
| `MuonCSCLayer` | per stored muon, per CSC chamber (with ring) where it has a chamber match or a segment, per segment there (or none), layers 1-6: three sources of a position on the layer, each with x, y (layer frame), `fs` (`CSCLayerGeometry::strip`, clamped to [0, N], with `InStrips`, `InStripPlane`) and `wg` (0-based): the segment's rechit (`recHit*`, wire group `hitWire() - 1` and its contiguous cluster size, strips, time, errors, flags), the segment line (`line*`) and reco's propagated state carried to the layer (`prop*`); for the lines also `fs`, `wg` at the gas-gap surfaces (`*Lo`, `*Hi`) | `CSCLCTRecoMuonTableProducer` |
| `MuonCSCRecHitADC` | the strip ADC samples of the `MuonCSCLayer` rechits: one row per (rechit, strip, time bin) | `CSCLCTRecoMuonTableProducer` |
| `SimTrackCSCLayer` | MC: per stored SimTrack, per CSC chamber (with ring) where it has own-type simhits, layers 1-6: the first own-type simhit on the layer (`simHit*`: midpoint x, y, `fs`, `wg`, entry and exit `fs`, `wg`, tof, energy loss, momentum, process) and the least-squares line through the chamber's own-type simhits (`line*`, with gas-gap surfaces) | `CSCLCTSimTrackTableProducer` |
| `MuonGEMChamber` | one row per GEM chamber match (`reco::Muon::matches()`, GEM subdetector) of each stored muon, in muon order: `muonIdx`, region, station (0 = GE0), ring, layer (1-2; 0 = a super-chamber: GE0, or any station in RECO before CMSSW_12_5_0), chamber; the track propagated by reco on the chamber surface: `x`, `y`, `dXdZ`, `dYdZ` and their errors, `edgeX`, `edgeY` (< 0: inside); the eta partition (`roll`) whose bounds hold the point where the straight line through that state crosses the partition plane, with its `strip` and `pad` coordinates (clamped to 0..nstrips, 0..npads; roll 0 and −1 if none: in a gap between partitions, outside the active area, no direction, or a super-chamber row); `eta`, `phi` of the surface point; `nHitMatches` (rows of `MuonGEMHitMatch`; only if the reco step had GEM rechits), `nSegmentMatches` (GEM segments: super-chamber rows only); ME0-subdetector matches are not stored | `CSCLCTRecoMuonTableProducer` |
| `MuonGEMHitMatch` | one row per GEM rechit reco matched to a `MuonGEMChamber` row (`reco::MuonGEMHitMatch`), grouped by chamber row: `muonGEMChamberIdx`, `x` (rechit local x in its eta partition), `dx` (x − the propagated x: reco's matching variable), `bx` | `CSCLCTRecoMuonTableProducer` |
| `TrigObj` | HLT muon objects with `filterBits`, as in the central NanoAOD (Run 3: its muon selection; Phase 2: `hltPhase2L3MuonCandidates`, bits 1 = Iso, 3 = 1mu, 10 = Mu50 set by the final filters of `HLT_IsoMu24_FromL1TkMuon` / `HLT_Mu50_FromL1TkMuon`, only objects with a bit are kept), plus bit 14 = the object passed the last filter of `HLT_IsoMu24_v*` (Phase 2: `HLT_IsoMu24_FromL1TkMuon`); the central bits 1 and 3 also match the final filters of IsoMu20 and IsoMu27 | PAT trigger chain + `triggerObjectTable` |
| `L1Mu` | L1 muons (legacy GMT, BX −2..2) | central `l1MuTable` |
| `L1gmtMuon`, `L1gmtTkMuon` | Phase 2 only: Phase-2 GMT standalone and tracker muons | `DPGAnalysis/Phase2L1TNanoAOD` |
| `HLT_*`, `L1_*` | HLT path bits of the HLT process; L1 bits (Run 3 from `gtStage2Digis`, Phase 2 from the Phase-2 GT) | `NanoAODOutputModule`, `l1bits` |
| `GenPart` | MC: generator muons and all their ancestors, central GenPart columns at full precision | `GenParticlePruner` + central `genParticleTable` |
| `SimTrack` | MC: muon SimTracks and every SimTrack with a CSC simhit (signal only): kinematics, vertex, `genPartIdx` (row in GenPart), `parentIdx`, `nCSCSimHits` | `CSCLCTSimTrackTableProducer` |

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

The producer reads the GEM geometry, so the era must have GEM (Run 3, Phase 2), and the geometry must be the RECO's: a GEM chamber missing from it throws, a moved one goes unnoticed. The header of the cfg shows how to write the ntuple from emulator LCTs re-emulated from data. The cff keeps the
central labels of the PAT trigger and NanoAOD modules (`patTrigger`, `triggerObjectTable`, `l1bits`, ...), so it cannot be
loaded into a job that also runs the PAT or NANO steps. Table names have no
underscore and the link columns follow the NanoAOD `<target>Idx` convention; coffea's `NanoAODSchema` does not know the
link tables, so they are read as plain collections.

## Known limitations

- A crossed chamber is described two ways: by the edges at the chamber surface (`edgeX < 0 && edgeY < 0`, the odd-layer
  wire plane) and per layer by the propagated state (`MuonCSCLayer_propInStrips`, `propInStripPlane`, `propWg`); they do
  not always agree. In ME1/1 the edges and the propagated state are nearly the same in the ring-1 and ring-4 rows.
- Phase 2: `L1Mu` and the `l1*` columns of `TrigObj` come from the legacy GMT emulation; the Phase-2 L1 muons are in
  `L1gmtMuon` / `L1gmtTkMuon`.
