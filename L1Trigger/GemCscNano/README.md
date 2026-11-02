# L1Trigger/GemCscNano

A NanoAOD-format ntuple of CSC LCTs, reco muons, the muons' CSC and GEM information, trigger information and
simulation truth, for the R&D of the LCT and GEM-CSC trigger algorithms. No LCT matching and no tag-and-probe
selection is applied.

Inputs need reco and raw in the same event (data LCTs are unpacked from RAW): e.g. the ZMu RAW-RECO skim, or AOD with
its RAW parent (not yet run). MiniAOD cannot be used (no CSC segments; the segment refs in `slimmedMuons` are not usable).

Library: `interface/CSCLayerPoint.h`, a point on a CSC layer in strip and wire-group coordinates (`fs`, `wg`) and a
straight line evaluated at a layer, used by the producers.

## NanoAOD ntuple (`python/GemCscNano_cff.py`, `test/runGemCscNano_cfg.py`)

One `Events` entry per selected event, written by `NanoAODOutputModule`. Events: those with a stored muon (`muonCut`)
with pT > 10 GeV and 0.9 < |eta| < 2.4, or in MC a final-state generator muon with pT > 10 GeV in that range
(`gemCscNanoMuonFilter`, `gemCscNanoGenFilter`). The reco muons are preselected by `muonCut` (cfg option of the same
name; an empty cut keeps all muons); every LCT digi is stored. After changing `muonCut` or `gemCscNanoSelection`
in a cfg, call `setGemCscNanoSelection(process)` of the cff, so that the event preselection follows. Run 3 only (the cfg
refuses other eras and requires `era=`). Row indices are 0-based, −1 = none.

| Collection | Content | Producer |
|---|---|---|
| `LCT` | every digi of the LCT collection, in collection order: chamber, `valid` (`isValid()`), hs, es, wg, quality, patterns, slope (`slope6` = `getSlopeEx(6)`), bend, bx, type, gemLayer | `CSCLCTRecoMuonTableProducer` |
| `Muon` | the reco muons passing `muonCut` (default `pt > 10 && (passed("CutBasedIdMedium") || passed("CutBasedIdTight"))`), in collection order (the doc string names the cut used): central NanoAOD names where the quantity is the same (pt, eta, phi, charge, isGlobal, isTracker, isStandalone, isPFcand, looseId, mediumId, tightId, pfRelIso04_all, dxy, dz, nStations), standalone-track kinematics (0 if none: `standalonePt == 0`) | `CSCLCTRecoMuonTableProducer` |
| `MuonCSCChamber` | one row per CSC chamber match (`reco::Muon::matches()`) of each stored muon, in muon order: `muonIdx`, chamber (ME1/1: normally a ring-1 = ME1b and a ring-4 = ME1a row for each crossing, with nearly the same propagated state and edges; `MuonCSCLayer_propInStripPlane` tells which ring's strip plane the track crosses); the track propagated by reco (inner track if any, else standalone; `TrackDetectorAssociator`) on the chamber surface: `x`, `y`, `dXdZ`, `dYdZ` and their errors, `edgeX`, `edgeY` (< 0: inside), `nDigisInRange`, `eta`, `phi` of that point; `nSegmentMatches` | `CSCLCTRecoMuonTableProducer` |
| `MuonCSCSegment` | one row per CSC segment of a stored muon: source 0 = a valid 4D CSC hit of its standalone track (`outerTrack()->recHits()` order), source 1 = a segment match of one of its CSC chamber matches (`matches()`, `segmentMatches` order) with reco's arbitration `mask`, link `muonCSCChamberIdx`; chamber ids (ring), nHits, chi2, ndof, time, position and direction (chamber frame) with errors | `CSCLCTRecoMuonTableProducer` |
| `MuonCSCLayer` | per stored muon, per CSC chamber (with ring) where it has a chamber match or a segment, per segment there (or none), layers 1-6: three sources of a position on the layer, each with x, y (layer frame), `fs` (`CSCLayerGeometry::strip`, clamped to [0, N], with `InStrips`, `InStripPlane`) and `wg` (0-based): the segment's rechit (`recHit*`, wire group `hitWire() - 1` and its contiguous cluster size, strips, time, errors, flags), the segment line (`line*`) and reco's propagated state carried to the layer (`prop*`); for the lines also `fs`, `wg` at the gas-gap surfaces (`*Lo`, `*Hi`) | `CSCLCTRecoMuonTableProducer` |
| `MuonCSCRecHitADC` | the strip ADC samples of the `MuonCSCLayer` rechits (`CSCRecHit2D::adcs`): one row per (rechit, strip, sample); samples 0-3 = SCA time bins tmax−1 … tmax+2 of the cluster's central strip (tmax not stored; 0.1 = a sample reco did not fill), ADC pedestal-subtracted and gain-corrected, shared strips split between clusters | `CSCLCTRecoMuonTableProducer` |
| `SimTrackCSCLayer` | MC: per stored SimTrack, per CSC chamber (with ring) where it has own-type simhits, layers 1-6: the first own-type simhit on the layer (`simHit*`: midpoint x, y, `fs`, `wg`, entry and exit `fs`, `wg`, tof, energy loss, momentum, process) and the least-squares line through the chamber's own-type simhits (`line*`, with gas-gap surfaces) | `CSCLCTSimTrackTableProducer` |
| `MuonGEMChamber` | one row per GEM chamber match (`reco::Muon::matches()`, GEM subdetector) of each stored muon, in muon order: `muonIdx`, region, station (0 = GE0), ring, layer (1-2; 0 = a super-chamber: GE0, or any station in RECO before CMSSW_12_5_0), chamber; the track propagated by reco on the chamber surface: `x`, `y`, `dXdZ`, `dYdZ` and their errors, `edgeX`, `edgeY` (< 0: inside); the eta partition (`roll`) whose bounds hold the point where the straight line through that state crosses the partition plane, with its `strip` and `pad` coordinates (clamped to 0..nstrips, 0..npads; roll 0 and −1 if none: in a gap between partitions, outside the active area, no direction, or a super-chamber row); `eta`, `phi` of the surface point; `nHitMatches` (rows of `MuonGEMHitMatch`; only if the reco step had GEM rechits), `nSegmentMatches` (GEM segments: super-chamber rows only); ME0-subdetector matches are not stored | `CSCLCTRecoMuonTableProducer` |
| `MuonGEMHitMatch` | one row per GEM rechit reco matched to a `MuonGEMChamber` row (`reco::MuonGEMHitMatch`), grouped by chamber row: `muonGEMChamberIdx`, `x` (rechit local x in its eta partition), `dx` (x − the propagated x: reco's matching variable), `bx` | `CSCLCTRecoMuonTableProducer` |
| `TrigObj` | HLT objects with `filterBits`, the central NanoAOD table unmodified (all its object types and bits) | PAT trigger chain + central `triggerObjectTable` |
| `L1Mu` | L1 muons of the uGMT readout (`gmtStage2Digis:Muon`), every BX of the readout | central `l1MuTable` |
| `HLT_*`, `L1_*` | HLT path bits of the HLT process; L1 bits from `gtStage2Digis`, with the central unprefireable-event bits (TCDS record unpacked from RAW; `L1_UnprefireableEvent_TriggerRules` is false in MC, where `L1TExtCondProducer` skips the TCDS check) | `NanoAODOutputModule`, central `l1bits` |
| `GenPart` | MC: the central NanoAOD GenPart (miniAOD pruning, then `finalGenParticles`, with `iso`) | central chain, `genParticleTable` |
| `SimTrack` | MC: muon SimTracks and every SimTrack with a CSC simhit (signal only): kinematics, vertex, `genPartIdx` (row in GenPart), `parentIdx`, `nCSCSimHits` | `CSCLCTSimTrackTableProducer` |

```bash
# data: ZMu RAW-RECO (LCTs unpacked here)
cmsRun L1Trigger/GemCscNano/test/runGemCscNano_cfg.py inputFiles=file:<file> \
  globalTag=<GT of the runs> era=<era of the runs, e.g. Run3_2026>
# Run-3 MC: GEN-SIM-RECO with its GEN-SIM-(DIGI-)RAW parent (LCTs unpacked from its RAW)
cmsRun ... mc=True inputFiles=file:<GEN-SIM-RECO> secondaryInputFiles=file:<RAW parent> \
  globalTag=<GT of the sample> era=Run3_2026
```

The producer reads the GEM geometry, and the geometry must be the RECO's: a GEM chamber missing from it throws, a moved one goes unnoticed. The header of the cfg shows how to write the ntuple from emulator LCTs re-emulated from data. The cff keeps the
central labels of the PAT trigger and NanoAOD modules (`patTrigger`, `triggerObjectTable`, `l1bits`, ...), so it cannot be
loaded into a job that also runs the PAT or NANO steps. Table names have no
underscore and the link columns follow the NanoAOD `<target>Idx` convention; coffea's `NanoAODSchema` does not know the
link tables, so they are read as plain collections.

## Known limitations

- A crossed chamber is described two ways: by the edges at the chamber surface (`edgeX < 0 && edgeY < 0`, the odd-layer
  wire plane) and per layer by the propagated state (`MuonCSCLayer_propInStrips`, `propInStripPlane`, `propWg`); they do
  not always agree. In ME1/1 the edges and the propagated state are nearly the same in the ring-1 and ring-4 rows.
