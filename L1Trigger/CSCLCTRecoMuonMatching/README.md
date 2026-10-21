# L1Trigger/CSCLCTRecoMuonMatching

Association of CSC LCTs with reco muons through the muons' own CSC segments: no propagation, no ΔR.

## Library

```cpp
#include "L1Trigger/CSCLCTRecoMuonMatching/interface/CSCLCTRecoMuonMatcher.h"

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

BuildFile: `<use name="L1Trigger/CSCLCTRecoMuonMatching"/>`.

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
- Muon selection is left to the caller (all muons of the collection are matched).
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

## Test analyzer

`CSCLCTRecoMuonMatchAnalyzer` (cfi `L1Trigger.CSCLCTRecoMuonMatching.cscLCTRecoMuonMatchAnalyzer_cfi`) writes the
tree `lctMuon`: one entry per LCT (`lct_*`, `nMuons`) with vector branches for its matched muons (`mu_*`, `seg_*`,
`key_*`, `dHS`, `dES`, `dWG`).

```bash
cmsRun L1Trigger/CSCLCTRecoMuonMatching/test/runCSCLCTRecoMuonMatch_cfg.py \
  inputFiles=file:<ZMu RAW-RECO file> globalTag=<GT of the runs> era=Run3_2026 maxEvents=1000
```

The header of the cfg shows how to run it on emulator LCTs (GT override and emulator options included).

## Known limitations (smoke test on one ZMu RAW-RECO file of run 403567 with ±2 HS / ±2 WG windows, provisional)

- `standalone` and `matches` pick the same segment in 95.8% of the chambers both have.
- In ME+1/1a and ME−1/1b mean ΔES is ≈ 0.5 eighth-strip lower than in ME+1/1b and ME−1/1a (unexplained).
- About 0.35% of segments have no valid key position (no layer-3 rechit and a projection just outside the wire
  plane); they cannot match.
- GEM-only LCT types (CLCT-2GEM, ALCT-2GEM) take a coordinate from GEM: expect larger residuals for them.
- A crossed chamber without a segment gives no candidate: an LCT-efficiency denominator needs the chamber matches
  without segments, which this matcher does not provide.
