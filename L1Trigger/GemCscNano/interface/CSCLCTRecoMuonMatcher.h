#ifndef L1Trigger_GemCscNano_CSCLCTRecoMuonMatcher_h
#define L1Trigger_GemCscNano_CSCLCTRecoMuonMatcher_h

// Association of CSC LCTs with reco muons through the muons' own CSC segments: no propagation, no deltaR.
//
// The position of each CSC segment of each muon on the key layer (layer 3) is expressed in LCT channel space (key
// half-strip, eighth-strip, key wire group) and compared with the LCTs of the segment's LCT chamber (ME1/1: ring-1
// and ring-4 segments both map to the ring-1 LCT chamber). The position comes from the segment's own layer-3
// rechit when it has a valid one, otherwise from the segment projected onto layer 3. An LCT is matched to every
// muon with a segment within |dHS| <= maxDeltaHalfStrip and |dWG| <= maxDeltaWireGroup (the muon's closest such
// segment is kept), so an LCT can be matched to several muons (e.g. muons sharing a segment) and a muon to several
// LCTs (also in the same chamber).
// In ME1/1 the half-strip window spans the ME1b/ME1a numbering boundary (127/128), so an ME1b segment near half-strip
// 127 can match an ME1a LCT near 128 and vice versa; compare MuonMatch::key.ring with the LCT half-strip (< 128: ME1b)
// to keep same-ring pairs only.
//
// Inputs need reco and raw in the same event (data LCTs are unpacked from RAW): e.g. the ZMu RAW-RECO skim, or AOD
// with its RAW parent. MiniAOD cannot be used (no CSC segments; slimmedMuons segment refs are not usable).
// Data LCTs carry a 1-bit BX: keep lctBXWindow < 0 for them; for emulator LCTs (readout window around BX 8)
// lctBXWindow selects |BX - 8| <= lctBXWindow.

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "L1Trigger/GemCscNano/interface/CSCLCTKeyLayer.h"

#include <utility>
#include <vector>

class CSCGeometry;
namespace edm {
  class ParameterSet;
  class ParameterSetDescription;
}  // namespace edm

class CSCLCTRecoMuonMatcher {
public:
  enum class SegmentSource {
    Standalone,  // the CSC segments used in the standalone-muon fit (outerTrack() 4D hits)
    Matches,     // the segment arbitrated to the muon in each crossed chamber (reco::Muon::matches())
    Auto         // Standalone when the standalone track has at least one CSC segment, otherwise Matches
  };

  // All values must be chosen by the caller. The defaults of fillPSetDescription: windows +-6 half-strips and
  // +-2 wire groups (deliberately loose, to be tightened), segmentSource auto (a study proposal), lctBXWindow -1
  // (no BX selection, as data LCTs need). The matcher rejects negative windows.
  struct Config {
    Config(SegmentSource segmentSource, int maxDeltaHalfStrip, int maxDeltaWireGroup, int lctBXWindow);
    SegmentSource segmentSource;
    int maxDeltaHalfStrip;
    int maxDeltaWireGroup;
    int lctBXWindow;  // < 0: no BX selection
  };

  struct MuonMatch {
    unsigned int muonIndex;     // index in the muon collection given to match()
    const CSCSegment* segment;  // the muon's segment in the LCT chamber (ME1/1: ring 1 or ring 4); points into
                                // the event product reached through the muon (standalone-muon rechits or cscSegments)
    SegmentSource source;       // Standalone or Matches: where the segment was taken from
    csclct::KeyLayerPosition key;  // segment position at the key layer, LCT units
    bool keyFromRecHit;            // key from the segment's layer-3 rechit, else from the segment projection
    int deltaHalfStrip;            // LCT - key
    float deltaEighthStrip;        // LCT - key
    int deltaWireGroup;            // LCT - key
    float distance;                // |deltaEighthStrip| / 4 + |deltaWireGroup|, orders the matches
  };

  struct LCTMatch {
    CSCDetId chamber;                 // LCT chamber (layer 0; ME1/1 LCTs are in ring 1)
    const CSCCorrelatedLCTDigi* lct;  // points into the LCT collection given to match()
    std::vector<MuonMatch> muons;     // matched muons, closest first; empty if none
  };

  explicit CSCLCTRecoMuonMatcher(const Config& config);
  explicit CSCLCTRecoMuonMatcher(const edm::ParameterSet& pset);
  static void fillPSetDescription(edm::ParameterSetDescription& desc);

  // One entry per valid LCT (inside the BX window when one is set), in the order of the collection. The pointers in
  // the result refer to the LCT collection and to the segment products reached through the muons; they are valid as
  // long as those event products are.
  std::vector<LCTMatch> match(const CSCCorrelatedLCTDigiCollection& lcts,
                              const reco::MuonCollection& muons,
                              const CSCGeometry& geometry) const;

  // The CSC segments of one muon according to segmentSource, each with the source it was taken from.
  std::vector<std::pair<const CSCSegment*, SegmentSource>> muonSegments(const reco::Muon& muon) const;

  const Config& config() const { return config_; }

private:
  Config config_;
};

#endif
