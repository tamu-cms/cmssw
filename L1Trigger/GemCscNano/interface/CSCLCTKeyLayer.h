#ifndef L1Trigger_GemCscNano_CSCLCTKeyLayer_h
#define L1Trigger_GemCscNano_CSCLCTKeyLayer_h

// Conversion of offline CSC objects (rechits, segments) into LCT channel space at the key layer (layer 3):
//   key half-strip, 0-based (ME1/1: ME1b 0-127, ME1a 128-223),
//   eighth-strip in the CCLUT convention (es = 4*hs at the half-strip centre),
//   key wire group, 0-based.
// The ring of the offline object is always used for the geometry (a CLCT is built from one ring's comparators).

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include <optional>

class CSCGeometry;
class CSCLayer;
class CSCRecHit2D;
class CSCSegment;

namespace csclct {

  struct KeyLayerPosition {
    bool inStrips = false;  // strip coordinate inside (0, N) of the layer-3 geometry of the object's own ring
    int ring = 0;           // that ring (ME1/1: 1 = ME1b, 4 = ME1a)
    int hs = -1;            // key half-strip, LCT units
    float es = -1.f;        // eighth-strip, LCT units, CCLUT convention
    int wg = -1;            // key wire group, LCT units (0-based)
    bool valid() const { return inStrips && wg >= 0; }
  };

  bool isME11(const CSCDetId& id);

  // ME1/1 LCTs (ME1a and ME1b) are stored under ring 1, layer 0, in data and in emulation.
  CSCDetId lctChamberId(const CSCDetId& id);

  // ME1a half-strips follow the 128 ME1b half-strips in the ME1/1 LCT.
  int halfStripOffset(const CSCDetId& id);

  // lp in the frame of the layer-3 geometry L3.
  KeyLayerPosition toLCTSpace(const CSCLayer& L3, const LocalPoint& lp);

  // gp converted with the layer-3 geometry of the same ring as id (no ME1a <-> ME1b switching).
  KeyLayerPosition toLCTSpace(const CSCGeometry& geometry, const CSCDetId& id, const GlobalPoint& gp);

  // The segment's rechit in layer 3, nullptr if it has none.
  const CSCRecHit2D* keyLayerRecHit(const CSCSegment& segment);

  // Key-layer position measured by a layer-3 rechit; the wire group is the rechit's.
  KeyLayerPosition fromKeyLayerRecHit(const CSCGeometry& geometry, const CSCRecHit2D& rechit);

  // Where a straight line (chamber frame of chamberId: point p0, direction d) crosses the plane of layer 3 of that
  // chamber; none if the chamber is unknown or the line is parallel to the layer (|d.z| < 1e-6).
  std::optional<GlobalPoint> lineAtKeyLayer(const CSCGeometry& geometry,
                                            const CSCDetId& chamberId,
                                            const LocalPoint& p0,
                                            const LocalVector& d);

  // Segment (chamber frame, origin at chamber z = 0) projected onto the plane of layer 3.
  KeyLayerPosition fromSegmentProjection(const CSCGeometry& geometry, const CSCSegment& segment);

  // LCT - key, LCT units. distance = |dES| / 4 + |dWG| orders the matches.
  struct Residual {
    int dHS;
    float dES;
    int dWG;
    float distance;
  };
  Residual residual(int lctHalfStrip, int lctEighthStrip, int lctWireGroup, const KeyLayerPosition& key);

}  // namespace csclct

#endif
