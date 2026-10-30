#include "L1Trigger/CSCLCTRecoMuonMatching/interface/CSCLCTKeyLayer.h"

#include "DataFormats/CSCDigi/interface/CSCConstants.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"

#include <algorithm>
#include <cmath>

namespace csclct {

  bool isME11(const CSCDetId& id) { return id.station() == 1 && (id.ring() == 1 || id.ring() == 4); }

  CSCDetId lctChamberId(const CSCDetId& id) {
    return CSCDetId(id.endcap(), id.station(), isME11(id) ? 1 : id.ring(), id.chamber(), 0);
  }

  int halfStripOffset(const CSCDetId& id) {
    return (id.station() == 1 && id.ring() == 4) ? CSCConstants::NUM_HALF_STRIPS_ME1B : 0;
  }

  // Only the strip coordinate must be inside (0, N): near the ME1a/ME1b cut a point can lie a few cm outside its
  // ring's strip plane in y (RecoLocalMuon/CSCRecHitD/src/CSCRecHitDBuilder.cc:131-136) with a well-defined strip.
  KeyLayerPosition toLCTSpace(const CSCLayer& L3, const LocalPoint& lp) {
    const CSCLayerGeometry& g = *L3.geometry();
    KeyLayerPosition r;
    r.ring = L3.id().ring();
    const float fs = g.strip(lp);  // [0, N]; strip k covers [k-1, k); includes the layer-3 stagger
    r.inStrips = fs > 1.e-6f && fs < g.numberOfStrips() - 1.e-6f;
    const int o = halfStripOffset(L3.id());
    r.hs = std::clamp(static_cast<int>(std::floor(2.f * fs)), 0, 2 * g.numberOfStrips() - 1) + o;
    r.es = 8.f * fs - 2.f + 4.f * o;
    r.wg = g.wireGroup(g.nearestWire(lp)) - 1;  // wireGroup() is 1-based and 0 outside or dead
    return r;
  }

  KeyLayerPosition toLCTSpace(const CSCGeometry& geometry, const CSCDetId& id, const GlobalPoint& gp) {
    const CSCLayer* L3 =
        geometry.layer(CSCDetId(id.endcap(), id.station(), id.ring(), id.chamber(), CSCConstants::KEY_CLCT_LAYER));
    return L3 != nullptr ? toLCTSpace(*L3, L3->toLocal(gp)) : KeyLayerPosition();
  }

  const CSCRecHit2D* keyLayerRecHit(const CSCSegment& segment) {
    for (const auto& rh : segment.specificRecHits())
      if (rh.cscDetId().layer() == CSCConstants::KEY_CLCT_LAYER)
        return &rh;
    return nullptr;
  }

  KeyLayerPosition fromKeyLayerRecHit(const CSCGeometry& geometry, const CSCRecHit2D& rechit) {
    const CSCLayer* L = geometry.layer(rechit.cscDetId());
    if (L == nullptr)
      return KeyLayerPosition();
    KeyLayerPosition r = toLCTSpace(geometry, rechit.cscDetId(), L->toGlobal(rechit.localPosition()));
    r.wg = rechit.hitWire() - 1;  // middle wire group of the cluster, 1-based; ME1a uses the ME1/1 numbering
    return r;
  }

  KeyLayerPosition fromSegmentProjection(const CSCGeometry& geometry, const CSCSegment& segment) {
    const CSCChamber* ch = geometry.chamber(segment.cscDetId());
    if (ch == nullptr)
      return KeyLayerPosition();
    const CSCLayer* L3 = ch->layer(CSCConstants::KEY_CLCT_LAYER);
    const LocalVector d = segment.localDirection();
    if (L3 == nullptr || std::abs(d.z()) < 1e-6f)
      return KeyLayerPosition();
    const LocalPoint p0 = segment.localPosition();
    const float z3 = ch->toLocal(L3->toGlobal(LocalPoint(0.f, 0.f, 0.f))).z();
    const float t = (z3 - p0.z()) / d.z();
    const LocalPoint p3(p0.x() + t * d.x(), p0.y() + t * d.y(), z3);
    return toLCTSpace(geometry, segment.cscDetId(), ch->toGlobal(p3));
  }

  Residual residual(int lctHalfStrip, int lctEighthStrip, int lctWireGroup, const KeyLayerPosition& key) {
    Residual r;
    r.dHS = lctHalfStrip - key.hs;
    r.dES = float(lctEighthStrip) - key.es;
    r.dWG = lctWireGroup - key.wg;
    r.distance = std::abs(r.dES) / 4.f + std::abs(r.dWG);
    return r;
  }

}  // namespace csclct
