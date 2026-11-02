#include "L1Trigger/GemCscNano/interface/CSCLayerPoint.h"

#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"

#include <cmath>

namespace csclct {

  LayerPoint layerPoint(const CSCLayer& L, const LocalPoint& lp) {
    const CSCLayerGeometry& g = *L.geometry();
    LayerPoint r;
    r.has = true;
    r.x = lp.x();
    r.y = lp.y();
    r.fs = g.strip(lp);
    r.inStrips = r.fs > 1.e-6f && r.fs < g.numberOfStrips() - 1.e-6f;
    const auto limits = g.yLimitsOfStripPlane();
    r.inStripPlane = limits.first <= r.y && r.y <= limits.second;
    r.wg = g.wireGroup(g.nearestWire(lp)) - 1;
    return r;
  }

  LayerPoint lineAtLayer(const CSCChamber& ch, const CSCLayer& L, const LocalPoint& p0, const LocalVector& d, float dz) {
    if (std::abs(d.z()) < 1e-6f)
      return LayerPoint();
    const float zL = ch.toLocal(L.toGlobal(LocalPoint(0.f, 0.f, dz))).z();
    const float t = (zL - p0.z()) / d.z();
    const LocalPoint p(p0.x() + t * d.x(), p0.y() + t * d.y(), zL);
    return layerPoint(L, L.toLocal(ch.toGlobal(p)));
  }

  float halfGasGap(const CSCLayer& L) { return 0.5f * L.surface().bounds().thickness(); }

}  // namespace csclct
