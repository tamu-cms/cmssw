#ifndef L1Trigger_GemCscNano_CSCLayerPoint_h
#define L1Trigger_GemCscNano_CSCLayerPoint_h

// A point on one CSC layer in the trigger's coordinates, unrounded (GemCscNano MuonCSCLayer, SimTrackCSCLayer):
//   x, y   in the layer frame (cm);
//   fs     CSCLayerGeometry::strip(lp): strip k covers [k-1, k), with the layer's stagger; CMSSW clamps it to [0, N]
//          (OffsetRadialStripTopology), so inStrips = 1e-6 < fs < N - 1e-6 tells a real value;
//   inStripPlane  y inside the layer's strip-plane y range (CSCLayerGeometry::yLimitsOfStripPlane);
//   wg     the nearest wire's group, 0-based (CSCLayerGeometry::wireGroup is 1-based and 0 outside or dead: -1).
// The LCT units follow from fs by the LCT numbering only: half-strip floor(2 fs) (+128 in ME1a), eighth-strip
// 8 fs - 2 (+512 in ME1a; 4 times the half-strip at its centre, CSCCorrelatedLCTDigi::getStrip(8)).

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

class CSCChamber;
class CSCLayer;

namespace csclct {

  struct LayerPoint {
    bool has = false;
    float x = 0.f, y = 0.f;
    float fs = -1.f;
    bool inStrips = false;
    bool inStripPlane = false;
    int wg = -1;
  };

  // lp in the frame of layer L (only x and y are used).
  LayerPoint layerPoint(const CSCLayer& L, const LocalPoint& lp);

  // Where a straight line (frame of chamber ch: p0 + t d) crosses the plane z = dz of layer L of ch (layer frame;
  // dz = 0 is the layer plane, the middle of its gas gap); has = false if |d.z| < 1e-6.
  LayerPoint lineAtLayer(const CSCChamber& ch, const CSCLayer& L, const LocalPoint& p0, const LocalVector& d, float dz);

  // Half the gas gap: CMSSW builds a CSC layer with the gas gap as its thickness (CSCGeometryBuilder.cc:204, 268).
  float halfGasGap(const CSCLayer& L);

}  // namespace csclct

#endif
