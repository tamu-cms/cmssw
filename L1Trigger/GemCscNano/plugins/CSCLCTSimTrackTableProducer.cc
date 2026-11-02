// NanoAOD tables of the SimTracks that reach the CSCs (MC):
//   SimTrack          the muon SimTracks and every SimTrack with a CSC simhit, in SimTrack collection order; genPartIdx is
//                     the row in the GenPart table (genParticles pruned by the genPartAssociations chain), parentIdx the row of the
//                     parent SimTrack;
//   SimTrackCSCLayer  per SimTrack, CSC chamber (with ring) and layer: its simhit and the straight line through its
//                     simhits in the chamber, in strip and wire-group coordinates. Only the simhits whose particle type
//                     is the track's own are used (a simhit of a secondary that was not saved carries its ancestor's
//                     track id; those of a secondary of the same particle type, e.g. a delta ray of an electron track,
//                     are kept).

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/CSCDigi/interface/CSCConstants.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "L1Trigger/GemCscNano/interface/CSCLayerPoint.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

  template <typename T, typename Row, typename F>
  void addColumn(nanoaod::FlatTable& table,
                 const std::string& name,
                 const std::vector<Row>& rows,
                 F value,
                 const std::string& doc) {
    std::vector<T> v;
    v.reserve(rows.size());
    for (const Row& r : rows)
      v.push_back(value(r));
    table.addColumn<T>(name, v, doc);
  }

  void checkIndexRange(size_t n, const char* what) {
    if (n > size_t(std::numeric_limits<int16_t>::max()))
      throw cms::Exception("LogicError") << n << " " << what << " in one event: the int16 index columns would overflow";
  }

  // Least-squares straight line x(z) = ax + bx z, y(z) = ay + by z through the simhit positions, chamber frame of ch;
  // false with fewer than two simhits or without two different z (normally: simhits on fewer than two layers).
  bool fitLine(const CSCGeometry& geometry,
               const CSCChamber& ch,
               const std::vector<const PSimHit*>& hits,
               double (&a)[4]) {
    double n = 0, sz = 0, szz = 0, sx = 0, sxz = 0, sy = 0, syz = 0;
    for (const PSimHit* h : hits) {
      const CSCLayer* L = geometry.layer(CSCDetId(h->detUnitId()));
      if (L == nullptr)
        continue;
      const LocalPoint p = ch.toLocal(L->toGlobal(h->localPosition()));
      n += 1;
      sz += p.z();
      szz += p.z() * p.z();
      sx += p.x();
      sxz += p.x() * p.z();
      sy += p.y();
      syz += p.y() * p.z();
    }
    const double det = n * szz - sz * sz;
    if (n < 2 || std::abs(det) < 1e-6)
      return false;
    const double bx = (n * sxz - sz * sx) / det, by = (n * syz - sz * sy) / det;
    a[0] = (sx - bx * sz) / n;
    a[1] = bx;
    a[2] = (sy - by * sz) / n;
    a[3] = by;
    return true;
  }

  // One layer of one (SimTrack, CSC chamber with ring where it has simhits of its own particle type).
  struct SimLayerRow {
    int16_t simTrackIdx;
    CSCDetId id;
    int16_t nSimHits = 0;          // own-type simhits on this layer
    const PSimHit* hit = nullptr;  // the first of them (collection order)
    csclct::LayerPoint simHit, entry, exit, line, lineLo, lineHi;
  };

}  // namespace

class CSCLCTSimTrackTableProducer : public edm::stream::EDProducer<> {
public:
  explicit CSCLCTSimTrackTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;
  const edm::EDGetTokenT<edm::SimVertexContainer> simVertexToken_;
  const edm::EDGetTokenT<edm::PSimHitContainer> simHitToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  const edm::EDGetTokenT<std::vector<int>> genBarcodeToken_;
  std::vector<edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>>> genAssociationTokens_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
};

CSCLCTSimTrackTableProducer::CSCLCTSimTrackTableProducer(const edm::ParameterSet& ps)
    : simTrackToken_(consumes(ps.getParameter<edm::InputTag>("simTracks"))),
      simVertexToken_(consumes(ps.getParameter<edm::InputTag>("simVertices"))),
      simHitToken_(consumes(ps.getParameter<edm::InputTag>("cscSimHits"))),
      genParticleToken_(consumes(ps.getParameter<edm::InputTag>("genParticles"))),
      genBarcodeToken_(consumes(ps.getParameter<edm::InputTag>("genParticles"))),
      cscGeomToken_(esConsumes()) {
  for (const auto& tag : ps.getParameter<std::vector<edm::InputTag>>("genPartAssociations"))
    genAssociationTokens_.push_back(consumes<edm::Association<reco::GenParticleCollection>>(tag));
  if (genAssociationTokens_.empty())
    throw cms::Exception("Configuration") << "genPartAssociations is empty";
  produces<nanoaod::FlatTable>("simTrack");
  produces<nanoaod::FlatTable>("simTrackCSCLayer");
}

void CSCLCTSimTrackTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("simVertices", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("cscSimHits", edm::InputTag("g4SimHits", "MuonCSCHits"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"))
      ->setComment("reco::GenParticle collection and its HepMC barcode vector (same label)");
  desc.add<std::vector<edm::InputTag>>("genPartAssociations",
                                       {edm::InputTag("prunedGenParticlesWithStatusOne"),
                                        edm::InputTag("prunedGenParticles"),
                                        edm::InputTag("finalGenParticles")})
      ->setComment(
          "GenParticlePruner associations from genParticles to the collection written as GenPart, applied in order "
          "(default: the central miniAOD and NanoAOD chain)");
  descriptions.add("cscLCTSimTrackTableProducer", desc);
}

void CSCLCTSimTrackTableProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  const auto& simTracks = event.get(simTrackToken_);
  const auto& simVertices = event.get(simVertexToken_);
  const auto& simHits = event.get(simHitToken_);
  const auto genParticles = event.getHandle(genParticleToken_);
  const auto& barcodes = event.get(genBarcodeToken_);
  std::vector<const edm::Association<reco::GenParticleCollection>*> genAssociations;
  genAssociations.reserve(genAssociationTokens_.size());
  for (const auto& token : genAssociationTokens_)
    genAssociations.push_back(&event.get(token));
  const CSCGeometry& geometry = setup.getData(cscGeomToken_);
  if (barcodes.size() != genParticles->size())
    throw cms::Exception("ProductNotValid")
        << "genParticles: " << genParticles->size() << " particles but " << barcodes.size() << " barcodes";

  std::unordered_map<unsigned int, std::vector<const PSimHit*>> hitsOfTrack;
  for (const PSimHit& h : simHits)
    hitsOfTrack[h.trackId()].push_back(&h);

  std::vector<size_t> rows;  // SimTrack collection indices of the table rows
  std::unordered_map<unsigned int, int16_t> rowOfTrackId;
  for (size_t i = 0; i < simTracks.size(); ++i)
    if (std::abs(simTracks[i].type()) == 13 || hitsOfTrack.count(simTracks[i].trackId())) {
      rowOfTrackId[simTracks[i].trackId()] = rows.size();
      rows.push_back(i);
    }
  checkIndexRange(rows.size(), "SimTracks");

  std::unordered_map<int, size_t> genIndexOfBarcode;
  for (size_t g = 0; g < barcodes.size(); ++g)
    genIndexOfBarcode[barcodes[g]] = g;
  const auto genPartIdx = [&](size_t i) -> int16_t {
    const SimTrack& t = simTracks[i];
    if (t.noGenpart())
      return -1;
    const auto g = genIndexOfBarcode.find(t.genpartIndex());
    if (g == genIndexOfBarcode.end())
      throw cms::Exception("ProductNotValid")
          << "SimTrack " << t.trackId() << ": barcode " << t.genpartIndex() << " not in genParticles";
    reco::GenParticleRef kept(genParticles, g->second);
    for (const auto* association : genAssociations) {
      kept = (*association)[kept];
      if (kept.isNull())
        return -1;
    }
    checkIndexRange(kept.key() + 1, "GenPart rows");
    return kept.key();
  };
  const auto parentTrackId = [&](size_t i) -> int {
    const SimTrack& t = simTracks[i];
    return t.noVertex() ? -1 : simVertices.at(t.vertIndex()).parentIndex();
  };

  // SimTrackCSCLayer: per stored SimTrack (row order), every CSC chamber (with ring) where it has simhits of its own
  // particle type (CSCDetId order), layers 1-6.
  std::vector<SimLayerRow> simLayers;
  for (size_t s = 0; s < rows.size(); ++s) {
    const SimTrack& t = simTracks[rows[s]];
    const auto trackHits = hitsOfTrack.find(t.trackId());
    if (trackHits == hitsOfTrack.end())
      continue;
    std::map<CSCDetId, std::vector<const PSimHit*>> hitsOfChamber;  // own ring (ME1a: 4)
    for (const PSimHit* h : trackHits->second)
      if (h->particleType() == t.type())
        hitsOfChamber[CSCDetId(h->detUnitId()).chamberId()].push_back(h);
    for (const auto& [chamberId, hits] : hitsOfChamber) {
      const CSCChamber* ch = geometry.chamber(chamberId);
      if (ch == nullptr)
        throw cms::Exception("LogicError") << chamberId << " is not in the CSC geometry";
      double a[4];
      const bool hasLine = fitLine(geometry, *ch, hits, a);
      for (int l = 1; l <= CSCConstants::NUM_LAYERS; ++l) {
        const CSCLayer* L = ch->layer(l);
        if (L == nullptr)
          throw cms::Exception("LogicError") << "layer " << l << " of " << chamberId << " is not in the CSC geometry";
        SimLayerRow row{int16_t(s), L->id(), 0, nullptr, {}, {}, {}, {}, {}, {}};
        for (const PSimHit* h : hits)
          if (CSCDetId(h->detUnitId()).layer() == l) {
            if (row.hit == nullptr) {
              row.hit = h;
              row.simHit = csclct::layerPoint(*L, h->localPosition());
              row.entry = csclct::layerPoint(*L, h->entryPoint());
              row.exit = csclct::layerPoint(*L, h->exitPoint());
            }
            ++row.nSimHits;
          }
        if (hasLine) {
          const float halfGap = csclct::halfGasGap(*L);
          const auto at = [&](float dz) {
            const double z = ch->toLocal(L->toGlobal(LocalPoint(0.f, 0.f, dz))).z();
            return csclct::layerPoint(*L, L->toLocal(ch->toGlobal(LocalPoint(a[0] + a[1] * z, a[2] + a[3] * z, z))));
          };
          row.line = at(0.f);
          row.lineLo = at(-halfGap);
          row.lineHi = at(halfGap);
        }
        simLayers.push_back(row);
      }
    }
  }
  auto simLayerTable = std::make_unique<nanoaod::FlatTable>(simLayers.size(), "SimTrackCSCLayer", false);
  simLayerTable->setDoc(
      "MC: per stored SimTrack (row order), every CSC chamber (with ring) where it has simhits of its own particle "
      "type (CSCDetId order), layers 1-6. Two sources of a position on the layer: the first own-type simhit on it "
      "(simHit*: its entry-exit midpoint, entry and exit points) and the least-squares straight line through the "
      "own-type simhit midpoints of the chamber, chamber frame (line*; none without own-type simhits at two different "
      "z, normally two layers). Per point as "
      "in MuonCSCLayer: x, y in the layer frame, fs = CSCLayerGeometry::strip (clamped to [0, N]), wg 0-based; "
      "*Lo / *Hi at the gas-gap surfaces, layer-frame z = -g/2 and +g/2");
  addColumn<int16_t>(
      *simLayerTable, "simTrackIdx", simLayers, [](const SimLayerRow& r) { return r.simTrackIdx; }, "row in SimTrack");
  addColumn<int16_t>(
      *simLayerTable, "endcap", simLayers, [](const SimLayerRow& r) { return r.id.zendcap(); }, "+1 or -1");
  addColumn<int16_t>(
      *simLayerTable, "station", simLayers, [](const SimLayerRow& r) { return r.id.station(); }, "station");
  addColumn<int16_t>(
      *simLayerTable,
      "ring",
      simLayers,
      [](const SimLayerRow& r) { return r.id.ring(); },
      "ring of the simhits (ME1/1: 1 = ME1b, 4 = ME1a)");
  addColumn<int16_t>(
      *simLayerTable, "chamber", simLayers, [](const SimLayerRow& r) { return r.id.chamber(); }, "chamber");
  addColumn<int16_t>(
      *simLayerTable, "layer", simLayers, [](const SimLayerRow& r) { return r.id.layer(); }, "layer, 1-6");
  addColumn<int16_t>(
      *simLayerTable,
      "nSimHits",
      simLayers,
      [](const SimLayerRow& r) { return r.nSimHits; },
      "simhits of the track's own particle type on this layer");
  addColumn<bool>(
      *simLayerTable,
      "hasSimHit",
      simLayers,
      [](const SimLayerRow& r) { return r.simHit.has; },
      "an own-type simhit lies on this layer (the first is described)");
  addColumn<float>(
      *simLayerTable,
      "simHitX",
      simLayers,
      [](const SimLayerRow& r) { return r.simHit.x; },
      "simhit midpoint (PSimHit::localPosition) x, layer frame, cm");
  addColumn<float>(
      *simLayerTable,
      "simHitY",
      simLayers,
      [](const SimLayerRow& r) { return r.simHit.y; },
      "simhit midpoint y, layer frame, cm");
  addColumn<float>(
      *simLayerTable, "simHitFs", simLayers, [](const SimLayerRow& r) { return r.simHit.fs; }, "fs of the midpoint");
  addColumn<bool>(
      *simLayerTable,
      "simHitInStrips",
      simLayers,
      [](const SimLayerRow& r) { return r.simHit.inStrips; },
      "1e-6 < simHitFs < N - 1e-6 (not clamped)");
  addColumn<bool>(
      *simLayerTable,
      "simHitInStripPlane",
      simLayers,
      [](const SimLayerRow& r) { return r.simHit.inStripPlane; },
      "simHitY inside the strip plane");
  addColumn<int16_t>(
      *simLayerTable,
      "simHitWg",
      simLayers,
      [](const SimLayerRow& r) { return r.simHit.wg; },
      "wire group of the nearest wire at the midpoint, 0-based; -1 outside or dead");
  addColumn<float>(
      *simLayerTable,
      "simHitFsEntry",
      simLayers,
      [](const SimLayerRow& r) { return r.entry.fs; },
      "fs of the simhit entry point (clamped like fs)");
  addColumn<int16_t>(
      *simLayerTable, "simHitWgEntry", simLayers, [](const SimLayerRow& r) { return r.entry.wg; }, "wire group there");
  addColumn<float>(
      *simLayerTable,
      "simHitFsExit",
      simLayers,
      [](const SimLayerRow& r) { return r.exit.fs; },
      "fs of the simhit exit point (clamped like fs)");
  addColumn<int16_t>(
      *simLayerTable, "simHitWgExit", simLayers, [](const SimLayerRow& r) { return r.exit.wg; }, "wire group there");
  addColumn<float>(
      *simLayerTable,
      "simHitTof",
      simLayers,
      [](const SimLayerRow& r) { return r.hit != nullptr ? r.hit->timeOfFlight() : 0.f; },
      "simhit time of flight, ns");
  addColumn<float>(
      *simLayerTable,
      "simHitEnergyLoss",
      simLayers,
      [](const SimLayerRow& r) { return r.hit != nullptr ? r.hit->energyLoss() : 0.f; },
      "simhit energy loss, GeV");
  addColumn<float>(
      *simLayerTable,
      "simHitPabs",
      simLayers,
      [](const SimLayerRow& r) { return r.hit != nullptr ? r.hit->pabs() : 0.f; },
      "momentum at the simhit entry, GeV");
  addColumn<int16_t>(
      *simLayerTable,
      "simHitProcessType",
      simLayers,
      [](const SimLayerRow& r) { return r.hit != nullptr ? int(r.hit->processType()) : 0; },
      "Geant4 process type of the simhit");
  addColumn<bool>(
      *simLayerTable,
      "hasLine",
      simLayers,
      [](const SimLayerRow& r) { return r.line.has; },
      "the least-squares line exists (own-type simhits at two or more different chamber z, normally on two or more "
      "layers)");
  addColumn<float>(
      *simLayerTable,
      "lineX",
      simLayers,
      [](const SimLayerRow& r) { return r.line.x; },
      "x where the line crosses the layer plane, layer frame, cm");
  addColumn<float>(
      *simLayerTable,
      "lineY",
      simLayers,
      [](const SimLayerRow& r) { return r.line.y; },
      "y where the line crosses the layer plane, layer frame, cm");
  addColumn<float>(
      *simLayerTable, "lineFs", simLayers, [](const SimLayerRow& r) { return r.line.fs; }, "fs at that point");
  addColumn<bool>(
      *simLayerTable,
      "lineInStrips",
      simLayers,
      [](const SimLayerRow& r) { return r.line.inStrips; },
      "1e-6 < lineFs < N - 1e-6 (not clamped)");
  addColumn<bool>(
      *simLayerTable,
      "lineInStripPlane",
      simLayers,
      [](const SimLayerRow& r) { return r.line.inStripPlane; },
      "lineY inside the strip plane");
  addColumn<int16_t>(
      *simLayerTable,
      "lineWg",
      simLayers,
      [](const SimLayerRow& r) { return r.line.wg; },
      "wire group of the nearest wire at that point, 0-based; -1 outside or dead");
  addColumn<float>(
      *simLayerTable,
      "lineFsLo",
      simLayers,
      [](const SimLayerRow& r) { return r.lineLo.fs; },
      "fs where the line crosses the gas-gap surface at layer-frame z = -g/2 (clamped like fs)");
  addColumn<int16_t>(
      *simLayerTable, "lineWgLo", simLayers, [](const SimLayerRow& r) { return r.lineLo.wg; }, "wire group there");
  addColumn<float>(
      *simLayerTable,
      "lineFsHi",
      simLayers,
      [](const SimLayerRow& r) { return r.lineHi.fs; },
      "fs where the line crosses the gas-gap surface at layer-frame z = +g/2 (clamped like fs)");
  addColumn<int16_t>(
      *simLayerTable, "lineWgHi", simLayers, [](const SimLayerRow& r) { return r.lineHi.wg; }, "wire group there");
  event.put(std::move(simLayerTable), "simTrackCSCLayer");

  std::vector<size_t> rowNumbers(rows.size());
  for (size_t s = 0; s < rows.size(); ++s)
    rowNumbers[s] = s;
  const auto track = [&](size_t s) -> const SimTrack& { return simTracks[rows[s]]; };
  auto simTrackTable = std::make_unique<nanoaod::FlatTable>(rows.size(), "SimTrack", false);
  simTrackTable->setDoc("muon SimTracks and SimTracks with CSC simhits (g4SimHits, signal only), in collection order");
  addColumn<int32_t>(
      *simTrackTable, "trackId", rowNumbers, [&](size_t s) { return track(s).trackId(); }, "Geant4 track id");
  addColumn<int32_t>(*simTrackTable, "pdgId", rowNumbers, [&](size_t s) { return track(s).type(); }, "PDG id");
  addColumn<int16_t>(
      *simTrackTable, "charge", rowNumbers, [&](size_t s) { return std::lround(track(s).charge()); }, "electric charge");
  addColumn<float>(
      *simTrackTable, "pt", rowNumbers, [&](size_t s) { return track(s).momentum().pt(); }, "pt at production");
  addColumn<float>(
      *simTrackTable, "eta", rowNumbers, [&](size_t s) { return track(s).momentum().eta(); }, "eta at production");
  addColumn<float>(
      *simTrackTable, "phi", rowNumbers, [&](size_t s) { return track(s).momentum().phi(); }, "phi at production");
  const auto vertexCoordinate = [&](size_t s, int c) -> float {
    const SimTrack& t = track(s);
    if (t.noVertex())
      return 0.f;
    const auto& p = simVertices.at(t.vertIndex()).position();
    return c == 0 ? p.x() : (c == 1 ? p.y() : p.z());
  };
  addColumn<float>(
      *simTrackTable,
      "vx",
      rowNumbers,
      [&](size_t s) { return vertexCoordinate(s, 0); },
      "production vertex x, cm (0 if none)");
  addColumn<float>(
      *simTrackTable,
      "vy",
      rowNumbers,
      [&](size_t s) { return vertexCoordinate(s, 1); },
      "production vertex y, cm (0 if none)");
  addColumn<float>(
      *simTrackTable,
      "vz",
      rowNumbers,
      [&](size_t s) { return vertexCoordinate(s, 2); },
      "production vertex z, cm (0 if none)");
  addColumn<int16_t>(
      *simTrackTable,
      "genPartIdx",
      rowNumbers,
      [&](size_t s) { return genPartIdx(rows[s]); },
      "row in GenPart of the generator particle (primary tracks), -1 if none or not stored");
  addColumn<int32_t>(
      *simTrackTable,
      "parentTrackId",
      rowNumbers,
      [&](size_t s) { return parentTrackId(rows[s]); },
      "Geant4 track id of the parent track (SimVertex::parentIndex), -1 if none");
  addColumn<int16_t>(
      *simTrackTable,
      "parentIdx",
      rowNumbers,
      [&](size_t s) {
        const int p = parentTrackId(rows[s]);
        const auto r = p < 0 ? rowOfTrackId.end() : rowOfTrackId.find(p);
        return r == rowOfTrackId.end() ? int16_t(-1) : r->second;
      },
      "row in SimTrack of the parent track, -1 if none or not stored");
  addColumn<int16_t>(
      *simTrackTable,
      "nCSCSimHits",
      rowNumbers,
      [&](size_t s) {
        const auto h = hitsOfTrack.find(track(s).trackId());
        return h == hitsOfTrack.end() ? 0 : h->second.size();
      },
      "number of CSC simhits with this track id (including those of unsaved secondaries)");

  event.put(std::move(simTrackTable), "simTrack");
}

DEFINE_FWK_MODULE(CSCLCTSimTrackTableProducer);
