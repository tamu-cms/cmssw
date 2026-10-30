// NanoAOD tables of the SimTracks that reach the CSCs and of their association with the rows of the LCT table (MC):
//   SimTrack          the muon SimTracks and every SimTrack with a CSC simhit, in SimTrack collection order; genPartIdx is
//                     the row in the GenPart table (genParticles pruned by genPartAssociation), parentIdx the row of the
//                     parent SimTrack;
//   LCTSimTrackMatch  one row per LCT-SimTrack pair within the windows, grouped by LCT, closest SimTrack first, with the
//                     residuals LCT - simhit position at the key layer (layer 3), LCT units;
//   LCT (extension)   nSimTrackMatch.
// The key-layer position of a SimTrack in a chamber is its own layer-3 simhit (entry-exit midpoint), otherwise the
// straight line through its simhits in that chamber taken at layer 3. Only the simhits whose particle type is the
// track's own are used for the position (a simhit of a secondary that was not saved carries its ancestor's track id).
// As for segments, ME1a (ring 4) and ME1b (ring 1) simhits both map to the ring-1 LCT chamber, each with the geometry of
// its own ring, and per LCT a SimTrack keeps its closest position.
// The LCTs are read from the LCT table, so the lctIdx are its rows by construction.

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
#include "L1Trigger/CSCLCTRecoMuonMatching/interface/CSCLCTKeyLayer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include <algorithm>
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

  // Key-layer position of a track from its simhits in one chamber (rows of the track's own particle type).
  csclct::KeyLayerPosition keyPosition(const CSCGeometry& geometry,
                                       const CSCDetId& chamberId,
                                       const std::vector<const PSimHit*>& hits,
                                       bool& fromLayer3) {
    fromLayer3 = false;
    const CSCChamber* ch = geometry.chamber(chamberId);
    if (ch == nullptr)
      return csclct::KeyLayerPosition();
    const CSCLayer* L3 = ch->layer(CSCConstants::KEY_CLCT_LAYER);
    if (L3 == nullptr)
      return csclct::KeyLayerPosition();
    for (const PSimHit* h : hits)
      if (CSCDetId(h->detUnitId()).layer() == CSCConstants::KEY_CLCT_LAYER) {
        fromLayer3 = true;
        return csclct::toLCTSpace(geometry, chamberId, L3->toGlobal(h->localPosition()));
      }
    // Straight-line least-squares fit x(z), y(z) in the chamber frame, taken at the layer-3 plane.
    double n = 0, sz = 0, szz = 0, sx = 0, sxz = 0, sy = 0, syz = 0;
    for (const PSimHit* h : hits) {
      const CSCLayer* L = geometry.layer(CSCDetId(h->detUnitId()));
      if (L == nullptr)
        continue;
      const LocalPoint p = ch->toLocal(L->toGlobal(h->localPosition()));
      n += 1;
      sz += p.z();
      szz += p.z() * p.z();
      sx += p.x();
      sxz += p.x() * p.z();
      sy += p.y();
      syz += p.y() * p.z();
    }
    const double det = n * szz - sz * sz;
    if (n < 2 || std::abs(det) < 1e-6)  // fewer than two layers
      return csclct::KeyLayerPosition();
    const double bx = (n * sxz - sz * sx) / det, by = (n * syz - sz * sy) / det;
    const double ax = (sx - bx * sz) / n, ay = (sy - by * sz) / n;
    const double z3 = ch->toLocal(L3->toGlobal(LocalPoint(0.f, 0.f, 0.f))).z();
    return csclct::toLCTSpace(geometry, chamberId, ch->toGlobal(LocalPoint(ax + bx * z3, ay + by * z3, z3)));
  }

  struct LinkRow {
    int16_t lctIdx;
    int16_t simTrackIdx;
    bool fromLayer3;
    int16_t nSimHits;
    csclct::KeyLayerPosition key;
    csclct::Residual residual;
  };

}  // namespace

class CSCLCTSimTrackTableProducer : public edm::stream::EDProducer<> {
public:
  explicit CSCLCTSimTrackTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<nanoaod::FlatTable> lctTableToken_;
  const edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;
  const edm::EDGetTokenT<edm::SimVertexContainer> simVertexToken_;
  const edm::EDGetTokenT<edm::PSimHitContainer> simHitToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  const edm::EDGetTokenT<std::vector<int>> genBarcodeToken_;
  const edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> genAssociationToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const int maxDeltaHalfStrip_;
  const int maxDeltaWireGroup_;
};

CSCLCTSimTrackTableProducer::CSCLCTSimTrackTableProducer(const edm::ParameterSet& ps)
    : lctTableToken_(consumes(ps.getParameter<edm::InputTag>("lctTable"))),
      simTrackToken_(consumes(ps.getParameter<edm::InputTag>("simTracks"))),
      simVertexToken_(consumes(ps.getParameter<edm::InputTag>("simVertices"))),
      simHitToken_(consumes(ps.getParameter<edm::InputTag>("cscSimHits"))),
      genParticleToken_(consumes(ps.getParameter<edm::InputTag>("genParticles"))),
      genBarcodeToken_(consumes(ps.getParameter<edm::InputTag>("genParticles"))),
      genAssociationToken_(consumes(ps.getParameter<edm::InputTag>("genPartAssociation"))),
      cscGeomToken_(esConsumes()),
      maxDeltaHalfStrip_(ps.getParameter<int>("maxDeltaHalfStrip")),
      maxDeltaWireGroup_(ps.getParameter<int>("maxDeltaWireGroup")) {
  if (maxDeltaHalfStrip_ < 0 || maxDeltaWireGroup_ < 0)
    throw cms::Exception("Configuration") << "CSCLCTSimTrackTableProducer: negative match window (maxDeltaHalfStrip "
                                          << maxDeltaHalfStrip_ << ", maxDeltaWireGroup " << maxDeltaWireGroup_ << ")";
  produces<nanoaod::FlatTable>("simTrack");
  produces<nanoaod::FlatTable>("lctSimTrackMatch");
  produces<nanoaod::FlatTable>("lctExtension");
}

void CSCLCTSimTrackTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("lctTable", edm::InputTag("cscLCTRecoMuonTable", "lct"))
      ->setComment("the LCT FlatTable of CSCLCTRecoMuonTableProducer");
  desc.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("simVertices", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("cscSimHits", edm::InputTag("g4SimHits", "MuonCSCHits"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"))
      ->setComment("reco::GenParticle collection and its HepMC barcode vector (same label)");
  desc.add<edm::InputTag>("genPartAssociation", edm::InputTag("cscLCTGenParticles"))
      ->setComment("GenParticlePruner association genParticles -> the collection written as GenPart");
  desc.add<int>("maxDeltaHalfStrip", 6)->setComment("match window |LCT - simhit| in key half-strips");
  desc.add<int>("maxDeltaWireGroup", 2)->setComment("match window |LCT - simhit| in key wire groups");
  descriptions.add("cscLCTSimTrackTableProducer", desc);
}

void CSCLCTSimTrackTableProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  const auto& lctTable = event.get(lctTableToken_);
  const auto& simTracks = event.get(simTrackToken_);
  const auto& simVertices = event.get(simVertexToken_);
  const auto& simHits = event.get(simHitToken_);
  const auto genParticles = event.getHandle(genParticleToken_);
  const auto& barcodes = event.get(genBarcodeToken_);
  const auto& genAssociation = event.get(genAssociationToken_);
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
  checkIndexRange(lctTable.nRows(), "LCTs");

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
    const reco::GenParticleRef kept = genAssociation[reco::GenParticleRef(genParticles, g->second)];
    return kept.isNonnull() ? int16_t(kept.key()) : int16_t(-1);
  };
  const auto parentTrackId = [&](size_t i) -> int {
    const SimTrack& t = simTracks[i];
    return t.noVertex() ? -1 : simVertices.at(t.vertIndex()).parentIndex();
  };

  // LCT rows per LCT chamber, from the LCT table.
  const auto column = [&lctTable](const char* name) {
    const int c = lctTable.columnIndex(name);
    if (c < 0)
      throw cms::Exception("ProductNotValid") << "LCT table has no column " << name;
    return lctTable.columnData<int16_t>(c);
  };
  const auto endcap = column("endcap"), station = column("station"), ring = column("ring"), chamber = column("chamber"),
             hs = column("hs"), es = column("es"), wg = column("wg");
  std::map<CSCDetId, std::vector<int16_t>> lctRowsOfChamber;
  for (unsigned int r = 0; r < lctTable.nRows(); ++r)
    lctRowsOfChamber[CSCDetId(endcap[r] > 0 ? 1 : 2, station[r], ring[r], chamber[r], 0)].push_back(r);

  std::vector<LinkRow> links;
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
      const auto lcts = lctRowsOfChamber.find(csclct::lctChamberId(chamberId));
      if (lcts == lctRowsOfChamber.end())
        continue;
      bool fromLayer3 = false;
      const csclct::KeyLayerPosition key = keyPosition(geometry, chamberId, hits, fromLayer3);
      if (!key.valid())
        continue;
      for (const int16_t r : lcts->second) {
        const csclct::Residual res = csclct::residual(hs[r], es[r], wg[r], key);
        if (std::abs(res.dHS) > maxDeltaHalfStrip_ || std::abs(res.dWG) > maxDeltaWireGroup_)
          continue;
        links.push_back({r, int16_t(s), fromLayer3, int16_t(hits.size()), key, res});
      }
    }
  }
  // Grouped by LCT, closest first; per LCT a SimTrack keeps its closest position (ME1a and ME1b of one track).
  std::stable_sort(links.begin(), links.end(), [](const LinkRow& a, const LinkRow& b) {
    return a.lctIdx != b.lctIdx ? a.lctIdx < b.lctIdx : a.residual.distance < b.residual.distance;
  });
  std::vector<LinkRow> unique;
  for (const LinkRow& l : links)
    if (std::none_of(unique.begin(), unique.end(), [&l](const LinkRow& u) {
          return u.lctIdx == l.lctIdx && u.simTrackIdx == l.simTrackIdx;
        }))
      unique.push_back(l);
  links = std::move(unique);
  std::vector<int16_t> lctNMatch(lctTable.nRows(), 0), simTrackNMatch(rows.size(), 0);
  for (const LinkRow& l : links) {
    ++lctNMatch[l.lctIdx];
    ++simTrackNMatch[l.simTrackIdx];
  }

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
  addColumn<int16_t>(
      *simTrackTable,
      "nLCTMatch",
      rowNumbers,
      [&](size_t s) { return simTrackNMatch[s]; },
      "number of LCTSimTrackMatch rows of this track");

  auto linkTable = std::make_unique<nanoaod::FlatTable>(links.size(), "LCTSimTrackMatch", false);
  linkTable->setDoc(
      "LCT-SimTrack pairs within the windows (many-to-many), grouped by LCT, closest SimTrack first; residuals are "
      "LCT - simhit position at the key layer (layer 3), LCT units");
  addColumn<int16_t>(*linkTable, "lctIdx", links, [](const LinkRow& l) { return l.lctIdx; }, "row in LCT");
  addColumn<int16_t>(
      *linkTable, "simTrackIdx", links, [](const LinkRow& l) { return l.simTrackIdx; }, "row in SimTrack");
  addColumn<bool>(
      *linkTable,
      "keyFromLayer3",
      links,
      [](const LinkRow& l) { return l.fromLayer3; },
      "key position from the track's layer-3 simhit, else from the line through its simhits in the chamber");
  addColumn<int16_t>(
      *linkTable,
      "nSimHits",
      links,
      [](const LinkRow& l) { return l.nSimHits; },
      "simhits of the track (own particle type) in the chamber (ME1/1: in its ring)");
  addColumn<uint8_t>(
      *linkTable, "keyRing", links, [](const LinkRow& l) { return l.key.ring; }, "ring of the simhits (ME1a: 4)");
  addColumn<int16_t>(*linkTable, "keyHS", links, [](const LinkRow& l) { return l.key.hs; }, "simhit key half-strip");
  addColumn<float>(
      *linkTable, "keyES", links, [](const LinkRow& l) { return l.key.es; }, "simhit eighth-strip (CCLUT)");
  addColumn<int16_t>(*linkTable, "keyWG", links, [](const LinkRow& l) { return l.key.wg; }, "simhit key wire group");
  addColumn<int16_t>(
      *linkTable, "dHS", links, [](const LinkRow& l) { return l.residual.dHS; }, "LCT - simhit half-strip");
  addColumn<float>(
      *linkTable, "dES", links, [](const LinkRow& l) { return l.residual.dES; }, "LCT - simhit eighth-strip");
  addColumn<int16_t>(
      *linkTable, "dWG", links, [](const LinkRow& l) { return l.residual.dWG; }, "LCT - simhit wire group");

  auto lctExtension = std::make_unique<nanoaod::FlatTable>(lctTable.nRows(), "LCT", false, true);
  lctExtension->addColumn<int16_t>("nSimTrackMatch", lctNMatch, "number of LCTSimTrackMatch rows of this LCT");

  event.put(std::move(simTrackTable), "simTrack");
  event.put(std::move(linkTable), "lctSimTrackMatch");
  event.put(std::move(lctExtension), "lctExtension");
}

DEFINE_FWK_MODULE(CSCLCTSimTrackTableProducer);
