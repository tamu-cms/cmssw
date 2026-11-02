// NanoAOD tables of CSC LCTs, reco muons and their CSCLCTRecoMuonMatcher associations (one entry per event):
//   LCT           every LCT the matcher considers (valid, inside lctBXWindow), in collection order;
//   Muon          the reco muons passing muonCut (default pT > 10 GeV and the medium or tight cut-based ID), in
//                 collection order; only these muons are matched;
//   LCTMuonMatch  one row per matched LCT-muon pair (rows of the LCT and Muon tables, segment, key-layer residuals),
//                 grouped by LCT, closest muon first;
//   MuonCSCChamber one row per CSC chamber match (reco::MuonChamberMatch) of each stored muon: the track propagated
//                 by reco to the chamber, that position carried to the key layer in LCT units, and which segments
//                 the muon has in the chamber;
//   MuonGEMChamber one row per GEM chamber match of each stored muon: the propagated track on the chamber surface,
//                 the eta partition, strip and pad where its line crosses the partition plane; MuonGEMHitMatch: the
//                 GEM rechits reco matched there.
// No event selection is applied: tag and probe is left to the analysis of the ntuple.

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/CSCDigi/interface/CSCConstants.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonGEMHitMatch.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
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
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "L1Trigger/GemCscNano/interface/CSCLCTRecoMuonMatcher.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace {

  // One column computed from each row.
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

  // One CSC chamber match of a stored muon, with the propagated position carried to the key layer.
  struct ChamberRow {
    int16_t muonIdx;
    const reco::MuonChamberMatch* match;
    CSCDetId id;                        // the chamber of the match (ME1/1: ring 1 = ME1b or ring 4 = ME1a)
    csclct::KeyLayerPosition key;       // propagated position on layer 3, LCT units, with the geometry of this ring
    bool inStripPlane = false;          // that point lies inside the y range of this ring's layer-3 strip plane
    float eta = 0.f, phi = 0.f;         // of that point (0 if none)
    bool hasMatcherSegment = false;     // the matcher's segmentSource gives a segment in this chamber
    bool hasStandaloneSegment = false;  // the standalone track has a CSC segment in this chamber
    bool hasArbitratedSegment = false;  // the arbitrated chamber matches have a segment here (segmentSource matches)
  };

  // One CSC segment of a stored muon: a hit of its standalone track (source 0) or a segment match of one of its
  // CSC chamber matches (source 1).
  struct SegmentRow {
    int16_t muonIdx;
    uint8_t source;
    int16_t chamberIdx;  // MuonCSCChamber row holding the segment match (source 1), -1 for source 0
    uint32_t mask;       // reco::MuonSegmentMatch::mask (source 1), 0 for source 0
    const CSCSegment* segment;
  };

  // One layer of one (muon, CSC chamber with ring, segment or none): the three sources of a position on that layer.
  struct LayerRow {
    int16_t muonIdx;
    int16_t segmentIdx;  // MuonCSCSegment row, -1 if none
    int16_t chamberIdx;  // MuonCSCChamber row of the same chamber (ring), -1 if none
    CSCDetId id;         // the layer
    const CSCRecHit2D* rechit = nullptr;
    csclct::LayerPoint recHit, line, lineLo, lineHi, prop, propLo, propHi;
  };

  struct ADCRow {
    int16_t layerIdx;
    int16_t strip;
    int16_t timeBin;
    float adc;
  };

  // One GEM chamber match of a stored muon, with the eta partition, strip and pad where the track crosses it.
  struct GEMChamberRow {
    int16_t muonIdx;
    const reco::MuonChamberMatch* match;
    GEMDetId id;                 // the chamber of the match (layer 0: a super-chamber)
    int16_t roll = 0;            // eta partition whose bounds hold the crossing point, 0 if none
    float strip = -1.f;          // strip coordinate of the crossing point in that partition, -1 if none
    float pad = -1.f;            // pad coordinate of the crossing point in that partition, -1 if none
    float eta = 0.f, phi = 0.f;  // of the state on the chamber surface
  };

  struct GEMHitRow {
    int16_t chamberIdx;
    const reco::MuonGEMHitMatch* hit;
    float dx;
  };

  struct LinkRow {
    int16_t lctIdx;
    int16_t muonIdx;
    const CSCLCTRecoMuonMatcher::MuonMatch* match;
  };

}  // namespace

class CSCLCTRecoMuonTableProducer : public edm::stream::EDProducer<> {
public:
  explicit CSCLCTRecoMuonTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lctToken_;
  const edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  const StringCutObjectSelector<reco::Muon> muonCut_;
  const CSCLCTRecoMuonMatcher matcher_;
  // the two explicit segment sources, for the MuonCSCChamber flags (their windows are not used)
  const CSCLCTRecoMuonMatcher standaloneSource_;
  const CSCLCTRecoMuonMatcher matchesSource_;
};

CSCLCTRecoMuonTableProducer::CSCLCTRecoMuonTableProducer(const edm::ParameterSet& ps)
    : lctToken_(consumes(ps.getParameter<edm::InputTag>("lcts"))),
      muonToken_(consumes(ps.getParameter<edm::InputTag>("muons"))),
      vertexToken_(consumes(ps.getParameter<edm::InputTag>("vertices"))),
      cscGeomToken_(esConsumes()),
      gemGeomToken_(esConsumes()),
      muonCut_(ps.getParameter<std::string>("muonCut")),
      matcher_(ps.getParameter<edm::ParameterSet>("matcher")),
      standaloneSource_(CSCLCTRecoMuonMatcher::Config(CSCLCTRecoMuonMatcher::SegmentSource::Standalone,
                                                      matcher_.config().maxDeltaHalfStrip,
                                                      matcher_.config().maxDeltaWireGroup,
                                                      matcher_.config().lctBXWindow)),
      matchesSource_(CSCLCTRecoMuonMatcher::Config(CSCLCTRecoMuonMatcher::SegmentSource::Matches,
                                                   matcher_.config().maxDeltaHalfStrip,
                                                   matcher_.config().maxDeltaWireGroup,
                                                   matcher_.config().lctBXWindow)) {
  produces<nanoaod::FlatTable>("lct");
  produces<nanoaod::FlatTable>("muon");
  produces<nanoaod::FlatTable>("lctMuonMatch");
  produces<nanoaod::FlatTable>("muonCSCChamber");
  produces<nanoaod::FlatTable>("muonGEMChamber");
  produces<nanoaod::FlatTable>("muonGEMHitMatch");
  produces<nanoaod::FlatTable>("muonCSCSegment");
  produces<nanoaod::FlatTable>("muonCSCLayer");
  produces<nanoaod::FlatTable>("muonCSCRecHitADC");
}

void CSCLCTRecoMuonTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("lcts", edm::InputTag("muonCSCDigis", "MuonCSCCorrelatedLCTDigi"))
      ->setComment("LCT collection: unpacked data LCTs, or emulator LCTs (then consider matcher.lctBXWindow)");
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<std::string>("muonCut", "pt > 10 && (passed(\"CutBasedIdMedium\") || passed(\"CutBasedIdTight\"))")
      ->setComment("reco::Muon selection: only these muons are stored and matched");
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"))
      ->setComment("Muon_dxy and Muon_dz are computed w.r.t. the first vertex, as Muon_tightId in RECO");
  edm::ParameterSetDescription matcher;
  CSCLCTRecoMuonMatcher::fillPSetDescription(matcher);
  desc.add<edm::ParameterSetDescription>("matcher", matcher);
  descriptions.add("cscLCTRecoMuonTableProducer", desc);
}

void CSCLCTRecoMuonTableProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  reco::MuonCollection muons;  // the muons passing muonCut, in collection order: the Muon rows
  for (const auto& mu : event.get(muonToken_))
    if (muonCut_(mu))
      muons.push_back(mu);
  const auto& vertices = event.get(vertexToken_);
  if (vertices.empty())
    throw cms::Exception("ProductNotValid") << "empty vertex collection: Muon_dxy and Muon_dz need a primary vertex";
  const reco::Vertex::Point& pv = vertices.front().position();
  const auto matches = matcher_.match(event.get(lctToken_), muons, setup.getData(cscGeomToken_));
  checkIndexRange(matches.size(), "LCTs");
  checkIndexRange(muons.size(), "muons");

  std::vector<LinkRow> links;
  std::vector<int16_t> muonNLCTMatch(muons.size(), 0);
  for (size_t i = 0; i < matches.size(); ++i)
    for (const auto& mm : matches[i].muons) {
      links.push_back({int16_t(i), int16_t(mm.muonIndex), &mm});
      ++muonNLCTMatch[mm.muonIndex];
    }

  using LCTMatch = CSCLCTRecoMuonMatcher::LCTMatch;
  auto lctTable = std::make_unique<nanoaod::FlatTable>(matches.size(), "LCT", false);
  lctTable->setDoc("CSC LCTs (valid, inside the matcher's lctBXWindow), in the order of the LCT collection");
  addColumn<int16_t>(*lctTable, "endcap", matches, [](const LCTMatch& m) { return m.chamber.zendcap(); }, "+1 or -1");
  addColumn<int16_t>(*lctTable, "station", matches, [](const LCTMatch& m) { return m.chamber.station(); }, "station");
  addColumn<int16_t>(
      *lctTable, "ring", matches, [](const LCTMatch& m) { return m.chamber.ring(); }, "ring (ME1/1 LCTs: 1)");
  addColumn<int16_t>(*lctTable, "chamber", matches, [](const LCTMatch& m) { return m.chamber.chamber(); }, "chamber");
  addColumn<int16_t>(
      *lctTable,
      "hs",
      matches,
      [](const LCTMatch& m) { return m.lct->getStrip(); },
      "key half-strip, 0-based (ME1/1: ME1b 0-127, ME1a 128-223)");
  addColumn<int16_t>(
      *lctTable, "es", matches, [](const LCTMatch& m) { return m.lct->getStrip(8); }, "key eighth-strip, getStrip(8)");
  addColumn<int16_t>(
      *lctTable, "wg", matches, [](const LCTMatch& m) { return m.lct->getKeyWG(); }, "key wire group, 0-based");
  addColumn<int16_t>(*lctTable, "quality", matches, [](const LCTMatch& m) { return m.lct->getQuality(); }, "quality");
  addColumn<int16_t>(
      *lctTable, "pattern", matches, [](const LCTMatch& m) { return m.lct->getPattern(); }, "Run-2 pattern");
  addColumn<int16_t>(
      *lctTable,
      "run3Pattern",
      matches,
      [](const LCTMatch& m) { return m.lct->getRun3Pattern(); },
      "Run-3 pattern (0 for Run3-b unpacked data)");
  addColumn<int16_t>(
      *lctTable, "slope", matches, [](const LCTMatch& m) { return m.lct->getSlope(); }, "slope, getSlope() (4 bits)");
  addColumn<int16_t>(
      *lctTable,
      "slope6",
      matches,
      [](const LCTMatch& m) { return m.lct->getSlopeEx(6); },
      "slope at 6-bit resolution, getSlopeEx(6)");
  addColumn<int16_t>(*lctTable, "bend", matches, [](const LCTMatch& m) { return m.lct->getBend(); }, "bend (L/R)");
  addColumn<int16_t>(
      *lctTable, "bx", matches, [](const LCTMatch& m) { return m.lct->getBX(); }, "BX (1 bit in unpacked data)");
  addColumn<int16_t>(
      *lctTable, "trknmb", matches, [](const LCTMatch& m) { return m.lct->getTrknmb(); }, "track number (1 or 2)");
  addColumn<int16_t>(
      *lctTable,
      "type",
      matches,
      [](const LCTMatch& m) { return m.lct->getType(); },
      "LCT type (CSCCorrelatedLCTDigi::Type; emulator only)");
  addColumn<int16_t>(
      *lctTable,
      "gemLayer",
      matches,
      [](const LCTMatch& m) { return m.lct->getGemLayerUsedForSlopeComputation(); },
      "GEM layer used for the slope: emulator 1, 2, 0 = none; Run3-b unpacked data the raw 1-bit field");
  addColumn<int16_t>(
      *lctTable,
      "nMuonMatch",
      matches,
      [](const LCTMatch& m) { return m.muons.size(); },
      "number of LCTMuonMatch rows of this LCT");

  // Muon rows are the selected muons, so a row index is the muon's index in the selected (muonCut) collection.
  std::vector<size_t> muonRows(muons.size());
  for (size_t m = 0; m < muons.size(); ++m)
    muonRows[m] = m;
  const auto mu = [&muons](size_t m) -> const reco::Muon& { return muons[m]; };
  auto muonTable = std::make_unique<nanoaod::FlatTable>(muons.size(), "Muon", false);
  muonTable->setDoc(
      "reco muons passing muonCut (default pT > 10 GeV and medium or tight cut-based ID), in the order "
      "of the muon collection");
  addColumn<float>(*muonTable, "pt", muonRows, [&](size_t m) { return mu(m).pt(); }, "pt");
  addColumn<float>(*muonTable, "eta", muonRows, [&](size_t m) { return mu(m).eta(); }, "eta");
  addColumn<float>(*muonTable, "phi", muonRows, [&](size_t m) { return mu(m).phi(); }, "phi");
  addColumn<float>(*muonTable, "mass", muonRows, [&](size_t m) { return mu(m).mass(); }, "mass");
  addColumn<int16_t>(*muonTable, "charge", muonRows, [&](size_t m) { return mu(m).charge(); }, "electric charge");
  addColumn<bool>(*muonTable, "isGlobal", muonRows, [&](size_t m) { return mu(m).isGlobalMuon(); }, "global muon");
  addColumn<bool>(*muonTable, "isTracker", muonRows, [&](size_t m) { return mu(m).isTrackerMuon(); }, "tracker muon");
  addColumn<bool>(
      *muonTable, "isStandalone", muonRows, [&](size_t m) { return mu(m).isStandAloneMuon(); }, "standalone muon");
  addColumn<bool>(*muonTable, "isPFcand", muonRows, [&](size_t m) { return mu(m).isPFMuon(); }, "PF candidate");
  addColumn<bool>(
      *muonTable,
      "looseId",
      muonRows,
      [&](size_t m) { return mu(m).passed(reco::Muon::CutBasedIdLoose); },
      "cut-based ID, loose WP");
  addColumn<bool>(
      *muonTable,
      "mediumId",
      muonRows,
      [&](size_t m) { return mu(m).passed(reco::Muon::CutBasedIdMedium); },
      "cut-based ID, medium WP");
  addColumn<bool>(
      *muonTable,
      "tightId",
      muonRows,
      [&](size_t m) { return mu(m).passed(reco::Muon::CutBasedIdTight); },
      "cut-based ID, tight WP (set in RECO w.r.t. the first offlinePrimaryVertices vertex)");
  addColumn<float>(
      *muonTable,
      "pfRelIso04_all",
      muonRows,
      [&](size_t m) {
        const auto& iso = mu(m).pfIsolationR04();
        return (iso.sumChargedHadronPt + std::max(iso.sumNeutralHadronEt + iso.sumPhotonEt - iso.sumPUPt / 2.f, 0.f)) /
               mu(m).pt();
      },
      "PF relative isolation dR=0.4, total (deltaBeta corrections)");
  addColumn<float>(
      *muonTable,
      "dxy",
      muonRows,
      [&](size_t m) { return mu(m).muonBestTrack()->dxy(pv); },
      "dxy (with sign) of the best track wrt first PV, in cm");
  addColumn<float>(
      *muonTable,
      "dz",
      muonRows,
      [&](size_t m) { return mu(m).muonBestTrack()->dz(pv); },
      "dz (with sign) of the best track wrt first PV, in cm");
  addColumn<int16_t>(
      *muonTable,
      "nStations",
      muonRows,
      [&](size_t m) { return mu(m).numberOfMatchedStations(); },
      "number of matched stations with default arbitration (segment & track)");
  addColumn<float>(
      *muonTable,
      "standalonePt",
      muonRows,
      [&](size_t m) { return mu(m).outerTrack().isNonnull() ? mu(m).outerTrack()->pt() : 0.f; },
      "pt of the standalone track, 0 if none");
  addColumn<float>(
      *muonTable,
      "standaloneEta",
      muonRows,
      [&](size_t m) { return mu(m).outerTrack().isNonnull() ? mu(m).outerTrack()->eta() : 0.f; },
      "eta of the standalone track, 0 if none (standalonePt == 0)");
  addColumn<float>(
      *muonTable,
      "standalonePhi",
      muonRows,
      [&](size_t m) { return mu(m).outerTrack().isNonnull() ? mu(m).outerTrack()->phi() : 0.f; },
      "phi of the standalone track, 0 if none (standalonePt == 0)");
  addColumn<int16_t>(
      *muonTable,
      "nCSCSegments",
      muonRows,
      [&](size_t m) { return matcher_.muonSegments(mu(m)).size(); },
      "number of CSC segments of this muon for matcher.segmentSource, including those without a valid key position");
  addColumn<int16_t>(
      *muonTable,
      "nLCTMatch",
      muonRows,
      [&](size_t m) { return muonNLCTMatch[m]; },
      "number of LCTMuonMatch rows of this muon");

  using Source = CSCLCTRecoMuonMatcher::SegmentSource;
  auto linkTable = std::make_unique<nanoaod::FlatTable>(links.size(), "LCTMuonMatch", false);
  linkTable->setDoc(
      "LCT-muon pairs within the matcher windows (many-to-many), grouped by LCT, closest muon first; residuals are "
      "LCT - segment at the key layer (layer 3), LCT units");
  addColumn<int16_t>(*linkTable, "lctIdx", links, [](const LinkRow& l) { return l.lctIdx; }, "row in LCT");
  addColumn<int16_t>(*linkTable, "muonIdx", links, [](const LinkRow& l) { return l.muonIdx; }, "row in Muon");
  addColumn<uint8_t>(
      *linkTable,
      "segSource",
      links,
      [](const LinkRow& l) { return l.match->source == Source::Standalone ? 0 : 1; },
      "segment taken from 0 = the standalone-muon fit, 1 = reco::Muon::matches()");
  addColumn<uint8_t>(
      *linkTable,
      "segRing",
      links,
      [](const LinkRow& l) { return l.match->segment->cscDetId().ring(); },
      "ring of the segment (ME1/1: 1 = ME1b, 4 = ME1a)");
  addColumn<uint8_t>(
      *linkTable,
      "segNHits",
      links,
      [](const LinkRow& l) { return l.match->segment->nRecHits(); },
      "number of rechits of the segment");
  addColumn<float>(
      *linkTable,
      "segChi2ndof",
      links,
      [](const LinkRow& l) { return l.match->segment->chi2() / std::max(1, l.match->segment->degreesOfFreedom()); },
      "segment chi2 / max(1, ndof)");
  addColumn<float>(
      *linkTable, "segTime", links, [](const LinkRow& l) { return l.match->segment->time(); }, "segment time, ns");
  addColumn<bool>(
      *linkTable,
      "keyFromRecHit",
      links,
      [](const LinkRow& l) { return l.match->keyFromRecHit; },
      "key position from the segment's layer-3 rechit, else from the segment projected onto layer 3");
  addColumn<int16_t>(
      *linkTable, "keyHS", links, [](const LinkRow& l) { return l.match->key.hs; }, "segment key half-strip");
  addColumn<float>(
      *linkTable, "keyES", links, [](const LinkRow& l) { return l.match->key.es; }, "segment eighth-strip (CCLUT)");
  addColumn<int16_t>(
      *linkTable, "keyWG", links, [](const LinkRow& l) { return l.match->key.wg; }, "segment key wire group");
  addColumn<int16_t>(
      *linkTable, "dHS", links, [](const LinkRow& l) { return l.match->deltaHalfStrip; }, "LCT - segment half-strip");
  addColumn<float>(
      *linkTable,
      "dES",
      links,
      [](const LinkRow& l) { return l.match->deltaEighthStrip; },
      "LCT - segment eighth-strip");
  addColumn<int16_t>(
      *linkTable, "dWG", links, [](const LinkRow& l) { return l.match->deltaWireGroup; }, "LCT - segment wire group");

  // MuonCSCChamber: the CSC chamber matches of the stored muons, in muon order, then in the order of matches()
  const CSCGeometry& geometry = setup.getData(cscGeomToken_);
  std::vector<ChamberRow> chambers;
  for (size_t m = 0; m < muons.size(); ++m) {
    // chambers (with ring) where the muon has a segment, per source
    auto segmentChambers = [&](const CSCLCTRecoMuonMatcher& source) {
      std::set<CSCDetId> ids;
      for (const auto& [seg, from] : source.muonSegments(muons[m]))
        ids.insert(seg->cscDetId().chamberId());
      return ids;
    };
    const auto matcherChambers = segmentChambers(matcher_);
    const auto standaloneChambers = segmentChambers(standaloneSource_);
    const auto arbitratedChambers = segmentChambers(matchesSource_);
    for (const auto& cm : muons[m].matches()) {
      if (cm.detector() != MuonSubdetId::CSC)
        continue;
      ChamberRow row{int16_t(m), &cm, CSCDetId(cm.id.rawId()), {}};
      row.hasMatcherSegment = matcherChambers.count(row.id) > 0;
      row.hasStandaloneSegment = standaloneChambers.count(row.id) > 0;
      row.hasArbitratedSegment = arbitratedChambers.count(row.id) > 0;
      // the chamber match is the track state on the chamber surface (local z = 0): a straight line to layer 3
      const auto p3 =
          csclct::lineAtKeyLayer(geometry, row.id, LocalPoint(cm.x, cm.y, 0.f), LocalVector(cm.dXdZ, cm.dYdZ, 1.f));
      if (p3) {
        row.key = csclct::toLCTSpace(geometry, row.id, *p3);
        const CSCLayer* L3 = geometry.chamber(row.id)->layer(CSCConstants::KEY_CLCT_LAYER);
        const auto limits = L3->geometry()->yLimitsOfStripPlane();
        const float y = L3->toLocal(*p3).y();
        row.inStripPlane = limits.first <= y && y <= limits.second;
        row.eta = p3->eta();
        row.phi = p3->phi();
      }
      chambers.push_back(row);
    }
  }
  auto chamberTable = std::make_unique<nanoaod::FlatTable>(chambers.size(), "MuonCSCChamber", false);
  chamberTable->setDoc(
      "CSC chamber matches of the stored muons (reco::Muon::matches()): the muon's track (the inner track if any, "
      "else the standalone track) propagated by reco (TrackDetectorAssociator) to every CSC chamber whose edge "
      "distances are both below 5 cm (muons1stStep; none if no chamber within dR 0.2 of the track has a segment or "
      "hit, preselectMuonTracks); in muon order, then in the order of matches(). ME1/1: normally a ring-1 (ME1b) and "
      "a ring-4 (ME1a) row with nearly the same state and edges; keyInStripPlane tells which ring the track "
      "crosses");
  addColumn<int16_t>(*chamberTable, "muonIdx", chambers, [](const ChamberRow& r) { return r.muonIdx; }, "row in Muon");
  addColumn<int16_t>(*chamberTable, "endcap", chambers, [](const ChamberRow& r) { return r.id.zendcap(); }, "+1 or -1");
  addColumn<int16_t>(*chamberTable, "station", chambers, [](const ChamberRow& r) { return r.id.station(); }, "station");
  addColumn<int16_t>(
      *chamberTable,
      "ring",
      chambers,
      [](const ChamberRow& r) { return r.id.ring(); },
      "ring of the chamber match (ME1/1: 1 = ME1b and 4 = ME1a, normally both rows for a crossing, with nearly the "
      "same state and edges; ME1/1 LCTs are in ring 1)");
  addColumn<int16_t>(*chamberTable, "chamber", chambers, [](const ChamberRow& r) { return r.id.chamber(); }, "chamber");
  addColumn<float>(
      *chamberTable,
      "x",
      chambers,
      [](const ChamberRow& r) { return r.match->x; },
      "propagated track position x, chamber frame, on the chamber surface (local z = 0), cm");
  addColumn<float>(
      *chamberTable,
      "y",
      chambers,
      [](const ChamberRow& r) { return r.match->y; },
      "propagated track position y, chamber frame, on the chamber surface (local z = 0), cm");
  addColumn<float>(
      *chamberTable, "xErr", chambers, [](const ChamberRow& r) { return r.match->xErr; }, "propagation error on x, cm");
  addColumn<float>(
      *chamberTable, "yErr", chambers, [](const ChamberRow& r) { return r.match->yErr; }, "propagation error on y, cm");
  addColumn<float>(
      *chamberTable,
      "dXdZ",
      chambers,
      [](const ChamberRow& r) { return r.match->dXdZ; },
      "propagated track direction dx/dz, chamber frame (9999 if the direction has no z component)");
  addColumn<float>(
      *chamberTable,
      "dYdZ",
      chambers,
      [](const ChamberRow& r) { return r.match->dYdZ; },
      "propagated track direction dy/dz, chamber frame (9999 if the direction has no z component)");
  addColumn<float>(
      *chamberTable,
      "dXdZErr",
      chambers,
      [](const ChamberRow& r) { return r.match->dXdZErr; },
      "propagation error on dx/dz (0 if the variance is not positive)");
  addColumn<float>(
      *chamberTable,
      "dYdZErr",
      chambers,
      [](const ChamberRow& r) { return r.match->dYdZErr; },
      "propagation error on dy/dz (0 if the variance is not positive)");
  addColumn<float>(
      *chamberTable,
      "edgeX",
      chambers,
      [](const ChamberRow& r) { return r.match->edgeX; },
      "distance in x of the propagated position to the edge of the chamber type's odd-layer wire plane, cm (< 0: "
      "inside; ME1/1: the whole ME1/1 plane, the same for both rings)");
  addColumn<float>(
      *chamberTable,
      "edgeY",
      chambers,
      [](const ChamberRow& r) { return r.match->edgeY; },
      "distance in y of the propagated position to the edge of the chamber type's odd-layer wire plane, cm (< 0: "
      "inside; ME1/1: the whole ME1/1 plane, the same for both rings)");
  addColumn<int16_t>(
      *chamberTable,
      "nDigisInRange",
      chambers,
      [](const ChamberRow& r) { return r.match->nDigisInRange; },
      "fired strip digis (an ADC sample > pedestal + 13.3) in layers 1-6 of this ring within 25 cm in x of the "
      "propagated position, strip x at the track's y (reco MuonShowerDigiFiller; 0 also when reco did not fill it)");
  addColumn<int16_t>(
      *chamberTable,
      "keyHS",
      chambers,
      [](const ChamberRow& r) { return r.key.hs; },
      "propagated position carried to layer 3 (straight line, chamber frame), converted with this "
      "ring's layer-3 geometry: key half-strip, LCT units (ME1a: 128-223), clamped to the ring's range when the "
      "point is outside the strips");
  addColumn<float>(
      *chamberTable,
      "keyES",
      chambers,
      [](const ChamberRow& r) { return r.key.es; },
      "the same point: eighth-strip, CCLUT convention (4*keyHS at the half-strip centre); saturates at -2 and "
      "8N-2 (ME1a: 510 and 894) when the point is outside the strips (strip coordinate clamped to [0, N])");
  addColumn<int16_t>(
      *chamberTable,
      "keyWG",
      chambers,
      [](const ChamberRow& r) { return r.key.wg; },
      "the same point: wire group of the nearest wire, 0-based; -1 beyond the first or last wire or at a dead "
      "wire (not bounded along the wires)");
  addColumn<bool>(
      *chamberTable,
      "keyValid",
      chambers,
      [](const ChamberRow& r) { return r.key.valid(); },
      "the layer-3 point has a strip coordinate inside (0, N) and a wire group; the ring's y range is not "
      "tested (ME1/1: often true in both rows, use keyInStripPlane)");
  addColumn<bool>(
      *chamberTable,
      "keyInStripPlane",
      chambers,
      [](const ChamberRow& r) { return r.inStripPlane; },
      "the layer-3 point lies inside the y range of this ring's strip plane (ME1/1: the ring the track crosses; "
      "neither row beyond the outer ends of ME1/1)");
  addColumn<float>(
      *chamberTable, "eta", chambers, [](const ChamberRow& r) { return r.eta; }, "eta of the layer-3 point");
  addColumn<float>(
      *chamberTable, "phi", chambers, [](const ChamberRow& r) { return r.phi; }, "phi of the layer-3 point");
  addColumn<int16_t>(
      *chamberTable,
      "nSegmentMatches",
      chambers,
      [](const ChamberRow& r) { return r.match->segmentMatches.size(); },
      "CSC segments reco matched to the propagated track in this chamber, whatever their arbitration flags");
  addColumn<bool>(
      *chamberTable,
      "hasMatcherSegment",
      chambers,
      [](const ChamberRow& r) { return r.hasMatcherSegment; },
      "the matcher's segmentSource gives this muon a segment in this chamber (ring included), with or without a "
      "valid key position (the matcher skips those without)");
  addColumn<bool>(
      *chamberTable,
      "hasStandaloneSegment",
      chambers,
      [](const ChamberRow& r) { return r.hasStandaloneSegment; },
      "the muon's standalone track has a CSC segment in this chamber (segmentSource standalone)");
  addColumn<bool>(
      *chamberTable,
      "hasArbitratedSegment",
      chambers,
      [](const ChamberRow& r) { return r.hasArbitratedSegment; },
      "a segment with BestInChamberByDR and BelongsToTrackByDR in this chamber (segmentSource matches)");

  // MuonCSCSegment: the CSC segments of the stored muons, in muon order. Per muon, source 0 first: every valid 4D CSC
  // hit (CSCSegment) of the standalone track, in outerTrack()->recHits() order; then source 1: every segment match of
  // every CSC chamber match, in matches() then segmentMatches order, whatever its arbitration bits.
  std::vector<SegmentRow> segments;
  for (size_t m = 0; m < muons.size(); ++m) {
    const reco::TrackRef trk = muons[m].outerTrack();
    if (trk.isNonnull()) {
      if (!trk.isAvailable() || !trk->recHitsOk())
        throw cms::Exception("ProductNotFound")
            << "muon " << m << ": its standalone track or the track's rechits are not in the input";
      for (auto const* hit : trk->recHits()) {
        if (!hit->isValid() || hit->dimension() != 4)
          continue;
        const DetId id = hit->geographicalId();
        if (id.det() != DetId::Muon || id.subdetId() != MuonSubdetId::CSC)
          continue;
        auto const* seg = dynamic_cast<const CSCSegment*>(hit);
        if (seg == nullptr)
          throw cms::Exception("LogicError")
              << "muon " << m << ": a 4D CSC hit of the standalone track is not a CSCSegment";
        segments.push_back({int16_t(m), 0, -1, 0u, seg});
      }
    }
    for (size_t r = 0; r < chambers.size(); ++r) {
      if (chambers[r].muonIdx != int16_t(m))
        continue;
      for (const auto& sm : chambers[r].match->segmentMatches) {
        if (sm.cscSegmentRef.isNull() || !sm.cscSegmentRef.isAvailable())
          throw cms::Exception("ProductNotFound")
              << "muon " << m << ": a CSC segment match without its segment in the input (" << chambers[r].id << ")";
        segments.push_back({int16_t(m), 1, int16_t(r), sm.mask, &*sm.cscSegmentRef});
      }
    }
  }
  checkIndexRange(segments.size(), "muon CSC segments");

  // MuonCSCLayer: per muon, every CSC chamber (with ring) holding a chamber match or a segment of the muon, the
  // chamber matches' chambers first (MuonCSCChamber order), then the others in segment order; per chamber, every
  // segment of the muon there (MuonCSCSegment order), or none; per segment, layers 1-6.
  std::vector<LayerRow> layers;
  for (size_t m = 0; m < muons.size(); ++m) {
    std::vector<CSCDetId> keys;
    for (const auto& c : chambers)
      if (c.muonIdx == int16_t(m))
        keys.push_back(c.id);
    for (const auto& sr : segments) {
      const CSCDetId id = sr.segment->cscDetId().chamberId();
      if (sr.muonIdx == int16_t(m) && std::find(keys.begin(), keys.end(), id) == keys.end())
        keys.push_back(id);
    }
    for (const CSCDetId& key : keys) {
      int16_t chamberIdx = -1;
      for (size_t r = 0; r < chambers.size(); ++r)
        if (chambers[r].muonIdx == int16_t(m) && chambers[r].id == key)
          chamberIdx = int16_t(r);
      std::vector<int16_t> segIdx;
      for (size_t i = 0; i < segments.size(); ++i)
        if (segments[i].muonIdx == int16_t(m) && segments[i].segment->cscDetId().chamberId() == key)
          segIdx.push_back(int16_t(i));
      if (segIdx.empty())
        segIdx.push_back(-1);
      const CSCChamber* ch = geometry.chamber(key);
      if (ch == nullptr)
        throw cms::Exception("LogicError")
            << key << " is not in the CSC geometry: the job's geometry differs from the RECO's";
      const reco::MuonChamberMatch* cm = chamberIdx >= 0 ? chambers[chamberIdx].match : nullptr;
      for (const int16_t si : segIdx)
        for (int l = 1; l <= CSCConstants::NUM_LAYERS; ++l) {
          const CSCLayer* L = ch->layer(l);
          if (L == nullptr)
            throw cms::Exception("LogicError") << "layer " << l << " of " << key << " is not in the CSC geometry";
          const float halfGap = csclct::halfGasGap(*L);
          LayerRow row{int16_t(m), si, chamberIdx, L->id(), nullptr, {}, {}, {}, {}, {}, {}, {}};
          if (si >= 0) {
            const CSCSegment& seg = *segments[si].segment;
            const LocalPoint p0 = seg.localPosition();
            const LocalVector d = seg.localDirection();
            row.line = csclct::lineAtLayer(*ch, *L, p0, d, 0.f);
            row.lineLo = csclct::lineAtLayer(*ch, *L, p0, d, -halfGap);
            row.lineHi = csclct::lineAtLayer(*ch, *L, p0, d, halfGap);
            for (const auto& rh : seg.specificRecHits())
              if (rh.cscDetId().layer() == l) {
                row.rechit = &rh;
                row.recHit = csclct::layerPoint(*L, rh.localPosition());
                row.recHit.wg = rh.hitWire() - 1;
                break;
              }
          }
          // the chamber match is reco's track state on the chamber surface (local z = 0); none without a direction
          if (cm != nullptr && cm->dXdZ != 9999.f) {
            const LocalPoint p0(cm->x, cm->y, 0.f);
            const LocalVector d(cm->dXdZ, cm->dYdZ, 1.f);
            row.prop = csclct::lineAtLayer(*ch, *L, p0, d, 0.f);
            row.propLo = csclct::lineAtLayer(*ch, *L, p0, d, -halfGap);
            row.propHi = csclct::lineAtLayer(*ch, *L, p0, d, halfGap);
          }
          layers.push_back(row);
        }
    }
  }
  checkIndexRange(layers.size(), "muon CSC layers");
  std::vector<ADCRow> adcs;
  for (size_t i = 0; i < layers.size(); ++i)
    if (const CSCRecHit2D* rh = layers[i].rechit)
      for (unsigned int s = 0; s < rh->nStrips(); ++s)
        for (unsigned int t = 0; t < rh->nTimeBins(); ++t)
          adcs.push_back({int16_t(i), int16_t(rh->channels(s)), int16_t(t), rh->adcs(s, t)});

  auto segmentTable = std::make_unique<nanoaod::FlatTable>(segments.size(), "MuonCSCSegment", false);
  segmentTable->setDoc(
      "CSC segments of the stored muons, in muon order; per muon, source 0 (every valid 4D CSC hit of the standalone "
      "track, outerTrack()->recHits() order), then source 1 (every segment match of every CSC chamber match, "
      "matches() then segmentMatches order, whatever its arbitration bits). A segment can appear under both sources");
  addColumn<int16_t>(*segmentTable, "muonIdx", segments, [](const SegmentRow& r) { return r.muonIdx; }, "row in Muon");
  addColumn<uint8_t>(
      *segmentTable,
      "source",
      segments,
      [](const SegmentRow& r) { return r.source; },
      "0 = a hit of the muon's standalone track, 1 = a segment match of a CSC chamber match");
  addColumn<int16_t>(
      *segmentTable,
      "muonCSCChamberIdx",
      segments,
      [](const SegmentRow& r) { return r.chamberIdx; },
      "row in MuonCSCChamber of the chamber match holding the segment match (source 1); -1 for source 0");
  addColumn<uint32_t>(
      *segmentTable,
      "mask",
      segments,
      [](const SegmentRow& r) { return r.mask; },
      "reco::MuonSegmentMatch::mask, the arbitration bits (e.g. 1<<10 BestInChamberByDR, 1<<18 BelongsToTrackByDR; "
      "DataFormats/MuonReco/interface/MuonSegmentMatch.h); 0 for source 0");
  addColumn<int16_t>(
      *segmentTable,
      "endcap",
      segments,
      [](const SegmentRow& r) { return r.segment->cscDetId().zendcap(); },
      "+1 or -1");
  addColumn<int16_t>(
      *segmentTable,
      "station",
      segments,
      [](const SegmentRow& r) { return r.segment->cscDetId().station(); },
      "station");
  addColumn<int16_t>(
      *segmentTable,
      "ring",
      segments,
      [](const SegmentRow& r) { return r.segment->cscDetId().ring(); },
      "ring (ME1/1: 1 = ME1b, 4 = ME1a)");
  addColumn<int16_t>(
      *segmentTable,
      "chamber",
      segments,
      [](const SegmentRow& r) { return r.segment->cscDetId().chamber(); },
      "chamber");
  addColumn<int16_t>(
      *segmentTable,
      "nHits",
      segments,
      [](const SegmentRow& r) { return r.segment->nRecHits(); },
      "rechits of the segment (at most one per layer)");
  addColumn<float>(
      *segmentTable, "chi2", segments, [](const SegmentRow& r) { return r.segment->chi2(); }, "segment fit chi2");
  addColumn<int16_t>(
      *segmentTable,
      "ndof",
      segments,
      [](const SegmentRow& r) { return r.segment->degreesOfFreedom(); },
      "segment fit degrees of freedom (2 nHits - 4)");
  addColumn<float>(
      *segmentTable, "time", segments, [](const SegmentRow& r) { return r.segment->time(); }, "segment time, ns");
  addColumn<float>(
      *segmentTable,
      "x",
      segments,
      [](const SegmentRow& r) { return r.segment->localPosition().x(); },
      "segment position x, chamber frame (localPosition: at chamber z = 0), cm");
  addColumn<float>(
      *segmentTable,
      "y",
      segments,
      [](const SegmentRow& r) { return r.segment->localPosition().y(); },
      "segment position y, chamber frame (localPosition: at chamber z = 0), cm");
  addColumn<float>(
      *segmentTable,
      "xErr",
      segments,
      [](const SegmentRow& r) { return std::sqrt(r.segment->localPositionError().xx()); },
      "error on x, cm");
  addColumn<float>(
      *segmentTable,
      "yErr",
      segments,
      [](const SegmentRow& r) { return std::sqrt(r.segment->localPositionError().yy()); },
      "error on y, cm");
  addColumn<float>(
      *segmentTable,
      "dirX",
      segments,
      [](const SegmentRow& r) { return r.segment->localDirection().x(); },
      "segment direction (unit vector, localDirection), chamber frame, x");
  addColumn<float>(
      *segmentTable,
      "dirY",
      segments,
      [](const SegmentRow& r) { return r.segment->localDirection().y(); },
      "segment direction (unit vector, localDirection), chamber frame, y");
  addColumn<float>(
      *segmentTable,
      "dirZ",
      segments,
      [](const SegmentRow& r) { return r.segment->localDirection().z(); },
      "segment direction (unit vector, localDirection), chamber frame, z");
  addColumn<float>(
      *segmentTable,
      "dXdZErr",
      segments,
      [](const SegmentRow& r) { return std::sqrt(r.segment->localDirectionError().xx()); },
      "error on dx/dz (localDirectionError)");
  addColumn<float>(
      *segmentTable,
      "dYdZErr",
      segments,
      [](const SegmentRow& r) { return std::sqrt(r.segment->localDirectionError().yy()); },
      "error on dy/dz (localDirectionError)");

  auto layerTable = std::make_unique<nanoaod::FlatTable>(layers.size(), "MuonCSCLayer", false);
  layerTable->setDoc(
      "Per stored muon, every CSC chamber (with ring) holding a chamber match (MuonCSCChamber) or a segment "
      "(MuonCSCSegment) of the muon, chamber matches first, then in segment order; per chamber, every segment of the "
      "muon there, or none; per segment, layers 1-6. Three sources of a position on the layer: the segment's rechit "
      "on it (recHit*), the segment line at the layer (line*), and reco's propagated state carried along its straight "
      "line to the layer (prop*; repeated on the rows of every segment of the chamber). Per point: x, y in the layer "
      "frame; fs = CSCLayerGeometry::strip (strip k covers [k-1, k), with the layer's stagger; CMSSW clamps it to "
      "[0, N]); wg = wire group, 0-based. *Lo / *Hi: where the line crosses the gas-gap surfaces, layer-frame z = "
      "-g/2 and +g/2 (g = the layer thickness, which CMSSW sets to the gas gap)");
  addColumn<int16_t>(*layerTable, "muonIdx", layers, [](const LayerRow& r) { return r.muonIdx; }, "row in Muon");
  addColumn<int16_t>(
      *layerTable,
      "muonCSCSegmentIdx",
      layers,
      [](const LayerRow& r) { return r.segmentIdx; },
      "row in MuonCSCSegment, -1 if none (the muon has a chamber match here but no segment)");
  addColumn<int16_t>(
      *layerTable,
      "muonCSCChamberIdx",
      layers,
      [](const LayerRow& r) { return r.chamberIdx; },
      "row in MuonCSCChamber of this chamber (ring), -1 if none (the muon has a segment here but no chamber match)");
  addColumn<int16_t>(*layerTable, "endcap", layers, [](const LayerRow& r) { return r.id.zendcap(); }, "+1 or -1");
  addColumn<int16_t>(*layerTable, "station", layers, [](const LayerRow& r) { return r.id.station(); }, "station");
  addColumn<int16_t>(
      *layerTable, "ring", layers, [](const LayerRow& r) { return r.id.ring(); }, "ring (ME1/1: 1 = ME1b, 4 = ME1a)");
  addColumn<int16_t>(*layerTable, "chamber", layers, [](const LayerRow& r) { return r.id.chamber(); }, "chamber");
  addColumn<int16_t>(*layerTable, "layer", layers, [](const LayerRow& r) { return r.id.layer(); }, "layer, 1-6");
  addColumn<bool>(
      *layerTable,
      "hasRecHit",
      layers,
      [](const LayerRow& r) { return r.recHit.has; },
      "the segment has a rechit on this layer");
  addColumn<float>(
      *layerTable, "recHitX", layers, [](const LayerRow& r) { return r.recHit.x; }, "rechit local x, layer frame, cm");
  addColumn<float>(
      *layerTable, "recHitY", layers, [](const LayerRow& r) { return r.recHit.y; }, "rechit local y, layer frame, cm");
  addColumn<float>(
      *layerTable, "recHitFs", layers, [](const LayerRow& r) { return r.recHit.fs; }, "fs of the rechit position");
  addColumn<bool>(
      *layerTable,
      "recHitInStrips",
      layers,
      [](const LayerRow& r) { return r.recHit.inStrips; },
      "1e-6 < recHitFs < N - 1e-6 (not clamped)");
  addColumn<bool>(
      *layerTable,
      "recHitInStripPlane",
      layers,
      [](const LayerRow& r) { return r.recHit.inStripPlane; },
      "recHitY inside the strip plane (CSCLayerGeometry::yLimitsOfStripPlane)");
  addColumn<int16_t>(
      *layerTable,
      "recHitWg",
      layers,
      [](const LayerRow& r) { return r.recHit.wg; },
      "the rechit's wire group, hitWire() - 1: the middle of its contiguous wire-group cluster");
  addColumn<int16_t>(
      *layerTable,
      "recHitNWireGroups",
      layers,
      [](const LayerRow& r) { return r.rechit != nullptr ? int(r.rechit->nWireGroups()) : 0; },
      "wire groups of the rechit's cluster, contiguous: recHitWg - n/2 ... recHitWg - n/2 + n - 1 (n/2 rounded down)");
  addColumn<int16_t>(
      *layerTable,
      "recHitNStrips",
      layers,
      [](const LayerRow& r) { return r.rechit != nullptr ? int(r.rechit->nStrips()) : 0; },
      "strips of the rechit (rows of MuonCSCRecHitADC per time bin)");
  addColumn<float>(
      *layerTable,
      "recHitTpeak",
      layers,
      [](const LayerRow& r) { return r.rechit != nullptr ? r.rechit->tpeak() : 0.f; },
      "rechit strip-pulse peak time (tpeak), ns");
  addColumn<int16_t>(
      *layerTable,
      "recHitWireBX",
      layers,
      [](const LayerRow& r) { return r.rechit != nullptr ? r.rechit->wgroupsBX() : 0; },
      "BX of the rechit's middle wire group (wgroupsBX)");
  addColumn<float>(
      *layerTable,
      "recHitXErr",
      layers,
      [](const LayerRow& r) { return r.rechit != nullptr ? std::sqrt(r.rechit->localPositionError().xx()) : 0.f; },
      "rechit error on x, cm");
  addColumn<float>(
      *layerTable,
      "recHitYErr",
      layers,
      [](const LayerRow& r) { return r.rechit != nullptr ? std::sqrt(r.rechit->localPositionError().yy()) : 0.f; },
      "rechit error on y, cm");
  addColumn<int16_t>(
      *layerTable,
      "recHitQuality",
      layers,
      [](const LayerRow& r) { return r.rechit != nullptr ? r.rechit->quality() : 0; },
      "rechit quality (CSCRecHit2D::quality)");
  addColumn<int16_t>(
      *layerTable,
      "recHitBadStrip",
      layers,
      [](const LayerRow& r) { return r.rechit != nullptr ? r.rechit->badStrip() : 0; },
      "rechit bad-strip flag (CSCRecHit2D::badStrip)");
  addColumn<int16_t>(
      *layerTable,
      "recHitBadWireGroup",
      layers,
      [](const LayerRow& r) { return r.rechit != nullptr ? r.rechit->badWireGroup() : 0; },
      "rechit bad-wire-group flag (CSCRecHit2D::badWireGroup)");
  // the line and the propagated state: the same columns
  for (const auto& [prefix, what] : std::vector<std::pair<std::string, std::string>>{
           {"line", "the segment line (localPosition, localDirection, chamber frame)"},
           {"prop", "reco's propagated state (MuonCSCChamber x, y at chamber z = 0, dXdZ, dYdZ)"}}) {
    const bool isLine = prefix == "line";
    auto pt = [isLine](const LayerRow& r) -> const csclct::LayerPoint& { return isLine ? r.line : r.prop; };
    auto lo = [isLine](const LayerRow& r) -> const csclct::LayerPoint& { return isLine ? r.lineLo : r.propLo; };
    auto hi = [isLine](const LayerRow& r) -> const csclct::LayerPoint& { return isLine ? r.lineHi : r.propHi; };
    addColumn<bool>(
        *layerTable,
        "has" + std::string(1, char(std::toupper(prefix[0]))) + prefix.substr(1),
        layers,
        [pt](const LayerRow& r) { return pt(r).has; },
        what + " crosses this layer (none: no " + (isLine ? "segment" : "chamber match, or no direction") +
            " here, or a line parallel to the layer)");
    addColumn<float>(
        *layerTable,
        prefix + "X",
        layers,
        [pt](const LayerRow& r) { return pt(r).x; },
        "x where " + what + " crosses the layer plane, layer frame, cm");
    addColumn<float>(
        *layerTable,
        prefix + "Y",
        layers,
        [pt](const LayerRow& r) { return pt(r).y; },
        "y where " + what + " crosses the layer plane, layer frame, cm");
    addColumn<float>(
        *layerTable, prefix + "Fs", layers, [pt](const LayerRow& r) { return pt(r).fs; }, "fs at that point");
    addColumn<bool>(
        *layerTable,
        prefix + "InStrips",
        layers,
        [pt](const LayerRow& r) { return pt(r).inStrips; },
        "1e-6 < " + prefix + "Fs < N - 1e-6 (not clamped)");
    addColumn<bool>(
        *layerTable,
        prefix + "InStripPlane",
        layers,
        [pt](const LayerRow& r) { return pt(r).inStripPlane; },
        prefix + "Y inside the strip plane (CSCLayerGeometry::yLimitsOfStripPlane)");
    addColumn<int16_t>(
        *layerTable,
        prefix + "Wg",
        layers,
        [pt](const LayerRow& r) { return pt(r).wg; },
        "wire group of the nearest wire at that point, 0-based; -1 outside the wires or at a dead wire");
    addColumn<float>(
        *layerTable,
        prefix + "FsLo",
        layers,
        [lo](const LayerRow& r) { return lo(r).fs; },
        "fs where the line crosses the gas-gap surface at layer-frame z = -g/2 (clamped like fs)");
    addColumn<int16_t>(
        *layerTable,
        prefix + "WgLo",
        layers,
        [lo](const LayerRow& r) { return lo(r).wg; },
        "wire group there (the line crosses the contiguous wire groups between WgLo and WgHi)");
    addColumn<float>(
        *layerTable,
        prefix + "FsHi",
        layers,
        [hi](const LayerRow& r) { return hi(r).fs; },
        "fs where the line crosses the gas-gap surface at layer-frame z = +g/2 (clamped like fs)");
    addColumn<int16_t>(
        *layerTable, prefix + "WgHi", layers, [hi](const LayerRow& r) { return hi(r).wg; }, "wire group there");
  }

  auto adcTable = std::make_unique<nanoaod::FlatTable>(adcs.size(), "MuonCSCRecHitADC", false);
  adcTable->setDoc(
      "Strip ADC samples of the rechits in MuonCSCLayer (CSCRecHit2D::adcs): one row per (MuonCSCLayer row with a "
      "rechit, strip of the rechit, time bin)");
  addColumn<int16_t>(
      *adcTable, "muonCSCLayerIdx", adcs, [](const ADCRow& r) { return r.layerIdx; }, "row in MuonCSCLayer");
  addColumn<int16_t>(
      *adcTable,
      "strip",
      adcs,
      [](const ADCRow& r) { return r.strip; },
      "strip number (CSCRecHit2D::channels), 1-based");
  addColumn<int16_t>(
      *adcTable, "timeBin", adcs, [](const ADCRow& r) { return r.timeBin; }, "time bin of the sample, 0-based");
  addColumn<float>(*adcTable, "adc", adcs, [](const ADCRow& r) { return r.adc; }, "ADC count (pedestal subtracted)");

  // MuonGEMChamber: the GEM chamber matches of the stored muons, in muon order, then in the order of matches()
  const GEMGeometry& gemGeometry = setup.getData(gemGeomToken_);
  std::vector<GEMChamberRow> gemChambers;
  for (size_t m = 0; m < muons.size(); ++m)
    for (const auto& cm : muons[m].matches()) {
      if (cm.detector() != MuonSubdetId::GEM)
        continue;
      GEMChamberRow row{int16_t(m), &cm, GEMDetId(cm.id.rawId())};
      const GeomDet* det = row.id.layer() > 0 ? static_cast<const GeomDet*>(gemGeometry.chamber(row.id))
                                              : static_cast<const GeomDet*>(gemGeometry.superChamber(row.id));
      if (det == nullptr)
        throw cms::Exception("LogicError") << "GEM chamber match " << row.id
                                           << " is not in the GEM geometry: the job's geometry differs from the RECO's";
      // the chamber match is the track state on the chamber surface (local z = 0)
      const GlobalPoint gp = det->toGlobal(LocalPoint(cm.x, cm.y, 0.f));
      row.eta = gp.eta();
      row.phi = gp.phi();
      // the eta partitions lie off the chamber surface: follow the straight line through that state to each
      // partition's plane (as MuonCSCChamber does to layer 3); no crossing if reco gives no direction (dXdZ = 9999,
      // the local direction has no z component)
      if (row.id.layer() > 0 && cm.dXdZ != 9999.f) {
        const GlobalVector dir = det->toGlobal(LocalVector(cm.dXdZ, cm.dYdZ, 1.f));
        for (const GEMEtaPartition* part : gemGeometry.chamber(row.id)->etaPartitions()) {
          const LocalPoint p0 = part->toLocal(gp);
          const LocalVector d = part->toLocal(dir);
          const LocalPoint lp = p0 - (p0.z() / d.z()) * d;
          if (part->surface().bounds().inside(Local2DPoint(lp.x(), lp.y()))) {
            row.roll = part->id().roll();
            row.strip = part->strip(lp);
            row.pad = part->pad(lp);
            break;
          }
        }
      }
      gemChambers.push_back(row);
    }
  checkIndexRange(gemChambers.size(), "GEM chamber matches");
  std::vector<GEMHitRow> gemHits;
  for (size_t r = 0; r < gemChambers.size(); ++r)
    for (const auto& hit : gemChambers[r].match->gemHitMatches)
      gemHits.push_back({int16_t(r), &hit, hit.x - gemChambers[r].match->x});

  auto gemTable = std::make_unique<nanoaod::FlatTable>(gemChambers.size(), "MuonGEMChamber", false);
  gemTable->setDoc(
      "GEM chamber matches of the stored muons (reco::Muon::matches(), GEM subdetector): the muon's track (the inner "
      "track if any, else the standalone track) propagated by reco (TrackDetectorAssociator) to every GEM chamber "
      "whose edge distances are both below 5 cm (muons1stStep; no chamber matches at all if no muon chamber within "
      "dR 0.2 of the track has a segment or hit, preselectMuonTracks, where GEM rechits do not count). RECO made with "
      "CMSSW_12_5_0 or later: GE1/1 and GE2/1 rows are per-layer chambers, with matched rechits (MuonGEMHitMatch) "
      "only if the reco step had GEM rechits, and GE0 rows are super-chambers (layer 0) with GEM segments; older "
      "RECO: super-chamber rows for every station, no rechits. In muon order, then in the order of matches()");
  addColumn<int16_t>(
      *gemTable, "muonIdx", gemChambers, [](const GEMChamberRow& r) { return r.muonIdx; }, "row in Muon");
  addColumn<int16_t>(
      *gemTable, "region", gemChambers, [](const GEMChamberRow& r) { return r.id.region(); }, "+1 or -1");
  addColumn<int16_t>(
      *gemTable,
      "station",
      gemChambers,
      [](const GEMChamberRow& r) { return r.id.station(); },
      "0 = GE0, 1 = GE1/1, 2 = GE2/1");
  addColumn<int16_t>(*gemTable, "ring", gemChambers, [](const GEMChamberRow& r) { return r.id.ring(); }, "ring");
  addColumn<int16_t>(
      *gemTable,
      "layer",
      gemChambers,
      [](const GEMChamberRow& r) { return r.id.layer(); },
      "layer of the chamber match (1-2; 0 = a super-chamber: GE0, or any station in RECO before CMSSW_12_5_0)");
  addColumn<int16_t>(
      *gemTable, "chamber", gemChambers, [](const GEMChamberRow& r) { return r.id.chamber(); }, "chamber");
  addColumn<float>(
      *gemTable,
      "x",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->x; },
      "propagated track position x, frame of the row's chamber or super-chamber, on its surface (local z = 0), cm");
  addColumn<float>(
      *gemTable,
      "y",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->y; },
      "propagated track position y, frame of the row's chamber or super-chamber, on its surface (local z = 0), cm");
  addColumn<float>(
      *gemTable,
      "xErr",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->xErr; },
      "propagation error on x, cm");
  addColumn<float>(
      *gemTable,
      "yErr",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->yErr; },
      "propagation error on y, cm");
  addColumn<float>(
      *gemTable,
      "dXdZ",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->dXdZ; },
      "propagated track direction dx/dz, chamber frame (9999 if the direction has no z component)");
  addColumn<float>(
      *gemTable,
      "dYdZ",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->dYdZ; },
      "propagated track direction dy/dz, chamber frame (9999 if the direction has no z component)");
  addColumn<float>(
      *gemTable,
      "dXdZErr",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->dXdZErr; },
      "propagation error on dx/dz (0 if the variance is not positive)");
  addColumn<float>(
      *gemTable,
      "dYdZErr",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->dYdZErr; },
      "propagation error on dy/dz (0 if the variance is not positive)");
  addColumn<float>(
      *gemTable,
      "edgeX",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->edgeX; },
      "distance in x of the propagated position to the edge of the chamber's trapezoid, cm (< 0: inside)");
  addColumn<float>(
      *gemTable,
      "edgeY",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->edgeY; },
      "distance in y of the propagated position to the edge of the chamber's trapezoid, cm (< 0: inside)");
  addColumn<int16_t>(
      *gemTable,
      "roll",
      gemChambers,
      [](const GEMChamberRow& r) { return r.roll; },
      "eta partition (GEMDetId roll) whose bounds hold the point where the straight line through the propagated "
      "state (x, y, dXdZ, dYdZ) crosses the partition plane (off the chamber surface; reco's rechit border cut "
      "projects the surface point orthogonally instead); 0 if none (in a gap between partitions, outside the active "
      "area, dXdZ = 9999, or a super-chamber row)");
  addColumn<float>(
      *gemTable,
      "strip",
      gemChambers,
      [](const GEMChamberRow& r) { return r.strip; },
      "strip coordinate of the crossing point in that partition (GEMEtaPartition::strip, clamped to 0 to nstrips: "
      "exactly 0 or nstrips outside the strips); -1 if roll is 0");
  addColumn<float>(
      *gemTable,
      "pad",
      gemChambers,
      [](const GEMChamberRow& r) { return r.pad; },
      "pad coordinate of the crossing point in that partition (GEMEtaPartition::pad, the strip coordinate on the "
      "pad topology, clamped to 0 to npads); -1 if roll is 0");
  addColumn<float>(
      *gemTable, "eta", gemChambers, [](const GEMChamberRow& r) { return r.eta; }, "eta of the point (x, y)");
  addColumn<float>(
      *gemTable, "phi", gemChambers, [](const GEMChamberRow& r) { return r.phi; }, "phi of the point (x, y)");
  addColumn<int16_t>(
      *gemTable,
      "nHitMatches",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->gemHitMatches.size(); },
      "GEM rechits reco matched in this chamber (rows of MuonGEMHitMatch): the propagated point within "
      "2 cm of the rechit's eta partition and |dx| <= 5 cm or dx^2 <= 16 xErr^2");
  addColumn<int16_t>(
      *gemTable,
      "nSegmentMatches",
      gemChambers,
      [](const GEMChamberRow& r) { return r.match->gemMatches.size(); },
      "GEM segments reco matched in this chamber: super-chamber rows only (reco attaches GEM segments to "
      "super-chambers; its muon detector associator builds them for GE0 only since CMSSW_12_5_0, for every station "
      "before)");

  auto gemHitTable = std::make_unique<nanoaod::FlatTable>(gemHits.size(), "MuonGEMHitMatch", false);
  gemHitTable->setDoc(
      "GEM rechits reco matched to the GEM chamber matches (reco::MuonGEMHitMatch), grouped by MuonGEMChamber row. "
      "Reco keeps x and bx only (the mask is always 0): the rechit's eta partition and y are not stored, and rechits "
      "of every partition within 2 cm of the point are matched");
  addColumn<int16_t>(
      *gemHitTable,
      "muonGEMChamberIdx",
      gemHits,
      [](const GEMHitRow& h) { return h.chamberIdx; },
      "row in MuonGEMChamber");
  addColumn<float>(
      *gemHitTable,
      "x",
      gemHits,
      [](const GEMHitRow& h) { return h.hit->x; },
      "local x of the rechit, cm, as reco stores it (frame of the rechit's eta partition)");
  addColumn<float>(
      *gemHitTable,
      "dx",
      gemHits,
      [](const GEMHitRow& h) { return h.dx; },
      "x - MuonGEMChamber_x, cm: reco's matching variable (rechit partition frame minus chamber frame)");
  addColumn<int16_t>(
      *gemHitTable, "bx", gemHits, [](const GEMHitRow& h) { return h.hit->bx; }, "rechit BX (GEMRecHit::BunchX)");

  event.put(std::move(lctTable), "lct");
  event.put(std::move(muonTable), "muon");
  event.put(std::move(linkTable), "lctMuonMatch");
  event.put(std::move(chamberTable), "muonCSCChamber");
  event.put(std::move(gemTable), "muonGEMChamber");
  event.put(std::move(gemHitTable), "muonGEMHitMatch");
  event.put(std::move(segmentTable), "muonCSCSegment");
  event.put(std::move(layerTable), "muonCSCLayer");
  event.put(std::move(adcTable), "muonCSCRecHitADC");
}

DEFINE_FWK_MODULE(CSCLCTRecoMuonTableProducer);
