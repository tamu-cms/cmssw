// NanoAOD tables of CSC LCTs, reco muons and their CSCLCTRecoMuonMatcher associations (one entry per event):
//   LCT           every LCT the matcher considers (valid, inside lctBXWindow), in collection order;
//   Muon          the reco muons passing muonCut (default pT > 10 GeV and the medium or tight cut-based ID), in
//                 collection order; only these muons are matched;
//   LCTMuonMatch  one row per matched LCT-muon pair (rows of the LCT and Muon tables, segment, key-layer residuals),
//                 grouped by LCT, closest muon first.
// No event selection is applied: tag and probe is left to the analysis of the ntuple.

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/MuonReco/interface/Muon.h"
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
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "L1Trigger/GemCscNano/interface/CSCLCTRecoMuonMatcher.h"

#include <algorithm>
#include <limits>
#include <memory>
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
  const StringCutObjectSelector<reco::Muon> muonCut_;
  const CSCLCTRecoMuonMatcher matcher_;
};

CSCLCTRecoMuonTableProducer::CSCLCTRecoMuonTableProducer(const edm::ParameterSet& ps)
    : lctToken_(consumes(ps.getParameter<edm::InputTag>("lcts"))),
      muonToken_(consumes(ps.getParameter<edm::InputTag>("muons"))),
      vertexToken_(consumes(ps.getParameter<edm::InputTag>("vertices"))),
      cscGeomToken_(esConsumes()),
      muonCut_(ps.getParameter<std::string>("muonCut")),
      matcher_(ps.getParameter<edm::ParameterSet>("matcher")) {
  produces<nanoaod::FlatTable>("lct");
  produces<nanoaod::FlatTable>("muon");
  produces<nanoaod::FlatTable>("lctMuonMatch");
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

  // Muon rows are the muon collection itself, so a row index is the muon's index in the collection.
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

  event.put(std::move(lctTable), "lct");
  event.put(std::move(muonTable), "muon");
  event.put(std::move(linkTable), "lctMuonMatch");
}

DEFINE_FWK_MODULE(CSCLCTRecoMuonTableProducer);
