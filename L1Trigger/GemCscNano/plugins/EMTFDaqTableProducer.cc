// Data only: the EMTF DAQ readout as the EMTF unpacker (emtfStage2Digis, EventFilter/L1TRawToDigi EMTFSetup) decodes
// it. MC RAW has no EMTF readout (L1Trigger/Configuration/python/L1TDigiToRaw_cff.py: "Missing: muon EMTF").
//   EMTFHit    every hit of the unpacker's EMTFHitCollection, in unpacking order: the CSC LCTs, GEM clusters and RPC
//              hits each sector received that the unpacker keeps (valid, BX within +-3, no exact duplicate: a GEM
//              cluster at the roll and pad of an earlier one counts as one, whatever its layer or size; RPC/GEM hits
//              only in sectors with a CSC hit, EMTFCollections.cc:12-26; neighbor copies included), in EMTF conventions;
//   EMTFTrack  every track of the EMTFTrackCollection, with the row in EMTFHit of its hit in each station (the unpacker
//              fills no hit index, EMTFTrack::HitIdx(): the rows are found by the hit's identity fields);
//   EMTFLCT    the CSCCorrelatedLCTDigis the unpacker builds from the CSC hits: one per distinct LCT, from its
//              first copy in unpacking order, home or neighbor sector (EMTFBlockME.cc:315-357), with the columns of the
//              LCT table (gemcscnano::makeLCTTable), and its match to the LCT table (below);
//   LCT        (extension of the LCT table) the match seen from the LCT rows.
// LCT <-> EMTFLCT: both are the LCT word the (O)TMB sent, decoded by two unpackers (CSC DAQ readout, EMTF input
// readout), so rows are matched by exact equality of their stored identity columns endcap, station, ring, chamber,
// es (the key eighth-strip, which contains the half-strip) and wg; no tolerance. Rows without a partner or with
// several are flagged, never dropped: nLCTMatch / nEMTFLCTMatch count the partners, the index is set only for a
// one-to-one pair, lctDiff lists the compared fields that differ in it.

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include "DataFormats/L1TMuon/interface/L1TMuonSubsystems.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "L1Trigger/GemCscNano/interface/CSCLCTTable.h"
#include "L1Trigger/GemCscNano/interface/TableColumns.h"

#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

using gemcscnano::addColumn;

namespace {

  // The identity fields of a hit, as the unpacker's duplicate checks use them (EMTFBlockME.cc:315-324,
  // EMTFBlockRPC.cc:253-263, EMTFBlockGEM.cc:251-261), plus the subsystem and the receiving sector.
  bool sameHit(const l1t::EMTFHit& a, const l1t::EMTFHit& b) {
    if (a.Subsystem() != b.Subsystem() || a.BX() != b.BX() || a.Endcap() != b.Endcap() || a.Station() != b.Station() ||
        a.Chamber() != b.Chamber() || a.Ring() != b.Ring() || a.Neighbor() != b.Neighbor() ||
        a.Sector_idx() != b.Sector_idx())
      return false;
    switch (a.Subsystem()) {
      case L1TMuon::kCSC:
        return a.Strip() == b.Strip() && a.Wire() == b.Wire();
      case L1TMuon::kRPC:
        return a.Theta_fp() == b.Theta_fp() && a.Phi_fp() == b.Phi_fp();
      case L1TMuon::kGEM:
        return a.Roll() == b.Roll() && a.Pad() == b.Pad();
      default:
        throw cms::Exception("InvalidInput") << "EMTFDaqTableProducer: hit of subsystem " << a.Subsystem();
    }
  }

  // a column of a table made by gemcscnano::makeLCTTable, as int (the identity and compared fields)
  std::vector<int> lctColumn(const nanoaod::FlatTable& table, const std::string& name) {
    const int i = table.columnIndex(name);
    if (i < 0)
      throw cms::Exception("LogicError") << "EMTFDaqTableProducer: table " << table.name() << " has no column " << name;
    std::vector<int> v;
    switch (table.columnType(i)) {
      case nanoaod::FlatTable::ColumnType::Int16:
        for (auto x : table.columnData<int16_t>(i))
          v.push_back(x);
        break;
      case nanoaod::FlatTable::ColumnType::Bool:
        for (auto x : table.columnData<bool>(i))
          v.push_back(x);
        break;
      default:
        throw cms::Exception("LogicError")
            << "EMTFDaqTableProducer: column " << name << " of " << table.name() << " is not int16 or bool";
    }
    return v;
  }

  constexpr std::array<const char*, 6> kIdentity = {{"endcap", "station", "ring", "chamber", "es", "wg"}};
  // compared fields of a one-to-one pair and their lctDiff bits
  constexpr std::array<std::pair<const char*, int>, 5> kCompared = {
      {{"quality", 1}, {"pattern", 2}, {"bend", 4}, {"slope", 8}, {"valid", 16}}};
  constexpr int kSlopeBit = 8;

}  // namespace

class EMTFDaqTableProducer : public edm::global::EDProducer<> {
public:
  explicit EMTFDaqTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  const edm::EDGetTokenT<l1t::EMTFHitCollection> hitToken_;
  const edm::EDGetTokenT<l1t::EMTFTrackCollection> trackToken_;
  const edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lctToken_;
  const edm::EDGetTokenT<nanoaod::FlatTable> lctTableToken_;
  const std::string label_;
};

EMTFDaqTableProducer::EMTFDaqTableProducer(const edm::ParameterSet& ps)
    : hitToken_(consumes(ps.getParameter<edm::InputTag>("src"))),
      trackToken_(consumes(ps.getParameter<edm::InputTag>("src"))),
      lctToken_(consumes(ps.getParameter<edm::InputTag>("src"))),
      lctTableToken_(consumes(ps.getParameter<edm::InputTag>("lctTable"))),
      label_(ps.getParameter<edm::InputTag>("src").encode()) {
  produces<nanoaod::FlatTable>("hit");
  produces<nanoaod::FlatTable>("track");
  produces<nanoaod::FlatTable>("lct");
  produces<nanoaod::FlatTable>("lctMatch");
}

void EMTFDaqTableProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const {
  const auto& hits = event.get(hitToken_);
  const auto& tracks = event.get(trackToken_);
  using Hit = l1t::EMTFHit;
  using Track = l1t::EMTFTrack;

  auto hitTable = std::make_unique<nanoaod::FlatTable>(hits.size(), "EMTFHit", false);
  hitTable->setDoc(
      "EMTF DAQ hits (data only): every hit of " + label_ +
      " (the CSC LCTs, GEM clusters, RPC hits received per sector that the unpacker keeps: valid, BX within "
      "+-3, no exact duplicate, RPC/GEM only in sectors with a CSC hit; neighbor copies included), in "
      "unpacking order; EMTF conventions");
  addColumn<int16_t>(
      *hitTable, "subsystem", hits, [](const Hit& h) { return h.Subsystem(); }, "1 CSC, 2 RPC, 3 GEM (L1TMuon)");
  addColumn<int16_t>(*hitTable, "endcap", hits, [](const Hit& h) { return h.Endcap(); }, "+1 or -1");
  addColumn<int16_t>(
      *hitTable,
      "sectorIdx",
      hits,
      [](const Hit& h) { return h.Sector_idx(); },
      "receiving sector index (0-5 ME+, 6-11 ME-)");
  addColumn<int16_t>(*hitTable, "sector", hits, [](const Hit& h) { return h.Sector(); }, "sector 1-6 of the hit");
  addColumn<bool>(
      *hitTable,
      "neighbor",
      hits,
      [](const Hit& h) { return h.Neighbor() == 1; },
      "received from the neighbor sector (a shared CSC LCT appears twice: neighbor 0 and 1)");
  addColumn<int16_t>(*hitTable, "station", hits, [](const Hit& h) { return h.Station(); }, "station 1-4 (GEM: 1)");
  addColumn<int16_t>(
      *hitTable,
      "ring",
      hits,
      [](const Hit& h) { return h.Ring(); },
      "ring, EMTF convention: CSC ME1/1a = 4 (the LCT table's ring 1), RPC 2-3, GEM 1");
  addColumn<int16_t>(
      *hitTable, "chamber", hits, [](const Hit& h) { return h.Chamber(); }, "chamber (CSCDetId numbering)");
  addColumn<int16_t>(
      *hitTable, "subsector", hits, [](const Hit& h) { return h.Subsector(); }, "subsector (ME2-4: -1 when unpacked)");
  addColumn<int16_t>(*hitTable, "cscId", hits, [](const Hit& h) { return h.CSC_ID(); }, "CSC ID in the sector (CSC)");
  addColumn<int16_t>(
      *hitTable,
      "bx",
      hits,
      [](const Hit& h) { return h.BX(); },
      "BX relative to the L1A: readout time bin - 3 (EMTF convention; not the LCT table's bx)");
  addColumn<int16_t>(*hitTable, "stubNum", hits, [](const Hit& h) { return h.Stub_num(); }, "stub number");
  addColumn<int16_t>(
      *hitTable,
      "strip",
      hits,
      [](const Hit& h) { return h.Strip(); },
      "CSC: key half-strip, 0-based (ME1/1 unified, ME1a >= 128) = the LCT table's hs; GEM: first pad of the cluster");
  addColumn<int16_t>(
      *hitTable, "quartStripBit", hits, [](const Hit& h) { return h.Strip_quart_bit(); }, "CSC quarter-strip bit");
  addColumn<int16_t>(
      *hitTable, "eighthStripBit", hits, [](const Hit& h) { return h.Strip_eighth_bit(); }, "CSC eighth-strip bit");
  addColumn<int16_t>(
      *hitTable,
      "wire",
      hits,
      [](const Hit& h) { return h.Wire(); },
      "CSC: key wire group, 0-based = the LCT table's wg");
  addColumn<int16_t>(
      *hitTable, "quality", hits, [](const Hit& h) { return h.Quality(); }, "CSC: LCT quality; GEM: cluster size");
  addColumn<int16_t>(*hitTable, "pattern", hits, [](const Hit& h) { return h.Pattern(); }, "CSC: Run-2 pattern");
  addColumn<int16_t>(*hitTable, "bend", hits, [](const Hit& h) { return h.Bend(); }, "CSC: bend (unpacked: +-1)");
  addColumn<int16_t>(
      *hitTable, "slope", hits, [](const Hit& h) { return h.Slope(); }, "CSC: slope (4 bits; OTMB chambers)");
  addColumn<int16_t>(*hitTable, "layer", hits, [](const Hit& h) { return h.Layer(); }, "GEM: layer 1-2");
  addColumn<int16_t>(*hitTable, "roll", hits, [](const Hit& h) { return h.Roll(); }, "GEM: eta partition");
  addColumn<int16_t>(*hitTable, "clusterId", hits, [](const Hit& h) { return h.ClusterID(); }, "GEM: cluster ID");
  addColumn<int16_t>(*hitTable, "phiFp", hits, [](const Hit& h) { return h.Phi_fp(); }, "RPC: fixed-point phi");
  addColumn<int16_t>(*hitTable, "thetaFp", hits, [](const Hit& h) { return h.Theta_fp(); }, "RPC: fixed-point theta");

  // each track's hit in each station: the row of the identical hit of the collection; -1 no hit in the station, -2
  // the hit is not in the collection
  std::vector<std::array<int, 4>> hitIdx;
  hitIdx.reserve(tracks.size());
  for (const Track& t : tracks) {
    std::array<int, 4> idx = {-1, -1, -1, -1};
    for (const Hit& th : t.Hits()) {
      const int st = th.Station();
      if (st < 1 || st > 4)
        throw cms::Exception("InvalidInput") << "EMTFDaqTableProducer: track hit in station " << st;
      if (idx[st - 1] != -1)
        throw cms::Exception("InvalidInput") << "EMTFDaqTableProducer: two track hits in station " << st;
      // -2: the hit is not in the collection (the unpacker's zero suppression of RPC/GEM hits, EMTFCollections.cc:12-26)
      int found = -2;
      for (size_t j = 0; j < hits.size(); ++j)
        if (sameHit(th, hits[j])) {
          if (found != -2)
            throw cms::Exception("InvalidInput") << "EMTFDaqTableProducer: a track hit matches two hits";
          found = j;
        }
      idx[st - 1] = found;
    }
    hitIdx.push_back(idx);
  }
  struct TrackRow {
    const Track* t;
    std::array<int, 4> idx;
  };
  std::vector<TrackRow> trackRows;
  trackRows.reserve(tracks.size());
  for (size_t i = 0; i < tracks.size(); ++i)
    trackRows.push_back({&tracks[i], hitIdx[i]});

  auto trackTable = std::make_unique<nanoaod::FlatTable>(trackRows.size(), "EMTFTrack", false);
  trackTable->setDoc("EMTF DAQ tracks (data only): every track of " + label_ + ", in unpacking order");
  using R = TrackRow;
  addColumn<int16_t>(*trackTable, "endcap", trackRows, [](const R& r) { return r.t->Endcap(); }, "+1 or -1");
  addColumn<int16_t>(*trackTable, "sector", trackRows, [](const R& r) { return r.t->Sector(); }, "sector 1-6");
  addColumn<int16_t>(
      *trackTable, "sectorIdx", trackRows, [](const R& r) { return r.t->Sector_idx(); }, "sector index (0-11)");
  addColumn<int16_t>(
      *trackTable, "bx", trackRows, [](const R& r) { return r.t->BX(); }, "BX relative to the L1A (EMTF convention)");
  addColumn<int16_t>(
      *trackTable, "mode", trackRows, [](const R& r) { return r.t->Mode(); }, "station mode (4 bits, ME1 = 8)");
  addColumn<int16_t>(
      *trackTable, "modeCSC", trackRows, [](const R& r) { return r.t->Mode_CSC(); }, "mode of the CSC hits");
  addColumn<int16_t>(
      *trackTable, "modeRPC", trackRows, [](const R& r) { return r.t->Mode_RPC(); }, "mode of the RPC hits");
  addColumn<int16_t>(
      *trackTable, "modeGEM", trackRows, [](const R& r) { return r.t->Mode_GEM(); }, "mode of the GEM hits");
  addColumn<int16_t>(
      *trackTable,
      "modeNeighbor",
      trackRows,
      [](const R& r) { return r.t->Mode_neighbor(); },
      "mode of the hits from the neighbor sector");
  addColumn<int16_t>(
      *trackTable,
      "hwPt",
      trackRows,
      [](const R& r) { return r.t->GMT_pt(); },
      "hardware pt sent to the uGMT: pt = (hwPt - 1) * 0.5 GeV");
  addColumn<int16_t>(
      *trackTable,
      "hwPtUnconstrained",
      trackRows,
      [](const R& r) { return r.t->GMT_pt_dxy(); },
      "hardware beam-spot-unconstrained pt");
  addColumn<int16_t>(
      *trackTable, "hwDXY", trackRows, [](const R& r) { return r.t->GMT_dxy(); }, "hardware impact-parameter class");
  addColumn<int16_t>(
      *trackTable,
      "hwEta",
      trackRows,
      [](const R& r) { return r.t->GMT_eta(); },
      "hardware eta: eta = hwEta * 0.010875");
  addColumn<int16_t>(
      *trackTable, "hwPhi", trackRows, [](const R& r) { return r.t->GMT_phi(); }, "hardware local phi (LSB 2pi/576)");
  addColumn<int16_t>(
      *trackTable, "hwQual", trackRows, [](const R& r) { return r.t->GMT_quality(); }, "hardware quality");
  addColumn<int16_t>(*trackTable, "hwSign", trackRows, [](const R& r) { return r.t->GMT_charge(); }, "charge sign bit");
  addColumn<int16_t>(
      *trackTable, "hwSignValid", trackRows, [](const R& r) { return r.t->GMT_charge_valid(); }, "charge valid bit");
  addColumn<int16_t>(
      *trackTable, "phiFp", trackRows, [](const R& r) { return r.t->Phi_fp(); }, "fixed-point phi in the sector");
  addColumn<int16_t>(
      *trackTable,
      "trackNum",
      trackRows,
      [](const R& r) { return r.t->Track_num(); },
      "track number (counts the sector's SP records, also of earlier BXs)");
  addColumn<int16_t>(
      *trackTable,
      "nHits",
      trackRows,
      [](const R& r) { return r.t->NumHits(); },
      "number of hits the unpacker matched (-99: none, EMTFTrack.h:85)");
  for (int st = 1; st <= 4; ++st)
    addColumn<int16_t>(
        *trackTable,
        "hitIdx" + std::to_string(st),
        trackRows,
        [st](const R& r) { return r.idx[st - 1]; },
        "row in EMTFHit of the station-" + std::to_string(st) +
            " hit; -1: no hit in the station; -2: the hit is "
            "not in EMTFHit (removed by the unpacker's zero suppression)");

  const auto& emtfLcts = event.get(lctToken_);
  auto lctTable = gemcscnano::makeLCTTable(
      emtfLcts,
      "EMTFLCT",
      "LCTs of the EMTF readout (data only): every CSCCorrelatedLCTDigi of " + label_ +
          ", built by the EMTF unpacker from the CSC hits: one per distinct LCT (its first copy, home or neighbor "
          "sector), in collection order; matched to LCT by exact equality of endcap, station, ring, chamber, es, wg",
      {{"bx",
        "BX: the EMTF hit bx + 8 (emulator convention, CSCConstants::LCT_CENTRAL_BX); not comparable with LCT_bx"},
       {"ring", "ring as CSCDetId (ME1/1a and ME1b: 1)"},
       {"trknmb", "track number: always 1 (set by the EMTF unpacker)"},
       {"run3Pattern", "Run-3 pattern as the EMTF unpacker sets it (a partial raw field)"},
       {"gemLayer", "not set by the EMTF unpacker (EMTFHit::CreateCSCCorrelatedLCTDigi)"}});

  // LCT <-> EMTFLCT on the stored values: the LCT table's rows, whatever collection it is built from
  const auto& lctRows = event.get(lctTableToken_);
  if (lctRows.name() != "LCT" || lctRows.extension() || lctRows.singleton())
    throw cms::Exception("Configuration")
        << "EMTFDaqTableProducer: lctTable is the table " << lctRows.name() << ", not the main LCT table";
  const unsigned nL = lctRows.nRows(), nE = lctTable->nRows();
  if (nL > unsigned(std::numeric_limits<int16_t>::max()) || nE > unsigned(std::numeric_limits<int16_t>::max()))
    throw cms::Exception("InvalidInput") << "EMTFDaqTableProducer: " << nL << " LCT or " << nE
                                         << " EMTFLCT rows do not fit an int16 index";
  // the EMTF unpacker sets the slope of the Run-3 LCTs only (OTMB chambers, Run-3 DAQ format, EMTFBlockME.cc:253-266);
  // its other LCTs keep EMTFHit's default (-99): slope is compared only for Run-3 EMTF LCTs
  std::vector<bool> emtfRun3;
  for (auto detUnit = emtfLcts.begin(); detUnit != emtfLcts.end(); ++detUnit)
    for (auto lct = (*detUnit).second.first; lct != (*detUnit).second.second; ++lct)
      emtfRun3.push_back(lct->isRun3());
  if (emtfRun3.size() != nE)
    throw cms::Exception("LogicError") << "EMTFDaqTableProducer: " << emtfRun3.size() << " EMTF LCT digis but " << nE
                                       << " EMTFLCT rows";
  using Key = std::array<int, kIdentity.size()>;
  auto keys = [](const nanoaod::FlatTable& t) {
    std::vector<Key> k(t.nRows());
    for (unsigned f = 0; f < kIdentity.size(); ++f) {
      const auto c = lctColumn(t, kIdentity[f]);
      for (unsigned r = 0; r < t.nRows(); ++r)
        k[r][f] = c[r];
    }
    return k;
  };
  const auto keyL = keys(lctRows), keyE = keys(*lctTable);
  std::map<Key, std::vector<int>> rowsL, rowsE;
  for (unsigned i = 0; i < nL; ++i)
    rowsL[keyL[i]].push_back(i);
  for (unsigned j = 0; j < nE; ++j)
    rowsE[keyE[j]].push_back(j);
  auto partners = [](const std::map<Key, std::vector<int>>& rows, const Key& k) -> const std::vector<int>* {
    const auto it = rows.find(k);
    return it == rows.end() ? nullptr : &it->second;
  };
  std::vector<std::vector<int>> comparedL, comparedE;
  for (const auto& [name, bit] : kCompared) {
    comparedL.push_back(lctColumn(lctRows, name));
    comparedE.push_back(lctColumn(*lctTable, name));
  }
  std::vector<int16_t> nLCTMatch(nE, 0), lctIdx(nE, -1), lctDiff(nE, -1), nEMTFLCTMatch(nL, 0), emtfLctIdx(nL, -1);
  for (unsigned i = 0; i < nL; ++i)
    if (const auto* p = partners(rowsE, keyL[i]))
      nEMTFLCTMatch[i] = p->size();
  for (unsigned j = 0; j < nE; ++j) {
    const auto* p = partners(rowsL, keyE[j]);
    if (!p)
      continue;
    nLCTMatch[j] = p->size();
    if (p->size() != 1 || rowsE.at(keyE[j]).size() != 1)
      continue;
    const int i = p->front();
    lctIdx[j] = i;
    emtfLctIdx[i] = j;
    int diff = 0;
    for (unsigned f = 0; f < kCompared.size(); ++f)
      if (comparedL[f][i] != comparedE[f][j] && (kCompared[f].second != kSlopeBit || emtfRun3[j]))
        diff |= kCompared[f].second;
    lctDiff[j] = diff;
  }
  lctTable->addColumn<int16_t>(
      "nLCTMatch",
      nLCTMatch,
      "number of LCT rows with the same endcap, station, ring, chamber, es and wg (exact equality): 0 = no LCT "
      "(e.g. another BX or a chamber the CSC DAQ did not read out), 2 or more = ambiguous");
  lctTable->addColumn<int16_t>(
      "lctIdx",
      lctIdx,
      "row in LCT of the one-to-one partner (nLCTMatch = 1 and that LCT's nEMTFLCTMatch = 1), else -1");
  lctTable->addColumn<int16_t>(
      "lctDiff",
      lctDiff,
      "one-to-one pairs: bits of the compared fields that differ, 1 quality, 2 pattern, 4 bend, 8 slope (only for "
      "Run-3 EMTF LCTs, the OTMB chambers: the EMTF unpacker sets no other slope), 16 valid; 0 = identical; -1 = no "
      "one-to-one partner. Not compared: bx (1 bit in LCT, EMTF bx + 8 here), run3Pattern (the EMTF unpacker keeps a "
      "partial field), trknmb, type, gemLayer, slope6");
  auto lctMatchTable = std::make_unique<nanoaod::FlatTable>(nL, "LCT", false, true);
  lctMatchTable->addColumn<int16_t>(
      "nEMTFLCTMatch",
      nEMTFLCTMatch,
      "data: number of EMTFLCT rows with the same endcap, station, ring, chamber, es and wg (exact equality): 0 = "
      "none, 2 or more = ambiguous");
  lctMatchTable->addColumn<int16_t>(
      "emtfLctIdx",
      emtfLctIdx,
      "data: row in EMTFLCT of the one-to-one partner (else -1); the field differences: EMTFLCT_lctDiff");

  event.put(std::move(hitTable), "hit");
  event.put(std::move(trackTable), "track");
  event.put(std::move(lctTable), "lct");
  event.put(std::move(lctMatchTable), "lctMatch");
}

void EMTFDaqTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("the EMTF unpacker (EMTFHitCollection, EMTFTrackCollection, LCT digis)");
  desc.add<edm::InputTag>("lctTable", edm::InputTag("cscLCTRecoMuonTable", "lct"))
      ->setComment("the LCT FlatTable (main table named LCT) the EMTFLCT rows are matched to");
  descriptions.add("emtfDaqTableProducer", desc);
}

DEFINE_FWK_MODULE(EMTFDaqTableProducer);
