// Data only: the EMTF DAQ readout as the EMTF unpacker (emtfStage2Digis, EventFilter/L1TRawToDigi EMTFSetup) decodes
// it. MC RAW has no EMTF readout (L1Trigger/Configuration/python/L1TDigiToRaw_cff.py: "Missing: muon EMTF").
//   EMTFHit    every hit of the unpacker's EMTFHitCollection, in unpacking order: the CSC LCTs, GEM clusters and RPC
//              hits each sector received (neighbor copies included; RPC/GEM hits only in sectors with a CSC hit,
//              EMTFCollections.cc:88-102), in EMTF conventions;
//   EMTFTrack  every track of the EMTFTrackCollection, with the row in EMTFHit of its hit in each station (the unpacker
//              fills no hit index, EMTFTrack::HitIdx(): the rows are found by the hit's identity fields);
//   EMTFLCT    the CSCCorrelatedLCTDigis the unpacker builds from the CSC hits (neighbor copies dropped), with the
//              columns of the LCT table (gemcscnano::makeLCTTable).

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
#include <memory>
#include <string>
#include <vector>

using gemcscnano::addColumn;

namespace {

  // The identity fields of a hit, as the unpacker's duplicate checks use them (EMTFBlockME.cc:315-324,
  // EMTFBlockRPC.cc:253-263, EMTFBlockGEM.cc:251-261), plus the subsystem and the receiving sector.
  bool sameHit(const l1t::EMTFHit& a, const l1t::EMTFHit& b) {
    if (a.Subsystem() != b.Subsystem() || a.BX() != b.BX() || a.Endcap() != b.Endcap() ||
        a.Station() != b.Station() || a.Chamber() != b.Chamber() || a.Ring() != b.Ring() ||
        a.Neighbor() != b.Neighbor() || a.Sector_idx() != b.Sector_idx())
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
  const std::string label_;
};

EMTFDaqTableProducer::EMTFDaqTableProducer(const edm::ParameterSet& ps)
    : hitToken_(consumes(ps.getParameter<edm::InputTag>("src"))),
      trackToken_(consumes(ps.getParameter<edm::InputTag>("src"))),
      lctToken_(consumes(ps.getParameter<edm::InputTag>("src"))),
      label_(ps.getParameter<edm::InputTag>("src").encode()) {
  produces<nanoaod::FlatTable>("hit");
  produces<nanoaod::FlatTable>("track");
  produces<nanoaod::FlatTable>("lct");
}

void EMTFDaqTableProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const {
  const auto& hits = event.get(hitToken_);
  const auto& tracks = event.get(trackToken_);
  using Hit = l1t::EMTFHit;
  using Track = l1t::EMTFTrack;

  auto hitTable = std::make_unique<nanoaod::FlatTable>(hits.size(), "EMTFHit", false);
  hitTable->setDoc("EMTF DAQ hits (data only): every hit of " + label_ +
                   " (CSC LCTs, GEM clusters, RPC hits received per sector, neighbor copies included), in unpacking "
                   "order; EMTF conventions");
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
  addColumn<int16_t>(*hitTable, "chamber", hits, [](const Hit& h) { return h.Chamber(); }, "chamber (CSCDetId numbering)");
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

  // each track's hit in each station: the row of the identical hit of the collection, -1 if none
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
      int found = -1;
      for (size_t j = 0; j < hits.size(); ++j)
        if (sameHit(th, hits[j])) {
          if (found != -1)
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
      *trackTable,
      "bx",
      trackRows,
      [](const R& r) { return r.t->BX(); },
      "BX relative to the L1A (EMTF convention)");
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
      *trackTable, "hwEta", trackRows, [](const R& r) { return r.t->GMT_eta(); }, "hardware eta: eta = hwEta * 0.010875");
  addColumn<int16_t>(
      *trackTable, "hwPhi", trackRows, [](const R& r) { return r.t->GMT_phi(); }, "hardware local phi (LSB 2pi/576)");
  addColumn<int16_t>(*trackTable, "hwQual", trackRows, [](const R& r) { return r.t->GMT_quality(); }, "hardware quality");
  addColumn<int16_t>(
      *trackTable, "hwSign", trackRows, [](const R& r) { return r.t->GMT_charge(); }, "charge sign bit");
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
  addColumn<int16_t>(*trackTable, "nHits", trackRows, [](const R& r) { return r.t->NumHits(); }, "number of hits");
  for (int st = 1; st <= 4; ++st)
    addColumn<int16_t>(
        *trackTable,
        "hitIdx" + std::to_string(st),
        trackRows,
        [st](const R& r) { return r.idx[st - 1]; },
        "row in EMTFHit of the station-" + std::to_string(st) + " hit, -1 if none");

  auto lctTable = gemcscnano::makeLCTTable(
      event.get(lctToken_),
      "EMTFLCT",
      "LCTs of the EMTF readout (data only): every CSCCorrelatedLCTDigi of " + label_ +
          ", built by the EMTF unpacker from the CSC hits (neighbor copies dropped), in collection order",
      {{"bx",
        "BX: the EMTF hit bx + 8 (emulator convention, CSCConstants::LCT_CENTRAL_BX); not comparable with LCT_bx"},
       {"ring", "ring as CSCDetId (ME1/1a and ME1b: 1)"},
       {"trknmb", "track number: always 1 (set by the EMTF unpacker)"},
       {"run3Pattern", "Run-3 pattern as the EMTF unpacker sets it (a partial raw field)"},
       {"gemLayer", "not set by the EMTF unpacker (EMTFHit::CreateCSCCorrelatedLCTDigi)"}});

  event.put(std::move(hitTable), "hit");
  event.put(std::move(trackTable), "track");
  event.put(std::move(lctTable), "lct");
}

void EMTFDaqTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("the EMTF unpacker (EMTFHitCollection, EMTFTrackCollection, LCT digis)");
  descriptions.add("emtfDaqTableProducer", desc);
}

DEFINE_FWK_MODULE(EMTFDaqTableProducer);
