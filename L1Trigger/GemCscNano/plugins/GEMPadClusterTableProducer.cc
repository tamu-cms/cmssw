// The GE1/1 pad clusters the ME1/1 OTMBs received and read out: every GEMPadDigiCluster of the CSC unpacker
// (muonCSCDigis:MuonGEMPadDigiCluster, decoded from the GEM block of the OTMB DAQ data, CSCDCCUnpacker.cc:633-651,
// CSCGEMData.cc:148-186), in collection order (detector unit, then digi). The unpacker writes each cluster as a
// first pad and the following ones (consecutive), its time bin as bx, and the TMB header's ALCT match time; it keeps
// only first pads below 192 (CSCGEMData.cc:166), so every cluster is valid; a readout that fails the TMB-data size
// check keeps its TMB header but loses its clusters (CSCDCCUnpacker.cc:627, 652-653). MC RAW is packed from the
// DIGI step's clusters (CSCEventData.cc:541-551, CSCGEMData.cc:192-209): the time bin is the simulated BX (only 0-11
// are packed) and the ALCT match time is the packer's -1 in a 4-bit field (15). Co-pad positions are built inside the
// OTMB and not read out (the header's co-pad counts are in TMBHeader).

#include "DataFormats/GEMDigi/interface/GEMPadDigiClusterCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "L1Trigger/GemCscNano/interface/TableColumns.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using gemcscnano::addColumn;

class GEMPadClusterTableProducer : public edm::global::EDProducer<> {
public:
  explicit GEMPadClusterTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  const edm::EDGetTokenT<GEMPadDigiClusterCollection> token_;
  const std::string label_;
};

GEMPadClusterTableProducer::GEMPadClusterTableProducer(const edm::ParameterSet& ps)
    : token_(consumes(ps.getParameter<edm::InputTag>("src"))), label_(ps.getParameter<edm::InputTag>("src").encode()) {
  produces<nanoaod::FlatTable>();
}

void GEMPadClusterTableProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const {
  const auto& clusters = event.get(token_);
  struct Row {
    GEMDetId id;
    const GEMPadDigiCluster* cluster;
  };
  std::vector<Row> rows;
  for (auto detUnit = clusters.begin(); detUnit != clusters.end(); ++detUnit) {
    const GEMDetId& id = (*detUnit).first;
    for (auto cluster = (*detUnit).second.first; cluster != (*detUnit).second.second; ++cluster) {
      // the columns hold the first pad and the size: the pads must be consecutive (as the unpacker builds them)
      const auto& pads = cluster->pads();
      if (pads.empty())
        throw cms::Exception("InvalidInput")
            << "GEMPadClusterTableProducer: a cluster of " << label_ << " in " << id << " has no pad";
      for (size_t i = 1; i < pads.size(); ++i)
        if (pads[i] != pads[0] + i)
          throw cms::Exception("InvalidInput")
              << "GEMPadClusterTableProducer: a cluster of " << label_ << " in " << id << " has non-consecutive pads";
      rows.push_back({id, &*cluster});
    }
  }

  auto table = std::make_unique<nanoaod::FlatTable>(rows.size(), "GEMPadCluster", false);
  table->setDoc("GE1/1 pad clusters read out by the ME1/1 OTMBs: every GEMPadDigiCluster of " + label_ +
                ", in collection order (detector unit, then digi); MC: the clusters the packer wrote into the RAW "
                "(simulated BX 0-11); co-pad positions are not read out");
  addColumn<int16_t>(*table, "region", rows, [](const Row& r) { return r.id.region(); }, "+1 or -1 (the CSC endcap)");
  addColumn<int16_t>(*table, "station", rows, [](const Row& r) { return r.id.station(); }, "GEM station");
  addColumn<int16_t>(*table, "ring", rows, [](const Row& r) { return r.id.ring(); }, "GEM ring");
  addColumn<int16_t>(
      *table,
      "chamber",
      rows,
      [](const Row& r) { return r.id.chamber(); },
      "chamber: the number of the ME1/1 chamber whose OTMB read it out");
  addColumn<int16_t>(
      *table, "layer", rows, [](const Row& r) { return r.id.layer(); }, "GEM layer 1 or 2 (the readout's GEM A or B)");
  addColumn<int16_t>(
      *table,
      "roll",
      rows,
      [](const Row& r) { return r.id.roll(); },
      "eta partition as GEMDetId numbers it (8 = closest to the beam pipe): the unpacker maps the readout eta e (0 = "
      "closest to the beam pipe) to roll 8 - e");
  addColumn<int16_t>(
      *table, "firstPad", rows, [](const Row& r) { return r.cluster->pads().front(); }, "first pad of the cluster");
  addColumn<int16_t>(
      *table,
      "nPads",
      rows,
      [](const Row& r) { return r.cluster->pads().size(); },
      "number of pads (consecutive from firstPad; the readout encodes up to 8)");
  addColumn<int16_t>(
      *table,
      "bx",
      rows,
      [](const Row& r) { return r.cluster->bx(); },
      "data: the readout time bin of the cluster (0 = the first time bin of the OTMB's GEM readout), not an "
      "L1A-relative BX; MC: the simulated BX the packer used as time bin (0 = in time; only 0-11 are packed)");
  addColumn<int16_t>(
      *table,
      "alctMatchTime",
      rows,
      [](const Row& r) { return r.cluster->alctMatchTime(); },
      "data: the TMB header's ALCT match time of the readout, as the unpacker attaches it; MC: 15, the packer's "
      "-1 (no match time) in the 4-bit header field");
  event.put(std::move(table));
}

void GEMPadClusterTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("GEMPadDigiClusterCollection (muonCSCDigis:MuonGEMPadDigiCluster)");
  descriptions.add("gemPadClusterTableProducer", desc);
}

DEFINE_FWK_MODULE(GEMPadClusterTableProducer);
