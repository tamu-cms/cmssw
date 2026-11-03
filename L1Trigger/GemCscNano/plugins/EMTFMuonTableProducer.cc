// The EMTF output tracks as the uGMT receives them: BXVector<l1t::RegionalMuonCand> (gmtStage2Digis:EMTF, the uGMT
// input readout, in data and MC), every BX of the readout, BX by BX in collection order. Per track: the hardware words
// as unpacked (RegionalMuonRawDigiTranslator), the EMTF track address, the uGMT input index of the central L1Mu table,
// and pt, eta, phi with the uGMT conversions (L1TMuonProducer.cc:329-337; MicroGMTConfiguration::calcGlobalPhi).

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
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
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using gemcscnano::addColumn;

class EMTFMuonTableProducer : public edm::global::EDProducer<> {
public:
  explicit EMTFMuonTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  const edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> token_;
  const std::string label_;
};

EMTFMuonTableProducer::EMTFMuonTableProducer(const edm::ParameterSet& ps)
    : token_(consumes(ps.getParameter<edm::InputTag>("src"))), label_(ps.getParameter<edm::InputTag>("src").encode()) {
  produces<nanoaod::FlatTable>();
}

void EMTFMuonTableProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const {
  const auto& cands = event.get(token_);
  struct Row {
    int bx;
    const l1t::RegionalMuonCand* mu;
    int endcap;
    int globalPhi;
  };
  std::vector<Row> rows;
  for (int bx = cands.getFirstBX(); bx <= cands.getLastBX(); ++bx)
    for (auto mu = cands.begin(bx); mu != cands.end(bx); ++mu) {
      const auto type = mu->trackFinderType();
      if (type != l1t::tftype::emtf_pos && type != l1t::tftype::emtf_neg)
        throw cms::Exception("InvalidInput") << "EMTFMuonTableProducer: a candidate of " << label_
                                             << " is not from EMTF (trackFinderType " << int(type) << ")";
      if (mu->muIdx() < 0)
        throw cms::Exception("InvalidInput")
            << "EMTFMuonTableProducer: a candidate of " << label_ << " has no muIdx (not an unpacked uGMT input)";
      rows.push_back({bx,
                      &*mu,
                      type == l1t::tftype::emtf_pos ? +1 : -1,
                      l1t::MicroGMTConfiguration::calcGlobalPhi(mu->hwPhi(), type, mu->processor())});
    }
  // the uGMT conversion to physical values (L1TMuonProducer.cc:329-330)
  auto p4 = [](const Row& r) {
    return math::PtEtaPhiMLorentzVector((r.mu->hwPt() - 1) * 0.5, r.mu->hwEta() * 0.010875, r.globalPhi * 0.010908, 0.);
  };
  using A = l1t::RegionalMuonCand;

  auto table = std::make_unique<nanoaod::FlatTable>(rows.size(), "EMTFMuon", false);
  table->setDoc("EMTF output tracks as received by the uGMT: every candidate of " + label_ +
                ", every BX, BX by BX in collection order");
  addColumn<int16_t>(*table, "bx", rows, [](const Row& r) { return r.bx; }, "BX of the L1 candidate");
  addColumn<int16_t>(
      *table, "endcap", rows, [](const Row& r) { return r.endcap; }, "+1 or -1 (trackFinderType emtf_pos, emtf_neg)");
  addColumn<int16_t>(
      *table,
      "processor",
      rows,
      [](const Row& r) { return r.mu->processor(); },
      "EMTF sector processor 0-5 (sector = processor + 1)");
  addColumn<int16_t>(
      *table, "muIdx", rows, [](const Row& r) { return r.mu->muIdx(); }, "position on its uGMT input link (0-2)");
  addColumn<uint16_t>(
      *table,
      "tfMuonIndex",
      rows,
      [](const Row& r) { return 3 * (r.mu->link() - 36) + r.mu->muIdx(); },
      "index at the uGMT input, 3 * (link - 36) + muIdx (EMTF+ 0-17, EMTF- 90-107): L1Mu_tfMuonIndex of a uGMT muon "
      "made from it, at the same bx");
  addColumn<int16_t>(
      *table,
      "hwPt",
      rows,
      [](const Row& r) { return r.mu->hwPt(); },
      "hardware pt (9 bits): pt = (hwPt - 1) * 0.5 GeV");
  addColumn<int16_t>(
      *table,
      "hwPtUnconstrained",
      rows,
      [](const Row& r) { return r.mu->hwPtUnconstrained(); },
      "hardware beam-spot-unconstrained pt (8 bits): hwPtUnconstrained - 1 GeV, 0 if 0 (0 when the uGMT firmware "
      "does not transport it)");
  addColumn<uint16_t>(
      *table,
      "hwDXY",
      rows,
      [](const Row& r) { return r.mu->hwDXY(); },
      "hardware impact-parameter class (2 bits; 0 when the uGMT firmware does not transport it)");
  addColumn<int16_t>(
      *table, "hwEta", rows, [](const Row& r) { return r.mu->hwEta(); }, "hardware eta: eta = hwEta * 0.010875");
  addColumn<int16_t>(
      *table,
      "hwPhi",
      rows,
      [](const Row& r) { return r.mu->hwPhi(); },
      "hardware local phi in the sector (LSB 2pi/576)");
  addColumn<int16_t>(
      *table,
      "hwGlobalPhi",
      rows,
      [](const Row& r) { return r.globalPhi; },
      "hardware global phi 0-575 (LSB 2pi/576), MicroGMTConfiguration::calcGlobalPhi; L1Mu_hwPhi of a uGMT muon made "
      "from it");
  addColumn<int16_t>(
      *table,
      "hwSign",
      rows,
      [](const Row& r) { return r.mu->hwSign(); },
      "charge sign bit: charge = (-1)^hwSign when hwSignValid");
  addColumn<int16_t>(*table, "hwSignValid", rows, [](const Row& r) { return r.mu->hwSignValid(); }, "charge valid bit");
  addColumn<int16_t>(
      *table, "hwQual", rows, [](const Row& r) { return r.mu->hwQual(); }, "hardware quality (4 bits, from the mode)");
  addColumn<bool>(*table, "hwHF", rows, [](const Row& r) { return r.mu->hwHF() != 0; }, "halo / HF bit");
  const std::vector<std::pair<std::string, A::emtfAddress>> address = {{"me1Seg", A::kME1Seg},
                                                                       {"me1Ch", A::kME1Ch},
                                                                       {"me2Seg", A::kME2Seg},
                                                                       {"me2Ch", A::kME2Ch},
                                                                       {"me3Seg", A::kME3Seg},
                                                                       {"me3Ch", A::kME3Ch},
                                                                       {"me4Seg", A::kME4Seg},
                                                                       {"me4Ch", A::kME4Ch}};
  for (const auto& [name, key] : address)
    addColumn<int16_t>(
        *table,
        name,
        rows,
        [key = key](const Row& r) { return r.mu->trackSubAddress(key); },
        "track address: " + name + " (trackSubAddress)");
  addColumn<uint64_t>(
      *table,
      "dataword",
      rows,
      [](const Row& r) { return r.mu->dataword(); },
      "the 64-bit uGMT input word as read out (bits the decoding zeroes stay here)");
  addColumn<float>(
      *table, "pt", rows, [&](const Row& r) { return float(p4(r).pt()); }, "pt, GeV: (hwPt - 1) * 0.5", 12);
  addColumn<float>(
      *table,
      "ptUnconstrained",
      rows,
      [](const Row& r) { return r.mu->hwPtUnconstrained() == 0 ? 0.f : float(r.mu->hwPtUnconstrained() - 1); },
      "beam-spot-unconstrained pt, GeV: hwPtUnconstrained - 1 (0 if 0)",
      12);
  addColumn<float>(*table, "eta", rows, [&](const Row& r) { return float(p4(r).eta()); }, "eta: hwEta * 0.010875", 12);
  addColumn<float>(
      *table,
      "phi",
      rows,
      [&](const Row& r) { return float(p4(r).phi()); },
      "global phi: hwGlobalPhi * 0.010908, in [-pi, pi]",
      12);
  event.put(std::move(table));
}

void EMTFMuonTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("BXVector<l1t::RegionalMuonCand> of EMTF candidates");
  descriptions.add("emtfMuonTableProducer", desc);
}

DEFINE_FWK_MODULE(EMTFMuonTableProducer);
