// The central NanoAOD PV columns that RECO input allows (PhysicsTools/NanoAOD/plugins/VertexTableProducer.cc without
// the miniAOD-only score and PF-candidate sums): the main primary vertex (the first of the collection) and the numbers
// of all and of good primary vertices. One row per event. Names, documentation, types and precisions as the central
// producer; the uint8 counts throw above 255 vertices instead of wrapping.

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>

class PrimaryVertexTableProducer : public edm::global::EDProducer<> {
public:
  explicit PrimaryVertexTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  const edm::EDGetTokenT<reco::VertexCollection> token_;
  const std::string goodPvCutString_;
  const StringCutObjectSelector<reco::Vertex> goodPvCut_;
  const std::string name_;
};

PrimaryVertexTableProducer::PrimaryVertexTableProducer(const edm::ParameterSet& ps)
    : token_(consumes(ps.getParameter<edm::InputTag>("src"))),
      goodPvCutString_(ps.getParameter<std::string>("goodPvCut")),
      goodPvCut_(goodPvCutString_),
      name_(ps.getParameter<std::string>("name")) {
  produces<nanoaod::FlatTable>();
}

void PrimaryVertexTableProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const {
  const auto& pvs = event.get(token_);
  // RECO's primary-vertex collection holds at least the beam-spot vertex
  if (pvs.empty())
    throw cms::Exception("InvalidInput") << "PrimaryVertexTableProducer: empty primary-vertex collection";
  if (pvs.size() > std::numeric_limits<uint8_t>::max())
    throw cms::Exception("InvalidInput") << "PrimaryVertexTableProducer: " << pvs.size() << " vertices exceed uint8";
  const reco::Vertex& pv = pvs.front();
  const auto good = std::count_if(pvs.begin(), pvs.end(), [this](const reco::Vertex& v) { return goodPvCut_(v); });

  auto table = std::make_unique<nanoaod::FlatTable>(1, name_, true);
  table->addColumnValue<float>("ndof", pv.ndof(), "main primary vertex number of degree of freedom", 8);
  table->addColumnValue<float>("x", pv.position().x(), "main primary vertex position x coordinate", 10);
  table->addColumnValue<float>("y", pv.position().y(), "main primary vertex position y coordinate", 10);
  table->addColumnValue<float>("z", pv.position().z(), "main primary vertex position z coordinate", 16);
  table->addColumnValue<float>("chi2", pv.normalizedChi2(), "main primary vertex reduced chi2", 8);
  table->addColumnValue<uint8_t>("npvs", pvs.size(), "total number of reconstructed primary vertices");
  table->addColumnValue<uint8_t>(
      "npvsGood", good, "number of good reconstructed primary vertices. selection:" + goodPvCutString_);
  event.put(std::move(table));
}

void PrimaryVertexTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("std::vector<reco::Vertex>: the primary vertices, main vertex first");
  desc.add<std::string>("goodPvCut")->setComment("selection of a good primary vertex (StringCutObjectSelector)");
  desc.add<std::string>("name")->setComment("name of the table");
  descriptions.add("primaryVertexTableProducer", desc);
}

DEFINE_FWK_MODULE(PrimaryVertexTableProducer);
