// A copy of an edm::TriggerResults with only the selected paths: those whose name starts with `prefix` and whose rest
// (after the prefix) matches the regular expression `select` and not `veto` (ECMAScript). The statuses are copied and
// the names are stored in the product, so NanoAODOutputModule writes them as <name> bits (version suffix dropped).

#include "DataFormats/Common/interface/HLTGlobalStatus.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <memory>
#include <regex>
#include <string>
#include <vector>

class TriggerResultsSelector : public edm::stream::EDProducer<> {
public:
  explicit TriggerResultsSelector(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<edm::TriggerResults> token_;
  const std::string prefix_;
  const std::regex select_;
  const std::regex veto_;
  // the selection of the last trigger menu seen (names change only with the menu)
  std::vector<std::string> lastNames_;
  std::vector<unsigned int> kept_;
};

TriggerResultsSelector::TriggerResultsSelector(const edm::ParameterSet& ps)
    : token_(consumes(ps.getParameter<edm::InputTag>("src"))),
      prefix_(ps.getParameter<std::string>("prefix")),
      select_(ps.getParameter<std::string>("select")),
      veto_(ps.getParameter<std::string>("veto")) {
  produces<edm::TriggerResults>();
}

void TriggerResultsSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("TriggerResults", "", "HLT"))->setComment("input edm::TriggerResults");
  desc.add<std::string>("prefix", "HLT_")->setComment("only names with this prefix are considered");
  desc.add<std::string>("select", "")->setComment("ECMAScript regex the name after the prefix must match");
  desc.add<std::string>("veto", "")->setComment("ECMAScript regex the name after the prefix must not match");
  descriptions.add("triggerResultsSelector", desc);
}

void TriggerResultsSelector::produce(edm::Event& event, const edm::EventSetup&) {
  const auto& results = event.get(token_);
  const edm::TriggerNames& names = event.triggerNames(results);
  if (names.triggerNames() != lastNames_) {
    lastNames_ = names.triggerNames();
    kept_.clear();
    for (unsigned int i = 0; i < lastNames_.size(); ++i) {
      const std::string& name = lastNames_[i];
      if (name.compare(0, prefix_.size(), prefix_) != 0)
        continue;
      const std::string rest = name.substr(prefix_.size());
      if (std::regex_search(rest, select_) && !std::regex_search(rest, veto_))
        kept_.push_back(i);
    }
  }
  edm::HLTGlobalStatus status(kept_.size());
  std::vector<std::string> keptNames;
  keptNames.reserve(kept_.size());
  for (unsigned int k = 0; k < kept_.size(); ++k) {
    status[k] = results[kept_[k]];
    keptNames.push_back(lastNames_[kept_[k]]);
  }
  event.put(std::make_unique<edm::TriggerResults>(status, keptNames));
}

DEFINE_FWK_MODULE(TriggerResultsSelector);
