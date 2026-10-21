// Test analyzer of CSCLCTRecoMuonMatcher: runs the matcher on one LCT collection and the reco muons, and writes
// one TTree entry per LCT with the muons matched to it (vector branches, closest first).

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "L1Trigger/CSCLCTRecoMuonMatching/interface/CSCLCTRecoMuonMatcher.h"

#include "TTree.h"

#include <algorithm>
#include <vector>

class CSCLCTRecoMuonMatchAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit CSCLCTRecoMuonMatchAnalyzer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  const edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lctToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const CSCLCTRecoMuonMatcher matcher_;

  TTree* tree_;
  unsigned int run_, lumi_;
  unsigned long long event_;
  int endcap_, station_, ring_, chamber_, hs_, es_, wg_, quality_, pattern_, run3Pattern_, slope_, slope6_, bend_, bx_,
      trknmb_, nMuons_;
  std::vector<int> muIndex_, muCharge_, muIsGlobal_, muIsTracker_, muIsStandAlone_, segSource_, segRing_, segNHits_,
      keyFromRecHit_, keyHS_, keyWG_, dHS_, dWG_;
  std::vector<float> muPt_, muEta_, muPhi_, segChi2ndof_, segTime_, keyES_, dES_;
};

CSCLCTRecoMuonMatchAnalyzer::CSCLCTRecoMuonMatchAnalyzer(const edm::ParameterSet& ps)
    : muonToken_(consumes(ps.getParameter<edm::InputTag>("muons"))),
      lctToken_(consumes(ps.getParameter<edm::InputTag>("lcts"))),
      cscGeomToken_(esConsumes()),
      matcher_(ps.getParameter<edm::ParameterSet>("matcher")) {
  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("lctMuon", "CSC LCT -> matched reco muons, one entry per LCT");
  tree_->Branch("run", &run_);
  tree_->Branch("lumi", &lumi_);
  tree_->Branch("event", &event_);
  tree_->Branch("lct_endcap", &endcap_);
  tree_->Branch("lct_station", &station_);
  tree_->Branch("lct_ring", &ring_);
  tree_->Branch("lct_chamber", &chamber_);
  tree_->Branch("lct_hs", &hs_);
  tree_->Branch("lct_es", &es_);
  tree_->Branch("lct_wg", &wg_);
  tree_->Branch("lct_quality", &quality_);
  tree_->Branch("lct_pattern", &pattern_);
  tree_->Branch("lct_run3Pattern", &run3Pattern_);
  tree_->Branch("lct_slope", &slope_);
  tree_->Branch("lct_slope6", &slope6_);
  tree_->Branch("lct_bend", &bend_);
  tree_->Branch("lct_bx", &bx_);
  tree_->Branch("lct_trknmb", &trknmb_);
  tree_->Branch("nMuons", &nMuons_);
  tree_->Branch("mu_index", &muIndex_);
  tree_->Branch("mu_pt", &muPt_);
  tree_->Branch("mu_eta", &muEta_);
  tree_->Branch("mu_phi", &muPhi_);
  tree_->Branch("mu_charge", &muCharge_);
  tree_->Branch("mu_isGlobal", &muIsGlobal_);
  tree_->Branch("mu_isTracker", &muIsTracker_);
  tree_->Branch("mu_isStandAlone", &muIsStandAlone_);
  tree_->Branch("seg_source", &segSource_);
  tree_->Branch("seg_ring", &segRing_);
  tree_->Branch("seg_nHits", &segNHits_);
  tree_->Branch("seg_chi2ndof", &segChi2ndof_);
  tree_->Branch("seg_time", &segTime_);
  tree_->Branch("key_fromRecHit", &keyFromRecHit_);
  tree_->Branch("key_hs", &keyHS_);
  tree_->Branch("key_es", &keyES_);
  tree_->Branch("key_wg", &keyWG_);
  tree_->Branch("dHS", &dHS_);
  tree_->Branch("dES", &dES_);
  tree_->Branch("dWG", &dWG_);
}

void CSCLCTRecoMuonMatchAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("lcts", edm::InputTag("muonCSCDigis", "MuonCSCCorrelatedLCTDigi"))
      ->setComment("LCT collection: unpacked data LCTs, or emulator LCTs (then set matcher.lctBXWindow)");
  edm::ParameterSetDescription matcher;
  CSCLCTRecoMuonMatcher::fillPSetDescription(matcher);
  desc.add<edm::ParameterSetDescription>("matcher", matcher);
  descriptions.add("cscLCTRecoMuonMatchAnalyzer", desc);
}

void CSCLCTRecoMuonMatchAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  const auto& muons = event.get(muonToken_);
  const auto matches = matcher_.match(event.get(lctToken_), muons, setup.getData(cscGeomToken_));

  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();
  event_ = event.id().event();
  for (const auto& lm : matches) {
    const CSCCorrelatedLCTDigi& lct = *lm.lct;
    endcap_ = lm.chamber.zendcap();
    station_ = lm.chamber.station();
    ring_ = lm.chamber.ring();
    chamber_ = lm.chamber.chamber();
    hs_ = lct.getStrip();
    es_ = lct.getStrip(8);
    wg_ = lct.getKeyWG();
    quality_ = lct.getQuality();
    pattern_ = lct.getPattern();
    run3Pattern_ = lct.getRun3Pattern();
    slope_ = lct.getSlope();
    slope6_ = lct.getSlopeEx(6);
    bend_ = lct.getBend();
    bx_ = lct.getBX();
    trknmb_ = lct.getTrknmb();
    nMuons_ = lm.muons.size();
    for (auto* v : {&muIndex_,
                    &muCharge_,
                    &muIsGlobal_,
                    &muIsTracker_,
                    &muIsStandAlone_,
                    &segSource_,
                    &segRing_,
                    &segNHits_,
                    &keyFromRecHit_,
                    &keyHS_,
                    &keyWG_,
                    &dHS_,
                    &dWG_})
      v->clear();
    for (auto* v : {&muPt_, &muEta_, &muPhi_, &segChi2ndof_, &segTime_, &keyES_, &dES_})
      v->clear();
    for (const auto& mm : lm.muons) {
      const reco::Muon& mu = muons[mm.muonIndex];
      muIndex_.push_back(mm.muonIndex);
      muPt_.push_back(mu.pt());
      muEta_.push_back(mu.eta());
      muPhi_.push_back(mu.phi());
      muCharge_.push_back(mu.charge());
      muIsGlobal_.push_back(mu.isGlobalMuon());
      muIsTracker_.push_back(mu.isTrackerMuon());
      muIsStandAlone_.push_back(mu.isStandAloneMuon());
      segSource_.push_back(mm.source == CSCLCTRecoMuonMatcher::SegmentSource::Standalone ? 0 : 1);
      segRing_.push_back(mm.segment->cscDetId().ring());
      segNHits_.push_back(mm.segment->nRecHits());
      segChi2ndof_.push_back(mm.segment->chi2() / std::max(1, mm.segment->degreesOfFreedom()));
      segTime_.push_back(mm.segment->time());
      keyFromRecHit_.push_back(mm.keyFromRecHit);
      keyHS_.push_back(mm.key.hs);
      keyES_.push_back(mm.key.es);
      keyWG_.push_back(mm.key.wg);
      dHS_.push_back(mm.deltaHalfStrip);
      dES_.push_back(mm.deltaEighthStrip);
      dWG_.push_back(mm.deltaWireGroup);
    }
    tree_->Fill();
  }
}

DEFINE_FWK_MODULE(CSCLCTRecoMuonMatchAnalyzer);
