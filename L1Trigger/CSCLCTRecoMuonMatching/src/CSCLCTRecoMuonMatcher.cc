#include "L1Trigger/CSCLCTRecoMuonMatching/interface/CSCLCTRecoMuonMatcher.h"

#include "DataFormats/CSCDigi/interface/CSCConstants.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <string>

namespace {

  CSCLCTRecoMuonMatcher::SegmentSource parseSegmentSource(const std::string& s) {
    if (s == "standalone")
      return CSCLCTRecoMuonMatcher::SegmentSource::Standalone;
    if (s == "matches")
      return CSCLCTRecoMuonMatcher::SegmentSource::Matches;
    if (s == "auto")
      return CSCLCTRecoMuonMatcher::SegmentSource::Auto;
    throw cms::Exception("Configuration")
        << "segmentSource must be 'standalone', 'matches' or 'auto', got '" << s << "'";
  }

  // The CSC segments used in the standalone fit: the 4D CSCSegment hits of the standalone track (kept in AOD).
  std::vector<const CSCSegment*> standaloneSegments(const reco::Muon& mu) {
    std::vector<const CSCSegment*> segs;
    const reco::TrackRef trk = mu.outerTrack();
    if (trk.isNull() || !trk.isAvailable() || !trk->recHitsOk())
      return segs;
    for (auto const* hit : trk->recHits()) {
      if (!hit->isValid() || hit->dimension() != 4)
        continue;
      const DetId id = hit->geographicalId();
      if (id.det() != DetId::Muon || id.subdetId() != MuonSubdetId::CSC)
        continue;
      if (auto const* seg = dynamic_cast<const CSCSegment*>(hit))
        segs.push_back(seg);
    }
    return segs;
  }

  // The segment arbitrated to the muon in each crossed chamber (SegmentAndTrackArbitration; standalone-only and
  // global-only muons receive both bits uncontested).
  std::vector<const CSCSegment*> matchedSegments(const reco::Muon& mu) {
    std::vector<const CSCSegment*> segs;
    for (const auto& chamber : mu.matches()) {
      if (chamber.detector() != MuonSubdetId::CSC)
        continue;
      for (const auto& sm : chamber.segmentMatches) {
        if (!sm.isMask(reco::MuonSegmentMatch::BestInChamberByDR) ||
            !sm.isMask(reco::MuonSegmentMatch::BelongsToTrackByDR))
          continue;
        if (sm.cscSegmentRef.isNull() || !sm.cscSegmentRef.isAvailable())
          continue;
        segs.push_back(&*sm.cscSegmentRef);
      }
    }
    return segs;
  }

  // A muon in one LCT chamber, represented by one segment.
  struct Candidate {
    unsigned int muonIndex;
    const CSCSegment* segment;
    CSCLCTRecoMuonMatcher::SegmentSource source;
    csclct::KeyLayerPosition key;
    bool keyFromRecHit;
  };

}  // namespace

CSCLCTRecoMuonMatcher::Config::Config(SegmentSource segmentSource,
                                      int maxDeltaHalfStrip,
                                      int maxDeltaWireGroup,
                                      int lctBXWindow)
    : segmentSource(segmentSource),
      maxDeltaHalfStrip(maxDeltaHalfStrip),
      maxDeltaWireGroup(maxDeltaWireGroup),
      lctBXWindow(lctBXWindow) {}

// The members of Config are public, so the windows are checked when the matcher is built.
CSCLCTRecoMuonMatcher::CSCLCTRecoMuonMatcher(const Config& config) : config_(config) {
  if (config_.maxDeltaHalfStrip < 0 || config_.maxDeltaWireGroup < 0)
    throw cms::Exception("Configuration")
        << "CSCLCTRecoMuonMatcher: negative match window (maxDeltaHalfStrip " << config_.maxDeltaHalfStrip
        << ", maxDeltaWireGroup " << config_.maxDeltaWireGroup << ")";
}

CSCLCTRecoMuonMatcher::CSCLCTRecoMuonMatcher(const edm::ParameterSet& pset)
    : CSCLCTRecoMuonMatcher(Config(parseSegmentSource(pset.getParameter<std::string>("segmentSource")),
                                   pset.getParameter<int>("maxDeltaHalfStrip"),
                                   pset.getParameter<int>("maxDeltaWireGroup"),
                                   pset.getParameter<int>("lctBXWindow"))) {}

// Default windows: deliberately loose (user, 2026-10-21), to be tightened from the residual distributions.
void CSCLCTRecoMuonMatcher::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<std::string>("segmentSource", "auto")
      ->setComment(
          "standalone: CSC segments of the standalone track; matches: arbitrated reco::Muon chamber matches; "
          "auto: standalone when the standalone track has at least one CSC segment, otherwise matches "
          "(default auto: a proposal of the matching study, not decided)");
  desc.add<int>("maxDeltaHalfStrip", 6)
      ->setComment("match window |LCT - segment| in key half-strips (default 6: loose, to be tightened)");
  desc.add<int>("maxDeltaWireGroup", 2)
      ->setComment("match window |LCT - segment| in key wire groups (default 2: loose, to be tightened)");
  desc.add<int>("lctBXWindow", -1)
      ->setComment("keep LCTs with |BX - 8| <= lctBXWindow; < 0: no BX selection (data LCTs carry a 1-bit BX)");
}

std::vector<std::pair<const CSCSegment*, CSCLCTRecoMuonMatcher::SegmentSource>> CSCLCTRecoMuonMatcher::muonSegments(
    const reco::Muon& muon) const {
  std::vector<std::pair<const CSCSegment*, SegmentSource>> result;
  if (config_.segmentSource != SegmentSource::Matches)
    for (auto const* seg : standaloneSegments(muon))
      result.emplace_back(seg, SegmentSource::Standalone);
  if (config_.segmentSource == SegmentSource::Matches ||
      (config_.segmentSource == SegmentSource::Auto && result.empty()))
    for (auto const* seg : matchedSegments(muon))
      result.emplace_back(seg, SegmentSource::Matches);
  return result;
}

std::vector<CSCLCTRecoMuonMatcher::LCTMatch> CSCLCTRecoMuonMatcher::match(const CSCCorrelatedLCTDigiCollection& lcts,
                                                                          const reco::MuonCollection& muons,
                                                                          const CSCGeometry& geometry) const {
  // Candidates per LCT chamber: every segment of every muon, with its key-layer position. ME1/1 ring-1 (ME1b) and
  // ring-4 (ME1a) segments share the ring-1 LCT chamber, so a muon crossing both can match an ME1b and an ME1a LCT.
  std::map<CSCDetId, std::vector<Candidate>> candidates;
  for (unsigned int m = 0; m < muons.size(); ++m) {
    for (const auto& [seg, source] : muonSegments(muons[m])) {
      Candidate c{m, seg, source, {}, false};
      if (const CSCRecHit2D* rh = csclct::keyLayerRecHit(*seg)) {
        c.key = csclct::fromKeyLayerRecHit(geometry, *rh);
        c.keyFromRecHit = c.key.valid();
      }
      if (!c.keyFromRecHit)
        c.key = csclct::fromSegmentProjection(geometry, *seg);
      if (c.key.valid())
        candidates[csclct::lctChamberId(seg->cscDetId())].push_back(c);
    }
  }

  std::vector<LCTMatch> result;
  for (auto detUnit = lcts.begin(); detUnit != lcts.end(); ++detUnit) {
    const CSCDetId chamber = csclct::lctChamberId((*detUnit).first);
    const auto cand = candidates.find(chamber);
    const auto& range = (*detUnit).second;
    for (auto lct = range.first; lct != range.second; ++lct) {
      if (!lct->isValid())
        continue;
      if (config_.lctBXWindow >= 0 &&
          std::abs(int(lct->getBX()) - int(CSCConstants::LCT_CENTRAL_BX)) > config_.lctBXWindow)
        continue;
      LCTMatch lm{chamber, &*lct, {}};
      if (cand != candidates.end()) {
        for (const Candidate& c : cand->second) {
          const int dHS = int(lct->getStrip()) - c.key.hs;
          const float dES = float(lct->getStrip(8)) - c.key.es;
          const int dWG = int(lct->getKeyWG()) - c.key.wg;
          if (std::abs(dHS) > config_.maxDeltaHalfStrip || std::abs(dWG) > config_.maxDeltaWireGroup)
            continue;
          lm.muons.push_back({c.muonIndex,
                              c.segment,
                              c.source,
                              c.key,
                              c.keyFromRecHit,
                              dHS,
                              dES,
                              dWG,
                              std::abs(dES) / 4.f + std::abs(dWG)});
        }
        // Closest first; a muon with several segments within the window keeps only its closest one.
        std::stable_sort(lm.muons.begin(), lm.muons.end(), [](const MuonMatch& a, const MuonMatch& b) {
          return a.distance < b.distance;
        });
        std::vector<MuonMatch> unique;
        for (const MuonMatch& mm : lm.muons)
          if (std::none_of(
                  unique.begin(), unique.end(), [&mm](const MuonMatch& u) { return u.muonIndex == mm.muonIndex; }))
            unique.push_back(mm);
        lm.muons = std::move(unique);
      }
      result.push_back(std::move(lm));
    }
  }
  return result;
}
