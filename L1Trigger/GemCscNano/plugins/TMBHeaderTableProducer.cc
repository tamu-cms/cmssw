// The (O)TMB DAQ header of each chamber readout, from the CSC unpacker's CSCTMBStatusDigi (muonCSCDigis with
// UnpackStatusDigis = True), one row per status digi in collection order (detector unit, then digi). CSCTMBStatusDigi
// keeps only header words 0-26 of the 43 (CSCTMBStatusDigi.cc:12-16): every column is a field that lies entirely in
// words 0-26; words 27-42 (ALCT words, MPC frames, GEM configuration, buffer status, HMT) are not available.
//   TMBHeader     chamber, decoded header format, firmware revision, and every field of words 0-26 of the format's
//                 CMSSW decoder (the CSCTMBHeader20xx struct `bits`, read by name); -1 where the format lacks a field;
//                 decoded for the 2013 and 2020 formats only (2006, 2007: identification and raw words);
//   TMBHeaderRaw  the same rows: header words 0-26 and the 7 trailer words the status digi holds.
// The format is chosen as CSCTMBHeader(const unsigned short*) chooses it (CSCTMBHeader.cc:140-236). A self-test at
// construction checks that no column reads beyond word 26.

#include "DataFormats/CSCDigi/interface/CSCTMBStatusDigi.h"
#include "DataFormats/CSCDigi/interface/CSCTMBStatusDigiCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2007.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2013.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_CCLUT.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_GEM.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_GEM_Run3b.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_Run2.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_TMB.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "L1Trigger/GemCscNano/interface/TableColumns.h"

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using gemcscnano::addColumn;

namespace {

  constexpr unsigned kStoredHeaderWords = 27;  // CSCTMBStatusDigi.cc:13 (54 bytes)
  constexpr unsigned kStoredTrailerWords = 7;  // CSCTMBTrailer::theData[7], CSCTMBTrailer.h:52 (16 bytes copied)
  constexpr unsigned kHeaderWords = 43;        // the 2013 and 2020 header formats

  // Header format codes (CSCTMBHeader.cc:140-236)
  const char* const kFormatDoc =
      "header format, as CSCTMBHeader chooses it: 1 = 2006, 2 = 2007, 3 = 2007_rev0x50c3, 4 = 2013, 5 = 2020_Run2, "
      "6 = 2020_TMB (copper TMB hybrid), 7 = 2020_CCLUT, 8 = 2020_GEM, 9 = 2020_GEM_Run3b, 0 = none";

  int formatCode(const unsigned short* buf) {
    if (buf[0] == 0x6B0C)
      return 1;
    if (buf[0] != 0xDB0C)
      return 0;
    const unsigned rev = CSCTMBHeader2007(buf).firmwareRevision();
    if (!(rev >= 0x50c3 || rev < 0x42D5))
      return 2;
    if (!(rev >= 0x6000 || rev < 0x42D5))
      return 3;
    if (!(rev < 0x4000 && rev > 0x0))
      return 4;
    const unsigned df = (rev >> 9) & 0xF, major = (rev >> 5) & 0xF;
    bool gem = false, cclut = false, otmbRun2 = false, tmbRun2 = false, tmbRun3 = false, hybrid = false, run3b = false;
    switch (df) {
      case 0x4:
        hybrid = true;
        break;
      case 0x3:
        gem = true;
        run3b = major == 2;
        break;
      case 0x2:
        cclut = true;
        break;
      case 0x1:
        otmbRun2 = true;
        break;
      case 0x0:
        (major == 1 ? tmbRun2 : tmbRun3) = true;
        break;
      default:
        gem = true;
    }
    const bool run2df = major == 1;
    if (gem)
      return run2df ? 5 : (run3b ? 9 : 8);
    if (cclut)
      return run2df ? 5 : 7;
    if (otmbRun2 || tmbRun2 || run2df)
      return 5;
    if (tmbRun3)
      return 7;
    if (hybrid)
      return 6;
    throw cms::Exception("LogicError") << "TMBHeaderTableProducer: no header format for revision " << rev;
  }

  // CSCTMBHeader::FirmwareVersion() of each format code
  int firmwareVersionOf(int code) { return code == 1 ? 2006 : code <= 3 ? 2007 : code == 4 ? 2013 : 2020; }

// The decoded columns: name, the field(s) of `bits` that hold it (the first one the format has), doc. Every field
// lies in header words 0-26 in the formats that have it (checked by the self-test).
#define GCN_TMB_COLUMNS(X, X2)                                                                                         \
  X(bxnCount, bxnCount, "BXN at which the L1A arrived (12 bits)")                                                      \
  X(l1aNumber, l1aNumber, "L1A counter (12 bits)")                                                                     \
  X(readoutCounter, readoutCounter, "counter of DAQ readouts sent (12 bits)")                                          \
  X(boardID, boardID, "board ID = VME slot (configuration)")                                                           \
  X(cscID, cscID, "chamber ID in the crate (configuration; MC: the DMB ID)")                                           \
  X(runID, runID, "run ID (4 bits, configuration)")                                                                    \
  X(stackOvf, stackOvf, "event-buffer stack overflow (at readout)")                                                    \
  X(syncError, syncError, "BXN sync error: BX0 did not match the BXN counter")                                         \
  X(nHeaderFrames, nHeaderFrames, "number of header words (42)")                                                       \
  X(fifoMode, fifoMode, "DAQ FIFO mode (configuration; 1 = full dump)")                                                \
  X(readoutType, r_type, "readout type: 0 header only, 1 full dump, 2 local dump, 3 short header")                     \
  X(l1aType, l1atype, "L1A pop type: 0 CLCT trigger + L1A in window, 1 ALCT only, 2 L1A without trigger, 3 no L1A")    \
  X(hasBuf, hasBuf, "the event has raw-hits buffer data")                                                              \
  X(bufFull, bufFull, "raw-hits buffer stalled at pre-trigger")                                                        \
  X(bdStatus, bd_status, "board status (15 bits, VME 0xEA) at readout")                                                \
  X(bxnPreTrigger, bxnPreTrigger, "BXN at the CLCT pre-trigger (the CLCT fullBX)")                                     \
  X(clct0Discard, tmb_clct0_discard, "CLCT0 discarded (in ME1a)")                                                      \
  X(clct1Discard, tmb_clct1_discard, "CLCT1 discarded (in ME1a)")                                                      \
  X(clockLockLost, clock_lock_lost, "40 MHz clock DLL lost lock (latched)")                                            \
  X2(preTrigCounter, preTrigCounterLow, preTrigCounter, "CLCT pre-trigger counter, low 15 bits")                       \
  X(preTrigCounterHigh, preTrigCounterHigh, "CLCT pre-trigger counter, bits 29:15 (formats 4-6)")                      \
  X(clctCounter, clctCounterLow, "counter of CLCTs sent to the TMB section, low 15 bits")                              \
  X(clctCounterHigh, clctCounterHigh, "CLCT counter, bits 29:15 (formats 4-7)")                                        \
  X2(trigCounter, trigCounterLow, trigCounter, "TMB trigger counter, low 15 bits")                                     \
  X(trigCounterHigh, trigCounterHigh, "TMB trigger counter, bits 29:15 (formats 4-6)")                                 \
  X(alctCounter, alctCounterLow, "counter of ALCTs received, low 15 bits")                                             \
  X(alctCounterHigh, alctCounterHigh, "ALCT counter, bits 29:15 (formats 4-7)")                                        \
  X(uptimeCounterLow, uptimeCounterLow, "OTMB: orbit counter at pre-trigger, bits 14:0 (CMSSW name uptimeCounterLow)") \
  X(uptimeCounterHigh, uptimeCounterHigh, "OTMB: orbit counter at pre-trigger, bits 29:15")                            \
  X(nCFEBs, nCFEBs, "number of CFEBs read out")                                                                        \
  X(nTBins, nTBins, "time bins per CFEB in the dump (configuration)")                                                  \
  X(fifoPretrig, fifoPretrig, "CFEB time bins before the pre-trigger (configuration)")                                 \
  X(scopeExists, scopeExists, "scope data in the readout (configuration)")                                             \
  X(vmeExists, vmeExists, "OTMB: miniscope in the readout (configuration; CMSSW name vmeExists)")                      \
  X(hitThresh, hitThresh, "pre-trigger threshold on layers hit (configuration)")                                       \
  X(pidThresh, pidThresh, "pre-trigger pattern-ID threshold (configuration)")                                          \
  X(nphThresh, nphThresh, "minimum layers hit after drift for a CLCT (configuration)")                                 \
  X(pidThreshPostdrift, pid_thresh_postdrift, "minimum pattern ID after drift (configuration)")                        \
  X(staggerCSC, staggerCSC, "staggered CSC (compile-time constant)")                                                   \
  X(triadPersist, triadPersist, "comparator triad persistence, BX (configuration)")                                    \
  X(dmbThresh, dmbThresh, "layer threshold for the active-CFEB flag to the DMB (configuration)")                       \
  X(alctDelay, alct_delay, "ALCT delay for the CLCT match window (configuration)")                                     \
  X(clctWidth, clct_width, "CLCT match window width (configuration)")                                                  \
  X(trigSourceVect, trigSourceVect, "pre-trigger source bits (9)")                                                     \
  X(nLayersHitVec, r_nlayers_hit_vec, "layers-hit vector (formats 4-6)")                                               \
  X(clct0Slope, clct0_slope, "CLCT0 CCLUT slope magnitude (formats 7-9)")                                              \
  X(clct0LR, clct0_LR_bend, "CLCT0 bend direction (formats 7-9)")                                                      \
  X(clct1LR, clct1_LR_bend, "CLCT1 bend direction (formats 7-9)")                                                      \
  X(activeCFEBsLow, activeCFEBs, "active-CFEB list sent to the DMB, bits 4:0 (bits 6:5 are in word 40, not kept)")     \
  X(readCFEBsLow, readCFEBs, "CFEBs read out, bits 4:0 (bits 6:5 are in word 40, not kept)")                           \
  X(l1aMatchWin, pop_l1a_match_win, "position of the L1A in the L1A-TMB match window")                                 \
  X(affSource, aff_source, "source of the active-CFEB list: 0 pre-trigger, 1 TMB match (configuration)")               \
  X(tmbMatch, tmbMatch, "ALCT and CLCT matched in time")                                                               \
  X(alctOnly, alctOnly, "only the ALCT triggered")                                                                     \
  X(clctOnly, clctOnly, "only the CLCT triggered")                                                                     \
  X(alctMatchWin, matchWin, "position of the ALCT in the CLCT match window (4 bits)")                                  \
  X(noALCT, noALCT, "no ALCT at the TMB match")                                                                        \
  X(oneALCT, oneALCT, "one ALCT at the TMB match")                                                                     \
  X(oneCLCT, oneCLCT, "one CLCT at the TMB match")                                                                     \
  X(twoALCT, twoALCT, "two ALCTs at the TMB match")                                                                    \
  X(twoCLCT, twoCLCT, "two CLCTs at the TMB match")                                                                    \
  X(dupeALCT, dupeALCT, "ALCT0 copied into ALCT1")                                                                     \
  X(dupeCLCT, dupeCLCT, "CLCT0 copied into CLCT1")                                                                     \
  X(lctRankErr, lctRankErr, "LCT1 of higher quality than LCT0")                                                        \
  X(clct0Valid, clct0_valid, "CLCT0 valid after drift (before the TMB match)")                                         \
  X(clct0Quality, clct0_quality, "CLCT0 layers hit")                                                                   \
  X(clct0Shape, clct0_shape, "CLCT0 pattern ID (4 bits)")                                                              \
  X(clct1Valid, clct1_valid, "CLCT1 valid after drift")                                                                \
  X(clct1Quality, clct1_quality, "CLCT1 layers hit")                                                                   \
  X(clct1Shape, clct1_shape, "CLCT1 pattern ID (4 bits)")                                                              \
  X(clct0CompCode, clct0_comparator_code, "CLCT0 comparator code (12 bits, formats 7-9)")                              \
  X(clct0Xky, clct0_xky, "CLCT0 quarter/eighth-strip bits (formats 7-9)")                                              \
  X(clct1CompCode, clct1_comparator_code, "CLCT1 comparator code (12 bits, formats 7-9)")                              \
  X(clct1Xky, clct1_xky, "CLCT1 quarter/eighth-strip bits (formats 7-9)")                                              \
  X(lct0NoGem, lct0_nogem, "LCT0 from ALCT+CLCT without GEM (formats 8, 9)")                                           \
  X(lct0WithGemA, lct0_with_gemA, "LCT0 from ALCT+CLCT+1 GEM cluster in layer 1 (formats 8, 9)")                       \
  X(lct0WithGemB, lct0_with_gemB, "LCT0 from ALCT+CLCT+1 GEM cluster in layer 2 (formats 8, 9)")                       \
  X(lct0WithCopad, lct0_with_copad, "LCT0 from ALCT+CLCT+copad (formats 8, 9)")                                        \
  X(lct1NoGem, lct1_nogem, "LCT1 from ALCT+CLCT without GEM (formats 8, 9)")                                           \
  X(lct1WithGemA, lct1_with_gemA, "LCT1 from ALCT+CLCT+1 GEM cluster in layer 1 (formats 8, 9)")                       \
  X(lct1WithGemB, lct1_with_gemB, "LCT1 from ALCT+CLCT+1 GEM cluster in layer 2 (formats 8, 9)")                       \
  X(lct1WithCopad, lct1_with_copad, "LCT1 from ALCT+CLCT+copad (formats 8, 9)")                                        \
  X(gemAvpf, gemA_vpf, "a GEM layer-1 cluster was available for the match (formats 8, 9)")                             \
  X(gemBvpf, gemB_vpf, "a GEM layer-2 cluster was available for the match (formats 8, 9)")                             \
  X(gemAOverflow, gemA_over_flow, "GEM layer-1 cluster overflow (formats 8, 9)")                                       \
  X(gemBOverflow, gemB_over_flow, "GEM layer-2 cluster overflow (formats 8, 9)")                                       \
  X(gemASyncErr, gemA_sync, "GEM layer-1 fibers not in sync: an error flag (CMSSW name gemA_sync; formats 8, 9)")      \
  X(gemBSyncErr, gemB_sync, "GEM layer-2 fibers not in sync: an error flag (formats 8, 9)")                            \
  X(gemsSyncErr, gems_sync, "GEM layer 1 vs layer 2 not in sync: an error flag (formats 8, 9)")                        \
  X(numCopad, num_copad, "number of GEM copads at the TMB match (formats 8, 9)")                                       \
  X(gemDelay, gem_delay, "GEM-ALCT match delay, low 4 bits of the register (configuration; formats 8, 9)")             \
  X(gemClctWin, gem_clct_win, "position of the CLCT in the GEM-CLCT match window (formats 8, 9)")                      \
  X(alctGemWin, alct_gem_win, "position of the GEM in the ALCT-GEM match window (formats 8, 9)")

  constexpr unsigned kNDecoded = 0
#define GCN_COUNT(name, field, doc) +1
#define GCN_COUNT2(name, a, b, doc) +1
      GCN_TMB_COLUMNS(GCN_COUNT, GCN_COUNT2);
#undef GCN_COUNT
#undef GCN_COUNT2

  struct ColumnSpec {
    const char* name;
    const char* doc;
  };
  const std::array<ColumnSpec, kNDecoded> kColumns = {{
#define GCN_SPEC(name, field, doc) {#name, doc},
#define GCN_SPEC2(name, a, b, doc) {#name, doc},
      GCN_TMB_COLUMNS(GCN_SPEC, GCN_SPEC2)
#undef GCN_SPEC
#undef GCN_SPEC2
  }};

  using Values = std::array<int, kNDecoded>;

  // the decoded columns of one header format: each field the format's `bits` has, -1 for the others
  template <typename Format>
  Values decode(const unsigned short* buf) {
    const Format header(buf);
    const auto& b = header.bits;
    Values v;
    unsigned i = 0;
#define GCN_GET(name, field, doc)      \
  if constexpr (requires { b.field; }) \
    v[i++] = b.field;                  \
  else                                 \
    v[i++] = -1;
#define GCN_GET2(name, a, c, doc)       \
  if constexpr (requires { b.a; })      \
    v[i++] = b.a;                       \
  else if constexpr (requires { b.c; }) \
    v[i++] = b.c;                       \
  else                                  \
    v[i++] = -1;
    GCN_TMB_COLUMNS(GCN_GET, GCN_GET2)
#undef GCN_GET
#undef GCN_GET2
    return v;
  }

  // which field names of each column a format's `bits` has: bit 0 the (first) name, bit 1 the second name of an X2
  // column; the self-test requires every name to exist in some format
  using Names = std::array<unsigned, kNDecoded>;
  template <typename Format>
  Names namesPresent() {
    using Bits = decltype(Format(nullptr).bits);
    Names n{};
    unsigned i = 0;
#define GCN_HAS(name, field, doc) n[i++] = requires(const Bits& b) { b.field; } ? 1u : 0u;
#define GCN_HAS2(name, a, c, doc) \
  n[i++] = (requires(const Bits& b) { b.a; } ? 1u : 0u) | (requires(const Bits& b) { b.c; } ? 2u : 0u);
    GCN_TMB_COLUMNS(GCN_HAS, GCN_HAS2)
#undef GCN_HAS
#undef GCN_HAS2
    return n;
  }
  const Names kExpectedNames = {{
#define GCN_EXP(name, field, doc) 1u,
#define GCN_EXP2(name, a, b, doc) 3u,
      GCN_TMB_COLUMNS(GCN_EXP, GCN_EXP2)
#undef GCN_EXP
#undef GCN_EXP2
  }};

  Values decodeFormat(int code, const unsigned short* buf) {
    switch (code) {
      case 4:
        return decode<CSCTMBHeader2013>(buf);
      case 5:
        return decode<CSCTMBHeader2020_Run2>(buf);
      case 6:
        return decode<CSCTMBHeader2020_TMB>(buf);
      case 7:
        return decode<CSCTMBHeader2020_CCLUT>(buf);
      case 8:
        return decode<CSCTMBHeader2020_GEM>(buf);
      case 9:
        return decode<CSCTMBHeader2020_GEM_Run3b>(buf);
      default: {
        Values none;
        none.fill(-1);
        return none;
      }
    }
  }

}  // namespace

class TMBHeaderTableProducer : public edm::global::EDProducer<> {
public:
  explicit TMBHeaderTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  const edm::EDGetTokenT<CSCTMBStatusDigiCollection> token_;
  const std::string label_;
};

TMBHeaderTableProducer::TMBHeaderTableProducer(const edm::ParameterSet& ps)
    : token_(consumes(ps.getParameter<edm::InputTag>("src"))), label_(ps.getParameter<edm::InputTag>("src").encode()) {
  // Self-test: with words 0-26 zero and words 27-42 all ones, every decoded field reads 0 (or -1 if absent), i.e.
  // no column reads a word the status digi does not keep; and every column exists in at least one format (a field
  // name missing from every format's `bits` would silently give -1).
  std::array<bool, kNDecoded> present{};
  for (int code = 4; code <= 9; ++code) {
    std::array<unsigned short, kHeaderWords> buf{};
    for (unsigned w = kStoredHeaderWords; w < kHeaderWords; ++w)
      buf[w] = 0xFFFF;
    const Values v = decodeFormat(code, buf.data());
    for (unsigned i = 0; i < kNDecoded; ++i) {
      if (v[i] != 0 && v[i] != -1)
        throw cms::Exception("LogicError") << "TMBHeaderTableProducer: column " << kColumns[i].name << " of format "
                                           << code << " reads beyond header word " << kStoredHeaderWords - 1;
      present[i] = present[i] || v[i] == 0;
    }
  }
  for (unsigned i = 0; i < kNDecoded; ++i)
    if (!present[i])
      throw cms::Exception("LogicError") << "TMBHeaderTableProducer: column " << kColumns[i].name
                                         << " is in no header format";
  Names names{};
  for (const Names& n : {namesPresent<CSCTMBHeader2013>(),
                         namesPresent<CSCTMBHeader2020_Run2>(),
                         namesPresent<CSCTMBHeader2020_TMB>(),
                         namesPresent<CSCTMBHeader2020_CCLUT>(),
                         namesPresent<CSCTMBHeader2020_GEM>(),
                         namesPresent<CSCTMBHeader2020_GEM_Run3b>()})
    for (unsigned i = 0; i < kNDecoded; ++i)
      names[i] |= n[i];
  for (unsigned i = 0; i < kNDecoded; ++i)
    if (names[i] != kExpectedNames[i])
      throw cms::Exception("LogicError") << "TMBHeaderTableProducer: a field name of column " << kColumns[i].name
                                         << " is in no header format";
  produces<nanoaod::FlatTable>("header");
  produces<nanoaod::FlatTable>("raw");
}

void TMBHeaderTableProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const {
  const auto& digis = event.get(token_);
  struct Row {
    CSCDetId chamber;
    std::array<unsigned short, kHeaderWords> header;  // words 0-26 as stored, 27-42 zero
    std::array<unsigned short, kStoredTrailerWords> trailer;
    int format;
    int revision;
    Values values;
  };
  std::vector<Row> rows;
  for (auto detUnit = digis.begin(); detUnit != digis.end(); ++detUnit) {
    const auto& range = (*detUnit).second;
    for (auto digi = range.first; digi != range.second; ++digi) {
      Row r;
      r.chamber = (*detUnit).first;
      r.header.fill(0);
      for (unsigned w = 0; w < kStoredHeaderWords; ++w)
        r.header[w] = digi->header()[w];
      for (unsigned w = 0; w < kStoredTrailerWords; ++w)
        r.trailer[w] = digi->trailer()[w];
      r.format = formatCode(r.header.data());
      // the dispatch must agree with CMSSW's own choice of the firmware version
      if (r.format != 0 && CSCTMBHeader(r.header.data()).FirmwareVersion() != firmwareVersionOf(r.format))
        throw cms::Exception("LogicError") << "TMBHeaderTableProducer: format " << r.format << " but CSCTMBHeader "
                                           << CSCTMBHeader(r.header.data()).FirmwareVersion();
      r.revision = r.header[0] == 0xDB0C ? int(CSCTMBHeader2007(r.header.data()).firmwareRevision()) : -1;
      r.values = decodeFormat(r.format, r.header.data());
      rows.push_back(r);
    }
  }

  auto table = std::make_unique<nanoaod::FlatTable>(rows.size(), "TMBHeader", false);
  table->setDoc("(O)TMB DAQ header of each chamber readout: every status digi of " + label_ +
                ", in collection order; only header words 0-26 are kept by CSCTMBStatusDigi; -1 where the format "
                "lacks a field; MC headers are written by the packer (most fields 0)");
  addColumn<int16_t>(*table, "endcap", rows, [](const Row& r) { return r.chamber.zendcap(); }, "+1 or -1");
  addColumn<int16_t>(*table, "station", rows, [](const Row& r) { return r.chamber.station(); }, "station");
  addColumn<int16_t>(
      *table, "ring", rows, [](const Row& r) { return r.chamber.ring(); }, "ring as stored (ME1/1 readouts: 1)");
  addColumn<int16_t>(*table, "chamber", rows, [](const Row& r) { return r.chamber.chamber(); }, "chamber");
  addColumn<int16_t>(*table, "format", rows, [](const Row& r) { return r.format; }, kFormatDoc);
  addColumn<int16_t>(
      *table,
      "firmwareRevision",
      rows,
      [](const Row& r) { return r.revision; },
      "firmware revision code (header word 7, 15 bits; formats 5-9: data format (rev>>9)&0xF, major (rev>>5)&0xF, "
      "minor rev&0x1F); -1 for the 2006 format");
  for (unsigned i = 0; i < kNDecoded; ++i)
    addColumn<int16_t>(*table, kColumns[i].name, rows, [i](const Row& r) { return r.values[i]; }, kColumns[i].doc);

  auto raw = std::make_unique<nanoaod::FlatTable>(rows.size(), "TMBHeaderRaw", false);
  raw->setDoc(
      "raw words of the TMBHeader rows (same rows): header words 0-26 and the 7 trailer words of the status digi");
  for (unsigned w = 0; w < kStoredHeaderWords; ++w)
    addColumn<uint16_t>(
        *raw,
        "header" + std::to_string(w),
        rows,
        [w](const Row& r) { return r.header[w]; },
        "header word " + std::to_string(w));
  for (unsigned w = 0; w < kStoredTrailerWords; ++w)
    addColumn<uint16_t>(
        *raw,
        "trailer" + std::to_string(w),
        rows,
        [w](const Row& r) { return r.trailer[w]; },
        "trailer word " + std::to_string(w) +
            " from 0x6E0C (optional 0x2AAA 0x5555 filler, 0xDE0F, 2 CRC words, word count; without filler the last 2 "
            "words lie beyond the TMB record)");

  event.put(std::move(table), "header");
  event.put(std::move(raw), "raw");
}

void TMBHeaderTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("CSCTMBStatusDigiCollection (muonCSCDigis:MuonCSCTMBStatusDigi)");
  descriptions.add("tmbHeaderTableProducer", desc);
}

DEFINE_FWK_MODULE(TMBHeaderTableProducer);
