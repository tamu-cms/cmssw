// The (O)TMB DAQ header of each chamber readout, one row per CSCTMBStatusDigi of the CSC unpacker (muonCSCDigis with
// UnpackStatusDigis = True) in collection order (detector unit, then digi). CSCTMBStatusDigi keeps only header words
// 0-26 of the 43 (CSCTMBStatusDigi.cc:12-16); words 27-42 (ALCT words, MPC frames, GEM configuration, buffer status,
// CFEB bits 6:5, HMT) come from a second decoding of the raw data with the unpacker's own classes and settings
// (tmbReadouts below), matched to the status digi by chamber and by exact equality of words 0-26 and the trailer.
//   TMBHeader          chamber, decoded header format, firmware revision, fullHeaderMatch, and the fields of the
//                      format's CMSSW decoder (the CSCTMBHeader20xx struct `bits`, read by name); -1 where the format
//                      lacks a field, -2 for a field of words 27-42 when fullHeaderMatch != 0; decoded for the 2013
//                      and 2020 formats only (2006, 2007: identification and raw words). The MPC frames (words 31-34,
//                      word 35 bits 3:0) are not decoded here: CMSSW decodes them into the LCT digis (LCT table);
//   TMBHeaderRaw       the same rows: header words 0-26 and the 7 trailer words of the status digi, words 27-42 of
//                      the second decoding (0 when fullHeaderMatch != 0);
//   TMBHeaderUnmatched readouts of the second decoding in chambers without a status digi (normally none).
// The format is chosen as CSCTMBHeader(const unsigned short*) chooses it (CSCTMBHeader.cc:140-236). Self-tests at
// construction check that the columns of words 0-26 read nothing beyond word 26 and the others read words 27-42.

#include "CondFormats/CSCObjects/interface/CSCChamberMap.h"
#include "CondFormats/CSCObjects/interface/CSCCrateMap.h"
#include "CondFormats/DataRecord/interface/CSCChamberMapRcd.h"
#include "CondFormats/DataRecord/interface/CSCCrateMapRcd.h"
#include "DataFormats/CSCDigi/interface/CSCTMBStatusDigi.h"
#include "DataFormats/CSCDigi/interface/CSCTMBStatusDigiCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDCCEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDCCExaminer.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDDUEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2007.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2013.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_CCLUT.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_GEM.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_GEM_Run3b.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_Run2.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader2020_TMB.h"
#include "FWCore/Common/interface/Provenance.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "L1Trigger/GemCscNano/interface/TableColumns.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
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
  X(activeCFEBsLow, activeCFEBs, "active-CFEB list sent to the DMB, bits 4:0 (bits 6:5: activeCFEBsHigh)")             \
  X(readCFEBsLow, readCFEBs, "CFEBs read out, bits 4:0 (bits 6:5: readCFEBsHigh)")                                     \
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
  X(clct0KeyLow, clct0_key_low, "CLCT0 key half-strip, bits 6:0 (key = clct0KeyLow + 128 * clct0KeyHigh)")             \
  X(clct1Valid, clct1_valid, "CLCT1 valid after drift")                                                                \
  X(clct1Quality, clct1_quality, "CLCT1 layers hit")                                                                   \
  X(clct1Shape, clct1_shape, "CLCT1 pattern ID (4 bits)")                                                              \
  X(clct1KeyLow, clct1_key_low, "CLCT1 key half-strip, bits 6:0 (key = clct1KeyLow + 128 * clct1KeyHigh)")             \
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

// The decoded columns of words 27-42 (the second decoding of the raw data). Meanings from the OTMB firmware source
// (otmb_fw_code evka-dev ce999aa, utils/sequencer.v:5514-5646; not the deployed builds, and the copper-TMB firmware of
// format 6 is not in that source: its fields follow the CMSSW struct): word 36 is the GEM word when GEM readout is enabled, the RPC word
// otherwise (:5592), decoded here by format; word 40 bits 5:4 are the CFEB RAM parity errors and bits 7:6 the bad-bits
// flags of CFEBs 6:5 (:5619-5620), the opposite of the CMSSW names; word 41 bits 14:13 in the Run-3 format are ALCT BXN
// bits 2:1, which carry the anode HMT in the Run-3 ALCT format (:5644-5645).
#define GCN_TMB_COLUMNS_27(X, X2)                                                                                      \
  X(clct0KeyHigh, clct0_key_high, "CLCT0 key half-strip bit 7")                                                        \
  X(clct1KeyHigh, clct1_key_high, "CLCT1 key half-strip bit 7")                                                        \
  X(clctBxn, clct_bxn, "CLCT0/1 BXN after drift (2 bits)")                                                             \
  X(clctSyncErr, clct_sync_err, "CLCT BX0 sync error")                                                                 \
  X(clct0InvalidPattern, clct0Invalid, "CLCT0 had an invalid pattern after drift")                                     \
  X(clct1InvalidPattern, clct1Invalid, "CLCT1 had an invalid pattern after drift")                                     \
  X(clct1Busy, clct1Busy, "second CLCT busy (logic error indicator)")                                                  \
  X(parityErrCfebRam, parity_err_cfeb_ram, "CFEB RAM parity error of CFEBs 4:0, latched (6:5: parityErrCfebRamHigh)")  \
  X(parityErrGemRpcMini,                                                                                               \
    parity_err_rpc,                                                                                                    \
    "GEM, RPC or miniscope RAM parity error, latched (CMSSW name parity_err_rpc)")                                     \
  X(parityErrSummary, parity_err_summary, "parity error summary, latched")                                             \
  X(alct0Valid, alct0Valid, "ALCT0 valid")                                                                             \
  X(alct0Quality, alct0Quality, "ALCT0 quality (2 bits)")                                                              \
  X(alct0Amu, alct0Amu, "ALCT0 accelerator-muon flag")                                                                 \
  X(alct0Key, alct0Key, "ALCT0 key wire group")                                                                        \
  X(alctPretrigWin, alct_pretrig_win, "position of the ALCT in the CLCT pre-trigger window (formats 4-6)")             \
  X(clct1Slope, clct1_slope, "CLCT1 CCLUT slope magnitude (formats 7-9)")                                              \
  X(alct1Valid, alct1Valid, "ALCT1 valid")                                                                             \
  X(alct1Quality, alct1Quality, "ALCT1 quality (2 bits)")                                                              \
  X(alct1Amu, alct1Amu, "ALCT1 accelerator-muon flag")                                                                 \
  X(alct1Key, alct1Key, "ALCT1 key wire group")                                                                        \
  X(driftDelay, drift_delay, "CLCT drift delay (configuration)")                                                       \
  X(bcbReadEnable, bcb_read_enable, "blocked CFEB bits in the readout (configuration)")                                \
  X(hsLayerTrig, hs_layer_trig, "layer-mode trigger")                                                                  \
  X(alctBxn, alctBXN, "ALCT0/1 BXN (5 bits; formats 4-6)")                                                             \
  X(hmtNhitsHigh, hmt_nhits_bits_high, "HMT number of hits, bits 6:2 (formats 7-9)")                                   \
  X(alctEccErr, alct_ecc_err, "ALCT trigger-path ECC error code (2 bits)")                                             \
  X(cfebBadbitsFound, cfeb_badbits_found, "CFEBs 4:0 with at least one bad bit (6:5: cfebBadbitsFoundHigh)")           \
  X(cfebBadbitsBlocked, cfeb_badbits_blocked, "a CFEB had bad bits that were blocked")                                 \
  X(alctCfgDone, alctCfg, "ALCT FPGA configuration done")                                                              \
  X(bx0Match, bx0_match, "ALCT BX0 and CLCT BX0 match")                                                                \
  X(mpcDelay, MPCDelay, "MPC transmit delay (configuration)")                                                          \
  X(mpcAccept, MPCAccept, "MPC accept response for the two LCTs (2 bits)")                                             \
  X(cfebsEnabled, CFEBsEnabled, "CFEBs 4:0 enabled for triggering (6:5: cfebsEnabledHigh)")                            \
  X(rpcList,                                                                                                           \
    RPCList,                                                                                                           \
    "RPCs in the readout (formats 4-7; word 36 holds the RPC fields only when the GEM readout is disabled, a "         \
    "register not in the header)")                                                                                     \
  X(nRPCs,                                                                                                             \
    NRPCs,                                                                                                             \
    "number of RPCs in the readout (formats 4-7; word 36 holds the RPC fields only when the GEM readout is disabled, " \
    "a register not in the header)")                                                                                   \
  X(rpcEnable,                                                                                                         \
    RPCEnable,                                                                                                         \
    "RPC readout enabled (configuration; formats 4-7; word 36 holds the RPC fields only when the GEM readout is "      \
    "disabled, a register not in the header)")                                                                         \
  X(fifoTbinsRPC,                                                                                                      \
    fifo_tbins_rpc,                                                                                                    \
    "RPC FIFO time bins read out (configuration; formats 4-7; word 36 holds the RPC fields only when the GEM readout " \
    "is disabled, a register not in the header)")                                                                      \
  X(fifoPretrigRPC,                                                                                                    \
    fifo_pretrig_rpc,                                                                                                  \
    "RPC FIFO time bins before the pre-trigger (configuration; formats 4-7; word 36 holds the RPC fields only when "   \
    "the GEM readout is disabled, a register not in the header)")                                                      \
  X(gemEnabledFibers,                                                                                                  \
    gem_enabled_fibers_,                                                                                               \
    "GEM fibers enabled for readout (4 bits, configuration; formats 8, 9; word 36 holds the GEM fields only when the " \
    "GEM readout is enabled, a register not in the header)")                                                           \
  X(gemZeroSuppression,                                                                                                \
    gem_zero_suppression_,                                                                                             \
    "GEM zero suppression enabled (configuration; formats 8, 9; word 36 holds the GEM fields only when the GEM "       \
    "readout is enabled, a register not in the header)")                                                               \
  X(fifoTbinsGEM,                                                                                                      \
    fifo_tbins_gem_,                                                                                                   \
    "GEM FIFO time bins read out (configuration; formats 8, 9; word 36 holds the GEM fields only when the GEM "        \
    "readout is enabled, a register not in the header)")                                                               \
  X(fifoPretrigGEM,                                                                                                    \
    fifo_pretrig_gem_,                                                                                                 \
    "GEM FIFO time bins before the pre-trigger (configuration; formats 8, 9; word 36 holds the GEM fields only when "  \
    "the GEM readout is enabled, a register not in the header)")                                                       \
  X(rWrBufAdr, r_wr_buf_adr, "buffer RAM write address at the pre-trigger (11 bits)")                                  \
  X(rWrBufReady, r_wr_buf_ready, "write buffer ready at the pre-trigger")                                              \
  X(wrBufReady, wr_buf_ready, "write buffer ready at readout")                                                         \
  X(bufQFull, buf_q_full, "all raw-hits RAM in use")                                                                   \
  X(bufQEmpty, buf_q_empty, "no fences left on the buffer stack")                                                      \
  X(rBufFenceDist, r_buf_fence_dist, "distance to the first fence address at the pre-trigger (11 bits)")               \
  X(bufQOvfErr, buf_q_ovf_err, "buffer stack: push when full")                                                         \
  X(bufQUdfErr, buf_q_udf_err, "buffer stack: pop when empty")                                                         \
  X(bufQAdrErr, buf_q_adr_err, "buffer stack: popped fence address differs from the release address")                  \
  X(bufStalled, buf_stalled, "buffer write pointer hit a fence and stalled (latched)")                                 \
  X(bufFenceCnt, buf_fence_cnt, "number of fences in the fence RAM (12 bits)")                                         \
  X(reverseHsCsc, reverse_hs_csc, "reversed staggered CSC, non-ME1 (configuration)")                                   \
  X(reverseHsMe1a, reverse_hs_me1a, "ME1a half-strip order reversed (configuration)")                                  \
  X(reverseHsMe1b, reverse_hs_me1b, "ME1b half-strip order reversed (configuration)")                                  \
  X(activeCFEBsHigh, activeCFEBs_2, "active-CFEB list sent to the DMB, bits 6:5")                                      \
  X(readCFEBsHigh, readCFEBs_2, "CFEBs read out, bits 6:5")                                                            \
  X(parityErrCfebRamHigh,                                                                                              \
    cfeb_badbits_found_2,                                                                                              \
    "CFEB RAM parity error of CFEBs 6:5, word 40 bits 5:4 (CMSSW name cfeb_badbits_found_2)")                          \
  X(cfebBadbitsFoundHigh,                                                                                              \
    parity_err_cfeb_ram_2,                                                                                             \
    "CFEBs 6:5 with at least one bad bit, word 40 bits 7:6 (CMSSW name parity_err_cfeb_ram_2)")                        \
  X(cfebsEnabledHigh, CFEBsEnabled_2, "CFEBs 6:5 enabled for triggering")                                              \
  X(bufFenceCntIsPeak, buf_fence_cnt_is_peak, "the current fence count is the peak")                                   \
  X(mxcfeb, mxcfeb, "board with 7 CFEBs (formats 4-7; OTMBs in the Run-3 format write the GEM-CSC bend enable)")       \
  X(gemCscBendEnable, gem_csc_bend_enable, "GEM-CSC bending enabled (configuration; formats 8, 9)")                    \
  X(trigSourceMe1aMe1b, trig_source_vec, "pre-trigger was ME1a only (bit 0) or ME1b only (bit 1)")                     \
  X(tmbTrigPulse, tmb_trig_pulse, "TMB trigger pulse coincident with the readout push")                                \
  X(tmbAllowAlct, tmb_allow_alct, "ALCT-only TMB trigger allowed (configuration; formats 4-6)")                        \
  X(tmbAllowClct, tmb_allow_clct, "CLCT-only TMB trigger allowed (configuration; formats 4-6)")                        \
  X(tmbAllowMatch, tmb_allow_match, "ALCT-CLCT match trigger allowed (configuration; formats 4-6)")                    \
  X(tmbAllowAlctRo, tmb_allow_alct_ro, "ALCT-only readout allowed (configuration; formats 4-6)")                       \
  X(tmbAllowClctRo, tmb_allow_clct_ro, "CLCT-only readout allowed (configuration; formats 4-6)")                       \
  X(tmbAllowMatchRo, tmb_allow_match_ro, "ALCT-CLCT match readout allowed (configuration; formats 4-6)")               \
  X(run3TrigDf, run3_trig_df, "Run-3 trigger data format (configuration; formats 7-9)")                                \
  X(gemEnable, gem_enable, "GEM enabled (configuration; formats 7-9)")                                                 \
  X(hmtMatchWin, hmt_match_win, "position of the HMT in its match window (4 bits; formats 7-9)")                       \
  X(tmbAlctOnlyRo, tmb_alct_only_ro, "only the ALCT triggered (non-triggering readout)")                               \
  X(tmbClctOnlyRo, tmb_clct_only_ro, "only the CLCT triggered (non-triggering readout)")                               \
  X(tmbMatchRo, tmb_match_ro, "ALCT and CLCT matched in time (non-triggering readout)")                                \
  X(tmbTrigKeep, tmb_trig_keep, "triggering readout event")                                                            \
  X(tmbNonTrigKeep, tmb_non_trig_keep, "non-triggering readout event")                                                 \
  X(lyrThreshPretrig, lyr_thresh_pretrig, "layer pre-trigger threshold (configuration; formats 4, 5)")                 \
  X(layerTrigEn, layer_trig_en, "layer trigger mode enabled (configuration; formats 4, 5)")                            \
  X(cathodeHmt, cathode_hmt, "cathode HMT bits (2; formats 6-9)")                                                      \
  X(anodeHmt, anode_hmt, "anode HMT bits (2; OTMBs: the ALCT BXN bits 2:1 of the Run-3 ALCT format; formats 6-9)")

#define GCN_TMB_ALL(X, X2) GCN_TMB_COLUMNS(X, X2) GCN_TMB_COLUMNS_27(X, X2)

  // kNWords0: the columns of words 0-26 (from the status digi), then those of words 27-42
  constexpr unsigned kNWords0 = 0
#define GCN_COUNT(name, field, doc) +1
#define GCN_COUNT2(name, a, b, doc) +1
      GCN_TMB_COLUMNS(GCN_COUNT, GCN_COUNT2);
  constexpr unsigned kNDecoded = 0 GCN_TMB_ALL(GCN_COUNT, GCN_COUNT2);
#undef GCN_COUNT
#undef GCN_COUNT2

  struct ColumnSpec {
    const char* name;
    const char* doc;
  };
  const std::array<ColumnSpec, kNDecoded> kColumns = {{
#define GCN_SPEC(name, field, doc) {#name, doc},
#define GCN_SPEC2(name, a, b, doc) {#name, doc},
      GCN_TMB_ALL(GCN_SPEC, GCN_SPEC2)
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
    GCN_TMB_ALL(GCN_GET, GCN_GET2)
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
    GCN_TMB_ALL(GCN_HAS, GCN_HAS2)
#undef GCN_HAS
#undef GCN_HAS2
    return n;
  }
  const Names kExpectedNames = {{
#define GCN_EXP(name, field, doc) 1u,
#define GCN_EXP2(name, a, b, doc) 3u,
      GCN_TMB_ALL(GCN_EXP, GCN_EXP2)
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

  // The unpacker settings that decide which readouts give a TMB status digi (CSCDCCUnpacker.cc:116-145), from the
  // configuration of the module that made the status digis (no default: a missing parameter throws).
  struct UnpackerSettings {
    edm::InputTag raw;
    bool useExaminer, useSelectiveUnpacking, disableMappingCheck, b904Setup;
    unsigned examinerMask, errorMask;
    int b904vmecrate, b904dmb;
  };
  UnpackerSettings unpackerSettings(const edm::ParameterSet& ps) {
    return {ps.getParameter<edm::InputTag>("InputObjects"),
            ps.getParameter<bool>("UseExaminer"),
            ps.getParameter<bool>("UseSelectiveUnpacking"),
            ps.getUntrackedParameter<bool>("DisableMappingCheck"),
            ps.getUntrackedParameter<bool>("B904Setup"),
            ps.getParameter<unsigned>("ExaminerMask"),
            ps.getParameter<unsigned>("ErrorMask"),
            ps.getUntrackedParameter<int>("B904vmecrate"),
            ps.getUntrackedParameter<int>("B904dmb")};
  }

  struct Readout {
    std::vector<unsigned short> header;  // CSCTMBHeader::sizeInWords() words, as read
    std::array<unsigned short, kStoredTrailerWords> trailer;
  };

  // The TMB headers of the chamber readouts that give a TMB status digi, by chamber, in unpacking order: the steps of
  // CSCDCCUnpacker::produce with the same CMSSW classes: FED list and size (CSCDCCUnpacker.cc:299-323), examiner
  // (:325-366), selective unpacking (:399-428), DDU error mask (:451-456), crate/DMB range and chamber id (:475-490),
  // FED/DDU mapping check (:494-509), goodTMB (:552-558), status digi (:667-670).
  std::map<CSCDetId, std::vector<Readout>> tmbReadouts(const FEDRawDataCollection& raw,
                                                       const CSCCrateMap& crateMap,
                                                       const CSCChamberMap& chamberMap,
                                                       const UnpackerSettings& u) {
    constexpr unsigned long kDccBinCheckMask = 0x06080016;  // CSCDCCUnpacker.cc:291
    constexpr std::array<unsigned, 36> kPostLS1Map = {{841, 842, 843, 844, 845, 846, 847, 848, 849, 831, 832, 833,
                                                       834, 835, 836, 837, 838, 839, 861, 862, 863, 864, 865, 866,
                                                       867, 868, 869, 851, 852, 853, 854, 855, 856, 857, 858, 859}};
    std::vector<unsigned> feds;
    for (unsigned id = FEDNumbering::MINCSCFEDID; id <= FEDNumbering::MAXCSCFEDID; ++id)
      feds.push_back(id);
    for (unsigned id = FEDNumbering::MINCSCDDUFEDID; id <= FEDNumbering::MAXCSCDDUFEDID; ++id)
      feds.push_back(id);
    std::map<CSCDetId, std::vector<Readout>> out;
    for (const unsigned id : feds) {
      const bool isDDU = id >= FEDNumbering::MINCSCDDUFEDID && id <= FEDNumbering::MAXCSCDDUFEDID;
      const FEDRawData& fed = raw.FEDData(id);
      if (fed.size() < 32)
        continue;
      std::unique_ptr<CSCDCCExaminer> examiner;
      if (u.useExaminer) {
        examiner = std::make_unique<CSCDCCExaminer>();
        if (u.examinerMask & 0x40000)
          examiner->crcCFEB(true);
        if (u.examinerMask & 0x8000)
          examiner->crcTMB(true);
        if (u.examinerMask & 0x0400)
          examiner->crcALCT(true);
        examiner->setMask(u.examinerMask);
        if (isDDU)
          examiner->modeDDU(true);
        const uint16_t* data = reinterpret_cast<const uint16_t*>(fed.data());  // check() advances it
        const int res = examiner->check(data, long(fed.size() / 2));
        if (res < 0 || (examiner->errors() & (u.useSelectiveUnpacking ? kDccBinCheckMask : u.examinerMask)))
          continue;
      }
      CSCDCCExaminer* ptrExaminer = u.useSelectiveUnpacking ? examiner.get() : nullptr;
      const uint16_t* buf = reinterpret_cast<const uint16_t*>(fed.data());
      std::vector<CSCDDUEventData> ddus;
      if (isDDU)
        ddus.push_back(CSCDDUEventData(buf, ptrExaminer));
      else
        ddus = CSCDCCEventData(buf, ptrExaminer).dduData();
      for (const CSCDDUEventData& ddu : ddus) {
        if (ddu.trailer().errorstat() & u.errorMask)
          continue;
        for (const CSCEventData& csc : ddu.cscData()) {
          const int vme = u.b904Setup ? u.b904vmecrate : csc.dmbHeader()->crateID();
          const int dmb = u.b904Setup ? u.b904dmb : csc.dmbHeader()->dmbID();
          if (!(vme >= 1 && vme <= 60 && dmb >= 1 && dmb <= 10 && dmb != 6))
            continue;
          const CSCDetId chamber = crateMap.detId(vme, dmb, 0, 0);
          if (isDDU) {
            unsigned dduid = chamberMap.ddu(chamber);
            if (dduid >= 1 && dduid <= 36)
              dduid = kPostLS1Map[dduid - 1];
            if (!u.disableMappingCheck && !u.b904Setup && id != dduid)
              continue;
          }
          if (!(csc.dmbHeader()->nclct() && csc.tmbData() && csc.tmbHeader()->check() && csc.comparatorData()->check()))
            continue;
          Readout r;
          CSCTMBHeader* header = csc.tmbHeader();
          r.header.assign(header->data(), header->data() + header->sizeInWords());
          const unsigned short* trailer = csc.tmbData()->tmbTrailer()->data();
          std::copy(trailer, trailer + kStoredTrailerWords, r.trailer.begin());
          out[chamber].push_back(r);
        }
      }
    }
    return out;
  }

  const char* const kFullHeaderMatchDoc =
      "words 27-42 of this row: 0 = the second decoding of the raw data has exactly one readout of the chamber, with "
      "words 0-26 and the trailer equal to the status digi's (words 27-42 filled); 1 = no readout of the chamber; 2 = "
      "more than one readout or status digi of the chamber; 3 = one, but words 0-26 or the trailer differ; 4 = the "
      "header has 27 words (2006 format)";

}  // namespace

class TMBHeaderTableProducer : public edm::global::EDProducer<> {
public:
  explicit TMBHeaderTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  const edm::EDGetTokenT<CSCTMBStatusDigiCollection> token_;
  const edm::EDGetTokenT<FEDRawDataCollection> rawToken_;
  const edm::ESGetToken<CSCCrateMap, CSCCrateMapRcd> crateMapToken_;
  const edm::ESGetToken<CSCChamberMap, CSCChamberMapRcd> chamberMapToken_;
  const std::string label_;
  const edm::InputTag rawTag_;
};

TMBHeaderTableProducer::TMBHeaderTableProducer(const edm::ParameterSet& ps)
    : token_(consumes(ps.getParameter<edm::InputTag>("src"))),
      rawToken_(consumes(ps.getParameter<edm::InputTag>("raw"))),
      crateMapToken_(esConsumes()),
      chamberMapToken_(esConsumes()),
      label_(ps.getParameter<edm::InputTag>("src").encode()),
      rawTag_(ps.getParameter<edm::InputTag>("raw")) {
  // Self-tests: with words 0-26 zero and words 27-42 all ones, every column of words 0-26 reads 0 (or -1 if absent),
  // i.e. reads no word the status digi does not keep, and every column of words 27-42 reads a non-zero value (or -1);
  // every column exists in at least one format (a field name missing from every format's `bits` would silently give
  // -1).
  std::array<bool, kNDecoded> present{};
  for (int code = 4; code <= 9; ++code) {
    std::array<unsigned short, kHeaderWords> buf{};
    for (unsigned w = kStoredHeaderWords; w < kHeaderWords; ++w)
      buf[w] = 0xFFFF;
    const Values v = decodeFormat(code, buf.data());
    for (unsigned i = 0; i < kNDecoded; ++i) {
      if (i < kNWords0 && v[i] != 0 && v[i] != -1)
        throw cms::Exception("LogicError") << "TMBHeaderTableProducer: column " << kColumns[i].name << " of format "
                                           << code << " reads beyond header word " << kStoredHeaderWords - 1;
      if (i >= kNWords0 && v[i] == 0)
        throw cms::Exception("LogicError")
            << "TMBHeaderTableProducer: column " << kColumns[i].name << " of format " << code
            << " reads no bit of header words " << kStoredHeaderWords << "-" << kHeaderWords - 1;
      present[i] = present[i] || v[i] != -1;
    }
  }
  for (unsigned i = 0; i < kNDecoded; ++i)
    if (!present[i])
      throw cms::Exception("LogicError") << "TMBHeaderTableProducer: column " << kColumns[i].name
                                         << " is in no header format";
  Names names{};
  const std::array<Names, 6> perFormat = {{namesPresent<CSCTMBHeader2013>(),
                                           namesPresent<CSCTMBHeader2020_Run2>(),
                                           namesPresent<CSCTMBHeader2020_TMB>(),
                                           namesPresent<CSCTMBHeader2020_CCLUT>(),
                                           namesPresent<CSCTMBHeader2020_GEM>(),
                                           namesPresent<CSCTMBHeader2020_GEM_Run3b>()}};
  for (int code = 4; code <= 9; ++code) {
    // the decoder gives a value exactly for the columns whose field name the format has
    const Names& n = perFormat[code - 4];
    std::array<unsigned short, kHeaderWords> buf{};
    const Values v = decodeFormat(code, buf.data());
    for (unsigned i = 0; i < kNDecoded; ++i) {
      if ((v[i] != -1) != (n[i] != 0))
        throw cms::Exception("LogicError") << "TMBHeaderTableProducer: column " << kColumns[i].name << " of format "
                                           << code << ": decoded presence differs from its field names";
      names[i] |= n[i];
    }
  }
  for (unsigned i = 0; i < kNDecoded; ++i)
    if (names[i] != kExpectedNames[i])
      throw cms::Exception("LogicError") << "TMBHeaderTableProducer: a field name of column " << kColumns[i].name
                                         << " is in no header format";
  produces<nanoaod::FlatTable>("header");
  produces<nanoaod::FlatTable>("raw");
  produces<nanoaod::FlatTable>("unmatched");
}

void TMBHeaderTableProducer::produce(edm::StreamID, edm::Event& event, const edm::EventSetup& setup) const {
  const auto digiHandle = event.getHandle(token_);
  const auto& digis = *digiHandle;
  // the unpacker's untracked settings are in the configuration of the current process only (files keep the tracked
  // part): the status digis must be made in this job
  if (digiHandle.provenance()->processName() != moduleDescription().processName())
    throw cms::Exception("Configuration")
        << "TMBHeaderTableProducer: " << label_ << " was made by process " << digiHandle.provenance()->processName()
        << ", not by this job: rerun the CSC unpacker (UnpackStatusDigis = True)";
  const UnpackerSettings unpacker =
      unpackerSettings(edm::parameterSet(digiHandle.provenance()->stable(), event.processHistory()));
  if (!(unpacker.raw == rawTag_))
    throw cms::Exception("Configuration") << "TMBHeaderTableProducer: raw = " << rawTag_.encode()
                                          << " but the module of " << label_ << " unpacks " << unpacker.raw.encode();
  const auto readouts =
      tmbReadouts(event.get(rawToken_), setup.getData(crateMapToken_), setup.getData(chamberMapToken_), unpacker);
  struct Row {
    CSCDetId chamber;
    std::array<unsigned short, kHeaderWords> header;  // words 0-26 of the status digi, 27-42 of the second decoding
    std::array<unsigned short, kStoredTrailerWords> trailer;
    int format;
    int revision;
    int fullHeaderMatch;
    Values values;
  };
  std::map<CSCDetId, unsigned> digisPerChamber;
  for (auto detUnit = digis.begin(); detUnit != digis.end(); ++detUnit)
    for (auto digi = (*detUnit).second.first; digi != (*detUnit).second.second; ++digi)
      ++digisPerChamber[(*detUnit).first];
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
      const auto readout = readouts.find(r.chamber);
      const unsigned nReadouts = readout == readouts.end() ? 0 : readout->second.size();
      if (nReadouts == 0)
        r.fullHeaderMatch = 1;
      else if (nReadouts > 1 || digisPerChamber.at(r.chamber) > 1)
        r.fullHeaderMatch = 2;
      else {
        const Readout& ro = readout->second.front();
        if (ro.header.size() < kStoredHeaderWords ||
            !std::equal(r.header.begin(), r.header.begin() + kStoredHeaderWords, ro.header.begin()) ||
            ro.trailer != r.trailer)
          r.fullHeaderMatch = 3;
        else if (ro.header.size() != kHeaderWords)
          r.fullHeaderMatch = 4;
        else {
          r.fullHeaderMatch = 0;
          std::copy(ro.header.begin() + kStoredHeaderWords, ro.header.end(), r.header.begin() + kStoredHeaderWords);
        }
      }
      r.format = formatCode(r.header.data());
      // the dispatch must agree with CMSSW's own choice of the firmware version
      if (r.format != 0 && CSCTMBHeader(r.header.data()).FirmwareVersion() != firmwareVersionOf(r.format))
        throw cms::Exception("LogicError") << "TMBHeaderTableProducer: format " << r.format << " but CSCTMBHeader "
                                           << CSCTMBHeader(r.header.data()).FirmwareVersion();
      r.revision = r.header[0] == 0xDB0C ? int(CSCTMBHeader2007(r.header.data()).firmwareRevision()) : -1;
      r.values = decodeFormat(r.format, r.header.data());
      if (r.fullHeaderMatch != 0)
        for (unsigned i = kNWords0; i < kNDecoded; ++i)
          if (r.values[i] != -1)
            r.values[i] = -2;
      rows.push_back(r);
    }
  }
  std::set<CSCDetId> withDigi;
  for (const Row& r : rows)
    withDigi.insert(r.chamber);
  std::vector<std::pair<CSCDetId, unsigned>> unmatched;
  for (const auto& [chamber, list] : readouts)
    if (!withDigi.count(chamber))
      unmatched.emplace_back(chamber, list.size());

  auto table = std::make_unique<nanoaod::FlatTable>(rows.size(), "TMBHeader", false);
  table->setDoc("(O)TMB DAQ header of each chamber readout: every status digi of " + label_ +
                ", in collection order; header words 0-26 from the status digi, words 27-42 from a second decoding "
                "of " +
                rawTag_.encode() +
                " (fullHeaderMatch); -1 where the format lacks a field; MC headers are "
                "written by the packer (most fields 0)");
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
      "minor rev&0x1F); -1 for the 2006 format and format 0");
  addColumn<int16_t>(
      *table, "fullHeaderMatch", rows, [](const Row& r) { return r.fullHeaderMatch; }, kFullHeaderMatchDoc);
  for (unsigned i = 0; i < kNDecoded; ++i)
    addColumn<int16_t>(
        *table,
        kColumns[i].name,
        rows,
        [i](const Row& r) { return r.values[i]; },
        i < kNWords0 ? std::string(kColumns[i].doc)
                     : std::string(kColumns[i].doc) + " [words 27-42: -2 unless fullHeaderMatch = 0]");

  auto raw = std::make_unique<nanoaod::FlatTable>(rows.size(), "TMBHeaderRaw", false);
  raw->setDoc(
      "raw words of the TMBHeader rows (same rows): header words 0-26 and the 7 trailer words of the status digi, "
      "header words 27-42 of the second decoding (0 unless TMBHeader_fullHeaderMatch = 0)");
  for (unsigned w = 0; w < kHeaderWords; ++w)
    addColumn<uint16_t>(
        *raw,
        "header" + std::to_string(w),
        rows,
        [w](const Row& r) { return r.header[w]; },
        "header word " + std::to_string(w) +
            (w < kStoredHeaderWords ? "" : " (second decoding; 0 unless TMBHeader_fullHeaderMatch = 0)"));
  for (unsigned w = 0; w < kStoredTrailerWords; ++w)
    addColumn<uint16_t>(
        *raw,
        "trailer" + std::to_string(w),
        rows,
        [w](const Row& r) { return r.trailer[w]; },
        "trailer word " + std::to_string(w) +
            " from 0x6E0C (optional 0x2AAA 0x5555 filler, 0xDE0F, 2 CRC words, word count; without filler the last 2 "
            "words lie beyond the TMB record)");

  auto unmatchedTable = std::make_unique<nanoaod::FlatTable>(unmatched.size(), "TMBHeaderUnmatched", false);
  unmatchedTable->setDoc("chambers whose TMB header the second decoding of " + rawTag_.encode() +
                         " finds but that have no status digi in " + label_ +
                         " (normally none: the decoding replicates the unpacker)");
  using U = std::pair<CSCDetId, unsigned>;
  addColumn<int16_t>(*unmatchedTable, "endcap", unmatched, [](const U& u) { return u.first.zendcap(); }, "+1 or -1");
  addColumn<int16_t>(*unmatchedTable, "station", unmatched, [](const U& u) { return u.first.station(); }, "station");
  addColumn<int16_t>(*unmatchedTable, "ring", unmatched, [](const U& u) { return u.first.ring(); }, "ring");
  addColumn<int16_t>(*unmatchedTable, "chamber", unmatched, [](const U& u) { return u.first.chamber(); }, "chamber");
  addColumn<int16_t>(
      *unmatchedTable, "nReadouts", unmatched, [](const U& u) { return u.second; }, "readouts of the chamber found");

  event.put(std::move(table), "header");
  event.put(std::move(raw), "raw");
  event.put(std::move(unmatchedTable), "unmatched");
}

void TMBHeaderTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("CSCTMBStatusDigiCollection (muonCSCDigis:MuonCSCTMBStatusDigi)");
  desc.add<edm::InputTag>("raw")->setComment(
      "FEDRawDataCollection for words 27-42: must be the InputObjects of the module that made src (checked)");
  descriptions.add("tmbHeaderTableProducer", desc);
}

DEFINE_FWK_MODULE(TMBHeaderTableProducer);
