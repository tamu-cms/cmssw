#include "L1Trigger/GemCscNano/interface/CSCLCTTable.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "L1Trigger/GemCscNano/interface/TableColumns.h"

#include <cstdint>
#include <set>
#include <vector>

namespace gemcscnano {

  std::unique_ptr<nanoaod::FlatTable> makeLCTTable(const CSCCorrelatedLCTDigiCollection& lctCollection,
                                                   const std::string& name,
                                                   const std::string& doc,
                                                   const std::map<std::string, std::string>& docOverrides) {
    struct LCTRow {
      CSCDetId chamber;
      const CSCCorrelatedLCTDigi* lct;
    };
    std::vector<LCTRow> lcts;
    for (auto detUnit = lctCollection.begin(); detUnit != lctCollection.end(); ++detUnit) {
      const auto& range = (*detUnit).second;
      for (auto lct = range.first; lct != range.second; ++lct)
        lcts.push_back({(*detUnit).first, &*lct});
    }
    std::set<std::string> used;
    auto d = [&](const std::string& column, const std::string& standard) {
      auto it = docOverrides.find(column);
      if (it == docOverrides.end())
        return standard;
      used.insert(column);
      return it->second;
    };

    auto lctTable = std::make_unique<nanoaod::FlatTable>(lcts.size(), name, false);
    lctTable->setDoc(doc);
    addColumn<int16_t>(
        *lctTable, "endcap", lcts, [](const LCTRow& m) { return m.chamber.zendcap(); }, d("endcap", "+1 or -1"));
    addColumn<bool>(
        *lctTable,
        "valid",
        lcts,
        [](const LCTRow& m) { return m.lct->isValid(); },
        d("valid", "CSCCorrelatedLCTDigi::isValid"));
    addColumn<int16_t>(
        *lctTable, "station", lcts, [](const LCTRow& m) { return m.chamber.station(); }, d("station", "station"));
    addColumn<int16_t>(
        *lctTable, "ring", lcts, [](const LCTRow& m) { return m.chamber.ring(); }, d("ring", "ring (ME1/1 LCTs: 1)"));
    addColumn<int16_t>(
        *lctTable, "chamber", lcts, [](const LCTRow& m) { return m.chamber.chamber(); }, d("chamber", "chamber"));
    addColumn<int16_t>(
        *lctTable,
        "hs",
        lcts,
        [](const LCTRow& m) { return m.lct->getStrip(); },
        d("hs", "key half-strip, 0-based (ME1/1: ME1b 0-127, ME1a 128-223)"));
    addColumn<int16_t>(
        *lctTable,
        "es",
        lcts,
        [](const LCTRow& m) { return m.lct->getStrip(8); },
        d("es", "key eighth-strip, getStrip(8)"));
    addColumn<int16_t>(
        *lctTable, "wg", lcts, [](const LCTRow& m) { return m.lct->getKeyWG(); }, d("wg", "key wire group, 0-based"));
    addColumn<int16_t>(
        *lctTable, "quality", lcts, [](const LCTRow& m) { return m.lct->getQuality(); }, d("quality", "quality"));
    addColumn<int16_t>(
        *lctTable, "pattern", lcts, [](const LCTRow& m) { return m.lct->getPattern(); }, d("pattern", "Run-2 pattern"));
    addColumn<int16_t>(
        *lctTable,
        "run3Pattern",
        lcts,
        [](const LCTRow& m) { return m.lct->getRun3Pattern(); },
        d("run3Pattern", "Run-3 pattern (0 for Run3-b unpacked data)"));
    addColumn<int16_t>(
        *lctTable,
        "slope",
        lcts,
        [](const LCTRow& m) { return m.lct->getSlope(); },
        d("slope", "slope, getSlope() (4 bits)"));
    addColumn<int16_t>(
        *lctTable,
        "slope6",
        lcts,
        [](const LCTRow& m) { return m.lct->getSlopeEx(6); },
        d("slope6", "slope at 6-bit resolution, getSlopeEx(6)"));
    addColumn<int16_t>(
        *lctTable, "bend", lcts, [](const LCTRow& m) { return m.lct->getBend(); }, d("bend", "bend (L/R)"));
    addColumn<int16_t>(
        *lctTable, "bx", lcts, [](const LCTRow& m) { return m.lct->getBX(); }, d("bx", "BX (1 bit in unpacked data)"));
    addColumn<int16_t>(
        *lctTable,
        "trknmb",
        lcts,
        [](const LCTRow& m) { return m.lct->getTrknmb(); },
        d("trknmb", "track number (1 or 2)"));
    addColumn<int16_t>(
        *lctTable,
        "type",
        lcts,
        [](const LCTRow& m) { return m.lct->getType(); },
        d("type", "LCT type (CSCCorrelatedLCTDigi::Type; emulator only)"));
    addColumn<int16_t>(
        *lctTable,
        "gemLayer",
        lcts,
        [](const LCTRow& m) { return m.lct->getGemLayerUsedForSlopeComputation(); },
        d("gemLayer", "GEM layer used for the slope: emulator 1, 2, 0 = none; Run3-b unpacked data the raw 1-bit field"));
    for (const auto& [column, text] : docOverrides)
      if (!used.count(column))
        throw cms::Exception("Configuration") << "makeLCTTable: doc override for unknown column " << column;
    return lctTable;
  }

}  // namespace gemcscnano
