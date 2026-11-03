#ifndef L1Trigger_GemCscNano_CSCLCTTable_h
#define L1Trigger_GemCscNano_CSCLCTTable_h

// The GemCscNano LCT columns: every CSCCorrelatedLCTDigi of a collection, in collection order (detector unit, then
// digi), with its accessors as stored. Used for the CSC-unpacked LCTs (table LCT) and for the LCTs the EMTF unpacker
// builds from the EMTF readout (table EMTFLCT), so that both tables have the same columns.

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <map>
#include <memory>
#include <string>

namespace gemcscnano {

  // docOverrides: column name -> doc, replacing the default doc of that column (the column must exist)
  std::unique_ptr<nanoaod::FlatTable> makeLCTTable(const CSCCorrelatedLCTDigiCollection& lcts,
                                                   const std::string& name,
                                                   const std::string& doc,
                                                   const std::map<std::string, std::string>& docOverrides = {});

}  // namespace gemcscnano

#endif
