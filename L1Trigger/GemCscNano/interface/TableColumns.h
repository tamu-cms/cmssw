#ifndef L1Trigger_GemCscNano_TableColumns_h
#define L1Trigger_GemCscNano_TableColumns_h

// One NanoAOD column computed from each row of a vector, shared by the GemCscNano table producers.

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <string>
#include <vector>

namespace gemcscnano {

  // mantissaBits: NanoAOD float precision (-1: full precision, the FlatTable default)
  template <typename T, typename Row, typename F>
  void addColumn(nanoaod::FlatTable& table,
                 const std::string& name,
                 const std::vector<Row>& rows,
                 F value,
                 const std::string& doc,
                 int mantissaBits = -1) {
    std::vector<T> v;
    v.reserve(rows.size());
    for (const Row& r : rows)
      v.push_back(value(r));
    table.addColumn<T>(name, v, doc, mantissaBits);
  }

}  // namespace gemcscnano

#endif
