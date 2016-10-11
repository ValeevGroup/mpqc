#pragma once
#ifndef MPQC_SCF_BUILDER_H
#define MPQC_SCF_BUILDER_H

#include "../../../../../include/tiledarray.h"
#include "../../../../../common/namespaces.h"

#include "../../../../../utility/json_handling.h"
#include <mpqc/chemistry/qc/expression/formula_registry.h>

#include <string>

namespace mpqc {
namespace scf {

class FockBuilder {
  public:
    using array_type = TA::TSpArrayD;
    virtual ~FockBuilder() = default;

    virtual array_type operator()(array_type const&, array_type const&) = 0;

    virtual void print_iter(std::string const &) = 0;

    virtual void register_fock(const TA::TSpArrayD &, FormulaRegistry<TA::TSpArrayD>& ) = 0;

    virtual rapidjson::Value results(rapidjson::Document &) = 0;
};

} // namespace scf
} // namespace mpqc


#endif // MPQC_SCF_BUILDER_H
