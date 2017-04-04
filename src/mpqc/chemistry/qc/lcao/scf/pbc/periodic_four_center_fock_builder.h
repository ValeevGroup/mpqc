#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicFourCenterFockBuilder : public PeriodicFockBuilder<Tile, Policy> {
 public:
	using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;

	PeriodicFourCenterFockBuilder(Factory &ao_factory)
			: ao_factory_(std::make_shared<Factory>(ao_factory)) {}

	~PeriodicFourCenterFockBuilder() {}

	array_type operator()(array_type const &D, array_type const &,
												double) override {
		// feed density matrix to Factory
		ao_factory_->set_density(D);

		array_type J, K, G;
		J = ao_factory_->compute_direct(L"(μ ν| J|κ λ)");
		K = ao_factory_->compute_direct(L"(μ ν| K|κ λ)");

		G("mu, nu") = 2.0 * J("mu, nu") - K("mu, nu");

		return G;
	}

	void register_fock(const array_type &fock,
										 FormulaRegistry<array_type> &registry) override {
		registry.insert(Formula(L"(κ|F|λ)"), fock);
	}

 private:
	std::shared_ptr<Factory> ao_factory_;
}

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_
