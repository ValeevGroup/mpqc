#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class ReferencePeriodicFourCenterFockBuilder : public PeriodicFockBuilder<Tile, Policy> {
 public:
	using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;

	ReferencePeriodicFourCenterFockBuilder(Factory &ao_factory)
			: ao_factory_(ao_factory) {}

	~ReferencePeriodicFourCenterFockBuilder() {}

	array_type operator()(array_type const &D, double) override {
		// feed density matrix to Factory
		ao_factory_.set_density(D);

		array_type J, K, G;
		J = ao_factory_.compute_direct(L"(μ ν| J|κ λ)");
		K = ao_factory_.compute_direct(L"(μ ν| K|κ λ)");

		G("mu, nu") = 2.0 * J("mu, nu") - K("mu, nu");

		return G;
	}

	void register_fock(const array_type &fock,
										 FormulaRegistry<array_type> &registry) override {
		registry.insert(Formula(L"(κ|F|λ)"), fock);
	}

 private:
	Factory &ao_factory_;
};

template <typename Tile, typename Policy>
class PeriodicFourCenterFockBuilder : public PeriodicFockBuilder<Tile, Policy> {
public:
	using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;

	PeriodicFourCenterFockBuilder(madness::World &world,
																std::shared_ptr<const Basis> basis,
																bool compute_J, bool compute_K,
																std::string screen = "schwarz",
																double screen_threshold = 1.0e-10)
		:  {}

	~ReferencePeriodicFourCenterFockBuilder() {}

	array_type operator()(array_type const &D, double) override {

	}

	void register_fock(const array_type &fock,
										 FormulaRegistry<array_type> &registry) override {
		registry.insert(Formula(L"(κ|F|λ)"), fock);
	}

private:


};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_
