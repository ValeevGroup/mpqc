#include <mpqc/chemistry/qc/scf/scf.h>
#include <mpqc/chemistry/qc/integrals/integrals.h>
#include "../../../../../utility/time.h"
#include <mpqc/chemistry/qc/scf/diagonalize_for_coffs.hpp>

namespace mpqc {
namespace scf {

double ClosedShellSCF::energy() const {
    return repulsion_
           + D_("i,j").dot(F_("i,j") + H_("i,j"), D_.get_world()).get();
}

bool ClosedShellSCF::solve(int64_t max_iters, double thresh) {
    auto &world = F_.get_world();

    if (world.rank() == 0) {
        std::cout << "Starting SCF:\n"
                  << "\tThreshold: " << thresh << "\n"
                  << "\tMaximum number of iterations: " << max_iters << "\n";
    }

    auto iter = 0;
    auto error = std::numeric_limits<double>::max();
    auto rms_error = std::numeric_limits<double>::max();
    auto old_energy = energy();
    const double volume = F_.trange().elements().volume();

    while (iter < max_iters && (thresh < (error / old_energy)
                                || thresh < (rms_error / volume))) {

        auto s0 = mpqc_time::fenced_now(world);

        build_F();
        auto b1 = mpqc_time::fenced_now(world);

        auto current_energy = energy();
        error = std::abs(old_energy - current_energy);
        old_energy = current_energy;

        array_type Grad;
        Grad("i,j") = F_("i,k") * D_("k,l") * S_("l,j")
                      - S_("i,k") * D_("k,l") * F_("l,j");

        rms_error = Grad("i,j").norm().get();

        F_diis_ = F_;
        diis_.extrapolate(F_diis_, Grad);

        auto d0 = mpqc_time::fenced_now(world);
        compute_density();
        auto s1 = mpqc_time::fenced_now(world);

        scf_times_.push_back(mpqc_time::duration_in_s(s0, s1));
        d_times_.push_back(mpqc_time::duration_in_s(d0, s1));
        build_times_.push_back(mpqc_time::duration_in_s(s0, b1));

        if (world.rank() == 0) {
            std::cout << "iteration: " << iter << "\n"
                      << "\tEnergy: " << old_energy << "\n"
                      << "\tabs(Energy Change)/energy: "
                      << (error / std::abs(old_energy)) << "\n"
                      << "\t(Gradient Norm)/n^2: " << (rms_error / volume)
                      << "\n"
                      << "\tScf Time: " << scf_times_.back() << "\n"
                      << "\t\tDensity Time: " << d_times_.back() << "\n"
                      << "\t\tFock Build Time: " << build_times_.back() << "\n";
        }
        f_builder_->print_iter("\t\t");
        d_builder_->print_iter("\t\t");
        ++iter;
    }

    if (iter == max_iters) {
        return false;
    } else {
        return true;
    }
}

void ClosedShellSCF::compute_density() {
    auto dc_pair = d_builder_->operator()(F_diis_);
    D_ = dc_pair.first;
    C_ = dc_pair.second;
}

void ClosedShellSCF::build_F() {
    F_("i,j") = H_("i,j") + f_builder_->operator()(D_, C_)("i,j");
}

rapidjson::Value ClosedShellSCF::results(rapidjson::Document &d) const {
    rapidjson::Value scf_object(rapidjson::kObjectType);
    scf_object.AddMember("Type", "ClosedShellSCF", d.GetAllocator());
    scf_object.AddMember("Energy", energy(), d.GetAllocator());

    auto avg = [](std::vector<double> const &v) {
        auto sum = 0.0;
        for (auto d : v) {
            sum += d;
        }
        return sum / double(v.size());
    };

    scf_object.AddMember("Avg Scf Time", avg(scf_times_), d.GetAllocator());
    scf_object.AddMember("Avg Density Build Time", avg(d_times_),
                         d.GetAllocator());
    scf_object.AddMember("Density Builder", d_builder_->results(d),
                         d.GetAllocator());
    scf_object.AddMember("Avg Fock Build Time", avg(build_times_),
                         d.GetAllocator());
    scf_object.AddMember("Fock Builder", f_builder_->results(d),
                         d.GetAllocator());

    return scf_object;
}

} // namespace scf
} // namespace mpqc
