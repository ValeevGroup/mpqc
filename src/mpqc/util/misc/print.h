/*
 * print.h
 *
 *  Created on: Nov 1, 2016
 *      Author: evaleev
 */

#ifndef MPQC4_SRC_MPQC_UTIL_MISC_PRINT_H_
#define MPQC4_SRC_MPQC_UTIL_MISC_PRINT_H_

#include <vector>

#include "mpqc/chemistry/units/units.h"
#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/util/misc/exenv.h"

namespace mpqc {
namespace util {

/// print progress
void print_progress(std::size_t lowprogress, std::size_t upprogress,
                    std::vector<std::size_t> &progress_points);

/// print excitation energy at each iteration
template <typename T>
inline void print_davidson_energy_iteration(std::size_t iter,
                                            const EigenVector<T> &delta_e,
                                            const EigenVector<T> &error,
                                            const EigenVector<T> &eig,
                                            double time1, double time2) {
  ExEnv::out0() << mpqc::printf("%4s %3i \t %s %10.1f \t %s %10.1f\n",
                                "iter=", iter, "total time/s=", time1 + time2,
                                "davidson time/s=", time2);

  ExEnv::out0() << mpqc::printf("%4s \t %10s \t %10s \t %15s \n", "root",
                                "deltaE", "error", "energy");
  const std::size_t size = eig.size();
  for (std::size_t i = 0; i < size; i++) {
    ExEnv::out0() << mpqc::printf("%4i \t %10.5e \t %10.5e \t %15.12f \n",
                                  i + 1, delta_e(i), error(i), eig(i));
  }

  ExEnv::out0() << "\n";
}

/// print excitation energy at each iteration
template <typename T>
inline void print_single_state_davidson_energy_iteration(
    std::size_t iter, const T &delta_e, const T &error, const T &eig,
    double time1, double time2) {
  if (iter == 0) {
    ExEnv::out0() << mpqc::printf(
        "%3s \t %10s \t %10s \t %15s \t %10s \t %10s \n", "iter", "deltaE",
        "residual", "energy", "total time/s", "davidson time/s");
  }
  ExEnv::out0() << mpqc::printf(
      "%3i \t %10.5e \t %10.5e \t %15.12f \t %10.1f \t %10.1f \n", iter,
      delta_e, error, eig, time1, time2);
}

/// print dominant index and element in cis vector, with index ordered as i,a
template <typename T>
inline void print_cis_dominant_elements(
    std::vector<std::pair<T, std::vector<std::size_t>>> &dominants) {
  for (const auto &x : dominants) {
    auto &value = x.first;
    auto &index = x.second;
    ExEnv::out0() << mpqc::printf("%5i \t > %5i \t : \t %15.12f \n", index[0],
                                  index[1], value);
  }
}

/// print dominant index and element in t1, with index ordered as a, i
template <typename T>
inline void print_t1_dominant_elements(
    std::vector<std::pair<T, std::vector<std::size_t>>> &dominants) {
  for (const auto &x : dominants) {
    auto &value = x.first;
    auto &index = x.second;
    ExEnv::out0() << mpqc::printf("%5i \t > %5i \t : \t %15.12f \n", index[1],
                                  index[0], value);
  }
}

/// print dominant index and element in t2, with index ordered as a,b,i,j
template <typename T>
inline void print_t2_dominant_elements(
    std::vector<std::pair<T, std::vector<std::size_t>>> &dominants) {
  for (const auto &x : dominants) {
    auto &value = x.first;
    auto &index = x.second;
    ExEnv::out0() << mpqc::printf(
        "%5i \t %5i \t > %5i \t %5i \t : \t %15.12f \n", index[2], index[3],
        index[0], index[1], value);
  }
}

/// print excitation energy with different unit
template <typename T>
inline void print_excitation_energy(const EigenVector<T> &eig, bool triplets) {
  const auto &unit_factory = UnitFactory::get_default();
  const auto Hartree_to_eV = unit_factory->make_unit("eV").from_atomic_units();
  const auto Hartree_to_wavenumber =
      unit_factory->make_unit("energy_wavenumber[cm]").from_atomic_units();

  ExEnv::out0() << "Excitation Energy: ( "
                << (triplets ? "Triplets" : "Singlets") << " )\n";

  ExEnv::out0() << mpqc::printf("%5s \t %10s \t %10s \t %10s \n", "state", "au",
                                "eV", "cm^-1");

  const std::size_t size = eig.size();

  for (std::size_t i = 1; i <= size; i++) {
    T e = eig[i - 1];
    ExEnv::out0() << mpqc::printf("%5i \t %10.8f \t %10.5f \t %10.2f \n", i, e,
                                  e * Hartree_to_eV, e * Hartree_to_wavenumber);
  }
  ExEnv::out0() << "\n";
}
}  // namespace  util
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_UTIL_MISC_PRINT_H_
