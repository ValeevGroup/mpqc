//
// Created by Chong Peng on 10/5/16.
//

#ifndef MPQC_CHEMISTRY_QC_SCF_RHF_H
#define MPQC_CHEMISTRY_QC_SCF_RHF_H


#include <tiledarray.h>
#include "mpqc/util/misc/json_handling.h"
#include <mpqc/chemistry/qc/wfn/ao_wfn.h>
#include <mpqc/chemistry/qc/scf/builder.h>
#include <mpqc/chemistry/qc/scf/density_builder.h>

/**
 *  RHF Class of AOWfn
 *
 */
namespace mpqc{
namespace scf{

class RHF : public qc::AOWavefunction<TA::TensorD,TA::SparsePolicy> {

public:
  using array_type = TA::TSpArrayD;

  RHF() = default;

  /**
   * KeyVal constructor for RHF
   * keywords
   * @param
   *
   */

  RHF(const KeyVal& kv);

  double value() override;
  void obsolete() override;

  double energy() const;

  inline array_type const &overlap() const { return S_; }
  inline array_type const &fock() const { return F_; }
  inline array_type const &density() const { return D_; }
  inline array_type const &coefficents() const { return C_; }
  inline double rhf_energy() const { return this->energy_; }

  /*! Function to compute the density to the desired accuracy.
   *
   * Takes some form of integral and does the rhf iterations.  The place to
   *specialized is in build_fock.
   *
   * returns true if the calculation converged to the desired threshold in
   *fewer than max_iters
   */
  bool solve(int64_t max_iters, double thresh);

  virtual rapidjson::Value results(rapidjson::Document &) const;


protected:
  double converge_;
  std::size_t max_iter_;
  double repulsion_;

  array_type H_;
  array_type S_;
  array_type F_;
  array_type F_diis_;
  array_type D_;
  array_type C_;

  std::unique_ptr<FockBuilder> f_builder_;
  std::unique_ptr<DensityBuilder> d_builder_;

  std::vector<double> rhf_times_;
  std::vector<double> d_times_;
  std::vector<double> build_times_;


private:
  void init(const KeyVal& kv);
  virtual void init_fock_builder();
  void compute_density();
  void build_F();

private:

  const KeyVal kv_;

};


/**
 *
 * RIRHF class, fock_builder is overide to use three center integral
 */

class RIRHF: public RHF{

public:
  RIRHF(const KeyVal& kv);

private:
  void init_fock_builder() override;

};


/**
 * DirectRIRHF, fock_builder is overide to use direct three center integral
 */

class DirectRIRHF : public RHF{

public:
  DirectRIRHF(const KeyVal& kv);

private:
  void init_fock_builder() override;
};

/**
 * DirectRHF, fock_builder is overide to use direct four center integral
 */
class DirectRHF : public RHF{

public:
  DirectRHF(const KeyVal& kv);

private:
  void init_fock_builder() override;
};

} // end of namespace scf
} // end of namespace mpqc

#endif //MPQC_CHEMISTRY_QC_SCF_RHF_H
