/*
 * wfn.h
 *
 *  Created on: Apr 27, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC_CHEMISTRY_QC_WFN_WFN_H_
#define MPQC_CHEMISTRY_QC_WFN_WFN_H_

#include <mpqc/chemistry/qc/wfn/wfn_world.h>
#include <mpqc/chemistry/qc/properties/propertybase.h>

#include <memory>
#include <functional>

namespace mpqc {
namespace qc {

class PropertyBase;

/**
 * Wfn is a wave function that lives in a WfnWorld
 */
class Wavefunction : public DescribedClass {
 private:
  /** Pointer to the WfnWorld
   *
   * \note No need to make this shared Wfn is just a member of the world it
   *lives in so no ownership here.
   *
   * \warning Wfn should never delete or allocate this pointer.
   *
   * \note by chong I changed this to shared pointer, for example, MP2 and HF
   *          will share the same wfn_world
   */
  std::shared_ptr<WavefunctionWorld> wfn_world_;

 protected:
  double energy_ = 0.0;

 public:
  /**
   *  \brief KeyVal constructor
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | wfn_world or $:wfn_world or ..:wfn_world | WavefunctionWorld | none | none |
   *
   *
   */
  Wavefunction(const KeyVal& kv);
  virtual ~Wavefunction();

  virtual void compute(PropertyBase* pb) = 0;
  virtual double value() = 0;
  virtual void obsolete() {
    energy_ = 0.0;
  };

  std::shared_ptr<WavefunctionWorld> wfn_world() { return wfn_world_; }
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_WFN_H_
