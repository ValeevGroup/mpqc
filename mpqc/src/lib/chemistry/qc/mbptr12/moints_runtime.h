//
// moints_runtime.h
//
// Copyright (C) 2008 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_mointsruntime_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_mointsruntime_h

#include <util/state/state.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/transform_factory.h>

namespace sc {

  /**
   *    Smart runtime support for computing MO-basis integrals.
   */
  class MOIntsRuntime : virtual public SavableState {
    public:
      // for now use the typedef, later replace R12IntsAcc with TwoBodyIntsAcc
      typedef R12IntsAcc TwoBodyIntsAcc;
      MOIntsRuntime(const Ref<MOIntsTransformFactory>& factory);
      MOIntsRuntime(StateIn& si);
      void save_data_state(StateOut& so);

      /** Returns the TwoBodyIntsAcc that contains the integrals described by key.
          The desired integrals may be a subset of the given TwoBodyIntsAcc.

          key must be in format recognized by ParsedTwoBodyIntKey, which matched the format used by key() or key_mulliken().
          If this key is not known, the integrals will be computed by an appropriate TwoBodyMOIntsTranform object.

          \sa key(), \sa key_mulliken()
        */
      Ref<TwoBodyIntsAcc> get(const std::string& key);   // non-const: can add transforms

      /** Returns key that describes params. If not registered, will do so using register_params(params)
       */
      std::string params_key(const Ref<IntParams>& params) const;
      /** Returns params that correspond to key. Returns null ptr if key is not known.
        */
      Ref<IntParams> params(const std::string& key) const;
      // register key->params mapping
      void register_params(const std::string& key, const Ref<IntParams>& params) const;
      // register key->params mapping. A unique random key is generated automatically and returned
      std::string register_params(const Ref<IntParams>& params) const;
      /** Returns key that corresponds to descr.
          \sa params_key
        */
      std::string descr_key(const Ref<TwoBodyIntDescr>& descr);

      /// returns the factory
      const Ref<MOIntsTransformFactory>& factory() const { return factory_; }

      /// describes the physical layout of the integrals in TwoBodyIntsAcc
      class Layout {
        public:
          Layout(const std::string& str);
          Layout(const Layout& other);
          Layout& operator=(const Layout& other);
          bool operator==(const Layout& other) const;

        private:
          typedef enum {
            b1b2_k1k2, // physicists layout
            b1k1_b2k2  // chemists layout
          } Type;
          Type type_;
      };
      static Layout Layout_b1b2_k1k2;  // physicists layout
      static Layout Layout_b1k1_b2k2;  // chemists layout

    private:
      Ref<MOIntsTransformFactory> factory_;  // that creates transforms

      // Maps key to IntParams
      typedef std::map< std::string, Ref<IntParams> > ParamsMap;
      mutable ParamsMap params_map_;  // key->params map can be modified silently, no need to observe const-ness

      // Map to TwoBodyMOIntsTransform objects that have been computed previously
      typedef std::map<std::string, Ref<TwoBodyMOIntsTransform> > TformMap;
      TformMap tform_map_;

      // creates a transform object for a given key
      const Ref<TwoBodyMOIntsTransform>& create_tform(const std::string& key);

      /** creates TwoBodyIntDescr given the operator and params keys.
          params_key must be in params_map_.
       */
      Ref<TwoBodyIntDescr> create_descr(const std::string& oper_key,
                                        const std::string& params_key);

  };

  /// Parsed representation of a string key that represents a set of 2-body integrals
  class ParsedTwoBodyIntKey {
    public:
      ParsedTwoBodyIntKey(const std::string& key);

      const std::string& key() const { return key_; }
      const std::string& bra1() const { return bra1_; }
      const std::string& bra2() const { return bra2_; }
      const std::string& ket1() const { return ket1_; }
      const std::string& ket2() const { return ket2_; }
      const std::string& oper() const { return oper_; }
      const std::string& params() const { return params_; }
      const std::string& layout() const { return layout_; }

      /// computes key from its components
      static std::string key(const std::string& bra1,
                             const std::string& bra2,
                             const std::string& ket1,
                             const std::string& ket2,
                             const std::string& oper,
                             const std::string& params,
                             const std::string& layout);
      /// computes operator part of the key given an TwoBodyIntDescr object
      static std::string key(const Ref<TwoBodyIntDescr>& descr);
      /// this factory method constructs a descriptor given operator key + IntParams object + Integrals object
      static Ref<TwoBodyIntDescr> create_descr(const std::string& oper_key,
                                               const Ref<IntParams>& p,
                                               const Ref<Integral>& integral);

    private:
      std::string key_;
      std::string bra1_, bra2_, ket1_, ket2_;
      std::string oper_;
      std::string params_;
      std::string layout_;
  };

  namespace detail {
    /// Convert 2 spaces to SpinCase2 assuming that lower-case symbols correspond to beta, upper -- to alpha
    SpinCase2 spincase2(const Ref<MOIndexSpace>& space1,
                        const Ref<MOIndexSpace>& space2);
    std::string id(SpinCase2 S);
  }

} // end of namespace sc

#endif /* end of header guard */

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
