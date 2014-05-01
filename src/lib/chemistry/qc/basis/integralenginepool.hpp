//
// integralenginepool.hpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef MPQC_CHEMISTRY_QC_BASIS_INTEGRALENGINEPOOL_HPP_
#define MPQC_CHEMISTRY_QC_BASIS_INTEGRALENGINEPOOL_HPP_

#include <pthread.h>
#include <mpqc/utility/mutex.hpp>
#include <chemistry/qc/basis/integral.h>
#include <util/misc/scexception.h>

namespace mpqc {
  namespace TA {
    /// @addtogroup ChemistryBasisIntegralTA
    /// @{

    /**
     * IntegralEnginePool is a class that will take a Ref<Engine> as a prototype
     * and then clone it multiple times such that each thread has its own
     * integral engine.  This is necesarry because the buffers in the
     * integral engines don't have thread safe access.
     */
    template<typename RefEngType>
    class IntegralEnginePool {
      /// The thread local storage key
      pthread_key_t key_;
      RefEngType prototype_;

    public:

      typedef RefEngType engine_type;

      /**
       * IntegralEnginePool constructor it takes a sc::Ref<sc::IntegralEngine>
       * and sets that as the prototype for all the thread local engines
       */
      IntegralEnginePool(const RefEngType engine) :
              prototype_(engine) {
        if (pthread_key_create(&key_, &destroy_thread_object) != 0)
          throw sc::SystemException(
                  "IntegralEnginePool::IntegralEnginePool() "
                  "Unable to register thread local storage key. "
                  "Likely due to a limited number of keys.",
                  __FILE__, __LINE__);
      }

      /**
       * Delete the key that corresponds to our thread local storage object
       */
      ~IntegralEnginePool() {
        std::cout << "Destructor Called" << std::endl;
        reinterpret_cast<RefEngType*>(pthread_getspecific(key_))->
            RefEngType::~RefEngType();
        delete reinterpret_cast<RefEngType*>(pthread_getspecific(key_));
        if(pthread_key_delete(key_) != 0){
          std::cout << "Warning pthread_key_delete failed" << std::endl;
        }
      }

      /**
       * Function that returns a clone of the prototype engine.  The clone
       * exist for only one thread to use and so cannot be written to or used
       * by any other thread.
       */
      RefEngType instance() const {

        // Declare a pointer to a sc::Ref<Engine>
        RefEngType *RefEngine =
                reinterpret_cast<RefEngType*>(pthread_getspecific(key_));

        if (RefEngine == NULL) {
          RefEngine = new RefEngType;

          // Get clone of prototype, must lock to ensure nobody else
          // clones at the same time. As long as number of threads isn't
          // Super high this should be ok.   Although in the future
          // clone might be thread safe.
          mutex::global::lock(); // <<< Begin Critical Section
          *RefEngine = prototype_->clone();
          mutex::global::unlock(); // <<< End Critical Section

          // Asign RefEngine to its thread specific partner.
          pthread_setspecific(key_, RefEngine);
        }

        return *RefEngine;

      }

    private:
      /// Function to destroy the thread objects
      static void destroy_thread_object(void* p) {
        std::cout << "I deleted the RefEngType" << std::endl;
        delete reinterpret_cast<RefEngType*>(p);
      }

      /**
       * Copy cosntruction and assignment are not allowed for
       * IntegralEnginePool.
       */
      IntegralEnginePool(const IntegralEnginePool &);
      IntegralEnginePool& operator=(const IntegralEnginePool &);

    };
// IntegralEnginePool

  }// namespace TA
} // namespace mpqc
/// @} // Integrals

#endif /* MPQC_CHEMISTRY_QC_BASIS_INTEGRALENGINEPOOL_HPP_ */
