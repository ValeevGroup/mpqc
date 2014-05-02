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
#include <vector>
#include <functional>
#include <util/misc/scexception.h>

namespace mpqc {
  namespace TA {
    /// @addtogroup ChemistryBasisIntegralTA
    /// @{

    /**
     * IntegralEnginePool is a class that will take a Ref<Engine> as a prototype
     * and then clone it multiple times such that each thread has its own
     * integral engine.  This is necessary because the buffers in the
     * integral engines don't have thread-safe access.
     */
    template<typename RefEngType>
    class IntegralEnginePool {
      /// The thread local storage key
      pthread_key_t engine_key_;
      RefEngType prototype_;
      // Pointers to the TLS pointers which will allows us to delete them in the constructor.
      std::vector<void*> tls_ptrs;

    public:

      typedef RefEngType engine_type;

      /**
       * IntegralEnginePool constructor it takes a sc::Ref<sc::IntegralEngine>
       * and sets that as the prototype for all the thread local engines
       */
      IntegralEnginePool(const RefEngType engine) :
              prototype_(engine) {
        if (pthread_key_create(&engine_key_, nullptr) != 0) // objects will be manually destroyed by the destructor
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
        /*
         * The following loop manually calls delete for each TLS RefEngine pointer
         * since they normally will not be deallocated until madness::Finalize()
         * is called effectively leading to a memory leak
         */
        for(void *ptr : tls_ptrs){ // For all TLS pointers
          std::cout << "deleting pointer " << ptr << std::endl;
          destroy_thread_object(ptr); // destory it
        }

        //// if main thread was used, call destructor for its TLS now
        //void* tls_ptr = pthread_getspecific(key_);
        //if (tls_ptr != 0) {
        //  destroy_thread_object(tls_ptr);
        //}

        if (pthread_key_delete(engine_key_) != 0) {
          std::cout << "WARNING: pthread_key_delete failed" << std::endl;
        }
      }

      /**
       * Function that returns a clone of the prototype engine.  The clone
       * exist for only one thread to use and so cannot be written to or used
       * by any other thread.
       */
      RefEngType instance() {

        // Declare a pointer to a sc::Ref<Engine>
        RefEngType *RefEngine =
                reinterpret_cast<RefEngType*>(pthread_getspecific(engine_key_));


        if (RefEngine == nullptr) {
          RefEngine = new RefEngType;

          // Get clone of prototype, must lock to ensure nobody else
          // clones at the same time. As long as number of threads isn't
          // Super high this should be ok.   Although in the future
          // clone might be thread safe.
          mutex::global::lock(); // <<< Begin Critical Section
          *RefEngine = prototype_->clone();
          mutex::global::unlock(); // <<< End Critical Section


          // Asign RefEngine to its thread specific partner.
          pthread_setspecific(engine_key_, RefEngine);

          tls_ptrs.push_back(pthread_getspecific(engine_key_));

        }

        return *RefEngine;

      }

    private:
      /// Function to destroy the thread objects
      static void destroy_thread_object(void* p) {
        RefEngType* ptr = reinterpret_cast<RefEngType*>(p);
        if(ptr){ // Make sure we havn't already deleted the data
          delete ptr;
          ptr = nullptr;
        }
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
