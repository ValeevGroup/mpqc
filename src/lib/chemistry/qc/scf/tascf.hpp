//
// tascf.hpp
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

#ifndef _MPQC_CHEMISTRY_QC_SCF_TASCF_HPP_
#define _MPQC_CHEMISTRY_QC_SCF_TASCF_HPP_

#include <chemistry/qc/wfn/tawfn.hpp>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/tbint.h>
#include <elemental-lite.hpp>

namespace mpqc{
  namespace TA {

    class SCF : public Wavefunction {
    public:
      typedef Wavefunction::TAMatrix TAMatrix;
      typedef Wavefunction::TAVector TAVector;
      //typedef Wavefunction::TAMatrixExpr TAMatrixExpr;
      typedef elem::DistMatrix<double, elem::VR, elem::STAR> ElemVector;
      typedef std::pair<ElemVector, TAMatrix> ElemTAEigenSystem;

      /** SCF KeyValue constructor
       *
       * */
      SCF(const sc::Ref<sc::KeyVal> &kval);
      virtual ~SCF();

      /// @return the number of electrons in the system
      virtual size_t nelectron() const override;

      /// @return the number of occupied orbitals in the system
      virtual size_t occupation() const {return occupation_;}

      /** @return the MO eigenvectors as a TiledArray::Array<double,2>
       * takes an accuracy, which supporst computing at a user defined
       * accuracy, but does not alter the classes internal desired_accuracy_
       * */
      TAMatrix MO_eigenvectors(double);

      /** @return the MO eigenvectors as a TiledArray::Array<double,2>
       * takes an accuracy, which supporst computing at a user defined
       * accuracy, but does not alter the classes internal desired_accuracy_
       * */
      TAMatrix
      MO_eigenvectors(){
        return MO_eigenvectors(MO_eigensystem_.desired_accuracy());
      }

      /** @return the MO eigenvalues as a TiledArray::Array<double,1>
       * takes an accuracy, which supporst computing at a user defined
       * accuracy, but does not alter the classes internal desired_accuracy_
       * */
      ElemVector
      MO_eigenvalues(double);

      /** @return the MO eigenvalues as a TiledArray::Array<double,1>
       * takes an accuracy, which supporst computing at a user defined
       * accuracy, but does not alter the classes internal desired_accuracy_
       * */
      ElemVector
      MO_eigenvalues(){
        return MO_eigenvalues(MO_eigensystem_.desired_accuracy());
      }

      /** @return the MO eigenvalues and eigenvectors as
       * a std::pair<TA::Array<double,1>, TA::Array<double,2>>
       * */
      ElemTAEigenSystem
      MO_eigensystem(double);

      /** @return the MO eigenvalues and eigenvectors as
       * a std::pair<TA::Array<double,1>, TA::Array<double,2>>
       * */
      ElemTAEigenSystem
      MO_eigensystem(){
        return MO_eigensystem(MO_eigensystem_.desired_accuracy());
      }

      /// @return the AO fock matrix computed to the desired accuracy
      TAMatrix& ao_fock(double);

      /// @return the AO fock matrix computed to the default accuracy
      TAMatrix& ao_fock(){return ao_fock(ao_fock_.desired_accuracy());}

      /** Returns an expression to ao_fock matrix.
       * If it has not been initialized or computed the it will compute
       * the matrix to the internal desired_accuracy
       * */
      //TAMatrixExpr ao_fock_expr(std::string);

      /// @return the converged scf energy
      virtual double scf_energy() = 0;

      void print(std::ostream &os = sc::ExEnv::out0()) const;

    protected:

      typedef Wavefunction::AccResultMatrix AccResultMatrix;
      typedef Wavefunction::AccResultVector AccResultVector;
      typedef sc::AccResult<ElemTAEigenSystem> AccResultEigenSystem;

      virtual void compute_ao_fock(double) = 0;

      // returns & to the current state of ao_fock_ used for computing.
      virtual TAMatrix& scf_ao_fock_(){return ao_fock_.result_noupdate();}

      unsigned int miniter() const {return miniter_; }
      unsigned int maxiter() const {return maxiter_; }

      void set_occupation(unsigned int i){occupation_ = i;}

    private:
      // default number of iterations to use
      unsigned int maxiter_ = 100;
      unsigned int miniter_ = 0;

      // The ao_fock matrix
      AccResultMatrix ao_fock_;
      // Holds the eigensystem of the ao_fock matrix
      AccResultEigenSystem MO_eigensystem_;

      // Number of electrons
      size_t occupation_ = 0;

      static sc::ClassDesc class_desc_;

    };
  } // namespace mpqc::TA
} //namespace mpqc


#endif /* _MPQC_CHEMISTRY_QC_SCF_TASCF_HPP_ */
