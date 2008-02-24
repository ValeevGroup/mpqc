//
// mbpt.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ibniels@kemi.aau.dk>
// Maintainer: LPS
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

#ifndef _chemistry_qc_mbpt_mbpt_h
#define _chemistry_qc_mbpt_mbpt_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/** The MBPT2 class implements several second-order perturbation theory
methods. */
class MBPT2: public Wavefunction {
  protected:
#define ref_to_mp2_acc 100.0

    Ref<SCF> reference_;
    Ref<MemoryGrp> mem;
    int nfzc, nfzv;
    size_t mem_alloc;

    double cphf_epsilon_;
    int eliminate_in_gmat_;
    const double *intbuf_;
    Ref<TwoBodyInt> tbint_;
    Ref<TwoBodyInt> *tbints_;
    Ref<TwoBodyDerivInt> *tbintder_;
    int nbasis;
    int noso;
    Ref<MessageGrp> msg_;
    int nvir, nocc, nsocc;

    Ref<ThreadGrp> thr_;

    // use a dynamic load balance algorithm if possible if true
    // (will not work if messagegrp not thread safe and
    // memorygrp needs catchup to work)
    int dynamic_;

    // control how frequently progress is printed
    double print_percent_;

    // The maximum number of orbitals in a pass.
    int max_norb_;

    // the irreps of the orbitals and the offset within the irrep
    int *symorb_irrep_;
    int *symorb_num_;

    std::string method_;
    std::string algorithm_;
    // if do_d1_ is true, D1(MP2) will be computed even if the gradient is not
    int do_d1_;
    // if do_d2_ is true, D2(MP1) will be computed 
    int do_d2_;
    
    int nfuncmax;

    double hf_energy_;
    RefSCVector hf_gradient_;

    double restart_ecorr_;
    int restart_orbital_v1_;
    int restart_orbital_memgrp_;

  protected:
    void init_variables();

    // implement the Compute::compute() function
    void compute();

    // Fill in the eigenvectors and eigenvalues (Guest & Saunders general
    // form is used for the Fock matrix in the open shell case).
    void eigen(RefDiagSCMatrix &vals, RefSCMatrix &vecs,
               RefDiagSCMatrix &occs);

    // calculate the opt2 energy using algorithm v1
    void compute_hsos_v1();

    // calculate the opt2 energy using algorithm v2
    distsize_t compute_v2_memory(int ni,
                             int nfuncmax, int nbfme, int nshell,
                             int ndocc, int nsocc, int nvir, int nproc);
    void compute_hsos_v2();

    // calculate the opt2 energy using the load balanced version of v2
    void compute_hsos_v2_lb();

    // calculate the closed shell mp2 energy and gradient
    int compute_cs_batchsize(size_t mem_static, int nocc_act);
    // distsize_t is used to allow memory requirements to be
    // estimated by starting the calculation on a single processor
    distsize_t compute_cs_dynamic_memory(int ni, int nocc_act);
    int make_cs_gmat(RefSymmSCMatrix& Gmat, double *DPmat);
    int make_cs_gmat_new(RefSymmSCMatrix& Gmat, const RefSymmSCMatrix& DPmat);
    void form_max_dens(double *DPmat, signed char *maxp);
    int init_cs_gmat();
    void done_cs_gmat();
    int make_g_d_nor(RefSymmSCMatrix& Gmat,
                     double *DPmat, const double *mgdbuff);
    void cs_cphf(double **scf_vector,
                 double *Laj, double *eigval, RefSCMatrix& P2aj);
    void s2pdm_contrib(const double *intderbuf, double *PHF,
                       double *P2AO, double **hf_ginter, double **ginter);
    void hcore_cs_grad(double *PHF, double *PMP2,
                       double **hf_ginter, double **ginter);
    void overlap_cs_grad(double *WHF, double *WMP2,
                         double **hf_ginter, double **ginter);
    void compute_cs_grad();
  public:
    MBPT2(StateIn&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>reference</tt><dd> This gives the reference wavefunction.
        It must be an object of type CLSCF for closed-shell molecules and
        HSOSSCF for open-shell molecules.  The is no default.

        <dt><tt>nfzc</tt><dd> The number of frozen core orbitals.  The
        default is 0.  If no atoms have an atomic number greater than 30,
        then the number of orbitals to be frozen can be automatically
        determined by specifying nfzc = auto.

        <dt><tt>nfzv</tt><dd> The number of frozen virtual orbitals.  The
        default is 0.

        <dt><tt>memory</tt><dd> The amount of memory, in bytes, that each
        processor may use.

        <dt><tt>method</tt><dd> This gives a string that must take on one
        of the values below.  The default is mp for closed-shell systems
        and zapt for open-shell systems.

        <dl>

          <dt><tt>mp</tt><dd> Use M&oslash;ller-Plesset perturbation theory.
          This is only valid for closed-shell systems.  Energies and
          gradients can be computed with this method.

          <dt><tt>opt1</tt><dd> Use the OPT1 variant of open-shell
          perturbation theory.  Only energies can be computed for
          open-shell systems.

          <dt><tt>opt2</tt><dd> Use the OPT2 variant of open-shell
          perturbation theory.  Only energies can be computed for
          open-shell systems.

          <dt><tt>zapt</tt><dd> Use the ZAPT variant of open-shell
          perturbation theory.  Only energies can be computed for
          open-shell systems.

        </dl>

        <dt><tt>algorithm</tt><dd> This gives a string that must take on
        one of the values given below.  The default is memgrp for
        closed-shell systems.  For open-shell systems v1 is used for a
        small number of processors and v2 is used otherwise.

        <dl>

          <dt><tt>memgrp</tt><dd> Use the distributed shared memory
          algorithm (which uses a MemoryGrp object).  This is only valid
          for MP2 energies and gradients.

          <dt><tt>v1</tt><dd> Use algorithm V1.  Only energies can be
          computed.  The maximum number of processors that can be utilized
          is the number of virtual orbitals.  This algorithm computes few
          integrals than the others, but has higher communication
          requirements.

          <dt><tt>v2</tt><dd> Use algorithm V2.  Only energies can be
          computed.  The maximum number of processors that can be utilized
          is the number of shells.

          <dt><tt>v2lb</tt><dd> Use a modified V2 algorithm that may
          compute more two electron integrals, but may get better load
          balance on the \f$O(n_\mathrm{basis}^5)\f$ part of the
          calculation.  Only energies can be computed.  This is recommended
          only for computations involving large molecules (where the
          transformation is dominant) on very many processors (approaching
          the number of shells).

        </dl>

        The v1 and v2 algorithms are discussed in Ida M. B. Nielsen and
        Edward T. Seidl, J. Comp. Chem. 16, 1301 (1995).  The memgrp
        algorithm is discussed in Ida M. B. Nielsen, Chem. Phys. Lett. 255,
        210 (1996).

        <dt><tt>memorygrp</tt><dd> A MemoryGrp object is used by the memgrp
        algorithm.  If this is not given the program will try to find an
        appropriate default.

        </dl> */
    MBPT2(const Ref<KeyVal>&);
    ~MBPT2();

    void save_data_state(StateOut&);

    Ref<SCF> ref() { return reference_; }
    double ref_energy();
    double corr_energy();
    RefSCVector ref_energy_gradient();
    RefSCVector corr_energy_gradient();

    int nelectron();

    int nfzcore() const { return nfzc; };
    int nfzvirt() const { return nfzv; };

    RefSymmSCMatrix density();
    int spin_polarized();

    int gradient_implemented() const;
    int value_implemented() const;
    /// set the value accuracy
    void set_desired_value_accuracy(double);
    void symmetry_changed();

    // override compute's obsolete so we can call the reference's obsolete
    void obsolete();

    void print(std::ostream&o=ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
