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
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>

////////////////////////////////////////////////////////////////////////////

class MBPT2: public Wavefunction {
#   define CLASSNAME MBPT2
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
#define ref_to_mp2_acc 100.0

    int debug_;
    RefSCF reference_;
    RefMemoryGrp mem;
    int nfzc, nfzv;
    unsigned int mem_alloc;

    int eliminate_in_gmat_;
    const double *intbuf_;
    RefTwoBodyInt tbint_;
    RefTwoBodyDerivInt tbintder_;
    int nbasis;
    RefMessageGrp msg_;
    int nvir, nocc, nsocc;

    // the irreps of the orbitals and the offset within the irrep
    int *symorb_irrep_;
    int *symorb_num_;

    char *method_;
    char *algorithm_;
    // if dos2_ is true, s2 will be computed even if the gradient is not
    int dos2_;

    int nfuncmax;

    double hf_energy_;
    RefSCVector hf_gradient_;

    double restart_ecorr_;
    int restart_orbital_v1_;

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
                             int nfuncmax, int nbasis, int nbfme, int nshell,
                             int ndocc, int nsocc, int nvir, int nproc);
    void compute_hsos_v2();

    // calculate the opt2 energy using the load balanced version of v2
    void compute_hsos_v2_lb();

    // calculate the closed shell mp2 energy and gradient
    int compute_cs_batchsize(int mem_static, int nocc_act);
    // distsize_t is used to allow memory requirements to be
    // estimated by starting the calculation on a single processor
    distsize_t compute_cs_dynamic_memory(int ni, int nocc_act);
    int make_cs_gmat(RefSymmSCMatrix& Gmat, double *DPmat);
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
    MBPT2(const RefKeyVal&);
    ~MBPT2();

    void save_data_state(StateOut&);

    RefSCF ref() { return reference_; }
    double ref_energy();
    double corr_energy();
    RefSCVector ref_energy_gradient();
    RefSCVector corr_energy_gradient();

    int nelectron();

    RefSymmSCMatrix density();
    int spin_polarized();

    int gradient_implemented() const;
    int value_implemented() const;

    void symmetry_changed();

    // override compute's obsolete so we can call the reference's obsolete
    void obsolete();

    void print(ostream&o=cout);
};
SavableState_REF_dec(MBPT2);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
