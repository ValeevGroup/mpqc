//
// g92.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifndef _chemistry_qc_g92_h
#define _chemistry_qc_g92_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>

#include <chemistry/qc/wfn/obwfn.h>
#include <math/scmat/matrix.h>

class Gaussian92: public OneBodyWavefunction
{
#   define CLASSNAME Gaussian92
#   include <util/state/stated.h>
#   include <util/class/classda.h>
protected:
  int charge_;
  int multiplicity_;
  int memory_;
  int use_ckpt_;
  char *name_;
  char *scr_dir_;
  char *g92_dir_;
  char *basis_;

  ResultRefSCMatrix _eigenvectors;

  void compute();

  virtual char * emethod() = 0;
  virtual char * gmethod() = 0;
  virtual char * hmethod() = 0;
  
  virtual int run_energy();
  virtual int parse_g92_energy();
  
  virtual int run_gradient();
  virtual int parse_g92_gradient();
  
  virtual int run_hessian();
  virtual int parse_g92_hessian();

  virtual int run_g92(const char *method);

public:
  Gaussian92(const RefKeyVal&);
  Gaussian92(StateIn&);
  virtual ~Gaussian92();
  
  void save_data_state(StateOut&);

  void print(ostream& =cout);
  
  RefSCMatrix eigenvectors();
  int do_eigenvectors(int);

  RefDiagSCMatrix eigenvalues();
  double occupation(int irrep, int vectornum);

  RefSCMatrix normal_modes();
  RefSCVector frequencies();
};

class Gaussian92SCF: public Gaussian92
{
#   define CLASSNAME Gaussian92SCF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
private:
  char * emethod();
  char * gmethod();
  char * hmethod();

public:
  Gaussian92SCF(const RefKeyVal&);
  Gaussian92SCF(StateIn&);
  virtual ~Gaussian92SCF();
  
  void save_data_state(StateOut&);

  double occupation(int f) { return 0; }

  int value_implemented() const { return 1; }
  int gradient_implemented() const { return 1; }
  int hessian_implemented() const { return 1; }
};

class Gaussian92UHF: public Gaussian92
{
#   define CLASSNAME Gaussian92UHF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
private:
  char * emethod();
  char * gmethod();
  char * hmethod();

public:
  Gaussian92UHF(const RefKeyVal&);
  Gaussian92UHF(StateIn&);
  virtual ~Gaussian92UHF();
  
  void save_data_state(StateOut&);

  double occupation(int f) { return 0; }

  int value_implemented() const { return 1; }
  int gradient_implemented() const { return 1; }
  int hessian_implemented() const { return 1; }
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
