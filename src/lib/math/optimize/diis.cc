//
// diis.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#include <math.h>

#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/cmatrix.h>
#include <math/optimize/diis.h>

using namespace sc;

static ClassDesc DIIS_cd(
  typeid(DIIS),"DIIS",4,"public SelfConsistentExtrapolation",
  0, create<DIIS>, create<DIIS>);

void
DIIS::init()
{
  int dim = ndiis+1;

  iter = 0;

  btemp = new double[dim];

  bmat = new double*[dim];
  bold = new double*[ndiis];

  if (!btemp || !bmat || !bold) {
    ExEnv::err0() << indent
         << "DIIS::init: alloc of bmat, bold, and btemp failed\n";
    abort();
  }

  for (int i=0; i < dim; i++) {
    bmat[i] = new double[dim];
    if (i < dim-1)
      bold[i] = new double[ndiis];
  }

  diism_data = new Ref<SCExtrapData>[ndiis];
  diism_error = new Ref<SCExtrapError>[ndiis];
  // diism_datain is on bigger than the other because
  // it must hold data associated with the next iteration
  diism_datain = new Ref<SCExtrapData>[ndiis+1];
}

DIIS::DIIS(int strt, int ndi, double dmp, int ngr, int ngrdiis, double mf) :
  btemp(0), bold(0), bmat(0), diism_data(0), diism_error(0), diism_datain(0),
  mixing_fraction(mf)
{
  start = strt;
  ndiis = ndi;
  damping_factor = dmp;
  mixing_fraction = mf;

  ngroup = ngr;
  ngroupdiis = ngrdiis;

  init();
}

DIIS::DIIS(StateIn& s) :
  SelfConsistentExtrapolation(s),
  btemp(0), bold(0), bmat(0), diism_data(0), diism_error(0), diism_datain(0)
{
  int i;

  s.get(start);
  s.get(ndiis);
  s.get(iter);
  if (s.version(::class_desc<DIIS>()) >= 2) {
    s.get(ngroup);
    s.get(ngroupdiis);
  }
  s.get(damping_factor);
  if (s.version(::class_desc<DIIS>()) >= 4) {
    s.get(mixing_fraction);
  }

  // alloc storage for arrays
  btemp = new double[ndiis+1];

  bold = new double*[ndiis];
  for (i=0; i < ndiis; i++)
    bold[i] = new double[ndiis];

  bmat = new double*[ndiis+1];
  for (i=0; i <= ndiis; i++)
    bmat[i] = new double[ndiis+1];

  diism_data = new Ref<SCExtrapData>[ndiis];
  diism_error = new Ref<SCExtrapError>[ndiis];
  // see DIIS::init() for why ndiis+1
  diism_datain = new Ref<SCExtrapData>[ndiis+1];

  // read arrays
  int ndat = iter;
  if (iter > ndiis) ndat = ndiis;
  if (s.version(::class_desc<DIIS>()) < 3) ndat = ndiis;

  s.get_array_double(btemp,ndat+1);

  for (i=0; i < ndat; i++)
    s.get_array_double(bold[i],ndat);

  for (i=0; i <= ndat; i++)
    s.get_array_double(bmat[i],ndat+1);

  for (i=0; i < ndat; i++) {
    diism_data[i] << SavableState::restore_state(s);
    diism_error[i] << SavableState::restore_state(s);
  }
  if (s.version(::class_desc<DIIS>()) >= 4) {
    for (i=0; i <= ndat; i++) {
      diism_datain[i] << SavableState::restore_state(s);
    }
  }
}

DIIS::DIIS(const Ref<KeyVal>& keyval):
  SelfConsistentExtrapolation(keyval),
  btemp(0), bold(0), bmat(0), diism_data(0), diism_error(0), diism_datain(0)
{
  ndiis = keyval->intvalue("n");
  if (keyval->error() != KeyVal::OK) ndiis = 5;

  start = keyval->intvalue("start");
  if (keyval->error() != KeyVal::OK) start = 1;

  ngroup = keyval->intvalue("ngroup");
  if (keyval->error() != KeyVal::OK) ngroup = 1;

  ngroupdiis = keyval->intvalue("ngroupdiis");
  if (keyval->error() != KeyVal::OK) ngroupdiis = 1;

  damping_factor = keyval->doublevalue("damping_factor");
  if (keyval->error() != KeyVal::OK) damping_factor = 0;

  mixing_fraction = keyval->doublevalue("mixing_fraction");
  if (keyval->error() != KeyVal::OK) mixing_fraction = 0;

  if (ndiis <= 0) {
    ExEnv::err0() << indent
         << "DIIS::DIIS(const Ref<KeyVal>& keyval): got ndiis = 0\n";
    abort();
  }

  init();
}

DIIS::~DIIS()
{
  if (btemp) {
    delete[] btemp;
    btemp=0;
  }

  if (bold) {
    for (int i=0; i < ndiis; i++) {
      if (bold[i])
        delete[] bold[i];
    }
    delete[] bold;
    bold=0;
  }

  if (bmat) {
    for (int i=0; i <= ndiis; i++) {
      if (bmat[i])
        delete[] bmat[i];
    }
    delete[] bmat;
    bmat=0;
  }

  if (diism_data) {
    delete[] diism_data;
    diism_data=0;
  }

  if (diism_datain) {
    delete[] diism_datain;
    diism_datain=0;
  }

  if (diism_error) {
    delete[] diism_error;
    diism_error=0;
  }
}

void
DIIS::save_data_state(StateOut& s)
{
  int i;

  SelfConsistentExtrapolation::save_data_state(s);
  s.put(start);
  s.put(ndiis);
  s.put(iter);
  s.put(ngroup);
  s.put(ngroupdiis);
  s.put(damping_factor);
  s.put(mixing_fraction);

  int ndat = iter;
  if (iter > ndiis) ndat = ndiis;

  s.put_array_double(btemp, ndat+1);

  for (i=0; i < ndat; i++)
    s.put_array_double(bold[i], ndat);

  for (i=0; i <= ndat; i++)
    s.put_array_double(bmat[i], ndat+1);

  for (i=0; i < ndat; i++) {
    SavableState::save_state(diism_data[i].pointer(),s);
    SavableState::save_state(diism_error[i].pointer(),s);
  }
  for (i=0; i <= ndat; i++) {
    SavableState::save_state(diism_datain[i].pointer(),s);
  }
}

void
DIIS::reinitialize(Ref<SCExtrapData> data)
{
  iter=0;
  if (data) {
    const bool do_mixing = (mixing_fraction != 0.0);
    if (do_mixing) diism_datain[0] = data->copy();
  }
  else {
    diism_datain[0] = 0;
  }
}

void
DIIS::start_extrapolation()
{
  if (start > iter) start = iter+1;
}

int
DIIS::extrapolate(const Ref<SCExtrapData>& data,
                  const Ref<SCExtrapError>& error)
{
  int i, j, k;
  int last = iter;
  int trial = 0;
  int col = iter + 2;
  double norm, determ;
  double scale;

  const bool do_mixing = (mixing_fraction != 0.0);

  iter++;

  scale = 1.0 + damping_factor;

  if (iter > ndiis) {
      last = ndiis-1;
      col = ndiis+1;
      for (i=0; i < last ; i++) {
          diism_data[i] = diism_data[i+1];
          diism_error[i] = diism_error[i+1];
        }
      for (i=0; i <= last ; i++) {
          diism_datain[i] = diism_datain[i+1];
      }
    }

  diism_data[last] = data->copy();
  diism_error[last] = error;

  set_error(error->error());

  // then set up B matrix, where B(i,j) = <ei|ej>

  // move bold(i+1,j+1) to bold(i,j)
  if (iter > ndiis) {
      for (i=0; i < last ; i++) {
          for (j=0; j <= i ; j++) {
              bold[i][j]=bold[j][i]=bold[i+1][j+1];
            }
        }
    }

  // and set the current rows of bold
  for (i=0; i <= last ; i++)
      bold[i][last]=bold[last][i] =
                      diism_error[i]->scalar_product(diism_error[last]);

  bmat[0][0] = 0.0;
  btemp[0] = -1.0;

  if (bold[0][0] > 1.e-10) {
      norm = 1.0/bold[0][0];
    }
  else {
      norm = 1.0;
    }

  for (i=1; i <= last+1 ; i++) {
      bmat[i][0]=bmat[0][i] = -1.0;
      btemp[i] = 0.0;
      for (j=1; j <= i ; j++) {
          bmat[i][j]=bmat[j][i] = bold[i-1][j-1]*norm;
          if (i==j) bmat[i][j] *= scale;
        }
    }

  // finally, solve the set of linear equations, obtain the coefficients,
  // and form the new fock matrix F= sum(i=1,n) ci*Fi

  if (iter-1) {
      determ = cmat_solve_lin(bmat,0,btemp,col);

      // test for poorly conditioned equations */
      while (fabs(determ) < 1.0e-19 && trial < last) {

          trial++;
          col--;

          bmat[0][0] = 0.0;
          btemp[0] = -1.0;

          if (bold[trial][trial] > 1.e-10) {
              norm=1.0/bold[trial][trial];
            }
          else {
              norm = 1.0;
            }

          for (i=1; i <= last-trial+1 ; i++) {
              bmat[i][0]=bmat[0][i] = -1.0;
              for (j=1; j <= i ; j++) {
                  bmat[i][j]=bmat[j][i]=bold[i+trial-1][j+trial-1]*norm;
                  if (i==j) bmat[i][j] *= scale;
                }
              btemp[i] = 0.0;
            }

          determ = cmat_solve_lin(bmat,0,btemp,col);
        }

      if (fabs(determ) < 10.0e-20) {
        ExEnv::err0() << indent
             << "DIIS::extrapolate:  trial " << trial << " no good\n";
        return -1;
      }

  //////////////////////////////////////////////////////////////
  // In the following, F_{in} is the F matrix we produce as the
  // input to the next iteration. F_{out} is the F matrix we
  // were given.

  // Charge mixing:
  // F_{in,i+1} = f F_{in,i} + (1-f) F_{out,i}
  // P_{i+1} <- F_{in,i+1}
  // F_{out,i+1} <- P_{i+1}

  // DIIS:
  // for (k=0; k<=i; k++) F_{in,i+1} = C_{k,i} F_{out,k}

  // DIIS + charge mixing:
  // for (k=0; k<=i; k++)
  //   F_{in,i+1} += C_{k,i} ( (1-f) F_{out,k} + f F_{in,k} )
  //
  // The problem with the above is that in the case of HF, we start with an
  // initial density and don't necessarily have anything like an initial
  // Fock matrix that generates that density. So we do not do density
  // mixing on the initial item:
  // F_{in,1} = C_{k,1} F_{out,0}
  // for (k=0; k<=i; k++)
  //   F_{in,i+1} += C_{k,i} ( (1-f) F_{out,k} + f F_{in,k} )

      if (iter >= start && (((iter-start)%ngroup) < ngroupdiis)) {
          int kk=1;

          data->zero();

          // This is the code before mixing_fraction was implemented.
//           for (k=trial; k < last+1 ; k++) {
//               data->accumulate_scaled(btemp[kk], diism_data[k]);
//               kk++;
//             }

          for (k=trial; k < last+1; k++) {
              if (diism_datain[k].null()) {
                data->accumulate_scaled(btemp[kk], diism_data[k]);
              }
              else {
                data->accumulate_scaled(btemp[kk]*(1.0-mixing_fraction),
                                        diism_data[k]);
                data->accumulate_scaled(btemp[kk]*mixing_fraction,
                                        diism_datain[k]);
              }
              kk++;
            }
        }
    }
  else {
      //first iteration
      data->zero();
      if (diism_datain[0].null()) {
          data->accumulate_scaled(1.0, diism_data[0]);
      }
      else {
          data->accumulate_scaled((1.0-mixing_fraction),
                                        diism_data[0]);
          data->accumulate_scaled(mixing_fraction,
                                  diism_datain[0]);
      }
  }

  // only need to keep data in diis_datain if doing mixing
  if (do_mixing) diism_datain[last+1] = data->copy();

  return 0;
}

void
DIIS::print(std::ostream& o) const
{
  o << indent
    << "DIIS: "
    << "n=" << ndiis
    << ", start=" << start
    << ", ngroup=" << ngroup
    << ", ngroupdiis=" << ngroupdiis
    << ", damping_factor=" << damping_factor
    << ", mixing_fraction=" << mixing_fraction
    << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
