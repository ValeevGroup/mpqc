//
// util.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#include <chemistry/qc/mbpt/util.h>

using namespace sc;

BiggestContribs::BiggestContribs(int nindex, int maxcontrib)
{
  maxcontrib_  = maxcontrib;
  nindex_ = nindex;

  vals_ = new double[maxcontrib_];

  indices_ = new int*[maxcontrib_];
  for (int i=0; i<maxcontrib_; i++) indices_[i] = new int[nindex_];

  ncontrib_ = 0;
}

BiggestContribs::~BiggestContribs()
{
  delete[] vals_;
  for (int i=0; i<maxcontrib_; i++) delete[] indices_[i];
  delete[] indices_;
}

void
BiggestContribs::insert(double val, int i0, int i1)
{
  int i[2];
  i[0] = i0;
  i[1] = i1;
  insert(val,i);
}

void
BiggestContribs::insert(double val, int i0, int i1, int i2, int i3)
{
  int i[4];
  i[0] = i0;
  i[1] = i1;
  i[2] = i2;
  i[3] = i3;
  insert(val,i);
}

void
BiggestContribs::insert(double val, int i0, int i1, int i2, int i3, int i4)
{
  int i[5];
  i[0] = i0;
  i[1] = i1;
  i[2] = i2;
  i[3] = i3;
  i[4] = i4;
  insert(val,i);
}

void
BiggestContribs::insert(double val, const int *ii)
{
  if (maxcontrib_ == 0) return;
  if (ncontrib_ == 0) {
      vals_[0] = val;
      memcpy(indices_[0], ii, nindex_*sizeof(int));
      ncontrib_ = 1;
      return;
    }
  int i;
  double fabsval = (val>=0.0?val:-val);
  int contrib = 0;
  for (i=ncontrib_-1; i>=0; i--) {
      double fabstmp = (vals_[i]>=0.0?vals_[i]:-vals_[i]);
      if (fabsval > fabstmp) {
          contrib = 1;
          if (i < maxcontrib_-1) {
              vals_[i+1] = vals_[i];
              memcpy(indices_[i+1], indices_[i], nindex_*sizeof(int));
            }
          vals_[i] = val;
          memcpy(indices_[i], ii, nindex_*sizeof(int));
        }
      else {
          break;
        }
    }
  if (ncontrib_ < maxcontrib_) {
      if (!contrib) {
          vals_[ncontrib_] = val;
          memcpy(indices_[ncontrib_], ii, nindex_*sizeof(int));
        }
      ncontrib_++;
    }
}

void
BiggestContribs::combine(const Ref<MessageGrp> &grp)
{
  int i;
  int n = grp->n();
  int me = grp->me();
  // allocate array to store the number of contribs from each node
  int *ncontrib_each = new int[n];
  memset(ncontrib_each, 0, n*sizeof(int));
  ncontrib_each[me] = ncontrib_;
  grp->sum(ncontrib_each,n);

  // compute the total number of contribs and my offset for contrib data
  int ncontrib_all = 0;
  int contrib_offset = 0;
  for (i=0; i<n; i++) {
      ncontrib_all += ncontrib_each[i];
      if (i<me) contrib_offset += ncontrib_each[i];
    }

  // allocate storage for contrib data from all nodes
  double *vals_all = new double[ncontrib_all];
  int **indices_all = new int*[ncontrib_all];
  indices_all[0] = new int[ncontrib_all*nindex_];
  for (i=1; i<ncontrib_all; i++) indices_all[i] = &indices_all[0][nindex_*i];
  memset(vals_all, 0, sizeof(double)*ncontrib_all);
  memset(indices_all[0], 0, sizeof(int)*ncontrib_all*nindex_);

  // send contrib data from all nodes to 0
  for (i=0; i<ncontrib_; i++) {
      vals_all[contrib_offset+i] = vals_[i];
      memcpy(indices_all[contrib_offset+i], indices_[i], sizeof(int)*nindex_);
    }
  grp->sum(vals_all, ncontrib_all);
  grp->sum(indices_all[0], ncontrib_all*nindex_);

  ncontrib_ = 0;
  for (i=0; i<ncontrib_all; i++) {
      insert(vals_all[i], indices_all[i]);
    }

  delete[] ncontrib_each;
  delete[] vals_all;
  delete[] indices_all[0];
  delete[] indices_all;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
