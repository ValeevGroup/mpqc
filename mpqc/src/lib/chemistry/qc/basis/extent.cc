//
// extent.cc
//
// Copyright (C) 2000 Sandia National Labs
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <chemistry/qc/basis/extent.h>

ARRAY_def(ExtentData);

ShellExtent::ShellExtent()
{
  contributing_shells_ = 0;
  n_[0] = n_[1] = n_[2] = 0;
}

ShellExtent::~ShellExtent()
{
  delete[] contributing_shells_;
}

ArrayExtentData &
ShellExtent::data(int x, int y, int z)
{
  if (x>=n_[0] || y>=n_[1] || z>= n_[2] || x<0 || y<0 || z<0) {
      ExEnv::out() << "ShellExtent::data: out of bounds" << endl;
      abort();
    }
  return contributing_shells_[z + n_[2]*(y + n_[1]*x)];
}

ArrayExtentData &
ShellExtent::data(int *b)
{
  return data(b[0],b[1],b[2]);
}

double
ShellExtent::distance(double loc, int axis, int origin, int point)
{
  if (point < origin)
      return loc - (lower_[axis]+resolution_*(point+1));
  else if (point == origin)
      return 0.0;
  else
      return loc - (lower_[axis]+resolution_*point);
}

void
ShellExtent::init(const RefGaussianBasisSet&gbs,
                  double resolution, double tolerance)
{
  int i,j,k,l,m,n;

  RefMolecule mol = gbs->molecule();
  resolution_ = resolution;

  delete[] contributing_shells_;
  n_[0] = n_[1] = n_[2] = 0;
  for (i=0; i<3; i++) lower_[i] = 0.0;
  contributing_shells_ = 0;
  if (mol->natom() == 0) return;

  double upper[3];

  for (i=0; i<3; i++) upper[i] = lower_[i] = mol->r(0,i);

  for (i=0; i<mol->natom(); i++) {
      for (j=0; j<gbs->nshell_on_center(i); j++) {
          double r = gbs->shell(gbs->shell_on_center(i, j)).extent(tolerance);
          for (k=0; k<3; k++) {
              if (lower_[k]>mol->r(i,k)-r) lower_[k] = mol->r(i,k)-r;
              if (upper [k]<mol->r(i,k)+r) upper [k] = mol->r(i,k)+r;
            }
        }
    }

  for (i=0; i<3; i++) {
      double l;
      l = upper[i]-lower_[i];
      n_[i] = int(l/resolution_);
      if (n_[i]*resolution_ + lower_[i] < upper[i]) n_[i]++;
    }

  contributing_shells_ = new ArrayExtentData[n_[0]*n_[1]*n_[2]];

  for (i=0; i<mol->natom(); i++) {
      //ExEnv::out() << indent << "working on atom " << i << endl;
      //ExEnv::out() << incindent;
      for (j=0; j<gbs->nshell_on_center(i); j++) {
          int ishell = gbs->shell_on_center(i,j);
          const GaussianShell &shell = gbs->shell(ishell);
          double r = shell.extent(tolerance);
          int ir = int(r/resolution_);
          if (ir*resolution_ < r) ir++;
          int atom_block[3];
          for (l=0; l<3; l++) {
              atom_block[l] = int((mol->r(i,l)-lower_[l])/resolution_);
            }
          int block[3];
          //ExEnv::out() << indent << "working on shell " << ishell << endl;
          //ExEnv::out() << incindent;
          for (k=atom_block[0]-ir; k<=atom_block[0]+ir; k++) {
              block[0] = k;
              for (l=atom_block[1]-ir; l<=atom_block[1]+ir; l++) {
                  block[1] = l;
                  for (m=atom_block[2]-ir; m<=atom_block[2]+ir; m++) {
                      block[2] = m;
                      //ExEnv::out() << indent
                      //     << "working on block " << block[0]
                      //     << " " << block[1]
                      //     << " " << block[2] << endl;
                      // find the distance to the atom from this block
                      double dist = 0.0;
                      for (n=0; n<3; n++) {
                          double r
                              = distance(mol->r(i,n),n,atom_block[n],block[n]);
                          dist += r*r;
                        }
                      dist = sqrt(dist);
                      double bound = shell.monobound(dist);
                      //ExEnv::out() << indent
                      //     << "dist = " << dist
                      //     << " bound = " << bound << endl;
                      if (bound < tolerance) continue;
                      data(block).push_back(ExtentData(ishell, bound));
                    }
                }
            }
          //ExEnv::out() << decindent;
        }
      //ExEnv::out() << decindent;
    }
}

const ArrayExtentData &
ShellExtent::contributing_shells(double x, double y, double z)
{
  int i, block[3];
  block[0] = int((x-lower_[0])/resolution_);
  block[1] = int((y-lower_[1])/resolution_);
  block[2] = int((z-lower_[2])/resolution_);
  for (i=0; i<3; i++) {
      if (block[i] < 0 || block[i] >= n_[i]) return null_;
    }
  return data(block);
}

void
ShellExtent::print(ostream &o)
{
  int i,j,k,l;

  o << node0
    << indent << "ShellExent:" << endl;

  o << node0 << incindent;

  o << node0
    << indent << "n = " << n_[0] << " " << n_[1] << " " << n_[2] << endl;;

  o << node0 << indent << "resolution = " << resolution_ << endl;

  o << node0 << indent
    << "lower = " << lower_[0]
    << " " << lower_[1]
    << " " << lower_[2] << endl;;

  o << node0 << indent
    << "upper = " << lower_[0] + n_[0] * resolution_
    << " " << lower_[1] + n_[1] * resolution_
    << " " << lower_[2] + n_[2] * resolution_ << endl;;

  for (i=0; i<n_[0]; i++) {
      for (j=0; j<n_[1]; j++) {
          for (k=0; k<n_[2]; k++) {
              const ArrayExtentData &d = data(i,j,k);
              if (d.size()) {
                  ExEnv::out() << node0 << indent
                       << i << " " << j << " " << k << ":" << endl;
                  for (l=0; l<d.size(); l++) {
                      ExEnv::out() << node0 << indent
                           << "  " << d[l].shell << " " << d[l].bound << endl;
                    }
                }
            }
        }
    }

  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
