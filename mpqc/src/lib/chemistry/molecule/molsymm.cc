//
// molsymm.cc
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

#include <util/misc/math.h>

#include <util/misc/formio.h>
#include <chemistry/molecule/molecule.h>

using namespace std;
using namespace sc;

#undef DEBUG

void
Molecule::clear_symmetry_info()
{
  for (int i=0; i<nuniq_; i++) {
    delete[] equiv_[i];
    }
  delete[] equiv_;
  delete[] nequiv_;
  delete[] atom_to_uniq_;
  nuniq_ = 0;
  equiv_ = 0;
  nequiv_ = 0;
  atom_to_uniq_ = 0;
}

void
Molecule::init_symmetry_info(double tol)
{
  if (equiv_)
    clear_symmetry_info();
  
  if (natom() == 0) {
    nuniq_ = 0;
    equiv_ = 0;
    nequiv_ = 0;
    atom_to_uniq_ = 0;
    return;
    }

  nequiv_ = new int[natom()];
  atom_to_uniq_ = new int[natom()];
  equiv_ = new int*[natom()];

  if (!strcmp(point_group()->symbol(),"c1")) {
    nuniq_ = natom();
    for (int i=0; i < natom(); i++) {
      nequiv_[i]=1;
      equiv_[i]=new int[1];
      equiv_[i][0]=i;
      atom_to_uniq_[i]=i;
      }
    return;
  }

  // the first atom is always unique
  nuniq_ = 1;
  nequiv_[0]=1;
  equiv_[0] = new int[1];
  equiv_[0][0]=0;
  atom_to_uniq_[0]=0;

  CharacterTable ct = point_group()->char_table();

  SCVector3 ac;
  SymmetryOperation so;
  SCVector3 np;

  // find the equivalent atoms
  int i;
  for (i=1; i < natom(); i++) {
    ac = r(i);
    int i_is_unique=1;
    int i_equiv=0;

    // apply all symmetry ops in the group to the atom
    for (int g=0; g < ct.order(); g++) {
      so = ct.symm_operation(g);
      for (int ii=0; ii < 3; ii++) {
        np[ii]=0;
        for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * ac[jj];
        }

      // see if the transformed atom is equivalent to a unique atom
      for (int j=0; j<nuniq_; j++) {
        int uniq = equiv_[j][0];
        SCVector3 aj(r(uniq));
        if (np.dist(aj) < tol
            && Z(uniq) == Z(i)
            && fabs(charge(uniq)-charge(i)) < tol
            && fabs(mass(uniq)-mass(i)) < tol) {
          i_is_unique = 0;
          i_equiv = j;
          break;
          }
        }
      }
    if (i_is_unique) {
      nequiv_[nuniq_]=1;
      equiv_[nuniq_]=new int[1];
      equiv_[nuniq_][0]=i;
      atom_to_uniq_[i] = nuniq_;
      nuniq_++;
      }
    else {
      int *tmp = new int[nequiv_[i_equiv]+1];
      memcpy(tmp,equiv_[i_equiv],nequiv_[i_equiv]*sizeof(int));
      delete[] equiv_[i_equiv];
      equiv_[i_equiv] = tmp;
      equiv_[i_equiv][nequiv_[i_equiv]] = i;
      nequiv_[i_equiv]++;
      atom_to_uniq_[i] = i_equiv;
      }
    }

  // The first atom in the equiv list is considered the primary unique
  // atom.  Just to make things look pretty, make the atom with the most
  // zeros in its x, y, z coordinate the unique atom.  Nothing else should
  // rely on this being done.
  double ztol=1.0e-5;
  for (i=0; i < nuniq_; i++) {
    int maxzero = 0;
    int jmaxzero = 0;
    for (int j=0; j<nequiv_[i]; j++) {
      int nzero = 0;
      for (int k=0; k<3; k++) if (fabs(r(equiv_[i][j],k)) < ztol) nzero++;
      if (nzero > maxzero) {
        maxzero = nzero;
        jmaxzero = j;
        }
      }
    int tmp = equiv_[i][jmaxzero];
    equiv_[i][jmaxzero] = equiv_[i][0];
    equiv_[i][0] = tmp;
    }
}

int
Molecule::has_inversion(SCVector3 &origin, double tol) const
{
  for (int i=0; i<natom(); i++) {
    SCVector3 inverted = origin-(SCVector3(r(i))-origin);
    int atom = atom_at_position(inverted.data(), tol);
    if (atom < 0
        || Z(atom) != Z(i)
        || fabs(charge(atom)-charge(i)) > tol
        || fabs(mass(atom)-mass(i)) > tol) {
      return 0;
      }
    }
  return 1;
}

int
Molecule::is_plane(SCVector3 &origin, SCVector3 &uperp, double tol) const
{
  for (int i=0; i<natom(); i++) {
    SCVector3 A = SCVector3(r(i))-origin;
    SCVector3 Apar = uperp.dot(A) * uperp;
    SCVector3 Aperp = A - Apar;
    A = (Aperp - Apar) + origin;
    int atom = atom_at_position(A.data(), tol);
    if (atom < 0
        || Z(atom) != Z(i)
        || fabs(charge(atom)-charge(i)) > tol
        || fabs(mass(atom)-mass(i)) > tol) {
      //ExEnv::outn() << "  is_plane: rejected (atom " << i << ")" << endl;
      return 0;
      }
    }
  return 1;
}

int
Molecule::is_axis(SCVector3 &origin, SCVector3 &axis,
                  int order, double tol) const
{
  // loop through atoms to see if axis is a c2 axis
  for (int i=0; i<natom(); i++) {
    SCVector3 A = SCVector3(r(i))-origin;
    for (int j=1; j<order; j++) {
      SCVector3 R = A;
      R.rotate(j*2.0*M_PI/order, axis);
      R += origin;
      int atom = atom_at_position(R.data(), tol);
      if (atom < 0
          || Z(atom) != Z(i)
          || fabs(charge(atom)-charge(i)) > tol
          || fabs(mass(atom)-mass(i)) > tol) {
        //ExEnv::outn() << "  is_axis: rejected (atom " << i << ")" << endl;
        return 0;
        }
      }
    }
  return 1;
}

enum AxisName { XAxis, YAxis, ZAxis };

static AxisName
like_world_axis(SCVector3 &axis,
                const SCVector3 &worldxaxis,
                const SCVector3 &worldyaxis,
                const SCVector3 &worldzaxis
                )
{
  AxisName like;
  double xlikeness = fabs(axis.dot(worldxaxis));
  double ylikeness = fabs(axis.dot(worldyaxis));
  double zlikeness = fabs(axis.dot(worldzaxis));
  if (xlikeness > ylikeness && xlikeness > zlikeness) {
    like = XAxis;
    if (axis.dot(worldxaxis) < 0) axis = - axis; 
    }
  else if (ylikeness > zlikeness) {
    like = YAxis;
    if (axis.dot(worldyaxis) < 0) axis = - axis; 
    }
  else {
    like = ZAxis;
    if (axis.dot(worldzaxis) < 0) axis = - axis; 
    }
  return like;
}

Ref<PointGroup>
Molecule::highest_point_group(double tol) const
{
  int i,j;

  SCVector3 com = center_of_mass();

  SCVector3 worldzaxis(0.,0.,1.);
  SCVector3 worldxaxis(1.,0.,0.);
  SCVector3 worldyaxis(0.,1.,0.);

  int linear,planar;
  is_linear_planar(linear,planar,tol);

  int have_inversion = has_inversion(com,tol);

  // check for C2 axis
  SCVector3 c2axis;
  int have_c2axis = 0;
  if (natom() < 2) {
    have_c2axis = 1;
    c2axis = SCVector3(0.0,0.0,1.0);
    }
  else if (linear) {
    have_c2axis = 1;
    c2axis = SCVector3(r(1)) - SCVector3(r(0));
    c2axis.normalize();
    }
  else if (planar && have_inversion) {
    // there is a c2 axis that won't be found using the usual algorithm.
    // find two noncolinear atom-atom vectors (we know that linear==0)
    SCVector3 BA = SCVector3(r(1))-SCVector3(r(0));
    BA.normalize();
    for (i=2; i<natom(); i++) {
      SCVector3 CA = SCVector3(r(i))-SCVector3(r(0));
      CA.normalize();
      SCVector3 BAxCA = BA.cross(CA);
      if (BAxCA.norm() > tol) {
        have_c2axis = 1;
        BAxCA.normalize();
        c2axis = BAxCA;
        break;
        }
      }
    }
  else {
    // loop through pairs of atoms to find C2 axis candidates
    for (i=0; i<natom(); i++) {
      SCVector3 A = SCVector3(r(i))-com;
      double AdotA = A.dot(A);
      for (j=0; j<=i; j++) {
        // the atoms must be identical
        if (Z(i) != Z(j) || fabs(mass(i)-mass(j)) > tol) continue;
        SCVector3 B = SCVector3(r(j))-com;
        // the atoms must be the same distance from the com
        if (fabs(AdotA - B.dot(B)) > tol) continue;
        SCVector3 axis = A+B;
        // atoms colinear with the com don't work
        if (axis.norm() < tol) continue;
        axis.normalize();
        if (is_axis(com,axis,2,tol)) {
          have_c2axis = 1;
          c2axis = axis;
          goto found_c2axis;
          }
        }
      }
    }
  found_c2axis:

  AxisName c2like = ZAxis;
  if (have_c2axis) {
    // try to make the sign of the axis correspond to one of the world axes
    c2like = like_world_axis(c2axis,worldxaxis,worldyaxis,worldzaxis);
    }

  // check for C2 axis perp to first C2 axis
  SCVector3 c2axisperp;
  int have_c2axisperp = 0;
  if (have_c2axis) {
    if (natom() < 2) {
      have_c2axisperp = 1;
      c2axisperp = SCVector3(1.0,0.0,0.0);
      }
    else if (linear) {
      if (have_inversion) {
        have_c2axisperp = 1;
        c2axisperp = c2axis.perp_unit(SCVector3(0.0,0.0,1.0));
        }
      }
    else {
      // loop through pairs of atoms to find C2 axis candidates
      for (i=0; i<natom(); i++) {
        SCVector3 A = SCVector3(r(i))-com;
        double AdotA = A.dot(A);
        for (j=0; j<i; j++) {
          // the atoms must be identical
          if (Z(i) != Z(j) || fabs(mass(i)-mass(j)) > tol) continue;
          SCVector3 B = SCVector3(r(j))-com;
          // the atoms must be the same distance from the com
          if (fabs(AdotA - B.dot(B)) > tol) continue;
          SCVector3 axis = A+B;
          // atoms colinear with the com don't work
          if (axis.norm() < tol) continue;
          axis.normalize();
          // if axis is not perp continue
          if (fabs(axis.dot(c2axis)) > tol) continue;
          if (is_axis(com,axis,2,tol)) {
            have_c2axisperp = 1;
            c2axisperp = axis;
            goto found_c2axisperp;
            }
          }
        }
      }
    }
  found_c2axisperp:

  AxisName c2perplike;
  if (have_c2axisperp) {
    // try to make the sign of the axis correspond to one of the world axes
    c2perplike = like_world_axis(c2axisperp,worldxaxis,worldyaxis,worldzaxis);

    // try to make c2axis the z axis
    if (c2perplike == ZAxis) {
      SCVector3 tmpv = c2axisperp;
      tmpv = c2axisperp; c2axisperp = c2axis; c2axis = tmpv;
      c2perplike = c2like;
      c2like = ZAxis;
      }
    if (c2like != ZAxis) {
      if (c2like == XAxis) c2axis = c2axis.cross(c2axisperp);
      else c2axis = c2axisperp.cross(c2axis);
      c2like = like_world_axis(c2axis,worldxaxis,worldyaxis,worldzaxis);
      }

    // try to make c2axisperp like the x axis
    if (c2perplike == YAxis) {
      c2axisperp = c2axisperp.cross(c2axis);
      c2perplike = like_world_axis(c2axisperp,
                                   worldxaxis,worldyaxis,worldzaxis);
      }
    }

  // check for vertical plane
  int have_sigmav = 0;
  SCVector3 sigmav;
  if (have_c2axis) {
    if (natom() < 2) {
      have_sigmav = 1;
      sigmav = c2axisperp;
      }
    else if (linear) {
      have_sigmav = 1;
      if (have_c2axisperp) {
        sigmav = c2axisperp;
        }
      else {
        sigmav = c2axis.perp_unit(SCVector3(0.0,0.0,1.0));
        }
      }
    else {
      // loop through pairs of atoms to find sigma v plane candidates
      for (i=0; i<natom(); i++) {
        SCVector3 A = SCVector3(r(i))-com;
        double AdotA = A.dot(A);
        // the second atom can equal i because i might be in the plane.
        for (j=0; j<=i; j++) {
          // the atoms must be identical
          if (Z(i) != Z(j) || fabs(mass(i)-mass(j)) > tol) continue;
          SCVector3 B = SCVector3(r(j))-com;
          // the atoms must be the same distance from the com
          if (fabs(AdotA - B.dot(B)) > tol) continue;
          SCVector3 inplane = B+A;
          double norm_inplane = inplane.norm();
          if (norm_inplane < tol) continue;
          inplane *= 1.0/norm_inplane;
          SCVector3 perp = c2axis.cross(inplane);
          double norm_perp = perp.norm();
          if (norm_perp < tol) continue;
          perp *= 1.0/norm_perp;
          if (is_plane(com,perp,tol)) {
            have_sigmav = 1;
            sigmav = perp;
            goto found_sigmav;
            }
          }
        }
      }
    }
  found_sigmav:

  if (have_sigmav) {
    // try to make the sign of the oop vec correspond to one of the world axes
    int sigmavlike = like_world_axis(sigmav,worldxaxis,worldyaxis,worldzaxis);

    // choose sigmav to be the world x axis, if possible
    if (c2like == ZAxis && sigmavlike == YAxis) {
      sigmav = sigmav.cross(c2axis);
      }
    else if (c2like == YAxis && sigmavlike == ZAxis) {
      sigmav = c2axis.cross(sigmav);
      }
    }

  // under certain conditions i need to know if there is any sigma plane
  int have_sigma = 0;
  SCVector3 sigma;
  if (!have_inversion && !have_c2axis) {
    if (planar) {
      // find two noncolinear atom-atom vectors
      // we know that linear==0 since !have_c2axis
      SCVector3 BA = SCVector3(r(1))-SCVector3(r(0));
      BA.normalize();
      for (i=2; i<natom(); i++) {
        SCVector3 CA = SCVector3(r(i))-SCVector3(r(0));
        CA.normalize();
        SCVector3 BAxCA = BA.cross(CA);
        if (BAxCA.norm() > tol) {
          have_sigma = 1;
          BAxCA.normalize();
          sigma = BAxCA;
          break;
          }
        }
      }
    else {
      // loop through pairs of atoms to construct trial planes
      for (i=0; i<natom(); i++) {
        SCVector3 A = SCVector3(r(i))-com;
        double AdotA = A.dot(A);
        for (j=0; j<i; j++) {
          //ExEnv::outn() << "sigma atoms = " << i << ", " << j << endl;
          // the atoms must be identical
          if (Z(i) != Z(j) || fabs(mass(i)-mass(j)) > tol) continue;
          SCVector3 B = SCVector3(r(j))-com;
          double BdotB = B.dot(B);
          // the atoms must be the same distance from the com
          if (fabs(AdotA - BdotB) > tol) continue;
          SCVector3 perp = B-A;
          double norm_perp = perp.norm();
          if (norm_perp < tol) {
            //ExEnv::outn() << "  rejected (atoms at same point?)" << endl;
            continue;
            }
          perp *= 1.0/norm_perp;
          if (is_plane(com,perp,tol)) {
            have_sigma = 1;
            sigma = perp;
            goto found_sigma;
            }
          }
        }
      }
    }
  found_sigma:

  if (have_sigma) {
    // try to make the sign of the oop vec correspond to one of the world axes
    double xlikeness = fabs(sigma.dot(worldxaxis));
    double ylikeness = fabs(sigma.dot(worldyaxis));
    double zlikeness = fabs(sigma.dot(worldzaxis));
    if (xlikeness > ylikeness && xlikeness > zlikeness) {
      if (sigma.dot(worldxaxis) < 0) sigma = - sigma; 
      }
    else if (ylikeness > zlikeness) {
      if (sigma.dot(worldyaxis) < 0) sigma = - sigma; 
      }
    else {
      if (sigma.dot(worldzaxis) < 0) sigma = - sigma; 
      }
    }

#ifdef DEBUG
  ExEnv::out0()
       << indent << "highest point group:" << endl
       << indent << "  linear          = " << linear << endl
       << indent << "  planar          = " << planar << endl
       << indent << "  have_inversion  = " << have_inversion << endl
       << indent << "  have_c2axis     = " << have_c2axis << endl
       << indent << "  have_c2axisperp = " << have_c2axisperp << endl
       << indent << "  have_sigmav     = " << have_sigmav << endl
       << indent << "  have_sigma      = " << have_sigma << endl;

  if (have_c2axis)
    ExEnv::out0() << indent << "  c2axis      = " << c2axis << endl;
  if (have_c2axisperp)
    ExEnv::out0() << indent << "  c2axisperp  = " << c2axisperp << endl;
  if (have_sigmav)
    ExEnv::out0() << indent << "  sigmav      = " << sigmav << endl;
  if (have_sigma)
    ExEnv::out0() << indent << "  sigma       = " << sigma << endl;
#endif

  // Find the three axes for the symmetry frame
  SCVector3 xaxis = worldxaxis;
  SCVector3 yaxis;
  SCVector3 zaxis = worldzaxis;;
  if (have_c2axis) {
    zaxis = c2axis;
    if (have_sigmav) {
      xaxis = sigmav;
      }
    else if (have_c2axisperp) {
      xaxis = c2axisperp;
      }
    else {
      // any axis orthogonal to the zaxis will do
      xaxis = zaxis.perp_unit(zaxis);
      }
    }
  else if (have_sigma) {
    zaxis = sigma;
    xaxis = zaxis.perp_unit(zaxis);
    }
  // the y axis is then -x cross z
  yaxis = - xaxis.cross(zaxis);

#ifdef DEBUG
  ExEnv::outn() << "X: " << xaxis << endl;
  ExEnv::outn() << "Y: " << yaxis << endl;
  ExEnv::outn() << "Z: " << zaxis << endl;
#endif

  SymmetryOperation frame;
  SCVector3 origin;
  for (i=0; i<3; i++) {
    frame(i,0) = xaxis[i];
    frame(i,1) = yaxis[i];
    frame(i,2) = zaxis[i];
    origin[i] = com[i];
    }

#ifdef DEBUG
  ExEnv::out0() << "frame:" << endl;
  frame.print(ExEnv::out0());

  ExEnv::out0() << "origin:" << endl;
  origin.print(ExEnv::out0());
#endif

  Ref<PointGroup> pg;
  if (have_inversion) {
    if (have_c2axis) {
      if (have_sigmav) {
        pg = new PointGroup("d2h",frame,origin);
        }
      else {
        pg = new PointGroup("c2h",frame,origin);
        }
      }
    else {
      pg = new PointGroup("ci",frame,origin);
      }
    }
  else {
    if (have_c2axis) {
      if (have_sigmav) {
        pg = new PointGroup("c2v",frame,origin);
        }
      else {
        if (have_c2axisperp) {
          pg = new PointGroup("d2",frame,origin);
          }
        else {
          pg = new PointGroup("c2",frame,origin);
          }
        }
      }
    else {
      if (have_sigma) {
        pg = new PointGroup("cs",frame,origin);
        }
      else {
        pg = new PointGroup("c1",frame,origin);
        }
      }
    }

  return pg;
}

int
Molecule::is_linear(double tol) const
{
  int linear, planar;

  is_linear_planar(linear,planar,tol);

  return linear;
}

int
Molecule::is_planar(double tol) const
{
  int linear, planar;

  is_linear_planar(linear,planar,tol);

  return planar;
}

void
Molecule::is_linear_planar(int &linear, int &planar, double tol) const
{
  if (natom() < 3) {
    linear = 1;
    planar = 1;
    return;
    }

  // find three atoms not on the same line
  SCVector3 A = r(0);
  SCVector3 B = r(1);
  SCVector3 BA = B-A;
  BA.normalize();
  SCVector3 CA;

  int i;
  double min_BAdotCA = 1.0;
  for (i=2; i<natom(); i++) {
    SCVector3 tmp = SCVector3(r(i))-A;
    tmp.normalize();
    if (fabs(BA.dot(tmp)) < min_BAdotCA) {
      CA = tmp;
      min_BAdotCA = fabs(BA.dot(tmp));
      }
    }
  if (min_BAdotCA >= 1.0 - tol) {
    linear = 1;
    planar = 1;
    return;
    }

  linear = 0;
  if (natom() < 4) {
    planar = 1;
    return;
    }

  // check for nontrivial planar molecules
  SCVector3 BAxCA = BA.cross(CA);
  BAxCA.normalize();
  for (i=2; i<natom(); i++) {
    SCVector3 tmp = SCVector3(r(i))-A;
    if (fabs(tmp.dot(BAxCA)) > tol) {
      planar = 0;
      return;
      }
    }
  planar = 1;
  return;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
