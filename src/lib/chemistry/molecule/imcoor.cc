//
// imcoor.cc
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

#include <cmath>
#include <memory>

#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <chemistry/molecule/localdef.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/simple.h>

#include <util/container/bitarray.h>

using namespace std;
using namespace sc;

#define DEFAULT_SIMPLE_TOLERANCE 1.0e-3

///////////////////////////////////////////////////////////////////////////
// members of IntMolecularCoor

static ClassDesc IntMolecularCoor_cd(
  typeid(IntMolecularCoor),"IntMolecularCoor",6,"public MolecularCoor",
  0, 0, 0);

IntMolecularCoor::IntMolecularCoor(Ref<Molecule>&mol):
  MolecularCoor(mol),
  update_bmat_(0),
  only_totally_symmetric_(1),
  symmetry_tolerance_(1.0e-5),
  simple_tolerance_(DEFAULT_SIMPLE_TOLERANCE),
  coordinate_tolerance_(1.0e-7),
  cartesian_tolerance_(1.0e-12),
  scale_bonds_(1.0),
  scale_bends_(1.0),
  scale_tors_(1.0),
  scale_outs_(1.0),
  given_fixed_values_(0),
  decouple_bonds_(0),
  decouple_bends_(0),
  max_update_steps_(100),
  max_update_disp_(0.5),
  form_print_simples_(0),
  form_print_variable_(0),
  form_print_constant_(0),
  form_print_molecule_(0)
{
  new_coords();
  generator_ = new IntCoorGen(mol);
}

IntMolecularCoor::IntMolecularCoor(const Ref<KeyVal>& keyval):
  MolecularCoor(keyval),
  update_bmat_(0),
  only_totally_symmetric_(1),
  symmetry_tolerance_(1.0e-5),
  simple_tolerance_(DEFAULT_SIMPLE_TOLERANCE),
  coordinate_tolerance_(1.0e-7),
  cartesian_tolerance_(1.0e-12),
  scale_bonds_(1.0),
  scale_bends_(1.0),
  scale_tors_(1.0),
  scale_outs_(1.0),
  decouple_bonds_(0),
  decouple_bends_(0)
{
  // intialize the coordinate sets
  new_coords();

  // actually read the keyval info
  read_keyval(keyval);
}

IntMolecularCoor::IntMolecularCoor(StateIn& s):
  MolecularCoor(s)
{
  generator_ << SavableState::restore_state(s);

  if (s.version(::class_desc<IntMolecularCoor>()) >= 3) {
      s.get(decouple_bonds_);
      s.get(decouple_bends_);
    }
  else {
      decouple_bonds_ = 0;
      decouple_bends_ = 0;
    }
  
  if (s.version(::class_desc<IntMolecularCoor>()) >= 2) {
    s.get(max_update_steps_);
    s.get(max_update_disp_);
    s.get(given_fixed_values_);
  } else {
    max_update_steps_ = 100;
    max_update_disp_ = 0.5;
    given_fixed_values_ = 0;
  }
  
  if (s.version(::class_desc<IntMolecularCoor>()) >= 4) {
    s.get(form_print_simples_);
    s.get(form_print_variable_);
    s.get(form_print_constant_);
  } else {
    form_print_simples_ = 0;
    form_print_variable_ = 0;
    form_print_constant_ = 0;
  }
  
  if (s.version(::class_desc<IntMolecularCoor>()) >= 5) {
    s.get(form_print_molecule_);
  } else {
    form_print_molecule_ = 0;
  }

  dim_ << SavableState::restore_state(s);
  dvc_ << SavableState::restore_state(s);

  all_ << SavableState::restore_state(s);

  variable_ << SavableState::restore_state(s);
  constant_ << SavableState::restore_state(s);

  fixed_ << SavableState::restore_state(s);
  followed_ << SavableState::restore_state(s);

  if (s.version(::class_desc<IntMolecularCoor>()) >= 6)
      watched_ << SavableState::restore_state(s);

  bonds_ << SavableState::restore_state(s);
  bends_ << SavableState::restore_state(s);
  tors_ << SavableState::restore_state(s);
  outs_ << SavableState::restore_state(s);
  extras_ << SavableState::restore_state(s);

  s.get(update_bmat_);
  s.get(only_totally_symmetric_);
  s.get(scale_bonds_);
  s.get(scale_bends_);
  s.get(scale_tors_);
  s.get(scale_outs_);
  s.get(simple_tolerance_);
  s.get(symmetry_tolerance_);
  s.get(coordinate_tolerance_);
  s.get(cartesian_tolerance_);
}

void
IntMolecularCoor::new_coords()
{
  // intialize the coordinate sets
  all_ = new SetIntCoor; // all redundant coors
  variable_ = new SetIntCoor; // internal coors to be varied
  constant_ = new SetIntCoor; // internal coors to be fixed
  bonds_ = new SetIntCoor;
  bends_ = new SetIntCoor;
  tors_ = new SetIntCoor;
  outs_ = new SetIntCoor;
  extras_ = new SetIntCoor;
  fixed_ = new SetIntCoor;
  followed_ = 0;
  watched_ = 0;
}

void
IntMolecularCoor::read_keyval(const Ref<KeyVal>& keyval)
{
  variable_ << keyval->describedclassvalue("variable");
  if (variable_.null()) variable_ = new SetIntCoor;
  fixed_ << keyval->describedclassvalue("fixed");
  if (fixed_.null()) fixed_ = new SetIntCoor;
  followed_ << keyval->describedclassvalue("followed");
  watched_ << keyval->describedclassvalue("watched");

  decouple_bonds_ = keyval->booleanvalue("decouple_bonds");
  decouple_bends_ = keyval->booleanvalue("decouple_bends");

  given_fixed_values_ = keyval->booleanvalue("have_fixed_values");

  max_update_steps_ = keyval->intvalue("max_update_steps");
  if (keyval->error() != KeyVal::OK) max_update_steps_ = 100;

  max_update_disp_ = keyval->doublevalue("max_update_disp");
  if (keyval->error() != KeyVal::OK) max_update_disp_ = 0.5;

  generator_ << keyval->describedclassvalue("generator");

  if (generator_.null()) {
      // the extra_bonds list is given as a vector of atom numbers
      // (atom numbering starts at 1)
      int nextra_bonds = keyval->count("extra_bonds");
      nextra_bonds /= 2;
      vector<int> extra_bonds;
      if (nextra_bonds) {
          extra_bonds.resize(nextra_bonds*2);
          for (int i=0; i<nextra_bonds*2; i++) {
              extra_bonds[i] = keyval->intvalue("extra_bonds",i);
              if (keyval->error() != KeyVal::OK) {
                  throw InputError("missing an expected integer value",
                                   __FILE__, __LINE__, "extra_bonds", 0,
                                   class_desc());
                }
            }
        }
      generator_ = new IntCoorGen(molecule_, nextra_bonds, extra_bonds.size() ? &(extra_bonds[0]) : 0);
    }
          

  update_bmat_ = keyval->booleanvalue("update_bmat");

  only_totally_symmetric_ = keyval->booleanvalue("only_totally_symmetric");
  if (keyval->error() != KeyVal::OK) only_totally_symmetric_ = 1;

  double tmp;
  tmp = keyval->doublevalue("scale_bonds");
  if (keyval->error() == KeyVal::OK) scale_bonds_ = tmp;
  tmp = keyval->doublevalue("scale_bends");
  if (keyval->error() == KeyVal::OK) scale_bends_ = tmp;
  tmp = keyval->doublevalue("scale_tors");
  if (keyval->error() == KeyVal::OK) scale_tors_ = tmp;
  tmp = keyval->doublevalue("scale_outs");
  if (keyval->error() == KeyVal::OK) scale_outs_ = tmp;
  tmp = keyval->doublevalue("symmetry_tolerance");
  if (keyval->error() == KeyVal::OK) symmetry_tolerance_ = tmp;
  tmp = keyval->doublevalue("simple_tolerance");
  if (keyval->error() == KeyVal::OK) simple_tolerance_ = tmp;
  tmp = keyval->doublevalue("coordinate_tolerance");
  if (keyval->error() == KeyVal::OK) coordinate_tolerance_ = tmp;
  tmp = keyval->doublevalue("cartesian_tolerance");
  if (keyval->error() == KeyVal::OK) cartesian_tolerance_ = tmp;

  form_print_simples_ = keyval->booleanvalue("form:print_simple");
  if (keyval->error() != KeyVal::OK) form_print_simples_ = 0;
  form_print_variable_ = keyval->booleanvalue("form:print_variable");
  if (keyval->error() != KeyVal::OK) form_print_variable_ = 0;
  form_print_constant_ = keyval->booleanvalue("form:print_constant");
  if (keyval->error() != KeyVal::OK) form_print_constant_ = 0;
  form_print_molecule_ = keyval->booleanvalue("form:print_molecule");
  if (keyval->error() != KeyVal::OK) form_print_molecule_ = 0;
}

void
IntMolecularCoor::init()
{
  Ref<SetIntCoor> redundant = new SetIntCoor;
  generator_->generate(redundant);

  // sort out the simple coordinates by type
  int i;
  for (i=0; i<redundant->n(); i++) {
      Ref<IntCoor> coor = redundant->coor(i);
      if (coor->class_desc()
          == ::class_desc<StreSimpleCo>()) {
          bonds_->add(coor);
        }
      else if (coor->class_desc() == ::class_desc<BendSimpleCo>()
               || coor->class_desc() == ::class_desc<LinIPSimpleCo>()
               || coor->class_desc() == ::class_desc<LinOPSimpleCo>()) {
          bends_->add(coor);
        }
      else if (coor->class_desc()
               == ::class_desc<TorsSimpleCo>()
               || coor->class_desc()
               == ::class_desc<ScaledTorsSimpleCo>()) {
          tors_->add(coor);
        }
      else if (coor->class_desc()
               == ::class_desc<OutSimpleCo>()) {
          outs_->add(coor);
        }
      else {
          extras_->add(coor);
        }
    }

  all_->add(bonds_);
  all_->add(bends_);
  all_->add(tors_);
  all_->add(outs_);
  all_->add(extras_);

  // don't let form_coordinates create new variables coordinates
  // if they were given by the user
  int keep_variable = (variable_->n() != 0);

  if (given_fixed_values_) {
      // save the given coordinate values
      RefSCDimension original_dfixed
          = new SCDimension(fixed_->n(),"Nfix");
      RefSCVector given_fixed_coords(original_dfixed,matrixkit_);
      for (i=0; i<original_dfixed.n(); i++) {
          given_fixed_coords(i) = fixed_->coor(i)->value();
        }

      // find the current fixed coordinates
      RefSCVector current_fixed_coords(original_dfixed,matrixkit_);
      fixed_->update_values(molecule_);
      for (i=0; i<original_dfixed.n(); i++) {
          current_fixed_coords(i) = fixed_->coor(i)->value();
        }

      // the difference between current fixed and desired fixed
      RefSCVector diff_fixed_coords = given_fixed_coords-current_fixed_coords;

      // break up the displacement into several manageable steps
      double maxabs = diff_fixed_coords.maxabs();
      int nstep = int(maxabs/max_update_disp_) + 1;
      diff_fixed_coords.scale(1.0/nstep);
      ExEnv::out0() << indent << "IntMolecularCoor: "
           << "displacing fixed coordinates to the requested values in "
           << nstep << " steps\n";
      for (int istep=1; istep<=nstep; istep++) {
          form_coordinates(keep_variable);

          dim_ = new SCDimension(variable_->n(), "Nvar");
          dvc_ = new SCDimension(variable_->n()+constant_->n(),
                             "Nvar+Nconst");

          RefSCVector new_internal_coordinates(dvc_,matrixkit_);
          for (i=0; i<variable_->n(); i++) {
              new_internal_coordinates(i) = variable_->coor(i)->value();
            }
          int j;
          for (j=0; j<original_dfixed.n(); j++,i++) {
              new_internal_coordinates(i)
                  = current_fixed_coords(j)+istep*double(diff_fixed_coords(j));
            }
          for (; j<constant_->n(); i++,j++) {
              new_internal_coordinates(i) = constant_->coor(j)->value();
            }

          all_to_cartesian(molecule_, new_internal_coordinates);
        }

      // make sure that the coordinates have exactly the
      // original values to avoid round-off error
      for (i=0; i<original_dfixed.n(); i++) {
          fixed_->coor(i)->set_value(given_fixed_coords(i));
        }
    }

  form_coordinates(keep_variable);

  dim_ = new SCDimension(variable_->n(), "Nvar");
  dvc_ = new SCDimension(variable_->n()+constant_->n(),
                         "Nvar+Nconst");

#if 0 // this will always think the rank has changed with redundant coordinates
    {
      const double epsilon = 0.001;

      // compute the condition number
      RefSCMatrix B(dim_, dnatom3_,matrixkit_);
      variable_->bmat(molecule_, B);

      // Compute the singular value decomposition of B
      RefSCMatrix U(dim_,dim_,matrixkit_);
      RefSCMatrix V(dnatom3_,dnatom3_,matrixkit_);
      RefSCDimension min;
      if (dnatom3_.n()<dim_.n()) min = dnatom3_;
      else min = dim_;
      int nmin = min.n();
      RefDiagSCMatrix sigma(min,matrixkit_);
      B.svd(U,sigma,V);

      // Compute the epsilon rank of B
      int i, rank = 0;
      for (i=0; i<nmin; i++) {
          if (sigma(i) > epsilon) rank++;
        }

      if (rank != dim_.n()) {
          ExEnv::out0() << indent << "IntMolecularCoor::init: rank changed\n";
          sigma.print("sigma");
        }

      double kappa2 = sigma(0)/sigma(dim_.n()-1);

      ExEnv::out0() << indent
           << scprintf("IntMolecularCoor::init: condition number = %14.8f\n",
                       kappa2);
    }
#endif

  if (watched_.nonnull()) {
      ExEnv::out0() << endl
           << indent << "Watched coordinate(s):\n" << incindent;
      watched_->update_values(molecule_);
      watched_->print_details(molecule_,ExEnv::out0());
      ExEnv::out0() << decindent;
    }
}

static int
count_nonzero(const RefSCVector &vec, double eps)
{
  int nz=0, i, n=vec.n();
  for (i=0; i<n; i++) {
      if (fabs(vec(i)) > eps) nz++;
    }
  return nz;
}

static RefSymmSCMatrix
form_partial_K(const Ref<SetIntCoor>& coor, Ref<Molecule>& molecule,
               const RefSCVector& geom,
               double epsilon,
               const RefSCDimension& dnatom3,
               const Ref<SCMatrixKit>& matrixkit,
               RefSCMatrix& projection,
               RefSCVector& totally_symmetric,
               RefSCMatrix& K,int debug)
{
  if (debug) {
      ExEnv::out0() << indent << "form_partial_K:" << endl;
      ExEnv::out0() << incindent;
    }

  // Compute the B matrix for the coordinates
  RefSCDimension dcoor = new SCDimension(coor->n());
  RefSCMatrix B(dcoor, dnatom3,matrixkit);
  coor->bmat(molecule, B);

  if (debug) B.print("B");

  // Project out the previously discovered internal coordinates
  if (projection.nonnull()) {
      B = B * projection;
      if (debug) B.print("Projected B");
    }

  // Compute the singular value decomposition of B
  RefSCMatrix U(dcoor,dcoor,matrixkit);
  RefSCMatrix V(dnatom3,dnatom3,matrixkit);
  RefSCDimension min;
  if (dnatom3.n()<dcoor.n()) min = dnatom3;
  else min = dcoor;
  int nmin = min.n();
  RefDiagSCMatrix sigma(min,matrixkit);
  B.svd(U,sigma,V);

  // Compute the epsilon rank of B
  int i, rank = 0;
  for (i=0; i<nmin; i++) {
      if (sigma(i) > epsilon) rank++;
    }

  if (debug)
      ExEnv::out0() << indent << "rank(" << epsilon << ",B) = " << rank
           << endl;

  RefSCMatrix SIGMA(dcoor, dnatom3,matrixkit);
  SIGMA.assign(0.0);
  for (i=0; i<nmin; i++) {
      SIGMA(i,i) = sigma(i);
    }

  // return if there are no new coordinates
  if (rank==0) {
      if (debug) ExEnv::out0() << decindent;
      return 0;
    }

  // Find an orthogonal matrix that spans the range of B
  RefSCMatrix Ur;
  RefSCDimension drank = new SCDimension(rank);
  if (rank) {
      Ur = matrixkit->matrix(dcoor,drank);
      Ur.assign_subblock(U,0, dcoor.n()-1, 0, drank.n()-1, 0, 0);
    }

  // Find an orthogonal matrix that spans the null space of B
  int rank_tilde = dnatom3.n() - rank;
  RefSCMatrix Vr_tilde;
  RefSCDimension drank_tilde = new SCDimension(rank_tilde);
  if (rank_tilde) {
      Vr_tilde = matrixkit->matrix(dnatom3,drank_tilde);
      Vr_tilde.assign_subblock(V,0, dnatom3.n()-1, 0, drank_tilde.n()-1,
                               0, drank.n());
    }

  // Find an orthogonal matrix that spans the null(B) perp
  RefSCMatrix Vr;
  if (rank) {
      Vr = matrixkit->matrix(dnatom3,drank);
      Vr.assign_subblock(V,0, dnatom3.n()-1, 0, drank.n()-1, 0, 0);
    }

  // compute the projection into the null space of B
  RefSymmSCMatrix proj_nullspace_B;
  if (rank_tilde) {
      proj_nullspace_B = matrixkit->symmmatrix(dnatom3);
      proj_nullspace_B.assign(0.0);
      proj_nullspace_B.accumulate_symmetric_product(Vr_tilde);
    }

  // compute the projection into the null(B) perp
  RefSymmSCMatrix proj_nullspace_B_perp;
  if (rank) {
      proj_nullspace_B_perp = matrixkit->symmmatrix(dnatom3);
      proj_nullspace_B_perp.assign(0.0);
      proj_nullspace_B_perp.accumulate_symmetric_product(Vr);
    }

  if (Ur.nonnull()) {
      // totally_symmetric will be nonzero for totally symmetric coordinates
      totally_symmetric = Ur.t() * B * geom;

      if (debug) {
          Ur.print("Ur");
          geom.print("geom");
          totally_symmetric.print("totally_symmetric = Ur.t()*B*geom");

          int ntotally_symmetric = count_nonzero(totally_symmetric,0.001);
          ExEnv::out0() << indent << "found " << ntotally_symmetric
               << " totally symmetric coordinates\n";
        }

      // compute the cumulative projection
      if (projection.null()) {
          projection = matrixkit->matrix(dnatom3,dnatom3);
          projection->unit();
        }
      projection = projection * proj_nullspace_B;
    }

  // give Ur to caller
  K = Ur;

  if (debug) ExEnv::out0() << decindent;

  return proj_nullspace_B_perp;
}

// this allocates storage for and computes K and is_totally_symmetric
void
IntMolecularCoor::form_K_matrix(RefSCDimension& dredundant,
                                RefSCDimension& dfixed,
                                RefSCMatrix& K,
                                int*& is_totally_symmetric)
{
  int i,j;

  // The cutoff for whether or not a coordinate is considered totally symmetric
  double ts_eps = 0.0001;

  // The geometry will be needed to check for totally symmetric
  // coordinates
  RefSCVector geom(dnatom3_,matrixkit_);
  for(i=0; i < geom.n()/3; i++) {
      geom(3*i  ) = molecule_->r(i,0);
      geom(3*i+1) = molecule_->r(i,1);
      geom(3*i+2) = molecule_->r(i,2);
    }

  RefSCDimension dcoor = new SCDimension(all_->n());

  // this keeps track of the total projection for the b matrices
  RefSCMatrix projection;
  if (dfixed.n()) {
      ExEnv::out0() << indent
           << "Forming fixed optimization coordinates:" << endl;
      RefSCMatrix Ktmp;
      RefSCVector totally_symmetric_fixed;
      RefSymmSCMatrix null_bfixed_perp
          = form_partial_K(fixed_, molecule_, geom, 0.001, dnatom3_,
                           matrixkit_, projection, totally_symmetric_fixed,
                           Ktmp,debug_);
      // require that the epsilon rank equal the number of fixed coordinates
      if (Ktmp.nrow() != dfixed.n()) {
          throw AlgorithmException("nfixed != rank",
                                   __FILE__, __LINE__,
                                   class_desc());
        }
      // check that fixed coordinates be totally symmetric
      //if (Ktmp.nrow() != count_nonzero(totally_symmetric_fixed, ts_eps)) {
      //    ExEnv::err0() << indent
      //         << scprintf("WARNING: only %d of %d fixed coordinates are"
      //                     " totally symmetric\n",
      //                     count_nonzero(totally_symmetric_fixed, ts_eps),
      //                     dfixed.n());
      //  }
    }

  ExEnv::out0() << indent << "Forming optimization coordinates:" << endl;

  int n_total = 0;

  RefSCVector totally_symmetric_bond;
  RefSCMatrix Kbond;
  if (decouple_bonds_) {
      ExEnv::out0() << indent << "looking for bonds" << endl;
      form_partial_K(bonds_, molecule_, geom, 0.1, dnatom3_, matrixkit_,
                     projection, totally_symmetric_bond, Kbond, debug_);
      if (Kbond.nonnull()) n_total += Kbond.ncol();
    }

  RefSCVector totally_symmetric_bend;
  RefSCMatrix Kbend;
  if (decouple_bends_) {
      ExEnv::out0() << indent << "looking for bends" << endl;
      form_partial_K(bends_, molecule_, geom, 0.1, dnatom3_, matrixkit_,
                     projection, totally_symmetric_bend, Kbend, debug_);
      if (Kbend.nonnull()) n_total += Kbend.ncol();
    }

  if (decouple_bonds_ || decouple_bends_) {
      ExEnv::out0() << indent << "looking for remaining coordinates" << endl;
    }
  RefSCVector totally_symmetric_all;
  RefSCMatrix Kall;
  // I hope the IntCoorSet keeps the ordering
  form_partial_K(all_, molecule_, geom, 0.001, dnatom3_, matrixkit_,
                 projection, totally_symmetric_all, Kall, debug_);
  if (Kall.nonnull()) n_total += Kall.ncol();

  // This requires that all_ coordinates is made up of first bonds,
  // bends, and finally the rest of the coordinates.
  RefSCDimension dtot = new SCDimension(n_total);
  K = matrixkit_->matrix(dcoor, dtot);
  K.assign(0.0);
  int istart=0, jstart=0;
  if (Kbond.nonnull()) {
      if (debug_) Kbond.print("Kbond");
      K.assign_subblock(Kbond, 0, Kbond.nrow()-1, 0, Kbond.ncol()-1, 0, 0);
      istart += Kbond.nrow();
      jstart += Kbond.ncol();
    }
  if (Kbend.nonnull()) {
      if (debug_) Kbend.print("Kbend");
      K.assign_subblock(Kbend, istart, istart+Kbend.nrow()-1,
                        jstart, jstart+Kbend.ncol()-1, 0, 0);
      istart += Kbend.nrow();
      jstart += Kbend.ncol();
    }
  if (Kall.nonnull()) {
      if (debug_) Kall.print("Kall");
      K.assign_subblock(Kall, 0, Kall.nrow()-1,
                        jstart, jstart+Kall.ncol()-1, 0, 0);
    }
  if (debug_) K.print("K");

  is_totally_symmetric = new int[K.ncol()];
  j=0;
  if (Kbond.nonnull()) {
      for (i=0; i<Kbond.ncol(); i++,j++) {
          if (fabs(totally_symmetric_bond(i)) > ts_eps)
              is_totally_symmetric[j] = 1;
          else is_totally_symmetric[j] = 0;
        }
    }
  if (Kbend.nonnull()) {
      for (i=0; i<Kbend.ncol(); i++,j++) {
          if (fabs(totally_symmetric_bend(i)) > ts_eps)
              is_totally_symmetric[j] = 1;
          else is_totally_symmetric[j] = 0;
        }
    }
  if (Kall.nonnull()) {
      for (i=0; i<Kall.ncol(); i++,j++) {
          if (fabs(totally_symmetric_all(i)) > ts_eps)
              is_totally_symmetric[j] = 1;
          else is_totally_symmetric[j] = 0;
        }
    }
}

IntMolecularCoor::~IntMolecularCoor()
{
}

void
IntMolecularCoor::save_data_state(StateOut&s)
{
  MolecularCoor::save_data_state(s);

  SavableState::save_state(generator_.pointer(),s);

  s.put(decouple_bonds_);
  s.put(decouple_bends_);

  s.put(max_update_steps_);
  s.put(max_update_disp_);
  s.put(given_fixed_values_);

  s.put(form_print_simples_);
  s.put(form_print_variable_);
  s.put(form_print_constant_);
  s.put(form_print_molecule_);

  SavableState::save_state(dim_.pointer(),s);
  SavableState::save_state(dvc_.pointer(),s);

  SavableState::save_state(all_.pointer(),s);
  
  SavableState::save_state(variable_.pointer(),s);
  SavableState::save_state(constant_.pointer(),s);

  SavableState::save_state(fixed_.pointer(),s);
  SavableState::save_state(followed_.pointer(),s);
  SavableState::save_state(watched_.pointer(),s);

  SavableState::save_state(bonds_.pointer(),s);
  SavableState::save_state(bends_.pointer(),s);
  SavableState::save_state(tors_.pointer(),s);
  SavableState::save_state(outs_.pointer(),s);
  SavableState::save_state(extras_.pointer(),s);

  s.put(update_bmat_);
  s.put(only_totally_symmetric_);
  s.put(scale_bonds_);
  s.put(scale_bends_);
  s.put(scale_tors_);
  s.put(scale_outs_);
  s.put(simple_tolerance_);
  s.put(symmetry_tolerance_);
  s.put(coordinate_tolerance_);
  s.put(cartesian_tolerance_);
}

RefSCDimension
IntMolecularCoor::dim()
{
  return dim_;
}

int
IntMolecularCoor::all_to_cartesian(const Ref<Molecule> &mol,
                                   RefSCVector&new_internal)
{
  // get a reference to Molecule for convenience
  Molecule& molecule = *(mol.pointer());

  // don't bother updating the bmatrix when the error is less than this
  const double update_tolerance = 1.0e-6;

  // compute the internal coordinate displacements
  RefSCVector old_internal(dvc_,matrixkit_);

  RefSCMatrix internal_to_cart_disp;
  double maxabs_cart_diff = 0.0;
  for (int step = 0; step < max_update_steps_; step++) {
      // compute the old internal coordinates
      all_to_internal(mol, old_internal);

      if (debug_) {
          ExEnv::out0()
               << indent << "Coordinates on step " << step << ":" << endl;
          variable_->print_details(0,ExEnv::out0());
        }

      // the displacements
      RefSCVector displacement = new_internal - old_internal;
      if (debug_ && step == 0) {
          displacement.print("Step 0 Internal Coordinate Displacement");
        }

      if ((update_bmat_ && maxabs_cart_diff>update_tolerance)
          || internal_to_cart_disp.null()) {
          if (debug_) {
              ExEnv::out0() << indent << "updating bmatrix" << endl;
            }

          int i;
          RefSCMatrix bmat(dvc_,dnatom3_,matrixkit_);

          // form the set of all coordinates
          Ref<SetIntCoor> variable_and_constant = new SetIntCoor();
          variable_and_constant->add(variable_);
          variable_and_constant->add(constant_);

          // form the bmatrix
          variable_and_constant->bmat(mol,bmat);

          // Compute the singular value decomposition of B
          RefSCMatrix U(dvc_,dvc_,matrixkit_);
          RefSCMatrix V(dnatom3_,dnatom3_,matrixkit_);
          RefSCDimension min;
          if (dnatom3_.n()<dvc_.n()) min = dnatom3_;
          else min = dvc_;
          int nmin = min.n();
          RefDiagSCMatrix sigma(min,matrixkit_);
          bmat.svd(U,sigma,V);

          // compute the epsilon rank of B
          int rank = 0;
          for (i=0; i<nmin; i++) {
              if (fabs(sigma(i)) > 0.0001) rank++;
            }

          RefSCDimension drank = new SCDimension(rank);
          RefDiagSCMatrix sigma_i(drank,matrixkit_);
          for (i=0; i<rank; i++) {
              sigma_i(i) = 1.0/sigma(i);
            }
          RefSCMatrix Ur(dvc_, drank, matrixkit_);
          RefSCMatrix Vr(dnatom3_, drank, matrixkit_);
          Ur.assign_subblock(U,0, dvc_.n()-1, 0, drank.n()-1, 0, 0);
          Vr.assign_subblock(V,0, dnatom3_.n()-1, 0, drank.n()-1, 0, 0);
          internal_to_cart_disp = Vr * sigma_i * Ur.t();

        }

      // compute the cartesian displacements
      RefSCVector cartesian_displacement = internal_to_cart_disp*displacement;
      if (debug_ && step == 0) {
          internal_to_cart_disp.print("Internal to Cartesian Transform");
          cartesian_displacement.print("Step 0 Cartesian Displacment");
        }
      // update the geometry
      for (int i=0; i < dnatom3_.n(); i++) {
#if OLD_BMAT
          molecule.r(i/3,i%3) += cartesian_displacement(i) * 1.88972666;
#else        
          molecule.r(i/3,i%3) += cartesian_displacement(i);
#endif          
        }

      // fix symmetry breaking due to numerical round-off
      molecule.cleanup_molecule();

      // check for convergence
      Ref<SCElementMaxAbs> maxabs = new SCElementMaxAbs();
      Ref<SCElementOp> op = maxabs.pointer();
      cartesian_displacement.element_op(op);
      maxabs_cart_diff = maxabs->result();
      if (maxabs_cart_diff < cartesian_tolerance_) {

          constant_->update_values(mol);
          variable_->update_values(mol);

          return 0;
        }
    }

  ExEnv::err0() << indent
       << "WARNING: IntMolecularCoor::all_to_cartesian(RefSCVector&):"
       << " too many iterations in geometry update" << endl;

  new_internal.print("desired internal coordinates");
  (new_internal
   - old_internal).print("difference of desired and actual coordinates");

  return -1;
}

int
IntMolecularCoor::to_cartesian(const Ref<Molecule> &mol,
                               const RefSCVector&new_internal)
{
  if (new_internal.dim().n() != dim_.n()
      || dvc_.n() != variable_->n() + constant_->n()
      || new_internal.dim().n() != variable_->n()) {
      throw ProgrammingError("to_cartesian: internal error in dim",
                             __FILE__, __LINE__, class_desc());
    }

  RefSCVector all_internal(dvc_,matrixkit_);

  int i,j;

  for (i=0; i<variable_->n(); i++) all_internal(i) = new_internal(i);
  for (j=0; j<constant_->n(); i++,j++) {
      all_internal(i) = constant_->coor(j)->value();
    }

  int ret = all_to_cartesian(mol, all_internal);

  if (watched_.nonnull()) {
      ExEnv::out0() << endl
           << indent << "Watched coordinate(s):\n" << incindent;
      watched_->update_values(mol);
      watched_->print_details(mol,ExEnv::out0());
      ExEnv::out0() << decindent;
    }
  
  return ret;
}

int
IntMolecularCoor::all_to_internal(const Ref<Molecule> &mol,RefSCVector&internal)
{
  if (internal.dim().n() != dvc_.n()
      || dim_.n() != variable_->n()
      || dvc_.n() != variable_->n() + constant_->n()) {
      throw ProgrammingError("all_to_internal: internal error in dim",
                             __FILE__, __LINE__, class_desc());
    }

  variable_->update_values(mol);
  constant_->update_values(mol);
   
  int n = dim_.n();
  int i;
  for (i=0; i<n; i++) {
      internal(i) = variable_->coor(i)->value();
    }
  n = dvc_.n();
  for (int j=0; i<n; i++,j++) {
      internal(i) = constant_->coor(j)->value();
    }

  return 0;
}

int
IntMolecularCoor::to_internal(RefSCVector&internal)
{
  if (internal.dim().n() != dim_.n()
      || dim_.n() != variable_->n()) {
      throw ProgrammingError("to_internal: internal error in dim",
                             __FILE__, __LINE__, class_desc());
    }

  variable_->update_values(molecule_);
   
  int n = dim_.n();
  for (int i=0; i<n; i++) {
      internal(i) = variable_->coor(i)->value();
    }

  return 0;
}

int
IntMolecularCoor::to_cartesian(RefSCVector&gradient,RefSCVector&internal)
{
  RefSCMatrix bmat(dim_,gradient.dim(),matrixkit_);
  variable_->bmat(molecule_,bmat);

  gradient = bmat.t() * internal;
  
  return 0;
}

// converts the gradient in cartesian coordinates to internal coordinates
int
IntMolecularCoor::to_internal(RefSCVector&internal,RefSCVector&gradient)
{
  RefSCMatrix bmat(dvc_,gradient.dim(),matrixkit_);
  RefSymmSCMatrix bmbt(dvc_,matrixkit_);

  Ref<SetIntCoor> variable_and_constant = new SetIntCoor();
  variable_and_constant->add(variable_);
  variable_and_constant->add(constant_);

  // form the bmatrix
  variable_and_constant->bmat(molecule_,bmat);
  
  // form the inverse of bmatrix * bmatrix_t
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(bmat);
  bmbt = bmbt.gi();

#if OLD_BMAT  
  RefSCVector all_internal = bmbt*bmat*(gradient*8.2388575);
#else
  RefSCVector all_internal = bmbt*bmat*gradient;
#endif  

  // put the variable coordinate gradients into internal
  int n = variable_->n();
  for (int i=0; i<n; i++) {
      internal.set_element(i,all_internal.get_element(i));
    }

  return 0;
}

int
IntMolecularCoor::to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal)
{
  cart.assign(0.0);
  RefSCMatrix bmat(dim_,cart.dim(),matrixkit_);
  variable_->bmat(molecule_,bmat);
  cart.accumulate_transform(bmat.t(),internal);
  return 0;
}

int
IntMolecularCoor::to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart)
{
  // form bmat
  RefSCMatrix bmat(dim_,cart.dim(),matrixkit_);
  variable_->bmat(molecule_,bmat);
  // and (B*B+)^-1
  RefSymmSCMatrix bmbt(dim_,matrixkit_);
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(bmat);
  bmbt = bmbt.gi();

  internal.assign(0.0);
  internal.accumulate_transform(bmbt*bmat,cart);
  return 0;
}

int
IntMolecularCoor::nconstrained()
{
  return fixed_->n();
}

void
IntMolecularCoor::print(ostream& os) const
{
  all_->update_values(molecule_);

  os << indent << "IntMolecularCoor Parameters:\n" << incindent
     << indent << "update_bmat = " << (update_bmat_?"yes":"no") << endl
     << indent << "scale_bonds = " << scale_bonds_ << endl
     << indent << "scale_bends = " << scale_bends_ << endl
     << indent << "scale_tors = " << scale_tors_ << endl
     << indent << "scale_outs = " << scale_outs_ << endl
     << indent << scprintf("symmetry_tolerance = %e\n",symmetry_tolerance_)
     << indent << scprintf("simple_tolerance = %e\n",simple_tolerance_)
     << indent << scprintf("coordinate_tolerance = %e\n",coordinate_tolerance_)
     << indent << "have_fixed_values = " << given_fixed_values_ << endl
     << indent << "max_update_steps = " << max_update_steps_ << endl
     << indent << scprintf("max_update_disp = %f\n",max_update_disp_)
     << indent << "have_fixed_values = " << given_fixed_values_ << endl
     << decindent << endl;

  molecule_->print(os);
  os << endl;

  print_simples(os);
  os << endl;

  if (form_print_variable_) {
      print_variable(os);
      os << endl;
    }

  if (form_print_constant_) {
      print_constant(os);
      os << endl;
    }
}

void
IntMolecularCoor::print_simples(ostream& os) const
{
  if (matrixkit()->messagegrp()->me()==0) {
    if (bonds_->n()) {
      os << indent << "Bonds:\n" << incindent;
      bonds_->print_details(molecule_,os);
      os << decindent;
    }
    if (bends_->n()) {
      os << indent << "Bends:\n" << incindent;
      bends_->print_details(molecule_,os);
      os << decindent;
    }
    if (tors_->n()) {
      os << indent << "Torsions:\n" << incindent;
      tors_->print_details(molecule_,os);
      os << decindent;
    }
    if (outs_->n()) {
      os << indent << "Out of Plane:\n" << incindent;
      outs_->print_details(molecule_,os);
      os << decindent;
    }
    if (extras_->n()) {
      os << indent << "Extras:\n" << incindent;
      extras_->print_details(molecule_,os);
      os << decindent;
    }
    if (fixed_->n()) {
      os << indent << "Fixed:\n" << incindent;
      fixed_->print_details(molecule_,os);
      os << decindent;
    }
    if (followed_.nonnull()) {
      os << indent << "Followed:\n" << incindent;
      followed_->print_details(molecule_,os);
      os << decindent;
    }
    if (watched_.nonnull()) {
      os << indent << "Watched:\n" << incindent;
      watched_->print_details(molecule_,os);
      os << decindent;
    }
  }
}

void
IntMolecularCoor::print_variable(ostream& os) const
{
  if (variable_->n() == 0) return;
  os << indent
     << "Variable Coordinates:" << endl;
  os << incindent;
  variable_->print_details(molecule_,os);
  os << decindent;
}

void
IntMolecularCoor::print_constant(ostream& os) const
{
  if (constant_->n() == 0) return;
  os << indent
     << "Constant Coordinates:" << endl;
  os << incindent;
  constant_->print_details(molecule_,os);
  os << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
