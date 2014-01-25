//
// symmcoor.cc
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

#define VERBOSE 0

namespace sc {

///////////////////////////////////////////////////////////////////////////
// SymmCoorTransform

class SymmCoorTransform: public NonlinearTransform {
  private:
    Ref<Molecule> molecule_;
    RefSCDimension dnatom3_;
    Ref<SetIntCoor> oldintcoor_;
    Ref<SetIntCoor> newintcoor_;
    Ref<SCMatrixKit> matrixkit_;
    int transform_hessian_;
  public:
    SymmCoorTransform(const Ref<Molecule>& molecule,
                      const RefSCDimension& dnatom3,
                      const Ref<SCMatrixKit>& kit,
                      const Ref<SetIntCoor>& oldintcoor,
                      const Ref<SetIntCoor>& newintcoor,
                      int transform_hessian);
    void to_cartesian(const RefSCVector& new_internal);
    void transform_coordinates(const RefSCVector& x);
    void transform_hessian(const RefSymmSCMatrix& h);
};

SymmCoorTransform::SymmCoorTransform(const Ref<Molecule>& molecule,
                                     const RefSCDimension& dnatom3,
                                     const Ref<SCMatrixKit>& kit,
                                     const Ref<SetIntCoor>& oldintcoor,
                                     const Ref<SetIntCoor>& newintcoor,
                                     int transform_hessian)
{
  molecule_ = new Molecule(*molecule.pointer());
  dnatom3_ = dnatom3;
  matrixkit_ = kit;
  oldintcoor_ = oldintcoor;
  newintcoor_ = newintcoor;
  transform_hessian_ = transform_hessian;
}

void
SymmCoorTransform::to_cartesian(const RefSCVector& new_internal)
{
  Ref<SCMatrixKit> kit = new_internal.kit();

  // get a reference to Molecule for convenience
  Molecule& molecule = *(molecule_.pointer());

  RefSCDimension dim = new_internal.dim();

  // don't bother updating the bmatrix when the error is less than this
  const double update_tolerance = 1.0e-3;
  const double cartesian_tolerance = 1.0e-8;

  // compute the internal coordinate displacements
  RefSCVector old_internal(new_internal.dim(),kit);

  RefSCMatrix internal_to_cart_disp;
  double maxabs_cart_diff = 0.0;
  const int maxiter = 100;
  for (int step = 0; step < maxiter; step++) {
      // compute the old internal coordinates
      oldintcoor_->update_values(molecule_);
      oldintcoor_->values_to_vector(old_internal);

      // the displacements
      RefSCVector displacement = new_internal - old_internal;

      if (maxabs_cart_diff>update_tolerance
          || internal_to_cart_disp.null()) {

          int i;
          RefSCMatrix bmat(dim,dnatom3_,kit);

          // form the bmatrix
          oldintcoor_->bmat(molecule_,bmat);

          // Compute the singular value decomposition of B
          RefSCMatrix U(dim,dim,kit);
          RefSCMatrix V(dnatom3_,dnatom3_,kit);
          RefSCDimension min;
          if (dnatom3_.n()<dim.n()) min = dnatom3_;
          else min = dim;
          int nmin = min.n();
          RefDiagSCMatrix sigma(min,kit);
          bmat.svd(U,sigma,V);

          // compute the epsilon rank of B
          int rank = 0;
          for (i=0; i<nmin; i++) {
              if (fabs(sigma(i)) > 0.0001) rank++;
            }

          RefSCDimension drank = new SCDimension(rank);
          RefDiagSCMatrix sigma_i(drank,kit);
          for (i=0; i<rank; i++) {
              sigma_i(i) = 1.0/sigma(i);
            }
          RefSCMatrix Ur(dim, drank, kit);
          RefSCMatrix Vr(dnatom3_, drank, kit);
          Ur.assign_subblock(U,0, dim.n()-1, 0, drank.n()-1, 0, 0);
          Vr.assign_subblock(V,0, dnatom3_.n()-1, 0, drank.n()-1, 0, 0);
          internal_to_cart_disp = Vr * sigma_i * Ur.t();
        }

      // compute the cartesian displacements
      RefSCVector cartesian_displacement = internal_to_cart_disp*displacement;
      // update the geometry
      for (int i=0; i < dnatom3_.n(); i++) {
          molecule.r(i/3,i%3) += cartesian_displacement(i);
        }

      // check for convergence
      Ref<SCElementMaxAbs> maxabs = new SCElementMaxAbs();
      Ref<SCElementOp> op = maxabs.pointer();
      cartesian_displacement.element_op(op);
      maxabs_cart_diff = maxabs->result();
      if (maxabs_cart_diff < cartesian_tolerance) {
          oldintcoor_->update_values(molecule_);
          return;
        }
    }

  throw MaxIterExceeded("too many iterations in geometry update",
                        __FILE__, __LINE__, maxiter);
}

void
SymmCoorTransform::transform_coordinates(const RefSCVector& x)
{
  if (x.null()) return;

  Ref<SCMatrixKit> kit = x.kit();
  RefSCDimension dim = x.dim();

  // using the old coordinates update molecule
  to_cartesian(x);

  // compute the new coordinates
  newintcoor_->update_values(molecule_);
  newintcoor_->values_to_vector(x);

  // compute the linear transformation information

  // the old B matrix
  RefSCMatrix B(dim, dnatom3_, kit);
  oldintcoor_->bmat(molecule_, B);

  // get the B matrix for the new coordinates
  RefSCMatrix Bnew(dim, dnatom3_, kit);
  newintcoor_->update_values(molecule_);
  newintcoor_->bmat(molecule_, Bnew);

  // the transform from cartesian to new internal coordinates
  RefSymmSCMatrix bmbt(dim,kit);
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(Bnew);
  RefSCMatrix cart_to_new_internal = bmbt.gi() * Bnew;
  Bnew = 0;
  bmbt = 0;

  linear_transform_ = cart_to_new_internal * B.t();
#if VERBOSE
  linear_transform_.print("old internal to new");
#endif
}

void
SymmCoorTransform::transform_hessian(const RefSymmSCMatrix& h)
{
  if (transform_hessian_) {
      NonlinearTransform::transform_hessian(h);
    }
  else {
      ExEnv::err0() << indent
           << "WARNING: SymmCoorTransform::transform_hessian: "
           << "skipping hessian transform";
    }
}

///////////////////////////////////////////////////////////////////////////
// members of SymmMolecularCoor

static ClassDesc SymmMolecularCoor_cd(
  typeid(SymmMolecularCoor),"SymmMolecularCoor",1,"public IntMolecularCoor",
  0, create<SymmMolecularCoor>, create<SymmMolecularCoor>);

SymmMolecularCoor::SymmMolecularCoor(Ref<Molecule>&mol):
  IntMolecularCoor(mol)
{
  init();
}

SymmMolecularCoor::SymmMolecularCoor(const Ref<KeyVal>& keyval):
  IntMolecularCoor(keyval)
{
  init();

  int itmp;
  double dtmp;
  itmp = keyval->booleanvalue("change_coordinates");
  if (keyval->error() == KeyVal::OK) change_coordinates_ = itmp;
  itmp = keyval->booleanvalue("transform_hessian");
  if (keyval->error() == KeyVal::OK) transform_hessian_ = itmp;
  dtmp = keyval->doublevalue("max_kappa2");
  if (keyval->error() == KeyVal::OK) max_kappa2_ = dtmp;
}

SymmMolecularCoor::SymmMolecularCoor(StateIn& s):
  IntMolecularCoor(s)
{
  s.get(change_coordinates_);
  s.get(transform_hessian_);
  s.get(max_kappa2_);
}

SymmMolecularCoor::~SymmMolecularCoor()
{
}

void
SymmMolecularCoor::save_data_state(StateOut&s)
{
  IntMolecularCoor::save_data_state(s);
  s.put(change_coordinates_);
  s.put(transform_hessian_);
  s.put(max_kappa2_);
}

void
SymmMolecularCoor::init()
{
  IntMolecularCoor::init();
  change_coordinates_ = 0;
  max_kappa2_ = 10.0;
  transform_hessian_ = 1;
}

void
SymmMolecularCoor::form_coordinates(int keep_variable)
{
  int i;
  int nbonds = bonds_->n();
  int nbends = bends_->n();
  int ntors = tors_->n();
  int nouts = outs_->n();
  int nextras = extras_->n();

  Ref<SetIntCoor> saved_fixed_ = fixed_;
  fixed_ = new SetIntCoor;
  fixed_->add(saved_fixed_);
  // if we're following coordinates, add them to the fixed list
  if (followed_.nonnull())
    fixed_->add(followed_);
  
  int nredundant = nbonds + nbends + ntors + nouts + nextras;
  int nfixed = fixed_->n();

  // see how many coords we expect
  int n3 = molecule_->natom()*3;
  int nunique = n3 - 6; // need to detect linear

  if (nredundant < nunique) {
      AlgorithmException ex("found too few redundant coordinates",
                            __FILE__, __LINE__, class_desc());
      try {
          ex.elaborate()
              << scprintf("nredundant = %d, 3n-6 = %d",
                          nredundant, nunique)
              << std::endl
              << "(the geometry is probably bad)"
              << std::endl;
        }
      catch (...) {}
      throw ex;
    }

  RefSCDimension dredundant = new SCDimension(nredundant, "Nredund");
  RefSCDimension dfixed = new SCDimension(nfixed, "Nfixed");
  RefSCMatrix K; // nredundant x nnonzero
  int* is_totally_symmetric; // nnonzero; if 1 coor has tot. symm. component

  form_K_matrix(dredundant, dfixed, K, is_totally_symmetric);

  RefSCDimension dnonzero = K.coldim();
  int nnonzero = dnonzero.n();

  if (!keep_variable) variable_->clear();
  constant_->clear();

  // now remove followed coords from the fixed list, and add to the
  // variable list
  if (followed_.nonnull()) {
    fixed_->pop();
    variable_->add(followed_);
  }
  
  // put the fixed coordinates into the constant list
  nfixed = fixed_->n();
  for (i=0; i<nfixed; i++) {
      constant_->add(fixed_->coor(i));
    }

  // ok, now we have the K matrix, the columns of which give us the
  // contribution from each red. coord to the ith non-red. coord.
  // this gets a little hairy since the red coords can themselves be
  // linear combinations of simple coords
  for(i=0; i < nnonzero; i++) {
      // construct the new linear combination coordinate
      char label[80];
      if (is_totally_symmetric[i]) {
          sprintf(label,"symm_coord_%03d",i+1);
        }
      else {
          sprintf(label,"asymm_coord_%03d",i+1);
        }
      SumIntCoor* coordinate = new SumIntCoor(label);

      int j;
      for(j=0; j < nredundant; j++) {
          if(pow(K(j,i),2.0) > simple_tolerance_) {
              Ref<IntCoor> c = all_->coor(j);
              coordinate->add(c,K(j,i));
              if (debug_) {
                  ExEnv::out0() << "added redund coor "
                       << j << " to coor " << i << ":" << endl;
                  c->print();
                }
            }
        }

      // normalize the coordinate
      coordinate->normalize();

      if (only_totally_symmetric_ && !is_totally_symmetric[i]) {
          // Don't put nonsymmetric coordinates into the
          // constant_ coordinate set.  This causes problems
          // when coordinates with small coefficients are eliminated
          // since they can then acquire symmetric components.
          // constant_->add(coordinate);
          delete coordinate;
        }
      else {
          if (!keep_variable) variable_->add(coordinate);
        }
    }

  constant_->update_values(molecule_);
  variable_->update_values(molecule_);

  ExEnv::out0() << incindent << indent
       << "SymmMolecularCoor::form_variable_coordinates()\n" << incindent
       << indent << "expected " << nunique << " coordinates\n"
       << indent << "found " << variable_->n() << " variable coordinates\n"
       << indent << "found " << constant_->n() << " constant coordinates\n"
       << decindent << decindent << flush;

  delete[] is_totally_symmetric;
  fixed_ = saved_fixed_;

  if (form_print_molecule_) molecule_->print();
  if (form_print_simples_) print_simples(ExEnv::out0());
  if (form_print_variable_) print_variable(ExEnv::out0());
  if (form_print_constant_) print_constant(ExEnv::out0());
}

void
SymmMolecularCoor::guess_hessian(RefSymmSCMatrix&hessian)
{
  // first form diagonal hessian in redundant internal coordinates
  RefSCDimension rdim = new SCDimension(all_->n(), "Nall");
  RefSymmSCMatrix rhessian(rdim,matrixkit_);
  rhessian.assign(0.0);
  all_->guess_hessian(molecule_,rhessian);

  // create redundant coordinate bmat
  RefSCDimension dn3 = dnatom3_;
  RefSCMatrix bmatr(rdim,dn3,matrixkit_);
  all_->bmat(molecule_,bmatr);

  // then form the variable coordinate bmat
  RefSCDimension dredundant = new SCDimension(variable_->n(), "Nvar");
  RefSCMatrix bmat(dredundant,dn3,matrixkit_);
  variable_->bmat(molecule_,bmat);

  // and (B*B+)^-1
  RefSymmSCMatrix bmbt(dredundant,matrixkit_);
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(bmat);
  bmbt = bmbt.gi();

  // now transform redundant hessian to internal coordinates
  // Hc = Br+ * Hr * Br
  // Hi = (B*B+)^-1 * B * Hc * B+ * (B*B+)^-1+
  //    = bmbt_inv*B*Br+ * Hr * Br*B+*bmbt_inv+
  //    = b * Hr * b+  (b = (B*B+)^-1 * B * Br+)
  RefSCMatrix b = bmbt * bmat * bmatr.t();
  
  hessian.assign(0.0);
  hessian.accumulate_transform(b,rhessian);
}

RefSymmSCMatrix
SymmMolecularCoor::inverse_hessian(RefSymmSCMatrix& hessian)
{
  return hessian.gi();
}

// Possibly change to a new coordinate system
Ref<NonlinearTransform>
SymmMolecularCoor::change_coordinates()
{
  if (dim_.n() == 0 || !change_coordinates_) return 0;

  const double epsilon = 0.001;

  // compute the condition number of the old coordinate system at the
  // current point
  RefSCMatrix B(dim_, dnatom3_, matrixkit_);
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
  // the rank could get bigger if there is a fixed coordinate
  if (rank < dim_.n() || ((fixed_.null()
                           || fixed_->n() == 0) && rank != dim_.n())) {
      throw AlgorithmException("disallowed rank change",
                               __FILE__, __LINE__, class_desc());
    }
  if (rank != dim_.n()) {
      ExEnv::out0() << indent
           << "SymmMolecularCoor::change_coordinates: rank changed\n";
    }

  double kappa2 = sigma(0)/sigma(dim_.n()-1);

  ExEnv::out0() << indent
       << scprintf(
           "SymmMolecularCoor: condition number = %14.8f (max = %14.8f)\n",
           kappa2, max_kappa2_);

  if (kappa2 > max_kappa2_) {
      Ref<SetIntCoor> oldvariable = new SetIntCoor;
      oldvariable->add(variable_);

      // form the new variable coordinates
      form_coordinates();

      SymmCoorTransform *trans = new SymmCoorTransform(molecule_,
                                                       dnatom3_,
                                                       matrixkit_,
                                                       oldvariable,
                                                       variable_,
                                                       transform_hessian_);
      return trans;
    }

  return 0;
}

void
SymmMolecularCoor::print(ostream& os) const
{
  IntMolecularCoor::print(os);
  
  os << indent << "SymmMolecularCoor Parameters:\n" << incindent
     << indent << "change_coordinates = "
     << (change_coordinates_?"yes":"no") << endl
     << indent << "transform_hessian = "
     << (transform_hessian_?"yes":"no") << endl
     << indent << scprintf("max_kappa2 = %f",max_kappa2_) << endl
     << decindent << endl;
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
