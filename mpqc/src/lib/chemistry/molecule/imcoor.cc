
extern "C" {
#include <math.h>
};

#include <util/misc/formio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <chemistry/molecule/localdef.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/simple.h>

#include <math/topology/bitarray.h>

#define DEFAULT_SIMPLE_TOLERANCE 1.0e-3

#define USE_SVD 1
#define VERBOSE 0

///////////////////////////////////////////////////////////////////////////
// members of IntMolecularCoor

#define CLASSNAME IntMolecularCoor
#define VERSION 3
#define PARENTS public MolecularCoor
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
IntMolecularCoor::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MolecularCoor::_castdown(cd);
  return do_castdowns(casts,cd);
}

IntMolecularCoor::IntMolecularCoor(RefMolecule&mol):
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
  max_update_steps_(100),
  max_update_disp_(0.5),
  given_fixed_values_(0),
  decouple_bonds_(0),
  decouple_bends_(0)
{
  new_coords();
  generator_ = new IntCoorGen(mol);
}

IntMolecularCoor::IntMolecularCoor(const RefKeyVal& keyval):
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
  generator_.restore_state(s);

  if (s.version(static_class_desc()) >= 3) {
      s.get(decouple_bonds_);
      s.get(decouple_bends_);
    }
  else {
      decouple_bonds_ = 0;
      decouple_bends_ = 0;
    }
  
  if (s.version(static_class_desc()) >= 2) {
    //printf("IntMolecularCoor::IntMolecularCoor(StateIn& s): new version\n");
    s.get(max_update_steps_);
    s.get(max_update_disp_);
    s.get(given_fixed_values_);
  } else {
    //printf("IntMolecularCoor::IntMolecularCoor(StateIn& s): old version\n");
    max_update_steps_ = 100;
    max_update_disp_ = 0.5;
    given_fixed_values_ = 0;
  }

  dim_.restore_state(s);
  dvc_.restore_state(s);

  all_.restore_state(s);

  variable_.restore_state(s);
  constant_.restore_state(s);

  fixed_.restore_state(s);
  followed_.restore_state(s);

  bonds_.restore_state(s);
  bends_.restore_state(s);
  tors_.restore_state(s);
  outs_.restore_state(s);
  extras_.restore_state(s);

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
  constant_ = new SetIntCoor; // internal coors to be head fixed
  bonds_ = new SetIntCoor;
  bends_ = new SetIntCoor;
  tors_ = new SetIntCoor;
  outs_ = new SetIntCoor;
  extras_ = new SetIntCoor;
  fixed_ = new SetIntCoor;
  followed_ = 0;
}

void
IntMolecularCoor::read_keyval(const RefKeyVal& keyval)
{
  fixed_ = keyval->describedclassvalue("fixed");
  if (fixed_.null()) fixed_ = new SetIntCoor;
  followed_ = keyval->describedclassvalue("followed");

  decouple_bonds_ = keyval->booleanvalue("decouple_bonds");
  decouple_bends_ = keyval->booleanvalue("decouple_bends");

  given_fixed_values_ = keyval->booleanvalue("have_fixed_values");

  max_update_steps_ = keyval->intvalue("max_update_steps");
  if (keyval->error() != KeyVal::OK) max_update_steps_ = 100;

  max_update_disp_ = keyval->doublevalue("max_update_disp");
  if (keyval->error() != KeyVal::OK) max_update_disp_ = 0.5;

  generator_ = keyval->describedclassvalue("generator");

  if (generator_.null()) {
      // the extra_bonds list is given as a vector of atom numbers
      // (atom numbering starts at 1)
      int nextra_bonds = keyval->count("extra_bonds");
      nextra_bonds /= 2;
      int *extra_bonds;
      if (nextra_bonds) {
          extra_bonds = new int[nextra_bonds*2];
          for (int i=0; i<nextra_bonds*2; i++) {
              extra_bonds[i] = keyval->intvalue("extra_bonds",i);
              if (keyval->error() != KeyVal::OK) {
                  fprintf(stderr,"  IntMolecularCoor:: keyval CTOR: "
                          "problem reading \"extra_bonds:%d\"\n",i);
                  abort();
                }
            }
        }
      else {
          extra_bonds = 0;
        }
      generator_ = new IntCoorGen(molecule_, nextra_bonds, extra_bonds);
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
}

void
IntMolecularCoor::init()
{
  RefSetIntCoor redundant = new SetIntCoor;
  generator_->generate(redundant);

  // sort out the simple coordinates by type
  int i;
  for (i=0; i<redundant->n(); i++) {
      RefIntCoor coor = redundant->coor(i);
      if (coor->class_desc()
          == StreSimpleCo::static_class_desc()) {
          bonds_->add(coor);
        }
      else if (coor->class_desc()
               == BendSimpleCo::static_class_desc()) {
          bends_->add(coor);
        }
      else if (coor->class_desc()
               == TorsSimpleCo::static_class_desc()
               || coor->class_desc()
               == ScaledTorsSimpleCo::static_class_desc()) {
          tors_->add(coor);
        }
      else if (coor->class_desc()
               == OutSimpleCo::static_class_desc()) {
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

  if (given_fixed_values_) {
      // save the given coordinate values
      RefSCDimension original_dfixed
          = new SCDimension(fixed_->n(),"Nfix");
      RefSCVector given_fixed_coords(original_dfixed,matrixkit_);
      for (i=0; i<original_dfixed.n(); i++) {
          given_fixed_coords(i) = fixed_->coor(i)->value();
        }

      // break up the displacement into several manageable steps
      double maxabs = given_fixed_coords.maxabs();
      int nstep = int(maxabs/max_update_disp_) + 1;
      given_fixed_coords.scale(1.0/nstep);
      printf("  IntMolecularCoor: displacing fixed coordinates to the"
             " requested values in %d steps\n", nstep);
      for (int istep=1; istep<=nstep; istep++) {
          form_coordinates();

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
                  = istep * double(given_fixed_coords(j));
            }
          for (; j<constant_->n(); i++,j++) {
              new_internal_coordinates(i) = constant_->coor(j)->value();
            }

          all_to_cartesian(new_internal_coordinates);
        }
    }

  form_coordinates();

  dim_ = new SCDimension(variable_->n(), "Nvar");
  dvc_ = new SCDimension(variable_->n()+constant_->n(),
                         "Nvar+Nconst");
#if 1
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
          printf("  IntMolecularCoor::init: rank changed\n");
        }

      double kappa2 = sigma(0)/sigma(dim_.n()-1);

      printf("  IntMolecularCoor::init: condition number = %14.8f\n", kappa2);
    }
#endif
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

#if USE_SVD
static RefSymmSCMatrix
form_partial_K(const RefSetIntCoor& coor, RefMolecule& molecule,
               const RefSCVector& geom,
               double epsilon,
               const RefSCDimension& dnatom3,
               const RefSCMatrixKit& matrixkit,
               RefSCMatrix& projection,
               RefSCVector& totally_symmetric,
               RefSCMatrix& K)
{
  // Compute the B matrix for the coordinates
  RefSCDimension dcoor = new SCDimension(coor->n());
  RefSCMatrix B(dcoor, dnatom3,matrixkit);
  coor->bmat(molecule, B);

  // Project out the previously discovered internal coordinates
  if (projection.nonnull()) {
      B = B * projection;
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

  RefSCMatrix SIGMA(dcoor, dnatom3,matrixkit);
  SIGMA.assign(0.0);
  for (i=0; i<nmin; i++) {
      SIGMA(i,i) = sigma(i);
    }

  // return if there are no new coordinates
  if (rank==0) return 0;

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

      int ntotally_symmetric = count_nonzero(totally_symmetric,0.001);
      printf("  found %d totally symmetric coordinates\n", ntotally_symmetric);

      // compute the cumulative projection
      if (projection.null()) {
          projection = matrixkit->matrix(dnatom3,dnatom3);
          projection->unit();
        }
      projection = projection * proj_nullspace_B;
    }

  // give Ur to caller
  K = Ur;

  return proj_nullspace_B_perp;
}
#endif // USE_SVD

// this allocates storage for and computes K, Kfixed, and is_totally_symmetric
#if USE_SVD
void
IntMolecularCoor::form_K_matrices(RefSCDimension& dredundant,
                                  RefSCDimension& dfixed,
                                  RefSCMatrix& K,
                                  RefSCMatrix& Kfixed,
                                  int*& is_totally_symmetric)
{
  int i,j;

  // The cutoff for whether or not a coordinate is considered totally symmetric
  double ts_eps = 0.0001;

  // The geometry will be needed to check for totally symmetric
  // coordinates
  RefSCVector geom(dnatom3_,matrixkit_);
  for(i=0; i < geom.n()/3; i++) {
      geom(3*i) = molecule_->operator[](i).point()[0];
      geom(3*i+1) = molecule_->operator[](i).point()[1];
      geom(3*i+2) = molecule_->operator[](i).point()[2];
    }

  RefSCDimension dcoor = new SCDimension(all_->n());

  // this keeps track of the total projection for the b matrices
  RefSCMatrix projection;
  if (dfixed.n()) {
      RefSCMatrix Ktmp;
      RefSCVector totally_symmetric_fixed;
      RefSymmSCMatrix null_bfixed_perp
          = form_partial_K(fixed_, molecule_, geom, 0.001, dnatom3_,
                           matrixkit_, projection, totally_symmetric_fixed,
                           Ktmp);
      // require that the epsilon rank equal the number of fixed coordinates
      if (Ktmp.nrow() != dfixed.n()) {
          fprintf(stderr, "  ERROR: IntMolecularCoor: nfixed = %d rank = %d\n",
                  dfixed.n(), Ktmp.ncol());
          abort();
        }
      // check that fixed coordinates be totally symmetric
      if (Ktmp.nrow() != count_nonzero(totally_symmetric_fixed, ts_eps)) {
          fprintf(stderr, "  WARNING: only %d of %d fixed coordinates are"
                  " totally symmetric\n",
                  count_nonzero(totally_symmetric_fixed, ts_eps), dfixed.n());
        }

      // Compute Kfixed
      RefSCMatrix B(dcoor, dnatom3_, matrixkit_);
      all_->bmat(molecule_, B);
      Kfixed = B * null_bfixed_perp;
    }

#define DECOUPLE 1
#if DECOUPLE
  int n_total = 0;

  RefSCVector totally_symmetric_bond;
  RefSCMatrix Kbond;
  if (decouple_bonds_) {
      cout << "looking for bonds" << endl;
      form_partial_K(bonds_, molecule_, geom, 0.1, dnatom3_, matrixkit_,
                     projection, totally_symmetric_bond, Kbond);
      if (Kbond.nonnull()) n_total += Kbond.ncol();
    }

  RefSCVector totally_symmetric_bend;
  RefSCMatrix Kbend;
  if (decouple_bends_) {
      cout << "looking for bends" << endl;
      form_partial_K(bends_, molecule_, geom, 0.1, dnatom3_, matrixkit_,
                     projection, totally_symmetric_bend, Kbend);
      if (Kbend.nonnull()) n_total += Kbend.ncol();
    }

  if (decouple_bonds_ || decouple_bends_) {
      cout << "looking for remaining coordinates" << endl;
    }
  RefSCVector totally_symmetric_all;
  RefSCMatrix Kall;
  // I hope the IntCoorSet keeps the ordering
  form_partial_K(all_, molecule_, geom, 0.001, dnatom3_, matrixkit_,
                 projection, totally_symmetric_all, Kall);
  int n_totally_symmetric_all = count_nonzero(totally_symmetric_all, ts_eps);
  if (Kall.nonnull()) n_total += Kall.ncol();

  // This requires that all_ coordinates is made up of first bonds,
  // bends, and finally the rest of the coordinates.
  RefSCDimension dtot = new SCDimension(n_total);
  K = matrixkit_->matrix(dcoor, dtot);
  K.assign(0.0);
  if (Kbond.nonnull()) {
      K.assign_subblock(Kbond, 0, Kbond.nrow()-1, 0, Kbond.ncol()-1, 0, 0);
    }
  if (Kbend.nonnull()) {
      K.assign_subblock(Kbend, Kbond.nrow(), Kbond.nrow()+Kbend.nrow()-1,
                        Kbond.ncol(), Kbond.ncol()+Kbend.ncol()-1, 0, 0);
    }
  if (Kall.nonnull()) {
      K.assign_subblock(Kall, 0, dcoor->n()-1,
                        Kbond.ncol()+Kbend.ncol(),
                        Kbond.ncol()+Kbend.ncol()+Kall.ncol()-1, 0, 0);
    }

  is_totally_symmetric = new int[K.ncol()];
  j=0;
  if (Kbond.nonnull()) {
      for (i=0; i<Kbond.ncol(); i++,j++) {
          if (fabs(totally_symmetric_bond(i)) > ts_eps)
              is_totally_symmetric[j] = 1;
          else is_totally_symmetric[j] = 0;
        }
    }
  if (Kbond.nonnull()) {
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
#else
  RefSCVector totally_symmetric_all;
  form_partial_K(all_, molecule_, geom, 0.001, dnatom3_, matrixkit_,
                 projection, totally_symmetric_all, K);
  int n_totally_symmetric_all = count_nonzero(totally_symmetric_all, ts_eps);

  is_totally_symmetric = new int[K.ncol()];
  for (i=0; i<K.ncol(); i++) {
      if (fabs(totally_symmetric_all(i)) > ts_eps) is_totally_symmetric[i] = 1;
      else is_totally_symmetric[i] = 0;
    }
#endif
}
#else // USE_SVD
void
IntMolecularCoor::form_K_matrices(RefSCDimension& dredundant,
                                  RefSCDimension& dfixed,
                                  RefSCMatrix& K,
                                  RefSCMatrix& Kfixed,
                                  int*& is_totally_symmetric)
{
  int i,j;
  int natom3 = dnatom3_.n();

  // form bmat for the set of redundant coordinates
  RefSCMatrix bmat(dredundant,dnatom3_,matrixkit_);
  all_->bmat(molecule_,bmat);
  int nredundant = dredundant.n();

  // scale the coordinates in the bmatrix
  i=0;
  int nbonds = bonds_->n();
  int nbends = bends_->n();
  int ntors = tors_->n();
  int nouts = outs_->n();
  // bonds
  if (scale_bonds_ != 1.0) {
      for (j=0; j<nbonds; i++,j++) {
          for (int k=0; k<natom3; k++) {
              bmat(i,k) = bmat(i,k) * scale_bonds_;
            }
        }
    }
  else {
      i += nbonds;
    }
  // bends
  if (scale_bends_ != 1.0) {
      for (j=0; j<nbends; i++,j++) {
          for (int k=0; k<natom3; k++) {
              bmat(i,k) = bmat(i,k) * scale_bends_;
            }
        }
    }
  else {
      i += nbends;
    }
  // torsions
  if (scale_tors_ != 1.0) {
      for (j=0; j<ntors; i++,j++) {
          for (int k=0; k<natom3; k++) {
              bmat(i,k) = bmat(i,k) * scale_tors_;
            }
        }
    }
  else {
      i += ntors;
    }
  // out of plane
  if (scale_outs_ != 1.0) {
      for (j=0; j<nouts; i++,j++) {
          for (int k=0; k<natom3; k++) {
              bmat(i,k) = bmat(i,k) * scale_outs_;
            }
        }
    }
  else {
      i += nouts;
    }

  // and form b*b~
  RefSymmSCMatrix bmbt(dredundant,matrixkit_);
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(bmat);

  // form the fixed part of the b matrix
  RefSCMatrix bmat_fixed(dfixed,dnatom3_,matrixkit_);
  fixed_->bmat(molecule_,bmat_fixed);
  int nfixed = dfixed.n();

  // and form the fixed b*b~
  RefSymmSCMatrix bmbt_fixed(dfixed,matrixkit_);
  bmbt_fixed.assign(0.0);
  bmbt_fixed.accumulate_symmetric_product(bmat_fixed);

  // need the cross terms in bmbt also
  RefSCMatrix bmbt_fix_red;
  if (nfixed != 0) bmbt_fix_red = bmat_fixed * bmat.t();

  // orthogonalize the redundant coordinates to the fixed coordinates
  RefSCMatrix redundant_ortho(dfixed,dredundant,matrixkit_);
  for (i=0; i<nredundant; i++) {
      for (j=0; j<nfixed; j++) {
          redundant_ortho(j,i) = - (bmbt_fix_red(j,i)/bmbt_fixed(j,j));
        }
    }

  // convert bmbt to the new coordinate system
  if (nfixed != 0) {
      bmbt.accumulate_transform(redundant_ortho.t(), bmbt_fixed);
      bmbt.accumulate_symmetric_sum(redundant_ortho.t() * bmbt_fix_red);
//       bmbt = bmbt
//         + redundant_ortho.t() * bmbt_fixed * redundant_ortho
//         + redundant_ortho.t() * bmbt_fix_red
//         + bmbt_fix_red.t() * redundant_ortho;
    }

  // now diagonalize bmbt, this should give you the 3n-6(5) symmetrized
  // internal coordinates
  RefSCMatrix vecs(dredundant,dredundant,matrixkit_);
  RefDiagSCMatrix vals(dredundant,matrixkit_);

  bmbt.diagonalize(vals,vecs);

  // ok, hopefully multiplying bmat*cart_coords will tell me which
  // coordinates are have totally symmetric components
  RefSCMatrix vecst = vecs.t();
  RefSCMatrix bm = vecst * bmat;
  // i'm done with bmat, get rid of it
  bmat = 0;
  RefSCVector geom(dnatom3_);
  for(i=0; i < geom.n()/3; i++) {
      geom(3*i) = molecule_->operator[](i).point()[0];
      geom(3*i+1) = molecule_->operator[](i).point()[1];
      geom(3*i+2) = molecule_->operator[](i).point()[2];
    }
  RefSCVector coords = bm * geom;

  // release unneeded memory
  geom = 0;
  bm = 0;

  // count the number of coordinates that have nonzero eigenvalues
  int nonzero=0;
  for(i=0; i < nredundant; i++) {
      if(fabs(vals(i)) > coordinate_tolerance_) nonzero++;
    }

  // size K and Kfixed
  RefSCDimension dnonzero = new SCDimension(nonzero, "Nnonzero");
  K = dredundant->create_matrix(dnonzero); // nredundant x nonzero
  K.assign(0.0);
  Kfixed = dfixed->create_matrix(dnonzero); // nfixed x nonzero
  Kfixed.assign(0.0);
  is_totally_symmetric = new int[nonzero];

  // generate K
  int coordno=0;
  for(i=0; i < nredundant; i++) {
      // nonzero eigenvalues are the non-redundant coordinates
      if(fabs(vals(i)) > coordinate_tolerance_) {

          int nonzero=0;
          for(j=0; j < nredundant; j++) {
              if(pow(vecs(j,i),2.0) > simple_tolerance_) nonzero++;
            }
          if(!nonzero) {
              fprintf(stderr,"  Geom_form_K: no nonzero simple coordinates");
              abort();
            }

          // bmat*cart_coords tells if a coordinate has a tot. symm. comp.
          if (fabs(coords(i)) > symmetry_tolerance_)
              is_totally_symmetric[coordno] = 1;
          else is_totally_symmetric[coordno] = 0;

          // construct the K arrays
          int ii;
          for(ii=0; ii < nredundant; ii++)
            K(ii,coordno) = vecs(ii,i);
          for(ii=0; ii < nfixed; ii++)
            Kfixed(ii,coordno) = redundant_ortho(ii,i);
          coordno++;
        }
    }
}
#endif // USE_SVD

IntMolecularCoor::~IntMolecularCoor()
{
}

void
IntMolecularCoor::save_data_state(StateOut&s)
{
  MolecularCoor::save_data_state(s);

  generator_.save_state(s);

  s.put(decouple_bonds_);
  s.put(decouple_bends_);

  s.put(max_update_steps_);
  s.put(max_update_disp_);
  s.put(given_fixed_values_);

  dim_.save_state(s);
  dvc_.save_state(s);

  all_.save_state(s);
  
  variable_.save_state(s);
  constant_.save_state(s);

  fixed_.save_state(s);
  followed_.save_state(s);

  bonds_.save_state(s);
  bends_.save_state(s);
  tors_.save_state(s);
  outs_.save_state(s);
  extras_.save_state(s);

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
IntMolecularCoor::all_to_cartesian(RefSCVector&new_internal)
{
  // get a reference to Molecule for convenience
  Molecule& molecule = *(molecule_.pointer());

  // don't bother updating the bmatrix when the error is less than this
  const double update_tolerance = 1.0e-6;

  // compute the internal coordinate displacements
  RefSCVector old_internal(dvc_,matrixkit_);

  RefSCMatrix internal_to_cart_disp;
  double maxabs_cart_diff = 0.0;
  for (int step = 0; step < max_update_steps_; step++) {
      // compute the old internal coordinates
      all_to_internal(old_internal);

#if VERBOSE
      cout << "Coordinates on step " << step << ":" << endl;
      variable_->print();
#endif

      // the displacements
      RefSCVector displacement = new_internal - old_internal;

      if ((update_bmat_ && maxabs_cart_diff>update_tolerance)
          || internal_to_cart_disp.null()) {
#if VERBOSE
          cout << "updating bmatrix" << endl;
#endif

          int i;
          RefSCMatrix bmat(dvc_,dnatom3_,matrixkit_);

          // form the set of all coordinates
          RefSetIntCoor variable_and_constant = new SetIntCoor();
          variable_and_constant->add(variable_);
          variable_and_constant->add(constant_);

          // form the bmatrix
          variable_and_constant->bmat(molecule_,bmat);

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

#if !USE_SVD
          RefSymmSCMatrix bmbt(dvc_,matrixkit_);
          RefSymmSCMatrix bmbt_i;

          // form the initial inverse of bmatrix * bmatrix_t
          bmbt.assign(0.0);
          bmbt.accumulate_symmetric_product(bmat);
          bmbt_i = bmbt.gi();

          internal_to_cart_disp = bmat.t() * bmbt_i;
#endif // !USE_SVD
        }

      // compute the cartesian displacements
      RefSCVector cartesian_displacement = internal_to_cart_disp*displacement;
      // update the geometry
      for (int i=0; i < dnatom3_.n(); i++) {
#if OLD_BMAT
          molecule[i/3][i%3] += cartesian_displacement(i) * 1.88972666;
#else        
          molecule[i/3][i%3] += cartesian_displacement(i);
#endif          
        }

      // check for convergence
      RefSCElementMaxAbs maxabs = new SCElementMaxAbs();
      RefSCElementOp op = maxabs;
      cartesian_displacement.element_op(op);
      maxabs_cart_diff = maxabs->result();
      if (maxabs_cart_diff < cartesian_tolerance_) {

          constant_->update_values(molecule_);
          variable_->update_values(molecule_);

          return 0;
        }
    }

  cout.flush();
  cerr.flush();
  fflush(stdout);
  fflush(stderr);

  fprintf(stderr,"  WARNING: IntMolecularCoor"
          "::all_to_cartesian(RefSCVector&):"
          " too many iterations in geometry update\n");

  fflush(stderr);

  new_internal.print("desired internal coordinates");
  (new_internal
   - old_internal).print("difference of desired and actual coordinates");

  cout.flush();

  return -1;
}

int
IntMolecularCoor::to_cartesian(const RefSCVector&new_internal)
{
  if (new_internal.dim().n() != dim_.n()
      || dvc_.n() != variable_->n() + constant_->n()
      || new_internal.dim().n() != variable_->n()) {
      fprintf(stderr,"  to_cartesian: internal error in dim\n");
      abort();
    }

  RefSCVector all_internal(dvc_,matrixkit_);

  int i,j;

  for (i=0; i<variable_->n(); i++) all_internal(i) = new_internal(i);
  for (j=0; j<constant_->n(); i++,j++) {
      all_internal(i) = constant_->coor(j)->value();
    }

  return all_to_cartesian(all_internal);
}

int
IntMolecularCoor::all_to_internal(RefSCVector&internal)
{
  if (internal.dim().n() != dvc_.n()
      || dim_.n() != variable_->n()
      || dvc_.n() != variable_->n() + constant_->n()) {
      fprintf(stderr,"  all_to_internal: internal error in dim\n");
      abort();
    }

  variable_->update_values(molecule_);
  constant_->update_values(molecule_);
   
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
      fprintf(stderr,"  to_internal: internal error in dim\n");
      abort();
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

  RefSetIntCoor variable_and_constant = new SetIntCoor();
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
IntMolecularCoor::print(ostream& os)
{
  os << indent << "IntMolecularCoor Parameters:\n";
  os << incindent;
  os << indent << "update_bmat = " << (update_bmat_?"yes":"no") << endl;
  os << indent << "scale_bonds = " << scale_bonds_ << endl;
  os << indent << "scale_bends = " << scale_bends_ << endl;
  os << indent << "scale_tors = " << scale_tors_ << endl;
  os << indent << "scale_outs = " << scale_outs_ << endl;
  os << indent << "symmetry_tolerance = " << symmetry_tolerance_ << endl;
  os << indent << "simple_tolerance = " << simple_tolerance_ << endl;
  os << indent << "coordinate_tolerance = " << coordinate_tolerance_ << endl;
  os << indent << "have_fixed_values = " << given_fixed_values_ << endl;
  os << indent << "max_update_steps = " << max_update_steps_ << endl;
  os << indent << "max_update_disp = " << max_update_disp_ << endl;
  os << indent << "have_fixed_values = " << given_fixed_values_ << endl;
  os << decindent;
  
  os << indent << "Molecule:\n";
  os << incindent;
  molecule_->print(os);
  os << decindent;

  print_simples(os);

  os << indent << "Variables:\n";
  os << incindent;
  variable_->print(molecule_,os);
  os << decindent;

  os << indent << "Constants:\n";
  os << incindent;
  constant_->print(molecule_,os);
  os << decindent;
}

void
IntMolecularCoor::print_simples(ostream& os)
{
  if (bonds_->n()) {
    os << indent << "Bonds:\n";
    os << incindent; bonds_->print(molecule_,os); os << decindent;
  }
  if (bends_->n()) {
    os << indent << "Bends:\n";
    os << incindent; bends_->print(molecule_,os); os << decindent;
  }
  if (tors_->n()) {
    os << indent << "Torsions:\n";
    os << incindent; tors_->print(molecule_,os); os << decindent;
  }
  if (outs_->n()) {
    os << indent << "Out of Plane:\n";
    os << incindent; outs_->print(molecule_,os); os << decindent;
  }
  if (extras_->n()) {
    os << indent << "Extras:\n";
    os << incindent; extras_->print(molecule_,os); os << decindent;
  }
  if (fixed_->n()) {
    os << indent << "Fixed:\n";
    os << incindent; fixed_->print(molecule_,os); os << decindent;
  }
  if (followed_.nonnull()) {
    os << indent << "Followed:\n";
    os << incindent; followed_->print(molecule_,os); os << decindent;
  }
}
