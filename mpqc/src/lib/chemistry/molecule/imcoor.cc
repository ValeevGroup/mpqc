
extern "C" {
#include <math.h>
};

#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <chemistry/molecule/localdef.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/simple.h>

#include <math/topology/bitarray.h>

#define DEFAULT_SIMPLE_TOLERANCE 1.0e-3

#define USE_SVD 1

static void add_bonds(RefSetIntCoor&, BitArray&, Molecule&);
static void add_bends(RefSetIntCoor&, BitArray&, Molecule&);
static void add_tors(RefSetIntCoor&, BitArray&, Molecule&);
static void add_out(RefSetIntCoor&, BitArray&, Molecule&);

static int linear(Molecule&,int,int,int);
static int hterminal(Molecule&, BitArray&, int);
static int nearest_contact(int,Molecule&);

///////////////////////////////////////////////////////////////////////////
// members of IntMolecularCoor

#define CLASSNAME IntMolecularCoor
#define VERSION 2
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
  nextra_bonds_(0),
  extra_bonds_(0),
  max_update_steps_(100),
  max_update_disp_(0.5),
  given_fixed_values_(0)
{
  new_coords();
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
  scale_outs_(1.0)
{
  // intialize the coordinate sets
  new_coords();

  // actually read the keyval info
  read_keyval(keyval);
}

IntMolecularCoor::IntMolecularCoor(StateIn& s):
  MolecularCoor(s)
{
  s.get(nextra_bonds_);
  s.get(extra_bonds_);

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

  given_fixed_values_ = keyval->booleanvalue("have_fixed_values");

  max_update_steps_ = keyval->intvalue("max_update_steps");
  if (keyval->error() != KeyVal::OK) max_update_steps_ = 100;

  max_update_disp_ = keyval->doublevalue("max_update_disp");
  if (keyval->error() != KeyVal::OK) max_update_disp_ = 0.5;

  // the extra_bonds list is given as a vector of atom numbers
  // (atom numbering starts at 1)
  nextra_bonds_ = keyval->count("extra_bonds");
  nextra_bonds_ /= 2;
  if (nextra_bonds_) {
      extra_bonds_ = new int[nextra_bonds_*2];
      for (int i=0; i<nextra_bonds_*2; i++) {
          extra_bonds_[i] = keyval->intvalue("extra_bonds",i);
          if (keyval->error() != KeyVal::OK) {
              fprintf(stderr,"IntMolecularCoor:: keyval CTOR: "
                      "problem reading \"extra_bonds:%d\"\n",i);
              abort();
            }
        }
    }
  else {
      extra_bonds_ = 0;
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
  Molecule& m = *molecule_.pointer();

  // let's go through the geometry and find all the close contacts
  // bonds is a lower triangle matrix of 1's and 0's indicating whether
  // there is a bond between atoms i and j

  BitArray bonds(m.natom(),m.natom());

  int i;
  for(i=0; i < m.natom(); i++) {
      double at_rad_i = m[i].element().atomic_radius();

      for(int j=0; j < i; j++) {
          double at_rad_j = m[j].element().atomic_radius();

          if (dist(m[i].point(),m[j].point()) < 1.1*(at_rad_i+at_rad_j))
            bonds.set(i,j);
        }
    }

  for (i=0; i<nextra_bonds_; i++) {
      bonds.set(extra_bonds_[i*2]-1,extra_bonds_[i*2+1]-1);
    }

  // check for atoms bound to nothing
  for (i=0; i < m.natom(); i++) {
    int bound=0;
    for (int j=0; j < m.natom(); j++) {
      if (bonds(i,j)) {
        bound=1;
        break;
      }
    }
    if (!bound) {
      int j = nearest_contact(i,m);
      fprintf(stderr,"\n  Warning!:  atom %d is not bound to anything.\n",i+1);
      fprintf(stderr,
              "             You may wish to add an entry to extra_bonds.\n");
      fprintf(stderr,
              "             Atom %d is only %f angstroms away...\n\n",j+1,
              bohr*dist(m[i].point(),m[j].point()));
    }
  }
      
  // compute the simple internal coordinates by type
  add_bonds(bonds_,bonds,m);
  add_bends(bends_,bonds,m);
  add_tors(tors_,bonds,m);
  add_out(outs_,bonds,m);

  if (given_fixed_values_) {
      // save the given coordinate values
      RefSCDimension original_dfixed
          = matrixkit_->dimension(fixed_->n(),"Nfix");
      RefSCVector given_fixed_coords(original_dfixed);
      for (i=0; i<original_dfixed.n(); i++) {
          given_fixed_coords(i) = fixed_->coor(i)->value();
        }

      // break up the displacement into several manageable steps
      double maxabs = given_fixed_coords.maxabs();
      int nstep = int(maxabs/max_update_disp_) + 1;
      given_fixed_coords.scale(1.0/nstep);
      printf("IntMolecularCoor: displacing fixed coordinates to the"
             " requested values in %d steps\n", nstep);
      for (int istep=1; istep<=nstep; istep++) {
          form_coordinates();

          dim_ = matrixkit_->dimension(variable_->n(), "Nvar");
          dvc_ = matrixkit_->dimension(variable_->n()+constant_->n(),
                                       "Nvar+Nconst");

          RefSCVector new_internal_coordinates(dvc_);
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

  dim_ = matrixkit_->dimension(variable_->n(), "Nvar");
  dvc_ = matrixkit_->dimension(variable_->n()+constant_->n(),
                               "Nvar+Nconst");
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
  RefSCDimension dcoor = matrixkit->dimension(coor->n());
  RefSCMatrix B(dcoor, dnatom3);
  coor->bmat(molecule, B);

  // Project out the previously discovered internal coordinates
  if (projection.nonnull()) {
      B = B * projection;
    }

  // Compute the singular value decomposition of B
  RefSCMatrix U(dcoor,dcoor);
  RefSCMatrix V(dnatom3,dnatom3);
  RefSCDimension min;
  if (dnatom3.n()<dcoor.n()) min = dnatom3;
  else min = dcoor;
  int nmin = min.n();
  RefDiagSCMatrix sigma(min);
  B.svd(U,sigma,V);

  // Compute the epsilon rank of B
  int i, rank = 0;
  for (i=0; i<nmin; i++) {
      if (sigma(i) > epsilon) rank++;
    }

  RefSCMatrix SIGMA(dcoor, dnatom3);
  SIGMA.assign(0.0);
  for (i=0; i<nmin; i++) {
      SIGMA(i,i) = sigma(i);
    }

  // return if there are no new coordinates
  if (rank==0) return 0;

  // Find an orthogonal matrix that spans the range of B
  RefSCMatrix Ur;
  RefSCDimension drank = matrixkit->dimension(rank);
  if (rank) {
      Ur = matrixkit->matrix(dcoor,drank);
      Ur.assign_subblock(U,0, dcoor.n()-1, 0, drank.n()-1, 0, 0);
    }

  // Find an orthogonal matrix that spans the null space of B
  int rank_tilde = dnatom3.n() - rank;
  RefSCMatrix Vr_tilde;
  RefSCDimension drank_tilde = matrixkit->dimension(rank_tilde);
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
      printf("found %d totally symmetric coordinates\n", ntotally_symmetric);

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
  int i;

  // The cutoff for whether or not a coordinate is considered totally symmetric
  double ts_eps = 0.0001;

  // The geometry will be needed to check for totally symmetric
  // coordinates
  RefSCVector geom(dnatom3_);
  for(i=0; i < geom.n()/3; i++) {
      geom(3*i) = molecule_->operator[](i).point()[0];
      geom(3*i+1) = molecule_->operator[](i).point()[1];
      geom(3*i+2) = molecule_->operator[](i).point()[2];
    }

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
          fprintf(stderr, "ERROR: IntMolecularCoor: nfixed = %d rank = %d\n",
                  dfixed.n(), Ktmp.ncol());
          abort();
        }
      // check that fixed coordinates be totally symmetric
      if (Ktmp.nrow() != count_nonzero(totally_symmetric_fixed, ts_eps)) {
          fprintf(stderr, "WARNING: only %d of %d fixed coordinates are"
                  "totally symmetric\n",
                  count_nonzero(totally_symmetric_fixed, ts_eps), dfixed.n());
        }

      // Compute Kfixed
      RefSCDimension dcoor = matrixkit_->dimension(all_->n());
      RefSCMatrix B(dcoor, dnatom3_);
      all_->bmat(molecule_, B);
      Kfixed = B * null_bfixed_perp;
    }

  RefSCVector totally_symmetric_all;
  form_partial_K(all_, molecule_, geom, 0.001, dnatom3_, matrixkit_,
                 projection, totally_symmetric_all, K);
  int n_totally_symmetric_all = count_nonzero(totally_symmetric_all, ts_eps);
  is_totally_symmetric = new int[K.ncol()];
  for (i=0; i<K.ncol(); i++) {
      if (fabs(totally_symmetric_all(i)) > ts_eps) is_totally_symmetric[i] = 1;
      else is_totally_symmetric[i] = 0;
    }
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
  RefSCMatrix bmat(dredundant,dnatom3_);
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
  RefSymmSCMatrix bmbt(dredundant);
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(bmat);

  // form the fixed part of the b matrix
  RefSCMatrix bmat_fixed(dfixed,dnatom3_);
  fixed_->bmat(molecule_,bmat_fixed);
  int nfixed = dfixed.n();

  // and form the fixed b*b~
  RefSymmSCMatrix bmbt_fixed(dfixed);
  bmbt_fixed.assign(0.0);
  bmbt_fixed.accumulate_symmetric_product(bmat_fixed);

  // need the cross terms in bmbt also
  RefSCMatrix bmbt_fix_red;
  if (nfixed != 0) bmbt_fix_red = bmat_fixed * bmat.t();

  // orthogonalize the redundant coordinates to the fixed coordinates
  RefSCMatrix redundant_ortho(dfixed,dredundant);
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
  RefSCMatrix vecs(dredundant,dredundant);
  RefDiagSCMatrix vals(dredundant);

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
  RefSCDimension dnonzero = matrixkit_->dimension(nonzero, "Nnonzero");
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
              fprintf(stderr,"Geom_form_K: no nonzero simple coordinates");
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
  if (extra_bonds_) delete[] extra_bonds_;
}

void
IntMolecularCoor::save_data_state(StateOut&s)
{
  MolecularCoor::save_data_state(s);

  s.put(nextra_bonds_);
  s.put(extra_bonds_,2*nextra_bonds_);

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
  RefSCVector old_internal(dvc_);

  all_to_internal(old_internal);

  RefSCMatrix internal_to_cart_disp;
  for (int step = 0; step < max_update_steps_; step++) {
      // compute the old internal coordinates
      all_to_internal(old_internal);

      // the displacements
      RefSCVector displacement = new_internal - old_internal;
      RefSCElementMaxAbs maxabs = new SCElementMaxAbs();
      RefSCElementOp op = maxabs;
      displacement.element_op(op);
      if (maxabs->result() < cartesian_tolerance_) {

          constant_->update_values(molecule_);
          variable_->update_values(molecule_);

          return 0;
        }

      if ((update_bmat_ && (maxabs->result()>update_tolerance))
          || internal_to_cart_disp.null()) {
          int i;
          RefSCMatrix bmat(dvc_,dnatom3_);

          // form the set of all coordinates
          RefSetIntCoor variable_and_constant = new SetIntCoor();
          variable_and_constant->add(variable_);
          variable_and_constant->add(constant_);

          // form the bmatrix
          variable_and_constant->bmat(molecule_,bmat);

          // Compute the singular value decomposition of B
          RefSCMatrix U(dvc_,dvc_);
          RefSCMatrix V(dnatom3_,dnatom3_);
          RefSCDimension min;
          if (dnatom3_.n()<dvc_.n()) min = dnatom3_;
          else min = dvc_;
          int nmin = min.n();
          RefDiagSCMatrix sigma(min);
          bmat.svd(U,sigma,V);

          // compute the epsilon rank of B
          int rank = 0;
          for (i=0; i<nmin; i++) {
              if (fabs(sigma(i)) > 0.0001) rank++;
            }

          RefSCDimension drank = matrixkit_->dimension(rank);
          RefDiagSCMatrix sigma_i(drank);
          for (i=0; i<rank; i++) {
              sigma_i(i) = 1.0/sigma(i);
            }
          RefSCMatrix Ur(dvc_, drank);
          RefSCMatrix Vr(dnatom3_, drank);
          Ur.assign_subblock(U,0, dvc_.n()-1, 0, drank.n()-1, 0, 0);
          Vr.assign_subblock(V,0, dnatom3_.n()-1, 0, drank.n()-1, 0, 0);
          internal_to_cart_disp = Vr * sigma_i * Ur.t();

#if !USE_SVD
          RefSymmSCMatrix bmbt(dvc_);
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
    }

  cout.flush();
  cerr.flush();
  fflush(stdout);
  fflush(stderr);

  fprintf(stderr,"WARNING: IntMolecularCoor"
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
IntMolecularCoor::to_cartesian(RefSCVector&new_internal)
{
  RefSCVector all_internal(dvc_);

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
  RefSCMatrix bmat(dim_,gradient.dim());
  variable_->bmat(molecule_,bmat);

  gradient = bmat.t() * internal;
  
  return 0;
}

// converts the gradient in cartesian coordinates to internal coordinates
int
IntMolecularCoor::to_internal(RefSCVector&internal,RefSCVector&gradient)
{
  RefSCMatrix bmat(dvc_,gradient.dim());
  RefSymmSCMatrix bmbt(dvc_);

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
  RefSCMatrix bmat(dim_,cart.dim());
  variable_->bmat(molecule_,bmat);
  cart.accumulate_transform(bmat.t(),internal);
  return 0;
}

int
IntMolecularCoor::to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart)
{
  // form bmat
  RefSCMatrix bmat(dim_,cart.dim());
  variable_->bmat(molecule_,bmat);
  // and (B*B+)^-1
  RefSymmSCMatrix bmbt(dim_);
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

///////////////////////////////////////////////////////////////////////////
// auxillary functions of IntMolecularCoor

/*
 * the following are translations of functions written by Gregory Humphreys
 * at the NIH
 */

/*
 * for each bonded pair, add an entry to the simple coord list
 */

static void
add_bonds(RefSetIntCoor& list, BitArray& bonds, Molecule& m)
{
  int i,j,ij;
  int labelc=0;
  char label[80];

  for(i=ij=0; i < m.natom(); i++) {
    for(j=0; j <= i; j++,ij++) {
      if(bonds[ij]) {
        labelc++;
        sprintf(label,"s%d",labelc);
        list->add(new Stre(label,j+1,i+1));
        }
      }
    }
  }

/*
 * for each pair of bonds sharing an atom, add a bend to the simple
 * coord list
 *
 * check each bend to see if it is linear.  if so, then we'll have to add
 * in-plane and out-of-plane linear bends as well
 *
 * let's do this later...I think I only want to do this for symmetric
 * tops, but I'm not sure...anyway, for now, let's just eliminate linear
 * bends so things like co2 won't blow up
 */

static int
linear(Molecule& m, int i, int j, int k)
{
  double dij = dist(m[i].point(),m[j].point());
  double dik = dist(m[i].point(),m[k].point());
  double djk = dist(m[j].point(),m[k].point());

 // if dij+djk==dik, then this bug is linear
  if((dij+djk - dik) < 1.0e-5) return 1;
  
  return 0;
  }

static void
add_bends(RefSetIntCoor& list, BitArray& bonds, Molecule& m)
{
  int i,j,k;
  int labelc=0;
  char label[80];

  int n = m.natom();

  for(i=0; i < n; i++) {
    for(j=0; j < n; j++) {
      if(bonds(i,j)) {
        for(k=0; k < i; k++) {
          if(bonds(j,k)) {
            if(linear(m,i,j,k)) continue;

	    labelc++;
	    sprintf(label,"b%d",labelc);
	    list->add(new Bend(label,k+1,j+1,i+1));
	    }
          }
	}
      }
    }
  }

/*
 * for each pair of bends which share a common bond, add a torsion
 */

/*
 * just look at the heavy-atom skeleton. return true if i is a terminal
 * atom.
 */

static int
hterminal(Molecule& m, BitArray& bonds, int i)
{
  int nh=0;
  for (int j=0; j < m.natom(); j++)
    if (bonds(i,j) && m[j].element().mass() > 1.1) nh++;
  return (nh==1);
}

static void
add_tors(RefSetIntCoor& list, BitArray& bonds, Molecule& m)
{
  int i,j,k,l;
  int labelc=0;
  char label[80];

  int n = m.natom();

  for(j=0; j < n; j++) {
    for(k=0; k < j; k++) {
      if(bonds(j,k)) {
        for(i=0; i < n; i++) {
          if(k==i) continue;

         // no hydrogen torsions, ok?
	  if (m[i].element().mass() < 1.1 && !hterminal(m,bonds,j)) continue;

          if (bonds(j,i)) {
	    if (linear(m,i,j,k)) continue;

            for (l=0; l < n; l++) {
              if (l==j || l==i) continue;

             // no hydrogen torsions, ok?
	      if (m[l].element().mass() < 1.1 && !hterminal(m,bonds,k))
                continue;

              if (bonds(k,l)) {
		if(linear(m,j,k,l)) continue;

		labelc++;
		sprintf(label,"t%d",labelc);
		list->add(new Tors(label,l+1,k+1,j+1,i+1));
		}
	      }
	    }
          }
	}
      }
    }
  }

static void
add_out(RefSetIntCoor& list, BitArray& bonds, Molecule& m)
{
  int i,j,k,l;
  int labelc=0;
  char label[80];

  int n = m.natom();

 // first find all tri-coordinate atoms
  for(i=0; i < n; i++) {
    if(bonds.degree(i)!=3) continue;

   // then look for terminal atoms connected to i
    for(j=0; j < n; j++) {
      if(bonds(i,j) && bonds.degree(j)==1) {

        for(k=0; k < n; k++) {
          if(k!=j && bonds(i,k)) {
            for(l=0; l < k; l++) {
              if(l!=j && bonds(i,l)) {
		labelc++;
		sprintf(label,"o%d",labelc);
		list->add(new Out(label,j+1,i+1,l+1,k+1));
		}
	      }
	    }
          }
	}
      }
    }
  }

void
IntMolecularCoor::print(SCostream& os)
{
  os.indent() << "Parameters:\n";
  os++;
  os.indent() << "update_bmat = " << (update_bmat_?"yes":"no") << endl;
  os.indent() << "scale_bonds = " << scale_bonds_ << endl;
  os.indent() << "scale_bends = " << scale_bends_ << endl;
  os.indent() << "scale_tors = " << scale_tors_ << endl;
  os.indent() << "scale_outs = " << scale_outs_ << endl;
  os.indent() << "symmetry_tolerance = " << symmetry_tolerance_ << endl;
  os.indent() << "simple_tolerance = " << simple_tolerance_ << endl;
  os.indent() << "coordinate_tolerance = " << coordinate_tolerance_ << endl;
  os.indent() << "have_fixed_values = " << given_fixed_values_ << endl;
  os.indent() << "max_update_steps = " << max_update_steps_ << endl;
  os.indent() << "max_update_disp = " << max_update_disp_ << endl;
  os.indent() << "have_fixed_values = " << given_fixed_values_ << endl;
  os--;
  
  os.indent() << "Molecule:\n"; os++; molecule_->print(os); os--;

  print_simples(os);

  os.indent() << "Variables:\n"; os++; variable_->print(molecule_,os); os--;
  os.indent() << "Constants:\n"; os++; constant_->print(molecule_,os); os--;
}

void
IntMolecularCoor::print_simples(SCostream& os)
{
  if (bonds_->n()) {
    os.indent() << "Bonds:\n"; os++; bonds_->print(molecule_,os); os--;
  }
  if (bends_->n()) {
    os.indent() << "Bends:\n";  os++; bends_->print(molecule_,os); os--;
  }
  if (tors_->n()) {
    os.indent() << "Torsions:\n";  os++; tors_->print(molecule_,os); os--;
  }
  if (outs_->n()) {
    os.indent() << "Out of Plane:\n";  os++; outs_->print(molecule_,os); os--;
  }
  if (extras_->n()) {
    os.indent() << "Extras:\n";  os++; extras_->print(molecule_,os); os--;
  }
  if (fixed_->n()) {
    os.indent() << "Fixed:\n";  os++; fixed_->print(molecule_,os); os--;
  }
  if (followed_.nonnull()) {
    os.indent() << "Followed:\n"; os++; followed_->print(molecule_,os); os--;
  }
}

static int
nearest_contact(int i, Molecule& m)
{
  double d=-1.0;
  int n=0;
  
  for (int j=0; j < m.natom(); j++) {
    double td = dist(m[i].point(),m[j].point());
    if (j==i)
      continue;
    else if (d < 0 || td < d) {
      d = td;
      n = j;
    }
  }
  
  return n;
}
