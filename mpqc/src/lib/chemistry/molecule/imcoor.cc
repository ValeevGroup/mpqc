
extern "C" {
#include <math.h>
};

#include <math/scmat/matrix.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/simple.h>

#include <math/topology/bitarray.h>

static void add_bonds(RefSetIntCoor&, BitArray&, Molecule&);
static void add_bends(RefSetIntCoor&, BitArray&, Molecule&);
static void add_tors(RefSetIntCoor&, BitArray&, Molecule&);
static void add_out(RefSetIntCoor&, BitArray&, Molecule&);

static int linear(Molecule&,int,int,int);
static int hterminal(Molecule&, BitArray&, int);

///////////////////////////////////////////////////////////////////////////
// utility functions

// this handles the inverse of matrices even if they are singular
static RefSymmSCMatrix
gen_inverse(RefSymmSCMatrix&mat)
{
  RefSCMatrix vecs(mat.dim(),mat.dim());
  RefDiagSCMatrix vals(mat.dim());

  mat.diagonalize(vals,vecs);


  RefSymmSCMatrix lamd(mat.dim());
  lamd.assign(0.0);
  for(int i=0; i < lamd.n(); i++)
    if(vals(i) > 1.0e-8) lamd(i,i) = 1.0/vals(i);

  RefSymmSCMatrix lam(mat.dim());
  lam.assign(0.0);
  lam.accumulate_transform(vecs,lamd);

  return lam;
}

static inline double
dist(Point& a, Point& b)
{
  return (sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+
               (a[2]-b[2])*(a[2]-b[2])));
}

///////////////////////////////////////////////////////////////////////////
// members of IntMolecularCoor

#define CLASSNAME IntMolecularCoor
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public MolecularCoor
#include <util/state/statei.h>
#include <util/class/classi.h>
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
  simple_tolerance_(1.0e-3),
  coordinate_tolerance_(1.0e-7),
  scale_bonds_(1.0),
  scale_bends_(1.0),
  scale_tors_(1.0),
  scale_outs_(1.0),
  nextra_bonds_(0),
  extra_bonds_(0)
{
  init();
}

IntMolecularCoor::IntMolecularCoor(KeyVal& keyval):
  MolecularCoor(keyval),
  update_bmat_(0),
  only_totally_symmetric_(1),
  symmetry_tolerance_(1.0e-5),
  simple_tolerance_(1.0e-3),
  coordinate_tolerance_(1.0e-7),
  scale_bonds_(1.0),
  scale_bends_(1.0),
  scale_tors_(1.0),
  scale_outs_(1.0)
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
  fixed_ = keyval.describedclassvalue("fixed");
  if (fixed_.null()) fixed_ = new SetIntCoor;

  // the extra_bonds list is given as a vector of atom numbers
  // (atom numbering starts at 1)
  nextra_bonds_ = keyval.count("extra_bonds");
  nextra_bonds_ /= 2;
  if (nextra_bonds_) {
      extra_bonds_ = new int[nextra_bonds_*2];
      for (int i=0; i<nextra_bonds_*2; i++) {
          extra_bonds_[i] = keyval.intvalue("extra_bonds",i);
          if (keyval.error() != KeyVal::OK) {
              fprintf(stderr,"IntMolecularCoor:: keyval CTOR: "
                      "problem reading \"extra_bonds:%d\"\n",i);
              abort();
            }
        }
    }
  else {
      extra_bonds_ = 0;
    }
          

  update_bmat_ = keyval.booleanvalue("update_bmat");

  only_totally_symmetric_ = keyval.booleanvalue("only_totally_symmetric");
  if (keyval.error() != KeyVal::OK) only_totally_symmetric_ = 1;

  double tmp;
  tmp = keyval.doublevalue("scale_bonds");
  if (keyval.error() == KeyVal::OK) scale_bonds_ = tmp;
  tmp = keyval.doublevalue("scale_bends");
  if (keyval.error() == KeyVal::OK) scale_bends_ = tmp;
  tmp = keyval.doublevalue("scale_tors");
  if (keyval.error() == KeyVal::OK) scale_tors_ = tmp;
  tmp = keyval.doublevalue("scale_outs");
  if (keyval.error() == KeyVal::OK) scale_outs_ = tmp;
  tmp = keyval.doublevalue("symmetry_tolerance");
  if (keyval.error() == KeyVal::OK) symmetry_tolerance_ = tmp;
  tmp = keyval.doublevalue("simple_tolerance");
  if (keyval.error() == KeyVal::OK) simple_tolerance_ = tmp;
  tmp = keyval.doublevalue("coordinate_tolerance");
  if (keyval.error() == KeyVal::OK) coordinate_tolerance_ = tmp;
  
  init();
}

IntMolecularCoor::IntMolecularCoor(StateIn& s):
  SavableState(s,IntMolecularCoor::class_desc_),
  MolecularCoor(s)
{
  s.get(nextra_bonds_);
  s.get(extra_bonds_);

  dim_.restore_state(s);
  dnatom3_.restore_state(s);
  dvc_.restore_state(s);

  all_.restore_state(s);

  variable_.restore_state(s);
  constant_.restore_state(s);

  fixed_.restore_state(s);

  bonds_.restore_state(s);
  bends_.restore_state(s);
  tors_.restore_state(s);
  outs_.restore_state(s);
  extras_.restore_state(s);

  s.get(update_bmat_);
  s.get(scale_bonds_);
  s.get(scale_bends_);
  s.get(scale_tors_);
  s.get(scale_outs_);
  s.get(simple_tolerance_);
  s.get(symmetry_tolerance_);
  s.get(coordinate_tolerance_);

  form_coordinates();
}

void
IntMolecularCoor::init()
{
  Molecule& m = *molecule_.pointer();

  // compute needed dimensions
  dnatom3_ = m.dim_natom3();

  // let's go through the geometry and find all the close contacts
  // bonds is a lower triangle matrix of 1's and 0's indicating whether
  // there is a bond between atoms i and j

  BitArray bonds(m.natom(),m.natom());

  for(int i=0; i < m.natom(); i++) {
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

  // compute the simple internal coordinates by type
  add_bonds(bonds_,bonds,m);
  add_bends(bends_,bonds,m);
  add_tors(tors_,bonds,m);
  add_out(outs_,bonds,m);

  form_coordinates();

  // finish computing needed dimensions
  dim_ = new LocalSCDimension(variable_->n());
  dvc_ = new LocalSCDimension(variable_->n()+constant_->n());
}

void
IntMolecularCoor::form_coordinates()
{
  int i;
  int nbonds = bonds_->n();
  int nbends = bends_->n();
  int ntors = tors_->n();
  int nouts = outs_->n();
  int nextras = extras_->n();
  for (i=0; i<nbonds; i++) all_->add(bonds_->coor(i));
  for (i=0; i<nbends; i++) all_->add(bends_->coor(i));
  for (i=0; i<ntors; i++) all_->add(tors_->coor(i));
  for (i=0; i<nouts; i++) all_->add(outs_->coor(i));
  for (i=0; i<nextras; i++) all_->add(extras_->coor(i));

  int nredundant = nbonds + nbends + ntors + nouts + nextras;
  int nfixed = fixed_->n();

  // see how many coords we expect
  int n3 = molecule_->natom()*3;
  int nunique = n3 - 6; // need to detect linear

  if (nredundant < nunique) {
      fprintf(stderr,"IntMolecularCoor::form_coordinates: "
              "found too few redundant coordinates\n");
      fprintf(stderr,"  (the geometry is probably bad)\n");
      abort();
    }

  RefSCDimension dredundant = new LocalSCDimension(nredundant);
  RefSCDimension dfixed = new LocalSCDimension(nfixed);
  RefSCMatrix K; // nredundant x nnonzero
  RefSCMatrix Kfixed; // nfixed x nnonzero
  int* is_totally_symmetric; // nnonzero; if 1 coor has tot. symm. component

  form_K_matrices(dredundant,
                  dfixed,
                  K,
                  Kfixed,
                  is_totally_symmetric);

  RefSCDimension dnonzero = K.coldim();
  int nnonzero = dnonzero.n();

  variable_->clear();
  constant_->clear();

  // put the fixed coordinates into the constant list
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

      for(int j=0; j < nredundant; j++) {
          if(fabs(K(j,i)) > simple_tolerance_) {
              coordinate->add(all_->coor(j),K(j,i));
            }
        }
      // now put the contribution from the fixed coordinates in the
      // coordinate list
      for(j=0; j < nfixed; j++) {
          if(fabs(Kfixed(j,i)) > simple_tolerance_) {
              coordinate->add(fixed_->coor(j),Kfixed(j,i));
            }
        }

      // normalize the coordinate
      coordinate->normalize();

      if (only_totally_symmetric_ && !is_totally_symmetric[i]) {
          constant_->add(coordinate);
        }
      else {
          variable_->add(coordinate);
        }
    }

  constant_->update_values(molecule_);
  variable_->update_values(molecule_);

  fflush(stdout);
  cout << "IntMolecularCoor::form_variable_coordinates()\n"
    << " expected " << nunique << " coordinates\n"
    << " found " << variable_->n() << " variable coordinates\n"
    << " found " << constant_->n() << " constant coordinates\n";
  cout.flush();


  delete[] is_totally_symmetric;
}

// this allocates storage for and computes K, Kfixed, and is_totally_symmetric
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
  RefSCDimension dnonzero = new LocalSCDimension(nonzero);
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
              if(fabs(vecs(j,i)) > simple_tolerance_) nonzero++;
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

  dim_.save_state(s);
  dnatom3_.save_state(s);
  dvc_.save_state(s);

  variable_.save_state(s);
  constant_.save_state(s);

  fixed_.save_state(s);

  bonds_.save_state(s);
  bends_.save_state(s);
  tors_.save_state(s);
  outs_.save_state(s);
  extras_.save_state(s);

  s.put(update_bmat_);
  s.put(scale_bonds_);
  s.put(scale_bends_);
  s.put(scale_tors_);
  s.put(scale_outs_);
  s.put(simple_tolerance_);
  s.put(symmetry_tolerance_);
  s.put(coordinate_tolerance_);
}

RefSCDimension
IntMolecularCoor::dim()
{
  return dim_;
}

void
IntMolecularCoor::to_cartesian(RefSCVector&new_internal)
{
  // get a reference to Molecule for convenience
  Molecule& molecule = *(molecule_.pointer());

  // convergence parameters
  const int maxstep = 20;
  // need to convergence tightly in case a tight optimization is being run
  const double cartesian_tolerance = 1.0e-12;
  // don't bother updating the bmatrix when the error is less than this
  const double update_tolerance = 1.0e-6;

  // compute the internal coordinate displacements
  RefSCVector old_internal(dim_);

  to_internal(old_internal);

  // form the set of all coordinates
  RefSetIntCoor variable_and_constant = new SetIntCoor();
  variable_and_constant->add(variable_);
  variable_and_constant->add(constant_);

  RefSCMatrix bmat(dvc_,dnatom3_);
  RefSymmSCMatrix bmbt(dvc_);
  RefSymmSCMatrix bmbt_i;
  for (int step = 0; step < maxstep; step++) {
      // compute the old internal coordinates
      to_internal(old_internal);

      // the displacements
      RefSCVector displacement = new_internal - old_internal;
      RefSCElementMaxAbs maxabs = new SCElementMaxAbs();
      RefSCElementOp op = maxabs;
      displacement.element_op(op);
      if (maxabs->result() < cartesian_tolerance) {
          return;
        }

      if ((update_bmat_ && (maxabs->result()>update_tolerance)) || step == 0) {
          // form the bmatrix
          variable_and_constant->bmat(molecule_,bmat);
          // form the initial inverse of bmatrix * bmatrix_t
          bmbt.assign(0.0);
          bmbt.accumulate_symmetric_product(bmat);
          bmbt_i = gen_inverse(bmbt);
        }

      // convert displacement to a displacement over all coordinates
      RefSCVector vc_displacement(dvc_);
      vc_displacement.assign(0.0);
      for (int i=0; i<variable_->n(); i++) {
          vc_displacement(i) = displacement(i);
        }

      // compute the cartesian displacements
      RefSCVector cartesian_displacement = bmat.t() * bmbt_i * vc_displacement;

      // update the geometry
      for(i=0; i < dnatom3_.n(); i++) {
          molecule[i/3][i%3] += cartesian_displacement(i);
        }
    }

  fprintf(stderr,"IntMolecularCoor"
          "::to_cartesian(RefSCVector&new_internal):"
          " too many iterations in geometry update\n");
  abort();
}

void
IntMolecularCoor::to_internal(RefSCVector&internal)
{
  variable_->update_values(molecule_);
   
  int n = dim_.n();
  for (int i=0; i<n; i++) {
      internal(i) = variable_->coor(i)->value();
    }
}

void
IntMolecularCoor::to_cartesian(RefSCVector&gradient,RefSCVector&internal)
{
  fprintf(stderr, "IntMolecularCoor::to_cartesian(RefSCVector&,RefSCVector&):"
          " not available\n");
  abort();
}

// converts the gradient in cartesian coordinates to internal coordinates
void
IntMolecularCoor::to_internal(RefSCVector&internal,RefSCVector&gradient)
{
  RefSCMatrix bmat(dvc_,gradient.dim());
  RefSymmSCMatrix bmbt(dvc_);
  RefSymmSCMatrix bmbt_i;

  RefSetIntCoor variable_and_constant = new SetIntCoor();
  variable_and_constant->add(variable_);
  variable_and_constant->add(constant_);

  // form the bmatrix
  variable_and_constant->bmat(molecule_,bmat);
  // form the inverse of bmatrix * bmatrix_t
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(bmat);
  bmbt_i = gen_inverse(bmbt);

  RefSCVector all_internal = bmbt_i*bmat*gradient;

  // put the variable coordinate gradients into internal
  int n = variable_->n();
  for (int i=0; i<n; i++) {
      internal.set_element(i,all_internal.get_element(i));
    }
}

void
IntMolecularCoor::to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal)
{
  fprintf(stderr, "IntMolecularCoor::"
          "to_cartesian(RefSymmSCMatrix&,RefSymmSCMatrix&): not available\n");
  abort();
}

void
IntMolecularCoor::to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart)
{
  fprintf(stderr, "IntMolecularCoor::"
          "to_internal(RefSymmSCMatrix&,RefSymmSCMatrix&): not ready\n");
  abort();
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
  os--;
  
  os.indent() << "Molecule:\n"; os++; molecule_->print(os); os--;

  os.indent() << "Bonds:\n"; os++; bonds_->print(molecule_,os); os--;
  os.indent() << "Bends:\n";  os++; bends_->print(molecule_,os); os--;
  os.indent() << "Torsions:\n";  os++; tors_->print(molecule_,os); os--;
  os.indent() << "Out of Plane:\n";  os++; outs_->print(molecule_,os); os--;
  os.indent() << "Extras:\n";  os++; extras_->print(molecule_,os); os--;
  os.indent() << "Fixed:\n";  os++; fixed_->print(molecule_,os); os--;

  os.indent() << "Variables:\n"; os++; variable_->print(molecule_,os); os--;
  os.indent() << "Constants:\n"; os++; constant_->print(molecule_,os); os--;

}

void
IntMolecularCoor::guess_hessian(RefSymmSCMatrix&hessian)
{
  hessian.assign(0.0);

  for (int i=0; i<variable_->n(); i++) {
      hessian(i,i) = variable_->coor(i)->force_constant(molecule_);
    }
}
