
extern "C" {
#include <math.h>
};

#include <math/scmat/matrix.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/localdef.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/simple.h>

#include <math/topology/bitarray.h>


///////////////////////////////////////////////////////////////////////////
// members of SymmMolecularCoor

#define CLASSNAME SymmMolecularCoor
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public IntMolecularCoor
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SymmMolecularCoor::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = IntMolecularCoor::_castdown(cd);
  return do_castdowns(casts,cd);
}

SymmMolecularCoor::SymmMolecularCoor(RefMolecule&mol):
  IntMolecularCoor(mol)
{
  init();
}

SymmMolecularCoor::SymmMolecularCoor(const RefKeyVal& keyval):
  IntMolecularCoor(keyval)
{
  init();
}

SymmMolecularCoor::SymmMolecularCoor(StateIn& s):
  IntMolecularCoor(s)
{
}

SymmMolecularCoor::~SymmMolecularCoor()
{
}

void
SymmMolecularCoor::save_data_state(StateOut&s)
{
  IntMolecularCoor::save_data_state(s);
}

void
SymmMolecularCoor::form_coordinates()
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

  // if we're following coordinates, add them to the fixed list
  if (followed_.nonnull())
    fixed_->add(followed_);
  
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

  // now remove followed coords from the fixed list, and add to the
  // variable list
  if (followed_.nonnull()) {
    fixed_->del(followed_);
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
              RefIntCoor c = all_->coor(j);
              coordinate->add(c,K(j,i));
            }
        }
      // now put the contribution from the fixed coordinates in the
      // coordinate list
      for(j=0; j < nfixed; j++) {
          if(pow(Kfixed(j,i),2.0) > simple_tolerance_) {
              RefIntCoor c = fixed_->coor(j);
              coordinate->add(c,Kfixed(j,i));
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
  cout << "  IntMolecularCoor::form_variable_coordinates()\n"
    << "    expected " << nunique << " coordinates\n"
    << "    found " << variable_->n() << " variable coordinates\n"
    << "    found " << constant_->n() << " constant coordinates\n";
  cout.flush();


  delete[] is_totally_symmetric;
}

void
SymmMolecularCoor::guess_hessian(RefSymmSCMatrix&hessian)
{
  // first form diagonal hessian in redundant internal coordinates
  RefSCDimension rdim = new LocalSCDimension(all_->n());
  RefSymmSCMatrix rhessian(rdim);
  rhessian.assign(0.0);
  all_->guess_hessian(molecule_,rhessian);

  // create redundant coordinate bmat
  RefSCDimension dn3 = molecule_->dim_natom3();
  RefSCMatrix bmatr(rdim,dn3);
  all_->bmat(molecule_,bmatr);

  // then form the variable coordinate bmat
  RefSCDimension dredundant = new LocalSCDimension(variable_->n());
  RefSCMatrix bmat(dredundant,dn3);
  variable_->bmat(molecule_,bmat);

  // and (B*B+)^-1
  RefSymmSCMatrix bmbt(dredundant);
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
