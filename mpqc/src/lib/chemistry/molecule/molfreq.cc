
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molfreq.h>

#define CLASSNAME MolecularFrequencies
#define PARENTS public SavableState
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
MolecularFrequencies::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

MolecularFrequencies::MolecularFrequencies(const RefKeyVal& keyval)
{
  mole_ = keyval->describedclassvalue("energy");

  debug_ = keyval->booleanvalue("debug");

  if (mole_.null()) {
      mol_ = keyval->describedclassvalue("molecule");
      kit_ = SCMatrixKit::default_matrixkit();
    }
  else {
      mol_ = mole_->molecule();
      kit_ = mole_->matrixkit();
      d3natom_ = mole_->moldim();
    }
  if (d3natom_.null()) d3natom_ = new SCDimension(3*mol_->natom());

  nirrep_ = mol_->point_group().char_table().nirrep();
  displacements_ = new RefSCMatrix[nirrep_];

  disp_ = keyval->doublevalue("displacement");
  if (keyval->error() != KeyVal::OK) disp_ = 0.001;

  gradients_ = 0;
}

MolecularFrequencies::~MolecularFrequencies()
{
  delete[] displacements_;
  delete[] gradients_;
}

MolecularFrequencies::MolecularFrequencies(StateIn& si):
  SavableState(si),
  original_point_group_(si)
{
  int i;

  mol_.restore_state(si);
  mole_.restore_state(si);

  if (mole_.null()) {
      kit_ = SCMatrixKit::default_matrixkit();
    }
  else {
      kit_ = mole_->matrixkit();
      d3natom_ = mole_->moldim();
    }
  if (d3natom_.null()) d3natom_ = new SCDimension(3*mol_->natom());

  si.get(disp_);
  si.get(ndisp_);
  si.get(nirrep_);
  displacements_ = new RefSCMatrix[nirrep_];

  for (i=0; i < nirrep_; i++) {
      int ndisp;
      si.get(ndisp);
      RefSCDimension ddisp = new SCDimension(ndisp);
      displacements_[i] = matrixkit()->restore_matrix(si,d3natom_,ddisp);
    }

  gradients_ = 0;
}

void
MolecularFrequencies::save_data_state(StateOut& so)
{
  original_point_group_.save_object_state(so);
  mol_.save_state(so);
  mole_.save_state(so);
  so.put(disp_);
  so.put(ndisp_);
  so.put(nirrep_);
  for (int i=0; i < nirrep_; i++) {
      so.put(displacements_[i].ncol());
      displacements_[i].save(so);
    }
}

void
MolecularFrequencies::compute_displacements()
{
  // create the character table for the point group
  CharacterTable ct = mol_->point_group().char_table();

  int ng = ct.order();
  int natom = mol_->natom();
  int nirrep = ct.nirrep();

  original_point_group_ = mol_->point_group();
  original_geometry_ = matrixkit()->vector(d3natom_);

  int i, coor;
  for (i=0, coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          original_geometry_(coor) = mol_->atom(i)[j];
        }
    }

  // Form the matrix of displacements in cartesian coordinates
  RefSCMatrix cartdisp(d3natom_,d3natom_,matrixkit());
  cartdisp.assign(0.0);
  for (i=0; i<3*natom; i++) {
      cartdisp(i,i) = 1.0;
    }

  // Project out translations and rotations
  RefSCDimension dext(new SCDimension(6));
  // form a basis for the translation and rotation coordinates
  RefSCMatrix externaldisp(d3natom_,dext,matrixkit());
  externaldisp.assign(0.0);
  for (i=0; i<natom; i++) {
      AtomicCenter &atom = mol_->atom(i);
      for (int j=0; j<3; j++) {
          externaldisp(i*3 + j,j) = 1.0;
        }
      externaldisp(i*3 + 1, 3 + 0) =  atom[2];
      externaldisp(i*3 + 2, 3 + 0) = -atom[1];
      externaldisp(i*3 + 0, 3 + 1) =  atom[2];
      externaldisp(i*3 + 2, 3 + 1) = -atom[0];
      externaldisp(i*3 + 0, 3 + 2) =  atom[1];
      externaldisp(i*3 + 1, 3 + 2) = -atom[0];
    }
  // do an SVD on the external displacements
  RefSCMatrix Uext(d3natom_,d3natom_,matrixkit());
  RefSCMatrix Vext(dext,dext,matrixkit());
  RefSCDimension min;
  if (d3natom_.n()<dext.n()) min = d3natom_;
  else min = dext;
  int nmin = min.n();
  RefDiagSCMatrix sigmaext(min,matrixkit());
  externaldisp.svd(Uext,sigmaext,Vext);
  // find the epsilon rank
  const double epsilonext = 1.0e-4;
  int rankext = 0;
  for (i=0; i<nmin; i++) {
      if (sigmaext(i) > epsilonext) rankext++;
    }
  cout << "The external rank is " << rankext << endl;
  // find the projection onto the externaldisp perp space
  if (rankext) {
      RefSCDimension drankext_tilde = new SCDimension(d3natom_.n() - rankext);
      RefSCMatrix Uextr_tilde(d3natom_,drankext_tilde,matrixkit());
      Uextr_tilde.assign_subblock(Uext,
                                  0, d3natom_.n()-1,
                                  0, drankext_tilde.n()-1,
                                  0, rankext);
      RefSymmSCMatrix projext_perp(d3natom_, matrixkit());
      projext_perp.assign(0.0);
      projext_perp.accumulate_symmetric_product(Uextr_tilde);
      cartdisp = projext_perp * cartdisp;
    }

  // Form the mapping of atom numbers to transformed atom number
  int **atom_map = new int*[natom];
  for (i=0; i < natom; i++) atom_map[i] = new int[ng];
  // loop over all centers
  for (i=0; i < natom; i++) {
      AtomicCenter ac = mol_->atom(i);
      // then for each symop in the pointgroup, transform the coordinates of
      // center "i" and see which atom it maps into
      for (int g=0; g < ng; g++) {
          double np[3];
          SymmetryOperation so = ct.symm_operation(g);
          for (int ii=0; ii < 3; ii++) {
              np[ii] = 0;
              for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * ac[jj];
            }
          atom_map[i][g] = mol_->atom_at_position(np, 0.05);
          if (atom_map[i][g] < 0) {
              cerr << "MolecularFrequencies: atom mapping bad" << endl;
              abort();
            }
        }
    }

  // Project the cartesian displacements into each irrep
  SymmetryOperation so;
  for (i=0; i<nirrep; i++) {
      IrreducibleRepresentation irrep = ct.gamma(i);
      RefSCMatrix *components = new RefSCMatrix[irrep.degeneracy()];
      // loop over the components of this irrep
      int j;
      for (j=0; j<irrep.degeneracy(); j++) {
          // form the projection matrix for this component of this irrep
          RefSCMatrix projmat(d3natom_,d3natom_,matrixkit());
          projmat.assign(0.0);
          // form the projection matrix for irrep i component j
          // loop over the symmetry operators
          for (int g=0; g < ng; g++) {
              double coef = ((double)irrep.character(g)*irrep.degeneracy())/ng;
              so = ct.symm_operation(g);
              for (int atom=0; atom<natom; atom++) {
                  for (int ii=0; ii < 3; ii++) {
                      for (int jj=0; jj < 3; jj++) {
                          projmat.accumulate_element(atom_map[atom][g]*3+ii,
                                                     atom*3 + jj,
                                                     coef * so(ii,jj));
                        }
                    }
                }
            }
          // projection matrix for irrep i, component j is formed
          RefSCMatrix cartdisp_ij = projmat * cartdisp;
          RefSCMatrix U(d3natom_, d3natom_, matrixkit());
          RefSCMatrix V(d3natom_, d3natom_, matrixkit());
          RefDiagSCMatrix sigma(d3natom_, matrixkit());
          cartdisp_ij.svd(U, sigma, V);
          // Compute the epsilon rank of cartdisp ij
          const double epsilon = 1.0e-3;
          int k, rank = 0;
          for (k=0; k<3*natom; k++) {
              if (sigma(k) > epsilon) rank++;
            }
          if (!rank) continue;
          // Find an orthogonal matrix that spans the range of cartdisp ij
          RefSCDimension drank = new SCDimension(rank);
          RefSCMatrix Ur(d3natom_,drank,matrixkit());
          Ur.assign_subblock(U,0, d3natom_.n()-1, 0, drank.n()-1, 0, 0);
          // Reassign cartdisp_ij to the orthonormal displacement
          cartdisp_ij = Ur;
          if (debug_) {
              cout << "Irrep " << irrep.symbol() << " component " << j << endl;
              cartdisp_ij.print("cartdisp:",cout);
            }
          components[j] = cartdisp_ij;
        }
      int ndisp = 0;
      for (j=0; j<irrep.degeneracy(); j++) ndisp += components[j].ncol();
      RefSCDimension ddisp = new SCDimension(ndisp);
      displacements_[i] = matrixkit()->matrix(d3natom_,ddisp);
      int offset = 0;
      for (j=0; j<irrep.degeneracy(); j++) {
          displacements_[i]->assign_subblock(
              components[j],
              0, d3natom_.n()-1,
              offset, offset+components[j].ncol()-1,
              0, 0);
          offset += components[j].ncol();
        }
      delete[] components;
    }

  for (i=0; i<natom; i++) delete[] atom_map[i];
  delete[] atom_map;

  gradients_ = new RefSCVector[ndisplace()];
}

void
MolecularFrequencies::get_disp(int disp, int &irrep, int &index, double &coef)
{
  int disp_offset = 0;

  // check for +ve totally symmetric displacements
  if (disp < disp_offset + displacements_[0].ncol()) {
      irrep = 0;
      coef = 1.0;
      index = disp - disp_offset;
      return;
    }
  disp_offset += displacements_[0].ncol();
  // check for -ve totally symmetric displacements
  if (disp < disp_offset + displacements_[0].ncol()) {
      irrep = 0;
      coef = -1.0;
      index = disp - disp_offset;
      return;
    }
  disp_offset += displacements_[0].ncol();
  for (int i=1; i<nirrep_; i++) {
      if (disp < disp_offset + displacements_[i].ncol()) {
          irrep = i;
          coef = 1.0;
          index = disp - disp_offset;
          return;
        }
      disp_offset += displacements_[i].ncol();
    }
  cerr << "MolecularFrequencies::get_disp: bad disp number" << endl;
  abort();
}

int
MolecularFrequencies::ndisplace() const
{
  int ndisp = 2 * displacements_[0].ncol();
  for (int i=1; i<nirrep_; i++) {
      ndisp += displacements_[i].ncol();
    }
  return ndisp;
}

void
MolecularFrequencies::displace(int disp)
{
  int irrep, index;
  double coef;
  get_disp(disp, irrep, index, coef);

  if (mole_.nonnull()) mole_->obsolete();

  for (int i=0, coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          mol_->atom(i)[j] = original_geometry_(coor)
                           + coef * disp_ * displacements_[irrep](coor,index);
        }
    }

  if (irrep == 0) {
      mol_->set_point_group(original_point_group_);
    }
  else {
      // Future work: doesn't need to be reduced to c1 symmetry here
      PointGroup pg("c1");
      mol_->set_point_group(pg);
    }
}

void
MolecularFrequencies::set_gradient(int disp, const RefSCVector &grad)
{
  int irrep, index;
  double coef;
  get_disp(disp, irrep, index, coef);

  // transform the gradient into symmetrized coordinates
  gradients_[disp] = displacements_[irrep].t() * grad;
  if (debug_) {
      grad.print("cartesian gradient");
      gradients_[disp].print("internal gradient");
    }
}

void
MolecularFrequencies::compute_frequencies_from_gradients()
{
  int i, coor;

  cout << "Frequencies (cm-1; negative is imaginary):";

  // find the inverse sqrt mass matrix
  RefDiagSCMatrix m(d3natom_, matrixkit());
  for (i=0,coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          m(coor) = 1.0/sqrt(mol_->atom(i).mass()*(1.0/5.48579903e-4));
        }
    }

  RefSymmSCMatrix dhessian;

  RefSymmSCMatrix xhessian;
  if (debug_) {
      xhessian = matrixkit()->symmmatrix(d3natom_);
      xhessian.assign(0.0);
    }

  // start with the totally symmetry frequencies
  RefSCMatrix dtrans = displacements_[0];
  RefSCDimension ddim = dtrans.coldim();
  dhessian = matrixkit()->symmmatrix(ddim);
  for (i=0; i<ddim.n(); i++) {
      for (int j=0; j<=i; j++) {
          dhessian(i,j) = (gradients_[i](j) - gradients_[i+ddim.n()](j)
                         + gradients_[j](i) - gradients_[j+ddim.n()](i))
                         /(4.0*disp_);
        }
    }
  do_freq_for_irrep(0, m, dhessian, xhessian);

  int offset = 2*ddim.n();
  for (int irrep=1; irrep<nirrep_; irrep++) {
      dtrans = displacements_[irrep];
      ddim = dtrans.coldim();
      if (ddim.n() == 0) continue;
      dhessian = matrixkit()->symmmatrix(ddim);
      for (i=0; i<ddim.n(); i++) {
          for (int j=0; j<=i; j++) {
              dhessian(i,j) = (gradients_[i+offset](j)
                             + gradients_[j+offset](i))
                             /(2.0*disp_);
            }
        }
      do_freq_for_irrep(irrep, m, dhessian, xhessian);
      offset += ddim.n();
    }

  if (debug_) {
      xhessian.print("xhessian");

      RefSCMatrix mrect(d3natom_,d3natom_,matrixkit());
      mrect.assign(0.0);
      mrect->accumulate(m.pointer());

      RefSymmSCMatrix mxhessian(d3natom_,matrixkit());
      mxhessian.assign(0.0);
      mxhessian.accumulate_transform(mrect,xhessian);
      mxhessian.print("mass weighted cartesian hessian");
      RefDiagSCMatrix freqs = mxhessian.eigvals();
      // convert the eigvals to frequencies in wavenumbers
      for (i=0; i<freqs.n(); i++) {
          if (freqs(i) >=0.0) freqs(i) = sqrt(freqs(i));
          else freqs(i) = -sqrt(-freqs(i));
          freqs(i) = freqs(i) * 219474.63;
        }
      freqs.print("Frequencies from cartesian hessian");
    }
}

void
MolecularFrequencies::do_freq_for_irrep(int irrep,
                                        const RefDiagSCMatrix &m,
                                        const RefSymmSCMatrix &dhessian,
                                        const RefSymmSCMatrix &xhessian)
{
  int i;
  RefSCMatrix dtrans = displacements_[irrep];
  RefSCDimension ddim = dtrans.coldim();
  if (ddim.n() == 0) return;
  if (debug_) {
      dhessian.print("dhessian");
      dtrans.print("dtrans");
      xhessian.accumulate_transform(dtrans, dhessian);
    }
  // find the basis for the normal coordinates
  RefSCMatrix ncbasis = m * dtrans;
  // use the SVD to orthogonalize and check this basis
  RefSCMatrix basU(d3natom_, d3natom_, matrixkit());
  RefSCMatrix basV(ddim, ddim, matrixkit());
  RefDiagSCMatrix bassigma(ddim, matrixkit());
  ncbasis.svd(basU, bassigma, basV);
  for (i=0; i<ddim.n(); i++) {
      if (bassigma(i) < 1.e0-3) {
          cerr << "MolecularFrequencies: displacements don't span"
               << " normal coordinates"
               << endl;
          abort();
        }
    }
  ncbasis.assign_subblock(basU, 0, d3natom_.n()-1, 0, ddim.n()-1, 0, 0);
  // a transform from disp to x to q (mass weighted x) to disp
  RefSCMatrix dxqd = ncbasis.t() * m * dtrans;
  // transform the dhessian to the mass weighted dhessian
  RefSymmSCMatrix mdhessian = matrixkit()->symmmatrix(dxqd.rowdim());
  mdhessian.assign(0.0);
  mdhessian.accumulate_transform(dxqd, dhessian);
  if (debug_) {
      mdhessian.print("mass weighted dhessian");
    }
  // diagonalize the hessian
  RefDiagSCMatrix freqs = mdhessian.eigvals();
  // convert the eigvals to frequencies in wavenumbers
  for (i=0; i<freqs.n(); i++) {
      if (freqs(i) >=0.0) freqs(i) = sqrt(freqs(i));
      else freqs(i) = -sqrt(-freqs(i));
      freqs(i) = freqs(i) * 219474.63;
    }
  freqs.print(original_point_group_.char_table().gamma(irrep).symbol());
}
