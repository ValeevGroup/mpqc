
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/optimize/diis.h>
#include <chemistry/qc/scf/grscf.h>
#include <chemistry/qc/integral/integralv2.h>

///////////////////////////////////////////////////////////////////////////
// GRSCF

#define CLASSNAME GRSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public OneBodyWavefunction
#include <util/class/classi.h>
void *
GRSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

GRSCF::GRSCF(StateIn& s) :
  OneBodyWavefunction(s)
{
  _extrap.restore_state(s);
  _data.restore_state(s);
  _error.restore_state(s);

  _accumdih.restore_state(s);
  _accumddh.restore_state(s);
  _accumeffh.restore_state(s);

  s.get(_ndocc);
  s.get(_nsocc);
  s.get(_density_reset_freq);
  s.get(_maxiter);
  s.get(_eliminate);

  s.getstring(ckptdir);
  s.getstring(fname);
}

GRSCF::GRSCF(const RefKeyVal& keyval) :
  OneBodyWavefunction(keyval)
{
  _maxiter = keyval->intvalue("maxiter");
  if (keyval->error() != KeyVal::OK) {
    _maxiter = 100;
  }

  _extrap = keyval->describedclassvalue("extrap");
  if (_extrap.null()) {
    _extrap = new DIIS;
  }

  _ndocc = keyval->intvalue("ndocc");
  if (keyval->error() != KeyVal::OK) {
    double charge = 0.0;
    for (int i=0; i<molecule()->natom(); i++) {
      charge += molecule()->atom(i).element().charge();
    }
    _ndocc = (int) (charge + 0.5);
  }

  _nsocc = keyval->intvalue("nsocc");

  init();
}

GRSCF::~GRSCF()
{
}

void
GRSCF::init()
{
}

void
GRSCF::compute()
{
  fprintf(stderr,"GRSCF::compute(): don't call me yet\n");
  abort();
}

RefSCMatrix
GRSCF::eigenvectors()
{
  return _eigenvectors;
}

void
GRSCF::save_data_state(StateOut& s)
{
  _extrap.save_state(s);
  _data.save_state(s);
  _error.save_state(s);

  _accumdih.save_state(s);
  _accumddh.save_state(s);
  _accumeffh.save_state(s);

  s.put(_ndocc);
  s.put(_nsocc);
  s.put(_density_reset_freq);
  s.put(_maxiter);
  s.put(_eliminate);

  s.putstring(ckptdir);
  s.putstring(fname);
}

double
GRSCF::occupation(int i)
{
  if (i < _ndocc) return 2.0;
  if (i < _ndocc + _nsocc) return 1.0;
  return 0.0;
}

#if 0
void
GRSCF::form_density_independent_h()
{
  RefSymmSCMatrix& h = _density_independent_h;
  if (h.null()) {
    h = _orth_dim->create_symmmatrix();
  }

  if (_density_independent_h.needed()) {
    h.assign(0.0);
    RefSCElementOp op = new GaussianKineticIntv2(basis(), molecule());
    h.element_op(op);
    op = new GaussianNuclearIntv2(basis(), molecule());
    h.element_op(op);
  }
}
#endif

void
GRSCF::converge_eigenvectors()
{
#if 0
  double etot, edif, neelec, delta;
  double plimit = 0.1 * _energy.desired_accuracy();

  _extrap->set_tolerance(plimit);

  RefDiagSCMatrix val(basis_dimension());
  RefSCMatrix vec(basis_dimension(), basis_dimension());
  RefSymmSCMatrix h(basis_dimension());
  RefSymmSCMatrix dih(basis_dimension());
  RefSymmSCMatrix h_open;
  RefSymmSCMatrix P(basis_dimension());
  RefSymmSCMatrix P_open;
  RefSymmSCMatrix DP(basis_dimension());
  RefSymmSCMatrix DP_open;

  if (_nsocc) {
    P_open = basis_dimension()->create_symmmatrix();
    DP_open = basis_dimension()->create_symmmatrix();
    h_open = basis_dimension()->create_symmmatrix();
  }

  _accumdih->init(basis(),_mol);
  _accumdih->accum(dih,h);
  _accumdih->done();

  fprintf(outfile, "GRSCF: AccumDIH = %s\n", _accumdih->class_name());
  fprintf(outfile, "GRSCF: AccumDDH = %s\n", _accumddh->class_name());

  fprintf(outfile,"\n  iter       total energy       "
          " delta E         delta P          error\n");


  double E_nuc = nuclear_energy();

  _accumeffh->docc(0, _ndocc);
  _accumeffh->socc(_ndocc, _ndocc + _nsocc);

  _accumddh->init(basis(),_mol);
  int iter = 0;
  double old_E_elec = 0.0;
  while (++iter < _maxiter) {
    form_density(iter, vec, P, P_open, DP, DP_open);

    h.assign(dih);
    if (h_open.nonnull()) h_open.assign(dih);

    _accumddh->accum(plimit, P, P_open, DP, DP_open, h, h_open);

    double E_elec = electronic_energy(dih, P, P_open, h, h_open);
    double E_dif = E_elec - old_E_elec;
    old_E_elec = E_elec;

    delta = rms_delta_density(DP, DP_open);

    fprintf(outfile, "%5d %20.10f %15.6e %15.6e %15.6e\n",
            iter+1, E_nuc + E_elec, E_dif, delta, _extrap->error());
    fflush(outfile);

    if (delta < plimit) {
      fprintf(outfile,"\n  converged scf energy is %20.10f au\n",
              E_tot);
      break;
    }

    if (((iter+1) % _ckpt_freq) == 0) checkpoint(vec);

    RefSCMatrix aovec(_orth_dim);
    aovec.assign(0.0);
    aovec.accumulate_product(vec, smhalf);

    RefSymmSCMatrix hmo(_orth_dim);
    hmo.assign(0.0);
    hmo.accumulate_transform(aovec, h);

    // convert hmo to an effective hmo
    if (hmo_open.nonnull()) hmo.element_op(_accumeffh, hmo_open);

    // convert the effective hmo back to the orthogonal ao basis
    RefSymmSCMatrix hoao(_orth_dim);
    hoao.assign(0.0);
    hoao.accum_transform(vec.t(), hmo);

    RefSCExtrapError error = form_error_matrix(hmo);
    RefSCExtrapData data = new SymmSCMatrixExtrapData(hoao);
    _extrap->extrapolate(data, error);

    // back to the mo basis so a level shift can be applied
    hmo.assign(0.0);
    hmo.accumulate_tranform(vec, hoao);
    level_shift(hmo);
    RefSCMatrix vec2(_orth_dim, _orth_dim);
    hmo.diagonalize(val, vec2);
    RefSCMatrix old_vec = vec.copy();
    vec.assign(0.0);
    vec.accumulate_transform(old_vec, vec2);
    vec.orthogonalize_this();
  }
  _accumddh->done();
#endif
}

void
GRSCF::form_density(int iter,
                    const RefSCMatrix& vec,
                    const RefSymmSCMatrix& P,
                    const RefSymmSCMatrix& P_open,
                    const RefSymmSCMatrix& DP,
                    const RefSymmSCMatrix& DP_open)
{
#if 0
  RefDiagSCMatrix occ(_orth_dim);
  RefSCElementAssignBlock occop = new SCElementAssignBlock();
  occop->value(2.0);
  occop->block(0, _ndocc, 0, _ndocc);
  occ.element_op(occup);
  occop->value(1.0);
  occop->block(_ndocc, _ndocc+_nsocc, _ndocc, _ndocc+_nsocc);
  occ.element_op(occup);
  occop->value(0.0);
  occop->block(_ndocc+_nsocc, _orth_dim.n(),
               _ndocc+_nsocc, _orth_dim.n());
  occ.element_op(occup);

  // P is the total density
  DP.assign(P);
  DP.scale(-1.0);
  P.assign(0.0);
  P.accumulate_transform(vec, occ);

  DP.accumulate(P);

  // P_open is the open shell density
  if (P_open.nonnull()) {
    DP_open.assign(P_open);
    DP_open.scale(-1.0);
    P_open.assign(0.0);
    occop->value(0.0);
    occop->block(0, _ndocc, 0, _ndocc);
    occ.element_op(occup);
    P_open.accumulate_transform(vec, occ);

    DP_open.accumulate(P_open);
  }

  if (iter && _eliminate && ((iter)%_density_reset_freq == 0)) {
    fprintf(outfile,"  GRSCF: resetting density matrices\n");
    DP.assign(P);
    P.assign(0.0);
    h.assign(dih);
    if (P_open.nonnull()) {
      DP_open.assign(P_open);
      P_open.assign(0.0);
      h_open.assign(0.0);
    }
  }
#endif
}

#if 0
RefSymmSCMatrix
GRSCF::form_error_matrix(const RefSymmSCMatrix& hmo,
                         const RefSCMatrix& vec)
{
  RefSymmSCMatrix errormat(_orth_dim);
  errormat.assign(hmo);

  // zero out arbitrary blocks of the error matrix
  RefSCElementAssignBlock zero = new SCElementAssignBlock();
  zero->value(0.0);
  zero->block(0, _ndocc, 0, _ndocc);
  errormat.element_op(zero);
  zero->block(_ndocc, _ndocc+_nsocc, _ndocc, _ndocc+_nsocc);
  errormat.element_op(zero);
  zero->block(_ndocc+_nsocc, _orth_dim.n(),
              _ndocc+_nsocc, _orth_dim.n());
  errormat.element_op(zero);

  // form the error matrix in the orthogonal AO basis
  RefSymmSCMatrix oaoerrormat(_orth_dim);
  oaoerrormat.assign(0.0);
  oaoerrormat.accumulate_transform(errormat, vec);

  RefSCExtrapError error = new SymmSCMatrixExtraData(aoaerrormat);

  return error;
}
#endif

#if 0
void
GRSCF::checkpoint(const RefSCMatrix& vec)
{
  char ckptfile[512];
  sprintf(ckptfile,"%s%s.scfvec",ckptdir,fname);
  char ckpttmpfile[512];
  sprintf(ckpttmpfile,"%s%s.scfvec.tmp",ckptdir,fname);
  fprintf(outfile,
          "  GRSCF: checkpointing vector\n");
  StateBinFileOut out(ckpttmpfile);
  vec.save_state(out);
  out.close();
  rename(ckpttmpfile,ckptfile);
}
#endif

#if 0
double
GRSCF::rms_delta_density(const RefSymmSCMatrix& DP,
                         const RefSymmSCMatrix& DP_open)
{
  RefSCElementScalarProduct op(new SCElementScalarProduct);
  DP->element_op(op,DP);
  double delta = op->result;
  if (DP_open.nonnull()) {
    op->init();
    DP_open->element_op(op,DP_open);
    delta += op->result();
  }
  int ntri = DP.dim().n() * (1 + DP.dim().n()) / 2;
  return delta/ntri;
}
#endif

void
GRSCF::print(SCostream&o)
{
  OneBodyWavefunction::print(o);
}

