
#ifdef __GNUC__
#pragma implementation
#endif

#include "mpqc.h"

extern "C" {

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
#include <util/misc/libmisc.h>

#if defined(I860) || defined(SGI)
void bzero(void*,int);
#endif
};

#include <util/keyval/keyval.h>
#include <util/sgen/sgen.h>
#include <util/misc/libmisc.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/force/libforce.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>

////////////////////////////////////////////////////////////////////////////

static void init_dmt(centers_t*,scf_struct_t*,sym_struct_t*);
static void clean_and_exit(int);
static void mkcostvec(centers_t*, sym_struct_t*, dmt_cost_t*);

///////////////////////////////////////////////////////////////////////////
int host;
char* argv0;
int MPSCF::active = 0;

#define CLASSNAME MPSCF
#define PARENTS public OneBodyWavefunction
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
MPSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

MPSCF::MPSCF(KeyVal&keyval):
  OneBodyWavefunction(keyval),
  _exchange_energy(this),
  _eigenvectors(this),
  _scf(this)
{
  RefGaussianBasisSet gbs = keyval.describedclassvalue("basis");
  centers_t *tcenters = gbs->convert_to_centers_t(_mol.pointer());

  if (!tcenters) {
    exit(3);
  }

  init_centers(&centers);
  init_centers(&oldcenters);
  assign_centers(&centers,tcenters);

  free_centers(tcenters);

  if (sym_struct_from_pg(_mol->point_group(),centers,sym_info) < 0) {
    fprintf(stderr,"MPSCF::MPSCF(KeyVal&):  could not form sym_info\n");
    exit(1);
  }

  if (scf_init_scf_struct(keyval,centers,scf_info) < 0) {
    fprintf(stderr,"MPSCF::MPSCF(KeyVal&):  could not form scf_info\n");
    exit(1);
  }

  _ndocc = scf_info.nclosed;
  _nsocc = scf_info.nopen;

 // override the default thresholds
  if (!keyval.exists("value_accuracy")) {
    set_desired_value_accuracy(1.0e-10);
  }
  if (!keyval.exists("gradient_accuracy")) {
    set_desired_gradient_accuracy(1.0e-9);
  }
  if (!keyval.exists("hessian_accuracy")) {
    set_desired_hessian_accuracy(1.0e-8);
  }

  // make sure only one MPSCF object exists
  if (active) {
    fprintf(stderr,"only one MPSCF allowed\n");
    abort();
  }
  active = 1;
  
 // some PICL emulators need argv0
  argv0 = "bad news: argv0 not available -- MPSCF";

 // get some info from the host program
  int nproc,me,top,ord,dir;

  open0(&nproc,&me,&host);
  setarc0(&nproc,&top,&ord,&dir);

 // initialize timing for mpqc 
  tim_enter("mpqcnode");
  tim_enter("input");

  outfile = stdout;

  int throttle,sync_loop;

  if (me == 0) {

    fprintf(outfile,
              "\n      MPSCF: Massively Parallel Quantum Chemistry\n\n\n");
    fprintf(outfile,"  Running on a %s with %d nodes.\n",machine_type(),nproc);
    fflush(outfile);

   // read input, and initialize various structs

    node_timings = keyval.booleanvalue("node_timings");

    throttle = keyval.intvalue("throttle");

    sync_loop = keyval.intvalue("sync_loop");
    if (keyval.error() != KeyVal::OK) sync_loop = 1;

    save_vector = keyval.booleanvalue("save_vector");
    if (keyval.error() != KeyVal::OK) save_vector=1;

    fprintf(outfile,"\n  mpqc options:\n");
    fprintf(outfile,"    node_timings       = %s\n",(node_timings)?"YES":"NO");
    fprintf(outfile,"    throttle           = %d\n",throttle);
    fprintf(outfile,"    sync_loop          = %d\n",sync_loop);
    fprintf(outfile,"    save_vector        = %s\n",(save_vector)?"YES":"NO");

    sprintf(vecfile,"%s.scfvec",scf_info.fname);

   // pretty print the scf struct
    scf_print_options(outfile,scf_info);
  }

  sgen_reset_bcast0();

  bcast0_scf_struct(&scf_info,0,0);
  bcast0_sym_struct(&sym_info,0,0);
  bcast0_centers(&centers,0,0);

  bcast0(&save_vector,sizeof(int),mtype_get(),0);
  bcast0(&throttle,sizeof(int),mtype_get(),0);
  bcast0(&sync_loop,sizeof(int),mtype_get(),0);
  bcast0(&node_timings,sizeof(int),mtype_get(),0);

  int size;
  if (me==0) size=strlen(vecfile)+1;
  bcast0(&size,sizeof(int),mtype_get(),0);
  bcast0(vecfile,size,mtype_get(),0);

 // if we're using a projected guess vector, then initialize oldcenters
  if (scf_info.proj_vector) {
    if (me==0) {
      RefGaussianBasisSet gbs = keyval.describedclassvalue("pbasis");
      tcenters = gbs->convert_to_centers_t(_mol.pointer());

      assign_centers(&oldcenters,tcenters);
      free_centers(tcenters);

      int_normalize_centers(&oldcenters);
    }

    bcast0_centers(&oldcenters,0,0);
  }

 // initialize the dmt routines
  init_dmt(&centers,&scf_info,&sym_info);

 // initialize force and geometry routines
  if (me==0) fprintf(outfile,"\n");
  if (scf_info.iopen)
    dmt_force_osscf_keyval_init(&keyval,outfile);
  else
    dmt_force_csscf_keyval_init(&keyval,outfile);
  if (me==0) fprintf(outfile,"\n");

 // set the throttle for libdmt loops
  dmt_set_throttle(throttle);

 // set the sync_loop for libdmt loops
  dmt_set_sync_loop(sync_loop);

 // allocate memory for vector and fock matrices

  Scf_Vec = dmt_create("scf vector",scf_info.nbfao,COLUMNS);
  Fock = dmt_create("fock matrix",scf_info.nbfao,SCATTERED);
  if (scf_info.iopen)
    FockO = dmt_create("open fock matrix",scf_info.nbfao,SCATTERED);
  else
    FockO = dmt_nil();

  if (mynode0() == 0 && save_vector)
      fprintf(outfile,"  scf vector will be written to file %s\n\n",vecfile);

 // if restart, then read in old scf vector if it exists

  FILE* test_vec = fopen(vecfile,"r");
  if (test_vec && scf_info.restart) {
    fclose(test_vec);
    dmt_read(vecfile,Scf_Vec);
    if (me==0) fprintf(outfile,"\n  read vector from file %s\n\n",vecfile);
    scf_info.restart=1;
  } else if (test_vec) {
    fclose(test_vec);
  }

  tim_exit("input");
}

MPSCF::~MPSCF()
{
  if (scf_info.iopen)
    dmt_force_osscf_done();
  else
    dmt_force_csscf_done();

  tim_print(node_timings);
  clean_and_exit(host0());
  active = 0;
}

MPSCF::MPSCF(StateIn&s):
  SavableState(s,MPSCF::class_desc_),
  OneBodyWavefunction(s),
  _exchange_energy(this),
  _eigenvectors(this),
  _scf(this)
{
  // make sure only one MPSCF object exists
  if (active) {
      fprintf(stderr,"only one MPSCF allowed\n");
      abort();
    }
  active = 1;

  abort();
}

void
MPSCF::save_data_state(StateOut&s)
{
  OneBodyWavefunction::save_data_state(s);
  abort();
}


int
MPSCF::do_exchange_energy(int f)
{
  int old = _exchange_energy.compute();
  _exchange_energy.compute() = f;
  return old;
}

int
MPSCF::do_eigenvectors(int f)
{
  int old = _eigenvectors.compute();
  _eigenvectors.compute() = f;
  return old;
}

void
MPSCF::compute()
{
  int i,j;

  // Adjust the value accuracy if gradients are needed and set up
  // minimal accuracies.
  if (desired_value_accuracy() > 1.0e-4)
    set_desired_value_accuracy(1.0e-4);

  if (_gradient.compute()) {
    if (desired_gradient_accuracy() > 1.0e-4)
      set_desired_gradient_accuracy(1.0e-4);
    if (desired_value_accuracy() > 0.1 * desired_gradient_accuracy()) {
      set_desired_value_accuracy(0.1 * desired_gradient_accuracy());
    }
  }

  if (mynode0()==0) fprintf(outfile,"\n");

 // copy geometry from Molecule to centers
  for (i=0; i < centers.n; i++) {
    for (j=0; j < 3; j++) {
      centers.center[i].r[j] = _mol->operator[](i)[j];
    }
  }

 // broadcast new geometry information
  for (i=0; i < centers.n; i++) {
    bcast0(centers.center[i].r,sizeof(double)*3,mtype_get(),0);
  }

  // make sure everybody knows if we're to compute a vector
  bcast0(&_scf.compute(),sizeof(int),mtype_get(),0);
  bcast0(&_scf.computed(),sizeof(int),mtype_get(),0);
  bcast0(&_energy.compute(),sizeof(int),mtype_get(),0);
  bcast0(&_energy.computed(),sizeof(int),mtype_get(),0);

  // calculate new scf_vector
  if (_scf.needed() || _energy.needed()) {
    tim_enter("scf_vect");

    scf_info.convergence = ((int) -log10(_energy.desired_accuracy())) + 1;
    scf_info.intcut = scf_info.convergence + 1;
    fprintf(outfile,"\n  MPSCF: computing energy with integral cutoff 10^-%d"
              " and convergence 10^-%d\n",
              scf_info.intcut,scf_info.convergence);

    if (scf_vector(&scf_info,&sym_info,&centers,
                     Fock,FockO,Scf_Vec,&oldcenters,outfile) < 0) {
      fprintf(stderr,"MPSCF::compute():  trouble forming scf vector\n");
      clean_and_exit(host);
    }

    if (save_vector) dmt_write(vecfile,Scf_Vec);

    scf_info.restart=1;
    _scf.computed() = 1;
    set_energy(scf_info.nuc_rep+scf_info.e_elec);
    _energy.set_actual_accuracy(_energy.desired_accuracy());

    tim_exit("scf_vect");
  }

  // make sure that everybody knows whether or not a gradient is needed.
  bcast0(&_gradient.compute(),sizeof(int),mtype_get(),0);
  bcast0(&_gradient.computed(),sizeof(int),mtype_get(),0);

  // compute the gradient if needed
  double_matrix_t grad;
  allocbn_double_matrix(&grad,"n1 n2",3,centers.n);
  if (_gradient.needed()) {
    int cutoff = ((int) -log10(desired_gradient_accuracy())) + 1;
    fprintf(outfile,"\n");
    fprintf(outfile,"  MPSCF: computing gradient with cutoff 10^-%d\n",cutoff);
    if (!scf_info.iopen) {
      dmt_force_csscf_threshold_10(cutoff);
      dmt_force_csscf(outfile,Fock,Scf_Vec,
                          &centers,&sym_info,scf_info.nclosed,&grad);
    } else {
      int ndoc=scf_info.nclosed;
      int nsoc=scf_info.nopen;
      dmt_force_osscf(outfile,Fock,FockO,Scf_Vec,
                      &centers,&sym_info,ndoc,nsoc,&grad);
    }

    // convert the gradient to a SCVector
    RefSCVector g(_moldim);
    for (int ii=0,i=0; i<centers.n; i++) {
      for (int j=0; j<3; j++,ii++) {
        g(ii) = grad.d[j][i];
      }
    }

    // update the gradient, converting to internal coordinates if needed
    set_gradient(g);
    _gradient.set_actual_accuracy(desired_gradient_accuracy());
  }

 // compute the hessian if needed
  if (_hessian.needed()) {
    if (mynode0()==0) {
      fprintf(stderr,"A hessian was requested, but cannot be computed\n");
      abort();
    }
  }
}

double
MPSCF::exchange_energy()
{
  return _exchange_energy;
}

RefSCMatrix
MPSCF::eigenvectors()
{
  return _eigenvectors;
}

double
MPSCF::occupation(int i)
{
  if (i < _ndocc) return 2.0;
  else if (i < _ndocc+_nsocc) return 1.0;
  else return 0.0;
}

void
MPSCF::print(SCostream&o)
{
  if (outfile) fflush(outfile);
  o.flush();
  OneBodyWavefunction::print(o);
  o.flush();
}

///////////////////////////////////////////////////////////////////////////
//
//  static functions used to initialize mpscf
//

static void
clean_and_exit(int host)
{
  int junk;

#if !defined(I860)
  picl_prober();

  /* tell host that we're done */
  if(mynode0()==0 && host0()) send0(&junk,sizeof(int),1,host);
#endif

 /* this one too */
  close0(0);
  exit(0);
}

static void
mkcostvec(centers_t *centers,sym_struct_t *sym_info,dmt_cost_t *costvec)
{
  int flags;
  int i,j,ij;
  double *intbuf;
  extern signed char *scf_bnd_Qvec;

 /* free these up for now */
  int_done_offsets1(centers,centers);
  int_done_1e();

  int_initialize_offsets2(centers,centers,centers,centers);

  flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;

  intbuf =
    int_initialize_erep(flags,0,centers,centers,centers,centers);

  scf_init_bounds(centers,intbuf);

  for (i=ij=0; i<centers->nshell; i++) {
    int nconi = centers->center[centers->center_num[i]].
                       basis.shell[centers->shell_num[i]].ncon;
    int nprimi = centers->center[centers->center_num[i]].
                       basis.shell[centers->shell_num[i]].nprim;
    int ami = centers->center[centers->center_num[i]].
                       basis.shell[centers->shell_num[i]].type[0].am;
    for (j=0; j<=i; j++,ij++) {
      int nconj = centers->center[centers->center_num[j]].
                         basis.shell[centers->shell_num[j]].ncon;
      int nprimj = centers->center[centers->center_num[j]].
                         basis.shell[centers->shell_num[j]].nprim;
      int amj = centers->center[centers->center_num[j]].
                         basis.shell[centers->shell_num[j]].type[0].am;

      costvec[ij].i = i;
      costvec[ij].j = j;
      costvec[ij].magnitude = (int) scf_bnd_Qvec[ij];
      costvec[ij].ami = ami;
      costvec[ij].amj = amj;
      costvec[ij].nconi = nconi;
      costvec[ij].nconj = nconj;
      costvec[ij].nprimi = nprimi;
      costvec[ij].nprimj = nprimj;
      costvec[ij].dimi = INT_SH_NFUNC(centers,i);;
      costvec[ij].dimj = INT_SH_NFUNC(centers,j);;
      }
    }

  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);
  scf_done_bounds();
}

static void
init_dmt(centers_t *centers, scf_struct_t *scf_info, sym_struct_t *sym_info)
{
 // initialize the centers
  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);

 // set up shell map for libdmt
  int nshell=centers->nshell;
  int nshtr=nshell*(nshell+1)/2;

  int *shellmap = new int[nshell];
  bzero(shellmap,sizeof(int)*nshell);

  for (int i=0; i < centers->nshell ; i++)
    shellmap[i] = INT_SH_NFUNC(centers,i);

  dmt_cost_t *costvec = (dmt_cost_t*) malloc(sizeof(dmt_cost_t)*nshtr);

  if (!costvec) {
    dmt_def_map2(scf_info->nbfao,centers->nshell,shellmap,costvec,1);
  } else {
    mkcostvec(centers,sym_info,costvec);
    dmt_def_map2(scf_info->nbfao,centers->nshell,shellmap,costvec,0);
    free(costvec);
  }

  free(shellmap);
  if (mynode0()==0) printf("\n");
  dmt_map_examine();

  int_done_offsets1(centers,centers);
  int_done_1e();
}
