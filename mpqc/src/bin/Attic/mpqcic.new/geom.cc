
/* geom.cc -- geometry optimization code
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      March, 1993
 */

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/symm.h>
#include <chemistry/molecule/simpleQCList.h>
#include <chemistry/molecule/symmQCList.h>
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/lmath.h>

#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>

#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "input.h"
#include "mpqcstate.h"

#define STATEOUT StateOutBinXDR
#define STATEIN StateInBinXDR
//#define STATEOUT StateOutBin
//#define STATEIN StateInBin
//#define STATEOUT StateOutText
//#define STATEIN StateInText

/////////////////////////////////////////////////////////////////

static FILE *outfp=stdout;
static FILE *errfp=stderr;

static centers_t *centers;
static Molecule mol;

static RefSymmCoList symm_coords;

static int iter=1;
static int justa1=1;

static double maxforce=0;
static double rmsforce=0;
static double maxdisp=0;
static double rmsdisp=0;

static int init=0;
static DMatrix oldhessian;
static DVector oldforce;
static DVector oldcoord;
static DVector olddeltax;

static int internal_forces(DMatrix&,Molecule&,double_matrix_t*);
static int update_hessian();
static int cartesian_disp(DMatrix&);
static void gdiis_disp();
static void guess_hessian(SimpleCoList*);
static int get_symmco();
static DMatrix gen_inverse(DMatrix&);
static int setup_cart(centers_t *cs);

/////////////////////////////////////////////////////////////////

void
init_geom(MPQC_state& mystate, MPQC_input& mpqcin, MPQC_geom_input& geomin,
          centers_t& cs)
{
#if 0
  int i;
  init=1;

  centers=&cs;

  struct stat buf;

 // well, if restart is 0, or this is the first go around, then initialize
 // all the junk we need
  if (!geomin.restart() || stat("mpqcstate.dat",&buf) < 0 || buf.st_size==0) {

    if (mpqcin.optimize_geometry()) {
     // ok, first let's create the molecule object using the coords in centers
      ParsedKeyVal rawin("mpqc.in");
      PrefixKeyVal pref1(":intco",rawin);
      PrefixKeyVal pref2(":default",rawin);
      AggregateKeyVal keyval(pref1,pref2);

      for (i=0; i < centers->n; i++) {
        char* label = keyval.pcharvalue("atom_labels",i);
        AtomicCenter ac(centers->center[i].atom,
                        centers->center[i].r[0],
                        centers->center[i].r[1],
                        centers->center[i].r[2],
                        label
                        );
        delete[] label;
        mol.add_atom(i,ac);
        }

     // now that we have a geometry, let's make the coordinates
      int ncoord=0;
      if(!geomin.cartesians()) {
        if(get_symmco() < 0) {
          err_msg("init_geom: yikes!  can't make symm coords");
          exit(1);
          }

        for(SymmCoListIter scp=symm_coords; scp; scp++,ncoord++) ;
        }
      else
         ncoord = mol.natom()*3;

      coord.resize(ncoord);
      coord.zero();

      force.resize(ncoord);
      force.zero();

      deltax.resize(ncoord);
      deltax.zero();

      if (!geomin.cartesians()) {
        i=0;
        for(SymmCoListIter scp=symm_coords; scp; scp++,i++) {
          coord[i] = scp->value();
          }
        }
      }

    mystate.checkpoint(mpqcin,geomin);
    }
  else {
    mystate.restore();
    }

  if (mpqcin.optimize_geometry()) geomin.print();

  fprintf(outfp,"\n initial geometry\n");
  mol.print(outfp);
#endif
  }

void
checkpoint_geom(StateOut& so, MPQC_geom_input& geomin)
{
#if 0
  so.put(init);
  so.put(iter);

  if (!init) return;

  int cart=geomin.cartesians();
  so.put(cart);

  if(!cart) symm_coords->save_state(so);

  mol.save_object_state(so);
  coord.save_object_state(so);
  force.save_object_state(so);
  deltax.save_object_state(so);
  hessian.save_object_state(so);
#endif
}

void
restore_geom_state(StateIn& si)
{
#if 0
  si.get(init);
  si.get(iter);

  if (!init) return;

  int cart; si.get(cart);

  if (!cart) symm_coords = SymmCoList::restore_state(si);

  Molecule mol2(si); mol=mol2;
  DVector ic(si);    coord=ic;
  DVector ifc(si);   force=ifc;
  DVector id(si);    deltax=id;
  DMatrix h(si);     hessian=h;

  if (centers && centers->n == mol.natom()) {
    for(int i=0; i < centers->n; i++) {
      centers->center[i].r[0] = mol[i][0];
      centers->center[i].r[1] = mol[i][1];
      centers->center[i].r[2] = mol[i][2];
      }
    }
#endif
}

#if 0
void
Geom_done_mpqc()
{
  symm_coords = 0;

  SimpleCoList * list = Geom_form_simples(mol);

  ParsedKeyVal rawin("mpqc.in");
  PrefixKeyVal pref1(":intco",rawin);
  PrefixKeyVal pref2(":default",rawin);
  AggregateKeyVal keyval(pref1,pref2);

  int nadd;
  if(nadd=keyval.count("add_simp")) {
    for(int i=0; i < nadd; i++) {
      char *val = keyval.pcharvalue("add_simp",i,0);

      if (!strcmp("stre",val))
        list->add(new Stre(&keyval,"add_simp",i));
      else if (!strcmp("bend",val))
        list->add(new Bend(&keyval,"add_simp",i));
      else if (!strcmp("tors",val))
        list->add(new Tors(&keyval,"add_simp",i));
      else if (!strcmp("out",val))
        list->add(new Out(&keyval,"add_simp",i));
      else if (!strcmp("linip",val))
        list->add(new LinIP(&keyval,"add_simp",i));
      else if (!strcmp("linop",val))
        list->add(new LinOP(&keyval,"add_simp",i));
      delete[] val;
      }
    }

  Geom_calc_simples(list,mol);

  fprintf(outfp,"\n  Final internal coordinates\n");
  Geom_print_pretty(list);
  delete list;
}
#endif

static int
get_symmco()
{
#if 0
  static int firsttime=1;

 // first let's screw up the input
  ParsedKeyVal rawin("mpqc.in");
  PrefixKeyVal pref1(":intco",rawin);
  PrefixKeyVal pref2(":default",rawin);
  AggregateKeyVal keyval(pref1,pref2);

  SimpleCoList *list=0;

 // if the simples are defined in the input, read them, otherwise generate
 // them based on the geometry

  if(keyval.count("simp"))
    list = Geom_read_simples(&keyval);
  else
    list = Geom_form_simples(mol);

  int nadd;
  if(nadd=keyval.count("add_simp")) {
    for(int i=0; i < nadd; i++) {
      char *val = keyval.pcharvalue("add_simp",i,0);

      if (!strcmp("stre",val))
        list->add(new Stre(&keyval,"add_simp",i));
      else if (!strcmp("bend",val))
        list->add(new Bend(&keyval,"add_simp",i));
      else if (!strcmp("tors",val))
        list->add(new Tors(&keyval,"add_simp",i));
      else if (!strcmp("out",val))
        list->add(new Out(&keyval,"add_simp",i));
      else if (!strcmp("linip",val))
        list->add(new LinIP(&keyval,"add_simp",i));
      else if (!strcmp("linop",val))
        list->add(new LinOP(&keyval,"add_simp",i));
      delete[] val;
      }
    }

  if(keyval.count("symm")) {
    symm_coords = Geom_read_symm(&keyval,"symm",list);
    if(!redundant) symm_coords = Geom_form_symm(mol,symm_coords,justa1);
    }
  else if(redundant)
    symm_coords = Geom_symm_from_simple(list,1);
  else
    symm_coords = Geom_form_symm(mol,list,justa1);

  if(!symm_coords) {
    err_msg("get_symmco: yikes!  no symmetrized internal coords");
    return -1;
    }

  int csize;
  if ((csize=keyval.count("chessian"))) {
    csize = (int) sqrt((double)csize);
    DMatrix chessian(csize,csize);
    int ij=0;
    for (int i=0; i < csize; i++)
      for (int j=0; j < csize; j++,ij++)
        chessian(i,j) = keyval.doublevalue("chessian",ij);

    chessian *= 4.359813653/(0.52917706*0.52917706);

    DMatrix bmat = Geom_make_bmat(symm_coords,mol);
    DMatrix bmbt = gen_inverse(bmat * bmat.transpose());

    hessian = bmbt * bmat * chessian * bmat.transpose() * bmbt;
    for (i=0; i < hessian.nrow(); i++) {
      for (int j=0; j < i; j++) {
        hessian(i,j) += hessian(j,i);
        hessian(j,i) = hessian(i,j) = 0.5*hessian(i,j);
        }
      }

    bmbt = bmat*bmat.transpose();
    DMatrix p = bmbt * gen_inverse(bmbt);
    hessian = p * gen_inverse(p*hessian*p) * p;
    }
  else
    guess_hessian(list);

  if(list) {
    if(firsttime) {
      fprintf(outfp,"  Simple Internal Coordinates\n");
      Geom_calc_simples(list,mol);
      Geom_print_pretty(list);
      firsttime=0;
      }
    delete list;
    }

  fprintf(outfp,"\n  in get_symmco:  %d symmetrized internal coordinates\n",
                hessian.nrow());

#endif
  return 0;
  }

#if 0
int
Geom_update_mpqc(double energy, double_matrix_t *grad, double_matrix_t *hess)
{

  oldforce = force;
  oldintco = intco;
  olddeltax = deltax;
  oldhessian = hessian;
    
  int i;
  DMatrix bmat;

  // find rms and max force
  rmsforce=0;
  maxforce=0;
  for (i=0; i < grad->n1; i++) {
    for (int j=0; j < grad->n2; j++) {
      double d = fabs(grad->d[i][j]);
      rmsforce += d*d;
      maxforce = (d>maxforce) ? d : maxforce;
      }
    }
  rmsforce = sqrt(rmsforce/(grad->n1*grad->n2));
      
  if (!cartesians) {

   // calculate new internal coordinates and the b matrix
    bmat = Geom_make_bmat(symm_coords,mol);

   // and now the new internal forces
    internal_forces(bmat,mol,grad);

    fprintf(outfp,"\n internal coordinates and forces\n");

    i=0;
    for(SymmCoListIter scp=symm_coords; scp; scp++,i++) {
      intco[i] = scp->value();
      fprintf(outfp," %5d %15.10f %15.10f\n",i+1,intco[i],force[i]);
      }
    }
  else {
    for (int i=0; i < force.dim(); i++) {
      force[i] = grad->d[i%3][i/3]*8.2388575;
      intco[i] = mol[i/3][i%3] * 0.52917706;
      }
    }

  if (update_hess)
    update_hessian();

  if (use_gdiis) {
    gdiis_disp();
    }
  else { // the NR step
    deltax = -1*(hessian * force);
    }

 // scale the displacement vector if too large
  {
     double tot=0,scal=1.0;
     tot = sqrt(deltax.dot(deltax));
     //for(i=0; i < deltax.dim(); i++) tot += fabs(deltax[i]);
     if(tot>maxstepsize) {
       scal=maxstepsize/tot;
       fprintf(outfp,"\n stepsize of %f is too big, scaling by %f\n",tot,scal);
       }
     fprintf(outfp,"\n taking step of size %f\n",tot*scal);
     for(i=0; i < deltax.dim(); i++) deltax[i] *= scal;
     }

  if (cartesians) {
    rmsdisp=0;
    for (i=0; i < deltax.dim(); i++) {
      mol[i/3][i%3] += deltax[i] / 0.52917706;
      rmsdisp += (deltax[i] / 0.52917706)*(deltax[i] / 0.52917706);
      }
    
    rmsdisp = sqrt(rmsdisp/deltax.dim());
    maxdisp = deltax.maxval()/0.52917706;
    }
      
  if (!cartesians && cartesian_disp(bmat) < 0) {
    fprintf(outfp,"\n could not update geometry, saving state and exiting\n");
    //StateOutBinXDR so("geom.dat");
    STATEOUT so("geom.dat");
    so.put(iter);
    symm_coords->save_state(so);
    mol.save_object_state(so);
    intco.save_object_state(so);
    force.save_object_state(so);
    deltax.save_object_state(so);
    hessian.save_object_state(so);
    return GEOM_ABORT;
    }
  else {
    fprintf(outfp,"\n displacements and new internal coordinates\n");
    for(i=0; i < intco.dim(); i++)
      fprintf(outfp," %5d %15.10f %15.10f\n",i+1,deltax[i],intco[i]);

    fprintf(outfp,"\n max of 1/2 deltax*force = %15.10g\n",
			      (deltax * force).maxval()/2);
    fprintf(outfp," max force               = %15.10g\n",maxforce);
    fprintf(outfp," rms force               = %15.10g\n",rmsforce);
    fprintf(outfp," max disp                = %15.10g\n",maxdisp);
    fprintf(outfp," rms disp                = %15.10g\n",rmsdisp);

    fprintf(outfp,"\n updated geometry\n");
    mol.print(outfp);
    }

  //StateOutBinXDR so("geom.dat");
  STATEOUT so("geom.dat");

  iter++;

  so.put(iter);
  symm_coords->save_state(so);
  mol.save_object_state(so);
  intco.save_object_state(so);
  force.save_object_state(so);
  deltax.save_object_state(so);
  hessian.save_object_state(so);

  for(i=0; i < centers->n; i++) {
    centers->center[i].r[0] = mol[i][0];
    centers->center[i].r[1] = mol[i][1];
    centers->center[i].r[2] = mol[i][2];
    }

  if (((deltax * force).maxval()/2) < conv_crit)
    return GEOM_DONE;
  else
    return GEOM_COMPUTE_GRADIENT;
  }


static int
internal_forces(DMatrix& bmat, Molecule& m, double_matrix_t *grad)
{
  int i,j,k;
  int nsym=bmat.nrow();

  DMatrix gminus = gen_inverse(bmat * bmat.transpose());

  DVector gradient(grad->n1*grad->n2);
  double sum=0;

  for(i=0; i < grad->n2; i++) {
    gradient[3*i] =   grad->d[0][i];
    gradient[3*i+1] = grad->d[1][i];
    gradient[3*i+2] = grad->d[2][i];
    }

 // convert gradient from au/bohr to mdyne
  gradient *= 8.2388575;

  force = gminus*bmat*gradient;

  return 0;
  }

static int
update_hessian()
{
  if(iter==1) return 0;

  if (!recalc_hess) {
    DVector delta = force - oldforce;
    DVector v = oldhessian * delta;
    double alpha = 1.0/olddeltax.dot(delta);
    double beta = (1.0 + alpha*delta.dot(v))*alpha;

    hessian = oldhessian + beta*(olddeltax.ccross(olddeltax)) -
             alpha*(olddeltax.ccross(v) + v.ccross(olddeltax));
    }
  else {
    SimpleCoList *list = Geom_form_simples(mol);
    Geom_calc_simples(list,mol);
    guess_hessian(list);
    fprintf(outfp," recalc hessian\n");
    }

  //hessian->inverse().print("\n improved hessian");

  return 0;
  }

static int
cartesian_disp(DMatrix& bmat)
{
  int i;
  int ntry=0;
  int step;
  int nsym=intco.dim();
  int n3 = bmat.ncol();
  SymmCoListIter p;

  DMatrix gminus,nbmat;

  if (!recalc_bmat) {
    gminus = gen_inverse(bmat*bmat.transpose());
    nbmat = bmat;
    }

  Molecule newg;
  DVector ndisp(nsym);
  DVector tintco(nsym);
  DVector tv1(nsym);

Foo:
  tintco = intco + deltax;
  ndisp= deltax;
  newg = mol;
  step=0;

  do {
    if (recalc_bmat) {
      nbmat = Geom_make_bmat(symm_coords,newg);
      gminus = gen_inverse(nbmat*nbmat.transpose());
      }

    DVector tv2 = gminus*ndisp;

    for(i=0; i < n3; i++) {
      double t=0;
      for(int j=0; j < nsym; j++)
        t += nbmat[j][i]*tv2[j];

      newg[i/3][i%3] += t*1.88972666;
      }

    for(p=symm_coords,i=0; p; p++,i++) {
      p->calc_intco(newg);
      tv1[i] = p->value();
      }

    if(tv2.maxval() < 1.0e-10) break;

    step++;

    ndisp = tintco - tv1;
    } while(step < maxstep);

  if(step==maxstep && ntry < 4) {
    fprintf(outfp," too many steps, scaling displacements by 0.5\n");
    deltax *= 0.5;
    ntry++;
    goto Foo;
    }
  else if (step==maxstep) {
    fprintf(outfp," just would not converge\n");
    deltax *= 16.0;
    return -1;
    }

  fprintf(outfp," converged in %d steps\n",step);

 // figure out max and rms cartesian displacements
  maxdisp=0;
  rmsdisp=0;
  for (i=0; i < mol.natom(); i++) {
    double dx = fabs(mol[i][0]-newg[i][0]);
    double dy = fabs(mol[i][1]-newg[i][1]);
    double dz = fabs(mol[i][2]-newg[i][2]);

    rmsdisp += dx*dx + dy*dy + dz*dz;

    maxdisp = (dx > maxdisp) ? dx : maxdisp;
    maxdisp = (dy > maxdisp) ? dy : maxdisp;
    maxdisp = (dz > maxdisp) ? dz : maxdisp;
    }
  rmsdisp = sqrt(rmsdisp/(mol.natom()*3));

  mol=newg;

  deltax = tv1 - intco;
  intco = tv1;

  return 0;
  }

// hack number one, let's see if we can get this to work

static int diis_iter=0;

static void
gdiis_disp()
{
  int i,j,ii,jj;

  if (rmsforce > gdiis_begin) {
    deltax = -1*(hessian * force);
    return;
    }

  diis_iter++;

  int howmany = (diis_iter < ngdiis) ? diis_iter : ngdiis;

  if(diis_iter==1) {
    save_diis.nsave=ngdiis;
    save_diis.coords = new DVector[ngdiis];
    save_diis.grads = new DVector[ngdiis];
    save_diis.error = new DVector[ngdiis];
    }

 // save the current gradient and coordinates
  if(diis_iter <= ngdiis) {
    save_diis.coords[diis_iter-1] = intco;
    save_diis.grads[diis_iter-1] = force;
    }
  else {
    for(i=0; i < ngdiis-1; i++) {
      save_diis.coords[i] = save_diis.coords[i+1];
      save_diis.grads[i] = save_diis.grads[i+1];
      }
    save_diis.coords[ngdiis-1] = intco;
    save_diis.grads[ngdiis-1] = force;
    }

  if (diis_iter==1) {
    deltax = -1*(hessian * force);
    return;
    }

 // now let's form the error matrices
  for(i=0; i < howmany; i++)
    save_diis.error[i] = -1*(hessian * save_diis.grads[i]); 

 // and now the A matrix

  DMatrix A;
  DVector coeff;
  double determ;
  int ntry=0;

  do {
    int num=howmany-ntry;

    A.resize(num+1,num+1);
    coeff.resize(num+1);

    for(ii=0,i=ntry; i < howmany; i++,ii++) {
      coeff[ii]=0;
      for(j=ntry,jj=0; j <= i; j++,jj++)
	A(ii,jj) = A(jj,ii) = save_diis.error[i].dot(save_diis.error[j]);
      }
    A *= 1/A(0,0);

    coeff[num]=1;
    for(i=0; i < num; i++)
      A(num,i) = A(i,num) = 1;
    A(num,num) = 0;
    
    ntry++;
    } while((determ=fabs(A.solve_lin(coeff))) < 1.0e-12);

  DVector delstar(intco.dim());
  DVector xstar(intco.dim());

  delstar.zero();
  xstar.zero();

  for(i=0,ii=ntry-1; ii < howmany; i++,ii++) {
    delstar += save_diis.grads[ii]*coeff[i];
    xstar += save_diis.coords[ii]*coeff[i];
    }

   deltax = xstar - intco - hessian*delstar;
  }

static void
guess_hessian(SimpleCoList *list)
{
  ParsedKeyVal rawin("mpqc.in");
  PrefixKeyVal pref1(":intco",rawin);
  PrefixKeyVal pref2(":default",rawin);
  AggregateKeyVal keyval(pref1,pref2);

 // see if there is a hessian in the input
  int size= keyval.count("hessian");

  if(size==0 || size != hessian.nrow()*(hessian.nrow()+1)/2) {
    hessian = Geom_form_hessian(mol,list,symm_coords);
    }
  else {
    int i,j,ij;
    for(i=ij=0; i < hessian.nrow(); i++)
      for(j=0; j <= i; j++,ij++)
        hessian(i,j) = hessian(j,i) = keyval.doublevalue("hessian",ij);
    }

  DMatrix bmat = Geom_make_bmat(symm_coords,mol);
  DMatrix bmbt = bmat*bmat.transpose();
  DMatrix p = bmbt * gen_inverse(bmbt);
  hessian = p * gen_inverse(p*hessian*p) * p;

  }

static DMatrix
gen_inverse(DMatrix& mat)
{
  DMatrix vecs(mat.nrow(),mat.nrow());
  DVector vals(mat.nrow());

  mat.diagonalize(vals,vecs);

  DMatrix lam(mat.nrow(),mat.nrow());
  lam.zero();

  for(int i=0; i < lam.nrow(); i++)
    if(vals[i] > 1.0e-8) lam(i,i) = 1.0/vals[i];

  lam = vecs * lam * vecs.transpose();

  return lam;
  }


static int
setup_cart(centers_t *cs)
{
  int i,j;

 // ok, first let's create the molecule object using the coords in centers
  ParsedKeyVal rawin("mpqc.in");
  PrefixKeyVal pref1(":intco",rawin);
  PrefixKeyVal pref2(":default",rawin);
  AggregateKeyVal keyval(pref1,pref2);

  for (i=0; i < centers->n; i++) {
    char* label = keyval.pcharvalue("atom_labels",i);
    AtomicCenter ac(centers->center[i].atom,
                    centers->center[i].r[0],
                    centers->center[i].r[1],
                    centers->center[i].r[2],
                    label
                    );
    delete[] label;
    mol.add_atom(i,ac);
    }

 // form simple diagonal hessian, and transform to cartesian coords
  SimpleCoList*  list = Geom_form_simples(mol);

  DMatrix ihessian = Geom_form_hessian(mol,list);
  SymmCoList* slist = Geom_symm_from_simple(list);
  DMatrix bmat = Geom_make_bmat(slist,mol);

  DMatrix thessian = bmat.transpose() * ihessian * bmat;

  delete list;

 // now remove the 3 translations and 3 rotations and invert the hessian
  int dim=thessian.nrow();
  DMatrix evecs(dim,dim);
  DVector evals(dim);

  thessian.diagonalize(evals,evecs);

  hessian.resize(dim,dim);
  hessian.zero();
  
  int ii,jj,kk;
  i=0;
  for (ii=0; ii < dim; ii++) {
    if (evals[ii] < 1e-8 && evals[ii] > -1e-8) {
      for (kk=0; kk < dim; kk++) evecs(kk,ii)=0;
      }
    else {
      double sum=0;
      for (kk=0; kk < dim; kk++) sum += mol[kk/3][kk%3]*evecs(kk,ii);
      if (sum < 1e-8 && sum > -1e-8) {
        for (kk=0; kk < dim; kk++) evecs(kk,ii)=0;
        }
      else
        i++;
      }
    }

  fprintf(outfp,
    "\n %d totally symmetric eigenvectors in cartesian hessian\n",i);

  for (ii=0; ii < dim; ii++) {
    for (jj=0; jj < dim; jj++) {
      for (kk=0; kk < dim; kk++) {
        if (evals[kk] > 1e-10 || evals[kk] < -1e-10)
          hessian(ii,jj) += evecs(ii,kk) * evecs(jj,kk) / evals[kk];
        }
      }
    }

  intco.resize(dim);
  intco.zero();

  force.resize(dim);
  force.zero();

  deltax.resize(dim);
  deltax.zero();

  return GEOM_COMPUTE_GRADIENT;
  }
#endif
