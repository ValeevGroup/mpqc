
/* mpqc_int.cc -- functions to allow use of libGeom with old mpqc code
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

#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <util/unix/cct_cprot.h>
#include <util/misc/libmisc.h>

#include <math/optimize/opt.h>
#include <math/scmat/local.h>

#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/localdef.h>

#include "mpqc_int.h"

#define STATEOUT StateOutBinXDR
#define STATEIN StateInBinXDR
// #define STATEOUT StateOutBin
// #define STATEIN StateInBin
// #define STATEOUT StateOutText
// #define STATEIN StateInText

/////////////////////////////////////////////////////////////////

static FILE *outfp=stdout;
static FILE *errfp=stderr;

static RefMolecule mol;
static RefMolecularCoor coor;
static RefHessianUpdate update;

static RefSymmSCMatrix hessian;
static RefSCVector xn;
static RefSCVector gn;
static RefSCVector cart_grad;

RefLocalSCDimension di, dc;

static int iter=1;

static int print_hessian=0;
static int print_internal=0;
static int cartesians=0;
static double conv_crit=1.0e-6;
static double conv_rmsf=1.0e-5;
static double conv_maxf=4.0e-5;
static double maxstepsize=0.3;

static double cart_tol=1.0e-10;

static double maxforce=0;
static double rmsforce=0;
static double maxdisp=0;
static double rmsdisp=0;


/////////////////////////////////////////////////////////////////

static void
get_input(const RefKeyVal& keyval)
{
  if (keyval->exists("rms_force")) {
    int conv = keyval->intvalue("rms_force");
    conv_rmsf = pow(10.0,(double)-conv);
  }

  if (keyval->exists("max_force")) {
    int conv = keyval->intvalue("max_force");
    conv_maxf = pow(10.0,(double)-conv);
  }

  if (keyval->exists("convergence")) {
    int conv = keyval->intvalue("convergence");
    conv_crit = pow(10.0,(double)-conv);
  }

  if (keyval->exists("maxstepsize")) {
    maxstepsize = keyval->doublevalue("maxstepsize");
  }

  if (keyval->exists("cartesians")) {
    cartesians = keyval->booleanvalue("cartesians");
  }

  if (keyval->exists("print_hessian")) {
    print_hessian = keyval->booleanvalue("print_hessian");
  }

  if (keyval->exists("print_internal")) {
    print_internal = keyval->booleanvalue("print_internal");
  }

  fprintf(outfp,"  intco:print_hessian       = %d\n",print_hessian);
  fprintf(outfp,"  intco:print_internal      = %d\n",print_internal);
  fprintf(outfp,"  intco:cartesians          = %d\n",cartesians);
  fprintf(outfp,"  intco:maxstepsize         = %g\n",maxstepsize);
  fprintf(outfp,"  intco:cartesian_tolerance = %g\n",cart_tol);
  fprintf(outfp,"  intco:convergence         = %g\n",conv_crit);
  fprintf(outfp,"  intco:max_force           = %g\n",conv_maxf);
  fprintf(outfp,"  intco:rms_force           = %g\n",conv_rmsf);
  fprintf(outfp,"\n");
}

////////////////////////////////////////////////////////////////////////////

int
Geom_init_mpqc(RefMolecule& molecule, const RefKeyVal& topkeyval)
{
  int i;

 // create a new keyval which adds :intco to the search list 
  RefKeyVal keyval = new PrefixKeyVal(":intco :default",topkeyval);
  
  get_input(keyval);

 // set mol = molecule
  mol = molecule;

 // if geom.dat exists, then read in all of the hooey stored there,
 // otherwise read intco stuff from mpqc.in

  struct stat buf;

  if (stat("geom.dat",&buf) < 0 || buf.st_size==0) {
    STATEOUT so("geom.dat","w+");

   // read coor and update from the input
    coor = keyval->describedclassvalue("coor");
    if (coor.null())
      coor = new IntMolecularCoor(mol);
    
   // for now, always read the update from the input. default to bfgs
    update = keyval->describedclassvalue("update");
    if (update.null())
      update = new BFGSUpdate;

    dc = mol->dim_natom3();
    di = coor->dim();

   // create the inverted hessian 
    hessian = new LocalSymmSCMatrix(di.pointer());
    coor->guess_hessian(hessian);
    hessian = hessian.i();

    xn = new LocalSCVector(di.pointer());
    gn = new LocalSCVector(di.pointer());
    cart_grad = new LocalSCVector(dc.pointer());

    coor->to_internal(xn);
    gn->assign(0.0);
    cart_grad->assign(0.0);
    
   // save it all to disk
    so.put(iter);
    mol.save_state(so);
    coor.save_state(so);
    update.save_state(so);
    dc.save_state(so);
    di.save_state(so);
    xn.save_state(so);
    gn.save_state(so);
    cart_grad.save_state(so);
    hessian.save_state(so);
  }  else {
    STATEIN si("geom.dat","r+");

    si.get(iter);
    mol.restore_state(si);
    coor.restore_state(si);
    update.restore_state(si);
    dc.restore_state(si);
    di.restore_state(si);
    xn.restore_state(si);
    gn.restore_state(si);
    cart_grad.restore_state(si);
    hessian.restore_state(si);

   // make sure molecule and mol refer to the same object
    molecule = mol;
    
    fprintf(outfp,
            "\n restarting geometry optimization at iteration %d\n",iter);
  }

  fprintf(outfp,"\n Initial geometry in Geom_init_mpqc\n");
  fprintf(outfp,"Molecule:\n");
  fflush(outfp);
  mol->print();
  
  fprintf(outfp,"\n Initial simple internal coordinates\n\n");
  fflush(outfp);
  coor->print_simples();
  fprintf(outfp,"\n");
  fflush(outfp);

  return GEOM_COMPUTE_GRADIENT;
}

void
Geom_done_mpqc(const RefKeyVal& keyval, int converged)
{
  const char *msg;

  if (converged)
    fprintf(outfp,"\nConverged Simple Internal Coordinates\n");
  else
    fprintf(outfp,"\nNonconverged Simple Internal Coordinates\n");

  coor->print_simples();
  fprintf(outfp,"\n");

  fprintf(outfp,"Final cartesian coordinates\n\nMolecule:\n");
  mol->print();
  fprintf(outfp,"\n");

  coor=0;
  mol=0;
}


int
Geom_update_mpqc(double_matrix_t *grad, const RefKeyVal& keyval)
{
  int i,j,ij;
  
  // find rms and max force
  rmsforce=0;
  maxforce=0;
  for (j=ij=0; j < grad->n2; j++) {
    for (i=0; i < grad->n1; i++,ij++) {
      double d = fabs(grad->d[i][j]);
      rmsforce += d*d;
      maxforce = (d>maxforce) ? d : maxforce;
      cart_grad->set_element(ij,grad->d[i][j]);
      }
    }
  rmsforce = sqrt(rmsforce/(grad->n1*grad->n2));

  // transform cartesian coordinates to internal coordinates 
  coor->to_internal(xn);
  xn.print("internal coordinates");

  // transform cartesian gradient to internal coordinates
  coor->to_internal(gn,cart_grad);
  gn.print("internal coordinate gradients");

  // update the inverse hessian
  RefNLP2 nlp = 0;
  update->update(hessian,nlp,xn,gn);

  // take the Newton-Raphson step
  RefSCVector xdisp = -1.0*(hessian * gn);
  xdisp.print("internal coordinate displacements");

  // calculate max and rms stepsize
  maxdisp = xdisp.maxabs();
  rmsdisp = 0;
  for (i=0; i < xdisp->n(); i++)
    rmsdisp += xdisp(i) * xdisp(i);
  rmsdisp = sqrt(rmsdisp/xdisp->n());

  // pulay's convergance criterion
  //
  // this should probably be done with element ops, but I'm too lazy to
  // do that now
  double iconv=0;
  for (i=0; i < xdisp->n(); i++) {
    double xg = fabs(xdisp(i)*gn(i));
    iconv = (xg>iconv) ? xg : iconv;
  }
  iconv /= 2.0;
  
  if (iconv < conv_crit && rmsforce < conv_rmsf && maxforce < conv_maxf) {
    fprintf(outfp,"\n max of 1/2 idisp*iforce = %15.10g\n",iconv);
    fprintf(outfp," max force               = %15.10g\n",maxforce);
    fprintf(outfp," rms force               = %15.10g\n",rmsforce);
    fprintf(outfp,"\n the geometry is converged\n");

    fprintf(outfp,"\n converged geometry\n");
    mol->print();

    return GEOM_DONE;
  }

 // scale the displacement vector if it's too large
  double tot = xdisp.scalar_product(xdisp);
  if (tot > maxstepsize) {
    double scal = maxstepsize/tot;
    fprintf(outfp,"\n stepsize of %f is too big, scaling by %f\n",tot,scal);
    xdisp.scale(scal);
    tot *= scal;
  }
  fprintf(outfp,"\n taking step of size %f\n",tot);

  // displace internal coordinates
  xn.accumulate(xdisp);
  xn.print("new internal coordinates");
  
  fprintf(outfp,"\n max of 1/2 idisp*iforce = %15.10g\n",iconv);
  fprintf(outfp," max force               = %15.10g\n",maxforce);
  fprintf(outfp," rms force               = %15.10g\n",rmsforce);
  fprintf(outfp," max disp                = %15.10g\n",maxdisp);
  fprintf(outfp," rms disp                = %15.10g\n",rmsdisp);

  // now transform new internal coords back to cartesian coordinates
  Molecule foo_save = *(mol.pointer());
  if (coor->to_cartesian(xn) < 0) {
    // that failed, try steepest descent and get out of here
    fprintf(stderr,"\n  trying cartesian steepest descent..."
                   "cross your fingers\n");
    for (int i=0; i < cart_grad.n(); i++)
      (*mol.pointer())[i/3][i%3] = foo_save[i/3][i%3]-cart_grad[i];
  }

  printf("\nnew molecular coordinates\n");
  mol->print();

  // checkpoint
  STATEOUT so("geom.dat","w+");
  iter++;
  so.put(iter);
  mol.save_state(so);
  coor.save_state(so);
  update.save_state(so);
  dc.save_state(so);
  di.save_state(so);
  xn.save_state(so);
  gn.save_state(so);
  cart_grad.save_state(so);
  hessian.save_state(so);

  return GEOM_COMPUTE_GRADIENT;
}

void
Geom_write_pdb(const RefKeyVal& keyval, RefMolecule& mol, char *title)
{
  double bohr = 0.52917706;

  FILE *pdbf;

  if (keyval->exists("filename")) {
    char path[512];
    char *fn = keyval->pcharvalue("filename");
    sprintf(path,"%s.pdb",fn);
    pdbf = fopen(path,"a+");
    delete[] fn;
  }  else
    pdbf = fopen("mpqc.pdb","a+");

  if (!pdbf) return;

  if (keyval->exists("title")) {
    char *title = keyval->pcharvalue("title");
    fprintf(pdbf,"%-10s%-60s\n","COMPND",title);
    delete[] title;
  }  else {
    fprintf(pdbf,"%-10s%-60s\n","COMPND","Title");
  }

  if (title)
    fprintf(pdbf,"REMARK   %s\n",title);

  for (int i=0; i < mol->natom(); i++) {
    char symb[4];
    sprintf(symb,"%s1",(*mol)[i].element().symbol());

    fprintf(pdbf,
            "HETATM%5d  %-3s UNK %5d    %8.3f%8.3f%8.3f  0.00  0.00   0\n",
            i+1, symb, 0, (*mol)[i][0]*bohr, (*mol)[i][1]*bohr,
            (*mol)[i][2]*bohr);
  }

  for (i=0; i < mol->natom(); i++) {
    double at_rad_i = (*mol)[i].element().atomic_radius();

    fprintf(pdbf,"CONECT%5d",i+1);

    for (int j=0; j < mol->natom(); j++) {

      if (j==i) continue;

      double at_rad_j = (*mol)[j].element().atomic_radius();

      if (dist((*mol)[i].point(),(*mol)[j].point()) < 1.1*(at_rad_i+at_rad_j))
        fprintf(pdbf,"%5d",j+1);
    }

    fprintf(pdbf,"\n");
  }

  fprintf(pdbf,"END\n");
  fclose(pdbf);
}
