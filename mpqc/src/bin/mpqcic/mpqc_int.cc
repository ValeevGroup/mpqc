
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

#include "mpqc_int.h"

#define STATEOUT StateOutBinXDR
#define STATEIN StateInBinXDR
//#define STATEOUT StateOutBin
//#define STATEIN StateInBin
//#define STATEOUT StateOutText
//#define STATEIN StateInText

/////////////////////////////////////////////////////////////////

static FILE *outfp=stdout;
static FILE *errfp=stderr;

static RefMolecule mol;
static RefMolecularCoor coor;
static RefIHessianUpdate update;

static RefSymmSCMatrix hessian;
static RefSCVector xn, xprev;
static RefSCVector gn, gprev;

RefLocalSCDimension di, dc;

static int iter=1;

static double gdiis_begin=1.0e-2;
static int ngdiis=5;
static int gdiis=1;
static int bfgs=1;
static int redundant=0;
static int justa1=1;
static int update_hess=1;
static int recalc_hess=0;
static int use_gdiis=1;
static int print_hessian=0;
static int recompute_internal=0;
static int print_internal=0;
static double conv_crit=1.0e-6;
static double conv_rmsf=1.0e-5;
static double conv_maxf=4.0e-5;
static int maxstep=50;
static double maxstepsize=0.3;
static int recalc_bmat=0;
static int cartesians=0;

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

  if (keyval->exists("update_hessian")) {
    update_hess = keyval->booleanvalue("update_hessian");
    }

  if (keyval->exists("maxstep")) {
    maxstep = keyval->intvalue("maxstep");
    }

  if (keyval->exists("recalc_bmat")) {
    recalc_bmat = keyval->booleanvalue("recalc_bmat");
    }

  if (keyval->exists("recalc_hessian")) {
    recalc_hess = keyval->booleanvalue("recalc_hessian");
    }

  if (keyval->exists("ngdiis")) {
    ngdiis = keyval->intvalue("ngdiis");
    }

  if (keyval->exists("use_gdiis")) {
    use_gdiis = keyval->booleanvalue("use_gdiis");
    }

  if (keyval->exists("print_hessian")) {
    print_hessian = keyval->booleanvalue("print_hessian");
    }

  if (keyval->exists("recompute_internal")) {
    recompute_internal = keyval->booleanvalue("recompute_internal");
    }

  if (keyval->exists("print_internal")) {
    print_internal = keyval->booleanvalue("print_internal");
    }

  if (keyval->exists("gdiis_begin")) {
    int tmp = keyval->intvalue("gdiis_begin");
    gdiis_begin = pow(10.0,(double)-tmp);
    }

  if (keyval->exists("cartesian_tolerance")) {
    int tmp = keyval->intvalue("cartesian_tolerance");
    cart_tol = pow(10.0,(double)-tmp);
    }

  redundant = keyval->booleanvalue("redundant");
  if (redundant > 0) justa1=0;

  fprintf(outfp,"\n  intco:use_gdiis           = %d\n",use_gdiis);
  fprintf(outfp,"  intco:recompute_internal  = %d\n",recompute_internal);
  fprintf(outfp,"  intco:print_hessian       = %d\n",print_hessian);
  fprintf(outfp,"  intco:print_internal      = %d\n",print_internal);
  fprintf(outfp,"  intco:redundant           = %d\n",redundant);
  fprintf(outfp,"  intco:cartesians          = %d\n",cartesians);
  fprintf(outfp,"  intco:update_hessian      = %d\n",update_hess);
  fprintf(outfp,"  intco:recalc_hessian      = %d\n",recalc_hess);
  fprintf(outfp,"  intco:recalc_bmat         = %d\n",recalc_bmat);
  fprintf(outfp,"  intco:ngdiis              = %d\n",ngdiis);
  fprintf(outfp,"  intco:maxstep             = %d\n",maxstep);
  fprintf(outfp,"  intco:maxstepsize         = %g\n",maxstepsize);
  fprintf(outfp,"  intco:cartesian_tolerance = %g\n",cart_tol);
  fprintf(outfp,"  intco:convergence         = %g\n",conv_crit);
  fprintf(outfp,"  intco:max_force           = %g\n",conv_maxf);
  fprintf(outfp,"  intco:rms_force           = %g\n",conv_rmsf);
  fprintf(outfp,"  intco:gdiis_begin         = %g\n\n",gdiis_begin);
}

////////////////////////////////////////////////////////////////////////////

int
Geom_init_mpqc(RefMolecule& molecule, const RefKeyVal& keyval)
{
  int i;

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
    update = keyval->describedclassvalue("update");
    
    dc = mol->dim_natom3();
    di = coor->dim();

    hessian = new LocalSymmSCMatrix(di.pointer());
    coor->guess_hessian(hessian);

    xn = new LocalSCVector(di.pointer());
    gn = new LocalSCVector(di.pointer());

    coor->to_internal(xn);
    gn->assign(0.0);
    
   // save it all to disk
    so.put(iter);
    mol.save_state(so);
    coor.save_state(so);
  }  else {
    STATEIN si("geom.dat","r+");

    si.get(iter);
    mol.restore_state(si); 
    coor.restore_state(si);

   // make sure molecule and mol refer to the same object
    molecule = mol;
    
    fprintf(outfp,
            "\n restarting geometry optimization at iteration %d\n",iter);
  }

  fprintf(outfp,"\n Initial geometry in Geom_init_mpqc\n");
  fprintf(outfp,"Molecule:\n");
  mol->print();
  
  fprintf(outfp,"\n Initial simple internal coordinates\n\n");
  coor->print_simples();
  fprintf(outfp,"\n");

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
  
  RefSCVector cgrad = new LocalSCVector(dc.pointer());
  
  // find rms and max force
  rmsforce=0;
  maxforce=0;
  for (j=ij=0; j < grad->n2; j++) {
    for (i=0; i < grad->n1; i++,ij++) {
      double d = fabs(grad->d[i][j]);
      rmsforce += d*d;
      maxforce = (d>maxforce) ? d : maxforce;
      cgrad->set_element(ij,grad->d[i][j]);
      }
    }
  rmsforce = sqrt(rmsforce/(grad->n1*grad->n2));

  cgrad.print("cgrad");
  coor->to_internal(gn,cgrad);

  gn.print("gn");
  exit(0);
#if 0  
  int nfixed = (fixed?fixed->length():0);

  oldiforce = iforce;
  oldintco = intco;
  oldidisp = idisp;
  oldhessian = hessian;
    
  int i;
  DMatrix bmat;

      
  if (!cartesians) {

   // calculate new internal coordinates and the b matrix
    bmat = Geom_make_bmat(symm_coords,mol);

   // and now the new internal forces
    internal_forces(bmat,mol,grad);

    fprintf(outfp,"\n internal coordinates and forces\n");

    i=0;
    for(SymmCoListIter scp=symm_coords; scp; scp++,i++) {
      intco[i] = scp->value();
      fprintf(outfp," %5d %15.10f %15.10f\n",i+1,intco[i],iforce[i]);
      }
    }
  else {
    for (int i=0; i < iforce.dim(); i++) {
      iforce[i] = grad->d[i%3][i/3]*8.2388575;
      intco[i] = mol[i/3][i%3] * 0.52917706;
      }
    }

  // if there are fixed coordinates, set their force to zero
  for (i=0; i<nfixed; i++) {
      iforce[i] = 0.0;
      // use a nonzero diagonal
      hessian(i,i) = 1.0;
      oldhessian(i,i) = 1.0;
      // and zero the coupling terms
      for (int j=i+1;j<hessian.nrow(); j++) {
          hessian(i,j) = 0.0;
          hessian(j,i) = 0.0;
          oldhessian(i,j) = 0.0;
          oldhessian(j,i) = 0.0;
        }
    }

  if (update_hess)
    update_hessian(keyval);

  if (print_hessian) {
      printf("The hessian is (after udpate):\n");
      hessian.print();
    }

  if (use_gdiis) {
    gdiis_disp();
    }
  else { // the NR step
    idisp = -1*(hessian * iforce);
    }

  // if there are fixed coordinates, set their displacement to zero
  for (i=0; i<nfixed; i++) {
      idisp[i] = 0.0;
    }

  // Compute the rms and max force in internal coordinates if
  // there are fixed coordinates.  Otherwise rms and max force
  // correspond to cartesian coords.  Really the forces with
  // the fixed coordinate forces set to zero should be backtransformed
  // to cartesian coordinates to compute rms and max force.
  if (fixed) {
      double imaxforce_all = 0.0;
      double irmsforce_all = 0.0;
      double imaxforce_nonfixed = 0.0;
      double irmsforce_nonfixed = 0.0;
      int i;
      for (i=0; i<nfixed; i++) {
          if (fabs(iforce[i]) > imaxforce_all)
            imaxforce_all = fabs(iforce[i]);
          irmsforce_all += iforce[i]*iforce[i];
        }
      for (; i<iforce.dim(); i++) {
          if (fabs(iforce[i]) > imaxforce_all)
            imaxforce_all = fabs(iforce[i]);
          irmsforce_all += iforce[i]*iforce[i];
          if (fabs(iforce[i]) > imaxforce_nonfixed)
            imaxforce_nonfixed = fabs(iforce[i]);
          irmsforce_nonfixed += iforce[i]*iforce[i];
        }
      maxforce = irmsforce_nonfixed;
      rmsforce = irmsforce_nonfixed;
    }

  // test the size of the displacement, and the other convergence criteria
  // return if converged without changing the geometry...that's just plain
  // stupid, ed
  if (((idisp * iforce).maxval()/2) < conv_crit &&
      rmsforce < conv_rmsf &&
      maxforce < conv_maxf) {

    fprintf(outfp,"\n max of 1/2 idisp*iforce = %15.10g\n",
			      (idisp * iforce).maxval()/2);
    fprintf(outfp," max force               = %15.10g\n",maxforce);
    fprintf(outfp," rms force               = %15.10g\n",rmsforce);
    fprintf(outfp,"\n the geometry is converged\n");

    fprintf(outfp,"\n converged geometry\n");
    mol.print(outfp);

    if (keyval->booleanvalue("output_pdb")) write_pdb(keyval);

    return GEOM_DONE;
    }

 // scale the displacement vector if too large
  {
     double tot=0,scal=1.0;
     tot = sqrt(idisp.dot(idisp));
     //for(i=0; i < idisp.dim(); i++) tot += fabs(idisp[i]);
     if(tot>maxstepsize) {
       scal=maxstepsize/tot;
       fprintf(outfp,"\n stepsize of %f is too big, scaling by %f\n",tot,scal);
       }
     fprintf(outfp,"\n taking step of size %f\n",tot*scal);
     for(i=0; i < idisp.dim(); i++) idisp[i] *= scal;
     }

  if (cartesians) {
    rmsdisp=0;
    for (i=0; i < idisp.dim(); i++) {
      mol[i/3][i%3] += idisp[i] / 0.52917706;
      rmsdisp += (idisp[i] / 0.52917706)*(idisp[i] / 0.52917706);
      }
    
    rmsdisp = sqrt(rmsdisp/idisp.dim());
    maxdisp = idisp.maxval()/0.52917706;
    }

  // abort if the displacement fails
  if (!cartesians && cartesian_disp(bmat) < 0) {
    fprintf(outfp,"\n could not update geometry, saving state and exiting\n");
    STATEOUT so("geom.dat");
    so.put(iter);
    symm_coords->save_state(so);
    so.put(fixed.nonnull());
    if (fixed.nonnull()) fixed->save_state(so);
    mol.save_object_state(so);
    intco.save_object_state(so);
    iforce.save_object_state(so);
    idisp.save_object_state(so);
    hessian.save_object_state(so);
    return GEOM_ABORT;
    }
  fprintf(outfp,"\n displacements and new internal coordinates\n");
  for(i=0; i < intco.dim(); i++)
    fprintf(outfp," %5d %15.10f %15.10f\n",i+1,idisp[i],intco[i]);
  
  fprintf(outfp,"\n max of 1/2 idisp*iforce = %15.10g\n",
          (idisp * iforce).maxval()/2);
  fprintf(outfp," max force               = %15.10g\n",maxforce);
  fprintf(outfp," rms force               = %15.10g\n",rmsforce);
  fprintf(outfp," max disp                = %15.10g\n",maxdisp);
  fprintf(outfp," rms disp                = %15.10g\n",rmsdisp);

  fprintf(outfp,"\n updated geometry\n");
  mol.print(outfp);

  if (keyval->booleanvalue("output_pdb")) write_pdb(keyval);

  if (print_internal) {
      Geom_form_and_print_simples(keyval,
                                  "Updated Simple Internal Coordinates");
    }

  // now that we have a geometry, let's make the list of symm coordinates
  // if the user asks for it or fixed is nonnull
  if(recompute_internal || fixed) {
     if(get_symmco(keyval) < 0) {
        err_msg("Geom_init_mpqc: yikes!  can't make symm coords");
        return GEOM_ABORT;
        }
     }

  //StateOutBinXDR so("geom.dat");
  STATEOUT so("geom.dat");

  iter++;

  so.put(iter);
  symm_coords->save_state(so);
  so.put(fixed.nonnull());
  if (fixed.nonnull()) fixed->save_state(so);
  mol.save_object_state(so);
  intco.save_object_state(so);
  iforce.save_object_state(so);
  idisp.save_object_state(so);
  hessian.save_object_state(so);

  for(i=0; i < centers->n; i++) {
    centers->center[i].r[0] = mol[i][0];
    centers->center[i].r[1] = mol[i][1];
    centers->center[i].r[2] = mol[i][2];
    }

#endif
  return GEOM_COMPUTE_GRADIENT;
}
