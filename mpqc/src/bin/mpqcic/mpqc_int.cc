
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

#include <util/misc/formio.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <util/unix/cct_cprot.h>
#include <util/misc/libmisc.h>

#include <math/optimize/update.h>
#include <math/scmat/local.h>

#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/localdef.h>

#include "mpqc_int.h"

#define VERBOSE 0

#define STATEOUT StateOutBinXDR
#define STATEIN StateInBinXDR
// #define STATEOUT StateOutBin
// #define STATEIN StateInBin
// #define STATEOUT StateOutText
// #define STATEIN StateInText

/////////////////////////////////////////////////////////////////

static RefMolecule mol;
static RefMolecularCoor coor;
static RefHessianUpdate update;

static RefSCMatrixKit matrixkit;

static RefSymmSCMatrix hessian;
static RefSCVector xn;
static RefSCVector gn;
static RefSCVector cart_grad;
static RefSCVector last_mode;
static double oldenergy = 0.0;

RefSCDimension di, dc;

static int iter=1;

static int checkpoint_geom=1;
static int print_hessian=0;
static int print_internal=0;
static int tstate=0;
static int efc=0;
static int modef=0;
static int use_internal_forces=-1;

static double conv_energy=1.0e-6;
static double conv_crit=1.0e-6;
static double conv_rmsf=1.0e-5;
static double conv_maxf=4.0e-5;
static double maxstepsize=0.6;

static double cart_tol=1.0e-10;

static double maxforce=0;
static double rmsforce=0;
static double maxdisp=0;
static double rmsdisp=0;


/////////////////////////////////////////////////////////////////

static void
get_input(const RefKeyVal& keyval)
{
  if (keyval->exists("delta_energy")) {
    double conv = fabs(keyval->doublevalue("delta_energy"));
    if (conv >= 1.0) conv_energy = pow(10.0,-conv);
    else conv_energy = conv;
  }

  if (keyval->exists("rms_force")) {
    double conv = fabs(keyval->doublevalue("rms_force"));
    if (conv >= 1.0) conv_rmsf = pow(10.0,-conv);
    else conv_rmsf = conv;
  }

  if (keyval->exists("max_force")) {
    double conv = fabs(keyval->doublevalue("max_force"));
    if (conv >= 1.0) conv_maxf = pow(10.0,-conv);
    else conv_maxf = conv;
  }

  if (keyval->exists("convergence")) {
    double conv = fabs(keyval->doublevalue("convergence"));
    if (conv >= 1.0) conv_crit = pow(10.0,-conv);
    else conv_crit = conv;
  }

  if (keyval->exists("maxstepsize")) {
    maxstepsize = keyval->doublevalue("maxstepsize");
  }

  if (keyval->exists("checkpoint")) {
    checkpoint_geom = keyval->booleanvalue("checkpoint");
  }

  if (keyval->exists("print_hessian")) {
    print_hessian = keyval->booleanvalue("print_hessian");
  }

  if (keyval->exists("print_internal")) {
    print_internal = keyval->booleanvalue("print_internal");
  }

  if (keyval->exists("eigenvector_following")) {
    efc = keyval->booleanvalue("eigenvector_following");
  }

  if (keyval->exists("use_internal_forces")) {
    use_internal_forces = keyval->booleanvalue("use_internal_forces");
  }
  
  if (keyval->exists("transition_state")) {
    efc = 1;
    tstate = keyval->booleanvalue("transition_state");
    if (keyval->exists("mode_following")) {
      modef = keyval->booleanvalue("mode_following");
    }
  }

  cout << scprintf("  intco:print_hessian       = %d\n",print_hessian);
  cout << scprintf("  intco:print_internal      = %d\n",print_internal);
  cout << scprintf("  intco:transition_state    = %d\n",tstate);
  cout << scprintf("  intco:mode_following      = %d\n",modef);
  cout << scprintf("  intco:maxstepsize         = %g\n",maxstepsize);
  cout << scprintf("  intco:cartesian_tolerance = %g\n",cart_tol);
  cout << scprintf("  intco:convergence         = %g\n",conv_crit);
  cout << scprintf("  intco:max_force           = %g\n",conv_maxf);
  cout << scprintf("  intco:rms_force           = %g\n",conv_rmsf);
  cout << scprintf("  intco:delta_energy        = %g\n",conv_energy);
  cout << endl;
}

////////////////////////////////////////////////////////////////////////////

int
Geom_init_mpqc(RefMolecule& molecule, const RefKeyVal& topkeyval,
               const char *geomfile)
{
  int i;

  matrixkit = SCMatrixKit::default_matrixkit();

 // create a new keyval which adds :intco to the search list 
  RefKeyVal keyval = new AggregateKeyVal(
      new PrefixKeyVal(":intco",topkeyval),
      new PrefixKeyVal(":default",topkeyval));
  
  get_input(keyval);

 // set mol = molecule
  mol = molecule;

 // if geom.dat exists, then read in all of the hooey stored there,
 // otherwise read intco stuff from mpqc.in

  struct stat buf;

  if (stat(geomfile,&buf) < 0 || buf.st_size==0) {
   // read coor and update from the input
    coor = keyval->describedclassvalue("coor");
    if (coor.null())
      coor = new SymmMolecularCoor(mol);
    
   // for now, always read the update from the input. default to bfgs
    update = keyval->describedclassvalue("update");
    if (update.null()) {
      if (tstate)
        update = new PowellUpdate;
      else
        update = new BFGSUpdate;
    }
    
   // use inverse hessian for QN, hessian for EFC
    if (!efc) update->set_inverse();

    dc = coor->dim_natom3();
    di = coor->dim();

   // create the inverted hessian 
    hessian = matrixkit->symmmatrix(di);
    coor->guess_hessian(hessian);
    if (!efc) hessian = hessian.gi();

    xn = matrixkit->vector(di);
    gn = matrixkit->vector(di);
    cart_grad = matrixkit->vector(dc);

    coor->to_internal(xn);
    gn->assign(0.0);
    cart_grad->assign(0.0);
    
   // save it all to disk
    if (checkpoint_geom) {
        STATEOUT so(geomfile,"w+");
        so.put(iter);
        mol.save_state(so);
        coor.save_state(so);
        update.save_state(so);
        dc.save_state(so);
        di.save_state(so);
        xn.save(so);
        gn.save(so);
        cart_grad.save(so);
        hessian.save(so);
      }
  }  else {
    STATEIN si(geomfile,"r+");

    si.get(iter);
    mol.restore_state(si);
    coor.restore_state(si);
    update.restore_state(si);
    dc.restore_state(si);
    di.restore_state(si);
    xn = matrixkit->vector(di);
    xn.restore(si);
    gn = matrixkit->vector(di);
    gn.restore(si);
    cart_grad = matrixkit->vector(dc);
    cart_grad.restore(si);
    hessian = matrixkit->symmmatrix(di);
    hessian.restore(si);

   // make sure molecule and mol refer to the same object
    molecule = mol;
    
    cout << scprintf(
            "\n restarting geometry optimization at iteration %d\n",iter);
  }

  if (use_internal_forces == -1) {
      if (coor->nconstrained()) use_internal_forces = 1;
      else use_internal_forces = 0;
    }
  cout << scprintf("  intco:use_internal_forces = %d\n", use_internal_forces);

  cout << scprintf("\n Initial geometry in Geom_init_mpqc\n");
  cout << scprintf("Molecule:") << endl;
  mol->print();
  
  if (print_internal) {
    cout << scprintf("\n Initial internal coordinates")
         << endl << endl;
    coor->print();
  } else {
    cout << scprintf("\n Initial simple internal coordinates");
    cout << endl << endl;
    coor->print_simples();
  }
  cout << endl;

  return GEOM_COMPUTE_GRADIENT;
}

void
Geom_done_mpqc(const RefKeyVal& keyval, int converged,
               const char *geomfile)
{
  const char *msg;

  if (converged)
    cout << scprintf("\nConverged Simple Internal Coordinates\n");
  else
    cout << scprintf("\nNonconverged Simple Internal Coordinates\n");

  cout << endl;
  
  coor->print_simples();

  cout << scprintf("\nFinal cartesian coordinates after %d iterations",iter);
  cout << endl << endl;
  cout << scprintf("Molecule:") << endl;;
  mol->print();
  cout << endl;

  coor=0;
  mol=0;
}

RefSCDimension
Geom_dim_natom3()
{
  return coor->dim_natom3();
}

int
Geom_update_mpqc(double energy, RefSCVector& grad, const RefKeyVal& keyval,
                 const char *geomfile)
{
  int i,j,ij;

  cart_grad.assign(grad);

  // transform cartesian coordinates to internal coordinates 
  coor->to_internal(xn);
  if (print_internal)
    xn.print("internal coordinates");

  // transform cartesian gradient to internal coordinates
  coor->to_internal(gn,cart_grad);
  if (print_internal)
    gn.print("internal coordinate gradients");

  // find rms and max force
  rmsforce=0;
  maxforce=0;
  RefSCVector conv_grad;
  if (use_internal_forces) conv_grad = gn;
  else conv_grad = cart_grad;
  for (i=0; i < conv_grad.n(); i++) {
    double d = fabs(conv_grad(i));
    rmsforce += d*d;
    maxforce = (d>maxforce) ? d : maxforce;
  }
  rmsforce = sqrt(rmsforce/conv_grad.n());

  // update the inverse hessian
  RefFunction nlp = 0;
  update->update(hessian,nlp,xn,gn);

  // possibly transform to a new coordinate system
  RefNonlinearTransform trans = coor->change_coordinates();

#if VERBOSE
  xn.print("xn before");
  gn.print("gn before");
  hessian.print("hessian before");
#endif

  update->apply_transform(trans);
  trans->transform_coordinates(xn);
  trans->transform_gradient(gn);
  if (efc) {
      trans->transform_hessian(hessian);
    }
  else {
      trans->transform_ihessian(hessian);
    }

#if VERBOSE
  xn.print("xn after");
  gn.print("gn after");
  hessian.print("hessian after");
#endif

  RefSCVector xdisp(hessian.dim(),matrixkit);
  
  if (!efc) {
    // take the Newton-Raphson step
    xdisp = -1.0*(hessian * gn);
  } else {
    // begin efc junk
    // first diagonalize hessian
    RefSCMatrix evecs(hessian.dim(),hessian.dim(),matrixkit);
    RefDiagSCMatrix evals(hessian.dim(),matrixkit);

    hessian.diagonalize(evals,evecs);
    //evals.print("hessian eigenvalues");
    //evecs.print("hessian eigenvectors");

    // form gradient to local hessian modes F = Ug
    RefSCVector F = evecs.t() * gn;
    //F.print("F");

    // figure out if hessian has the right number of negative eigenvalues
    int ncoord = evals.n();
    int npos=0,nneg=0;
    for (i=0; i < ncoord; i++) {
      if (evals.get_element(i) >= 0.0) npos++;
      else nneg++;
    }

    xdisp.assign(0.0);
  
    // for now, we always take the P-RFO for tstate (could take NR if
    // nneg==1, but we won't make that an option yet)
    if (tstate) {
      // no mode following yet, just follow lowest mode
      int mode = 0;

      if (modef) {
        printf("\n following mode %d\n",mode);
      }
      
      
      double bk = evals(mode);
      double Fk = F(mode);
      double lambda_p = 0.5*bk + 0.5*sqrt(bk*bk + 4*Fk*Fk);
    
      double lambda_n;
      double nlambda=1.0;
      do {
        lambda_n=nlambda;
        nlambda=0;
        for (i=0; i < ncoord; i++) {
          if (i==mode) continue;
          
          nlambda += F.get_element(i)*F.get_element(i) /
                    (lambda_n - evals.get_element(i));
        }
      } while(fabs(nlambda-lambda_n) > 1.0e-8);

      printf("\n  lambda_p = %g\n",lambda_p);
      printf("  lambda_n = %g\n",lambda_n);

      // form Xk
      double Fkobkl = F(mode)/(evals(mode)-lambda_p);
      for (j=0; j < F.n(); j++)
        xdisp.accumulate_element(j, -(evecs.get_element(j,mode) * Fkobkl));
    
      // form displacement x = sum -Fi*Vi/(bi-lam)
      for (i=0; i < F.n(); i++) {
        if (i==mode) continue;
      
        double Fiobil = F(i) / (evals(i)-lambda_n);
        for (j=0; j < F.n(); j++)
          xdisp.accumulate_element(j, -(evecs.get_element(j,i) * Fiobil));
      }
    
    // minimum search
    } else {
      // evaluate lambda
      double lambda;
      double nlambda=1.0;
      do {
        lambda=nlambda;
        nlambda=0;
        for (i=0; i < F.n(); i++) {
          double Fi = F(i);
          nlambda += Fi*Fi / (lambda - evals.get_element(i));
        }
      } while(fabs(nlambda-lambda) > 1.0e-8);
      printf("\n  lambda = %g\n",lambda);

      // form displacement x = sum -Fi*Vi/(bi-lam)
      for (i=0; i < F.n(); i++) {
        double Fiobil = F(i) / (evals(i)-lambda);
        for (j=0; j < F.n(); j++)
          xdisp.accumulate_element(j, -(evecs.get_element(j,i) * Fiobil));
      }
    }
  }
  
  //xdisp.print("internal coordinate displacements");

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

  double delta_energy = energy - oldenergy;
  oldenergy = energy;
  
  if (fabs(delta_energy) < conv_energy
      && iconv < conv_crit && rmsforce < conv_rmsf && maxforce < conv_maxf) {
    cout << scprintf("\n max of 1/2 idisp*iforce = %15.10g (%5.2g)\n",
            iconv, conv_crit);
    cout << scprintf(" max force               = %15.10g (%5.2g)\n",
            maxforce, conv_maxf);
    cout << scprintf(" rms force               = %15.10g (%5.2g)\n",
            rmsforce, conv_rmsf);
    cout << scprintf(" delta energy            = %15.10g (%5.2g)\n",
            delta_energy, conv_energy);
    cout << scprintf("\n the geometry is converged\n");

    cout << scprintf("\n converged geometry") << endl;
    mol->print();

    return GEOM_DONE;
  }

 // scale the displacement vector if it's too large
  double tot = sqrt(xdisp.scalar_product(xdisp));
  if (tot > maxstepsize) {
    double scal = maxstepsize/tot;
    cout << scprintf("\n stepsize of %f is too big, scaling by %f\n",tot,scal);
    xdisp.scale(scal);
    tot *= scal;
  }
  cout << scprintf("\n taking step of size %f\n",tot);

  // displace internal coordinates
  xn.accumulate(xdisp);
  //xn.print("new internal coordinates");
  
  cout << scprintf("\n max of 1/2 idisp*iforce = %15.10g\n",iconv);
  cout << scprintf(" max force               = %15.10g\n",maxforce);
  cout << scprintf(" rms force               = %15.10g\n",rmsforce);
  cout << scprintf(" max disp                = %15.10g\n",maxdisp);
  cout << scprintf(" rms disp                = %15.10g\n",rmsdisp);
  cout.flush();

  // now transform new internal coords back to cartesian coordinates
  Molecule foo_save = *(mol.pointer());
  if (coor->to_cartesian(xn) < 0) {
    // that failed, try steepest descent and get out of here
    cerr << scprintf("\n  trying cartesian steepest descent..."
                     "cross your fingers\n");
    for (i=0; i < cart_grad.n(); i++)
      (*mol.pointer())[i/3][i%3] = foo_save[i/3][i%3]-cart_grad[i];
  }

  cout << scprintf("\nnew molecular coordinates") << endl;
  mol->print();

  // checkpoint
  iter++;
  if (checkpoint_geom) {
      STATEOUT so(geomfile,"w+");
      so.put(iter);
      mol.save_state(so);
      coor.save_state(so);
      update.save_state(so);
      dc.save_state(so);
      di.save_state(so);
      xn.save(so);
      gn.save(so);
      cart_grad.save(so);
      hessian.save(so);
    }

  return GEOM_COMPUTE_GRADIENT;
}

void
Geom_write_pdb(const RefKeyVal& keyval, RefMolecule& mol,
               const char *pdbfile, char *title)
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
    pdbf = fopen(pdbfile,"a+");

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

  int i;
  for (i=0; i < mol->natom(); i++) {
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
