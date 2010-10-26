//
// dfttest.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//
#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

#include <signal.h>

#include <util/misc/math.h>

#include <util/misc/formio.h>
#include <util/group/pregtime.h>
#include <chemistry/qc/dft/functional.h>
#include <chemistry/qc/dft/am05.h>
#include <chemistry/qc/dft/integrator.h>

#include <chemistry/qc/dft/linkage.h>
#include <chemistry/qc/scf/linkage.h>

#ifdef HAVE_MPI
#include <mpi.h>
#include <util/group/messmpi.h>
#endif

#ifdef HAVE_FENV_H
#  include <fenv.h>
#endif

using namespace std;
using namespace sc;

inline static double
norm(double v[3])
{
  double x,y,z;
  return sqrt((x=v[0])*x + (y=v[1])*y + (z=v[2])*z);
}

inline static double
dot(double v[3], double w[3])
{
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

double *
density_matrix(const Ref<Wavefunction> &wfn)
{
  int nbasis = wfn->basis()->nbasis();
  RefSymmSCMatrix adens = wfn->alpha_ao_density();
  double * alpha_dmat = new double[(nbasis*(nbasis+1))/2];
  adens->convert(alpha_dmat);
  return alpha_dmat;
}

void
get_density(PointInputData::SpinData &d, const SCVector3 &r,
            const Ref<Wavefunction> &wfn, double *pdmat = 0)
{
  double *dmat;
  if (pdmat) dmat = pdmat;
  else dmat = density_matrix(wfn);

  double * bsg_values_ = new double[3*wfn->basis()->nbasis()];
  double * bs_values_ = new double[wfn->basis()->nbasis()];
  GaussianBasisSet::ValueData vdat(wfn->basis(), wfn->integral());
  wfn->basis()->grad_values(r,&vdat,bsg_values_,bs_values_);

  int i, j;

  double tmp = 0.0;
  double densij;
  double bvi, bvix, bviy, bviz;
  double bvixx, bviyx, bviyy, bvizx, bvizy, bvizz;

  const int X = PointInputData::X;
  const int Y = PointInputData::Y;
  const int Z = PointInputData::Z;
  const int XX = PointInputData::XX;
  const int YX = PointInputData::YX;
  const int YY = PointInputData::YY;
  const int ZX = PointInputData::ZX;
  const int ZY = PointInputData::ZY;
  const int ZZ = PointInputData::ZZ;

  double grad[3], hess[6];
  for (i=0; i<3; i++) grad[i] = 0.0;
  for (i=0; i<6; i++) hess[i] = 0.0;

  int nbasis_ = wfn->basis()->nbasis();
  int ij=0;
  for (i=0; i < nbasis_; i++) {
    bvi = bs_values_[i];
    bvix = bsg_values_[i*3+X];
    bviy = bsg_values_[i*3+Y];
    bviz = bsg_values_[i*3+Z];
    for (j=0; j < i; j++,ij++) {
      densij = dmat[ij];
      double bvj = bs_values_[j];
      double bvjx = bsg_values_[j*3+X];
      double bvjy = bsg_values_[j*3+Y];
      double bvjz = bsg_values_[j*3+Z];

      tmp += 2.0*densij*bvi*bvj;

      grad[X] += densij*(bvi*bvjx + bvj*bvix);
      grad[Y] += densij*(bvi*bvjy + bvj*bviy);
      grad[Z] += densij*(bvi*bvjz + bvj*bviz);
      }

    densij = dmat[ij]*bvi;
    tmp += densij*bvi;
    grad[X] += densij*bvix;
    grad[Y] += densij*bviy;
    grad[Z] += densij*bviz;
    ij++;
    }


  d.rho = tmp;
  for (i=0; i<3; i++) d.del_rho[i] = 2.0 * grad[i];
  for (i=0; i<6; i++) d.hes_rho[i] = 2.0 * hess[i];

  d.lap_rho = d.hes_rho[XX] + d.hes_rho[YY] + d.hes_rho[ZZ];

  d.gamma = dot(d.del_rho,d.del_rho);

  delete[] bsg_values_;
  delete[] bs_values_;
  if (!pdmat) delete[] dmat;
}

double
fd_test_do_point(const SCVector3 &point,
                 const Ref<DenFunctional> &func, const Ref<Wavefunction> &wfn,
                 double *frozen_dmat = 0)
{
  PointInputData id(point);
  get_density(id.a, point, wfn, frozen_dmat);
  id.compute_derived(0,func->need_density_gradient(),false);
  PointOutputData od;
  if ( (id.a.rho + id.b.rho) > 1e2*DBL_EPSILON) func->point(id, od);
  else return 0.0;
  return od.energy;
}

void
fd_test_point(int acenter, const SCVector3 &tpoint,
              const Ref<DenFunctional> &functional, const Ref<Wavefunction> &wfn)
{
  SCVector3 point(tpoint);
  Ref<Molecule> mol = wfn->molecule();

  double *fd_grad_f = new double[mol->natom()*3];
  memset(fd_grad_f,0, 3*mol->natom() * sizeof(double));

  double *dmat = density_matrix(wfn);

  // frozen_dmat = dmat makes the overlap derivative contribution 0
  double *frozen_dmat = dmat;

  double delta = 0.0001;
  for (int i=0; i<mol->natom(); i++) {
    for (int j=0; j<3; j++) {
      if (acenter == i) point[j] += delta;
      mol->r(i,j) += delta;
      wfn->obsolete();
      double f_plus = fd_test_do_point(point, functional, wfn, frozen_dmat);
      if (acenter == i) point[j] -= 2*delta;
      mol->r(i,j) -= 2*delta;
      wfn->obsolete();
      double f_minus = fd_test_do_point(point, functional, wfn, frozen_dmat);
      if (acenter == i) point[j] += delta;
      mol->r(i,j) += delta;
      wfn->obsolete();
      fd_grad_f[i*3+j] = (f_plus-f_minus)/(2.0*delta);
      }
    }

  double * bsh_values_ = new double[6*wfn->basis()->nbasis()];
  double * bsg_values_ = new double[3*wfn->basis()->nbasis()];
  double * bs_values_ = new double[wfn->basis()->nbasis()];
  GaussianBasisSet::ValueData vdat(wfn->basis(), wfn->integral());
  wfn->basis()->hessian_values(point,&vdat,bsh_values_,bsg_values_,bs_values_);

  PointInputData id(point);
  get_density(id.a, point, wfn);
  id.compute_derived(0,functional->need_density_gradient(),false);

  PointOutputData od;
  functional->set_compute_potential(1);
  if ( (id.a.rho + id.b.rho) > 1e2*DBL_EPSILON) functional->point(id, od);
  else od.zero();

  double *an_grad_f = new double[mol->natom()*3];
  memset(an_grad_f,0, 3*mol->natom() * sizeof(double));

  int ncontrib = wfn->basis()->nshell();
  int ncontrib_bf = wfn->basis()->nbasis();
  int *contrib = new int[ncontrib];
  int *contrib_bf = new int[ncontrib_bf];
  for (int i=0; i<ncontrib; i++) contrib[i] = i;
  for (int i=0; i<ncontrib_bf; i++) contrib_bf[i] = i;
  functional->gradient(id, od, an_grad_f, acenter, wfn->basis(),
                       dmat, dmat,
                       ncontrib, contrib,
                       ncontrib_bf, contrib_bf,
                       bs_values_, bsg_values_, bsh_values_);
  delete[] contrib;
  delete[] contrib_bf;

  cout << " acenter = " << acenter << " point = " << point << endl;
  cout << "FD df/dx:" << endl;
  for (int i=0; i<mol->natom(); i++) {
    cout << scprintf(" % 16.12f % 16.12f % 16.12f",
                     fd_grad_f[3*i+0],
                     fd_grad_f[3*i+1],
                     fd_grad_f[3*i+2])
                     << endl;
    }

  cout << "AN df/dx:" << endl;
  for (int i=0; i<mol->natom(); i++) {
    cout << scprintf(" % 16.12f % 16.12f % 16.12f",
                     an_grad_f[3*i+0],
                     an_grad_f[3*i+1],
                     an_grad_f[3*i+2])
                     << endl;
    }

  delete[] bsh_values_;
  delete[] bsg_values_;
  delete[] bs_values_;
  delete[] dmat;
  delete[] fd_grad_f;
  delete[] an_grad_f;
}

void
fd_test(const Ref<DenFunctional> &functional, const Ref<Wavefunction> &wfn)
{
  cout << "fd_test with functional:" << endl;
  cout << functional;
  for (int i=0; i<wfn->molecule()->natom(); i++) {
    for (double x = -1.0; x <= 1.0; x++) {
      for (double y = -1.0; y <= 1.0; y++) {
        for (double z = -1.0; z <= 1.0; z++) {
          fd_test_point(i, SCVector3(x,y,z), functional, wfn);
          }
        }
      }
    }
}

void
fd_e_test(const Ref<Wavefunction> &wfn)
{
  Ref<Molecule> mol = wfn->molecule();

  cout << "Testing dE/dx with:" << endl;
  cout << incindent;
  wfn->print();
  cout << decindent;

  double *fd_grad_e = new double[mol->natom()*3];
  memset(fd_grad_e,0, 3*mol->natom() * sizeof(double));

  double delta = 0.0001;
  for (int i=0; i<mol->natom(); i++) {
    for (int j=0; j<3; j++) {
      mol->r(i,j) += delta;
      wfn->obsolete();
      double e_plus = wfn->energy();
      mol->r(i,j) -= 2*delta;
      wfn->obsolete();
      double e_minus = wfn->energy();
      mol->r(i,j) += delta;
      wfn->obsolete();
      fd_grad_e[i*3+j] = (e_plus-e_minus)/(2.0*delta);
      }
    }

  double *an_grad_e = new double[mol->natom()*3];
  memset(an_grad_e,0, 3*mol->natom() * sizeof(double));

  RefSCVector grad = wfn->get_cartesian_gradient();
  grad->convert(an_grad_e);

  cout << "FD dE/dx:" << endl;
  for (int i=0; i<mol->natom(); i++) {
    cout << scprintf(" % 16.12f % 16.12f % 16.12f",
                     fd_grad_e[3*i+0],
                     fd_grad_e[3*i+1],
                     fd_grad_e[3*i+2])
                     << endl;
    }

  cout << "AN dE/dx:" << endl;
  for (int i=0; i<mol->natom(); i++) {
    cout << scprintf(" % 16.12f % 16.12f % 16.12f",
                     an_grad_e[3*i+0],
                     an_grad_e[3*i+1],
                     an_grad_e[3*i+2])
                     << endl;
    }

  delete[] fd_grad_e;
  delete[] an_grad_e;
}

#define USE_ORIGINAL_CODE 0
#define USE_MPQC_CODE 1
#if USE_ORIGINAL_CODE
extern "C" easypbe_(double *up,double *agrup,double *delgrup,double *uplap,
                    double *dn,double *agrdn,double *delgrdn,double *dnlap,
                    double *agr,double *delgr,int *lcor,int *lpot,
                    double *exlsd,double *vxuplsd,double *vxdnlsd,
                    double *eclsd,double *vcuplsd,double *vcdnlsd,
                    double *expw91,double *vxuppw91,double *vxdnpw91,
                    double *ecpw91,double *vcuppw91,double *vcdnpw91,
                    double *expbe,double *vxuppbe,double *vxdnpbe,
                    double *ecpbe,double *vcuppbe,double *vcdnpbe);
#endif

void
do_valtest(const Ref<DenFunctional> &valtest)
{
  valtest->set_spin_polarized(1);

  SCVector3 point = 0.0;
  PointInputData id(point);
  // zero out data that is never used
  int ii=0;
  for (ii=0; ii<3; ii++) id.a.del_rho[ii] = 0.0;
  for (ii=0; ii<3; ii++) id.b.del_rho[ii] = 0.0;
  id.a.lap_rho = 0.0;
  id.b.lap_rho = 0.0;
  for (ii=0; ii<6; ii++) id.a.hes_rho[ii] = 0.0;
  for (ii=0; ii<6; ii++) id.b.hes_rho[ii] = 0.0;

  // taken from PBE.f (PBE alpha2.1)
  const double thrd = 1.0/3.0;
  const double thrd2 = 2.0/3.0;
  const double pi = M_PI;
  double conf=pow((3.0*pi*pi),thrd);
  double conrs=pow((3.0/(4.0*pi)),thrd);
  cout << "     Fup Fdn Zup Zdn             Exc" << endl;
  // BEGIN THE LOOP THAT SELECTS A TRIAL DENSITY
  // spin-densities are of the form
  //          rho(r)=f*(Z**3/pi)*dexp(-2*Z*r)
  // delzdn=small change in zdn to test potentials
  // jdens=counter for which density
  for (int jdens = 1; jdens <= 10; jdens++) {
    double fup=1.0;
    double fdn=0.2*(jdens-1);
    double zup=1.0;
    double zdn=0.5;
    if (jdens > 6) {
      fdn=1.0;
      zup=0.5+0.5*(jdens-7);
      zdn=zup;
      }
    double delzdn=1e-5;
    double sumexc, mpqc_sumexc;
    double sumexco;
    // BEGIN THE LOOP THAT INCREMENTS THE DENSITY DIFFERENTIALLY
    // kdif=1=>density as above
    // kdif=2=>Zdn changed by DELZDN
    for (int kdif=1; kdif<=2; kdif++) {
      if (kdif == 2) zdn=zdn+delzdn;
      // BEGIN THE RADIAL LOOP
      // sumexc=integrated exchange-correlation energy 
      // chng1=integrated xc energy change, based on vxc
      // nr=number of points in radial loop
      // rf=final value of r in integrals
      // dr=change in r
      // wt=weight of r in trapezoidal rule
      // dup=up density
      // agrup=|grad up|
      // delgrup=(grad up).(grad |grad up|) 
      // uplap=grad^2 up=Laplacian of up
      // dn,agrdn,delgrdn,dnlap=corresponding down quantities
      // d=up+dn
      // agrad=|grad rho|
      // delgrad=(grad rho).(grad |grad rho|) 
      sumexc=0.0;
      mpqc_sumexc = 0.0;
      sumexco=0.0;
      double chng1=0.0;
      int nr=10000;
      double rf=20.0;
      double dr=rf/nr;
      for (int i=1; i<=nr; i++) {
        double r=i*dr;
        double wt=4.*pi*r*r*dr;
        double dup=fup*(zup*zup*zup/pi)*exp(-2.0*zup*r);
        double ddn=fdn*(zdn*zdn*zdn/pi)*exp(-2.0*zdn*r);
        if (dup+ddn < DBL_EPSILON) continue;
        double zdnnu=zdn+delzdn;
        double delddn=fdn*(zdnnu*zdnnu*zdnnu/pi)*exp(-2.0*zdnnu*r)-ddn;
        double agrup=2.0*zup*dup;
        double delgrup=8.0*(zup*zup*zup)*dup*dup;
        double uplap=4.0*zup*dup*(zup-1.0/r);
        double agrdn=2.0*zdn*ddn;
        double delgrdn=8.0*(zdn*zdn*zdn)*ddn*ddn;
        double dnlap=4.0*zdn*ddn*(zdn-1.0/r);
        double d=dup+ddn;
        double agrad=2.0*(zup*dup+zdn*ddn);
        double delgrad=4.0*agrad*(zup*zup*dup+zdn*zdn*ddn);
#if USE_ORIGINAL_CODE
        double exlsd;
        double vxuplsd;
        double vxdnlsd;
        double exclsd;
        double vxcuplsd;
        double vxcdnlsd;
        double expw91,vxuppw91,vxdnpw91,ecpw91;
        double expbe,vxuppbe,vxdnpbe,ecpbe;
        double eclsd, vcuplsd, vcdnlsd, vcuppw91, vcdnpw91, vcuppbe, vcdnpbe;
        int ione=1;
        easypbe_(&dup,&agrup,&delgrup,&uplap,&ddn,&agrdn,&delgrdn,
                 &dnlap,&agrad,&delgrad,&ione,&ione,
                 &exlsd,&vxuplsd,&vxdnlsd,&eclsd,&vcuplsd,&vcdnlsd,
                 &expw91,&vxuppw91,&vxdnpw91,&ecpw91,&vcuppw91,&vcdnpw91,
                 &expbe,&vxuppbe,&vxdnpbe,&ecpbe,&vcuppbe,&vcdnpbe);
        //sumexc=sumexc+d*(expbe+ecpbe)*wt;
        sumexc=sumexc+d*(expbe+ecpbe)*wt;
        // CHNG1=CHNG1+(vxdnpbe+vcdnpbe)*DELDDN*WT/DELZDN
#endif
#if USE_MPQC_CODE
        PointOutputData od;
        id.a.rho = dup;
        id.a.gamma = agrup*agrup;
        id.b.rho = ddn;
        id.b.gamma = agrdn*agrdn;
        id.gamma_ab = 0.5*(agrad*agrad-id.a.gamma-id.b.gamma);
        if (id.gamma_ab > sqrt(id.a.gamma*id.b.gamma))
          id.gamma_ab = sqrt(id.a.gamma*id.b.gamma);
        if (id.gamma_ab < -sqrt(id.a.gamma*id.b.gamma))
          id.gamma_ab = -sqrt(id.a.gamma*id.b.gamma);
        if (id.gamma_ab < -0.5*(id.a.gamma*id.b.gamma))
          id.gamma_ab = -0.5*(id.a.gamma*id.b.gamma);
        id.compute_derived(1, valtest->need_density_gradient(),false);
        valtest->point(id,od);
        mpqc_sumexc += od.energy*wt;
#endif
//          cout << scprintf("d = %12.8f wt = %12.8f OK = %12.8f MPQC = %12.8f",
//                           d, wt, expbe, od.energy/d) << endl;
        }
      if(kdif==1) {
        sumexco=sumexc;
        }
      }
    // CHNG: DIRECT XC ENERGY INCREMENT
    // IF THE FUNCTIONAL DERIVATIVE IS CORRECT, THEN CHNG1=CHNG
//      CHNG=(sumEXC-sumEXCO)/DELZDN
//      PRINT 200,FUP,FDN,ZUP,ZDN,sumEXC,CHNG1,chng
#if USE_ORIGINAL_CODE
    cout << scprintf("orig %3.1f %3.1f %3.1f %3.1f %16.12f",
                     fup,fdn,zup,zdn,sumexc) << endl;
#endif
#if USE_MPQC_CODE
    cout << scprintf("mpqc %3.1f %3.1f %3.1f %3.1f %16.12f",
                     fup,fdn,zup,zdn,mpqc_sumexc) << endl;
#endif
    }
}

int
main(int argc, char**argv)
{


  ExEnv::init(argc, argv);

  Ref<MessageGrp> grp;
#if defined(HAVE_MPI) && defined(ALWAYS_USE_MPI)
  grp = new MPIMessageGrp(&argc, &argv);
  MessageGrp::set_default_messagegrp(grp);
#endif

#if 0
#ifdef HAVE_FEENABLEEXCEPT
// this uses a glibc extension to trap on individual exceptions
# ifdef FE_DIVBYZERO
  feenableexcept(FE_DIVBYZERO);
# endif
# ifdef FE_INVALID
  feenableexcept(FE_INVALID);
# endif
# ifdef FE_OVERFLOW
  feenableexcept(FE_OVERFLOW);
# endif
#endif
#endif

#ifdef HAVE_FEDISABLEEXCEPT
// this uses a glibc extension to not trap on individual exceptions
# ifdef FE_UNDERFLOW
  fedisableexcept(FE_UNDERFLOW);
# endif
# ifdef FE_INEXACT
  fedisableexcept(FE_INEXACT);
# endif
#endif

  int i;
  const char *input = (argc > 1)? argv[1] : SRCDIR "/dfttest.in";

  // open keyval input
  Ref<KeyVal> keyval(new ParsedKeyVal(input));

  cout << "=========== Value f Tests ===========" << endl;
  int nvaltest = keyval->count("valtest");
  for (i=0; i<nvaltest; i++) {
    Ref<DenFunctional> valtest;
    valtest << keyval->describedclassvalue("valtest", i);
    if (valtest.nonnull()) valtest->print();
    do_valtest(valtest);
    }

  Ref<Wavefunction> dft; dft << keyval->describedclassvalue("dft");
  if (dft.nonnull()) {
    cout << "=========== FD dE/dx Tests ===========" << endl;
    fd_e_test(dft);
    }

  cout << "=========== FD df/drho Tests ===========" << endl;
  Ref<DenFunctional> funcs[] = {
    new AM05Functional,
    new PBECFunctional,
    new PW91CFunctional,
    new PW91XFunctional,
    new PBEXFunctional,
    new PW92LCFunctional,
    new mPW91XFunctional(mPW91XFunctional::B88),
    new mPW91XFunctional(mPW91XFunctional::PW91),
    new mPW91XFunctional(mPW91XFunctional::mPW91),
    new SlaterXFunctional,
    new Becke88XFunctional,
    new VWN1LCFunctional(1),
    new VWN1LCFunctional,
    new VWN2LCFunctional,
    new VWN3LCFunctional,
    new VWN4LCFunctional,
    new VWN5LCFunctional,
    new PZ81LCFunctional,
    new P86CFunctional,
//     new NewP86CFunctional,
    new XalphaFunctional,
    new LYPCFunctional,
    new PW86XFunctional,
    new G96XFunctional,
    0
  };
  const int maxerr = 1000;
  int errcount[maxerr];
  for (i=0; funcs[i]; i++) {
    cout << "-----------------"
         << funcs[i]->class_name()
         << "-----------------"
         << endl;
    int nerr = funcs[i]->test();
    if (i<maxerr) errcount[i] = nerr;
    else { cout << "dfttest: maxerr exceeded" << endl; abort(); }
    }
  cout << "-------------- ERROR RESULTS --------------" << endl;
  for (int i=0; funcs[i]; i++) {
    cout << funcs[i]->class_name() << ": " << errcount[i];
    if (errcount[i] == 0) cout << " (OK)";
    cout << endl;
    }

  Ref<Molecule> mol; mol << keyval->describedclassvalue("molecule");
  if (mol.nonnull()) {
    cout << "=========== FD Weight Tests ===========" << endl;
    Ref<IntegrationWeight> weights[] = {
      new BeckeIntegrationWeight,
      0
    };
    for (i=0; weights[i]; i++) {
      cout << "-----------------"
           << weights[i]->class_name()
           << "-----------------"
           << endl;
      weights[i]->init(mol,1.0e-8);
      weights[i]->test();
      }
    }

  Ref<DenFunctional> functional;
  functional << keyval->describedclassvalue("functional");
  Ref<Wavefunction>  wfn;
  wfn << keyval->describedclassvalue("wfn");
  if (functional.nonnull() && wfn.nonnull()) {
    cout << "=========== FD df/dx Tests ===========" << endl;
    fd_test(functional, wfn);
    }

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
