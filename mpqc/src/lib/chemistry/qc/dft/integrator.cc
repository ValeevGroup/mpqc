//
// integrator.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/formio.h>
#include <chemistry/qc/dft/integrator.h>

///////////////////////////////////////////////////////////////////////////
// DenIntegrator

#define CLASSNAME DenIntegrator
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
DenIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

DenIntegrator::DenIntegrator(StateIn& s):
  SavableState(s)
{
  abort();
}

DenIntegrator::DenIntegrator()
{
  bs_values_ = 0;
  bsg_values_ = 0;
  alpha_dmat_ = 0;
  beta_dmat_ = 0;
  alpha_vmat_ = 0;
  beta_vmat_ = 0;
  compute_potential_integrals_ = 0;
}

DenIntegrator::DenIntegrator(const RefKeyVal& keyval)
{
  bs_values_ = 0;
  bsg_values_ = 0;
  alpha_dmat_ = 0;
  beta_dmat_ = 0;
  alpha_vmat_ = 0;
  beta_vmat_ = 0;
  compute_potential_integrals_ = 0;
}

DenIntegrator::~DenIntegrator()
{
  delete[] bs_values_;
  delete[] bsg_values_;
  delete[] alpha_dmat_;
  delete[] beta_dmat_;
  delete[] alpha_vmat_;
  delete[] beta_vmat_;
}

void
DenIntegrator::save_data_state(StateOut& s)
{
  cout << class_name() << ": cannot save state" << endl;
  abort();
}

void
DenIntegrator::set_wavefunction(const RefWavefunction &wfn)
{
  wfn_ = wfn;
}

void
DenIntegrator::init_integration(const RefDenFunctional &func)
{
  value_ = 0.0;

  wfn_->basis()->set_integral(wfn_->integral());

  func->set_compute_potential(compute_potential_integrals_);

  need_gradient_ = func->need_density_gradient();

  spin_polarized_ = wfn_->spin_polarized();

  nbasis_ = wfn_->basis()->nbasis();
  delete[] bs_values_;
  bs_values_ = new double[nbasis_];
  delete[] bsg_values_;
  bsg_values_ = new double[3*nbasis_];

  delete[] alpha_dmat_;
  alpha_dmat_ = new double[(nbasis_*(nbasis_+1))/2];
  RefSymmSCMatrix admat = wfn_->alpha_density();
  admat->convert(alpha_dmat_);
  admat = 0;

  delete[] beta_dmat_;
  beta_dmat_ = 0;
  if (wfn_->spin_polarized()) {
      beta_dmat_ = new double[(nbasis_*(nbasis_+1))/2];
      RefSymmSCMatrix bdmat = wfn_->beta_density();
      bdmat->convert(beta_dmat_);
      bdmat = 0;
    }

  delete[] alpha_vmat_;
  delete[] beta_vmat_;
  alpha_vmat_ = 0;
  beta_vmat_ = 0;
  if (compute_potential_integrals_) {
      int ntri = (nbasis_*(nbasis_+1))/2;
      alpha_vmat_ = new double[ntri];
      beta_vmat_ = new double[ntri];
      for (int i=0; i<ntri; i++) alpha_vmat_[i] = beta_vmat_[i] = 0.0;
    }
}

void
DenIntegrator::get_density(double *dmat, double &den, double &den_grad_mag)
{
  int i, j;
  double bvi, bvix, bviy, bviz;
  double densij;
  
  den = 0.0;
  den_grad_mag = 0.0;

  double tmp = 0.0;

  if (need_gradient_) {
      double grad[3];
      grad[0] = grad[1] = grad[2] = 0.0;

      int ij=0;
      for (i=0; i<nbasis_; i++) {
          bvi = bs_values_[i];
          bvix = bsg_values_[i*3];
          bviy = bsg_values_[i*3+1];
          bviz = bsg_values_[i*3+2];

          for (j=0; j < i; j++,ij++) {
              densij = dmat[ij];
              double bvj = bs_values_[j];

              tmp += 2.0*densij*bvi*bvj;
              grad[0] += densij*(bvi*bsg_values_[j*3+0] + bvj*bvix);
              grad[1] += densij*(bvi*bsg_values_[j*3+1] + bvj*bviy);
              grad[2] += densij*(bvi*bsg_values_[j*3+2] + bvj*bviz);
            }

          densij = dmat[ij]*bvi;
          tmp += densij*bvi;
          grad[0] += densij*bvix;
          grad[1] += densij*bviy;
          grad[2] += densij*bviz;

          ij++;
        }


      den = tmp;
  
      double x = grad[0] * 2.0;
      double y = grad[1] * 2.0;
      double z = grad[2] * 2.0;
      den_grad_mag = sqrt(x*x+y*y+z*z);
    }
  else {
      int ij=0;
      for (i=0; i<nbasis_; i++) {
          bvi = 2.0*bs_values_[i];

          for (j=0; j < i; j++,ij++)
              tmp += dmat[ij]*bvi*bs_values_[j];

          tmp += 0.25*dmat[ij]*bvi*bvi;
          ij++;
        }

      den = tmp;
    }
}

double
DenIntegrator::do_point(const SCVector3 &r,
                        const RefDenFunctional &func,
                        double weight)
{
  int i,j,k;

  // compute the basis set values
  if (need_gradient_) {
      wfn_->basis()->grad_values(r,bsg_values_,bs_values_);
    }
  else {
      wfn_->basis()->values(r,bs_values_);
    }

  // loop over basis functions adding contributions to the density
  PointInputData id;

  get_density(alpha_dmat_, id.dens_alpha, id.dens_grad_alpha);
  id.dens_alpha13 = pow(id.dens_alpha,1./3.);

  if (!spin_polarized_) {
      id.dens_beta = id.dens_alpha;
      id.dens_grad_beta = id.dens_grad_alpha;
      id.dens_beta13 = id.dens_alpha13;
    }
  else {
      get_density(beta_dmat_, id.dens_beta, id.dens_grad_beta);
      id.dens_beta13 = pow(id.dens_alpha,1./3.);
    }
  //id.dens13 = pow(id.dens_alpha + id.dens_beta, 1./3.);

  PointOutputData od;
  func->point(id, od);

  value_ += od.energy * weight;

  if (compute_potential_integrals_) {
      // the contribution to the potential integrals
      for (j=0; j<nbasis_; j++) {
          int joff = (j*(j+1))/2;
          double tmp = bs_values_[j] * weight;
          double tmpa = tmp * od.alpha_pot;
          double tmpb = tmp * od.beta_pot;
          for (k=0; k<=j; k++) {
              alpha_vmat_[joff+k] += tmpa * bs_values_[k];
              beta_vmat_[joff+k] += tmpb * bs_values_[k];
            }
        }
    }
  return id.dens_alpha + id.dens_beta;
}

///////////////////////////////////////////////////////////////////////////
// Murray93Integrator

// utility functions

static void
gauleg(double x1, double x2, double x[], double w[], int n)
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  const double EPS = 10.0 * DBL_EPSILON;

  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (i=1;i<=m;i++)  {
      z=cos(M_PI*(i-0.25)/(n+0.5));
      do {
          p1=1.0;
          p2=0.0;
          for (j=1;j<=n;j++) {
              p3=p2;
              p2=p1;
              p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	    }
          pp=n*(z*p1-p2)/(z*z-1.0);
          z1=z;
          z=z1-p1/pp;
	} while (fabs(z-z1) > EPS);
      x[i-1]=xm-xl*z;
      x[n-i]=xm+xl*z;
      w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
      w[n-i]=w[i-1];
    }
}

static double
calc_s(double m)
{
#if 1
  double value, value2;
  double a= -969969./262144.;
  double m2=m*m;
  value2= -3.814697265625e-06*m*(m2*(m2*(m2*(m2*(m2*(m2*(m2*(m2*(m2*(46189.*m2-510510.)+2567565.)-7759752.)+15668730.)-22221108.)+22632610.)-16628040.)+8729721.)-3233230.)+969969.);

  value = (1.+value2)/2.0;
    
  return value;
#else
  double m1 = 1.5*m - 0.5*m*m*m;
  double m2 = 1.5*m1 - 0.5*m1*m1*m1;
  double m3 = 1.5*m2 - 0.5*m2*m2*m2;
  return 0.5*(1.0-m3);
#endif  
}

static double **a_mat;
static double **oorab;

static double
calc_w(int this_center, SCVector3 &point, int ncenters,
       SCVector3 *centers, double *bragg_radius)
{
  int icenter, jcenter;
  double p_sum=0, p_point, p_tmp;
    
  for (icenter=0; icenter<ncenters; icenter++) {
      double ra = point.dist(centers[icenter]);
      double *ooraba = oorab[icenter];
      double *aa = a_mat[icenter];

      p_tmp = 1.0;
      for (jcenter=0; jcenter < ncenters; jcenter++) {
	if (icenter != jcenter) {
            double mu = (ra-point.dist(centers[jcenter]))*ooraba[jcenter];

            if (mu < -1.)
                continue; // s(-1) == 1.0
            else if (mu > 1.) {
                p_tmp = 0; // s(1) == 0.0
                break;
              }
            else
                p_tmp *= calc_s(mu + aa[jcenter]*(1.-mu*mu));
          }
      }

      if (icenter==this_center)
          p_point=p_tmp; 

      p_sum += p_tmp;
    }

  return p_point/p_sum;
}

#define CLASSNAME Murray93Integrator
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenIntegrator
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Murray93Integrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenIntegrator::_castdown(cd);
  return do_castdowns(casts,cd);
}

Murray93Integrator::Murray93Integrator(StateIn& s):
  DenIntegrator(s)
  maybe_SavableState(s)
{
  abort();
}

Murray93Integrator::Murray93Integrator()
{
  nr_ = 64;
  ntheta_ = 16;
  nphi_ = 32;
  Ktheta_ = 5;
}

Murray93Integrator::Murray93Integrator(const RefKeyVal& keyval):
  DenIntegrator(keyval)
{
  nr_ = keyval->intvalue("nr");
  if (nr_ == 0) nr_ = 64;
  ntheta_ = keyval->intvalue("ntheta");
  if (ntheta_ == 0) ntheta_ = 16;
  nphi_ = keyval->intvalue("nphi");
  if (nphi_ == 0) nphi_ = 2*ntheta_;
  Ktheta_ = keyval->intvalue("Ktheta");
  if (keyval->error() != KeyVal::OK)
      Ktheta_ = 5;
}

Murray93Integrator::~Murray93Integrator()
{
}

void
Murray93Integrator::save_data_state(StateOut& s)
{
  cout << ": cannot save state" << endl;
  abort();
}


// taken from Mike Colvin 1997/07/21 and modified
void
Murray93Integrator::integrate(const RefDenFunctional &denfunc)
{
  init_integration(denfunc);

  RefMolecule mol = wavefunction()->molecule();

  int ncenters=mol->natom();   // number of centers
  int icenter;                 // Loop index over centers

  int *nr = new int[ncenters];
  int ntheta = ntheta_;
  int nphi = nphi_;
  for (icenter=0; icenter<ncenters; icenter++) nr[icenter] = nr_;

  SCVector3 *centers = new SCVector3[ncenters];
  for (icenter=0; icenter<ncenters; icenter++) {
      centers[icenter].x() = mol->atom(icenter)[0];
      centers[icenter].y() = mol->atom(icenter)[1];
      centers[icenter].z() = mol->atom(icenter)[2];
    }

  double *bragg_radius = new double[ncenters];
  for (icenter=0; icenter<ncenters; icenter++) {
      bragg_radius[icenter] = mol->atom(icenter).element().bragg_radius();
    }

  if (a_mat) {
      delete[] a_mat[0];
      delete[] a_mat;
    }

  if (oorab) {
      delete[] oorab[0];
      delete[] oorab;
    }
  
  a_mat = new double*[ncenters];
  a_mat[0] = new double[ncenters*ncenters];
  oorab = new double*[ncenters];
  oorab[0] = new double[ncenters*ncenters];

  for (icenter=0; icenter < ncenters; icenter++) {
      if (icenter) {
          a_mat[icenter] = &a_mat[icenter-1][ncenters];
          oorab[icenter] = &oorab[icenter-1][ncenters];
        }

      double bragg_radius_a = bragg_radius[icenter];
      
      for (int jcenter=0; jcenter < ncenters; jcenter++) {
          double chi=bragg_radius_a/bragg_radius[jcenter];
          double uab=(chi-1.)/(chi+1.);
          a_mat[icenter][jcenter] = uab/(uab*uab-1.);
          oorab[icenter][jcenter] = 1./centers[icenter].dist(centers[jcenter]);
        }
    }

  int ir, itheta, iphi;       // Loop indices for diff. integration dim
  int point_count;            // Counter for # integration points per center
  int point_count_total=0;    // Counter for # integration points

  SCVector3 point;            // sph coord of current integration point
  SCVector3 center;           // Cartesian position of center
  SCVector3 integration_point;

  double w,q,int_volume,sin_theta;
        
  // Determine maximium # grid points
  int nr_max=0;
  for (icenter=0;icenter<ncenters;icenter++) {
      if (nr[icenter]>nr_max) nr_max=nr[icenter];
    }

  // allocate memory for integration sub-expressions
  double *theta_quad_points =  new double[ntheta];
  double *theta_quad_weights = new double[ntheta];
  double *r_values = new double[nr_max];
  double *dr_dq = new double[nr_max];

#if 0
  // Precalcute theta Guass-Legendre abcissas and weights
  gauleg(0.0, M_PI, theta_quad_points, theta_quad_weights, ntheta);

  for (icenter=0; icenter < ncenters; icenter++) {
      point_count=0;

      // Precalculate radial integration points
      for (ir=0; ir < nr[icenter]; ir++) {
          // Mike Colvin's interpretation of Murray's radial grid
          q=(double) ir/nr[icenter];
          double value=q/(1-q);
          r_values[ir]=bragg_radius[icenter]*value*value;
          dr_dq[ir]=2.0*bragg_radius[icenter]*q*pow(1.0-q,-3.);
        }

      center = centers[icenter];

      for (iphi=0; iphi < nphi; iphi++) {
          point.phi() = (double)iphi/(double)nphi * 2.0 * M_PI;

          for (itheta=0; itheta < ntheta; itheta++) {
              point.theta() = theta_quad_points[itheta];
              sin_theta = sin(point.theta());

              for (ir=0; ir < nr[icenter]; ir++) {
                  point.r() = r_values[ir];
                    
                  point.spherical_to_cartesian(integration_point);
                  integration_point += center;

                  //cout << "  " << integration_point << ": ";

                  // calculate weighting factor
                  w=calc_w(icenter, integration_point, ncenters,
                           centers, bragg_radius);
                    
                  // calculate integration volume 
                  int_volume=dr_dq[ir]*point.r()*point.r()*sin_theta;

                  // Update center point count
                  point_count++;

                  // Make contribution to Euler-Maclaurin formula
                  double multiplier = w*int_volume
                                    * theta_quad_weights[itheta]/nr[icenter]
                                    * 2.0 * M_PI / ((double)nphi);

                  //cout << scprintf("% 10.6f ", int_volume);
                  //cout << scprintf("% 10.6f ", multiplier);
                  if (do_point(integration_point, denfunc, multiplier)
                      * int_volume < DBL_EPSILON
                      && int_volume > DBL_EPSILON) break;
                }
            }
        }
      point_count_total+=point_count;
    }
#else

  for (icenter=0; icenter < ncenters; icenter++) {
      point_count=0;
      center = centers[icenter];
      int r_done = 0;

      for (ir=0; ir < nr[icenter]; ir++) {
          // Mike Colvin's interpretation of Murray's radial grid
          q=(double) ir/nr[icenter];
          double value=q/(1-q);
          double r = bragg_radius[icenter]*value*value;
          double dr_dq = 2.0*bragg_radius[icenter]*q*pow(1.0-q,-3.);
          double drdqr2 = dr_dq*r*r;
          point.r() = r;

          // Precalcute theta Guass-Legendre abcissas and weights
          if (ir==0) {
              ntheta=1;
              nphi=1;
            }
          else {
              ntheta = (int) (r/bragg_radius[icenter] * Ktheta_*ntheta_);
              if (ntheta > ntheta_)
                  ntheta = ntheta_;
              if (ntheta < 6)
                  ntheta=6;
              nphi = 2*ntheta;
            }

          gauleg(0.0, M_PI, theta_quad_points, theta_quad_weights, ntheta);

          for (itheta=0; itheta < ntheta; itheta++) {
              point.theta() = theta_quad_points[itheta];
              sin_theta = sin(point.theta());

              // calculate integration volume 
              int_volume=drdqr2*sin_theta;

              for (iphi=0; iphi < nphi; iphi++) {
                  point.phi() = (double)iphi/(double)nphi * 2.0 * M_PI;

                  point.spherical_to_cartesian(integration_point);
                  integration_point += center;

                  // calculate weighting factor
                  w=calc_w(icenter, integration_point, ncenters,
                           centers, bragg_radius);
                    
                  // Update center point count
                  point_count++;

                  // Make contribution to Euler-Maclaurin formula
                  double multiplier = w*int_volume
                                    * theta_quad_weights[itheta]/nr[icenter]
                                    * 2.0 * M_PI / ((double)nphi);

                  if (do_point(integration_point, denfunc, multiplier)
                      * int_volume < DBL_EPSILON
                      && int_volume > DBL_EPSILON) {
                      r_done=1;
                      break;
                    }
                }

              if (r_done)
                  break;
            }

          if (r_done)
              break;
        }
      point_count_total+=point_count;
    }
#endif  

    cout << " Total integration points = " << point_count_total << endl;
    cout << scprintf(" Value of integral = %16.14f", value()) << endl;

    delete[] theta_quad_points;
    delete[] theta_quad_weights;
    delete[] r_values;
    delete[] dr_dq;
    delete[] bragg_radius;
    delete[] nr;
    delete[] centers;
    delete[] a_mat[0];
    delete[] a_mat;
    a_mat=0;
    delete[] oorab[0];
    delete[] oorab;
    oorab=0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
