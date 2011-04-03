//
// efc.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#include <math.h>
#include <float.h>

#include <util/state/stateio.h>
#include <math/optimize/efc.h>
#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// EFCOpt

static ClassDesc EFCOpt_cd(
  typeid(EFCOpt),"EFCOpt",2,"public Optimize",
  0, create<EFCOpt>, create<EFCOpt>);

EFCOpt::EFCOpt(const Ref<KeyVal>&keyval):
  Optimize(keyval),
  maxabs_gradient(-1.0)
{
  update_ << keyval->describedclassvalue("update");
  
  accuracy_ = keyval->doublevalue("accuracy");
  if (keyval->error() != KeyVal::OK) accuracy_ = 0.0001;

  tstate = keyval->booleanvalue("transition_state");
  if (keyval->error() != KeyVal::OK) tstate = 0;

  modef = keyval->booleanvalue("mode_following");
  if (keyval->error() != KeyVal::OK) modef = 0;

  if (tstate)
    ExEnv::out0() << endl << indent
         << "performing a transition state search\n\n";
  
  RefSymmSCMatrix hessian(dimension(),matrixkit());
  // get a guess hessian from function
  function()->guess_hessian(hessian);
  
  // see if any hessian matrix elements have been given in the input
  if (keyval->exists("hessian")) {
    int n = hessian.n();
    for (int i=0; i<n; i++) {
      if (keyval->exists("hessian",i)) {
        for (int j=0; j<=i; j++) {
          double tmp = keyval->doublevalue("hessian",i,j);
          if (keyval->error() == KeyVal::OK) hessian(i,j) = tmp;
        }
      }
    }
  }
  hessian_ = hessian;
  last_mode_ = 0;
}

EFCOpt::EFCOpt(StateIn&s):
  SavableState(s),
  Optimize(s)
{
  s.get(tstate);
  s.get(modef);
  hessian_ = matrixkit()->symmmatrix(dimension());
  hessian_.restore(s);
  update_ << SavableState::restore_state(s);
  last_mode_ = matrixkit()->vector(dimension());
  last_mode_.restore(s);
  if (s.version(::class_desc<EFCOpt>()) < 2) {
    double convergence;
    s.get(convergence);
  }
  s.get(accuracy_);
  s.get(maxabs_gradient);
}

EFCOpt::~EFCOpt()
{
}

void
EFCOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  s.put(tstate);
  s.put(modef);
  hessian_.save(s);
  SavableState::save_state(update_.pointer(),s);
  last_mode_.save(s);
  s.put(accuracy_);
  s.put(maxabs_gradient);
}

void
EFCOpt::init()
{
  Optimize::init();
  maxabs_gradient = -1.0;
}

int
EFCOpt::update()
{
  int i,j;
  
  // these are good candidates to be input options
  const double maxabs_gradient_to_desired_accuracy = 0.05;
  const double maxabs_gradient_to_next_desired_accuracy = 0.005;
  const double roundoff_error_factor = 1.1;

  // the gradient convergence criterion.
  double old_maxabs_gradient = maxabs_gradient;
  RefSCVector xcurrent;
  RefSCVector gcurrent;

  ExEnv::out0().flush();
    
  // get the next gradient at the required level of accuracy.
  // usually only one pass is needed, unless we happen to find
  // that the accuracy was set too low.
  int accurate_enough;
  do {
    // compute the current point
    function()->set_desired_gradient_accuracy(accuracy_);
    
    xcurrent = function()->get_x();
    gcurrent = function()->gradient().copy();

    // compute the gradient convergence criterion now so i can see if
    // the accuracy needs to be tighter
    maxabs_gradient = gcurrent.maxabs();
    // compute the required accuracy
    accuracy_ = maxabs_gradient * maxabs_gradient_to_desired_accuracy;

    if (accuracy_ < DBL_EPSILON) accuracy_ = DBL_EPSILON;

    // The roundoff_error_factor is thrown in to allow for round off making
    // the current gcurrent.maxabs() a bit smaller than the previous,
    // which would make the current required accuracy less than the
    // gradient's actual accuracy and cause everything to be recomputed.
    accurate_enough = (function()->actual_gradient_accuracy() <=
                       accuracy_*roundoff_error_factor);

    if (!accurate_enough) {
      ExEnv::out0() << indent
           << "NOTICE: function()->actual_gradient_accuracy() > accuracy_:\n"
           << indent << scprintf(
             "        function()->actual_gradient_accuracy() = %15.8e",
             function()->actual_gradient_accuracy()) << endl
           << scprintf(
             "                                     accuracy_ = %15.8e",
             accuracy_) << endl;
    }
  } while(!accurate_enough);

  if (old_maxabs_gradient >= 0.0 && old_maxabs_gradient < maxabs_gradient) {
    ExEnv::out0() << indent
         << scprintf("NOTICE: maxabs_gradient increased from %8.4e to %8.4e",
                     old_maxabs_gradient, maxabs_gradient) << endl;
  }

  // update the hessian
  if (update_.nonnull()) {
    update_->update(hessian_,function(),xcurrent,gcurrent);
  }

  // begin efc junk
  // first diagonalize hessian
  RefSCMatrix evecs(dimension(),dimension(),matrixkit());
  RefDiagSCMatrix evals(dimension(),matrixkit());

  hessian_.diagonalize(evals,evecs);
  //evals.print("hessian eigenvalues");
  //evecs.print("hessian eigenvectors");

  // form gradient to local hessian modes F = Ug
  RefSCVector F = evecs.t() * gcurrent;
  //F.print("F");

  // figure out if hessian has the right number of negative eigenvalues
  int ncoord = evals.n();
  int npos=0,nneg=0;
  for (i=0; i < ncoord; i++) {
    if (evals.get_element(i) >= 0.0) npos++;
    else nneg++;
  }

  RefSCVector xdisp(dimension(),matrixkit());
  xdisp.assign(0.0);
  
  // for now, we always take the P-RFO for tstate (could take NR if
  // nneg==1, but we won't make that an option yet)
  if (tstate) {
    int mode = 0;

    if (modef) {
      // which mode are we following.  find mode with maximum overlap with
      // last mode followed
      if (last_mode_.nonnull()) {
        double overlap=0;
        for (i=0; i < ncoord; i++) {
          double S=0;
          for (j=0; j < ncoord; j++) {
            S += last_mode_.get_element(j)*evecs.get_element(j,i);
          }
          S = fabs(S);
          if (S > overlap) {
            mode = i;
            overlap = S;
          }
        }
      } else {
        last_mode_ = matrixkit()->vector(dimension());
      
        // find mode with max component = coord 0 which should be the
        // mode being followed
        double comp=0;
        for (i=0; i < ncoord; i++) {
          double S = fabs(evecs.get_element(0,i));
          if (S>comp) {
            mode=i;
            comp=S;
          }
        }
      }
    
      for (i=0; i < ncoord; i++)
        last_mode_(i) = evecs(i,mode);

      ExEnv::out0() << endl << indent << "\n following mode " << mode << endl;
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

    ExEnv::out0()
         << indent << scprintf("lambda_p = %8.5g",lambda_p) << endl
         << indent << scprintf("lambda_n = %8.5g",lambda_n) << endl;

    // form Xk
    double Fkobkl = F(mode)/(evals(mode)-lambda_p);
    for (j=0; j < F.n(); j++)
      xdisp(j) = xdisp(j) - evecs(j,mode) * Fkobkl;
    
    // form displacement x = sum -Fi*Vi/(bi-lam)
    for (i=0; i < F.n(); i++) {
      if (i==mode) continue;
      
      double Fiobil = F(i) / (evals(i)-lambda_n);
      for (j=0; j < F.n(); j++) {
        xdisp(j) = xdisp(j) - evecs(j,i) * Fiobil;
      }
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

    ExEnv::out0() << indent << scprintf("lambda = %8.5g", lambda) << endl;

  // form displacement x = sum -Fi*Vi/(bi-lam)
    for (i=0; i < F.n(); i++) {
      double Fiobil = F(i) / (evals(i)-lambda);
      for (j=0; j < F.n(); j++) {
        xdisp(j) = xdisp(j) - evecs(j,i) * Fiobil;
      }
    }
  }

  // scale the displacement vector if it's too large
  double tot = sqrt(xdisp.scalar_product(xdisp));
  if (tot > max_stepsize_) {
    double scal = max_stepsize_/tot;
    ExEnv::out0() << endl << indent
         << scprintf("stepsize of %f is too big, scaling by %f",tot,scal)
         << endl;
    xdisp.scale(scal);
    tot *= scal;
  }

  //xdisp.print("xdisp");

  // try steepest descent
  // RefSCVector xdisp = -1.0*gcurrent;
  RefSCVector xnext = xcurrent + xdisp;

  conv_->reset();
  conv_->get_grad(function());
  conv_->get_x(function());
  conv_->set_nextx(xnext);

  // check for conergence before resetting the geometry
  int converged = conv_->converged();
  if (converged)
    return converged;

  ExEnv::out0() << endl
       << indent << scprintf("taking step of size %f",tot) << endl;
                    
  function()->set_x(xnext);
  Ref<NonlinearTransform> t = function()->change_coordinates();
  apply_transform(t);

  // make the next gradient computed more accurate, since it will
  // be smaller
  accuracy_ = maxabs_gradient * maxabs_gradient_to_next_desired_accuracy;
  
  return converged;
}

void
EFCOpt::apply_transform(const Ref<NonlinearTransform> &t)
{
  if (t.null()) return;
  Optimize::apply_transform(t);
  if (last_mode_.nonnull()) t->transform_gradient(last_mode_);
  if (hessian_.nonnull()) t->transform_hessian(hessian_);
  if (update_.nonnull()) update_->apply_transform(t);
}

void
EFCOpt::print(std::ostream&o) const
{
  o << indent
    << "EFCOpt:"
    << std::endl
    << incindent
    << indent << "accuracy         = " << accuracy_
    << std::endl
    << indent << "transition_state = " << (tstate?"yes":"no")
    << std::endl
    << indent << "mode_following   = " << (modef?"yes":"no")
    << std::endl;

  if (update_.null()) {
    o << indent << "update           = 0 (hessian updates will not be performed)"
      << std::endl;
  }
  else {
    o << indent << "update           =" << std::endl;
    o << incindent;
    update_->print(o);
    o << decindent;
  }

  Optimize::print(o);
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
