static char rcsid[] = "$Header$";

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#define WANT_MATH
#define WANT_STREAM
#include <stdio.h>
#include "opt.h"

int mcsrch(NLP1* nlp, ColumnVector s, double &stp, int& nfev,
	   int maxfev, double ftol, double xtol, double gtol, 
	   double stpmax, double stpmin)
{
/****************************************************************************
 *   subroutine mcsrch
 *   Purpose
 *   find a step which satisfies 
 *   a sufficient decrease condition and a curvature condition. 
 *
 *   at each stage the subroutine updates an interval of 
 *   uncertainty with endpoints stx and sty. the interval of 
 *   uncertainty is initially chosen so that it contains a 
 *   minimizer of the modified function 
 *        f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s). 
 *
 *   if a step is obtained for which the modified function 
 *   has a nonpositive function value and nonnegative derivative, 
 *   then the interval of uncertainty is chosen so that it 
 *   contains a minimizer of f(x+stp*s). 
 *   the algorithm is designed to find a step which satisfies 
 *   the sufficient decrease condition 
 *         f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s), 
 *   and the curvature condition 
 *         abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s). 
 *   if ftol is less than gtol and if, for example, the function 
 *   is bounded below, then there is always a step which satisfies 
 *   both conditions. if no step can be found which satisfies both 
 *   conditions, then the algorithm usually stops when rounding 
 *   errors prevent further progress. in this case stp only 
 *   satisfies the sufficient decrease condition. 
 *
 *   Parameters
 *     n is a positive int input variable set to the number 
 *       of variables. 
 *     x is an array of length n. on input it must contain the 
 *       base point for the line search. on output it contains 
 *       x + stp*s. 
 *     f is a variable. on input it must contain the value of f 
 *       at x. on output it contains the value of f at x + stp*s. 
 *     g is an array of length n. on input it must contain the 
 *       gradient of f at x. on output it contains the gradient 
 *       of f at x + stp*s. 
 *     s is an input array of length n which specifies the 
 *       search direction. 
 *     stp is a nonnegative variable. on input stp contains an 
 *       initial estimate of a satisfactory step. on output 
 *       stp contains the final estimate. 
 *     ftol and gtol are nonnegative input variables. 
 *       termination occurs when the sufficient decrease 
 *       condition and the directional derivative condition are 
 *       satisfied. 
 *       ftol should be smaller than 5.e-1
 *       suggested value = 1.e-4 for newton methods
 *                       = 1.e-1 for more exact line searches
 *       Default Value = 1.e-4
 *
 *       gtol should be greater than 1.e-4 
 *       Default Value = 0.9
 *     xtol is a nonnegative input variable. termination occurs 
 *       when the relative width of the interval of uncertainty 
 *       is at most xtol. 
 *       Default Value = 2.2e-16
 *     stpmin and stpmax are nonnegative input variables which 
 *       specify lower and upper bounds for the step. (in this reverse 
 *       communication implementatin they are defined in a common 
 *       statement). 
 *       stpmin Default Value = 1.e-9
 *       stpmax Default Value = 1.e3
 *     maxfev is a positive int input variable. termination 
 *       occurs when the number of calls to fcn is at least 
 *       maxfev by the end of an iteration. 
 *     info is an integer output variable set as follows: 
 *       info =-1 improper input parameters. 
 *       info = 1  the sufficient decrease condition and the 
 *                 directional derivative condition hold. 
 *       info = 2  relative width of the interval of uncertainty 
 *                 is at most xtol. 
 *       info = 3  number of calls to fcn has reached maxfev. 
 *       info = 4  the step is at the lower bound stpmin. 
 *       info = 5  the step is at the upper bound stpmax. 
 *       info = 6  rounding errors prevent further progress. 
 *                 there may not be a step which satisfies the 
 *                 sufficient decrease and curvature conditions. 
 *                 tolerances may be too small. 
 *     nfev is an integer output variable set to the number of 
 *       calls to fcn. 
 *   subprograms called 
 *     mcstep 
 *     fortran-supplied...abs,max,min 
 *   argonne national laboratory. minpack project. june 1983 
 *   jorge j. more', david j. thuente 
 *
 *
 *   Recoded in C++ by Juan Meza December 1992
 *
 *
 *****************************************************************************/

  /* initialized data */

  static double half = .5;
  static double p66 = .66;
  static double xtrapf = 4.;
  static double zero = 0.;
  
  /* local variables */
  static double dgxm, dgym;
  static int j, info, infoc;
  static double finit, width, stmin, stmax;
  static long int stage1;
  static double width1, ftest1, dg, fm, fx, fy;
  static long int brackt;
  static double dginit, dgtest;
  static double  dgm, dgx, dgy, fxm, fym, stx, sty;

  int    siter;
  int    maxiter = 10;
  int    n = nlp->GetDim();

  double fvalue;
  ColumnVector work(n), xc(n), grad(n);
  
  infoc = 1;
  
  /*   check the input parameters for errors. */
  
  if (n <= 0 || stp <= zero || ftol < zero || gtol < zero ||
      xtol < zero || stpmin < zero || stpmax < stpmin || 
      maxfev <= 0) {
    return -1;
  }
  
  /* compute the initial gradient in the search direction */
  /* and check that s is a descent direction. */
  
  dginit = zero;
  grad = nlp->GetGrad();
  for (j = 1; j <= n; ++j) {
    dginit += grad(j) * s(j);
  }
  if (dginit >= zero) {
    printf("\nmcsrch: Initial search direction is not a descent direction\n");
    printf("gradient:\n");
    Print(grad);
    printf("displacement:\n");
    Print(s);
    return -1;
  }
  
  /* initialize local variables. */
  
  brackt = FALSE;
  stage1 = TRUE;
  nfev = 0;
  finit = nlp->GetF();
  dgtest = ftol * dginit;
  width = stpmax - stpmin;
  width1 = width / half;

  work = nlp->GetXc();

  /* the variables stx, fx, dgx contain the values of the step, */
  /* function, and directional derivative at the best step. */
  /* the variables sty, fy, dgy contain the value of the step, */
  /* function, and derivative at the other endpoint of */
  /* the interval of uncertainty. */
  /* the variables stp, f, dg contain the values of the step, */
  /* function, and derivative at the current step. */
  
  stx = zero;
  fx = finit;
  dgx = dginit;
  sty = zero;
  fy = finit;
  dgy = dginit;
  
  /* start of iteration. */
  
  for(siter=1; siter < maxiter; ++siter) {

    /*  set the minimum and maximum steps to correspond */
    /*  to the present interval of uncertainty. */
    
    if (brackt) {
      stmin = min(stx,sty);
      stmax = max(stx,sty);
    } else {
      stmin = stx;
      stmax = stp + xtrapf * (stp - stx);
    }
    
    /*    force the step to be within the bounds stpmax and stpmin. */
    
    stp = max(stp,stpmin);
    stp = min(stp,stpmax);
    
    /*    if an unusual termination is to occur then let */
    /*    stp be the lowest point obtained so far. */
    
    if (brackt && (stp <= stmin || stp >= stmax) || nfev >= maxfev - 1 || 
	infoc == 0 || brackt && stmax - stmin <= xtol * stmax) {
      stp = stx;
    }
    
    /*    evaluate the function and gradient at stp */
    /*    and compute the directional derivative. */
    
    xc = work + s*(stp); 
    nlp->SetX(xc);
    fvalue = nlp->EvalF(); 
    grad = nlp->EvalG(); 
    
    info = 0;
    ++(nfev);
    dg = zero;
    for (j = 1; j <= n; ++j) {
      dg += grad(j) * s(j);
    }
    ftest1 = finit + stp * dgtest;
    
    /*    test for convergence. */

    if (brackt && (stp <= stmin || stp >= stmax) || infoc == 0) {
      info = 6;
    }
    if (stp == stpmax && fvalue <= ftest1 && dg <= dgtest) {
      info = 5;
    }
    if (stp == stpmin && (fvalue > ftest1 || dg >= dgtest)) {
      info = 4;
    }
    if (nfev >= maxfev) {
      info = 3;
    }
    if (brackt && stmax - stmin <= xtol * stmax) {
      info = 2;
    }
    if (fvalue <= ftest1 && fabs(dg) <= gtol * (-dginit)) {
      info = 1;
    }
    
    /*    check for termination. */
    
    if (info != 0) {
      return info;
    }
    
    /*  in the first stage we seek a step for which the modified */
    /*  function has a nonpositive value and nonnegative derivative. */
    
    if (stage1 && fvalue <= ftest1 && dg >= min(ftol,gtol) * dginit) {
      stage1 = FALSE;
    }
    
    /*  a modified function is used to predict the step only if */
    /*  we have not obtained a step for which the modified */
    /*  function has a nonpositive function value and nonnegative */
    /*  derivative, and if a lower function value has been */
    /*  obtained but the decrease is not sufficient. */
    
    if (stage1 && fvalue <= fx && fvalue > ftest1) {
      
      /* define the modified function and derivative values. */
      
      fm = fvalue - stp * dgtest;
      fxm = fx - stx * dgtest;
      fym = fy - sty * dgtest;
      dgm = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;
      
      /* call cstep to update the interval of uncertainty */
      /* and to compute the new step. */
      
      mcstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm,
	     dgm, brackt, stmin, stmax, infoc);
      
      /* reset the function and gradient values for f. */
      
      fx = fxm + stx * dgtest;
      fy = fym + sty * dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;
    } else {
      
      /* call mcstep to update the interval of uncertainty */
      /* and to compute the new step. */
      
      mcstep(stx, fx, dgx, sty, fy, dgy, stp, fvalue,
	     dg, brackt, stmin, stmax, infoc);
    }
    
    /*  force a sufficient decrease in the size of the */
    /*  interval of uncertainty. */
    
    if (brackt) {
      if (fabs(sty - stx) >= p66 * width1) {
	stp = stx + (sty - stx)/2.;
      }
      width1 = width;
      width = fabs(sty - stx);
    }
    
  }    /*  end of iteration. */
  
  return info;

}  /* last line of subroutine mcsrch. */




int mcstep(double& stx, double& fx, double& dx, double& sty, double& fy, 
	   double& dy, double& stp, double& fp, double& dp, long int& brackt, 
	   double& stpmin, double& stpmax, int& info)
{
    /* system generated locals */
    double d_1, d_2, d_3;

    /* local variables */
    static double sgnd, stpc, stpf, stpq, p, q, gamma, r, s, theta;
    static long int bound;


/******************************************************************************
 *     subroutine mcstep 

 *     the purpose of mcstep is to compute a safeguarded step for 
 *     a linesearch and to update an interval of uncertainty for 
 *     a minimizer of the function. 

 *     the parameter stx contains the step with the least function 
 *     value. the parameter stp contains the current step. it is 
 *     assumed that the derivative at stx is negative in the 
 *     direction of the step. if brackt is set true then a 
 *     minimizer has been bracketed in an interval of uncertainty 
 *     with endpoints stx and sty. 

 *     the subroutine statement is 

 *       subroutine mcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt, 
 *                        stpmin,stpmax,info) 

 *     where 

 *       stx, fx, and dx are variables which specify the step, 
 *         the function, and the derivative at the best step obtained 
 *         so far. the derivative must be negative in the direction 
 *         of the step, that is, dx and stp-stx must have opposite 
 *         signs. on output these parameters are updated appropriately. 


 *       sty, fy, and dy are variables which specify the step, 
 *         the function, and the derivative at the other endpoint of 
 *         the interval of uncertainty. on output these parameters are 
 *         updated appropriately. 

 *       stp, fp, and dp are variables which specify the step, 
 *         the function, and the derivative at the current step. 
 *         if brackt is set true then on input stp must be 
 *         between stx and sty. on output stp is set to the new step. 

 *       brackt is a long int variable which specifies if a minimizer 
 *         has been bracketed. if the minimizer has not been bracketed 
 *         then on input brackt must be set false. if the minimizer 
 *         is bracketed then on output brackt is set true. 

 *       stpmin and stpmax are input variables which specify lower 
 *         and upper bounds for the step. 

 *       info is an int output variable set as follows: 
 *         if info = 1,2,3,4,5, then the step has been computed 
 *         according to one of the five cases below. otherwise 
 *         info = 0, and this indicates improper input parameters. 

 *     subprograms called 

 *       fortran-supplied ... abs,max,min,sqrt 

 *     argonne national laboratory. minpack project. june 1983 
 *     jorge j. more', david j. thuente 

 ******************************************************************************/

    info = 0;

    /* check the input parameters for errors. */

    if (brackt && 
	(stp <= min(stx,sty) || stp >= max(stx,sty)) || 
	dx * (stp - stx) >= 0. || stpmax < stpmin) {
	return 0;
      }

    /* determine if the derivatives have opposite sign. */
    
    sgnd = dp * (dx / fabs(dx));
    
    /* first case. a higher function value. */
    /* the minimum is bracketed. if the cubic step is closer */
    /* to stx than the quadratic step, the cubic step is taken, */
    /* else the average of the cubic and quadratic steps is taken. */
    
    if (fp > fx) {
      info = 1;
      bound = TRUE;
      theta = (fx - fp) * 3 / (stp - stx) + dx + dp;
      /* computing max */
      d_1 = fabs(theta), d_2 = fabs(dx), d_1 = max(d_2,d_1), d_2 = fabs(dp);
      s = max(d_2,d_1);
      /* computing 2nd power */
      d_1 = theta / s;
      gamma = s * sqrt(d_1 * d_1 - dx / s * (dp / s));
      if (stp < stx) {
	gamma = -gamma;
      }
      p = gamma - dx + theta;
      q = gamma - dx + gamma + dp;
      r = p / q;
      stpc = stx + r * (stp - stx);
      stpq = stx + dx / ((fx - fp) / (stp - stx) + dx) / 2 * (stp - 
								     stx);
      if ((d_1 = stpc - stx, fabs(d_1)) < (d_2 = stpq - stx, fabs(d_2))) {
	stpf = stpc;
      } else {
	stpf = stpc + (stpq - stpc) / 2;
      }
      brackt = TRUE;
      
      /* second case. a lower function value and derivatives of */
      /* opposite sign. the minimum is bracketed. if the cubic */
      /* step is closer to stx than the quadratic (secant) step, */
      /* the cubic step is taken, else the quadratic step is taken. */
      
    } else if (sgnd < 0.) {
	info = 2;
	bound = FALSE;
	theta = (fx - fp) * 3 / (stp - stx) + dx + dp;
	/* computing max */
	d_1 = fabs(theta), d_2 = fabs(dx), d_1 = max(d_2,d_1), d_2 = fabs(dp);
	s = max(d_2,d_1);
	/* computing 2nd power */
	d_1 = theta / s;
	gamma = s * sqrt(d_1 * d_1 - dx / s * (dp / s));
	if (stp > stx) {
	  gamma = -gamma;
	}
	p = gamma - dp + theta;
	q = gamma - dp + gamma + dx;
	r = p / q;
	stpc = stp + r * (stx - stp);
	stpq = stp + dp / (dp - dx) * (stx - stp);
	if ((d_1 = stpc - stp, fabs(d_1)) > (d_2 = stpq - stp, fabs(d_2))) {
	  stpf = stpc;
	} else {
	  stpf = stpq;
	}
	brackt = TRUE;
	
	/* third case. a lower function value, derivatives of the */
	/* same sign, and the magnitude of the derivative decreases. */
	/* the cubic step is only used if the cubic tends to infinity */
	/* in the direction of the step or if the minimum of the cubic */
	/* is beyond stp. otherwise the cubic step is defined to be */
	/* either stpmin or stpmax. the quadratic (secant) step is also */
	
	/* computed and if the minimum is bracketed then the the step */
	/* closest to stx is taken, else the step farthest away is taken. 
	 */
	
      } else if (fabs(dp) < fabs(dx)) {
	info = 3;
	bound = TRUE;
	theta = (fx - fp) * 3 / (stp - stx) + dx + dp;
	/* computing max */
	d_1 = fabs(theta), d_2 = fabs(dx), d_1 = max(d_2,d_1), d_2 = fabs(dp);
	s = max(d_2,d_1);
	
	/* the case gamma = 0 only arises if the cubic does not tend */
	/* to infinity in the direction of the step. */
	
	/* computing max */
	/* computing 2nd power */
	d_3 = theta / s;
	d_1 = 0., d_2 = d_3 * d_3 - dx / s * (dp / s);
	gamma = s * sqrt((max(d_2,d_1)));
	if (stp > stx) {
	  gamma = -gamma;
	}
	p = gamma - dp + theta;
	q = gamma + (dx - dp) + gamma;
	r = p / q;
	if (r < 0. && gamma != 0.) {
	  stpc = stp + r * (stx - stp);
	} else if (stp > stx) {
	  stpc = stpmax;
	} else {
	  stpc = stpmin;
	}
	stpq = stp + dp / (dp - dx) * (stx - stp);
	if (brackt) {
	  if ((d_1 = stp - stpc, fabs(d_1)) < (d_2 = stp - stpq, fabs(d_2)))
	    {
	      stpf = stpc;
	    } else {
	      stpf = stpq;
	    }
	} else {
	  if ((d_1 = stp - stpc, fabs(d_1)) > (d_2 = stp - stpq, fabs(d_2)))
	    {
	      stpf = stpc;
	    } else {
	      stpf = stpq;
	    }
	}
	
	/* fourth case. a lower function value, derivatives of the */
	/* same sign, and the magnitude of the derivative does */
	/* not decrease. if the minimum is not bracketed, the step */
	/* is either stpmin or stpmax, else the cubic step is taken. */
	
      } else {
	info = 4;
	bound = FALSE;
	if (brackt) {
	  theta = (fp - fy) * 3 / (sty - stp) + dy + dp;
	  /* computing max */
	  d_1 = fabs(theta), d_2 = fabs(dy), d_1 = max(d_2,d_1), d_2 = fabs(dp);
	  s = max(d_2,d_1);
	  /* computing 2nd power */
	  d_1 = theta / s;
	  gamma = s * sqrt(d_1 * d_1 - dy / s * (dp / s));
	  if (stp > sty) {
	    gamma = -gamma;
	  }
	  p = gamma - dp + theta;
	  q = gamma - dp + gamma + dy;
	  r = p / q;
	  stpc = stp + r * (sty - stp);
	  stpf = stpc;
	} else if (stp > stx) {
	  stpf = stpmax;
	} else {
	  stpf = stpmin;
	}
      }
    
    /* update the interval of uncertainty. this update does not */
    /* depend on the new step or the case analysis above. */
    
    if (fp > fx) {
      sty = stp;
      fy = fp;
      dy = dp;
    } else {
      if (sgnd < 0.) {
	sty = stx;
	fy = fx;
	dy = dx;
      }
      stx = stp;
      fx = fp;
      dx = dp;
    }
    
    /*  compute the new step and safeguard it. */
    
    stpf = min(stpmax,stpf);
    stpf = max(stpmin,stpf);
    stp = stpf;
    if (brackt && bound) {
      if (sty > stx) {
	/* computing max */
	d_1 = stx + (sty - stx) * (float).66;
	stp = min(stp,d_1);
      } else {
	/* computing max */
	d_1 = stx + (sty - stx) * (float).66;
	stp = max(stp,d_1);
      }
    }
    return 0;
    
  }
