static char rcsid[] = "$Header$";

#include <stdio.h>

#define WANT_MATH
#define WANT_STREAM
#include "opt.h"

/*-----------------------------------------------------------------*/
/*                                                                 */
/* test OPT                                                        */
/*                                                                 */
/* Copyright (C) 1993:                                             */
/* Juan C. Meza                                                    */
/* Center for Computational Engineering                            */
/* Sandia National Laboratories                                    */
/* (415)294-2425, (415)294-1016                                    */
/* email: meza@sandia.llnl.gov                                     */
/*-----------------------------------------------------------------*/

#define DEBUG
#define NMAX 10

main ()
{

  int i, n;
  int testcase, method;
  int deriv, maxiter;
  
  TOLS tol;

  srand48(1);

  ProgramException::SetAction(1);

  printf("Enter test case:");
  scanf("%d", &testcase);

  printf("Enter dimension:");
  scanf("%d",&n);

// Construct problem 

  ColumnVector x(n),  g(n);

  USERFCN2 tstf;


  Exception::PrintTrace(TRUE);

  printf("Derivative Level?:");
  scanf("%d",&deriv);

  printf( "Optimization Method ?\n");
  printf( "0 (Direct)       \n");
  printf( "1 (CG)           \n");
  printf( "2 (Quasi-Newton) \n");
  printf( "3 (FD-Newton)    \n");
  printf( "4 (Newton)       \n");
  scanf("%d", &method);
  
  printf("Number of Iterations?:");
  scanf("%d",&maxiter);

  switch (testcase) {

  case 1:
    printf( "\nQuadratic function test case\n");
    for (i=1; i<=n; ++i) {x(i) =  1.0;}
    tstf = &quad;
    break;
  case 2:
    printf( "\nDennis and Schnabel Ex 9.2.2, p 201\n");
    x(1) =  1.0;
    x(2) =  1.0;
    tstf = &ds_ex9;
    break;
  case 3:
    printf( "\nRosenbrock function test case\n");
    x(1) = -1.2;
    x(2)   =  1.0;
    tstf = &rosen;
    break;
  case 4:
    printf( "\nExtended Rosenbrock function test case\n");
    for (i=1; i<=n/2; ++i) {
      x(2*i-1) = -1.2;
      x(2*i)   =  1.0;
    }
    tstf = &erosen;
    break;
  case 5:
    printf( "\nWood test function\n");
    x(1) =  -3.0;
    x(2) =   1.0;
    x(3) =  -3.0;
    x(4) =   1.0;
    tstf = &wood;
    break;
  default:
    cerr << "Something is wrong\n";
    return -1;
  }

  printf("main: test case %d:",testcase);

  switch (method) {

  case 0:
// Direct Method

#ifdef LATER /* Thu May 20 16:42:16 1993 */
/**/    tol.DefaultTol();
/**/    tol.SetFtol(1.e-6);
/**/    tol.SetMaxIter(maxiter);
/**/    objfcn.OptPDS();
/**/    nlp.PrintState("Solution from nlcg");
/**/    objfcn.PrintStatus(); 
#endif /* LATER Thu May 20 16:42:16 1993 */
	break;

      case 1: // Nonlinear CG 
	{
	  NLF2 nlp(n,tstf);
	  nlp.SetX(x);
	  nlp.Eval();
	  nlp.PrintState("Initial Guess");

	  tol.SetDefaultTol();
	  tol.SetFtol(1.e-6);
	  tol.SetMaxIter(maxiter);
	  tol.PrintTol();

	  OptCG objfcn(&nlp,&tol);

	  objfcn.optimize();

	  nlp.PrintState("Solution from nlcg");
	  
	  objfcn.PrintStatus(); 
	}
	break;

      case 2:  // Quasi-Newton Method
	{
	  NLF2 nlp(n,tstf);
	  nlp.SetX(x);
	  nlp.Eval();
	  nlp.PrintState("Initial Guess");

	  tol.SetDefaultTol();
	  tol.SetFtol(1.e-9);
	  tol.SetMaxIter(maxiter);
	  tol.PrintTol();

	  OptQNewton objfcn(&nlp,&tol);

	  objfcn.optimize();

	  nlp.PrintState("Solution from quasi-newton");
	  objfcn.PrintStatus();
	}
	break;
      case 3:  // Finite-Difference Newton Method
	{
	  NLF2 nlp(n,tstf);
	  nlp.SetX(x);
	  nlp.Eval();
	  nlp.PrintState("Initial Guess");

	  tol.SetDefaultTol();
	  tol.SetFtol(1.e-9);
	  tol.SetMaxIter(maxiter);
	  tol.PrintTol();

          
          FDupdate fdu;
	  OptQNewton objfcn(&nlp,&tol,&fdu,0);

	  objfcn.optimize();

	  nlp.PrintState("Solution from finite-difference newton");
	  objfcn.PrintStatus();
	}
	break;
      case 4: // Newton Method
	{
	  NLF2 nlp(n,tstf);
	  nlp.SetX(x);
	  nlp.Eval();
	  nlp.PrintState("Initial Guess");

	  tol.SetDefaultTol();
	  tol.SetFtol(1.e-9);
	  tol.SetMaxIter(maxiter);
	  tol.PrintTol();

	  OptNewton objfcn(&nlp,&tol);

// Check Derivatives

	  int deriv_ok = objfcn.CheckDeriv();
	  if (deriv_ok) {
	    printf( "main: Deriv O.K.\n");
	  }
	  else {
	    printf("main: Derivatives may be wrong\n");
	  }

	  objfcn.optimize();

	  nlp.PrintState("Solution from newton");
	  
	  objfcn.PrintStatus();
	}
	break;
      default:
	break;
      }
  
#ifdef ATandT
  cout.flush();
#endif
#ifdef DEBUG_MALLOC
  malloc_dump(2);
#endif
}
