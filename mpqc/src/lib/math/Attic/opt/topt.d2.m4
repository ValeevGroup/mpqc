Enter test case:Enter dimension:
Derivative Level?:Optimization Method ?
0 (Direct)       
1 (CG)           
2 (Quasi-Newton) 
3 (FD-Newton)    
4 (Newton)       
Number of Iterations?:
Rosenbrock function test case
main: test case 3:
NLP2: Initial Guess

    i	    xc 		 grad 
     1  -1.2000e+00  -2.1560e+02
     2   1.0000e+00  -8.8000e+01
f(x) =   2.4200e+01 

Tolerances
machine epsilon    = 2.220446e-16
maximum step       = 1.000000e+02
maximum iter       = 100
maximum fcn eval   = 250
step tolerance     = 2.220446e-16
function tolerance = 1.000000e-09
gradient tolerance = 6.055454e-06
Check_Deriv: Checking gradients versus finite-differences
    i    gradient     fd grad       error
    1  -2.1560e+02  -2.1560e+02   1.0506e-05
    2  -8.8000e+01  -8.8000e+01   1.4305e-06
maxerror =   1.0506e-05, tolerance =    1.3057e-03

Check_Deriv: Checking Hessian versus finite-differences

Matrix type: Sym   (2, 2)

1.90735e-05	3.97364e-07	
3.97364e-07	0	
maxerror =   1.9471e-05, tolerance =    1.3057e-03
main: Deriv O.K.

		Newton Method
Iter	  ||step||	 F(x) 	   ||grad|| 
     0   0.0000e+00   2.4200e+01   2.3287e+02
     1   1.0000e+00   4.7319e+00   4.6395e+00
     2   4.9677e-02   4.3319e+00   8.6375e+00
     3   2.1730e-01   3.8081e+00   1.6487e+01
     4   1.0000e+00   3.2714e+00   2.6603e+01
     5   1.0000e+00   2.4182e+00   6.0913e+00
     6   2.8506e-01   2.0761e+00   8.8933e+00
     7   1.0000e+00   1.6827e+00   1.2119e+01
     8   1.0000e+00   1.1886e+00   3.3586e+00
     9   3.7171e-01   9.7082e-01   4.8151e+00
    10   1.0000e+00   7.0629e-01   6.0165e+00
    11   1.0000e+00   4.6928e-01   2.5897e+00
    12   1.0000e+00   3.6105e-01   9.5472e+00
    13   1.0000e+00   1.8157e-01   7.1325e-01
    14   3.5467e-01   1.2342e-01   2.8246e+00
    15   1.0000e+00   6.6217e-02   4.2577e+00
    16   1.0000e+00   2.7490e-02   1.3369e+00
    17   1.0000e+00   1.1308e-02   3.1564e+00
    18   1.0000e+00   2.0083e-03   2.7065e-01
    19   1.0000e+00   2.4564e-04   6.3182e-01
    20   1.0000e+00   1.9134e-06   7.6330e-03
    21   1.0000e+00   3.7901e-10   8.1506e-04
    22   1.0000e+00   7.3530e-18   1.3083e-08
CheckConvg: deltaf =   3.7901e-10, ftol =   1.0000e-09

NLP2: Solution from newton

    i	    xc 		 grad 
     1   1.0000e+00   1.0430e-08
     2   1.0000e+00  -7.8978e-09
f(x) =   7.3530e-18 

Function tolerance test passed
Dimension of the problem = 2
Optimization method      = Newton
Return code = 2
Number of iterations taken    = 22

Matrix type: Sym   (2, 2)

All elements are zero
Number of matrices tested   = 2
Number of non-zero matrices = 1
