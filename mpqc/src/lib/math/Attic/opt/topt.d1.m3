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

		Finite-Difference Newton Method
Iter	  ||step||	 F(x) 	   ||grad|| 
     0   0.0000e+00   2.4200e+01   2.3287e+02
     1   1.0000e+00   4.7319e+00   4.6395e+00
     2   4.9676e-02   4.3319e+00   8.6374e+00
     3   2.1730e-01   3.8081e+00   1.6487e+01
     4   1.0000e+00   3.2713e+00   2.6603e+01
     5   1.0000e+00   2.4182e+00   6.0914e+00
     6   2.8506e-01   2.0761e+00   8.8933e+00
     7   1.0000e+00   1.6827e+00   1.2119e+01
     8   1.0000e+00   1.1886e+00   3.3587e+00
     9   3.7172e-01   9.7082e-01   4.8152e+00
    10   1.0000e+00   7.0629e-01   6.0164e+00
    11   1.0000e+00   4.6928e-01   2.5898e+00
    12   1.0000e+00   3.6104e-01   9.5467e+00
    13   1.0000e+00   1.8157e-01   7.1332e-01
    14   3.5469e-01   1.2342e-01   2.8247e+00
    15   1.0000e+00   6.6217e-02   4.2574e+00
    16   1.0000e+00   2.7491e-02   1.3371e+00
    17   1.0000e+00   1.1307e-02   3.1559e+00
    18   1.0000e+00   2.0084e-03   2.7075e-01
    19   1.0000e+00   2.4564e-04   6.3179e-01
    20   1.0000e+00   1.9146e-06   7.6451e-03
    21   1.0000e+00   3.8513e-10   8.2214e-04
    22   1.0000e+00   2.2330e-14   6.6768e-06
CheckConvg: deltaf =   3.8511e-10, ftol =   1.0000e-09

NLP2: Solution from finite-difference newton

    i	    xc 		 grad 
     1   1.0000e+00   5.9708e-06
     2   1.0000e+00  -2.9882e-06
f(x) =   2.2330e-14 

Function tolerance test passed
Dimension of the problem = 2
Optimization method      = Finite-Difference Newton
Return code = 2
Number of iterations taken    = 22

Matrix type: Sym   (2, 2)

801.991	-399.998	
-399.998	200	
Number of matrices tested   = 1
Number of non-zero matrices = 1
