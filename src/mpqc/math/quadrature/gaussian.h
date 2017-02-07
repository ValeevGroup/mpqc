/*
 * Created By Fabijan Pavosovic on 02/07/2017
 */
#ifndef MPQC_MATH_QUADRATURE_GAUSSIAN_H
#define MPQC_MATH_QUADRATURE_GAUSSIAN_H

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace math {

	/*! \brief this function uses orthogonal polynomial approach for obtaining
	 * quadrature roots and weights.  Taken from paper: Gaussian Quadrature
	 * and the Eigenvalue Problem by John A.  Gubner.
	 * http://gubner.ece.wisc.edu/gaussquad.pdf The code is outlined at
	 * Example 15: Legendre polynomials
	 * 
	 * \param N is the number of quadrature points 
	 * 
	 * \param w will return the weights 
	 * 
	 * \param x will return the roots
	 */ 
	void gauss_legendre(int N, Eigen::VectorXd &w, Eigen::VectorXd &x);

} // namespace math
} // namespace mpqc

#endif // MPQC_MATH_QUADRATURE_GAUSSIAN_H
