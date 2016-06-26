//
// Created by Chong Peng on 10/21/15.
//

#include "f12_utility.h"

#include <boost/algorithm/string.hpp>


namespace mpqc{
namespace f12{

    double basis_to_f12exponent(const std::string &basis_name) {

        std::string basis = basis_name;
        boost::trim(basis);
        boost::to_lower(basis);

        double f12_gamma;

        if(basis == "cc-pvdz-f12") {
            f12_gamma = 0.9;
        }
        else if(basis == "cc-pvtz-f12") {
            f12_gamma = 1.0;
        }
        else if(basis == "cc-pvqz-f12") {
            f12_gamma = 1.1;
        }
        else if(basis == "aug-cc-pvdz") {
            f12_gamma = 1.1;
        }
        else if(basis == "aug-cc-pvtz") {
            f12_gamma = 1.2;
        }
        else if(basis == "aug-cc-pvqz") {
            f12_gamma = 1.4;
        }
        else if(basis == "aug-cc-pv5z") {
            f12_gamma = 1.4;
        }
        else {
            f12_gamma = 1.2;
        }

        return f12_gamma;
    }

    namespace detail{

        double fstg(double zeta, double x){
            return -std::exp(-zeta*x)/zeta;
        }

        double fngtg(const std::vector<double>& cc,
                     const std::vector<double>& aa,
                     double x)
        {
            double value = 0.0;
            const double x2 = x*x;
            const auto n = cc.size();

            for (auto i = 0ul; i < n; ++i){
                value += cc[i] * std::exp(-aa[i]*x2);
            }

            return value;
        }

        double norm(const std::vector<double>& vec){
            double value = 0.0;
            const auto n = vec.size();
            for (auto i = 0ul; i < n; ++i){
                value += vec[i]*vec[i];
            }
            return value;
        }

        void LinearSolveDamped(int n, Eigen::MatrixXd& A, const double* b, double lambda, double* x) {
            Eigen::MatrixXd Acopy(A);
            for(int m=0; m<n; ++m) Acopy(m,m)  *= (1 + lambda);
            Eigen::VectorXd e(n);
            for( int i=0; i<n; i++ ) e(i) = b[i];
//            LINEQ_Gauss_Jordan( Acopy, e, 1);
            Eigen::VectorXd result = Acopy.colPivHouseholderQr().solve(e);
            for( int i=0; i<n; i++ ) x[i] = result(i);
        }

        // --- weighting functions ---
        // L2 error is weighted by ww(x)
        // hence error is weighted by sqrt(ww(x))
        double wwtewklopper(double x) {
            const double x2 = x * x;
            return x2 * std::exp(-2 * x2);
        }
        double wwcusp(double x) {
            const double x2 = x * x;
            const double x6 = x2 * x2 * x2;
            return std::exp(-0.005 * x6);
        }

        // default is Tew-Klopper
        double ww(double x) {
            return wwtewklopper(x);
        }

    }

    std::vector<std::pair<double,double>> stg_ng_fit(std::size_t n, double zeta) {
        const std::size_t nparams = 2* n;
        std::vector<double> cc(n,1.0);
        std::vector<double> aa(n);
        for(auto i=0; i<n; ++i){
            aa[i] = std::pow(3.0, (i + 2 - (n + 1)/2.0));
        }

        // first rescale cc for ff[x] to match the norm of f[x]
        double ffnormfac = 0.0;
        for(auto i=0; i<n; ++i){
            for(auto j=0; j<n; ++j){
                ffnormfac += cc[i] * cc[j]/std::sqrt(aa[i] + aa[j]);
            }
        }

        const double Nf = std::sqrt(2.0 * zeta) * zeta;
        const double Nff = std::sqrt(2.0) / (std::sqrt(ffnormfac) * std::sqrt(std::sqrt(M_PI)));
        for(auto i=0; i<n; ++i) cc[i] *= -Nff/Nf;

        const unsigned int npts = 1001;
        const double xmin = 0.0;
        const double xmax = 10.0;
        double lambda0 = 1000; // damping factor is initially set to 1000, eventually should end up at 0
        const double nu = 3.0; // increase/decrease the damping factor scale it by this
        const double epsilon = 1e-15;
        const unsigned int maxniter = 200;

        // grid points on which we will fit
        std::vector<double> xi(npts);
        for(unsigned int i=0; i<npts; ++i) xi[i] = xmin + (xmax - xmin)*i/(npts - 1);

        Eigen::VectorXd err(npts);
        Eigen::MatrixXd J(npts,nparams);

        std::vector<double> delta(nparams);

        double errnormI;
        double errnormIm1 = 1e3;
        bool converged = false;
        unsigned int iter = 0;
        while (!converged && iter < maxniter) {
            //std::cout << "Iteration " << ++iter << ": lambda = " << lambda0/nu << std::endl;

            for(unsigned int i=0; i<npts; ++i) {
                const double x = xi[i];
                err[i] = (detail::fstg(zeta, x) - detail::fngtg(cc, aa, x)) * std::sqrt(detail::ww(x));
            }
            errnormI = err.norm()/std::sqrt((double)npts);

            //std::cout << "|err|=" << errnormI << std::endl;
            converged = std::abs((errnormI - errnormIm1)/errnormIm1) <= epsilon;
            if (converged) break;
            errnormIm1 = errnormI;

            for(unsigned int i=0; i<npts; ++i) {
                const double x2 = xi[i] * xi[i];
                const double sqrt_ww_x = std::sqrt(detail::ww(xi[i]));
                for(auto j=0; j<n; ++j)
                    J(i,j) = (std::exp(-aa[j] * x2)) * sqrt_ww_x;
                for(auto j=0; j<n; ++j)
                    J(i,n+j) = -  sqrt_ww_x * x2 * cc[j] * std::exp(-aa[j] * x2);
            }

            Eigen::MatrixXd A(nparams,nparams);

            A = J.transpose()*J;
//            BLAS_Mat_x_Mat( A, J, J, true, false );

            Eigen::VectorXd b(nparams);

            b = J.transpose()*err;
//            BLAS_Mat_x_Vec( b, J, err, true );

            // try decreasing damping first
            // if not successful try increasing damping until it results in a decrease in the error
            lambda0 /= nu;
            for(int l=-1; l<1000; ++l) {
                detail::LinearSolveDamped(nparams, A, &(b[0]), lambda0, &(delta[0]) );
                std::vector<double> cc_0(cc); for(auto i=0; i<n; ++i) cc_0[i] += delta[i];
                std::vector<double> aa_0(aa); for(auto i=0; i<n; ++i) aa_0[i] += delta[i+n];

                // if any of the exponents are negative the step is too large and need to increase damping
                bool step_too_large = false;
                for(auto i=0; i<n; ++i)
                    if (aa_0[i] < 0.0) {
                        step_too_large = true;
                        break;
                    }
                if (!step_too_large) {
                    std::vector<double> err_0(npts);
                    for(unsigned int i=0; i<npts; ++i) {
                        const double x = xi[i];
                        err_0[i] = (detail::fstg(zeta, x) - detail::fngtg(cc_0, aa_0, x)) * std::sqrt(detail::ww(x));
                    }
                    const double errnorm_0 = detail::norm(err_0)/std::sqrt((double)npts);
                    if (errnorm_0 < errnormI) {
                        cc = cc_0;
                        aa = aa_0;
                        break;
                    }
                    else // step lead to increase of the error -- try dampening a bit more
                        lambda0 *= nu;
                }
                else // too large of a step
                    lambda0 *= nu;
            } // done adjusting the damping factor

        } // end of iterative minimization

        // if reached max # of iterations throw if the error is too terrible
        if (iter == maxniter && errnormI > 1e-10)
            throw std::runtime_error("after the max # of iterations Levenberg-Marquardt failed to converge to better than 1e-10");

        std::vector<std::pair<double,double>> result(n);
        for(auto i=0; i<n; ++i) {
            result[i].second = cc[i];
            result[i].first = aa[i];
        }
        return result;

    }


TiledArray::SparseShape<float> make_ijij_ijji_shape(const TiledArray::TiledRange& trange){


    // number of occ tiles
    auto n_occ =  trange.data()[0].tiles().second;

    TiledArray::Tensor<float> tile_norms(trange.tiles());

    auto max = std::numeric_limits<float>::max();

    // set sparse tile
    auto total = 0;
    for(auto i = 0; i < n_occ; i++){
        for(auto j=0; j < n_occ; j++){
            for(auto k=0; k < n_occ; k++){
                for(auto l=0; l < n_occ; l++, total++){
                    if(i==k && j==l){
                        tile_norms[total] = max;
                    }
                    else if(i==l && j==k){
                        tile_norms[total] = max;
                    }
                    else{
                        tile_norms[total] = 0.0;
                    }
                }
            }
        }
    }

    TiledArray::SparseShape<float> shape(tile_norms,trange);
    return shape;
}

}// end of f12 namespace
}// end of mpqc namespace
