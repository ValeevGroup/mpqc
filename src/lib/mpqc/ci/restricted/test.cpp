#define BOOST_TEST_MODULE CI Restricted

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "mpqc/ci/restricted/vector.hpp"

BOOST_AUTO_TEST_CASE(vector) {

    using namespace mpqc;

    size_t N = 3;
    Matrix matrices[N][N];
    std::vector<mpqc::range> spaces{range(0,1), range(1,10), range(10,20)};

    {
        int rank = 0;
        for (auto r : std::vector<mpqc::range>{range(0,1), range(1,10), range(10,20)}) {
        }
    }

    std::vector< ci::Space<ci::Alpha> > alpha;
    std::vector< ci::Space<ci::Beta> > beta;
    mpqc::matrix<bool> allowed(N,N);
    for (int j = 0; j < N; ++j) {
        alpha.push_back(ci::Space<ci::Alpha>(spaces.at(j), j));
        beta.push_back(ci::Space<ci::Beta>(spaces.at(j), j));
        for (int i = 0; i < N; ++i) {
            allowed(i,j) = (i < N-j);
        }
    }
    ci::SubspaceGrid grid(alpha, beta, allowed);

    ci::Vector<ci::Restricted> V("", grid);

    BOOST_TEST_MESSAGE("forbidden subspaces");
    {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                if (allowed(i,j)) continue;
                BOOST_TEST_MESSAGE("forbidden subspace (" << i << "," << j << ") should throw");
                range a = spaces[i];
                range b = spaces[j];
                matrices[i][j] = Matrix::Random(a.size(), b.size());
                BOOST_CHECK_THROW(V(a,b) = matrices[i][j], mpqc::Exception);
            }
        }
    }

    BOOST_TEST_MESSAGE("set/get");
    {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N-j; ++i) {
                range a = spaces[i];
                range b = spaces[j];
                matrices[i][j] = Matrix::Random(a.size(), b.size());
                V(a,b) = matrices[i][j];
            }
        }
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N-j; ++i) {
                range a = spaces[i];
                range b = spaces[j];
                BOOST_CHECK(Matrix(V(a,b)) == matrices[i][j]);
            }
        }
    }

}
