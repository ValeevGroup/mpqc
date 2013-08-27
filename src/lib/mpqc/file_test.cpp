#define BOOST_TEST_MODULE File

#include "mpqc/file.hpp"
#include "mpqc/math/matrix.hpp"

#include <boost/test/unit_test.hpp>
#include <vector>

BOOST_AUTO_TEST_CASE(FileTest) {

    using mpqc::File;
    using namespace mpqc;

    size_t m = 100, n = 10;
    std::vector<size_t> dims;
    dims.push_back(m);
    dims.push_back(n);

    BOOST_TEST_MESSAGE("File create");
    File file = File("file.h5");
    BOOST_REQUIRE(file);

    BOOST_TEST_MESSAGE("File::Dataset create");
    File::Dataset<double> ds(file, "my dataset", dims);
    BOOST_REQUIRE(ds);

    Matrix a = Matrix::Random(m,n);
    BOOST_TEST_MESSAGE("File::Dataset write");
    ds << a;

    Matrix b = Matrix::Random(m,n);
    BOOST_TEST_MESSAGE("File::Dataset read");
    ds >> b;
    BOOST_CHECK(a == b);

}


