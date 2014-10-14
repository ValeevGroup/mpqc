#include "../matrix/tile_algebra.h"
#include "../include/eigen.h"
#include "create_low_rank_array.h"
#include "gtest.h"

using namespace algebra;
using namespace Eigen;

TEST(EigenTileAlgebraTest, CreateCMatrixSquareGEMM) {
    const auto rows = 10;
    const auto cols = 10;
    const auto inner_dim = 10;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = eigen_version::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, CreateCMatrixMoreColsGEMM) {
    const auto rows = 4;
    const auto cols = 11;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = eigen_version::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, CreateCMatrixMoreRowsGEMM) {
    const auto rows = 11;
    const auto cols = 4;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = eigen_version::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, CreateCMatrixBisColVectorGEMM) {
    const auto rows = 11;
    const auto cols = 1;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = eigen_version::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, CreateCMatrixAisRowVectorGEMM) {
    const auto rows = 1;
    const auto cols = 11;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = eigen_version::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, CreateCMatrixInnerDimIsOneGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 1;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = eigen_version::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, CreateCMatrixInnerDimIsZeroGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 0;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = eigen_version::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}


TEST(EigenTileAlgebraTest, SquareInPlaceGEMM) {
    const auto rows = 11;
    const auto cols = 11;
    const auto inner_dim = 11;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    eigen_version::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, MoreColsThanRowsInPlaceGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    eigen_version::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, MoreRowsThanColsInPlaceGEMM) {
    const auto rows = 11;
    const auto cols = 9;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    eigen_version::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, AisRowVectorInPlaceGEMM) {
    const auto rows = 1;
    const auto cols = 11;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    eigen_version::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, BisColVectorInPlaceGEMM) {
    const auto rows = 11;
    const auto cols = 1;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    eigen_version::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, InnerDimIsOneInPlaceGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 1;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    eigen_version::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, InnerDimIsZeroInPlaceGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 0;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    eigen_version::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(EigenTileAlgebraTest, FullRankColPivQR) {
    const auto rows = 9;
    const auto cols = 11;
    const auto rank = std::min(rows, cols);
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_TRUE(eigen_version::Decompose_Matrix(C, L, R, cut))
        << "Matrix C = \n" << C << "\n";
    EXPECT_EQ(L.size(), 0ul); // Shouldn't decompose if full rank
    EXPECT_EQ(R.size(), 0ul);
}

TEST(EigenTileAlgebraTest, LowRankSquareColPivQR) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_FALSE(eigen_version::Decompose_Matrix(C, L, R, cut))
        << "Matrix C = \n" << C << "\n";
    EXPECT_EQ(L.cols(), rank); // Shouldn't decompose if full rank
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R))
        << "2 norm of diff = " << (C - (L * R)).lpNorm<2>() << "\nC = \n" << C
        << "\nL = \n" << L << "\nR = \n" << R << "\n";
}

TEST(EigenTileAlgebraTest, LowRankMoreRowsThanColsColPivQR) {
    const auto rows = 13;
    const auto cols = 9;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_FALSE(eigen_version::Decompose_Matrix(C, L, R, cut))
        << "Matrix C = \n" << C << "\n";
    EXPECT_EQ(L.cols(), rank); // Shouldn't decompose if full rank
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R))
        << "2 norm of diff = " << (C - (L * R)).lpNorm<2>() << "\nC = \n" << C
        << "\nL = \n" << L << "\nR = \n" << R << "\n";
}

TEST(EigenTileAlgebraTest, LowRankMoreColsThanRowsColPivQR) {
    const auto rows = 9;
    const auto cols = 13;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_FALSE(eigen_version::Decompose_Matrix(C, L, R, cut))
        << "Matrix C = \n" << C << "\n";
    EXPECT_EQ(L.cols(), rank); // Shouldn't decompose if full rank
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R))
        << "2 norm of diff = " << (C - (L * R)).lpNorm<2>() << "\nC = \n" << C
        << "\nL = \n" << L << "\nR = \n" << R << "\n";
}

TEST(EigenTileAlgebraTest, RankOneColPivQR) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 1;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_FALSE(eigen_version::Decompose_Matrix(C, L, R, cut))
        << "Matrix C = \n" << C << "\n";
    EXPECT_EQ(L.cols(), rank); // Shouldn't decompose if full rank
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R))
        << "2 norm of diff = " << (C - (L * R)).lpNorm<2>() << "\nC = \n" << C
        << "\nL = \n" << L << "\nR = \n" << R << "\n";
}


/*****************************************************
 * BEGIN LAPACK TESTS
 *****************************************************/

TEST(LapackTileAlgebraTest, CreateCMatrixSquareGEMM) {
    const auto rows = 10;
    const auto cols = 10;
    const auto inner_dim = 10;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = lapack::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, CreateCMatrixMoreColsGEMM) {
    const auto rows = 4;
    const auto cols = 11;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = lapack::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, CreateCMatrixMoreRowsGEMM) {
    const auto rows = 11;
    const auto cols = 4;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = lapack::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, CreateCMatrixBisColVectorGEMM) {
    const auto rows = 11;
    const auto cols = 1;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = lapack::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, CreateCMatrixAisRowVectorGEMM) {
    const auto rows = 1;
    const auto cols = 11;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = lapack::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, CreateCMatrixInnerDimIsOneGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 1;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = lapack::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, CreateCMatrixInnerDimIsZeroGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 0;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);

    auto alpha = 2.0;

    auto C = lapack::cblas_gemm(A, B, alpha);

    auto Ccorrect = alpha * A * B;

    EXPECT_TRUE(C.isApprox(Ccorrect))
        << "C = \n" << C << "\nCorrect Answer = \n" << Ccorrect << "\n";
}
TEST(LapackTileAlgebraTest, SquareInPlaceGEMM) {
    const auto rows = 11;
    const auto cols = 11;
    const auto inner_dim = 11;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    lapack::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, MoreColsThanRowsInPlaceGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    lapack::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, MoreRowsThanColsInPlaceGEMM) {
    const auto rows = 11;
    const auto cols = 9;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    lapack::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, AisRowVectorInPlaceGEMM) {
    const auto rows = 1;
    const auto cols = 11;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    lapack::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, BisColVectorInPlaceGEMM) {
    const auto rows = 11;
    const auto cols = 1;
    const auto inner_dim = 7;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    lapack::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;


    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, InnerDimIsOneInPlaceGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 1;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    lapack::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, InnerDimIsZeroInPlaceGEMM) {
    const auto rows = 9;
    const auto cols = 11;
    const auto inner_dim = 0;
    MatrixXd A = MatrixXd::Random(rows, inner_dim);
    MatrixXd B = MatrixXd::Random(inner_dim, cols);
    MatrixXd C = MatrixXd::Random(rows, cols);
    MatrixXd Ccorrect = C;

    auto alpha = 2.0;
    auto beta = 3.0;

    lapack::cblas_gemm_inplace(A, B, C, alpha, beta);

    Ccorrect = alpha * A * B + beta * Ccorrect;

    EXPECT_TRUE(C.isApprox(Ccorrect));
}

TEST(LapackTileAlgebraTest, FullRankColPivQR) {
    const auto rows = 9;
    const auto cols = 11;
    const auto rank = std::min(rows, cols);
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_TRUE(lapack::Decompose_Matrix(C, L, R, cut)) << "Matrix C = \n" << C
                                                        << "\n";
    EXPECT_EQ(L.size(), 0ul); // Shouldn't decompose if full rank
    EXPECT_EQ(R.size(), 0ul);
}

TEST(LapackTileAlgebraTest, LowRankSquareColPivQR) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_FALSE(lapack::Decompose_Matrix(C, L, R, cut)) << "Matrix C = \n" << C
                                                         << "\n";
    EXPECT_EQ(L.cols(), rank); // Shouldn't decompose if full rank
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R))
        << "2 norm of diff = " << (C - (L * R)).lpNorm<2>() << "\nC = \n" << C
        << "\nL = \n" << L << "\nR = \n" << R << "\n";
}

TEST(LapackTileAlgebraTest, LowRankMoreRowsThanColsColPivQR) {
    const auto rows = 13;
    const auto cols = 9;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_FALSE(lapack::Decompose_Matrix(C, L, R, cut)) << "Matrix C = \n" << C
                                                         << "\n";
    EXPECT_EQ(L.cols(), rank); // Shouldn't decompose if full rank
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R))
        << "2 norm of diff = " << (C - (L * R)).lpNorm<2>() << "\nC = \n" << C
        << "\nL = \n" << L << "\nR = \n" << R << "\n";
}

TEST(LapackTileAlgebraTest, LowRankMoreColsThanRowsColPivQR) {
    const auto rows = 9;
    const auto cols = 13;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_FALSE(lapack::Decompose_Matrix(C, L, R, cut)) << "Matrix C = \n" << C
                                                         << "\n";
    EXPECT_EQ(L.cols(), rank); // Shouldn't decompose if full rank
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R))
        << "2 norm of diff = " << (C - (L * R)).lpNorm<2>() << "\nC = \n" << C
        << "\nL = \n" << L << "\nR = \n" << R << "\n";
}

TEST(LapackTileAlgebraTest, RankOneColPivQR) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 1;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_FALSE(lapack::Decompose_Matrix(C, L, R, cut)) << "Matrix C = \n" << C
                                                         << "\n";
    EXPECT_EQ(L.cols(), rank); // Shouldn't decompose if full rank
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R))
        << "2 norm of diff = " << (C - (L * R)).lpNorm<2>() << "\nC = \n" << C
        << "\nL = \n" << L << "\nR = \n" << R << "\n";
}

TEST(LapackTileAlgebraTest, RankZeroColPivQR) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 0;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    EXPECT_FALSE(lapack::Decompose_Matrix(C, L, R, cut)) << "Matrix C = \n" << C
                                                         << "\n";
    EXPECT_EQ(L.cols(), rank); // Shouldn't decompose if full rank
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R))
        << "2 norm of diff = " << (C - (L * R)).lpNorm<2>() << "\nC = \n" << C
        << "\nL = \n" << L << "\nR = \n" << R << "\n";
}

TEST(LapackTileAlgebraTest, FullRankColPivotedQr) {
    const auto rows = 9;
    const auto cols = 11;
    const auto rank = std::min(rows, cols);
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    algebra::ColPivotedQr(C, L, R, cut);
    EXPECT_EQ(L.cols(), rank);
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, LowRankSquareQRInit) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    EXPECT_EQ(L.cols(), rank);
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, LowRankMoreRowsThanColsColPivotedQr) {
    const auto rows = 13;
    const auto cols = 9;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    EXPECT_EQ(L.cols(), rank);
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, LowRankMoreColsThanRowsColPivotedQr) {
    const auto rows = 9;
    const auto cols = 13;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    EXPECT_EQ(L.cols(), rank);
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, RankOneColPivotedQr) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 1;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    EXPECT_EQ(L.cols(), rank);
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, RankZeroColPivotedQr) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 0;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    EXPECT_EQ(L.cols(), rank);
    EXPECT_EQ(R.rows(), rank);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, FullRankCompressLeft) {
    const auto rows = 9;
    const auto cols = 11;
    const auto rank = std::min(rows, cols);
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressLeft(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, LowRankSquareCompressLeft) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressLeft(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, LowRankMoreRowsThanColsCompressLeft) {
    const auto rows = 13;
    const auto cols = 9;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressLeft(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, LowRankMoreColsThanRowsCompressLeft) {
    const auto rows = 9;
    const auto cols = 13;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressLeft(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, RankOneCompressLeft) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 1;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressLeft(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, FullRankCompressRight) {
    const auto rows = 9;
    const auto cols = 11;
    const auto rank = std::min(rows, cols);
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressRight(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, LowRankSquareCompressRight) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    // Return true if mat is full rank
    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressRight(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, LowRankMoreRowsThanColsCompressRight) {
    const auto rows = 13;
    const auto cols = 9;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressRight(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, LowRankMoreColsThanRowsCompressRight) {
    const auto rows = 9;
    const auto cols = 13;
    const auto rank = 3;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressRight(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}

TEST(LapackTileAlgebraTest, RankOneCompressRight) {
    const auto rows = 10;
    const auto cols = 10;
    const auto rank = 1;
    const auto cut = 1e-07;

    MatrixXd C = TCC::test::low_rank_matrix<double>(rows, cols, rank);
    MatrixXd L, R;

    algebra::ColPivotedQr(C, L, R, cut);
    algebra::CompressRight(L, R, cut, true);
    EXPECT_TRUE(C.isApprox(L * R));
}
