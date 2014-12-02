#include "../tensor/full_rank_tile.h"
#include "../tensor/low_rank_tile.h"
#include "../tensor/tile_mutations.h"
#include "create_low_rank_array.h"
#include "gtest.h"

template <typename T>
class TileMutationsTest : public ::testing::Test {
  public:
    TileMutationsTest() = default;
};

typedef ::testing::Types<float, double> TileTypes;
TYPED_TEST_CASE(TileMutationsTest, TileTypes);

TYPED_TEST(TileMutationsTest, FullSquareGemm) {
    const int rows = 10;
    const int cols = 10;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    tcc::tensor::FullRankTile<TypeParam> left(mat_type::Random(rows, cols));
    tcc::tensor::FullRankTile<TypeParam> right(mat_type::Random(rows, cols));
    tcc::tensor::FullRankTile<TypeParam> tile(mat_type::Random(rows, cols));

    double alpha = 3.0;
    double beta = 2.0;

    mat_type result = alpha * left.matrix() * right.matrix() + beta
                                                               * tile.matrix();
    tile = tcc::tensor::ternary_mutations::gemm(std::move(tile), left, right, alpha, beta);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, FullMoreColsGemm) {
    const int rows = 3;
    const int cols = 10;
    const int inner_index = 7;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    tcc::tensor::FullRankTile<TypeParam> left(mat_type::Random(rows, inner_index));
    tcc::tensor::FullRankTile<TypeParam> right(mat_type::Random(inner_index, cols));
    tcc::tensor::FullRankTile<TypeParam> tile(mat_type::Random(rows, cols));

    double alpha = 3.0;
    double beta = 2.0;

    mat_type result = alpha * left.matrix() * right.matrix() + beta
                                                               * tile.matrix();
    tile = tcc::tensor::ternary_mutations::gemm(std::move(tile), left, right, alpha, beta);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, FullMoreRowsGemm) {
    const int rows = 10;
    const int cols = 3;
    const int inner_index = 7;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    tcc::tensor::FullRankTile<TypeParam> left(mat_type::Random(rows, inner_index));
    tcc::tensor::FullRankTile<TypeParam> right(mat_type::Random(inner_index, cols));
    tcc::tensor::FullRankTile<TypeParam> tile(mat_type::Random(rows, cols));

    double alpha = 3.0;
    double beta = 2.0;

    mat_type result = alpha * left.matrix() * right.matrix() + beta * tile.matrix();
    tile = tcc::tensor::ternary_mutations::gemm(std::move(tile), left, right, alpha, beta);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, LowSquareGemmRightLargerRank) {
    const int rows = 10;
    const int cols = 10;
    const int rank_left = 3;
    const int rank_right = 4;
    const int tile_rank = 2;
    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_right);
    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, tile_rank);

    double alpha = 3.0;
    double beta = 2.0;

    mat_type result = alpha * left.matrix() * right.matrix() + beta * tile.matrix();
    tile = tcc::tensor::ternary_mutations::gemm(std::move(tile), left, right, alpha, beta);

    auto out_rank = tile_rank + std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, LowMoreColsGemmLeftLargerRank) {
    const int rows = 5;
    const int cols = 10;
    const int k = 7;
    const int rank_left = 4;
    const int rank_right = 3;
    const int tile_rank = 2;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);
    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, tile_rank);

    double alpha = 3.0;
    double beta = 2.0;

    mat_type result = alpha * left.matrix() * right.matrix() + beta * tile.matrix();
    tile = tcc::tensor::ternary_mutations::gemm(std::move(tile), left, right, alpha, beta);

    auto out_rank = tile_rank + std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, LowMoreRowsGemmABLargerRank) {
    const int rows = 10;
    const int cols = 5;
    const int k = 7;
    const int rank_left = 3;
    const int rank_right = 4;
    const int tile_rank = 2;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);
    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, tile_rank);

    double alpha = 3.0;
    double beta = 2.0;

    mat_type result = alpha * left.matrix() * right.matrix() + beta * tile.matrix();
    tile = tcc::tensor::ternary_mutations::gemm(std::move(tile), left, right, alpha, beta);

    auto out_rank = tile_rank + std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, LowSquareGemmResultLargerRankInPlace) {
    const int rows = 10;
    const int cols = 10;
    const int rank_left = 3;
    const int rank_right = 2;
    const int tile_rank = 4;
    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_right);
    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, tile_rank);

    double alpha = 3.0;
    double beta = 2.0;

    mat_type result = alpha * left.matrix() * right.matrix() + beta * tile.matrix();
    tile = tcc::tensor::ternary_mutations::gemm(std::move(tile), left, right, alpha, beta);

    auto out_rank = tile_rank + std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, LowMoreColsGemmResultLargerRankInPlace) {
    const int rows = 5;
    const int cols = 10;
    const int k = 7;
    const int rank_left = 2;
    const int rank_right = 3;
    const int tile_rank = 3;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);
    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, tile_rank);

    double alpha = 3.0;
    double beta = 2.0;

    mat_type result = alpha * left.matrix() * right.matrix() + beta * tile.matrix();
    tile = tcc::tensor::ternary_mutations::gemm(std::move(tile), left, right, alpha, beta);

    auto out_rank = tile_rank + std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, LowMoreRowsGemmResultLargerRank) {
    const int rows = 10;
    const int cols = 5;
    const int k = 7;
    const int rank_left = 3;
    const int rank_right = 2;
    const int tile_rank = 3;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);
    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, tile_rank);

    double alpha = 3.0;
    double beta = 2.0;

    mat_type result = alpha * left.matrix() * right.matrix() + beta * tile.matrix();
    tile = tcc::tensor::ternary_mutations::gemm(std::move(tile), left, right, alpha, beta);

    auto out_rank = tile_rank + std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, SubtLowLow) {
    const int rows = 10;
    const int cols = 10;

    const int rank_right = 2;
    const int tile_rank = 3;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto right = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_right);
    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, tile_rank);

    double beta = 2.0;

    mat_type result = tile.matrix() - beta * right.matrix();

    tile = tcc::tensor::binary_mutations::subt(std::move(tile), right, beta);

    auto out_rank = tile_rank + rank_right;

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}


TYPED_TEST(TileMutationsTest, LowCompressSquare) {
    const int rows = 10;
    const int cols = 10;
    const int rank = 3;

    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank);

    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::template Matrix
        <TypeParam>;
    mat_type result = tile.matrix();

    tcc::tensor::unary_mutations::compress(tile,1e-07);

    EXPECT_GE(rank, tile.rank());
    EXPECT_EQ(rows * rank + rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, LowCompressRightLarger) {
    const int rows = 5;
    const int cols = 10;
    const int rank = 3;

    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank);

    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::template Matrix
        <TypeParam>;
    mat_type result = tile.matrix();
    tcc::tensor::unary_mutations::compress(tile,1e-07);

    EXPECT_GE(rank, tile.rank());
    EXPECT_EQ(rows * rank + rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileMutationsTest, LowCompressLeftLarger) {
    const int rows = 10;
    const int cols = 5;
    const int rank = 3;

    auto tile = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank);

    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::template Matrix
        <TypeParam>;
    mat_type result = tile.matrix();
    tcc::tensor::unary_mutations::compress(tile, 1e-07);

    EXPECT_GE(rank, tile.rank());
    EXPECT_EQ(rows * rank + rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}
