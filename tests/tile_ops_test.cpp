#include "../tensor/full_rank_tile.h"
#include "../tensor/low_rank_tile.h"
#include "../tensor/tile_ops.h"
#include "create_low_rank_array.h"
#include "gtest.h"

template <typename T>
class TileOpsTest : public ::testing::Test {
  public:
    TileOpsTest() = default;
};

typedef ::testing::Types<float, double> TileTypes;
TYPED_TEST_CASE(TileOpsTest, TileTypes);

TYPED_TEST(TileOpsTest, FullSquareAdd) {
    const int rows = 10;
    const int cols = 10;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::
        template Matrix<TypeParam>;

    tcc::tensor::FullRankTile<TypeParam> left(mat_type::Random(rows, cols));
    tcc::tensor::FullRankTile<TypeParam> right(mat_type::Random(rows, cols));

    double beta = 3.0;

    mat_type result = beta * left.matrix() + right.matrix();
    auto tile = tcc::tensor::binary_ops::add(left, right, beta);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, FullSkinnyAdd) {
    const int rows = 10;
    const int cols = 6;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::
        template Matrix<TypeParam>;

    tcc::tensor::FullRankTile<TypeParam> left_tall(
        mat_type::Random(rows, cols));
    tcc::tensor::FullRankTile<TypeParam> right_tall(
        mat_type::Random(rows, cols));

    tcc::tensor::FullRankTile<TypeParam> left_wide(
        mat_type::Random(cols, rows));
    tcc::tensor::FullRankTile<TypeParam> right_wide(
        mat_type::Random(cols, rows));

    double beta = 3.0;

    mat_type result_tall = beta * left_tall.matrix() + right_tall.matrix();
    mat_type result_wide = beta * left_wide.matrix() + right_wide.matrix();
    auto tile_tall = tcc::tensor::binary_ops::add(left_tall, right_tall, beta);
    auto tile_wide = tcc::tensor::binary_ops::add(left_wide, right_wide, beta);

    EXPECT_EQ(std::min(rows, cols), tile_tall.rank());
    EXPECT_EQ(rows * cols, tile_tall.size());
    EXPECT_EQ(cols, tile_tall.Cols());
    EXPECT_EQ(rows, tile_tall.Rows());
    EXPECT_TRUE(tile_tall.matrix().isApprox(result_tall));

    EXPECT_EQ(std::min(rows, cols), tile_wide.rank());
    EXPECT_EQ(rows * cols, tile_wide.size());
    EXPECT_EQ(cols, tile_wide.Rows()); // Reversed for wide
    EXPECT_EQ(rows, tile_wide.Cols());
    EXPECT_TRUE(tile_wide.matrix().isApprox(result_wide));
}

TYPED_TEST(TileOpsTest, FullSquareGemm) {
    const int rows = 10;
    const int cols = 10;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::
        template Matrix<TypeParam>;

    tcc::tensor::FullRankTile<TypeParam> left(mat_type::Random(rows, cols));
    tcc::tensor::FullRankTile<TypeParam> right(mat_type::Random(rows, cols));

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, FullMoreColsGemm) {
    const int rows = 3;
    const int cols = 10;
    const int inner_index = 7;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::
        template Matrix<TypeParam>;

    tcc::tensor::FullRankTile<TypeParam> left(
        mat_type::Random(rows, inner_index));
    tcc::tensor::FullRankTile<TypeParam> right(
        mat_type::Random(inner_index, cols));

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, FullMoreRowsGemm) {
    const int rows = 10;
    const int cols = 3;
    const int inner_index = 7;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::
        template Matrix<TypeParam>;

    tcc::tensor::FullRankTile<TypeParam> left(
        mat_type::Random(rows, inner_index));
    tcc::tensor::FullRankTile<TypeParam> right(
        mat_type::Random(inner_index, cols));

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowSquareAdd) {
    const int rows = 10;
    const int cols = 10;
    const int rankA = 2;
    const int rankB = 3;
    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, cols, rankA);
    auto right = TCC::test::low_rank_tile<TypeParam>(rows, cols, rankB);

    double beta = 3.0;

    mat_type result = beta * left.matrix() + right.matrix();
    auto tile = tcc::tensor::binary_ops::add(left, right, beta);

    EXPECT_EQ(rankA + rankB, tile.rank());
    EXPECT_EQ(rows * (rankA + rankB) + cols * (rankA + rankB), tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowSkinnyAdd) {
    const int rows = 10;
    const int cols = 6;
    const int rankA = 2;
    const int rankB = 1;
    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left_tall = TCC::test::low_rank_tile<TypeParam>(rows, cols, rankA);
    auto right_tall = TCC::test::low_rank_tile<TypeParam>(rows, cols, rankB);

    auto left_wide = TCC::test::low_rank_tile<TypeParam>(cols, rows, rankA);
    auto right_wide = TCC::test::low_rank_tile<TypeParam>(cols, rows, rankB);


    double beta = 3.0;

    mat_type result_tall = beta * left_tall.matrix() + right_tall.matrix();
    mat_type result_wide = beta * left_wide.matrix() + right_wide.matrix();
    auto tile_tall = tcc::tensor::binary_ops::add(left_tall, right_tall, beta);
    auto tile_wide = tcc::tensor::binary_ops::add(left_wide, right_wide, beta);

    EXPECT_EQ(rankA + rankB, tile_tall.rank());
    EXPECT_EQ(rows * (rankA + rankB) + cols * (rankA + rankB),
              tile_tall.size());
    EXPECT_EQ(cols, tile_tall.Cols());
    EXPECT_EQ(rows, tile_tall.Rows());
    EXPECT_TRUE(tile_tall.matrix().isApprox(result_tall));

    EXPECT_EQ(rankA + rankB, tile_wide.rank());
    EXPECT_EQ(rows * (rankA + rankB) + cols * (rankA + rankB),
              tile_wide.size());
    EXPECT_EQ(cols, tile_wide.Rows()); // Reversed for wide
    EXPECT_EQ(rows, tile_wide.Cols());
    EXPECT_TRUE(tile_wide.matrix().isApprox(result_wide));
}

TYPED_TEST(TileOpsTest, LowSquareGemmSameRank) {
    const int rows = 10;
    const int cols = 10;
    const int rank = 3;
    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank);
    auto right = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank);

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    EXPECT_EQ(rank, tile.rank());
    EXPECT_EQ(rows * rank + rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreColsGemmSameRank) {
    const int rows = 5;
    const int cols = 10;
    const int inner_index = 7;
    const int rank = 3;
    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, inner_index, rank);
    auto right = TCC::test::low_rank_tile<TypeParam>(inner_index, cols, rank);

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    EXPECT_EQ(rank, tile.rank());
    EXPECT_EQ(rows * rank + rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreRowsGemmSameRank) {
    const int rows = 10;
    const int cols = 5;
    const int inner_index = 7;
    const int rank = 3;
    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, inner_index, rank);
    auto right = TCC::test::low_rank_tile<TypeParam>(inner_index, cols, rank);

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    EXPECT_EQ(rank, tile.rank());
    EXPECT_EQ(rows * rank + rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowSquareGemmLeftLargerRank) {
    const int rows = 10;
    const int cols = 10;
    const int rank_left = 4;
    const int rank_right = 3;
    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_right);

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    auto out_rank = std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreColsGemmLeftLargerRank) {
    const int rows = 5;
    const int cols = 10;
    const int k = 7;
    const int rank_left = 3;
    const int rank_right = 2;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    auto out_rank = std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreRowsGemmLeftLargerRank) {
    const int rows = 10;
    const int cols = 5;
    const int k = 7;
    const int rank_left = 3;
    const int rank_right = 2;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    auto out_rank = std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}


TYPED_TEST(TileOpsTest, LowSquareGemmRightLargerRank) {
    const int rows = 10;
    const int cols = 10;
    const int rank_left = 3;
    const int rank_right = 4;
    using mat_type = typename tcc::tensor::LowRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_right);

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    auto out_rank = std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreColsGemmRightLargerRank) {
    const int rows = 5;
    const int cols = 10;
    const int k = 7;
    const int rank_left = 2;
    const int rank_right = 3;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    auto out_rank = std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreRowsGemmRightLargerRank) {
    const int rows = 10;
    const int cols = 5;
    const int k = 7;
    const int rank_left = 2;
    const int rank_right = 3;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::
        template Matrix<TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);

    double alpha = 3.0;

    mat_type result = alpha * left.matrix() * right.matrix();
    auto tile = tcc::tensor::binary_ops::gemm(left, right, alpha);

    auto out_rank = std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrix().isApprox(result));
}
