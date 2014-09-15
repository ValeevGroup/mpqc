#include "../matrix/full_rank_tile.h"
#include "../matrix/low_rank_tile.h"
#include "../matrix/tile_ops.h"
#include "create_low_rank_array.h"
#include "gtest.h"

template <typename T>
class TileOpsTest : public ::testing::Test {
  public:
    TileOpsTest() = default;
};

typedef ::testing::Types<float, double> TileTypes;
TYPED_TEST_CASE(TileOpsTest, TileTypes);

TYPED_TEST(TileOpsTest, FullSquareGemm) {
    const int rows = 10;
    const int cols = 10;
    using mat_type = typename FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    FullRankTile<TypeParam> left(mat_type::Random(rows, cols));
    FullRankTile<TypeParam> right(mat_type::Random(rows, cols));

    double alpha = 3.0;

    mat_type result = alpha * left.data() * right.data();
    auto tile = tile_ops::gemm(left, right, alpha);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.data().isApprox(result));
}

TYPED_TEST(TileOpsTest, FullMoreColsGemm) {
    const int rows = 3;
    const int cols = 10;
    const int inner_index = 7;
    using mat_type = typename FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    FullRankTile<TypeParam> left(mat_type::Random(rows, inner_index));
    FullRankTile<TypeParam> right(mat_type::Random(inner_index, cols));

    double alpha = 3.0;

    mat_type result = alpha * left.data() * right.data();
    auto tile = tile_ops::gemm(left, right, alpha);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.data().isApprox(result));
}

TYPED_TEST(TileOpsTest, FullMoreRowsGemm) {
    const int rows = 10;
    const int cols = 3;
    const int inner_index = 7;
    using mat_type = typename FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    FullRankTile<TypeParam> left(mat_type::Random(rows, inner_index));
    FullRankTile<TypeParam> right(mat_type::Random(inner_index, cols));

    double alpha = 3.0;

    mat_type result = alpha * left.data() * right.data();
    auto tile = tile_ops::gemm(left, right, alpha);

    EXPECT_EQ(std::min(rows, cols), tile.rank());
    EXPECT_EQ(rows * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.data().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowSquareGemmSameRank) {
    const int rows = 10;
    const int cols = 10;
    const int rank = 3;
    using mat_type = typename LowRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank);
    auto right = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank);

    double alpha = 3.0;

    mat_type result = alpha * left.matrixLR() * right.matrixLR();
    auto tile = tile_ops::gemm(left, right, alpha);

    EXPECT_EQ(rank, tile.rank());
    EXPECT_EQ(rows * rank + rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrixLR().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreColsGemmSameRank) {
    const int rows = 5;
    const int cols = 10;
    const int inner_index = 7;
    const int rank = 3;
    using mat_type = typename LowRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, inner_index, rank);
    auto right = TCC::test::low_rank_tile<TypeParam>(inner_index, cols, rank);

    double alpha = 3.0;

    mat_type result = alpha * left.matrixLR() * right.matrixLR();
    auto tile = tile_ops::gemm(left, right, alpha);

    EXPECT_EQ(rank, tile.rank());
    EXPECT_EQ(rows * rank + rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrixLR().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreRowsGemmSameRank) {
    const int rows = 10;
    const int cols = 5;
    const int inner_index = 7;
    const int rank = 3;
    using mat_type = typename LowRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, inner_index, rank);
    auto right = TCC::test::low_rank_tile<TypeParam>(inner_index, cols, rank);

    double alpha = 3.0;

    mat_type result = alpha * left.matrixLR() * right.matrixLR();
    auto tile = tile_ops::gemm(left, right, alpha);

    EXPECT_EQ(rank, tile.rank());
    EXPECT_EQ(rows * rank + rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrixLR().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowSquareGemmDiffRank) {
    const int rows = 10;
    const int cols = 10;
    const int rank_left = 3;
    const int rank_right = 4;
    using mat_type = typename LowRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(rows, cols, rank_right);

    double alpha = 3.0;

    mat_type result = alpha * left.matrixLR() * right.matrixLR();
    auto tile = tile_ops::gemm(left, right, alpha);

    auto out_rank = std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrixLR().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreColsGemmDiffRank) {
    const int rows = 5;
    const int cols = 10;
    const int k = 7;
    const int rank_left = 3;
    const int rank_right = 2;
    using mat_type = typename FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);

    double alpha = 3.0;

    mat_type result = alpha * left.matrixLR() * right.matrixLR();
    auto tile = tile_ops::gemm(left, right, alpha);

    auto out_rank = std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrixLR().isApprox(result));
}

TYPED_TEST(TileOpsTest, LowMoreRowsGemmDiffRank) {
    const int rows = 10;
    const int cols = 5;
    const int k = 7;
    const int rank_left = 3;
    const int rank_right = 2;
    using mat_type = typename FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    auto left = TCC::test::low_rank_tile<TypeParam>(rows, k, rank_left);
    auto right = TCC::test::low_rank_tile<TypeParam>(k, cols, rank_right);

    double alpha = 3.0;

    mat_type result = alpha * left.matrixLR() * right.matrixLR();
    auto tile = tile_ops::gemm(left, right, alpha);

    auto out_rank = std::min(rank_left, rank_right);

    EXPECT_EQ(out_rank, tile.rank());
    EXPECT_EQ(rows * out_rank + out_rank * cols, tile.size());
    EXPECT_EQ(cols, tile.Cols());
    EXPECT_EQ(rows, tile.Rows());
    EXPECT_TRUE(tile.matrixLR().isApprox(result));
}
