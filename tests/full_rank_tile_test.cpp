#include "../tensor/full_rank_tile.h"
#include "gtest.h"

template <typename T>
class FullRankTileTest : public ::testing::Test {
  public:
    FullRankTileTest() = default;
    tcc::tensor::FullRankTile<T> tile;
};

typedef ::testing::Types<int, long, float, double> FRTileTypes;
TYPED_TEST_CASE(FullRankTileTest, FRTileTypes);

TYPED_TEST(FullRankTileTest, DefaultConstructor) {
    EXPECT_EQ(0ul, this->tile.rank());
    EXPECT_EQ(0ul, this->tile.matrix().size());
    EXPECT_EQ(0ul, this->tile.matrix().cols());
    EXPECT_EQ(0ul, this->tile.matrix().rows());
}

TYPED_TEST(FullRankTileTest, AssignmentTest) {
    const int rows = 10;
    const int cols = 10;
    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;
    mat_type mat = mat_type::Random(rows, cols);
    tcc::tensor::FullRankTile<TypeParam> from_mat(mat);
    this->tile = from_mat;

    EXPECT_EQ(std::min(rows, cols), this->tile.rank());
    EXPECT_EQ(rows * cols, this->tile.matrix().size());
    EXPECT_EQ(cols, this->tile.matrix().cols());
    EXPECT_EQ(rows, this->tile.matrix().rows());
    EXPECT_EQ(this->tile.size(), this->tile.matrix().size());
    EXPECT_EQ(this->tile.Cols(), this->tile.matrix().cols());
    EXPECT_EQ(this->tile.Rows(), this->tile.matrix().rows());
    EXPECT_TRUE(this->tile.matrix().isApprox(mat));
}

TYPED_TEST(FullRankTileTest, CopyConstructor) {
    const int rows = 10;
    const int cols = 10;

    using mat_type = typename tcc::tensor::FullRankTile<TypeParam>::template Matrix
        <TypeParam>;

    mat_type mat = mat_type::Random(rows, cols);
    tcc::tensor::FullRankTile<TypeParam> from_mat(mat);
    this->tile = from_mat;

    // Actual copy
    tcc::tensor::FullRankTile<TypeParam> tile_copy(this->tile);

    EXPECT_EQ(tile_copy.rank(), this->tile.rank());
    EXPECT_EQ(tile_copy.size(), this->tile.size());
    EXPECT_EQ(tile_copy.Cols(), this->tile.Cols());
    EXPECT_EQ(tile_copy.Rows(), this->tile.Rows());
    EXPECT_TRUE(tile_copy.matrix().isApprox(this->tile.matrix()));
}
