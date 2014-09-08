#include "../matrix/full_rank_tile.h"
#include <gtest/gtest.h>

template <typename T>
class FullRankTileTest : public ::testing::Test {
  public:
    FullRankTileTest() = default;
    FullRankTile<T> tile;
};

typedef ::testing::Types<int, long, float, double> FRTileTypes;
TYPED_TEST_CASE(FullRankTileTest, FRTileTypes);

TYPED_TEST(FullRankTileTest, DefaultConstructor) {
    EXPECT_EQ(0ul, this->tile.rank());
    EXPECT_EQ(0ul, this->tile.data().size());
    EXPECT_EQ(0ul, this->tile.data().cols());
    EXPECT_EQ(0ul, this->tile.data().rows());
}

TYPED_TEST(FullRankTileTest, AssignmentTest) {
    const int rows = 10;
    const int cols = 10;
    using mat_type = typename FullRankTile<TypeParam>::template Matrix<TypeParam>;
    mat_type mat = mat_type::Random(rows, cols);
    FullRankTile<TypeParam> from_mat(mat);
    this->tile = from_mat;

    EXPECT_EQ(std::min(rows,cols), this->tile.rank());
    EXPECT_EQ(rows * cols, this->tile.data().size());
    EXPECT_EQ(cols, this->tile.data().cols());
    EXPECT_EQ(rows, this->tile.data().rows());
}
