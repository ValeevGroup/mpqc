#include "../matrix/low_rank_tile.h"
#include "gtest.h"


template <typename T>
class LowRankTileTest : public ::testing::Test {
  public:
    LowRankTileTest() = default;
    LowRankTile<T> tile;
};

typedef ::testing::Types<int, long, float, double> LRTileTypes;
TYPED_TEST_CASE(LowRankTileTest, LRTileTypes);

TYPED_TEST(LowRankTileTest, DefaultConstructor) {
    EXPECT_EQ(0ul, this->tile.rank());
    EXPECT_EQ(0ul, this->tile.Rows());
    EXPECT_EQ(0ul, this->tile.Cols());
    EXPECT_EQ(0ul, this->tile.size());
    EXPECT_EQ(this->tile.matrixL().rows(), this->tile.Rows());
    EXPECT_EQ(this->tile.matrixR().cols(), this->tile.Cols());
    EXPECT_EQ(this->tile.matrixL().size() + this->tile.matrixR().size(),
              this->tile.size());
}
