#include "../matrix/tile_pimpl.h"
#include "gtest.h"

template <typename T>
class TilePimplTest : public ::testing::Test {
  public:
    TileVariant<T> tile;
    TiledArray::Range range;
    TiledArray::math::GemmHelper gemm_helper;

    TilePimplTest()
        : tile(FullRankTile
               <T>(Eigen::Matrix
                   <T, Eigen::Dynamic, Eigen::Dynamic>::Random(10, 10))),
          range(),
          gemm_helper(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                      madness::cblas::CBLAS_TRANSPOSE::NoTrans, 10, 10, 10) {}
};

typedef ::testing::Types<float, double> TilePimplTypes;
TYPED_TEST_CASE(TilePimplTest, TilePimplTypes);

TYPED_TEST(TilePimplTest, DefaultConstructor) {
    TilePimpl<TypeParam> tile;
    EXPECT_EQ(1e-07, tile.cut());
    EXPECT_EQ(tile.range().volume(), 0);
    EXPECT_TRUE(tile.empty());
}

TYPED_TEST(TilePimplTest, RangeConstructor) {
    TilePimpl<TypeParam> tile(this->range);
    EXPECT_EQ(1e-07, tile.cut());
    EXPECT_EQ(tile.range(), this->range);
    EXPECT_EQ(tile.range().volume(), 0);
    EXPECT_TRUE(tile.empty());
}

TYPED_TEST(TilePimplTest, RangeCutConstructor) {
    auto cut = 1e-04;
    TilePimpl<TypeParam> tile(this->range, cut);
    EXPECT_EQ(cut, tile.cut());
    EXPECT_EQ(tile.range(), this->range);
    EXPECT_EQ(tile.range().volume(), 0);
    EXPECT_TRUE(tile.empty());
}

TYPED_TEST(TilePimplTest, TileVariantConstructor) {
    TilePimpl<TypeParam> tile(this->range, this->tile);
    EXPECT_EQ(1e-07, tile.cut());
    EXPECT_EQ(tile.range(), this->range);
    EXPECT_EQ(tile.range().volume(), 0);
    EXPECT_TRUE(tile.isFull());
    EXPECT_EQ(10, tile.rank());
    EXPECT_FALSE(tile.empty());
}

TYPED_TEST(TilePimplTest, Gemm_AB) {
    double alpha = 2.0;
    TilePimpl<TypeParam> tileA(this->range, this->tile);
    TilePimpl<TypeParam> tileB(tileA);
    auto tileC = tileA.gemm(tileB, alpha, this->gemm_helper);

    EXPECT_EQ(1e-07, tileC.cut());
}
