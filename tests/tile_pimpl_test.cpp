#include "../tensor/tile_pimpl.h"
#include "gtest.h"

template <typename T>
class TilePimplTest : public ::testing::Test {
  public:
    const int rows = 10;
    const int cols = 10;

    template <typename U>
    using Matrix = Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>;
    Matrix<T> matrix_;
    TileVariant<T> tile;
    TiledArray::Range range;
    TiledArray::math::GemmHelper gemm_helper;

    TilePimplTest()
        : matrix_(Matrix<T>::Random(rows, cols)), tile(matrix_), range(),
          gemm_helper(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                      madness::cblas::CBLAS_TRANSPOSE::NoTrans, rows, cols,
                      std::min(rows, cols)) {}
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
    EXPECT_TRUE(tile.empty());

    EXPECT_EQ(cut, tile.cut());
    EXPECT_EQ(tile.range().volume(), 0);
    EXPECT_EQ(tile.range(), this->range);
}

TYPED_TEST(TilePimplTest, TileVariantConstructor) {
    TilePimpl<TypeParam> tile(this->range, this->tile);
    EXPECT_FALSE(tile.empty());

    EXPECT_EQ(1e-07, tile.cut());
    EXPECT_EQ(tile.range(), this->range);
    EXPECT_EQ(tile.range().volume(), 0);

    EXPECT_EQ(10, tile.rank());
    EXPECT_TRUE(tile.isFull());
    EXPECT_TRUE(tile.tile().ftile().matrix().isApprox(this->matrix_));
}

TYPED_TEST(TilePimplTest, CopyConstructor) {
    TilePimpl<TypeParam> tile(this->range, this->tile);
    TilePimpl<TypeParam> copied = tile;

    EXPECT_FALSE(tile.empty());
    EXPECT_FALSE(copied.empty());

    EXPECT_EQ(1e-07, tile.cut());
    EXPECT_EQ(copied.cut(), tile.cut());
    EXPECT_EQ(copied.range(), tile.range());

    EXPECT_TRUE(tile.isFull());
    EXPECT_TRUE(copied.isFull());

    // Check ptr's should be same for copy
    EXPECT_TRUE(&(tile.tile()) == &(copied.tile()));

    EXPECT_TRUE(
        tile.tile().ftile().matrix().isApprox(copied.tile().ftile().matrix()));
    EXPECT_EQ(10, tile.rank());
    EXPECT_EQ(copied.rank(), tile.rank());
}

TYPED_TEST(TilePimplTest, TileClone) {
    TilePimpl<TypeParam> tile(this->range, this->tile);
    TilePimpl<TypeParam> cloned = tile.clone();
    EXPECT_FALSE(tile.empty());
    EXPECT_FALSE(cloned.empty());

    EXPECT_EQ(1e-07, tile.cut());
    EXPECT_EQ(cloned.cut(), tile.cut());
    EXPECT_EQ(cloned.range(), tile.range());

    EXPECT_TRUE(tile.isFull());
    EXPECT_TRUE(cloned.isFull());

    // Check ptr's should be different for clone
    EXPECT_FALSE(&(tile.tile()) == &(cloned.tile()));

    EXPECT_TRUE(
        tile.tile().ftile().matrix().isApprox(cloned.tile().ftile().matrix()));
    EXPECT_EQ(10, tile.rank());
    EXPECT_EQ(cloned.rank(), tile.rank());
}

TYPED_TEST(TilePimplTest, Gemm_AB) {
    double alpha = 2.0;
    TilePimpl<TypeParam> tileA(this->range, this->tile);
    TilePimpl<TypeParam> tileB(tileA);
    decltype(this->matrix_) correct_matrix = alpha * this->matrix_
                                             * this->matrix_;
    auto tileC = tileA.gemm(tileB, alpha, this->gemm_helper);

    EXPECT_FALSE(tileC.empty());

    EXPECT_EQ(1e-07, tileC.cut());
    EXPECT_TRUE(tileC.isFull());
    EXPECT_TRUE(tileC.tile().ftile().matrix().isApprox(correct_matrix));
}
