#include "../matrix/tile_variant.h"
#include "gtest.h"

template <typename T>
class TileVariantTest : public ::testing::Test{
public:
  using Lmat = typename LowRankTile<T>::template Matrix<T>;
  using Fmat = typename FullRankTile<T>::template Matrix<T>;

  const int rows = 10;
  const int cols = 10;
  const int rank = 3;

  Lmat L;
  Lmat R;

  LowRankTile<T> lr_tile;
  FullRankTile<T> f_tile;

  TileVariantTest() : L(Lmat::Random(rows,rank)), R(Lmat::Random(rank,cols)),
    lr_tile(L,R), f_tile(lr_tile.matrixLR()) {}
};

typedef ::testing::Types<float, double> TileVariantTypes;
TYPED_TEST_CASE(TileVariantTest, TileVariantTypes);

TYPED_TEST(TileVariantTest, DefaultConstructor){
  TileVariant<TypeParam> tile;
}

TYPED_TEST(TileVariantTest, LowRankTileConstructor){
  TileVariant<TypeParam> tile(this->lr_tile);
  EXPECT_EQ(TileVariant<TypeParam>::TileType::LowRank, tile.tag());
  EXPECT_TRUE(this->lr_tile.matrixLR().isApprox(tile.lrtile().matrixLR()));
}

TYPED_TEST(TileVariantTest, LowRankTileMoveConstructor){
  TileVariant<TypeParam> tile(std::move(this->lr_tile));
  EXPECT_EQ(TileVariant<TypeParam>::TileType::LowRank, tile.tag());
  EXPECT_EQ(this->lr_tile.size(),0);
  EXPECT_TRUE((this->L * this->R).isApprox(tile.lrtile().matrixLR()));
}

TYPED_TEST(TileVariantTest, FullRankTileConstructor){
  TileVariant<TypeParam> tile(this->f_tile);
  EXPECT_EQ(TileVariant<TypeParam>::TileType::FullRank, tile.tag());
  EXPECT_TRUE(this->f_tile.matrix().isApprox(tile.ftile().matrix()));
}

TYPED_TEST(TileVariantTest, FullRankTileMoveConstructor){
  TileVariant<TypeParam> tile(std::move(this->f_tile));
  EXPECT_EQ(TileVariant<TypeParam>::TileType::FullRank, tile.tag());
  EXPECT_EQ(this->f_tile.size(),0);
  EXPECT_TRUE((this->L * this->R).isApprox(tile.ftile().matrix()));
}

TYPED_TEST(TileVariantTest, LowRankFuncRank){
  TileVariant<TypeParam> tile(std::move(this->lr_tile));
  EXPECT_EQ(this->rank, tile.rank());
}

TYPED_TEST(TileVariantTest, FullRankFuncRank){
  TileVariant<TypeParam> tile(std::move(this->f_tile));
  EXPECT_EQ(this->cols, tile.rank());
}
