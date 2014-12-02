#include "../tensor/tile_variant.h"
#include "gtest.h"

template <typename T>
class TileVariantTest : public ::testing::Test{
public:
  using Lmat = typename tcc::tensor::LowRankTile<T>::template Matrix<T>;
  using Fmat = typename tcc::tensor::FullRankTile<T>::template Matrix<T>;

  const int rows = 10;
  const int cols = 10;
  const int rank = 3;

  Lmat L;
  Lmat R;

  tcc::tensor::LowRankTile<T> lr_tile;
  tcc::tensor::FullRankTile<T> f_tile;

  TileVariantTest() : L(Lmat::Random(rows,rank)), R(Lmat::Random(rank,cols)),
    lr_tile(L,R), f_tile(lr_tile.matrix()) {}
};

typedef ::testing::Types<float, double> TileVariantTypes;
TYPED_TEST_CASE(TileVariantTest, TileVariantTypes);

TYPED_TEST(TileVariantTest, DefaultConstructor){
  tcc::tensor::TileVariant<TypeParam> tile;
}

TYPED_TEST(TileVariantTest, LowRankTileConstructor){
  tcc::tensor::TileVariant<TypeParam> tile(this->lr_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile.tag());
  EXPECT_TRUE(this->lr_tile.matrix().isApprox(tile.lrtile().matrix()));
}

TYPED_TEST(TileVariantTest, LowRankTileMoveConstructor){
  tcc::tensor::TileVariant<TypeParam> tile(std::move(this->lr_tile));
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile.tag());
  EXPECT_EQ(this->lr_tile.size(),0);
  EXPECT_TRUE((this->L * this->R).isApprox(tile.lrtile().matrix()));
}

TYPED_TEST(TileVariantTest, FullRankTileConstructor){
  tcc::tensor::TileVariant<TypeParam> tile(this->f_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile.tag());
  EXPECT_TRUE(this->f_tile.matrix().isApprox(tile.ftile().matrix()));
}

TYPED_TEST(TileVariantTest, FullRankTileMoveConstructor){
  tcc::tensor::TileVariant<TypeParam> tile(std::move(this->f_tile));
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile.tag());
  EXPECT_EQ(this->f_tile.size(),0);
  EXPECT_TRUE((this->L * this->R).isApprox(tile.ftile().matrix()));
}

TYPED_TEST(TileVariantTest, FullRankTileToFullRankTileAssignment){
  tcc::tensor::TileVariant<TypeParam> tile(this->f_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile.tag());

  tcc::tensor::TileVariant<TypeParam> tile2(this->f_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile2.tag());

  tile = tile2;

  EXPECT_EQ(tile2.tag(),tile.tag());
  EXPECT_TRUE(tile2.matrix().isApprox(tile.matrix()));
}

TYPED_TEST(TileVariantTest, FullRankTileToFullRankTileMoveAssignment){
  tcc::tensor::TileVariant<TypeParam> tile(this->f_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile.tag());

  tcc::tensor::TileVariant<TypeParam> tile2(this->f_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile2.tag());

  auto tile3 = tile2;
  tile = std::move(tile2);

  EXPECT_EQ(tile3.tag(),tile.tag());
  EXPECT_TRUE(tile3.matrix().isApprox(tile.matrix()));
}

TYPED_TEST(TileVariantTest, FullRankTileToLowRankTileAssignment){
  tcc::tensor::TileVariant<TypeParam> tile(this->f_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile.tag());

  tcc::tensor::TileVariant<TypeParam> tile2(this->lr_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile2.tag());

  tile = tile2;

  EXPECT_EQ(tile2.tag(),tile.tag());
  EXPECT_TRUE(tile2.matrix().isApprox(tile.matrix()));
}

TYPED_TEST(TileVariantTest, FullRankTileToLowRankTileMoveAssignment){
  tcc::tensor::TileVariant<TypeParam> tile(this->f_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile.tag());

  tcc::tensor::TileVariant<TypeParam> tile2(this->lr_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile2.tag());

  auto tile3 = tile2;
  tile = std::move(tile2);

  EXPECT_EQ(tile3.tag(),tile.tag());
  EXPECT_TRUE(tile3.matrix().isApprox(tile.matrix()));
}

TYPED_TEST(TileVariantTest, LowRankTileToLowRankTileAssignment){
  tcc::tensor::TileVariant<TypeParam> tile(this->lr_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile.tag());

  tcc::tensor::TileVariant<TypeParam> tile2(this->lr_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile2.tag());

  tile = tile2;

  EXPECT_EQ(tile2.tag(),tile.tag());
  EXPECT_TRUE(tile2.matrix().isApprox(tile.matrix()));
}

TYPED_TEST(TileVariantTest, LowRankTileToLowRankTileMoveAssignment){
  tcc::tensor::TileVariant<TypeParam> tile(this->lr_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile.tag());

  tcc::tensor::TileVariant<TypeParam> tile2(this->lr_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile2.tag());

  auto tile3 = tile2;
  tile = std::move(tile2);

  EXPECT_EQ(tile3.tag(),tile.tag());
  EXPECT_TRUE(tile3.matrix().isApprox(tile.matrix()));
}

TYPED_TEST(TileVariantTest, LowRankTileToFullRankTileAssignment){
  tcc::tensor::TileVariant<TypeParam> tile(this->lr_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile.tag());

  tcc::tensor::TileVariant<TypeParam> tile2(this->f_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile2.tag());

  tile = tile2;

  EXPECT_EQ(tile2.tag(),tile.tag());
  EXPECT_TRUE(tile2.matrix().isApprox(tile.matrix()));
}

TYPED_TEST(TileVariantTest, LowRankTileToFullRankTileMoveAssignment){
  tcc::tensor::TileVariant<TypeParam> tile(this->lr_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::LowRank, tile.tag());

  tcc::tensor::TileVariant<TypeParam> tile2(this->f_tile);
  EXPECT_EQ(tcc::tensor::TileVariant<TypeParam>::TileType::FullRank, tile2.tag());

  auto tile3 = tile2;
  tile = std::move(tile2);

  EXPECT_EQ(tile3.tag(),tile.tag());
  EXPECT_TRUE(tile3.matrix().isApprox(tile.matrix()));
}

TYPED_TEST(TileVariantTest, LowRankFuncRank){
  tcc::tensor::TileVariant<TypeParam> tile(std::move(this->lr_tile));
  EXPECT_EQ(this->rank, tile.rank());
}

TYPED_TEST(TileVariantTest, FullRankFuncRank){
  tcc::tensor::TileVariant<TypeParam> tile(std::move(this->f_tile));
  EXPECT_EQ(this->cols, tile.rank());
}
