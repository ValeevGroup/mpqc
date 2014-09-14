#include "../matrix/tile_concept.h"
#include "gtest.h"

TEST(TileConceptTest, DefaultConstructor){
  EXPECT_NO_THROW(TileAble(3.0));
  TileAble tile(3.0);
  EXPECT_TRUE((tile.type() == typeid(double)));
  EXPECT_EQ(3.0, *(tile.cast_data<double>()));
}



