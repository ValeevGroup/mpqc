#include "../include/tcc_all_same.h"
#include "gtest.h"

TEST(AllSameTest, AllSameTrue){
  bool same = tcc::all_same<double,double,double>::value;
  EXPECT_TRUE(same);
}

TEST(AllSameTest, AllSameFalse){
  bool same = tcc::all_same<double,float,double>::value;
  EXPECT_FALSE(same);
}


