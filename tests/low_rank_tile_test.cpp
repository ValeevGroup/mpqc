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

TYPED_TEST(LowRankTileTest, LRConstructor) {
    const int rows = 11;
    const int cols = 9;
    const int rank = 3;
    auto L = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rows, rank).eval();
    auto R = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rank, cols).eval();
    auto mat = (L * R).eval();

    LowRankTile<TypeParam> lr_tile(L, R);

    EXPECT_EQ(mat.rows(), lr_tile.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(mat.cols(), lr_tile.Cols()) << "Cols were not equal\n";
    EXPECT_EQ(rows, lr_tile.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(cols, lr_tile.Cols()) << "Cols were not equal\n";
    EXPECT_EQ(rank, lr_tile.rank()) << "Rank was not equal\n";
    EXPECT_EQ(L.size() + R.size(), lr_tile.size()) << "Size was not equal";

    EXPECT_TRUE(L.isApprox(lr_tile.matrixL()))
        << "L was not equal to matrixL()\nL = \n" << L << "\nmatrixL() = \n"
        << lr_tile.matrixL() << "\n";

    EXPECT_TRUE(R.isApprox(lr_tile.matrixR()))
        << "R was not equal to matrixR()\nR = \n" << R << "\nmatrixR() = \n"
        << lr_tile.matrixR() << "\n";

    EXPECT_TRUE(mat.isApprox(lr_tile.matrixLR()))
        << "L * R was not equal to matrixLR()\nL * R = \n" << mat
        << "\nmatrixLR() = \n" << lr_tile.matrixLR() << "\n";
}

TYPED_TEST(LowRankTileTest, LRMoveConstructor) {
    const int rows = 11;
    const int cols = 9;
    const int rank = 3;
    auto L = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rows, rank).eval();
    auto R = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rank, cols).eval();
    auto L_moved_from = L;
    auto R_moved_from = R;
    auto mat = (L * R).eval();

    LowRankTile
        <TypeParam> lr_tile(std::move(L_moved_from), std::move(R_moved_from));

    EXPECT_EQ(0ul, L_moved_from.rows());
    EXPECT_EQ(0ul, L_moved_from.cols());
    EXPECT_EQ(0ul, L_moved_from.size());
    EXPECT_EQ(0ul, R_moved_from.rows());
    EXPECT_EQ(0ul, R_moved_from.cols());
    EXPECT_EQ(0ul, R_moved_from.size());

    EXPECT_EQ(mat.rows(), lr_tile.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(mat.cols(), lr_tile.Cols()) << "Cols were not equal\n";
    EXPECT_EQ(rows, lr_tile.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(cols, lr_tile.Cols()) << "Cols were not equal\n";
    EXPECT_EQ(rank, lr_tile.rank()) << "Rank was not equal\n";
    EXPECT_EQ(L.size() + R.size(), lr_tile.size()) << "Size was not equal";

    EXPECT_TRUE(L.isApprox(lr_tile.matrixL()))
        << "L was not equal to matrixL()\nL = \n" << L << "\nmatrixL() = \n"
        << lr_tile.matrixL() << "\n";

    EXPECT_TRUE(R.isApprox(lr_tile.matrixR()))
        << "R was not equal to matrixR()\nR = \n" << R << "\nmatrixR() = \n"
        << lr_tile.matrixR() << "\n";

    EXPECT_TRUE(mat.isApprox(lr_tile.matrixLR()))
        << "L * R was not equal to matrixLR()\nL * R = \n" << mat
        << "\nmatrixLR() = \n" << lr_tile.matrixLR() << "\n";
}


TYPED_TEST(LowRankTileTest, CopyConstructor) {
    const int rows = 9;
    const int cols = 11;
    const int rank = 3;

    auto L = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rows, rank).eval();
    auto R = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rank, cols).eval();
    auto mat = (L * R).eval();

    LowRankTile<TypeParam> lr_tile(L, R);
    LowRankTile<TypeParam> lr_tile_copy(lr_tile);

    EXPECT_EQ(rows, lr_tile_copy.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(cols, lr_tile_copy.Cols()) << "Cols were not equal\n";
    EXPECT_EQ(rank, lr_tile_copy.rank()) << "Rank was not equal\n";
    EXPECT_EQ(lr_tile_copy.size(), lr_tile.size()) << "Size was not equal";
    EXPECT_EQ(lr_tile_copy.Rows(), lr_tile.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(lr_tile_copy.Cols(), lr_tile.Cols()) << "Cols were not equal\n";

    EXPECT_TRUE(lr_tile_copy.matrixL().isApprox(lr_tile.matrixL()));
    EXPECT_TRUE(lr_tile_copy.matrixR().isApprox(lr_tile.matrixR()));
    EXPECT_TRUE(lr_tile_copy.matrixLR().isApprox(lr_tile.matrixLR()));
}

TYPED_TEST(LowRankTileTest, MoveConstructor) {
    const int rows = 9;
    const int cols = 11;
    const int rank = 3;

    auto L = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rows, rank).eval();
    auto R = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rank, cols).eval();
    auto mat = (L * R).eval();

    LowRankTile<TypeParam> lr_tile(L, R);
    LowRankTile<TypeParam> lr_tile_moved_from(lr_tile);
    LowRankTile<TypeParam> lr_tile_moved_into(std::move(lr_tile_moved_from));

    // Inspect moved from value
    EXPECT_EQ(0ul, lr_tile_moved_from.size());
    // EXPECT_EQ(0ul, lr_tile_moved_from.rank()); // Can't test compiler
    // determined
    EXPECT_EQ(0ul, lr_tile_moved_from.Rows());
    EXPECT_EQ(0ul, lr_tile_moved_from.Cols());

    EXPECT_EQ(rows, lr_tile_moved_into.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(cols, lr_tile_moved_into.Cols()) << "Cols were not equal\n";
    EXPECT_EQ(rank, lr_tile_moved_into.rank()) << "Rank was not equal\n";
    EXPECT_EQ(lr_tile_moved_into.size(), lr_tile.size())
        << "Size was not equal";
    EXPECT_EQ(lr_tile_moved_into.Rows(), lr_tile.Rows())
        << "Rows were not equal\n";
    EXPECT_EQ(lr_tile_moved_into.Cols(), lr_tile.Cols())
        << "Cols were not equal\n";

    EXPECT_TRUE(lr_tile_moved_into.matrixL().isApprox(lr_tile.matrixL()));
    EXPECT_TRUE(lr_tile_moved_into.matrixR().isApprox(lr_tile.matrixR()));
    EXPECT_TRUE(lr_tile_moved_into.matrixLR().isApprox(lr_tile.matrixLR()));
}

TYPED_TEST(LowRankTileTest, AssignmentOperator) {
    const int rows = 9;
    const int cols = 11;
    const int rank = 3;

    auto L = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rows, rank).eval();
    auto R = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rank, cols).eval();
    auto mat = (L * R).eval();

    LowRankTile<TypeParam> lr_tile(L, R);
    this->tile = lr_tile;

    EXPECT_EQ(rows, this->tile.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(cols, this->tile.Cols()) << "Cols were not equal\n";
    EXPECT_EQ(rank, this->tile.rank()) << "Rank was not equal\n";
    EXPECT_EQ(this->tile.size(), lr_tile.size()) << "Size was not equal";
    EXPECT_EQ(this->tile.Rows(), lr_tile.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(this->tile.Cols(), lr_tile.Cols()) << "Cols were not equal\n";

    EXPECT_TRUE(this->tile.matrixL().isApprox(lr_tile.matrixL()));
    EXPECT_TRUE(this->tile.matrixR().isApprox(lr_tile.matrixR()));
    EXPECT_TRUE(this->tile.matrixLR().isApprox(lr_tile.matrixLR()));
}

TYPED_TEST(LowRankTileTest, MoveAssignmentOperator) {
    const int rows = 9;
    const int cols = 11;
    const int rank = 3;

    auto L = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rows, rank).eval();
    auto R = LowRankTile<TypeParam>::template Matrix
             <TypeParam>::Random(rank, cols).eval();
    auto mat = (L * R).eval();

    LowRankTile<TypeParam> lr_tile(L, R);
    LowRankTile<TypeParam> lr_tile_moved_from(lr_tile);
    this->tile = std::move(lr_tile_moved_from);

    // Inspect moved from value
    EXPECT_EQ(0ul, lr_tile_moved_from.size());
    // EXPECT_EQ(0ul, lr_tile_moved_from.rank()); // Can't test compiler
    // determined
    EXPECT_EQ(0ul, lr_tile_moved_from.Rows());
    EXPECT_EQ(0ul, lr_tile_moved_from.Cols());


    EXPECT_EQ(rows, this->tile.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(cols, this->tile.Cols()) << "Cols were not equal\n";
    EXPECT_EQ(rank, this->tile.rank()) << "Rank was not equal\n";
    EXPECT_EQ(this->tile.size(), lr_tile.size()) << "Size was not equal";
    EXPECT_EQ(this->tile.Rows(), lr_tile.Rows()) << "Rows were not equal\n";
    EXPECT_EQ(this->tile.Cols(), lr_tile.Cols()) << "Cols were not equal\n";

    EXPECT_TRUE(this->tile.matrixL().isApprox(lr_tile.matrixL()));
    EXPECT_TRUE(this->tile.matrixR().isApprox(lr_tile.matrixR()));
    EXPECT_TRUE(this->tile.matrixLR().isApprox(lr_tile.matrixLR()));
}
