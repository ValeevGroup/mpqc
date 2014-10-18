#include "../../include/eigen.h"
#include <iostream>
#include <string>
#include <chrono>

using namespace Eigen;

MatrixXd formQ(ColPivHouseholderQR<MatrixXd> const &qr) {
    MatrixXd Q = MatrixXd::Identity(qr.matrixQR().rows(), qr.rank());
    qr.householderQ().applyThisOnTheLeft(Q);
    return Q;
}

MatrixXd formR(ColPivHouseholderQR<MatrixXd> const &qr) {
    MatrixXd R = MatrixXd(qr.matrixR()
                              .topLeftCorner(qr.rank(), qr.matrixQR().cols())
                              .template triangularView<Upper>())
                 * qr.colsPermutation().transpose();
    return R;
}

MatrixXd makeLR_mat(int rows, int rank) {
    MatrixXd mat = MatrixXd::Random(rows, rows);
    mat = mat.transpose() * mat;
    SelfAdjointEigenSolver<MatrixXd> es(mat);
    MatrixXd C = es.eigenvectors();
    VectorXd V = es.eigenvalues();

    // Create Low Rank Matrix
    std::for_each(V.data(), V.data() + (rows - rank),
                  [](double &x) { x = 0.0; });
    mat = C * V.asDiagonal() * C.transpose();

    return mat;
}

int main(int argc, char **argv) {
    int input = (argc > 1) ? std::stoi(argv[1]) : 500;
    int rank = (argc > 2) ? std::stoi(argv[2]) : input / 2;

    MatrixXd mat = makeLR_mat(input, rank);
    MatrixXd mat2 = makeLR_mat(input, rank);

    ColPivHouseholderQR<MatrixXd> qr(mat);
    MatrixXd Q = formQ(qr);
    MatrixXd R = formR(qr);
    auto rank1 = qr.rank();
    std::cout << "Rank of mat1 = " << rank1 << std::endl;

    ColPivHouseholderQR<MatrixXd> qr2(mat2);
    MatrixXd Q2 = formQ(qr2);
    MatrixXd R2 = formR(qr2);
    auto rank2 = qr2.rank();
    std::cout << "Rank of mat2 = " << rank2 << std::endl;

    // Create 2 * rank mats for add
    MatrixXd La(mat.rows(), rank1 + rank2);
    MatrixXd Ra(rank1 + rank2, mat.cols());

    // Fill matrices with data
    La.leftCols(rank1) = Q;
    La.rightCols(rank2) = Q2;
    Ra.topRows(rank1) = R;
    Ra.bottomRows(rank2) = R2;

    // Assuming C = A + B, C currently has rank(C) = rank(A) + rank(B), we
    // should be able to reduce this at least a little.

    // Time Eigen add for reference
    auto start = std::chrono::steady_clock::now();
    MatrixXd Correct = mat + mat2;
    auto end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast
        <std::chrono::duration<double>>(end - start);

    std::cout << "Eigen add took " << time.count() << " s\n\n";

    // Original Method
    {
        // Lets try using ColPivQr to reduce the rank of L^C and R^C
        auto start = std::chrono::steady_clock::now();
        qr.compute(La);
        ColPivHouseholderQR<MatrixXd> qrR(Ra);

        // At this point just do the same thing we did for multiplication
        MatrixXd LaQ = formQ(qr);
        MatrixXd RaQ = formQ(qrR);
        MatrixXd LaR = formR(qr);
        MatrixXd RaR = formR(qrR);

        // New low rank L and low rank R.  Arbitrarrialy picked L to absorb
        // the small matrix.
        MatrixXd FinalL = LaQ * (LaR * RaQ);
        MatrixXd FinalR = RaR;
        auto end = std::chrono::steady_clock::now();
        auto time = std::chrono::duration_cast
            <std::chrono::duration<double>>(end - start);

        // Calculate Full added matrix and check if correct!
        MatrixXd approx = MatrixXd(FinalL * FinalR);
        std::cout << "Original add method took " << time.count() << " s\n";
        std::cout << "Is original add correct ? " << (Correct).isApprox(approx)
                  << "\n\n";
    }

    // Trying something new.
    {
        auto start = std::chrono::steady_clock::now();
        qr.compute(La);
        MatrixXd LaQ = formQ(qr);
        MatrixXd LaR = formR(qr);
        MatrixXd FinalL = LaQ;
        MatrixXd FinalR = LaR * Ra;
        auto end = std::chrono::steady_clock::now();
        auto time = std::chrono::duration_cast
            <std::chrono::duration<double>>(end - start);

        // Calculate Full added matrix and check if correct!
        MatrixXd approx = MatrixXd(FinalL * FinalR);
        std::cout << "New add method took " << time.count() << " s\n";
        std::cout << "Is new add correct ? " << (Correct).isApprox(approx)
                  << "\n\n";
    }

    return 0;
}
