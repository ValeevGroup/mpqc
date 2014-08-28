#include "lr_tile.h"
#include <iostream>
#include <vector>
#include "../include/tiledarray.h"

int main(int argc, char **argv) {
    madness::World &world = madness::initialize(argc, argv);
    world.rank(); // Silence warning

    int input = (argc > 1) ? std::stoi(argv[1]) : 50;
    int job_rank = (argc > 2) ? std::stoi(argv[2]) : input / 2;

    // Test that multiplication works
    Eigen::MatrixXd Z = Eigen::MatrixXd::Random(input, input);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(input, input);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(input, input);
    {
        Z = Z.transpose() * Z;
        Q = Q.transpose() * Q;
        C = C.transpose() * C;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esZ(Z);
        Eigen::MatrixXd Cz = esZ.eigenvectors();
        Eigen::VectorXd Vz = esZ.eigenvalues();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esQ(Q);
        Eigen::MatrixXd Cq = esQ.eigenvectors();
        Eigen::VectorXd Vq = esQ.eigenvalues();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esC(Q);
        Eigen::MatrixXd Cc = esC.eigenvectors();
        Eigen::VectorXd Vc = esC.eigenvalues();

        for (auto i = 0; i < (input - job_rank); ++i) {
            Vz[i] = 0.0;
            Vq[i] = 0.0;
            Vc[i] = 0.0;
        }

        Z = Cz * Vz.asDiagonal() * Cz.transpose();
        Q = Cq * Vq.asDiagonal() * Cq.transpose();
        C = Cc * Vc.asDiagonal() * Cc.transpose();
    }

    LRTile<double> B(Z);
    LRTile<double> F(Q);
    LRTile<double> X(Q);

    /*
     * Test Add and Subt
     */

    double LR_add_start = madness::wall_time();
    LRTile<double> BpF = B.add(F);
    double LR_add_end = madness::wall_time();
    Eigen::MatrixXd ZpQ = Z + Q;
    double Reg_add_end = madness::wall_time();
    std::cout << "Does add work (1:yes,0:no)? " << BpF.matrixLR().isApprox(ZpQ)
              << "\n"
              << "LR add took " << LR_add_end - LR_add_start << " s\n"
              << "Reg add took " << Reg_add_end - LR_add_end << " s\n\n";


    double LR_subt_start = madness::wall_time();
    LRTile<double> BmF = B.subt(F);
    double LR_subt_end = madness::wall_time();
    Eigen::MatrixXd ZmQ = Z - Q;
    double Reg_subt_end = madness::wall_time();
    std::cout << "Does subt work (1:yes,0:no)? " << BmF.matrixLR().isApprox(ZmQ)
              << "\n"
              << "LR subt took " << LR_subt_end - LR_subt_start << " s\n"
              << "Reg subt took " << Reg_subt_end - LR_subt_end << " s\n\n";

    double LR_subt_to_start = madness::wall_time();
    B.subt_to(F);
    double LR_subt_to_end = madness::wall_time();
    Z -= Q;
    double Reg_subt_to_end = madness::wall_time();
    std::cout << "Does subt_to work (1:yes,0:no)? " << B.matrixLR().isApprox(Z)
              << "\n"
              << "LR subt_to took " << LR_subt_to_end - LR_subt_to_start
              << " s\n"
              << "Reg subt_to took " << Reg_subt_to_end - LR_subt_to_end
              << " s\n\n";

    /*
     * Test Mult and Gemm
     */
    double LR_mult_start = madness::wall_time();
    LRTile<double> BF = B * F;
    double LR_mult_end = madness::wall_time();
    Eigen::MatrixXd ZQ = Z * Q;
    double Reg_mult_end = madness::wall_time();

    std::cout << "Does multiply work (1:yes,0:no)? "
              << BF.matrixLR().isApprox(Z * Q) << "\n"
              << "LR mult took " << LR_mult_end - LR_mult_start << " s\n"
              << "Reg mult took " << Reg_mult_end - LR_mult_end << " s\n\n";

    double LR_gemm_start = madness::wall_time();
    LRTile<double> temp = B * F;
    LRTile<double> Gemm = temp.add(X);
    double LR_gemm_end = madness::wall_time();
    C = Z * Q + C;
    double Reg_gemm_end = madness::wall_time();
    std::cout << "Does gemm work (1:yes,0:no)? " << Gemm.matrixLR().isApprox(C)
              << "\n"
              << "LR gemm took " << LR_gemm_end - LR_gemm_start << " s\n"
              << "Reg gemm took " << Reg_gemm_end - LR_gemm_end << " s\n\n";


    madness::finalize();
    return 0;
}
