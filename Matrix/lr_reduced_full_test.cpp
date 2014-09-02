
#include "lr_tile.h"
#include <iostream>
#include "../include/tiledarray.h"

int main(int argc, char **argv) {
    int input = (argc > 1) ? std::stoi(argv[1]) : 50;
    int job_rank = (argc > 2) ? std::stoi(argv[2]) : input / 2;

    std::cout << "=========================================================\n"
              << "Performing test of functions for tiles with reduced rank.\n"
              << "=========================================================\n";

    // Test that multiplication works
    Eigen::MatrixXd Z = Eigen::MatrixXd::Random(input, input);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(input, input);
    {
        Z = Z.transpose() * Z;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esZ(Z);
        Eigen::MatrixXd Cz = esZ.eigenvectors();
        Eigen::VectorXd Vz = esZ.eigenvalues();


        for (auto i = 0; i < (input - job_rank); ++i) {
            Vz[i] = 0.0;
        }

        Z = Cz * Vz.asDiagonal() * Cz.transpose();
    }

    LRTile<double> B(Z);
    LRTile<double> F(Q);

    std::cout << "Test that rank of tile equals input rank.\n";
    std::cout << "\tB rank = " << job_rank << "? "
              << (B.rank() == std::size_t(job_rank)) << "\n";
    std::cout << "\tF is full ? " << (F.is_full()) << "\n";


    /*
     * Test scale
     */
    std::cout << "Testing scaling funcitons.\n";
    double LR_scale = madness::wall_time();
    LRTile<double> Bscale = B.scale(5.2);
    LRTile<double> Fscale = F.scale(5.2);
    double LR_scale_end = madness::wall_time();
    Eigen::MatrixXd Zscale = 5.2 * Z;
    Eigen::MatrixXd Qscale = 5.2 * Q;
    double Reg_scale_end = madness::wall_time();
    std::cout << "\tDoes scale work (1:yes,0:no)? "
              << (Bscale.matrixLR().isApprox(Zscale)
                  && Fscale.matrixLR().isApprox(Qscale)) << "\n"
              << "\tLR scale took " << LR_scale_end - LR_scale << " s\n"
              << "\tReg scale took " << Reg_scale_end - LR_scale_end
              << " s\n\n";


    /*
     * Test Add and Subt
     */
    std::cout << "Testing add and subt functions.\n";
    double LR_add_start = madness::wall_time();
    LRTile<double> BpF = B.add(F);
    LRTile<double> FpB = F.add(B);
    double LR_add_end = madness::wall_time();
    Eigen::MatrixXd ZpQ = Z + Q;
    Eigen::MatrixXd QpZ = Q + Z;
    double Reg_add_end = madness::wall_time();
    std::cout << "\tDoes add work (1:yes,0:no)? "
              << (BpF.matrixLR().isApprox(ZpQ) && FpB.matrixLR().isApprox(QpZ))
              << "\n"
              << "\tLR add took " << LR_add_end - LR_add_start << " s\n"
              << "\tReg add took " << Reg_add_end - LR_add_end << " s\n"
              << "Was the output low rank? " << (!BpF.is_full()) << "\n\n";


    double LR_subt_start = madness::wall_time();
    LRTile<double> BmF = B.subt(F);
    LRTile<double> FmB = F.subt(B);
    double LR_subt_end = madness::wall_time();
    Eigen::MatrixXd ZmQ = Z - Q;
    Eigen::MatrixXd QmZ = Q - Z;
    double Reg_subt_end = madness::wall_time();
    std::cout << "\tDoes subt work (1:yes,0:no)? "
              << (BmF.matrixLR().isApprox(ZmQ) && FmB.matrixLR().isApprox(QmZ))
              << "\n"
              << "\tLR subt took " << LR_subt_end - LR_subt_start << " s\n"
              << "\tReg subt took " << Reg_subt_end - LR_subt_end << " s\n"
              << "Was the output low rank? " << (!BmF.is_full()) << "\n\n";


    /*
     * Test Mult and Gemm
     */
    std::cout << "Testing mult and Gemm functions.\n";
    double LR_mult_start = madness::wall_time();
    LRTile<double> BF = B * F;
    LRTile<double> FB = F * B;
    double LR_mult_end = madness::wall_time();
    Eigen::MatrixXd ZQ = Z * Q;
    Eigen::MatrixXd QZ = Q * Z;
    double Reg_mult_end = madness::wall_time();

    std::cout << "\tDoes multiply work (1:yes,0:no)? "
              << (BF.matrixLR().isApprox(ZQ) && FB.matrixLR().isApprox(QZ))
              << "\n"
              << "\tLR mult took " << LR_mult_end - LR_mult_start << " s\n"
              << "\tReg mult took " << Reg_mult_end - LR_mult_end << " s\n"
              << "\toutput is low rank? " << (!BF.is_full() && !FB.is_full())
              << " s\n\n";


    double LR_gemm_start = madness::wall_time();
    LRTile<double> temp = B.add(F);
    LRTile<double> Gemm = (B * F).add(temp);
    double LR_gemm_end = madness::wall_time();
    Eigen::MatrixXd C = Z + Q;
    C = Z * Q + C;
    double Reg_gemm_end = madness::wall_time();
    std::cout << "\tDoes gemm work (1:yes,0:no)? "
              << Gemm.matrixLR().isApprox(C) << "\n"
              << "\tLR gemm took " << LR_gemm_end - LR_gemm_start << " s\n"
              << "\tReg gemm took " << Reg_gemm_end - LR_gemm_end << " s\n"
              << "\tout is low rank? " << (!Gemm.is_full()) << "\n\n";



    /*
     * Test the XXX_to methods.
     */
    std::cout << "Testing xxx_to methods.\n";
    double LR_subt_to_start = madness::wall_time();
    B.subt_to(F);
    double LR_subt_to_end = madness::wall_time();
    Z -= Q;
    double Reg_subt_to_end = madness::wall_time();
    std::cout << "\tDoes subt_to work (1:yes,0:no)? "
              << B.matrixLR().isApprox(Z) << "\n"
              << "\tLR subt_to took " << LR_subt_to_end - LR_subt_to_start
              << " s\n"
              << "\tReg subt_to took " << Reg_subt_to_end - LR_subt_to_end
              << " s\n\n";

    double LR_subt_to_factor_start = madness::wall_time();
    B.subt_to(F, 4.7);
    double LR_subt_to_factor_end = madness::wall_time();
    (Z -= Q) *= 4.7;
    double Reg_subt_to_factor_end = madness::wall_time();
    std::cout << "\tDoes subt_to with factor work (1:yes,0:no)? "
              << B.matrixLR().isApprox(Z) << "\n"
              << "\tLR subt_to factor took "
              << LR_subt_to_factor_end - LR_subt_to_factor_start << " s\n"
              << "\tReg subt_to factor took "
              << Reg_subt_to_factor_end - LR_subt_to_factor_end << " s\n\n";

    double LR_scale_to = madness::wall_time();
    B.scale_to(5.2);
    double LR_scale_to_end = madness::wall_time();
    Z *= 5.2;
    double Reg_scale_to_end = madness::wall_time();
    std::cout << "\tDoes scale_to work (1:yes,0:no)? "
              << B.matrixLR().isApprox(Z) << "\n"
              << "\tLR scale took " << LR_scale_to_end - LR_scale_to << " s\n"
              << "\tReg scale took " << Reg_scale_to_end - LR_scale_to_end
              << " s\n\n";

#if 0
#endif

    return 0;
}
