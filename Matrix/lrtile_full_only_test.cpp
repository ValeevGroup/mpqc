#include "lr_tile.h"
#include "../include/tiledarray.h"

#include <iostream>

int main(int argc, char *argv[]) {
    int input = (argc > 1) ? std::stoi(argv[1]) : 50;

    std::cout << "======================================================\n"
              << "Performing test of functions for tiles with full rank.\n"
              << "======================================================\n";

    // Test that multiplication works
    Eigen::MatrixXd Z = Eigen::MatrixXd::Random(input, input);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(input, input);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(input, input);

    LRTile<double> B(Z, false);
    LRTile<double> F(Q, false);
    LRTile<double> X(C, false);

    std::cout << "Testing that all tiles are full rank.\n";
    std::cout << "\tB is full? " << B.is_full() << "\n";
    std::cout << "\tF is full? " << F.is_full() << "\n";
    std::cout << "\tX is full? " << X.is_full() << "\n\n";


    /*
     * Test scale
     */
    std::cout << "Testing scaling funcitons.\n";
    double LR_scale = madness::wall_time();
    LRTile<double> Bscale = B.scale(5.2);
    double LR_scale_end = madness::wall_time();
    Eigen::MatrixXd Zscale = 5.2 * Z;
    double Reg_scale_end = madness::wall_time();
    std::cout << "\tDoes scale work (1:yes,0:no)? "
              << (Bscale.is_full() && Bscale.matrixLR().isApprox(Zscale))
              << "\n"
              << "\tLR scale took " << LR_scale_end - LR_scale << " s\n"
              << "\tReg scale took " << Reg_scale_end - LR_scale_end
              << " s\n\n";

    /*
     * Test Add and Subt
     */
    std::cout << "Testing add and subt functions.\n";
    double LR_add_start = madness::wall_time();
    LRTile<double> BpF = B.add(F);
    double LR_add_end = madness::wall_time();
    Eigen::MatrixXd ZpQ = Z + Q;
    double Reg_add_end = madness::wall_time();
    std::cout << "\tDoes add work (1:yes,0:no)? "
              << BpF.matrixLR().isApprox(ZpQ) << "\n"
              << "\tLR add took " << LR_add_end - LR_add_start << " s\n"
              << "\tReg add took " << Reg_add_end - LR_add_end << " s\n\n";


    double LR_subt_start = madness::wall_time();
    LRTile<double> BmF = B.subt(F);
    double LR_subt_end = madness::wall_time();
    Eigen::MatrixXd ZmQ = Z - Q;
    double Reg_subt_end = madness::wall_time();
    std::cout << "\tDoes subt work (1:yes,0:no)? "
              << BmF.matrixLR().isApprox(ZmQ) << "\n"
              << "\tLR subt took " << LR_subt_end - LR_subt_start << " s\n"
              << "\tReg subt took " << Reg_subt_end - LR_subt_end << " s\n\n";


    /*
     * Test Mult and Gemm
     */
    std::cout << "Testing mult and Gemm functions.\n";
    double LR_mult_start = madness::wall_time();
    LRTile<double> BF = B * F;
    double LR_mult_end = madness::wall_time();
    Eigen::MatrixXd ZQ = Z * Q;
    double Reg_mult_end = madness::wall_time();

    std::cout << "\tDoes multiply work (1:yes,0:no)? "
              << BF.matrixLR().isApprox(Z * Q) << "\n"
              << "\tLR mult took " << LR_mult_end - LR_mult_start << " s\n"
              << "\tReg mult took " << Reg_mult_end - LR_mult_end << " s\n\n";

    double LR_gemm_start = madness::wall_time();
    LRTile<double> temp = B * F;
    LRTile<double> Gemm = temp.add(X);
    double LR_gemm_end = madness::wall_time();
    C = Z * Q + C;
    double Reg_gemm_end = madness::wall_time();
    std::cout << "\tDoes gemm work (1:yes,0:no)? "
              << Gemm.matrixLR().isApprox(C) << "\n"
              << "\tLR gemm took " << LR_gemm_end - LR_gemm_start << " s\n"
              << "\tReg gemm took " << Reg_gemm_end - LR_gemm_end << " s\n\n";

    /*
     * Test the XXX_to methods.
     */
    std::cout << "Testing xxx_to methods.\n";

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




#if 0
#endif
    return 0;
}
