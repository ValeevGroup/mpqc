#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>

int main(int argc, char* argv[]){
  Eigen::MatrixXd m = Eigen::MatrixXd::Random(5, 5);
  m = m.transpose() + m;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(m);
  Eigen::MatrixXd m_invsqrt = eig.operatorInverseSqrt();
  std::cout << m_invsqrt << std::endl;
  return 0;
}
