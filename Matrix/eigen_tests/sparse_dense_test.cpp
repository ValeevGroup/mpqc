#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <random>
#include <chrono>
#include <string>

using namespace Eigen;

int main(int argc, char** argv){
  int input = (argc > 1) ? std::stoi(argv[1]) : 1000;
  int num_insig = (argc > 2) ? std::stoi(argv[2]) : 990;

  MatrixXd mat = MatrixXd::Random(input, input);
  std::mt19937 gen(100);
  std::uniform_int_distribution<int> dist(0,input-1);

  for(auto j = 0; j < mat.cols(); ++j){
    for(auto i = 0; i < num_insig; ++i){
        mat(j, dist(gen)) = 0.0; 
    }
  }

  if(mat.size() < 100){
    std::cout << "mat = \n" << mat << std::endl;
  }

  long num_zeros = 0;
  for(auto i = 0; i < mat.cols(); ++i){
    for(auto j = 0; j < mat.rows(); ++j){
      num_zeros += (mat(i,j) == 0) ? 1 : 0;
    }
  }

  double percent_zero =  double(num_zeros)/double(mat.size()) * 100;
    

  SparseMatrix<double> smat = mat.sparseView();
  smat.makeCompressed();

  std::chrono::steady_clock::time_point sparse_mult0 = std::chrono::steady_clock::now(); 
  smat = smat * smat;
  std::chrono::steady_clock::time_point sparse_mult1 = std::chrono::steady_clock::now();
  double sparse_time = std::chrono::duration_cast<std::chrono::duration<double>>(sparse_mult1-sparse_mult0).count();

  std::chrono::steady_clock::time_point dense_mult0 = std::chrono::steady_clock::now(); 
  mat = mat * mat;
  std::chrono::steady_clock::time_point dense_mult1 = std::chrono::steady_clock::now();
  double dense_time = std::chrono::duration_cast<std::chrono::duration<double>>(dense_mult1-dense_mult0).count();

  std::cout << input << ", " << percent_zero << ", " << dense_time << ", " << sparse_time << "\n";

  return 0;
}
