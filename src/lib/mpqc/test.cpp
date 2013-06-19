#include "mpqc/comm.hpp"
#include "mpqc/math.hpp"

using namespace mpqc;

void openmp() {
  int N = 10000;
  Vector v(N), u(N);
  double e = 0;
#pragma omp parallel for private(e)
  for (int i = 0; i < v.size(); ++i) {
    for (int j = 0; j < v.size(); ++j) {
      e += v[i]*u[j];
    }
  }
  printf("%e\n", e);
}

int main() {

    {
	Comm comm(MPI_COMM_WORLD);
	openmp();
    }

    mpi::finalize();

    return 0;

}
