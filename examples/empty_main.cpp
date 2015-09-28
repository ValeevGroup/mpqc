#include <memory>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include "../include/libint.h"
#include "../include/tiledarray.h"
#include "../include/btas.h"

#include "../utility/make_array.h"
#include "../utility/parallel_print.h"
#include "../utility/parallel_break_point.h"
#include "../utility/array_storage.h"
#include "../utility/time.h"
#include "../utility/json_input.h"

#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "../integrals/btas_to_ta_tensor.h"
#include "../integrals/make_engine.h"
#include "../integrals/integral_engine_pool.h"
#include "../integrals/sparse_task_integrals.h"

#include "../scf/soad.h"
#include "../scf/diagonalize_for_coffs.hpp"
#include "../scf/clusterd_coeffs.h"

#include "../ta_routines/array_to_eigen.h"

int main(int argc, char** argv){
    auto &world = madness::initialize(argc, argv);
    return 0;
}
