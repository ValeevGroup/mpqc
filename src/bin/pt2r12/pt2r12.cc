
#include <iostream>
#include <exception>

// Force linkages:
#include <scdirlist.h>
#include <util/group/linkage.h>
#include <chemistry/qc/basis/linkage.h>
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_PSI
#include <chemistry/qc/psi/linkage.h>
#endif
#include <chemistry/qc/wfn/linkage.h>
#include <chemistry/qc/scf/linkage.h>
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_MBPTR12
#  include <chemistry/qc/mbptr12/linkage.h>
#endif
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_LIBINT2
#  include <chemistry/qc/libint2/linkage.h>
#endif
#include <util/state/linkage.h>

#ifdef HAVE_MPI
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <util/group/messmpi.h>
#endif

using std::cout;
using std::endl;

extern int main_gamess(int argc, char *argv[]);
extern int main_molcas(int argc, char *argv[]);

int
main(int argc, char *argv[])
{
  try {
    //main_gamess(argc, argv);
    main_molcas(argc, argv);
  }
  catch (std::bad_alloc &e) {
      cout << argv[0] << ": ERROR: MEMORY ALLOCATION FAILED:" << endl
           << e.what()
           << endl;
      exit(1);
    }
  catch (std::exception &e) {
      cout << argv[0] << ": ERROR: EXCEPTION RAISED:" << endl
           << e.what()
           << endl;
      exit(1);
    }
  catch (...) {
      cout << argv[0] << ": ERROR: UNKNOWN EXCEPTION RAISED" << endl;
      exit(1);
    }
  return 0;
}
