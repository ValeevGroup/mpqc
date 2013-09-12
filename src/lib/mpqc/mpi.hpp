#ifndef MPQC_MPI_HPP
#define MPQC_MPI_HPP

/// @defgroup MPI mpqc.Core.MPI
/// <a href=http://www.mpi-forum.org>MPI</a> wrappers/stubs

#include "mpqc_config.h"

#ifndef HAVE_MPI
#warning MPQC will use serial MPI stub
#endif // HAVE_MPI

#include "mpqc/mpi/base.hpp"
#include "mpqc/mpi/comm.hpp"

// N.B. may require ARMCI
//#include "mpqc/mpi/task.hpp"

#endif // MPQC_MPI_HPP
