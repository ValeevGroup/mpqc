#ifndef CLUSTERABLE_SHELL_H
#define CLUSTERABLE_SHELL_H

//#include <chemistry/qc/basis/shell.h>
#include <memory>
#include "../molecule/molecule_fwd.h"

/**
 * @brief The ClusterableShell class is used to reorder shells in given angular
 * momentum blocks.
 */
class ClusterableShell {
private :
  std::shared_ptr<const Atom> shell_owner;
};

#endif // CLUSTERABLE_SHELL_H
