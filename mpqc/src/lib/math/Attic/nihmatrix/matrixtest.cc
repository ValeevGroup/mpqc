
#include <ostream.h>
#include "nihmatrix.h"
#include <util/keyval/keyval.h>

main()
{
  ParsedKeyVal kv("matrixtest.in");

  // read in the matrix
  RefDMatrix matrix(kv.describedclassvalue("matrix"));
  cout << "matrix:\n";
  if (matrix.nonnull()) matrix->print();

  // another way to read in the matrix
  RefDMatrix matrix2(new DMatrix(PrefixKeyVal("matrix",kv)));
  cout << "matrix2:\n";
  if (matrix2.nonnull()) matrix2->print();

  // reading references
  RefDMatrix refmatrix(kv.describedclassvalue("reference_to_matrix"));

  // show that the reference to matrix and matrix are identical
  cout << "matrix pointers: "
       << matrix.pointer()
       << " should be the same as "
       << refmatrix.pointer()
       << "\n";

  // try reading vectors too
  RefDVector v(kv.describedclassvalue("vector"));
  cout << "vector:\n";
  if (v.nonnull()) v->print();

}
