//
// unittest.cc
//

#include <util/misc/units.h>
#include <util/misc/formio.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Unit conversion test program 

int main(int argc, char **argv) {
  if (argc != 2) {
      cerr << "One argument, the unit to be converted, must be given" << endl;
      abort();
    }
  
  Ref<Units> unit = new Units(argv[1]);
  cout << indent << "Conversion between " << unit->string_rep() << " and atomic units:" << endl;
  cout << setprecision(10);
  cout << indent << "From atomic units: " << unit->from_atomic_units() << endl;
  cout << indent << "To atomic units: " << unit->to_atomic_units() << endl;

  return 0;

}


// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
