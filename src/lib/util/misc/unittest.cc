//
// unittest.cc
//

#include <util/misc/units.h>
#include <util/misc/formio.h>

#include <util/state/state.h>
#include <util/state/stateout.h>
static void (sc::SavableState::*force_state_link)(sc::StateOut&)
    = &sc::SavableState::save_state;
#include <util/state/linkage.h>

using namespace std;
using namespace sc;

//////////////////////////////////////////////////////////////////////
// Unit conversion test program 

int main(int argc, char **argv) {
  const char *unitstr = "kcal/mol";

  if (argc == 2) {
      unitstr = argv[1];
    }
  else if (argc != 1) {
      cerr << "One argument, the unit to be converted, must be given" << endl;
      abort();
    }
  
  Ref<Units> unit = new Units(unitstr);
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
