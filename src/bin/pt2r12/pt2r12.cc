
#include <iostream>
#include <exception>

using std::cout;
using std::endl;

extern int main_gamess(int argc, char *argv[]);

int
main(int argc, char *argv[])
{
  try {
    main_gamess(argc, argv);
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
