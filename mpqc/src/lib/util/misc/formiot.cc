
#include <util/misc/formio.h>

int
main()
{
  int elem = ios::xalloc();
  cout << " elem = " << elem << endl;

  cout << indent << "l0" << endl;
  cout << incindent;
  cout << indent << "l1" << endl;
  cout << incindent;
  cout << indent << "l2" << endl;
  cout << indent << "l2" << endl;
  long ind = SCFormIO::getindent(cout);
  cout << indent << "xyz = " << skipnextindent;
  SCFormIO::setindent(cout,SCFormIO::getindent(cout) + 6);
  cout << indent << "lxyz0" << endl;
  cout << indent << "lxyz1" << endl;
  cout << indent << "lxyz2" << endl;
  SCFormIO::setindent(cout,ind);
  cout << decindent;
  cout << indent << "l1" << endl;
  cout << decindent;
  cout << indent << "l0" << endl;

  cout << node0 << indent << scprintf("%3d %10.5f",10,3.14) << endl;
}

