#include "../molecule/atom.h"
BOOST_TEST_MODULE("atom_test");
#include <boost/test/included/unit_test.hpp>

using position_t = Atom::position_t;

BOOST_AUTO_TEST_CASE(atom_constructors){
  // Default
  BOOST_REQUIRE_NO_THROW(Atom a);
  Atom a;
  BOOST_CHECK_CLOSE(a.charge(),0, 1e-16);
  BOOST_CHECK_CLOSE(a.mass(),0, 1e-16);
  //BOOST_CHECK_CLOSE(a.center(),position_t({0,0,0}));

  //Copy
  BOOST_REQUIRE_NO_THROW(Atom a_copy(a));
  Atom a_copy(a);
  BOOST_CHECK_CLOSE(a_copy.charge(),0, 1e-16);
  BOOST_CHECK_CLOSE(a_copy.mass(),0, 1e-16);
  //BOOST_CHECK_CLOSE(a_copy.center(),position_t({0,0,0}));

  // Normal ctor
  BOOST_REQUIRE_NO_THROW(Atom n({0,0,1},1.0,1.0));
  Atom n({0,0,1},1.0,1.0);
  BOOST_CHECK_CLOSE(n.charge(),1.0, 1e-16);
  BOOST_CHECK_CLOSE(n.mass(),1.0, 1e-16);
  //BOOST_CHECK_CLOSE(n.center(),position_t({0,0,0}));

}

