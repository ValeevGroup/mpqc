#include <iostream>
#include <sstream>
#include <iomanip>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/molecule/atominfo.h>
#include <Chemistry_QC_GaussianBasis_MolecularInterface.hxx>
#include <Chemistry_MoleculeInterface.hxx>
#include "basis_cca_to_sc.h"

using namespace std;
using namespace sc;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;

namespace MPQC {
  
Ref<GaussianBasisSet> basis_cca_to_sc( MolecularInterface &cca_basis ) {

  const char* am_to_symbol[] = {"s", "p", "d", "f", "g", "h", "i", "k", "l",
			       "m", "n", "o", "p", "q", "r", "s", "t", "u",
			       "v", "w", "x", "y", "z"};

  //cca_basis.print_molecular();

  Chemistry::MoleculeInterface cca_mol = cca_basis.get_molecule();

  ostringstream input;

  // form molecule keyval
  double conv = cca_mol.get_units().convert_to("bohr");
  input
    << "molecule<Molecule>: (\n"
    << "  unit = bohr\n"
    << "  symmetry = C1\n"
    << "  {n atoms geometry } = {\n";
  for( int i=0; i<cca_mol.get_n_atom(); ++i ) {
    input << setprecision(16);
    input << "\t" << i << "\t" << cca_mol.get_atomic_number(i)
      << "\t[  " << cca_mol.get_cart_coor(i,0)*conv
      << "  " << cca_mol.get_cart_coor(i,1)*conv
      << "  " << cca_mol.get_cart_coor(i,2)*conv << "  ]\n";
  }
  input << "  }\n" << ")\n";

  // form basis keyval
  input.precision(18);
  input << "scbasis<GaussianBasisSet>:(\n" 
	<< "  molecule = $:molecule\n";
  if (cca_basis.get_label().size() > 0) {
      input << "  name = \"CCA(" << cca_basis.get_label() << ")\"\n";
    }
  input << "  basis = [";
  for(int i=0; i<cca_mol.get_n_atom(); ++i) 
    input << " basis" << i;
  input << " ]\n" << ")\n";

  input << "basis:(\n";

  // form atomic set for each individual center (possibly redundant)
  // <atomname>: <basisname>: [
  AtomInfo empty_info;
  for(int i=0; i<cca_mol.get_n_atom(); ++i) {
    AtomicInterface atomic = cca_basis.get_atomic(i);
    input << " " << empty_info.name( cca_mol.get_atomic_number(i) ) << ": "
	  << " basis" << i << ": [\n";
    
    // form shells
    for(int ishell=0; ishell<atomic.get_n_shell(); ++ishell) {
      ShellInterface shell = atomic.get_shell(ishell);
      // (type: [ am = <symbol> ...]
      input << "  (normalized = 0\n";
      //input << "  (normalized = 1\n";
      input << "   type: [";
      for(int icon=0; icon<shell.get_n_contraction(); ++icon) {
	input << " (am = " << am_to_symbol[shell.get_angular_momentum(icon)];
        if( shell.get_max_angular_momentum() > 1 ) {
          if( shell.get_angular_type() == AngularType_CARTESIAN )
            input << " puream = 0)";
          else if( shell.get_angular_type() == AngularType_SPHERICAL )
            input << " puream = 1)";
          else if( shell.get_angular_type() == AngularType_MIXED )
            std::cerr << " mixed angular types?";
        }
        else input << ")";
      }
      input << "]\n";
      // {exp coef:<am> ...} = {
      input << "   {exp";
      for(int icon=0; icon<shell.get_n_contraction(); ++icon)
	input << " coef:" << icon;
      input << "} = {\n";
      // <exp> <coef> ...
      for(int iprim=0; iprim<shell.get_n_primitive(); ++iprim) {
	input << "\t" << shell.get_exponent(iprim);
	for(int icon=0; icon<shell.get_n_contraction(); ++icon)
	  input << "\t" << shell.get_contraction_coef(icon, iprim);
	input << endl;
      }
      input << "\n   })\n";
    }
    input << " ]\n";
  }
  input << ")\n";

  //cout << "  basis input:\n" << input.str() << endl;

  Ref<ParsedKeyVal> kv = new ParsedKeyVal();
  kv->parse_string( input.str().c_str() );
  Ref<DescribedClass> dc = kv->describedclassvalue("scbasis");
  GaussianBasisSet *sc_basis = 
    dynamic_cast< GaussianBasisSet* >( dc.pointer() );

  Ref<GaussianBasisSet> gbs;
  gbs.assign_pointer(sc_basis);

  //std::cerr << "basis converter's molcule:\n";
  //gbs->molecule()->print();
  
  //for(int i=0; i<gbs->nshell(); ++i)
  //  gbs->shell(i).print();

  return gbs;
}

} // end of namespace
