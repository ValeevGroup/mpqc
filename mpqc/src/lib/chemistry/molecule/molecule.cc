
#include "molecule.h"

DescribedClass_REF_def(Molecule);

#define CLASSNAME Molecule
#define PARENTS virtual public SavableState
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
Molecule::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

Molecule::Molecule() :
  atoms(0), natoms(0)
{
}

Molecule::Molecule(Molecule& mol) :
  atoms(0), natoms(0)
{
  *this=mol;
}

Molecule::~Molecule()
{
  if(atoms) delete[] atoms; atoms=0;
  natoms=0;
}

Molecule::Molecule(KeyVal&input):
//atoms(0,10)
atoms(0),
natoms(0)
{

  // get the number of atoms and make sure that the geometry and the
  // atoms array have the same number of atoms.
  int natom = input.count("geometry");
  if (natom != input.count("atoms")) {
      fprintf(stderr,"Molecule: size of \"geometry\" != size of \"atoms\"\n");
      return;
    }

  int aangstroms = input.booleanvalue("angstrom");
  if (input.error() != KeyVal::OK) {
    aangstroms = input.booleanvalue("aangstrom");
    }
  double conv = 1.0;
  if (aangstroms) {
      conv = 1.0/0.52917706;
    }
  int i;
  for (i=0; i<natom; i++) {
      char* name;
      char* labels;
      AtomicCenter ac(name = input.pcharvalue("atoms",i),
		      input.doublevalue("geometry",i,0)*conv,
		      input.doublevalue("geometry",i,1)*conv,
		      input.doublevalue("geometry",i,2)*conv,
                      labels = input.pcharvalue("atom_labels",i)
		      );
      delete[] name;
      delete[] labels;
      add_atom(i,ac);
    }
}

Molecule& Molecule::operator=(Molecule& mol)
{
  if(atoms) delete[] atoms; atoms=0;
  natoms=0;

  for(int i=0; i < mol.natom(); i++) {
       AtomicCenter ac = mol[i];
       add_atom(i,ac);
    }

  return *this;
}

int
Molecule::natom() const
{
  return natoms;
}

const AtomicCenter&
Molecule::get_atom(int i) const
{
  return atoms[i];
}

AtomicCenter&
Molecule::get_atom(int i)
{
  return atoms[i];
}

Pix Molecule::first()
{
  if (natoms) return (Pix) 1;
  else return (Pix) 0;
}

void Molecule::next(Pix& i)
{
  if ((int)i < natoms) ((int)i)++;
  else i = 0;
}

int Molecule::owns(Pix i)
{
  if ((int)i > 0 && (int)i <= natoms) return 1;
  else return 0;
}

void
Molecule::add_atom(int i,AtomicCenter& ac)
{
  if (i>=natoms) {
      AtomicCenter* new_atoms = new AtomicCenter[i+1];
      for (int j=0; j<natoms; j++) {
	  new_atoms[j] = atoms[j];
	}
      if (atoms) delete[] atoms;
      atoms = new_atoms;
      natoms = i+1;
    }
  atoms[i] = ac;
}

void
Molecule::print(FILE*fp)
{
  fprintf(fp,"Molecule:\n");
  fprintf(fp,"    n  atom  label          x               y               z"
             "          mass\n");

  int i;
  for (i=0; i<natom(); i++) {
      fprintf(fp," %5d%5s%8s%16.10f%16.10f%16.10f%10.5f\n",
              i+1,
	      get_atom(i).element().symbol(),
	      (get_atom(i).label()) ? get_atom(i).label(): " ",
	      get_atom(i)[0],
	      get_atom(i)[1],
	      get_atom(i)[2],
              get_atom(i).element().mass()
              );
    }
}

AtomicCenter& Molecule::operator[](int ind)
{
  return get_atom(ind);
}

AtomicCenter& Molecule::operator()(Pix pix)
{
  return get_atom(((int)pix)-1);
}

const AtomicCenter& Molecule::operator[](int ind) const
{
  return get_atom(ind);
}

const AtomicCenter& Molecule::operator()(Pix pix) const
{
  return get_atom(((int)pix)-1);
}

PointBag_double* Molecule::charges() const
{
  PointBag_double*result = new PointBag_double;
  int i;
  for (i=0; i<natom(); i++) {
      result->add(get_atom(i).point(),get_atom(i).element().charge());
    }
  return result;
}

void Molecule::save_data_state(StateOut& so)
{
  so.put(natoms);
  for (int i=0; i < natoms; i++) {
      get_atom(i).save_object_state(so);
    }
}

Molecule::Molecule(StateIn& si):
  atoms(0),
  natoms(0),
  SavableState(si,class_desc_)
{
  int natom;
  si.get(natom);
  natoms=0;
  for (int i=0; i < natom; i++) {
      AtomicCenter ac(si);
      add_atom(i,ac);
    }
}


//QCMolecule::QCMolecule(KeyVal&input):
//Molecule(input)
//{
//  // get the name of the basis set
//  basissetname = input.pcharvalue("basis");
//
//  // generate the GaussianBasisSet
//  //basisset = 0;
//}
//
//QCMolecule::~QCMolecule()
//{
//  //if (basisset) delete basisset;
//  if (basissetname) delete basissetname;
//}
//
//void QCMolecule::print(FILE*fp)
//{
//  Molecule::print(fp);
//  if (basissetname) fprintf(fp,"basis = %s\n",basissetname);
//}
