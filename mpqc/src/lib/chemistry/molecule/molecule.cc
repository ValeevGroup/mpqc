
#include <math.h>
#include "molecule.h"

SavableState_REF_def(Molecule);

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
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
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
  const double ang_to_bohr = 1.0/0.52917706;

  if (input.exists("pdb_file")) {
      const int LineLength = 85;
      char line[LineLength];
      char* filename = input.pcharvalue("pdb_file");
      FILE*fp = fopen(filename,"r");
      if (!fp) {
          fprintf(stderr,
                  "Molecule::Molecule(KeyVal&input): pdb file not found\n");
          abort();
        }
      int i=0;
      while(fgets(line,LineLength,fp)) {
          if (!strncmp(line,"HETA",4) || !strncmp(line,"ATOM",4)) {
              char atomsym[3];
              // find the atomic symbol
              if (line[12] == ' ') {
                  atomsym[0] = line[13];
                  atomsym[1] = '\0';
                }
              else {
                  atomsym[0] = line[12];
                  atomsym[1] = line[13];
                  atomsym[2] = '\0';
                }
              // skip dummy atoms
              if (!strcmp(atomsym,"Q")) continue;
              char position[9];
              position[8] = '\0';
              // x
              strncpy(position,&line[30],8);
              double x = atof(position);
              // y
              strncpy(position,&line[38],8);
              double y = atof(position);
              // z
              strncpy(position,&line[46],8);
              double z = atof(position);
              AtomicCenter ac(atomsym,
                              x*ang_to_bohr,
                              y*ang_to_bohr,
                              z*ang_to_bohr,
                              0);
              add_atom(i,ac);
              i++;
            }
        }
      fclose(fp);
    }
  else {
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
          conv = ang_to_bohr;
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
#ifdef __GNUC__
  if ((int)i < natoms) ((int)i)++;
  else i = 0;
#else
  if ((int)i < natoms) {
    int ii = (int) i;
    ii++;
    i = (Pix)ii;
  }
  else i = 0;
#endif
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
Molecule::print(SCostream& os)
{
  int i;

  os.indent();
  os << "    n  atom  label          x               y               z"
      "          mass\n";

  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);

  for (i=0; i<natom(); i++) {
      os.indent();
      os.precision(10);
      os << ' ';
      os.width(5);
      os << i+1;
      os.width(5);
      os << get_atom(i).element().symbol();
      os.width(8);
      if (get_atom(i).label()) os << get_atom(i).label();
      else os << " ";
      for (int j=0; j<3; j++) {
          os.width(16);
          os << get_atom(i)[j];
        }
      os.precision(5);
      os.width(10);
      os << get_atom(i).element().mass();
      os << endl;
    }
}

void
Molecule::print(FILE*fp)
{
  fprintf(fp,"Molecule:\n");
  fprintf(fp,"    n  atom  label          x               y               z"
             "          mass\n");

  int i;
#if defined(I860) && !defined(PARAGON)
 // there seems to  be a bug with somewhere in the iPSC version of the
 // code which causes negative numbers greater than -0.1 to be written
 // as positive, so on the iPSC write everything in sci notation
  double x;
  for (i=0; i<natom(); i++) {
      fprintf(fp," %5d%5s%8s%16.9e %16.9e %16.9e%10.5f\n",
              i+1,
              get_atom(i).element().symbol(),
              (get_atom(i).label()) ? get_atom(i).label(): " ",
              (((x=get_atom(i)[0]) < 1e-12 && x > -1e-12) ? 0.0 : x),
              (((x=get_atom(i)[1]) < 1e-12 && x > -1e-12) ? 0.0 : x),
              (((x=get_atom(i)[2]) < 1e-12 && x > -1e-12) ? 0.0 : x),
              get_atom(i).element().mass()
              );
    }
#else
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
#endif
}

AtomicCenter& Molecule::operator[](int ind)
{
  return get_atom(ind);
}

AtomicCenter& Molecule::atom(int ind)
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

const AtomicCenter& Molecule::atom(int ind) const
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
      result->add(get_atom(i).point(), get_atom(i).element().charge());
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
  atoms = 0;
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
