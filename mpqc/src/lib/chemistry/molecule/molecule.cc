
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <string.h>

#include <util/misc/formio.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/localdef.h>
#include <math/scmat/cmatrix.h>

SavableState_REF_def(Molecule);

#define CLASSNAME Molecule
#define PARENTS public SavableState
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
  if (atoms) delete[] atoms;
  atoms=0;
  natoms=0;
}

Molecule::Molecule(const RefKeyVal&input) :
  pg(input), atoms(0), natoms(0)
{
  const double ang_to_bohr = 1.0/0.52917706;

  if (input->exists("pdb_file")) {
      const int LineLength = 85;
      char line[LineLength];
      char* filename = input->pcharvalue("pdb_file");
      FILE*fp = fopen(filename,"r");
      if (!fp) {
          fprintf(stderr,
                  "Molecule::Molecule(const RefKeyVal&input): "
                  "pdb file not found: \"%s\"\n", filename);
          abort();
        }
      int i=0;
      while(fgets(line,LineLength,fp)) {
          if (!strncmp(line,"HETA",4) || !strncmp(line,"ATOM",4)) {
              char atomsym[3];
              // find the atomic symbol
              int symletter=0, offset;
              for (offset=12; offset<16 && symletter<2; offset++) {
                  if (line[offset] != ' '
                      && (line[offset] < '0' || line[offset] > '9')) {
                      atomsym[symletter] = line[offset];
                      symletter++;
                    }
                }
              atomsym[symletter] = '\0';
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
    // first let's see if the input is in bohr or angstrom units
      int aangstroms = input->booleanvalue("angstrom");
      if (input->error() != KeyVal::OK) {
        aangstroms = input->booleanvalue("aangstrom");
      }
      if (input->error() != KeyVal::OK) {
        aangstroms = input->booleanvalue("angstroms");
      }
      if (input->error() != KeyVal::OK) {
        aangstroms = input->booleanvalue("aangstroms");
      }
        
      double conv = 1.0;
      if (aangstroms) {
          conv = ang_to_bohr;
        }

      // get the number of atoms and make sure that the geometry and the
      // atoms array have the same number of atoms.
      // right now we read in the unique atoms...then we will symmetrize.
      // the length of atoms must still equal the length of geometry, but
      // we'll try to set up atom_labels such that different lengths are
      // possible
      int natom = input->count("geometry");
      if (natom != input->count("atoms")) {
        fprintf(stderr,"Molecule: size of \"geometry\" != size of \"atoms\"\n");
        return;
      }

      int i;
      for (i=0; i<natom; i++) {
          char* name;
          char* labels;
          AtomicCenter ac(name = input->pcharvalue("atoms",i),
                          input->doublevalue("geometry",i,0)*conv,
                          input->doublevalue("geometry",i,1)*conv,
                          input->doublevalue("geometry",i,2)*conv,
                          labels = input->pcharvalue("atom_labels",i)
                          );
          if (name) delete[] name;
          if (labels) delete[] labels;
          add_atom(i,ac);
        }
    }
    
  // we'll assume that only unique atoms are given in the input unless
  // told otherwise
  if (!input->booleanvalue("redundant_atoms"))
    symmetrize();
}

Molecule& Molecule::operator=(Molecule& mol)
{
  if(atoms) delete[] atoms; atoms=0;
  natoms=0;

  pg = mol.pg;

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
  return (Pix)atoms;
}

void Molecule::next(Pix& i)
{
  if ((AtomicCenter*)i < &atoms[natoms-1]) {
      i = (Pix) &((AtomicCenter*)i)[1];
    }
  else i = 0;
}

int Molecule::owns(Pix i)
{
  if (i >= (Pix)atoms && i <= (Pix)&atoms[natoms-1]) return 1;
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
Molecule::print(ostream& os)
{
  int i;

  os << indent;
  os << "    n  atom  label          x               y               z"
      "          mass\n";

  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);

  for (i=0; i<natom(); i++) {
      os << indent;
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

AtomicCenter& Molecule::atom(int ind)
{
  return get_atom(ind);
}

AtomicCenter& Molecule::operator()(Pix pix)
{
  return *(AtomicCenter*)pix;
}

const AtomicCenter& Molecule::operator[](int ind) const
{
  return get_atom(ind);
}

const AtomicCenter& Molecule::atom(int ind) const
{
  return get_atom(ind);
}

int
Molecule::atom_label_to_index(const char *label) const
{
  int i;
  for (i=0; i<natom(); i++) {
      if (!strcmp(label,atom(i).label())) return i;
    }
  return -1;
}

const AtomicCenter& Molecule::operator()(Pix pix) const
{
  return *(AtomicCenter*)pix;
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
  pg.save_object_state(so);
  so.put(natoms);
  for (int i=0; i < natoms; i++) {
      get_atom(i).save_object_state(so);
    }
}

Molecule::Molecule(StateIn& si):
  atoms(0),
  natoms(0),
  SavableState(si)
{
  PointGroup tpg(si);
  pg=tpg;

  int natom;
  si.get(natom);
  for (int i=0; i < natom; i++) {
      AtomicCenter ac(si);
      add_atom(i,ac);
    }
}

void
Molecule::set_point_group(const PointGroup&ppg)
{
  pg = ppg;
  symmetrize();
}

const PointGroup& Molecule::point_group() const
{
  return pg;
}

RefPoint Molecule::center_of_mass()
#ifdef __GNUC__
  return ret;
#endif
{
#ifndef __GNUC__
  RefPoint ret;
#endif
  double X,Y,Z,M;

  X = Y = Z = M = 0;

  for (int i=0; i < natom(); i++) {
    double m = atom(i).element().mass();
    X += m * atom(i)[0];
    Y += m * atom(i)[1];
    Z += m * atom(i)[2];
    M += m;
  }

  X /= M;
  Y /= M;
  Z /= M;

  //printf("center of mass = %f %f %f\n",X,Y,Z);

  ret = new Point;
  ret->operator[](0) = X;
  ret->operator[](1) = Y;
  ret->operator[](2) = Z;

#ifndef __GNUC__
  return ret;
#endif
}

double
Molecule::nuclear_repulsion_energy()
{
  int i, j;
  double r, e=0;

  for (i=1; i < natoms; i++) {
    AtomicCenter& ai = get_atom(i);
    double Zi = ai.element().charge();
    
    for (j=0; j < i; j++) {
      AtomicCenter& aj = get_atom(j);
      e += Zi * aj.element().charge() / dist(ai.point(), aj.point());
    }
  }

  return e;
}

void
Molecule::nuclear_repulsion_1der(int center, double xyz[3])
{
  int i,j,k;
  double r[3],r2;
  double factor;

  xyz[0] = 0.0;
  xyz[1] = 0.0;
  xyz[2] = 0.0;
  for (i=0; i < natoms; i++) {
    AtomicCenter& ai = get_atom(i);
    double Zi = ai.element().charge();

    for (j=0; j < i; j++) {
      if (center==i || center==j) {
        AtomicCenter& aj = get_atom(j);

        r2 = 0.0;
        for (k=0; k < 3; k++) {
          r[k] = ai[k] - aj[k];
          r2 += r[k]*r[k];
        }
        
        factor = - Zi * aj.element().charge() * pow(r2,-1.5);
        if (center==j) factor = -factor;
        for (k=0; k<3; k++) {
          xyz[k] += factor * r[k];
        }
      }
    }
  }
}

int
Molecule::atom_at_position(double *v, double tol)
{
  Point p(v,3);
  for (int i=0; i < natom(); i++) {
    if (dist(p,atom(i).point()) < tol)
      return i;
  }
  return -1;
}

///////////////////////////////////////////////////////////////////////////

// pass in natoms if you don't want to search through the entire molecule
static int
is_unique(Point& p, Molecule *mol, int natoms)
{
  for (int i=natoms-1; i >= 0; i--) {
    if (dist(p,mol->atom(i).point()) < 0.5) {
      return 0;
    }
  }
  return 1;
}

static int
is_unique(Point& p, Molecule *mol)
{
  for (int i=0; i < mol->natom(); i++) {
    if (dist(p,mol->atom(i).point()) < 0.5) {
      return 0;
    }
  }

  return 1;
}

// We are given a molecule which may or may not have just the symmetry
// distinct atoms in it.  We have to go through the existing set of atoms,
// perform each symmetry operation in the point group on each of them, and
// then add the new atom if it isn't in the list already

void
Molecule::symmetrize()
{
  // if molecule is c1, don't do anything
  if (!strcmp(this->point_group().symbol(),"c1")) {
    return;
  }

  Molecule *newmol = new Molecule;
  newmol->pg = this->pg;
  
  CharacterTable ct = this->point_group().char_table();

  Point np;
  SymmetryOperation so;
  int nnew=0;

  // first off, copy the un-symmetrized molecule into the new one
  int i;
  for (i=0; i < this->natom(); i++) {
    newmol->add_atom(nnew,this->atom(i));
    nnew++;
  }
  
  for (i=0; i < this->natom(); i++) {
    AtomicCenter ac = this->atom(i);

    for (int g=0; g < ct.order(); g++) {
      so = ct.symm_operation(g);
      for (int ii=0; ii < 3; ii++) {
        np[ii]=0;
        for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * ac[jj];
      }

      if (is_unique(np,newmol)) {
        AtomicCenter nac(ac.element().symbol(),np[0],np[1],np[2],ac.label());
        newmol->add_atom(nnew,nac);
        nnew++;
      }
    }
  }
  
  *this = *newmol;
  delete newmol;
}

// move the molecule to the center of mass
void
Molecule::move_to_com()
{
  RefPoint com = center_of_mass();

  double X = com->operator[](0);
  double Y = com->operator[](1);
  double Z = com->operator[](2);

  for (int i=0; i < natom(); i++) {
    atom(i)[0] -= X;
    atom(i)[1] -= Y;
    atom(i)[2] -= Z;
  }
}

// find the 3 principal coordinate axes, and rotate the molecule to be 
// aligned along them.  also rotate the symmetry frame contained in point_group
void
Molecule::transform_to_principal_axes(int trans_frame)
{
  const double au_to_angs = 0.2800283608302436;

  // mol_move_to_com(mol);

  double *inert[3], *evecs[3];
  double evals[3];

  int i;
  for (i=0; i < 3; i++) {
    inert[i] = new double[3];
    evecs[i] = new double[3];
    memset(inert[i],'\0',sizeof(double)*3);
    memset(evecs[i],'\0',sizeof(double)*3);
  }

  AtomicCenter ac;
  for (i=0; i < natom(); i++) {
    ac = atom(i);
    double m=au_to_angs*ac.element().mass();
    inert[0][0] += m * (ac[1]*ac[1] + ac[2]*ac[2]);
    inert[1][0] -= m * ac[0]*ac[1];
    inert[1][1] += m * (ac[0]*ac[0] + ac[2]*ac[2]);
    inert[2][0] -= m * ac[0]*ac[2];
    inert[2][1] -= m * ac[1]*ac[2];
    inert[2][2] += m * (ac[0]*ac[0] + ac[1]*ac[1]);
  }

  inert[0][1] = inert[1][0];
  inert[0][2] = inert[2][0];
  inert[1][2] = inert[2][1];

 // cleanup inert
  for (i=0; i < 3; i++) {
    for (int j=0; j <= i; j++) {
      if (fabs(inert[i][j]) < 1.0e-5) {
        inert[i][j]=inert[j][i]=0.0;
      }
    }
  }

  cmat_diag(inert, evals, evecs, 3, 1, 1e-14);

 // cleanup evecs
  for (i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      if (fabs(evecs[i][j]) < 1.0e-5) {
        evecs[i][j]=0.0;
      }
    }
  }

  double x,y,z;
  for (i=0; i < natom(); i++) {
    x = atom(i)[0]; y = atom(i)[1]; z = atom(i)[2];

    atom(i)[0] = evecs[0][0]*x + evecs[1][0]*y + evecs[2][0]*z;
    atom(i)[1] = evecs[0][1]*x + evecs[1][1]*y + evecs[2][1]*z;
    atom(i)[2] = evecs[0][2]*x + evecs[1][2]*y + evecs[2][2]*z;
  }

  if (!trans_frame) return;
  
  SymmetryOperation tso=point_group().symm_frame();

  for (i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      double t=0;
      for (int k=0; k < 3; k++) t += tso[i][k]*evecs[k][j];
      if (fabs(t) < 0.5)
        t = 0;
      else if (fabs(t) >= .5)
        t = 1;
      
      pg.symm_frame()[i][j] = t;
    }
  }
}

// returns an array containing indices of the unique atoms
int *
Molecule::find_unique_atoms()
{
  // if this is a c1 molecule, then return all indices
  if (!strcmp(point_group().symbol(),"c1")) {
    int * ret = new int[natom()];
    for (int i=0; i < natom(); i++) ret[i]=i;
    return ret;
  }

  int nuniq = num_unique_atoms();
  int * ret = new int[nuniq];

  // so that we don't have side effects, copy this to mol
  RefMolecule mol=new Molecule(*this);

  // the first atom is always unique
  ret[0]=0;

//   mol->transform_to_principal_axes(0);
//   mol->symmetrize();

  CharacterTable ct = mol->point_group().char_table();

  AtomicCenter ac;
  SymmetryOperation so;
  Point np;

  int nu=1;
  int i;
  for (i=1; i < mol->natom(); i++) {
    ac = mol->atom(i);
    int i_is_unique=1;

    // subject i to all symmetry ops...if one of these maps into an atom
    // further down in atoms, then break out
    for (int g=0; g < ct.order(); g++) {
      so = ct.symm_operation(g);
      for (int ii=0; ii < 3; ii++) {
        np[ii]=0;
        for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * ac[jj];
      }

      if (!is_unique(np,mol.pointer(),i)) {
        i_is_unique=0;
        break;
      }
    }
    if (i_is_unique) {
      ret[nu]=i;
      nu++;
    }
  }

 // now let's go through the unique atoms and see if they have any zero
 // coordinates.  if not, then see if any of the equivalent atoms do.
 // change the coordinate number if so
  double tol=1.0e-5;
  for (i=0; i < nuniq; i++) {
    ac = mol->atom(ret[i]);
    if (fabs(ac[0]) < tol || fabs(ac[1]) < tol || fabs(ac[2]) < tol)
      continue;

    for (int g=1; g < ct.order(); g++) {
      int breakg=0;
      so = ct.symm_operation(g);
      for (int ii=0; ii < 3; ii++) {
        np[ii]=0;
        for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * ac[jj];
      }

     // if np has a zero coordinate then find it in atoms
      if (fabs(np[0]) < tol || fabs(np[1]) < tol || fabs(np[2]) < tol) {
        for (int j=0; j < mol->natom(); j++) {
          if (dist(np,mol->atom(j).point()) < 0.1) {
            ret[i] = j;
            breakg=1;
            break;
          }
        }
      }
    if (breakg) break;
    }
  }

  return ret;
}

int
Molecule::num_unique_atoms()
{
 // if this is a c1 molecule, then return natom
  if (!strcmp(point_group().symbol(),"c1"))
    return natom();

 // so that we don't have side effects, copy this to mol
  RefMolecule mol = new Molecule(*this);

//   mol_transform_to_principal_axes(mol,0);
//   mol->symmetrize();

  AtomicCenter ac;
  SymmetryOperation so;
  Point np;

  CharacterTable ct = mol->point_group().char_table();

  int nu=1;

  for (int i=1; i < mol->natom(); i++) {
    ac = mol->atom(i);
    int i_is_unique=1;

    // subject i to all symmetry ops...if one of these maps into an atom
    // further down in atoms, then break out
    for (int g=0; g < ct.order(); g++) {
      so = ct.symm_operation(g);
      for (int ii=0; ii < 3; ii++) {
        np[ii]=0;
        for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * ac[jj];
      }

      if (!is_unique(np,mol.pointer(),i)) {
        i_is_unique=0;
        break;
      }
    }
    if (i_is_unique) nu++;
  }

  return nu;
}

// I hate small numbers...they make my teeth itch.  Let's go through
// all the atoms, and find coordinates which are close to zero ( < 1.0e-5)
// and set these to zero.  this is cheating, but who cares

static void
get_rid_of_annoying_numbers(Molecule* mol)
{
  for (int i=0; i < mol->natom(); i++) {
    for (int j=0; j < 3; j++) {
      if (fabs(mol->atom(i)[j]) < 5.0e-2) mol->atom(i).point()[j]=0;
    }
  }
}

// given a molecule, make sure that equivalent centers have coordinates
// that really map into each other

void
Molecule::cleanup_molecule()
{
  // this may have already been done, but let's first move to the principal
  // axes
//   transform_to_principal_axes(0);
//   symmetrize();

  // if symmetry is c1, do nothing else
  if (!strcmp(point_group().symbol(),"c1")) return;

  get_rid_of_annoying_numbers(this);

  // now let's find out how many unique atoms there are and who they are
  int nuniq = num_unique_atoms();
  int *uniq = find_unique_atoms();

  Point up,np;
  SymmetryOperation so;
  CharacterTable ct = point_group().char_table();

  // grab the coordinates of each unique atom and stuff into up
  for (int i=0; i < nuniq; i++) {
    up = atom(uniq[i]).point();

   // subject up to all symmetry ops...find the atom this maps to and
   // reset it's coordinates. skip E.
    for (int g=1; g < ct.order(); g++) {
      so = ct.symm_operation(g);
      for (int ii=0; ii < 3; ii++) {
        np[ii]=0;
        for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * up[jj];
      }

      for (int j=0; j < natom(); j++) {
        if (dist(np,atom(j).point()) < 0.1) {
          //printf("gamma(%d)*atom(%d)(%f,%f,%f) = atom(%d)(%f,%f,%f)\n",
          //                  g,uniq[i],up[0],up[1],up[2],j,np[0],np[1],np[2]);

          atom(j)[0] = np[0];
          atom(j)[1] = np[1];
          atom(j)[2] = np[2];
          break;
        }
      }
    }
  }

 // one last pass to make me happy
  get_rid_of_annoying_numbers(this);
}

///////////////////////////////////////////////////////////////////
// Compute the principal axes and the principal moments of inertia
///////////////////////////////////////////////////////////////////

void
Molecule::principal_moments_of_inertia(double *evals, double **evecs)
{

  // The principal moments of inertia are computed in amu*angstrom^2
  // evals: principal moments of inertia
  // evecs: principal axes (optional argument)

  const double au_to_angs = 0.2800283608302436; // for moments of inertia

  double *inert[3];  // inertia tensor

  int i, j;
  int delete_evecs = 0;

  // (allocate and) initialize evecs, evals, and inert
  if (!evecs) {
    evecs = new double*[3];
    for (i=0; i<3; i++) evecs[i] = new double[3];
    delete_evecs = 1;
    }
  for (i=0; i<3; i++) {
    inert[i] = new double[3];
    memset(inert[i],'\0',sizeof(double)*3);
    memset(evecs[i],'\0',sizeof(double)*3);
    }
  memset(evals,'\0',sizeof(double)*3);

  // translate molecule so origin becomes the center of mass
  center_of_mass();
  move_to_com();

  // compute inertia tensor
  AtomicCenter ac;
  for (i=0; i<natom(); i++) {
    ac = atom(i);
    double m=au_to_angs*ac.element().mass();
    inert[0][0] += m * (ac[1]*ac[1] + ac[2]*ac[2]);
    inert[1][0] -= m * ac[0]*ac[1];
    inert[1][1] += m * (ac[0]*ac[0] + ac[2]*ac[2]);
    inert[2][0] -= m * ac[0]*ac[2];
    inert[2][1] -= m * ac[1]*ac[2];
    inert[2][2] += m * (ac[0]*ac[0] + ac[1]*ac[1]);
    }
  inert[0][1] = inert[1][0];
  inert[0][2] = inert[2][0];
  inert[1][2] = inert[2][1];

  cmat_diag(inert, evals, evecs, 3, 1, 1e-14);

  if (delete_evecs) {
    for (i=0; i<3; i++) delete[] evecs[i];
    delete[] evecs;
    }
  for (i=0; i<3; i++) {
    delete[] inert[i];
    }
}
