
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/local.h>
#include <chemistry/molecule/localdef.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/petite.h>

////////////////////////////////////////////////////////////////////////////

class Bfn {
  private:
    int _num;
    Point a;

  public:
    int& num() { return _num; }
    double& operator[](int i) { return a[i]; }
    Point transform(const SymmetryOperation& so) const;
    void print() { printf("%5d ",_num); a.print(); }
};
    
inline Point
Bfn::transform(const SymmetryOperation& so) const
{
  Point ret;
  for (int ii=0; ii < 3; ii++) {
    ret[ii] = 0;
    for (int jj=0; jj < 3; jj++)
      ret[ii] += so(ii,jj) * a[jj];
  }

  return ret;
}

class Shell {
  private:
    int _num;
    int _am;
    int _nbf;
    Bfn *bfns;
  public:
    Shell(int m, int a, int n) : _num(m), _am(a), _nbf(n) { bfns = new Bfn[_nbf]; }
    ~Shell() { delete[] bfns; }

    int am() const { return _am; }
    int nbf() const { return _nbf; }
    int num() const { return _num; }
    Bfn& operator[](int i) { return bfns[i]; }
    void print() {
      printf("  %5d am = %d:\n",_num,_am);
      for (int i=0; i < _nbf; i++) {
        printf("    ");
        bfns[i].print();
      }
    }
};
    
class Atom {
  private:
    int _nsh;
    Shell **_shells;
  public:
    Atom(int n) : _nsh(n) { _shells = new Shell*[_nsh]; }

    int nsh() { return _nsh; }
    Shell& operator[](int i) { return *_shells[i]; }
    Shell*& operator()(int i) { return _shells[i]; }

    void print() {
      for (int i=0; i < _nsh; i++)
        _shells[i]->print();
    }
};
    

////////////////////////////////////////////////////////////////////////////

PetiteList::PetiteList()
{
  _molecule=0;
  _gbs=0;
  _natom=0;
  _nshell=0;
  _ng=0;
  _p1=0;
  _atom_map=0;
  _shell_map=0;
  _lamij=0;
}

PetiteList::PetiteList(const RefMolecule& mol, const RefGaussianBasisSet &gbs)
{
  init(mol,gbs);
}

PetiteList::~PetiteList()
{
  if (_p1)
    delete[] _p1;

  if (_lamij)
    delete[] _lamij;

  if (_atom_map) {
    for (int i=0; i < _natom; i++)
      delete[] _atom_map[i];
    delete[] _atom_map;
  }

  if (_shell_map) {
    for (int i=0; i < _nshell; i++)
      delete[] _shell_map[i];
    delete[] _shell_map;
  }

  _molecule=0;
  _gbs=0;
  _natom=0;
  _nshell=0;
  _ng=0;
  _p1=0;
  _atom_map=0;
  _shell_map=0;
  _lamij=0;
}

static int
atom_num(Point& p, Molecule& mol)
{
  for (int i=0; i < mol.natom(); i++) {
    if (dist(p,mol.atom(i).point()) < 0.05)
      return i;
  }
  return -1;
}

void
PetiteList::init(const RefMolecule& _mol, const RefGaussianBasisSet &_gb)
{
  _molecule = _mol;
  _gbs = _gb;
    
  CharacterTable ct = _molecule->point_group().char_table();
  Molecule& mol = *_molecule.pointer();
  GaussianBasisSet& gbs = *_gbs.pointer();
  
  _ng = ct.order();
  
  _natom = mol.natom();
  _nshell = gbs.nshell();

  _p1 = new char[_nshell];
  _lamij = new char[ioff(_nshell)];

  _atom_map = new int*[_natom];
  for (int i=0; i < _natom; i++)
    _atom_map[i] = new int[_ng];
  
  _shell_map = new int*[_nshell];
  for (int i=0; i < _nshell; i++)
    _shell_map[i] = new int[_ng];
  
  // set up atom and shell mappings
  Point np;
  SymmetryOperation so;
  
  // loop over all centers
  for (int i=0; i < _natom; i++) {
    AtomicCenter ac = mol.atom(i);

    // then for each symop in the pointgroup, transform the coordinates of
    // center "i" and see which atom it maps into
    for (int g=0; g < _ng; g++) {
      so = ct.symm_operation(g);

      for (int ii=0; ii < 3; ii++) {
        np[ii] = 0;
        for (int jj=0; jj < 3; jj++)
          np[ii] += so(ii,jj) * ac[jj];
      }

      _atom_map[i][g] = atom_num(np,mol);
    }

    // hopefully, shells on equivalent centers will be numbered in the same
    // order
    for (int s=0; s < gbs.nshell_on_center(i); s++) {
      int shellnum = gbs.shell_on_center(i,s);
      for (int g=0; g < _ng; g++) {
        _shell_map[shellnum][g] = gbs.shell_on_center(_atom_map[i][g],s);
      }
    }
  }

  memset(_p1,0,_nshell);
  memset(_lamij,0,ioff(_nshell));
  
  // now we do _p1 and _lamij
  for (int i=0; i < _nshell; i++) {
    int leave=0;

    // we want the highest numbered shell in a group of equivalent shells
    for (int g=0; g < _ng; g++)
      if (_shell_map[i][g] > i)
        leave=1;
    
    if (leave)
      continue;
    
    _p1[i] = 1;

    for (int j=0; j <= i; j++) {
      int ij = ioff(i)+j;
      int nij = 0;

      leave=0;
      for (int g=0; g < _ng; g++) {
        int gi = _shell_map[i][g];
        int gj = _shell_map[j][g];
        int gij = ioff(gi,gj);
        if (gij > ij) {
          leave=1;
          break;
        }
        else if (gij == ij)
          nij++;
      }

      if (leave)
        continue;

      _lamij[ij] = (char) (_ng/nij);
    }
  }

  _function_map = new int*[gbs.nbasis()];
  for (int i=0; i < gbs.nbasis(); i++)
    _function_map[i] = new int[_ng];
  
  for (int i=0; i < mol.natom(); i++) {
    for (int si=0; si < gbs.nshell_on_center(i); si++) {
      int shell_i = gbs.shell_on_center(i,si);
      
      for (int fi=0; fi < gbs(i,si).nfunction(); fi++) {
        int func_i = gbs.shell_to_function(shell_i)+fi;

        for (int g=0; g < _ng; g++) {
          int j=_atom_map[i][g];
          _function_map[func_i][g] =
            gbs.shell_to_function(gbs.shell_on_center(j,si))+fi;
        }
      }
    }
  }

  // try to see how functions will transform under symops
  Atom **atoms = new Atom*[_natom];
  
  int gns=0;
  for (int i=0; i < _natom; i++) {
    int ns=0;
    for (int si=0; si < gbs.nshell_on_center(i); si++)
      ns += gbs(i,si).ncontraction();

    atoms[i] = new Atom(ns);
    Atom& ati = *atoms[i];

    ns=0;
    for (int si=0; si < gbs.nshell_on_center(i); si++) {
      int shell_i = gbs.shell_on_center(i,si);
      GaussianShell& gsi = gbs(i,si);
      int func_i = gbs.shell_to_function(shell_i);
      
      int nf=0;
      for (int nc=0; nc < gsi.ncontraction(); nc++) {
        int amc = gsi.am(nc);
        CartesianIter j(amc);
        
        ati(ns) = new Shell(shell_i,amc,gsi.nfunction(nc));
                                  
        Shell& sh = ati[ns];
        
        int f=0;
        for (j.start(); j; j.next()) {
          sh[f].num() = func_i+nf;
          sh[f][0] = (double) j.a();
          sh[f][1] = (double) j.b();
          sh[f][2] = (double) j.c();
          f++;
          nf++;
        }
        ns++;
        gns++;
      }
    }

    printf("atom %d:\n",i+1);
    atoms[i]->print();
  }
  printf("\n");
}

void
PetiteList::print(FILE *o)
{
  fprintf(o,"PetiteList:\n");
  fprintf(o,"  _natom = %d\n",_natom);
  fprintf(o,"  _nshell = %d\n",_nshell);
  fprintf(o,"  _ng = %d\n",_ng);

  fprintf(o,"\n");
  fprintf(o,"  _atom_map = \n");
  for (int i=0; i < _natom; i++) {
    fprintf(o,"    ");
    for (int g=0; g < _ng; g++)
      fprintf(o,"%5d ",_atom_map[i][g]);
    fprintf(o,"\n");
  }

  fprintf(o,"\n");
  fprintf(o,"  _shell_map = \n");
  for (int i=0; i < _nshell; i++) {
    fprintf(o,"    ");
    for (int g=0; g < _ng; g++)
      fprintf(o,"%5d ",_shell_map[i][g]);
    fprintf(o,"\n");
  }

  fprintf(o,"\n");
  fprintf(o,"  _function_map = \n");
  for (int i=0; i < _gbs->nbasis(); i++) {
    fprintf(o,"    ");
    for (int g=0; g < _ng; g++)
      fprintf(o,"%5d ",_function_map[i][g]);
    fprintf(o,"\n");
  }

  fprintf(o,"\n");
  fprintf(o,"  _p1 = \n");
  for (int i=0; i < _nshell; i++)
    fprintf(o,"    %5d\n",_p1[i]);
    
  fprintf(o,"  _lamij = \n");
  for (int i=0; i < _nshell; i++) {
    fprintf(o,"    ");
    for (int j=0; j <= i; j++)
      fprintf(o,"%5d ",_lamij[ioff(i)+j]);
    fprintf(o,"\n");
  }
  
}
