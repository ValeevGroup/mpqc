
#include <math.h>
#include <moleculei.h>

#include <chemistry/molecule/molecule.h>

C_MoleculeImpl::C_MoleculeImpl()
{
}

C_MoleculeImpl::C_MoleculeImpl(Molecule *m):
  C_KeyValCreatableImpl(m)
{
}

C_MoleculeImpl::~C_MoleculeImpl()
{
}

Molecule *
C_MoleculeImpl::mol()
{
  Molecule *ret;
  ret = Molecule::castdown(dc_);
  return ret;
}

long
C_MoleculeImpl::natom(CORBA_Environment &)
{
  Molecule *m = mol();
  if (m) return m->natom();
  else return 0;
}

double
C_MoleculeImpl::x(long i, CORBA_Environment &)
{
  Molecule *m = mol();
  if (m) return m->atom(i)[0];
  else return 0.0;
}

double
C_MoleculeImpl::y(long i, CORBA_Environment &)
{
  Molecule *m = mol();
  if (m) return m->atom(i)[1];
  else return 0.0;
}

double
C_MoleculeImpl::z(long i, CORBA_Environment &)
{
  Molecule *m = mol();
  if (m) return m->atom(i)[2];
  else return 0.0;
}

double
C_MoleculeImpl::r(long i, long j, CORBA_Environment &e)
{
  double dx = x(i,e) - x(j,e);
  double dy = y(i,e) - y(j,e);
  double dz = z(i,e) - z(j,e);
  return sqrt(dx*dx + dy*dy + dz*dz);
}

long
C_MoleculeImpl::atomic_number(long i, CORBA_Environment &)
{
  Molecule *m = mol();
  if (m) return m->atom(i).element().number();
  else return 0;
}

