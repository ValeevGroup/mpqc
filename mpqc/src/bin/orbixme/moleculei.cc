
#include <strstream.h>
#include <iostream.h>
#include <moleculei.h>

#include <chemistry/molecule/molecule.h>

C_MoleculeImpl::C_MoleculeImpl()
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

