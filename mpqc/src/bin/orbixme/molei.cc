
#include <strstream.h>
#include <iostream.h>
#include <molei.h>

#include <chemistry/molecule/energy.h>

// force linkage of MPSCF
#include <chemistry/qc/mpqc/mpqc.h>
const ClassDesc &fl0 = MPSCF::class_desc_;

MolEImpl::MolEImpl(const char *name, CORBA(LoaderClass)* ld):
  MolEBOAImpl(name, ld)
{
  mole_ = 0;
}

void
MolEImpl::create(const char *s, CORBA_Environment &IT_env)
{
  cout << "create: " << s << endl;

  istrstream in(s);

  RefKeyVal keyval = new ParsedKeyVal(in);

  RefMolecularEnergy mole = keyval->describedclassvalue("mole");
  if (mole.null()) {
      cout << "create failed" << endl;
      return;
    }
  else {
      cout << "created a \"" << mole->class_name() << "\"" << endl;
    }

  mole_ = mole.pointer();
  mole_->reference();
}

MolEImpl::~MolEImpl()
{
  mole_->dereference();
  if (mole_->nreference() == 0) delete mole_;
}

double
MolEImpl::energy(CORBA_Environment &)
{
  if (!mole_) return 0.0;
  return mole_->energy();
}

