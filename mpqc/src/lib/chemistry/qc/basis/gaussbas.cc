
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/files.h>

SavableState_REF_def(GaussianBasisSet);

#define CLASSNAME GaussianBasisSet
#define PARENTS public SavableState
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define VERSION 2
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
GaussianBasisSet::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

GaussianBasisSet::GaussianBasisSet(const RefKeyVal&topkeyval)
{
  molecule_ =
    Molecule::require_castdown(topkeyval
                               ->describedclassvalue("molecule").pointer(),
                               "molecule of wrong type");

  // see if the user requests pure am or cartesian functions
  int pure;
  pure = topkeyval->booleanvalue("puream");
  if (topkeyval->error() != KeyVal::OK) pure = -1;

  // construct a keyval that contains the basis library

  // this ParsedKeyVal CTOR looks at the basisdir and basisfiles
  // variables to find out what basis set files are to be read in
  RefKeyVal libkeyval = new ParsedKeyVal("basis",topkeyval);
  RefKeyVal keyval = new AggregateKeyVal(topkeyval,libkeyval);

  // Bases keeps track of what basis set data bases have already
  // been read in.  It also handles the conversion of basis
  // names to file names.
  BasisFileSet bases(keyval);
  init(molecule_,keyval,bases,1,pure);
}

GaussianBasisSet::GaussianBasisSet(StateIn&s):
  SavableState(s),
  center_to_nshell_(s)
{
  molecule_.restore_state(s);

  ncenter_ = center_to_nshell_.length();
  s.getstring(name_);

  nshell_ = 0;
  int i;
  for (i=0; i<ncenter_; i++) {
      nshell_ += center_to_nshell_(i);
    }
  
  shell = new GaussianShell*[nshell_];
  for (i=0; i<nshell_; i++) {
      shell[i] = new GaussianShell(s);
    }

  init2();
}

void
GaussianBasisSet::save_data_state(StateOut&s)
{
  center_to_nshell_.save_object_state(s);

  molecule_.save_state(s);
  s.putstring(name_);
  for (int i=0; i<nshell_; i++) {
      shell[i]->save_object_state(s);
    }
}

void
GaussianBasisSet::init(RefMolecule&molecule,
                       RefKeyVal&keyval,
                       BasisFileSet& bases,
                       int have_userkeyval,
                       int pur)
{
  int pure, havepure;
  pure = pur;
  if (pur == -1) {
      havepure = 0;
    }
  else {
      havepure = 1;
    }

  name_ = keyval->pcharvalue("name");
  int have_custom, nelement;

  if (keyval->exists("basis")) {
      have_custom = 1;
      nelement = keyval->count("element");
    }
  else {
      have_custom = 0;
      nelement = 0;
      if (!name_) {
          fprintf(stderr,"GaussianBasisSet: No name given for basis set\n");
          abort();
        }
    }

  // construct prefixes for each atom: :basis:<atom>:<basisname>:<shell#>
  // and read in the shell
  nbasis_ = 0;
  int ishell = 0;
  ncenter_ = molecule->natom();
  int iatom;
  for (iatom=0; iatom < ncenter_; iatom++) {
      ChemicalElement currentelement(molecule->operator[](iatom).element());
      // see if there is a specific basisname for this atom
      char* sbasisname = 0;
      if (have_custom && !nelement) {
          sbasisname = keyval->pcharvalue("basis",iatom);
        }
      else if (nelement) {
          int i;
          for (i=0; i<nelement; i++) {
              char *tmpelementname = keyval->pcharvalue("element", i);
              ChemicalElement tmpelement(tmpelementname);
              if (tmpelement == currentelement) {
                  sbasisname = keyval->pcharvalue("basis", i);
                  break;
                }
              delete[] tmpelementname;
            }
        }
      if (!sbasisname) {
          if (!name_) {
              fprintf(stderr,
                      "GaussianBasisSet: no basis name for atom %d (%s)\n",
                      iatom, currentelement.name());
              abort();
            }
          sbasisname = strcpy(new char[strlen(name_)+1],name_);
        }
      recursively_get_shell(ishell,keyval,
                            currentelement.name(),
                            sbasisname,bases,havepure,pure,0);
      delete[] sbasisname;
    }
  nshell_ = ishell;
  shell = new GaussianShell*[nshell_];
  ishell = 0;
  center_to_nshell_.set_length(ncenter_);
  for (iatom=0; iatom<ncenter_; iatom++) {
      ChemicalElement currentelement(molecule->operator[](iatom).element());
      // see if there is a specific basisname for this atom
      char* sbasisname = 0;
      if (have_custom && !nelement) {
          sbasisname = keyval->pcharvalue("basis",iatom);
        }
      else if (nelement) {
          int i;
          for (i=0; i<nelement; i++) {
              char *tmpelementname = keyval->pcharvalue("element", i);
              ChemicalElement tmpelement(tmpelementname);
              if (tmpelement == currentelement) {
                  sbasisname = keyval->pcharvalue("basis", i);
                  break;
                }
              delete[] tmpelementname;
            }
        }
      if (!sbasisname) {
          if (!name_) {
              fprintf(stderr,
                      "GaussianBasisSet: no basis name for atom %d (%s)\n",
                      iatom, currentelement.name());
              abort();
            }
          sbasisname = strcpy(new char[strlen(name_)+1],name_);
        }

      int ishell_old = ishell;
      recursively_get_shell(ishell,keyval,
                            currentelement.name(),
                            sbasisname,bases,havepure,pure,1);

      center_to_nshell_[iatom] = ishell - ishell_old;

      delete[] sbasisname;
     }

  // finish with the initialization steps that don't require any
  // external information
  init2();
}

double
GaussianBasisSet::r(int icenter, int xyz) const
{
  return molecule_->atom(icenter).point()[xyz];
}

void
GaussianBasisSet::init2()
{
  // center_to_shell_ and shell_to_center_
  shell_to_center_.set_length(nshell_);
  center_to_shell_.set_length(ncenter_);
  int ishell = 0;
  for (int icenter=0; icenter<ncenter_; icenter++) {
      center_to_shell_[icenter] = ishell;
      ishell += center_to_nshell_[icenter];
      for (int j = center_to_shell_[icenter]; j<icenter; j++) {
	  shell_to_center_[j] = icenter;
	}
     }

  // compute nbasis_ and shell_to_function_[]
  shell_to_function_.set_length(nshell_);
  nbasis_ = 0;
  nprim_ = 0;
  for (ishell=0; ishell<nshell_; ishell++) {
      shell_to_function_[ishell] = nbasis_;
      nbasis_ += shell[ishell]->nfunction();
      nprim_ += shell[ishell]->nprimitive();
    }

  // would like to do this in function_to_shell(), but it is const
  int n = nbasis();
  int nsh = nshell();
  function_to_shell_.set_length(n);
  int ifunc = 0;
  for (int i=0; i<nsh; i++) {
      int nfun = operator()(i).nfunction();
      for (int j=0; j<nfun; j++) {
          function_to_shell_[ifunc] = i;
          ifunc++;
        }
    }
}

void
GaussianBasisSet::
  recursively_get_shell(int&ishell,RefKeyVal&keyval,
			const char*element,
			const char*basisname,
                        BasisFileSet &bases,
			int havepure,int pure,
			int get)
{
  char keyword[KeyVal::MaxKeywordLength],prefix[KeyVal::MaxKeywordLength];

  sprintf(keyword,":basis:%s:%s",
	  element,basisname);
  int count = keyval->count(keyword);
  if (keyval->error() != KeyVal::OK) {
      keyval = bases.keyval(keyval, basisname);
    }
  count = keyval->count(keyword);
  if (keyval->error() != KeyVal::OK) {
      fprintf(stderr,"GaussianBasisSet:: couldn't find \"%s\":\n",
	      keyword);
      keyval->errortrace(cerr);
      exit(1);
    }
  if (!count) return;
  for (int j=0; j<count; j++) {
      sprintf(prefix,":basis:%s:%s",
	      element,basisname);
      RefKeyVal prefixkeyval = new PrefixKeyVal(prefix,keyval,j);
      if (prefixkeyval->exists("get")) {
          char* newbasis = prefixkeyval->pcharvalue("get");
          if (!newbasis) {
	      fprintf(stderr,"GaussianBasisSet: error processing get for \"%s\"\n",
		      prefix);
              keyval->errortrace(cerr);
	      exit(1);
	    }
	  recursively_get_shell(ishell,keyval,element,newbasis,bases,
                                havepure,pure,get);
          delete[] newbasis;
	}
      else {
          if (get) {
	      if (havepure) shell[ishell] = new GaussianShell(prefixkeyval,pure);
	      else shell[ishell] = new GaussianShell(prefixkeyval);
	    }
	  ishell++;
	}
    }
}

GaussianBasisSet::~GaussianBasisSet()
{
  delete[] name_;

  int ii;
  for (ii=0; ii<nshell_; ii++) {
      delete shell[ii];
    }
  delete[] shell;
}

int
GaussianBasisSet::max_nfunction_in_shell() const
{
  int i;
  int max = 0;
  for (i=0; i<nshell_; i++) {
      if (max < shell[i]->nfunction()) max = shell[i]->nfunction();
    }
  return max;
}

int
GaussianBasisSet::function_to_shell(int func) const
{
  return function_to_shell_[func];
}

int
GaussianBasisSet::ncenter() const
{
  return ncenter_;
}

int
GaussianBasisSet::nshell_on_center(int icenter) const
{
  return center_to_nshell_[icenter];
}

const GaussianShell&
GaussianBasisSet::operator()(int icenter,int ishell) const
{
  return *shell[center_to_shell_(icenter) + ishell];
}

GaussianShell&
GaussianBasisSet::operator()(int icenter,int ishell)
{
  return *shell[center_to_shell_(icenter) + ishell];
}

centers_t*
GaussianBasisSet::convert_to_centers_t()
{
  centers_t* c;
  c = (centers_t*)malloc(sizeof(centers_t));
  c->n = ncenter();
  c->center = (center_t*) malloc(sizeof(center_t)*ncenter());
  c->center_num = 0;
  c->shell_num = 0;
  c->func_num = 0;
  c->nfunc = 0;
  c->nshell = 0;
  c->nprim = 0;
  c->func_offset = 0;
  c->prim_offset = 0;
  c->shell_offset = 0;

  Molecule *mol = molecule_.pointer();
  
  for (int icenter = 0; icenter < ncenter_; icenter++) {
      if (mol) {
          const Molecule& molr = *mol;
          c->center[icenter].atom = strdup(molr[icenter].element().symbol());
          c->center[icenter].charge = molr[icenter].element().charge();
        }
      else {
          c->center[icenter].atom = strdup("dummy");
          c->center[icenter].charge = 0.0;
        }
      c->center[icenter].r = (double*)malloc(sizeof(double)*3);
      for (int xyz=0; xyz<3; xyz++)
        c->center[icenter].r[xyz]=r(icenter,xyz);
      c->center[icenter].basis.n = nshell_on_center(icenter);
      c->center[icenter].basis.name = strdup(name());
      c->center[icenter].basis.shell =
	(shell_t*)malloc(sizeof(shell_t)*nshell_on_center(icenter));
      shell_t*shell = c->center[icenter].basis.shell;
      int ishell;
      for (ishell = 0; ishell < center_to_nshell_[icenter]; ishell++) {
          const GaussianShell &igshell = operator()(icenter, ishell);
          int nprim = shell[ishell].nprim = igshell.nprimitive();
          int ncon = shell[ishell].ncon = igshell.ncontraction();
          shell[ishell].exp = (double*)malloc(sizeof(double)*nprim);
          shell[ishell].type
            = (shell_type_t*)malloc(sizeof(shell_type_t)*ncon);
	  shell[ishell].coef = (double**)malloc(sizeof(double*)*ncon);
          shell[ishell].norm = 0;
	  int i;
	  for (i=0; i<ncon; i++) {
	      shell[ishell].coef[i] = (double*)malloc(sizeof(double)*nprim);
	      shell[ishell].type[i].am
                = operator()(icenter,ishell).am(i);
	      shell[ishell].type[i].puream
                = operator()(icenter,ishell).is_pure(i);
	    }
	  for (i=0; i<nprim; i++) {
	      shell[ishell].exp[i]
                = operator()(icenter,ishell).exponent(i);
	      int j;
	      for (j=0; j<ncon; j++) {
		  shell[ishell].coef[j][i]
		    = operator()(icenter,ishell).coefficient_norm(j,i);
		}
	    }
	}
    }
  return c;
}

void GaussianBasisSet::print(FILE*fp) const
{
  fprintf(fp,"GaussianBasisSet: nshell = %d, nbasis = %d\n",nshell_,nbasis_);

  // Loop over centers
  int icenter = 0;
  int ioshell = 0;
  for (icenter=0; icenter < ncenter_; icenter++) {
      fprintf(fp,"center %d: %12.8f %12.8f %12.8f, nshell = %d, shellnum = %d\n",
	      icenter,
	      r(icenter,0),
	      r(icenter,1),
	      r(icenter,2),
	      center_to_nshell_[icenter],
	      center_to_shell_[icenter]);
      for (int ishell=0; ishell < center_to_nshell_[icenter]; ishell++) {
	  fprintf(fp,"Shell %d: functionnum = %d\n",
		  ishell,shell_to_function_[ioshell]);
	  operator()(icenter,ishell).print(fp);
          ioshell++;
	}
    }

}
