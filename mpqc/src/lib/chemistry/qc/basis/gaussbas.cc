
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include "gaussshell.h"
#include "gaussbas.h"

SavableState_REF_def(GaussianBasisSet);

#define CLASSNAME GaussianBasisSet
#define PARENTS virtual public SavableState
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
GaussianBasisSet::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

GaussianBasisSet::GaussianBasisSet(KeyVal&topkeyval)
{
  RefMolecule molecule =
    Molecule::require_castdown(topkeyval
                               .describedclassvalue("molecule").pointer(),
                               "molecule of wrong type");
  
  char* basisname = topkeyval.pcharvalue("name");

  // see if the user requests pure am or cartesian functions
  int pure;
  pure = topkeyval.booleanvalue("puream");
  if (topkeyval.error() != KeyVal::OK) pure = -1;

  // construct a keyval that contains the basis library

  // this ParsedKeyVal CTOR looks at the basisdir and basisfiles
  // variables to find out what basis set files are to be read in
  ParsedKeyVal libkeyval("basis",topkeyval); libkeyval.unmanage();
  AggregateKeyVal keyval(topkeyval,libkeyval); keyval.unmanage();

  init(molecule,keyval,basisname,1,pure);

  delete[] basisname;
}

GaussianBasisSet::GaussianBasisSet(RefMolecule&molecule,
                                   const char*basisname,
                                   int pure)
{
  // construct a keyval that contains the basis library
  ParsedKeyVal libkeyval(BASISDIR "/v2g90.ipv2"); libkeyval.unmanage();
  libkeyval.read(BASISDIR "/v2g90supp.ipv2");

  init(molecule,libkeyval,basisname,0,pure);
}

GaussianBasisSet::GaussianBasisSet(StateIn&s):
  SavableState(s,class_desc_),
  center_to_r_(s),
  center_to_nshell_(s)
{
  ncenter_ = center_to_nshell_.length();
  s.getstring(name_);

  nshell_ = 0;
  for (int i=0; i<ncenter_; i++) {
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
  center_to_r_.save_object_state(s);
  center_to_nshell_.save_object_state(s);
  s.putstring(name_);
  for (int i=0; i<nshell_; i++) {
      shell[i]->save_object_state(s);
    }
}

void
GaussianBasisSet::init(RefMolecule&molecule,
                       KeyVal&keyval,
                       const char*basisname,
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

  name_ = strcpy(new char[strlen(basisname)+1], basisname);

  // construct prefixes for each atom: :basis:<atom>:<basisname>:<shell#>
  // and read in the shell
  nbasis_ = 0;
  int ishell = 0;
  ncenter_ = molecule->natom();
  for (int iatom=0; iatom < ncenter_; iatom++) {
      // see if there is a specific basisname for this atom
      char* sbasisname = 0;
      if (have_userkeyval) {
          sbasisname = keyval.pcharvalue("basis",iatom);
        }
      if (!sbasisname) {
          sbasisname = strcpy(new char[strlen(basisname)+1],basisname);
        }
      recursively_get_shell(ishell,keyval,
                            molecule->operator[](iatom).element().name(),
			    sbasisname,havepure,pure,0);
      delete[] sbasisname;
    }
  nshell_ = ishell;
  shell = new GaussianShell*[nshell_];
  ishell = 0;
  center_to_nshell_.set_length(ncenter_);
  for (iatom=0; iatom<ncenter_; iatom++) {
      // see if there is a specific basisname for this atom
      char* sbasisname = 0;
      if (have_userkeyval) {
          sbasisname = keyval.pcharvalue("basis",iatom);
        }
      if (!sbasisname) {
          sbasisname = strcpy(new char[strlen(basisname)+1],basisname);
        }

      int ishell_old = ishell;
      recursively_get_shell(ishell,keyval,
                            molecule->operator[](iatom).element().name(),
			    sbasisname,havepure,pure,1);

      center_to_nshell_[iatom] = ishell - ishell_old;

      delete[] sbasisname;
     }

  // compute center_to_r_
  center_to_r_.set_lengths(ncenter_,3);
  ishell = 0;
  for (iatom=0; iatom<ncenter_; iatom++) {
      center_to_r_(iatom,0) = molecule->operator[](iatom)[0];
      center_to_r_(iatom,1) = molecule->operator[](iatom)[1];
      center_to_r_(iatom,2) = molecule->operator[](iatom)[2];
     }

  // finish with the initialization steps that don't require any
  // external information
  init2();
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
  recursively_get_shell(int&ishell,KeyVal&keyval,
			const char*element,
			const char*basisname,
			int havepure,int pure,
			int get)
{
  char keyword[KeyVal::MaxKeywordLength],prefix[KeyVal::MaxKeywordLength];

  sprintf(keyword,":basis:%s:%s",
	  element,basisname);
  int count = keyval.count(keyword);
  if (keyval.error() != KeyVal::OK) {
      fprintf(stderr,"GaussianBasisSet:: couldn't find \"%s\":\n",
	      keyword);
      keyval.errortrace(stderr);
      exit(1);
    }
  if (!count) return;
  for (int j=0; j<count; j++) {
      sprintf(prefix,":basis:%s:%s",
	      element,basisname);
      PrefixKeyVal prefixkeyval(prefix,keyval,j); prefixkeyval.unmanage();
      if (prefixkeyval.exists("get")) {
          char* newbasis = prefixkeyval.pcharvalue("get");
          if (!newbasis) {
	      fprintf(stderr,"GaussianBasisSet: error processing get for \"%s\"\n",
		      prefix);
              keyval.errortrace(stderr);
	      exit(1);
	    }
	  recursively_get_shell(ishell,keyval,element,newbasis,havepure,pure,get);
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
GaussianBasisSet::convert_to_centers_t(const Molecule*mol) const
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
        c->center[icenter].r[xyz]=center_to_r_(icenter,xyz);
      c->center[icenter].basis.n = nshell_on_center(icenter);
      c->center[icenter].basis.name = strdup(name());
      c->center[icenter].basis.shell =
	(shell_t*)malloc(sizeof(shell_t)*nshell_on_center(icenter));
      shell_t*shell = c->center[icenter].basis.shell;
      int ishell = 0;
      for (ishell = 0; ishell < center_to_nshell_[icenter]; ishell++) {
	  int nprim = shell[ishell].nprim
            = operator()(icenter,ishell).nprimitive();
          int ncon = shell[ishell].ncon
            = operator()(icenter,ishell).ncontraction();
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
	      center_to_r_[icenter][0],
	      center_to_r_[icenter][1],
	      center_to_r_[icenter][2],
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
