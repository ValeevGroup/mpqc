
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include "gaussshell.h"
#include "gaussbas.h"

GaussianBasisSet::GaussianBasisSet(KeyVal&topkeyval,Molecule&molecule)
{
  char* basisname = topkeyval.pcharvalue("basis");
  char keyword[KeyVal::MaxKeywordLength],prefix[KeyVal::MaxKeywordLength];
  int havepure,pure;

  name_ = new char[strlen(basisname)+1];
  strcpy(name_,basisname);

  // see if the user requests pure am or cartesian functions
  pure = topkeyval.booleanvalue("puream");
  if (topkeyval.error() == KeyVal::OK) havepure = 1;

  // construct a keyval that contains the basis library
  ParsedKeyVal lib1keyval("/net/munin/usr/people/cljanss/src/SC/lib/v2g90.ipv2");
  ParsedKeyVal lib2keyval("/net/munin/usr/people/cljanss/src/SC/lib/v2g90supp.ipv2");
  AggregateKeyVal keyval(topkeyval,lib1keyval,lib2keyval);

  // construct prefixes for each atom: :basis:<atom>:<basisname>:<shell#>
  // and read in the shell
  nbasis_ = 0;
  Pix i;
  int ishell = 0;
  for (i=molecule.first(); i != 0; molecule.next(i)) {
      AtomicCenter ac = molecule(i);
      ChemicalElement e = ac.element();
      const char* name = e.name();
      sprintf(keyword,":basis:%s:%s",
              molecule(i).element().name(),basisname);
      int count = keyval.count(keyword);
      nshell_ += count;
      recursively_get_shell(ishell,keyval,molecule(i).element().name(),
			    basisname,havepure,pure,0);
    }
  nshell_ = ishell;
  shell = new GaussianShell*[nshell_];
  shell_to_centerpix = new Pix[nshell_];
  shell_to_function_ = new int[nshell_];
  ishell = 0;
  int ifunction = 0;
  for (i=molecule.first(); i != 0; molecule.next(i)) {
      centerpix_to_r[i] = new double[3];
      centerpix_to_r[i][0] = molecule(i)[0];
      centerpix_to_r[i][1] = molecule(i)[1];
      centerpix_to_r[i][2] = molecule(i)[2];
      centerpix_to_shellnum[i] = ishell;
      recursively_get_shell(ishell,keyval,molecule(i).element().name(),
			    basisname,havepure,pure,1);
      centerpix_to_nshell[i] = ishell - centerpix_to_shellnum[i];
      for (int j = centerpix_to_shellnum[i]; j<ishell; j++) {
	  shell_to_centerpix[j] = i;
	}
     }

  // compute nbasis_ and shell_to_function_[]
  nbasis_ = 0;
  for (ishell=0; ishell<nshell_; ishell++) {
      shell_to_function_[ishell] = nbasis_;
      nbasis_ += shell[ishell]->nfunction();
    }

  delete[] basisname;
}

void GaussianBasisSet::
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
      PrefixKeyVal prefixkeyval(prefix,keyval,j);
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
  delete name_;
  delete shell_to_centerpix;
  delete shell_to_function_;

  int ii;
  for (ii=0; ii<nshell_; ii++) {
      delete shell[ii];
    }
  delete[] shell;

  Pix i;
  for (i=centerpix_to_r.first(); i!=0; centerpix_to_r.next(i)) {
      delete[] centerpix_to_r.contents(i);
    }
}

int GaussianBasisSet::max_nfunction_in_shell() const
{
  int i;
  int max = 0;
  for (i=0; i<nshell_; i++) {
      if (max < shell[i]->nfunction()) max = shell[i]->nfunction();
    }
  return max;
}

int GaussianBasisSet::shell_to_function(Pix shellpix) const
{
  return shell_to_function_[((int)shellpix)-1];
}

int GaussianBasisSet::ncenter() const
{
  return centerpix_to_shellnum.length();
}

// The first_center and next_center members use the keys of the
// centerpix_to_shellnum PixMap to iterator over centerpixes.
Pix GaussianBasisSet::first_center() const
{
  return centerpix_to_shellnum.key(centerpix_to_shellnum.first());
}

void GaussianBasisSet::next_center(Pix& centerpix) const
{
  Pix pix = centerpix_to_shellnum.seek(centerpix);
  centerpix_to_shellnum.next(pix);
  if (pix != 0) centerpix = centerpix_to_shellnum.key(pix);
  else centerpix = 0;
}

int GaussianBasisSet::nshell_on_center(Pix centerpix) const
{
  return centerpix_to_nshell[centerpix];
}

Pix GaussianBasisSet::first_shell_on_center(Pix centerpix) const
{
  if (centerpix_to_nshell[centerpix]) {
      return (Pix) (centerpix_to_shellnum[centerpix] + 1);
    }
  else {
      return 0;
    }
}

void GaussianBasisSet::next_shell_on_center(Pix centerpix,Pix& shellpix) const
{
  if ((int)shellpix
      < centerpix_to_shellnum[centerpix]+centerpix_to_nshell[centerpix]) {
      ((int)shellpix)++;
    }
  else {
      shellpix = 0;
    }
}

Pix GaussianBasisSet::first() const
{
  if (nshell_) return (Pix)1;
  else return 0;
}

void GaussianBasisSet::next(Pix&shellpix) const
{
  if ((int)shellpix < nshell_) ((int)shellpix)++;
  else shellpix = 0;
}

const GaussianShell& GaussianBasisSet::operator[](Pix shellpix) const
{
  return *(shell[((int)shellpix)-1]);
}

const double& GaussianBasisSet::r_center(Pix centerpix,int xyz) const
{
  return centerpix_to_r[centerpix][xyz];
}

const double& GaussianBasisSet::r_shell(Pix shellpix,int xyz) const
{
  return centerpix_to_r[shell_to_centerpix[((int)shellpix)-1]][xyz];
}

GaussianShell& GaussianBasisSet::operator[](Pix shellpix)
{
  return *(shell[((int)shellpix)-1]);
}

double& GaussianBasisSet::r_center(Pix centerpix,int xyz)
{
  return centerpix_to_r[centerpix][xyz];
}

double& GaussianBasisSet::r_shell(Pix shellpix,int xyz)
{
  return centerpix_to_r[shell_to_centerpix[((int)shellpix)-1]][xyz];
}

centers_t* GaussianBasisSet::convert_to_centers_t(const Molecule*mol) const
{
  centers_t* c;
  c = (centers_t*)malloc(sizeof(centers_t));
  c->n = ncenter();
  c->center = (center_t*) malloc(sizeof(center_t)*ncenter());
  c->center_num = 0;
  c->shell_num = 0;
  c->func_num = 0;
  int icenter = 0;
  for (Pix center = first_center(); center; next_center(center)) {
      if (mol) {
          const Molecule& molr = *mol;
          c->center[icenter].atom = strdup(molr(center).element().symbol());
          c->center[icenter].charge = molr(center).element().charge();
        }
      else {
          c->center[icenter].atom = strdup("dummy");
          c->center[icenter].charge = 0.0;
        }
      c->center[icenter].r = (double*)malloc(sizeof(double)*3);
      for (int xyz=0; xyz<3; xyz++) c->center[icenter].r[xyz]=r_center(center,xyz);
      c->center[icenter].basis.n = nshell_on_center(center);
      c->center[icenter].basis.name = strdup(name());
      c->center[icenter].basis.shell =
	(shell_t*)malloc(sizeof(shell_t)*nshell_on_center(center));
      shell_t*shell = c->center[icenter].basis.shell;
      int ishell = 0;
      for (Pix pshell = first_shell_on_center(center);
	   pshell;
	   next_shell_on_center(center,pshell)) {
	  int nprim = shell[ishell].nprim = operator[](pshell).nprimitive();
          int ncon = shell[ishell].ncon  = operator[](pshell).ncontraction();
          int nfunc = shell[ishell].nfunc = operator[](pshell).nfunction();
          shell[ishell].exp = (double*)malloc(sizeof(double)*nprim);
          shell[ishell].type = (shell_type_t*)malloc(sizeof(shell_type_t)*ncon);
	  shell[ishell].coef = (double**)malloc(sizeof(double*)*ncon);
          shell[ishell].norm = 0;
	  int i;
	  for (i=0; i<ncon; i++) {
	      shell[ishell].coef[i] = (double*)malloc(sizeof(double)*nprim);
	      shell[ishell].type[i].am = operator[](pshell).am(i);
	      shell[ishell].type[i].puream = operator[](pshell).is_pure(i);
	    }
	  for (i=0; i<nprim; i++) {
	      shell[ishell].exp[i] = operator[](pshell).exponent(i);
	      int j;
	      for (j=0; j<ncon; j++) {
		  shell[ishell].coef[j][i]
		    = operator[](pshell).coefficient_norm(j,i);
		}
	    }
	  ishell++;
	}
      icenter++;
    }
  return c;
}

void GaussianBasisSet::print(FILE*fp) const
{
  fprintf(fp,"GaussianBasisSet: nshell = %d, nbasis = %d\n",nshell_,nbasis_);

  // Loop over centers
  int icenter = 0;
  int ishell = 0;
  for (Pix center=first_center(); center != 0; next_center(center)) {
      fprintf(fp,"center %d: %12.8f %12.8f %12.8f, nshell = %d, shellnum = %d\n",
	      icenter,
	      centerpix_to_r[center][0],
	      centerpix_to_r[center][1],
	      centerpix_to_r[center][2],
	      centerpix_to_nshell[center],
	      centerpix_to_shellnum[center]);
      for (Pix shell=first_shell_on_center(center);
	   shell != 0;
	   next_shell_on_center(center,shell)) {
	  fprintf(fp,"Shell %d: functionnum = %d, centerpix = 0x%x\n",
		  ishell,shell_to_function_[ishell],shell_to_centerpix[ishell]);
	  operator[](shell).print(fp);
          ishell++;
	}
      icenter++;
    }

}
