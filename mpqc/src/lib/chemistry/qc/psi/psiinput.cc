/*
**
** PSI Input Class
**
** This helper class will set up input decks for the PSI suite of
** ab initio quantum chemistry programs.
**
** David Sherrill & Justin Fermann
** Center for Computational Quantum Chemistry, University of Georgia
**
*/

#ifdef __GNUG__
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <math/symmetry/corrtab.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/psi/psiinput.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
extern "C" {
#include <stdio.h>
#include <math.h>
}


PSI_Input::PSI_Input()
{
  printf("in default constructor\n");
}
PSI_Input::PSI_Input(const RefKeyVal&keyval)
{
  int i,n;
  int tmp;
  char *ts;

  indentation = 0;
  _mol = keyval->describedclassvalue("molecule");
  _gbs = keyval->describedclassvalue("psibasis");
  _gbs->print();
  // _gbs = keyval->describedclassvalue("basis");
  _origpg = new PointGroup(*_mol->point_group().pointer());
  nirrep = _mol->point_group()->char_table().nirrep();
  docc = new int[nirrep];
  socc = new int[nirrep];
  frozen_docc = new int[nirrep];
  frozen_uocc = new int[nirrep];
  RefKeyVal psifiles = new PrefixKeyVal("psifiles",keyval);

  memory = keyval->intvalue("memory");
  if (keyval->error() != KeyVal::OK) {
      memory = 10;
    }
  name =  psifiles->pcharvalue("name");
  nunit = psifiles->count("unit");
  nvolume = new int[nunit];
  unit = new char*[nunit];
  volumes = new char**[nunit];
  for(n=0; n<nunit; n++){
    ts = psifiles->pcharvalue("unit",n);
    if(strcmp(ts, "default")){
      unit[n] = new char[strlen(ts)+5];
      strcpy(unit[n], "file");
      strcat(unit[n], ts);
      delete[] ts;
      }
    else
      unit[n] = ts;
    nvolume[n] = psifiles->intvalue("nvolume", n);
    tmp = psifiles->count("volumes", n);
    if (tmp != nvolume[n])
        cerr << "bad nvolumes" << endl;
    volumes[n] = new char*[tmp];
    for(int j=0; j<tmp; j++)
      volumes[n][j] = psifiles->pcharvalue("volumes", n, j);
    }

  opentype = keyval->pcharvalue("opentype");
  label = keyval->pcharvalue("label");
  _test = keyval->booleanvalue("test");

  RefOneBodyWavefunction obwfn = keyval->describedclassvalue("obwfn");
  const double epsilon = 0.001;
  RefPetiteList pl;
  if (obwfn.nonnull()) pl = obwfn->integral()->petite_list(obwfn->basis());

  n = keyval->count("docc");
  if (keyval->error() != KeyVal::OK || n != nirrep) {
      if (obwfn.nonnull()) {
          for (i=0; i<nirrep; i++) {
              docc[i] = 0;
              for (int j=0; j<pl->nfunction(i); j++) {
                  if (obwfn->occupation(i,j) > 2.0-epsilon) docc[i]++;
                }
            }
        }
      else {
          cerr << "change size of docc array or give obwfn" << endl;
          abort();
        }
    }
  else {
      for (i=0; i<nirrep; i++) 
          docc[i] = keyval->intvalue("docc",i);
    }
  
  if (strcmp(opentype, "NONE")) {
    n = keyval->count("socc");
    if (keyval->error() != KeyVal::OK || n != nirrep) {
      if (obwfn.nonnull()) {
          for (i=0; i<nirrep; i++) {
              socc[i] = 0;
              for (int j=0; j<pl->nfunction(i); j++) {
                  if (obwfn->occupation(i,j) > 1.0-epsilon
                      && obwfn->occupation(i,j) < 1.0+epsilon) socc[i]++;
                }
            }
        }
      else {
          cerr << "change size of socc array or give obwfn" << endl;
          abort();
        }
      }
    else {
        for (i=0; i<nirrep; i++) {
            socc[i] = keyval->intvalue("socc",i);
          }
      }
    }
  else {
      for (i=0; i<nirrep; i++) 
          socc[i] = 0;
    }


  if (keyval->exists("frozen_docc")) {
    n = keyval->count("frozen_docc");
    if (keyval->error() != KeyVal::OK || n != nirrep) {
        char *tmp = keyval->pcharvalue("frozen_docc");
        if (tmp) {
            if (!strcmp(tmp,"auto") && obwfn.nonnull()) {
                int nfzc = _mol->n_core_electrons()/2;
                cout << node0 << indent
                     << "PSI: auto-freezing "<<nfzc<<" core orbitals" << endl;
                RefDiagSCMatrix eigvals = obwfn->eigenvalues().copy();
                for (i=0; i<nirrep; i++) frozen_docc[i] = 0;
                while (nfzc) {
                    double smallest = 0.0;
                    int smallesti=0;
                    for (i=0; i<obwfn->basis()->nbasis(); i++) {
                        if (smallest > eigvals(i)) {
                            smallest = eigvals(i);
                            smallesti = i;
                          }
                      }
                    eigvals(smallesti) = 0.0;
                    int orbnum = 0;
                    for (i=0; i<nirrep; i++) {
                        orbnum += pl->nfunction(i);
                        if (smallesti < orbnum) {
                            frozen_docc[i]++;
                            break;
                          }
                      }
                    nfzc--;
                  }
              }
            else {
                cerr << "bad value for frozen_docc or missing obwfn" << endl;
                abort();
              }
            delete[] tmp;
          }
        else {
            cerr << "change size of frozen_docc array" << endl;
            abort();
          }
      }
    else {
        for (i=0; i<nirrep; i++) 
            frozen_docc[i] = keyval->intvalue("frozen_docc",i);
      }
    }
  else for (i=0; i<nirrep; i++) 
    frozen_docc[i] = 0;

  if (keyval->exists("frozen_uocc")) {
    n = keyval->count("frozen_uocc");
    if (keyval->error() != KeyVal::OK || n != nirrep) {
        cerr << "change size of frozen_uocc array" << endl;
        abort();
      }
    for (i=0; i<nirrep; i++) 
        frozen_uocc[i] = keyval->intvalue("frozen_uocc",i);
    }
  else for (i=0; i<nirrep; i++) 
    frozen_uocc[i] = 0;

  if (keyval->exists("ex_lvl")) ex_lvl = keyval->intvalue("ex_lvl");
  else ex_lvl = 0;

  cout << node0 << indent << "docc = [";
  for (i=0; i<nirrep; i++) cout << node0 << " " << docc[i];
  cout << node0 << " ]" << endl;

  cout << node0 << indent << "socc = [";
  for (i=0; i<nirrep; i++) cout << node0 << " " << socc[i];
  cout << node0 << " ]" << endl;

  cout << node0 << indent << "frozen_docc = [";
  for (i=0; i<nirrep; i++) cout << node0 << " " << frozen_docc[i];
  cout << node0 << " ]" << endl;

  cout << node0 << indent << "frozen_uocc = [";
  for (i=0; i<nirrep; i++) cout << node0 << " " << frozen_uocc[i];
  cout << node0 << " ]" << endl;
}


PSI_Input::~PSI_Input()
{
  delete[] opentype;
  delete[] docc;
  delete[] socc;
  delete[] frozen_docc;
  delete[] frozen_uocc;
  delete[] label;
  delete[] name;
  for (int i=0; i<nunit; i++) {
    delete[] unit[i];
    for (int j=0; j<nvolume[i]; j++) {
      delete[] volumes[i][j];
      }
    delete[] volumes[i];
    }
  delete[] volumes;
  delete[] unit;
  delete[] nvolume;
}

void
PSI_Input::open(const char *fname)
{
  fp = fopen(fname, "w");
  if (fp == NULL) {
    cerr << "(PSI_Input_CI::write_input_file): Can't open "
         << fname << endl;
    abort();
    }
}

void
PSI_Input::close()
{
  fclose(fp);
  fp = 0;
}

void
PSI_Input::begin_section(const char * s)
{
   write_indent();
   indentation += 2;
   fprintf(fp, "%s: (\n", s);
}


void
PSI_Input::end_section(void)
{
   indentation -= 2;
   write_indent();
   fprintf(fp, ")\n");
}


int
PSI_Input::write_keyword(const char *s, const char *t)
{
   write_indent();
   if (fprintf(fp, "%s = %s\n", s, t) < 0) {
      fprintf(stderr,"(PSI_Input::write_keyword): trouble writing %s = %s\n",
         s, t);
      return(0);
      }
   else return(1);
}

int
PSI_Input::write_keyword(const char *s, int t)
{
   write_indent();
   if (fprintf(fp, "%s = %d\n", s, t) < 0) {
      fprintf(stderr,"(PSI_Input::write_keyword): trouble writing %s = %d\n",
         s, t);
      return(0);
      }
   else return(1);
}

int
PSI_Input::write_keyword(const char *s, double t)
{
   write_indent();
   if (fprintf(fp, "%s = %f\n", s, t) < 0) {
      fprintf(stderr,"(PSI_Input::write_keyword): trouble writing %s = %f\n",
         s, t);
      return(0);
      }
   else return(1);
}

int
PSI_Input::write_keyword(const char *s, int num, int * t)
{
int errcod=0;

   write_indent();
   fprintf(fp, "%s = (", s);
   for (int i=0; i<num; i++) {
      if (fprintf(fp, "%d ", t[i]) < 0) errcod=1;
      }
   fprintf(fp, ")\n");
   if (errcod) {
      fprintf(stderr,"(PSI_Input::write_keyword): trouble writing array ");
      fprintf(stderr,"for %s\n",s);
      return(0);
      }
   else return(1);
}

int
PSI_Input::write_keyword(const char *s, int num, double * t)
{
int errcod=0;

   write_indent(); 
   fprintf(fp, "%s = (", s);
   for (int i=0; i<num; i++) {
      if (fprintf(fp, "%f ", t[i]) < 0) errcod=1;
      }
   fprintf(fp, ")\n");
   if (errcod) {
      fprintf(stderr,"(PSI_Input::write_keyword): trouble writing array ");  
      fprintf(stderr,"for %s\n",s);
      return(0);
      }
   else return(1);
}

int
PSI_Input::write_geom()
{
int i;
char ts[133];

    write_string("% original MPQC geometry:\n");
    for (i=0; i < _mol->natom(); i++) {
        sprintf(ts, "%% %3d % 14.12f % 14.12f % 14.12f\n",
                _mol->Z(i),
                _mol->r(i,0),
                _mol->r(i,1),
                _mol->r(i,2));
        write_string(ts);
      }
    write_string("geometry = (\n");
    for (i=0; i < _mol->nunique(); i++) {
        sprintf(ts, "  (% 14.12f % 14.12f % 14.12f)\n",
                _mol->r(_mol->unique(i),0),
                _mol->r(_mol->unique(i),1),
                _mol->r(_mol->unique(i),2));
        write_string(ts);
        } 
    write_string("    )\n");
    return 0;
}


int
PSI_Input::write_basis(void)
{
  int i, j, k, l, am;
  char ts[133];

  begin_section("basis");
  RefAtomInfo atominfo = _mol->atominfo();
  for (i=0; i<_mol->nunique(); i++) {
    int uniquei = _mol->unique(i);
    sprintf(ts, "%s:SCdefined = (\n", 
	    atominfo->name(_mol->Z(uniquei)));
    write_string(ts);
    for(am=0; am<6; am++){
      for (j = 0; j < _gbs->nshell_on_center(uniquei); j++) {
	for (l=(*_gbs)(uniquei,j).ncontraction()-1; l>-1;l--) {
          if((*_gbs)(uniquei,j).am(l)==am){
            const char *purestring = "";
            if (am==2&&(*_gbs)(uniquei,j).is_pure(l)) purestring = "5";
            if (am==3&&(*_gbs)(uniquei,j).is_pure(l)) purestring = "7";
            if (am==4&&(*_gbs)(uniquei,j).is_pure(l)) purestring = "9";
            if (am==5&&(*_gbs)(uniquei,j).is_pure(l)) purestring = "11";
	    sprintf(ts, "  (%c%s\n", (*_gbs)(uniquei,j).amchar(l),
                    purestring);
	    write_string(ts);
	    for(k=0; k<(*_gbs)(uniquei,j).nprimitive(); k++){
	      sprintf(ts, "    (%22.16f   % 18.16f)\n",
                      (*_gbs)(uniquei,j).exponent(k),
		      (*_gbs)(uniquei,j).coefficient_norm(l,k));
	      write_string(ts);
	      }
	    write_string("  )\n");
	    }
          }
	}
      }
    write_string("  )\n");
    }
  end_section();
  return 0;
}

void
PSI_Input::write_orbvec(const CorrelationTable &corrtab,
                        const char *orbvec_name,
                        const int *orbvec)
{
  int *orbvecnew = new int[corrtab.subn()];
  memset(orbvecnew,0,sizeof(int)*corrtab.subn());

  for (int i=0; i<corrtab.n(); i++) {
      for (int j=0; j<corrtab.ngamma(i); j++) {
          int gam = corrtab.gamma(i,j);
          orbvecnew[gam] += (corrtab.subdegen(gam)*orbvec[i])/corrtab.degen(i);
        }
    }

  write_keyword(orbvec_name, corrtab.subn(), orbvecnew);

  delete[] orbvecnew;
}

int
PSI_Input::write_defaults(const char *dertype, const char *wavefn)
{
   char ts[133];
   int i;
   double *x_vec, *y_vec, *z_vec;

   // PSI doesn't do d2 right. The X and Y axes are wrong and
   // sometimes intsth fails.  So d2 gets lowered to c2.
   if (!strcmp(_mol->point_group()->symbol(),"d2")) {
       cout << node0 << indent
            << "DOING D2 calc in C2 because of bugs in PSI"
            << endl;
       RefPointGroup newgrp(new PointGroup("c2",
                                           _mol->point_group()->symm_frame(),
                                           _mol->point_group()->origin()));
       _mol->set_point_group(newgrp);
     }

   int dograd = !strcmp(dertype,"FIRST");

   begin_section("default");
   // Gradient calcs write the energy to file11 so the output
   // can go to the terminal.  Otherwise, the output must go
   // to output.dat from which the energy will be read.
   if (dograd) write_keyword("output", "terminal");
   write_keyword("wfn", wavefn);
   write_keyword("dertype", dertype);
   write_keyword("opentype", opentype);
   write_key_wq("label",label);
   sprintf(ts, "memory = (%d MB)\n", memory/1000000);
   write_string(ts);

   write_keyword("symmetry",_mol->point_group()->symbol());

   x_vec = new double[3];
   y_vec = new double[3];
   z_vec = new double[3];
   for (i=0; i<3; i++) {
     x_vec[i] = _mol->point_group()->symm_frame()[i][0];
     y_vec[i] = _mol->point_group()->symm_frame()[i][1];
     z_vec[i] = _mol->point_group()->symm_frame()[i][2];
     }
 
   write_keyword("x_axis", 3, x_vec);
   write_keyword("y_axis", 3, y_vec);
   write_keyword("z_axis", 3, z_vec);
   delete[] x_vec;
   delete[] y_vec;
   delete[] z_vec;

   // perhaps the symmetry has been lowered
   // make sure that the occupation vectors are still correct
   CorrelationTable corrtab;
   int rc;
   if ((rc = corrtab.initialize_table(_origpg, _mol->point_group()))) {
       cerr << node0
            << "ERROR: couldn't initialize correlation table:" << endl
            << "  " << corrtab.error(rc) << endl;
       abort();
     }
   write_orbvec(corrtab, "docc", docc);
   write_orbvec(corrtab, "socc", socc);
   write_orbvec(corrtab, "frozen_docc", frozen_docc);
   write_orbvec(corrtab, "frozen_uocc", frozen_uocc);

   if (ex_lvl) write_keyword("ex_lvl", ex_lvl);
   begin_section("files");
   for (i=0; i<nunit; i++) {
     begin_section(unit[i]);
     if (strcmp(unit[i], "default")==0) write_key_wq("name",name);
     write_keyword("nvolume", nvolume[i]);
     for (int j=0; j<nvolume[i]; j++) {
        sprintf(ts, "volume%d = \"%s\"\n", j+1, volumes[i][j]); 
        write_string(ts);
        }
     end_section();
     }
   end_section();
   end_section(); 
   return 0;
}


void 
PSI_Input::write_string(const char *s)
{
   write_indent();
   fprintf(fp, "%s", s);
}


void 
PSI_Input::write_indent()
{
   for (int i=0; i<indentation; i++) fprintf(fp, " ");
}

void
PSI_Input::print(ostream& o)
{
  int i;

  o << indent << "opentype = " << opentype << endl;
  o << indent << "label = " << label << endl;
  o << indent << "test = " << _test << endl;
  o << indent << "memory = " << memory << endl;
             
  o << indent << "docc:";
  for (i=0; i<nirrep; i++) {
    o << " " << docc[i];
    }
  o << endl;

  o << indent << "socc:";
  for (i=0; i<nirrep; i++) {
    o << " " << socc[i];
    }
  o << endl;

  o << indent << "frozen_docc:";
  for (i=0; i<nirrep; i++) {
     o << " " << frozen_docc[i];
     }
  o << endl;

  o << indent << "frozen_uocc:";
  for (i=0; i<nirrep; i++) {
    o << " " << frozen_uocc[i];
    }
  o << endl;

}

void
PSI_Input::write_input_file(const char *dertype, const char *wavefn, 
    const int convergence, const char *fname )
{

  fp = fopen(fname, "w");
  if (fp == NULL) {
    cerr << "(PSI_Input::write_input_file): Can't open " << fname << endl;
    abort();
    }
  write_defaults(dertype, wavefn);
  write_input();
  fclose(fp);
}

void
PSI_Input::write_input(void)
{
  int i;
  char t1[133];
  char t2[133];

  begin_section("input");
  RefAtomInfo atominfo = _mol->atominfo();
  sprintf(t1, "atoms = (");
  for (i=0; i < _mol->nunique(); i++) {
    sprintf(t2, "%s ", atominfo->symbol(_mol->Z(_mol->unique(i))));
    strcat(t1, t2);
    }
  strcat(t1, ")\n");
  write_string(t1);

  write_geom();
  write_keyword("basis", "SCdefined");
  end_section();
}

int
PSI_Input::write_key_wq(const char *s, const char *t)
{
   write_indent();
   if (fprintf(fp, "%s = \"%s\"\n", s, t) < 0) {
      cerr << "(PSI_Input::write_key_wq): trouble writing"
           << s << "= \"" << t << "\"" << endl;
      return(0);
      }
   else return(1);
}

