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
#include <chemistry/qc/psi/psiinput.h>
#include <chemistry/qc/basis/basis.h>
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
  printf("here it comes\n");
  _gbs->print();
  // _gbs = keyval->describedclassvalue("basis");
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
  n = keyval->count("docc");
  if (keyval->error() != KeyVal::OK || n != nirrep) {
      cerr << "change size of docc array" << endl;
      abort();
    }
  for (i=0; i<nirrep; i++) 
      docc[i] = keyval->intvalue("docc",i);
  
  if (!strcmp(opentype, "NONE")) {
    n = keyval->count("socc");
    if (keyval->error() != KeyVal::OK || n != nirrep) {
        cerr << "change size of socc array" << endl;
        abort();
        }
    for (i=0; i<nirrep; i++) {
      socc[i] = keyval->intvalue("socc",i);
      }
    }
  else for (i=0; i<nirrep; i++) 
    socc[i] = 0;


  if (keyval->exists("frozen_docc")) {
    n = keyval->count("frozen_docc");
    if (keyval->error() != KeyVal::OK || n != nirrep) {
        cerr << "change size of frozen_docc array" << endl;
        abort();
      }
    for (i=0; i<nirrep; i++) 
        frozen_docc[i] = keyval->intvalue("frozen_docc",i);
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
   if (fprintf(fp, "%s = %lf\n", s, t) < 0) {
      fprintf(stderr,"(PSI_Input::write_keyword): trouble writing %s = %lf\n",
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
      if (fprintf(fp, "%lf ", t[i]) < 0) errcod=1;
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
int errcod;
char ts[133];

    int *unique ;
    unique = _mol->find_unique_atoms();
    write_string("geometry = (\n");
    for (int i=0; i < _mol->num_unique_atoms(); i++) {
        sprintf(ts, "  (%f %f %f)\n", _mol->r(unique[i],0),
           _mol->r(unique[i],1), _mol->r(unique[i],2));
        write_string(ts);
        } 
    write_string("    )\n");
    delete[] unique;
    if (errcod) return(0);
    else return(1);
}


int
PSI_Input::write_basis(void)
{
  int i, j, k, l, am;
  char ts[133];

  begin_section("basis");
  RefAtomInfo atominfo = _mol->atominfo();
  int *unique ;
  unique = _mol->find_unique_atoms();
  for (i=0; i<_mol->num_unique_atoms(); i++) {
    sprintf(ts, "%s:SCdefined = (\n", 
	    atominfo->name(_mol->Z(unique[i])));
    write_string(ts);
    for(am=0; am<6; am++){
      for (j = 0; j < _gbs->nshell_on_center(unique[i]); j++) {
	for (l=(*_gbs)(unique[i],j).ncontraction()-1; l>-1;l--) {
          if((*_gbs)(unique[i],j).am(l)==am){
	    sprintf(ts, "  (%c\n", (*_gbs)(unique[i],j).amchar(l));
	    write_string(ts);
	    for(k=0; k<(*_gbs)(unique[i],j).nprimitive(); k++){
	      sprintf(ts, "    (%f   %f)\n", (*_gbs)(unique[i],j).exponent(k),
		      (*_gbs)(unique[i],j).coefficient_norm(l,k));
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
  delete[] unique;
          
}


int
PSI_Input::write_defaults(const char *dertype, const char *wavefn)
{
   char ts[133];
   int i;
   double *x_vec, *y_vec, *z_vec;

   begin_section("default");
   write_keyword("output", "terminal");
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
   
   write_keyword("docc", nirrep, docc);
   write_keyword("socc", nirrep, socc);
   write_keyword("frozen_docc", nirrep, frozen_docc);
   write_keyword("frozen_uocc", nirrep, frozen_uocc);

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
  int *unique ;
  char t1[133];
  char t2[133];

  begin_section("input");
  RefAtomInfo atominfo = _mol->atominfo();
  unique = _mol->find_unique_atoms();
  sprintf(t1, "atoms = (");
  for (i=0; i < _mol->num_unique_atoms(); i++) {
    sprintf(t2, "%s ", atominfo->symbol(_mol->Z(unique[i])));
    strcat(t1, t2);
    }
  strcat(t1, ")\n");
  write_string(t1);
  delete[] unique;

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


void
PSI_Input_SCF::write_input_file(const char *dertype, const char *wavefn, 
    const int convergence, const char *fname )
{
  fp = fopen(fname, "w");
  if (fp == NULL) {
    cerr << "(PSI_Input_SCF::write_input_file): Can't open " << fname << endl;
    abort();
    }
  write_defaults(dertype, wavefn);
  write_input();
  write_basis();
  if(convergence != 0){
    begin_section("scf");
    write_keyword("convergence", convergence);
    end_section();
    }

  fclose(fp);
}

void
PSI_Input_CI::write_input_file(const char *dertype, const char *wavefn,
    const int convergence, const char *fname )
{
  fp = fopen(fname, "w");
  if (fp == NULL) {
    cerr << "(PSI_Input_CI::write_input_file): Can't open "
         << fname << endl;
    abort();
    }
  write_defaults(dertype, wavefn);
  write_input();
  write_basis();
  if(convergence != 0){
    begin_section("gugaci");
    write_keyword("convergence", convergence);
    end_section();
    }

  fclose(fp);
}

void
PSI_Input_CC::write_input_file(const char *dertype, const char *wavefn,
    const int convergence, const char *fname )
{
  fp = fopen(fname, "w");
  if (fp == NULL) {
    cerr << "(PSI_Input_CC::write_input_file): Can't open " << fname << endl;
    abort();
    }
  write_defaults(dertype, wavefn);
  write_input();
  write_basis();
  if(convergence != 0){
    begin_section("cceg");
    write_keyword("convergence", convergence);
    end_section();
    }

  fclose(fp);
}

