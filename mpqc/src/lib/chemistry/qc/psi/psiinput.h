/*
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
#pragma interface
#endif

#ifndef _CHEMISTRY_QC_PSI_PSI_INPUT_H
#define _CHEMISTRY_QC_PSI_PSI_INPUT_H

#include<chemistry/molecule/molecule.h>
#include<chemistry/qc/basis/basis.h>
extern "C" {
   #include <stdio.h>
}

class PSI_Input {

   private:
      int indentation;
      int memory; // the memory in megabytes
      char * opentype;
      int nirrep;
      int * docc;
      int * socc;
      int * frozen_docc;
      int * frozen_uocc;
      int ex_lvl;
      char * label;
      char * name;
      int nunit;
      char **unit;
      int *nvolume;
      char ***volumes;
      int _test;

   protected:
      RefMolecule _mol;
      RefGaussianBasisSet _gbs;
      FILE *fp;
      void begin_section(const char * s);
      void end_section();
      void write_indent();
      int write_keyword(const char *, const char *);
      int write_keyword(const char *, int);
      int write_keyword(const char *, double);
      int write_keyword(const char *, int, int *);
      int write_keyword(const char *, int, double *);
      int write_geom();
      void write_string(const char *);
//      int write_basis(RefGaussianBasisSet&);
      int write_basis(void);
      int write_defaults(const char *, const char *);
      void write_input();
      int write_key_wq(const char *, const char *);

   public:
      PSI_Input(const RefKeyVal&);
      PSI_Input();
      virtual ~PSI_Input();
      void print(ostream&);
      virtual void write_input_file(const char *,const char *,
               const int convergence = 0, const char *s = "input.dat");
      int test() { return _test; }
};


class PSI_Input_SCF : public PSI_Input
{

   public:
      PSI_Input_SCF(const RefKeyVal&keyval) : PSI_Input(keyval) {}
      PSI_Input_SCF() : PSI_Input() {}
      ~PSI_Input_SCF() {}
      virtual void write_input_file(const char *,const char *,
               const int convergence = 0, const char *s = "input.dat");

};


class PSI_Input_CI : public PSI_Input
{

   public:
      PSI_Input_CI(const RefKeyVal&keyval) : PSI_Input(keyval) {}
      PSI_Input_CI() : PSI_Input() {}
      ~PSI_Input_CI() {}
      virtual void write_input_file(const char *,const char *,
               const int convergence = 0, const char *s = "input.dat");

};

class PSI_Input_CC : public PSI_Input
{

   public:
      PSI_Input_CC(const RefKeyVal&keyval) : PSI_Input(keyval) {}
      PSI_Input_CC() : PSI_Input() {}
      ~PSI_Input_CC() {}
      virtual void write_input_file(const char *,const char *,
               const int convergence = 0, const char *s = "input.dat");
};

#endif
