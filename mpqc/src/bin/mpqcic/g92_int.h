
#ifndef _mpqcic_g92_int_h
#define _mpqcic_g92_int_h

typedef struct g92_calc_type {
  int number;
  int basis_req;
  char *name;
  char *command;
  char *parse_string;
} g92_calc_t;

// don't forget units=au, Freq or Force
// Initialize the list of runtypes and their corresponding names
static int n_g92_calc_types=7;
static g92_calc_t g92_calc[] = {
  {1, 1, "SCF", "SCF=DIRECT RHF", "\\HF="},
  {2, 1, "MP2", "SCF=DIRECT MP2=fulldirect","\\MP2="},
  {3, 0, "AM1", "MNDO", "\\HF="},
  {4, 0, "PM3", "PM3", "\\HF="},
  {5, 0, "MNDO", "MNDO", "\\HF="},
  {6, 1, "UHF", "SCF=DIRECT UHF", "\\HF="},
  {7, 1, "ROHF", "SCF=DIRECT ROHF", "\\HF="}
};

//////////////////////////////////////////////////////////////////////////////

int run_g92(char *name_in, const RefKeyVal &, RefMolecule&,
            double &energy, RefSCVector& gradient);

int run_g92_calc(char *prefix, int runtype, char *basis, int memory,
                 int chk_guess, char *scratch, char *g92_dir, RefMolecule &,
                 int charge, int multiplicity);

int parse_g92(char *prefix, char *parse_string, int natoms, double & energy,
              RefSCVector &gradient);

/////////////////////////////////////////////////////////////////////////////

int g92_freq_driver(char *name_in, RefMolecule&, const RefKeyVal&,
                    char *method, double &energy,
                    RefSCVector& gradient, RefSCVector &frequencies,
                    RefSCVector &normalmodes, int &nmodes, int &nimag);

int run_g92_freq(char *prefix, int runtype, char *basis, int memory,
                 int chk_guess,char *scratch, char *g92_dir, RefMolecule &mole,
                 int charge, int multiplicity);

int parse_g92_freq(char *prefix,char *parse_string, RefMolecule mole,
                   double &energy, RefSCVector& gradient,
                   RefSCVector& frequencies, RefSCVector &normalmodes,
                   int &nmodes, int &nimag);

#endif
