/*
** Yet another read11 routine
** For use in FILE11 class
**
** David Sherrill, Justin Fermann  -  August 1994
**
*/

/* include's */
#include <stdio.h>

/* define's */
#define MAX_LN_11 133

/*
** READ_FILE11(): This function reads 'file11.dat' and prints out again
**      some of the info to an output file.
**
** Arguments: 
**      natom   = number of atoms
**      readto  = which entry of file11 to read until
**      label   = character string to hold file11 label field
**      theory  = level of theory (in label line)
**      dertype = derivative type (in label line)
**      energy  = pointer to double to hold energy
**      X, Y, Z = pointers to arrays of cartesian coordinates (assume bohr)
**      AN      = pointer to atomic number array 
**      grad    = matrix (3 x natom) to hold gradients
**      ngrad   = pointer to hold number of entries in file11.dat
**
** Returns: success (1) or failure (0) 
**
*/
int read_file11(int natom, int readto,
      char *label, char *theory, char *dertype, 
      double *energy, double *X, double *Y, double *Z, int *AN, 
      double **grad, int *ngrad) 
{
   int i,j ;                /* loop variables */
   FILE *fpi ;              /* pointer for input file file11.dat */
   double g1, g2, g3;       /* junk variables for gradients */
   char line_in[MAX_LN_11]; /* input buffer */
   int count = 0 ;          /* which entry of file11 are we on? */
   int entrynum;
   double tmp;

   /* open the file */
   fpi = fopen("file11.dat", "r");
   if (!fpi) {
      fprintf(stderr, "(read_file11): Can't open file11.dat\n");
      return(0) ;
      }

   if (readto != 0) entrynum = readto;
   else entrynum = 999;

   while (fgets(line_in, MAX_LN_11,fpi)!=NULL && count!=entrynum) {

      sscanf(line_in, "%s %s %s", label, theory, dertype);

      /* now read the number of atoms and the energy */
      fgets(line_in, MAX_LN_11, fpi);
      if (sscanf(line_in, "%d %lf", &j, energy) != 2) {
         fprintf(stderr, "(read_file11): Trouble reading natoms and energy\n");
         return(0);
         }
      if (j != natom) {
         fprintf(stderr, "(read_file11): Read %d atoms, expecting %d\n",
            j, natom);
         return(0);
         }

      /* read Cartesians */
      for (i=0; i<natom; i++) {
         if (fscanf(fpi, " %lf %lf %lf %lf ", &tmp, X+i, Y+i, Z+i) != 4){
            fprintf(stderr, "(read_file11): Trouble reading Cartesians\n");
	    }
         AN[i] = (int)tmp;
         }

      /* read Gradients */
      for (i=0; i<natom; i++) {
         if (fscanf(fpi, "%lf %lf %lf", &g1, &g2, &g3) != 3)
            fprintf(stderr, "(read_file11): Trouble reading gradients\n");

         else {
            grad[0][i] = g1;
            grad[1][i] = g2;
            grad[2][i] = g3;
            }
         }

      fgets(line_in, MAX_LN_11, fpi) ;   /* go to end of line */
      }

   fclose(fpi) ;

   *ngrad = count;
   if (readto == 0) return(1);
   else if (count == entrynum) return(1);
   else return(0);
}



/*
** PRINT_FILE11(): Function prints out the information from file11.dat
** obtained from the read_file11() function.
**
** Arguments: 
**      label   = character string to hold file11 label field
**      natom   = number of atoms
**      energy  = energy from file11
**      X, Y, Z = arrays of cartesian coordinates (assume bohr)
**      AN      = atomic number array (allocated here)
**      fpo     = file pointer for output
*/
void print_file11(char *label, int natom, double energy, double *X, 
      double *Y, double *Z, int *AN, FILE *fpo) 
{
int i ;

   fprintf(fpo, "DATA FROM FILE11.DAT\n") ;
   fprintf(fpo, "Label :\n%s\n", label) ;
   fprintf(fpo, "Number of atoms = %d\n", natom) ;
   fprintf(fpo, "Energy = %.10f\n", energy) ;
   fprintf(fpo, "Cartesian coordinates (bohr) :\n") ;
   for (i=0; i<natom; i++) {
      fprintf(fpo, "     %4d    %12.7f    %12.7f    %12.7f\n",
            AN[i], X[i], Y[i], Z[i]) ;
      }
   fprintf(fpo, "\n") ;
   fflush(fpo) ;
}


