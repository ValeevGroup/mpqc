
#include <stdio.h>
#include <string.h>
#include <tmpl.h>
//#include <math_lib.h>
//#include <int_libv2.h>
#include <chemistry/molecule/chemelem.h>

double get_fortran_double(char*line,int position);
void printbasissp(FILE* fp,int n,double* exp,double* con1,double* con2,char* prec);
void printbasis(FILE* fp,char* type,int n,double* exp,double* con,char* prec);
double scan_fortran_double(FILE*);
double sscan_fortran_double(char*);

main(int argc,char** argv)
{
  int i;
  char shelltype[40];
  int ngauss;
  double scale;
  double exp[100];
  double con1[100];
  double con2[100];
  char form2[50],form3[50];
  char *prec = "% .9e";
  int nconv;
  FILE *input;
  const int MaxLineLength = 200;
  char line[MaxLineLength];
  char atomsymbol[MaxLineLength],*basisname;
  const char *atom;
  char tokens[200];
  int nconv0;
  char strexp[100];
  int usestrexp;

  if (argc != 2) {
    fprintf(stderr,"the basis set name is a required argument\n");
    return 1;
    }

  strcpy(tokens,argv[1]);
  input = fopen(argv[1],"r");
  if (!input) {
    fprintf(stderr,"open failed for basis file %s\n",argv[1]);
    return 1;
    }

  printf("basis: (\n");

  strcpy(form2,prec);
  strcat(form2," ");
  strcat(form2,prec);
  strcat(form2,"\n");

  strcpy(form3,prec);
  strcat(form3," ");
  strcat(form3,prec);
  strcat(form3," ");
  strcat(form3,prec);
  strcat(form3,"\n");

  /* Get the basis name from the filename. */
  strcpy(tokens,argv[1]);
  basisname  = strtok(tokens,".");
  /* fprintf(stderr,"basis = %s\n",basisname); */

  while (fgets(line,MaxLineLength,input)) {
      // ignore comments
      if (line[0] == '!') {
          //fprintf(stderr,"skipping comment\n");
          continue;
        }

      if ((nconv0 = sscanf(line,"%s",atomsymbol)) != 1)  {
          fprintf(stderr,"conversion to atomsymbol failed\n");
          break;
        }

      char * tatomsymbol;

      if (atomsymbol[0] == '-') {
          tatomsymbol = &atomsymbol[1];
        }
      else {
          tatomsymbol = &atomsymbol[0];
        }

      ChemicalElement element(tatomsymbol);
      atom = element.name();

    /* fprintf(stdout," atomsymbol = %s\n",atomsymbol); */
    /* fprintf(stdout," atom = %s\n",atom); */
    fprintf(stdout,"  %s: %s: [\n",atom,basisname);
    fgets(line,MaxLineLength,input);
    sscanf(line,"%s",shelltype);
    while (1) {

      usestrexp = 0;
      if (shelltype[0]<'A' || shelltype[0]>'Z') {
          /* We didn't actually get what we wanted. So we use default
           * shelltype, ngauss and scale info -- or at least what is the
           * apparent default. */
          strcpy(strexp,shelltype);
          strcpy(shelltype,"P");
          ngauss = 1;
          scale = 1.0;
          usestrexp = 1;
        }
      else {
          sscanf(line,"%*s %d %lf",&ngauss,&scale);
        }
      for (i=0; i<strlen(shelltype); i++) shelltype[i] -= 'A'-'a';
      for (i=0; i<ngauss; i++) {
        if (!usestrexp) fgets(line,MaxLineLength,input);
        if (!strcmp(shelltype,"sp")) {
          exp[i] = get_fortran_double(line,1);
          con1[i] = get_fortran_double(line,2);
          con2[i] = get_fortran_double(line,3);
          }
        else {
          if (usestrexp) {
            exp[i] = sscan_fortran_double(strexp);
	    }
	  else {
            exp[i] = get_fortran_double(line,1);
	    }
          con1[i] = get_fortran_double(line,2);
          }
        if (nconv == EOF) {
          fprintf(stderr,"unexpected end of file\n");
          return 1;
          }
        exp[i] *= scale*scale;
        if (!strcmp(shelltype,"sp")) {
          /* fprintf(stderr,form3,exp[i],con1[i],con2[i]); */
          }
        else {
          /* fprintf(stderr,form2,exp[i],con1[i]); */
          }
        }

      /* Print out the basis in Psi v2 format: */
#define USE_GC
#ifdef USE_GC
      if (!strcmp(shelltype,"sp")) {
        printbasissp(stdout,ngauss,exp,con1,con2,prec);
        }
#else
      if (!strcmp(shelltype,"sp")) {
        printbasis(stdout,"p",ngauss,exp,con2,prec);
        printbasis(stdout,"s",ngauss,exp,con1,prec);
        }
#endif
      else {
        printbasis(stdout,shelltype,ngauss,exp,con1,prec);
        }

      if (!fgets(line,MaxLineLength,input)) {
          fprintf(stderr,"couldn't find ending ****\n");
          fprintf(stderr,"got eof\n",line);
          break;
        }
      if (!strncmp(line,"****",4)) {
          break;
        }

      sscanf(line,"%s",shelltype);
      }

    fprintf(stdout,"    ]\n");

    }

  printf("  )\n");
  return 0;
  }

void
printbasissp(FILE* fp,int n,double* exp,double* con1,double* con2,char* prec)
{
  int i;
  char *end;
  char fmt[200];

  sprintf(fmt,"      %s  %s  %s%%s\n",prec,prec,prec);

  fprintf(fp,"    (type: [am = p  am = s] \n");
  fprintf(fp,"            {exp               coef:0        coef:1} = {\n");
  end = "";
  for (i=0; i<n; i++) {
    if (i==n-1) end = "})";
    fprintf(fp,fmt,exp[i],con2[i],con1[i],end);
    }
  }

void
printbasis(FILE* fp,char* type,int n,double* exp,double* con,char* prec)
{
  int i;
  char *end;
  char fmt[200];

  sprintf(fmt,"      %s  %s%%s\n",prec,prec);

  fprintf(fp,"    (type: [am = %s]\n",type);
  fprintf(fp,"            {exp               coef:0} = {\n",type);
  end = "";
  for (i=0; i<n; i++) {
    if (i==n-1) end = "})";
    fprintf(fp,fmt,exp[i],con[i],end);
    }
  }

double
get_fortran_double(char*line,int position)
{
  char buff[100];
  char* fmt;

  if (position == 1) fmt = "%s";
  else if (position == 2) fmt = "%*s %s";
  else if (position == 3) fmt = "%*s %*s %s";
  
  sscanf(line,fmt,buff);
  return sscan_fortran_double(buff);
}

double
scan_fortran_double(FILE* fp)
{
  char buff[100];
  double result;

  fscanf(fp,"%s",buff);

  return sscan_fortran_double(buff);
  }

double
sscan_fortran_double(char* str)
{
  double result;
  char *ch;

  for (ch=str; *ch!='\0'; ch++) {
    if (*ch == 'D') *ch = 'e';
    }

  sscanf(str,"%lf",&result);

  return result;
  }
