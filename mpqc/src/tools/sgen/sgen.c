
/* $Log$
 * Revision 1.3  1994/10/18 23:04:01  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.2  1994/06/08  01:14:15  cljanss
 * Many changes.  These include: newmat7 and nihmatrix -> scmat
 * and mpqcic -> MPSCF and updated optimize stuff.
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1993/04/28  17:29:07  jannsen
 * add -ifndef option
 *
 * Revision 1.1.1.1  1992/03/17  18:11:03  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:11:02  seidl
 * Initial revision
 *
 * Revision 1.3  1992/01/09  12:15:16  seidl
 * add send0 and recv0 generators
 *
 * Revision 1.2  1991/12/23  22:23:48  seidl
 * add gen functions iseq, assign, and zero
 *
 * Revision 1.1  1991/12/20  16:18:23  seidl
 * Initial revision
 *
 * Revision 1.8  91/11/18  17:24:04  cljanss
 * added fwrite and fread
 * 
 * Revision 1.7  91/11/06  14:20:40  cljanss
 * Modified so flex/bison can be used instead of lex/yacc.
 * 
 * Revision 1.6  91/09/28  16:40:39  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.5  91/07/19  19:25:47  cljanss
 * added bcast0
 * 
 * Revision 1.4  1991/07/19  14:42:33  cljanss
 * Changed the way that generation modules are selected.
 *
 * Revision 1.3  1991/06/17  17:25:57  seidl
 * add bread_gen
 *
 * Revision 1.2  1991/06/16  20:33:16  seidl
 * added bwrite_gen
 *
 * Revision 1.1  1991/06/15  21:13:57  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <tmpl.h>
#include "types.h"
#define  _ALLOCATE_GLOBAL_
#include "global.h"

#include "sgen.gbl"
#include "sgen.lcl"

#include "sgen_util.gbl"

int
main(argc,argv)
int argc;
char **argv;
{
  char *tmpBaseName;
  int i;
  char *ch;
  int ifndef;
  char *ifndefname;
  int first = 1;

  progname = argv[0];

  for (i=1; i<argc; i++) {
    ifndef = 0;

    if (!strcmp(argv[i],"-ifndef")) {
      ifndef = 1;
      i++;
      ifndefname = argv[i];
      if (i == argc) {
        fprintf(stderr,"sgen: The -ifndef option requires an argument.\n");
        exit(1);
        }
      i++;
      }
    input = fopen(argv[i],"r");
    if (!input) {
      fprintf(stderr,"%s: couldn't open the input file %s\n",progname,argv[i]);
      exit(1);
      }

    yyin = input;
    init_declarations();

    /* Compute a basename which is used to determine the output file
     * name. */
    tmpBaseName = strdup(argv[i]);
    /* remove the suffix */
    ch = strrchr(tmpBaseName,'.');
    if (ch) *ch = '\0';
    /* remove directory names */
    if (ch = strrchr(tmpBaseName,'/')) {
        strcpy(BaseName,&ch[1]);
      }
    else {
        strcpy(BaseName,tmpBaseName);
      }
    free(tmpBaseName);

    /* This opens the structure output file, since the parsing process
     * can cause some stuff to be written to it. */
    struct_open();

    if (ifndef) {
      fprintf(includeoutput,"#ifndef %s\n",ifndefname);
      fprintf(includeoutput,"#define %s\n",ifndefname);
      }

    /* Parse the sgen file. */
#ifdef FLEX
    if (!first) yyrestart();
#endif /* FLEX */
    yyparse();
    fclose(input);
    input = NULL;

    /* The structure file is always generated. */
    printf("Generating for %s: struct",argv[i]); struct_gen();
    if (ifndef) {
      fprintf(includeoutput,"#endif /* %s */\n",ifndefname);
      }
    struct_close();

    /* The other files are optionally generated. */
    if (!is_entirely_excluded("print")) { printf(" print"); print_gen(); }
    if (!is_entirely_excluded("alloc")) { printf(" alloc"); alloc_gen(); }
    if (!is_entirely_excluded("free")) { printf(" free"); free_gen(); }
    if (!is_entirely_excluded("init")) { printf(" init"); init_gen(); }
    /* keyvalip and ip are mutually exclusive */
    if (!is_entirely_excluded("ip"))
      { printf(" ip"); ip_gen(0); }
    else if (!is_entirely_excluded("keyvalip"))
      { printf(" keyvalip"); ip_gen(1); }
    if (!is_entirely_excluded("iseq")) { printf(" iseq"); iseq_gen(); }
    if (!is_entirely_excluded("assign")) { printf(" assign"); assign_gen(); }
    if (!is_entirely_excluded("bwrite")) { printf(" bwrite"); bwrite_gen(); }
    if (!is_entirely_excluded("bread")) { printf(" bread"); bread_gen(); }
    if (!is_entirely_excluded("fwrite")) { printf(" fwrite"); fwrite_gen(); }
    if (!is_entirely_excluded("fread")) { printf(" fread"); fread_gen(); }
    if (!is_entirely_excluded("bcast0")) { printf(" bcast0"); bcast0_gen(); }
    if (!is_entirely_excluded("sbcast0")) { printf(" sbcast0"); sbcast0_gen(); }
    if (!is_entirely_excluded("rbcast0")) { printf(" rbcast0"); rbcast0_gen(); }
    if (!is_entirely_excluded("recv0")) { printf(" recv0"); recv0_gen(); }
    if (!is_entirely_excluded("send0")) { printf(" send0"); send0_gen(); }
    if (!is_entirely_excluded("zero")) { printf(" zero"); zero_gen(); }
    printf("\n");
    first = 0;
    }
  return 0;
  }

