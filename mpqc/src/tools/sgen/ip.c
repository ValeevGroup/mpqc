
/* $Log$
 * Revision 1.3  1994/10/18 23:03:54  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.2  1994/06/08  01:14:14  cljanss
 * Many changes.  These include: newmat7 and nihmatrix -> scmat
 * and mpqcic -> MPSCF and updated optimize stuff.
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1992/06/17  23:07:29  jannsen
 * modified to generate clean code
 *
 * Revision 1.1.1.1  1992/03/17  18:10:36  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:10:34  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/15  17:19:20  seidl
 * apply clj's patches
 *
 * Revision 1.5  1991/12/14  00:13:08  cljanss
 * fix bug affecting the determination of array dimensions
 *
 * Revision 1.4  1991/09/28  18:13:37  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 *
 * Revision 1.3  91/09/28  16:40:31  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.2  91/06/16  15:11:15  janssen
 * added support for boolean
 * 
 * Revision 1.1  1991/06/15  21:13:57  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <string.h>
#include <tmpl.h>
#include "types.h"
#include "global.h"

#include "ip.gbl"
#include "ip.lcl"

#include "sgen_util.gbl"
#include "error.gbl"

#define F fprintf
#define O output

GLOBAL_FUNCTION VOID
ip_gen(use_keyval_ipv2)
    int use_keyval_ipv2;
{
  declaration_t *I;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"ip");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"ip");
  strcat(outfile,".h");

  /* Open the output file. */
  include = fopen(outfile,"w");
  if (!include) {
    fprintf(stderr,"%s: couldn't open the include file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Include the following files. */
  fprintf(output,"\n");
  fprintf(output,"#include <stdio.h>\n");
  fprintf(output,"#include <stdlib.h>\n");
  fprintf(output,"#include <stdarg.h>\n");
  fprintf(output,"#include <util/sgen/sgen.h>\n");
  if (use_keyval_ipv2) {
      fprintf(output,"#include <util/keyval/ipv2c.h>\n");
    }
  else {
      fprintf(output,"#include <util/ipv2/ip_libv2.h>\n");
    }
  fprintf(output,"#include \"%s.h\"\n",BaseName);
  fprintf(output,"#include \"%sinit.h\"\n",BaseName);
  fprintf(output,"#include \"%sip.h\"\n",BaseName);

  /* Go thru the list of declarations and generate the ip functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  int ip_read_%s();\n",I->name);
    fprintf(include,"  int ip_read_%s_v();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  int ip_read_%s(char *keyword, %s_t *_%s, int n, ...);\n",
            I->name,I->name,I->name);
    fprintf(include,"  int ip_read_%s_v(char *keyword, %s_t *_%s, int n, int *v);\n",
            I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (use_keyval_ipv2 && is_excluded("keyvalip",I)) continue;
    else if (!use_keyval_ipv2 && is_excluded("ip",I)) continue;
    /* Generate the vararg version of the vector function. */
    ip_vararg(I);

    /* Write to the C source file. */
    ip_declaration(I);
    }

  fclose(output);
  }

LOCAL_FUNCTION VOID
ip_vararg(dec)
declaration_t *dec;
{
  /* The RS/6000 xlc compiler requires the following distinction for
   * var arg function declarations. */
  fprintf(output,"int\n");
  fprintf(output,"#if defined(NO_PROTO)\n");
    fprintf(output,"ip_read_%s(keyword, _%s, n)\n",dec->name,dec->name);
    fprintf(output,"char *keyword;\n");
    fprintf(output,"%s_t *_%s;\n",dec->name,dec->name);
    fprintf(output,"int n;\n");
  fprintf(output,"#else\n");
    fprintf(output,"ip_read_%s(",dec->name);
    fprintf(output,"char *keyword,");
    fprintf(output," %s_t *_%s,",dec->name,dec->name);
    fprintf(output," int n, ...)\n");
  fprintf(output,"#endif\n");

  fprintf(output,"{\n");

  fprintf(output,"  va_list args;\n");
  fprintf(output,"  int i;\n");
  fprintf(output,"  int *v;\n");

  fprintf(output,"  if (n==0) {\n");
  fprintf(output,"    return ip_read_%s_v(keyword,_%s,n,NULL);\n",
                 dec->name,dec->name);
  fprintf(output,"    }\n");
  fprintf(output,"  else {\n");
  fprintf(output,"    v = (int *) malloc(sizeof(int)*n);\n");
  fprintf(output,"    if (!v) {\n");
  fprintf(output,"      ip_warn(\"ip_read_%s: problem with malloc of %%d integers\",n);\n",dec->name);
  fprintf(output,"      return IPE_MALLOC;\n");
  fprintf(output,"      }\n");
  fprintf(output,"    va_start(args, n);\n");
  fprintf(output,"    for (i=0; i<n; i++) {\n");
  fprintf(output,"      v[i] = va_arg(args,int);\n");
  fprintf(output,"      }\n");
  fprintf(output,"    va_end(args);\n");
  fprintf(output,"    free(v);\n");
  fprintf(output,"    return ip_read_%s_v(keyword,_%s,n,v);\n",
                 dec->name,dec->name);
  fprintf(output,"    }\n");
  fprintf(output,"  }\n");

  }

LOCAL_FUNCTION VOID
ip_declaration(dec)
declaration_t *dec;
{
  member_list_t *J;

  fprintf(output,"\n/* Generated ip function for %s: */\n",dec->name);
  fprintf(output,"int\n");

  fprintf(output,"ip_read_%s_v(keyword, _%s, n, v)\n",dec->name,dec->name);
  fprintf(output,"char *keyword;\n");
  fprintf(output,"%s_t *_%s;\n",dec->name,dec->name);
  fprintf(output,"int n;\n");
  fprintf(output,"int *v;\n");

  fprintf(output,"{\n");
  declare_indices(dec->members,dec->name);
  fprintf(output,"  int errcod;\n");
  fprintf(output,"  char newkey[KEYWORD_LENGTH],basekey[KEYWORD_LENGTH];\n");

  /* Initialize the data in this structure. */
  fprintf(output,"  init_%s(_%s);\n",
          dec->name,dec->name);

  /* Construct the base keyword. */
  fprintf(output,"  ip_construct_key_v(keyword,basekey,n,v);\n");

  /* Loop over the union selectors of the structure. */
  for (J=dec->members; J!=NULL; J=J->p) {
    if (!is_union_selector(J->member,dec->members)) continue;
    ip_member(J->member,dec->name);
    }

  /* Loop over the dimensions of the structure. */
  for (J=dec->members; J!=NULL; J=J->p) {
    if (!is_array_index(J->member,dec->members)) continue;
    ip_dimension(J->member,dec);
    }

  /* Allocate storage for the elements of the structure. */
  fprintf(output,"  if (alloc_%s(_%s)) {\n",
                  dec->name,dec->name);
  fprintf(output,"    ip_warn(\"ip_read_%s_v: problem with alloc_%s\");\n",
                 dec->name,dec->name);
  fprintf(output,"    return IPE_MALLOC;\n");
  fprintf(output,"    }\n");

  /* Loop over the other elements of the structure. */
  for (J=dec->members; J!=NULL; J=J->p) {
    if ((  is_array_index(J->member,dec->members)
         ||is_union_selector(J->member,dec->members))) continue;
    ip_member(J->member,dec->name);
    }

  /* If we get here then return that everything is OK. */
  fprintf(output,"  return IPE_OK;\n");
  fprintf(output,"  }\n");
  }

/* This reads in dimensions.  If dimensions are not explicitly given
 * in the input, then the matrices with the variable dimensions are
 * examined to determine how much storage will be needed. */
LOCAL_FUNCTION VOID
ip_dimension(member,dec)
member_t *member;
declaration_t *dec;
{
  member_list_t *I;
  index_list_t *J;
  char spaces[STRING_LENGTH];
  char indices[STRING_LENGTH];

  spaces[0] = '\0';

  /* Construct the new keyword. */
  F(O,"%s  sprintf(newkey,\"%%s:%%s",spaces);
  F(O,"\",basekey,\"%s\"",member->name);
  F(O,");\n");

  F(O,"%s  if (ip_exist(newkey,0)) {\n",spaces);
  strcat(spaces,"  ");
  indices[0] = '\0';
  ip_elementary(member,dec->name,spaces,indices,0);
  F(O,"%s  }\n",spaces);
  spaces[strlen(spaces)-2] = '\0';

  /* If the keyword doesn't exist, then must examine the arrays using
   * this dimension to get the size. */
  F(O,"%s  else {\n",spaces);
  strcat(spaces,"  ");
  F(O,"%s  int count;\n",spaces);

  for (I=dec->members; I!=NULL; I=I->p) {
    indices[0] = '\0';
    for (J=I->member->indices; J!=NULL; J=J->p) {
      /* The index must be the same as the name of the integer we seek. */
      if (  (J->index.type == TYPE_CONSTANT)
          ||strcmp(J->index.v.v,member->name)) {
        strcat(indices,":0");
        continue;
        }
      F(O,"%s  sprintf(newkey,\"%%s:%%s%s\",basekey,\"%s\");\n",
          spaces,indices,I->member->name);
      F(O,"%s  if (ip_count_v(newkey,&count,0,NULL) == IPE_OK) {\n",spaces);
      strcat(spaces,"  ");
      F(O,"%s  if (count > %s) %s = count;\n",
           spaces,member_name(dec->name,member),member_name(dec->name,member));
      F(O,"%s  }\n",spaces);
      spaces[strlen(spaces)-2] = '\0';
      strcat(indices,":0");
      }
    }

  F(O,"%s  }\n",spaces);
  spaces[strlen(spaces)-2] = '\0';

  }

LOCAL_FUNCTION VOID
ip_member(member,structname)
member_t *member;
char *structname;
{
  char indices[STRING_LENGTH];
  char spaces[STRING_LENGTH];
  index_list_t *I;
  int i;
  int n_indices;

  spaces[0] = '\0';

  /* If this is a union then we conditional execute the following code. */
  if (member->uname) {
    fprintf(output,"  if (_%s->%s==%s) {\n",
            structname,member->uselname,member->uselval);
    strcat(spaces,"  ");
    }

  /* Compute the number of indices on the member. */
  for ((I=member->indices),(n_indices=0);  I!=NULL; (I=I->p),(n_indices++));

  /* We only deal with the member if all of its indices have nonzero
   * dimension. */
  if (member->indices) {
    fprintf(output,"%s  if (",spaces);
    for (I=member->indices; I!=NULL; I=I->p) {
      fprintf(output,"%s",index_dimension(structname,&I->index));
      if (I->p) fprintf(output," && ");
      }
    fprintf(output,") {\n");
    strcat(spaces,"  ");
    }

  /* Generate code to read in the data. */
  indices[0] = '\0';
  if (member->indices) strcat(spaces,"  ");
  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    strcat(indices,"[");
    indices[strlen(indices)+1] = '\0';
    indices[strlen(indices)] = 'i' + i;
    strcat(indices,"]");
    strcat(spaces,"  ");
    fprintf(output,"%sfor (%c=0; %c<%s; %c++)",spaces,
            'i'+i,'i'+i,index_dimension(structname,&I->index),'i'+i);
    fprintf(output,"  {\n");
    }

  /* Construct the new keyword. */
  fprintf(output,"%s  sprintf(newkey,\"%%s:%%s",spaces);
  for (I=member->indices; I!=NULL; I=I->p) {
    fprintf(output,":%%d");
    }
  fprintf(output,"\",basekey,\"%s\"",member->name);
  for (I=member->indices,i=0; I!=NULL; I=I->p,i++) {
    fprintf(output,",%c",'i'+i);
    }
  fprintf(output,");\n");

  ip_elementary(member,structname,spaces,indices,n_indices);

  for (I=member->indices;  I!=NULL; I=I->p) {
    fprintf(output,"%s  }\n",spaces);
    spaces[strlen(spaces)-2] = '\0';
    }

  /* Finish the test for nonzero dimensions. */
  if (member->indices) {
    spaces[strlen(spaces)-2] = '\0';
    fprintf(output,"%s    }\n",spaces);
    }

  if (member->uname) {
    fprintf(output,"  }\n");
    }

  }

LOCAL_FUNCTION VOID
ip_elementary(member,structname,spaces,indices,n_indices)
member_t *member;
char *structname;
char *spaces;
char *indices;
int n_indices;
{
  int i;
  char pointer[STRING_LENGTH];
  char pointermalloc[STRING_LENGTH];

  /* If member is a pointer, then storage must be allocated for that
   * to which it points.  If the keyword does not exist, then storage
   * is not allocated. */
  if (member->pointer) {
    /* See if the keyword exists. */
    fprintf(output,"%s  if (ip_exist(newkey,0)) {\n",spaces);
    strcat(spaces,"  ");
    pointer[0] = '\0';
    for (i=0; i<member->pointer; i++) strcat(pointer,"*");
    strcpy(pointermalloc,pointer);
    pointermalloc[member->pointer-1] = '\0';
    for (i=0; i<member->pointer; i++) {
      fprintf(output,"%s  %s%s = (%s%s)malloc(sizeof(%s%s));\n",
              spaces,
              member_name(structname,member),
              indices,
              type_name(member),
              pointer,
              type_name(member),
              pointermalloc);
      fprintf(output,"%s  if (!%s%s) return IPE_MALLOC;\n",
              spaces,
              member_name(structname,member),
              indices);
      if (strlen(pointer)) pointer[strlen(pointer)-1] = '\0';
      if (strlen(pointermalloc)) pointermalloc[strlen(pointermalloc)-1] = '\0';
      strcat(indices,"[0]");
      }
    }

  if (!strcmp(member->type,"int")) {
    if (member->qualifier == Q_UNSIGNED) {
      fprintf(output,"%s  errcod = ip_read_unsigned_%s_v(newkey,&(%s%s),0,NULL);\n",
              spaces,member->type,member_name(structname,member),indices);
      }
    else {
      fprintf(output,"%s  errcod = ip_read_%s_v(newkey,&(%s%s),0,NULL);\n",
              spaces,member->type,member_name(structname,member),indices);
      }
    }
  else if (!strcmp(member->type,"long")) {
    if (member->qualifier == Q_UNSIGNED) {
      fprintf(output,"%s  errcod = ip_read_unsigned_%s_v(newkey,&(%s%s),0,NULL);\n",
              spaces,member->type,member_name(structname,member),indices);
      }
    else {
      fprintf(output,"%s  errcod = ip_read_%s_v(newkey,&(%s%s),0,NULL);\n",
              spaces,member->type,member_name(structname,member),indices);
      }
    }
  else if (!strcmp(member->type,"string")) {
    fprintf(output,"%s  errcod = ip_read_string_v(newkey,&(%s%s),0,NULL);\n",
            spaces,member_name(structname,member),indices);
    }
  else {
    /* Read the data. */
    fprintf(output,"%s  errcod = ip_read_%s_v(newkey,&(%s%s),0,NULL);\n",
            spaces,member->type,member_name(structname,member),indices);
    }

  /* Return for unexpected errors. */
  fprintf(output,"%s  if (errcod!=IPE_OK && errcod!=IPE_KEY_NOT_FOUND)",
          spaces);
  fprintf(output," return errcod;\n");

  /* Cleanup for members that are pointers. */
  if (member->pointer) {
    fprintf(output,"%s  }\n",spaces);
    spaces[strlen(spaces)-2] = '\0';
    }
  }

