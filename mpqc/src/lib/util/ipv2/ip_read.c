
/* This file provides the routines to do the initial parse of the input
 * file. */

#include <stdio.h>
#include <string.h>
#ifdef DEC
#include <math.h>
#else
#include <stdlib.h>
#endif
#include <tmpl.h>
#include "ip_types.h"
#define _IP_ALLOCATE_GLOBAL_
#include "ip_global.h"
#include "parse.h"

/* The input for yacc. */
extern FILE *yyin;
#ifndef FLEX
extern int yylineno;
#endif

#include "ip_read.gbl"
#include "ip_read.lcl"

#include "ip_data.gbl"
#include "ip_error.gbl"
#include "ip_print.gbl"
#include "ip_alloc.gbl"
#include "ip_cwk.gbl"

#include "ip_error.h"

static ip_string_list_t* table_keywords = NULL;
static ip_string_list_t* current_table_keyword = NULL;
static ip_keyword_tree_t* table_sub_tree = NULL;
static int table_row_number;
static int table_array_depth;

/* This integer list is used to keep track of the karray index. */
struct intlist_struct {
  int i;
  struct intlist_struct *p;
  };
typedef struct intlist_struct intlist_t;

static intlist_t *karray_indices = NULL;

/* Set up static variables. */
static ip_keyword_tree_t *sub_tree = NULL;

#ifdef FLEX
static int ip_done_called=0;
#endif

/* Set the ip_uppercase global. */
GLOBAL_FUNCTION VOID
ip_set_uppercase(uc)
int uc;
{
  if (uc) ip_uppercase = 1;
  else ip_uppercase = 0;
  }

/* Initialize the ip routines.  This involves parsing the entire file and
 * converting it into an internal representation. */
/* in = the input file. */
/* out = the output file. */
GLOBAL_FUNCTION VOID
ip_initialize(in,out)
FILE *in;
FILE *out;
{

  if (in)  ip_in = yyin = in;
  else     ip_in = yyin = stdin;

  if (out) ip_out = out;
  else     ip_out = stdout;

#ifdef FLEX
  /* Just in case this is the second call to ip_initialize. */
  if (ip_done_called) yyrestart(in);
#endif

  /* If ip_tree is not NULL, then ip_initialize has been called twice,
   * with a ip_done inbetween. Call ip_done now. */
  if (ip_tree) {
    ip_warn("ip_initialize has been called twice without an ip_done");
    ip_done();
    }

  sub_tree = ip_tree;

  yyparse();

  /* The initial cwk list is nothing. */
  ip_cwk_clear();

  ip_internal_values();
  }

/* Continue adding to the ip_tree with, presumably, another input file. 
 * This should be called after ip_initialize has been called with different
 * input file.  Multiple calls to ip_append, with different input files,
 * are allowed. */
/* in = the input file. */
/* out = the output file. */
GLOBAL_FUNCTION VOID
ip_append(in,out)
FILE *in;
FILE *out;
{

#ifdef FLEX
  if (in)  {
    ip_in = in;
    yyrestart(in);
    }
#else
  if (in)  ip_in = yyin = in;
#endif

  if (out) ip_out = out;

  if (sub_tree != NULL) {
    ip_error("ip_append: sub_tree != NULL - impossible");
    }

#ifndef FLEX
  yylineno = 0;
#endif
  yyparse();

  ip_internal_values();
  }

/* This routine can be called by the user after the ip routines have been
 * initialized.  It look for a prefix"dir" in the cwk list and a prefix"files"
 * array.  If prefix"dir" is found this concatenated with each of the
 * prefix"files" that does not begin with a '/'.  Each of these files is
 * ip_append'ed to the set of inputs. */
GLOBAL_FUNCTION VOID
ip_append_from_input(prefix,outfile)
char *prefix;
FILE *outfile;
{
  char keyword[KEYWORD_LENGTH];
  char *dir;
  char *file;
  char dirfile[512];
  int i,nfile;
  FILE *infile;

  /* Get the prefix. */
  strcpy(keyword,prefix);
  strcat(keyword,"dir");
  if (ip_value(keyword,&dir,0) != IPE_OK) dir = NULL;

  strcpy(keyword,prefix);
  strcat(keyword,"files");
  if (ip_count(keyword,&nfile,0)!=IPE_OK) return;
  for (i=0; i<nfile; i++) {
    if (ip_value_v(keyword,&file,1,&i) == IPE_OK) {
      if (dir && (file[0] != '/')) strcpy(dirfile,dir);
      else dirfile[0] = '\0';
      strcat(dirfile,file);
      infile = fopen(dirfile,"r");
      if (!infile) {
        ip_warn("ip_append_from_input: couldn't open the file %s",dirfile);
        }
      else {
        fprintf(outfile,"appending %s to input\n",dirfile);
        ip_append(infile,outfile);
        fclose(infile);
        }
      }
    }
  }


/* Set up internal ip variables based on the input that has been read in. */
GLOBAL_FUNCTION VOID
ip_internal_values()
{
  int errcod;

  errcod = ip_boolean(":ip:keyword",&ip_keyword,0);
  if (errcod) ip_keyword = 0;

  }

/* Free all of the data. */
GLOBAL_FUNCTION VOID
ip_done()
{
  ip_free_keyword_tree(ip_tree);
  ip_tree = NULL;
  sub_tree = NULL;
  ip_in = NULL;
  ip_out = NULL;
#ifdef FLEX
  ip_done_called = 1;
#endif
  }

GLOBAL_FUNCTION VOID
ip_push_keyword(keyword)
char *keyword;
{
  ip_push_keyclass(keyword,0,0);
  }

GLOBAL_FUNCTION VOID
ip_push_keyclass(keyword,classname,parentlist)
char *keyword;
char *classname;
ip_string_list_t *parentlist;
{
  ip_keyword_tree_t *I, *new_keyword;

  /* Make the parentlist a part of the keyword. */
  if (parentlist) {
    int newkeysize = strlen(keyword) + 4 + 1;
    ip_string_list_t *pl;
    char* newkey;

    for (pl=parentlist; pl != NULL; pl=pl->p) {
      newkeysize += strlen(pl->string);
      if (pl->p) newkeysize++;
      }

    newkey = (char*)malloc(newkeysize);
    strcpy(newkey,keyword);
    strcat(newkey,"<<");

    for (pl=parentlist; pl != NULL; pl=pl->p) {
      strcat(newkey,pl->string);
      if (pl->p) strcat(newkey,",");
      }
    strcat(newkey,">>");

    free(keyword);
    keyword = newkey;
    }

  /* If this is the first keyword, then create the tree. */
  if (!ip_tree) {
    sub_tree = ip_tree = ip_alloc_keyword_tree();
    sub_tree->across = sub_tree;
    sub_tree->keyword = keyword;
    sub_tree->classname = classname;
    return;
    }

  /* This is not the first keyword, so descend the tree. */

  /* If sub_tree is at the top (NULL), then move to ip_tree. */
  if (!sub_tree) {
    sub_tree = ip_tree;
    }
  /* If there is not already a sub_tree->down, then create it. */
  else if (!sub_tree->down) {
    sub_tree->down = ip_alloc_keyword_tree();
    sub_tree->down->across = sub_tree->down;
    sub_tree->down->up = sub_tree;

    sub_tree = sub_tree->down;
    sub_tree->keyword = keyword;
    sub_tree->classname = classname;
    return;
    }
  /* Descend the tree, but keep track of where we were. */
  else {
    sub_tree = sub_tree->down;
    }

  /* Does the keyword exist in the current sub tree? */
  I=sub_tree;
  do {

    if (!strcmp(I->keyword,keyword)) {
      /* We found it. */
      sub_tree = I;
      if ((classname && (!I->classname))||((!classname) && I->classname)) {
        ip_error("A class specification was given for only one keyword: %k\n");
        }
      else if (classname) {
        if (strcmp(classname,I->classname)) {
          ip_error("Class specifications differ for keyword %k\n");
          }
        }
      free(keyword);
      if (classname) free(classname);
      return;
      }

    } while ((I = I->across) != sub_tree);

  /* We could not find it -- create a new entry. */

  new_keyword = ip_alloc_keyword_tree();
  new_keyword->across = sub_tree->across;
  new_keyword->keyword = keyword;
  new_keyword->classname = classname;
  sub_tree->across = new_keyword;

  new_keyword->up = sub_tree->up;

  /* Move us down to the new keyword. */
  sub_tree = new_keyword;
  }

GLOBAL_FUNCTION VOID
ip_pop_keyword()
{
  /* Make sure we aren't already on top. */
  if (!sub_tree) {
    ip_error("ip_pop_keyword: tried to pop above top");
    }
  sub_tree = sub_tree->up;
  }

GLOBAL_FUNCTION VOID
ip_begin_table(keywords)
ip_string_list_t *keywords;
{
  current_table_keyword = table_keywords = keywords;
  table_sub_tree = sub_tree;
  table_row_number = 0;
  }

/* Given a string containing keywords separated by ':', push the
 * keywords. */
LOCAL_FUNCTION VOID
ip_push_table_col(keys)
char *keys;
{
  char cindex[10];
  char * tmp = dup_string(keys);
  char * keyword = strtok(tmp,":");
  int n = 0;
  do {
    ip_push_keyword(dup_string(keyword));
    n++;
    } while((keyword = strtok(NULL,":")) != NULL);
  free(tmp);
  sprintf(cindex,"%d",table_row_number);
  ip_push_keyword(dup_string(cindex));
  }

LOCAL_FUNCTION VOID
ip_next_table_entry()
{

  if (table_array_depth>0) return;

  sub_tree = table_sub_tree;
  ip_push_table_col(current_table_keyword->string);

  /* Advance the current_table_keyword pointer */
  if (current_table_keyword->p == NULL) {
    current_table_keyword = table_keywords;
    table_row_number++;
    }
  else {
    current_table_keyword = current_table_keyword->p;
    }
  }

GLOBAL_FUNCTION VOID
ip_done_table()
{
  ip_string_list_t *I,*J;
  /* Free the keywords strings and string list */
  for (I=table_keywords; I!=NULL; ) {
    free(I->string);
    J = I->p;
    free(I);
    I = J;
    }
  table_keywords = NULL;
  current_table_keyword = NULL;
  sub_tree = table_sub_tree;
  table_sub_tree = NULL;
  }

/* This adds the string, s, to the string list linked list, sl. */
GLOBAL_FUNCTION ip_string_list_t *
ip_add_string_list(sl,s)
ip_string_list_t *sl;
char *s;
{
  ip_string_list_t *I;

  if (!sl) return ip_string_to_string_list(s);

  for (I=sl; I->p!=NULL; I=I->p);
  I->p = ip_string_to_string_list(s);
  return sl;
  }

GLOBAL_FUNCTION ip_string_list_t *
ip_string_to_string_list(s)
char *s;
{
  ip_string_list_t *r;
  r = (ip_string_list_t *) malloc(sizeof(ip_string_list_t));
  r->string = s;
  r->p = NULL;
  return r;
  }

LOCAL_FUNCTION char *
dup_string(s)
char *s;
{
  char *r;
  r = (char *) malloc(strlen(s)+1);
  strcpy(r,s);
  return r;
  }

LOCAL_FUNCTION ip_keyword_tree_t *
ip_get_variable_kt(variable)
char *variable;
{
  ip_keyword_tree_t *kt;

  /* If sub_tree is still NULL then we have a problem. */
  if (!sub_tree) {
    ip_error("tried to put a keyword at the top level - impossible");
    }

  /* Push the old path and set up the new one. */
  ip_cwk_push();
  ip_cwk_clear();
  ip_cwk_add_kt(sub_tree);

  /* Descend the keyword tree. */
  kt = ip_cwk_descend_tree(variable);

  if (!kt) {
    ip_warn("couldn't find the variable %s",variable);
    return NULL;
    }

  /* Release storage for the variable. */
  free(variable);

  /* Restore the old path. */
  ip_cwk_pop();

  return(kt);
  }

GLOBAL_FUNCTION VOID
ip_assign_variable(variable)
char *variable;
{
  ip_keyword_tree_t *kt;

  if (table_keywords) ip_next_table_entry();

  /* Get the keyword tree associated with the variable. */
  kt = ip_get_variable_kt(variable);

  /* Copy kt to the current subtree. */
  if (kt) ip_copy_keyword_tree(sub_tree,kt);
  }

LOCAL_FUNCTION char *
ip_get_variable_value(variable)
char *variable;
{
  ip_keyword_tree_t *kt;

  /* Get the keyword tree associated with the variable. */
  kt = ip_get_variable_kt(variable);

  /* Return the value associated with the keyword. */
  if (kt) return(kt->value);
  else return NULL;
  }

GLOBAL_FUNCTION double
ip_get_variable_double(variable)
char *variable;
{
  char *value;

  value = ip_get_variable_value(variable);

  if (value == NULL) return 0.0;
  else return atof(value);
  }

GLOBAL_FUNCTION char *
ip_double_to_string(val)
double val;
{
  char *result;

  result = (char *) malloc(64);

  sprintf(result,"%22.15e",val);

  return result;
  }

/* kt1 is made to look like kt2. */
LOCAL_FUNCTION VOID
ip_copy_keyword_tree(kt1,kt2)
ip_keyword_tree_t *kt1;
ip_keyword_tree_t *kt2;
{
  ip_keyword_tree_t *new_kt,*I;
  int first = 1;

  if (!kt2) return;
  if (!kt1) return;

  if (kt2->value) {
    if (kt1->value) return; /* Duplicate value--return. */
    kt1->value = (char *) malloc(strlen(kt2->value)+1);
    strcpy(kt1->value,kt2->value);
    return;
    }

  if (!kt2->down) return;

  I = kt2->down;
  do {
    if (first) {
      first = 0;
      kt1->down = ip_alloc_keyword_tree();
      kt1->down->up = kt1;
      kt1->down->across = kt1->down;
      kt1->down->down = NULL;
      kt1->down->value = NULL;
      kt1 = kt1->down;
      new_kt = kt1;
      }
    else {
      new_kt = ip_alloc_keyword_tree();
      new_kt->across = kt1->across;
      new_kt->up = kt1->up;
      new_kt->down = NULL;
      kt1->across = new_kt;
      }
    new_kt->keyword = (char *) malloc(strlen(I->keyword)+1);
    strcpy(new_kt->keyword,I->keyword);
    ip_copy_keyword_tree(new_kt,I);
    } while ((I = I->across) != kt2->down);
  }

GLOBAL_FUNCTION VOID
ip_assign_value(value)
char *value;
{

  if (table_keywords) ip_next_table_entry();

  /* If sub_tree is still NULL then we have a problem. */
  if (!sub_tree) {
    ip_error("tried to put a keyword at the top level - impossible");
    }

  /* Check for duplicate definitions. */
  if (sub_tree->value) {
#   ifdef DUP_WARN
    /* Warn the user about duplicate definitions. */
    ip_warn("duplicate definition of the following keyword:");
    ip_print_keyword(ip_out,sub_tree);
    fprintf(ip_out,"\n");
    ip_warn("the new value will be ignored");
#   endif /* DUP_WARN */
    free(value);
    }
  else sub_tree->value = value;
  }

/* These routines maintain the keyword index stack. */
static int init_karray;

GLOBAL_FUNCTION VOID
ip_start_karray()
{
  if (table_keywords && table_array_depth == 0) ip_next_table_entry(); 
  if (table_keywords) table_array_depth++;
  if (!karray_indices) {
    karray_indices = (intlist_t *) malloc(sizeof(intlist_t));
    karray_indices->p = NULL;
    }
  else {
    intlist_t *tmp;
    tmp = (intlist_t *) malloc(sizeof(intlist_t));
    tmp->p = karray_indices;
    karray_indices = tmp;
    }
  }

GLOBAL_FUNCTION VOID
ip_init_karray()
{
  init_karray = 1;
  karray_indices->i = 0;
  }

GLOBAL_FUNCTION VOID
ip_incr_karray()
{
  char *key;

  if (init_karray) {
    init_karray = 0;
    }
  else {
    ip_pop_keyword();
    }

  /* Construct a keyword to push. */
    /* A cheap, yet inaccurate estimate of the string size needed. */
  key = (char *) malloc(karray_indices->i/1000 + 4);
  sprintf(key,"%d",karray_indices->i);

  ip_push_keyword(key);

  /* Increment the current karray index. */
  karray_indices->i++;
  }

GLOBAL_FUNCTION VOID
ip_pop_karray()
{
  intlist_t *tmp;
  if (table_keywords) table_array_depth--;
  if (!init_karray) ip_pop_keyword();
  tmp = karray_indices;
  karray_indices = karray_indices->p;
  if (tmp) free(tmp);
  }

GLOBAL_FUNCTION char *
ip_append_keystrings(s1,s2)
char *s1;
char *s2;
{
  char *r;

  if (s1) r = (char *) malloc(strlen(s1)+strlen(s2)+2);
  else    r = (char *) malloc(strlen(s2)+2);
  r[0] = '\0';
  if (s1) strcat(r,s1);
  strcat(r,":");
  strcat(r,s2);
  return r;
  }
