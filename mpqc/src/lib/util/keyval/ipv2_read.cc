
/* This file provides the routines to do the initial parse of the input
 * file. */

#include <stdio.h>
#include <string.h>
#ifdef DEC
#include <math.h>
#else
#include <stdlib.h>
#endif
#include "ipv2.h"

/* Set the ip_uppercase global. */
void
IPV2::set_uppercase(int uc)
{
  if (uc) ip_uppercase = 1;
  else ip_uppercase = 0;
}

/* Initialize the ip routines.  This involves parsing the entire file and
 * converting it into an internal representation. */
/* in = the input file. */
/* out = the output file. */
void
IPV2::ip_initialize(FILE*in,FILE*out)
{
  ip_initialized = 1;
  
  if (in)  ip_in = in;
  else     ip_in = stdin;
  
  if (out) ip_out = out;
  else     ip_out = stdout;

  // initialize the line number
  lineno = 1;
  
  /* Just in case a scanner has already been running. */
  if (yin) yrestart(ip_in);
  else yin = ip_in;
  
  /* If ip_tree is not NULL, then ip_initialize has been called twice,
   * with a done inbetween. Call done now. */
  if (ip_tree) {
      warn("ip_initialize has been called twice without an IPV2::done");
      done();
    }
  
  sub_tree = ip_tree;
  
  yparse();
  
  /* The initial cwk list is nothing. */
  cwk_clear();
  
  ip_internal_values();
}

/* Continue adding to the ip_tree with, presumably, another input file. 
 * This should be called after ip_initialize has been called with different
 * input file.  Multiple calls to ip_append, with different input files,
 * are allowed. */
/* in = the input file. */
/* out = the output file. */
void
IPV2::ip_append(FILE*in,FILE*out)
{
  
  if (in)  {
      ip_in = in;
      lineno = 1;
      yrestart(in);
    }
  
  if (out) ip_out = out;
  
  if (sub_tree != NULL) {
      error("ip_append: sub_tree != NULL - impossible");
    }
  
#ifndef FLEX
  yylineno = 0;
#endif
  yparse();
  
  ip_internal_values();
}

/* This routine can be called by the user after the ip routines have been
 * initialized.  It look for a prefix"dir" in the cwk list and a prefix"files"
 * array.  If prefix"dir" is found this concatenated with each of the
 * prefix"files" that does not begin with a '/'.  Each of these files is
 * ip_append'ed to the set of inputs. */
void
IPV2::append_from_input(const char*prefix,FILE*outfile)
{
  char keyword[KEYWORD_LENGTH];
  const char *dir;
  const char *file;
  char dirfile[512];
  int i,nfile;
  FILE *infile;
  
  /* Get the prefix. */
  strcpy(keyword,prefix);
  strcat(keyword,"dir");
  if (value(keyword,&dir,0) != OK) dir = NULL;
  
  strcpy(keyword,prefix);
  strcat(keyword,"files");
  if (count(keyword,&nfile,0)!=OK) return;
  for (i=0; i<nfile; i++) {
      if (value_v(keyword,&file,1,&i) == OK) {
	  if (dir && (file[0] != '/')) strcpy(dirfile,dir);
	  else dirfile[0] = '\0';
	  strcat(dirfile,file);
	  infile = fopen(dirfile,"r");
	  if (!infile) {
	      warn("IPV2::append_from_input: couldn't open the file %s",dirfile);
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
void
IPV2::ip_internal_values()
{
  int errcod;
  
  errcod = boolean(":ip:keyword",&ip_keyword,0);
  if (errcod) ip_keyword = 0;
  
}

/* Free all of the data. */
void
IPV2::done()
{
  ip_free_keyword_tree(ip_tree);
  ip_tree = NULL;
  sub_tree = NULL;
  ip_in = NULL;
  ip_out = NULL;
}

void
IPV2::ip_push_keyword(char*keyword)
{
  ip_push_keyclass(keyword,0,0);
}

void
IPV2::ip_push_keyclass(char*keyword,char*classname,ip_string_list_t*parentlist)
{
  ip_keyword_tree_t *I, *new_keyword;

  // if there is no keyword, then set the classname of the current
  // sub_tree
  if (!keyword) {
      if (classname && sub_tree && sub_tree->keyword && !sub_tree->classname) {
          sub_tree->classname = classname;
        }
      else {
          error("got a classname only in invalid context: %k");
        }
      return;
    }
  
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
	      error("A class specification was given for only one keyword: %k\n");
	    }
	  else if (classname) {
	      if (strcmp(classname,I->classname)) {
		  error("Class specifications differ for keyword %k\n");
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

void
IPV2::ip_pop_keyword()
{
  /* Make sure we aren\'t already on top. */
  if (!sub_tree) {
      error("ip_pop_keyword: tried to pop above top");
    }
  sub_tree = sub_tree->up;
}

void
IPV2::ip_begin_table(ip_string_list_t*keywords)
{
  //fprintf(stderr,"ip_begin_table()\n");
  current_table_keyword = table_keywords = keywords;
  table_sub_tree = sub_tree;
  table_row_number = 0;
}

/* Given a string containing keywords separated by ':', push the
 * keywords. */
void
IPV2::ip_push_table_col(char*keys)
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

void
IPV2::ip_next_table_entry()
{
  
  //fprintf(stderr,"ip_next_table_entry: table_array_depth = %d\n",table_array_depth);

  if (table_array_depth>0) return;

  //fprintf(stderr,"ip_next_table_entry: table keyword = \"");
  //ip_print_keyword(stderr,table_sub_tree);
  //fprintf(stderr,"\"\n");
  
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

void
IPV2::ip_done_table()
{
  ip_string_list_t *I,*J;

  //fprintf(stderr,"ip_done_table()\n");

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
ip_string_list_t *
IPV2::ip_add_string_list(ip_string_list_t*sl,char*s)
{
  ip_string_list_t *I;
  
  if (!sl) return ip_string_to_string_list(s);
  
  for (I=sl; I->p!=NULL; I=I->p);
  I->p = ip_string_to_string_list(s);
  return sl;
}

ip_string_list_t *
IPV2::ip_string_to_string_list(char*s)
{
  ip_string_list_t *r;
  r = (ip_string_list_t *) malloc(sizeof(ip_string_list_t));
  r->string = s;
  r->p = NULL;
  return r;
}

char *
IPV2::dup_string(const char*s)
{
  char *r;
  r = (char *) malloc(strlen(s)+1);
  strcpy(r,s);
  return r;
}

ip_keyword_tree_t *
IPV2::ip_get_variable_kt(char* variable)
{
  char* passed_variable = variable;
  ip_keyword_tree_t *kt;
  ip_keyword_tree_t *top;

  top = sub_tree;

  /* One or more occurrences of "..:" at the beginning of the keyword
   * move us up the keyword tree */
  while(top && !strncmp(variable,"..:",3)) {
      variable = &variable[3];
      top = top->up;
    }
  
  /* If top is still then we have a problem. */
  if (!top) {
      error("tried to get a variable above the top level - impossible");
    }

  /* Descend the keyword tree, creating nodes if needed. */
  if (variable[0] == ':') {
      kt = ip_descend_tree(ip_tree,variable);
    }
  else {
      kt = ip_descend_tree(top,variable);
    }

  /* This should never be the case since variable keyword trees are
   * created as needed. */
  if (!kt) {
      warn("couldn't find the variable %s",variable);
      return NULL;
    }
  
  /* Release storage for the variable. */
  free(passed_variable);
  
  return(kt);
}

void
IPV2::ip_assign_variable(char* variable)
{
  if (table_keywords) ip_next_table_entry();

  /* Note that the subtree is really a reference to another subtree. */
  sub_tree->variable = (char*)malloc(strlen(variable)+1);
  strcpy(sub_tree->variable,variable);
}

char *
IPV2::ip_get_variable_value(char*variable)
{
  ip_keyword_tree_t *kt;
  
  /* Get the keyword tree associated with the variable. */
  kt = ip_get_variable_kt(variable);
  
  /* Return the value associated with the keyword. */
  if (kt) return(kt->value);
  else return NULL;
}

double
IPV2::ip_get_variable_double(char*variable)
{
  char *value;
  
  value = ip_get_variable_value(variable);
  
  if (value == NULL) return 0.0;
  else return atof(value);
}

char *
IPV2::ip_double_to_string(double val)
{
  char *result;
  
  result = (char *) malloc(64);
  
  sprintf(result,"%22.15e",val);
  
  return result;
}

void
IPV2::ip_assign_value(char*value)
{

  if (table_keywords) ip_next_table_entry();

  /* If sub_tree is still NULL then we have a problem. */
  if (!sub_tree) {
    error("tried to put a keyword at the top level - impossible");
    }

  /* Check for duplicate definitions. */
  if (sub_tree->value) {
#   ifdef DUP_WARN
      /* Warn the user about duplicate definitions. */
      warn("duplicate definition of the following keyword:");
      ip_print_keyword(ip_out,sub_tree);
      fprintf(ip_out,"\n");
      warn("the new value will be ignored");
#   endif /* DUP_WARN */
      free(value);
    }
  else sub_tree->value = value;
  //ip_print_keyword(stderr,sub_tree);
  //fprintf(stderr," assigned to \"%s\"\n",value);
  }

void
IPV2::ip_start_karray()
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

void
IPV2::ip_init_karray()
{
  init_karray = 1;
  karray_indices->i = 0;
}

void
IPV2::ip_incr_karray()
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

void
IPV2::ip_pop_karray()
{
  intlist_t *tmp;
  if (table_keywords) table_array_depth--;
  if (!init_karray) ip_pop_keyword();
  tmp = karray_indices;
  karray_indices = karray_indices->p;
  if (tmp) free(tmp);
}

char *
IPV2::ip_append_keystrings(char*s1,char*s2)
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
