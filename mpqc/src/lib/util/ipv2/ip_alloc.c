
#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include "ip_types.h"
#include "ip_global.h"

#include "ip_alloc.gbl"
#include "ip_alloc.lcl"

#include "ip_error.gbl"

GLOBAL_FUNCTION ip_keyword_tree_t *
ip_alloc_keyword_tree()
{
  ip_keyword_tree_t *result;

  result = (ip_keyword_tree_t *) malloc(sizeof(ip_keyword_tree_t));
  if (!result) {
    perror("ip_alloc_keyword_tree: malloc failed");
    ip_error(NULL);
    }

  result->up = NULL;
  result->down = NULL;
  result->across = NULL;
  result->keyword = NULL;
  result->classname = NULL;
  result->value = NULL;

  return result;
  }

GLOBAL_FUNCTION VOID
ip_free_keyword_tree(tree)
ip_keyword_tree_t *tree;
{
  ip_keyword_tree_t *I,*start,*nextI;

  if (!tree) return;

  /* Convert the circular list into a standard linked list (to
   * avoid saber-c error messages) */
  start = tree->across;
  tree->across = NULL;
  for (I=start; I!=NULL; I=nextI) {
    ip_free_keyword_tree(I->down);
    if (I->keyword) free(I->keyword);
    if (I->value) free(I->value);
    nextI = I->across;
    free(I);
    }
  }

