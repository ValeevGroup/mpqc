
#include <stdlib.h>
#include <util/keyval/ipv2.h>

ip_keyword_tree_t *
IPV2::ip_alloc_keyword_tree()
{
  ip_keyword_tree_t *result;

  result = (ip_keyword_tree_t *) malloc(sizeof(ip_keyword_tree_t));
  if (!result) {
    cerr << "ip_alloc_keyword_tree: malloc failed";
    error(0);
    }

  result->up = 0;
  result->down = 0;
  result->across = 0;
  result->keyword = 0;
  result->classname = 0;
  result->truename = 0;
  result->value = 0;
  result->variable = 0;

  return result;
  }

void
IPV2::ip_free_keyword_tree(ip_keyword_tree_t* tree)
{
  ip_keyword_tree_t *I,*start,*nextI;

  if (!tree) return;

  /* Convert the circular list into a standard linked list (to
   * avoid saber-c error messages) */
  start = tree->across;
  tree->across = 0;
  for (I=start; I!=0; I=nextI) {
    ip_free_keyword_tree(I->down);
    if (I->keyword) free(I->keyword);
    if (I->classname) free(I->classname);
    if (I->truename) free(I->truename);
    if (I->value) free(I->value);
    if (I->variable) free(I->variable);
    nextI = I->across;
    free(I);
    }
  }

