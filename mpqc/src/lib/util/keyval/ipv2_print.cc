
#include <stdio.h>
#include "ipv2.h"

#define N_INDENT 2

void
IPV2::ip_print_keyword(FILE*fp,ip_keyword_tree_t*st)
{
  if (st->up) ip_print_keyword(fp,st->up);
  fprintf(fp,"%s:",st->keyword);
  }

/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out. */
void
IPV2::ip_print_tree(FILE*fp,ip_keyword_tree_t*tree)
{
  if (!fp) fp = stdout;
  if (!tree) tree = ip_tree;
  if (!tree) return;

  ip_print_tree_(fp,tree,0);
  }


/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out.  Indent is used to record how deep in the tree we
 * are, so we know how far to indent things. */
void
IPV2::ip_print_tree_(FILE*fp,ip_keyword_tree_t*tree,int indent)
{
  ip_keyword_tree_t *I;

  I=tree;
  do {
    //if (I->value && I->down) {
    //  warn("ip_print_tree: tree has both value and subtrees - can't print");
    //  warn("keyword is %s, value is %s, subtree key is %s\n",
    //       I->keyword,I->value,I->down->keyword);
    //  }

    if (!I->keyword) {
      warn("ip_print_tree: tree has no keyword - impossible");
      }

    ip_indent(fp,indent);
    if (ip_special_characters(I->keyword)) {
      fprintf(fp,"\"%s\"",I->keyword);
      }
    else {
      fprintf(fp,"%s",I->keyword);
      }

    if (I->classname) {
      fprintf(fp,"<%s>",I->classname);
      }

    if (!(I->value || I->down || I->variable)) {
      printf(": (\n");
      }

    if (I->variable) {
        fprintf(fp," = $%s\n",I->variable);
      }
    if (I->truename) {
        fprintf(fp,"\"%s\"",I->truename);
      }

    if (I->value) {
      if (I->down) fprintf(fp," (= %s)",I->value);
      else fprintf(fp," = %s\n",I->value);
      }
    if (I->down) {
      fprintf(fp,": (\n");
      ip_print_tree_(fp,I->down,indent + N_INDENT);
      ip_indent(fp,indent + N_INDENT);
      fprintf(fp,")\n");
      }

    } while ((I = I->across) != tree);

  }

void
IPV2::ip_indent(FILE*fp,int n)
{
  int i;

  for (i=0; i<n; i++) fprintf(fp," ");
  }

int
IPV2::ip_special_characters(char*keyword)
{
  char *ch=keyword;

  if (!keyword) return 0;
  while (*ch) {
    if (!(  (*ch >= 'a' && *ch <= 'z')
          ||(*ch >= 'A' && *ch <= 'Z')
          ||(*ch >= '0' && *ch <= '9')
          ||(*ch == '<')
          ||(*ch == '>')
          ||(*ch == '+')
          ||(*ch == '-')
          ||(*ch == '.')
          ||(*ch == '_'))) return 1;

    ch++;
    }

  return 0;
  }
