
#include <iostream.h>
#include <util/keyval/ipv2.h>

#define N_INDENT 2

void
IPV2::ip_print_keyword(ostream&fp,ip_keyword_tree_t*st)
{
  if (st->up) ip_print_keyword(fp,st->up);
  fp << st->keyword << ":";
  }

/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out. */
void
IPV2::ip_print_tree(ostream&fp,ip_keyword_tree_t*tree)
{
  if (!tree) tree = ip_tree;
  if (!tree) return;

  ip_print_tree_(fp,tree,0);
  }


/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out.  Indent is used to record how deep in the tree we
 * are, so we know how far to indent things. */
void
IPV2::ip_print_tree_(ostream&fp,ip_keyword_tree_t*tree,int indent)
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
      fp << "\"" << I->keyword << "\"" << endl;
      }
    else {
      fp << I->keyword << endl;
      }

    if (I->classname) {
      fp << "<" << I->keyword << ">" << endl;
      }

    if (!(I->value || I->down || I->variable)) {
      fp << ": (" << endl;
      }

    if (I->variable) {
      fp << " = $" << I->variable << endl;
      }
    if (I->truename) {
      fp << "\"" << I->truename << "\"";
      }

    if (I->value) {
      if (I->down) fp << " (= " << I->value << ")";
      else fp << " = " << I->value << endl;
      }
    if (I->down) {
      fp << ": (" << endl;
      ip_print_tree_(fp,I->down,indent + N_INDENT);
      ip_indent(fp,indent + N_INDENT);
      fp << ")" << endl;
      }

    } while ((I = I->across) != tree);

  }

void
IPV2::ip_indent(ostream&fp,int n)
{
  int i;

  for (i=0; i<n; i++) fp << " ";
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
