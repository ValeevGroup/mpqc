#ifdef __GNUG__
#pragma implementation
#endif

#include <stdlib.h>
#include <string.h>

#include <util/keyval/ipv2.h>
#include <util/keyval/ipv2_parse.h>

// IPV2's static member
IPV2* IPV2::global_ = 0;

IPV2::IPV2():
table_keywords(0),
current_table_keyword(0),
table_sub_tree(0),
table_array_depth(0),
karray_indices(0),
sub_tree(0),
cwkstack(0),
ip_initialized(0),
ip_in(0),
ip_out(0),
ip_tree(0),
ip_cwk(0),
ip_keyword(0),
filename_(0)
{
  lastkeyword[0] = '\0';
  lexer = new IPV2FlexLexer;
}

IPV2::~IPV2()
{
  if (ip_tree) ip_free_keyword_tree(ip_tree);
  ip_tree = 0;
  delete lexer;
  delete[] filename_;
}

void IPV2::read(istream&in,ostream&out,const char *filename)
{
  delete[] filename_;
  if (filename) {
      filename_ = strcpy(new char[strlen(filename)+1], filename);
    }
  else filename_ = 0;
  if (ip_initialized) {
    ip_append(in,out);
    }
  else {
    ip_initialize(in,out);
    cwk_root();
    }
}

void IPV2::yerror(const char* s)
{
  error(s);
}

void
IPV2::set_global(IPV2* i)
{
  global_ = i;
}

int
IPV2::have_global()
{
  return global_ != 0;
}

IPV2*
IPV2::global()
{
  if (!global_) {
      cerr << "IPV2::global: global not set" << endl;
      abort();
    }
  return global_;
}

  // some routines for debugging
void
IPV2::print_keyword(ostream&f,ip_keyword_tree_t*k)
{
  ip_print_keyword(f,k);
}

void
IPV2::print_tree(ostream&f,ip_keyword_tree_t*k)
{
  ip_print_tree(f,k);
}

/* Show position. */
void
IPV2::showpos()
{
  cerr << "error occurred at line number "
       << lexer->lineno() << " (roughly)" << endl;
  if (filename_) {
      cerr << "in file \"" << filename_ << "\"";
    }
  cerr << endl;
}

int
IPV2wrap()
{
  return 1;
}
