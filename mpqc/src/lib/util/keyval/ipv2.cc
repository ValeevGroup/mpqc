#ifdef __GNUG__
#pragma implementation
#endif

#include <stdlib.h>

#include "ipv2.h"
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
ip_uppercase(0),
ip_keyword(0)
{
  lastkeyword[0] = '\0';
}

IPV2::~IPV2()
{
  if (ip_tree) ip_free_keyword_tree(ip_tree);
  ip_tree = 0;
}

void IPV2::read(FILE*in,FILE*out)
{
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

int IPV2::ywrap()
{
  return 1;
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
      fprintf(stderr,"IPV2::global: global not set\n");
      abort();
    }
  return global_;
}

  // some routines for debugging
void
IPV2::print_keyword(FILE*f,ip_keyword_tree_t*k)
{
  ip_print_keyword(f,k);
}

void
IPV2::print_tree(FILE*f,ip_keyword_tree_t*k)
{
  ip_print_tree(f,k);
}
