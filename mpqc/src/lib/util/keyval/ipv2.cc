
#include "ipv2.h"
#include <util/keyval/ipv2_parse.h>

IPV2::IPV2():
ip_in(0),
ip_out(0),
ip_tree(0),
ip_cwk(0),
ip_uppercase(0),
ip_keyword(0),
table_keywords(0),
current_table_keyword(0),
table_sub_tree(0),
karray_indices(0),
sub_tree(0),
ip_initialized(0),
cwkstack(0),
table_array_depth(0)
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
