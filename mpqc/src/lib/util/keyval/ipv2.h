
#ifndef _util_keyval_ipv2_ipv2_h
#define _util_keyval_ipv2_ipv2_h
#ifdef __GNUG__
#pragma interface
#endif

#include <iostream.h>
#include <util/keyval/ipv2_scan.h>

#undef yyFlexLexer
#define yyFlexLexer IPV2FlexLexer
#include <FlexLexer.h>

// For temporary data (only used while parsing)
/* This integer list is used to keep track of the karray index. */
struct intlist_struct {
  int i;
  struct intlist_struct *p;
  };
typedef struct intlist_struct intlist_t;

// For permanent data
struct ip_keyword_tree_struct {
  char *keyword;
  char *classname;
  char *truename;
  struct ip_keyword_tree_struct *across; /* Circular list. */
  struct ip_keyword_tree_struct *up;    /* Terminated by NULL. */
  struct ip_keyword_tree_struct *down;  /* Terminated by NULL. */
  char *variable;  /* If this node points to another name, this
                    * is the name, otherwise NULL. */
  char *value;
  };

struct ip_keyword_tree_list_struct {
  struct ip_keyword_tree_struct *kt;
  struct ip_keyword_tree_list_struct *p;
  };

struct ip_cwk_stack_struct {
  struct ip_keyword_tree_list_struct *ktl;
  struct ip_cwk_stack_struct *p;
  };
typedef struct ip_cwk_stack_struct ip_cwk_stack_t;

typedef struct ip_keyword_tree_struct ip_keyword_tree_t;
typedef struct ip_keyword_tree_list_struct ip_keyword_tree_list_t;

class IPV2
{
 public:
  enum Status {
      OK=0          ,  /* No problem. */
      KeyNotFound=1 ,  /* The keyword was not found. */
      OutOfBounds=2 ,  /* An array subscript was out of bounds. */
      Malloc=3      ,  /* Memory allocation failed. */
      NotAnArray=4  ,  /* Gave index for data which isn't an array */
      NotAScalar=5  ,  /* Didn't give index for data which is an array */
      Type=6        ,  /* The datum is not of the appropiate type. */
      HasNoValue=7  ,  /* The keyword has no value. */
      ValNotExpd=8     /* A value was not expected for the keyword. */
      };
  enum { KEYWORD_LENGTH=256 };
  
 private:
  // this is what is read by the ip_data, etc, func calls
  // no IPV2 members allocate or delete storage for global_
  static IPV2* global_;

  char *filename_;
    
  // These are needed only when the input is being read in:
  ip_string_list_t* table_keywords;
  ip_string_list_t* current_table_keyword;
  ip_keyword_tree_t* table_sub_tree;
  int table_row_number;
  int table_array_depth;
  intlist_t *karray_indices;
  ip_keyword_tree_t *sub_tree;
  int init_karray;

  // this maintains a list of current working keyword lists (for cwk_push
  // and cwk_pop)
  ip_cwk_stack_t *cwkstack;

  // This keeps track of whether or not we've been initialized
  int ip_initialized;

  // This is used for error processing
  char lastkeyword[KEYWORD_LENGTH];
  
  // These are needed always:
  istream* ip_in;
  ostream* ip_out;
  ip_keyword_tree_t* ip_tree;
  ip_keyword_tree_list_t* ip_cwk;
  int ip_keyword;

  // private routines mainly used for parsing the input
  void ip_push_table_col(char*);
  void ip_next_table_entry();
  char* dup_string(const char*);
  ip_keyword_tree_t* ip_get_variable_kt(char*);
  char* ip_get_variable_value(char*);
  void ip_internal_values();
  void ip_push_keyword(char*);
  void ip_push_keyclass(char*,char*,ip_string_list_t*);
  void ip_pop_keyword();
  void ip_begin_table(ip_string_list_t*);
  void ip_done_table();
  ip_string_list_t* ip_add_string_list(ip_string_list_t*,char*);
  ip_string_list_t* ip_string_to_string_list(char*);
  void ip_assign_variable(char*);
  double ip_get_variable_double(char*);
  char* ip_double_to_string(double);
  void ip_assign_value(char*value);
  void ip_start_karray();
  void ip_init_karray();
  void ip_incr_karray();
  void ip_lastkeyword(const char*);
  void ip_lastkeywordtree(ip_keyword_tree_t*);
  void ip_lastkeyword_(ip_keyword_tree_t*);
  ip_keyword_tree_t* ip_alloc_keyword_tree();
  void ip_free_keyword_tree(ip_keyword_tree_t*);
  void ip_cwk_add_kt(ip_keyword_tree_t*);
  ip_keyword_tree_t* ip_cwk_descend_tree(const char*);
  ip_keyword_tree_t* ip_descend_tree(ip_keyword_tree_t*,const char*);
  char* ip_key_value(const char*);
  void free_keyword_tree_list(ip_keyword_tree_list_t*);
  ip_keyword_tree_list_t* splice_keyword_tree_list(ip_keyword_tree_t*,
                                                   ip_keyword_tree_list_t*);
  void ip_cwk_karray_add_v(int,int*);
  void ip_cwk_karray_add(int,...);
  ip_keyword_tree_t* ip_karray_descend_v(ip_keyword_tree_t*,int,int*);
  ip_keyword_tree_t* ip_karray_descend(ip_keyword_tree_t*,int,...);
  void ip_print_keyword(ostream&,ip_keyword_tree_t*);
  void ip_print_tree(ostream&,ip_keyword_tree_t*);
  void ip_print_tree_(ostream&,ip_keyword_tree_t*,int);
  void ip_indent(ostream&,int);
  int ip_special_characters(char*);
  char* ip_append_keystrings(char*,char*);
  void ip_pop_karray();
  void ip_initialize(istream&,ostream&);
  void ip_append(istream&,ostream&);
  char* get_truename(ip_keyword_tree_t*kt);

  void showpos();

  IPV2FlexLexer *lexer;

  int ylex() { return lexer->yylex(); }
  int yparse();
  void yerror(const char* s);

 public:
  IPV2();
  virtual ~IPV2();
  static int have_global();
  static void set_global(IPV2*);
  static IPV2* global();
  // calls either ip_append or ip_initialize based on ip_initialized
  void read(istream&,ostream&,const char *filename=0);
  void append_from_input(const char*,ostream&);
  void done();
  const char* error_message(IPV2::Status);
  void error(const char*);
  void warn(const char*);
  void cwk_root();
  void cwk_clear();
  void cwk_add(const char*);
  void cwk_push();
  void cwk_pop();
  IPV2::Status boolean(const char*,int*,int,...);
  IPV2::Status boolean_v(const char*,int*,int,int*);
  int exist(const char*,int,...);
  int exist_v(const char*,int,int*);
  IPV2::Status data(const char*,const char*,void*,int,...);
  IPV2::Status data_v(const char*,const char*,void*,int,int*);
    // the character string produced by classname must not be delete[]'ed
  IPV2::Status classname(const char*,const char**,int,...);
  IPV2::Status classname_v(const char*,const char**,int,int*);
    // the character string produced by truekeyword must not be delete[]'ed
    // if there is no alias for the keyword the string pointer is set to
    // null and if the keyword exists OK is returned
  IPV2::Status truekeyword(const char*,const char**,int,...);
  IPV2::Status truekeyword_v(const char*,const char**,int,int*);
  IPV2::Status string(const char*,char**,int,...);
  IPV2::Status string_v(const char*,char**,int,int*);
    // the character string produced by value must not be delete[]'ed
    // or free'ed.
  IPV2::Status value(const char*,const char**,int,...);
  IPV2::Status value_v(const char*,const char**,int,int*);

  IPV2::Status construct_key_v(const char*,char*,int,int*);
  IPV2::Status count(const char*,int*,int,...);
  IPV2::Status count_v(const char*,int*,int,int*);

  // some routines for debugging
  void print_keyword(ostream&f=cerr,ip_keyword_tree_t*k=0);
  void print_tree(ostream&f=cerr,ip_keyword_tree_t*k=0);
};

#endif
