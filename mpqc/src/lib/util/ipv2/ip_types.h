
/* This is a tree with all of the keywords and values. */
struct ip_keyword_tree_struct {
  char *keyword;
  char *classname;
  struct ip_keyword_tree_struct *across; /* Circular list. */
  struct ip_keyword_tree_struct *up;    /* Terminated by NULL. */
  struct ip_keyword_tree_struct *down;  /* Terminated by NULL. */
  char *value;
  };

struct ip_keyword_tree_list_struct {
  struct ip_keyword_tree_struct *kt;
  struct ip_keyword_tree_list_struct *p;
  };

struct ip_string_list_struct {
  char *string;
  struct ip_string_list_struct *p;
  };

typedef struct ip_keyword_tree_struct ip_keyword_tree_t;
typedef struct ip_keyword_tree_list_struct ip_keyword_tree_list_t;
typedef struct ip_string_list_struct ip_string_list_t;

