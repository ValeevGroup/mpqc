
/* $Log$
 * Revision 1.1  1993/12/29 12:53:57  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  18:11:18  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:11:16  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/20  16:19:40  seidl
 * Initial revision
 *
 * Revision 1.2  1991/07/19  14:42:33  cljanss
 * Changed the way that generation modules are selected.
 *
 * Revision 1.1  1991/06/15  21:15:38  janssen
 * Initial revision
 * */

#define STRING_LENGTH 512

#define TYPE_CONSTANT 1
#define TYPE_VARIABLE 2
struct index_struct {
  int type;
  union {
    int c;
    char *v;
    } v;
  };

struct index_list_struct {
  struct index_struct index;
  struct index_list_struct *p;
  };

  /* These are the possible values for the qualifier. */
#define Q_NONE 1
#define Q_STRUCT 2
#define Q_SIGNED 3
#define Q_UNSIGNED 4
#define Q_UNION 5
struct member_struct {
  int qualifier;
  char *type;
  int pointer;  /* The number of levels of pointers (number of '*'s). */
  char *name;
  struct index_list_struct *indices;
  char *uname;    /* union name */
  char *uselname; /* union select name */
  char *uselval;  /* union select val */
  int mark; /* A generic int. */
  };

struct member_list_struct {
  struct member_struct *member;
  struct member_list_struct *p;
  };

struct name_list_struct {
  char *name;
  struct name_list_struct *p;
  };

struct declaration_struct {
  char *name;
  struct member_list_struct *members;
  struct name_list_struct *modules;
  struct declaration_struct *p;
  };

typedef struct index_struct index_t;
typedef struct index_list_struct index_list_t;
typedef struct member_struct member_t;
typedef struct member_list_struct member_list_t;
typedef struct declaration_struct declaration_t;
typedef struct name_list_struct name_list_t;

