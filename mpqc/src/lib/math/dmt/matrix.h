
typedef int dmt_matrix;
typedef int mdescr;
typedef int bool;

#define MAXARGS 50

struct ngl_descr {
  int type;
  mdescr m;
  double *a;
  double *t;
  int readonly;
  };
typedef struct ngl_descr ngl_descr_t;

struct ngl_list {
  ngl_descr_t args[MAXARGS];
  int nargs;
  int rmsgtype;
  int smsgtype;
  void *buf;
  int current_bufsize;
  int total_bufsize;
  int loop_count;
  int readonly;
  int nextblock;
  int innerloopmat;
  struct ngl_list *p;
  };
typedef struct ngl_list loop_t;

#define COLUMNS   1
#define SCATTERED 2

struct dmt_block_info {
  int i;
  int j;
  int magnitude; /* the value of Qvec for this block */
  int ami;
  int amj;
  int nconi;
  int nconj;
  int nprimi;
  int nprimj;
  int dimi;
  int dimj;
};
typedef struct dmt_block_info dmt_cost_t;

