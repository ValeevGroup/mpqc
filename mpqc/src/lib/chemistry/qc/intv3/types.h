
#ifndef _chemistry_qc_intv3_types_h
#define _chemistry_qc_intv3_types_h

/* Types that are used for integrals, but for which we don't need all
 * of the sgen utilities, are defined here. */

struct struct_der_centers {
  int n;
  centers_t *cs[4];
  int num[4];
  centers_t *ocs; /* The omitted center's centers_t. */
  int onum;        /* The omitted center's number. */
  };

typedef struct struct_der_centers der_centers_t;

#endif
