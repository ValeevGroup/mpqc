
/* $Id$ */
/* $Log$
 * Revision 1.1  1993/12/29 12:53:41  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/04/06  12:28:40  seidl
 * merge in sandia changes
 *
 * Revision 1.2  1992/03/30  22:45:53  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  91/06/16  20:28:50  seidl
 * Initial revision
 *  */

typedef
struct {
  int junk;
  } r_async_t;

typedef
struct {
  int junk;
  } s_async_t;

typedef
struct {
  char *path;
  int stream;
  } sequential_volume_t;

typedef
struct {
  int n;
  int blocksize;
  int last_ioop;
  int verbose;
  long next;
  int unit; /* The unit number. */
  unsigned incount;  /* The number of bytes read. */
  unsigned outcount;  /* The number of bytes written. */
  sequential_volume_t v[MAX_VOLUME];
  } sequential_t;

typedef
struct {
  int junk;
  } ram_t;
  

typedef
struct {
  int method;
  int status;
  union {
    r_async_t *r_async;
    s_async_t *s_async;
    sequential_t *sequential;
    ram_t *ram;
    } ptr;
  } ioFILE_t;
