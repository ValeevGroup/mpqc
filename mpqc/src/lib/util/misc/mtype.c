
#define MTYPE_START 10000
#define MTYPE_END   16384 /* = 0x4000, this constraint is due to the nCUBE */

static int next_mtype = MTYPE_START;


int
mtype_get()
{
  int r = next_mtype;

  if (next_mtype==MTYPE_END) next_mtype = MTYPE_START;
  else next_mtype++;

  return r;
  }

