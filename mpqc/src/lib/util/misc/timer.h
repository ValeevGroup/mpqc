
#include <util/misc/timer.gbl>

#ifdef TIMING

#define TENTER(x) tim_enter(x)
#define TEXIT(x) tim_exit(x)
#define TCHANGE(x) tim_change(x)
#define TPRINT(x) tim_print(x)
#else
#ifdef PRINTING
#define TENTER(x) fprintf(stdout,"TENTER: %s\n",x);
#define TEXIT(x) fprintf(stdout,"TEXIT: %s\n",x);
#define TCHANGE(x) fprintf(stdout,"TCHANGE: %s\n",x);
#define TPRINT(x) fprintf(stdout,"TPRINT\n");
#else
#define TENTER(x)
#define TEXIT(x) 
#define TCHANGE(x)
#define TPRINT(x)
#endif
#endif
