#define NO_MEM_DEBUG

char *mem_malloc();
void mem_free();
#ifdef MEM_DEBUG

#define Malloc(x) mem_malloc((int)(x), __FILE__, __LINE__)
#define Free(x) mem_free((char *)(x))

#else

#ifdef NCUBE
#define Malloc(x) (((x) > 0) ? malloc((unsigned) (x)) : malloc(sizeof(double)))
#else
#define Malloc(x) malloc((unsigned) (x))
#endif

#define Free(x) free((char *)(x))

#define mem_dump()
#define mem_set_max(x) 
#define mem_used() 0

#endif
