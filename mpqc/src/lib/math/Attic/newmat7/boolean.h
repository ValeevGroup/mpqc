//$$ boolean.h                       Boolean class

#ifndef Boolean_LIB
#define Boolean_LIB 0

class Boolean
{
   int value;
public:
   Boolean(const int b);
   Boolean(const void* b);
   Boolean();
   operator int() const;
   int operator!() const;
   FREE_CHECK(Boolean);
};

#ifndef __GNUG__
   const Boolean TRUE = 1;
   const Boolean FALSE = 0;
#else
#define FALSE 0
#define TRUE 1
#endif

#endif
