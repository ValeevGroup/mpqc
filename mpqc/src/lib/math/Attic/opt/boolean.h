//$$ boolean.h                       Boolean class

#ifndef Boolean_LIB
#define Boolean_LIB 0

class Boolean
{
   int value;
public:
   Boolean(const int b) { value = b ? 1 : 0; }
   Boolean(const void* b) { value = b ? 1 : 0; }
   Boolean() {}
   operator int() const { return value; }
   int operator!() const { return !value; }
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
