//$$ precisio.h                          floating point constants

#ifndef PRECISION_LIB
#define PRECISION_LIB 0

#ifndef SystemV                    // if there is float.h


#ifdef USING_FLOAT


class FloatingPointPrecision
{
public:
   static int Dig();              // number of decimal digits or precision
   static Real Epsilon();          // smallest number such that 1+Eps!=Eps
   static int Mantissa();         // bits in mantisa
   static Real Maximum();              // maximum value
   static int MaximumDecimalExponent();       // maximum decimal exponent
   static int MaximumExponent();          // maximum binary exponent
   static Real Minimum();              // minimum positive value
   static int MinimumDecimalExponent();       // minimum decimal exponent
   static int MinimumExponent();          // minimum binary exponent
   static int Radix();            // exponent radix
   static int Rounds();           // addition rounding (1 = does round)
   FREE_CHECK(FloatingPointPrecision)
};

#endif


#ifdef USING_DOUBLE

class FloatingPointPrecision
{
public:
   static int Dig();
   static Real Epsilon();
   static int Mantissa();
   static Real Maximum();
   static int MaximumDecimalExponent();
   static int MaximumExponent();
   static Real Minimum();
   static int MinimumDecimalExponent();
   static int MinimumExponent();
   static int Radix();
   static int Rounds();
   FREE_CHECK(FloatingPointPrecision)
};

#endif

#endif

#ifdef SystemV                    // if there is no float.h

#ifdef USING_FLOAT

class FloatingPointPrecision
{
public:
   static Real Epsilon();
   static Real Maximum();
   static Real Minimum();
   FREE_CHECK(FloatingPointPrecision)
};

#endif


#ifdef USING_DOUBLE

class FloatingPointPrecision
{
public:
   static Real Epsilon();
   static Real Maximum();
   static Real Minimum();
   FREE_CHECK(FloatingPointPrecision)
};

#endif

#endif




#endif
