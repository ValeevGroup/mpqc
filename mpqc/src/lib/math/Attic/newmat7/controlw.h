//$$ controlw.h                Control word class

#ifndef CONTROL_WORD_LIB
#define CONTROL_WORD_LIB 0

// for organising an int as a series of bits which indicate whether an
// option is on or off.

class ControlWord
{
protected:
   int cw;                                      // the control word
public:
   ControlWord();                     // do nothing
   ControlWord(int i);                // load an integer

      // select specific bits (for testing at least one set)
   ControlWord operator*(ControlWord i) const;
   void operator*=(ControlWord i);

      // set bits
   ControlWord operator+(ControlWord i) const;
   void operator+=(ControlWord i);

      // reset bits
   ControlWord operator-(ControlWord i) const;
   void operator-=(ControlWord i);

      // check if all of selected bits set or reset
   Boolean operator>=(ControlWord i) const;
   Boolean operator<=(ControlWord i) const;

      // flip selected bits
   ControlWord operator^(ControlWord i) const;
   ControlWord operator~() const;

      // convert to integer
   int operator+() const;
   int operator!() const;
   FREE_CHECK(ControlWord)
};


#endif
