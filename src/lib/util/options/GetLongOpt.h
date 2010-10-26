/* $Id$ */
/* S Manoharan. Advanced Computer Research Institute. Lyon. France */

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _GetLongOpt_h_
#define _GetLongOpt_h_

#include <iostream>
#include <string.h>

namespace sc {

/// Parse command line options.
class GetLongOpt {
public:
  /// Used by the enroll member to specify whether or not a value is
  /// expected or optional.
   enum OptType { 
      NoValue, OptionalValue, MandatoryValue
   };
private:
   struct Cell {
      const char *option;	// option name
      OptType type;		// option type
      const char *description;	// a description of option
      const char *value;	// value of option (string)
      Cell *next;		// pointer to the next cell

      Cell() { option = description = value = 0; next = 0; }
   };
private:
  Cell *table;				// option table
  const char *ustring;			// usage message
  char *pname;				// program basename
  char optmarker;			// option marker

  int enroll_done;			// finished enrolling
  Cell *last;				// last entry in option table 

private:
  char *basename(char * const p) const;
  int setcell(Cell *c, const char *valtoken, const char *nexttoken, const char *p);
public:
   /// Initialize the object.
   /// @param optmark the option flag marker (default is <tt>-</tt>).
   GetLongOpt(const char optmark = '-');
   ~GetLongOpt();

   /// Parse command line options.
   /// @param argc the number of arguments, as passed to <tt>main</tt>
   /// @param argv the arguments, as passed to <tt>main</tt>
   /// @return the index to the start of arguments that were not
   ///         processed (an error occurred if the return value is < 1)
   int parse(int argc, char * const *argv);
   /// Parse options in a string.
   /// @param str the string to be parsed
   /// @param p a prefix that will be prefixed to error messages
   /// @return the index to the start of arguments that were not
   ///         processed (an error occurred if the return value is < 1)
   int parse(char * const str, char * const p);

   /// Enroll an option.
   /// @param opt the option name
   /// @param t whether or not a value is expected
   /// @param desc a description of the option
   /// @param val a default value for the option with an optional value
   int enroll(const char * const opt, const OptType t,
      const char * const desc, const char * const val);
   /// Retrieve an option.
   /// @param opt the name of the option
   const char *retrieve(const char * const opt) const;

   /// Print usage information.
   /// @param outfile stream to use for printing (default: <tt>std::cout</tt>)
   void usage(std::ostream &outfile = std::cout) const;
   /// Initialize usage synopsis.
   /// @param str the usage synopsis
   void usage(const char *str)		{ ustring = str; }
};

}

#endif /* _GetLongOpt_h_ */
