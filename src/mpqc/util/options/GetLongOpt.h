
#ifndef SRC_MPQC_UTIL_OPTIONS_GETLONGOPT_H_
#define SRC_MPQC_UTIL_OPTIONS_GETLONGOPT_H_

#include <string>
#include <iostream>
#include <list>

namespace mpqc {

/// @addtogroup Init
/// @{

/// Parse command line options.
class GetLongOpt {
 public:
  /// Used by the enroll member to specify whether or not a value is
  /// expected or optional.
  enum OptType { Invalid, NoValue, OptionalValue, MandatoryValue };

 private:
  struct Cell {
    std::string option;       // option name
    OptType type;             // option type
    std::string description;  // a description of option
    std::string value;        // value of option (string)
  };

 private:
  std::list<Cell> table;          // option table
  std::string ustring;  // usage message
  std::string pname;          // program basename
  char optmarker;       // option marker

  bool finalized;  // finished enrolling

 private:
  std::string basename(const std::string& p) const;
  int setcell(Cell& c, const char *valtoken, const char *nexttoken,
              const std::string &name);

 public:
  /// Initialize the object.
  /// @param optmark the option flag marker (default is <tt>-</tt>).
  GetLongOpt(const char optmark = '-');
  ~GetLongOpt() = default;

  /// Parse command line options.
  /// @param argc the number of arguments, as passed to <tt>main</tt>
  /// @param argv the arguments, as passed to <tt>main</tt>
  /// @return the index to the start of arguments that were not
  ///         processed (an error occurred if the return value is < 1)
  int parse(int argc, char *const *argv);

  /// Enroll an option.
  /// @param opt the option name
  /// @param t whether or not a value is expected
  /// @param desc a description of the option
  /// @param val a default value for the option with an optional value
  int enroll(const std::string &opt, const OptType t,
             const std::string &desc, const std::string &val);
  /// Retrieve an option.
  /// @param opt the name of the option
  std::string retrieve(const std::string& opt) const;

  /// Print usage information.
  /// @param outfile stream to use for printing (default: <tt>std::cout</tt>)
  void usage(std::ostream &outfile = std::cout) const;
  /// Initialize usage synopsis.
  /// @param str the usage synopsis
  void usage(std::string&& str) { ustring = std::move(str); }
};

/// @}
// end of addtogroup Init
}

#endif  // SRC_MPQC_UTIL_OPTIONS_GETLONGOPT_H_
