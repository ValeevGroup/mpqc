
#ifndef MPQC4_SRC_MPQC_UTIL_OPTIONS_GETLONGOPT_H_
#define MPQC4_SRC_MPQC_UTIL_OPTIONS_GETLONGOPT_H_

#include <iostream>
#include <list>
#include <memory>
#include <string>

#include <boost/optional.hpp>

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
    std::string option;                  // option name
    OptType type;                        // option type
    std::string description;             // a description of option
    std::unique_ptr<std::string> value;  // ptr to value (nullptr, if not given)

    Cell() = default;
    Cell(const Cell&) = default;
    Cell(Cell&& other)
        : option(std::move(other.option)),
          type(other.type),
          description(std::move(other.description)),
          value(std::move(other.value)) {}
    Cell(std::string opt, OptType t, std::string descr,
         std::unique_ptr<std::string> val)
        : option(std::move(opt)),
          type(t),
          description(std::move(descr)),
          value(std::move(val)) {}
    ~Cell() { }
  };

 private:
  std::list<Cell> table;  // option table
  std::string ustring;    // usage message
  std::string pname;      // program basename
  char optmarker;         // option marker

  bool finalized;  // finished enrolling
  int first_unprocessed_arg_;  //!< index of the first argument that was not processed

 private:
  std::string basename(const std::string& p) const;
  int setcell(Cell& c, const char* valtoken, const char* nexttoken,
              const std::string& name);

 public:
  /// Initialize the object.
  /// @param optmark the option flag marker (default is <tt>-</tt>).
  GetLongOpt(const char optmark = '-');
  ~GetLongOpt() { }

  /// Parse command line options.
  /// @note call this once, after all options have been enrolled
  /// @warning this object becomes finalized, additional options cannot be
  /// enrolled
  /// @param argc the number of arguments, as passed to <tt>main</tt>
  /// @param argv the arguments, as passed to <tt>main</tt>
  void parse(int argc, char* const* argv);
  /// Parse options in a string.
  /// @note call this once, after all options have been enrolled
  /// @warning this object becomes finalized, additional options cannot be
  /// enrolled
  /// @param str the string to be parsed
  /// @param p a prefix that will be prefixed to error messages
  void parse(const std::string& str, const std::string& p);

  /// After calling parse() this will return the index of the first argument that was not processed
  /// @return the index to the start of arguments that were not
  ///         processed (an error occurred if the return value is < 1)
  int first_unprocessed_arg() const {
    return first_unprocessed_arg_;
  }

  /// Enroll an option.
  /// @param opt the option name
  /// @param t whether or not a value is expected
  /// @param desc a description of the option
  int enroll(std::string opt, const OptType t, std::string desc);

  /// Enroll an option, with the default value provided.
  /// @param opt the option name
  /// @param t whether or not a value is expected
  /// @param desc a description of the option
  /// @param val the default value for the option with an optional value
  int enroll(std::string opt, const OptType t, std::string desc,
             std::string default_value);

  /// Retrieve the value of the option.
  /// @param opt the name of the option
  /// @return if \c opt was given, return \c boost::optional<std::string>
  ///         initialized with the value of the option (empty string for \c NoValue option types),
  ///         otherwise a default-initialized \c boost::optional<std::string>
  boost::optional<std::string> retrieve(const std::string& opt) const;

  /// Print usage information.
  /// @param outfile stream to use for printing (default: <tt>std::cout</tt>)
  void usage(std::ostream& outfile = std::cout) const;
  /// Initialize usage synopsis.
  /// @param str the usage synopsis
  void usage(std::string&& str) { ustring = std::move(str); }
};

/// @}
// end of addtogroup Init
}

#endif  // MPQC4_SRC_MPQC_UTIL_OPTIONS_GETLONGOPT_H_
