
#ifndef MPQC4_SRC_BIN_MPQC_MPQC_INIT_H_
#define MPQC4_SRC_BIN_MPQC_MPQC_INIT_H_

#include <cstdlib>
#include <string>

#include <madness/world/world.h>

#include "mpqc/mpqc_config.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/options/GetLongOpt.h"

namespace mpqc {

/// @addtogroup Init
/// @{

/// \brief Static MPQC initializer. Must be called from the main thread of every process
///        in the default MADWorld World before doing any MPQC-specific computation
///        (e.g. before creating mpqc::MPQC objects).
///
/// madness::initialize() must have been called before this.
/// \throw mpqc::ProgrammingError if the MADWorld runtime is not initialized.
/// \throw mpqc::ProgrammingError if already been called and mpqc::finalize() has not been
///        called since the previous call to mpqc::initialize()
/// \note The initialize()/finalize() sequence can occur more than once.
/// \param[in/out] argc the argument count
/// \param[in/out] argv the sequence of argument strings
/// \param[in] world the top World object in which MPQC will execute (MPQCTask object will
///                  live in subworlds of this)
/// \param[in/out] opt command-line options parser, if non-null
///                additional MPQC-specific options will be enrolled
void initialize(int &argc, char **argv,
                const madness::World& top_world,
    std::shared_ptr<GetLongOpt> opt = std::shared_ptr<GetLongOpt>());

/// Finalize MPQC.
void finalize();

/// @brief Constructs a KeyVal object on every rank of @c world by reading file @c filename on rank 0.

/// Will try every known file format from which KeyVal can be constructed (see KeyVal::InputFormat ).
/// @note This is a collective operation.
/// @param[in] world the World object
/// @param[in] filename the file name
/// @return a tuple consisting of a pointer to the KeyVal object and the input format identifier
std::tuple<std::shared_ptr<mpqc::KeyVal>, KeyVal::InputFormat>
make_keyval(madness::World &world, const std::string &filename);

/** \brief This helper singleton class simplifies initialization of MPQC.
 *
 * An object of this type is created on every process of the default MADWordl World
 * by calling mpqc::initialize() and can
 * be accessed for additional customization via mpqc::MPQCInit::instance(). The object
 * is destroyed by calling mpqc::finalize() .
 * \note this object is meant to be manipulated by 1 (usually, main) thread.
 */
class MPQCInit {
 private:
  struct singleton_ctor_tag {};

 public:
  using InputFormat = KeyVal::InputFormat;

  ~MPQCInit();

  /// \return the reference to the only instance of this object
  /// \throw mpqc::ProgrammingError if mpqc::initialize() had not been called
  static MPQCInit &instance();

  /// \return the command-line argument count
  const int &argc() const { return argc_; }
  /// \return the command-line arguments
  char **argv() const { return argv_; }
  /// \return the command-line options parser
  std::shared_ptr<const GetLongOpt> opt() const {
    return std::const_pointer_cast<const GetLongOpt, GetLongOpt>(opt_);
  }
  InputFormat input_format() const { return input_format_; }

  /// Creates the KeyVal object from the contents of file \c filename .
  /// The file will be read by one of the processes in \c world and broadcast to
  /// every other process.
  /// \note Must be called on every process in \c world
  std::shared_ptr<mpqc::KeyVal> make_keyval(madness::World &world,
                                            const std::string &filename);

  /// Set the name used to construct data file names. A noncollective operation.
  void set_basename(const std::string &input_filename,
                    const std::string &output_filename = "");

  /// Create the initializer. Only one object of this time can be created.
  /// Needed options will be enrolled
  /// in the opt object. The parse member of opt must be called
  /// after this constructor completes, but before any of the
  /// other members of MPQCInit are called.
  ///
  /// N.B. This is not explicitly implemented as a Singleton for syntactic
  /// reasons.
  ///
  /// \param[in] world the top World object in which MPQC will execute
  MPQCInit(int &argc, char **argv, std::shared_ptr<GetLongOpt> opt,
           const madness::World &world, singleton_ctor_tag);

 private:
  std::shared_ptr<GetLongOpt> opt_;
  char **argv_;
  int &argc_;
  InputFormat input_format_;

  /// the unique instance
  static std::unique_ptr<MPQCInit> instance_;

  friend void ::mpqc::initialize(int &argc, char **argv,
                                 const madness::World &world,
                                 std::shared_ptr<GetLongOpt> opt);
  friend void ::mpqc::finalize();

  /// Initialize the floating point control word.
  void init_fp();
  /// Initialize the resource limits.
  void init_limits();
  /// Initialize the default ConsumableResources object.
  //  void init_resources(std::shared_ptr<mpqc::KeyVal> keyval);

  /// Initialize the path to directory used for POSIX I/O of large text/binary files.
  /// Input read from the environment variable MPQC_WORK_DIR . The default is the current working
  /// directory of this process.
  /// @note The existence of this path will be checked in every rank of the World object used to initialize MPQCInit .
  void init_work_dir();

  /// Initialize formatted I/O. This is a collective operation.
  /// \param[in] world the top World object in which MPQC will execute
  void init_io(const madness::World &top_world);
  /// Calls all of the initialize routines in the proper sequence.
  /// The parse member for the GetLongOpt object given to the
  /// constructor must have been called before this is called.
  void init(const std::string &input_filename,
            const std::string &output_filename = "");
};

inline std::string
to_string(MPQCInit::InputFormat f) {
  switch (f) {
    case MPQCInit::InputFormat::json:
      return "JSON";
    case MPQCInit::InputFormat::xml:
      return "XML";
    case MPQCInit::InputFormat::info:
      return "INFO";
    default:
      return "invalid";
  }
}

// clang-format off
/** \brief Creates a default options parser object for an MPQC executable

    Creates a new options parser object and enrolls standard MPQC options:
    | option          | accept value?| description                                           |
    |-----------------|--------------|-------------------------------------------------------|
    | <tt>-i</tt>     | mandatory    | The name of the input file. MPQC will attempt to parse the given input file using the JSON, XML, and INFO formats (in that order). If <tt>-i</tt> is not given, and no options with omitted optional values are used (e.g., <tt>-D</tt>), the last command-line argument that does not specify an option will be assumed to give the input file name.|
    | <tt>-o</tt>     | mandatory    | The name of the output file.  The default is to send output to the console. |
    | <tt>-p</tt>     | mandatory    | The prefix for all relative file paths in the input file               |
    | <tt>-W</tt>     | mandatory    | the working directory in which to compute             |
    | <tt>-D</tt>     | optional     | unless "debugger" keyword is given in input KeyVal, create a debugger at start, with the optional argument as its JSON KeyVal input |
    | <tt>-v</tt>     | no           | print the version number and exit                     |
    | <tt>-w</tt>     | no           | print the warranty and exit                           |
    | <tt>-L</tt>     | no           | print the license and exit                            |
    | <tt>-k</tt>     | no           | print all registered (KeyVal-constructible) DescribedClass classes |
    | <tt>-h</tt>     | no           | print the usage info and exit                         |
    | <tt>-d</tt>     | no           | start the program and attach a debugger               |
    | <tt>-t</tt>     | no           | throw if a deprecated keyword is read                 |

    \return a smart pointer to the newly created options parser object (see mpqc::GetLongOpt )
    */
// clang-format on
std::shared_ptr<GetLongOpt> make_options();

/// \brief Processes command-line options parsed by the options parser.
///
/// \param options the options parser object whose \c parse method has already
/// been called.
/// \return the {input,output} file name tuple.
/// \note will exit via \c std::exit(0) if any of these options are given: -h, -v, -w, or -L
std::tuple<std::string, std::string>
process_options(const std::shared_ptr<GetLongOpt>& options);

/// @}
// end of addtogroup Init
}  // namespace mpqc

#endif  // MPQC4_SRC_BIN_MPQC_MPQC_INIT_H_
