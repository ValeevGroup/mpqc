
#ifndef SRC_BIN_MPQC_MPQCINIT_H_
#define SRC_BIN_MPQC_MPQCINIT_H_

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
/// \param[in/out] opt command-line options parser, if non-null
///                additional MPQC-specific options will be enrolled
void initialize(int &argc, char **argv, std::shared_ptr<GetLongOpt> opt = nullptr);

/// Finalize MPQC.
void finalize();

/** \brief This helper singleton class simplifies initialization of MPQC.
 *
 * An object of this type is created on every process of the default MADWordl World
 * by calling mpqc::initialize() and can
 * be accessed for additional customization via mpqc::MPQCInit::instance(). The object
 * is destroyed by calling mpqc::finalize() .
 * \note this object is meant to be manipulated by 1 (usually, main) thread.
 */
class MPQCInit {
  std::shared_ptr<GetLongOpt> opt_;
  char **argv_;
  int &argc_;

  /// the unique instance
  static std::unique_ptr<MPQCInit> instance_;

 public:
  ~MPQCInit();

  /// \return the reference to the only instance of this object
  /// \throw mpqc::ProgrammingError if mpqc::initialize() had not been called
  static MPQCInit& instance();

  /// \return the command-line argument count
  const int& argc() const { return argc_; }
  /// \return the command-line arguments
  char** argv() const { return argv_; }
  /// \return the command-line options parser
  std::shared_ptr<const GetLongOpt> opt() const {
    return std::const_pointer_cast<const GetLongOpt, GetLongOpt>(opt_);
  }

  /// Creates the KeyVal object from the contents of file \c filename .
  /// The file will be read by one of the processes in \c world and broadcast to
  /// every other process.
  /// \note Must be called on every process in \c world
  std::shared_ptr<mpqc::KeyVal> make_keyval(madness::World& world,
                                            const std::string &filename);

  /// Set the name used to construct data file names. A noncollective operation.
  void set_basename(const std::string &input_filename,
                    const std::string &output_filename = "");

 private:

  friend void ::mpqc::initialize(int &argc, char **argv, std::shared_ptr<GetLongOpt> opt);
  friend void ::mpqc::finalize();

  /// Create the initializer. Only one object of this time can be created.
  /// Needed options will be enrolled
  /// in the opt object. The parse member of opt must be called
  /// after this constructor completes, but before any of the
  /// other members of MPQCInit are called.
  ///
  /// N.B. This is not explicitly implemented as a Singleton for syntactic
  /// reasons.
  MPQCInit(int &argc, char **argv, std::shared_ptr<GetLongOpt> opt);

  /// Initialize the floating point control word.
  void init_fp();
  /// Initialize the resource limits.
  void init_limits();
  /// Initialize the default ConsumableResources object.
  void init_resources(std::shared_ptr<mpqc::KeyVal> keyval);

  /// Initialize formatted I/O.
  void init_io();
  /// Calls all of the initialize routines in the proper sequence.
  /// The parse member for the GetLongOpt object given to the
  /// constructor must have been called before this is called.
  void init(const std::string &input_filename,
            const std::string &output_filename = "");
};

/// @}
// end of addtogroup Init
}  // namespace mpqc

#endif
