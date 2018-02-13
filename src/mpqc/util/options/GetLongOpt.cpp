#include "mpqc/util/options/GetLongOpt.h"

#include <cstring>
#include <memory>

using namespace std;
using namespace mpqc;

GetLongOpt::GetLongOpt(const char optmark)
    : ustring("[valid options and arguments]"),
      optmarker(optmark),
      finalized(false), first_unprocessed_arg_(0) {}

std::string GetLongOpt::basename(const std::string &pname) const {
  auto pos = pname.find_last_of('/');
  if (pos == std::string::npos)
    return pname;
  else
    return pname.substr(pos + 1);
}

int GetLongOpt::enroll(std::string opt, const OptType t,
                       std::string desc) {
  if (finalized) return 0;
  table.emplace_back(Cell{
      opt, t, (!desc.empty() ? desc : std::string("no description available")),
      nullptr});
  return 1;
}

int GetLongOpt::enroll(std::string opt, const OptType t,
                       std::string desc, std::string val) {
  if (finalized) return 0;
  table.emplace_back(Cell{
      opt, t, (!desc.empty() ? desc : std::string("no description available")),
      std::make_unique<std::string>(std::move(val))});
  return 1;
}

boost::optional<std::string> GetLongOpt::retrieve(
    const std::string &opt) const {
  using result_t = boost::optional<std::string>;
  if (finalized == true) {
    for (const auto &t : table) {
      if (opt == t.option)
        return t.value ? result_t(std::string(*t.value)) : result_t();
    }
    std::cerr << "GetLongOpt::retrieve - unenrolled option ";
    std::cerr << optmarker << opt << "\n";
  }
  return result_t();
}

void GetLongOpt::parse(int argc, char *const *argv) {
  int optind = 1;

  pname = basename(*argv);
  finalized = true;
  if (argc-- <= 1) { first_unprocessed_arg_ = optind; return; }

  while (argc >= 1) {
    char *token = *++argv;
    --argc;

    if (token[0] != optmarker || token[1] == optmarker)
      break; /* end of options */

    ++optind;
    char *tmptoken = ++token;
    while (*tmptoken && *tmptoken != '=') ++tmptoken;
    /* (tmptoken - token) is now equal to the command line option
       length. */

    enum { NoMatch, ExactMatch, PartialMatch } matchStatus = NoMatch;
    Cell *pc = 0;  // pointer to the partially-matched cell
    for (auto &t : table) {
      if (strncmp(t.option.c_str(), token, (tmptoken - token)) == 0) {
        if (t.option.size() == (tmptoken - token)) {
          /* an exact match found */
          int stat = setcell(t, tmptoken, *(argv + 1), pname);
          if (stat == 1) {
            ++argv;
            --argc;
            ++optind;
          }
          matchStatus = ExactMatch;
          break;
        } else {
          /* partial match found */
          matchStatus = PartialMatch;
          pc = &t;
        }
      } /* end if */
    }   /* end for */

    if (matchStatus == PartialMatch) {
      int stat = setcell(*pc, tmptoken, *(argv + 1), pname);
      if (stat == 1) {
        ++argv;
        --argc;
        ++optind;
      }
    } else if (matchStatus == NoMatch) {
      cerr << pname << ": unrecognized option ";
      cerr << optmarker << strtok(token, "= ") << "\n";
      throw std::runtime_error("unrecognized option " + std::string(1, optmarker) +
                               std::string(strtok(token, "= ")));
    }

  } /* end while */

 first_unprocessed_arg_ = optind;
}

void GetLongOpt::parse(const std::string &cppstr, const std::string &p) {
  finalized = true;
  std::unique_ptr<char[]> str_ptr(strdup(cppstr.c_str()));
  char* str = str_ptr.get();
  char *token = strtok(str, " \t");
  const char *name = p.c_str();

  while (token) {
    if (token[0] != optmarker || token[1] == optmarker) {
      cerr << name << ": nonoptions not allowed\n";
      throw std::runtime_error("nonoptions not allowed");
    }

    char *ladtoken = 0; /* lookahead token */
    char *tmptoken = ++token;
    while (*tmptoken && *tmptoken != '=') ++tmptoken;
    /* (tmptoken - token) is now equal to the command line option
       length. */

    enum { NoMatch, ExactMatch, PartialMatch } matchStatus = NoMatch;
    Cell *pc = 0;  // pointer to the partially-matched cell
    for (auto &t : table) {
      if (strncmp(t.option.c_str(), token, (tmptoken - token)) == 0) {
        if (t.option.size() == (tmptoken - token)) {
          /* an exact match found */
          ladtoken = strtok(0, " \t");
          int stat = setcell(t, tmptoken, ladtoken, name);
          if (stat == 1) {
            ladtoken = 0;
          }
          matchStatus = ExactMatch;
          break;
        } else {
          /* partial match found */
          matchStatus = PartialMatch;
          pc = &t;
        }
      } /* end if */
    }   /* end for */

    if (matchStatus == PartialMatch) {
      ladtoken = strtok(0, " \t");
      int stat = setcell(*pc, tmptoken, ladtoken, name);
      if (stat == 1) {
        ladtoken = 0;
      }
    } else if (matchStatus == NoMatch) {
      cerr << name << ": unrecognized option ";
      cerr << optmarker << strtok(token, "= ") << "\n";
      throw std::runtime_error("unrecognized option " + std::string(1, optmarker) +
                               std::string(strtok(token, "= ")));
    }

    token = ladtoken ? ladtoken : strtok(0, " \t");
  } /* end while */

  first_unprocessed_arg_ = 1;
}

/* ----------------------------------------------------------------
GetLongOpt::setcell returns
    0	if the nexttoken was not consumed
    1	if the nexttoken was consumed
------------------------------------------------------------------- */

int GetLongOpt::setcell(Cell &c, const char *valtoken, const char *nexttoken,
                        const std::string &name) {
  switch (c.type) {
    case GetLongOpt::NoValue:
      if (*valtoken == '=') {
        throw std::runtime_error("unsolicited value for program option " +
                                 std::string(1, optmarker) + c.option);
      }
      c.value = std::make_unique<std::string>("");
      return 0;
    case GetLongOpt::OptionalValue:
      if (*valtoken == '=') {
        c.value = std::make_unique<std::string>(++valtoken);
      } else if (nexttoken != 0 && nexttoken[0] != optmarker) {
        c.value = std::make_unique<std::string>(nexttoken);
        return 1;
      }
      else
        c.value = std::make_unique<std::string>("");
      return 0;
    case GetLongOpt::MandatoryValue:
      if (*valtoken == '=') {
        c.value = std::make_unique<std::string>(++valtoken);
        return 0;
      } else if (nexttoken != 0 && nexttoken[0] != optmarker) {
        c.value = std::make_unique<std::string>(nexttoken);
        return 1;
      }

      throw std::runtime_error("mandatory value for program option " +
                               std::string(1, optmarker) + c.option + " not specified");
    default:
      break;
  }
  return -1;
}

void GetLongOpt::usage(ostream &outfile) const {
  outfile << "usage: " << pname << " " << ustring << "\n";
  for (const auto &t : table) {
    outfile << "\t" << optmarker << t.option;
    if (t.type == GetLongOpt::MandatoryValue)
      outfile << " <$val>";
    else if (t.type == GetLongOpt::OptionalValue)
      outfile << " [$val]";
    outfile << " (" << t.description << ")\n";
  }
  outfile.flush();
}
