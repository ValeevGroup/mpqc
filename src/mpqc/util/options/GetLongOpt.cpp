#include "mpqc/util/options/GetLongOpt.h"

using namespace std;
using namespace mpqc;

GetLongOpt::GetLongOpt(const char optmark)
    : ustring("[valid options and arguments]"),
      optmarker(optmark),
      finalized(false) {}

std::string GetLongOpt::basename(const std::string &pname) const {
  auto pos = pname.find_last_of('/');
  if (pos == std::string::npos)
    return pname;
  else
    return pname.substr(pos + 1);
}

int GetLongOpt::enroll(const std::string &opt, const OptType t,
                       const std::string &desc, const std::string &val) {
  if (finalized) return 0;
  table.emplace_back(Cell{
      opt, t, (!desc.empty() ? desc : std::string("no description available")),
      val});
  return 1;
}

std::string GetLongOpt::retrieve(const std::string &opt) const {
  for (const auto &t : table) {
    if (opt == t.option) return t.value;
  }
  std::cerr << "GetLongOpt::retrieve - unenrolled option ";
  std::cerr << optmarker << opt << "\n";
  return 0;
}

int GetLongOpt::parse(int argc, char *const *argv) {
  int optind = 1;

  pname = basename(*argv);
  finalized = true;
  if (argc-- <= 1) return optind;

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
      throw std::runtime_error("unrecognized option " +
                               std::to_string(optmarker) +
                               std::string(strtok(token, "= ")));
    }

  } /* end while */

  return optind;
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
                                 std::to_string(optmarker) + c.option);
      }
      return 0;
    case GetLongOpt::OptionalValue:
      if (*valtoken == '=') {
        c.value = ++valtoken;
      } else if (nexttoken != 0 && nexttoken[0] != optmarker) {
        c.value = nexttoken;
        return 1;
      }
      return 0;
    case GetLongOpt::MandatoryValue:
      if (*valtoken == '=') {
        c.value = ++valtoken;
        return 0;
      } else if (nexttoken != 0 && nexttoken[0] != optmarker) {
        c.value = nexttoken;
        return 1;
      }

      throw std::runtime_error("mandatory value for program option " +
                               std::to_string(optmarker) + c.option +
                               " not specified");
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
