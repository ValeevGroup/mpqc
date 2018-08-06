#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/core/exception.h"

namespace mpqc {

namespace {
  template <typename T>
  T
  min (const T& t1,
       const T& t2,
       const T& t3) {
    return std::min(std::min(t1,t2),t3);
  }

std::string::size_type
string_distance(const std::string& str1,
                const std::string& str2) {

  typedef std::string::size_type intsize;
  const intsize lenStr1 = str1.length();
  const intsize lenStr2 = str2.length();

  std::vector< std::vector<intsize> > d(lenStr1+1);
  for(intsize i=0; i<=lenStr1; ++i)
    d[i].resize(lenStr2+1);

  for(intsize i=0; i<=lenStr1; ++i)
    d[i][0] = i;
  for(intsize j=1; j<=lenStr2; ++j)
    d[0][j] = j;

  for(intsize i=1; i<=lenStr1; ++i) {
    for(intsize j=1; j<=lenStr2; ++j) {

      intsize cost;
      if (str1[i-1] == str2[j-1])
        cost = 0;
      else
        cost = 1;

      d[i][j] = min(
                   d[i-1][j] + 1,     // deletion
                   d[i][j-1] + 1,     // insertion
                   d[i-1][j-1] + cost   // substitution
                );

      if (i > 1 and j > 1 and str1[i] == str2[j-1] and str1[i-1] == str2[j])
        d[i][j] = std::min(
            d[i][j],
            d[i-2][j-2] + cost   // transposition
        );
    }
  }

  return d[lenStr1][lenStr2];
}
}  // anonymous namespace

DescribedClass::keyval_ctor_wrapper_type DescribedClass::type_to_keyval_ctor(
    const std::string& type_name) {
  auto& registry = keyval_ctor_registry_instance();
  auto iter = registry.find(type_name);
  if (iter == registry.end()) {
    // check if the name was simply misspelled
    std::vector<std::string> candidates;
    // suggest possible misspelling if distance between type_name and any key
    // less than 4 (MPQC3 default)
    for (const auto& e : registry) {
      if (string_distance(type_name, e.first) <= std::size_t(4)) {
        candidates.push_back(e.first);
      }
    }
    std::ostringstream oss;
    oss << "DescribedClass::type_to_keyval_ctor -- type \"" << type_name
        << "\" not registered\n";

    if (!candidates.empty()) {
      for (const auto& candidate : candidates)
        oss << "misspelling of the registered class name \""
            << candidate << "\"?\n";
    }
    throw mpqc::InputError(oss.str().c_str(), __FILE__, __LINE__, "type",
                           type_name.c_str());
  }
  return iter->second.first;
}

///////////////////////////////////////////////

KeyVal::KeyVal()
    : top_tree_(std::make_shared<ptree>()),
      dc_registry_(std::make_shared<dc_registry_type>()), path_(""),
      default_class_key_(std::make_shared<dck_registry_type>())
      {}

KeyVal KeyVal::clone() const {
  return KeyVal(std::make_shared<ptree>(*this->top_tree()),
                std::make_shared<dc_registry_type>(), std::make_shared<dck_registry_type>(), std::string());
}

std::shared_ptr<KeyVal::ptree> KeyVal::tree() const {
  std::shared_ptr<ptree> result(
      top_tree_, &top_tree_->get_child(ptree::path_type{path_, separator}));
  return result;
}

const KeyVal::is_nonnegative_t KeyVal::is_nonnegative{};
const KeyVal::is_positive_t KeyVal::is_positive{};

}  // namespace mpqc
