
#include <mpqc/util/keyval/keyval.hpp>

namespace mpqc {

  DescribedClass::keyval_ctor_wrapper_type
  DescribedClass::type_to_keyval_ctor(const std::string& type_name) {
    auto& registry = keyval_ctor_registry();
    if (registry.find(type_name) == registry.end())
      throw std::runtime_error("DescribedClass::type_to_keyval_ctor -- type not registered");
    return registry[type_name];
  }

  ///////////////////////////////////////////////

  KeyVal::KeyVal() : top_tree_(std::make_shared<ptree>()),
      class_registry_(std::make_shared<class_registry_type>()),
      path_("") {}

  std::shared_ptr<KeyVal::ptree>
  KeyVal::tree() const {
    std::shared_ptr<ptree> result(top_tree_,
                                  &top_tree_->get_child(ptree::path_type{path_, separator}));
    return result;
  }

} // namespace mpqc
