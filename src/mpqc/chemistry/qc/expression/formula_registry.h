//
// Created by pchong on 11/10/15.
//

#ifndef MPQC_FORMULA_MAP_H
#define MPQC_FORMULA_MAP_H

#include "mpqc/chemistry/qc/expression/formula.h"
#include "mpqc/math/external/tiledarray/array_info.h"
#include "mpqc/util/external/madworld/parallel_print.h"
#include <map>

namespace mpqc {

/**
 *
 * \brief a wrapper of std::map, not thread_safe for writing
 *
 *
 *
 */
template <typename Key, typename Value>
class Registry {
 public:
  using value_type = std::pair<Key, std::shared_ptr<Value>>;
  using element_type = std::map<Key, std::shared_ptr<Value>>;
  using iterator = typename element_type::iterator;
  using const_iterator = typename element_type::const_iterator;

  Registry() = default;

  Registry(const element_type& map) : registry_(map) {}

  Registry(Registry const&) = default;
  Registry& operator=(Registry const&) = default;

  Registry(Registry&&) = default;
  Registry& operator=(Registry&&) = default;

  virtual ~Registry() = default;

  /// return the registry
  const element_type& registry() const { return registry_; }

  /// insert to registry by std::pair<Key, Value>
  void insert(const value_type& val) {
    auto insert_result = registry_.insert(val);
    if (insert_result.second == false) {
      throw std::runtime_error("Key Already Exist!!!");
    }
  }

  /// insert to registry by Key and Value
  void insert(const Key& key, const Value& val) {
    insert(std::make_pair(key, std::make_shared<Value>(val)));
  }

  /// remove Value by Key
  void remove(const Key& key) { registry_.erase(key); }

  /// clear the registry
  void clear() { registry_.clear(); }

  /// find item, return iterator
  iterator find(const Key& key) { return registry_.find(key); }

  /// find item, return const iterator
  const_iterator find(const Key& key) const { return registry_.find(key); }

  /// check if have key in registry
  bool have(const Key& key) const {
    auto iter = registry_.find(key);
    return iter != registry_.end();
  }

  /// find item, return iterator, if not found, throw
  Value& retrieve(const Key& key) {
    auto iter = registry_.find(key);
    if (iter == registry_.end()) {
      throw std::runtime_error("Key not found!");
    }
    return *(iter->second);
  }

  /// find item, return const iterator, if not found, throw
  const Value& retrieve(const Key& key) const {
    auto iter = registry_.find(key);
    if (iter == registry_.cend()) {
      throw std::runtime_error("Key not found!");
    }
    return *(iter->second);
  }

  /// return begin of iterator
  iterator begin() { return registry_.begin(); }

  /// return end of iterator
  iterator end() { return registry_.end(); }

  /// return begin of const_iterator
  const_iterator cbegin() const { return registry_.cbegin(); }

  /// return end of const_iterator
  const_iterator cend() const { return registry_.cend(); }

  /// purges all objects if p(key) == true
  template <typename Pred>
  void purge_if(const Pred& p) {
    auto i = registry_.begin();
    for (; i != registry_.end();) {
      if (p(*i)) {
        registry_.erase(i++);
      } else {
        ++i;
      }
    }
  }

 protected:
  element_type registry_;
};

/**
 *
 *  \brief map Formula to Value object
 *
 */
template <typename Value>
class FormulaRegistry : public Registry<Formula, Value> {
 public:
  using Key = Formula;
  using value_type = typename Registry<Key, Value>::value_type;
  using element_type = typename Registry<Key, Value>::element_type;
  using iterator = typename Registry<Key, Value>::iterator;
  using const_iterator = typename Registry<Key, Value>::const_iterator;

  FormulaRegistry() = default;
  FormulaRegistry(const element_type& map) : Registry<Key, Value>(map) {}

  /// prevent from copy and assign of FormulaRegistry
  FormulaRegistry(FormulaRegistry const&) = delete;
  FormulaRegistry& operator=(FormulaRegistry const&) = delete;

  FormulaRegistry(FormulaRegistry&&) = default;
  FormulaRegistry& operator=(FormulaRegistry&&) = default;

  /// print out formula that stored in registry
  void print_formula(madness::World& world) const {
    if (world.rank() == 0) {
      for (const auto& item : this->registry_) {
        mpqc::detail::print_size_info(item.second, item.first.string());
      }
      std::cout << std::endl;
    }
  }

  /// purges all objects if p(key) == true
  template <typename Pred>
  void purge_if(madness::World& world, const Pred& p) {
    //    world.gop.fence();
    auto i = this->registry_.begin();
    for (; i != this->registry_.end();) {
      if (p(*i)) {
        utility::print_par(world, "Removed from Registry: ");
        utility::print_par(world, utility::to_string(i->first.string()));
        utility::print_par(world, "\n");
        this->registry_.erase(i++);
      } else {
        ++i;
      }
    }
  }

  /// purges formulae that contain Operator whose type matches \c optype
  void purge_operator(madness::World& world, const Operator::Type& optype) {
    auto pred = [optype](const value_type& item) {
      return item.first.oper().type() == optype;
    };

    this->purge_if(world, pred);
  }

  /// purges formulae that contain Operator described by string \c opstr
  void purge_operator(madness::World& world, const std::wstring& opstr) {
    Operator oper(opstr);
    Operator::Type oper_type = oper.type();
    purge_operator(world, oper_type);
  }

  /// purges the Formula object that equals \c formula from the registry
  void purge_formula(madness::World& world, const Formula& formula) {
    auto pred = [&formula](const value_type& item) {
      return item.first == formula;
    };

    this->purge_if(world, pred);
  }

  /// purges the formula that that corresponds to string \c str
  void purge_formula(madness::World& world, const std::wstring& str) {
    Formula formula(str);
    purge_formula(world, formula);
  }

  /// purges formulae that contain index \c idx
  void purge_index(madness::World& world, const OrbitalIndex& idx) {
    auto pred = [&idx](const value_type& item) {
      return item.first.has_index(idx);
    };

    this->purge_if(world, pred);
  }

  /// purges all formula in registry
  void purge(madness::World& world) {
    auto pred = [](const value_type& item) { return true; };

    this->purge_if(world, pred);
  }
};
}  // end of namespace mpqc

#endif  // MPQC_FORMULA_MAP_H
