//
// Created by pchong on 11/10/15.
//

#ifndef TILECLUSTERCHEM_FORMULA_MAP_H
#define TILECLUSTERCHEM_FORMULA_MAP_H

#include <map>
#include "formula.h"
#include "../utility/print_size_info.h"

namespace mpqc{

    // a wrapper of std::map, not thread_safe for writing
    template<typename Key, typename Value>
    class Registry {
    public:

        using value_type = std::pair<Key,Value>;
        using element_type = std::map<Key,Value>;
        using iterator = typename element_type::iterator;
        using const_iterator = typename element_type::const_iterator;

        Registry() = default;

        Registry(const element_type &map) : registry_(map) { }

        Registry(Registry const &) = default;
        Registry& operator=(Registry const &) = default;

        Registry(Registry&&) = default;
        Registry& operator=(Registry &&) = default;

        virtual ~Registry()= default;

        const element_type &registry() const {
            return registry_;
        }

        void insert(const value_type& val){
            registry_.insert(val);
        }

        void insert(const Key& key, const Value& val){
            registry_[key] = val;
        }

        void remove(const Key& key){
            registry_.erase(key);
        }

        void clear(){
            registry_.clear();
        }

        // find item, return iterator
        iterator find(const Key& key) {
            return registry_.find(key);
        }

        // find item, return const iterator
        const_iterator find(const Key& key) const{
            return registry_.find(key);
        }

        // find item, return iterator, if not found, throw
        Value retrieve(const Key& key){
            auto iter = registry_.find(key);
            if(iter==registry_.end()){
                throw std::runtime_error("Key not found!");
            }
            return iter->second;
        }

        // return item, return const iterator, if not found, throw
        const Value retrieve(const Key& key) const{
            auto iter = registry_.find(key);
            if(iter==registry_.cend()){
                throw std::runtime_error("Key not found!");
            }
            return iter->second;
        }

        iterator begin() {
            return registry_.begin();
        }

        iterator end()  {
            return registry_.end();
        }

        const_iterator cbegin() const {
            return registry_.cbegin();
        }

        const_iterator cend() const {
            return registry_.cend();
        }

        /// removes all objects if p(key) == true
        template<typename Pred>
        void remove_if(const Pred& p){
            auto i = registry_.begin();
            for(; i != registry_.end(); ){
                if (p(*i)){
                    registry_.erase(i++);
                }else{
                    ++i;
                }
            }
        }

    protected:
        element_type registry_;
    };

    // map formula to template object
    template<typename Value>
    class FormulaRegistry : public Registry<Formula,Value>{

    public:

        using Key = Formula;
        using value_type = typename Registry<Key,Value>::value_type;
        using element_type = typename Registry<Key,Value>::element_type;
        using iterator = typename Registry<Key,Value>::iterator;
        using const_iterator = typename Registry<Key,Value>::const_iterator;

        FormulaRegistry()= default;
        FormulaRegistry(const element_type& map) : Registry<Key,Value>(map){}

        FormulaRegistry(FormulaRegistry const &) = delete;
        FormulaRegistry& operator=(FormulaRegistry const &) = delete;

        FormulaRegistry(FormulaRegistry&&) = default;
        FormulaRegistry& operator=(FormulaRegistry &&) = default;

        void print_formula(madness::World& world) const {

            if(world.rank()==0){
                for(const auto& item : this->registry_){
                    utility::wprint_size_info(item.second, item.first.formula_string());
                }
                std::cout << std::endl;
            }
        }
    };
} // end of namespace mpqc



#endif //TILECLUSTERCHEM_FORMULA_MAP_H
