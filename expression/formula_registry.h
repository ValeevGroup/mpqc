//
// Created by pchong on 11/10/15.
//

#ifndef TILECLUSTERCHEM_FORMULA_MAP_H
#define TILECLUSTERCHEM_FORMULA_MAP_H

#include <map>
#include "formula.h"

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


        virtual ~Registry()= default;

        const element_type & map() const {
            return registry_;
        }

        void insert(const value_type& val){
            registry_.insert(val);
        }

        void insert(const Formula& formula, const Value& val){
            registry_[formula] = val;
        }

        void remove(const Formula& formula){
            registry_.erase(formula);
        }

        void clear(){
            registry_.clear();
        }

        iterator find(const Formula& formula){
            return registry_.find(formula);
        }

        const_iterator find(const Formula& formula) const{
            return registry_.find(formula);
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

    };
} // end of namespace mpqc



#endif //TILECLUSTERCHEM_FORMULA_MAP_H
