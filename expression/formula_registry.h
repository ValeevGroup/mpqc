//
// Created by pchong on 11/10/15.
//

#ifndef TILECLUSTERCHEM_FORMULA_MAP_H
#define TILECLUSTERCHEM_FORMULA_MAP_H

#include <map>
#include "formula.h"

namespace mpqc{

    // map formula to template object
    // a wrapper of std::map, not thread_safe for writing
    template<typename Value>
    class FormulaRegistry {
    public:

        using value_type = std::pair<Formula,Value>;
        using element_type = std::map<Formula,Value>;
        using iterator = typename element_type::iterator;
        using const_iterator = typename element_type::const_iterator;

        FormulaRegistry() = default;

        FormulaRegistry(const std::map<Formula, Value> &formula_map) : formula_registry_(formula_map) { }

        const std::map<Formula, Value>& map() const {
            return formula_registry_;
        }

        void insert(const value_type& val){
            formula_registry_.insert(val);
        }

        void insert(const Formula& formula, const Value& val){
            formula_registry_[formula] = val;
        }

        void remove(const Formula& formula){
            formula_registry_.erase(formula);
        }

        void clear(){
            formula_registry_.clear();
        }

        iterator find(const Formula& formula){
            return formula_registry_.find(formula);
        }

        const_iterator find(const Formula& formula) const{
            return formula_registry_.find(formula);
        }

        iterator begin() {
            return formula_registry_.begin();
        }

        iterator end()  {
            return formula_registry_.end();
        }

        const_iterator cbegin() const {
            return formula_registry_.cbegin();
        }

        const_iterator cend() const {
            return formula_registry_.cend();
        }

    private:

        /// removes all objects if p(key) == true
        template<typename Pred>
        void remove_if(const Pred& p){
            auto i = formula_registry_.begin();
            for(; i != formula_registry_.end(); ){
                if (p(*i)){
                    formula_registry_.erase(i++);
                }else{
                    ++i;
                }
            }
        }
        std::map<Formula, Value> formula_registry_;
    };
}



#endif //TILECLUSTERCHEM_FORMULA_MAP_H
