//
// Created by pchong on 11/10/15.
//

#ifndef TILECLUSTERCHEM_FORMULA_MAP_H
#define TILECLUSTERCHEM_FORMULA_MAP_H

#include <bits/unordered_map.h>
#include "formula.h"

namespace mpqc{

    // map formula to template object
    // a wrapper of std::unordered_map, not thread_safe for writing
    template<typename Value>
    class FormulaMap {
    public:

        using value_type = std::pair<Formula,Value>;
        using element_type = std::unordered_map<Formula,Value>;
        using iterator = typename element_type::iterator;
        using const_iterator = typename element_type::const_iterator;

        FormulaMap() = default;

        FormulaMap(const std::unordered_map<Formula, Value> &formula_map) : formula_map_(formula_map) { }

        const std::unordered_map<Formula, Value>& map() const {
            return formula_map_;
        }

        void insert(const value_type& val){
            formula_map_.insert(val);
        }

        void insert(const Formula& formula, const Value& val){
            formula_map_[formula] = val;
        }

        void remove(const Formula& formula){
            formula_map_.erase(formula);
        }

        /// removes all objects if p(key) == true
        template<typename Pred>
        void remove_if(const Pred& p){
            auto i = formula_map_.begin();
            for(; i != formula_map_.end(); ){
                if (p(*i)){
                    formula_map_.erase(i++);
                }else{
                    ++i;
                }
            }
        }

        void clear(){
            formula_map_.clear();
        }

        iterator find(const Formula& formula){
            return formula_map_.find(formula);
        }

        const_iterator find(const Formula& formula) const{
            return formula_map_.find(formula);
        }

    private:
        std::unordered_map<Formula, Value> formula_map_;
    };
}



#endif //TILECLUSTERCHEM_FORMULA_MAP_H
