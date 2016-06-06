//
// Created by Chong Peng on 10/15/15.
//

#include <algorithm>
#include <iostream>


#include <mpqc/util/expression/formula.h>

#include <boost/algorithm/string.hpp>
#include <TiledArray/error.h>

namespace mpqc{


    Formula::Formula(std::wstring formula_string) {
        if (formula_string.empty()) {
            throw std::runtime_error("Empty Formula!");
        }

        // detect the brackets <> or ()
        auto left_mark = [](const wchar_t letter) {
            return letter == L'<' || letter == L'(';
        };

        // find formula between < > or ( )
        auto bra_symbol = std::find_if(formula_string.cbegin(), formula_string.cend(), left_mark);

        auto ket_symbol = formula_string.cbegin();

        if(*bra_symbol==L'<'){
            auto pos = formula_string.find_last_of(L'>');
            ket_symbol += pos;
        }
        else{
            auto pos = formula_string.find_last_of(L')');
            ket_symbol += pos;
        }

        TA_ASSERT(bra_symbol < ket_symbol);

        if (bra_symbol == formula_string.cend() || ket_symbol == formula_string.cend()) {
            throw std::runtime_error("Formula should start with < or ( and end with > or ). \n");
        }

        TA_ASSERT( (*bra_symbol==L'<'&&*ket_symbol==L'>') || (*bra_symbol==L'(')&&*ket_symbol==L')');

        if(*bra_symbol==L'<'){
            notation_ = Notation::Physical;
        }else{
            notation_ = Notation::Chemical;
        }

        std::wstring formula(bra_symbol + 1, ket_symbol);

//        std::wcout << formula_string;


        // find option between [ and ]
        auto option_left = std::find(formula_string.cbegin(), formula_string.cend(), L'[');
        auto option_right = std::find(formula_string.cbegin(), formula_string.cend(), L']');

        std::wstring option;

        if( (option_left!=formula_string.cend()) && (option_right!=formula_string.cend())){
            TA_ASSERT(option_left < option_right);
            option = std::wstring(option_left+1,option_right);
        }


        // split the formula by |
        std::vector<std::wstring> split_formula;
        boost::split(split_formula, formula, boost::is_any_of(L"|"), boost::token_compress_on);

        if (split_formula.size() == 2) {
            std::wstring left_formula = std::move(split_formula[0]);
            std::wstring right_formula = std::move(split_formula[1]);

            operation_ = Operation(L" ",option);
            left_index_ = check_orbital_index(left_formula);
            right_index_ = check_orbital_index(right_formula);

        }
        else if (split_formula.size() == 3) {
            std::wstring left_formula = std::move(split_formula[0]);
            std::wstring operation = std::move(split_formula[1]);
            std::wstring right_formula = std::move(split_formula[2]);

            operation_ = Operation(operation, option);
            left_index_ = check_orbital_index(left_formula);
            right_index_ = check_orbital_index(right_formula);

        }
        else {
            throw std::runtime_error("Formula in wrong length. \n");
        }

        // error detecting


        // one body operation
        if(operation_.is_onebody() && (left_index_.size()!=1) && (right_index_.size()!=1)){
            throw std::runtime_error("One body Operator with Wrong Index Size!");
        }

        // more than three index
        if( (left_index_.size()>=3) || (right_index_.size()>=3)){
            throw std::runtime_error("Wrong Number of Index!");
        }

    }

    std::size_t Formula::rank() const {
        return (left_index_.size() + right_index_.size());
    }

    std::vector<OrbitalIndex> Formula::check_orbital_index(std::wstring index_array) {

        std::vector<std::wstring> split_index;
        boost::split(split_index,index_array,boost::is_any_of(L" "),boost::token_compress_on);

        if(split_index.empty()){
            return std::vector<OrbitalIndex>();
        }else{
            std::vector<OrbitalIndex> result;

            for(auto index : split_index){
                if(!index.empty()){
                    result.emplace_back(index);
                }
            }
            return result;
        }
    }

    bool Formula::operator<(const Formula &other) const {

        if(operation()!=other.operation()){
            return operation() < other.operation();
        }
        else if(notation_ != other.notation()){
            return notation_ < other.notation();
        }
        else if(left_index() != other.left_index()){
            return left_index() < other.left_index();
        }
        else{
            return right_index() < other.right_index();
        }
    }

    bool Formula::operator==(const Formula &other) const{
        return (operation_== other.operation_) && (left_index_ == other.left_index_) && (right_index_ == other.right_index_) && (notation_ == other.notation_);
    }

    std::wstring Formula::formula_string() const {

        std::wstring left, right, result;
        for (const auto& index : left_index_){
            left += index.name() + L" ";
        }

        for (const auto& index : right_index_){
            right += L" " + index.name();
        }

        // add operation
        auto oper_str = operation_.oper_string();
        if(!oper_str.empty()){
            result = left + L"|" + operation_.oper_string() + L"|" + right;
        }
        else{
            result = left + L"|" + right;
        }

        if(notation_==Notation::Chemical){
            result = L"( " + result + L" )";
        }
        else{
            result = L"< " + result + L" >";
        }
        // add option
        result = result + operation_.option_string();

        return result;

    }

    std::string Formula::to_ta_expression() const {

        std::string ta_expression;
        std::size_t rank = this->rank();
        std::size_t count = 0;

        // add left index
        for (const auto & index : left_index_){
            std::string index_expression = index.to_ta_expression();
            ta_expression.append(index_expression.begin(),index_expression.end());
            ++count;
            ta_expression.append(", ");
        }

        // add right index
        for (const auto & index : right_index_){
            std::string index_expression = index.to_ta_expression();
            ta_expression.append(index_expression.begin(),index_expression.end());
            ++count;
            if(count < rank){
                ta_expression.append(", ");
            }
        }

        return ta_expression;
    }

bool Formula::has_index(const OrbitalIndex &ind) const {

    for(auto& index : left_index_){
        if(ind == index){
            return true;
        }
    }

    for(auto& index : right_index_){
        if(ind == index){
            return true;
        }
    }

    return false;
}


}
