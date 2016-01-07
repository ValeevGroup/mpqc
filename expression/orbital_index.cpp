//
// Created by Chong Peng on 10/14/15.
//


#include <cwchar>

#include <stdexcept>
#include <TiledArray/error.h>

#include "orbital_index.h"


namespace mpqc{

// set the range of key index to use
const wchar_t OrbitalIndex::inactocc_wchar[2] = {L'm',L'n'};
const wchar_t OrbitalIndex::occ_wchar[2] = {L'm',L'n'};
const wchar_t OrbitalIndex::actocc_wchar[2] = {L'i',L'l'};
const wchar_t OrbitalIndex::virt_wchar[2] = {L'a',L'd'};
const wchar_t OrbitalIndex::active_wchar[2] = {L'x',L'y'};
const wchar_t OrbitalIndex::any_wchar[2] = {L'p',L's'};
const wchar_t OrbitalIndex::othervirt_wchar[2] = {L'a',L'd'};
const wchar_t OrbitalIndex::allvirt_wchar[2] = {L'A',L'D'};
const wchar_t OrbitalIndex::allany_wchar[2] = {L'P',L'S'};
const wchar_t OrbitalIndex::obs_wchar[4] = {L'κ',L'λ',L'μ',L'ν'};
const wchar_t OrbitalIndex::dfbs_wchar[4] = {L'Κ',L'Λ',L'Μ', L'Ν'};
const wchar_t OrbitalIndex::abs_wchar[4] = {L'α', L'β',L'γ',L'δ'};
const wchar_t OrbitalIndex::ribs_wchar[4] = {L'ρ',L'σ',L'τ',L'υ'};

OrbitalIndex::OrbitalIndex(const wchar_t *letter) {
    init(letter);
}

OrbitalIndex::OrbitalIndex(std::wstring letter) {
    init(letter.c_str());
}

//TODO better error handling in throw
void OrbitalIndex::init(const wchar_t* letter) {

    std::string error_message("Wrong Key Index!");

    auto length = wcslen(letter);
        
    // single letter case, for example i, a , p

    const auto first = letter[0];
    if(length == 1) {
        index_ = wchar_to_index(first);
    }
    // letter with ' or letter with number, number is ignored
    // i1, a2, or A', a', P'
    else if(length == 2){
        if(letter[1]==L'\'') {
            index_ = wchar_with_prime_to_index(first);
        }
        else if(isdigit(letter[1])){
            index_ = wchar_to_index(first);
        }
        else{
            throw std::runtime_error(error_message);
        }
    }
    // letter with number and prime
    // A'1, a'1, P'12, i12
    else if(wcslen(letter) >= 3){
        // with prime
        if(letter[1]==L'\'') {
            
            // check number
            for(auto i = 2ul; i < length; i++ ){
                TA_ASSERT(isdigit(letter[i]));
            }
            
            index_ = wchar_with_prime_to_index(first);
        }
        //without prime
        else if(isdigit(letter[1])) {
            // check number
            for (auto i = 1ul; i < length; i++) {
                TA_ASSERT(isdigit(letter[i]));
            }

            index_ = wchar_to_index(first);
        }
        else{
            throw std::runtime_error(error_message);
        }
    }
    std::wstring tmp(letter);
    name_ = tmp;
}

bool OrbitalIndex::operator==(const OrbitalIndex &other) const{
    return this->index_ == other.index_;
}

bool OrbitalIndex::operator==(const OrbitalIndex::Index i) const{
    return this->index_ == i;
}

bool OrbitalIndex::same(const OrbitalIndex &other) const{
    return (index_ == other.index()) && (name_ == other.name());
}

bool OrbitalIndex::is_ao() const {
    int index = static_cast<int> (index_);
    return index < 0;
}

bool OrbitalIndex::is_mo() const {
    int index = static_cast<int> (index_);
    return index > 0;
}

OrbitalIndex::Index OrbitalIndex::wchar_to_index(const wchar_t first) {
    if (first >= occ_wchar[0] && first <= occ_wchar[1]) {
        return Index::occ;
    }
    else if (first >= actocc_wchar[0] && first <= actocc_wchar[1]){
        return Index::actocc;
    }
    else if ( first >= active_wchar[0] && first <= active_wchar[1]){
        return Index::active;
    }
    else if (first >= virt_wchar[0] && first <= virt_wchar[1]) {
        return Index::virt;
    }
    else if (first >= any_wchar[0] && first <= any_wchar[1]) {
        return Index::any;
    }
    else if (first >= obs_wchar[0] && first <= obs_wchar[3]) {
        return Index::obs;
    }
    else if (first >= abs_wchar[0] && first <= abs_wchar[3]){
        return Index::abs;
    }
    else if (first >= dfbs_wchar[0] && first <= dfbs_wchar[3]){
        return Index::dfbs;
    }
    else if (first >= ribs_wchar[0] && first <= ribs_wchar[3]){
        return Index::ribs;
    }
    else{
        throw std::runtime_error("Wrong Key Index!");
        return Index();
    }
}

OrbitalIndex::Index OrbitalIndex::wchar_with_prime_to_index(const wchar_t first) {
    if (first >= othervirt_wchar[0] && first <= othervirt_wchar[1]) {
        return Index::othervirt;
    }
    else if (first >= inactocc_wchar[0] && first <= inactocc_wchar[1]){
        return Index::inactocc;
    }
    else if (first >= allvirt_wchar[0] && first <= allvirt_wchar[1]) {
        return Index::allvirt;
    }
    else if (first >= allany_wchar[0] && first <= allany_wchar[1]) {
        return Index::allany;
    }
    else{
        throw std::runtime_error("Wrong Key Index!");
        return Index();
    }
}

bool OrbitalIndex::is_mo_in_obs() const {
    int index = static_cast<int> (index_);
    return (index > 0 ) && (index <= 9);
}

bool OrbitalIndex::is_mo_in_abs() const {
    return index_ == Index::othervirt;
}

bool OrbitalIndex::is_mo_in_ribs() const {
    int index = static_cast<int> (index_);
    return index >= 15;
}

OrbitalIndex OrbitalIndex::mo_to_ao() {

    std::wstring new_string;

    if(this->is_mo_in_abs()){
        new_string.push_back(abs_wchar[0]);
    }
    else if(this->is_mo_in_obs()){
        new_string.push_back(obs_wchar[0]);
    }
    else if(this->is_mo_in_ribs()){
        new_string.push_back(ribs_wchar[0]);
    }

    if(name_[1] == L'\''){
        new_string.insert(new_string.end(),name_.begin()+2,name_.end());
    }
    else{
        new_string.insert(new_string.end(),name_.begin()+1,name_.end());
    }

    return OrbitalIndex(new_string);
}
}
