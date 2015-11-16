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
const wchar_t OrbitalIndex::obs_wchar[2] = {L'κ', L'ν'};
const wchar_t OrbitalIndex::dfbs_wchar[2] = {L'Κ', L'Ν'};
const wchar_t OrbitalIndex::abs_wchar[2] = {L'α', L'δ'};

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
        if (first >= occ_wchar[0] && first <= occ_wchar[1]) {
            index_ = Index::occ;
        }
        else if (first >= actocc_wchar[0] && first <= actocc_wchar[1]){
            index_ = Index::actocc;
        }
        else if ( first >= active_wchar[0] && first <= active_wchar[1]){
            index_ = Index::active;
        }
        else if (first >= virt_wchar[0] && first <= virt_wchar[1]) {
            index_ = Index::virt;
        }
        else if (first >= any_wchar[0] && first <= any_wchar[1]) {
            index_ = Index::any;
        }
        else if (first >= obs_wchar[0] && first <= obs_wchar[1]) {
            index_ = Index::obs;
        }
        else if (first >= abs_wchar[0] && first <= abs_wchar[1]){
            index_ = Index::abs;
        }
        else if (first >= dfbs_wchar[0] && first <= dfbs_wchar[1]){
            index_ = Index::dfbs;
        }
        else{
            throw std::runtime_error(error_message);
        }

    }
    // letter with ' or letter with number, number is ignored
    // i1, a2, or A', a', P'
    else if(length == 2){
        if(letter[1]==L'\'') {
            if (first >= othervirt_wchar[0] && first <= othervirt_wchar[1]) {
                index_ = Index::othervirt;
            }
            else if (first >= inactocc_wchar[0] && first <= inactocc_wchar[1]){
                index_ = Index::inactocc;
            }
            else if (first >= allvirt_wchar[0] && first <= allvirt_wchar[1]) {
                index_ = Index::allvirt;
            }
            else if (first >= allany_wchar[0] && first <= allany_wchar[1]) {
                index_ = Index::allany;
            }
            else{
                throw std::runtime_error(error_message);
            }
        }
        else if(isdigit(letter[1])){
            if (first >= occ_wchar[0] && first <= occ_wchar[1]) {
                index_ = Index::occ;
            }
            else if (first >= actocc_wchar[0] && first <= actocc_wchar[1]){
                index_ = Index::actocc;
            }
            else if ( first >= active_wchar[0] && first <= active_wchar[1]){
                index_ = Index::active;
            }
            else if (first >= virt_wchar[0] && first <= virt_wchar[1]) {
                index_ = Index::virt;
            }
            else if (first >= any_wchar[0] && first <= any_wchar[1]) {
                index_ = Index::any;
            }
            else if (first >= obs_wchar[0] && first <= obs_wchar[1]) {
                index_ = Index::obs;
            }
            else if (first >= abs_wchar[0] && first <= abs_wchar[1]){
                index_ = Index::abs;
            }
            else if (first >= dfbs_wchar[0] && first <= dfbs_wchar[1]){
                index_ = Index::dfbs;
            }
            else{
                throw std::runtime_error(error_message);
            }
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
            
            if (first >= othervirt_wchar[0] && first <= othervirt_wchar[1]) {
                index_ = Index::othervirt;
            }
            else if (first >= inactocc_wchar[0] && first <= inactocc_wchar[1]){
                index_ = Index::inactocc;
            }
            else if (first >= allvirt_wchar[0] && first <= allvirt_wchar[1]) {
                index_ = Index::allvirt;
            }
            else if (first >= allany_wchar[0] && first <= allany_wchar[1]) {
                index_ = Index::allany;
            }
            else{
                throw std::runtime_error(error_message);
            }
        }
        //without prime
        else if(isdigit(letter[1])) {
            // check number
            for (auto i = 1ul; i < length; i++) {
                TA_ASSERT(isdigit(letter[i]));
            }

            if (first >= occ_wchar[0] && first <= occ_wchar[1]) {
                index_ = Index::occ;
            }
            else if (first >= actocc_wchar[0] && first <= actocc_wchar[1]) {
                index_ = Index::actocc;
            }
            else if (first >= active_wchar[0] && first <= active_wchar[1]) {
                index_ = Index::active;
            }
            else if (first >= virt_wchar[0] && first <= virt_wchar[1]) {
                index_ = Index::virt;
            }
            else if (first >= any_wchar[0] && first <= any_wchar[1]) {
                index_ = Index::any;
            }
            else if (first >= obs_wchar[0] && first <= obs_wchar[1]) {
                index_ = Index::obs;
            }
            else if (first >= abs_wchar[0] && first <= abs_wchar[1]) {
                index_ = Index::abs;
            }
            else if (first >= dfbs_wchar[0] && first <= dfbs_wchar[1]) {
                index_ = Index::dfbs;
            }
            else {
                throw std::runtime_error(error_message);
            }
        }
        else{
            throw std::runtime_error(error_message);
        }
    }
    std::wstring tmp(letter);
    name_ = tmp;
}

bool OrbitalIndex::operator==(const OrbitalIndex &other) {
    return this->index_ == other.index_;
}

bool OrbitalIndex::operator==(const OrbitalIndex::Index i) {
    return this->index_ == i;
}

bool OrbitalIndex::same(const OrbitalIndex &other) {
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
}
