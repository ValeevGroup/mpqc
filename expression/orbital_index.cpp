//
// Created by Chong Peng on 10/14/15.
//


#include <cwchar>

#include <stdexcept>

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

    // single letter case, for example i, a , p
    if(wcslen(letter) == 1) {
        if (letter[0] >= occ_wchar[0] && letter[0] <= occ_wchar[1]) {
            index_ = Index::occ;
        }
        else if (letter[0] >= actocc_wchar[0] && letter[0] <= actocc_wchar[1]){
            index_ = Index::actocc;
        }
        else if ( letter[0] >= active_wchar[0] && letter[0] <= active_wchar[1]){
            index_ = Index::active;
        }
        else if (letter[0] >= virt_wchar[0] && letter[0] <= virt_wchar[1]) {
            index_ = Index::virt;
        }
        else if (letter[0] >= any_wchar[0] && letter[0] <= any_wchar[1]) {
            index_ = Index::any;
        }
        else if (letter[0] >= obs_wchar[0] && letter[0] <= obs_wchar[1]) {
            index_ = Index::obs;
        }
        else if (letter[0] >= abs_wchar[0] && letter[0] <= abs_wchar[1]){
            index_ = Index::abs;
        }
        else if (letter[0] >= dfbs_wchar[0] && letter[0] <= dfbs_wchar[1]){
            index_ = Index::dfbs;
        }
        else{
            throw std::runtime_error(error_message);
        }

    }
    // letter with ' or letter with number, number is ignored
    // i1, a2, or A', a', P'
    else if(wcslen(letter) == 2){
        if(letter[1]=='\'') {
            if (letter[0] >= othervirt_wchar[0] && letter[0] <= othervirt_wchar[1]) {
                index_ = Index::othervirt;
            }
            else if (letter[0] >= inactocc_wchar[0] && letter[0] <= inactocc_wchar[1]){
                index_ = Index::inactocc;
            }
            else if (letter[0] >= allvirt_wchar[0] && letter[0] <= allvirt_wchar[1]) {
                index_ = Index::allvirt;
            }
            else if (letter[0] >= allany_wchar[0] && letter[0] <= allany_wchar[1]) {
                index_ = Index::allany;
            }
            else{
                throw std::runtime_error(error_message);
            }
        }
        else if(isdigit(letter[1])){
            if (letter[0] >= occ_wchar[0] && letter[0] <= occ_wchar[1]) {
                index_ = Index::occ;
            }
            else if (letter[0] >= actocc_wchar[0] && letter[0] <= actocc_wchar[1]){
                index_ = Index::actocc;
            }
            else if ( letter[0] >= active_wchar[0] && letter[0] <= active_wchar[1]){
                index_ = Index::active;
            }
            else if (letter[0] >= virt_wchar[0] && letter[0] <= virt_wchar[1]) {
                index_ = Index::virt;
            }
            else if (letter[0] >= any_wchar[0] && letter[0] <= any_wchar[1]) {
                index_ = Index::any;
            }
            else if (letter[0] >= obs_wchar[0] && letter[0] <= obs_wchar[1]) {
                index_ = Index::obs;
            }
            else if (letter[0] >= abs_wchar[0] && letter[0] <= abs_wchar[1]){
                index_ = Index::abs;
            }
            else if (letter[0] >= dfbs_wchar[0] && letter[0] <= dfbs_wchar[1]){
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
    // A1', a2', P3'
    else if(wcslen(letter) == 3){
        if(letter[2]=='\'' && isdigit(letter[1])) {
            if (letter[0] >= othervirt_wchar[0] && letter[0] <= othervirt_wchar[1]) {
                index_ = Index::othervirt;
            }
            else if (letter[0] >= inactocc_wchar[0] && letter[0] <= inactocc_wchar[1]){
                index_ = Index::inactocc;
            }
            else if (letter[0] >= allvirt_wchar[0] && letter[0] <= allvirt_wchar[1]) {
                index_ = Index::allvirt;
            }
            else if (letter[0] >= allany_wchar[0] && letter[0] <= allany_wchar[1]) {
                index_ = Index::allany;
            }
            else{
                throw std::runtime_error(error_message);
            }
        }
        else{
            throw std::runtime_error(error_message);
        }
    }
    else{
        throw std::runtime_error(error_message);
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
