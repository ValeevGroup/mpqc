//
// Created by Chong Peng on 10/14/15.
//

#include "orbital_index.h"

namespace mpqc{


// set the range of key index to use
const char OrbitalIndex::occ_char[2] = {'i','l'};
const char OrbitalIndex::virt_char[2] = {'a','d'};
const char OrbitalIndex::any_char[2] = {'p','s'};
const char OrbitalIndex::othervirt_char[2] = {'a','d'};
const char OrbitalIndex::allvirt_char[2] = {'A','D'};
const char OrbitalIndex::allany_char[2] = {'P','S'};

OrbitalIndex::OrbitalIndex(const char *letter) {
    init(letter);
}

OrbitalIndex::OrbitalIndex(std::string letter) {
    init(letter.c_str());
}

void OrbitalIndex::init(const char* letter) {

    std::string error_message(letter);
    error_message = "Wrong Key Index: " + error_message;

    // single letter case, for example i, a , p
    if(strlen(letter) == 1) {
        if (letter[0] >= occ_char[0] && letter[0] <= occ_char[1]) {
            index_ = Index::occ;
        }
        else if (letter[0] >= virt_char[0] && letter[0] <= virt_char[1]) {
            index_ = Index::virt;
        }
        else if (letter[0] >= any_char[0] && letter[0] <= any_char[1]) {
            index_ = Index::any;
        }
        else{
            throw std::runtime_error(error_message);
        }

    }
    // letter with ' or letter with number, number is ignored
    // i1, a2, or A', a', P'
    else if(strlen(letter) == 2){
        if(letter[1]=='\'') {
            if (letter[0] >= othervirt_char[0] && letter[0] <= othervirt_char[1]) {
                index_ = Index::othervirt;
            }
            else if (letter[0] >= allvirt_char[0] && letter[0] <= allvirt_char[1]) {
                index_ = Index::allvirt;
            }
            else if (letter[0] >= allany_char[0] && letter[0] <= allany_char[1]) {
                index_ = Index::allany;
            }
            else{
                throw std::runtime_error(error_message);
            }
        }
        else if(isdigit(letter[1])){
            if (letter[0] >= occ_char[0] && letter[0] <= occ_char[1]) {
                index_ = Index::occ;
            }
            else if (letter[0] >= virt_char[0] && letter[0] <= virt_char[1]) {
                index_ = Index::virt;
            }
            else if (letter[0] >= any_char[0] && letter[0] <= any_char[1]) {
                index_ = Index::any;
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
    else if(strlen(letter) == 3){
        if(letter[2]=='\'' && isdigit(letter[1])) {
            if (letter[0] >= othervirt_char[0] && letter[0] <= othervirt_char[1]) {
                index_ = Index::othervirt;
            }
            else if (letter[0] >= allvirt_char[0] && letter[0] <= allvirt_char[1]) {
                index_ = Index::allvirt;
            }
            else if (letter[0] >= allany_char[0] && letter[0] <= allany_char[1]) {
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
    std::string tmp(letter);
    name_ = tmp;
}

bool OrbitalIndex::operator==(const OrbitalIndex &other) {
    return this->index_ == other.index_;
}


bool OrbitalIndex::operator==(const OrbitalIndex::Index i) {
    return this->index_ == i;
}

}
