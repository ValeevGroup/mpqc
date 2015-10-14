//
// Created by Chong Peng on 10/14/15.
//

#include "orbital_index.h"

namespace mpqc{


OrbitalIndex::OrbitalIndex(const char* letter) {

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
            throw std::runtime_error("Wrong Key Index");
        }

    }
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
                throw std::runtime_error("Wrong Key Index");
            }
        }else{
            throw std::runtime_error("Wrong Key Index");
        }
    }
    else{
        throw std::runtime_error("Wrong Key Length");
    }
}

bool OrbitalIndex::operator==(const OrbitalIndex &other) {
    return this->index_ == other.index_;
}


}
