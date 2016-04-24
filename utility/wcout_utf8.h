//
// Created by Chong Peng on 11/3/15.
//

#ifndef MPQC_WCOUT_UTF8_H
#define MPQC_WCOUT_UTF8_H

#include <iostream>
#include <string>

//TODO wcout on linux system
inline void wcout_utf8(const std::wstring& s){

    std::wcout.imbue(std::locale("en_US.UTF-8"));
    std::wcout << s;
}

#endif //MPQC_WCOUT_UTF8_H
