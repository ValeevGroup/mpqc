//
// Created by Chong Peng on 11/3/15.
//

#ifndef TILECLUSTERCHEM_WCOUT_UTF8_H
#define TILECLUSTERCHEM_WCOUT_UTF8_H

#include <iostream>
#include <string>

void wcout_utf8(const std::wstring& s){

    std::wcout.imbue(std::locale("en_US.UTF-8"));
    std::wcout << s;
}

#endif //TILECLUSTERCHEM_WCOUT_UTF8_H
