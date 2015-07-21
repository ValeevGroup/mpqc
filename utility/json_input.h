#pragma once

#ifndef TCC_UTILITY_JSONINPUT_H
#define TCC_UTILITY_JSONINPUT_H

#include <rapidjson/document.h>
#include <iostream>
#include <fstream>

/*! \brief looks for input json file and parses it to a document
 *
 * Parse an input json file into a document.  You must create a document 
 * and pass it in yourself because rapidjson::Document does not have 
 * a copy constructor.
 */
void parse_input(int argc, char *argv[], rapidjson::Document &d) {
    if (argc == 1) {
        std::cout << "No input file" << std::endl;
        std::abort();
    }

    std::string input_file_name = argv[1];
    std::ifstream input_file(input_file_name, std::ifstream::in);

    std::string contents((std::istreambuf_iterator<char>(input_file)),
                         std::istreambuf_iterator<char>());

    const char *json = contents.c_str();
    d.Parse(json);
}

#endif // TCC_UTILITY_JSONINPUT_H
