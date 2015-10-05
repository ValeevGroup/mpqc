#pragma once

#ifndef TCC_UTILITY_JSONINPUT_H
#define TCC_UTILITY_JSONINPUT_H

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <iostream>
#include <fstream>
#include <stdexcept>


using namespace rapidjson;

/*! \brief looks for input json file and parses it to a document
 *
 * Parse an input json file into a document.  You must create a document 
 * and pass it in yourself because rapidjson::Document does not have 
 * a copy constructor.
 */
void parse_input(int argc, char *argv[], rapidjson::Document &d) {
    if (argc != 2) {
      std::cout << "usage: " << argv[0] << " <input_file.json>" << std::endl;
      throw std::invalid_argument("no input file given");
    }

    std::string input_file_name = argv[1];
    std::ifstream input_file(input_file_name, std::ifstream::in);

    std::string contents((std::istreambuf_iterator<char>(input_file)),
                         std::istreambuf_iterator<char>());

    const char *json = contents.c_str();
    d.Parse(json);
}

Document get_nested(Document &d, std::string key){
    rapidjson::StringBuffer buffer;
    const char *key_ctr = key.c_str();

    assert(d[key_ctr].IsObject());

    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    d[key_ctr].Accept(writer);

    rapidjson::Document result;
    rapidjson::StringStream s(buffer.GetString());
    result.ParseStream(s);

    return result;

}

#endif // TCC_UTILITY_JSONINPUT_H
