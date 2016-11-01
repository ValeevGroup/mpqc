#include "json_handling.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>

#include "mpqc/util/misc/time.h"

namespace mpqc {
namespace json {

using namespace rapidjson;

OutputWriter::OutputWriter() : writer_(buffer_) { writer_.StartObject();
    out_doc_.SetObject();
}

std::ostream &OutputWriter::finalize_and_print(std::ostream &os) {
    writer_.String("Results");
    out_doc_.Accept(writer_);
    writer_.EndObject();
    os << buffer_.GetString() << std::endl;
    return os;
}

void parse_input(int argc, char *argv[], Document &d) {
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

std::unique_ptr<OutputWriter> init_json_writer() {
    std::unique_ptr<OutputWriter> owriter{new OutputWriter};
    auto &writer = owriter->writer();

    auto now = mpqc::system_now();
    auto date = std::chrono::system_clock::to_time_t(now);
    char date_string[512];
    std::strftime(date_string, 512, "%F %T", std::localtime(&date));

    writer.String("Date");
    writer.String(date_string);

    return owriter;
}

std::unique_ptr<OutputWriter> init_json_writer(Document &input) {
    std::unique_ptr<OutputWriter> owriter{new OutputWriter};
    auto &writer = owriter->writer();

    auto now = mpqc::system_now();
    auto date = std::chrono::system_clock::to_time_t(now);
    char date_string[512];
    std::strftime(date_string, 512, "%F %T", std::localtime(&date));

    writer.String("Date");
    writer.String(date_string);
    writer.String("Input");
    input.Accept(writer);

    return owriter;
}

Document get_nested(Document &d, std::string key) {
    const char *key_ctr = key.c_str();

    // if key find
    if(d.HasMember(key_ctr)){

        // if key is object
        if(d[key_ctr].IsObject()){
            rapidjson::StringBuffer buffer;
            rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
            d[key_ctr].Accept(writer);

            rapidjson::Document result;
            rapidjson::StringStream s(buffer.GetString());
            result.ParseStream(s);

            return result;
        }
        else{
            return rapidjson::Document();
        }
    }
    // if key not find, return empty object
    else{
        return rapidjson::Document();
    }
}
} //  namespace json
} //  namespace mpqc
