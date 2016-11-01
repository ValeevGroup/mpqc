#pragma once

#ifndef MPQC_UTILITY_JSONHANDLER_H
#define MPQC_UTILITY_JSONHANDLER_H

#include <rapidjson/rapidjson.h>
#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>

#include <ostream>
#include <memory>

namespace mpqc {
namespace json {

class OutputWriter {
  private:
    using buffer_type = rapidjson::StringBuffer;
    buffer_type buffer_;

    using writer_type = rapidjson::PrettyWriter<buffer_type>;
    writer_type writer_;

    rapidjson::Document out_doc_;

  public:
    OutputWriter();
    OutputWriter(OutputWriter const &) = delete;
    OutputWriter(OutputWriter &&) = default;
    OutputWriter &operator=(OutputWriter const &) = delete;
    OutputWriter &operator=(OutputWriter &&) = default;
    ~OutputWriter() = default;

    inline writer_type &writer() { return writer_; }
    inline rapidjson::Document &doc() { return out_doc_; }

    std::ostream &finalize_and_print(std::ostream &os);
};

/*! \brief looks for input json file and parses it to a document
 *
 * Parse an input json file into a document.  You must create a document
 * and pass it in yourself because rapidjson::Document does not have
 * a copy constructor.
 */
void parse_input(int argc, char *argv[], rapidjson::Document &d);

std::unique_ptr<OutputWriter> init_json_writer();

std::unique_ptr<OutputWriter> init_json_writer(rapidjson::Document &);

rapidjson::Document get_nested(rapidjson::Document &d, std::string key);

} //  namespace json
} //  namespace mpqc

#endif // MPQC_UTILITY_JSONHANDLER_H
