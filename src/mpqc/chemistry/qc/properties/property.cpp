

#include "property.h"

namespace mpqc {

std::atomic<typename TimestampFactory::timestamp_type>
    TimestampFactory::current_timestamp (typename TimestampFactory::timestamp_type{0});
}
