#pragma once
#ifndef TCC_UTILITY_MEMORYUSAGE_H
#define TCC_UTILITY_MEMORYUSAGE_H

namespace mpqc {
namespace utility {

void process_linux_mem_usage(double &vm_usage, double &resident_set);

} // namespace utility
} // namespace mpqc 

#endif // TCC_UTILITY_MEMORYUSAGE_H
