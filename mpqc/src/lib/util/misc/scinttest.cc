
#include <iostream>

#include <util/misc/scint.h>

using namespace sc;

main()
{
  std::cout << "sizeof(sc_int8_t) = " << sizeof(sc_int8_t) << std::endl;
  std::cout << "sizeof(sc_int16_t) = " << sizeof(sc_int16_t) << std::endl;
  std::cout << "sizeof(sc_int32_t) = " << sizeof(sc_int32_t) << std::endl;
  std::cout << "sizeof(sc_int64_t) = " << sizeof(sc_int64_t) << std::endl;

  return 0;
}
