
#include "color.h"

Color::Color(const RefKeyVal& keyval)
{
  const char* rgb = "rgb";
  red_ = keyval->doublevalue(rgb,0);
  green_ = keyval->doublevalue(rgb,1);
  blue_ = keyval->doublevalue(rgb,2);
}
