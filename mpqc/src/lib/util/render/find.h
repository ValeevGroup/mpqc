
#ifndef _util_render_find_h
#define _util_render_find_h

#include <util/render/stack.h>
#include <util/render/parameter.h>

// cannot be used with g++ 2.6-94q4 and has other bugs anyway
template <class T1, class T2>
void
find_parameter_in_stack(Stack<T1>& stack,
                        Parameter<T2>& (T1::*access)(),
                        T2& result
                        )
{
  int have_result = 0;
  for (int i=stack.n()-1; i>=0; i--) {
      if ((stack[i]->*access)().is_set()) {
          if (!have_result || (stack[i]->*access)().overrides()) {
              result = (stack[i]->*access)().value();
              have_result = 1;
            }
        }
    }
}

inline void
find_int_parameter_in_appearance_stack(Stack<RefAppearance>& stack,
                        Parameter<int>& (Appearance::*access)(),
                        int& result
                        )
{
  int have_result = 0;
  for (int i=stack.n()-1; i>=0; i--) {
      if ((stack[i].pointer()->*access)().is_set()) {
          if (!have_result || (stack[i].pointer()->*access)().overrides()) {
              result = (stack[i].pointer()->*access)().value();
              have_result = 1;
            }
        }
    }
}

#endif
