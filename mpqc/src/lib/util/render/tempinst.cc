
#include <util/render/material.h>
#include <util/render/appearance.h>
#include <util/render/transform.h>
#include <util/render/stack.h>
#include <util/render/parameter.h>

#ifdef __GNUC__
template class Stack<DCRefMaterial>;
template class Stack<DCRefAppearance>;
template class Stack<DCRefTransform>;
template class Parameter<Color>;
template class Parameter<int>;
#endif
