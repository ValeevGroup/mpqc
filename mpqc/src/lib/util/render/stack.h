//
// stack.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _util_render_stack_h
#define _util_render_stack_h

#include <iostream>

namespace sc {

#define STACK_MAX_STACK_SIZE 20
template <class T>
class Stack {
  private:
    T objects[STACK_MAX_STACK_SIZE];
    int nobjects;
  public:
    Stack(): nobjects(0) {}
    void push(const T&a) {
        if (nobjects >= STACK_MAX_STACK_SIZE) {
            ExEnv::errn() << "Stack: overflow" << std::endl;
            abort();
          }
        objects[nobjects++] = a;
      }
    T pop() {
        if (!nobjects) {
            ExEnv::errn() << "Stack: underflow" << std::endl;
            abort();
          }
        nobjects -= 1;
        return objects[nobjects];
      }
    T top() const {
        if (!nobjects) {
            ExEnv::errn() << "Stack: underflow" << std::endl;
            abort();
          }
        return objects[nobjects - 1];
      }
    int n() const { return nobjects; }
    T operator[](int i) { return objects[i]; }
    void print(std::ostream& os = ExEnv::out0()) {
        os << "Stack (depth = " << nobjects << "):" << std::endl;
        for (int i=0; i<nobjects; i++) {
            os << "  object " << i << ":" << std::endl;
            objects[i]->print(os);
          }
      }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
