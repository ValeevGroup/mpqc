
#ifndef _util_render_stack_h
#define _util_render_stack_h

#include <stdio.h>

template <class T>
class Stack {
  private:
    const int max_stack_size = 20;
    T objects[max_stack_size];
    int nobjects;
  public:
    Stack(): nobjects(0) {}
    void push(const T&a) {
        if (nobjects >= max_stack_size) {
            fprintf(stderr,"Stack: overflow\n");
            abort();
          }
        objects[nobjects++] = a;
      }
    T pop() {
        if (!nobjects) {
            fprintf(stderr,"Stack: underflow\n");
            abort();
          }
        nobjects -= 1;
        return objects[nobjects];
      }
    T top() const {
        if (!nobjects) {
            fprintf(stderr,"Stack: underflow\n");
            abort();
          }
        return objects[nobjects - 1];
      }
    int n() const { return nobjects; }
    T operator[](int i) { return objects[i]; }
    void print(FILE*fp = stdout) {
        fprintf(fp, "Stack (depth = %d):\n", nobjects);
        for (int i=0; i<nobjects; i++) {
            fprintf(fp, "  object %d:\n", i);
            objects[i]->print(fp);
          }
      }
};

#endif
