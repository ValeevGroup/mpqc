
#ifndef _util_misc_formio_h
#define _util_misc_formio_h

#include <iostream.h>

class SCFormIO {
  private:
    static bool ready_;
    static long nindent_;
    static long indent_size_;
    static long skip_indent_;
    static void init();
  public:
    static ios& indent(ios&o);
    static ios& decindent(ios&o);
    static ios& incindent(ios&o);
    static ios& skipnextindent(ios&o);

    static void setindent(ios&o, long column);
    static long getindent(ios&o);
};

ios& indent(ios&);

ios& decindent(ios&);

ios& incindent(ios&);

ios& skipnextindent(ios&);

#endif
