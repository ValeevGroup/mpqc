
#ifndef _util_misc_formio_h
#define _util_misc_formio_h

#include <iostream.h>
#include <fstream.h>

#include <util/class/class.h>

DescribedClass_REF_fwddec(MessageGrp);

class SCFormIO {
  private:
    static char *default_basename_;
    static int  ready_;
    static long nindent_;
    static long indent_size_;
    static long skip_indent_;
    static int node_to_print_;
    static int debug_;
    static ofstream nullstream_;
    static RefMessageGrp grp_;
    static void init();
  public:
    static ios& indent(ios&o);
    static ios& decindent(ios&o);
    static ios& incindent(ios&o);
    static ios& skipnextindent(ios&o);
    static ostream& node0(ostream&o);

    static void setindent(ios&o, long column);
    static long getindent(ios&o);
    static void set_printnode(int);
    static void set_debug(int);
    static void set_messagegrp(const RefMessageGrp& grp_);

    static void set_default_basename(const char *);
    static const char *default_basename();
    static char *fileext_to_filename(const char *extension);
};

ios& indent(ios&);

ios& decindent(ios&);

ios& incindent(ios&);

ios& skipnextindent(ios&);

ostream& node0(ostream&);

/////////////////////////////////////////////////////////////////////////////

class scprintf {
  private:
    char str[1024];

  public:
    scprintf(const char*,...);
    friend ostream& operator<<(ostream&, const scprintf&);
};

ostream& operator<<(ostream&, const scprintf&);

#endif
