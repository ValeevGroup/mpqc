
#ifndef _util_misc_scostream_h
#define _util_misc_scostream_h

#include <stdio.h>
#include <string.h>
#include <iostream.h>

class SCostream: public ostream {
  private:
    int nskip;
    int indentation;
    int indentation_increment;
    int indentation_maximum;
  public:
    static SCostream cout;
    static SCostream cerr;

    SCostream();
    SCostream(FILE*);
    SCostream(int fd);
    //SCostream(ostream&);
    ~SCostream();

    void set_indent_to_column();
    int get_column();
    void skip_next_indent();
    void set_indent(int);
    int get_indent();
    void set_indent_inc(int);
    int get_indent_inc();
    ostream& indent();
    ostream& indent(int);
    SCostream& operator--();
    SCostream& operator++();
    SCostream& operator--(int);
    SCostream& operator++(int);
};

#endif
