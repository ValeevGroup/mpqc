
#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mkclasses_h
#define _mkclasses_h

#include <stdio.h>
#include <string>
#include <vector>

#ifndef in_mkclasses_scanner
#undef yyFlexLexer
#define yyFlexLexer MkClassesFlexLexer
#include <FlexLexer.h>
#undef yyFlexLexer
#endif

string to_outfilename(const string &infilename, const string &suffix);

class Out {
  private:
    FILE *out_;
    int indent_;
  public:
    // output file manipulation members
    void open(const string &filename);
    void close();
    void operator()(const char *fmt, ...);
    //void operator()(void) { operator(""); }
    void operator()(bool indent, bool newline, const char *fmt, ...);
    int &depth() { return indent_; }
    void operator ++() { depth()++; }
    void operator --() { depth()--; }
    void operator ++(int) { depth()++; }
    void operator --(int) { depth()--; }
};

class Parent {
  public:
    enum Access { Private, Protected, Public };
  private:
    Access access_;
    bool virtual_;
    string name_;
  public:
    Parent();
    Parent(const string &);
    ~Parent();

    Access access() const { return access_; }
    bool is_virtual() const { return virtual_; }
    const string &name() const { return name_; }

    string stringrep() const;
};

class DataMember {
  private:
    string name_;
    string type_;
    string long_name_;
    string description_;
  public:
    DataMember();
    ~DataMember();

    void set_name(const string &s) { name_ = s; }
    void set_type(const string &s) { type_ = s; }
    void set_long_name(const string &s) { long_name_ = s; }
    void set_description(const string &s) { description_ = s; }

    const string &name() const { return name_; }
    const string &type() const { return type_; }
    const string &long_name() const { return long_name_; }
    const string &description() const { return description_; }
};

class Class {
  private:
    string name_;
    vector<Parent> parents_;
    vector<DataMember> member_data_;
    int version_;
    string description_;
    bool castdown_;
    bool require_castdown_;
    bool static_class_desc_;
    bool ctor_void_;
    bool ctor_keyval_;
    bool ctor_statein_;

    vector<DataMember>::iterator current_member_datum_;
  public:
    Class();
    ~Class();

    void set_name(const string &n) { name_ = n; }
    void set_parents(const vector<string> &p);
    void set_version(int v) { version_ = v; }
    void set_description(const string &s) { description_ = s; }
    void set_castdown(bool b) { castdown_ = b; }
    void set_require_castdown(bool b) { require_castdown_ = b; }
    void set_static_class_desc(bool b) { static_class_desc_ = b; }
    void set_ctor_void(bool b) { ctor_void_ = b; }
    void set_ctor_keyval(bool b) { ctor_keyval_ = b; }
    void set_ctor_statein(bool b) { ctor_statein_ = b; }

    const string &name() const { return name_; }
    const vector<Parent> &parents() const { return parents_; }
    int version() const { return version_; }
    const string &description() const { return description_; }
    bool castdown() const { return castdown_; }
    bool require_castdown() const { return require_castdown_; }
    bool static_class_desc() const { return static_class_desc_; }
    bool ctor_void() const { return ctor_void_; }
    bool ctor_keyval() const { return ctor_keyval_; }
    bool ctor_statein() const { return ctor_statein_; }
    bool void_ctor() const { return ctor_void_; }
    bool keyval_ctor() const { return ctor_keyval_; }
    bool statein_ctor() const { return ctor_statein_; }

    void write_declaration(Out &);
    void write_definition(Out &);
    void write_member_declaration(Out &, const DataMember &);

    void write_ctor(Out &);
    void write_members(Out &);

    void write_castdown(Out &);

    void write_keyvalin(Out &);
    void write_keyvalout(Out &);

    void write_statein(Out &);
    void write_stateout(Out &);

    string parentstring() const;

    const vector<DataMember> &member_data() { return member_data_; }

    DataMember &current_member_datum() { return *current_member_datum_; }
    void start_member_datum();
    void end_member_datum();
};

class MkClassesOptions {
  private:
    string libname_;
    string dbname_;

    string hdrname_;
    string srcname_;

    bool gensrc_;
    bool geninc_;
  public:
    MkClassesOptions();
    ~MkClassesOptions();

    void set_libname(const string & s) { libname_ = s; }
    void set_dbname(const string & s) { dbname_ = s; }
    void set_hdrname(const string & s) { hdrname_ = s; }
    void set_srcname(const string & s) { srcname_ = s; }
    void set_gensrc(bool b) { gensrc_ = b; }
    void set_geninc(bool b) { geninc_ = b; }

    const string &libname() const { return libname_; }
    const string &dbname() const { return dbname_; }
    const string &hdrname() const { return hdrname_; }
    const string &srcname() const { return srcname_; }
    bool gensrc() const { return gensrc_; }
    bool geninc() const { return geninc_; }
};

class MkClasses {
  private:
    string unique_name_;
    vector<string> header_includes_;
    vector<string> source_includes_;

    vector<Class> classes_;

    MkClassesFlexLexer *lexer_;

    MkClassesOptions options_;

    Out p;
    void open(const string &filename) { p.open(filename); }
    void close() { p.close(); }
    int &depth() { return p.depth(); }

    int ylex() { return lexer_->yylex(); }
    int yparse();
    void yerror(const char* s);

    void write_source();
    void write_header();

    vector<Class>::iterator current_class_;

    int nerror_;
    int nwarning_;
  public:
    MkClasses(const MkClassesOptions &options);
    ~MkClasses();

    void set_unique_name(const string &s) { unique_name_ = s; }
    void set_header_inc(const vector<string> &v) { header_includes_ = v; }
    void set_source_inc(const vector<string> &v) { source_includes_ = v; }
    void start_class();
    void end_class();
    Class& current_class() { return *current_class_; }

    const string &unique_name() { return unique_name_; }

    void writedb();
    void process(istream&);
};

#endif
