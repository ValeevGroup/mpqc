
#ifdef __GNUG__
#pragma implementation
#endif

#include <mkclasses.h>

MkClasses::MkClasses(const MkClassesOptions &options):
  options_(options)
{
  nerror_ = 0;
  nwarning_ = 0;
  lexer_ = new MkClassesFlexLexer;
}

MkClasses::~MkClasses()
{
  delete lexer_;
}

void
MkClasses::process(istream &input)
{
  lexer_->switch_streams(&input, &cout);
  yparse();

  if (nerror_ || nwarning_) {
      cerr << "MkClasses: " << nerror_ << " error(s) and "
           << nwarning_ << " warning(s)" << endl;
      if (nerror_) {
          exit(1);
        }
    }

  if (options_.geninc()) {
      open(options_.hdrname());
      write_header();
      close();
    }

  if (options_.gensrc()) {
      open(options_.srcname());
      write_source();
      close();
    }

  writedb();
}

void
MkClasses::write_header()
{
  p("#ifndef %s\n", unique_name().c_str());
  p("#define %s\n", unique_name().c_str());

  p("#include <util/class/class.h>");

  vector<string>::iterator i;
  for (i=header_includes_.begin(); i != header_includes_.end();  i++) {
      const char *inc = (*i).c_str();
      if (inc) {
          if (*inc == '<') {
              p("#include %s", inc);
            }
          else {
              p("#include \"%s\"", inc);
            }
        }
    }

  for (vector<Class>::iterator ci = classes_.begin();
       ci != classes_.end();
       ci++) {
      (*ci).write_declaration(p);
    }

  p("#endif");
}

void
MkClasses::write_source()
{
  vector<string>::iterator i;
  for (i=header_includes_.begin(); i != header_includes_.end();  i++) {
      const char *inc = (*i).c_str();
      if (inc) {
          if (*inc == '<') {
              p("#include %s", inc);
            }
          else {
              p("#include \"%s\"", inc);
            }
        }
    }

  for (vector<Class>::iterator ci = classes_.begin();
       ci != classes_.end();
       ci++) {
      (*ci).write_ctor(p);
      (*ci).write_members(p);
      (*ci).write_keyvalin(p);
      (*ci).write_keyvalout(p);
      (*ci).write_statein(p);
      (*ci).write_stateout(p);
    }
}

void
MkClasses::yerror(const char* s)
{
  cerr << "MkClasses: "
       << s << ": line number "
       << (lexer_->lineno() - 1) << endl;
  nerror_++;
  if (nerror_ > 10) {
      cerr << "MkClasses: too many errors" << endl;
      exit(1);
    }
}

void
MkClasses::start_class()
{
  Class c;
  current_class_ = classes_.insert(classes_.end(), c);
}

void
MkClasses::end_class()
{
}

extern "C" {
  int
  MkClasseswrap()
  {
   return 1;
  }
}
