
#include <util/misc/formio.h>
#include <util/group/message.h>

#include <stdio.h> // for vsprintf
#include <string.h>
#include <stdarg.h>

char *SCFormIO::default_basename_ = 0;
int  SCFormIO::ready_ = 0;
long SCFormIO::nindent_ = 0;
long SCFormIO::indent_size_ = 0;
long SCFormIO::skip_indent_ = 0;
int SCFormIO::node_to_print_ = 0;
int SCFormIO::debug_ = 0;
ofstream SCFormIO::nullstream_;
RefMessageGrp SCFormIO::grp_;

char *
SCFormIO::fileext_to_filename(const char *ext)
{
  const char *basename;

  if (default_basename_) basename = default_basename_;
  else basename = "SC";

  char * res = new char[strlen(basename) + strlen(ext) + 1];
  strcpy(res, basename);
  strcat(res, ext);

  return res;
}

void
SCFormIO::set_default_basename(const char *basename)
{
  if (default_basename_) delete[] default_basename_;

  if (basename)
      default_basename_ = strcpy(new char[strlen(basename)+1], basename);
  else
      default_basename_ = 0;
}

const char *
SCFormIO::default_basename()
{
  return default_basename_;
}

void
SCFormIO::set_printnode(int n)
{
  node_to_print_ = n;
}

void
SCFormIO::set_debug(int n)
{
  debug_ = n;
}

void
SCFormIO::set_messagegrp(const RefMessageGrp& g)
{
  grp_ = g;
}
  
void
SCFormIO::init()
{
  ready_ = 1;
  nindent_ = ios::xalloc();
  indent_size_ = ios::xalloc();
  skip_indent_ = ios::xalloc();

  if (grp_.null())
    set_messagegrp(MessageGrp::get_default_messagegrp());

  if (nullstream_.bad() || nullstream_.fail())
    nullstream_.open("/dev/null");
}

ios&
SCFormIO::indent(ios&o)
{
  if (!ready_) init();
  if (debug_) {
      char nn[24];
      sprintf(nn,"node %5d:",grp_->me());
      for (int i=0; i < strlen(nn); i++) o.rdbuf()->sputc(nn[i]);
    }
  long &skip = o.iword(skip_indent_);
  if (skip) {
      skip--;
      return o;
    }
  long n = o.iword(nindent_);
  for (int i=0; i<n; i++) o.rdbuf()->sputc(' ');
  return o;
}

ios&
SCFormIO::incindent(ios&o)
{
  if (!ready_) init();
  long &n = o.iword(nindent_);
  long size = o.iword(indent_size_);
  if (size == 0) size = 2;
  else if (size < 0) size = 0;
  n += size;
  return o;
}

ios&
SCFormIO::decindent(ios&o)
{
  if (!ready_) init();
  long &n = o.iword(nindent_);
  long size = o.iword(indent_size_);
  if (size == 0) size = 2;
  else if (size < 0) size = 0;
  n -= size;
  if (n<0) n=0;
  return o;
}

long
SCFormIO::getindent(ios&o)
{
  if (!ready_) init();
  return o.iword(nindent_);
}

void
SCFormIO::setindent(ios&o, long n)
{
  if (!ready_) init();
  o.iword(nindent_) = n;
}

ios&
SCFormIO::skipnextindent(ios&o)
{
  if (!ready_) init();
  o.iword(skip_indent_)++;
  return o;
}

ostream&
SCFormIO::node0(ostream& o)
{
  if (!ready_) init();
  
  if (!debug_ && node_to_print_ >= 0 && node_to_print_ != grp_->me())
    return nullstream_;

  return o;
}

ios&
indent(ios& o)
{
  return SCFormIO::indent(o);
}

ios&
decindent(ios& o)
{
  return SCFormIO::decindent(o);
}

ios&
incindent(ios& o)
{
  return SCFormIO::incindent(o);
}

ios&
skipnextindent(ios& o)
{
  return SCFormIO::skipnextindent(o);
}

ostream&
node0(ostream& o)
{
  return SCFormIO::node0(o);
}

/////////////////////////////////////////////////////////////////////////////

scprintf::scprintf(const char *fmt, ...)
{
  va_list args;
  
  va_start(args, fmt);

  str[0] = '\0';
  
  // hopefully this won't overflow
  if (fmt && fmt[0]!='\0') {
    if (vsprintf(str, fmt, args) > 1023) {
      cerr << indent << "scprintf overflow\n";
      abort();
    }
  }

  va_end(args);
}

ostream&
operator<<(ostream& o, const scprintf& s)
{
  o << s.str << flush;
  return o;
}
