
#include <util/misc/formio.h>

const char *SCFormIO::default_basename_ = 0;
int  SCFormIO::ready_ = 0;
long SCFormIO::nindent_ = 0;
long SCFormIO::indent_size_ = 0;
long SCFormIO::skip_indent_ = 0;

void
SCFormIO::set_default_basename(const char *basename)
{
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
SCFormIO::init()
{
  ready_ = 1;
  nindent_ = ios::xalloc();
  indent_size_ = ios::xalloc();
  skip_indent_ = ios::xalloc();
}

ios&
SCFormIO::indent(ios&o)
{
  if (!ready_) init();
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

