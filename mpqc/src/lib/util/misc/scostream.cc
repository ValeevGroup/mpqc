
#include <string.h>
#ifdef __GNUC__
#pragma implementation
#include <streambuf.h>
#else
#include <fstream.h>
typedef int streamsize;
#endif

#include <stdiostream.h>
#include <util/misc/scostream.h>

// static FILE* debug = fopen("debug.out","w");

// create a standard SCostream output and error stream using file
// descriptors
SCostream SCostream::cout(1);
SCostream SCostream::cerr(2);

#ifdef __GNUC__
#define do_system_write sys_write
#else
#define do_system_write xsputn
#endif

// see if a sensible filebuf can be created
class SCfilebuf: public filebuf {
  private:
    int _column;
  protected:
    streamsize do_system_write(const char*, streamsize);
  public:
    SCfilebuf(int);
    SCfilebuf(filebuf*);
    ~SCfilebuf();
    int get_column();
};

SCfilebuf::SCfilebuf(int fd):
  filebuf(fd),
  _column(0)
{
}

SCfilebuf::SCfilebuf(filebuf * fb):
#if defined(SGI) && !defined(__GNUC__)
  filebuf(fb->fd()),
#else  
  filebuf(*fb),
#endif  
  _column(0)
{
}

SCfilebuf::~SCfilebuf()
{
}

streamsize
SCfilebuf::do_system_write(const char* s, streamsize n)
{
//   fprintf(debug, "sys_write called n = %d\n",n);
  streamsize ret = filebuf::do_system_write(s, n);

  // find the last newline in s
  int i;
  for (i=n-1; i>=0; i--) {
      if (s[i] == '\n') break;
    }
  if (i == -1) _column += n;
  else _column = n - i - 1;
//   fprintf(stderr,"i = %d, n = %d, _column = %d\n",i,n,_column);
//   fprintf(stderr,"\"");
//   for (i=0; i<n; i++) fprintf(stderr,"%c",s[i]);
//   fprintf(stderr,"\"\n");

  return ret;
}

int
SCfilebuf::get_column()
{
//   fprintf(debug,"get_column got %d\n",_column);
  // flush updates _column:
  return _column;
}

SCostream::SCostream():
  nskip(0),
  indentation(0),
  indentation_increment(2),
  indentation_maximum(20)
{
}

SCostream::SCostream(FILE*fp):
  nskip(0),
  indentation(0),
  indentation_increment(2),
  indentation_maximum(20)
{
  init(new stdiobuf(fp));
}

SCostream::SCostream(int fd):
  nskip(0),
  indentation(0),
  indentation_increment(2),
  indentation_maximum(20)
{
#ifdef linux
  init(new filebuf(fd));
#else
  init(new SCfilebuf(fd));
#endif
}

// this causes very strange problems
// SCostream::SCostream(ostream&o):
//   ostream(o.rdbuf()), // the streambuf is not destroyed if the ctor is used
//   indentation(0),
//   indentation_increment(2),
//   indentation_maximum(20),
//   nskip(0)
// {
//   // the streambuf is destroyed if init is used.
//   //init(o.rdbuf());
// }

SCostream::~SCostream()
{
}

void
SCostream::skip_next_indent()
{
  nskip++;
}

void
SCostream::set_indent_to_column()
{
  int tmp = get_column();
  if (tmp > 0) indentation += tmp;
}

int
SCostream::get_column()
{
  // flush updates SCfilebuf's _column
  flush();
#ifdef __GNUC__
  return rdbuf()->get_column();
#else
  return 0;
#endif
}

void
SCostream::set_indent(int i)
{
  indentation = i;
}

int
SCostream::get_indent()
{
  return indentation;
}

void
SCostream::set_indent_inc(int i)
{
  indentation_increment = i;
}

int
SCostream::get_indent_inc()
{
  return indentation_increment;
}

ostream&
SCostream::indent()
{
  return indent(indentation);
}

ostream&
SCostream::indent(int ni)
{
  if (nskip) {
      nskip--;
    }
  else {
      for (int i=0; i<ni; i++) {
          (*this) << ' ';
          //put(' ');
        }
    }
  return *this;
}

SCostream&
SCostream::operator--()
{
  indentation -= indentation_increment;
  return *this;
}

SCostream&
SCostream::operator++()
{
  indentation += indentation_increment;
  return *this;
}

SCostream&
SCostream::operator--(int)
{
  indentation -= indentation_increment;
  return *this;
}

SCostream&
SCostream::operator++(int)
{
  indentation += indentation_increment;
  return *this;
}
