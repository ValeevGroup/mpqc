
#include <stdio.h>
#include <stdarg.h>

#include "mkclasses.h"

string
to_outfilename(const string &infilenamestr, const string &suffixstr)
{
  const char *infilename = infilenamestr.c_str();
  const char *suffix = suffixstr.c_str();
  char *slash = strrchr((char*)infilename, '/');
  if (slash) infilename = &slash[1];
  char *outfilename = new char[strlen(infilename) + strlen(suffix) + 1];
  strcpy(outfilename, infilename);
  char *dot = strrchr(outfilename, '.');
  if (dot) *dot = '\0';
  strcat(outfilename, suffix);
  string result(outfilename);
  delete[] outfilename;
  return result;
}

void
Out::open(const string &filename)
{
  close();
  indent_ = 0;
  out_ = fopen(filename.c_str(), "w");
  if (!out_) {
      fprintf(stderr, "MkClasses::open(\"%s\"): failed\n", filename.c_str());
      exit(1);
    }
}

void
Out::close()
{
  if (out_) {
      fclose(out_);
      out_ = 0;
    }
}

void
Out::operator ()(const char *format, ...)
{
  int i;
  for (i=0; i<indent_; i++) {
      fprintf(out_, "  ");
    }

  va_list args;
  va_start(args, format);
  vfprintf(out_, format, args);
  fprintf(out_, "\n");
  va_end(args);
}

void
Out::operator()(bool indent, bool newline, const char* format, ...)
{
  if (indent) {
      for (int i = -1; i<indent_; i++) {
          fprintf(out_, "  ");
        }
    }

  va_list args;
  va_start(args, format);
  vfprintf(out_, format, args);
  if (newline) fprintf(out_, "\n");
  va_end(args);
}
