
#include <stdio.h>
#include <util/class/class.h>
#include <util/state/state.h>

StateOutBin::StateOutBin() :
  StateOutFile()
{
}

StateOutBin::StateOutBin(FILE* fp) :
  StateOutFile(fp)
{
}

StateOutBin::StateOutBin(const char *path, const char * mode) :
  StateOutFile(path,mode)
{
}

StateOutBin::~StateOutBin()
{
}

StateInBin::StateInBin() :
  StateInFile()
{
}

StateInBin::StateInBin(FILE* fp) :
  StateInFile(fp)
{
}

StateInBin::StateInBin(const char *path, const char * mode) :
  StateInFile(path,mode)
{
}

StateInBin::~StateInBin()
{
}

////////////////////////////////////////////////////////////////

int StateOutBin::put_array_void(const void*p,int size)
{
  return fwrite(p,1,size,fp_);
}

int StateInBin::get_array_void(void*p,int size)
{
  return fread(p,1,size,fp_);
}
