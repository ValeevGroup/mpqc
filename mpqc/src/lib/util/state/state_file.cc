
#include <stdio.h>
#include <util/class/class.h>
#include "state.h"

#include "statenumSet.h"
#include "classdImplMap.h"

StateOutFile::StateOutFile() :
  fp_(stdout), opened_(0)
{
}

StateOutFile::StateOutFile(FILE* fp) :
  fp_(fp), opened_(0)
{
}

StateOutFile::StateOutFile(const char * path, const char * mode) :
  opened_(1)
{
  fp_ = fopen(path,mode);
}

StateOutFile::~StateOutFile()
{
  if (opened_) close();
}

void StateOutFile::flush() { fflush(fp_); }
void StateOutFile::close()
{
  if(opened_) fclose(fp_);
  opened_=0; fp_=0;

  _classidmap->clear(); _nextclassid=0;
  forget();
}

void StateOutFile::rewind() { if(fp_) fseek(fp_,0,0); }

int StateOutFile::open(const char *path, const char * mode)
{
  if (opened_) close();

  if ((fp_ = fopen(path,mode))==0) {
    fprintf(stderr,"StateOutFile::open(%s,%s) failed\n",path,mode);
    return -1;
    }

  opened_ = 1;
  return 0;
}

////////////////////////////////////

StateInFile::StateInFile() :
  fp_(stdin), opened_(0)
{
}

StateInFile::StateInFile(FILE* fp) :
  fp_(fp), opened_(0)
{
}

StateInFile::StateInFile(const char * path, const char * mode) :
  opened_(1)
{
  fp_ = fopen(path,mode);
}

StateInFile::~StateInFile()
{
  if (opened_) close();
}

void StateInFile::flush() { fflush(fp_); }
void StateInFile::close()
{
  if(opened_) fclose(fp_);
  opened_=0; fp_=0;

  _cd.clear();
  forget();
}
void StateInFile::rewind() { if(fp_) fseek(fp_,0,0); }

int StateInFile::open(const char *path, const char * mode)
{
  if (opened_) close();

  if ((fp_ = fopen(path,mode))==0) {
    fprintf(stderr,"StateInFile::open(%s,%s) failed\n",path,mode);
    return -1;
    }

  opened_ = 1;
  return 0;
}

